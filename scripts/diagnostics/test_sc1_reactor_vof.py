"""Reference isolation, observer-only source changes, budgets and failure gates."""
from contextlib import redirect_stdout
import difflib
import hashlib
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import sc1_reactor_vof as ref
import sc1_reactor_vof_build as build
import sc1_reactor_vof_case as cases
import sc1_reactor_vof_report as report
import sc1_reactor_vof_instrument as instrument
import sc1_reactor_watch as watch
import sc1_reactor_failure as failure
from sc1_reactor_startup_case import field

INITIAL={'liquid_volume_m3':1.,'gas_volume_m3':.642,'interface_midheight_m':1.0186}


def good_log():
    lines=['SC1_VOF_INITIAL liquid_m3=1 gas_m3=.642 version=1']
    for i in range(1,201):
        t=i*report.DT;lines.append(f'Time = {t:.12g}')
        lines.extend([f'SC1_VOF_PRESSURE t_s={t:.17g} min_Pa=99983 max_Pa=111800 nonfinite=0']*4)
        # Known throughflow: 1e-5 m3/s liquid leaves, gas replaces it.
        vl=1-t*1e-5;vg=.642+t*1e-5
        lines.append(f'SC1_VOF_STEP t_s={t:.17g} dt_s=5e-5 liquid_m3={vl:.17g} gas_m3={vg:.17g} '
            'liquidOut_m3_s=1e-5 gasOut_m3_s=-1e-5 postPressureOut_m3_s=0 '
            'liquidResidual_kg=0 gasResidual_kg=0 liquidL1_kg=0 gasL1_kg=0 '
            'alphaMin=0 alphaMax=1 complement=0 liquidMax_m_s=.0001 gasMax_m_s=.0001 '
            'alphaChangeL1_m3=1e-7 Co=.01 alphaCalls=2')
    return '\n'.join(lines)+'\nEnd\n'


class Budgets(unittest.TestCase):
    def test_known_export_and_replacement_close_both_phases(self):
        r,_=report.budget(good_log(),INITIAL)
        self.assertFalse(r['errors']);self.assertAlmostEqual(r['phases']['liquid']['boundaryOut_kg'],.0001)
        self.assertLess(abs(r['phases']['liquid']['physicalResidual_kg']),1e-12)
        self.assertLess(abs(r['phases']['gas']['physicalResidual_kg']),1e-15)

    def test_incomplete_corrupt_and_nonfinite_budgets_rejected(self):
        for changed in (good_log().replace('End\n',''),good_log().replace('nonfinite=0','nonfinite=1',1),
            good_log().replace('dt_s=5e-5','dt_s=nan',1),good_log().replace('alphaCalls=2','alphaCalls=4',1),
            good_log().replace('liquidOut_m3_s=1e-5','liquidOut_m3_s=2e-5',1),
            good_log().replace('liquidL1_kg=0','liquidL1_kg=.001',1)):
            with self.assertRaises(ValueError):report.budget(changed,INITIAL)

    def test_quiescent_failure_cannot_be_accepted(self):
        r,_=report.budget(good_log().replace('liquidMax_m_s=.0001','liquidMax_m_s=.002',1),INITIAL)
        self.assertIn('VOF liquid speed exceeded 1 mm/s',r['errors'])

    def test_native_snapshots_match_logged_inventory_and_shared_velocity(self):
        with tempfile.TemporaryDirectory() as d:
            case=Path(d);(case/'log.solver').write_text(good_log())
            volumes=np.array([1.,.642]);centres=np.array([[0,0,0],[0,1.5,0]])
            def read(p,n,vector=False):
                t=float(p.parent.name)
                if p.name=='C':return centres
                if p.name=='V':return volumes
                if p.name=='alpha.liquid':return np.array([1-t*1e-5,0.])
                if p.name=='U':return np.tile([.0001,0,0],(2,1))
                if p.name=='p_rgh':return np.ones(2)*110000
                raise AssertionError(p)
            masks={'walls':np.array([True,True]),'propeller':np.array([False,False]),'outlet':np.array([False,True])}
            with patch.object(report,'internal',side_effect=read),patch.object(report,'patch_cells',return_value=masks):
                r=report.analyze(case,INITIAL,2)
                self.assertTrue(r['reference_pass']);self.assertEqual(len(r['snapshots']),40)
                self.assertEqual(r['snapshots'][-1]['phases']['liquid']['present']['max_m_s'],.0001)


class SourceAndCase(unittest.TestCase):
    def test_exact_upstream_copy_is_observer_only_and_private(self):
        fixture=json.loads(build.FIXTURE.read_text())['sources']
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);project=root/'OpenFOAM-9'
            for name,text in fixture.items():
                p=project/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text)
            before=build.native_sources(project);out=root/'build';private=build.private_sources(project,out)
            self.assertEqual(before,build.native_sources(project))
            self.assertEqual(private['alphaEqn.H'],fixture['src/twoPhaseModels/twoPhaseMixture/VoF/alphaEqn.H'])
            original=fixture['applications/solvers/multiphase/interFoam/pEqn.H']
            self.assertEqual(private['pEqn.H'].replace('\n        sc1Audit.pressure(p_rgh,rho,gh);',''),original)
            main=private['interFoam.C']
            main_original=fixture['applications/solvers/multiphase/interFoam/interFoam.C']
            self.assertFalse([line for line in difflib.ndiff(main_original.splitlines(),main.splitlines()) if line.startswith('- ')])
            self.assertLess(main.index('sc1Audit.step('),main.index('runTime.write();'))
            self.assertLess(main.index('sc1Audit.alphaFlux('),main.index('#include "UEqn.H"'))
            self.assertIn(str(out/'bin'/build.APP),private['Make/files'])
            self.assertNotIn('FOAM_APPBIN',private['Make/files'])
            (project/'applications/solvers/multiphase/interFoam/pEqn.H').write_text('changed')
            with self.assertRaisesRegex(ValueError,'differs'):build.native_sources(project)
        with self.assertRaises(ValueError):instrument.pressure(original+'\n        p_rghEqn.solve();')

    def test_private_build_verification_rejects_changed_sources(self):
        fixture=json.loads(build.FIXTURE.read_text())['sources']
        with tempfile.TemporaryDirectory() as d,redirect_stdout(io.StringIO()):
            root=Path(d);project=root/'OpenFOAM-9'
            for name,text in fixture.items():
                p=project/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text)
            out=root/'build';identity={'binary':str(out/'bin'/build.APP),'sha256':'binary','assets':{}}
            before=build.native_sources(project)
            with (patch.dict(os.environ,WM_PROJECT_VERSION='9',WM_PROJECT_DIR=str(project)),
                  patch.object(build.startup,'require_allocation'),patch.object(build.shutil,'which',return_value='wmake'),
                  patch.object(build.base,'execute') as execute,patch.object(build,'identity',return_value=identity)):
                self.assertEqual(build.build(out),out);self.assertEqual(build.verify(out),identity)
                execute.assert_called_once();self.assertEqual(execute.call_args.args[0],['wmake'])
                self.assertEqual(build.native_sources(project),before)
                (out/'src/pEqn.H').write_text('modified')
                with self.assertRaisesRegex(ValueError,'source identity changed'):build.verify(out)

    def test_boundaries_properties_and_alpha_contract(self):
        files=cases.case_files();u=files['0/U']
        self.assertIn('object U;',u);self.assertNotIn('object U.liquid;',u)
        self.assertIn('outlet { type pressureInletOutletVelocity;',u)
        self.assertIn('propeller { type noSlip;',u);self.assertIn('walls { type noSlip;',u)
        self.assertIn('inside_to_hub { type empty;',u)
        self.assertIn('nu 1.8e-5;',files['constant/transportProperties'])
        self.assertIn('sigma 0.07;',files['constant/transportProperties'])
        self.assertIn('nAlphaSubCycles 1;',files['system/fvSolution'])
        self.assertIn('nOuterCorrectors 2;',files['system/fvSolution'])
        self.assertNotIn('limitPressure',files['system/fvConstraints'])


class Runner(unittest.TestCase):
    def parent(self,p):
        old=p/'control';case=old/'quiescent'
        for name in ref.COPIED+('constant/polyMesh/points','constant/polyMesh/faces'):
            q=case/name;q.parent.mkdir(parents=True,exist_ok=True);q.write_text('preserved '+name)
        for name in ('log.checkMesh','log.writeCellCentres','log.writeCellVolumes'):(case/name).write_text('preserved')
        m=dict(root=str(old),inputs=ref.base.case_input_hashes(case),initial=INITIAL,ranks=8,
            archive=str(p/'archive'),regression_root=str(p/'rates'),build=str(p/'oldbuild'),mesh_hashes={'points':'p'})
        for key in ('baseline','failed_control','face_control','coupled_control','pressure_probe','velocity_control'):
            m[key]={'root':str(p/key)}
        buildroot=p/'build';buildroot.mkdir();(buildroot/'vof_build_manifest.json').write_text('{}')
        identity=dict(root=str(old),manifest_sha256='m',log_sha256=ref.CONTROL_LOG,result_sha256=ref.CONTROL_RESULT)
        binary=dict(binary=str(buildroot/'bin'/build.APP),assets={'test':'hash'},sha256='hash',dependencies={})
        return old,m,identity,buildroot,binary

    def test_fresh_control_preserved_and_tampering_rejected(self):
        with tempfile.TemporaryDirectory() as d,redirect_stdout(io.StringIO()):
            p=Path(d);old,parent,identity,buildroot,binary=self.parent(p);before=ref.startup.tree_hashes(old)
            with patch.object(ref,'control',return_value=(parent,identity)),patch.object(build,'verify',return_value=binary),patch.object(ref.startup,'require_allocation'):
                root=ref.prepare(old,buildroot,p/'new');m=ref.verify(root)
                self.assertEqual(before,ref.startup.tree_hashes(old))
                self.assertEqual(m['inputs']['0/alpha.liquid'],parent['inputs']['0/alpha.liquid'])
                self.assertNotIn('constant/phaseProperties',m['inputs']);self.assertNotIn('0/U.gas',m['inputs'])
                (root/'quiescent/attempt.json').write_text('attempt')
                with self.assertRaisesRegex(ValueError,'attempted'):ref.run(root)
                (root/'quiescent/0/p_rgh').write_text('tampered')
                with self.assertRaisesRegex(ValueError,'inputs changed'):ref.verify(root)

    def test_failed_solver_collected_without_reconstruction(self):
        with tempfile.TemporaryDirectory() as d,redirect_stdout(io.StringIO()):
            p=Path(d);old,parent,identity,buildroot,binary=self.parent(p)
            with patch.object(ref,'control',return_value=(parent,identity)),patch.object(build,'verify',return_value=binary),patch.object(ref.startup,'require_allocation'):
                root=ref.prepare(old,buildroot,p/'new')
                def solver(command,case,log,env):log.write_text('FOAM FATAL ERROR\n');raise ValueError('failed')
                with patch.object(watch,'execute',side_effect=solver),patch.object(ref.base,'execute') as utility,patch.object(failure,'collect') as collect:
                    with self.assertRaisesRegex(ValueError,'failed'):ref.run(root)
                    self.assertEqual(utility.call_count,1);collect.assert_called_once()

    def test_success_reconstructs_only_complete_snapshots(self):
        with tempfile.TemporaryDirectory() as d,redirect_stdout(io.StringIO()):
            p=Path(d);old,parent,identity,buildroot,binary=self.parent(p)
            with patch.object(ref,'control',return_value=(parent,identity)),patch.object(build,'verify',return_value=binary),patch.object(ref.startup,'require_allocation'):
                root=ref.prepare(old,buildroot,p/'new')
                def solver(command,case,log,env):
                    log.write_text(good_log());self.assertIn('--time=00:30:00',command)
                    for rank in range(8):
                        for i in range(5,201,5):
                            q=case/f'processor{rank}'/f'{i*report.DT:.12g}';q.mkdir(parents=True)
                            for f in ('U','alpha.liquid','p_rgh'):(q/f).write_text('output')
                result={'reference_pass':True,'snapshots':[{'time_s':.01}]}
                with patch.object(watch,'execute',side_effect=solver),patch.object(ref.base,'execute') as utility,patch.object(report,'analyze',return_value=result),patch.object(failure,'collect') as collect:
                    self.assertTrue(ref.run(root)['reference_pass'])
                    self.assertEqual([x.args[0][0] for x in utility.call_args_list],['decomposePar','reconstructPar'])
                    collect.assert_called_once()


if __name__=='__main__':unittest.main()
