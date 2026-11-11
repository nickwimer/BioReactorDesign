"""Hydrostatic identity, active buoyancy, immutable preparation and acceptance gates."""
from contextlib import ExitStack,redirect_stdout
import io
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import sc1_reactor_vof_hydro as candidate
import sc1_reactor_vof_hydro_build as build
import sc1_reactor_vof_hydro_instrument as instrument
import sc1_reactor_vof_hydro_report as report
import sc1_reactor_watch as watch
import test_sc1_reactor_vof_force as force_tests
import test_sc1_reactor_vof_nonorth as old

HEADER='SC1_HYDROSTATIC_REFERENCE version=1 frozen=1 rhoError_kg_m3=0 pressureError_Pa=0 initialSpeed_m_s=0\n'
def log():return HEADER+old.corrected_log().replace('Solving for p_rgh,','Solving for pPert,')


class Physics(unittest.TestCase):
    def test_compiled_profile_and_full_density_buoyancy_identity(self):
        # Compile the actual independent helper; this is not an OpenFOAM build.
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);cpp=p/'profile.C';exe=p/'profile'
            cpp.write_text('#include <iostream>\n#include <iomanip>\n#include "sc1HydrostaticProfile.H"\nint main(){double y;while(std::cin>>y) std::cout<<std::setprecision(17)<<sc1Hydro::alpha(y,1.0185981622921294,.0125)<<" "<<sc1Hydro::pressure(y,1.0185981622921294,.0125)<<"\\n";}\n')
            subprocess.run(['g++','-std=c++14','-Wall','-Wextra','-Werror','-I',str(build.base.HERE/'include'),str(cpp),'-o',str(exe)],check=True,capture_output=True)
            y=np.array([-.2,.5,1.0125,1.015,1.0185981622921294,1.02,1.0247,1.04,1.74]);eps=1e-5
            data=np.concatenate((y,y-eps,y+eps))
            values=np.loadtxt(io.StringIO(subprocess.check_output([str(exe)],input='\n'.join(map(str,data)),text=True)))
            alpha,ph=values[:len(y)].T;gradient=(values[2*len(y):,1]-values[len(y):2*len(y),1])/(2*eps)
            rho0=1+999*alpha
            np.testing.assert_allclose(gradient,-9.81*rho0,rtol=2e-7,atol=1e-6)
            delta=np.linspace(-.4,.7,len(y));qgrad=np.linspace(-20,30,len(y));rho=rho0+delta
            np.testing.assert_allclose(-gradient-qgrad-9.81*rho,-qgrad-9.81*delta,rtol=1e-6,atol=2e-6)
            self.assertTrue(np.any(np.abs(9.81*delta)>1))

    def test_rotated_face_projection_keeps_nonzero_buoyancy(self):
        normals=np.array([[.8,.6,0],[-.6,.8,0],[0,0,1.]])
        normals=np.concatenate((normals,-normals));areas=np.array([1,2,3,1,2,3.])
        gravity=np.array([0,-9.81,0]);delta=7.5;force=delta*(normals@gravity)
        recovered=np.linalg.solve(np.einsum('i,ij,ik->jk',areas,normals,normals),np.einsum('i,i,ij->j',areas,force,normals))
        np.testing.assert_allclose(recovered,delta*gravity,atol=1e-13)
        self.assertNotEqual(np.linalg.norm(recovered),0)

    def test_initial_fill_rejects_nonhydrostatic_pressure(self):
        f=candidate.fields;c=np.array([[0,.9,0],[0,1.2,0.]]);v=np.array([1.,.642])
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);files,initial=f.initialized_fields(c,v)
            for rel,text in files.items():
                p=root/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text)
            original=f.internal;perturb=0
            def read(p,n,vector=False):
                if p.name in ('C','V','U'):return {'C':c,'V':v,'U':np.zeros_like(c)}[p.name]
                value=original(p,n,vector)
                return value+perturb if p.name=='p_rgh' else value
            with patch.object(f,'internal',side_effect=read):
                candidate.check_initial(root,initial);perturb=1
                with self.assertRaises(ValueError):candidate.check_initial(root,initial)


class Build(unittest.TestCase):
    def test_private_build_keeps_transport_audit_inertia_and_surface_tension(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);project=p/'OpenFOAM-9';base=p/'base'
            for rel,text in json.loads(build.reference.FIXTURE.read_text())['sources'].items():
                q=project/rel;q.parent.mkdir(parents=True,exist_ok=True);q.write_text(text)
            original=build.reference.private_sources(project,base)
            for rel,text in original.items():
                q=base/'src'/rel;q.parent.mkdir(parents=True,exist_ok=True);q.write_text(text)
            source=build.private_sources(base,p/'new')
            self.assertEqual({k for k in original if source[k]!=original[k]}, {'interFoam.C','pEqn.H','UEqn.H','Make/files'})
            for k in ('alphaEqn.H','alphaEqnSubCycle.H','sc1VofAudit.H','createFields.H','Make/options'):self.assertEqual(source[k],original[k])
            for k in ('pEqn.H','UEqn.H'):
                self.assertIn(instrument.BUOYANCY,source[k]);self.assertIn('mixture.surfaceTensionForce()',source[k]);self.assertNotIn('ghf*fvc::snGrad(rho)',source[k])
            self.assertIn('fvm::ddt(rho, U)',source['UEqn.H']);self.assertIn('fvm::laplacian(rAUf, pPert)',source['pEqn.H'])
            self.assertIn('constrainPressure(pPert, U',source['pEqn.H']);self.assertIn('p_rgh == p - rho*gh;',source['pEqn.H'])
            self.assertIn('const volScalarField sc1RhoHydro',source['sc1CreateHydrostatic.H'])
            self.assertEqual(source['interFoam.C'].count('#include "sc1CreateHydrostatic.H"'),1)
            self.assertLess(source['interFoam.C'].index('#include "sc1CreateHydrostatic.H"'),source['interFoam.C'].index('while (pimple.run'))
            self.assertIn('runTime.value()>0.010000000001',source['interFoam.C'])
            before=candidate.base.source_hashes(base/'src')
            old_id=dict(binary=str(base/'bin/sc1InterFoamReference'),assets={'lib':'hash'},dependencies={'lib':'/lib'})
            new_id=dict(binary=str(p/'new/bin'/build.APP),assets={'lib':'hash'},sha256='hash',dependencies={'lib':'/lib'})
            with (patch.object(build.reference,'verify',return_value=old_id),patch.object(build,'identity',return_value=new_id),
                patch.object(build.startup,'require_allocation'),patch.object(build.shutil,'which',return_value='wmake'),
                patch.dict(os.environ,WM_PROJECT_VERSION='9',WM_PROJECT_DIR=str(project)),
                patch.object(build.base,'execute'),redirect_stdout(io.StringIO())):
                out=build.build(base,p/'new');self.assertEqual(build.verify(out),new_id);self.assertEqual(before,candidate.base.source_hashes(base/'src'))
                (out/'src/pEqn.H').write_text('changed')
                with self.assertRaises(ValueError):build.verify(out)


class Reports(unittest.TestCase):
    def test_requires_reference_physical_pressure_convergence_and_complete_budgets(self):
        text=log();r,_=report.budget(text,old.INITIAL);self.assertTrue(r['complete']);self.assertTrue(r['pressureIteration']['pressure_iteration_pass'])
        for t in (text.replace(HEADER,''),text.replace('frozen=1','frozen=0'),text.replace('rhoError_kg_m3=0','rhoError_kg_m3=nan'),text.replace('End\n',''),text.replace('Solving for pPert,','Solving for p_rgh,'),text.replace('min_Pa=99983','min_Pa=999')):
            with self.assertRaises(ValueError):report.budget(t,old.INITIAL)

    def test_snapshots_require_pressure_identity_and_both_phase_speeds(self):
        with tempfile.TemporaryDirectory() as d:
            case=Path(d);(case/'log.solver').write_text(log());corrupt=False;gas_speed=.0001
            def read(p,n,vector=False):
                t=float(p.parent.name);a=np.array([1-t*1e-5,0.]);c=np.array([[0,0,0],[0,1.5,0]])
                ph=np.array([110000,110000-1.5*9.81])
                return {'C':c,'V':np.array([1.,.642]),'alpha.liquid':a,'U':np.array([[.0001,0,0],[gas_speed,0,0]]),
                    'p':ph,'p_rgh':ph+(1+999*a)*9.81*c[:,1],'pPert':np.ones(2) if corrupt else np.zeros(2)}[p.name]
            masks={'walls':np.array([True,True]),'propeller':np.array([False,False])}
            with patch.object(report,'internal',side_effect=read),patch.object(report,'patch_cells',return_value=masks):
                r=report.analyze(case,old.INITIAL,2);self.assertTrue(r['hydrostatic_pass']);self.assertEqual(len(r['snapshots']),40)
                gas_speed=.002;(case/'log.solver').write_text(log().replace('gasMax_m_s=.0001','gasMax_m_s=.002'))
                self.assertFalse(report.analyze(case,old.INITIAL,2)['hydrostatic_pass'])
                gas_speed=.0001;(case/'log.solver').write_text(log().replace('Initial residual = 5e-13','Initial residual = 2e-12'))
                self.assertFalse(report.analyze(case,old.INITIAL,2)['hydrostatic_pass'])
                (case/'log.solver').write_text(log())
                corrupt=True
                with self.assertRaisesRegex(ValueError,'pressure identity'):report.analyze(case,old.INITIAL,2)


class Cases(unittest.TestCase):
    def fixture(self,p):
        root,m,_,_=force_tests.Cases().fixture(p)
        m['initial'].update(cells=2,interface_midheight_m=1.0185981622921294,transition_width_m=.0125)
        fr=p/'force';fr.mkdir();identity=dict(root=str(fr),build=str(p/'force_build'),manifest_sha256='manifest',log_sha256='log',nonorth=dict(root=str(root)))
        br=p/'hydro_build';br.mkdir();(br/'hydro_build_manifest.json').write_text(json.dumps({'base_build':m['build']}))
        binary=dict(binary=str(br/'bin'/build.APP),assets={'hydro solver':'hash'})
        return root,m,fr,identity,br,binary

    def patches(self,m,identity,binary):
        s=ExitStack();s.enter_context(patch.object(candidate,'control',return_value=(m,identity)))
        s.enter_context(patch.object(candidate,'check_initial'));s.enter_context(patch.object(build,'verify',return_value=binary))
        s.enter_context(patch.object(candidate.startup,'require_allocation'));s.enter_context(redirect_stdout(io.StringIO()))
        return s

    def test_fresh_preparation_preserves_existing_fields(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);root,m,fr,identity,br,binary=self.fixture(p)
            with self.patches(m,identity,binary):
                before=candidate.startup.tree_hashes(root);out=candidate.prepare(fr,br,p/'candidate');n=candidate.verify(out)
                self.assertEqual(before,candidate.startup.tree_hashes(root))
                self.assertEqual({k for k in m['inputs'] if m['inputs'][k]!=n['inputs'][k]}, {'system/controlDict','system/fvSchemes','system/fvSolution'})
                self.assertEqual(set(n['inputs'])-set(m['inputs']),{'0/pPert','constant/sc1HydrostaticReference'})
                self.assertNotIn('100000',(out/'quiescent/0/pPert').read_text());self.assertFalse((out/'quiescent/processor0').exists())
                for target in (out,root/'nested',fr/'nested',br/'nested',Path(identity['build'])/'nested'):
                    with self.assertRaises(ValueError):candidate.prepare(fr,br,target)
                (out/'quiescent/0/p_rgh').write_text('tamper')
                with self.assertRaises(ValueError):candidate.verify(out)

    def test_run_requires_written_pressures_and_bundles_failures(self):
        for mode in ('complete','failure','missing'):
            with self.subTest(mode=mode),tempfile.TemporaryDirectory() as d:
                p=Path(d);root,m,fr,identity,br,binary=self.fixture(p)
                with self.patches(m,identity,binary):
                    out=candidate.prepare(fr,br,p/'candidate');before=candidate.startup.tree_hashes(root)
                    def execute(command,case,logfile,env):
                        self.assertIn(binary['binary'],command);logfile.write_text(log())
                        if mode=='failure':raise ValueError('controlled stop')
                        for rank in range(8):
                            for i in range(5,201,5):
                                for name in ('U','alpha.liquid','p_rgh','p','pPert'):
                                    if mode=='missing' and rank==7 and i==200 and name=='pPert':continue
                                    q=case/f'processor{rank}'/f'{i*5e-5:.12g}'/name;q.parent.mkdir(parents=True,exist_ok=True);q.write_text('saved')
                    with patch.object(watch,'execute',side_effect=execute),patch.object(candidate.base,'execute') as shell,patch.object(report,'analyze',return_value={'hydrostatic_pass':True,'snapshots':[{}]}) as analyze:
                        if mode=='complete':
                            self.assertTrue(candidate.run(out)['hydrostatic_pass']);self.assertEqual(shell.call_count,2)
                        else:
                            with self.assertRaises(ValueError):candidate.run(out)
                            self.assertEqual(shell.call_count,1);analyze.assert_not_called()
                    self.assertEqual(before,candidate.startup.tree_hashes(root));self.assertTrue((p/'candidate_evidence/reactor_failure_evidence.tar.gz').is_file())
                    with self.assertRaisesRegex(ValueError,'attempted'):candidate.run(out)


if __name__=='__main__':unittest.main()
