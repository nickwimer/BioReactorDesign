"""Spatial diagnostic isolation, prefix checks and native-field peak selection."""
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import sc1_reactor_velocity_probe as probe
import sc1_reactor_velocity_fields as fields
import sc1_reactor_failure as failure
import sc1_reactor_watch as watch
from sc1_reactor_startup_case import case_files, field


def sample_log():
    lines=[]
    for t in probe.TIMES:
        lines.append('Time = '+t)
        for i in range(4):
            lines.append(f'SC1_PRESSURE_FIELD sequence={i} term=pUnclipped min=99983 max=111800 nonfinite=0')
            lines.append('SC1_COUPLED_FLUX maxAbsError_m3_s=1e-18 maxAbsMixture_m3_s=1e-6')
        for phase in ('gas','liquid'):lines.append(f'SC1_TWO_STEP phase={phase} t_s={t} physicalResidual_kg=1e-17')
        lines.append(f'SC1_QUIESCENT_STEP t_s={t} maxLiquidPresent_m_s=0.0008')
    return '\n'.join(lines)+'\nEnd\n'


class Prefix(unittest.TestCase):
    def test_three_steps_match_but_do_not_accept_the_parent_failed_step(self):
        text=sample_log();parent=text.replace('End\n','')+'Time = 0.0002\nSC1_QUIESCENT_STEP t_s=0.0002 maxLiquidPresent_m_s=0.00106\nFOAM FATAL ERROR\n'
        self.assertEqual(len(probe.compare_prefix(text,parent)),33)
        with self.assertRaisesRegex(ValueError,'trajectory'):
            probe.compare_prefix(text.replace('0.0008','0.0009',1),parent)
        with self.assertRaisesRegex(ValueError,'Incomplete'):
            probe.compare_prefix(text.replace('SC1_COUPLED_FLUX','omitted',1),parent)
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'log.solver').write_text(parent)
            with self.assertRaisesRegex(ValueError,'exactly three'):probe.analyze(p,{},{} )

    def test_output_recipe_changes_only_duration_and_write_interval(self):
        old=case_files(8,segregated_drag=True,face_momentum=True,coupled_drag=True)['system/controlDict']
        new=probe.output_control(old)
        self.assertEqual(new.replace('endTime 0.00015;','endTime 0.01;').replace('writeInterval 1;','writeInterval 5;'),old)
        self.assertIn('deltaT 5e-05;',new)
        with self.assertRaises(ValueError):probe.output_control(new)


class NativeFields(unittest.TestCase):
    def test_boundary_adjacency_uses_face_owner_cells(self):
        with tempfile.TemporaryDirectory() as d:
            case=Path(d);mesh=case/'constant/polyMesh';mesh.mkdir(parents=True)
            (mesh/'owner').write_text('FoamFile { format ascii; class labelList; object owner; }\n6\n(0 1 2 3 0 2)\n')
            (mesh/'boundary').write_text('FoamFile { format ascii; class polyBoundaryMesh; object boundary; }\n3\n(\n'
                'walls { type wall; nFaces 2; startFace 2; }\n'
                'propeller { type wall; nFaces 1; startFace 4; }\n'
                'outlet { type patch; nFaces 1; startFace 5; }\n)\n')
            masks=fields.patch_cells(case,4)
            np.testing.assert_array_equal(np.flatnonzero(masks['walls']),[2,3])
            np.testing.assert_array_equal(np.flatnonzero(masks['propeller']),[0])
            np.testing.assert_array_equal(np.flatnonzero(masks['outlet']),[2])

    def write_field(self,p,name,data):
        data=np.asarray(data);vector=data.ndim==2
        payload='\n'.join('('+' '.join(str(v) for v in row)+')' for row in data) if vector else '\n'.join(str(v) for v in data)
        value=f'nonuniform List<{"vector" if vector else "scalar"}>\n{len(data)}\n(\n{payload}\n)'
        (p/name).write_text(field(name,'[0 0 0 0 0 0 0]',value,vector))

    def test_peak_masks_locations_counterflow_and_pressure(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'0.00015';p.mkdir()
            c=np.array([[.4,.2,0],[.1,1.018,0],[.2,1.2,0],[0,.7,.2]])
            a=np.array([1.,.5,0.,.99]);v=np.ones(4)*.25;p0=np.ones(4)*100000
            ul=np.array([[.0001,0,0],[0,.0008,0],[3,0,0],[0,.0003,0]])
            ug=np.array([[9,0,0],[0,-.04,0],[.001,0,0],[.01,0,0]])
            for name,data in {'alpha.liquid':a,'alpha.gas':1-a,'U.liquid':ul,'U.gas':ug,
                'p_rgh':p0+(1+999*a)*fields.G*c[:,1]}.items():self.write_field(p,name,data)
            patches={'walls':np.array([True,False,False,False]),'propeller':np.array([False,False,False,True]),'outlet':np.array([False,False,True,False])}
            report=fields.measure(p,c,v,a,p0,patches,1.018)
            liquid=report['phases']['liquid'];self.assertEqual(liquid['present']['peak']['cell'],1)
            self.assertEqual(liquid['bulk']['peak']['cell'],3)
            self.assertEqual(liquid['wallAdjacent']['peak']['cell'],0)
            self.assertEqual(liquid['present']['peak']['adjacentPatches'],[])
            self.assertAlmostEqual(liquid['present']['peak']['pressureChangeFromInitial_Pa'],0)
            self.assertEqual(liquid['present']['peak']['gasVelocity_m_s'],[0,-.04,0])
            self.assertAlmostEqual(report['interfaceSlip']['max_m_s'],.0408)
            self.assertAlmostEqual(report['alphaChangeL1_m3'],0)
            self.write_field(p,'alpha.gas',np.array([-.1,.5,1,.01]))
            with self.assertRaisesRegex(ValueError,'fractions'):fields.measure(p,c,v,a,p0,patches,1.018)


class Cases(unittest.TestCase):
    def setup_control(self,p):
        old=p/'failed';case=old/'quiescent'
        for folder in ('0','constant','system'):(case/folder).mkdir(parents=True)
        (case/'0/alpha.liquid').write_text('original')
        (case/'system/controlDict').write_text(case_files(8,segregated_drag=True,face_momentum=True,coupled_drag=True)['system/controlDict'])
        for name in ('log.checkMesh','log.writeCellCentres','log.writeCellVolumes'):(case/name).write_text('original')
        (case/'log.solver').write_text('preserved failure');(case/'processor0').mkdir()
        build=p/'build';build.mkdir();(build/'coupling_build_manifest.json').write_text(json.dumps({'base_build':str(p/'rates')}))
        m=dict(root=str(old),variant=probe.coupled.TAG,inputs=probe.base.case_input_hashes(case),
            build=str(build),assets={'binary':'sha'},ranks=8,archive=str(p/'archive'),regression_root=str(p/'regression'),
            candidate_identity={'binary':str(build/'bin/birdPressureCoupled')},initial={'interface_midheight_m':1.018})
        for k in ('baseline','failed_control','face_control','coupled_control','pressure_probe'):m[k]={'root':str(p/k)}
        identity=dict(root=str(old),log_sha256=probe.FAILED_LOG,manifest_sha256='control')
        return old,m,identity

    def test_prepare_uses_same_binary_preserves_control_and_rejects_reuse_or_other_changes(self):
        with tempfile.TemporaryDirectory() as d,redirect_stdout(io.StringIO()):
            p=Path(d);old,parent,identity=self.setup_control(p);before=probe.startup.tree_hashes(old)
            with patch.object(probe,'control',return_value=(parent,identity)),patch.object(probe.startup,'require_allocation'):
                out=probe.prepare(old,p/'new');m,_=probe.verify(out)
                self.assertEqual(m['build'],parent['build']);self.assertEqual(m['assets'],parent['assets'])
                self.assertEqual([k for k in m['inputs'] if m['inputs'][k]!=parent['inputs'][k]],['system/controlDict'])
                self.assertEqual(probe.startup.tree_hashes(old),before)
                self.assertFalse((out/'quiescent/processor0').exists())
                (out/'quiescent/attempt.json').write_text('attempt')
                with self.assertRaisesRegex(ValueError,'attempted'):probe.run(out)
                (out/'quiescent/0/alpha.liquid').write_text('different')
                with self.assertRaisesRegex(ValueError,'inputs changed'):probe.verify(out)

    def test_success_reconstructs_only_three_steps_and_never_builds(self):
        with tempfile.TemporaryDirectory() as d,redirect_stdout(io.StringIO()):
            p=Path(d);old,parent,identity=self.setup_control(p)
            with patch.object(probe,'control',return_value=(parent,identity)),patch.object(probe.startup,'require_allocation'):
                out=probe.prepare(old,p/'new')
                def solver(command,cwd,log,env,**kwargs):
                    self.assertIn('--time=00:10:00',command);self.assertEqual(kwargs['wall_seconds'],600)
                    log.write_text(sample_log())
                    for rank in range(8):
                        for time in probe.TIMES:
                            folder=cwd/f'processor{rank}'/time;folder.mkdir(parents=True)
                            for name in fields.FIELDS:(folder/name).write_text('mock native field')
                result=dict(prefixAgreement={'records':33},snapshots=[],meaning='diagnostic only')
                with (patch.object(watch,'execute',side_effect=solver),patch.object(probe.base,'execute') as utilities,
                      patch.object(probe,'analyze',return_value=result),patch.object(failure,'collect') as collect):
                    self.assertEqual(probe.run(out),result)
                    self.assertEqual([c.args[0][0] for c in utilities.call_args_list],['decomposePar','reconstructPar'])
                    self.assertIn('5e-05:0.00015',utilities.call_args_list[-1].args[0]);collect.assert_called_once()

    def test_failed_solver_is_collected_without_reconstruction(self):
        with tempfile.TemporaryDirectory() as d,redirect_stdout(io.StringIO()):
            p=Path(d);old,parent,identity=self.setup_control(p)
            with patch.object(probe,'control',return_value=(parent,identity)),patch.object(probe.startup,'require_allocation'):
                out=probe.prepare(old,p/'new')
                def solver(command,cwd,log,env,**kwargs):log.write_text('FOAM FATAL ERROR\n');raise ValueError('failed')
                with (patch.object(watch,'execute',side_effect=solver),patch.object(probe.base,'execute') as utilities,
                      patch.object(failure,'collect') as collect):
                    with self.assertRaisesRegex(ValueError,'failed'):probe.run(out)
                    self.assertEqual(utilities.call_count,1);collect.assert_called_once()
                self.assertIn(probe.TAG,failure.VARIANTS)


if __name__=='__main__':unittest.main()
