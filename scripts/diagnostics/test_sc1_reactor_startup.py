"""Native inventory/pressure checks, complete budget reports and mocked OF workflow."""
from contextlib import redirect_stdout
import copy
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import sc1_reactor_startup as run
import sc1_reactor_startup_case as case
from sc1_two_field_cases import foam
from test_sc1_two_field_bounds import records
from test_sc1_two_field import render


def array_field(obj, values):
    a=np.asarray(values); vector=a.ndim==2
    payload='\n'.join('('+' '.join(str(z) for z in row)+')' if vector else str(row) for row in a)
    return foam(obj,'internalField nonuniform List<'+('vector' if vector else 'scalar')+f'>\n{len(a)}\n(\n{payload}\n);',
        'volVectorField' if vector else 'volScalarField')


def log():
    rows=records('bound_fine_gas_reference')
    for kind,r in rows:
        r['t_s']=round(r['t_s']/10,12)
        if 'dt_s' in r:r['dt_s']=case.DT
        if kind in ('SC1_TWO_INITIAL','SC1_TWO_STEP'):
            r['volume_m3']=1;r['mass_kg']=1 if r['phase']=='gas' else 1000
        if kind=='SC1_TWO_STEP':
            r['optionSource_kg']=r['boundaryOut_kg']=0
            r['alphaMin']=0;r['alphaMax']=1
        if kind=='SC1_PHASE_TERM':r['volume_m3']=r['frozenRhoMass_kg']=0
    states=[]
    for s in range(case.STEPS):
        for c in (1,2):
            states.append(f'SC1_RATE_STATE t_s={round((s+1)*case.DT,12)} call={c} valid=1 applied=2 '
                'qGasMin=0 qGasMax=0 qLiquidMin=0 qLiquidMax=0 alphaSumMaxError=0 rawIdentityScaled=0')
    return render(rows).replace('SC1_TWO_LIBRARY version=1','SC1_RATE_LIBRARY version=1\nSC1_TWO_LIBRARY version=1').replace(
        '\nEnd\n','\n'+'\n'.join(states)+'\nCourant Number mean: 0 max: 0\nEnd\n')


class Initializer(unittest.TestCase):
    def test_native_volume_and_hydrostatic_gradient(self):
        y=np.linspace(0,2,20001); c=np.column_stack((y*0,y,y*0));v=np.full(len(y),2/len(y))
        files,meta=case.initialized_fields(c,v)
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            for key,value in files.items():(p/Path(key).name).write_text(value)
            a=case.internal(p/'alpha.liquid',len(y));b=case.internal(p/'alpha.gas',len(y))
            pressure=case.internal(p/'p',len(y));prgh=case.internal(p/'p_rgh',len(y))
        self.assertAlmostEqual(float(np.dot(a,v)),1,13)
        np.testing.assert_allclose(a+b,1,atol=0)
        np.testing.assert_allclose(prgh,pressure+(1+999*a)*case.G*y,rtol=1e-14)
        mid=(y>.995)&(y<1.005)
        np.testing.assert_allclose(np.gradient(pressure,y)[mid],-(1+999*a[mid])*case.G,rtol=1e-7)
        np.testing.assert_allclose(prgh[y>1.1],100000,atol=1e-10)
        self.assertTrue(np.all(np.diff(pressure)<=0))
        self.assertAlmostEqual(meta['interface_midheight_m'],1,12)
        with self.assertRaises(ValueError):case.initialized_fields(c,v,target=3)

    def test_reader_rejects_corruption_and_nonfinite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'V';valid=array_field('V',[1,2])
            for bad in (valid.replace('ascii','binary'),valid.replace('\n2\n(','\n3\n('),
                        valid.replace('\n1\n','\nnan\n'),valid.replace('\n2\n);','\nbad\n);')):
                p.write_text(bad)
                with self.assertRaises(ValueError):case.internal(p,2)
            p.write_text(foam('V','internalField uniform 2;'))
            np.testing.assert_array_equal(case.internal(p,3),[2,2,2])

    def test_physics_and_boundary_recipe(self):
        f=case.case_files(8)
        self.assertNotIn('constant/MRFProperties',f)
        self.assertNotIn('system/blockMeshDict',f)
        self.assertNotIn('massSource',f['constant/fvModels'])
        self.assertIn('type segregated;',f['constant/phaseProperties'])
        self.assertNotIn('SchillerNaumann',f['constant/phaseProperties'])
        self.assertIn('interfaceCompression ((gas and liquid) 1)',f['constant/phaseProperties'])
        for p in ('gas','liquid'):
            self.assertEqual(f['0/U.'+p].count('type noSlip;'),2)
            self.assertEqual(f['0/alpha.'+p].count('type empty;'),4)
        self.assertIn('inletValue uniform 0',f['0/alpha.liquid'])
        self.assertIn('inletValue uniform 1',f['0/alpha.gas'])
        self.assertIn('numberOfSubdomains 8;',f['system/decomposeParDict'])
        with self.assertRaises(ValueError):case.case_files(0)


class Reports(unittest.TestCase):
    def setup_case(self,p):
        for step in (0,*range(case.WRITE,case.STEPS+1,case.WRITE)):
            folder=p/f'{step*case.DT:.12g}';folder.mkdir()
            for name,a in {'alpha.liquid':[1,1,0,0],'alpha.gas':[0,0,1,1],
                           'U.liquid':np.zeros((4,3)),'U.gas':np.zeros((4,3))}.items():
                (folder/name).write_text(array_field(name,a))
        (p/'0/V').write_text(array_field('V',[.5]*4));(p/'log.solver').write_text(log())
        return {'initial':{'gas_volume_m3':1,'liquid_volume_m3':1}}

    def test_complete_zero_feed_report_and_quiet_failure(self):
        with tempfile.TemporaryDirectory() as d,patch.object(run,'CELLS',4):
            p=Path(d);m=self.setup_case(p)
            r=run.analyze(p,m)
            self.assertTrue(r['startup_pass'],r['errors'])
            self.assertEqual(r['positiveRateCalls'],0)
            self.assertEqual(len(r['snapshots']),5)
            self.assertEqual(r['phases']['liquid']['optionSource_kg'],0)
            u=np.zeros((4,3));u[0,0]=.002
            (p/'0.004/U.liquid').write_text(array_field('U.liquid',u))
            r=run.analyze(p,m)
            self.assertTrue(r['screen_pass']);self.assertFalse(r['startup_pass'])
            self.assertIn('Quiescent liquid speed exceeded 1 mm/s',r['quiescent_errors'])

    def test_incomplete_audit_and_final_mass_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as d,patch.object(run,'CELLS',4):
            p=Path(d);m=self.setup_case(p)
            (p/'log.solver').write_text(log().replace('SC1_TWO_STEP phase=gas t_s=0.01 ','REMOVED phase=gas t_s=0.01 '))
            with self.assertRaises(ValueError):run.analyze(p,m)
            (p/'log.solver').write_text(log())
            (p/'0.01/alpha.liquid').write_text(array_field('alpha.liquid',[.5,1,0,0]))
            (p/'0.01/alpha.gas').write_text(array_field('alpha.gas',[.5,0,1,1]))
            with self.assertRaisesRegex(ValueError,'inventory disagrees'):run.analyze(p,m)

    def test_complete_mass_failure_is_reported_without_relaxing_threshold(self):
        with tempfile.TemporaryDirectory() as d,patch.object(run,'CELLS',4):
            p=Path(d);m=self.setup_case(p);lines=log().splitlines()
            for i,line in enumerate(lines):
                if line.startswith('SC1_TWO_STEP phase=liquid '):
                    lines[i]=line.replace('boundaryOut_kg=0 ','boundaryOut_kg=0.01 ').replace(
                        'physicalResidual_kg=0 ','physicalResidual_kg=0.01 ')
                    break
            (p/'log.solver').write_text('\n'.join(lines)+'\n')
            r=run.analyze(p,m)
            self.assertFalse(r['startup_pass']);self.assertFalse(r['screen_pass'])
            self.assertTrue(r['diagnostic_pass'],r['errors'])
            self.assertEqual(r['errors'],['liquid: phase mass budget exceeds screening tolerance'])


class Workflow(unittest.TestCase):
    def test_prepare_verify_preserves_mesh_and_detects_edits(self):
        with tempfile.TemporaryDirectory() as d,redirect_stdout(io.StringIO()):
            p=Path(d);archive=p/'archive';mesh=archive/'constant/polyMesh';mesh.mkdir(parents=True)
            for name in (*run.PINNED,'points','faces'):(mesh/name).write_text('fixture '+name)
            hashes=run.tree_hashes(mesh)
            adj=dict(inputs={k:{'sha256':hashes[k]} for k in run.PINNED},
                patches=[dict(patch=k,type=v) for k,v in case.PATCHES.items()])
            reg=p/'regression';reg.mkdir();build=p/'build';build.mkdir()
            c=np.array([[0,y,0] for y in (.999,1.01,1.022,1.033)])
            def execute(cmd,cwd,log,env):
                log.write_text('Mesh OK.\n')
                if 'writeCellCentres' in cmd:(cwd/'0/C').write_text(array_field('C',c))
                if 'writeCellVolumes' in cmd:(cwd/'0/V').write_text(array_field('V',[.5]*4))
            with (patch.object(run,'regression',return_value=({'build':str(build),'assets':{}},'proof')),
                 patch.object(run,'zone_adjacency',return_value=adj),patch.object(run,'PINNED',{k:hashes[k] for k in run.PINNED}),
                 patch.object(run,'CELLS',4),patch.object(run.base,'environment',return_value={}),
                 patch.object(run.base,'execute',side_effect=execute),patch.object(run.base,'verify_build',return_value={'assets':{}})):
                root=run.prepare(reg,archive,p/'out',8)
                m=run.verify(root)
                self.assertAlmostEqual(m['initial']['liquid_volume_m3'],1,12)
                self.assertEqual(run.tree_hashes(mesh),hashes)
                with self.assertRaises(ValueError):run.prepare(reg,archive,p/'out',8)
                (mesh/'points').write_text('changed')
                with self.assertRaisesRegex(ValueError,'Archived mesh changed'):run.verify(root)

    def test_parallel_launch_reconstruct_and_no_reuse(self):
        with tempfile.TemporaryDirectory() as d,redirect_stdout(io.StringIO()):
            root=Path(d);p=root/'quiescent';p.mkdir()
            for n in ('0','constant','system'):(p/n).mkdir()
            for n in ('checkMesh','writeCellCentres','writeCellVolumes'):(p/('log.'+n)).write_text('ok')
            m={'build':'/private','assets':{},'ranks':8,'inputs':{}}
            calls=[]
            def execute(cmd,cwd,log,env):calls.append(cmd);log.write_text('End\n')
            with (patch.object(run,'verify',return_value=m),patch.object(run,'require_allocation'),
                 patch.object(run.base,'environment',return_value={}),patch.object(run.base,'execute',side_effect=execute),
                 patch.object(run,'analyze',return_value={'startup_pass':True}),patch.object(run.base,'case_input_hashes',return_value={})):
                self.assertTrue(run.run(root)['startup_pass'])
                self.assertEqual([c[0] for c in calls],['decomposePar','srun','reconstructPar'])
                self.assertEqual(calls[1][1:3],['-n','8'])
                self.assertIn('0.002:0.01',calls[-1])
                with self.assertRaisesRegex(ValueError,'attempted'):run.run(root)
                self.assertEqual(len(calls),3)

    def test_allocation_and_cli_gate(self):
        with patch.dict(os.environ,{},clear=True):
            with self.assertRaises(ValueError):run.require_allocation(8)
        with patch.dict(os.environ,{'SLURM_JOB_ID':'1','SLURM_JOB_CPUS_PER_NODE':'4'},clear=True):
            with self.assertRaises(ValueError):run.require_allocation(8)
        with patch.dict(os.environ,{'SLURM_JOB_ID':'1','SLURM_JOB_CPUS_PER_NODE':'4(x2)'},clear=True):
            run.require_allocation(8)
        with self.assertRaises(SystemExit):run.main(['run','--root','/tmp/unused'])


if __name__=='__main__':unittest.main()
