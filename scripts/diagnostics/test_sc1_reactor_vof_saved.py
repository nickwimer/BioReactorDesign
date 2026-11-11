"""Saved-field mapping, failed-prefix budgets and preservation without OpenFOAM."""
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import sc1_reactor_vof_saved as saved
from sc1_reactor_startup_case import field
from test_sc1_reactor_failure import vof_step


def log():
    lines=['SC1_VOF_INITIAL liquid_m3=1 gas_m3=.642 version=1']
    for i in range(1,169):
        t=i*5e-5;lines.append(f'Time = {t:.12g}')
        lines += [f'SC1_VOF_PRESSURE t_s={t:.17g} min_Pa=100000 max_Pa=120000 nonfinite=0']*4
        lines.append(vof_step(t_s=t,gas_m3=.642,liquidMax_m_s=.001003 if i==168 else .0001))
    return '\n'.join(lines)+'\n[1] --> FOAM FATAL ERROR:\nVOF reference failed unchanged startup screens\n'


def write_labels(p,name,values):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(f'FoamFile {{ format ascii; class labelList; object {name}; }}\n{len(values)}\n('
        +' '.join(map(str,values))+')\n')


def write_field(p,name,values):
    a=np.asarray(values);vector=a.ndim==2;kind='vector' if vector else 'scalar'
    payload='nonuniform List<'+kind+'>\n'+str(len(a))+'\n(\n'
    payload+='\n'.join('('+' '.join(f'{x:.17g}' for x in row)+')' if vector else f'{row:.17g}' for row in a)+'\n)'
    p.parent.mkdir(parents=True,exist_ok=True);p.write_text(field(name,'[0 0 0 0 0 0 0]',payload,vector))


class SavedVof(unittest.TestCase):
    def test_failed_prefix_reconciles_without_accepting_startup(self):
        initial={'liquid_volume_m3':1,'gas_volume_m3':.642}
        r,steps=saved.budget(log(),initial)
        self.assertEqual(r['observedSteps'],168);self.assertEqual(r['pressureRecords'],672)
        self.assertEqual(r['phases']['liquid']['residualL1_kg'],0)
        self.assertTrue(r['phases']['gas']['withinBudgetTolerance'])
        self.assertGreater(steps[-1]['liquidMax_m_s'],.001)
        self.assertNotIn('reference_pass',r)
        for text in (log().replace('alphaCalls=2','alphaCalls=3',1),
                     log().replace('gas_m3=0.642','gas_m3=nan',1),
                     log().replace('gasOut_m3_s=0','gasOut_m3_s=0.1',1),
                     log().replace('liquidL1_kg=0','liquidL1_kg=.1',1),
                     log().replace('nonfinite=0','nonfinite=1',1),
                     log().replace('dt_s=5e-05','dt_s=0.0001',1),
                     log().replace('Time = 0.0084','Time = 0.0084\n'+vof_step(t_s=.0084))):
            with self.subTest(text=text[-80:]),self.assertRaises(ValueError):saved.budget(text,initial)

    def test_known_nonzero_throughput_is_reconciled_independently(self):
        from test_sc1_reactor_vof import good_log,INITIAL
        text=good_log().split('Time = 0.00845\n')[0]
        r,_=saved.budget(text,INITIAL)
        self.assertAlmostEqual(r['phases']['liquid']['boundaryOut_kg'],8.4e-5)
        self.assertLess(abs(r['phases']['gas']['physicalResidual_kg']),1e-15)

    def fixture(self,root):
        case=root/'quiescent';c=np.array([[0,0,0],[0,1,0],[.1,2,0.]])
        a=np.array([1.,1.,0.]);u=np.array([[.0001,0,0],[.0001,0,0],[.1,0,0]])
        originals={'C':c,'V':np.array([.6,.4,.642]),'alpha.liquid':a,
            'p':120000-(1+999*a)*9.81*c[:,1],'p_rgh':np.full(3,120000.),'U':np.zeros((3,3))}
        for name,values in originals.items():write_field(case/'0'/name,name,values)
        write_labels(case/'constant/polyMesh/owner','owner',[0,1,0,1,2])
        write_labels(case/'constant/polyMesh/neighbour','neighbour',[1,2])
        (case/'constant/polyMesh/boundary').write_text('FoamFile { format ascii; class polyBoundaryMesh; object boundary; }\n3\n(\n'
            'propeller { type wall; nFaces 1; startFace 2; }\nwalls { type wall; nFaces 1; startFace 3; }\n'
            'outlet { type patch; nFaces 1; startFace 4; }\n)\n')
        inputs=saved.base.case_input_hashes(case)
        for rank,ids in enumerate((np.array([2,0]),np.array([1]))):
            folder=case/f'processor{rank}'
            write_labels(folder/'constant/polyMesh/cellProcAddressing','cellProcAddressing',ids)
            for name in ('C','alpha.liquid','U','p_rgh'):write_field(folder/'0'/name,name,originals[name][ids])
            for time in saved.TIMES:
                for name,values in (('alpha.liquid',a),('U',u),('p_rgh',originals['p_rgh'])):
                    write_field(folder/time/name,name,values[ids])
        (case/'log.solver').write_text(log())
        m=dict(root=str(root),variant='sc1_meml_vof_reference_v1',ranks=2,
            initial=dict(cells=3,liquid_volume_m3=1,gas_volume_m3=.642,interface_midheight_m=1.5),
            build=str(root.parent/'build'),archive=str(root.parent/'archive'),regression_root=str(root.parent/'rates'),
            baseline={'root':str(root.parent/'baseline')},control={'root':str(root.parent/'control')},
            inputs=inputs,assets={'solver':'fixture'},source_commit='test')
        (root/'startup_manifest.json').write_text(json.dumps(m))
        (case/'attempt.json').write_text(json.dumps({k:m[k] for k in ('build','assets','ranks')}))
        return case,m

    def run_fixture(self,root,out):
        with patch.multiple(saved,CELLS=3,RANKS=2,MANIFEST_SHA=saved.sha(root/'startup_manifest.json'),
                LOG_SHA=saved.sha(root/'quiescent/log.solver')),redirect_stdout(io.StringIO()):
            return saved.inspect(root,out)

    def test_maps_partition_order_matches_audits_and_preserves_every_source_file(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)/'source';case,m=self.fixture(root);out=Path(d)/'report'
            before={str(p):saved.sha(p) for p in root.rglob('*') if p.is_file()}
            r=self.run_fixture(root,out)
            self.assertEqual(before,{str(p):saved.sha(p) for p in root.rglob('*') if p.is_file()})
            last=r['snapshots'][-1];peak=last['phases']['gas']['present']['peak']
            self.assertEqual(peak['cell'],2);self.assertEqual(peak['velocity_m_s'],[.1,0,0])
            self.assertEqual(peak['adjacentPatches'],['outlet'])
            self.assertAlmostEqual(last['liquidVolume_L'],1000)
            stencil=next(s for s in last['peakNeighbourStencils'] if s['peakCell']==2)
            self.assertEqual([r['cell'] for r in stencil['cells']],[2,1])
            with tarfile.open(out/'vof_spatial_evidence.tar.gz') as t:
                self.assertEqual(t.extractfile('source/quiescent/log.solver').read().decode(),log())
                self.assertEqual(len(t.getmembers()),4)
            with self.assertRaises(ValueError):self.run_fixture(root,out)
            with self.assertRaises(ValueError):self.run_fixture(root,root/'nested')

    def test_rejects_missing_duplicate_or_out_of_range_processor_addressing(self):
        for ids in ([2,1],[2],[3,0]):
            with tempfile.TemporaryDirectory() as d:
                root=Path(d)/'source';case,_=self.fixture(root)
                write_labels(case/'processor0/constant/polyMesh/cellProcAddressing','cellProcAddressing',ids)
                with self.assertRaisesRegex(ValueError,'permutation'):self.run_fixture(root,Path(d)/'report')

    def test_rejects_wrong_mapping_native_field_mismatch_and_redirection(self):
        for kind in ('mapping','velocity','missing','symlink','initial','during_read'):
            with self.subTest(kind=kind),tempfile.TemporaryDirectory() as d:
                root=Path(d)/'source';case,_=self.fixture(root)
                if kind=='mapping':write_labels(case/'processor0/constant/polyMesh/cellProcAddressing','cellProcAddressing',[0,2])
                elif kind=='velocity':write_field(case/'processor0/0.00825/U','U',[[.2,0,0],[.0001,0,0]])
                elif kind=='initial':(case/'0/C').write_text('changed')
                elif kind in ('missing','symlink'):
                    p=case/'processor0/0.00825/U';p.unlink()
                    if kind=='symlink':p.symlink_to(case/'processor0/0.004/U')
                if kind=='during_read':
                    original=saved.measure
                    def change(*args,**kwargs):
                        r=original(*args,**kwargs)
                        with (case/'processor0/0.00825/U').open('a') as f:f.write('\n// changed\n')
                        return r
                    with patch.object(saved,'measure',side_effect=change),self.assertRaises(ValueError):
                        self.run_fixture(root,Path(d)/'report')
                else:
                    with self.assertRaises(ValueError):self.run_fixture(root,Path(d)/'report')


if __name__=='__main__':unittest.main()
