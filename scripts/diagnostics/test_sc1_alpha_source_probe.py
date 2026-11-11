#!/usr/bin/env python3
"""Counterfactual bounds, exact diagnostic identity and isolated preparation."""
import contextlib
import difflib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import sc1_alpha_source_probe as probe
import test_sc1_alpha_preparation as fixtures

CPP=r'''
#include "sc1AlphaSourceProbe.H"
#include <cassert>
#include <string>
using namespace sc1SourceProbe;
void near(double a,double b) { assert(std::abs(a-b)<1e-12); }
int main(int argc,char** argv)
{
    assert(argc==2); const std::string test(argv[1]);
    if (test=="feasible_floor")
    {
        // Target removes the tiny donor exactly; legacy floor removes less.
        double a=1e-8,dt=1e-4,g=-1e-4,sp=-1,after=a/(1-dt*sp);
        auto r=cell(a,after,dt,0,0,sp,a,0,g,g,1,1000,1);
        near(r[rawBelow],0); near(r[rawAbove],0);
        near(r[rawExcessMass],0); assert(r[splitMass]>0);
        assert(r[rawCorrectionMass]<0); near(r[identityRemainderMass],0);
        near(r[maxDtSp],1e-4);
    }
    else if (test=="impossible_source")
    {
        auto r=cell(0,0,.1,0,0,0,0,0,-1,-1,2,1000,1);
        near(r[rawBelow],1); near(r[pressureBelow],1);
        near(r[rawExcessMass],200); near(r[rawCorrectionMass],-200);
        // Bounds violation is reported, never clipped away.
        near(r[identityRemainderMass],0);
    }
    else if (test=="transport_rescues_source")
    {
        auto r=cell(0,.01,.1,-1.1,0,0,0,0,-1,-1,1,1,1);
        near(r[rawBelow],0); near(r[rawAbove],0); near(r[rawExcessMass],0);
        near(r[rawCorrectionMass],0); near(r[identityRemainderMass],0);
    }
    else if (test=="pressure_and_raw_differ")
    {
        auto r=cell(.9,.9,.1,0,0,0,.9,0,2,0,1,1,1);
        near(r[rawAbove],1); near(r[pressureAbove],0);
        near(r[rawExcessMass],.1); near(r[pressureExcessMass],0);
    }
    else if (test=="alpha_max_and_tolerance")
    {
        auto r=cell(.4,.4,.1,0,0,0,.4,0,1,1,1,1,.45);
        near(r[rawAbove],1); near(r[rawExcessMass],.05);
        r=cell(1,1,.1,0,0,0,1,0,1e-13,1e-13,1,1,1);
        near(r[rawAbove],0); assert(r[rawExcessMass]>0);
    }
    else if (test=="bins")
    {
        assert(bin(0)==0 && bin(1e-5)==0 && bin(1e-4)==2);
        assert(bin(1)==1 && bin(1-1e-5)==1 && bin(1-1e-4)==2);
        assert(bin(-1e-16)==3 && bin(1+1e-6)==3);
    }
    else if (test=="accumulation")
    {
        Row a{}, b{}; a[samples]=3; b[samples]=2;
        a[splitMass]=.3; b[splitMass]=-.1;
        a[maxDtSp]=5; b[maxDtSp]=2; add(a,b);
        near(a[samples],5); near(a[splitMass],.2); near(a[maxDtSp],5);
        Row rank2{}; rank2[maxDtSp]=8; rank2[samples]=11; add(a,rank2);
        near(a[samples],16); near(a[maxDtSp],8);
    }
    else if (test=="independent_equation_error")
    {
        auto r=cell(.4,.7,.02,.3,.5,-.2,.4,.1,.2,.3,.02,1000,1);
        near(r[identityRemainderMass],0);
        assert(std::abs(r[rawCorrectionMass]+r[splitMass]+r[implicitMass])>1);
    }
    else return 2;
}
'''


class KernelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory(prefix='sc1-source-probe-')
        cls.addClassCleanup(cls.temp.cleanup)
        source=Path(cls.temp.name)/'test.C'; source.write_text(CPP)
        cls.binary=source.with_suffix('')
        compiler=shutil.which('g++') or shutil.which('c++')
        if not compiler:
            raise RuntimeError('C++ compiler required')
        subprocess.run([compiler,'-std=c++11','-Wall','-Wextra','-Werror',
                        '-I',str(probe.HERE/'include'),str(source),'-o',str(cls.binary)],
                       check=True,capture_output=True,text=True)

    def check(self,name):
        subprocess.run([str(self.binary),name],check=True,capture_output=True,text=True)

    def test_feasible_floor_correction(self): self.check('feasible_floor')
    def test_impossible_raw_source_is_not_clipped(self): self.check('impossible_source')
    def test_flux_can_make_a_negative_source_feasible(self): self.check('transport_rescues_source')
    def test_raw_and_pressure_targets_are_distinct(self): self.check('pressure_and_raw_differ')
    def test_phase_maximum_and_diagnostic_tolerance(self): self.check('alpha_max_and_tolerance')
    def test_assembly_alpha_bins(self): self.check('bins')
    def test_counts_sum_but_stiffness_uses_maximum(self): self.check('accumulation')
    def test_identity_includes_independent_equation_residual(self): self.check('independent_equation_error')


class PreparationTests(unittest.TestCase):
    def setUp(self):
        self.fixture=fixtures.PreparationTests()
        self.fixture.setUp(); self.addCleanup(self.fixture.doCleanups)
        self.out=self.fixture.out

    def prepare(self):
        return probe.prepare(self.fixture.source,self.out,self.fixture.installed,
                             base_prepare=self.fixture.base_prepare)

    def test_observation_preserves_every_original_header_line(self):
        original=(probe.alpha.PHASE_HEADERS/'sc1AlphaSubcycleAudit.H').read_text()
        transformed=probe.instrument_header(original)
        self.assertFalse([line for line in difflib.ndiff(original.splitlines(),transformed.splitlines())
                          if line.startswith('- ')])
        self.assertIn('("sc1AlphaSourceProbe",false)',transformed)
        self.assertIn('const auto& cumulative=totals_.sourceProbeByCall[callInStep_]',transformed)
        self.assertIn('listCombineGather(maxima,maxEqOp<scalar>())',transformed)
        self.assertLess(transformed.index('sc1SourceProbe::cell'),transformed.index('void finish()'))
        with self.assertRaises(ValueError): probe.instrument_header(transformed)

    def test_only_private_observer_and_mandatory_binding_change(self):
        original=probe.alpha.source_hashes(self.fixture.phase)
        live=probe.alpha.source_hashes(self.fixture.source)
        installed=probe.alpha.sha256(self.fixture.installed)
        manifest=self.prepare()
        self.assertEqual(probe.alpha.source_hashes(self.fixture.phase),original)
        self.assertEqual(probe.alpha.source_hashes(self.fixture.source),live)
        self.assertEqual(probe.alpha.sha256(self.fixture.installed),installed)
        a=json.loads((self.out/'alpha_audit_manifest.json').read_text())
        self.assertEqual(a['private_phase_source_sha256_after'],probe.alpha.source_hashes(self.out/'phaseSystems'))
        self.assertEqual(a['solver_source_sha256_after_alpha'],probe.alpha.source_hashes(self.out/'solver'))
        self.assertEqual(manifest['alpha_manifest_sha256'],probe.alpha.sha256(self.out/'alpha_audit_manifest.json'))
        main=(self.out/'solver/putidaFoam.C').read_text()
        self.assertIn('extern "C" int '+probe.MARKER+'();',main)
        self.assertIn('if (sc1AlphaSourceProbeLibraryVersion()!=1)',main)
        self.assertIn('LOCAL_BIOLOGY_SENTINEL',main)
        phase=(self.out/'phaseSystems/phaseSystem/phaseSystemSolve.C').read_text()
        baseline=(self.fixture.phase/'phaseSystem/phaseSystemSolve.C').read_text()
        self.assertFalse([line for line in difflib.ndiff(baseline.splitlines(),phase.splitlines())
                          if line.startswith('- ')])

    def test_existing_output_and_changed_header_are_rejected(self):
        self.prepare(); before=probe.alpha.source_hashes(self.out)
        with self.assertRaises(ValueError): self.prepare()
        self.assertEqual(probe.alpha.source_hashes(self.out),before)
        (self.out/'phaseSystems/phaseSystem'/probe.HEADER).write_text('changed')
        with self.assertRaisesRegex(ValueError,'changed source-probe'):
            probe.verify_probe(self.out,[],[])

    def test_missing_new_symbol_is_rejected_after_base_verification(self):
        self.prepare()
        def runner(command,env):
            return {'returncode':0,'stdout':'0001 T other\n','stderr':''}
        with patch.object(probe,'verify',return_value={'ok':True,'LD_LIBRARY_PATH':'/private','checks':[]}):
            with self.assertRaisesRegex(ValueError,'mandatory source-probe marker'):
                probe.verify_probe(self.out,[],[],runner=runner)

    def test_private_definition_and_executable_reference_are_required(self):
        self.prepare()
        def runner(command,env):
            kind='T' if command[-1].endswith('.so') else 'U'
            return {'returncode':0,'stdout':'0001 '+kind+' '+probe.MARKER+'\n','stderr':''}
        with patch.object(probe,'verify',return_value={'ok':True,'LD_LIBRARY_PATH':'/private','checks':[]}):
            self.assertTrue(probe.verify_probe(self.out,[],[],runner=runner)['ok'])

    def test_cli_cannot_redirect_output_through_symlink(self):
        self.out.symlink_to(self.fixture.root/'absent')
        with patch.dict(os.environ,{'FOAM_USER_LIBBIN':'/dummy'}),contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(probe.main(['prepare','--out',str(self.out)]),1)
        self.assertFalse((self.fixture.root/'absent').exists())

    def test_build_wrapper_stops_after_failed_phase_build(self):
        stub=self.fixture.root/'stub'; stub.mkdir()
        python=stub/'python'
        python.write_text('''#!/usr/bin/env python3
import pathlib,sys
if sys.argv[2]=='prepare':
    root=pathlib.Path(sys.argv[4]); root.mkdir()
    (root/'phaseSystems').mkdir(); (root/'solver').mkdir()
else:
    raise RuntimeError('Must not verify a failed build')
''')
        python.chmod(0o755)
        wmake=stub/'wmake'
        wmake.write_text('#!/usr/bin/env bash\nprintf "%s\\n" "$PWD"\nfalse\n')
        wmake.chmod(0o755)
        env=os.environ.copy(); env.update(SLURM_JOB_ID='test',SC1_PROBE_PY=str(python),
                                          PATH=str(stub)+':'+env['PATH'])
        result=subprocess.run(['bash',str(probe.HERE/'build_sc1_alpha_source_probe.sh'),
                               str(self.out)],env=env,text=True,capture_output=True)
        self.assertEqual(result.returncode,1)
        self.assertTrue((self.out/'log.build.phase').exists())
        self.assertFalse((self.out/'log.build.solver').exists())
        self.assertIn('Build stopped',result.stdout)

    def test_build_wrapper_requires_allocation_before_creating_anything(self):
        env=os.environ.copy(); env.pop('SLURM_JOB_ID',None)
        result=subprocess.run(['bash',str(probe.HERE/'build_sc1_alpha_source_probe.sh'),
                               str(self.out)],env=env,text=True,capture_output=True)
        self.assertEqual(result.returncode,1)
        self.assertFalse(self.out.exists())


if __name__=='__main__':
    unittest.main(verbosity=2)
