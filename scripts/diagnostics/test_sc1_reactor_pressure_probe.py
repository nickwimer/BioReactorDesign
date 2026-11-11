"""Pressure probe isolation, observation ordering, dependency identity and stop handling."""
from contextlib import redirect_stdout
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import sc1_reactor_pressure_probe as probe
import sc1_reactor_pressure_instrument as instrument
import sc1_reactor_failure as failure
import sc1_reactor_watch as watch
from sc1_two_field_rates_instrument import solver_source, RECIPE, MARKER


class Instrumentation(unittest.TestCase):
    def test_source_hooks_are_observers_with_prelimit_stop_and_two_step_bound(self):
        source=probe.base.SOURCE/'multiphaseEulerFoam'
        main=solver_source((source/'multiphaseEulerFoam.C').read_text())
        result=instrument.solver(main)
        self.assertIn('if (!faceMomentum || !partialElimination)',result)
        self.assertLess(result.index('sc1Stages.endStep();'),result.index('sc1Probe.endStep()'))
        with self.assertRaises(ValueError):instrument.solver(result)
        old=(source/'pUf/pEqn.H').read_text();new=instrument.pressure(old)
        self.assertEqual(new.count('pEqn.solve();'),1)
        self.assertLess(new.index('pEqn.solve();'),new.index('sc1Probe.pressure('))
        self.assertLess(new.index('sc1Probe.pressure('),new.index('fluid.partialEliminationf('))
        self.assertLess(new.index('sc1Probe.pressure('),new.index('fvConstraints.constrain(p);'))
        # Every original line remains in the same order: only observations inserted.
        pending=iter(new.splitlines())
        for line in old.splitlines():
            self.assertTrue(any(candidate==line for candidate in pending),line)
        with self.assertRaises(ValueError):instrument.pressure(new)
        self.assertIn('pEqn.source()',new)

    def test_actual_range_and_step_bound_cpp(self):
        compiler=shutil.which('g++')
        if not compiler:self.skipTest('g++ unavailable')
        source=probe.HEADER.read_text()
        range_code=source.split('    struct Range',1)[1].split('\n    void emit',1)[0]
        end_code=source.split('    bool endStep()',1)[1].split('\n};',1)[0]
        code='''#include <algorithm>
#include <cmath>
#include <cassert>
#include <limits>
#include <sstream>
using scalar=double; using label=int; using std::min; using std::max;
constexpr double GREAT=1e30; std::ostringstream Info; constexpr char endl='\\n';
struct Range'''+range_code+'''
struct Steps { int steps_=0; bool endStep()'''+end_code+'''\n};
int main() {
 Range r; r.add(7);r.add(-3);r.add(std::numeric_limits<double>::infinity());
 r.add(std::numeric_limits<double>::quiet_NaN());
 assert(r.lo==-3 && r.hi==7 && r.bad==2);
 Steps s;assert(!s.endStep());assert(s.endStep());
 assert(Info.str().find("SC1_PRESSURE_DONE completedSteps=2")!=std::string::npos);
}
'''
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'test.C').write_text(code)
            subprocess.run([compiler,'-std=c++14',str(p/'test.C'),'-o',str(p/'test')],check=True,capture_output=True,text=True)
            subprocess.run([str(p/'test')],check=True)


class Build(unittest.TestCase):
    def test_only_new_application_built_and_base_preserved(self):
        with tempfile.TemporaryDirectory() as d,redirect_stdout(io.StringIO()):
            p=Path(d);installed=p/'installed';installed.mkdir()
            parent=probe.base.prepare_build(p/'rates',installed=installed,recipe=RECIPE)
            before=probe.base.selected_hashes(parent/'src');out=p/'probe'
            def execute(command,cwd,log,env):
                self.assertEqual(command,['wmake']);self.assertEqual(cwd,out/'src')
                self.assertEqual(env['FOAM_USER_LIBBIN'],str(parent/'lib'))
                log.write_text('mock');(out/'bin'/probe.APP).write_text('mock binary')
            identity=dict(binary=str(out/'bin'/probe.APP),sha256='mock',dependencies={})
            with (patch.object(probe.startup,'require_allocation'),patch.dict(os.environ,WM_PROJECT_VERSION='9'),
                  patch.object(probe.shutil,'which',return_value='/mock/wmake'),
                  patch.object(probe.base,'verify_build',return_value={'assets':{}}),
                  patch.object(probe.base,'execute',side_effect=execute) as build,
                  patch.object(probe,'executable_identity',return_value=identity)):
                probe.build(parent,out);build.assert_called_once()
                self.assertEqual(probe.verify_build(out),identity)
                self.assertEqual(probe.base.selected_hashes(parent/'src'),before)
                self.assertFalse((out/'lib').exists())
                self.assertIn('-I'+str(parent/'src/phaseSystems/lnInclude'),(out/'src/Make/options').read_text())
                self.assertIn('EXE = '+str(out/'bin'/probe.APP),(out/'src/Make/files').read_text())
                with self.assertRaises(ValueError):probe.build(parent,out)
                (out/'src/pUf/pEqn.H').write_text('changed')
                with self.assertRaisesRegex(ValueError,'observation recipe'):probe.verify_build(out)

    def test_probe_runtime_dependencies_are_pinned(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'bin').mkdir();binary=root/'bin'/probe.APP;binary.write_text('binary');binary.chmod(0o755)
            m={'base_build':'/private/rates'}
            lines='\n'.join(f'{name} => /private/rates/lib/{name} (0x0)' for _,name in probe.base.COMPONENTS)
            def command(cmd,**kwargs):
                return subprocess.CompletedProcess(cmd,0,lines if cmd[0]=='ldd' else ' U '+MARKER+'\n','')
            with patch.object(probe.subprocess,'run',side_effect=command):
                identity=probe.executable_identity(root,m)
                self.assertEqual(identity['sha256'],probe.base.sha(binary))
                lines=lines.replace('/private/rates/lib/libbirdphaseSystem.so','/installed/libbirdphaseSystem.so')
                with self.assertRaisesRegex(ValueError,'project libraries'):probe.executable_identity(root,m)


class Cases(unittest.TestCase):
    def control(self,p):
        root=p/'coupled';case=root/'quiescent'
        for folder in ('0','constant','system'):(case/folder).mkdir(parents=True)
        for name in ('0/alpha.gas','constant/g','system/controlDict'):(case/name).write_text('original')
        for name in ('log.checkMesh','log.writeCellCentres','log.writeCellVolumes'):(case/name).write_text('saved')
        (case/'processor0').mkdir();(case/'log.solver').write_text('preserve failed log')
        (p/'probe').mkdir()
        parent=dict(variant=probe.startup.COUPLED_TAG,root=str(root),inputs=probe.base.case_input_hashes(case),
            ranks=8,build=str(p/'rates'),assets={'oldbinary':'oldsha'},initial={},archive=str(p/'archive'),
            mesh_hashes={},regression_root=str(p/'regression'),regression_sha256='proof',
            baseline={'root':str(p/'v1')},failed_control={'root':str(p/'v2')},face_control={'root':str(p/'v3')})
        control=dict(root=str(root),log_sha256=probe.base.sha(case/'log.solver'),manifest_sha256='manifest')
        build_meta=dict(base_build=parent['build'],base_assets=parent['assets'])
        (p/'probe/probe_build_manifest.json').write_text(json.dumps(build_meta))
        identity=dict(binary=str(p/'probe/bin'/probe.APP),sha256='probe',dependencies={})
        return root,parent,control,identity

    def test_prepare_identical_inputs_and_refuse_changed_or_attempted_cases(self):
        with tempfile.TemporaryDirectory() as d,redirect_stdout(io.StringIO()):
            p=Path(d);old,parent,control,identity=self.control(p)
            before=probe.startup.tree_hashes(old)
            with (patch.object(probe,'failed_control',return_value=(parent,control)),
                  patch.object(probe,'verify_build',return_value=identity),patch.object(probe.startup,'require_allocation')):
                out=probe.prepare(old,p/'probe',p/'new');m,_=probe.verify(out)
                self.assertEqual(m['inputs'],parent['inputs'])
                self.assertFalse((out/'quiescent/processor0').exists())
                self.assertFalse((out/'quiescent/log.solver').exists())
                self.assertEqual(probe.startup.tree_hashes(old),before)
                with self.assertRaises(ValueError):probe.prepare(old,p/'probe',out)
                (out/'quiescent/attempt.json').write_text('attempt')
                with self.assertRaisesRegex(ValueError,'attempted'):probe.run(out)
                (out/'quiescent/system/controlDict').write_text('changed')
                with self.assertRaisesRegex(ValueError,'inputs changed'):probe.verify(out)

    def test_received_control_hash_and_variant_required(self):
        with patch.object(probe.startup,'failed_case',return_value=('m','identity')) as call:
            self.assertEqual(probe.failed_control(Path('/control')),('m','identity'))
            self.assertEqual(call.call_args.args[1:3],(probe.startup.COUPLED_TAG,probe.FAILED_LOG))

    def test_monitored_run_bundles_prelimit_failure_and_never_reconstructs(self):
        with tempfile.TemporaryDirectory() as d,redirect_stdout(io.StringIO()):
            p=Path(d);old,parent,control,identity=self.control(p)
            with (patch.object(probe,'failed_control',return_value=(parent,control)),
                  patch.object(probe,'verify_build',return_value=identity),patch.object(probe.startup,'require_allocation')):
                out=probe.prepare(old,p/'probe',p/'new')
                def execute(command,cwd,log,env):
                    self.assertIn('--kill-on-bad-exit=1',command)
                    log.write_text('Time = 5e-05\nSC1_PRESSURE_PROBE version=1 maxCompletedSteps=2\nSC1_PRESSURE_STAGE sequence=1 t_s=5e-05 stage=pressureSolved\nSC1_PRESSURE_STOP reason=unclippedPressureBelow1000Pa min_Pa=-2499398\nFOAM FATAL ERROR\n')
                    raise ValueError('solver failure')
                with (patch.object(probe.base,'execute') as utility,patch.object(watch,'execute',side_effect=execute),
                      patch.object(failure,'collect') as collect):
                    with self.assertRaisesRegex(ValueError,'preserve all'):probe.run(out)
                    self.assertEqual(utility.call_count,1);self.assertEqual(utility.call_args.args[0][0],'decomposePar')
                    collect.assert_called_once_with(out,out.with_name(out.name+'_evidence'))
                result=json.loads((out/'quiescent/result.json').read_text())
                self.assertEqual(result['stops'][0]['min_Pa'],'-2499398')
                self.assertNotIn('startup_pass',result)
                self.assertIn(probe.TAG,failure.VARIANTS)


if __name__=='__main__':unittest.main()
