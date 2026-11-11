"""Passive prefix, face reconstruction, immutable build/case and failure evidence."""
from contextlib import ExitStack,redirect_stdout
import difflib
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import sc1_reactor_vof_force as probe
import sc1_reactor_vof_force_build as build
import sc1_reactor_vof_force_report as reporting
import sc1_reactor_vof_force_instrument as instrument
import sc1_reactor_watch as watch
import test_sc1_reactor_vof_nonorth as previous


def observed_log():
    control=previous.corrected_log().split('Time = 0.00015\n')[0]+'End\n'
    lines=['SC1_VOF_FORCE_PROBE version=1 maxSteps=2'];sequence=0;pressures=0
    normals=np.array([[.8,.6,0],[-.6,.8,0],[0,0,1]])
    normals=np.concatenate((normals,-normals))
    def emit(kind,row):lines.append('SC1_VOF_FORCE_'+kind+' '+' '.join(f'{k}={v}' for k,v in row.items()))
    for line in control.splitlines():
        lines.append(line)
        if not line.startswith('SC1_VOF_PRESSURE '):continue
        pressures+=1
        if pressures%10:continue
        sequence+=1;t=((sequence-1)//4+1)*5e-5
        emit('CHECK',dict(sequence=sequence,t_s=t,forceClosure_m_s=0,updateClosure_m_s=0,nonfinite=0))
        for phase,a in (('gas',0),('liquid',1)):
            rho=1+999*a;rau=5e-5/rho
            terms=dict(U=np.array([.0001,0,0]),H=np.zeros(3),gravity=np.array([.25,-.1,.05]),capillary=np.array([1e-5,2e-5,-3e-5]),net=np.array([.0001,0,0]))
            terms['pressure']=terms['net']-terms['gravity']-terms['capillary']
            terms['freshPressure']=terms['pressure'];terms['geometryDelta']=-.1*terms['gravity']
            cell=dict(sequence=sequence,t_s=t,phase=phase,proc=0,localCell=a,x_m=.1,y_m=1,z_m=.2,alphaLiquid=a,rho_kg_m3=rho,rAU=rau,p_rgh_Pa=100000,gh_m2_s2=-10)
            for key,value in terms.items():
                cell.update({key+k:float(x) for k,x in zip('xyz' if key in ('U','H') else 'XYZ',value)})
            emit('CELL',cell)
            for i,n in enumerate(normals):
                f=dict(sequence=sequence,phase=phase,proc=0,localFace=i,outwardSign=1,patch=-1,coupled=0,
                    ghFace=-10,ghInterpolated=-9,rAUf=rau)
                for name in ('Sf','Cf','dC'):f.update({name+k:float(x) for k,x in zip('XYZ',n)})
                for name in reporting.VECTORS[2:]:f[name+'_Pa_m']=float(np.dot(terms[name]/rau,n))
                f['snGradRho']=f['gravity_Pa_m']/10
                emit('FACE',f)
            emit('FACES',dict(sequence=sequence,phase=phase,count=6))
    return '\n'.join(lines)+'\n',control


class Reports(unittest.TestCase):
    def test_independent_rotated_face_reconstruction_and_passive_prefix(self):
        text,control=observed_log();r=reporting.analyze(text,control)
        self.assertTrue(r['diagnostic_complete']);self.assertFalse(r['accepted_startup'])
        self.assertEqual(len(r['faces']),96);self.assertEqual(r['prefixAgreement']['records'],83)
        self.assertAlmostEqual(r['cells'][0]['magnitudes_m_s']['net'],.0001)

    def test_prefix_change_or_incomplete_run_rejected(self):
        text,control=observed_log()
        for t in (text.replace('End\n',''),text.replace('gasMax_m_s=.0001','gasMax_m_s=.0002',1),
            text.replace('SC1_VOF_FORCE_PROBE version=1 maxSteps=2\n',''),text+'Time = 0.00015\n'):
            with self.assertRaises(ValueError):reporting.analyze(t,control)

    def test_missing_duplicate_nonfinite_and_inconsistent_force_data_rejected(self):
        text,control=observed_log();face=next(l for l in text.splitlines() if l.startswith('SC1_VOF_FORCE_FACE '))
        for t in (text.replace(face+'\n','',1),text.replace(face,face+'\n'+face,1),
            text.replace('forceClosure_m_s=0','forceClosure_m_s=0.1',1),
            text.replace('rAU=5e-05','rAU=nan',1),text.replace('ghFace=-10','ghFace=-20',1),
            text.replace('netX=0.0001','netX=0.001',1),text.replace('count=6','count=7',1)):
            with self.assertRaises(ValueError):reporting.analyze(t,control)


class Build(unittest.TestCase):
    def test_only_observer_additions_and_private_application_destination(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);project=p/'OpenFOAM-9';base=p/'base'
            for rel,text in json.loads(build.reference.FIXTURE.read_text())['sources'].items():
                q=project/rel;q.parent.mkdir(parents=True,exist_ok=True);q.write_text(text)
            source=build.reference.private_sources(project,base)
            for rel,text in source.items():
                q=base/'src'/rel;q.parent.mkdir(parents=True,exist_ok=True);q.write_text(text)
            actual=build.private_sources(base,p/'observer')
            changed={k for k in source if actual[k]!=source[k]}
            self.assertEqual(changed,{'interFoam.C','pEqn.H','Make/files'})
            for name in ('interFoam.C','pEqn.H'):
                self.assertFalse([l for l in difflib.ndiff(source[name].splitlines(),actual[name].splitlines()) if l.startswith('- ')])
            self.assertIn('sc1InterFoamForces',actual['Make/files'])
            self.assertIn('runTime.value()>0.000100000001',actual['interFoam.C'])
            self.assertIn('sequence_>8',actual['sc1VofForceAudit.H'])
            self.assertEqual(actual['sc1VofAudit.H'],source['sc1VofAudit.H'])
            with self.assertRaises(ValueError):instrument.pressure(source['pEqn.H']*2)
            before=probe.base.source_hashes(base/'src')
            old=dict(binary=str(base/'bin/sc1InterFoamReference'),assets={'library':'hash'},dependencies={'lib':'/lib'})
            new=dict(binary=str(p/'observer/bin/sc1InterFoamForces'),sha256='executable',assets={'library':'hash'},dependencies=old['dependencies'])
            with (patch.object(build.reference,'verify',return_value=old),patch.object(build,'identity',return_value=new),
                patch.object(build.startup,'require_allocation'),patch.object(build.shutil,'which',return_value='wmake'),
                patch.dict(os.environ,WM_PROJECT_VERSION='9',WM_PROJECT_DIR=str(project)),
                patch.object(build.base,'execute') as execute,redirect_stdout(io.StringIO())):
                out=build.build(base,p/'observer');self.assertEqual(build.verify(out),new)
                execute.assert_called_once();self.assertEqual(before,probe.base.source_hashes(base/'src'))
                (out/'src/pEqn.H').write_text('changed')
                with self.assertRaises(ValueError):build.verify(out)


class Cases(unittest.TestCase):
    def fixture(self,p):
        root,spatial,m=previous.Cases().fixture(p);case=root/'quiescent'
        sol=case/'system/fvSolution';sol.write_text(probe.nonorth.solution(sol.read_text()))
        log=previous.corrected_log().split('Time = 0.00845\n')[0]+'FOAM FATAL ERROR:\n'
        (case/'log.solver').write_text(log)
        m.update(variant=probe.nonorth.TAG,inputs=probe.base.case_input_hashes(case),nonorth_recipe=probe.nonorth.recipe(),nonorth_control={'root':str(p/'gradient')})
        (root/'startup_manifest.json').write_text(json.dumps(m))
        out=p/'force_build';out.mkdir();(out/'force_build_manifest.json').write_text(json.dumps({'base_build':m['build']}))
        binary=dict(binary=str(out/'bin/sc1InterFoamForces'),assets={'force solver':'hash'})
        return root,m,out,binary

    def patches(self,root,m,binary):
        stack=ExitStack();stack.enter_context(patch.object(probe.nonorth,'verify',return_value=m))
        stack.enter_context(patch.object(probe,'CONTROL_LOG',probe.base.sha(root/'quiescent/log.solver')))
        stack.enter_context(patch.object(build,'verify',return_value=binary))
        stack.enter_context(patch.object(probe.startup,'require_allocation'));stack.enter_context(redirect_stdout(io.StringIO()))
        return stack

    def test_fresh_case_only_changes_application_duration_and_preserves_all_controls(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);root,m,binary_root,binary=self.fixture(p)
            with self.patches(root,m,binary):
                before=probe.startup.tree_hashes(root);out=probe.prepare(root,binary_root,p/'probe');n=probe.verify(out)
                self.assertEqual([k for k in m['inputs'] if m['inputs'][k]!=n['inputs'][k]],['system/controlDict'])
                self.assertEqual(probe.startup.tree_hashes(root),before)
                self.assertFalse((out/'quiescent/processor0').exists())
                with self.assertRaises(ValueError):probe.prepare(root,binary_root,out)
                with self.assertRaises(ValueError):probe.prepare(root,binary_root,root/'nested')
                (out/'quiescent/0/p_rgh').write_text('changed')
                with self.assertRaises(ValueError):probe.verify(out)

    def test_run_collects_full_force_evidence_without_reconstruction_and_rejects_reuse(self):
        for failed in (False,True):
            with self.subTest(failed=failed),tempfile.TemporaryDirectory() as d:
                p=Path(d);root,m,binary_root,binary=self.fixture(p)
                with self.patches(root,m,binary):
                    out=probe.prepare(root,binary_root,p/'probe');before=probe.startup.tree_hashes(root)
                    def run(command,case,log,env,wall_seconds):
                        self.assertEqual(wall_seconds,600);self.assertIn(binary['binary'],command)
                        log.write_text(observed_log()[0])
                        if failed:raise ValueError('controlled stop')
                    with patch.object(watch,'execute',side_effect=run),patch.object(probe.base,'execute') as shell:
                        if failed:
                            with self.assertRaises(ValueError):probe.run(out)
                        else:self.assertTrue(probe.run(out)['diagnostic_complete'])
                        shell.assert_called_once();self.assertEqual(shell.call_args.args[0][0],'decomposePar')
                    self.assertTrue((p/'probe_evidence/reactor_failure_evidence.tar.gz').is_file())
                    self.assertEqual(probe.startup.tree_hashes(root),before)
                    with self.assertRaisesRegex(ValueError,'attempted'):probe.run(out)


if __name__=='__main__':unittest.main()
