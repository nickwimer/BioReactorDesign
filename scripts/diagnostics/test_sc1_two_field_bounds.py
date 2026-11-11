"""Near-bound fixture, report and isolated runner tests. CFD is mocked."""
import copy
from contextlib import redirect_stdout
import hashlib
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import sc1_two_field as base
import sc1_two_field_bounds as runner
from sc1_two_field_bounds_cases import ALPHA, CASES, LEVELS, case_files, contract
from sc1_two_field_bounds_report import summarize
from sc1_two_field_cases import CASES as BASELINE_CASES, case_files as baseline_files, OUTERS, SUBCYCLES, SOURCE_KG_S
from sc1_two_field_report import TERMS, STAGES
import test_sc1_two_field as baseline_tests
from test_sc1_two_field import records as baseline_records, render


def records(case):
    """Analytic budget/identity fixture; its split and raw terms cancel exactly."""
    spec=contract(case); rows=[]
    reference=case.endswith('gas_reference')
    for i,p in enumerate(('gas','liquid')):
        density=1 if p=='gas' else 1000
        volume=spec.volumes[i]; mass=volume*density
        rows.append(('SC1_TWO_INITIAL',dict(phase=p,t_s=0,volume_m3=volume,mass_kg=mass)))
        for step in range(spec.steps):
            t=round((step+1)*spec.dt,12)
            for outer in range(1,OUTERS+1):
                solved=not (reference and p=='gas')
                rows.append(('SC1_PHASE_SOLVE',dict(phase=p,t_s=t,dt_s=spec.dt,call=outer,
                    solved=int(solved),mulesCalls=SUBCYCLES if solved else 0,
                    deltaVolume_m3=0,deltaFrozenRhoMass_kg=0,changeL1Volume_m3=0,
                    unassignedVolume_m3=0,unassignedFrozenRhoMass_kg=0,
                    mulesBoundaryVolume_m3=0,divBoundaryDifference_m3=0)))
                for term in TERMS if solved else ('referenceSubcycle','referenceFinal'):
                    value=0
                    if p=='liquid' and term in ('rawDgdt','splitAtSourceAlpha'):
                        value=spec.dt*1e-11*(1 if term=='rawDgdt' else -1)
                    rows.append(('SC1_PHASE_TERM',dict(phase=p,t_s=t,call=outer,term=term,
                        volume_m3=value,frozenRhoMass_kg=value*density)))
                for stage in STAGES:
                    rows.append(('SC1_TWO_STAGE',dict(phase=p,t_s=t,outer=outer,stage=stage,
                        deltaVolume_m3=0,deltaMass_kg=0,alphaPart_kg=0,densityPart_kg=0,identityResidual_kg=0)))
            rows.append(('SC1_TWO_STAGE',dict(phase=p,t_s=t,outer=OUTERS,stage='endStep',
                deltaVolume_m3=0,deltaMass_kg=0,alphaPart_kg=0,densityPart_kg=0,identityResidual_kg=0)))
            supply=SOURCE_KG_S*spec.dt if p=='liquid' else 0
            rows.append(('SC1_TWO_STEP',dict(phase=p,t_s=t,dt_s=spec.dt,outers=OUTERS,
                volume_m3=volume,mass_kg=mass,deltaMass_kg=0,boundaryOut_kg=supply,
                optionSource_kg=supply,transferSource_kg=0,physicalResidual_kg=0,
                alphaMin=float(ALPHA[p]),alphaMax=float(ALPHA[p]),phaseSumMaxError=0)))
    return rows


class Fixtures(unittest.TestCase):
    def test_reference_and_timestep_comparisons_are_controlled(self):
        for level in LEVELS:
            a=case_files(f'bound_{level}_no_reference');b=case_files(f'bound_{level}_gas_reference')
            self.assertEqual([p for p in a if a[p]!=b[p]],['constant/phaseProperties'])
            self.assertEqual(a['constant/phaseProperties'],b['constant/phaseProperties'].replace('referencePhase gas;\n',''))
        for ref in ('no_reference','gas_reference'):
            a=case_files('bound_coarse_'+ref);b=case_files('bound_fine_'+ref)
            self.assertEqual([p for p in a if a[p]!=b[p]],['system/controlDict'])
            expected=a['system/controlDict'].replace('deltaT 0.001;','deltaT 0.0005;').replace(
                'maxDeltaT 0.001;','maxDeltaT 0.0005;').replace('writeInterval 100;','writeInterval 200;')
            self.assertEqual(expected,b['system/controlDict'])

    def test_initial_and_outlet_fractions_match_without_other_physics_changes(self):
        for name in CASES:
            files=case_files(name); spec=contract(name)
            base_files=baseline_files('source_gas_reference' if name.endswith('gas_reference') else 'source_no_reference')
            self.assertEqual({p for p in files if files[p]!=base_files[p]},
                {'0/alpha.gas','0/alpha.liquid','system/setFieldsDict'} |
                ({'system/controlDict'} if name.startswith('bound_fine_') else set()))
            for phase,alpha in ALPHA.items():
                self.assertEqual(files['0/alpha.'+phase].count('uniform '+alpha+';'),3)
                self.assertIn('alpha.'+phase+' '+alpha,files['system/setFieldsDict'])
            self.assertEqual(spec.dt*spec.steps,.1)
            self.assertAlmostEqual(sum(spec.volumes),.001)
            self.assertEqual(spec.volumes[0],1e-8)
            self.assertAlmostEqual(spec.volumes[1]*1000,.99999)

    def test_original_report_contracts_are_unchanged(self):
        # Canonical JSON digests captured from the committed baseline analyzer
        # before extracting its configurable inventory/timestep contract.
        digests=('599f3ff92e0d41d025cc223f8074c8d9d6c93c71fef5cc6dad1a1526b81ae7c8',
                 'b9ccbfa7c7f6c02badce9fc269aee9b844e7ecb64eb33c4bf99724c3046c50a5',
                 '9f2ed23f6b9ea415aa3c60133f189684a54d259045c3c826133305ed7f4176c2',
                 '305e9d45ef02c4ed266b8a6271a3aa0575dc7869cbb98787dda3f710d7ccd8fe')
        for name,digest in zip(BASELINE_CASES,digests):
            report=base.summarize(render(baseline_records(name)),name)
            self.assertEqual(hashlib.sha256(json.dumps(report,sort_keys=True).encode()).hexdigest(),digest)


class Reporting(unittest.TestCase):
    def test_both_timestep_levels_cover_full_same_duration(self):
        for name in CASES:
            report=summarize(render(records(name)),name)
            self.assertTrue(report['screen_pass'],report['errors'])
            self.assertEqual(report['steps'],contract(name).steps)
            self.assertEqual(report['endTime_s'],.1)
            self.assertTrue(report['sourceSplitObserved'])
            self.assertEqual(report['phases']['gas']['nearBoundSteps'],contract(name).steps)
            self.assertEqual(report['phases']['gas']['screenTolerance_kg'],1e-15)
            self.assertAlmostEqual(report['phases']['liquid']['screenTolerance_kg'],1.01e-10)

    def test_missing_duplicate_off_grid_and_wrong_case_records_rejected(self):
        name=CASES[2]; original=records(name)
        bad=[original[:-1],original+[original[-1]],copy.deepcopy(original)]
        next(r for k,r in bad[-1] if k=='SC1_TWO_STEP')['t_s']+=.00001
        for rows in bad:
            with self.assertRaises(ValueError): summarize(render(rows),name)
        with self.assertRaises(ValueError): summarize(render(original),CASES[0])
        with self.assertRaises(ValueError): summarize(render(original),'source_no_reference')

    def test_fine_time_roundoff_and_canonical_duplicates(self):
        name=CASES[2]; rows=records(name)
        for _,r in rows:
            if r['t_s']>0:r['t_s']+=1e-16
        self.assertTrue(summarize(render(rows),name)['screen_pass'])
        duplicate=copy.deepcopy(next(r for k,r in rows if k=='SC1_TWO_STEP'))
        duplicate['t_s']-=1e-16
        with self.assertRaisesRegex(ValueError,'Duplicate'):
            summarize(render(rows+[('SC1_TWO_STEP',duplicate)]),name)

    def test_wrong_inventory_source_and_unexercised_regime_fail(self):
        changes=(('SC1_TWO_INITIAL','gas','mass_kg',.0005),
                 ('SC1_TWO_STEP','liquid','optionSource_kg',0),
                 ('SC1_TWO_STEP','gas','transferSource_kg',1e-16),
                 ('SC1_TWO_STEP','gas','alphaMax',.001),
                 ('SC1_TWO_STEP','liquid','alphaMin',.5),
                 ('SC1_TWO_STEP','gas','alphaMin',0))
        for kind,phase,key,value in changes:
            with self.subTest(key=key,phase=phase):
                rows=records(CASES[0]);next(r for k,r in rows if k==kind and r['phase']==phase)[key]=value
                self.assertFalse(summarize(render(rows),CASES[0])['screen_pass'])
        rows=records(CASES[0])
        for k,r in rows:
            if k=='SC1_PHASE_TERM':r['volume_m3']=r['frozenRhoMass_kg']=0
        result=summarize(render(rows),CASES[0])
        self.assertFalse(result['sourceSplitObserved']);self.assertFalse(result['screen_pass'])

    def test_small_gas_mass_errors_cannot_use_baseline_absolute_allowance(self):
        rows=records(CASES[0])
        gas=next(r for k,r in rows if k=='SC1_TWO_STEP' and r['phase']=='gas')
        gas['boundaryOut_kg']=gas['physicalResidual_kg']=1e-13
        result=summarize(render(rows),CASES[0])
        self.assertFalse(result['screen_pass'])
        self.assertIn('gas: phase mass budget exceeds screening tolerance',result['errors'])
        self.assertGreater(result['phases']['gas']['residualL1_kg'],result['phases']['gas']['screenTolerance_kg'])


class Isolation(unittest.TestCase):
    def seed_baseline(self, root):
        root.mkdir(); m=baseline_tests.ClosedContinuation().seed(root)
        first=base.completed_report(root,BASELINE_CASES[0],m)
        base.write_json(root/BASELINE_CASES[0]/'result.json',first)
        with patch.object(base,'execute',side_effect=baseline_tests.ClosedContinuation.fake_execute),redirect_stdout(io.StringIO()):
            base.run_fresh_cases(root,m,BASELINE_CASES[1:])
        return m

    @staticmethod
    def fake_execute(command,case,log,env):
        if Path(command[0]).name==base.APP:
            log.write_text(render(records(case.name)))
        else: baseline_tests.ClosedContinuation.fake_execute(command,case,log,env)

    def test_complete_lifecycle_preserves_baseline_and_uses_existing_build(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); baseline=root/'baseline';m=self.seed_baseline(baseline)
            before=baseline_tests.ClosedContinuation.hashes(baseline)
            with patch.object(base,'verify_build',return_value={'assets':m['assets']}),redirect_stdout(io.StringIO()):
                out=runner.prepare(baseline,root/'near')
                with patch.object(base,'execute') as execute:
                    with self.assertRaisesRegex(ValueError,'completed saved'):
                        runner.run(out,'fine')
                    execute.assert_not_called()
                with patch.object(base,'execute',side_effect=self.fake_execute) as execute:
                    runner.run(out,'coarse');runner.run(out,'fine')
                self.assertEqual(len(execute.call_args_list),16)
                self.assertEqual({c.args[1].name for c in execute.call_args_list},set(CASES))
                solver_calls=[c for c in execute.call_args_list if Path(c.args[0][0]).name==base.APP]
                self.assertTrue(all(c.args[0][0]==Path(m['build'])/'bin'/base.APP for c in solver_calls))
                out_before=baseline_tests.ClosedContinuation.hashes(out)
                result=runner.compare(out)
                self.assertTrue(result['complete'])
                self.assertTrue(all(r['screen_pass'] for r in result['cases'].values()))
                self.assertEqual(out_before,baseline_tests.ClosedContinuation.hashes(out))
                with patch.object(base,'execute') as execute:
                    with self.assertRaisesRegex(ValueError,'attempted/modified'):runner.run(out,'coarse')
                    execute.assert_not_called()
            self.assertEqual(before,baseline_tests.ClosedContinuation.hashes(baseline))

    def test_fine_accepts_complete_failed_coarse_screen_without_overwriting_it(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);m=self.seed_baseline(root/'baseline')
            def failed_coarse(command,case,log,env):
                if Path(command[0]).name==base.APP:
                    rows=records(case.name)
                    step=next(r for k,r in rows if k=='SC1_TWO_STEP' and r['phase']=='gas')
                    step['boundaryOut_kg']=step['physicalResidual_kg']=1e-13
                    log.write_text(render(rows))
                else:self.fake_execute(command,case,log,env)
            with patch.object(base,'verify_build',return_value={'assets':m['assets']}),redirect_stdout(io.StringIO()):
                out=runner.prepare(root/'baseline',root/'near')
                with patch.object(base,'execute',side_effect=failed_coarse):runner.run(out,'coarse')
                coarse_before={n:baseline_tests.ClosedContinuation.hashes(out/n) for n in CASES[:2]}
                self.assertFalse(json.loads((out/CASES[0]/'result.json').read_text())['screen_pass'])
                with patch.object(base,'execute',side_effect=self.fake_execute):runner.run(out,'fine')
                self.assertEqual(coarse_before,{n:baseline_tests.ClosedContinuation.hashes(out/n) for n in CASES[:2]})

    def test_fine_blocks_coarse_regime_or_accounting_failures(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);m=self.seed_baseline(root/'baseline')
            def bad_coarse(command,case,log,env):
                if Path(command[0]).name==base.APP:
                    rows=records(case.name)
                    next(r for k,r in rows if k=='SC1_TWO_STEP' and r['phase']=='gas')['alphaMax']=.5
                    log.write_text(render(rows))
                else:self.fake_execute(command,case,log,env)
            with patch.object(base,'verify_build',return_value={'assets':m['assets']}),redirect_stdout(io.StringIO()):
                out=runner.prepare(root/'baseline',root/'near')
                with patch.object(base,'execute',side_effect=bad_coarse):runner.run(out,'coarse')
                before=baseline_tests.ClosedContinuation.hashes(out)
                with patch.object(base,'execute') as execute:
                    with self.assertRaisesRegex(ValueError,'regime/accounting'):
                        runner.run(out,'fine')
                    execute.assert_not_called()
                self.assertEqual(before,baseline_tests.ClosedContinuation.hashes(out))

    def test_modified_partner_and_changed_baseline_block_before_execution(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);m=self.seed_baseline(root/'baseline')
            with patch.object(base,'verify_build',return_value={'assets':m['assets']}),redirect_stdout(io.StringIO()):
                out=runner.prepare(root/'baseline',root/'near')
                target=out/CASES[1]/'system/fvSolution'; original=target.read_text();target.write_text('changed')
                before=baseline_tests.ClosedContinuation.hashes(out)
                with patch.object(base,'execute') as execute:
                    with self.assertRaisesRegex(ValueError,'inputs changed'):runner.run(out,'coarse')
                    execute.assert_not_called()
                self.assertEqual(before,baseline_tests.ClosedContinuation.hashes(out));target.write_text(original)
                (root/'baseline'/BASELINE_CASES[0]/'log.solver').write_text('changed')
                with patch.object(base,'execute') as execute:
                    with self.assertRaises(ValueError):runner.run(out,'coarse')
                    execute.assert_not_called()

    def test_prepare_refuses_unfinished_baseline_and_existing_destination(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);m=self.seed_baseline(root/'baseline')
            with patch.object(base,'verify_build',return_value={'assets':m['assets']}),redirect_stdout(io.StringIO()):
                runner.prepare(root/'baseline',root/'near')
                with self.assertRaisesRegex(ValueError,'already exists'):runner.prepare(root/'baseline',root/'near')
                (root/'baseline'/BASELINE_CASES[3]/'result.json').unlink()
                with self.assertRaisesRegex(ValueError,'Missing saved baseline'):runner.prepare(root/'baseline',root/'other')
                self.assertFalse((root/'other').exists())

    def test_run_requires_explicit_flag_and_allocation(self):
        with patch('sys.stderr',new=io.StringIO()),patch.object(runner,'run') as execute:
            with self.assertRaises(SystemExit):runner.main(['run','--root','/unused','--level','coarse'])
            with patch.dict(os.environ,{},clear=True):
                with self.assertRaisesRegex(ValueError,'allocation'):
                    runner.main(['run','--root','/unused','--level','coarse','--run'])
            execute.assert_not_called()


if __name__=='__main__': unittest.main()
