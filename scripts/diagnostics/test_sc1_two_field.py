"""Local checks for observer placement, private preparation and evidence validation.

No OpenFOAM compilation or CFD is claimed by these tests.
"""
import copy
from contextlib import redirect_stdout
import io
import difflib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import sc1_two_field as harness
from sc1_two_field_cases import CASES, DT, OUTERS, SOURCE_KG_S, STEPS, SUBCYCLES, case_files, foam
from sc1_two_field_instrument import phase_source, pressure_source, solver_source
from sc1_two_field_report import TERMS, STAGES, TIME_KEY_TOLERANCE, parse, summarize, unique
from sc1_wall_mesh import Reader


def records(case):
    """Analytic zero-net-change, finite-transport log for report contract tests."""
    reference=case.endswith('gas_reference')
    rows=[]
    def row(kind, **values): rows.append((kind,values))
    for p in ('gas','liquid'):
        row('SC1_TWO_INITIAL',phase=p,t_s=0,volume_m3=.0005,mass_kg=.0005 if p=='gas' else .5)
    for step in range(STEPS):
        t=round((step+1)*DT,12)
        for p in ('gas','liquid'):
            for outer in range(1,OUTERS+1):
                solved=not(reference and p=='gas')
                row('SC1_PHASE_SOLVE',phase=p,t_s=t,dt_s=DT,call=outer,solved=int(solved),
                    mulesCalls=SUBCYCLES if solved else 0,deltaVolume_m3=0,deltaFrozenRhoMass_kg=0,
                    changeL1Volume_m3=1e-6,unassignedVolume_m3=0,unassignedFrozenRhoMass_kg=0,
                    mulesBoundaryVolume_m3=0,divBoundaryDifference_m3=0)
                for term in TERMS if solved else ('referenceSubcycle','referenceFinal'):
                    row('SC1_PHASE_TERM',phase=p,t_s=t,call=outer,term=term,volume_m3=0,frozenRhoMass_kg=0)
                for stage in STAGES:
                    row('SC1_TWO_STAGE',phase=p,t_s=t,outer=outer,stage=stage,
                        deltaVolume_m3=0,deltaMass_kg=0,alphaPart_kg=0,densityPart_kg=0,identityResidual_kg=0)
            row('SC1_TWO_STAGE',phase=p,t_s=t,outer=OUTERS,stage='endStep',deltaVolume_m3=0,
                deltaMass_kg=0,alphaPart_kg=0,densityPart_kg=0,identityResidual_kg=0)
            feed=SOURCE_KG_S*DT if case.startswith('source_') and p=='liquid' else 0
            row('SC1_TWO_STEP',phase=p,t_s=t,dt_s=DT,outers=OUTERS,volume_m3=.0005,
                mass_kg=.0005 if p=='gas' else .5,deltaMass_kg=0,boundaryOut_kg=feed,
                optionSource_kg=feed,transferSource_kg=0,physicalResidual_kg=0,
                alphaMin=.2,alphaMax=.8,phaseSumMaxError=0)
    return rows


def render(rows):
    return 'SC1_TWO_LIBRARY version=1\n'+'\n'.join(k+' '+' '.join(f'{x}={v}' for x,v in r.items()) for k,r in rows)+'\nEnd\n'


class Instrumentation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.phase=(harness.SOURCE/'phaseSystems/phaseSystem/phaseSystemSolve.C').read_text()
        cls.main=(harness.SOURCE/'multiphaseEulerFoam/multiphaseEulerFoam.C').read_text()
        cls.pressure=(harness.SOURCE/'multiphaseEulerFoam/pEqnComps.H').read_text()

    def test_phase_equations_retained_and_rescale_bracketed(self):
        transformed=phase_source(self.phase)
        removed=[x for x in difflib.ndiff(self.phase.splitlines(),transformed.splitlines()) if x.startswith('- ')]
        self.assertEqual(removed,[])
        rescale=transformed.index('movingPhases()[movingPhasei] *= alphaVoid/sumAlphaMoving;')
        self.assertGreater(rescale,transformed.index('sc1Observer->beforeChange(movingPhases()'))
        self.assertLess(rescale,transformed.index('sc1Observer->afterChange(movingPhases()'))
        self.assertIn('sc1PhaseAudit',transformed)
        self.assertNotIn('sc1AlphaDiagnostic',transformed)
        self.assertIn('sc1AlphaBudget::referenceFinal',transformed)

    def test_source_evaluated_once_with_same_pressure_use(self):
        transformed=pressure_source(self.pressure)
        self.assertEqual(transformed.count('fvModels.source(alpha, rho)'),1)
        self.assertIn('pEqnComp -= sc1OptionSource()/rho;',transformed)
        self.assertIn('sc1Stages.transferSource(phasei, dmdts[phasei]);',transformed)
        self.assertLess(transformed.index('sourceReset'),transformed.index('optionSource'))

    def test_main_has_complete_stage_order_and_no_biology(self):
        transformed=solver_source(self.main)
        self.assertLess(transformed.index('sc1Stages.beginStep'),transformed.index('runTime++;'))
        self.assertLess(transformed.index('sc1Stages.endStep'),transformed.index('runTime.write();'))
        self.assertEqual(transformed.count('sc1Stages.capture("alpha")'),2)
        self.assertNotIn('putida',transformed)
        self.assertIn('sc1TwoFieldAuditLibraryVersion()!=1',transformed)

    def test_incompatible_or_repeat_instrumentation_rejected(self):
        for transform,source in ((phase_source,self.phase),(solver_source,self.main),(pressure_source,self.pressure)):
            with self.assertRaises(ValueError): transform(transform(source))
        with self.assertRaises(ValueError): phase_source(self.phase.replace(' *= alphaVoid/sumAlphaMoving;', ' *= 1;'))


class FieldAssignmentCompatibility(unittest.TestCase):
    """Compile actual initialization statements against OF9's void-assignment contract.

    This targeted regression does not substitute for an OpenFOAM header build.
    """
    @classmethod
    def setUpClass(cls):
        cls.compiler=shutil.which('g++') or shutil.which('c++')
        if not cls.compiler:
            raise unittest.SkipTest('A C++ compiler is required for the assignment regression')
        header=(harness.HERE/'include/sc1PhaseFractionAudit.H').read_text()
        cls.initialization=header.split('sourceAlpha_[i]=entry_[i];\n',1)[1].split('\n        }',1)[0]

    def compile(self, initialization, run=False):
        code=r'''
#include <cassert>
#include <vector>
struct scalarField
{
    std::vector<double> cells;
    scalarField() = default;
    scalarField(unsigned n, double value) : cells(n,value) {}
    void operator=(const scalarField& other) { cells=other.cells; }
};
struct Phase { unsigned size() const { return 3; } };
int main()
{
    const Phase phase;
    std::vector<scalarField> base_(2),dgdt_(2),requested_(2);
    for (unsigned i=0;i<2;++i)
    {
INITIALIZATION
    }
    for (const auto* fields : {&base_,&dgdt_,&requested_})
        for (const auto& field : *fields)
        {
            assert(field.cells.size()==phase.size());
            for (double value : field.cells) assert(value==0);
        }
}
'''.replace('INITIALIZATION',initialization)
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/'assignment.C'; source.write_text(code)
            executable=Path(tmp)/'assignment'
            result=subprocess.run([self.compiler,'-std=c++14','-Wall','-Wextra','-Werror',
                                   str(source),'-o',str(executable)],text=True,capture_output=True)
            if run and result.returncode==0:
                subprocess.run([str(executable)],check=True)
            return result

    def test_observer_initialization_compiles_and_sizes_all_fields(self):
        result=self.compile(self.initialization,run=True)
        self.assertEqual(result.returncode,0,result.stderr)

    def test_regression_reproduces_rejected_chained_assignment(self):
        result=self.compile('base_[i]=dgdt_[i]=requested_[i]=scalarField(phase.size(),0);')
        self.assertNotEqual(result.returncode,0)
        self.assertIn('void',result.stderr)

    def test_new_observer_headers_have_no_chained_assignments(self):
        for name in ('sc1PhaseFractionAudit.H','sc1TwoFieldStages.H'):
            header=(harness.HERE/'include'/name).read_text()
            # Ignore comments/strings; detect assignment to another assignment's result.
            clean=re.sub(r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"','',header,flags=re.S)
            chained=re.search(r'(?<![=!<>])=\s*[A-Za-z_]\w*(?:\[[^\]\n]+\])?\s*=(?!=)',clean)
            self.assertIsNone(chained,name)


class Fixtures(unittest.TestCase):
    @staticmethod
    def alpha_controls(solution, field):
        # Restricted lookup for these flat generated solver entries: literal
        # keys match exactly; only quoted keys are regular expressions.
        # This is a fixture-contract check, not an OpenFOAM parser/runtime test.
        solvers=solution.split('solvers\n{',1)[1].split('\nPIMPLE',1)[0]
        entries=re.findall(r'("[^"\n]+"|\w+)\s*(\{[^{}]*\})',solvers)
        for key,body in entries:
            if key==field or (key.startswith('"') and re.fullmatch(key[1:-1],field)):
                reader=Reader(foam('controls',body).encode())
                controls=reader.dictionary(); reader.finish()
                return controls
        raise KeyError(field)

    def test_alpha_controls_cover_system_and_mules_phase_lookups(self):
        # phaseSystemSolve.C requests "alpha"; the momentum-transfer implicit
        # pressure check and OF9 MULES::limiter request alpha.gas/alpha.liquid,
        # even without an implicit alpha equation.
        for name in CASES:
            with self.subTest(case=name):
                solution=case_files(name)['system/fvSolution']
                controls=[self.alpha_controls(solution,field)
                          for field in ('alpha','alpha.gas','alpha.liquid')]
                self.assertTrue(all(c==controls[0] for c in controls))
                self.assertEqual(controls[0]['nAlphaCorr'],[b'1'])
                self.assertEqual(controls[0]['nAlphaSubCycles'],[str(SUBCYCLES).encode()])
                self.assertEqual(controls[0]['sc1PhaseAudit'],[b'true'])
                for option in ('meanFluxReference','splitPhaseFlux'):
                    self.assertEqual(controls[0][option],[b'false'])

    def test_literal_alpha_key_reproduces_missing_phase_controls(self):
        legacy=case_files(CASES[0])['system/fvSolution'].replace('"alpha.*"','alpha')
        self.assertIn('nAlphaSubCycles',self.alpha_controls(legacy,'alpha'))
        for field in ('alpha.gas','alpha.liquid'):
            with self.assertRaises(KeyError): self.alpha_controls(legacy,field)

    def test_of9_pressure_reference_and_limit_are_in_separate_dictionaries(self):
        # OF9 pressureReference rejects these legacy keys in PIMPLE. Its
        # fvConstraints reader instead loads system/fvConstraints/limitp.
        rejected={'pMin','pMax','pMinFactor','pMaxFactor','rhoMin','rhoMax'}
        for name in CASES:
            with self.subTest(case=name):
                files=case_files(name)
                body=re.search(r'\bPIMPLE\s*(\{[^{}]*\})',files['system/fvSolution']).group(1)
                reader=Reader(foam('fvSolution',body).encode())
                pimple=reader.dictionary(); reader.finish()
                self.assertFalse(rejected.intersection(pimple))
                self.assertEqual(pimple['pRefCell'],[b'0'])
                self.assertEqual(pimple['pRefValue'],[b'100000'])
                reader=Reader(files['system/fvConstraints'].encode())
                reader.expect(b'limitp')
                limit=reader.dictionary(); reader.finish()
                self.assertEqual(limit,{'type':[b'limitPressure'],'min':[b'1000']})
                self.assertNotIn('constant/fvConstraints',files)

    def test_each_pair_only_changes_reference(self):
        for pair in ('closed','source'):
            a=case_files(pair+'_no_reference'); b=case_files(pair+'_gas_reference')
            self.assertEqual(set(a),set(b))
            self.assertEqual([p for p in a if a[p]!=b[p]],['constant/phaseProperties'])
            self.assertEqual(b['constant/phaseProperties'].replace('referencePhase gas;\n',''),a['constant/phaseProperties'])

    def test_known_source_is_vented_closed_box_is_periodic(self):
        closed=case_files(CASES[0]); source=case_files(CASES[2])
        self.assertEqual(closed['system/blockMeshDict'].count('neighbourPatch'),6)
        self.assertNotIn('massSource',closed['constant/fvModels'])
        self.assertIn('selectionMode all',source['constant/fvModels'])
        self.assertIn('xMax { type fixedValue; value uniform 100000;',source['0/p_rgh'])
        self.assertIn('xMax { type pressureInletOutletVelocity;',source['0/U.liquid'])
        for files in (closed,source):
            self.assertNotIn('populationBalance',files['constant/phaseProperties'])
            self.assertNotIn('biology',''.join(files.values()))
            self.assertIn('phaseTransfer ();',files['constant/phaseProperties'])
            self.assertIn('equationOfState rhoConst',files['constant/thermophysicalProperties.gas'])

    def test_reject_unknown_case(self):
        with self.assertRaises(ValueError): case_files('production')


class Reporting(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.baseline=records(CASES[0])

    def test_all_four_complete_contracts(self):
        for case in CASES:
            report=summarize(render(records(case)),case)
            self.assertTrue(report['screen_pass'],report['errors'])

    def test_printed_time_roundoff_preserves_all_four_budgets(self):
        for case in CASES:
            with self.subTest(case=case):
                rows=records(case)
                baseline=summarize(render(rows),case)
                changed=0
                for kind,row in rows:
                    if row['t_s']>0:
                        # Different observers can print slightly different times.
                        row['t_s']+=1e-16 if kind=='SC1_PHASE_TERM' else -1e-16
                        changed+=1
                report=summarize(render(rows),case)
                self.assertTrue(report['screen_pass'],report['errors'])
                self.assertEqual(report['phases'],baseline['phases'])
                self.assertEqual(report['timeKeyNormalization']['adjustedRecords'],changed)
                self.assertLess(report['timeKeyNormalization']['maxAdjustment_s'],TIME_KEY_TOLERANCE)

    def test_time_key_normalization_retains_raw_record(self):
        rows=parse('SC1_TWO_STEP t_s=0.0720000000000001 phase=gas mass_kg=0.0005\n')
        original=copy.deepcopy(rows['SC1_TWO_STEP'])
        indexed=unique(rows['SC1_TWO_STEP'],('t_s','phase'))
        self.assertEqual(rows['SC1_TWO_STEP'],original)
        self.assertEqual(indexed[(0.072,'gas')]['t_s'],0.0720000000000001)

    def test_near_duplicate_time_keys_are_rejected(self):
        rows=copy.deepcopy(self.baseline)
        duplicate=copy.deepcopy(next(r for k,r in rows if k=='SC1_TWO_STEP'))
        duplicate['t_s']+=1e-16
        rows.append(('SC1_TWO_STEP',duplicate))
        with self.assertRaisesRegex(ValueError,'Duplicate record'):
            summarize(render(rows),CASES[0])

    def test_off_grid_and_out_of_range_times_are_rejected_with_evidence(self):
        for time in (0,0.001+4*TIME_KEY_TOLERANCE,0.0015,0.101,1e300):
            with self.subTest(time=time):
                rows=copy.deepcopy(self.baseline)
                next(r for k,r in rows if k=='SC1_TWO_STEP')['t_s']=time
                with self.assertRaisesRegex(ValueError,
                        r'steps: expected=200 observed=200 missing=1 unexpected=1'):
                    summarize(render(rows),CASES[0])

    def test_coverage_failure_names_missing_steps_and_alpha_calls(self):
        rows=[(k,r) for k,r in self.baseline if not(
            k in ('SC1_TWO_STEP','SC1_PHASE_SOLVE') and r['phase']=='gas' and r['t_s']==0.001)]
        with self.assertRaises(ValueError) as caught:
            summarize(render(rows),CASES[0])
        message=str(caught.exception)
        self.assertIn('steps: expected=200 observed=199 missing=1 unexpected=0',message)
        self.assertIn('alpha_calls: expected=400 observed=398 missing=2 unexpected=0',message)
        self.assertIn("missing_examples=[(0.001, 'gas')]",message)

    def test_missing_or_duplicate_records_rejected(self):
        with self.assertRaises(ValueError): summarize(render(self.baseline[:-1]),CASES[0])
        with self.assertRaises(ValueError): summarize(render(self.baseline+[self.baseline[-1]]),CASES[0])

    def test_nonfinite_and_no_end_rejected(self):
        text=render(self.baseline)
        with self.assertRaises(ValueError): summarize(text.replace('mass_kg=0.0005','mass_kg=nan',1),CASES[0])
        with self.assertRaises(ValueError): summarize(text.replace('\nEnd\n',''),CASES[0])

    def mutate(self,kind,key,value):
        rows=copy.deepcopy(self.baseline)
        next(r for k,r in rows if k==kind)[key]=value
        return summarize(render(rows),CASES[0])

    def test_nonconservation_is_failure_despite_solver_end(self):
        report=self.mutate('SC1_TWO_STEP','physicalResidual_kg',1e-5)
        self.assertFalse(report['screen_pass'])
        self.assertIn('gas: phase mass budget exceeds screening tolerance',report['errors'])

    def test_unassigned_alpha_defect_is_failure(self):
        self.assertFalse(self.mutate('SC1_PHASE_SOLVE','unassignedVolume_m3',1e-6)['screen_pass'])

    def test_bounds_and_wrong_source_are_failures(self):
        self.assertFalse(self.mutate('SC1_TWO_STEP','alphaMax',1.01)['screen_pass'])
        self.assertFalse(self.mutate('SC1_TWO_STEP','optionSource_kg',1e-4)['screen_pass'])

    def test_stage_decomposition_and_wrong_outers_fail(self):
        self.assertFalse(self.mutate('SC1_TWO_STAGE','densityPart_kg',.001)['screen_pass'])
        self.assertFalse(self.mutate('SC1_TWO_STEP','outers',4)['screen_pass'])

    def test_cancellation_cannot_hide_mass_error(self):
        rows=copy.deepcopy(self.baseline)
        selected=[r for k,r in rows if k=='SC1_TWO_STEP' and r['phase']=='gas'][:2]
        for r,sign in zip(selected,(1,-1)):
            r['deltaMass_kg']=sign*1e-5; r['physicalResidual_kg']=sign*1e-5
        report=summarize(render(rows),CASES[0])
        self.assertEqual(report['phases']['gas']['physicalResidual_kg'],0)
        self.assertGreater(report['phases']['gas']['residualL1_kg'],1e-5)
        self.assertFalse(report['screen_pass'])

    def test_missing_term_or_outer_fails(self):
        rows=[(k,r) for k,r in self.baseline if not(k=='SC1_PHASE_TERM' and r['term']=='rescale')]
        with self.assertRaises(ValueError): summarize(render(rows),CASES[0])


class ClosedContinuation(unittest.TestCase):
    """Synthetic completed evidence and mocked CFD; no OpenFOAM run is claimed."""
    def seed(self, root):
        manifest=dict(version=1,root=str(root),build=str(root/'private_build'),
                      assets={'private_solver':'verified_hash'},cases={})
        for name in CASES:
            for rel,text in case_files(name).items():
                path=root/name/rel; path.parent.mkdir(parents=True,exist_ok=True)
                path.write_text(text)
            manifest['cases'][name]=harness.case_input_hashes(root/name)
        harness.write_json(root/'cases_manifest.json',manifest)
        first=root/CASES[0]
        harness.write_json(first/'attempt.json',dict(case=CASES[0],
            build=manifest['build'],assets=manifest['assets']))
        mesh=first/'constant/polyMesh'; mesh.mkdir(); (mesh/'points').write_text('mesh sentinel')
        harness.write_json(first/'initialized.json',harness.case_input_hashes(first))
        for name in ('blockMesh','checkMesh','setFields'):
            (first/('log.'+name)).write_text('Mesh OK.\n' if name=='checkMesh' else 'End\n')
        (first/'log.solver').write_text(render(records(CASES[0])))
        return manifest

    @staticmethod
    def fake_execute(command, case, log, env):
        name=Path(command[0]).name
        if name=='blockMesh':
            mesh=case/'constant/polyMesh'; mesh.mkdir(); (mesh/'points').write_text('mesh sentinel')
        elif name not in ('checkMesh','setFields',harness.APP):
            raise AssertionError('Unexpected command: '+str(command))
        log.write_text(render(records(case.name)) if name==harness.APP else
                       ('Mesh OK.\n' if name=='checkMesh' else 'End\n'))

    @staticmethod
    def hashes(root):
        return {str(p.relative_to(root)):harness.sha(p) for p in root.rglob('*') if p.is_file()}

    def test_recovers_report_and_runs_only_fresh_partner(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); m=self.seed(root)
            first_before=self.hashes(root/CASES[0])
            source_before=[self.hashes(root/name) for name in CASES[2:]]
            with patch.object(harness,'verify_build',return_value={'assets':m['assets']}), \
                    patch.object(harness,'execute',side_effect=self.fake_execute) as execute, \
                    redirect_stdout(io.StringIO()):
                harness.continue_closed(root)
            self.assertEqual(len(execute.call_args_list),4)
            self.assertTrue(all(call.args[1]==root/CASES[1] for call in execute.call_args_list))
            self.assertEqual([Path(call.args[0][0]).name for call in execute.call_args_list],
                             ['blockMesh','checkMesh','setFields',harness.APP])
            first_after=self.hashes(root/CASES[0]); first_after.pop('result.json')
            self.assertEqual(first_before,first_after)
            self.assertEqual(source_before,[self.hashes(root/name) for name in CASES[2:]])
            for name in CASES[:2]:
                result=json.loads((root/name/'result.json').read_text())
                self.assertTrue(result['screen_pass'])
                self.assertEqual(result['log_sha256'],harness.sha(root/name/'log.solver'))
            # A second invocation must refuse the now-attempted partner.
            before=self.hashes(root)
            with patch.object(harness,'verify_build',return_value={'assets':m['assets']}), \
                    patch.object(harness,'execute') as execute:
                with self.assertRaisesRegex(ValueError,'attempted/modified'):
                    harness.continue_closed(root)
                execute.assert_not_called()
            self.assertEqual(before,self.hashes(root))

    def test_reuses_matching_existing_report_without_overwriting(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); m=self.seed(root)
            report=harness.completed_report(root,CASES[0],m)
            destination=root/CASES[0]/'result.json'
            # Different formatting makes an accidental rewrite observable.
            destination.write_text(json.dumps(report,separators=(',',':')))
            before=destination.read_bytes()
            with patch.object(harness,'verify_build',return_value={'assets':m['assets']}), \
                    patch.object(harness,'execute',side_effect=self.fake_execute), \
                    redirect_stdout(io.StringIO()):
                harness.continue_closed(root)
            self.assertEqual(destination.read_bytes(),before)

    def test_invalid_evidence_or_partner_blocks_before_any_mutation(self):
        def wrong_attempt(root):
            p=root/CASES[0]/'attempt.json'; data=json.loads(p.read_text())
            data['assets']={}; p.write_text(json.dumps(data))
        def failed_screen(root):
            rows=records(CASES[0])
            next(r for k,r in rows if k=='SC1_TWO_STEP')['physicalResidual_kg']=1e-5
            (root/CASES[0]/'log.solver').write_text(render(rows))
        def changed_recipe_with_updated_snapshot(root):
            (root/CASES[0]/'system/fvSolution').write_text('changed solver controls')
            (root/CASES[0]/'initialized.json').write_text(json.dumps(harness.case_input_hashes(root/CASES[0])))
        mutations={
            'incomplete_log':lambda r:(r/CASES[0]/'log.solver').write_text('SC1_TWO_LIBRARY version=1\nEnd\n'),
            'failed_screen':failed_screen,
            'wrong_attempt_build':wrong_attempt,
            'missing_snapshot':lambda r:(r/CASES[0]/'initialized.json').unlink(),
            'changed_initial_field':lambda r:(r/CASES[0]/'0/alpha.gas').write_text('changed field'),
            'changed_prepared_recipe':changed_recipe_with_updated_snapshot,
            'failed_mesh':lambda r:(r/CASES[0]/'log.checkMesh').write_text('failed mesh'),
            'stale_result':lambda r:(r/CASES[0]/'result.json').write_text('{}'),
            'attempted_partner':lambda r:(r/CASES[1]/'attempt.json').write_text('{}'),
            'changed_partner':lambda r:(r/CASES[1]/'system/fvSolution').write_text('changed'),
        }
        for label,mutate in mutations.items():
            with self.subTest(failure=label), tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp); m=self.seed(root); mutate(root); before=self.hashes(root)
                with patch.object(harness,'verify_build',return_value={'assets':m['assets']}), \
                        patch.object(harness,'execute') as execute:
                    with self.assertRaises((ValueError,OSError)):
                        harness.continue_closed(root)
                    execute.assert_not_called()
                self.assertEqual(self.hashes(root),before)

    def test_changed_build_blocks_continuation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); self.seed(root); before=self.hashes(root)
            with patch.object(harness,'verify_build',return_value={'assets':{'changed':'binary'}}), \
                    patch.object(harness,'execute') as execute:
                with self.assertRaisesRegex(ValueError,'different build'):
                    harness.continue_closed(root)
                execute.assert_not_called()
            self.assertEqual(self.hashes(root),before)

    def test_continuation_requires_explicit_run_and_allocation(self):
        with redirect_stdout(io.StringIO()), patch('sys.stderr',new=io.StringIO()), \
                patch.object(harness,'continue_closed') as launch:
            with self.assertRaises(SystemExit):
                harness.main(['continue-closed','--root','/unused'])
            with patch.dict(os.environ,{},clear=True):
                with self.assertRaisesRegex(ValueError,'active Kestrel allocation'):
                    harness.main(['continue-closed','--root','/unused','--run'])
            launch.assert_not_called()


class Isolation(unittest.TestCase):
    def test_private_preparation_preserves_sources_and_installed_library(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp); installed=tmp/'installed'; installed.mkdir()
            original=installed/'libbirdphaseSystem.so'; original.write_bytes(b'installed sentinel')
            before=harness.selected_hashes(harness.SOURCE)
            out=harness.prepare_build(tmp/'new',installed=installed)
            self.assertEqual(original.read_bytes(),b'installed sentinel')
            self.assertEqual(before,harness.selected_hashes(harness.SOURCE))
            for directory,_ in harness.COMPONENTS:
                text=(out/'src'/directory/'Make/files').read_text()
                self.assertIn(str(out/'lib'),text)
                self.assertNotIn('LIB = $(FOAM_USER_LIBBIN)',text)
            with self.assertRaises(ValueError): harness.prepare_build(out,installed=installed)

    def test_overlap_and_symlink_destinations_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); original=root/'original'; original.mkdir()
            with self.assertRaises(ValueError): harness.fresh(original/'new',[original])
            link=root/'link'; link.symlink_to(root/'missing')
            with self.assertRaises(ValueError): harness.fresh(link,[])

    def test_changed_case_or_attempt_blocks_before_any_execution(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            manifest={'build':str(root/'build'),'cases':{},'assets':{}}
            for name in CASES[:2]:
                case=root/name; case.mkdir()
                for rel,content in case_files(name).items():
                    p=case/rel; p.parent.mkdir(exist_ok=True); p.write_text(content)
                manifest['cases'][name]=harness.case_input_hashes(case)
            (root/CASES[1]/'attempt.json').write_text('{}')
            with patch.object(harness,'verify_cases',return_value=manifest), patch.object(harness,'execute') as execute:
                with self.assertRaises(ValueError): harness.run_pair(root,'closed')
                execute.assert_not_called()
            (root/CASES[1]/'attempt.json').unlink()
            (root/CASES[1]/'0/alpha.liquid').write_text('changed field')
            with patch.object(harness,'verify_cases',return_value=manifest), patch.object(harness,'execute') as execute:
                with self.assertRaises(ValueError): harness.run_pair(root,'closed')
                execute.assert_not_called()

    def test_source_pair_blocked_without_completed_closed_results(self):
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(harness,'verify_cases',return_value={}), patch.object(harness,'execute') as execute:
                with self.assertRaises(FileNotFoundError): harness.run_pair(Path(tmp),'source')
                execute.assert_not_called()

    def test_runtime_verification_rejects_fallback_and_changed_artifacts(self):
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); installed=root/'installed'; installed.mkdir()
            installed_lib=installed/'libbirdphaseSystem.so'; installed_lib.write_bytes(b'original')
            out=harness.prepare_build(root/'build',installed=installed)
            for path in [out/'bin'/harness.APP]+[out/'lib'/n for _,n in harness.COMPONENTS]:
                path.write_bytes(b'fake binary'); path.chmod(0o755)
            def runner(command,**kwargs):
                if command[0]=='nm':
                    kind='U' if Path(command[-1]).name==harness.APP else 'T'
                    return SimpleNamespace(stdout='0000 '+kind+' '+harness.MARKER+'\n',stderr='')
                return SimpleNamespace(stdout='libbirdphaseSystem.so => '+str(out/'lib/libbirdphaseSystem.so')+' (0x1)\n',stderr='')
            with patch.object(harness.subprocess,'run',side_effect=runner):
                verified=harness.verify_build(out,require_record=False)
                harness.write_json(out/'build_verified.json',verified)
                self.assertEqual(harness.verify_build(out),verified)
                def fallback(command,**kwargs):
                    result=runner(command,**kwargs)
                    if command[0]=='ldd': result.stdout=result.stdout.replace(str(out/'lib'),str(installed))
                    return result
                with patch.object(harness.subprocess,'run',side_effect=fallback):
                    with self.assertRaisesRegex(ValueError,'escaped private'): harness.verify_build(out)
                (out/'bin'/harness.APP).write_bytes(b'changed')
                with self.assertRaisesRegex(ValueError,'binaries changed'): harness.verify_build(out)
                (out/'bin'/harness.APP).write_bytes(b'fake binary')
                installed_lib.write_bytes(b'changed installed')
                with self.assertRaisesRegex(ValueError,'Installed binary changed'): harness.verify_build(out)

    def test_failure_output_includes_first_error_before_long_note_tail(self):
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); output=io.StringIO()
            def fail(command,**kwargs):
                kwargs['stdout'].write('observer.H:70: error: assignment from void\n'+
                                       'note: candidate declaration\n'*100)
                return SimpleNamespace(returncode=2)
            with patch.object(harness.subprocess,'run',side_effect=fail), redirect_stdout(output):
                with self.assertRaises(ValueError): harness.execute(['wmake','libso'],root,root/'build.log',{})
            self.assertIn('First compiler errors:',output.getvalue())
            self.assertIn('observer.H:70: error: assignment from void',output.getvalue())
            self.assertTrue((root/'build.log').is_file())

    def test_failed_solver_is_not_labelled_a_compiler_error(self):
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);output=io.StringIO()
            def fail(command,**kwargs):
                kwargs['stdout'].write('srun: error: task 2: Floating point exception\n')
                return SimpleNamespace(returncode=143)
            with patch.object(harness.subprocess,'run',side_effect=fail),redirect_stdout(output):
                with self.assertRaises(ValueError):harness.execute(['srun','solver'],root,root/'log.solver',{})
            self.assertIn('First command errors:',output.getvalue())
            self.assertNotIn('compiler',output.getvalue())

    def test_runtime_environment_does_not_change_parent(self):
        before=dict(os.environ)
        env=harness.environment(Path('/scratch/fresh'))
        self.assertEqual(dict(os.environ),before)
        self.assertEqual(env['FOAM_USER_LIBBIN'],'/scratch/fresh/lib')
        self.assertTrue(env['LD_LIBRARY_PATH'].startswith('/scratch/fresh/lib'))


if __name__=='__main__': unittest.main()
