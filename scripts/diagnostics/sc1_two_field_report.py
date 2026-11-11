"""Validate complete observer streams; never mistake process success for conservation."""
from collections import defaultdict
from dataclasses import dataclass
import math
from sc1_two_field_cases import CASES, DT, OUTERS, SOURCE_KG_S, STEPS, SUBCYCLES

TERMS = ('mulesChange oldTimeReset explicitSu implicitSp divergence equationResidual '
         'equationResidualL1 originalSu rawDgdt pressureDilatation couplingDifference '
         'splitAtSourceAlpha implicitSourceChange sourceDecompositionRemainder '
         'implicitCorrection referenceSubcycle clamp referenceFinal implicitDivergence '
         'implicitEquationResidual implicitEquationResidualL1 finalTransportDivergence '
         'fluxAveragingDefect rescale').split()
STAGES = ('alpha', 'modelCorrect', 'transportPressure', 'endOuter')
KINDS = ('SC1_TWO_INITIAL', 'SC1_PHASE_SOLVE', 'SC1_PHASE_TERM', 'SC1_TWO_STAGE', 'SC1_TWO_STEP')
TEXT_KEYS = {'phase', 'term', 'stage'}
# Match only printed-time roundoff around the contract's prescribed grid.
# This is not permission to omit or merge records.
TIME_KEY_TOLERANCE = 1e-12


@dataclass(frozen=True)
class AuditContract:
    """Explicit fixed-mesh two-phase reporting contract; densities are 1/1000."""
    dt: float = DT
    steps: int = STEPS
    volumes: tuple = (.0005, .0005)
    volume_initial_tol: tuple = (1e-11, 1e-11)
    mass_initial_tol: tuple = (1e-8, 1e-8)
    volume_identity_tol: tuple = (1e-13, 1e-13)
    mass_identity_tol: tuple = (1e-11, 1e-11)
    frozen_mass_tol: tuple = (1e-10, 1e-10)
    equation_l1_tol: tuple = (1e-12, 1e-12)
    budget_arithmetic_tol: tuple = (1e-12, 1e-12)
    budget_absolute: tuple = (1e-10, 1e-10)
    budget_relative: tuple = (1e-3, 1e-3)
    option_source_kg_s: tuple = (0, SOURCE_KG_S)


BASELINE_CONTRACT = AuditContract()


def parse(text):
    rows = defaultdict(list)
    for line in text.splitlines():
        if not line.startswith(KINDS):
            continue
        tokens = line.split()
        if tokens[0] not in KINDS:
            raise ValueError('Unknown audit record: '+tokens[0])
        r = {}
        for token in tokens[1:]:
            if '=' not in token:
                raise ValueError('Malformed audit token: '+token)
            key, value = token.split('=', 1)
            if key in r:
                raise ValueError('Duplicate audit key: '+key)
            if key not in TEXT_KEYS:
                value = float(value)
                if not math.isfinite(value):
                    raise ValueError('Nonfinite audit value')
            r[key] = value
        rows[tokens[0]].append(r)
    return rows


def time_key(value, contract=BASELINE_CONTRACT):
    """Canonicalize near-grid times for indexing; retain all raw record values."""
    if contract.dt-TIME_KEY_TOLERANCE <= value <= contract.steps*contract.dt+TIME_KEY_TOLERANCE:
        step = round(value/contract.dt)
        target = round(step*contract.dt, 12)
        if 1 <= step <= contract.steps and abs(value-target) <= TIME_KEY_TOLERANCE:
            return target
    return value


def unique(rows, fields, contract=BASELINE_CONTRACT):
    result = {}
    for row in rows:
        key = tuple(time_key(row[k],contract) if k=='t_s' else row[k] for k in fields)
        if key in result:
            raise ValueError('Duplicate record: '+repr(key))
        result[key] = row
    return result


def require_coverage(groups):
    failures = []
    for name, observed, expected in groups:
        observed = set(observed)
        missing, unexpected = sorted(expected-observed), sorted(observed-expected)
        if missing or unexpected:
            failures.append(
                f'{name}: expected={len(expected)} observed={len(observed)} '
                f'missing={len(missing)} unexpected={len(unexpected)} '
                f'missing_examples={missing[:3]!r} unexpected_examples={unexpected[:3]!r}')
    if failures:
        raise ValueError('Incomplete or unexpected audit records; '+'; '.join(failures))


def close(a, b, tol=1e-12):
    return abs(a-b) <= tol


def summarize(text, case):
    if case not in CASES:
        raise ValueError('Unknown case')
    return summarize_contract(text,case,BASELINE_CONTRACT,
        closed=case.startswith('closed_'),reference=case.endswith('gas_reference'))


def summarize_contract(text, case, contract, *, closed, reference):
    if text.count('SC1_TWO_LIBRARY version=1') != 1 or '\nEnd\n' not in text:
        raise ValueError('Missing unique library marker or successful solver End')
    if any(s in text for s in ('FOAM FATAL', 'Floating point exception', 'Segmentation fault')):
        raise ValueError('Solver failure in log')
    rows = parse(text)
    initial = unique(rows['SC1_TWO_INITIAL'], ('phase',),contract)
    steps = unique(rows['SC1_TWO_STEP'], ('t_s', 'phase'),contract)
    solves = unique(rows['SC1_PHASE_SOLVE'], ('t_s', 'phase', 'call'),contract)
    terms = unique(rows['SC1_PHASE_TERM'], ('t_s', 'phase', 'call', 'term'),contract)
    stages = unique(rows['SC1_TWO_STAGE'], ('t_s', 'phase', 'outer', 'stage'),contract)
    times = [round((i+1)*contract.dt, 12) for i in range(contract.steps)]
    phases = ('gas', 'liquid')
    phase_index = {p:i for i,p in enumerate(phases)}
    expected_steps = {(t,p) for t in times for p in phases}
    expected_solves = {(t,p,o) for t in times for p in phases for o in range(1,OUTERS+1)}
    expected_terms = {(t,p,o,j) for t,p,o in expected_solves
                      for j in (('referenceSubcycle','referenceFinal') if reference and p=='gas' else TERMS)}
    expected_stages = {(t,p,o,s) for t,p,o in expected_solves for s in STAGES}
    expected_stages |= {(t,p,OUTERS,'endStep') for t in times for p in phases}
    require_coverage((
        ('initial', initial, {(p,) for p in phases}),
        ('steps', steps, expected_steps),
        ('alpha_calls', solves, expected_solves),
        ('alpha_terms', terms, expected_terms),
        ('stages', stages, expected_stages),
    ))
    errors = []
    def require(condition, message):
        if not condition and message not in errors:
            errors.append(message)
    for (t,p,o), row in solves.items():
        i=phase_index[p]
        solved = not (reference and p=='gas')
        require(row['solved']==int(solved) and row['mulesCalls']==(SUBCYCLES if solved else 0), 'Wrong alpha-call coverage')
        require(close(row['dt_s'],contract.dt), 'Unexpected alpha dt')
        require(abs(row['unassignedVolume_m3'])<contract.volume_identity_tol[i]
                and abs(row['unassignedFrozenRhoMass_kg'])<contract.frozen_mass_tol[i],
                'Alpha change ledger does not close')
        require(abs(row['divBoundaryDifference_m3'])<contract.volume_identity_tol[i], 'MULES boundary/divergence mismatch')
        stage=stages[(t,p,o,'alpha')]
        require(close(row['deltaVolume_m3'],stage['deltaVolume_m3'],contract.volume_identity_tol[i])
                and close(row['deltaFrozenRhoMass_kg'],stage['alphaPart_kg'],contract.frozen_mass_tol[i]),
                'Alpha-solve and solver-stage observers disagree')
    # Validate the observer identities independently from the reported remainders.
    for t,p,o in expected_solves:
        i=phase_index[p]
        get=lambda j: terms[(t,p,o,j)]['volume_m3']
        if reference and p=='gas':
            change=get('referenceSubcycle')+get('referenceFinal')
        else:
            change=sum(get(j) for j in ('mulesChange','implicitCorrection','referenceSubcycle','referenceFinal','clamp','rescale'))
            require(close(get('mulesChange'),get('oldTimeReset')+get('explicitSu')+get('implicitSp')-get('divergence')+get('equationResidual'),contract.volume_identity_tol[i]), 'MULES identity inconsistent')
            require(abs(get('equationResidualL1'))<contract.equation_l1_tol[i], 'MULES equation residual exceeds diagnostic tolerance')
            require(close(get('explicitSu')+get('implicitSp'),sum(get(j) for j in ('originalSu','rawDgdt','splitAtSourceAlpha','implicitSourceChange','sourceDecompositionRemainder')),contract.volume_identity_tol[i]), 'Source decomposition inconsistent')
        require(close(change,solves[(t,p,o)]['deltaVolume_m3'],contract.volume_identity_tol[i]), 'Reported alpha changes do not reproduce observed change')
    summary = {}
    for p in phases:
        i=phase_index[p]
        init = initial[(p,)]
        require(close(init['t_s'],0) and close(init['volume_m3'],contract.volumes[i],contract.volume_initial_tol[i]), 'Incorrect initial time or prescribed inventory')
        require(close(init['mass_kg'],contract.volumes[i]*(1 if p=='gas' else 1000),contract.mass_initial_tol[i]), 'Incorrect initial phase mass')
        phase_steps = [steps[(t,p)] for t in times]
        sums = {key: math.fsum(r[key] for r in phase_steps) for key in ('deltaMass_kg','boundaryOut_kg','optionSource_kg','transferSource_kg','physicalResidual_kg')}
        expected_source = contract.option_source_kg_s[i]*contract.dt*contract.steps if not closed else 0
        require(close(sums['optionSource_kg'],expected_source,1e-12), 'Actual assembled source differs from prescribed source')
        require(abs(sums['transferSource_kg'])<1e-13, 'Unexpected interphase mass transfer')
        require(close(phase_steps[-1]['mass_kg']-init['mass_kg'],sums['deltaMass_kg'],contract.mass_identity_tol[i]), 'Accepted-step deltas do not reproduce inventory change')
        for t,r in zip(times,phase_steps):
            require(close(r['dt_s'],contract.dt) and r['outers']==OUTERS, 'Unexpected timestep or outer-iteration count')
            require(r['alphaMin']>=-1e-10 and r['alphaMax']<=1+1e-10 and r['phaseSumMaxError']<=1e-10,
                    'Phase fractions violate bounds or sum constraint')
            require(close(r['physicalResidual_kg'],r['deltaMass_kg']+r['boundaryOut_kg']-r['optionSource_kg']-r['transferSource_kg'],contract.budget_arithmetic_tol[i]), 'Physical budget arithmetic inconsistent')
            sr=[stages[(t,p,o,s)] for o in range(1,OUTERS+1) for s in STAGES]+[stages[(t,p,OUTERS,'endStep')]]
            require(close(math.fsum(s['deltaMass_kg'] for s in sr),r['deltaMass_kg'],contract.mass_identity_tol[i]), 'Stage deltas do not reproduce accepted timestep')
            for s in sr:
                require(close(s['deltaMass_kg'],s['alphaPart_kg']+s['densityPart_kg'],contract.mass_identity_tol[i]) and abs(s['identityResidual_kg'])<contract.mass_identity_tol[i],
                        'Alpha/density stage identity inconsistent')
        # Numerical-screen tolerances, not experimental validation tolerances.
        tolerance = 1e-12+1e-9*init['mass_kg'] if closed else contract.budget_absolute[i]+contract.budget_relative[i]*max(abs(sums['optionSource_kg']),abs(sums['boundaryOut_kg']))
        residual_l1=math.fsum(abs(r['physicalResidual_kg']) for r in phase_steps)
        require(residual_l1 <= tolerance, p+': phase mass budget exceeds screening tolerance')
        if closed:
            require(abs(sums['boundaryOut_kg'])<1e-13, 'Unexpected external flux in closed periodic box')
            require(sum(solves[(t,p,o)]['changeL1Volume_m3'] for t in times for o in range(1,OUTERS+1))>1e-8, 'Closed transport test did not move the phase fractions')
        summary[p] = dict(initialMass_kg=init['mass_kg'],finalMass_kg=phase_steps[-1]['mass_kg'],
            **sums, residualL1_kg=residual_l1,screenTolerance_kg=tolerance,
            finalOuterRescaleVolume_m3=sum(terms.get((t,p,OUTERS,'rescale'),{}).get('volume_m3',0) for t in times),
            finalOuterSplitVolume_m3=sum(terms.get((t,p,OUTERS,'splitAtSourceAlpha'),{}).get('volume_m3',0) for t in times))
    adjustments = [abs(r['t_s']-time_key(r['t_s'],contract))
                   for kind in KINDS[1:] for r in rows[kind]]
    return dict(case=case,complete=True,screen_pass=not errors,errors=errors,steps=contract.steps,phases=summary,
                timeKeyNormalization=dict(tolerance_s=TIME_KEY_TOLERANCE,
                    adjustedRecords=sum(d>0 for d in adjustments),
                    maxAdjustment_s=max(adjustments,default=0)),
                meaning='Numerical screening only; not SC1 reactor, free-surface, oxygen or biology validation.')
