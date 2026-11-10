"""Validate complete observer streams; never mistake process success for conservation."""
from collections import defaultdict
import math
from pathlib import Path
from sc1_two_field_cases import CASES, DT, END, OUTERS, SOURCE_KG_S, STEPS, SUBCYCLES

TERMS = ('mulesChange oldTimeReset explicitSu implicitSp divergence equationResidual '
         'equationResidualL1 originalSu rawDgdt pressureDilatation couplingDifference '
         'splitAtSourceAlpha implicitSourceChange sourceDecompositionRemainder '
         'implicitCorrection referenceSubcycle clamp referenceFinal implicitDivergence '
         'implicitEquationResidual implicitEquationResidualL1 finalTransportDivergence '
         'fluxAveragingDefect rescale').split()
STAGES = ('alpha', 'modelCorrect', 'transportPressure', 'endOuter')
KINDS = ('SC1_TWO_INITIAL', 'SC1_PHASE_SOLVE', 'SC1_PHASE_TERM', 'SC1_TWO_STAGE', 'SC1_TWO_STEP')
TEXT_KEYS = {'phase', 'term', 'stage'}


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


def unique(rows, fields):
    result = {}
    for row in rows:
        key = tuple(row[k] for k in fields)
        if key in result:
            raise ValueError('Duplicate record: '+repr(key))
        result[key] = row
    return result


def close(a, b, tol=1e-12):
    return abs(a-b) <= tol


def summarize(text, case):
    if case not in CASES:
        raise ValueError('Unknown case')
    if text.count('SC1_TWO_LIBRARY version=1') != 1 or '\nEnd\n' not in text:
        raise ValueError('Missing unique library marker or successful solver End')
    if any(s in text for s in ('FOAM FATAL', 'Floating point exception', 'Segmentation fault')):
        raise ValueError('Solver failure in log')
    rows = parse(text)
    initial = unique(rows['SC1_TWO_INITIAL'], ('phase',))
    steps = unique(rows['SC1_TWO_STEP'], ('t_s', 'phase'))
    solves = unique(rows['SC1_PHASE_SOLVE'], ('t_s', 'phase', 'call'))
    terms = unique(rows['SC1_PHASE_TERM'], ('t_s', 'phase', 'call', 'term'))
    stages = unique(rows['SC1_TWO_STAGE'], ('t_s', 'phase', 'outer', 'stage'))
    times = [round((i+1)*DT, 12) for i in range(STEPS)]
    phases = ('gas', 'liquid')
    expected_steps = {(t,p) for t in times for p in phases}
    expected_solves = {(t,p,o) for t in times for p in phases for o in range(1,OUTERS+1)}
    if set(initial) != {(p,) for p in phases} or set(steps) != expected_steps or set(solves) != expected_solves:
        raise ValueError('Incomplete or unexpected initial/step/alpha-call records')
    reference = case.endswith('gas_reference')
    expected_terms = {(t,p,o,j) for t,p,o in expected_solves
                      for j in (('referenceSubcycle','referenceFinal') if reference and p=='gas' else TERMS)}
    expected_stages = {(t,p,o,s) for t,p,o in expected_solves for s in STAGES}
    expected_stages |= {(t,p,OUTERS,'endStep') for t in times for p in phases}
    if set(terms) != expected_terms or set(stages) != expected_stages:
        raise ValueError('Incomplete or unexpected alpha-term/stage records')
    errors = []
    def require(condition, message):
        if not condition and message not in errors:
            errors.append(message)
    for (t,p,o), row in solves.items():
        solved = not (reference and p=='gas')
        require(row['solved']==int(solved) and row['mulesCalls']==(SUBCYCLES if solved else 0), 'Wrong alpha-call coverage')
        require(close(row['dt_s'],DT), 'Unexpected alpha dt')
        require(abs(row['unassignedVolume_m3'])<1e-13 and abs(row['unassignedFrozenRhoMass_kg'])<1e-10,
                'Alpha change ledger does not close')
        require(abs(row['divBoundaryDifference_m3'])<1e-13, 'MULES boundary/divergence mismatch')
        stage=stages[(t,p,o,'alpha')]
        require(close(row['deltaVolume_m3'],stage['deltaVolume_m3'],1e-13)
                and close(row['deltaFrozenRhoMass_kg'],stage['alphaPart_kg'],1e-10),
                'Alpha-solve and solver-stage observers disagree')
    # Validate the observer identities independently from the reported remainders.
    for t,p,o in expected_solves:
        get=lambda j: terms[(t,p,o,j)]['volume_m3']
        if reference and p=='gas':
            change=get('referenceSubcycle')+get('referenceFinal')
        else:
            change=sum(get(j) for j in ('mulesChange','implicitCorrection','referenceSubcycle','referenceFinal','clamp','rescale'))
            require(close(get('mulesChange'),get('oldTimeReset')+get('explicitSu')+get('implicitSp')-get('divergence')+get('equationResidual'),1e-13), 'MULES identity inconsistent')
            require(abs(get('equationResidualL1'))<1e-12, 'MULES equation residual exceeds diagnostic tolerance')
            require(close(get('explicitSu')+get('implicitSp'),sum(get(j) for j in ('originalSu','rawDgdt','splitAtSourceAlpha','implicitSourceChange','sourceDecompositionRemainder')),1e-13), 'Source decomposition inconsistent')
        require(close(change,solves[(t,p,o)]['deltaVolume_m3'],1e-13), 'Reported alpha changes do not reproduce observed change')
    summary = {}
    for p in phases:
        init = initial[(p,)]
        require(close(init['t_s'],0) and close(init['volume_m3'],.0005,1e-11), 'Incorrect initial time or prescribed inventory')
        require(close(init['mass_kg'], .0005 if p=='gas' else .5,1e-8), 'Incorrect initial phase mass')
        phase_steps = [steps[(t,p)] for t in times]
        sums = {key: math.fsum(r[key] for r in phase_steps) for key in ('deltaMass_kg','boundaryOut_kg','optionSource_kg','transferSource_kg','physicalResidual_kg')}
        expected_source = SOURCE_KG_S*END if case.startswith('source_') and p=='liquid' else 0
        require(close(sums['optionSource_kg'],expected_source,1e-12), 'Actual assembled source differs from prescribed source')
        require(abs(sums['transferSource_kg'])<1e-13, 'Unexpected interphase mass transfer')
        require(close(phase_steps[-1]['mass_kg']-init['mass_kg'],sums['deltaMass_kg'],1e-11), 'Accepted-step deltas do not reproduce inventory change')
        for t,r in zip(times,phase_steps):
            require(close(r['dt_s'],DT) and r['outers']==OUTERS, 'Unexpected timestep or outer-iteration count')
            require(r['alphaMin']>=-1e-10 and r['alphaMax']<=1+1e-10 and r['phaseSumMaxError']<=1e-10,
                    'Phase fractions violate bounds or sum constraint')
            require(close(r['physicalResidual_kg'],r['deltaMass_kg']+r['boundaryOut_kg']-r['optionSource_kg']-r['transferSource_kg'],1e-12), 'Physical budget arithmetic inconsistent')
            sr=[stages[(t,p,o,s)] for o in range(1,OUTERS+1) for s in STAGES]+[stages[(t,p,OUTERS,'endStep')]]
            require(close(math.fsum(s['deltaMass_kg'] for s in sr),r['deltaMass_kg'],1e-11), 'Stage deltas do not reproduce accepted timestep')
            for s in sr:
                require(close(s['deltaMass_kg'],s['alphaPart_kg']+s['densityPart_kg'],1e-11) and abs(s['identityResidual_kg'])<1e-11,
                        'Alpha/density stage identity inconsistent')
        closed=case.startswith('closed_')
        # Numerical-screen tolerances, not experimental validation tolerances.
        tolerance = 1e-12+1e-9*init['mass_kg'] if closed else 1e-10+1e-3*max(abs(sums['optionSource_kg']),abs(sums['boundaryOut_kg']))
        residual_l1=math.fsum(abs(r['physicalResidual_kg']) for r in phase_steps)
        require(residual_l1 <= tolerance, p+': phase mass budget exceeds screening tolerance')
        if closed:
            require(abs(sums['boundaryOut_kg'])<1e-13, 'Unexpected external flux in closed periodic box')
            require(sum(solves[(t,p,o)]['changeL1Volume_m3'] for t in times for o in range(1,OUTERS+1))>1e-8, 'Closed transport test did not move the phase fractions')
        summary[p] = dict(initialMass_kg=init['mass_kg'],finalMass_kg=phase_steps[-1]['mass_kg'],
            **sums, residualL1_kg=residual_l1,screenTolerance_kg=tolerance,
            finalOuterRescaleVolume_m3=sum(terms.get((t,p,OUTERS,'rescale'),{}).get('volume_m3',0) for t in times),
            finalOuterSplitVolume_m3=sum(terms.get((t,p,OUTERS,'splitAtSourceAlpha'),{}).get('volume_m3',0) for t in times))
    return dict(case=case,complete=True,screen_pass=not errors,errors=errors,steps=STEPS,phases=summary,
                meaning='Numerical screening only; not SC1 reactor, free-surface, oxygen or biology validation.')
