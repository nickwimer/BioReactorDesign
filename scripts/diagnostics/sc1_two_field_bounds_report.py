"""Near-bound phase budgets and evidence of floor activity; no tolerance fitting."""
import math
from sc1_two_field_cases import OUTERS, SOURCE_KG_S
from sc1_two_field_report import parse, summarize_contract
from sc1_two_field_bounds_cases import ALPHA, FLOOR, VERSION, contract


def summarize(text, case):
    spec=contract(case)
    report=summarize_contract(text,case,spec,closed=False,
                              reference=case.endswith('gas_reference'))
    rows=parse(text)
    errors=report['errors']
    def require(condition,message):
        if not condition and message not in errors: errors.append(message)
    for i,p in enumerate(('gas','liquid')):
        steps=[r for r in rows['SC1_TWO_STEP'] if r['phase']==p]
        terms=[r for r in rows['SC1_PHASE_TERM'] if r['phase']==p and r['call']==OUTERS]
        selected=lambda term:[r for r in terms if r['term']==term]
        summary=report['phases'][p]
        minimum=min(r['alphaMin'] for r in steps)
        maximum=max(r['alphaMax'] for r in steps)
        exposed=sum((r['alphaMax']<FLOOR if p=='gas' else r['alphaMin']>1-FLOOR)
                    for r in steps)
        require(exposed==spec.steps,p+': accepted steps did not all stay in the target near-bound regime')
        if p=='gas':
            require(minimum>1e-8,'Gas reached the residual-alpha floor; source-floor isolation failed')
        for r in steps:
            require(r['alphaMin']<=r['alphaMax'],p+': reversed alpha extrema')
            require(r['phaseSumMaxError']>=0,p+': negative phase-sum error')
            require(abs(r['mass_kg']-r['volume_m3']*(1 if p=='gas' else 1000))
                    <=spec.mass_identity_tol[i],p+': constant-density mass/volume mismatch')
            supply=SOURCE_KG_S*spec.dt if p=='liquid' else 0
            require(abs(r['optionSource_kg']-supply)<=spec.budget_arithmetic_tol[i],
                    p+': unexpected per-step option source')
            require(abs(r['transferSource_kg'])<=spec.budget_arithmetic_tol[i],
                    p+': unexpected per-step transfer source')
        for r in terms:
            require(abs(r['frozenRhoMass_kg']-r['volume_m3']*(1 if p=='gas' else 1000))
                    <=spec.frozen_mass_tol[i],p+': constant-density term mass/volume mismatch')
        throughput=max(abs(summary['optionSource_kg']),abs(summary['boundaryOut_kg']))
        summary.update(
            alphaMinimum=minimum,alphaMaximum=maximum,nearBoundSteps=exposed,
            maxAbsStepResidual_kg=max(abs(r['physicalResidual_kg']) for r in steps),
            residualL1OverInitialMass=summary['residualL1_kg']/summary['initialMass_kg'],
            residualL1OverThroughput=summary['residualL1_kg']/throughput if throughput else None,
            finalOuterSplitL1Volume_m3=math.fsum(abs(r['volume_m3']) for r in selected('splitAtSourceAlpha')),
            finalOuterSplitFrozenMass_kg=math.fsum(r['frozenRhoMass_kg'] for r in selected('splitAtSourceAlpha')),
            finalOuterImplicitSourceFrozenMass_kg=math.fsum(r['frozenRhoMass_kg'] for r in selected('implicitSourceChange')),
            finalOuterRescaleFrozenMass_kg=math.fsum(r['frozenRhoMass_kg'] for r in selected('rescale')))
    floor_activity=report['phases']['liquid']['finalOuterSplitL1Volume_m3']>1e-18
    require(floor_activity,'Liquid source-split effect was not measurably exercised')
    errors.sort()
    mass_failures={p+': phase mass budget exceeds screening tolerance' for p in ('gas','liquid')}
    report.update(suite='sc1_two_field_bounds',contractVersion=VERSION,dt_s=spec.dt,
                  endTime_s=spec.dt*spec.steps,
                  initialFractions={p:float(a) for p,a in ALPHA.items()},
                  sourceFloor=FLOOR,sourceSplitObserved=floor_activity,
                  diagnostic_pass=not any(e not in mass_failures for e in errors),
                  screen_pass=not errors,
                  meaning='Near-bound numerical diagnostic; a completed failed screen remains evidence. '
                          'No reactor, interface or interphase-transfer validation.')
    return report
