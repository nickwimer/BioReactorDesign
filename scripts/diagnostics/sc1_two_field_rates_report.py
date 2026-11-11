"""Same physical budgets and tolerances, plus mandatory rate-split evidence."""
import math
from sc1_two_field_cases import CASES as BASELINE_CASES, OUTERS
from sc1_two_field_bounds_cases import CASES as BOUND_CASES, contract
from sc1_two_field_report import BASELINE_CONTRACT, parse, unique, require_coverage
from sc1_two_field_report import summarize as baseline_report
from sc1_two_field_bounds_report import summarize as bound_report
from sc1_two_field_rates_instrument import TAG


def summarize(text, case):
    if text.splitlines().count('SC1_RATE_LIBRARY version=1')!=1:
        raise ValueError('Missing unique rate-split runtime marker')
    if case in BASELINE_CASES:
        report=baseline_report(text,case); spec=BASELINE_CONTRACT
    elif case in BOUND_CASES:
        report=bound_report(text,case,require_split_effect=False); spec=contract(case)
    else: raise ValueError('Unknown rate-split case')
    return validate_rates(text,report,spec,source_expected=not case.startswith('closed_'))


def validate_rates(text, report, spec, *, source_expected):
    """Apply the same runtime evidence checks to an explicit source contract."""
    if text.splitlines().count('SC1_RATE_LIBRARY version=1')!=1:
        raise ValueError('Missing unique rate-split runtime marker')
    lines=[line for line in text.splitlines() if line.startswith('SC1_RATE_STATE')]
    rows=parse('\n'.join(line.replace('SC1_RATE_STATE','SC1_PHASE_SOLVE',1) for line in lines))['SC1_PHASE_SOLVE']
    keys={'t_s','call','valid','applied','qGasMin','qGasMax','qLiquidMin','qLiquidMax',
          'alphaSumMaxError','rawIdentityScaled'}
    if any(set(r)!=keys for r in rows):raise ValueError('Unexpected rate-state record fields')
    indexed=unique(rows,('t_s','call'),spec)
    expected={(round((s+1)*spec.dt,12),o) for s in range(spec.steps) for o in range(1,OUTERS+1)}
    require_coverage((('rate_states',indexed,expected),))
    errors=report['errors']
    def require(condition,message):
        if not condition and message not in errors:errors.append(message)
    for r in rows:
        require(r['valid'] in (0,1) and r['applied']==2*r['valid'],'Incomplete rate-split application')
        require(0<=r['alphaSumMaxError']<=1e-12 and 0<=r['rawIdentityScaled']<=1e-12,
                'Rate-split complement/raw-source identity failed')
        require(r['qGasMin']==r['qGasMax']==0,'Unexpected gas pressure-source rate')
        require(0<=r['qLiquidMin']<=r['qLiquidMax'],'Negative/reversed liquid pressure-source rates')
        require(r['valid']==1 or r['qLiquidMax']==0,'Unavailable rate was used')
    if not source_expected:
        require(all(r['qLiquidMax']==0 for r in rows),'Unexpected closed-case pressure source')
    else:
        require(any(r['qLiquidMax']>0 for r in rows),'Liquid pressure-source rate was never exercised')
    terms=parse(text)['SC1_PHASE_TERM']
    for phase in ('gas','liquid'):
        split=math.fsum(abs(r['volume_m3']) for r in terms
            if r['phase']==phase and r['term']=='splitAtSourceAlpha')
        report['phases'][phase]['allOuterSplitL1Volume_m3']=split
        require(split<=1e-18,phase+': rate split does not preserve the assembled coupled source')
    errors.sort()
    mass_errors={p+': phase mass budget exceeds screening tolerance' for p in ('gas','liquid')}
    report.update(suite=TAG,sourceFloor=None,screen_pass=not errors,
        diagnostic_pass=not any(e not in mass_errors for e in errors),
        rateStates=len(rows),positiveRateCalls=sum(r['qLiquidMax']>0 for r in rows),
        meaning='Private two-field nonnegative-rate experiment. Original mass thresholds retained; '
                'no reactor, signed-source or interphase-transfer validation.')
    return report
