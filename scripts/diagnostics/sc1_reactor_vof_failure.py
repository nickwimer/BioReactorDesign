"""Summarize incomplete VOF audit logs without accepting a run or changing files."""
import math
import re

FIELDS = {
    'INITIAL': 'liquid_m3 gas_m3 version',
    'PRESSURE': 't_s min_Pa max_Pa nonfinite',
    'STEP': ('t_s dt_s liquid_m3 gas_m3 liquidOut_m3_s gasOut_m3_s '
             'postPressureOut_m3_s liquidResidual_kg gasResidual_kg liquidL1_kg '
             'gasL1_kg alphaMin alphaMax complement liquidMax_m_s gasMax_m_s '
             'alphaChangeL1_m3 Co alphaCalls'),
}


def scan(path):
    counts={key:0 for key in FIELDS}; rows={key:[] for key in FIELDS}; malformed=[]
    with path.open(errors='replace') as source:
        for line,text in enumerate(source,1):
            match=re.match(r'^SC1_VOF_(INITIAL|PRESSURE|STEP)\s+',text.strip())
            if not match:continue
            kind=match[1];counts[kind]+=1
            pairs=re.findall(r'(\w+)=([^\s]+)',text); raw=dict(pairs)
            try:
                values={key:float(raw[key]) for key in FIELDS[kind].split()}
                if len(pairs)!=len(raw) or not all(math.isfinite(x) for x in values.values()):
                    raise ValueError('Duplicate key or nonfinite value')
                if kind=='INITIAL' and (values['version']!=1 or min(values['liquid_m3'],values['gas_m3'])<0):
                    raise ValueError('Unsupported initial record')
                if kind=='PRESSURE' and (values['t_s']<=0 or values['min_Pa']>values['max_Pa']
                        or values['nonfinite']<0 or not values['nonfinite'].is_integer()):
                    raise ValueError('Invalid pressure record')
                if kind=='STEP' and (values['t_s']<=0 or values['dt_s']<=0
                        or values['alphaMin']>values['alphaMax']
                        or not values['alphaCalls'].is_integer()
                        or any(values[k]<0 for k in ('liquid_m3','gas_m3','liquidL1_kg','gasL1_kg',
                            'complement','liquidMax_m_s','gasMax_m_s','alphaChangeL1_m3','Co','alphaCalls'))):
                    raise ValueError('Invalid step record')
            except (KeyError,ValueError):
                malformed.append(dict(line=line,kind=kind,record=raw));continue
            rows[kind].append(dict(line=line,record=values))
    if not sum(counts.values()):return None
    initial=rows['INITIAL'][0]['record'] if counts['INITIAL']==1 and len(rows['INITIAL'])==1 else None
    steps=rows['STEP'];pressures=rows['PRESSURE'];seen=set();duplicates=0
    first={};events={}

    def event(kind,row,value,limit):
        events[kind]=events.get(kind,0)+1
        first.setdefault(kind,dict(line=row['line'],t_s=row['record']['t_s'],value=value,limit=limit))

    for row in pressures:
        r=row['record']
        if r['nonfinite']:event('nonfinite_pressure',row,r['nonfinite'],0)
        if r['min_Pa']<1000:event('pressure_below_1000_Pa',row,r['min_Pa'],1000)
    for row in steps:
        r=row['record'];key=round(r['t_s'],12)
        if key in seen:duplicates+=1
        seen.add(key)
        for key,limit in (('alphaMax',1+1e-10),('complement',1e-10),('Co',.2),
                          ('liquidMax_m_s',1e-3),('alphaChangeL1_m3',1e-5)):
            if r[key]>limit:event(key,row,r[key],limit)
        if r['alphaMin']< -1e-10:event('alphaMin',row,r['alphaMin'],-1e-10)
        if r['alphaCalls']!=2:event('alphaCalls',row,r['alphaCalls'],2)
        if initial:
            for phase,rho in (('liquid',1000),('gas',1)):
                limit=1e-10+1e-8*rho*initial[phase+'_m3'];key=phase+'L1_kg'
                if r[key]>limit:event(key,row,r[key],limit)
    return dict(recordCounts=counts,validStepRecordCount=len(steps),uniqueStepTimeCount=len(seen),
        duplicateStepRecords=duplicates,malformedRecordCount=len(malformed),
        firstMalformedRecord=malformed[0] if malformed else None,initial=initial,
        lastValidStep=steps[-1] if steps else None,
        minimumFinitePressure_Pa=min((p['record']['min_Pa'] for p in pressures),default=None),
        maximumFinitePressure_Pa=max((p['record']['max_Pa'] for p in pressures),default=None),
        maximumLoggedCo=max((s['record']['Co'] for s in steps),default=None),
        maximumLiquidSpeed_m_s=max((s['record']['liquidMax_m_s'] for s in steps),default=None),
        maximumGasSpeed_m_s=max((s['record']['gasMax_m_s'] for s in steps),default=None),
        largestFiniteAlphaMagnitude=max((max(abs(s['record']['alphaMin']),abs(s['record']['alphaMax']))
            for s in steps),default=None),firstLimitEvents=first,limitEventCounts=events,
        meaning='Logged diagnostics only. A step is printed before its guard; record counts are not accepted steps. '
                'No independent budget reconciliation, field validation or startup acceptance.',
        budgetLimitMeaning='Initial-inventory limits are evaluated only with one valid version-1 initial record.',
        legacyCounterMeaning='Top-level pairedStepRecordCount counts SC1_TWO_STEP only, not SC1_VOF_STEP.')
