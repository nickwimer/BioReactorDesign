"""Validate every nonorthogonal pressure pass before reusing the VOF mass audit."""
import math
import re
import numpy as np
import sc1_reactor_vof_report as reference
from sc1_reactor_startup_case import internal, DT, STEPS, G
from sc1_reactor_startup_flow import statistics
from sc1_reactor_velocity_fields import patch_cells

PASSES=10
TOLERANCE=1e-12
SOLVE=re.compile(r'^\S+:\s+Solving for p_rgh, Initial residual = (\S+), Final residual = (\S+), No Iterations (\d+)\s*$')


def pressure_coverage(text,require_complete=True):
    """Keep final pressure observations for the old four-per-step budget parser.

    All 40 observations and solve records per step are validated first. This is
    an in-memory reporting view; the solver log and CFD are never modified.
    """
    view=[];groups=[];pending=None;time=None;step=0;count=0;total=0;minimum=math.inf
    for line in text.splitlines():
        if line.startswith('Time = '):
            if time is not None or pending is not None:raise ValueError('Incomplete pressure step')
            time=float(line.split('=',1)[1]);count=0
            if not math.isfinite(time) or abs(time-(step+1)*DT)>1e-12:raise ValueError('Unexpected pressure time')
        if 'Solving for p_rgh,' in line:
            match=SOLVE.fullmatch(line)
            if not match or pending is not None or time is None:raise ValueError('Malformed/duplicate pressure solve')
            pending=dict(initial=float(match[1]),final=float(match[2]),iterations=int(match[3]))
            if not all(math.isfinite(v) and v>=0 for v in pending.values()) or pending['final']>TOLERANCE:
                raise ValueError('Invalid/unconverged linear pressure solve')
        if line.startswith('SC1_VOF_PRESSURE '):
            pairs=re.findall(r'(\w+)=([^\s]+)',line);p={k:float(v) for k,v in pairs}
            if pending is None or time is None or len(pairs)!=len(p) or set(p)!={'t_s','min_Pa','max_Pa','nonfinite'}:
                raise ValueError('Unpaired/malformed pressure audit')
            if not all(math.isfinite(v) for v in p.values()) or abs(p['t_s']-time)>1e-12 or p['nonfinite']!=0 or p['min_Pa']<1000 or p['max_Pa']<p['min_Pa']:
                raise ValueError('Invalid pressure audit')
            minimum=min(minimum,p['min_Pa']);count+=1;total+=1
            if count>4*PASSES:raise ValueError('Extra pressure passes')
            if count%PASSES==0:
                groups.append(dict(t_s=time,correction=count//PASSES,**pending))
                view.append(line)
            pending=None
            continue
        if line.startswith('SC1_VOF_STEP '):
            if time is None or count!=4*PASSES or pending is not None:raise ValueError('Incomplete pressure passes')
            step+=1;time=None
        view.append(line)
    if step>STEPS:raise ValueError('Extra pressure steps')
    complete=(time is None and pending is None and step==STEPS and total==STEPS*4*PASSES
        and text.splitlines().count('End')==1)
    if require_complete and not complete:
        raise ValueError('Incomplete nonorthogonal pressure coverage')
    failures=[g for g in groups if g['initial']>TOLERANCE or g['iterations']!=0]
    return '\n'.join(view)+'\n',dict(complete=complete,observedStepRecords=step,
        unfinishedTime_s=time,pressureRecords=total,passesPerCorrection=PASSES,
        corrections=len(groups),convergedCorrections=len(groups)-len(failures),
        maxFinalPassInitialResidual=max((g['initial'] for g in groups),default=None),
        finalPassTolerance=TOLERANCE,pressure_iteration_pass=complete and not failures,
        firstUnconvergedCorrection=failures[0] if failures else None,
        minimumPressureAllPasses_Pa=minimum if total else None,
        convergenceMeaning='Final pass in every fixed ten-pass group must enter at <=1e-12 and perform zero linear iterations. This tests explicit-correction convergence, not facewise force balance.')


def budget(text,initial):
    view,coverage=pressure_coverage(text)
    report,steps=reference.budget(view,initial)
    report.update(pressureIteration=coverage)
    report['minimumPressure_Pa']=coverage['minimumPressureAllPasses_Pa']
    return report,steps


def analyze(case,initial,cells):
    result,steps=budget((case/'log.solver').read_text(),initial)
    c=internal(case/'0/C',cells,True);v=internal(case/'0/V',cells);a0=internal(case/'0/alpha.liquid',cells)
    if np.any(v<=0):raise ValueError('Invalid native volumes')
    patches=patch_cells(case,cells);snapshots=[]
    for index in range(5,STEPS+1,5):
        t=index*DT;folder=case/f'{t:.12g}';a=internal(folder/'alpha.liquid',cells);u=internal(folder/'U',cells,True)
        pressure=internal(folder/'p_rgh',cells)-(1+999*a)*G*c[:,1]
        if a.min()< -1e-10 or a.max()>1+1e-10 or pressure.min()<1000:raise ValueError('Invalid native state')
        inventory=float(np.dot(a,v));row=steps[index-1]
        if abs(inventory-row['liquid_m3'])>1e-11:raise ValueError('Native inventory disagrees with audit')
        phase_rows={}
        for phase,alpha in (('liquid',a),('gas',1-a)):
            speed=float(np.linalg.norm(u[alpha>.01],axis=1).max())
            if abs(speed-row[phase+'Max_m_s'])>1e-10:raise ValueError('Native phase speed disagrees with audit')
            phase_rows[phase]={name:statistics(u,v*np.maximum(alpha,0),mask,c,a) for name,mask in
                (('present',alpha>.01),('bulk',alpha>=.99),('interface',(a>.01)&(a<.99)),
                 ('wallAdjacent',(alpha>.01)&patches['walls']),('propellerAdjacent',(alpha>.01)&patches['propeller']))}
        snapshots.append(dict(time_s=t,liquidVolume_L=1000*inventory,
            alphaChangeL1_m3=float(np.dot(abs(a-a0),v)),pressureMinimum_Pa=float(pressure.min()),
            pressureMaximum_Pa=float(pressure.max()),phases=phase_rows))
    result.update(snapshots=snapshots,reference_pass=not result['errors'])
    return result
