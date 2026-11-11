"""Independently reconcile the observed VOF phase fluxes and native snapshots."""
import math
import re
import numpy as np
from sc1_reactor_startup_case import internal,DT,END,STEPS,G
from sc1_reactor_startup_flow import statistics
from sc1_reactor_velocity_fields import patch_cells


def rows(text,marker):
    return [{k:float(v) for k,v in re.findall(r'(\w+)=([^\s]+)',line)}
        for line in text.splitlines() if line.startswith(marker+' ')]


def budget(text,initial):
    starts=rows(text,'SC1_VOF_INITIAL');steps=rows(text,'SC1_VOF_STEP');ps=rows(text,'SC1_VOF_PRESSURE')
    if len(starts)!=1 or len(steps)!=STEPS or len(ps)!=4*STEPS or text.splitlines().count('End')!=1:
        raise ValueError('Incomplete VOF reference record coverage')
    if not all(math.isfinite(v) for r in starts+steps+ps for v in r.values()):raise ValueError('Nonfinite VOF records')
    from sc1_reactor_watch import Monitor
    monitor=Monitor()
    if any(monitor.feed(line) for line in text.splitlines()):raise ValueError('VOF log contains a numerical stop')
    for i,r in enumerate(ps):
        if abs(r['t_s']-(i//4+1)*DT)>1e-12 or r['nonfinite']!=0 or r['min_Pa']<1000 or r['max_Pa']<r['min_Pa']:
            raise ValueError('Invalid VOF pressure coverage')
    phases={};errors=[]
    for name,rho,key in (('liquid',1000,'liquid_volume_m3'),('gas',1,'gas_volume_m3')):
        start=starts[0][name+'_m3'];previous=start;out=0.;residuals=[];logged=[]
        if abs(start-initial[key])>1e-12:raise ValueError('Initial VOF inventory differs from preserved fill')
        for i,r in enumerate(steps,1):
            if abs(r['t_s']-i*DT)>1e-12 or abs(r['dt_s']-DT)>1e-12 or r['alphaCalls']!=2:
                raise ValueError('Unexpected VOF time or alpha-call coverage')
            v=r[name+'_m3'];flux=r[name+'Out_m3_s'];residual=rho*(v-previous+DT*flux)
            if v<0 or abs(residual-r[name+'Residual_kg'])>2e-11*rho/1000:
                raise ValueError('VOF printed budget does not reconcile')
            residuals.append(residual);logged.append(r[name+'Residual_kg']);out+=rho*DT*flux;previous=v
            if abs(sum(abs(x) for x in logged)-r[name+'L1_kg'])>1e-12*max(rho/1000,1e-3):
                raise ValueError('VOF cumulative residual disagrees')
        tolerance=1e-10+1e-8*rho*start;l1=sum(abs(x) for x in residuals)
        if l1>tolerance:errors.append(name+' phase budget exceeds unchanged screen')
        phases[name]=dict(initialMass_kg=rho*start,finalMass_kg=rho*previous,boundaryOut_kg=out,
            physicalResidual_kg=rho*(previous-start)+out,residualL1_kg=l1,screenTolerance_kg=tolerance)
    for r in steps:
        if r['alphaMin']< -1e-10 or r['alphaMax']>1+1e-10 or r['complement']>1e-10:
            errors.append('VOF fractions failed unchanged bounds')
        if r['Co']<0 or r['Co']>.2:errors.append('VOF Courant number exceeded unchanged limit')
        if r['liquidMax_m_s']<0 or r['liquidMax_m_s']>1e-3:errors.append('VOF liquid speed exceeded 1 mm/s')
        if r['alphaChangeL1_m3']<0 or r['alphaChangeL1_m3']>1e-5:errors.append('VOF redistribution exceeded 0.01 L')
    return dict(complete=True,errors=sorted(set(errors)),phases=phases,steps=STEPS,
        maxCourant=max(r['Co'] for r in steps),maxLiquid_m_s=max(r['liquidMax_m_s'] for r in steps),
        maxGas_m_s=max(r['gasMax_m_s'] for r in steps),minimumPressure_Pa=min(r['min_Pa'] for r in ps)),steps


def analyze(case,initial,cells):
    result,steps=budget((case/'log.solver').read_text(),initial)
    c=internal(case/'0/C',cells,True);v=internal(case/'0/V',cells);a0=internal(case/'0/alpha.liquid',cells)
    if np.any(v<=0):raise ValueError('Invalid VOF native volumes')
    patches=patch_cells(case,cells);snapshots=[]
    for index in range(5,STEPS+1,5):
        t=index*DT;folder=case/f'{t:.12g}';a=internal(folder/'alpha.liquid',cells);u=internal(folder/'U',cells,True)
        pressure=internal(folder/'p_rgh',cells)-(1+999*a)*G*c[:,1]
        if a.min()< -1e-10 or a.max()>1+1e-10 or pressure.min()<1000:raise ValueError('Invalid VOF native state')
        inventory=float(np.dot(a,v));speed=float(np.linalg.norm(u[a>.01],axis=1).max())
        row=steps[index-1]
        if abs(inventory-row['liquid_m3'])>1e-11 or abs(speed-row['liquidMax_m_s'])>1e-10:
            raise ValueError('VOF native inventory/speed disagrees with logged state')
        phase_rows={}
        for phase,alpha in (('liquid',a),('gas',1-a)):
            phase_rows[phase]={name:statistics(u,v*np.maximum(alpha,0),mask,c,a) for name,mask in
                (('present',alpha>.01),('bulk',alpha>=.99),('interface',(a>.01)&(a<.99)),
                 ('wallAdjacent',(alpha>.01)&patches['walls']),('propellerAdjacent',(alpha>.01)&patches['propeller']))}
        snapshots.append(dict(time_s=t,liquidVolume_L=1000*inventory,
            alphaChangeL1_m3=float(np.dot(abs(a-a0),v)),pressureMinimum_Pa=float(pressure.min()),
            pressureMaximum_Pa=float(pressure.max()),phases=phase_rows))
    result.update(snapshots=snapshots,reference_pass=not result['errors'],
        meaning='10-ms unaerated single-velocity VOF reference only. No BiRD, dispersed-bubble, oxygen-transfer, biology or operating-checkpoint acceptance.',
        fluxMeaning='Last outer Euler alpha transport flux; gas flux is its complement using the pre-pressure mixture flux employed by that alpha solve. Post-pressure mixture boundary flux is logged separately.')
    return result
