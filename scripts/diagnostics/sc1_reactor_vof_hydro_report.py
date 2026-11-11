"""Original budgets and all 40 snapshots, with physical pressure reconstruction."""
import numpy as np
import sc1_reactor_vof_nonorth_report as previous
from sc1_reactor_startup_case import internal, DT, STEPS, G
from sc1_reactor_startup_flow import statistics
from sc1_reactor_velocity_fields import patch_cells


def budget(text,initial):
    # Translate only the pressure unknown's name for the existing pass parser.
    # The stored solver log remains unchanged; all audits use physical pressure.
    if 'Solving for p_rgh,' in text:raise ValueError('Wrong pressure variable')
    records=previous.reference.rows(text,'SC1_HYDROSTATIC_REFERENCE')
    if len(records)!=1 or set(records[0])!={'version','frozen','rhoError_kg_m3','pressureError_Pa','initialSpeed_m_s'}:
        raise ValueError('Missing/malformed hydrostatic reference check')
    r=records[0]
    if not all(np.isfinite(v) for v in r.values()) or r['version']!=1 or r['frozen']!=1 or not 0<=r['rhoError_kg_m3']<=1e-9 or not 0<=r['pressureError_Pa']<=1e-7 or r['initialSpeed_m_s']!=0:
        raise ValueError('Invalid hydrostatic reference')
    report,steps=previous.budget(text.replace('Solving for pPert,','Solving for p_rgh,'),initial)
    report['hydrostaticReference']=r
    return report,steps


def analyze(case,initial,cells):
    result,steps=budget((case/'log.solver').read_text(),initial)
    c=internal(case/'0/C',cells,True);v=internal(case/'0/V',cells);a0=internal(case/'0/alpha.liquid',cells)
    if np.any(v<=0):raise ValueError('Invalid native volumes')
    ph=internal(case/'0/p_rgh',cells)-(1+999*a0)*G*c[:,1]
    patches=patch_cells(case,cells);snapshots=[]
    for index in range(5,STEPS+1,5):
        t=index*DT;folder=case/f'{t:.12g}';a=internal(folder/'alpha.liquid',cells);u=internal(folder/'U',cells,True)
        pressure=internal(folder/'p_rgh',cells)-(1+999*a)*G*c[:,1]
        if a.min()< -1e-10 or a.max()>1+1e-10 or pressure.min()<1000:raise ValueError('Invalid native state')
        inventory=float(np.dot(a,v));row=steps[index-1]
        if abs(inventory-row['liquid_m3'])>1e-11:raise ValueError('Native inventory disagrees with audit')
        perturbation=internal(folder/'pPert',cells)
        absolute=internal(folder/'p',cells)
        if not np.allclose(pressure,ph+perturbation,rtol=0,atol=1e-7) or not np.allclose(absolute,pressure,rtol=0,atol=1e-7):
            raise ValueError('Written physical/perturbation pressure identity failed')
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
    if result['maxGas_m_s']>1e-3:result['errors'].append('Gas speed exceeded 1 mm/s')
    if not result['pressureIteration']['pressure_iteration_pass']:
        result['errors'].append('Nonorthogonal pressure corrections did not converge')
    result.update(snapshots=snapshots,hydrostatic_pass=not result['errors'],
        meaning='10-ms fixed-reference hydrostatic startup only. No stirred/aerated flow, operating checkpoint or biology acceptance.')
    return result
