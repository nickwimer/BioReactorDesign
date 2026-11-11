"""Read saved end fields to compare counterflow, bulk motion and peak locations."""
from pathlib import Path
import numpy as np
from sc1_reactor_startup_case import END, internal


def statistics(velocity, weights, mask, centres, alpha):
    count=int(mask.sum())
    if not count:return dict(cells=0)
    ids=np.flatnonzero(mask);u=velocity[mask];w=weights[mask]
    speed=np.linalg.norm(u,axis=1);cell=int(ids[np.argmax(speed)])
    return dict(cells=count,rms_m_s=float(np.sqrt(np.dot(w,speed*speed)/w.sum())),
        meanVertical_m_s=float(np.dot(w,u[:,1])/w.sum()),max_m_s=float(speed.max()),
        peak=dict(cell=cell,centre_m=centres[cell].tolist(),alphaLiquid=float(alpha[cell]),
            velocity_m_s=velocity[cell].tolist()))


def flow(case, cells):
    c=internal(case/'0/C',cells,True);v=internal(case/'0/V',cells)
    end=case/f'{END:.12g}'
    a=internal(end/'alpha.liquid',cells);b=internal(end/'alpha.gas',cells)
    if np.max(np.abs(a+b-1))>1e-10 or np.min(a)<-1e-10 or np.max(a)>1+1e-10:
        raise ValueError('Invalid fractions in saved flow comparison')
    ug=internal(end/'U.gas',cells,True);ul=internal(end/'U.liquid',cells,True)
    mixed=(a>.01)&(b>.01)
    phases={}
    for phase,alpha,u in (('gas',b,ug),('liquid',a,ul)):
        phases[phase]={name:statistics(u,v*alpha,mask,c,a) for name,mask in
            (('present',alpha>.01),('bulk',alpha>=.99),('interface',mixed))}
    rho=b+1000*a
    volume_velocity=b[:,None]*ug+a[:,None]*ul
    mass_velocity=(b[:,None]*ug+1000*a[:,None]*ul)/rho[:,None]
    return dict(phases=phases,
        interfaceSlip=statistics(ug-ul,v,mixed,c,a),
        volumeVelocityRms_m_s=float(np.sqrt(np.dot(v,np.sum(volume_velocity**2,axis=1))/v.sum())),
        massVelocityRms_m_s=float(np.sqrt(np.dot(v*rho,np.sum(mass_velocity**2,axis=1))/np.dot(v,rho))))


def compare_flow(root,manifest,candidate):
    import sc1_reactor_startup as startup
    baseline_root=Path(manifest['baseline']['root'])
    _,baseline=startup.saved_startup(baseline_root)
    if startup.baseline_identity(baseline_root,baseline)!=manifest['baseline']:
        raise ValueError('Control changed during comparison')
    measures={}
    for label,folder,report in (('control',baseline_root,baseline),('drag_enabled',root,candidate)):
        measures[label]=dict(log_sha256=report['log_sha256'],mass_screen_pass=report['screen_pass'],
            startup_pass=report['startup_pass'],maxCourant=report['maxCourant'],
            finalLiquidVolume_L=report['snapshots'][-1]['liquidVolume_L'],
            finalAlphaChangeL1_m3=report['snapshots'][-1]['alphaChangeL1_m3'],
            finalFlow=flow(folder/'quiescent',startup.CELLS))
    return dict(suite='sc1_startup_drag_comparison',changed_input='constant/phaseProperties: blending/drag only',
        time_s=END,results=measures,
        meaning='Same thresholds and initial fields. A listed drag model with canceled blending is not active coupling; '
                'restoring its weight does not itself establish quiescent or reactor validity.')
