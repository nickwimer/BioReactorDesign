"""Locate early quiescent velocities in reconstructed native fields."""
import numpy as np
from sc1_reactor_startup_case import internal, G
from sc1_reactor_startup_flow import statistics
from sc1_wall_mesh import Reader, read_file

FIELDS=('alpha.gas','alpha.liquid','U.gas','U.liquid','p_rgh')


def patch_cells(case,cells):
    mesh=case/'constant/polyMesh'
    owner_reader=Reader(read_file(mesh/'owner')[0]);owner=owner_reader.label_list();owner_reader.finish()
    boundary_reader=Reader(read_file(mesh/'boundary')[0]);patches=boundary_reader.named_list();boundary_reader.finish()
    if np.any(owner<0) or np.any(owner>=cells):raise ValueError('Invalid boundary owner labels')
    masks={}
    for name,p in patches.items():
        start=int(p['startFace'][0]);count=int(p['nFaces'][0])
        if start<0 or count<0 or start+count>len(owner):raise ValueError('Boundary faces outside owner array')
        mask=np.zeros(cells,dtype=bool);mask[owner[start:start+count]]=True;masks[name]=mask
    for name in ('walls','propeller','outlet'):
        if name not in masks:raise ValueError('Missing required reactor patch '+name)
    return masks


def measure(folder,centres,volumes,a0,p0,patches,midheight):
    n=len(volumes);a=internal(folder/'alpha.liquid',n);b=internal(folder/'alpha.gas',n)
    if min(a.min(),b.min())<-1e-10 or max(a.max(),b.max())>1+1e-10 or np.max(abs(a+b-1))>1e-10:
        raise ValueError('Invalid reconstructed fractions')
    # Validate then remove only tolerated floating-point negatives from weights.
    aw=np.maximum(a,0);bw=np.maximum(b,0)
    ug=internal(folder/'U.gas',n,True);ul=internal(folder/'U.liquid',n,True)
    prgh=internal(folder/'p_rgh',n)
    pressure=prgh-(b+1000*a)*G*centres[:,1]
    if pressure.min()<1000:raise ValueError('Invalid reconstructed static pressure')
    mixed=(a>.01)&(b>.01);rows={}
    for name,alpha,weight,u in (('gas',b,bw,ug),('liquid',a,aw,ul)):
        present=alpha>.01
        masks=dict(present=present,bulk=alpha>=.99,interface=mixed,
            wallAdjacent=present&patches['walls'],propellerAdjacent=present&patches['propeller'],
            outletAdjacent=present&patches['outlet'])
        rows[name]={}
        for region,mask in masks.items():
            result=statistics(u,volumes*weight,mask,centres,a)
            if result['cells']:
                peak=result['peak'];i=peak['cell']
                peak.update(radius_m=float(np.hypot(centres[i,0],centres[i,2])),
                    heightFromInitialMidplane_m=float(centres[i,1]-midheight),
                    adjacentPatches=[k for k,v in patches.items() if v[i]],
                    alphaGas=float(b[i]),p_Pa=float(pressure[i]),p_rgh_Pa=float(prgh[i]),
                    pressureChangeFromInitial_Pa=float(pressure[i]-p0[i]),
                    gasVelocity_m_s=ug[i].tolist(),liquidVelocity_m_s=ul[i].tolist(),
                    slipVelocity_m_s=(ug[i]-ul[i]).tolist(),
                    volumeMixtureVelocity_m_s=(b[i]*ug[i]+a[i]*ul[i]).tolist())
            rows[name][region]=result
    mixture=b[:,None]*ug+a[:,None]*ul
    return dict(time_s=float(folder.name),liquidVolume_L=float(1000*np.dot(a,volumes)),
        alphaChangeL1_m3=float(np.dot(abs(a-a0),volumes)),
        pressureMinimum_Pa=float(pressure.min()),pressureMaximum_Pa=float(pressure.max()),
        volumeMixtureRms_m_s=float(np.sqrt(np.dot(volumes,np.sum(mixture*mixture,axis=1))/volumes.sum())),
        phases=rows,interfaceSlip=statistics(ug-ul,volumes,mixed,centres,a))


def snapshots(case,times,cells,midheight):
    c=internal(case/'0/C',cells,True);v=internal(case/'0/V',cells)
    if np.any(v<=0):raise ValueError('Invalid native cell volumes')
    a0=internal(case/'0/alpha.liquid',cells);p0=internal(case/'0/p',cells)
    masks=patch_cells(case,cells)
    return [measure(case/time,c,v,a0,p0,masks,midheight) for time in times]
