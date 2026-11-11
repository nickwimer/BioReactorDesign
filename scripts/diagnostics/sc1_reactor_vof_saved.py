#!/usr/bin/env python3
"""Inspect existing VOF processor fields in memory; no OpenFOAM command or CFD."""
import argparse
import json
import math
from pathlib import Path
import re
import sys
import tarfile
import numpy as np
import sc1_two_field as base
from sc1_reactor_failure import regular, sha
from sc1_reactor_vof_failure import scan
from sc1_reactor_startup_case import internal, G
from sc1_reactor_startup_flow import statistics
from sc1_reactor_velocity_fields import patch_cells
from sc1_wall_mesh import Reader, read_file

LOG_SHA='2811c733dc1fe00499022c7d1b447d26ff6ec45facf145261f5f5dfae9190a9b'
MANIFEST_SHA='12418a7b61431d0fbdacd6e4178e7bacc6e6a070912ce0a86f58ec01bacb2fc9'
CELLS=281664
RANKS=8
TIMES=('0.00025','0.001','0.004','0.00825')
FIELDS=('alpha.liquid','U','p_rgh')


def budget(text,initial):
    """Reconcile the observed prefix, including the final guard-rejected record."""
    def rows(marker):
        result=[]
        for line in text.splitlines():
            if not line.startswith(marker+' '):continue
            pairs=re.findall(r'(\w+)=([^\s]+)',line)
            if len(dict(pairs))!=len(pairs):raise ValueError('Duplicate audit key')
            result.append({k:float(v) for k,v in pairs})
        return result
    starts=rows('SC1_VOF_INITIAL');steps=rows('SC1_VOF_STEP');ps=rows('SC1_VOF_PRESSURE')
    if len(starts)!=1 or len(steps)!=168 or len(ps)!=4*len(steps):
        raise ValueError('Unexpected stopped VOF record coverage')
    if not all(math.isfinite(v) for r in starts+steps+ps for v in r.values()):
        raise ValueError('Nonfinite VOF record')
    if starts[0]['version']!=1:raise ValueError('Unknown VOF audit version')
    for i,r in enumerate(steps,1):
        if abs(r['t_s']-i*5e-5)>1e-12 or abs(r['dt_s']-5e-5)>1e-12 or r['alphaCalls']!=2:
            raise ValueError('Unexpected VOF time/alpha coverage')
    for i,r in enumerate(ps):
        if abs(r['t_s']-(i//4+1)*5e-5)>1e-12 or r['nonfinite']!=0 or r['min_Pa']<1000 or r['max_Pa']<r['min_Pa']:
            raise ValueError('Unexpected VOF pressure coverage')
    phases={}
    for phase,rho in (('liquid',1000),('gas',1)):
        start=starts[0][phase+'_m3'];previous=start;out=0.;residuals=[];logged_l1=0.;max_error=0.
        if abs(start-initial[phase+'_volume_m3'])>1e-12:raise ValueError('VOF initial inventory mismatch')
        for r in steps:
            value=r[phase+'_m3'];flux=r[phase+'Out_m3_s'];res=rho*(value-previous+r['dt_s']*flux)
            difference=abs(res-r[phase+'Residual_kg']);max_error=max(max_error,difference)
            if value<0 or difference>2e-11*rho/1000:raise ValueError('VOF phase budget does not reconcile')
            residuals.append(res);logged_l1+=abs(r[phase+'Residual_kg']);out+=rho*r['dt_s']*flux;previous=value
            if abs(logged_l1-r[phase+'L1_kg'])>1e-12*max(rho/1000,1e-3):
                raise ValueError('VOF cumulative residual mismatch')
        tolerance=1e-10+1e-8*rho*start;l1=sum(abs(r) for r in residuals)
        phases[phase]=dict(initialMass_kg=rho*start,finalMass_kg=rho*previous,boundaryOut_kg=out,
            physicalResidual_kg=rho*(previous-start)+out,residualL1_kg=l1,
            screenTolerance_kg=tolerance,withinBudgetTolerance=l1<=tolerance,
            maxPrintedResidualDifference_kg=max_error)
    return dict(observedSteps=len(steps),pressureRecords=len(ps),phases=phases,
        meaning='Observed prefix only, including its final rejected state; no completed startup acceptance.'),steps


def labels(path,obj):
    reader=Reader(read_file(path)[0])
    if reader.header.get('object')!=[obj.encode()] or reader.header.get('class')!=[b'labelList']:
        raise ValueError('Wrong addressing header: '+str(path))
    values=reader.label_list();reader.finish();return values


def addressing(case,ranks,cells,track):
    maps=[labels(track(f'processor{rank}/constant/polyMesh/cellProcAddressing'),'cellProcAddressing')
          for rank in range(ranks)]
    all_ids=np.concatenate(maps)
    if not np.array_equal(np.sort(all_ids),np.arange(cells)):
        raise ValueError('Processor addressing is not a permutation of global cells')
    return maps


def field(path,count,name):
    vector=name in ('U','C');header=Reader(path.read_bytes()).header
    if header.get('object')!=[name.encode()] or header.get('class')!=[b'volVectorField' if vector else b'volScalarField']:
        raise ValueError('Wrong saved field header: '+str(path))
    return internal(path,count,vector)


def assemble(time,name,maps,cells,track):
    result=np.empty((cells,3) if name in ('U','C') else cells)
    for rank,ids in enumerate(maps):result[ids]=field(track(f'processor{rank}/{time}/{name}'),len(ids),name)
    return result


def measure(time,a,u,prgh,c,v,a0,p0,patches,midheight,owner,neighbour):
    if a.min()< -1e-10 or a.max()>1+1e-10:raise ValueError('Saved fractions outside tolerance')
    pressure=prgh-(1+999*a)*G*c[:,1];speed=np.linalg.norm(u,axis=1)
    if pressure.min()<1000:raise ValueError('Saved pressure below floor')
    phases={};peaks=set()
    for phase,alpha in (('liquid',a),('gas',1-a)):
        present=alpha>.01
        regions=dict(present=present,bulk=alpha>=.99,interface=(a>.01)&(a<.99),
            wallAdjacent=present&patches['walls'],propellerAdjacent=present&patches['propeller'],
            outletAdjacent=present&patches['outlet'])
        phases[phase]={}
        for name,mask in regions.items():
            stat=statistics(u,v*np.maximum(alpha,0),mask,c,a)
            if stat['cells']:
                peak=stat['peak'];i=peak['cell'];peaks.add(i)
                peak.update(radius_m=float(np.hypot(c[i,0],c[i,2])),
                    heightFromInitialMidplane_m=float(c[i,1]-midheight),alphaInitial=float(a0[i]),
                    adjacentPatches=[k for k,m in patches.items() if m[i]],
                    p_Pa=float(pressure[i]),p_rgh_Pa=float(prgh[i]),
                    pressureChangeFromInitial_Pa=float(pressure[i]-p0[i]),
                    horizontalSpeed_m_s=float(np.hypot(u[i,0],u[i,2])))
            phases[phase][name]=stat
    # Exact native values around each peak, not an approximation to OF force operators.
    stencils=[];internal_owner=owner[:len(neighbour)]
    for i in sorted(peaks):
        adjacent=np.unique(np.r_[neighbour[internal_owner==i],internal_owner[neighbour==i]])
        ids=np.r_[i,adjacent]
        stencils.append(dict(peakCell=i,cells=[dict(cell=int(j),centre_m=c[j].tolist(),
            volume_m3=float(v[j]),alphaLiquid=float(a[j]),alphaInitial=float(a0[j]),
            p_rgh_Pa=float(prgh[j]),p_Pa=float(pressure[j]),initialPressure_Pa=float(p0[j]),
            velocity_m_s=u[j].tolist()) for j in ids]))
    return dict(time_s=float(time),liquidVolume_L=float(1000*np.dot(a,v)),
        alphaChangeL1_m3=float(np.dot(abs(a-a0),v)),pressureMinimum_Pa=float(pressure.min()),
        pressureMaximum_Pa=float(pressure.max()),volumeVelocityRms_m_s=float(np.sqrt(np.dot(v,speed**2)/v.sum())),
        phases=phases,peakNeighbourStencils=stencils)


def inspect(root,out):
    root=root.resolve(strict=True);mp=regular(root,'startup_manifest.json');m=json.loads(mp.read_text())
    if sha(mp)!=MANIFEST_SHA:raise ValueError('Requires the received VOF manifest')
    if m['root']!=str(root) or m['variant']!='sc1_meml_vof_reference_v1' or m['ranks']!=RANKS or m['initial']['cells']!=CELLS:
        raise ValueError('Wrong VOF source case')
    log=regular(root,'quiescent/log.solver')
    if sha(log)!=LOG_SHA:raise ValueError('Requires the received stopped VOF log')
    attempt=regular(root,'quiescent/attempt.json')
    if json.loads(attempt.read_text())!={k:m[k] for k in ('build','assets','ranks')}:
        raise ValueError('Attempt identity mismatch')
    protected=[root,base.REPO]+[Path(m[k]) for k in ('build','archive','regression_root')]
    protected += [Path(m[k]['root']) for k in ('baseline','control')]
    out=base.fresh(out,protected);case=root/'quiescent';hashes={}
    for rel,digest in m['inputs'].items():
        p=regular(case,rel)
        if sha(p)!=digest:raise ValueError('Initial input changed: '+rel)
        hashes[rel]=digest
    def track(rel):
        p=regular(case,rel);digest=sha(p)
        if rel in hashes and hashes[rel]!=digest:raise ValueError('Input changed during inspection: '+rel)
        hashes[rel]=digest;return p
    source_hashes={name:sha(regular(root,name)) for name in ('startup_manifest.json','quiescent/attempt.json','quiescent/log.solver')}
    budget_result,steps=budget(log.read_text(),m['initial']);diagnostic=scan(log)
    cells=m['initial']['cells'];maps=addressing(case,m['ranks'],cells,track)
    c=field(track('0/C'),cells,'C');v=field(track('0/V'),cells,'V')
    if np.any(v<=0):raise ValueError('Invalid cell volume')
    a0=field(track('0/alpha.liquid'),cells,'alpha.liquid');p0=field(track('0/p'),cells,'p')
    for name in ('C','alpha.liquid','U','p_rgh'):
        original=field(track('0/'+name),cells,name);mapped=assemble('0',name,maps,cells,track)
        if not np.allclose(original,mapped,rtol=0,atol=1e-12):
            raise ValueError('Processor mapping does not reproduce initial '+name)
    patches=patch_cells(case,cells)
    owner=labels(track('constant/polyMesh/owner'),'owner');neighbour=labels(track('constant/polyMesh/neighbour'),'neighbour')
    if len(neighbour)>len(owner) or np.any(owner>=cells) or np.any(neighbour>=cells):
        raise ValueError('Invalid neighbour addressing')
    snapshots=[]
    for time in TIMES:
        print('Reading saved processor fields at '+time+' s',flush=True)
        a,u,p=[assemble(time,name,maps,cells,track) for name in FIELDS]
        snapshot=measure(time,a,u,p,c,v,a0,p0,patches,m['initial']['interface_midheight_m'],owner,neighbour)
        row=next(r for r in steps if abs(r['t_s']-float(time))<1e-12)
        for phase in ('liquid','gas'):
            if abs(snapshot['phases'][phase]['present']['max_m_s']-row[phase+'Max_m_s'])>1e-10:
                raise ValueError('Saved velocity disagrees with audit: '+phase)
        if abs(snapshot['liquidVolume_L']/1000-row['liquid_m3'])>1e-11 or abs(snapshot['alphaChangeL1_m3']-row['alphaChangeL1_m3'])>1e-11:
            raise ValueError('Saved inventory/redistribution disagrees with audit')
        snapshots.append(snapshot)
    if any(sha(regular(case,rel))!=digest for rel,digest in hashes.items()) or any(
            sha(regular(root,rel))!=digest for rel,digest in source_hashes.items()):
        raise ValueError('Source changed during inspection')
    report=dict(suite='sc1_vof_saved_spatial_v1',sourceRoot=str(root),sourceHashes=source_hashes,
        sourceCommit=m['source_commit'],loggedPrefix=budget_result,loggedDiagnostics=diagnostic,
        snapshots=snapshots,fieldHashes=hashes,initial=m['initial'],
        provenance='Native processor internal fields mapped by cellProcAddressing; full global-cell permutation and mapped initial fields checked. No external build/control-chain revalidation.',
        meaning='Saved-output diagnostic only. Source startup failed. No CFD, reconstruction, restart, field edits or acceptance.',
        patchMeaning='Adjacent means owner of a boundary face, not distance to a wall.',
        stencilMeaning='Raw cell/neighbour values only; no discrete force balance or cause inferred.')
    out.mkdir(parents=True,exist_ok=False);base.write_json(out/'report.json',report)
    bundle=out/'vof_spatial_evidence.tar.gz'
    with tarfile.open(bundle,'w:gz') as archive:
        archive.add(out/'report.json',arcname='report.json')
        for rel in source_hashes:archive.add(regular(root,rel),arcname='source/'+rel)
    if any(sha(regular(root,rel))!=digest for rel,digest in source_hashes.items()):
        raise ValueError('Source changed while bundling; preserve output')
    print(json.dumps(dict(loggedPrefix=budget_result,snapshots=[dict(time_s=s['time_s'],
        liquidVolume_L=s['liquidVolume_L'],volumeVelocityRms_m_s=s['volumeVelocityRms_m_s'],
        liquidPresent=s['phases']['liquid']['present'],gasPresent=s['phases']['gas']['present'])
        for s in snapshots]),indent=2))
    print('REPORT='+str(out/'report.json'));print('BUNDLE='+str(bundle));print('BUNDLE_SHA256='+sha(bundle))
    print('Saved VOF inspection complete; no CFD or case changes.')
    return report


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args(argv);inspect(a.root,a.out);return 0


if __name__=='__main__':
    try:raise SystemExit(main())
    except (OSError,ValueError,KeyError,TypeError,StopIteration,tarfile.TarError) as error:
        print('Saved VOF inspection stopped: '+str(error)+'; preserve the source and any output',file=sys.stderr)
        raise SystemExit(1)
