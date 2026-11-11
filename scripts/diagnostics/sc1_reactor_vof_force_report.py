"""Check passive two-step agreement and independently reconstruct peak-cell forces."""
import math
import re
import numpy as np
from sc1_reactor_watch import Monitor

END=1e-4
VECTORS=('U','H','gravity','capillary','pressure','net','freshPressure','geometryDelta')


def vector(row,name):
    return np.array([row[name+k] for k in ('x','y','z')]) if name in ('U','H') else np.array([row[name+k] for k in 'XYZ'])


def audit_prefix(text):
    rows=[]
    for line in text.splitlines():
        if not re.match(r'^SC1_VOF_(INITIAL|PRESSURE|STEP) ',line):continue
        pairs=re.findall(r'(\w+)=([^\s]+)',line);row={k:float(v) for k,v in pairs}
        if len(row)!=len(pairs) or not all(math.isfinite(v) for v in row.values()):raise ValueError('Malformed audit prefix')
        if row.get('t_s',0)<=END+1e-12:rows.append(dict(marker=line.split()[0],**row))
    return rows


def compare_prefix(text,control):
    left=audit_prefix(text);right=audit_prefix(control)
    counts={m:sum(r['marker']==m for r in left) for m in ('SC1_VOF_INITIAL','SC1_VOF_PRESSURE','SC1_VOF_STEP')}
    if counts!={'SC1_VOF_INITIAL':1,'SC1_VOF_PRESSURE':80,'SC1_VOF_STEP':2} or len(left)!=len(right):
        raise ValueError('Incomplete two-step passive-observer prefix')
    for a,b in zip(left,right):
        if a.keys()!=b.keys():raise ValueError('Audit prefix shape changed')
        for k,value in a.items():
            if k=='marker':
                if value!=b[k]:raise ValueError('Audit prefix order changed')
            elif not math.isclose(value,b[k],rel_tol=1e-9,abs_tol=1e-12):
                raise ValueError('Force observer changed control trajectory: '+k)
    return left


def rows(text,kind):
    result=[]
    for line in text.splitlines():
        if not line.startswith('SC1_VOF_FORCE_'+kind+' '):continue
        pairs=re.findall(r'(\w+)=([^\s]+)',line)
        row={k:v if k=='phase' else float(v) for k,v in pairs}
        if len(row)!=len(pairs) or not all(math.isfinite(v) for k,v in row.items() if k!='phase'):
            raise ValueError('Malformed force record')
        result.append(row)
    return result


def close(a,b,description):
    if not np.allclose(a,b,rtol=1e-9,atol=1e-10):raise ValueError(description)


def analyze(text,control):
    monitor=Monitor()
    if any(monitor.feed(l) for l in text.splitlines()):raise ValueError('Force probe contains a solver stop')
    times=re.findall(r'^Time = (\S+)$',text,re.M)
    if len(times)!=2 or any(abs(float(t)-v)>1e-12 for t,v in zip(times,(5e-5,1e-4))) or text.splitlines().count('End')!=1:
        raise ValueError('Force probe did not complete exactly two steps')
    if text.splitlines().count('SC1_VOF_FORCE_PROBE version=1 maxSteps=2')!=1:raise ValueError('Missing force observer marker')
    prefix=compare_prefix(text,control)
    checks=rows(text,'CHECK');cells=rows(text,'CELL');faces=rows(text,'FACE');coverage=rows(text,'FACES')
    if len(checks)!=8 or len(cells)!=16 or len(coverage)!=16:raise ValueError('Incomplete force observation coverage')
    expected={(i,p) for i in range(1,9) for p in ('gas','liquid')}
    if {(r['sequence'],r['phase']) for r in cells}!=expected or {(r['sequence'],r['phase']) for r in coverage}!=expected:
        raise ValueError('Duplicate/missing force observations')
    if not {(r['sequence'],r['phase']) for r in faces}.issubset(expected):raise ValueError('Unexpected force face group')
    for i,r in enumerate(checks,1):
        if r['sequence']!=i or abs(r['t_s']-((i-1)//4+1)*5e-5)>1e-12 or r['nonfinite']!=0:
            raise ValueError('Wrong force-check sequence/time')
        if not 0<=r['forceClosure_m_s']<=1e-9 or not 0<=r['updateClosure_m_s']<=1e-10:
            raise ValueError('Force decomposition or velocity update does not close')
    for cell in cells:
        seq=cell['sequence'];phase=cell['phase']
        if abs(cell['t_s']-((seq-1)//4+1)*5e-5)>1e-12 or cell['proc']!=int(cell['proc']) or not 0<=cell['proc']<8 or cell['localCell']<0 or cell['localCell']!=int(cell['localCell']) or cell['rAU']<=0:
            raise ValueError('Invalid force-cell identity')
        a=cell['alphaLiquid']
        if not -1e-10<=a<=1+1e-10 or (a if phase=='liquid' else 1-a)<=.01:raise ValueError('Wrong present-phase peak')
        close(cell['rho_kg_m3'],1+999*a,'Force-cell density disagrees with fraction')
        ff=[r for r in faces if r['sequence']==seq and r['phase']==phase]
        count=next(r['count'] for r in coverage if r['sequence']==seq and r['phase']==phase)
        if len(ff)!=count or count<4 or len({r['localFace'] for r in ff})!=count:raise ValueError('Missing/duplicate peak faces')
        sf=np.array([[r['Sf'+k] for k in 'XYZ'] for r in ff]);area=np.linalg.norm(sf,axis=1)
        if np.any(area<=0):raise ValueError('Invalid face area')
        nh=sf/area[:,None];matrix=np.einsum('fi,fj,f->ij',nh,nh,area)
        for f in ff:
            if f['proc']!=cell['proc'] or f['outwardSign'] not in (-1,1) or f['rAUf']<=0:raise ValueError('Wrong face identity')
            close(f['gravity_Pa_m'],-f['ghFace']*f['snGradRho'],'Gravity face term disagrees')
            close(f['geometryDelta_Pa_m'],(f['ghFace']-f['ghInterpolated'])*f['snGradRho'],'Geometry delta disagrees')
            parts=[f[k+'_Pa_m'] for k in ('gravity','capillary','pressure')]
            # Large pressure/gravity terms can cancel to a near-zero net force.
            roundoff=128*np.finfo(float).eps*sum(abs(x) for x in parts)
            if abs(f['net_Pa_m']-sum(parts))>roundoff+1e-10+1e-9*abs(f['net_Pa_m']):
                raise ValueError('Face force sum does not close')
        for term in VECTORS[2:]:
            try:reconstructed=cell['rAU']*np.linalg.solve(matrix,np.sum(nh*np.array([r[term+'_Pa_m'] for r in ff])[:,None]*area[:,None],axis=0))
            except np.linalg.LinAlgError as error:raise ValueError('Singular peak-cell reconstruction') from error
            close(reconstructed,vector(cell,term),'Independent face reconstruction disagrees: '+term)
        close(vector(cell,'net'),sum((vector(cell,k) for k in ('gravity','capillary','pressure'))),'Cell force sum does not close')
        close(vector(cell,'U'),vector(cell,'H')+vector(cell,'net'),'Cell velocity update does not close')
        if seq in (4,8):
            step=next(r for r in prefix if r['marker']=='SC1_VOF_STEP' and abs(r['t_s']-cell['t_s'])<1e-12)
            close(np.linalg.norm(vector(cell,'U')),step[phase+'Max_m_s'],'Peak speed disagrees with step audit')
        cell['magnitudes_m_s']={k:float(np.linalg.norm(vector(cell,k))) for k in VECTORS}
        cell['pressureFreshDifference_m_s']=float(np.linalg.norm(vector(cell,'pressure')-vector(cell,'freshPressure')))
    return dict(diagnostic_complete=True,accepted_startup=False,steps=2,prefixAgreement=dict(records=len(prefix),relativeTolerance=1e-9,absoluteTolerance=1e-12),
        checks=checks,cells=cells,faces=faces,faceCoverage=coverage,
        meaning='Two-step passive force decomposition, not a corrected simulation or startup acceptance.',
        units='U/H and reconstructed force contributions are m/s; each contribution includes the same cell rAU as the native velocity correction. Face force terms are Pa/m.',
        geometryDeltaMeaning='Frozen-pressure diagnostic difference if ghFace were replaced by interpolated cell gh. It is not applied, reprojected or established as a correction.',
        faceMeaning='Native local face orientation; outwardSign relates it to the selected cell. Processor copies are local evidence, not a global face integral. dC is owner-to-neighbour for internal faces and zero for boundary faces.')
