"""Fresh 1000-L, unaerated, stationary-vessel startup on the archived MEML mesh.

This is a quiescent water/gas diagnostic, not a copy of the reacting case.
The verified rate-split binary's pure-phase density guards remain unchanged.
"""
import re
import numpy as np
from sc1_two_field_cases import case_files as fixture, foam
from sc1_two_field_rates_instrument import APP

DT = 0.00005
STEPS = 200
WRITE = 40
FACE_WRITE = 5
END = DT*STEPS
WIDTH = 0.0125
TARGET_M3 = 1.0
G = 9.81
PATCHES = {'outlet':'patch', 'propeller':'wall', 'walls':'wall',
    'inside_to_hub':'empty', 'inside_to_hub_copy':'empty',
    'hub_to_rotor':'empty', 'hub_to_rotor_copy':'empty'}


def field(obj, dimensions, internal, vector=False):
    boundaries=[]
    for patch,kind in PATCHES.items():
        if kind=='empty': bc='type empty;'
        elif obj.startswith('U.'):
            bc=('type pressureInletOutletVelocity; value uniform (0 0 0);'
                if patch=='outlet' else 'type noSlip;')
        elif obj.startswith('alpha.'):
            incoming=1 if obj=='alpha.gas' else 0
            bc=(f'type inletOutlet; inletValue uniform {incoming}; value uniform {incoming};'
                if patch=='outlet' else 'type zeroGradient;')
        elif obj=='p_rgh':
            bc=('type fixedValue; value uniform 100000;' if patch=='outlet'
                else 'type fixedFluxPressure; value uniform 100000;')
        elif obj=='p': bc='type calculated; value uniform 100000;'
        else: bc='type zeroGradient;'
        boundaries.append(f'{patch} {{ {bc} }}')
    body=f'dimensions {dimensions};\ninternalField {internal};\nboundaryField\n{{\n'+ '\n'.join(boundaries)+'\n}'
    return foam(obj,body,'volVectorField' if vector else 'volScalarField')


def case_files(ranks, *, segregated_drag=False, face_momentum=False):
    if not isinstance(ranks,int) or not 1<=ranks<=52:
        raise ValueError('Startup requires 1 to 52 ranks')
    if face_momentum and not segregated_drag:raise ValueError('Face comparison requires restored segregated drag')
    old=fixture('source_gas_reference')
    files={k:v for k,v in old.items() if k.startswith('constant/') or k in ('system/fvSchemes','system/fvSolution','system/fvConstraints')}
    files['constant/g']=foam('g','dimensions [0 1 -2 0 0 0 0]; value (0 -9.81 0);','uniformDimensionedVectorField')
    files['constant/fvModels']=foam('fvModels','')
    files['constant/phaseProperties']=foam('phaseProperties','''type basicMultiphaseSystem;
phases (gas liquid); referencePhase gas;
gas { type pureIsothermalPhaseModel; diameterModel constant;
      constantCoeffs { d 0.01; } residualAlpha 1e-6; }
liquid { type pureIsothermalPhaseModel; diameterModel constant;
         constantCoeffs { d 0.001; } residualAlpha 1e-6; }
blending { default { type linear;
    minFullyContinuousAlpha.gas 0.7; minPartlyContinuousAlpha.gas 0.3;
    minFullyContinuousAlpha.liquid 0.7; minPartlyContinuousAlpha.liquid 0.3; } }
surfaceTension ((gas and liquid) { type constant; sigma 0.07; });
interfaceCompression ((gas and liquid) 1);
aspectRatio ();
drag ((gas and liquid) { type segregated; m 0.5; n 8; });
virtualMass (); heatTransfer (); phaseTransfer (); lift (); wallLubrication (); turbulentDispersion ();
''')
    if segregated_drag:
        # The unordered segregated model is multiplied by (1-f1-f2).
        # The legacy symmetric 0.3/0.7 linear blend makes that nearly zero
        # for complementary fractions. OF9's none/none blend gives weight 1.
        anchor='minFullyContinuousAlpha.liquid 0.7; minPartlyContinuousAlpha.liquid 0.3; } }'
        assert files['constant/phaseProperties'].count(anchor)==1
        files['constant/phaseProperties']=files['constant/phaseProperties'].replace(
            anchor,anchor[:-2]+'\n    drag { type none; continuousPhase none; }\n}')
    # MEML globalVars specifies this water-viscosity law at T=303.15 K.
    mu=2.414e-5*10**(247.8/(303.15-140))
    files['constant/thermophysicalProperties.liquid']=files['constant/thermophysicalProperties.liquid'].replace('mu 0.001;',f'mu {mu:.17g};')
    # No MRF dictionary: both walls and propeller are stationary in this test.
    files['system/controlDict']=foam('controlDict',f'''application {APP};
startFrom startTime; startTime 0; stopAt endTime; endTime {END};
deltaT {DT}; adjustTimeStep false; maxCo 0.2; maxDeltaT {DT};
writeControl timeStep; writeInterval {WRITE}; purgeWrite 0;
writeFormat ascii; writePrecision 17; writeCompression off;
timeFormat general; timePrecision 12; runTimeModifiable false;''')
    if face_momentum:
        files['system/fvSolution']=files['system/fvSolution'].replace('PIMPLE {','PIMPLE { faceMomentum true;')
        files['system/controlDict']=files['system/controlDict'].replace(f'writeInterval {WRITE};',f'writeInterval {FACE_WRITE};')
    files['system/decomposeParDict']=foam('decomposeParDict',f'numberOfSubdomains {ranks}; method scotch;')
    for phase in ('gas','liquid'):
        files['0/alpha.'+phase]=field('alpha.'+phase,'[0 0 0 0 0 0 0]','uniform 0.5')
        files['0/U.'+phase]=field('U.'+phase,'[0 1 -1 0 0 0 0]','uniform (0 0 0)',True)
        files['0/T.'+phase]=field('T.'+phase,'[0 0 0 1 0 0 0]','uniform 303.15')
    for name in ('p','p_rgh'):
        files['0/'+name]=field(name,'[1 -1 -2 0 0 0 0]','uniform 100000')
    return files


def internal(path, count, vector=False):
    """Strict reader for our ASCII uniform/nonuniform scalar/vector fields."""
    text=path.read_text()
    text=re.sub(r'/\*.*?\*/|//[^\n]*',' ',text,flags=re.S)
    if not re.search(r'\bformat\s+ascii\s*;',text): raise ValueError('Expected ASCII field: '+str(path))
    width=3 if vector else 1
    m=re.search(r'\binternalField\s+uniform\s+([^;]+);',text)
    if m:
        values=np.fromstring(m[1].replace('(',' ').replace(')',' '),sep=' ')
        if len(values)!=width:raise ValueError('Invalid uniform field')
        data=np.tile(values,(count,1))
    else:
        kind='vector' if vector else 'scalar'
        m=re.search(r'\binternalField\s+nonuniform\s+List<'+kind+r'>\s+(\d+)\s*\((.*?)\)\s*;',text,re.S)
        if not m or int(m[1])!=count:raise ValueError('Missing/mismatched internal field: '+str(path))
        payload=m[2].replace('(',' ').replace(')',' ')
        tokens=payload.split()
        if len(tokens)!=count*width:raise ValueError('Field length mismatch')
        data=np.array([float(s) for s in tokens]).reshape(count,width)
    if not np.isfinite(data).all():raise ValueError('Nonfinite field')
    return data if vector else data[:,0]


def initialized_fields(centres, volumes, target=TARGET_M3):
    """Set exactly the prescribed inventory with a planar transition layer.

    Integrate the same continuous density profile analytically for p and p_rgh.
    Discrete force balance remains an outcome measured by the startup run.
    """
    y=np.asarray(centres)[:,1]; v=np.asarray(volumes)
    if len(v)!=len(y) or not np.isfinite(y).all() or not np.isfinite(v).all() or np.any(v<=0):
        raise ValueError('Invalid native mesh geometry')
    if not 0<target<float(v.sum()):raise ValueError('Target fill does not fit mesh')
    lo,hi=float(y.min()-WIDTH),float(y.max()+WIDTH)
    for _ in range(80):
        h=(lo+hi)/2
        if np.dot(np.clip((h+WIDTH/2-y)/WIDTH,0,1),v)<target:lo=h
        else:hi=h
    h=(lo+hi)/2
    a=np.clip((h+WIDTH/2-y)/WIDTH,0,1); x=y-h
    integral=np.where(x<-WIDTH/2,x+WIDTH/8,np.where(x>WIDTH/2,WIDTH/8,x/2-x*x/(2*WIDTH)))
    p=100000-G*y+999*G*(WIDTH/8-integral)
    prgh=p+(1+999*a)*G*y
    if abs(float(np.dot(a,v))-target)>1e-12:raise ValueError('Fill inversion failed')
    values={'alpha.liquid':a,'alpha.gas':1-a,'p':p,'p_rgh':prgh}
    files={}
    for name,data in values.items():
        payload='nonuniform List<scalar>\n'+str(len(data))+'\n(\n'+'\n'.join(f'{z:.17g}' for z in data)+'\n)'
        files['0/'+name]=field(name,'[0 0 0 0 0 0 0]' if name.startswith('alpha.') else '[1 -1 -2 0 0 0 0]',payload)
    return files,dict(cells=len(v),volume_m3=float(v.sum()),liquid_volume_m3=float(np.dot(a,v)),
        gas_volume_m3=float(np.dot(1-a,v)),interface_midheight_m=h,transition_width_m=WIDTH,
        pressure='Hydrostatic continuum profile matching alpha; vent p_rgh=100000 Pa',
        density_kg_m3={'gas':1,'liquid':1000})
