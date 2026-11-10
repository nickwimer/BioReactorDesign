"""Exact small-case fixtures; these are numerical tests, not SC1 reactor settings."""
from pathlib import Path

CASES = ('closed_no_reference', 'closed_gas_reference',
         'source_no_reference', 'source_gas_reference')
END = 0.1
DT = 0.001
STEPS = 100
OUTERS = 2
SUBCYCLES = 6
SOURCE_KG_S = 0.001


def foam(obj, body, cls='dictionary'):
    return f'FoamFile\n{{ format ascii; class {cls}; object {obj}; }}\n'+body+'\n'


def case_files(name):
    if name not in CASES:
        raise ValueError('Unknown case: '+name)
    closed = name.startswith('closed_')
    reference = name.endswith('gas_reference')
    patches = ('xMin', 'xMax', 'yMin', 'yMax', 'zMin', 'zMax')
    faces = ('(0 4 7 3)', '(1 2 6 5)', '(0 1 5 4)', '(3 7 6 2)', '(0 3 2 1)', '(4 5 6 7)')
    boundary = ''
    for i, p in enumerate(patches):
        typ = f'type cyclic; neighbourPatch {patches[i ^ 1]};' if closed else ('type patch;' if p == 'xMax' else 'type wall;')
        boundary += f'{p} {{ {typ} faces ({faces[i]}); }}\n'
    files = {}
    def put(path, body, cls='dictionary'):
        files[path] = foam(Path(path).name, body, cls)
    put('system/blockMeshDict', '''convertToMeters 1;
vertices ((0 0 0) (0.1 0 0) (0.1 0.1 0) (0 0.1 0)
          (0 0 0.1) (0.1 0 0.1) (0.1 0.1 0.1) (0 0.1 0.1));
blocks (hex (0 1 2 3 4 5 6 7) (16 4 4) simpleGrading (1 1 1));
edges ();
boundary (\n'''+boundary+');\nmergePatchPairs ();')
    put('system/controlDict', f'''application birdTwoFieldConservation;
startFrom startTime; startTime 0; stopAt endTime; endTime {END};
deltaT {DT}; adjustTimeStep false; maxCo 0.2; maxDeltaT {DT};
writeControl timeStep; writeInterval {STEPS}; purgeWrite 0;
writeFormat ascii; writePrecision 15; writeCompression off;
timeFormat general; timePrecision 12; runTimeModifiable false;
''')
    put('system/fvSchemes', r'''ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; limited cellLimited Gauss linear 1; }
divSchemes
{
    default Gauss linear;
    "div\(phi,alpha.*\)" Gauss upwind;
    "div\(phir,alpha.*,alpha.*\)" Gauss upwind;
    "div\(alphaRhoPhi.*,U.*\)" Gauss upwind;
    "div\(phi.*,U.*\)" Gauss upwind;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; bounded corrected; }
wallDist { method meshWave; }
''')
    put('system/fvSolution', f'''solvers
{{
    alpha {{ nAlphaCorr 1; nAlphaSubCycles {SUBCYCLES}; sc1PhaseAudit true;
             meanFluxReference false; splitPhaseFlux false; }}
    p_rgh {{ solver PCG; preconditioner DIC; tolerance 1e-12; relTol 0; }}
    p_rghFinal {{ $p_rgh; }}
    "U.*" {{ solver smoothSolver; smoother symGaussSeidel; tolerance 1e-12; relTol 0; }}
}}
PIMPLE {{ nOuterCorrectors {OUTERS}; nCorrectors 2; nNonOrthogonalCorrectors 0;
          pRefCell 0; pRefValue 100000; }}
relaxationFactors {{ equations {{ ".*" 1; }} }}
''')
    put('constant/g', 'dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 0);', 'uniformDimensionedVectorField')
    put('constant/dynamicMeshDict', 'dynamicFvMesh staticFvMesh;')
    put('constant/phaseProperties', '''type basicMultiphaseSystem;
phases (gas liquid);
'''+('referencePhase gas;\n' if reference else '')+'''
gas { type pureIsothermalPhaseModel; diameterModel constant;
      constantCoeffs { d 0.003; } residualAlpha 1e-8; }
liquid { type pureIsothermalPhaseModel; diameterModel constant;
         constantCoeffs { d 0.003; } residualAlpha 1e-8; }
blending { default { type linear;
    minFullyContinuousAlpha.gas 0.7; minPartlyContinuousAlpha.gas 0.3;
    minFullyContinuousAlpha.liquid 0.7; minPartlyContinuousAlpha.liquid 0.3; } }
surfaceTension ((gas and liquid) { type constant; sigma 0; });
interfaceCompression ();
aspectRatio ();
drag ((gas in liquid) { type SchillerNaumann; residualRe 1e-3; swarmCorrection { type none; } }
      (liquid in gas) { type SchillerNaumann; residualRe 1e-3; swarmCorrection { type none; } });
virtualMass (); heatTransfer (); phaseTransfer (); lift (); wallLubrication (); turbulentDispersion ();
''')
    put('constant/fvModels', '' if closed else f'''liquidFeed
{{ type massSource; selectionMode all; massFlowRate {SOURCE_KG_S};
   phase liquid; rho thermo:rho.liquid;
   fieldValues {{ U.liquid (0 0 0); T.liquid 300; }} }}''')
    put('system/fvConstraints', '''limitp
{
    type limitPressure;
    min 1000;
}''')
    for phase, density, cp, mu, mw in [('gas', 1, 1000, 1.8e-5, 28.97), ('liquid', 1000, 4180, .001, 18.015)]:
        put(f'constant/thermophysicalProperties.{phase}', f'''thermoType
{{ type heRhoThermo; mixture pureMixture; transport const; thermo hConst;
   equationOfState rhoConst; specie specie; energy sensibleInternalEnergy; }}
mixture {{ specie {{ molWeight {mw}; }} equationOfState {{ rho {density}; }}
    thermodynamics {{ Cp {cp}; Hf 0; }} transport {{ mu {mu}; Pr 1; }} }}''')
        put(f'constant/momentumTransport.{phase}', 'simulationType laminar;')
    def field(obj, dimensions, internal, kind):
        b = ''
        for p in patches:
            if closed:
                bc = 'type cyclic;'
            elif obj.startswith('U.'):
                bc = 'type pressureInletOutletVelocity; value uniform (0 0 0);' if p == 'xMax' else 'type noSlip;'
            elif obj.startswith('alpha.'):
                bc = 'type inletOutlet; inletValue uniform 0.5; value uniform 0.5;' if p == 'xMax' else 'type zeroGradient;'
            elif obj == 'p_rgh':
                bc = 'type fixedValue; value uniform 100000;' if p == 'xMax' else 'type fixedFluxPressure; value uniform 100000;'
            elif obj == 'p':
                bc = 'type calculated; value uniform 100000;'
            else:
                bc = 'type zeroGradient;'
            b += f'{p} {{ {bc} }}\n'
        put('0/'+obj, f'dimensions {dimensions};\ninternalField uniform {internal};\nboundaryField {{\n{b}}}', kind)
    for phase in ('gas','liquid'):
        field('alpha.'+phase, '[0 0 0 0 0 0 0]', '0.5', 'volScalarField')
        field('U.'+phase, '[0 1 -1 0 0 0 0]', '(0.1 0 0)' if closed else '(0 0 0)', 'volVectorField')
        field('T.'+phase, '[0 0 0 1 0 0 0]', '300', 'volScalarField')
    for obj in ('p','p_rgh'):
        field(obj, '[1 -1 -2 0 0 0 0]', '100000', 'volScalarField')
    put('system/setFieldsDict', '''defaultFieldValues (volScalarFieldValue alpha.liquid 0.2 volScalarFieldValue alpha.gas 0.8);
regions (boxToCell { box (0 0 0) (0.05 0.1 0.1);
    fieldValues (volScalarFieldValue alpha.liquid 0.8 volScalarFieldValue alpha.gas 0.2); });'''
        if closed else 'defaultFieldValues (volScalarFieldValue alpha.liquid 0.5 volScalarFieldValue alpha.gas 0.5);\nregions ();')
    return files
