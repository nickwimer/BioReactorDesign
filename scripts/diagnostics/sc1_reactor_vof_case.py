"""Single-velocity resolved-interface control: same mesh, fill and fluid properties."""
from sc1_reactor_startup_case import field
from sc1_two_field_cases import foam


def case_files():
    mu=2.414e-5*10**(247.8/(303.15-140))
    return {
        '0/U':velocity_field(),
        'constant/transportProperties':foam('transportProperties',f'''phases (liquid gas);
liquid {{ transportModel Newtonian; nu {mu/1000:.17g}; rho 1000; }}
gas {{ transportModel Newtonian; nu 1.8e-5; rho 1; }}
sigma 0.07;'''),
        'constant/momentumTransport':foam('momentumTransport','simulationType laminar;'),
        'constant/dynamicMeshDict':foam('dynamicMeshDict','dynamicFvMesh staticFvMesh;'),
        'constant/fvModels':foam('fvModels',''),
        'system/fvConstraints':foam('fvConstraints',''),
        'system/controlDict':foam('controlDict', '''application sc1InterFoamReference;
startFrom startTime; startTime 0; stopAt endTime; endTime 0.01;
deltaT 5e-05; adjustTimeStep false; maxCo 0.2; maxAlphaCo 0.2; maxDeltaT 5e-05;
writeControl timeStep; writeInterval 5; purgeWrite 0;
writeFormat ascii; writePrecision 17; writeCompression off;
timeFormat general; timePrecision 12; runTimeModifiable false;'''),
        'system/fvSchemes':foam('fvSchemes','''ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes {
    default none;
    div(rhoPhi,U) Gauss linearUpwind grad(U);
    div(phi,alpha) Gauss interfaceCompression vanLeer 1;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }'''),
        'system/fvSolution':foam('fvSolution','''solvers {
    "alpha.liquid.*" {
        nAlphaCorr 1; nAlphaSubCycles 1; MULESCorr yes; nLimiterIter 5;
        alphaApplyPrevCorr false;
        solver smoothSolver; smoother symGaussSeidel; tolerance 1e-12; relTol 0;
    }
    "pcorr.*" { solver PCG; preconditioner DIC; tolerance 1e-12; relTol 0; }
    p_rgh { solver PCG; preconditioner DIC; tolerance 1e-12; relTol 0; }
    p_rghFinal { $p_rgh; }
    "U.*" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-12; relTol 0; }
}
PIMPLE { momentumPredictor no; nOuterCorrectors 2; nCorrectors 2;
    nNonOrthogonalCorrectors 0; pRefCell 0; pRefValue 100000; }
relaxationFactors { equations { ".*" 1; } }'''),
    }


def velocity_field():
    # Reuse the exact phase-velocity boundary recipe, changing only object name.
    return field('U.liquid','[0 1 -1 0 0 0 0]','uniform (0 0 0)',True).replace('object U.liquid;', 'object U;')
