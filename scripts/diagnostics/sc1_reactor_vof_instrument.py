"""Only insert read-only audit calls in the unmodified Foundation OF9 solver."""
from instrument_sc1_alpha_subcycles import _once


def solver(text):
    text=_once(text,'#include "fvcSmooth.H"', '#include "fvcSmooth.H"\n#include "sc1VofAudit.H"')
    text=_once(text,'    turbulence->validate();', '''    turbulence->validate();
    if (mesh.dynamic() || LTS || rho1.value()!=1000 || rho2.value()!=1
     || alpha1.name()!="alpha.liquid" || alpha2.name()!="alpha.gas"
     || !isType<twoPhaseChangeModels::noPhaseChange>(phaseChange))
        FatalErrorInFunction << "Unsupported VOF reference contract" << exit(FatalError);
    sc1VofAudit sc1Audit(alpha1,alpha2);''')
    text=_once(text,'            #include "alphaEqnSubCycle.H"', '''            if (nAlphaSubCycles!=1 || nAlphaCorr!=1 || alphaApplyPrevCorr)
                FatalErrorInFunction << "Unsupported VOF alpha settings" << exit(FatalError);
            #include "alphaEqnSubCycle.H"
            sc1Audit.alphaFlux(alphaPhi10,phi);''')
    return _once(text,'        runTime.write();',
        '        sc1Audit.step(alpha1,alpha2,U,phi);\n        runTime.write();')


def pressure(text):
    return _once(text,'        p_rghEqn.solve();',
        '        p_rghEqn.solve();\n        sc1Audit.pressure(p_rgh,rho,gh);')
