"""Pressure perturbation about a fixed, independently checked hydrostatic profile."""
import re
from instrument_sc1_alpha_subcycles import _once

BUOYANCY='fvc::interpolate(rho-sc1RhoHydro)*(g & mesh.Sf())/mesh.magSf()'


def solver(text):
    text=_once(text,'#include "sc1VofAudit.H"',
        '#include "sc1VofAudit.H"\n#include "sc1HydrostaticProfile.H"')
    text=_once(text,'    sc1VofAudit sc1Audit(alpha1,alpha2);',
        '    sc1VofAudit sc1Audit(alpha1,alpha2);\n    #include "sc1CreateHydrostatic.H"')
    return _once(text,'        runTime++;','''        runTime++;
        if (runTime.value()>0.010000000001)
            FatalErrorInFunction << "Hydrostatic candidate limited to 10 ms" << exit(FatalError);''')


def momentum(text):
    text=_once(text,'- ghf*fvc::snGrad(rho)','+ '+BUOYANCY)
    return _once(text,'fvc::snGrad(p_rgh)','fvc::snGrad(pPert)')


def pressure(text):
    text=_once(text,'- ghf*fvc::snGrad(rho)','+ '+BUOYANCY)
    text=re.sub(r'\bp_rgh\b','pPert',text)
    text=_once(text,'        sc1Audit.pressure(pPert,rho,gh);','''        p == sc1PHydro + pPert;
        p_rgh == p - rho*gh;
        sc1Audit.pressure(p_rgh,rho,gh);''')
    text=_once(text,'    p == pPert + rho*gh;',
        '    p == sc1PHydro + pPert;\n    p_rgh == p - rho*gh;')
    # The vent fixes perturbation pressure; a reference-cell shift is forbidden.
    start=text.index('    if (pPert.needReference())',text.index('    p == sc1PHydro'))
    end=text.index('    if (!correctPhi)',start)
    text=text[:start]+'''    if (pPert.needReference())
        FatalErrorInFunction << "Lost fixed perturbation-pressure vent" << exit(FatalError);

'''+text[end:]
    return text
