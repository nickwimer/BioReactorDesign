"""Insert observation calls into the verified reference; retain native equations."""
from instrument_sc1_alpha_subcycles import _once


def solver(text):
    text=_once(text,'#include "sc1VofAudit.H"','#include "sc1VofAudit.H"\n#include "sc1VofForceAudit.H"')
    text=_once(text,'    sc1VofAudit sc1Audit(alpha1,alpha2);',
        '''    sc1VofAudit sc1Audit(alpha1,alpha2);
    if (mesh.nGeometricD()!=3)
        FatalErrorInFunction << "Force probe requires the 3D reactor mesh" << exit(FatalError);
    sc1VofForceAudit sc1ForceAudit;''')
    return _once(text,'        runTime++;','''        runTime++;
        if (runTime.value()>0.000100000001)
            FatalErrorInFunction << "Force probe is limited to two steps" << exit(FatalError);''')


def pressure(text):
    anchor='            U = HbyA + rAU()*fvc::reconstruct((phig + p_rghEqn.flux())/rAUf);'
    return _once(text,anchor,anchor+'''
            sc1ForceAudit.sample(alpha1,U,HbyA,rAU(),rAUf,gh,ghf,rho,
                mixture.surfaceTensionForce()(),phig,p_rghEqn.flux()(),p_rgh);''')
