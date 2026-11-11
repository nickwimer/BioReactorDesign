"""Observe the existing private face-pressure solve; do not change its equations."""
from instrument_sc1_alpha_subcycles import _once


def solver(source):
    if 'sc1Probe' in source:raise ValueError('Pressure probe already present')
    source=_once(source,'#include "phaseSystem.H"',
        '#include "phaseSystem.H"\n#include "sc1PressureProbe.H"')
    source=_once(source,'    sc1TwoFieldStages sc1Stages(fluid, LTS);',
        '    sc1TwoFieldStages sc1Stages(fluid, LTS);\n'
        '    if (!faceMomentum || !partialElimination)\n'
        '        FatalErrorInFunction << "Pressure probe requires coupled face momentum" << exit(FatalError);\n'
        '    sc1PressureProbe sc1Probe(fluid);')
    return _once(source,'        sc1Stages.endStep();',
        '        sc1Stages.endStep();\n        if (sc1Probe.endStep()) return 0;')


def pressure(source):
    if 'sc1Probe' in source:raise ValueError('Pressure probe already present')
    source=_once(source,'// Face volume fractions', '''sc1Probe.stage("momentum");
sc1Probe.phases();
forAll(fluid.movingPhases(), sc1i)
{
    const phaseModel& sc1p=fluid.movingPhases()[sc1i];
    sc1Probe.phase(sc1p.name());
    sc1Probe.field("momentumA",UEqns[sc1p.index()].A()());
    sc1Probe.field("momentumHmag",mag(UEqns[sc1p.index()].H()())());
}
sc1Probe.phase("mixture");
// Face volume fractions''')
    source=_once(source,'fluid.fillFields("rAUf", dimTime/dimDensity, rAUfs);', '''fluid.fillFields("rAUf", dimTime/dimDensity, rAUfs);
forAll(phases,sc1i)
{
    sc1Probe.phase(phases[sc1i].name());
    sc1Probe.field("alphaRhoOldFace",alphaRho0fs[sc1i]);
    sc1Probe.field("rAUfPhase",rAUfs[sc1i]);
}
sc1Probe.phase("mixture");''')
    source=_once(source,'        rAUfs.set\n', '''        sc1Probe.phase(phase.name());
        sc1Probe.field("faceDiagonalExtra",AFfs[phase.index()]);
        rAUfs.set
''')
    source=_once(source,'    volScalarField rho("rho", fluid.rho());', '''    sc1Probe.stage("pressureBegin");
    sc1Probe.field("pBefore",p); sc1Probe.field("p_rghBefore",p_rgh);
    volScalarField rho("rho", fluid.rho());
    sc1Probe.field("mixtureRho",rho);''')
    source=_once(source,'    // Compressible pressure equations', '''    sc1Probe.stage("pressureInputs");
    sc1Probe.field("rAUfMixture",rAUf);
    sc1Probe.field("phiHbyAMixture",phiHbyA);
    sc1Probe.field("divPhiHbyA",fvc::div(phiHbyA)());
    forAll(phases,sc1i)
    {
        sc1Probe.phase(phases[sc1i].name());
        sc1Probe.field("alphaFace",alphafs[sc1i]);
        sc1Probe.field("alphaRAUf",alpharAUfs[sc1i]);
        sc1Probe.field("forcePhi",phigFs[sc1i]);
        sc1Probe.field("phiHbyAPhase",phiHbyAs[sc1i]);
        if (phiFfs.set(sc1i)) sc1Probe.field("extraForcePhi",phiFfs[sc1i]);
        if (phiKdPhifs.set(sc1i)) sc1Probe.field("explicitDragPhi",phiKdPhifs[sc1i]);
    }
    sc1Probe.phase("mixture");
    // Compressible pressure equations''')
    source=_once(source,'            pEqn.solve();', '''            sc1Probe.array("pressureDiagonal",pEqn.diag());
            sc1Probe.array("pressureSource",pEqn.source());
            pEqn.solve();
            sc1Probe.stage("pressureSolved");
            sc1Probe.pressure(p_rgh,rho,gh);''')
    return _once(source,'    // Update and limit the static pressure', '''    sc1Probe.stage("fluxCorrected");
    sc1Probe.phases();
    // Update and limit the static pressure''')
