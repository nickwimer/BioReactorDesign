"""Exact two-phase face-drag elimination before pressure assembly, private only."""
from instrument_sc1_alpha_subcycles import _once
import sc1_reactor_pressure_instrument as observation


def solver(source):
    source=observation.solver(source)
    source=_once(source,'#include "sc1PressureProbe.H"',
        '#include "sc1PressureProbe.H"\n#include "sc1FaceDragCoupling.H"')
    source=_once(source,'    sc1PressureProbe sc1Probe(fluid);','''    if (phases.size()!=2 || fluid.movingPhases().size()!=2
        || phases[0].name()!="gas" || phases[1].name()!="liquid"
        || !fluid.incompressible() || mesh.dynamic() || LTS)
        FatalErrorInFunction << "Unsupported pressure-coupling candidate case" << exit(FatalError);
    sc1PressureProbe sc1Probe(fluid);''')
    return _once(source,'        if (sc1Probe.endStep()) return 0;',
        '        sc1QuiescentStep(fluid);\n        sc1Probe.endStep();')


def pressure(source):
    source=observation.pressure(source)
    source=_once(source,'rAUfs.setSize(phases.size());\n{\n    PtrList<surfaceScalarField> AFfs(fluid.AFfs());',
        'rAUfs.setSize(phases.size());\n'
        'PtrList<surfaceScalarField> sc1As(phases.size());\n'
        'PtrList<surfaceScalarField> AFfs(fluid.AFfs());\n{')
    source=_once(source,'        sc1Probe.phase(phase.name());', '''        sc1As.set(phase.index(),
            (byDt(alphaRho0fs[phase.index()])
             + fvc::interpolate(UEqns[phase.index()].A())).ptr());
        sc1Probe.phase(phase.name());''')
    source=_once(source,'''                    byDt(alphaRho0fs[phase.index()])
                  + fvc::interpolate(UEqns[phase.index()].A())
                  + AFfs[phase.index()]''',
        '                    sc1As[phase.index()] + AFfs[phase.index()]')
    source=_once(source,'fluid.fillFields("rAUf", dimTime/dimDensity, rAUfs);',
        'fluid.fillFields("rAUf", dimTime/dimDensity, rAUfs);\n'
        'sc1FaceDragCoupling sc1Coupling(sc1As,AFfs,rAUfs);')
    source=_once(source,'// Explicit force fluxes', '''// Apply the same inverse to the pressure response as to the explicit RHS.
// Keep the diagonal alpha*rAUf for the original phase-force construction.
PtrList<surfaceScalarField> sc1PressureAlpha(phases.size());
forAll(phases,i)
    sc1PressureAlpha.set(i,new surfaceScalarField
        (IOobject::groupName("sc1PressureAlpha",phases[i].name()),alpharAUfs[i]));
sc1Coupling.apply(sc1PressureAlpha[0],sc1PressureAlpha[1]);

// Explicit force fluxes''')
    start=source.index('    // Predicted fluxes for each phase')
    end=source.index('    // Total predicted flux',start)
    source=source[:start]+'''    // Diagonal-normalized RHS, excluding every lagged other-phase drag flux.
    PtrList<surfaceScalarField> phiHbyAs(phases.size());
    forAll(phases,i)
    {
        phaseModel& phase=phases[i];
        phiHbyAs.set(i,
            (rAUfs[i]*(fvc::flux(UEqns[i].H())
              + alphaRho0fs[i]*byDt(phase.Uf().valid()
                  ? (mesh.Sf() & phase.Uf()().oldTime())
                  : MRF.absolute(phase.phi()().oldTime())))
              - phigFs[i]).ptr());
    }
    sc1Coupling.apply(phiHbyAs[0],phiHbyAs[1]);
    // Apply velocity boundary constraints to the coupled predictor.
    forAll(phases,i)
        phiHbyAs[i]=constrainPhiHbyA(tmp<surfaceScalarField>(phiHbyAs[i]),phases[i].U(),p_rgh);
    // Empty observer slot: no explicit drag is added to the pressure RHS.
    PtrList<surfaceScalarField> phiKdPhifs(phases.size());

'''+source[end:]
    source=_once(source,'rAUf += alphafs[phasei]*alpharAUfs[phasei];',
        'rAUf += alphafs[phasei]*sc1PressureAlpha[phasei];')
    source=_once(source,'+ alpharAUfs[phase.index()]*mSfGradp;',
        '+ sc1PressureAlpha[phase.index()]*mSfGradp;')
    start=source.index('            if (partialElimination)')
    end=source.index('            forAll(fluid.movingPhases(), movingPhasei)',
        source.index('            else',start))
    # Keep the existing relative-frame conversion loop, dropping its else scope.
    loop_end=source.index('\n            forAll(fluid.movingPhases(), movingPhasei)',end+1)
    source=source[:start]+'''            sc1CoupledFluxCheck(phi,
                (alphafs[0]*phases[0].phi()+alphafs[1]*phases[1].phi()-phi)());
            // Drag was already eliminated in BOTH pressure coefficients and RHS.
            // A second post-pressure elimination would apply it twice.
            forAll(phases,i)
            {
                MRF.makeRelative(phases[i].phiRef());
                fvc::makeRelative(phases[i].phiRef(),phases[i].U());
            }
'''+source[loop_end:]
    source=_once(source,'sc1Probe.field("alphaRAUf",alpharAUfs[sc1i]);',
        'sc1Probe.field("alphaRAUf",alpharAUfs[sc1i]);\n'
        '        sc1Probe.field("coupledPressureAlpha",sc1PressureAlpha[sc1i]);')
    return source


def observer_header(source):
    source=_once(source,'maxCompletedSteps=2','maxCompletedSteps=200')
    source=_once(source,'if (++steps_<2) return false;', 'if (++steps_<200) return false;')
    return _once(source,'completedSteps=2 meaning=diagnosticOnly',
        'completedSteps=200 meaning=requiresStartupAnalysis')
