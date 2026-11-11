"""Observer transformations for private copies only. Never edit production sources."""
from instrument_sc1_alpha_subcycles import _once, instrument_phase_system_solve

MARKER = 'sc1TwoFieldAuditLibraryVersion'


def phase_source(source):
    result = instrument_phase_system_solve(source)
    # Reuse validated hook locations, but supply a new phase-count-independent observer.
    result = result.replace('sc1AlphaSubcycleAudit.H', 'sc1PhaseFractionAudit.H')
    result = result.replace('sc1AlphaDiagnostic', 'sc1PhaseDiagnostic')
    result = result.replace('sc1AlphaAuditLibraryVersion', MARKER)
    result = result.replace('sc1AlphaTotals', 'sc1PhaseTotals')
    result = result.replace('sc1AlphaAudit', 'sc1PhaseAudit')
    result = result.replace('sc1Alpha', 'sc1Observer')
    # Keep the shared, independently tested arithmetic enum name.
    result = result.replace('sc1ObserverBudget', 'sc1AlphaBudget')
    anchor = '                    movingPhases()[movingPhasei] *= alphaVoid/sumAlphaMoving;'
    result = _once(result, anchor, '''                    if (sc1Observer.valid())
                        sc1Observer->beforeChange(movingPhases()[movingPhasei], sc1PhaseDiagnostic::budget::rescale);
''' + anchor + '''
                    if (sc1Observer.valid())
                        sc1Observer->afterChange(movingPhases()[movingPhasei], sc1PhaseDiagnostic::budget::rescale);''')
    return result


def solver_source(source):
    if 'sc1TwoField' in source:
        raise ValueError('Two-field solver observer already present')
    source = _once(source, '#include "phaseSystem.H"', '#include "phaseSystem.H"\n'
                   '#include "sc1TwoFieldStages.H"\nextern "C" int '+MARKER+'();')
    source = _once(source, '    while (pimple.run(runTime))', '''    if (sc1TwoFieldAuditLibraryVersion()!=1)
        FatalErrorInFunction << "Wrong private two-field library" << exit(FatalError);
    Info<< "SC1_TWO_LIBRARY version=1" << endl;
    sc1TwoFieldStages sc1Stages(fluid, LTS);

    while (pimple.run(runTime))''')
    source = _once(source, '        runTime++;', '        sc1Stages.beginStep();\n        runTime++;')
    source = _once(source, '            if (!pimple.flow())', '''            sc1Stages.beginOuter();
            if (!pimple.flow())
                FatalErrorInFunction << "Two-field harness requires flow PIMPLE" << exit(FatalError);
            if (!pimple.flow())''')
    for anchor, suffix in (
        ('fluid.solve(rAUs, rAUfs);', '\n                sc1Stages.capture("alpha");'),
        ('fluid.correct();', '\n                sc1Stages.capture("modelCorrect");'),
    ):
        if source.count(anchor) != 2:
            raise ValueError('Unexpected hydrodynamic solver structure: '+anchor)
        source = source.replace(anchor, anchor+suffix)
    source = _once(source, '                fluid.correctKinematics();',
                   '                sc1Stages.capture("transportPressure");\n                fluid.correctKinematics();')
    source = _once(source, '        }\n\n        runTime.write();',
                   '            sc1Stages.capture("endOuter");\n        }\n\n'
                   '        sc1Stages.endStep();\n        runTime.write();')
    return source


def pressure_source(source):
    if 'sc1Stages' in source:
        raise ValueError('Pressure observer already present')
    source = _once(source, '        // Density variation',
                   '        sc1Stages.sourceReset(phasei);\n\n        // Density variation')
    source = _once(source,
        '            pEqnComp -= (fvModels.source(alpha, rho) & rho)/rho;',
        '''            // Evaluate the original source expression exactly once.
            const tmp<volScalarField> sc1OptionSource = fvModels.source(alpha, rho) & rho;
            sc1Stages.optionSource(phasei, sc1OptionSource());
            pEqnComp -= sc1OptionSource()/rho;''')
    source = _once(source, '            pEqnComp -= dmdts[phasei]/rho;',
                   '            sc1Stages.transferSource(phasei, dmdts[phasei]);\n'
                   '            pEqnComp -= dmdts[phasei]/rho;')
    return source
