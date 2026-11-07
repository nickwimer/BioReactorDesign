#!/usr/bin/env python3
"""Pure source transformation for a PRIVATE copy of phaseSystemSolve.C.

This module never writes files or invokes a build. The private-library preparer
owns destinations, dependency checks and library linkage. Original equations
and source expressions are retained byte-for-byte between inserted hooks.
"""


def _once(source, anchor, replacement):
    if source.count(anchor) != 1:
        raise ValueError(f"SC1 alpha instrumentation expected one anchor: {anchor!r}")
    return source.replace(anchor, replacement)


def instrument_phase_system_solve(source):
    if "sc1Alpha" in source:
        raise ValueError("SC1 alpha audit already/partially installed")
    source = _once(source, '#include "phaseSystem.H"',
                   '#include "phaseSystem.H"\n#include "sc1AlphaSubcycleAudit.H"')
    source = _once(source, 'void Foam::phaseSystem::solve',
                   'extern "C" int sc1AlphaAuditLibraryVersion() { return 1; }\n\n'
                   'void Foam::phaseSystem::solve')
    anchor = '    const bool LTS = fv::localEulerDdt::enabled(mesh_);'
    source = _once(source, anchor, anchor + '''

    static sc1AlphaDiagnostic::Totals sc1AlphaTotals;
    autoPtr<sc1AlphaDiagnostic::Audit> sc1Alpha;
    if (alphaControls.lookupOrDefault<Switch>("sc1AlphaAudit", false))
    {
        sc1Alpha.reset(new sc1AlphaDiagnostic::Audit
            (*this, sc1AlphaTotals, LTS, nAlphaCorr));
    }''')
    source = _once(source, '            if (dilatation)', '''            if (sc1Alpha.valid())
                sc1Alpha->originalSource(phase, Sus[phasei]);

            if (dilatation)''')
    source = _once(source, '                volScalarField::Internal& Sp = Sps[phasei];',
                   '''                if (sc1Alpha.valid())
                    sc1Alpha->rawDilatation(phase, dgdt);

                volScalarField::Internal& Sp = Sps[phasei];''')
    anchor = '''                MULES::explicitSolve
                (
                    geometricOneField(),
                    alpha,
                    alphaPhi,
                    Sps[phase.index()],
                    Sus[phase.index()]
                );'''
    source = _once(source, anchor, '''                if (sc1Alpha.valid()) sc1Alpha->beforeMules(phase);
''' + anchor + '''
                if (sc1Alpha.valid())
                    sc1Alpha->afterMules
                        (phase, alphaPhi, Sps[phase.index()], Sus[phase.index()]);''')
    anchor = '                    alphaEqn.solve();'
    source = _once(source, anchor, '''                    if (sc1Alpha.valid())
                        sc1Alpha->beforeChange(phase, sc1AlphaBudget::implicitCorrection);
''' + anchor + '''
                    if (sc1Alpha.valid())
                        sc1Alpha->afterImplicit(phase, alphaEqn);''')
    anchor = '''                volScalarField& referenceAlpha = *referencePhasePtr;
                referenceAlpha = alphaVoid;

                forAll(solvePhases, solvePhasei)
                {
                    referenceAlpha -= solvePhases[solvePhasei];
                }'''
    source = _once(source, anchor, '''                if (sc1Alpha.valid())
                    sc1Alpha->beforeChange(*referencePhasePtr, sc1AlphaBudget::referenceSubcycle);
''' + anchor + '''
                if (sc1Alpha.valid())
                    sc1Alpha->afterChange(*referencePhasePtr, sc1AlphaBudget::referenceSubcycle);''')
    anchor = '            phase.maxMin(0, 1);'
    source = _once(source, anchor, '''            if (sc1Alpha.valid())
                sc1Alpha->beforeChange(phase, sc1AlphaBudget::clamp);
''' + anchor + '''
            if (sc1Alpha.valid())
                sc1Alpha->afterChange(phase, sc1AlphaBudget::clamp);''')
    anchor = '''            volScalarField& referenceAlpha = referencePhase;
            referenceAlpha = alphaVoid;

            forAll(solvePhases, solvePhasei)
            {
                referenceAlpha -= solvePhases[solvePhasei];
            }'''
    source = _once(source, anchor, '''            if (sc1Alpha.valid())
                sc1Alpha->beforeChange(referencePhase, sc1AlphaBudget::referenceFinal);
''' + anchor + '''
            if (sc1Alpha.valid())
                sc1Alpha->afterChange(referencePhase, sc1AlphaBudget::referenceFinal);''')
    anchor = '\n}\n\n\n// ************************************************************************* //'
    source = _once(source, anchor,
                   '\n    if (sc1Alpha.valid()) sc1Alpha->finish();' + anchor)
    return source
