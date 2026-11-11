"""Replace only the private dilatation split; retain all existing audit hooks."""
from instrument_sc1_alpha_subcycles import _once
import sc1_two_field as base
import sc1_two_field_instrument as observer

APP='birdTwoFieldRates'
TAG='sc1_two_field_nonnegative_rates_v1'
MARKER='sc1TwoFieldRatesLibraryVersion'
OLD_SPLIT='''                forAll(dgdt, celli)
                {
                    if (dgdt[celli] > 0)
                    {
                        Sp[celli] -= dgdt[celli]/max(1 - alpha[celli], 1e-4);
                        Su[celli] += dgdt[celli]/max(1 - alpha[celli], 1e-4);
                    }
                    else if (dgdt[celli] < 0)
                    {
                        Sp[celli] += dgdt[celli]/max(alpha[celli], 1e-4);
                    }
                }'''


def phase_source(source):
    result=observer.phase_source(source).replace(observer.MARKER,MARKER)
    result=_once(result,'#include "sc1PhaseFractionAudit.H"',
                 '#include "sc1PhaseFractionAudit.H"\n#include "sc1TwoFieldRateSplit.H"')
    result=_once(result,'    bool dilatation = false;',
        '    sc1TwoFieldRateSplit sc1Rates(*this, sc1PhaseTotals.call, LTS, nAlphaCorr, sc1Observer.valid());\n\n'
        '    bool dilatation = false;')
    result=_once(result,OLD_SPLIT,'                sc1Rates.apply(phase, Sp, Su, dgdt);')
    return _once(result,'        tmp<volScalarField> trSubDeltaT;',
        '        sc1Rates.report();\n\n        tmp<volScalarField> trSubDeltaT;')


def solver_source(source):
    result=observer.solver_source(source).replace(observer.MARKER,MARKER)
    result=_once(result,'#include "sc1TwoFieldStages.H"',
                 '#include "sc1TwoFieldStages.H"\n#include "sc1TwoFieldRateSplit.H"')
    return _once(result,'    Info<< "SC1_TWO_LIBRARY version=1" << endl;',
        '    Info<< "SC1_RATE_LIBRARY version=1" << endl;\n'
        '    Info<< "SC1_TWO_LIBRARY version=1" << endl;')


def pressure_source(source):
    result=observer.pressure_source(source)
    result=_once(result,'            sc1Stages.optionSource(phasei, sc1OptionSource());',
        '            sc1TwoFieldRateSplit::checkOptionSource(sc1OptionSource());\n'
        '            sc1Stages.optionSource(phasei, sc1OptionSource());')
    return _once(result,'            sc1Stages.transferSource(phasei, dmdts[phasei]);',
        '            sc1TwoFieldRateSplit::checkTransfer(dmdts[phasei]);\n'
        '            sc1Stages.transferSource(phasei, dmdts[phasei]);')


RECIPE=base.BuildRecipe(TAG,APP,MARKER,phase_source,solver_source,pressure_source,
    tuple(('phaseSystems/phaseSystem/'+name,base.HERE/'include'/name)
          for name in ('sc1TwoFieldRateAlgebra.H','sc1TwoFieldRateSplit.H')))
