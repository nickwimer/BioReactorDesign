/*---------------------------------------------------------------------------*\
  Optional phase-targeted alpha-aware Courant number for multiphaseEulerFoam.

  With no phaseCourantControl dictionary, this reproduces the stock
  multiphaseEulerFoam Courant-number calculation.

  When phaseCourantControl is enabled:
    - mixture phi may be included independently;
    - only phases named in "phases" contribute to timestep control;
    - each selected phase can have an alpha cutoff;
    - raw phase fluxes are never modified;
    - an unfiltered Courant maximum remains available for diagnostics.

  Example:

  phaseCourantControl
  {
      enabled             true;
      includeMixture      true;
      diagnostics         true;

      phases
      (
          bubbleGas
          headspaceGas
          liquid
      );

      defaultAlphaCutoff  1e-10;

      alphaCutoff
      {
          bubbleGas       1e-10;
          headspaceGas    1e-2;
          liquid          1e-10;
      }
  }

\*---------------------------------------------------------------------------*/

scalar CoNum = 0.0;
scalar meanCoNum = 0.0;
scalar rawCoNum = 0.0;

const dictionary& controlDict = runTime.controlDict();

const dictionary* phaseCoDictPtr = nullptr;

if (controlDict.found("phaseCourantControl"))
{
    phaseCoDictPtr = &controlDict.subDict("phaseCourantControl");
}

bool targetedCoEnabled = false;
bool includeMixture = true;
bool diagnostics = true;

wordList targetedPhaseNames;

scalar defaultAlphaCoCutoff = 1e-10;

const dictionary* alphaCutoffDictPtr = nullptr;


if (phaseCoDictPtr)
{
    targetedCoEnabled =
        phaseCoDictPtr->lookupOrDefault<Switch>("enabled", true);

    if (targetedCoEnabled)
    {
        includeMixture =
            phaseCoDictPtr->lookupOrDefault<Switch>
            (
                "includeMixture",
                true
            );

        diagnostics =
            phaseCoDictPtr->lookupOrDefault<Switch>
            (
                "diagnostics",
                true
            );

        defaultAlphaCoCutoff =
            phaseCoDictPtr->lookupOrDefault<scalar>
            (
                "defaultAlphaCutoff",
                1e-10
            );

        if (defaultAlphaCoCutoff < 0)
        {
            FatalErrorInFunction
                << "phaseCourantControl defaultAlphaCutoff must be >= 0; got "
                << defaultAlphaCoCutoff
                << exit(FatalError);
        }

        if (!phaseCoDictPtr->found("phases"))
        {
            FatalErrorInFunction
                << "phaseCourantControl is enabled but no 'phases' list "
                << "was supplied"
                << exit(FatalError);
        }

        phaseCoDictPtr->lookup("phases") >> targetedPhaseNames;

        if (phaseCoDictPtr->found("alphaCutoff"))
        {
            alphaCutoffDictPtr =
                &phaseCoDictPtr->subDict("alphaCutoff");
        }

        // Validate names and reject duplicates before calculating Co.
        forAll(targetedPhaseNames, namei)
        {
            const word& requestedName = targetedPhaseNames[namei];

            bool phaseFound = false;

            forAll(phases, phasei)
            {
                if (phases[phasei].name() == requestedName)
                {
                    phaseFound = true;
                    break;
                }
            }

            if (!phaseFound)
            {
                FatalErrorInFunction
                    << "phaseCourantControl requests unknown phase '"
                    << requestedName << "'"
                    << exit(FatalError);
            }

            for (label otheri = 0; otheri < namei; ++otheri)
            {
                if (targetedPhaseNames[otheri] == requestedName)
                {
                    FatalErrorInFunction
                        << "phaseCourantControl phase '" << requestedName
                        << "' appears more than once"
                        << exit(FatalError);
                }
            }
        }

        if (!includeMixture && targetedPhaseNames.empty())
        {
            FatalErrorInFunction
                << "phaseCourantControl would include neither mixture phi "
                << "nor any phase fluxes"
                << exit(FatalError);
        }
    }
}


if (mesh.nInternalFaces())
{
    scalarField mixtureSumPhi
    (
        fvc::surfaceSum(mag(phi))().primitiveField()
    );

    // Raw diagnostic always contains mixture and every phase.
    scalarField rawSumPhi(mixtureSumPhi);

    // Stock behavior when targeted control is disabled.
    // Targeted behavior starts from either mixture phi or zero.
    scalarField sumPhi(mesh.nCells(), 0.0);

    if (!targetedCoEnabled || includeMixture)
    {
        sumPhi = mixtureSumPhi;
    }

    if (targetedCoEnabled && diagnostics)
    {
        const scalar mixtureCoNum =
            0.5
           *gMax(mixtureSumPhi/mesh.V().field())
           *runTime.deltaTValue();

        Info<< "Active mixture Courant max: "
            << (includeMixture ? mixtureCoNum : scalar(0))
            << " included: " << includeMixture
            << endl;
    }

    forAll(phases, phasei)
    {
        const phaseModel& phase = phases[phasei];
        const volScalarField& alpha = phase;

        scalarField phaseSumPhi
        (
            fvc::surfaceSum(mag(phase.phi()))().primitiveField()
        );

        rawSumPhi = max(rawSumPhi, phaseSumPhi);

        bool usePhase = !targetedCoEnabled;

        if (targetedCoEnabled)
        {
            usePhase = false;

            forAll(targetedPhaseNames, namei)
            {
                if (targetedPhaseNames[namei] == phase.name())
                {
                    usePhase = true;
                    break;
                }
            }
        }

        if (!usePhase)
        {
            if (targetedCoEnabled && diagnostics)
            {
                Info<< "Active phase Courant " << phase.name()
                    << " excluded"
                    << endl;
            }

            continue;
        }

        scalar phaseAlphaCoCutoff = 0.0;

        if (targetedCoEnabled)
        {
            phaseAlphaCoCutoff = defaultAlphaCoCutoff;

            if (alphaCutoffDictPtr)
            {
                phaseAlphaCoCutoff =
                    alphaCutoffDictPtr->lookupOrDefault<scalar>
                    (
                        phase.name(),
                        defaultAlphaCoCutoff
                    );
            }

            if (phaseAlphaCoCutoff < 0)
            {
                FatalErrorInFunction
                    << "phaseCourantControl alphaCutoff for phase '"
                    << phase.name() << "' must be >= 0; got "
                    << phaseAlphaCoCutoff
                    << exit(FatalError);
            }

            const scalarField& alphaI = alpha.primitiveField();

            forAll(phaseSumPhi, celli)
            {
                if (alphaI[celli] <= phaseAlphaCoCutoff)
                {
                    phaseSumPhi[celli] = 0.0;
                }
            }
        }

        scalarField phaseRate
        (
            phaseSumPhi/mesh.V().field()
        );

        const scalar phaseCoNum =
            0.5
           *gMax(phaseRate)
           *runTime.deltaTValue();

        if (targetedCoEnabled && diagnostics)
        {
            Info<< "Active phase Courant " << phase.name()
                << " max: " << phaseCoNum
                << " alphaCutoff: " << phaseAlphaCoCutoff
                << endl;

            // Print detailed data only on a processor containing a
            // global-max controlling cell, and only for non-zero Co.
            const label localCelli = findMax(phaseRate);

            const scalar localPhaseCo =
                0.5
               *phaseRate[localCelli]
               *runTime.deltaTValue();

            if
            (
                phaseCoNum > SMALL
             && localPhaseCo >= phaseCoNum*(1.0 - 1e-12)
            )
            {
                tmp<volVectorField> tU = phase.U();
                const volVectorField& U = tU();

                Pout<< "ACTIVE_CFL_OWNER"
                    << " phase=" << phase.name()
                    << " Co=" << localPhaseCo
                    << " proc=" << Pstream::myProcNo()
                    << " cell=" << localCelli
                    << " C=" << mesh.C().primitiveField()[localCelli]
                    << " alpha=" << alpha.primitiveField()[localCelli]
                    << " magU=" << mag(U.primitiveField()[localCelli])
                    << " V=" << mesh.V().field()[localCelli]
                    << " sumAbsPhi=" << phaseSumPhi[localCelli]
                    << " alphaCutoff=" << phaseAlphaCoCutoff
                    << endl;
            }
        }

        sumPhi = max(sumPhi, phaseSumPhi);
    }

    rawCoNum =
        0.5
       *gMax(rawSumPhi/mesh.V().field())
       *runTime.deltaTValue();

    CoNum =
        0.5
       *gMax(sumPhi/mesh.V().field())
       *runTime.deltaTValue();

    meanCoNum =
        0.5
       *(gSum(sumPhi)/gSum(mesh.V().field()))
       *runTime.deltaTValue();
}


Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum << endl;

if (targetedCoEnabled)
{
    Info<< "Raw unfiltered Courant Number max: " << rawCoNum
        << " targetedPhaseCount: " << targetedPhaseNames.size()
        << " includeMixture: " << includeMixture
        << endl;
}

// ************************************************************************* //
