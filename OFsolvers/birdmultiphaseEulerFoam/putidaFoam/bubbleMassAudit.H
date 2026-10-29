// Diagnostic-only bubble-gas mass budget for a three-phase putidaFoam run.
//
// Included after the PIMPLE loop and before runTime.write(), only when
// PIMPLE/bubbleMassAudit is enabled and runTime.writeTime() is true.
//
// All rates are signed in kg/s.  Positive boundary flux is outward.
// Positive transfer losses remove mass from bubbleGas.  This is a
// write-interval, endpoint-rate budget, NOT a discrete solver residual.
{
    const phaseModel& auditBubble = phases["bubbleGas"];
    const phaseModel& auditHeadspace = phases["headspaceGas"];
    const volScalarField& auditAlphaB = auditBubble;

    const tmp<volScalarField> auditRhoBTmp = auditBubble.rho();
    const volScalarField& auditRhoB = auditRhoBTmp();
    const scalarField& auditV = mesh.V();

    const tmp<volScalarField> auditDepTmp
    (
        fluid.dmdtf(phasePairKey("bubbleGas", "headspaceGas", true))
    );
    const volScalarField& auditDep = auditDepTmp();

    const tmp<volScalarField> auditO2Tmp
    (
        fluid.dmdtf(phasePairKey("bubbleGas", "liquid"))
    );
    const volScalarField& auditO2 = auditO2Tmp();

    scalar auditBubbleMass = 0;
    scalar auditDepositionLoss = 0;
    scalar auditDissolutionLoss = 0;

    forAll(auditV, celli)
    {
        auditBubbleMass +=
            auditAlphaB[celli]*auditRhoB[celli]*auditV[celli];

        // The signed pair source is negative for bubbleGas mass loss.
        auditDepositionLoss -= auditDep[celli]*auditV[celli];
        auditDissolutionLoss -= auditO2[celli]*auditV[celli];
    }

    reduce(auditBubbleMass, sumOp<scalar>());
    reduce(auditDepositionLoss, sumOp<scalar>());
    reduce(auditDissolutionLoss, sumOp<scalar>());

    // Re-evaluate the fvModels source to the bubbleGas density equation.
    // For the SC1 case this should be the sum of 15 sparger mass sources.
    // It is an endpoint source, not a recorded time-integral.
    scalar auditInjection = 0;
    if (fvModels.addsSupToField(auditRhoB.name()))
    {
        const tmp<volScalarField> auditInjectionTmp
        (
            fvModels.source(auditAlphaB, auditRhoB) & auditRhoB
        );
        const volScalarField& auditInjectionField = auditInjectionTmp();

        forAll(auditV, celli)
        {
            auditInjection += auditInjectionField[celli]*auditV[celli];
        }
    }
    reduce(auditInjection, sumOp<scalar>());

    // Use the phase mass flux actually assembled by the phase solver,
    // excluding coupled processor/cyclic patch faces.
    const tmp<surfaceScalarField> auditBubbleFluxTmp
    (
        auditBubble.alphaRhoPhi()
    );
    const surfaceScalarField& auditBubbleFlux = auditBubbleFluxTmp();

    const tmp<surfaceScalarField> auditHeadFluxTmp
    (
        auditHeadspace.alphaRhoPhi()
    );
    const surfaceScalarField& auditHeadFlux = auditHeadFluxTmp();

    scalar auditBubbleBoundaryOut = 0;
    scalar auditHeadBoundaryOut = 0;

    forAll(auditBubbleFlux.boundaryField(), patchi)
    {
        const fvsPatchScalarField& bubblePatch =
            auditBubbleFlux.boundaryField()[patchi];
        const fvsPatchScalarField& headPatch =
            auditHeadFlux.boundaryField()[patchi];

        if (bubblePatch.coupled())
        {
            continue;
        }

        forAll(bubblePatch, facei)
        {
            auditBubbleBoundaryOut += bubblePatch[facei];
            auditHeadBoundaryOut += headPatch[facei];
        }
    }

    reduce(auditBubbleBoundaryOut, sumOp<scalar>());
    reduce(auditHeadBoundaryOut, sumOp<scalar>());

    // Preserve a history across included calls, so the inventory rate
    // is the finite difference between consecutive audit write times.
    static bool auditHavePrevious = false;
    static scalar auditPreviousTime = 0;
    static scalar auditPreviousMass = 0;

    const scalar auditTime = runTime.value();
    scalar auditDMdt = 0;
    scalar auditClosure = 0;

    if (auditHavePrevious && auditTime > auditPreviousTime)
    {
        auditDMdt =
            (auditBubbleMass - auditPreviousMass)
           /(auditTime - auditPreviousTime);

        // + injection, - phase transfer, - outward boundary transport.
        auditClosure =
            auditDMdt
          - (
                auditInjection
              - auditDepositionLoss
              - auditDissolutionLoss
              - auditBubbleBoundaryOut
            );
    }

    if (Pstream::master())
    {
        Info<< "BUBBLE_MASS_AUDIT"
            << " t_s=" << auditTime
            << " inventory_kg=" << auditBubbleMass
            << " injection_kg_s=" << auditInjection
            << " depositionLoss_kg_s=" << auditDepositionLoss
            << " oxygenDissolutionLoss_kg_s=" << auditDissolutionLoss
            << " bubbleBoundaryOut_kg_s=" << auditBubbleBoundaryOut
            << " headspaceBoundaryOut_kg_s=" << auditHeadBoundaryOut;

        if (auditHavePrevious && auditTime > auditPreviousTime)
        {
            Info<< " inventoryRate_kg_s=" << auditDMdt
                << " endpointClosure_kg_s=" << auditClosure;
        }
        else
        {
            Info<< " inventoryRate_kg_s=NA"
                << " endpointClosure_kg_s=NA";
        }

        Info<< endl;
    }

    // Audit the pressure boundary on each physical patch, rather than
    // assuming that the top pressure patch is called "outlet".
    const volVectorField& auditHeadU =
        mesh.lookupObject<volVectorField>("U.headspaceGas");

    forAll(mesh.boundary(), patchi)
    {
        const fvPatchScalarField& auditPrghPatch =
            p_rgh.boundaryField()[patchi];

        if (auditPrghPatch.coupled())
        {
            continue;
        }

        const fvPatchScalarField& auditPPatch =
            p.boundaryField()[patchi];

        const fvsPatchScalarField& auditAreas =
            mesh.magSf().boundaryField()[patchi];

        const fvsPatchScalarField& auditBFluxPatch =
            auditBubbleFlux.boundaryField()[patchi];
        const fvsPatchScalarField& auditHFluxPatch =
            auditHeadFlux.boundaryField()[patchi];

        scalar auditPatchArea = 0;
        scalar auditPressureAreaSum = 0;
        scalar auditPrghAreaSum = 0;
        scalar auditPatchBubbleOut = 0;
        scalar auditPatchHeadOut = 0;

        forAll(auditAreas, facei)
        {
            auditPatchArea += auditAreas[facei];
            auditPressureAreaSum +=
                auditAreas[facei]*auditPPatch[facei];
            auditPrghAreaSum +=
                auditAreas[facei]*auditPrghPatch[facei];
            auditPatchBubbleOut += auditBFluxPatch[facei];
            auditPatchHeadOut += auditHFluxPatch[facei];
        }

        reduce(auditPatchArea, sumOp<scalar>());
        reduce(auditPressureAreaSum, sumOp<scalar>());
        reduce(auditPrghAreaSum, sumOp<scalar>());
        reduce(auditPatchBubbleOut, sumOp<scalar>());
        reduce(auditPatchHeadOut, sumOp<scalar>());

        if (Pstream::master() && auditPatchArea > SMALL)
        {
            Info<< "HEADSPACE_PRESSURE_AUDIT"
                << " t_s=" << auditTime
                << " patch=" << mesh.boundary()[patchi].name()
                << " pRghBC=" << auditPrghPatch.type()
                << " headspaceUBC="
                << auditHeadU.boundaryField()[patchi].type()
                << " area_m2=" << auditPatchArea
                << " meanP_Pa="
                << auditPressureAreaSum/auditPatchArea
                << " meanPrgh_Pa="
                << auditPrghAreaSum/auditPatchArea
                << " bubbleOut_kg_s=" << auditPatchBubbleOut
                << " headspaceOut_kg_s=" << auditPatchHeadOut
                << endl;
        }
    }

    auditHavePrevious = true;
    auditPreviousTime = auditTime;
    auditPreviousMass = auditBubbleMass;
}
