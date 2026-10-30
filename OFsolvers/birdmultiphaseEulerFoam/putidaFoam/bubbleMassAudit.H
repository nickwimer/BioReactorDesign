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
    const phaseModel& auditLiquid = phases["liquid"];
    const volScalarField& auditAlphaB = auditBubble;
    const volScalarField& auditAlphaH = auditHeadspace;

    const tmp<volScalarField> auditRhoBTmp = auditBubble.rho();
    const volScalarField& auditRhoB = auditRhoBTmp();
    const scalarField& auditV = mesh.V();

    const tmp<volScalarField> auditDepTmp
    (
        fluid.dmdtf(phasePairKey("bubbleGas", "headspaceGas"))
    );
    const volScalarField& auditDep = auditDepTmp();
    const PtrList<volScalarField> auditActualSources(fluid.dmdts());

    const tmp<volScalarField> auditOrderedDepTmp
    (
        fluid.dmdtf(phasePairKey("bubbleGas", "headspaceGas", true))
    );
    const volScalarField& auditOrderedDep = auditOrderedDepTmp();

    // Independent evaluation of the configured SC1 deposition law
    // (efficiency = 1); this is NOT the solver's assembled source.
    const tmp<volScalarField> auditHeadDiameterTmp = auditHeadspace.d();
    const volScalarField& auditHeadDiameter = auditHeadDiameterTmp();
    const tmp<volVectorField> auditBubbleUTmp = auditBubble.U();
    const tmp<volVectorField> auditHeadUTmp = auditHeadspace.U();
    const volVectorField& auditBubbleVel = auditBubbleUTmp();
    const volVectorField& auditHeadVel = auditHeadUTmp();

    const tmp<volScalarField> auditO2Tmp
    (
        fluid.dmdtf(phasePairKey("bubbleGas", "liquid"))
    );
    const volScalarField& auditO2 = auditO2Tmp();

    scalar auditBubbleMass = 0;
    scalar auditDepositionLoss = 0;
    scalar auditOrderedDepLoss = 0;
    scalar auditBubblePhaseSource = 0;
    scalar auditCalculatedDepLoss = 0;
    scalar auditDissolutionLoss = 0;
    const bool auditHasBubbleSource =
        auditActualSources.set(auditBubble.index());

    forAll(auditV, celli)
    {
        auditBubbleMass +=
            auditAlphaB[celli]*auditRhoB[celli]*auditV[celli];

        // The signed pair source is negative for bubbleGas mass loss.
        auditDepositionLoss -= auditDep[celli]*auditV[celli];
        auditOrderedDepLoss -=
            auditOrderedDep[celli]*auditV[celli];
        if (auditHasBubbleSource)
        {
            auditBubblePhaseSource +=
                auditActualSources[auditBubble.index()][celli]*auditV[celli];
        }
        if (auditHeadDiameter[celli] > SMALL)
        {
            auditCalculatedDepLoss +=
                1.5*auditRhoB[celli]*auditAlphaB[celli]
               *auditAlphaH[celli]
               *mag(auditBubbleVel[celli] - auditHeadVel[celli])
               /auditHeadDiameter[celli]*auditV[celli];
        }

        auditDissolutionLoss -= auditO2[celli]*auditV[celli];
    }

    reduce(auditBubbleMass, sumOp<scalar>());
    reduce(auditDepositionLoss, sumOp<scalar>());
    reduce(auditOrderedDepLoss, sumOp<scalar>());
    reduce(auditBubblePhaseSource, sumOp<scalar>());
    reduce(auditCalculatedDepLoss, sumOp<scalar>());
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

    const tmp<surfaceScalarField> auditLiquidFluxTmp
    (
        auditLiquid.alphaRhoPhi()
    );
    const surfaceScalarField& auditLiquidFlux = auditLiquidFluxTmp();

    scalar auditBubbleBoundaryOut = 0;
    scalar auditHeadBoundaryOut = 0;
    scalar auditLiquidBoundaryOut = 0;
    scalar auditBubbleBoundaryIn = 0;
    scalar auditHeadBoundaryIn = 0;
    scalar auditLiquidBoundaryIn = 0;

    forAll(auditBubbleFlux.boundaryField(), patchi)
    {
        const fvsPatchScalarField& bubblePatch =
            auditBubbleFlux.boundaryField()[patchi];
        const fvsPatchScalarField& headPatch =
            auditHeadFlux.boundaryField()[patchi];
        const fvsPatchScalarField& liquidPatch =
            auditLiquidFlux.boundaryField()[patchi];

        if (bubblePatch.coupled())
        {
            continue;
        }

        forAll(bubblePatch, facei)
        {
            auditBubbleBoundaryOut += bubblePatch[facei];
            auditHeadBoundaryOut += headPatch[facei];
            auditLiquidBoundaryOut += liquidPatch[facei];

            auditBubbleBoundaryIn += min(bubblePatch[facei], scalar(0));
            auditHeadBoundaryIn += min(headPatch[facei], scalar(0));
            auditLiquidBoundaryIn += min(liquidPatch[facei], scalar(0));
        }
    }

    reduce(auditBubbleBoundaryOut, sumOp<scalar>());
    reduce(auditHeadBoundaryOut, sumOp<scalar>());
    reduce(auditLiquidBoundaryOut, sumOp<scalar>());
    reduce(auditBubbleBoundaryIn, sumOp<scalar>());
    reduce(auditHeadBoundaryIn, sumOp<scalar>());
    reduce(auditLiquidBoundaryIn, sumOp<scalar>());

    // Preserve a history across included calls, so the inventory rate
    // is the finite difference between consecutive audit write times.
    static bool auditHavePrevious = false;
    static scalar auditPreviousTime = 0;
    static scalar auditPreviousMass = 0;

    const scalar auditTime = runTime.value();
    scalar auditDMdt = 0;
    scalar auditClosure = 0;
    scalar auditAssembledClosure = 0;

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
        // This uses the actual sum of all the phase-system's
        // interphase transfer terms, bypassing pair-key lookup errors.
        auditAssembledClosure =
            auditDMdt
          - (
                auditInjection
              + auditBubblePhaseSource
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
            << " orderedDepLoss_kg_s=" << auditOrderedDepLoss
            << " assembledBubbleDmdt_kg_s=" << auditBubblePhaseSource
            << " assembledBubbleDmdtPresent=" << auditHasBubbleSource
            << " correlationDepLoss_kg_s=" << auditCalculatedDepLoss
            << " oxygenDissolutionLoss_kg_s=" << auditDissolutionLoss
            << " bubbleBoundaryOut_kg_s=" << auditBubbleBoundaryOut
            << " headspaceBoundaryOut_kg_s=" << auditHeadBoundaryOut
            << " liquidBoundaryOut_kg_s=" << auditLiquidBoundaryOut
            << " bubbleBoundaryIn_kg_s=" << auditBubbleBoundaryIn
            << " headspaceBoundaryIn_kg_s=" << auditHeadBoundaryIn
            << " liquidBoundaryIn_kg_s=" << auditLiquidBoundaryIn;

        if (auditHavePrevious && auditTime > auditPreviousTime)
        {
            Info<< " inventoryRate_kg_s=" << auditDMdt
                << " endpointClosure_kg_s=" << auditClosure
                << " assembledEndpointClosure_kg_s="
                << auditAssembledClosure;
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
    const volVectorField& auditBubbleU =
        mesh.lookupObject<volVectorField>("U.bubbleGas");

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
        const fvsPatchScalarField& auditLFluxPatch =
            auditLiquidFlux.boundaryField()[patchi];

        const fvPatchScalarField& auditLAlphaPatch =
            auditLiquid.boundaryField()[patchi];
        const fvPatchScalarField& auditBAlphaPatch =
            auditBubble.boundaryField()[patchi];
        const fvPatchScalarField& auditHAlphaPatch =
            auditHeadspace.boundaryField()[patchi];

        scalar auditPatchArea = 0;
        scalar auditPressureAreaSum = 0;
        scalar auditPrghAreaSum = 0;
        scalar auditPatchBubbleOut = 0;
        scalar auditPatchHeadOut = 0;
        scalar auditPatchLiquidOut = 0;
        scalar auditPatchBubbleIn = 0;
        scalar auditPatchHeadIn = 0;
        scalar auditPatchLiquidIn = 0;
        scalar auditPatchLiquidAlphaArea = 0;
        scalar auditPatchBubbleAlphaArea = 0;
        scalar auditPatchHeadAlphaArea = 0;

        forAll(auditAreas, facei)
        {
            auditPatchArea += auditAreas[facei];
            auditPressureAreaSum +=
                auditAreas[facei]*auditPPatch[facei];
            auditPrghAreaSum +=
                auditAreas[facei]*auditPrghPatch[facei];
            auditPatchBubbleOut += auditBFluxPatch[facei];
            auditPatchHeadOut += auditHFluxPatch[facei];
            auditPatchLiquidOut += auditLFluxPatch[facei];
            auditPatchBubbleIn += min(auditBFluxPatch[facei], scalar(0));
            auditPatchHeadIn += min(auditHFluxPatch[facei], scalar(0));
            auditPatchLiquidIn += min(auditLFluxPatch[facei], scalar(0));
            auditPatchLiquidAlphaArea +=
                auditAreas[facei]*auditLAlphaPatch[facei];
            auditPatchBubbleAlphaArea +=
                auditAreas[facei]*auditBAlphaPatch[facei];
            auditPatchHeadAlphaArea +=
                auditAreas[facei]*auditHAlphaPatch[facei];
        }

        reduce(auditPatchArea, sumOp<scalar>());
        reduce(auditPressureAreaSum, sumOp<scalar>());
        reduce(auditPrghAreaSum, sumOp<scalar>());
        reduce(auditPatchBubbleOut, sumOp<scalar>());
        reduce(auditPatchHeadOut, sumOp<scalar>());
        reduce(auditPatchLiquidOut, sumOp<scalar>());
        reduce(auditPatchBubbleIn, sumOp<scalar>());
        reduce(auditPatchHeadIn, sumOp<scalar>());
        reduce(auditPatchLiquidIn, sumOp<scalar>());
        reduce(auditPatchLiquidAlphaArea, sumOp<scalar>());
        reduce(auditPatchBubbleAlphaArea, sumOp<scalar>());
        reduce(auditPatchHeadAlphaArea, sumOp<scalar>());

        if (Pstream::master() && auditPatchArea > SMALL)
        {
            Info<< "HEADSPACE_PRESSURE_AUDIT"
                << " t_s=" << auditTime
                << " patch=" << mesh.boundary()[patchi].name()
                << " pRghBC=" << auditPrghPatch.type()
                << " headspaceUBC="
                << auditHeadU.boundaryField()[patchi].type()
                << " bubbleUBC="
                << auditBubbleU.boundaryField()[patchi].type()
                << " liquidAlphaBC=" << auditLAlphaPatch.type()
                << " bubbleAlphaBC=" << auditBAlphaPatch.type()
                << " headspaceAlphaBC=" << auditHAlphaPatch.type()
                << " area_m2=" << auditPatchArea
                << " meanP_Pa="
                << auditPressureAreaSum/auditPatchArea
                << " meanPrgh_Pa="
                << auditPrghAreaSum/auditPatchArea
                << " meanAlphaLiquid="
                << auditPatchLiquidAlphaArea/auditPatchArea
                << " meanAlphaBubble="
                << auditPatchBubbleAlphaArea/auditPatchArea
                << " meanAlphaHeadspace="
                << auditPatchHeadAlphaArea/auditPatchArea
                << " bubbleOut_kg_s=" << auditPatchBubbleOut
                << " headspaceOut_kg_s=" << auditPatchHeadOut
                << " liquidOut_kg_s=" << auditPatchLiquidOut
                << " bubbleIn_kg_s=" << auditPatchBubbleIn
                << " headspaceIn_kg_s=" << auditPatchHeadIn
                << " liquidIn_kg_s=" << auditPatchLiquidIn
                << endl;
        }
    }

    auditHavePrevious = true;
    auditPreviousTime = auditTime;
    auditPreviousMass = auditBubbleMass;
}
