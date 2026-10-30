// Diagnostic-only: compare the assembled phase-system deposition source
// with the deposition correlation at the SAME PIMPLE state.
// Include immediately after fluid.correct(), before pressure/velocity update.
// Do not modify any state fields or call fluid.correct() again.
{
    const phaseModel& checkBubble = phases["bubbleGas"];
    const phaseModel& checkHead = phases["headspaceGas"];
    const volScalarField& checkAlphaB = checkBubble;
    const volScalarField& checkAlphaH = checkHead;
    const scalarField& checkV = mesh.V();

    const tmp<volScalarField> checkRhoTmp = checkBubble.rho();
    const volScalarField& checkRho = checkRhoTmp();
    const tmp<volScalarField> checkDTmp = checkHead.d();
    const volScalarField& checkD = checkDTmp();
    const tmp<volVectorField> checkUbTmp = checkBubble.U();
    const tmp<volVectorField> checkUhTmp = checkHead.U();
    const volVectorField& checkUb = checkUbTmp();
    const volVectorField& checkUh = checkUhTmp();

    const tmp<volScalarField> checkDmdtfTmp
    (
        fluid.dmdtf(phasePairKey("bubbleGas", "headspaceGas"))
    );
    const volScalarField& checkDmdtf = checkDmdtfTmp();

    // This is the reference formula in the checked-in BiRD deposition.C:
    // efficiency = 1; headspace diameter from the active phase model.
    scalar checkModelDepLoss = 0;
    scalar checkCorrelationDepLoss = 0;
    scalar checkWeightedSlip = 0;
    scalar checkOverlapWeight = 0;

    forAll(checkV, celli)
    {
        checkModelDepLoss -= checkDmdtf[celli]*checkV[celli];

        if (checkD[celli] > SMALL)
        {
            checkCorrelationDepLoss +=
                1.5*checkRho[celli]*checkAlphaB[celli]
               *checkAlphaH[celli]
               *mag(checkUb[celli] - checkUh[celli])
               /checkD[celli]*checkV[celli];
        }
        const scalar overlap =
            max(checkAlphaB[celli], scalar(0))
           *max(checkAlphaH[celli], scalar(0))*checkV[celli];
        checkOverlapWeight += overlap;
        checkWeightedSlip += overlap*mag(checkUb[celli] - checkUh[celli]);
    }

    reduce(checkModelDepLoss, sumOp<scalar>());
    reduce(checkCorrelationDepLoss, sumOp<scalar>());
    reduce(checkWeightedSlip, sumOp<scalar>());
    reduce(checkOverlapWeight, sumOp<scalar>());

    static label checkStageNumber = 0;
    ++checkStageNumber;

    if (Pstream::master())
    {
        Info<< "BUBBLE_STAGE_AUDIT"
            << " t_s=" << runTime.value()
            << " stage=" << checkStageNumber
            << " depositionLoss_kg_s=" << checkModelDepLoss
            << " sameStateCorrelationLoss_kg_s="
            << checkCorrelationDepLoss
            << " overlapSlip_m_s="
            << (checkOverlapWeight > SMALL
                ? checkWeightedSlip/checkOverlapWeight : 0)
            << endl;
    }
}
