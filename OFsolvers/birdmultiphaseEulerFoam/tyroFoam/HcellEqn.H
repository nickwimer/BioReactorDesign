phaseModel* liquidPtr = nullptr;

forAll(phases, phasei)
{
    if (phases[phasei].name() == "liquid")
    {
        liquidPtr = &phases[phasei];
        break;
    }
}

if (!liquidPtr)
{
    FatalErrorInFunction
        << "Could not find phase named 'liquid'" << nl
        << "Available phases are: ";

    forAll(phases, phasei)
    {
        FatalError << phases[phasei].name() << " ";
    }

    FatalError << exit(FatalError);
}

phaseModel& liquid = *liquidPtr;
const volScalarField& alphaLiquid = liquid;
const surfaceScalarField& alphaPhiLiquid = liquid.alphaPhi();

fvScalarMatrix HcellEqn
(
    fvm::ddt(alphaLiquid, Hcell)
  + fvm::div(alphaPhiLiquid, Hcell)
  - fvm::laplacian(alphaLiquid*DHcell, Hcell)
);

HcellEqn.relax();
HcellEqn.solve();
Hcell.correctBoundaryConditions();