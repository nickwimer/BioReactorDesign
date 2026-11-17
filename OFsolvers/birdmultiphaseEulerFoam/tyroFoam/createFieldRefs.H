surfaceScalarField& phi = fluid.phi();

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

phaseModel& phase = *liquidPtr;
rhoThermo& thermo = phase.thermoRef();
volScalarField& rho = thermo.rho();

const IOMRFZoneList& MRF = fluid.MRF();
Foam::fvModels& fvModels = fluid.fvModels(mesh);
Foam::fvConstraints& fvConstraints = fluid.fvConstraints(mesh);
