surfaceScalarField& phi = fluid.phi();

phaseModel& phase = fluid.multiComponentPhases()[1];
rhoThermo& thermo = phase.thermoRef();
volScalarField& rho = thermo.rho();

const IOMRFZoneList& MRF = fluid.MRF();
Foam::fvModels& fvModels = fluid.fvModels(mesh);
Foam::fvConstraints& fvConstraints = fluid.fvConstraints(mesh);
