{
    autoPtr<phaseSystem::specieTransferTable>
        specieTransferPtr(fluid.specieTransfer());

    phaseSystem::specieTransferTable&
        specieTransfer(specieTransferPtr());

    fluid.correctReactions();

    forAll(fluid.multiComponentPhases(), multiComponentPhasei)
    {
        phaseModel& phase = fluid.multiComponentPhases()[multiComponentPhasei];

        UPtrList<volScalarField>& Y = phase.YActiveRef();
        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();

        forAll(Y, i)
        {
	  // Info<< "Accessing species:\t" << Y[i].name() << "\n";
	  // Info<< "where mass fracs SHOULD be given by:\t" << Y[i] << "\n";
            const tmp<fvScalarMatrix> sc1ExternalSource
            (
                fvModels.source(alpha, rho, Y[i])
            );
            fvScalarMatrix YiEqn
            (
                phase.YiEqn(Y[i])
             ==
               *specieTransfer[Y[i].name()]
              + sc1ExternalSource()
            );

            if
            (
                sc1Conservation.valid() && pimple.finalPimpleIter()
             && Y[i].member() == "O2"
            )
            {
                sc1Conservation->equation
                (
                    phase, Y[i], YiEqn,
                    *specieTransfer[Y[i].name()], sc1ExternalSource()
                );
            }
            YiEqn.relax();
            YiEqn.solve("Yi");
        }
    }

    if (sc1Conservation.valid() && pimple.finalPimpleIter())
        sc1Conservation->beforeSpeciesCorrection();
    fluid.correctSpecies();
    if (sc1Conservation.valid() && pimple.finalPimpleIter())
        sc1Conservation->afterSpeciesCorrection();
}
