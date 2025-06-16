{
    const vectorField& C = mesh.C();
    std::vector<double> solnvec(microbemodel::nvars);
    int O2id=0;
    double liqvol;
    double dt_react=0.01; //hours (need to have it as an input)
    double mwt_O2=32.0;
    double our_avg;

    
    //in liquid phase
    //hardcoded under assumption that i=1 references liquid.
    //absolutely no clue if that's correct or how to pull that properly
    //phaseModel& phase = fluid.multiComponentPhases()[1];
    UPtrList<volScalarField>& Y = phase.YActiveRef();
    const volScalarField& alpha2 = phase;
    //const volScalarField& rho = phase.rho();
    
    //find final OUR
    forAll(C,i)
    {
        our[i] = 0.0; 
    }

    //exit(0);
}
