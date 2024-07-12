{
    const vectorField& C = mesh.C();
    std::vector<double> solnvec(microbemodel::nvars);
    std::vector<double> mwt(microbemodel::nvars);
    std::vector<int> species_id(microbemodel::nvars);
    int O2id=0;
    int id=0;
    double liqvol;
    double dt_react=0.01; //hours (need to have it as an input)
    double mwt_O2=32.0;
    double our_avg;

    
    //in liquid phase
    //hardcoded i=1 references liquid.
    phaseModel& phase = fluid.multiComponentPhases()[1];    
    UPtrList<volScalarField>& Y = phase.YActiveRef();
    const volScalarField& rho = phase.rho();
    const volScalarField& Wt = phase.thermoRef().W();

    Info << "Y? \t" << Y << "\n";
    Info << "weight? \t" << Wt << "\n";
    
    forAll(Y,i)
    {
        Info<<"species names:\t"<<Y[i].name()<<"\n";
	id = microbemodel::get_sp_id(Y[i].name());
	species_id[i] = id;
	Info<<"species id:\t"<<species_id[i]<<"\n";
	Info<<"molecular weight?\t"<< Wt[i] <<"\n";
	mwt[i] = 10.0+i*5;
        if(Y[i].name()=="O2.liquid")
        {
            O2id=i;
        }
    }

    Info << "Starting reaction update..\n";
    Info << "O2 id is = "<<O2id<<"\tspecname:"<<Y[O2id].name()<<"\n";


    //get averaged values as mixing time scales are small compared to 
    //reactions
/*
    liqvol = fvc::domainIntegrate(alpha2).value();
    microbemodel::X_avg  = fvc::domainIntegrate(alpha2*microbe_conc).value()/liqvol;
    microbemodel::G_avg  = fvc::domainIntegrate(alpha2*glucose).value()/liqvol;
    microbemodel::M_avg  = fvc::domainIntegrate(alpha2*muconate).value()/liqvol;
*/  
    
    //find initial OUR and advance solution
    forAll(C,i)
      {
      forAll(Y,j)
	{
	  solnvec[species_id[j]] = rho[i]*Y[j][i]*1000.0/Wt[j]; //mol/m3
	}
      our[i]=0.5*microbemodel::get_our(solnvec,microbemodel::nvars);
      microbemodel::advance(solnvec,microbemodel::nvars,0.0,reaction_update_time.value(),dt_react);
    }

    // set update species mass fractions and get final OUR
    forAll(C,i)
    {
      forAll(Y,j)
	{
	  
	  Y[j][i] = solnvec[species_id[j]]/rho[i]/1000.0*Wt[j];
	}
      our[i] += 0.5*microbemodel::get_our(solnvec,microbemodel::nvars);
    }
   
    //exit(0);
}
