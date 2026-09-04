
{
    const vectorField& C = mesh.C();
    const scalarField& Vcell = mesh.V();
    
    std::vector<double> abSystem(acidbasemodel::nvars); // 6 mixture components to pull; 2 phosphates, 2 ammonias, 1 butyric, 1 acetic
    std::vector<double> totalConc(acidbasemodel::nvars); // coincidentally also 6 total systems; buffer, base, 2 acids, spectators, bio-contribution
    std::vector<double> mwt(acidbasemodel::nvars);
    std::vector<double> ZConc(acidbasemodel::nvars);
    std::vector<double> activity(3);
    std::vector<int> species_id(acidbasemodel::nvars);
    int id=0;
    int i=0;
    double liqvol=0;
    double wAcidConc = 0.0;
    double wBaseConc = 0.0;
    double sAcidConc = 0.0;
    double sBaseConc = 0.0;
    double Kw = 1.0E-14;
    bool do_ml = bool(use_ml);

    mwt[0] = 95.98; //HPO4
    mwt[1] = 96.98; //H2PO4
    mwt[2] = 17.031; //NH3
    mwt[3] = 18.039; //NH4
    // mwt[4] = 59.045; //CH3COO
    mwt[4] = 60.05; //CH3COOH; acetic acid
    // mwt[6] = 87.11; //C3H7COO
    mwt[5] = 88.11; //C3H7COOH; butyric acid

    for(int i = 0; i < acidbasemodel::nvars; i++)
      {
	ZConc[i] = 0.0;
      }

    // phase, thermo, and rho are set in createFieldRefs.H
    UPtrList<volScalarField>& Y = phase.YActiveRef();
    const volScalarField& alpha2 = phase;

    // Map species names in Y[] to acid-base model indices via system_keys
    int mappedSpecies = 0;
    forAll(Y, id)
    {
        auto it = acidbasemodel::system_keys.find(Y[id].name());
        if (it == acidbasemodel::system_keys.end())
        {
            Info << "phReact: skipping non-acid-base species " << Y[id].name() << nl;
            continue;
        }

        const int acidIndex = it->second;
        species_id[acidIndex] = id;
        mappedSpecies++;

        Info << "phReact: mapped " << Y[id].name()
             << " -> modelIndex=" << acidIndex
             << " foamIndex=" << id << nl;
    }

    if (mappedSpecies != acidbasemodel::nvars)
    {
        FatalErrorInFunction
            << "Expected " << (acidbasemodel::nvars)
            << " acid-base species mapped, got " << mappedSpecies << nl
            << exit(FatalError);
    }

    for (int k = 0; k < acidbasemodel::nvars; ++k)
    {
        if (species_id[k] < 0)
        {
            FatalErrorInFunction
                << "Acid-base species index " << k << " is not mapped in YActiveRef()" << nl
                << exit(FatalError);
        }
    }
    
    Info << "Starting pH reaction update..\n";

    const scalar alphaLiquidMin = 1.0e-6;
    
    forAll(C,i)
      {

	if (alpha2[i] < alphaLiquidMin)
        {
            // No liquid in this cell, use neutral pH in gas region
	    // Neutral or ambient? Ambient/target pH for tyro is 6.3 so that's what I put here
            pH[i] = 6.3;
            continue;
        }
	
	// fill abSystem with mol/L = M amounts of each component present in cell i (kg/m3 = g/L)
	abSystem[acidbasemodel::PB] = rho[i]*Y[species_id[acidbasemodel::PB]][i]/mwt[acidbasemodel::PB]; //mol/L
	abSystem[acidbasemodel::PA] = rho[i]*Y[species_id[acidbasemodel::PA]][i]/mwt[acidbasemodel::PA]; //mol/L
	abSystem[acidbasemodel::NB] = rho[i]*Y[species_id[acidbasemodel::NB]][i]/mwt[acidbasemodel::NB]; //mol/L
	abSystem[acidbasemodel::NA] = rho[i]*Y[species_id[acidbasemodel::NA]][i]/mwt[acidbasemodel::NA]; //mol/L
	abSystem[acidbasemodel::AA] = rho[i]*Y[species_id[acidbasemodel::AA]][i]/mwt[acidbasemodel::AA]; //mol/L
	abSystem[acidbasemodel::BA] = rho[i]*Y[species_id[acidbasemodel::BA]][i]/mwt[acidbasemodel::BA]; //mol/L

	activity[0] = 1.0;
	activity[1] = 1.0;
	activity[2] = 1.0;
	
	// get conserved totals for system components
	acidbasemodel::getTotalConc(abSystem, totalConc);

	// initialize spectator ions with activity = uniform(1)
	ZConc[acidbasemodel::nvars/2 + 1] = acidbasemodel::getSpectatorsNonBio(totalConc, activity);
	double Z_hold = 1.0;
	int Z_count = 0;
	double Is = 0.0;
	// iterate over activity correction and Z+ calculation to hone in on Z+ value for starting pH = 6.3
	while( (std::abs( ZConc[acidbasemodel::nvars/2+1] - Z_hold ) >= 1.0E-06 ) && (Z_count <= 100))
	  {
	    if(Z_count > 0)
	      {
		ZConc[acidbasemodel::nvars/2 + 1] = Z_hold;
	      }
	    Is = acidbasemodel::ionicStrength(std::pow(10, -6.3), ZConc, activity, do_ml);
	    acidbasemodel::activity(Is, activity);
	    Z_hold = acidbasemodel::getSpectatorsNonBio(ZConc, activity);
	    Z_count++;
	  }

	// reset activity and fill totalConc Z+ with final value from iterative process
	for(int i = 0; i < 3; i++)
	  {
	    activity[i] = 1.0;
	  }
	totalConc[acidbasemodel::nvars/2 + 1] = Z_hold;

	if(Z_hold == 1.0)
	  {
	    Info << "Z_hold stuck at 1.0..." << "\n";
	  }
	
	// get spectator ions from bio components
	// ML proton-equivalent source from microbe model, stored in Hcell [mol/m3]
	// convert to mol/L for acid-base balance
	totalConc[acidbasemodel::nvars/2 + 2] = do_ml ? (Hcell[i]/1000.0) : 0.0;
	
	double aH_init = 1.0E-15;
	double aH_update = 1.0;
	double pH_init = -std::log10(aH_init);
	double pH_update = -std::log10(aH_update);
	int H_counter = 0;
	// double H = acidbasemodel::NewtonRaphson(H_set, totalConc, activity,  do_ml);

	while((std::abs(pH_init - pH_update) >= 1.0E-08) && (H_counter <= 100))
	  {
	    if(H_counter > 0)
	      {
		aH_init = aH_update;
		pH_init = pH_update;
	      }
	    if (aH_init < 1.0E-15)
	      {
		aH_init = 1.0E-14;
		pH_init = -std::log10(aH_init);
	      }
	    Is = acidbasemodel::ionicStrength(aH_init, totalConc, activity, do_ml);
	    acidbasemodel::activity(Is, activity);
	    aH_update = acidbasemodel::NewtonRaphson(aH_init, totalConc, activity, do_ml);
	    pH_update = -std::log10(aH_update);
	    H_counter++;
	  }
	
	pH[i] = pH_update;
	Info << "Ionic Strength is: " << Is << "\n";
	Info << "pH is: " << pH[i] << "\n";

	// update individual components for systems; future - might could only track the lumped components for each of these within the broth...
	// tracking the acids by lumped components; i.e., dissociation only considered for pH calc. 
	abSystem[acidbasemodel::PB] = totalConc[acidbasemodel::PhosphateBuffer] * acidbasemodel::Ka[acidbasemodel::PhosphateBuffer] * activity[1] \
	  / (acidbasemodel::Ka[acidbasemodel::PhosphateBuffer] * activity[1] + aH_update * activity[2]);
	abSystem[acidbasemodel::PA] = totalConc[acidbasemodel::PhosphateBuffer] * aH_update * activity[2] \
	  / (acidbasemodel::Ka[acidbasemodel::PhosphateBuffer]*activity[1] + aH_update * activity[2]);
	abSystem[acidbasemodel::NA] = totalConc[acidbasemodel::AmmoniaBase] * aH_update / (acidbasemodel::Ka[acidbasemodel::AmmoniaBase] * activity[1] + aH_update);
	abSystem[acidbasemodel::NB] = totalConc[acidbasemodel::AmmoniaBase] - abSystem[acidbasemodel::NA];

	
	// convert updated concentrations back into massfracs for CFD 
	// liqvol = Vcell[i]*alpha2[i];
	Y[species_id[acidbasemodel::PB]][i] = abSystem[acidbasemodel::PB]/rho[i]*mwt[acidbasemodel::PB];
	Y[species_id[acidbasemodel::PA]][i] = abSystem[acidbasemodel::PA]/rho[i]*mwt[acidbasemodel::PA];
	Y[species_id[acidbasemodel::NB]][i] = abSystem[acidbasemodel::NB]/rho[i]*mwt[acidbasemodel::NB];
	Y[species_id[acidbasemodel::NA]][i] = abSystem[acidbasemodel::NA]/rho[i]*mwt[acidbasemodel::NA];
      }

    Info << "pH update complete.\n";
    
    liqvol = fvc::domainIntegrate(alpha2).value();
    if (liqvol > 1e-6)
      {
	acidbasemodel::pH_avg  = fvc::domainIntegrate(alpha2*pH).value();
      }
    else
      {
	acidbasemodel::pH_avg = 6.3;
      }
    
    if(Pstream::master())
      {                                                                                                                                             
	os_timehist << " " << acidbasemodel::pH_avg << "\n";
	os_timehist.flush();
      }

    //exit(0);
}
