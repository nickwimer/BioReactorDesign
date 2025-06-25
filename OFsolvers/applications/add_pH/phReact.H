
{
    const vectorField& C = mesh.C();
    const scalarField& Vcell = mesh.V();
    
    std::vector<double> abSystem(acidbasemodel::nvars - 1);
    std::vector<double> totalConc(acidbasemodel::nvars/2 + 2);
    std::vector<double> mwt(acidbasemodel::nvars - 1);
    std::vector<int> species_id(acidbasemodel::nvars - 1);
    int id=0;
    int i=0;
    double liqvol=0;
    double wAcidConc = 0.0;
    double wBaseConc = 0.0;
    double sAcidConc = 0.0;
    double sBaseConc = 0.0;
    double Kw = 1.0E-14;

    mwt[0] = 95.98; //HPO4
    mwt[1] = 96.98; //H2PO4
    mwt[2] = 17.031; //NH3
    mwt[3] = 18.039; //NH4
    mwt[4] = 140.09; //C6H4O4
    // mwt[5] = 142.11; //C6H6O4

    // phase, thermo, and rho are set in createFieldRefs.H
    UPtrList<volScalarField>& Y = phase.YActiveRef();
    const volScalarField& alpha2 = phase;

    //dictionary dict(IFstream("thermophysicalProperties.liquid")());
    
    forAll(Y,i)
    {
        Info<<"species name:\t"<<Y[i].name()<<"\n";
	acidbasemodel::getSystemID(Y[i].name(), species_id, i);
	Info<<"OpenFOAM id is:\t"<<i<<"\n";

	//ThermoType s_hold(dict.subDict(Y[i].name()));
	//mwt[i] = s_hold.W();
    }

    for(i=0; i<acidbasemodel::nvars-1; i++)
      {
	Info << "species foam_id = " << species_id[i] << "\n";
	Info << "which indicates foam species: " << Y[species_id[i]].name() << "\n";
	Info << "this foam species should reflect the acid-base model species matching the call order of this loop \n";
	Info << "Loop index: " << i << "\t acid/base index: " << acidbasemodel::system_keys[Y[species_id[i]].name()] << "\n\n";
      }
    
    Info << "Starting reaction update..\n";
    
    // find initial OUR and advance solution
    forAll(C,i)
      {
	// fill abSystem with mol/L = M amounts of each component present in cell i (kg/m3 = g/L)
	abSystem[acidbasemodel::PB] = rho[i]*Y[species_id[acidbasemodel::PB]][i]/mwt[acidbasemodel::PB]; //mol/L
	abSystem[acidbasemodel::PA] = rho[i]*Y[species_id[acidbasemodel::PA]][i]/mwt[acidbasemodel::PA]; //mol/L
	abSystem[acidbasemodel::NB] = rho[i]*Y[species_id[acidbasemodel::NB]][i]/mwt[acidbasemodel::NB]; //mol/L
	abSystem[acidbasemodel::NA] = rho[i]*Y[species_id[acidbasemodel::NA]][i]/mwt[acidbasemodel::NA]; //mol/L
	abSystem[acidbasemodel::MB] = rho[i]*Y[species_id[acidbasemodel::MB]][i]/mwt[acidbasemodel::MB]; //mol/L

	// get conserved totals for system components
	acidbasemodel::getTotalConc(abSystem, totalConc);

	// get spectator ions from buffer
	acidbasemodel::getSpectatorsNonBio(totalConc);

	// get spectator ions from bio components
	// this would be stuff the bio spits out that contributes to charge balance but not pH;
	// not needed for this system so fill with 0.0 instead
	// getSpectatorsBio(totalConc);
	totalConc[acidbasemodel::nvars/2 + 1] = 0.0;
	double H_set = 1.0E-15;
	double H = acidbasemodel::NewtonRaphson(H_set, totalConc);
	pH[i] = -std::log10(H);
	Info << "pH is: " << pH[i] << "\n";

	abSystem[acidbasemodel::PB] = totalConc[acidbasemodel::PhosphateBuffer] * acidbasemodel::Ka[acidbasemodel::PhosphateBuffer] / (acidbasemodel::Ka[acidbasemodel::PhosphateBuffer] + H);
	abSystem[acidbasemodel::PA] = totalConc[acidbasemodel::PhosphateBuffer] * H / (acidbasemodel::Ka[acidbasemodel::PhosphateBuffer] + H);
	abSystem[acidbasemodel::NA] = totalConc[acidbasemodel::AmmoniaBase] * H / (acidbasemodel::Ka[acidbasemodel::AmmoniaBase] + H);
	abSystem[acidbasemodel::NB] = totalConc[acidbasemodel::AmmoniaBase] - abSystem[acidbasemodel::NA];
	
	// convert updated concentrations back into massfracs for CFD 
	// liqvol = Vcell[i]*alpha2[i];
	Y[species_id[acidbasemodel::PB]][i] = abSystem[acidbasemodel::PB]/rho[i]*mwt[acidbasemodel::PB];
	Y[species_id[acidbasemodel::PA]][i] = abSystem[acidbasemodel::PA]/rho[i]*mwt[acidbasemodel::PA];
	Y[species_id[acidbasemodel::NB]][i] = abSystem[acidbasemodel::NB]/rho[i]*mwt[acidbasemodel::NB];
	Y[species_id[acidbasemodel::NA]][i] = abSystem[acidbasemodel::NA]/rho[i]*mwt[acidbasemodel::NA];
      }

    liqvol = fvc::domainIntegrate(alpha2).value();
    acidbasemodel::pH_avg  = fvc::domainIntegrate(alpha2*pH).value();
    
    if(Pstream::master())
      {                                                                                                                                             
	os_timehist << " " << acidbasemodel::pH_avg << "\n";
	os_timehist.flush();
      }

    //exit(0);
}
