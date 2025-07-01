{
    std::vector<double> solnvec(microbemodel::nvars);
    std::vector<double> totalConc(acidbasemodel::nvars - 1);
    std::vector<double> baseConc(13);
    double advance_time=4.0; //hours
    double dt_react=advance_time/50.0; //hours
    double final_time=36.0; //hours
    double pH = 7.0;
    double H = 1.0E-14;
    double H_start = 1.0E-7;
    double H_new;

    // BR52 04112025 initial conditions (except base concentration, which is endpoint condition)
    solnvec[microbemodel::X]=0.350147052723467;  // kg/m^3
    solnvec[microbemodel::O2]=0.05;               // mol/m^3
    solnvec[microbemodel::G]=545.9946375;        // mol/m^3
    solnvec[microbemodel::M]=0.157763;           // mol/m^3

    totalConc[acidbasemodel::PhosphateBuffer]= 0.139608336; // mol/L
    totalConc[acidbasemodel::AmmoniaBase]= 0.0; //0.736052744;     // mol/L
    totalConc[acidbasemodel::MuconicAcid]= 0.0; // 0.000157763;     // mol/L
    totalConc[3]= 0.193039;                                 // mol/L
    totalConc[4]= 0.0;

    // Base conc roughly interpolated from BR52 04112025 (U.B.)
    baseConc[0] = 0.0;
    baseConc[1] = 0.0;
    baseConc[2] = 0.079665546;
    baseConc[3] = 0.083732977;
    baseConc[4] = 0.356565343;
    baseConc[5] = 0.4;
    baseConc[6] = 0.487798268;
    baseConc[7] = 0.615927712;
    baseConc[8] = 0.698626333;
    baseConc[9] = 0.702196826;
    baseConc[10] = 0.736052744;
    
    std::ofstream outfile("wellmixed.dat");
    
    double current_time=0.0;
    Info << "current time: " << current_time << "\n";
    for(int i=0;i<microbemodel::nvars;i++)
    {
        Info << solnvec[i] <<"\t";
    }
    Info << "\n";
    for(int i=0;i<acidbasemodel::nvars-1;i++)
    {
        Info <<totalConc[i]<<"\t";
    }
    Info << "\n";
    outfile<<current_time<<"\t";
    for(int i=0;i<microbemodel::nvars;i++)
    {
        outfile<<solnvec[i]<<"\t";
    }
    for(int i=0;i<acidbasemodel::nvars-1;i++)
    {
        outfile<<totalConc[i]<<"\t";
    }
    outfile<<-std::log10(H_start)<<"\t";
    outfile<<"\n";
    acidbasemodel::getSpectatorsNonBio(totalConc);
    // outfile << "spec ions: " << totalConc[acidbasemodel::nvars/2] << "\n";
    // outfile<<"\n";
    int baseTime = 0;
    while(current_time < final_time)
    {
        current_time += advance_time;
        microbemodel::advance(solnvec,microbemodel::nvars,current_time,advance_time,dt_react);
	double tauC = microbemodel::get_tauC(solnvec, microbemodel::nvars);
	totalConc[acidbasemodel::MuconicAcid] = solnvec[microbemodel::M]/1000;
	totalConc[acidbasemodel::AmmoniaBase] = baseConc[baseTime];
	H_new = acidbasemodel::NewtonRaphson(H, totalConc);
	pH = -std::log10(H_new);
	//H = H_new;
        outfile<<current_time<<"\t";
        for(int i=0;i<microbemodel::nvars;i++)
        {
            outfile<<solnvec[i]<<"\t";
        }
	for(int i=0;i<acidbasemodel::nvars-1;i++)
	{
	    outfile<<totalConc[i]<<"\t";
	}
	outfile<<pH<<"\t";
	outfile<<tauC<<"\t";
        outfile<<"\n";
	baseTime += 1;
    }
    outfile.close();
}
