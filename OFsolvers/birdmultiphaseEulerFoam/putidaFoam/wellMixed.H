{
    std::vector<double> solnvec(microbemodel::nvars);
    std::vector<double> totalConc(acidbasemodel::nvars - 1);
    std::vector<double> baseConc(6);
    double advance_time=5.0; //hours
    double dt_react=advance_time/50.0; //hours
    double final_time=30.0; //hours
    double pH = 7.0;
    double H = 1.0E-14;
    double H_start = 1.0E-7;
    double H_new;

    std::vector<double> MM_params(9);
    MM_params[microbemodel::MM_param_keys["y_xs"]] = y_xs.value();
    MM_params[microbemodel::MM_param_keys["y_ms"]] = y_ms.value();
    MM_params[microbemodel::MM_param_keys["y_os"]] = y_os.value();
    MM_params[microbemodel::MM_param_keys["x_max"]] = x_max.value();
    MM_params[microbemodel::MM_param_keys["q_max"]] = q_max.value();
    MM_params[microbemodel::MM_param_keys["o2_max"]] = o2_max.value();
    MM_params[microbemodel::MM_param_keys["K_o"]] = K_o.value();
    MM_params[microbemodel::MM_param_keys["K_s"]] = K_s.value();
    MM_params[microbemodel::MM_param_keys["kLa"]] = kLa.value();
    bool do_ml = bool(use_ml);

    // BR52 04112025 initial conditions (except base concentration, which is endpoint condition)
    solnvec[microbemodel::X]=0.350147052723467;  // kg/m^3
    solnvec[microbemodel::O2]=0.214;             // mol/m^3
    solnvec[microbemodel::G]=545.9946375;        // mol/m^3
    solnvec[microbemodel::M]=0.157763;           // mol/m^3

    totalConc[acidbasemodel::PhosphateBuffer]= 0.139608336; // mol/L
    totalConc[acidbasemodel::AmmoniaBase]= 0.0; //0.736052744;     // mol/L
    totalConc[acidbasemodel::MuconicAcid]= 0.0; // 0.000157763;     // mol/L
    totalConc[3]= 0.193039;                                 // mol/L
    totalConc[4]= 0.0;

    baseConc[0] = 0.0;
    baseConc[1] = 0.0;
    baseConc[2] = 0.083732977;
    baseConc[3] = 0.2;
    baseConc[4] = 0.487798268;
    baseConc[5] = 0.615927712;
    
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
    outfile << "spec ions: " << totalConc[acidbasemodel::nvars/2] << "\n";
    outfile<<"\n";
    int baseTime = 0;
    while(current_time < final_time)
    {
        current_time += advance_time;
        microbemodel::advance(solnvec,microbemodel::nvars,current_time,advance_time,dt_react, MM_params, do_ml);
        totalConc[acidbasemodel::MuconicAcid] = solnvec[microbemodel::M]/1000;
        totalConc[acidbasemodel::AmmoniaBase] = baseConc[baseTime];
        H_new = acidbasemodel::NewtonRaphson(H, totalConc, do_ml);
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
        outfile<<"\n";
	baseTime += 1;
    }
    outfile.close();
}
