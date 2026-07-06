{
    std::vector<double> solnvec(microbemodel::nvars);
    std::vector<double> totalConc(acidbasemodel::nvars);
    std::vector<double> baseConc(6);
    double advance_time=5.0; //hours
    double dt_react=advance_time/50.0; //hours
    double final_time=30.0; //hours
    double pH = 7.0;
    double H = 1.0E-14;
    double H_start = 1.0E-7;
    double H_new;

    std::vector<double> MM_params(13);
    MM_params[microbemodel::MM_param_keys["q_max"]] = q_max.value();
    MM_params[microbemodel::MM_param_keys["K_s"]] = K_s.value();
    MM_params[microbemodel::MM_param_keys["K_i"]] = K_i.value();
    MM_params[microbemodel::MM_param_keys["P_d"]] = P_d.value();
    MM_params[microbemodel::MM_param_keys["m_i"]] = m_i.value();
    MM_params[microbemodel::MM_param_keys["alpha_ba"]] = alpha_ba.value();
    MM_params[microbemodel::MM_param_keys["beta_ba"]] = beta_ba.value();
    MM_params[microbemodel::MM_param_keys["alpha_aa"]] = alpha_aa.value();
    MM_params[microbemodel::MM_param_keys["beta_aa"]] = beta_aa.value();
    MM_params[microbemodel::MM_param_keys["Y_x"]] = Y_x.value();
    MM_params[microbemodel::MM_param_keys["Y_ba"]] = Y_ba.value();
    MM_params[microbemodel::MM_param_keys["Y_aa"]] = Y_aa.value();
    MM_params[microbemodel::MM_param_keys["m_s"]] = m_s.value();
    bool do_ml = bool(use_ml);

    // BR52 04112025 initial conditions (except base concentration, which is endpoint condition)
    solnvec[microbemodel::X]=0.350147052723467;  // kg/m^3
    solnvec[microbemodel::G]=545.9946375;        // mol/m^3
    solnvec[microbemodel::B]=193.039;            // mol/m^3
    solnvec[microbemodel::A]=0.0;                // mol/m^3
    solnvec[microbemodel::H]=H_start*1000.0;     // mol/m^3

    totalConc[acidbasemodel::PhosphateBuffer]= 0.139608336; // mol/L
    totalConc[acidbasemodel::AmmoniaBase]= 0.0; //0.736052744;     // mol/L
    totalConc[acidbasemodel::AceticAcid]= solnvec[microbemodel::A]/1000.0; // mol/L
    totalConc[acidbasemodel::ButyricAcid]= solnvec[microbemodel::B]/1000.0; // mol/L
    totalConc[acidbasemodel::nvars/2 + 1]= 0.0; // spectator ions
    totalConc[acidbasemodel::nvars/2 + 2]= 0.0; // bio proton-equivalent source
    acidbasemodel::getSpectatorsNonBio(totalConc);

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
    for(int i=0;i<acidbasemodel::nvars;i++)
    {
        Info <<totalConc[i]<<"\t";
    }
    Info << "\n";
    outfile<<current_time<<"\t";
    for(int i=0;i<microbemodel::nvars;i++)
    {
        outfile<<solnvec[i]<<"\t";
    }
    for(int i=0;i<acidbasemodel::nvars;i++)
    {
        outfile<<totalConc[i]<<"\t";
    }
    outfile<<-std::log10(H_start)<<"\t";
    outfile<<"\n";
    outfile << "spec ions: " << totalConc[acidbasemodel::nvars/2] << "\n";
    outfile<<"\n";
    int baseTime = 0;
    while(current_time < final_time)
    {
        current_time += advance_time;
        solnvec[microbemodel::H] = H*1000.0;
        microbemodel::advance(solnvec,microbemodel::nvars,current_time,advance_time,dt_react, MM_params, do_ml);
        totalConc[acidbasemodel::AceticAcid] = solnvec[microbemodel::A]/1000.0;
        totalConc[acidbasemodel::ButyricAcid] = solnvec[microbemodel::B]/1000.0;
        totalConc[acidbasemodel::AmmoniaBase] = baseConc[baseTime];
        totalConc[acidbasemodel::nvars/2 + 2] = do_ml ? (solnvec[microbemodel::H]/1000.0) : 0.0;
        H_new = acidbasemodel::NewtonRaphson(H, totalConc, do_ml);
        H = H_new;
        pH = -std::log10(H_new);
        outfile<<current_time<<"\t";
        for(int i=0;i<microbemodel::nvars;i++)
        {
            outfile<<solnvec[i]<<"\t";
        }
	for(int i=0;i<acidbasemodel::nvars;i++)
	{
	    outfile<<totalConc[i]<<"\t";
	}
	outfile<<pH<<"\t";
        outfile<<"\n";
	baseTime += 1;
    }
    outfile.close();
}
