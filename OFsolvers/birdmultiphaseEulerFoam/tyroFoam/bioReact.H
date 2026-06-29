
{
    const vectorField& C = mesh.C();
    const scalarField& Vcell = mesh.V();
    std::vector<double> solnvec(microbemodel::nvars);
    std::vector<double> mwt(microbemodel::nvars);
    std::vector<int> species_id(microbemodel::nvars);
    std::vector<double> conc(microbemodel::nvars);
    // int O2id=0;
    int NH3id=-1;
    int NH4id=-1;
    int id=0;
    int i=0;
    double liqvol;
    double dt_react=0.01; //hours (need to have it as an input)
    // double mwt_O2=32.0;
    double molarity_NH4OH = 7.398273736; // mol/L
    const scalar baseFlowRate = base_flowrate.value(); // m3/s, from constant/microbeUpdateProperties
    std::vector<double> MM_params(13); // 13 parameters for the Monod Michaelis Menten kinetics
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

    double NH3_avg;
    double NH3_avg_hold;
    double NH4_avg;
    double NH4_avg_hold;
    double delNH3;
    double delNH4;
    double NH3_conc;
    double NH4_conc;

    double total_moles_added;
    double total_new_concentration;
    double cellVol;
    double NH3_moles;
    double new_moles;
    double new_concentration;
    
    mwt[0] = 663.43; // g/mol putida
    mwt[1] = 180.156; //g/mol glucose
    mwt[2] = 88.01; // g/mol butyric acid
    mwt[3] = 60.05; // g/mol acetic acid
    // mwt[3] = 142.11; // g/mol muconic acid
    // mwt[4] = 44.00995; // g/mol CO2

    
    // phase, thermo, and rho are set in createFieldRefs.H
    UPtrList<volScalarField>& Y = phase.YActiveRef();
    const volScalarField& alpha2 = phase;

    
    forAll(Y,i)
    {
        Info<<"species name:\t"<<Y[i].name()<<"\n";
        microbemodel::get_sp_id(Y[i].name(), species_id, i);
        Info<<"OpenFOAM id is:\t"<<i<<"\n";
            if(Y[i].name()=="NH3.liquid")
            {
                NH3id=i;
            }
        if(Y[i].name()=="NH4.liquid")
            {
                NH4id=i;
            }
    }

    for(i=0; i<microbemodel::nvars; i++)
      {
        Info << "species foam_id = " << species_id[i] << "\n";
        Info << "which indicates foam species: " << Y[species_id[i]].name() << "\n";
        Info << "this foam species should reflect the microbe model species matching the call order of this loop \n";
        Info << "Loop index: " << i << "\t microbe index: " << microbemodel::sp_keys[Y[species_id[i]].name()] << "\n\n";
      }
    
    Info << "Starting bioreaction update..\n";

    //get averaged values as mixing time scales are small compared to reactions
    liqvol = fvc::domainIntegrate(alpha2).value();
    microbemodel::X_avg  = fvc::domainIntegrate(alpha2*rho*Y[species_id[microbemodel::X]]).value()/liqvol;
    microbemodel::G_avg  = fvc::domainIntegrate(alpha2*1000*rho*Y[species_id[microbemodel::G]]/mwt[microbemodel::G]).value()/liqvol;
    microbemodel::B_avg  = fvc::domainIntegrate(alpha2*1000*rho*Y[species_id[microbemodel::B]]/mwt[microbemodel::B]).value()/liqvol;
    microbemodel::A_avg  = fvc::domainIntegrate(alpha2*1000*rho*Y[species_id[microbemodel::A]]/mwt[microbemodel::A]).value()/liqvol;
  
    NH3_avg = fvc::domainIntegrate(alpha2*rho*1000*Y[NH3id]/17.031).value()/liqvol/1000; // avg mol/L 

    // get total amount of NH3 added over reaction update time: m3/s * 1000L/m3 * reaction h * 3600s/h * moles/L = moles 
    total_moles_added = baseFlowRate * 1000 * reaction_update_time.value() * 3600 * molarity_NH4OH;
    total_new_concentration = total_moles_added / liqvol / 1000; //mol/L 
  
    Info<<"\n Base flow rate (m3/s) is: "<<base_flowrate<<"\n";
    Info<<"\n Use ML for bio kinetics? "<<do_ml<<"\n";
    Info<<"\n MM parameters are: \n";
    for (const auto& param : microbemodel::MM_param_keys) {
        Info << "\t" << param.first << ": " << MM_params[param.second] << "\n";
    }
    Info<<"\n Avg quantities before and after mixing window: "<<NH3_avg<<"\t"<<total_new_concentration<< "\n\n";
    
    forAll(C,i)
      {
        // liquid density * species massFrac * kg -> g * / g/mol weight
        solnvec[microbemodel::X] = rho[i]*Y[species_id[microbemodel::X]][i]; // kg/m3
        solnvec[microbemodel::G] = rho[i]*Y[species_id[microbemodel::G]][i]*1000/mwt[microbemodel::G]; //mol/m3
        solnvec[microbemodel::B] = rho[i]*Y[species_id[microbemodel::B]][i]*1000/mwt[microbemodel::B]; //mol/m3
        solnvec[microbemodel::A] = rho[i]*Y[species_id[microbemodel::A]][i]*1000/mwt[microbemodel::A]; //mol/m3
        // Update H+ concentration in solnvec for use in ML model
        solnvec[microbemodel::H] = Hcell[i]; // mol/m3
        
        microbemodel::advance(solnvec,microbemodel::nvars,0.0,reaction_update_time.value(),dt_react, MM_params, do_ml);
        // Update solution vector with new values
        Y[species_id[microbemodel::X]][i] = solnvec[microbemodel::X]/rho[i];
        Y[species_id[microbemodel::G]][i] = solnvec[microbemodel::G]/rho[i]/1000*mwt[microbemodel::G];
        Y[species_id[microbemodel::B]][i] = solnvec[microbemodel::B]/rho[i]/1000*mwt[microbemodel::B];
        Y[species_id[microbemodel::A]][i] = solnvec[microbemodel::A]/rho[i]/1000*mwt[microbemodel::A];
        // Update H+ concentration for ML model
        Hcell[i] = solnvec[microbemodel::H]; // mol/m3
	
        // Get NH3 concentration to scale up over temporal disparity btwn reaction time and mixing time
        NH3_conc = rho[i]*Y[NH3id][i]/17.031; // mol/L
        // NH4_conc = rho[i]*Y[NH4id][i]/18.04;  // mol/L

        cellVol = Vcell[i] * alpha2[i] * 1000; // L of liquid in Cell
        // get current moles of NH3 present in mixture:
        NH3_moles = NH3_conc * cellVol; // moles in this specific cell

        new_moles = total_new_concentration * cellVol; // moles added to this specific cell
        new_concentration = (new_moles + NH3_moles) / cellVol; // mol/L
	
        //delNH3 = NH3_conc - NH3_hold[i];
        //delNH4 = NH4_conc - NH4_hold[i];

        // NH3_hold[i] = NH3_moles + 3600/fluid_update_time.value()*reaction_update_time.value()*delNH3;
        // NH4_hold[i] = NH4_avg + 3600/fluid_update_time.value()*reaction_update_time.value()*delNH4;

        // TODO: Remove hard code
        Y[NH3id][i] = new_concentration*17.031/rho[i]; 
    }

    if(Pstream::master())
    {                                                                                                        
        os_timehist << reaction_time <<" " << microbemodel::X_avg <<
        " " << microbemodel::G_avg <<
        " " << microbemodel::B_avg <<
	    " " << microbemodel::A_avg;
	os_timehist.flush();
    }

    //check average values after reaction
    liqvol = fvc::domainIntegrate(alpha2).value();
    microbemodel::X_avg  = fvc::domainIntegrate(alpha2*rho*Y[species_id[microbemodel::X]]).value()/liqvol;
    microbemodel::G_avg  = fvc::domainIntegrate(alpha2*1000*rho*Y[species_id[microbemodel::G]]/mwt[microbemodel::G]).value()/liqvol;
    microbemodel::B_avg  = fvc::domainIntegrate(alpha2*1000*rho*Y[species_id[microbemodel::B]]/mwt[microbemodel::B]).value()/liqvol;
    microbemodel::A_avg  = fvc::domainIntegrate(alpha2*1000*rho*Y[species_id[microbemodel::A]]/mwt[microbemodel::A]).value()/liqvol;
    
    Info<<"Avg quantities after reaction:"<<microbemodel::X_avg<<"\t"<<microbemodel::G_avg<<"\t"<<microbemodel::B_avg<<"\t"<<
    microbemodel::A_avg<<"\n";
    //exit(0);
}
