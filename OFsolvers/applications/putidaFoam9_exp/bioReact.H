
{
    const vectorField& C = mesh.C();
    std::vector<double> solnvec(microbemodel::nvars);
    std::vector<double> mwt(microbemodel::nvars);
    std::vector<int> species_id(microbemodel::nvars);
    std::vector<double> conc(microbemodel::nvars);
    int O2id=0;
    int id=0;
    double liqvol;
    double dt_react=0.01; //hours (need to have it as an input)
    double mwt_O2=32.0;
    double our_avg;

    mwt[0] = 663.43; //putida
    mwt[1] = 31.9988; //O2
    mwt[2] = 180.156; //glucose
    // mwt[3] = 142.11; //muconic acid
    mwt[3] = 140.09; //muconate
    mwt[4] = 44.00995; //CO2

    
    // phase, thermo, and rho are set in createFieldRefs.H
    UPtrList<volScalarField>& Y = phase.YActiveRef();
    const volScalarField& alpha2 = phase;

    //dictionary dict(IFstream("thermophysicalProperties.liquid")());
    
    forAll(Y,i)
    {
        Info<<"species name:\t"<<Y[i].name()<<"\n";
	microbemodel::get_sp_id(Y[i].name(), species_id, i);
	Info<<"OpenFOAM id is:\t"<<i<<"\n";

	//ThermoType s_hold(dict.subDict(Y[i].name()));
	//mwt[i] = s_hold.W();
    }

    for(i=0, i<microbemodel::nvars; i++)
      {
	Info << "species foam_id = " << species_id[i] << "\n";
	Info << "which indicates foam species: " << Y[species_id[i]].name() << "\n";
	Info << "this foam species should reflect the microbe model species matching the call order of this loop \n";
	Info << "Loop index: " << i << "\t microbe index: " << microbemodel::sp_keys[Y[species_id[i]].name()] << "\n\n";
      }
    
    Info << "Starting reaction update..\n";


    //get averaged values as mixing time scales are small compared to 
    //reactions
    liqvol = fvc::domainIntegrate(alpha2).value();
    microbemodel::X_avg  = fvc::domainIntegrate(alpha2*1000*rho*Y[species_id[microbemodel::X]]/mwt[microbemodel::X]).value()/liqvol;
    microbemodel::G_avg  = fvc::domainIntegrate(alpha2*1000*rho*Y[species_id[microbemodel::G]]/mwt[microbemodel::G]).value()/liqvol;
    microbemodel::M_avg  = fvc::domainIntegrate(alpha2*1000*rho*Y[species_id[microbemodel::M]]/mwt[microbemodel::M]).value()/liqvol;
    microbemodel::O2_avg  = fvc::domainIntegrate(alpha2*1000*rho*Y[species_id[microbemodel::O2]]/mwt[microbemodel::O2]).value()/liqvol;
  
    Info<<"Avg quantities before reaction:"<<microbemodel::X_avg<<"\t"<<microbemodel::O2<<"\t"<<microbemodel::G_avg<<"\t"<<
    microbemodel::M_avg<<"\t"<<our_avg<<"\n";
    
    //find initial OUR and advance solution
    forAll(C,i)
      {
	//Info << "Are we looping over these cells or not?\n";
	solnvec[microbemodel::X] = rho[i]*Y[species_id[microbemodel::X]][i]*1000.0/mwt[microbemodel::X]; //mol/m3
	solnvec[microbemodel::O2] = rho[i]*Y[species_id[microbemodel::O2]][i]*1000.0/mwt[microbemodel::O2]; //mol/m3
	solnvec[microbemodel::G] = rho[i]*Y[species_id[microbemodel::G]][i]*1000.0/mwt[microbemodel::G]; //mol/m3
	solnvec[microbemodel::M] = rho[i]*Y[species_id[microbemodel::M]][i]*1000.0/mwt[microbemodel::M]; //mol/m3
	solnvec[microbemodel::CO2] = rho[i]*Y[species_id[microbemodel::CO2]][i]*1000.0/mwt[microbemodel::CO2]; //mol/m3

	//Info << "before putida concentration for this cell is:\t" << solnvec[microbemodel::X] << "\n";
	//Info << "before putida mass frac for this cell is:\t" << Y[species_id[microbemodel::X]][i] << "\n";
	//Info << "O2 concentration for this cell is:\t" << solnvec[microbemodel::O2] << "\n";
	/*if(Y[species_id[microbemodel::X]][i]>0.01){
	    Info << "Location\t" << C[i] << "\n";
	    Info << "For species\t" << Y[species_id[microbemodel::X]].name() << "\n";
	    Info << "\tMass frac here is\t" << Y[species_id[microbemodel::X]][i] << "\n";
	    Info << "\tFluid density  here is\t" << rho[i] << "\n";
	    Info << "\tMolWeight here is\t" << mwt[microbemodel::X] << "\n\n";
	    }*/
	our[i]=0.5*microbemodel::get_our(solnvec,microbemodel::nvars);
	microbemodel::advance(solnvec,microbemodel::nvars,0.0,reaction_update_time.value(),dt_react);

	Y[species_id[microbemodel::X]][i] = solnvec[microbemodel::X]/rho[i]/1000.0*mwt[microbemodel::X];
	Y[species_id[microbemodel::O2]][i] = solnvec[microbemodel::O2]/rho[i]/1000.0*mwt[microbemodel::O2];
	Y[species_id[microbemodel::G]][i] = solnvec[microbemodel::G]/rho[i]/1000.0*mwt[microbemodel::G];
	Y[species_id[microbemodel::M]][i] = solnvec[microbemodel::M]/rho[i]/1000.0*mwt[microbemodel::M];
	Y[species_id[microbemodel::CO2]][i] = solnvec[microbemodel::CO2]/rho[i]/1000.0*mwt[microbemodel::CO2];

	//Info << "after putida concentration for this cell is:\t" << solnvec[microbemodel::X] << "\n";
	//Info << "after putida mass frac for this cell is:\t" << Y[species_id[microbemodel::X]][i] << "\n";
	
	our[i] += 0.5*microbemodel::get_our(solnvec,microbemodel::nvars);
    }

    our_avg=fvc::domainIntegrate(alpha2*our).value()/liqvol;

    if(Pstream::master())
    {
        //this our is the old our that is printed                                                                                                        
        os_timehist << reaction_time <<" " << microbemodel::O2_avg <<
        " " << microbemodel::X_avg <<
        " " << microbemodel::G_avg <<
        " " << microbemodel::M_avg << "\n";
	os_timehist.flush();
    }

    // set update species mass fractions and get final OUR
    /*forAll(C,i)
    {	  
      Y[species_id[microbemodel::X]][i] = solnvec[species_id[microbemodel::X]]/rho[i]/1000.0*mwt[microbemodel::X];
      Y[species_id[microbemodel::O2]][i] = solnvec[species_id[microbemodel::O2]]/rho[i]/1000.0*mwt[microbemodel::O2];
      Y[species_id[microbemodel::G]][i] = solnvec[species_id[microbemodel::G]]/rho[i]/1000.0*mwt[microbemodel::G];
      Y[species_id[microbemodel::M]][i] = solnvec[species_id[microbemodel::M]]/rho[i]/1000.0*mwt[microbemodel::M];
      Y[species_id[microbemodel::CO2]][i] = solnvec[species_id[microbemodel::CO2]]/rho[i]/1000.0*mwt[microbemodel::CO2];
      
      our[i] += 0.5*microbemodel::get_our(solnvec,microbemodel::nvars);

      if(Y[species_id[microbemodel::X]][i]>1 || Y[species_id[microbemodel::X]][i]<0){
	    Info << "Location\t" << C[i] << "\n";
	    Info << "For species\t" << Y[species_id[microbemodel::X]].name() << "\n";
	    Info << "Value here is\t" << Y[species_id[microbemodel::X]][i] << "\n\n";
	    Info << "!!!!!!!!!!NO BUENO!!!!!!!!!!!" << "\n";
	    exit(0);
	  }
	  }*/

    //check average values after reaction
    liqvol = fvc::domainIntegrate(alpha2).value();
    microbemodel::X_avg  = fvc::domainIntegrate(alpha2*1000*rho*Y[species_id[microbemodel::X]]/mwt[microbemodel::X]).value()/liqvol;
    microbemodel::G_avg  = fvc::domainIntegrate(alpha2*1000*rho*Y[species_id[microbemodel::G]]/mwt[microbemodel::G]).value()/liqvol;
    microbemodel::M_avg  = fvc::domainIntegrate(alpha2*1000*rho*Y[species_id[microbemodel::M]]/mwt[microbemodel::M]).value()/liqvol;
    microbemodel::O2_avg  = fvc::domainIntegrate(alpha2*1000*rho*Y[species_id[microbemodel::O2]]/mwt[microbemodel::O2]).value()/liqvol;
    our_avg=fvc::domainIntegrate(alpha2*our).value()/liqvol;
    
    Info<<"Avg quantities after reaction:"<<microbemodel::X_avg<<"\t"<<microbemodel::O2_avg<<"\t"<<microbemodel::G_avg<<"\t"<<
    microbemodel::M_avg<<"\t"<<our_avg<<"\n";
    //exit(0);
}
