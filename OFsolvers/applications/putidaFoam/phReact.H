{
    const vectorField& C = mesh.C();
    const volScalarField& Vcell = mesh.V();
    
    std::vector<double> abSystem(acidbasemodel::nvars);
    std::vector<double> mwt(acidbasemodel::nvars);
    std::vector<int> species_id(acidbasemodel::nvars);
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
    mwt[5] = 142.11; //C6H6O4

    // phase, thermo, and rho are set in createFieldRefs.H
    UPtrList<volScalarField>& Y = phase.YActiveRef();
    const volScalarField& alpha2 = phase;

    //dictionary dict(IFstream("thermophysicalProperties.liquid")());
    
    forAll(Y,i)
    {
        Info<<"species name:\t"<<Y[i].name()<<"\n";
	acidbaseemodel::getSystemID(Y[i].name(), species_id, i);
	Info<<"OpenFOAM id is:\t"<<i<<"\n";

	//ThermoType s_hold(dict.subDict(Y[i].name()));
	//mwt[i] = s_hold.W();
    }

    for(i=0; i<acidbasemodel::nvars; i++)
      {
	Info << "species foam_id = " << species_id[i] << "\n";
	Info << "which indicates foam species: " << Y[species_id[i]].name() << "\n";
	Info << "this foam species should reflect the acid-base model species matching the call order of this loop \n";
	Info << "Loop index: " << i << "\t acid/base index: " << acidbasemodel::system_keys[Y[species_id[i]].name()] << "\n\n";
      }
    
    Info << "Starting reaction update..\n";
    
    //find initial OUR and advance solution
    forAll(C,i)
      {
	liqvol = Vcell[i]*alpha2[i];
	//fill abSystem with mol amounts of each component present in cell i
	abSystem[acidbasemodel::PB] = rho[i]*Y[species_id[acidbasemodel::PB]][i]*1000.0*liqvol/mwt[acidbasemodel::PB]; //mol
	abSystem[acidbasemodel::PA] = rho[i]*Y[species_id[acidbasemodel::PA]][i]*1000.0*liqvol/mwt[acidbasemodel::PA]; //mol
	abSystem[acidbasemodel::NB] = rho[i]*Y[species_id[acidbasemodel::NB]][i]*1000.0*liqvol/mwt[acidbasemodel::NB]; //mol
	abSystem[acidbasemodel::NA] = rho[i]*Y[species_id[acidbasemodel::NA]][i]*1000.0*liqvol/mwt[acidbasemodel::NA]; //mol
	abSystem[acidbasemodel::MB] = rho[i]*Y[species_id[acidbasemodel::MB]][i]*1000.0*liqvol/mwt[acidbasemodel::MB]; //mol
	abSystem[acidbasemodel::MA] = rho[i]*Y[species_id[acidbasemodel::MA]][i]*1000.0*liqvol/mwt[acidbasemodel::MA]; //mol

	// compare moles of base vs moles of acid present in cell, where each mole of muconic acid needs 2 moles of NH3 to fully react
	double leftoverCheck = 2.0*abSystem[acidbasemodel::MA] - abSystem[acidbasemodel::NB];

	// more muconic acid present than NH3 base (or full neutralization); update limited by amount of NH3 present
	if(leftoverCheck >= 0.0)
	  {
	    acidRemainderUpdate(abSystem, acidbasemodel::MuconicAcid, acidbasemodel::AmmoniaBase);
	    // react remaining acid with buffer system
	    // if no acid is leftover, pH is determined by buffer components using HH-approx.
	    if(leftoverCheck == 0.0)
	      {
		pH[i] = phBufferSystem(abSystem[acidbasemodel::PA], abSystem[acidbasemodel::PB], acidbasemodel::PhosphateBuffer);
	      }
	    // if acid is leftover, react with weak base in phosphate buffer
	    else
	      {
		double bufferCheck = 2.0*abSystem[acidbasemodel::MA] - abSystem[acidbasemodel::PB];

		// if more acid is present than the buffer can handle, reaction limited by buffer base
		// pH determined by remaining strong acid
		if(bufferCheck > 0.0)
		  {
		    acidRemainderUpdate(abSystem, acidbasemodel::MuconicAcid, acidbasemodel::PhosphateBuffer);
		    sAcidConc = abSystem[acidbasemodel::MA]/liqvol;
		    pH[i] = phStrongAcid(sAcidConc);
		  }
		// if enough buffer is present to fully neutralize muconica acid, pH determined by remaining weak acid in buffer
		else if(bufferCheck = 0.0)
		  {
		    acidRemainderUpdate(abSystem, acidbasemodel::MuconicAcid, acidbasemodel::PhosphateBuffer);
		    wAcidConc = abSystem[acidbasemodel::PA]/liqvol;
		    pH[i] = phWeakAcid(wAcidConc, acidbasemodel::PhosphateBuffer);
		  }
		// if minimal amount of muconic acid is present and weak base remains after neutralization,
		// pH approx with HH-equation for buffer systems
		else
		  {
		    baseRemainderUpdate(abSystem, acidbasemodel::MuconicAcid, acidbasemodel::PhosphateBuffer);
		    pH[i] = phBufferSystem(abSystem[acidbasemodel::PA], abSystem[acidbasemodel::PB], acidbasemodel::PhosphateBuffer);
		  }
	      }
	  }
	else
	  {
	    baseRemainderUpdate(abSystem, acidbasemodel::MuconicAcid, acidbasemodel::AmmoniaBase);
	    // react remaining base with weak acid in buffer system
	    double bufferCheck = abSystem[acidbasemodel::NB] - abSystem[acidbasemodel::PA];

	    // if more base is present than the buffer can handle, reaction limited by buffer acid
	    // pH determined by remaining "strong" base (NH3)
	    if(bufferCheck > 0.0)
	      {
		baseRemainderUpdate(abSystem, acidbasemodel::PhosphateBuffer, acidbasemodel::AmmoniaBase);
		sBaseConc = abSystem[acidbasemodel::NB]/liqvol;
		pH[i] = -log(Kw/sBaseConc);
	      }
	    // if enough buffer is present to fully neutralize NH3, pH determined by remaining weak base in buffer
	    else if(bufferCheck = 0.0)
	      {
		baseRemainderUpdate(abSystem, acidbasemodel::PhosphateBuffer, acidbasemodel::AmmoniaBase);
		wBaseConc = abSystem[acidbasemodel::PA]/liqvol;
		pH[i] = phWeakBase(wBaseConc, acidbasemodel::PhosphateBuffer);
	      }
	    // if minimal amount of NH3 is present and weak acid remains after neutralization,
	    // pH approx with HH-equation for buffer systems
	    else
	      {
		acidRemainderUpdate(abSystem, acidbasemodel::PhosphateBuffer, acidbasemodel::AmmoniaBase);
		pH[i] = phBufferSystem(abSystem[acidbasemodel::PA], abSystem[acidbasemodel::PB], acidbasemodel::PhosphateBuffer);
	      }
	  }
      }

    // update species mass fractions post-pH reactions
    forAll(C,i)
      {
	liqvol = Vcell[i]*alpha2[i];
	Y[species_id[acidbasemodel::PB]][i] = abSystem[acidbasemodel::PB]/rho[i]/1000.0/liqvol*mwt[microbemodel::PB];
	Y[species_id[acidbasemodel::PA]][i] = abSystem[acidbasemodel::PA]/rho[i]/1000.0/liqvol*mwt[microbemodel::PA];
	Y[species_id[acidbasemodel::NB]][i] = abSystem[acidbasemodel::NB]/rho[i]/1000.0/liqvol*mwt[microbemodel::NB];
	Y[species_id[acidbasemodel::NA]][i] = abSystem[acidbasemodel::NA]/rho[i]/1000.0/liqvol*mwt[microbemodel::NA];
	Y[species_id[acidbasemodel::MB]][i] = abSystem[acidbasemodel::MB]/rho[i]/1000.0/liqvol*mwt[microbemodel::MB];
	Y[species_id[acidbasemodel::MA]][i] = abSystem[acidbasemodel::MA]/rho[i]/1000.0/liqvol*mwt[microbemodel::MA];
      }

    //exit(0);
}
