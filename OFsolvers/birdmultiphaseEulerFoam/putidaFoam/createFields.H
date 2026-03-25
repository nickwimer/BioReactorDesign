#include "createRDeltaT.H"
#include "readGravitationalAcceleration.H"
#include "readhRef.H"

Info<< "Creating phaseSystem\n" << endl;

autoPtr<phaseSystem> fluidPtr
(
    phaseSystem::New(mesh)
);
phaseSystem& fluid = fluidPtr();
phaseSystem::phaseModelList& phases = fluid.phases();

#include "gh.H"

volScalarField& p = phases[0].thermoRef().p();

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

if (fluid.found("pMin"))
  {
    IOWarningInFunction(fluid)
      << "Pressure limits, pMin and pMax, are now read from "
      << pimple.dict().name() << endl;
  }

pressureReference pressureReference
(
    p,
    p_rgh,
    pimple.dict(),
    fluid.incompressible()
);

if (fluid.incompressible())
  {
    p = p_rgh + fluid.rho()*gh;
  }

if (p_rgh.needReference() && fluid.incompressible())
{
    p += dimensionedScalar
    (
        "p",
	p.dimensions(),
        pressureReference.refValue()
      - getRefCellValue(p, pressureReference.refCell())
    );
}

p_rgh = p - fluid.rho()*gh;

mesh.setFluxRequired(p_rgh.name());

PtrList<volScalarField> rAUs;
PtrList<surfaceScalarField> rAUfs;

//Read microbe update properties============================
IOdictionary reaction_props
(
 IOobject
 (
  "microbeUpdateProperties",
  runTime.constant(),
  mesh,
  IOobject::MUST_READ_IF_MODIFIED,
  IOobject::NO_WRITE
 )
);

/*
IOdictionary microbe_props
(
 IOobject
 (
  "globalTransport",
  runTime.constant(),
  mesh,
  IOobject::MUST_READ_IF_MODIFIED,
  IOobject::NO_WRITE
 )
);
*/

dimensionedScalar reaction_update_time
(   "reaction_update_time",
     dimTime,
     reaction_props.lookup("reaction_update_time"));

dimensionedScalar fluid_update_time
(   "fluid_update_time",
     dimTime,
     reaction_props.lookup("fluid_update_time"));

dimensionedScalar ph_update_time
(   "ph_update_time",
     dimTime,
     reaction_props.lookup("ph_update_time"));

// base flow rate into reactor (m3/s)
dimensionedScalar base_flowrate
(   "base_flowrate",
     dimVolume/dimTime,
     reaction_props.lookup("base_flowrate"));

// boolean for whether to use ML model for reaction updates or not
Switch use_ml(false);
reaction_props.readIfPresent("use_ml", use_ml);
// MM parameters to pass into microbe model
dimensionedScalar y_xs
(   "y_xs",
     dimless,
     reaction_props.lookup("y_xs"));
dimensionedScalar y_ms
(   "y_ms",
     dimless,
     reaction_props.lookup("y_ms"));
dimensionedScalar y_os
(   "y_os",
     dimless,     
     reaction_props.lookup("y_os"));
dimensionedScalar x_max
(   "x_max",
     dimMoles/dimMass,
     reaction_props.lookup("x_max"));
dimensionedScalar q_max
(   "q_max",
     dimMoles/dimMass/dimTime,
     reaction_props.lookup("q_max"));
dimensionedScalar o2_max
(   "o2_max",
     dimMoles/dimVolume,
     reaction_props.lookup("o2_max"));
dimensionedScalar K_o
(   "K_o",
     dimMoles/dimVolume,
     reaction_props.lookup("K_o"));
dimensionedScalar K_s
(   "K_s",
     dimMoles/dimVolume,
     reaction_props.lookup("K_s"));
dimensionedScalar kLa
(   "kLa",
     dimless/dimTime,
     reaction_props.lookup("kLa"));


/*
scalar initial_microbe_conc = microbe_props.lookup<scalar>("putidaLiquid");

scalar initial_glucose = microbe_props.lookup<scalar>("glLiquid");

scalar initial_muconate = microbe_props.lookup<scalar>("maLiquid");
*/
/*
dimensionedScalar init_glucose
(   "initial_glucose",
     dimMoles/dimVolume,
     microbe_props.lookup("glLiquid"));

dimensionedScalar init_muconate
(   "initial_muconate",
    dimMoles/dimVolume,
    microbe_props.lookup("maLiquid"));
*/
//=================================================================
//
//Create volume scalars for bioeractor fields======================
/*
volScalarField pH
(
 IOobject
 (
  "pH",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 
 );

volScalarField glucose
(
 IOobject
 (
  "glucose",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 initial_glucose
 );

volScalarField muconate
(
 IOobject
 (
  "muconate",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 initial_muconate
 );
*/
volScalarField our
(
 IOobject
 (
  "our",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("our", dimMoles/dimVolume/dimTime, 0)
);

volScalarField pH
(
 IOobject
 (
  "pH",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
 ),
 mesh,
 7.0
 );

volScalarField NH3_hold
(
 IOobject
 (
  "NH3_hold",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
 ),
 mesh,
 0.0
 );

volScalarField NH4_hold
(
 IOobject
 (
  "NH4_hold",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
 ),
 mesh,
 0.0
);

/*double liqvol = fvc::domainIntegrate(alpha2).value();
microbemodel::X_avg  = fvc::domainIntegrate(alpha2*microbe_conc).value()/liqvol;
microbemodel::G_avg  = fvc::domainIntegrate(alpha2*glucose).value()/liqvol;
microbemodel::M_avg  = fvc::domainIntegrate(alpha2*muconate).value()/liqvol;*/

//=================================================================
