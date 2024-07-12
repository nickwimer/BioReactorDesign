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
volScalarField microbe_conc
(
 IOobject
 (
  "microbe_conc",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 initial_microbe_conc
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
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
 ),
 mesh,
 dimensionedScalar("our", dimMoles/dimVolume/dimTime, 0)
);

/*double liqvol = fvc::domainIntegrate(alpha2).value();
microbemodel::X_avg  = fvc::domainIntegrate(alpha2*microbe_conc).value()/liqvol;
microbemodel::G_avg  = fvc::domainIntegrate(alpha2*glucose).value()/liqvol;
microbemodel::M_avg  = fvc::domainIntegrate(alpha2*muconate).value()/liqvol;*/

//=================================================================
