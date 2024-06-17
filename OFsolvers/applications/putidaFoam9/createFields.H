#include "createRDeltaT.H"
#include "readGravitationalAcceleration.H"
#include "readhRef.H"

Info<< "Creating phaseSystem\n" << endl;

autoPtr<twoPhaseSystem> fluidPtr
(
    twoPhaseSystem::New(mesh)
);
twoPhaseSystem& fluid = fluidPtr();

dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    fluid
);

#include "gh.H"

volScalarField& p = fluid.phase1().thermoRef().p();

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

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);
mesh.setFluxRequired(p_rgh.name());

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

dimensionedScalar reaction_update_time
(   "reaction_update_time",
     dimTime,
     reaction_props.lookup("reaction_update_time"));

dimensionedScalar fluid_update_time
(   "fluid_update_time",
     dimTime,
     reaction_props.lookup("fluid_update_time"));


dimensionedScalar init_microbe_conc
(   "initial_microbe_conc",
     dimMoles/dimVolume,
     microbe_props.lookup("putidaLiquid"));

dimensionedScalar init_glucose
(   "initial_glucose",
     dimMoles/dimVolume,
     microbe_props.lookup("glLiquid"));

dimensiondScalar init_muconate
(   "initial_muconate",
    dimMoles/dimVolume,
    microbe_props.lookup("maLiquid"));

//=================================================================
//
//Create volume scalars for bioeractor fields======================
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
 init_microbe_conc
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
 init_glucose
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
 init_muconate)
 );

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
