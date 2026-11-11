// Included once, before the time loop. This candidate is fresh-start only.
if (runTime.value()!=0 || mesh.dynamic() || mesh.nGeometricD()!=3
 || mag(g.value()-vector(0,-9.81,0))>1e-12 || p_rgh.needReference())
    FatalErrorInFunction << "Unsupported hydrostatic reference startup" << exit(FatalError);

IOdictionary sc1Reference
(
    IOobject("sc1HydrostaticReference",runTime.constant(),mesh,
        IOobject::MUST_READ,IOobject::NO_WRITE)
);
const scalar sc1Height=readScalar(sc1Reference.lookup("interfaceHeight"));
const scalar sc1Width=readScalar(sc1Reference.lookup("transitionWidth"));
if (!std::isfinite(sc1Height) || sc1Height<1 || sc1Height>1.04
 || sc1Width!=0.0125)
    FatalErrorInFunction << "Wrong SC1 hydrostatic profile" << exit(FatalError);

// The boundary extension is the original density boundary extension. Thus
// rho-rhoHydro is identically zero, including boundaries, at the reference.
const volScalarField sc1RhoHydro
(
    IOobject("rhoHydro",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),rho
);
volScalarField sc1PHydroStorage
(
    IOobject("pHydro",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),p
);
scalar sc1RhoError=0,sc1PressureError=0,sc1PotentialError=0,sc1InitialSpeed=0;
forAll(rho,i)
{
    const scalar y=mesh.C()[i].y();
    sc1RhoError=max(sc1RhoError,mag(rho[i]-(1+999*sc1Hydro::alpha(y,sc1Height,sc1Width))));
    sc1PressureError=max(sc1PressureError,mag(p[i]-sc1Hydro::pressure(y,sc1Height,sc1Width)));
    sc1PotentialError=max(sc1PotentialError,mag(gh[i]+9.81*y));
    sc1InitialSpeed=max(sc1InitialSpeed,mag(U[i]));
    if (!std::isfinite(rho[i]) || !std::isfinite(p[i]) || !std::isfinite(mag(U[i])))
        sc1PressureError=GREAT;
}
reduce(sc1RhoError,maxOp<scalar>());reduce(sc1PressureError,maxOp<scalar>());
reduce(sc1PotentialError,maxOp<scalar>());
reduce(sc1InitialSpeed,maxOp<scalar>());
if (sc1RhoError>1e-9 || sc1PressureError>1e-7 || sc1PotentialError>1e-12 || sc1InitialSpeed!=0)
    FatalErrorInFunction << "Initial fields are not the prescribed hydrostatic rest state" << exit(FatalError);
forAll(sc1PHydroStorage.boundaryField(),patchi)
    if (!mesh.boundary()[patchi].coupled())
        forAll(sc1PHydroStorage.boundaryField()[patchi],facei)
            sc1PHydroStorage.boundaryFieldRef()[patchi][facei]=sc1Hydro::pressure
            (mesh.Cf().boundaryField()[patchi][facei].y(),sc1Height,sc1Width);
const volScalarField& sc1PHydro=sc1PHydroStorage;

volScalarField pPert
(
    IOobject("pPert",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE),mesh
);
if (pPert.needReference() || gMax(mag(pPert.primitiveField()))!=0)
    FatalErrorInFunction << "Expected zero perturbation pressure with fixed vent" << exit(FatalError);
forAll(pPert.boundaryField(),patchi)
    forAll(pPert.boundaryField()[patchi],facei)
        if (pPert.boundaryField()[patchi][facei]!=0)
            FatalErrorInFunction << "Nonzero initial perturbation boundary" << exit(FatalError);
mesh.setFluxRequired(pPert.name());
Info<< "SC1_HYDROSTATIC_REFERENCE version=1 frozen=1 rhoError_kg_m3=" << sc1RhoError
    << " pressureError_Pa=" << sc1PressureError << " initialSpeed_m_s=" << sc1InitialSpeed << endl;
