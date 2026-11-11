#ifndef sc1VofForceAudit_H
#define sc1VofForceAudit_H
#include <cmath>

namespace Foam
{
// Observation only. All added fields are temporaries; no equation/input changes.
class sc1VofForceAudit
{
    label sequence_=0;

    static scalar faceValue(const surfaceScalarField& f, const label facei)
    {
        const fvMesh& mesh=f.mesh();
        if (facei<mesh.nInternalFaces()) return f[facei];
        const label patchi=mesh.boundaryMesh().whichPatch(facei);
        return f.boundaryField()[patchi][facei-mesh.boundaryMesh()[patchi].start()];
    }

    static void put(scalarField& values, const label offset, const vector& v)
    {
        values[offset]=v.x(); values[offset+1]=v.y(); values[offset+2]=v.z();
    }

    static void gather(scalarField& values)
    {
        // Exactly one elected rank contributes. Master prints complete records.
        forAll(values,i) reduce(values[i],sumOp<scalar>());
        forAll(values,i)
            if (!std::isfinite(values[i]))
                FatalErrorInFunction << "Nonfinite force-probe record" << exit(FatalError);
    }

public:
    sc1VofForceAudit()
    {
        Info<< "SC1_VOF_FORCE_PROBE version=1 maxSteps=2" << endl;
    }

    void sample
    (
        const volScalarField& alpha, const volVectorField& U,
        const volVectorField& HbyA, const volScalarField& rAU,
        const surfaceScalarField& rAUf, const volScalarField& gh,
        const surfaceScalarField& ghf, const volScalarField& rho,
        const surfaceScalarField& surfaceTension,
        const surfaceScalarField& phig, const surfaceScalarField& pressureFlux,
        const volScalarField& prgh
    )
    {
        ++sequence_;
        const fvMesh& mesh=alpha.mesh();
        if (sequence_>8)
            FatalErrorInFunction << "Force probe exceeded eight corrections" << exit(FatalError);
        const surfaceScalarField gradRho("sc1ForceGradRho",fvc::snGrad(rho));
        const surfaceScalarField interpGh("sc1ForceInterpGh",fvc::interpolate(gh));
        const surfaceScalarField gravity("sc1ForceGravity",-ghf*gradRho*mesh.magSf());
        const surfaceScalarField capillary("sc1ForceCapillary",surfaceTension*mesh.magSf());
        const surfaceScalarField pressure("sc1ForcePressure",pressureFlux/rAUf);
        const surfaceScalarField net("sc1ForceNet",(phig+pressureFlux)/rAUf);
        const surfaceScalarField freshPressure("sc1ForceFreshPressure",-fvc::snGrad(prgh)*mesh.magSf());
        // Diagnostic delta ONLY: replacing face-centre gh with interpolated gh.
        // No such change is made in the solver or pressure equation.
        const surfaceScalarField geometryDelta("sc1ForceGeometryDelta",(ghf-interpGh)*gradRho*mesh.magSf());
        const volVectorField dg("sc1DeltaGravity",rAU*fvc::reconstruct(gravity));
        const volVectorField dc("sc1DeltaCapillary",rAU*fvc::reconstruct(capillary));
        const volVectorField dp("sc1DeltaPressure",rAU*fvc::reconstruct(pressure));
        const volVectorField dn("sc1DeltaNet",rAU*fvc::reconstruct(net));
        const volVectorField df("sc1DeltaFreshPressure",rAU*fvc::reconstruct(freshPressure));
        const volVectorField dx("sc1DeltaGeometry",rAU*fvc::reconstruct(geometryDelta));
        scalar forceClosure=0, updateClosure=0; label bad=0;
        forAll(alpha,i)
        {
            const scalar a=mag(dn[i]-dg[i]-dc[i]-dp[i]);
            const scalar b=mag(U[i]-HbyA[i]-dn[i]);
            if (!std::isfinite(a) || !std::isfinite(b)) ++bad;
            else { forceClosure=max(forceClosure,a); updateClosure=max(updateClosure,b); }
        }
        reduce(bad,sumOp<label>()); reduce(forceClosure,maxOp<scalar>());
        reduce(updateClosure,maxOp<scalar>());
        Info<< "SC1_VOF_FORCE_CHECK sequence=" << sequence_ << " t_s=" << mesh.time().value()
            << " forceClosure_m_s=" << forceClosure << " updateClosure_m_s=" << updateClosure
            << " nonfinite=" << bad << endl;
        if (bad) FatalErrorInFunction << "Nonfinite force decomposition" << exit(FatalError);

        for (label phasei=0; phasei<2; ++phasei)
        {
            const word phase=phasei==0 ? "gas" : "liquid";
            scalar peak=-1; label celli=-1;
            forAll(alpha,i)
            {
                const scalar fraction=phasei==0 ? 1-alpha[i] : alpha[i];
                if (fraction>.01 && mag(U[i])>peak) { peak=mag(U[i]); celli=i; }
            }
            scalar globalPeak=peak; reduce(globalPeak,maxOp<scalar>());
            label proc=peak==globalPeak ? Pstream::myProcNo() : Pstream::nProcs();
            reduce(proc,minOp<label>());
            if (globalPeak<0) FatalErrorInFunction << "Missing present phase" << exit(FatalError);
            const bool own=Pstream::myProcNo()==proc;
            scalarField values(33,0);
            if (own)
            {
                values[0]=celli; put(values,1,mesh.C()[celli]);
                values[4]=alpha[celli]; values[5]=rho[celli]; values[6]=rAU[celli];
                put(values,7,U[celli]); put(values,10,HbyA[celli]);
                put(values,13,dg[celli]); put(values,16,dc[celli]); put(values,19,dp[celli]);
                put(values,22,dn[celli]); put(values,25,df[celli]); put(values,28,dx[celli]);
                values[31]=prgh[celli]; values[32]=gh[celli];
            }
            gather(values);
            static const char* keys[]={"localCell","x_m","y_m","z_m","alphaLiquid","rho_kg_m3","rAU",
                "Ux","Uy","Uz","Hx","Hy","Hz","gravityX","gravityY","gravityZ",
                "capillaryX","capillaryY","capillaryZ","pressureX","pressureY","pressureZ",
                "netX","netY","netZ","freshPressureX","freshPressureY","freshPressureZ",
                "geometryDeltaX","geometryDeltaY","geometryDeltaZ","p_rgh_Pa","gh_m2_s2"};
            Info<< "SC1_VOF_FORCE_CELL sequence=" << sequence_ << " t_s=" << mesh.time().value()
                << " phase=" << phase << " proc=" << proc;
            forAll(values,i) Info<< ' ' << keys[i] << '=' << values[i];
            Info<< endl;

            label nFaces=own ? mesh.cells()[celli].size() : 0;
            reduce(nFaces,sumOp<label>());
            label activeFaces=0;
            for (label j=0; j<nFaces; ++j)
            {
                scalarField v(23,0); label empty=0;
                if (own)
                {
                    const label facei=mesh.cells()[celli][j];
                    const label patchi=facei<mesh.nInternalFaces() ? -1 : mesh.boundaryMesh().whichPatch(facei);
                    empty=patchi>=0 && mesh.boundary()[patchi].type()=="empty";
                    if (!empty)
                    {
                        v[0]=facei; v[1]=mesh.faceOwner()[facei]==celli ? 1 : -1;
                        v[2]=patchi; v[3]=patchi>=0 && mesh.boundary()[patchi].coupled();
                        put(v,4,mesh.faceCentres()[facei]); put(v,7,mesh.faceAreas()[facei]);
                        if (patchi<0) put(v,10,mesh.C()[mesh.faceNeighbour()[facei]]-mesh.C()[mesh.faceOwner()[facei]]);
                        const scalar area=mag(mesh.faceAreas()[facei]);
                        v[13]=faceValue(ghf,facei); v[14]=faceValue(interpGh,facei);
                        v[15]=faceValue(gradRho,facei);
                        v[16]=faceValue(gravity,facei)/area; v[17]=faceValue(capillary,facei)/area;
                        v[18]=faceValue(pressure,facei)/area; v[19]=faceValue(net,facei)/area;
                        v[20]=faceValue(freshPressure,facei)/area; v[21]=faceValue(geometryDelta,facei)/area;
                        v[22]=faceValue(rAUf,facei);
                    }
                }
                reduce(empty,sumOp<label>());
                if (empty) continue;
                ++activeFaces;
                gather(v);
                static const char* faceKeys[]={"localFace","outwardSign","patch","coupled",
                    "CfX","CfY","CfZ","SfX","SfY","SfZ","dCX","dCY","dCZ",
                    "ghFace","ghInterpolated","snGradRho","gravity_Pa_m","capillary_Pa_m",
                    "pressure_Pa_m","net_Pa_m","freshPressure_Pa_m","geometryDelta_Pa_m","rAUf"};
                Info<< "SC1_VOF_FORCE_FACE sequence=" << sequence_ << " phase=" << phase << " proc=" << proc;
                forAll(v,i) Info<< ' ' << faceKeys[i] << '=' << v[i];
                Info<< endl;
            }
            Info<< "SC1_VOF_FORCE_FACES sequence=" << sequence_ << " phase=" << phase
                << " count=" << activeFaces << endl;
        }
    }
};
}
#endif
