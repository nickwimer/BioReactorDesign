// Whole-timestep and stage inventory observer for the isolated nonreacting harness.
#ifndef sc1TwoFieldStages_H
#define sc1TwoFieldStages_H
#include "fvCFD.H"
#include "phaseSystem.H"
#include <vector>
namespace Foam {
class sc1TwoFieldStages
{
    phaseSystem& fluid_;
    const fvMesh& mesh_;
    std::vector<scalarField> alpha_, rho_, stepAlpha_, stepRho_;
    scalarField source_, transfer_;
    label outer_=0;
    scalar boundaryRate(const surfaceScalarField& f) const
    {
        scalar v=0;
        forAll(f.boundaryField(),p) if (!f.boundaryField()[p].coupled())
            forAll(f.boundaryField()[p],i) v+=f.boundaryField()[p][i];
        reduce(v,sumOp<scalar>()); return v;
    }
    void snapshot(std::vector<scalarField>& a, std::vector<scalarField>& r) const
    {
        forAll(fluid_.phases(),i)
        { a[i]=fluid_.phases()[i].primitiveField(); r[i]=fluid_.phases()[i].rho()().primitiveField(); }
    }
    scalar integral(const volScalarField& field) const
    { return gSum(field.primitiveField()*mesh_.V().field()); }
public:
    sc1TwoFieldStages(phaseSystem& fluid, bool LTS)
    : fluid_(fluid), mesh_(fluid.mesh()), alpha_(fluid.phases().size()), rho_(alpha_.size()),
      stepAlpha_(alpha_.size()), stepRho_(alpha_.size()), source_(alpha_.size(),0), transfer_(alpha_.size(),0)
    {
        if (LTS || fluid.phases().size()!=2 || fluid.stationaryPhases().size()
            || mesh_.moving() || mesh_.changing()
            || !mesh_.solverDict("alpha").lookupOrDefault<Switch>("sc1PhaseAudit",false))
            FatalErrorInFunction << "Two-field harness needs two moving phases, fixed mesh, global dt and phase audit"
                << exit(FatalError);
        snapshot(alpha_,rho_); snapshot(stepAlpha_,stepRho_);
        const unsigned precision=Info().precision(15);
        forAll(fluid_.phases(),i)
            Info<< "SC1_TWO_INITIAL phase=" << fluid_.phases()[i].name()
                << " t_s=" << mesh_.time().value()
                << " volume_m3=" << gSum(alpha_[i]*mesh_.V().field())
                << " mass_kg=" << gSum(alpha_[i]*rho_[i]*mesh_.V().field()) << endl;
        Info().precision(precision);
    }
    void beginStep()
    { snapshot(stepAlpha_,stepRho_); snapshot(alpha_,rho_); outer_=0; source_=0; transfer_=0; }
    void beginOuter()
    { ++outer_; source_=0; transfer_=0; }
    void sourceReset(label i) { source_[i]=0; transfer_[i]=0; }
    void optionSource(label i, const volScalarField& source) { source_[i]=integral(source); }
    void transferSource(label i, const volScalarField& source) { transfer_[i]=integral(source); }
    void capture(const word& stage)
    {
        const unsigned precision=Info().precision(15);
        forAll(fluid_.phases(),i)
        {
            const phaseModel& phase=fluid_.phases()[i];
            const scalarField a(phase.primitiveField()), r(phase.rho()().primitiveField());
            const scalar mass=gSum((a*r-alpha_[i]*rho_[i])*mesh_.V().field());
            const scalar byAlpha=gSum((a-alpha_[i])*rho_[i]*mesh_.V().field());
            const scalar byRho=gSum(a*(r-rho_[i])*mesh_.V().field());
            Info<< "SC1_TWO_STAGE phase=" << phase.name() << " t_s=" << mesh_.time().value()
                << " outer=" << outer_ << " stage=" << stage
                << " deltaVolume_m3=" << gSum((a-alpha_[i])*mesh_.V().field())
                << " deltaMass_kg=" << mass << " alphaPart_kg=" << byAlpha
                << " densityPart_kg=" << byRho << " identityResidual_kg=" << mass-byAlpha-byRho << endl;
            alpha_[i]=a; rho_[i]=r;
        }
        Info().precision(precision);
    }
    void endStep()
    {
        if (outer_==0 || mesh_.moving() || mesh_.changing())
            FatalErrorInFunction << "Unsupported timestep in two-field harness" << exit(FatalError);
        capture("endStep");
        const unsigned precision=Info().precision(15);
        scalarField sumAlpha(mesh_.nCells(),0);
        forAll(fluid_.phases(),i) sumAlpha+=fluid_.phases()[i].primitiveField();
        const scalar sumError=gMax(mag(sumAlpha-scalar(1)));
        const scalar dt=mesh_.time().deltaTValue();
        forAll(fluid_.phases(),i)
        {
            const phaseModel& phase=fluid_.phases()[i];
            const scalarField a(phase.primitiveField()), r(phase.rho()().primitiveField());
            const scalar mass=gSum(a*r*mesh_.V().field());
            const scalar delta=gSum((a*r-stepAlpha_[i]*stepRho_[i])*mesh_.V().field());
            const scalar boundary=boundaryRate(phase.alphaRhoPhi()());
            Info<< "SC1_TWO_STEP phase=" << phase.name() << " t_s=" << mesh_.time().value()
                << " dt_s=" << dt << " outers=" << outer_
                << " volume_m3=" << gSum(a*mesh_.V().field()) << " mass_kg=" << mass
                << " deltaMass_kg=" << delta << " boundaryOut_kg=" << dt*boundary
                << " optionSource_kg=" << dt*source_[i] << " transferSource_kg=" << dt*transfer_[i]
                << " physicalResidual_kg=" << delta+dt*(boundary-source_[i]-transfer_[i])
                << " alphaMin=" << gMin(a) << " alphaMax=" << gMax(a)
                << " phaseSumMaxError=" << sumError << endl;
        }
        Info().precision(precision);
    }
};
}
#endif
