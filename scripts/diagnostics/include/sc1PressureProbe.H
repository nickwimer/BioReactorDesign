#ifndef sc1PressureProbe_H
#define sc1PressureProbe_H
#include <cmath>
#include <limits>

namespace Foam
{
// Read-only field reductions. Boundary face counts include processor duplicates.
class sc1PressureProbe
{
    phaseSystem& fluid_;
    label sequence_, steps_;
    word phase_;

    struct Range
    {
        scalar lo, hi;
        label bad;
        Range() : lo(std::numeric_limits<scalar>::max()),
                  hi(-std::numeric_limits<scalar>::max()), bad(0) {}
        void add(const scalar value)
        {
            if (!std::isfinite(value)) { ++bad; return; }
            lo=min(lo,value); hi=max(hi,value);
        }
    };

    void emit(const word& term, const char* scope, Range r) const
    {
        reduce(r.lo,minOp<scalar>()); reduce(r.hi,maxOp<scalar>());
        reduce(r.bad,sumOp<label>());
        Info<< "SC1_PRESSURE_FIELD sequence=" << sequence_
            << " phase=" << phase_ << " term=" << term << " scope=" << scope
            << " min=" << r.lo << " max=" << r.hi << " nonfinite=" << r.bad << endl;
        if (r.bad)
            FatalErrorInFunction << "Nonfinite pressure-probe field " << term << exit(FatalError);
    }

public:
    explicit sc1PressureProbe(phaseSystem& fluid)
    : fluid_(fluid), sequence_(0), steps_(0), phase_("mixture")
    {
        Info<< "SC1_PRESSURE_PROBE version=1 maxCompletedSteps=2" << endl;
    }

    void stage(const word& name)
    {
        ++sequence_; phase_="mixture";
        Info<< "SC1_PRESSURE_STAGE sequence=" << sequence_
            << " t_s=" << fluid_.mesh().time().value() << " stage=" << name << endl;
    }

    void phase(const word& name) { phase_=name; }

    void array(const word& term, const scalarField& f) const
    {
        Range r; forAll(f,i) r.add(f[i]); emit(term,"cells",r);
    }

    void field(const word& term, const volScalarField& f) const
    {
        array(term,f.primitiveField());
    }

    void field(const word& term, const surfaceScalarField& f) const
    {
        Range r;
        forAll(f.primitiveField(),i) r.add(f[i]);
        forAll(f.boundaryField(),patchi)
            forAll(f.boundaryField()[patchi],i) r.add(f.boundaryField()[patchi][i]);
        emit(term,"allFacesIncludingProcessorCopies",r);
    }

    void phases()
    {
        forAll(fluid_.phases(),i)
        {
            const phaseModel& p=fluid_.phases()[i]; phase(p.name());
            field("alpha",p); field("rho",p.rho()());
            const tmp<volScalarField> speed=mag(p.U()());
            field("speed",speed()); field("phasePhi",p.phi()());
            scalar present=0, absent=0;
            forAll(p,celli)
            {
                if (p[celli]>1e-6) present=max(present,speed()[celli]);
                else absent=max(absent,speed()[celli]);
            }
            reduce(present,maxOp<scalar>()); reduce(absent,maxOp<scalar>());
            Info<< "SC1_PRESSURE_SPEED sequence=" << sequence_ << " phase=" << phase_
                << " presentThreshold=1e-6 maxPresent_m_s=" << present
                << " maxAbsent_m_s=" << absent << endl;
        }
        phase("mixture");
    }

    void pressure(const volScalarField& prgh, const volScalarField& rho,
                  const volScalarField& gh)
    {
        phase("mixture"); field("p_rghSolved",prgh);
        const tmp<volScalarField> candidate=prgh+rho*gh;
        field("pUnclipped",candidate());
        scalar minimum=std::numeric_limits<scalar>::max();
        forAll(candidate(),celli) minimum=min(minimum,candidate()[celli]);
        reduce(minimum,minOp<scalar>());
        // Local cell labels/coordinates on every rank owning an exact minimum.
        // Pout keeps processor identity; these records are evidence, not restarts.
        forAll(candidate(),celli)
        {
            if (candidate()[celli]!=minimum) continue;
            const vector& c=fluid_.mesh().C()[celli];
            Pout<< "SC1_PRESSURE_LOCATION sequence=" << sequence_
                << " localCell=" << celli << " p_Pa=" << minimum
                << " x_m=" << c.x() << " y_m=" << c.y() << " z_m=" << c.z();
            forAll(fluid_.phases(),i)
            {
                const phaseModel& p=fluid_.phases()[i];
                Pout<< " alpha_" << p.name() << "=" << p[celli]
                    << " speed_" << p.name() << "_m_s=" << mag(p.U()()[celli]);
            }
            Pout<< endl; break;
        }
        if (minimum<1000)
        {
            Info<< "SC1_PRESSURE_STOP reason=unclippedPressureBelow1000Pa min_Pa=" << minimum << endl;
            FatalErrorInFunction << "Pressure probe stopped before limiting or advancing invalid pressure"
                << exit(FatalError);
        }
    }

    bool endStep()
    {
        if (++steps_<2) return false;
        Info<< "SC1_PRESSURE_DONE completedSteps=2 meaning=diagnosticOnly" << endl;
        return true;
    }
};
}
#endif
