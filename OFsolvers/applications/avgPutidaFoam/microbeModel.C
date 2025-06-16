#include"microbeModel.H"

namespace microbemodel
{
    const int X=0;
    const int O2=1;
    const int G=2;
    const int M=3;
    const int nvars=4;

    // Summary of parameters that need to be updated
    // y_os, Fo_max, Fs_max, K_o, kLa
    const double y_xs = 0.04804688;
    const double y_ms = 0.31582621;
    const double y_os = 0.19930839;

    const double x_max = 7.9;
    const double q_max = 4.50091315;
    const double o2_max = 0.214;

    const double K_o = 0.01384016; 
    const double K_s = 0.92049743;

    const double kLa = 60; // NEED UPDATE; approximated with current range provided 

    double X_avg=0.0;
    double G_avg=0.0;
    double M_avg=0.0;
    double O2_avg=0.0;

    void get_rhs(std::vector<double>& rhs,std::vector<double> solnvec,double t,int nvars)
    {
 
        // calculate q_s
        double F_s = solnvec[G]/(solnvec[G] + K_s);
        double F_e = solnvec[O2]/(solnvec[O2] + K_o);
        double q_s = q_max*F_s*F_e;

        // calculate final rates
        rhs[X] = y_xs*q_s*solnvec[X]*(1 - solnvec[X]/x_max);

        //set RHS of o2 to 0 as it is solved in CFD
        rhs[O2] = 0.0;
        rhs[G] = -q_s*solnvec[X];
        rhs[M] = y_ms*q_s*solnvec[X];

    }

    void advance(std::vector<double>& solnvec,int nvars,double t_now,double t_adv,double dt)
    {
        double current_time=t_now;
        double final_time=t_now+t_adv;


	//std::cout << "current time as read in by advance():\t" << current_time << "\n";
	
        std::vector<double> rhs(nvars);
        std::vector<double> solnvec_n(nvars);

        while(current_time < final_time)
        {
            current_time += dt;
	    //std::cout << "current time is:\t" << current_time << "\n";

            //at current time level n
            solnvec_n=solnvec;

            //Doing RK23

            //stage 1
            get_rhs(rhs,solnvec,current_time,nvars);
            for(int i=0;i<nvars;i++)
            {
                solnvec[i] = solnvec_n[i] + 0.5*rhs[i]*dt;
            }

            //stage 2
            get_rhs(rhs,solnvec,current_time,nvars);
            for(int i=0;i<nvars;i++)
            {
                solnvec[i] = solnvec_n[i] + rhs[i]*dt;
            }
	}
    }

    double get_our(std::vector<double> solnvec,int nvars)
    {
        // calculate q_s
        double F_s = solnvec[G]/(solnvec[G] + K_s);
        double F_e = solnvec[O2]/(solnvec[O2] + K_o);
        double q_s = q_max*F_s*F_e;

        double our = y_os*q_s*solnvec[X];

        return our;
    }
}
