#include"microbeModel.H"
#include<map>

namespace microbemodel
{
    const int X=0;
    const int O2=1;
    const int G=2;
    const int M=3;
    const int CO2=4;
    const int nvars=5;
    std::map<std::string, int> sp_keys = {{"putida.liquid", 0},
						  {"O2.liquid", 1},
						  {"glucose.liquid", 2},
						  {"muconate.liquid", 3},
						  {"CO2.liquid", 4}};

    // Summary of parameters that need to be updated
    // y_os, Fo_max, Fs_max, K_o, kLa
    const double y_xs = 0.109;
    const double y_ms = 0.3;
    const double y_os = 0.0467; // NEED UPDATE

    const double x_max = 7.9;
    const double Fo_max = 1; // NEED UPDATE
    const double Fs_max = 1.1; // NEED UPDATE
    const double o2_max = 0.214;

    const double K_o = 0.0214; // NEED UPDATE 
    const double K_s = 0.92;

    const double kLa = 50; // NEED UPDATE; approximated with current range provided 

    double X_avg=0.0;
    double G_avg=0.0;
    double M_avg=0.0;
    double O2_avg=0.0;

  int get_sp_id(std::string name)
  {
    return sp_keys[name];
  }

  
    void get_rhs(std::vector<double>& rhs,std::vector<double> solnvec,double t,int nvars)
    {
 
        // calculate q_s
        double F_s = Fs_max*solnvec[G]/(solnvec[G] + K_s);
        double F_e = Fo_max*solnvec[O2]/(solnvec[O2] + K_o);
        double q_s = F_s*F_e;

        // calculate final rates
        rhs[X] = y_xs*q_s*solnvec[X]*(1 - solnvec[X]/x_max);

        //set RHS of o2 to 0 as it is solved in CFD
        rhs[O2] = 0.0;
        rhs[G] = -q_s*solnvec[X];
        rhs[M] = y_ms*q_s*solnvec[X];
	rhs[CO2] = 0.0;

    }

    void advance(std::vector<double>& solnvec,int nvars,double t_now,double t_adv,double dt)
    {
        double current_time=t_now;
        double final_time=t_now+t_adv;

        std::vector<double> rhs(nvars);
        std::vector<double> solnvec_n(nvars);

        while(current_time < final_time)
        {
            current_time += dt;

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
        double F_s = Fs_max*solnvec[G]/(solnvec[G] + K_s);
        double F_e = Fo_max*solnvec[O2]/(solnvec[O2] + K_o);
        double q_s = F_s*F_e;

        double our = y_os*q_s*solnvec[X];

        return our;
    }
}
