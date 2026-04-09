#include<vector>
#include<scalar.H>
#include<array>
#include<map>
//#include <boost/math/special_functions/gamma.hpp>
//#include <gamma.hpp>

namespace microbemodel
{
    extern const int X;
    extern const int G;
    extern const int B;
    extern const int A;
    extern const int H;
    //extern const int CO2;
    extern const int nvars;
    extern std::map<std::string, int> sp_keys;
    // extern const double y_xs;
    // extern const double y_ms;
    // extern const double y_os;
    // extern const double x_max;
    // extern const double q_max;
    // extern const double o2_max;
    // extern const double K_o;
    // extern const double K_s;
    // extern const double kLa;
    extern std::map<std::string, int> MM_param_keys;

    extern double X_avg;
    extern double G_avg;
    extern double B_avg;
    extern double A_avg;

    void get_rhs(std::vector<double>& rhs,std::vector<double> solnvec,double t,int nvars, std::vector<double>& MM_params);
    void get_rhs_ml(std::vector<double>& rhs,std::vector<double> solnvec,double t,int nvars, std::vector<double>& MM_params);
    void advance(std::vector<double>& solnvec,int nvars,double t_now,double t_adv,double dt, std::vector<double>& MM_params, bool do_ml);
    double get_ratio(std::vector<double> solnvec,int nvars, std::vector<double>& MM_params);
    void get_sp_id(std::string name, std::vector<int>& id_map, int foam_id);

    class sp_info;
}
