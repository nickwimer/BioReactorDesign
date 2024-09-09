#include<vector>
#include<scalar.H>
#include<array>
#include<map>
//#include <boost/math/special_functions/gamma.hpp>
//#include <gamma.hpp>

namespace microbemodel
{
    extern const int X;
    extern const int O2;
    extern const int G;
    extern const int M;
    extern const int CO2;
    extern const int nvars;
    extern std::map<std::string, int> sp_list;
    extern const double y_xs;
    extern const double y_ms;
    extern const double y_os;
    extern const double x_max;
    extern const double Fo_max;
    extern const double Fs_max;
    extern const double o2_max;
    extern const double K_o;
    extern const double K_s;
    extern const double kLa;

    extern double X_avg;
    extern double G_avg;
    extern double M_avg;
    extern double O2_avg;

    void get_rhs(std::vector<double>& rhs,std::vector<double> solnvec,double t,int nvars);
    void advance(std::vector<double>& solnvec,int nvars,double t_now,double t_adv,double dt);
    double get_our(std::vector<double> solnvec,int nvars);
    void get_sp_id(std::string name, std::vector<int>& id_map, int foam_id);

    class sp_info;
}
