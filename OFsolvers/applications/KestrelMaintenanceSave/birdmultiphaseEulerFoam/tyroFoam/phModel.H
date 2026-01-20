#include<vector>
#include<scalar.H>
#include<array>
#include<map>
//#include <boost/math/special_functions/gamma.hpp>
//#include <gamma.hpp>

namespace acidbasemodel
{
    extern const int PB;
    extern const int PA;
    extern const int NB;
    extern const int NA;
    extern const int AB;
    extern const int AA;
    extern const int BB;
    extern const int BA;
    // extern const int MA;
    extern const int PhosphateBuffer;
    extern const int AmmoniaBase;
    extern const int AceticAcid;
    extern const int ButyricAcid;
    extern const int nvars;
    extern std::map<std::string, int> system_keys;

    extern const double Kw;
    extern std::vector<double> Ka;
    extern double pH_avg;

    void getTotalConc(std::vector<double>& abSystem, std::vector<double>& totalConc);
    void getSpectatorsNonBio(std::vector<double>& totalConc);
    void getSystemID(std::string name, std::vector<int>& id_map, int foam_id);
    double charge(double H, std::vector<double>& totalConc);
    double dchargedH(double H, std::vector<double>& totalConc);
    double NewtonRaphson(double H_init, std::vector<double>& totalConc);  
}
