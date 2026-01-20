#include "phModel.H"
#include<map>
//#include<vector>

namespace acidbasemodel
{
  const int PB = 0;
  const int PA = 1;
  const int NB = 2;
  const int NA = 3;
  const int AB = 4;
  const int AA = 5;
  const int BB = 6;
  const int BA = 7;
  
  const int PhosphateBuffer = 0;
  const int AmmoniaBase = 1;
  const int AceticAcid = 2;
  const int ButyricAcid = 3;
  const int nvars = 8;

  std::vector<double> Ka = {6.2E-8, 5.56E-10, 1.8E-05, 1.5E-05};
  //Ka[0] = 6.2E-8;
  //Ka[AmmoniaBase] = 5.56E-10;
  //Ka[MuconicAcid] = 0.0;
  const double Kw = 1.0E-14;

  double pH_avg = 0.0;
  
  std::map<std::string, int> system_keys = {{"HPO4.liquid", 0},
					    {"H2PO4.liquid", 1},
					    {"NH3.liquid", 2},
					    {"NH4.liquid", 3},
					    {"CH3COO.liquid", 4},
					    {"CH3COOH.liquid", 5},
                                            {"C3H7COO.liquid", 6},
					    {"C3H7COOH.liquid", 7}};
  
  void getSystemID(std::string name, std::vector<int>& id_map, int foam_id)
  {
    auto it = system_keys.find(name);
    if(it != system_keys.end())
      {
	id_map[it->second] = foam_id;
      }
  }

  void getTotalConc(std::vector<double>& abSystem, std::vector<double>& totalConc)
  {    
    for(int i = 0; i<nvars/2 - 1; i++)
      {
	totalConc[i] = abSystem[2*i] + abSystem[2*i+1];
      }
  }

  void getSpectatorsNonBio(std::vector<double>& totalConc)
  {
    double H_set = std::pow(10, -6.3);;
    double Z = H_set - H_set*totalConc[PhosphateBuffer]/(H_set + Ka[PhosphateBuffer]) \
      - 2.0 * totalConc[PhosphateBuffer] * Ka[PhosphateBuffer]/(H_set + Ka[PhosphateBuffer]) \
      - Kw/H_set;
      
    totalConc[nvars/2] = -Z;
  }

  double charge(double H, std::vector<double>& totalConc)
  {
    // getRatio; // figure this out later...

    double ratio = 2.0;
    double ch = H						\
      + totalConc[AmmoniaBase] * H / (H + Ka[AmmoniaBase])	\
      - totalConc[AceticAcid] * H / (H + Ka[AceticAcid]) \
      - totalConc[ButyricAcid] * H / (H + Ka[ButyricAcid]) \
      - totalConc[PhosphateBuffer] * H / (H + Ka[PhosphateBuffer]) \
      - 2.0 * totalConc[PhosphateBuffer] * Ka[PhosphateBuffer] / (H + Ka[PhosphateBuffer]) \
      - Kw / H \
      + totalConc[nvars/2];

    return ch;
  }

  double dchargedH(double H, std::vector<double>& totalConc)
  {
    double dchdH = 1.0 \
      + Ka[AmmoniaBase] * totalConc[AmmoniaBase] / ((Ka[AmmoniaBase] + H) * (Ka[AmmoniaBase] + H)) \
      - Ka[ButyricAcid] * totalConc[ButyricAcid] / ((Ka[ButyricAcid] + H) * (Ka[ButyricAcid] + H)) \
      - Ka[AceticAcid] * totalConc[AceticAcid] / ((Ka[AceticAcid] + H) * (Ka[AceticAcid] + H)) \
      - Ka[PhosphateBuffer] * totalConc[PhosphateBuffer] / ((Ka[PhosphateBuffer] + H) * (Ka[PhosphateBuffer] + H)) \
      + 2.0 * totalConc[PhosphateBuffer] * Ka[PhosphateBuffer] / ((Ka[PhosphateBuffer] + H) * (Ka[PhosphateBuffer] + H)) \
      + Kw/(H*H);

    return dchdH;
  }

  double NewtonRaphson(double H_init, std::vector<double>& totalConc)
  {
    double tol = 1.0E-14;
    double test = charge(H_init, totalConc);
    // std::cout << "Initialized Charge is: " << test << "\n";
    double H_next = 0.0;
    double H_now = H_init;
    int iter = 0;
    while(std::abs(test) > tol && iter<100)
      {
	double chargeCurrent = charge(H_now, totalConc);
	double dchargedHCurrent = dchargedH(H_now, totalConc);
	// std::cout << "current Charge: " << chargeCurrent << "\n";
	// std::cout << "current derivative: " << dchargedHCurrent << "\n";
	H_next =  H_now - charge(H_now, totalConc)/dchargedH(H_now, totalConc);
	test = charge(H_next, totalConc);
	H_now = H_next;
	// std::cout << "\n";
	// std::cout << "Updated Charge: " << test << "\n";
	// std::cout << "Updated H: " << H_now << "\n";
	// std::cout << "\n\n";
	iter += 1;
      }
    return H_now;
  }
}
