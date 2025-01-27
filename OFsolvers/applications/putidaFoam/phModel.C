#include "acidBaseModel.H"
#include<map>

namespace acidbasemodel
{

  const int PB = 0;
  const int PA = 1;
  const int NB = 2;
  const int NA = 3;
  const int MB = 4;
  const int MA = 5;
  
  const int PhosphateBuffer = 0;
  const int AmmoniaBase = 1;
  const int MuconicAcid = 2;
  const int nvars = 6;

  const std::vector<double> Ka(nvars/2);
  Ka[PhosphateBuffer] = 6.2E-8;
  Ka[AmmoniaBase] = 5.56E-10;
  Ka[MuconicAcid] = 0.0;
  const double Kw = 1.0E-14;
  
  std::map<std::string, int> system_keys = {{"HPO4.liquid", 0},
					    {"H2PO4.liquid", 1},
					    {"NH3.liquid", 2},
					    {"NH4.liquid", 3},
					    {"C6H4O4.liquid", 4},
					    {"C6H6O4.liquid", 5}};
  
  void getSystemID(std::string name, std::vector<int>& id_map, int foam_id)
  {
    auto it = system_keys.find(name);
    if(it != system_keys.end())
      {
	id_map[it->second] = foam_id;
      }
  }

  void acidRemainderUpdate(std::vector<double>& abSystem, int aRef, int bRef)
  {
    int aID = aRef*2 - 1;
    int conjBaseID = aID - 1;
    int bID = (bRef - 1)*2;
    int conjAcidID = bID + 1;
    
    double baseTotalUse = abSystem[bID];

    abSystem[aID] -= 0.5*baseTotalUse;
    abSystem[conjBaseID] += 0.5baseTotalUse;
    abSystem[bID] -= baseTotalUse;
    abSystem[conjAcidID] += baseTotalUse;
  }

  void baseRemainderUpdate(std::vector<double>& abSystem, int aRef, int bRef)
  {
    int aID = aRef*2 - 1;
    int conjBaseID = aID - 1;
    int bID = (bRef - 1)*2;
    int conjAcidID = bID + 1;
    
    double acidTotalUse = abSystem[aID];

    abSystem[aID] -= acidTotalUse;
    abSystem[conjBaseID] += acidTotalUse;
    abSystem[bID] -= 2.0*acidTotalUse;
    abSystem[conjAcidID] += 2.0*acidTotalUse;
  }
  
  double phBufferSystem(double wAcidConc, double wBaseConc, int abRef)
  {    
    double pH = -log(Ka[abRef-1]) + log(wAcidConc/wBaseConc);
    
    return pH;
  }

  double phWeakAcid(double wAcidConc, int abRef)
  { 
    double hConc = 0.5*(-Ka[abRef-1] + sqrt(pow(Ka[abRef-1],2) - 4.0*Ka[abRef-1]*wAcidConc));
    double pH = -log(hConc);
    
    return pH;
  }

  double phWeakBase(double wBaseConc, int abRef)
  { 
    double ohConc = 0.5*(-Kw/Ka[abRef-1] + sqrt(pow(Kw/Ka[abRef-1],2) - 4.0*Kw/Ka[abRef-1]*wBaseConc));
    double pH = -log(Kw/ohConc);
    
    return pH;
  }
  
  double phStrongAcid(double sAcidConc)
  {
    double pH = -log(sAcidConc);

    return pH;
  }
}
