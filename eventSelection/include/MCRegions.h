#ifndef H_MCREGIONS
#define H_MCREGIONS

#include <cstdlib>
#include <map>
#include <iostream>

namespace MCRegions {
  std::map<int, std::string> regionNames = {
    {1, "bulk_wellExcluded"},
    {2, "bulk_closeToContours"},
    {3, "bulk_notExcluded"},
    {4, "lowNeutralinoMass"},
    {5, "gluinoNeutralinoDegenerate"}
  };

  int getRegionIndex(const float& generated_gluinoMass, const float& generated_neutralinoMass) {
    if (generated_gluinoMass > 1200.) {
      if (generated_neutralinoMass < 150.) return 4;
      else if (generated_neutralinoMass > 250.) {
        if ((generated_gluinoMass - generated_neutralinoMass) < 150.) return 5;
        else if ((generated_gluinoMass - generated_neutralinoMass) > 250.) {
	  if (generated_gluinoMass < 2000) return 1;
	  else if (generated_gluinoMass < 2200) return 2;
	  else return 3;
	}
      }
    }
    return 0;
  }
}

#endif
