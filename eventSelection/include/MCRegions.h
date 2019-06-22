#ifndef H_MCREGIONS
#define H_MCREGIONS

#include <cstdlib>
#include <map>
#include <iostream>

namespace MCRegions {
  std::map<int, std::string> regionNames = {
    {1, "bulk"},
    {2, "lowNeutralinoMass"},
    {3, "gluinoNeutralinoDegenerate"}
  };

  int getRegionIndex(const float& generated_gluinoMass, const float& generated_neutralinoMass) {
    if (generated_gluinoMass > 1200.) {
      if (generated_neutralinoMass < 150.) return 2;
      else if (generated_neutralinoMass > 250.) {
        if ((generated_gluinoMass - generated_neutralinoMass) < 150.) return 3;
        else if ((generated_gluinoMass - generated_neutralinoMass) > 250.) return 1;
      }
    }
    return 0;
  }
}

#endif
