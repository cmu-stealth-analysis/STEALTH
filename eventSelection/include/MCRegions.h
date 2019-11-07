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
    {5, "eventProgenitorNeutralinoDegenerate"}
  };

  int getRegionIndex(const float& generated_eventProgenitorMass, const float& generated_neutralinoMass) {
    if (generated_eventProgenitorMass > 1200.) {
      if (generated_neutralinoMass < 150.) return 4;
      else if (generated_neutralinoMass > 250.) {
        if ((generated_eventProgenitorMass - generated_neutralinoMass) < 150.) return 5;
        else if ((generated_eventProgenitorMass - generated_neutralinoMass) > 250.) {
	  if (generated_eventProgenitorMass < 2000) return 1;
	  else if (generated_eventProgenitorMass < 2200) return 2;
	  else return 3;
	}
      }
    }
    return 0;
  }
}

#endif
