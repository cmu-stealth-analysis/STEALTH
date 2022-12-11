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
    {5, "eventProgenitorNeutralinoDegenerate"},
    {6, "eventProgenitor2100_neutralino1000"},
    {7, "eventProgenitor1950_neutralino1850"},
    {8, "eventProgenitor1850_neutralino200"}
  };

  int getRegionIndex(const float& generated_eventProgenitorMass, const float& generated_neutralinoMass) {
    if (((generated_eventProgenitorMass > 2075.) && (generated_eventProgenitorMass < 2125.)) &&
        ((generated_neutralinoMass > 993.75) && (generated_neutralinoMass < 1006.25))) {
      return 6;
    }
    if (((generated_eventProgenitorMass > 1925.) && (generated_eventProgenitorMass < 1975.)) &&
        ((generated_neutralinoMass > 1843.75) && (generated_neutralinoMass < 1856.25))) {
      return 7;
    }
    if (((generated_eventProgenitorMass > 1825.) && (generated_eventProgenitorMass < 1875.)) &&
        ((generated_neutralinoMass > 193.75) && (generated_neutralinoMass < 206.25))) {
      return 8;
    }
    if (generated_eventProgenitorMass > 1200.) {
      if (generated_neutralinoMass < 150.) {
        // generated_eventProgenitorMass > 1200. &&
        // generated_neutralinoMass < 150.
        return 4;
      }
      else if (generated_neutralinoMass > 250.) {
        if ((generated_eventProgenitorMass - generated_neutralinoMass) < 150.) {
          // generated_eventProgenitorMass > 1200. &&
          // generated_neutralinoMass > 250. &&
          // (generated_eventProgenitorMass - generated_neutralinoMass) < 150.
          return 5;
        }
        else if ((generated_eventProgenitorMass - generated_neutralinoMass) > 250.) {
	  if (generated_eventProgenitorMass < 2000) {
            // generated_eventProgenitorMass > 1200. &&
            // generated_neutralinoMass > 250. &&
            // (generated_eventProgenitorMass - generated_neutralinoMass) > 250. &&
            // (generated_eventProgenitorMass < 2000)
            return 1;
          }
	  else if (generated_eventProgenitorMass < 2200) {
            // generated_eventProgenitorMass > 1200. &&
            // generated_neutralinoMass > 250. &&
            // (generated_eventProgenitorMass - generated_neutralinoMass) > 250. &&
            // (generated_eventProgenitorMass >= 2000) &&
            // (generated_eventProgenitorMass < 2200)
            return 2;
          }
	  else {
            // generated_eventProgenitorMass > 1200. &&
            // generated_neutralinoMass > 250. &&
            // (generated_eventProgenitorMass - generated_neutralinoMass) > 250. &&
            // (generated_eventProgenitorMass >= 2000) &&
            // (generated_eventProgenitorMass >= 2200)
            return 3;
          }
	}
      }
    }
    return 0;
  }
}

#endif
