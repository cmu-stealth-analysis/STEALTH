#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

#include "tmArgumentParser.h"
#include "tmProgressBar.h"
#include "tmMiscellaneous.h"
#include "tmROOTSaverUtils.h"
#include "TROOT.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TFile.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TObjArray.h"
#include "TLegend.h"
#include "TLegendEntry.h"

#include "STRegionsStruct.h"

#define MCPID_PHOTON 22
#define MCPID_GLUINO 1000021
#define MCPID_NEUTRALINO 1000022

struct parameterSpaceRegion {
  double minGluinoMass;
  double maxGluinoMass;
  double minNeutralinoMass;
  double maxNeutralinoMass;
  parameterSpaceRegion() {
    minGluinoMass = -1.0;
    maxGluinoMass = -1.0;
    minNeutralinoMass = -1.0;
    maxNeutralinoMass = -1.0;
  }
  void setParameters(double minGluinoMass_, double maxGluinoMass_, double minNeutralinoMass_, double maxNeutralinoMass_) {
    minGluinoMass = minGluinoMass_;
    maxGluinoMass = maxGluinoMass_;
    minNeutralinoMass = minNeutralinoMass_;
    maxNeutralinoMass = maxNeutralinoMass_;
  }
  bool contains(double gluinoMass, double neutralinoMass) {
    return (((gluinoMass > minGluinoMass) && (gluinoMass < maxGluinoMass)) && ((neutralinoMass > minNeutralinoMass) && (neutralinoMass < maxNeutralinoMass)));
  }
};

struct argumentsStruct {
  std::string inputMCPath, inputMCPath_JECUp, inputMCPath_JECDown, crossSectionsFilePath, MCTemplate, outputDirectory, outputPrefix;
  int n_sTBinsToPlot, nGeneratedEventsPerBin, nGluinoMassBins, nNeutralinoMassBins;
  std::string inputFile_STRegionBoundaries;
  long maxMCEvents;
  double sTMax_toPlot, integratedLuminosity, minGluinoMass, maxGluinoMass, minNeutralinoMass, maxNeutralinoMass;
  std::map<int, parameterSpaceRegion> specialZonesFor_sTDistributions;
};

struct outputHistogramsStruct {
  // syntax: histograms[JEC][regionIndex][nJetsBin] where JEC belongs to allowedJECs and regionIndex ranges from 1 to (1 + number of ST signal bins), where regionIndex 1 corresponds to the normalization bin
  std::map< std::string, std::map< int, std::map< int, TH2F* > > > h_totalNEvents;
  std::map< std::string, std::map< int, std::map< int, TH2F* > > > h_weightedNEvents;
  // syntax: histograms[specialZoneIndex][JEC][nJetsBin] where JEC belongs to allowedJECs
  std::map<int, std::map<std::string, std::map<int, TH1F* > > > h_sTDistributions;
};
