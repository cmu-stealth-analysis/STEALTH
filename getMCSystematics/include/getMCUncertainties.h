#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

#include "tmArgumentParser.h"
#include "tmROOTSaverUtils.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH2F.h"
#include "TGraph2D.h"

struct optionsStruct{
  std::string inputPath, MCTemplate, outputDirectory, outputPrefix;
  int nGluinoMassBins, nNeutralinoMassBins;
  double minGluinoMass, maxGluinoMass, minNeutralinoMass, maxNeutralinoMass;
  double sTMin_normWindow, sTMax_normWindow, sTStartMainRegion, minGluinoMassToPlot;
};

struct outputHistogramsStruct{
  // syntax: outputHistogram[zone][nJetsBin]
  std::map< std::string, std::map< int, TH2F* > > h_JECUncertainty;
  std::map< std::string, std::map< int, TH2F* > > h_MCStatisticsFractionalError;
};

struct inputHistogramsStruct{
  // syntax: inputHistograms[JEC][zone][nJetsBin] where JEC belongs to allowedJECs and zone belongs to allowedZones
  std::map< std::string, std::map< std::string, std::map< int, TH2F* > > > h_totalNEvents;
  std::map< std::string, std::map< std::string, std::map< int, TH2F* > > > h_weightedNEvents;
};
