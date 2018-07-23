#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
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
#include "TH1F.h"
#include "TH2F.h"
#include "TObjArray.h"
#include "TLegend.h"
#include "TLegendEntry.h"

#define MCPID_PHOTON 22
#define MCPID_GLUINO 1000021
#define MCPID_NEUTRALINO 1000022

std::string inputMCPath, inputMCPath_JECUp, inputMCPath_JECDown, crossSectionsFilePath, MCTemplate, outputDirectory, outputPrefix;
int n_sTBinsToPlot, nGeneratedEventsPerBin, nGluinoMassBins, nNeutralinoMassBins;
long maxMCEvents;
double sTMin_normWindow, sTMax_normWindow, sTStartMainRegion, sTMax_toPlot, integratedLuminosity, minGluinoMass, maxGluinoMass, minNeutralinoMass, maxNeutralinoMass;

std::map<int, double> crossSections;
std::map<int, double> crossSectionsFractionalUncertainty;

std::vector<std::string> allowedJECs{"JECDown", "JECNominal", "JECUp"};
std::vector<std::string> allowedZones{"norm", "sub", "main"};

// syntax: histograms[JEC][zone][nJetsBin] where JEC belongs to allowedJECs and zone belongs to allowedZones
std::map< std::string, std::map< std::string, std::map< int, TH2F* > > > h_totalNEvents;
std::map< std::string, std::map< std::string, std::map< int, TH2F* > > > h_weightedNEvents;

// syntax: histograms[JEC][nJetsBin] where JEC belongs to allowedJECs
std::map< std::string, std::map< int, TH1F* > > h_sTDistributions;
