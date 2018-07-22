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

std::string inputPath_, MCTemplate_, outputDirectory_, outputPrefix_;
int nGluinoMassBins_, nNeutralinoMassBins_;
double minGluinoMass_, maxGluinoMass_, minNeutralinoMass_, maxNeutralinoMass_;
double sTMin_normWindow_, sTMax_normWindow_, sTStartMainRegion_, minGluinoMassToPlot_;

std::vector<std::string> allowedJECs{"JECDown", "JECNominal", "JECUp"};
std::vector<std::string> allowedZones{"sub", "main"}; // "norm" not needed

// syntax: inputHistograms[JEC][zone][nJetsBin] where JEC belongs to allowedJECs and zone belongs to allowedZones
std::map< std::string, std::map< std::string, std::map< int, TH2F* > > > h_totalNEvents;
std::map< std::string, std::map< std::string, std::map< int, TH2F* > > > h_weightedNEvents;

// syntax: outputHistogram[zone][nJetsBin]
std::map< std::string, std::map< int, TGraph2D* > > g_JECUncertainty;
std::map< std::string, std::map< int, TH2F* > >     h_JECUncertainty;
std::map< std::string, std::map< int, TGraph2D* > > g_JECUncertainty_fractionalError;
std::map< std::string, std::map< int, TH2F* > >     h_JECUncertainty_fractionalError;
