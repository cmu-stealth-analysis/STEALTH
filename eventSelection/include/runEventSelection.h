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
/* #include "TTreeReader.h" */
/* #include "TTreeReaderValue.h" */
/* #include "TTreeReaderArray.h" */
/* #include "TFile.h" */
/* #include "TH1F.h" */
/* #include "TH2F.h" */
/* #include "TObjArray.h" */
/* #include "TLegend.h" */
/* #include "TLegendEntry.h" */

#define MCPID_PHOTON 22
#define MCPID_GLUINO 1000021
#define MCPID_NEUTRALINO 1000022

struct PIDsStruct {
  const int photon = 22;
  const int gluino = 1000021;
  const int neutralino = 1000022;
};

struct rangeStruct{
  float rangeLower, rangeUpper;
rangeStruct(float rangeLower_, float rangeUpper_) : rangeLower(rangeLower_),
    rangeUpper(rangeUpper_){}
  bool isInRange(float candidate) {
    return (candidate >= rangeLower && candidate < rangeUpper);
  }
};

struct neutralIsolationCutStruct{
  float constCoefficient, linearCoefficient, squareCoefficient;
neutralIsolationCutStruct(float constCoefficient_, float linearCoefficient_, float squareCoefficient_) : constCoefficient(constCoefficient_),
    linearCoefficient(linearCoefficient_),
    squareCoefficient(squareCoefficient_) {}
  float polynomialValue(float pT) {
    return (constCoefficient + linearCoefficient*pT + squareCoefficient*pT*pT);
  }
};

struct photonIsolationCutStruct{
  float constCoefficient, linearCoefficient;
photonIsolationCutStruct(float constCoefficient_, float linearCoefficient_) : constCoefficient(constCoefficient_),
    linearCoefficient(linearCoefficient_) {}
  float polynomialValue(float pT) {
    return (constCoefficient + linearCoefficient*pT);
  }
};

struct EAValuesStruct{
  float regionUpperBound, chargedHadronsEA, neutralHadronsEA, photonsEA;
EAValuesStruct(float regionUpperBound_, float chargedHadronsEA_, float neutralHadronsEA_, float photonsEA_) : regionUpperBound(regionUpperBound_),
    chargedHadronsEA(chargedHadronsEA_),
    neutralHadronsEA(neutralHadronsEA_),
    photonsEA(photonsEA_) {}
};

struct parametersStruct {
  const float pTCutSubLeading = 25.0;
  const float pTCutLeading = 35.0;
  const float photonEtaCut = 1.442;
  const float R9Cut = 1.0;
  const int nSubLeadingPhotons = 2;
  const int nLeadingPhotons = 1;
  const float jetEtaCut = 2.4;
  const float jetpTCut = 30.;
  const float jetPUIDThreshold = 0.61;
  const int jetSelectionID = 6;
  const float minDeltaRCut = 0.4;
  const int nJetsCut = 2;
  const float HTCut = 60.;
  const float electronPtCut = 15.;
  const float electronEtaCut = 2.5;
  const float electronDzCut = 0.1;
  const float electronPFPUIsoCut = 0.1;
  const int nElectronsCut = 0;
  const float muonPtCut = 15.0;
  const float muonPFPUIsoCut = 0.12;
  const int nMuonsCut = 0;
  const float METThreshold = 15.0;
  const EAValuesStruct region1EAs = EAValuesStruct(1.0, 0.0385, 0.0636, 0.124);
  const EAValuesStruct region2EAs = EAValuesStruct(1.479, 0.0468, 0.1103, 0.1093);
  const PIDsStruct PIDs;
  const int nPhotonsWithNeutralinoMom = 2;
  const float towerHOverECut = 0.035;
  const rangeStruct sigmeietaietaRange = rangeStruct(0.0103, 0.015);
  const rangeStruct chargedIsolationRange = rangeStruct(1.416, 15.0);
  const neutralIsolationCutStruct neutralIsolationCut = neutralIsolationCutStruct(2.491, 0.0126, 0.000026);
  const photonIsolationCutStruct photonIsolationCut = photonIsolationCutStruct(2.952, 0.004);
  int HLTPhotonBit = -1;
  const int HLTPhotonBit2016 = -1;
  const int HLTPhotonBit2017 = -1;
  void setHLTPhotonBitForYear(int year) {
    if (year == 2017) HLTPhotonBit = HLTPhotonBit2017;
    else if (year == 2016) HLTPhotonBit = HLTPhotonBit2016;
  }
};

struct optionsStruct {
  std::string inputFilePath, outputFilePath, photonSelectionType;
  long counterStartInclusive, counterEndInclusive;
  int year, JECUncertainty;
};
