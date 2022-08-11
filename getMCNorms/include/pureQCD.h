#include <cstdlib>
#include <algorithm>
#include <initializer_list>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cassert>

#include "tmProgressBar.h"
#include "tmMiscellaneous.h"
#include "tmROOTSaverUtils.h"

#include "Rtypes.h"
#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TAxis.h"
#include "TH1.h"
#include "TH1D.h"

#include "common.h"

#define PROG_MASS_MIN 1975.0
#define PROG_MASS_MAX 2025.0
#define NEUT_MASS_MIN 993.75
#define NEUT_MASS_MAX 1006.25
#define LUMI_2016_PB 36303.1
#define LUMI_2017_PB 41527.3
#define LUMI_2018_PB 59736.0
#define XSEC_WGT 0.00101
#define N_RECORDED_EVTS 41496

struct eventDataStruct {
  float evtST;
  float pT_leadingJet;
  double MCXSecWeight;
  float MCGenWeight;
  float prefiringWeight;
  float photonMCScaleFactor;
  double MCPUWeight;
  std::vector<UShort_t> * phoID = nullptr;
  int photonIndex_leading;
  int photonIndex_subLeading;
  int nJetsDR;
  int nJetsAll;
  float event_progenitor_mass;
  float event_neutralino_mass;

  eventDataStruct() : evtST(-1.), pT_leadingJet(-1.), MCXSecWeight(-1.), MCGenWeight(-1.), prefiringWeight(-1.), photonMCScaleFactor(-1.), MCPUWeight(-1.), phoID(nullptr), photonIndex_leading(-1), photonIndex_subLeading(-1), nJetsDR(-1), nJetsAll(-1), event_progenitor_mass(-1), event_neutralino_mass(-1) {}

  eventDataStruct(const float & evtST_, const float & pT_leadingJet_, const double & MCXSecWeight_, const float & MCGenWeight_, const float & prefiringWeight_, const float & photonMCScaleFactor_, const float & MCPUWeight_, std::vector<UShort_t> * phoID_, const int & photonIndex_leading_, const int & photonIndex_subLeading_, const int & nJetsDR_, const int & nJetsAll_, const float & event_progenitor_mass_, const float & event_neutralino_mass_) : evtST(evtST_), pT_leadingJet(pT_leadingJet_), MCXSecWeight(MCXSecWeight_), MCGenWeight(MCGenWeight_), prefiringWeight(prefiringWeight_), photonMCScaleFactor(photonMCScaleFactor_), MCPUWeight(MCPUWeight_), phoID(phoID_), photonIndex_leading(photonIndex_leading_), photonIndex_subLeading(photonIndex_subLeading_), nJetsDR(nJetsDR_), nJetsAll(nJetsAll_), event_progenitor_mass(event_progenitor_mass_), event_neutralino_mass(event_neutralino_mass_) {}

};
