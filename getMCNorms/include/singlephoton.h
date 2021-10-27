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

struct eventDataStruct {
  float evtST;
  float pT_leadingPhoton;
  double MCXSecWeight;
  float MCGenWeight;
  float prefiringWeight;
  float photonMCScaleFactor;
  double MCPUWeight;
  int nJetsDR;
  int nJetsAll;

eventDataStruct() : evtST(-1.), pT_leadingPhoton(-1.), MCXSecWeight(-1.), MCGenWeight(-1.), prefiringWeight(-1.), photonMCScaleFactor(-1.), MCPUWeight(-1.), nJetsDR(-1), nJetsAll(-1) {}

eventDataStruct(const float & evtST_, const float & pT_leadingPhoton_, const double & MCXSecWeight_, const float & MCGenWeight_, const float & prefiringWeight_, const float & photonMCScaleFactor_, const float & MCPUWeight_, const int & nJetsDR_, const int & nJetsAll_) : evtST(evtST_), pT_leadingPhoton(pT_leadingPhoton_), MCXSecWeight(MCXSecWeight_), MCGenWeight(MCGenWeight_), prefiringWeight(prefiringWeight_), photonMCScaleFactor(photonMCScaleFactor_), MCPUWeight(MCPUWeight_), nJetsDR(nJetsDR_), nJetsAll(nJetsAll_) {}

};
