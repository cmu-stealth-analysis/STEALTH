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
  float pT_leadingJet;
  double MCXSecWeight;
  float prefiringWeight;
  float photonMCScaleFactor;
  int nJetsDR;

  eventDataStruct() : pT_leadingJet(-1.), MCXSecWeight(-1.), prefiringWeight(-1.), photonMCScaleFactor(-1.), nJetsDR(-1) {}

  eventDataStruct(const float & pT_leadingJet_, const double & MCXSecWeight_, const float & prefiringWeight_, const float & photonMCScaleFactor_, const int & nJetsDR_) : pT_leadingJet(pT_leadingJet_), MCXSecWeight(MCXSecWeight_), prefiringWeight(prefiringWeight_), photonMCScaleFactor(photonMCScaleFactor_), nJetsDR(nJetsDR_) {}

};