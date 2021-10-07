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
  double MCXSecWeight;
  float prefiringWeight;
  float photonMCScaleFactor;
  std::vector<UShort_t> * phoID = nullptr;
  int photonIndex_leading;
  int photonIndex_subLeading;
  int nJetsDR;

  eventDataStruct() : evtST(-1.), MCXSecWeight(-1.), prefiringWeight(-1.), photonMCScaleFactor(-1.), phoID(nullptr), photonIndex_leading(-1), photonIndex_subLeading(-1), nJetsDR(-1) {}

  eventDataStruct(const float & evtST_, const double & MCXSecWeight_, const float & prefiringWeight_, const float & photonMCScaleFactor_, std::vector<UShort_t> * phoID_, const int & photonIndex_leading_, const int & photonIndex_subLeading_, const int & nJetsDR_) : evtST(evtST_), MCXSecWeight(MCXSecWeight_), prefiringWeight(prefiringWeight_), photonMCScaleFactor(photonMCScaleFactor_), phoID(phoID_), photonIndex_leading(photonIndex_leading_), photonIndex_subLeading(photonIndex_subLeading_), nJetsDR(nJetsDR_) {}

};
