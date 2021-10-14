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

#define DIPH_STRINGIFY(x) #x
#define DIPH_NJETS_TH1_NAME DIPH_STRINGIFY(nJets_in_normST)
#define DIPH_STNORM_MIN 1200.0
#define DIPH_STNORM_MAX 1300.0
#define DIPH_INVMASS_TH1_NAME DIPH_STRINGIFY(invMass_zeroJets)
#define DIPH_INVMASS_MIN 110.0
#define DIPH_INVMASS_MAX 200.0
#define DIPH_INVMASS_NBINS 9
#define DIPH_TOLERANCE_SANITY_CHECK 0.0001

struct eventDataStruct {
  float evtST;
  float invMass;
  double MCXSecWeight;
  float MCGenWeight;
  float prefiringWeight;
  float photonMCScaleFactor;
  std::vector<UShort_t> * phoID = nullptr;
  int photonIndex_leading;
  int photonIndex_subLeading;
  int nJetsDR;

eventDataStruct() : evtST(-1.), invMass(-1.), MCXSecWeight(-1.), MCGenWeight(-1.), prefiringWeight(-1.), photonMCScaleFactor(-1.), phoID(nullptr), photonIndex_leading(-1), photonIndex_subLeading(-1), nJetsDR(-1) {}

eventDataStruct(const float & evtST_, const float & invMass_, const double & MCXSecWeight_, const float & MCGenWeight_, const float & prefiringWeight_, const float & photonMCScaleFactor_, std::vector<UShort_t> * phoID_, const int & photonIndex_leading_, const int & photonIndex_subLeading_, const int & nJetsDR_) : evtST(evtST_), invMass(invMass_), MCXSecWeight(MCXSecWeight_), MCGenWeight(MCGenWeight_), prefiringWeight(prefiringWeight_), photonMCScaleFactor(photonMCScaleFactor_), phoID(phoID_), photonIndex_leading(photonIndex_leading_), photonIndex_subLeading(photonIndex_subLeading_), nJetsDR(nJetsDR_) {}

};
