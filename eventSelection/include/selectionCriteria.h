#ifndef H_SELECTIONCRITERIA
#define H_SELECTIONCRITERIA

#include <cstdlib>
#include <iostream>
#include <map>
#include <cassert>

#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"

#include "miscUtils.h"

enum class selectionRegion{signal=0, signal_loose, control_fakefake, control_singlemedium, control_singleloose, control_singlefake, nSelectionRegions};
int selectionRegionFirst = static_cast<int>(selectionRegion::signal);
std::map<selectionRegion, std::string> selectionRegionNames = {
  {selectionRegion::signal, "signal"},
  {selectionRegion::signal_loose, "signal_loose"},
  {selectionRegion::control_fakefake, "control_fakefake"},
  {selectionRegion::control_singlemedium, "control_singlemedium"},
  {selectionRegion::control_singleloose, "control_singleloose"},
  {selectionRegion::control_singlefake, "control_singlefake"}
};

enum class photonType{medium=0, vetoed, fake, nPhotonTypes};
int photonTypeFirst = static_cast<int>(photonType::medium);
std::map<photonType, std::string> photonTypeNames = {
  {photonType::medium, "medium"},
  {photonType::vetoed, "vetoed"},
  {photonType::fake, "fake"}
};

enum class mediumPhotonCriterion{eta=0, pT, conversionSafeElectronVeto, hOverE, neutralIsolation, photonIsolation, chargedIsolation, sigmaIEtaIEta, nMediumPhotonCriteria};
int mediumPhotonCriterionFirst = static_cast<int>(mediumPhotonCriterion::eta);
std::map<mediumPhotonCriterion, std::string> mediumPhotonCriterionNames = {
  {mediumPhotonCriterion::eta, "eta"},
  {mediumPhotonCriterion::pT, "pT"},
  {mediumPhotonCriterion::conversionSafeElectronVeto, "conversionSafeElectronVeto"},
  {mediumPhotonCriterion::hOverE, "hOverE"},
  {mediumPhotonCriterion::neutralIsolation, "neutralIsolation"},
  {mediumPhotonCriterion::photonIsolation, "photonIsolation"},
  {mediumPhotonCriterion::chargedIsolation, "chargedIsolation"},
  {mediumPhotonCriterion::sigmaIEtaIEta, "sigmaIEtaIEta"}
};

enum class vetoedPhotonCriterion{eta=0, pT, conversionSafeElectronVeto, failsMediumID, hOverE, neutralIsolation, photonIsolation, chargedIsolation, sigmaIEtaIEta, nVetoedPhotonCriteria};
int vetoedPhotonCriterionFirst = static_cast<int>(vetoedPhotonCriterion::eta);
std::map<vetoedPhotonCriterion, std::string> vetoedPhotonCriterionNames = {
  {vetoedPhotonCriterion::eta, "eta"},
  {vetoedPhotonCriterion::pT, "pT"},
  {vetoedPhotonCriterion::conversionSafeElectronVeto, "conversionSafeElectronVeto"},
  {vetoedPhotonCriterion::failsMediumID, "failsMediumID"},
  {vetoedPhotonCriterion::hOverE, "hOverE"},
  {vetoedPhotonCriterion::neutralIsolation, "neutralIsolation"},
  {vetoedPhotonCriterion::photonIsolation, "photonIsolation"},
  {vetoedPhotonCriterion::chargedIsolation, "chargedIsolation"},
  {vetoedPhotonCriterion::sigmaIEtaIEta, "sigmaIEtaIEta"}
};

enum class fakePhotonCriterion{eta=0, pT, conversionSafeElectronVeto, failsMediumID, passesShowerShapeMedIDCuts, passesInvertedChIso, nFakePhotonCriteria};
int fakePhotonCriterionFirst = static_cast<int>(fakePhotonCriterion::eta);
std::map<fakePhotonCriterion, std::string> fakePhotonCriterionNames = {
  {fakePhotonCriterion::eta, "eta"},
  {fakePhotonCriterion::pT, "pT"},
  {fakePhotonCriterion::conversionSafeElectronVeto, "conversionSafeElectronVeto"},
  {fakePhotonCriterion::failsMediumID, "failsMediumID"},
  {fakePhotonCriterion::passesShowerShapeMedIDCuts, "passesShowerShapeMedIDCuts"},
  {fakePhotonCriterion::passesInvertedChIso, "passesInvertedChIso"}
};

enum class jetCriterion{eta=0, pT, puID, jetID, deltaR_photon, nJetCriteria};
int jetCriterionFirst = static_cast<int>(jetCriterion::eta);
std::map<jetCriterion, std::string> jetCriterionNames = {
  {jetCriterion::eta, "eta"},
  {jetCriterion::pT, "pT"},
  {jetCriterion::puID, "puID"},
  {jetCriterion::jetID, "jetID"},
  {jetCriterion::deltaR_photon, "deltaR_photon"}
};

enum class eventSelectionCriterion{HLTSelection=0, MCGenInformation, overlap, doublePhoton, invariantMass, NJets, ST, /* hTCut,  */nEventSelectionCriteria};
int eventSelectionCriterionFirst = static_cast<int>(eventSelectionCriterion::HLTSelection);
std::map<eventSelectionCriterion, std::string> eventSelectionCriterionNames = {
  {eventSelectionCriterion::HLTSelection, "HLTSelection"},
  {eventSelectionCriterion::MCGenInformation, "MCGenInformation"},
  {eventSelectionCriterion::overlap, "overlap"},
  {eventSelectionCriterion::doublePhoton, "doublePhoton"},
  {eventSelectionCriterion::invariantMass, "invariantMass"},
  {eventSelectionCriterion::NJets, "NJets"},
  {eventSelectionCriterion::ST, "ST"}/*,*/
  /* {eventSelectionCriterion::hT, "hT"} */
};

struct cutflowCountersStruct{
  Long64_t N_analyzed;
  std::map<eventSelectionCriterion, Long64_t> N_passing_cut;
  std::map<eventSelectionCriterion, Long64_t> N_passing_all_cuts_upto;
  std::map<eventSelectionCriterion, Long64_t> N_passing_all_cuts_besides;
  Long64_t N_selected;

  cutflowCountersStruct() {
    N_analyzed = 0;
    for (int criterionIndex = eventSelectionCriterionFirst; criterionIndex != static_cast<int>(eventSelectionCriterion::nEventSelectionCriteria); ++criterionIndex) {
      eventSelectionCriterion criterion = static_cast<eventSelectionCriterion>(criterionIndex);
      N_passing_cut[criterion] = 0;
      N_passing_all_cuts_upto[criterion] = 0;
      N_passing_all_cuts_besides[criterion] = 0;
    }
    N_selected = 0;
  }

  /* instantiate from already filled histogram */
  cutflowCountersStruct(TH1D * counts_TH1) {
    /*
      syntax:
      bin 0: underflow (nothing stored)
      bin 1: N_analyzed
      bin 2: N_selected
      M := number of selection criteria
      bin (2+1) to (2+M): N_passing_cut
      bin (2+M+1) to (2+M+M): N_passing_all_cuts_upto
      bin (2+2M+1) to (2+2M+M): N_passing_all_cuts_besides
      bin (2+3M+1): overflow (nothing stored)
    */
    assert(static_cast<int>(counts_TH1->GetXaxis()->GetNbins()) == (2+(3*(static_cast<int>(eventSelectionCriterion::nEventSelectionCriteria)))));
    N_analyzed = static_cast<Long64_t>(0.5 + counts_TH1->GetBinContent(1));
    N_selected = static_cast<Long64_t>(0.5 + counts_TH1->GetBinContent(2));
    for (int criterionIndex = eventSelectionCriterionFirst; criterionIndex != static_cast<int>(eventSelectionCriterion::nEventSelectionCriteria); ++criterionIndex) {
      eventSelectionCriterion criterion = static_cast<eventSelectionCriterion>(criterionIndex);
      N_passing_cut[criterion] = static_cast<Long64_t>(0.5 + counts_TH1->GetBinContent(3+criterionIndex));
      N_passing_all_cuts_upto[criterion] = static_cast<Long64_t>(0.5 + counts_TH1->GetBinContent(3+(static_cast<int>(eventSelectionCriterion::nEventSelectionCriteria))+criterionIndex));
      N_passing_all_cuts_besides[criterion] = static_cast<Long64_t>(0.5 + counts_TH1->GetBinContent(3+(2*(static_cast<int>(eventSelectionCriterion::nEventSelectionCriteria)))+criterionIndex));
    }
  }

  TH1D * getCountsTH1(const std::string & name) {
    TH1D* counts_TH1 = new TH1D(name.c_str(), name.c_str(), 2+(3*(static_cast<int>(eventSelectionCriterion::nEventSelectionCriteria))), 0.5, 1.0*(0.5+(2+(3*(static_cast<int>(eventSelectionCriterion::nEventSelectionCriteria))))));
    counts_TH1->SetBinContent(1, 1.0*N_analyzed); counts_TH1->SetBinError(1, 0.); counts_TH1->GetXaxis()->SetBinLabel(1, "analyzed");
    counts_TH1->SetBinContent(2, 1.0*N_selected); counts_TH1->SetBinError(2, 0.); counts_TH1->GetXaxis()->SetBinLabel(2, "selected");
    for (int criterionIndex = eventSelectionCriterionFirst; criterionIndex != static_cast<int>(eventSelectionCriterion::nEventSelectionCriteria); ++criterionIndex) {
      eventSelectionCriterion criterion = static_cast<eventSelectionCriterion>(criterionIndex);
      Int_t offset = 0;
      offset = 3;
      counts_TH1->SetBinContent(offset+criterionIndex, 1.0*N_passing_cut.at(criterion));              counts_TH1->SetBinError(offset+criterionIndex, 0.);
      offset = 3+(static_cast<int>(eventSelectionCriterion::nEventSelectionCriteria));
      counts_TH1->SetBinContent(offset+criterionIndex, 1.0*N_passing_all_cuts_upto.at(criterion));    counts_TH1->SetBinError(offset+criterionIndex, 0.);
      offset = 3+(2*(static_cast<int>(eventSelectionCriterion::nEventSelectionCriteria)));
      counts_TH1->SetBinContent(offset+criterionIndex, 1.0*N_passing_all_cuts_besides.at(criterion)); counts_TH1->SetBinError(offset+criterionIndex, 0.);
    }
    return counts_TH1;
  }

  void fill(std::map<eventSelectionCriterion, bool> &bits) {
    N_analyzed += 1;

    bool matches_up_to_cut = true;
    for (int criterionIndex = eventSelectionCriterionFirst; criterionIndex != static_cast<int>(eventSelectionCriterion::nEventSelectionCriteria); ++criterionIndex) {
      eventSelectionCriterion criterion = static_cast<eventSelectionCriterion>(criterionIndex);
      if (bits.at(criterion)) {
        N_passing_cut[criterion] += 1;
      }
      else {
        matches_up_to_cut = false;
      }
      if (matches_up_to_cut) {
        N_passing_all_cuts_upto[criterion] += 1;
      }
    }

    int n_false = miscUtils::getNFalseBits(bits);
    if (n_false == 0) {
      N_selected += 1;
      for (int criterionIndex = eventSelectionCriterionFirst; criterionIndex != static_cast<int>(eventSelectionCriterion::nEventSelectionCriteria); ++criterionIndex) {
        eventSelectionCriterion criterion = static_cast<eventSelectionCriterion>(criterionIndex);
        N_passing_all_cuts_besides[criterion] += 1;
      }
    }
    else if (n_false == 1) {
      eventSelectionCriterion criterion = miscUtils::getFirstFalseCriterion(bits);
      N_passing_all_cuts_besides[criterion] += 1;
    }
  }

  friend std::ostream& operator<< (std::ostream& out, const cutflowCountersStruct& counters) {
    out << "N_analyzed: " << counters.N_analyzed << std::endl;
    for (int criterionIndex = eventSelectionCriterionFirst; criterionIndex != static_cast<int>(eventSelectionCriterion::nEventSelectionCriteria); ++criterionIndex) {
      eventSelectionCriterion criterion = static_cast<eventSelectionCriterion>(criterionIndex);
      out << "For criterion: " << eventSelectionCriterionNames.at(criterion)
          << ", N_passing_cut: " << counters.N_passing_cut.at(criterion)
          << ", N_passing_all_cuts_upto: " << counters.N_passing_all_cuts_upto.at(criterion)
          << ", N_passing_all_cuts_besides: " << counters.N_passing_all_cuts_besides.at(criterion) << std::endl;
    }
    out << "N_selected: " << counters.N_selected;
    return out;
  }
};

void do_sanity_checks_selectionCriteria() {
  assert(static_cast<int>(selectionRegionNames.size()) == static_cast<int>(selectionRegion::nSelectionRegions));
  assert(static_cast<int>(mediumPhotonCriterionNames.size()) == static_cast<int>(mediumPhotonCriterion::nMediumPhotonCriteria));
  assert(static_cast<int>(vetoedPhotonCriterionNames.size()) == static_cast<int>(vetoedPhotonCriterion::nVetoedPhotonCriteria));
  assert(static_cast<int>(fakePhotonCriterionNames.size()) == static_cast<int>(fakePhotonCriterion::nFakePhotonCriteria));
  assert(static_cast<int>(jetCriterionNames.size()) == static_cast<int>(jetCriterion::nJetCriteria));
  assert(static_cast<int>(eventSelectionCriterionNames.size()) == static_cast<int>(eventSelectionCriterion::nEventSelectionCriteria));
}

#endif
