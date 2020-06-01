#ifndef H_SELECTIONCRITERIA
#define H_SELECTIONCRITERIA

#include <cstdlib>
#include <iostream>
#include <map>
#include <cassert>

#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"

enum class selectionRegion{signal=0, signal_loose, control_fakefake, control_singlemedium, nSelectionRegions};
int selectionRegionFirst = static_cast<int>(selectionRegion::signal);
std::map<selectionRegion, std::string> selectionRegionNames = {
  {selectionRegion::signal, "signal"},
  {selectionRegion::signal_loose, "signal_loose"},
  {selectionRegion::control_fakefake, "control_fakefake"},
  {selectionRegion::control_singlemedium, "control_singlemedium"}
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

enum class vetoedPhotonCriterion{eta=0, pT, conversionSafeElectronVeto, failsMediumID, passesShowerShapeLooseIDCuts, chIsoBetweenMedAndLoose, passesNeutIsoAndPhoIsoLooseCriteria, nVetoedPhotonCriteria};
int vetoedPhotonCriterionFirst = static_cast<int>(vetoedPhotonCriterion::eta);
std::map<vetoedPhotonCriterion, std::string> vetoedPhotonCriterionNames = {
  {vetoedPhotonCriterion::eta, "eta"},
  {vetoedPhotonCriterion::pT, "pT"},
  {vetoedPhotonCriterion::conversionSafeElectronVeto, "conversionSafeElectronVeto"},
  {vetoedPhotonCriterion::failsMediumID, "failsMediumID"},
  {vetoedPhotonCriterion::passesShowerShapeLooseIDCuts, "passesShowerShapeLooseIDCuts"},
  {vetoedPhotonCriterion::chIsoBetweenMedAndLoose, "chIsoBetweenMedAndLoose"},
  {vetoedPhotonCriterion::passesNeutIsoAndPhoIsoLooseCriteria, "passesNeutIsoAndPhoIsoLooseCriteria"}
};

enum class fakePhotonCriterion{eta=0, pT, conversionSafeElectronVeto, failsMediumID, passesShowerShapeMedIDCuts, chIsoBetweenLooseAndExtraLoose, nFakePhotonCriteria};
int fakePhotonCriterionFirst = static_cast<int>(fakePhotonCriterion::eta);
std::map<fakePhotonCriterion, std::string> fakePhotonCriterionNames = {
  {fakePhotonCriterion::eta, "eta"},
  {fakePhotonCriterion::pT, "pT"},
  {fakePhotonCriterion::conversionSafeElectronVeto, "conversionSafeElectronVeto"},
  {fakePhotonCriterion::failsMediumID, "failsMediumID"},
  {fakePhotonCriterion::passesShowerShapeMedIDCuts, "passesShowerShapeMedIDCuts"},
  {fakePhotonCriterion::chIsoBetweenLooseAndExtraLoose, "chIsoBetweenLooseAndExtraLoose"}
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

enum class eventSelectionCriterion{HLTSelection=0, MCGenInformation, doublePhoton, invariantMass, NJets, /* hTCut,  */nEventSelectionCriteria};
int eventSelectionCriterionFirst = static_cast<int>(eventSelectionCriterion::HLTSelection);
std::map<eventSelectionCriterion, std::string> eventSelectionCriterionNames = {
  {eventSelectionCriterion::HLTSelection, "HLTSelection"},
  {eventSelectionCriterion::MCGenInformation, "MCGenInformation"},
  {eventSelectionCriterion::doublePhoton, "doublePhoton"},
  {eventSelectionCriterion::invariantMass, "invariantMass"},
  {eventSelectionCriterion::NJets, "NJets"}/* , */
  /* {eventSelectionCriterion::hT, "hT"} */
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
