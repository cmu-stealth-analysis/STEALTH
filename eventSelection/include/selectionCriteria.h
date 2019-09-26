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

enum class selectionRegion{signal=0, control_fakefake, control_mediumfake, nSelectionRegions};
std::map<selectionRegion, std::string> selectionRegionNames = {
  {selectionRegion::signal, "signal"},
  {selectionRegion::control_fakefake, "control_fakefake"},
  {selectionRegion::control_mediumfake, "control_mediumfake"}
};

enum class mediumPhotonCriterion{eta=0, pT, hOverE, neutralIsolation, photonIsolation, conversionSafeElectronVeto, chargedIsolation, sigmaIEtaIEta, nMediumPhotonCriteria};
int mediumPhotonCriterionFirst = static_cast<int>(mediumPhotonCriterion::eta);
std::map<mediumPhotonCriterion, std::string> mediumPhotonCriterionNames = {
  {mediumPhotonCriterion::eta, "eta"},
  {mediumPhotonCriterion::pT, "pT"},
  {mediumPhotonCriterion::hOverE, "hOverE"},
  {mediumPhotonCriterion::neutralIsolation, "neutralIsolation"},
  {mediumPhotonCriterion::photonIsolation, "photonIsolation"},
  {mediumPhotonCriterion::conversionSafeElectronVeto, "conversionSafeElectronVeto"},
  {mediumPhotonCriterion::chargedIsolation, "chargedIsolation"},
  {mediumPhotonCriterion::sigmaIEtaIEta, "sigmaIEtaIEta"}
};

/* enum class fakePhotonCriterion{eta=0, pT, hOverE, neutralIsolation, photonIsolation, conversionSafeElectronVeto, notTightChIsoAndSigmaIEtaIEta, chargedIsolationLoose, sigmaIEtaIEtaLoose, nFakePhotonCriteria}; */
/* int fakePhotonCriterionFirst = static_cast<int>(fakePhotonCriterion::eta); */
/* std::map<fakePhotonCriterion, std::string> fakePhotonCriterionNames = { */
/*   {fakePhotonCriterion::eta, "eta"}, */
/*   {fakePhotonCriterion::pT, "pT"}, */
/*   {fakePhotonCriterion::hOverE, "hOverE"}, */
/*   {fakePhotonCriterion::neutralIsolation, "neutralIsolation"}, */
/*   {fakePhotonCriterion::photonIsolation, "photonIsolation"}, */
/*   {fakePhotonCriterion::conversionSafeElectronVeto, "conversionSafeElectronVeto"}, */
/*   {fakePhotonCriterion::notTightChIsoAndSigmaIEtaIEta, "notTightChIsoAndSigmaIEtaIEta"}, */
/*   {fakePhotonCriterion::chargedIsolationLoose, "chargedIsolationLoose"}, */
/*   {fakePhotonCriterion::sigmaIEtaIEtaLoose, "sigmaIEtaIEtaLoose"} */
/* }; */

enum class fakePhotonCriterion{eta=0, pT, failsMediumID, passesChIsoVeto, passesPhoIsoVeto, passesOtherLooseCuts, nFakePhotonCriteria};
int fakePhotonCriterionFirst = static_cast<int>(fakePhotonCriterion::eta);
std::map<fakePhotonCriterion, std::string> fakePhotonCriterionNames = {
  {fakePhotonCriterion::eta, "eta"},
  {fakePhotonCriterion::pT, "pT"},
  {fakePhotonCriterion::failsMediumID, "failsMediumID"},
  {fakePhotonCriterion::passesChIsoVeto, "passesChIsoVeto"},
  {fakePhotonCriterion::passesPhoIsoVeto, "passesPhoIsoVeto"},
  {fakePhotonCriterion::passesOtherLooseCuts, "passesOtherLooseCuts"}
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

enum class eventSelectionCriterion{HLTPhoton=0, MCGenInformation, photonEnergy, photonQuality, invariantMass, NJets, /* hTCut,  */nEventSelectionCriteria};
int eventSelectionCriterionFirst = static_cast<int>(eventSelectionCriterion::HLTPhoton);
std::map<eventSelectionCriterion, std::string> eventSelectionCriterionNames = {
  {eventSelectionCriterion::HLTPhoton, "HLTPhoton"},
  {eventSelectionCriterion::MCGenInformation, "MCGenInformation"},
  {eventSelectionCriterion::photonEnergy, "photonEnergy"},
  {eventSelectionCriterion::photonQuality, "photonQuality"},
  {eventSelectionCriterion::invariantMass, "invariantMass"},
  {eventSelectionCriterion::NJets, "NJets"}/* , */
  /* {eventSelectionCriterion::hT, "hT"} */
};

void do_sanity_checks_selectionCriteria() {
  assert(static_cast<int>(selectionRegionNames.size()) == static_cast<int>(selectionRegion::nSelectionRegions));
  assert(static_cast<int>(mediumPhotonCriterionNames.size()) == static_cast<int>(mediumPhotonCriterion::nMediumPhotonCriteria));
  assert(static_cast<int>(fakePhotonCriterionNames.size()) == static_cast<int>(fakePhotonCriterion::nFakePhotonCriteria));
  assert(static_cast<int>(jetCriterionNames.size()) == static_cast<int>(jetCriterion::nJetCriteria));
  assert(static_cast<int>(eventSelectionCriterionNames.size()) == static_cast<int>(eventSelectionCriterion::nEventSelectionCriteria));
}

#endif
