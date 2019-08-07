#ifndef H_OBJECTPROPERTIES
#define H_OBJECTPROPERTIES

#include <map>
#include <vector>
#include <iostream>
#include <sstream>
#include <cassert>

#include "selectionCriteria.h"
#include "constants.h"

#define NBINS_DEFAULT 500

struct propertyAttributes{
  std::string name;
  float defaultValue;
  int plot_nBins;
  float plot_minRange, plot_maxRange;

  propertyAttributes () : name(std::string("")),
    defaultValue(-0.5),
    plot_nBins(NBINS_DEFAULT),
    plot_minRange(0.),
    plot_maxRange(-1.) {}

  void runChecks(std::string name_, int nBins, float minRange, float maxRange) {
    assert(name_ != "");
    assert(nBins != 0);
    assert(minRange < maxRange);
  }

  propertyAttributes (std::string name_, int nBins, float minRange, float maxRange) : name(name_) {
    runChecks(name_, nBins, minRange, maxRange);
    float binWidth = std::fabs((maxRange - minRange)/nBins);
    plot_nBins = 1+nBins;
    defaultValue = minRange - 0.5*binWidth;
    plot_minRange = minRange - binWidth;
    plot_maxRange = maxRange;
  }

  propertyAttributes (std::string name_, float minRange, float maxRange) : name(name_) {
    runChecks(name_, NBINS_DEFAULT, minRange, maxRange);
    float binWidth = (maxRange - minRange)/NBINS_DEFAULT;
    plot_nBins = 1+NBINS_DEFAULT;
    defaultValue = minRange - 0.5*binWidth;
    plot_minRange = minRange - binWidth;
    plot_maxRange = maxRange;
  }
};

enum class eventProperty{invariantMass=0, hT, MC_nPhotonsWithNeutralinoMom, MC_nGenJets, MC_nNonPhotonGenJets, MC_nJetCandidatesWithStealthMom, MC_nJetCandidatesWithGluinoMom, MC_nJetCandidatesWithSingletMom, MC_nStealthJetsCloseToTruePhoton, nMediumPhotons, MC_nTruthMatchedMediumPhotons, nFakePhotons, MC_nTruthMatchedFakePhotons, nSelectedPhotonsPassingSubLeadingpTCut, nSelectedPhotonsPassingLeadingpTCut, nGoodJetsCloseToSelectedPhoton, nJetsDR, ST, selectionRegionIndex, nEventProperties};
int eventPropertyFirst = static_cast<int>(eventProperty::invariantMass);
std::map<eventProperty, propertyAttributes> eventPropertyAttributes = {
  {eventProperty::invariantMass, propertyAttributes(std::string("invariantMass"), 0., 1000.)},
  {eventProperty::hT, propertyAttributes(std::string("hT"), 0., 5000.)},
  {eventProperty::MC_nGenJets, propertyAttributes(std::string("MC_nGenJets"), 21, -0.5, 20.5)},
  {eventProperty::MC_nNonPhotonGenJets, propertyAttributes(std::string("MC_nNonPhotonGenJets"), 21, -0.5, 20.5)},
  {eventProperty::MC_nPhotonsWithNeutralinoMom, propertyAttributes(std::string("MC_nPhotonsWithNeutralinoMom"), 5, -0.5, 4.5)},
  {eventProperty::MC_nJetCandidatesWithStealthMom, propertyAttributes(std::string("MC_nJetCandidatesWithStealthMom"), 21, -0.5, 20.5)},
  {eventProperty::MC_nJetCandidatesWithGluinoMom, propertyAttributes(std::string("MC_nJetCandidatesWithGluinoMom"), 21, -0.5, 20.5)},
  {eventProperty::MC_nJetCandidatesWithSingletMom, propertyAttributes(std::string("MC_nJetCandidatesWithSingletMom"), 21, -0.5, 20.5)},
  {eventProperty::MC_nStealthJetsCloseToTruePhoton, propertyAttributes(std::string("MC_nStealthJetsCloseToTruePhoton"), 4, -0.5, 3.5)},
  {eventProperty::nMediumPhotons, propertyAttributes(std::string("nMediumPhotons"), 5, -0.5, 4.5)},
  {eventProperty::MC_nTruthMatchedMediumPhotons, propertyAttributes(std::string("MC_nTruthMatchedMediumPhotons"), 5, -0.5, 4.5)},
  {eventProperty::nFakePhotons, propertyAttributes(std::string("nFakePhotons"), 5, -0.5, 4.5)},
  {eventProperty::MC_nTruthMatchedFakePhotons, propertyAttributes(std::string("MC_nTruthMatchedFakePhotons"), 5, -0.5, 4.5)},
  {eventProperty::nSelectedPhotonsPassingSubLeadingpTCut, propertyAttributes(std::string("nSelectedPhotonsPassingSubLeadingpTCut"), 5, -0.5, 4.5)},
  {eventProperty::nSelectedPhotonsPassingLeadingpTCut, propertyAttributes(std::string("nSelectedPhotonsPassingLeadingpTCut"), 5, -0.5, 4.5)},
  {eventProperty::nGoodJetsCloseToSelectedPhoton, propertyAttributes(std::string("nGoodJetsCloseToSelectedPhoton"), 5, -0.5, 4.5)},
  {eventProperty::nJetsDR, propertyAttributes(std::string("nJetsDR"), 21, -0.5, 20.5)},
  {eventProperty::ST, propertyAttributes(std::string("ST"), 0., 5000.)},
  {eventProperty::selectionRegionIndex, propertyAttributes(std::string("selectionRegionIndex"), 1+static_cast<int>(selectionRegion::nSelectionRegions), -0.5, 0.5+static_cast<int>(selectionRegion::nSelectionRegions))}
};
typedef std::map<eventProperty, float> eventProperties;
eventProperties initialize_eventProperties_with_defaults() {
  eventProperties properties;
  for (int index = eventPropertyFirst; index != static_cast<int>(eventProperty::nEventProperties); ++index) {
    eventProperty propertiesIndex = static_cast<eventProperty>(index);
    properties[propertiesIndex] = (eventPropertyAttributes[propertiesIndex]).defaultValue;
  }
  return properties;
}
typedef std::vector<eventProperties> eventPropertiesCollection;
typedef std::pair<eventSelectionCriterion, eventProperties> unselectedEventProperties;
typedef std::vector<unselectedEventProperties> unselectedEventPropertiesCollection;

enum class truthPhotonProperty{eta=0, phi, pT, status, deltaR_nearestTruthJetCandidate, deltaR_nearestGenJet, deltaR_nearestNonPhotonGenJet, nTruthPhotonProperties};
int truthPhotonPropertyFirst = static_cast<int>(truthPhotonProperty::eta);
std::map<truthPhotonProperty, propertyAttributes> truthPhotonPropertyAttributes = {
  {truthPhotonProperty::eta, propertyAttributes(std::string("eta"), -5., 5.)},
  {truthPhotonProperty::phi, propertyAttributes(std::string("phi"), (-1.0*constants::PI), constants::PI)},
  {truthPhotonProperty::pT, propertyAttributes(std::string("pT"), 0., 1000.)},
  {truthPhotonProperty::status, propertyAttributes(std::string("status"), 201, -0.5, 200.5)},
  {truthPhotonProperty::deltaR_nearestTruthJetCandidate, propertyAttributes(std::string("deltaR_nearestTruthJetCandidate"), 0., constants::TWOPI)},
  {truthPhotonProperty::deltaR_nearestGenJet, propertyAttributes(std::string("deltaR_nearestGenJet"), 0., constants::TWOPI)},
  {truthPhotonProperty::deltaR_nearestNonPhotonGenJet, propertyAttributes(std::string("deltaR_nearestNonPhotonGenJet"), 0., constants::TWOPI)}
};
typedef std::map<truthPhotonProperty, float> truthPhotonProperties;
truthPhotonProperties initialize_truthPhotonProperties_with_defaults() {
  truthPhotonProperties properties;
  for (int index = truthPhotonPropertyFirst; index != static_cast<int>(truthPhotonProperty::nTruthPhotonProperties); ++index) {
    truthPhotonProperty propertiesIndex = static_cast<truthPhotonProperty>(index);
    properties[propertiesIndex] = (truthPhotonPropertyAttributes[propertiesIndex]).defaultValue;
  }
  return properties;
}
typedef std::vector<truthPhotonProperties> truthPhotonPropertiesCollection;

enum class truthJetCandidateProperty{eta=0, phi, pT, momID, status, statusFlag, deltaR_nearestTruePhoton, nTruthJetCandidateProperties};
int truthJetCandidatePropertyFirst = static_cast<int>(truthJetCandidateProperty::eta);
std::map<truthJetCandidateProperty, propertyAttributes> truthJetCandidatePropertyAttributes = {
  {truthJetCandidateProperty::eta, propertyAttributes(std::string("eta"), -5., 5.)},
  {truthJetCandidateProperty::phi, propertyAttributes(std::string("phi"), (-1.0*constants::PI), constants::PI)},
  {truthJetCandidateProperty::pT, propertyAttributes(std::string("pT"), 0., 1000.)},
  {truthJetCandidateProperty::momID, propertyAttributes(std::string("momID"), 10, 0.5, 10.5)},
  {truthJetCandidateProperty::status, propertyAttributes(std::string("status"), 201, -0.5, 200.5)},
  {truthJetCandidateProperty::statusFlag, propertyAttributes(std::string("statusFlag"), 9, -0.5, 8.5)},
  {truthJetCandidateProperty::deltaR_nearestTruePhoton, propertyAttributes(std::string("deltaR_nearestTruePhoton"), 0., 0.5)}
};
typedef std::map<truthJetCandidateProperty, float> truthJetCandidateProperties;
truthJetCandidateProperties initialize_truthJetCandidateProperties_with_defaults() {
  truthJetCandidateProperties properties;
  for (int index = truthJetCandidatePropertyFirst; index != static_cast<int>(truthJetCandidateProperty::nTruthJetCandidateProperties); ++index) {
    truthJetCandidateProperty propertiesIndex = static_cast<truthJetCandidateProperty>(index);
    properties[propertiesIndex] = (truthJetCandidatePropertyAttributes[propertiesIndex]).defaultValue;
  }
  return properties;
}
typedef std::vector<truthJetCandidateProperties> truthJetCandidatePropertiesCollection;

enum class photonProperty{eta=0, phi, pT, hOverE, energy, rhoCorrectedNeutralIsolation, rhoCorrectedPhotonIsolation, rawChargedIsolation, rhoCorrectedChargedIsolation, sigmaIEtaIEta, deltaR_nearestTruePhoton, R9, ecalClusIso, trkIso, deltaR_nearestGenJet, deltaR_nearestNonPhotonGenJet, nPhotonProperties};
int photonPropertyFirst = static_cast<int>(photonProperty::eta);
std::map<photonProperty, propertyAttributes> photonPropertyAttributes = {
  {photonProperty::eta, propertyAttributes(std::string("eta"), -5., 5.)},
  {photonProperty::phi, propertyAttributes(std::string("phi"), (-1.0*constants::PI), constants::PI)},
  {photonProperty::pT, propertyAttributes(std::string("pT"), 0., 1000.)},
  {photonProperty::energy, propertyAttributes(std::string("energy"), 0., 1000.)},
  {photonProperty::hOverE, propertyAttributes(std::string("hOverE"), 0., 0.1)},
  {photonProperty::rhoCorrectedNeutralIsolation, propertyAttributes(std::string("rhoCorrectedNeutralIsolation"), 0., 50.)},
  {photonProperty::rhoCorrectedPhotonIsolation, propertyAttributes(std::string("rhoCorrectedPhotonIsolation"), 0., 10.)},
  {photonProperty::rawChargedIsolation, propertyAttributes(std::string("rawChargedIsolation"), 0., 15.)},
  {photonProperty::rhoCorrectedChargedIsolation, propertyAttributes(std::string("rhoCorrectedChargedIsolation"), 0., 15.)},
  {photonProperty::sigmaIEtaIEta, propertyAttributes(std::string("sigmaIEtaIEta"), 0.008, 0.032)},
  {photonProperty::deltaR_nearestTruePhoton, propertyAttributes(std::string("deltaR_nearestTruePhoton"), 0., 0.5)},
  {photonProperty::R9, propertyAttributes(std::string("R9"), 0., 1.0)},
  {photonProperty::ecalClusIso, propertyAttributes(std::string("ecalClusIso"), 0., 15.)},
  {photonProperty::trkIso, propertyAttributes(std::string("trkIso"), 0., 15.)},
  {photonProperty::deltaR_nearestGenJet, propertyAttributes(std::string("deltaR_nearestGenJet"), 0., constants::TWOPI)},
  {photonProperty::deltaR_nearestNonPhotonGenJet, propertyAttributes(std::string("deltaR_nearestNonPhotonGenJet"), 0., constants::TWOPI)}
};
typedef std::map<photonProperty, float> photonProperties;
photonProperties initialize_photonProperties_with_defaults() {
  photonProperties properties;
  for (int index = photonPropertyFirst; index != static_cast<int>(photonProperty::nPhotonProperties); ++index) {
    photonProperty propertiesIndex = static_cast<photonProperty>(index);
    properties[propertiesIndex] = (photonPropertyAttributes[propertiesIndex]).defaultValue;
  }
  return properties;
}
typedef std::vector<photonProperties> photonPropertiesCollection;
typedef std::pair<mediumPhotonCriterion, photonProperties> unselectedMediumPhotonProperties;
typedef std::vector<unselectedMediumPhotonProperties> unselectedMediumPhotonPropertiesCollection;
typedef std::pair<fakePhotonCriterion, photonProperties> unselectedFakePhotonProperties;
typedef std::vector<unselectedFakePhotonProperties> unselectedFakePhotonPropertiesCollection;

enum class jetProperty{eta=0, phi, pT, PUID, jetID, deltaR_nearestCaloPhoton, deltaR_nearestTruePhoton, deltaR_nearestTrueJetCandidate, truthPTRatio, deltaR_genJet, nJetProperties};
int jetPropertyFirst = static_cast<int>(jetProperty::eta);
std::map<jetProperty, propertyAttributes> jetPropertyAttributes = {
  {jetProperty::eta, propertyAttributes(std::string("eta"), -5., 5.)},
  {jetProperty::phi, propertyAttributes(std::string("phi"), (-1.0*constants::PI), constants::PI)},
  {jetProperty::pT, propertyAttributes(std::string("pT"), 0., 1000.)},
  {jetProperty::PUID, propertyAttributes(std::string("PUID"), 0., 1.)},
  {jetProperty::jetID, propertyAttributes(std::string("jetID"), 0., 200.)},
  {jetProperty::deltaR_nearestCaloPhoton, propertyAttributes(std::string("deltaR_nearestCaloPhoton"), 0., constants::TWOPI)},
  {jetProperty::deltaR_nearestTruePhoton, propertyAttributes(std::string("deltaR_nearestTruePhoton"), 0., constants::TWOPI)},
  {jetProperty::deltaR_nearestTrueJetCandidate, propertyAttributes(std::string("deltaR_nearestTrueJetCandidate"), 0., 0.5)},
  {jetProperty::truthPTRatio, propertyAttributes(std::string("truthPTRatio"), 0., 2.)},
  {jetProperty::deltaR_genJet, propertyAttributes(std::string("deltaR_genJet"), 0., 0.5)}
};
typedef std::map<jetProperty, float> jetProperties;
jetProperties initialize_jetProperties_with_defaults() {
  jetProperties properties;
  for (int index = jetPropertyFirst; index != static_cast<int>(jetProperty::nJetProperties); ++index) {
    jetProperty propertiesIndex = static_cast<jetProperty>(index);
    properties[propertiesIndex] = (jetPropertyAttributes[propertiesIndex]).defaultValue;
  }
  return properties;
}
typedef std::vector<jetProperties> jetPropertiesCollection;
typedef std::pair<jetCriterion, jetProperties> unselectedJetProperties;
typedef std::vector<unselectedJetProperties> unselectedJetPropertiesCollection;

enum class genJetProperty{eta=0, phi, pT, partonID, partonMomID, nGenJetProperties};
int genJetPropertyFirst = static_cast<int>(genJetProperty::eta);
std::map<genJetProperty, propertyAttributes> genJetPropertyAttributes = {
  {genJetProperty::eta, propertyAttributes(std::string("eta"), -5., 5.)},
  {genJetProperty::phi, propertyAttributes(std::string("phi"), (-1.0*constants::PI), constants::PI)},
  {genJetProperty::pT, propertyAttributes(std::string("pT"), 0., 1000.)},
  {genJetProperty::partonID, propertyAttributes(std::string("partonID"), 16, -0.5, 15.5)},
  {genJetProperty::partonMomID, propertyAttributes(std::string("partonMomID"), 16, -0.5, 15.5)}
};
typedef std::map<genJetProperty, float> genJetProperties;
genJetProperties initialize_genJetProperties_with_defaults() {
  genJetProperties properties;
  for (int index = genJetPropertyFirst; index != static_cast<int>(genJetProperty::nGenJetProperties); ++index) {
    genJetProperty propertiesIndex = static_cast<genJetProperty>(index);
    properties[propertiesIndex] = (genJetPropertyAttributes[propertiesIndex]).defaultValue;
  }
  return properties;
}
typedef std::vector<genJetProperties> genJetPropertiesCollection;

void do_sanity_checks_objectProperties() {
  assert(static_cast<int>(eventPropertyAttributes.size()) == static_cast<int>(eventProperty::nEventProperties));
  assert(static_cast<int>(truthPhotonPropertyAttributes.size()) == static_cast<int>(truthPhotonProperty::nTruthPhotonProperties));
  assert(static_cast<int>(truthJetCandidatePropertyAttributes.size()) == static_cast<int>(truthJetCandidateProperty::nTruthJetCandidateProperties));
  assert(static_cast<int>(photonPropertyAttributes.size()) == static_cast<int>(photonProperty::nPhotonProperties));
  assert(static_cast<int>(jetPropertyAttributes.size()) == static_cast<int>(jetProperty::nJetProperties));
  assert(static_cast<int>(genJetPropertyAttributes.size()) == static_cast<int>(genJetProperty::nGenJetProperties));
}

#endif
