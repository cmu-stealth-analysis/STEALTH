#ifndef H_OBJECTPROPERTIES
#define H_OBJECTPROPERTIES

#include <map>
#include <vector>
#include <iostream>
#include <sstream>
#include <cassert>

#include "selectionCriteria.h"
#include "constants.h"

struct propertyAttributes{
  std::string name;
  float defaultValue;
  int plot_nBins;
  float plot_minRange, plot_maxRange;

  propertyAttributes () : name(std::string("")),
    plot_nBins(100),
    plot_minRange(0.),
    plot_maxRange(-1.) {}
  
  propertyAttributes (std::string name_, float defaultValue_, int plot_nBins_, float plot_minRange_, float plot_maxRange_) : name(name_),
    defaultValue(defaultValue_),
    plot_nBins(plot_nBins_),
    plot_minRange(plot_minRange_),
    plot_maxRange(plot_maxRange_) {}
};

enum class eventProperty{invariantMass=0, hT, MC_nPhotonsWithNeutralinoMom, nMediumPhotons, nFakePhotons, nSelectedPhotonsPassingSubLeadingpTCut, nSelectedPhotonsPassingLeadingpTCut, nJetsDR, ST, nEventProperties};
int eventPropertyFirst = static_cast<int>(eventProperty::invariantMass);
std::map<eventProperty, propertyAttributes> eventPropertyAttributes = {
  {eventProperty::invariantMass, propertyAttributes(std::string("invariantMass"), 0., 100, 50., 1050.)},
  {eventProperty::hT, propertyAttributes(std::string("hT"), 0., 100, -5., 4995.)},
  {eventProperty::MC_nPhotonsWithNeutralinoMom, propertyAttributes(std::string("MC_nPhotonsWithNeutralinoMom"), -1., 6, -1.5, 4.5)},
  {eventProperty::nMediumPhotons, propertyAttributes(std::string("nMediumPhotons"), -1., 6, -1.5, 4.5)},
  {eventProperty::nFakePhotons, propertyAttributes(std::string("nFakePhotons"), -1., 6, -1.5, 4.5)},
  {eventProperty::nSelectedPhotonsPassingSubLeadingpTCut, propertyAttributes(std::string("nSelectedPhotonsPassingSubLeadingpTCut"), -1., 6, -1.5, 4.5)},
  {eventProperty::nSelectedPhotonsPassingLeadingpTCut, propertyAttributes(std::string("nSelectedPhotonsPassingLeadingpTCut"), -1., 6, -1.5, 4.5)},
  {eventProperty::nJetsDR, propertyAttributes(std::string("nJetsDR"), -1., 22, -1.5, 20.5)},
  {eventProperty::ST, propertyAttributes(std::string("ST"), 0., 100, -5., 4995.)}
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

enum class truthPhotonProperty{eta=0, phi, pT, nTruthPhotonProperties};
int truthPhotonPropertyFirst = static_cast<int>(truthPhotonProperty::eta);
std::map<truthPhotonProperty, propertyAttributes> truthPhotonPropertyAttributes = {
  {truthPhotonProperty::eta, propertyAttributes(std::string("eta"), 0., 100, -3., 3.)},
  {truthPhotonProperty::phi, propertyAttributes(std::string("phi"), 0., 100, 0., constants::TWOPI)},
  {truthPhotonProperty::pT, propertyAttributes(std::string("pT"), -1., 100, -2., 398.)}
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

enum class photonProperty{eta=0, phi, pT, hOverE, rhoCorrectedNeutralIsolation, rhoCorrectedPhotonIsolation, rhoCorrectedChargedIsolation, sigmaIEtaIEta, deltaR_nearestTruePhoton, nPhotonProperties};
int photonPropertyFirst = static_cast<int>(photonProperty::eta);
std::map<photonProperty, propertyAttributes> photonPropertyAttributes = {
  {photonProperty::eta, propertyAttributes(std::string("eta"), 0., 100, -3., 3.)},
  {photonProperty::phi, propertyAttributes(std::string("phi"), 0., 100, 0., constants::TWOPI)},
  {photonProperty::pT, propertyAttributes(std::string("pT"), 0., 100, -2., 398.)},
  {photonProperty::hOverE, propertyAttributes(std::string("hOverE"), 0., 100, 0.01, 0.09)},
  {photonProperty::rhoCorrectedNeutralIsolation, propertyAttributes(std::string("rhoCorrectedNeutralIsolation"), 0., 100, 0.1, 50.1)},
  {photonProperty::rhoCorrectedPhotonIsolation, propertyAttributes(std::string("rhoCorrectedPhotonIsolation"), 0., 100, 1., 10.1)},
  {photonProperty::rhoCorrectedChargedIsolation, propertyAttributes(std::string("rhoCorrectedChargedIsolation"), 0., 100, 0.1, 15.1)},
  {photonProperty::sigmaIEtaIEta, propertyAttributes(std::string("sigmaIEtaIEta"), 0., 100, 0.008, 0.032)},
  {photonProperty::deltaR_nearestTruePhoton, propertyAttributes(std::string("deltaR_nearestTruePhoton"), 0., 100, 0., 1.41*constants::TWOPI)}
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

enum class jetProperty{eta=0, phi, pT, PUID, jetID, deltaR_nearestCaloPhoton, deltaR_nearestTruePhoton, nJetProperties};
int jetPropertyFirst = static_cast<int>(jetProperty::eta);
std::map<jetProperty, propertyAttributes> jetPropertyAttributes = {
  {jetProperty::eta, propertyAttributes(std::string("eta"), 0., 100, -3., 3.)},
  {jetProperty::phi, propertyAttributes(std::string("phi"), 0., 100, 0., constants::TWOPI)},
  {jetProperty::pT, propertyAttributes(std::string("pT"), 0., 100, 0., 1000.)},
  {jetProperty::PUID, propertyAttributes(std::string("PUID"), 0., 100, 0., 1.)},
  {jetProperty::jetID, propertyAttributes(std::string("jetID"), 0., 100, 0., -1.)},
  {jetProperty::deltaR_nearestCaloPhoton, propertyAttributes(std::string("deltaR_nearestCaloPhoton"), 0., 100, 0., 1.41*constants::TWOPI)},
  {jetProperty::deltaR_nearestTruePhoton, propertyAttributes(std::string("deltaR_nearestTruePhoton"), 0., 100, 0., 1.41*constants::TWOPI)}
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

void do_sanity_checks_objectProperties() {
  assert(static_cast<int>(eventPropertyAttributes.size()) == static_cast<int>(eventProperty::nEventProperties));
  assert(static_cast<int>(truthPhotonPropertyAttributes.size()) == static_cast<int>(truthPhotonProperty::nTruthPhotonProperties));
  assert(static_cast<int>(photonPropertyAttributes.size()) == static_cast<int>(photonProperty::nPhotonProperties));
  assert(static_cast<int>(jetPropertyAttributes.size()) == static_cast<int>(jetProperty::nJetProperties));
}

#endif
