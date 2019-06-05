#ifndef H_STATISTICSHISTOGRAMS
#define H_STATISTICSHISTOGRAMS

#include <iostream>
#include <cstdlib>
#include <map>

#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"

#include "objectProperties.h"
#include "selectionCriteria.h"

class statisticsHistograms {
 public:
  std::map<std::string, TH1F> stats;
  std::map<int, std::string> MCBinNames;

  void initializeWithCheck(std::string& name, int& nBins, float& xmin, float& xmax) {
    if (stats.find(name) == stats.end()) {
      stats.insert(std::make_pair(name, TH1F(name.c_str(), name.c_str(), nBins, xmin, xmax)));
      // stats[name].SetCanExtend(TH1::kAllAxes);
    }
    else {
      std::cout << "ERROR: tried to create new statistics histogram with name \"" << name << "\", but it already exists!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  statisticsHistograms(const bool& isMC, std::map<int, std::string>& MCBinNames_) {
    for (auto&& MCBinNamesElement: MCBinNames_) {
      MCBinNames.insert(std::make_pair(MCBinNamesElement.first, MCBinNamesElement.second));
    }
    for (auto&& eventPropertyAttributesElement: eventPropertyAttributes) {
      auto& attributesElement = (eventPropertyAttributesElement.second);
      std::string& eventPropertyName = attributesElement.name;
      for (auto&& selectionRegionNamesElement: selectionRegionNames) {
        std::string& selectionRegionName = selectionRegionNamesElement.second;
        std::string fullName = eventPropertyName + std::string("_") + selectionRegionName + std::string("_selectedEvents");
        initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
      }
      for (auto&& eventSelectionCriterionNamesElement: eventSelectionCriterionNames) {
        std::string& eventSelectionCriterionName = eventSelectionCriterionNamesElement.second;
        std::string fullName = eventPropertyName + std::string("_marginallyUnselectedEvents_") + eventSelectionCriterionName;
        initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
      }
    }

    if (isMC) {
      for (auto&& truthPhotonPropertyAttributesElement: truthPhotonPropertyAttributes) {
        auto& attributesElement = (truthPhotonPropertyAttributesElement.second);
        std::string& truthPhotonPropertyName = attributesElement.name;
        for (auto&& selectionRegionNamesElement: selectionRegionNames) {
          std::string& selectionRegionName = selectionRegionNamesElement.second;
          for (auto&& MCBinNamesElement: MCBinNames) {
            std::string& MCBinName = MCBinNamesElement.second;
            std::string fullName = truthPhotonPropertyName + "_" + selectionRegionName + "_truePhotons_MC_" + MCBinName;
            initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
          }
        }
      }
    }

    for (auto&& selectionRegionNamesElement: selectionRegionNames) {
      std::string& selectionRegionName = selectionRegionNamesElement.second;
      for (auto&& photonPropertyAttributesElement: photonPropertyAttributes) {
        auto& attributesElement = (photonPropertyAttributesElement.second);
        std::string& photonPropertyName = attributesElement.name;
        std::string fullName;
        fullName = photonPropertyName + "_" + selectionRegionName + "_selectedMediumCaloPhotons";
        initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
        fullName = photonPropertyName + "_" + selectionRegionName + "_selectedFakeCaloPhotons";
        initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
        if (isMC) {
          for (auto&& MCBinNamesElement: MCBinNames) {
            std::string& MCBinName = MCBinNamesElement.second;
            fullName = photonPropertyName + "_" + selectionRegionName + "_selectedMediumCaloPhotons_closeToTruePhoton_MC_" + MCBinName;
            initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
            fullName = photonPropertyName + "_" + selectionRegionName + "_selectedMediumCaloPhotons_awayFromTruePhotons_MC_" + MCBinName;
            initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
          }
        } // MC plots
        for (auto&& mediumPhotonCriterionNamesElement: mediumPhotonCriterionNames) {
          std::string& mediumPhotonCriterionName = mediumPhotonCriterionNamesElement.second;
          std::string fullName = photonPropertyName + "_" + selectionRegionName + "_marginallyUnselectedMediumCaloPhotons_" + mediumPhotonCriterionName;
          initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
          if (isMC) {
            for (auto&& MCBinNamesElement: MCBinNames) {
              std::string& MCBinName = MCBinNamesElement.second;
              fullName = photonPropertyName + "_" + selectionRegionName + "_marginallyUnselectedMediumCaloPhotons_" + mediumPhotonCriterionName + "_closeToTruePhoton_MC_" + MCBinName;
              initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
              fullName = photonPropertyName + "_" + selectionRegionName + "_marginallyUnselectedMediumCaloPhotons_" + mediumPhotonCriterionName + "_awayFromTruePhotons_MC_" + MCBinName;
              initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
            }
          } // MC plots
        } // medium photon criteria
        for (auto&& fakePhotonCriterionNamesElement: fakePhotonCriterionNames) {
          std::string& fakePhotonCriterionName = fakePhotonCriterionNamesElement.second;
          std::string fullName = photonPropertyName + "_" + selectionRegionName + "_marginallyUnselectedFakeCaloPhotons_" + fakePhotonCriterionName;
          initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
          if (isMC) {
            for (auto&& MCBinNamesElement: MCBinNames) {
              std::string& MCBinName = MCBinNamesElement.second;
              fullName = photonPropertyName + "_" + selectionRegionName + "_marginallyUnselectedFakeCaloPhotons_" + fakePhotonCriterionName + "_closeToTruePhoton_MC_" + MCBinName;
              initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
              fullName = photonPropertyName + "_" + selectionRegionName + "_marginallyUnselectedFakeCaloPhotons_" + fakePhotonCriterionName + "_awayFromTruePhotons_MC_" + MCBinName;
              initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
            }
          } // MC plots
        } // fake photon criteria
      } // photon plots
      for (auto&& jetPropertyAttributesElement: jetPropertyAttributes) {
        auto& attributesElement = (jetPropertyAttributesElement.second);
        std::string& jetPropertyName = attributesElement.name;
        std::string fullName = jetPropertyName + "_" + selectionRegionName + "_selectedCaloJets";
        initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
        if (isMC) {
          for (auto&& MCBinNamesElement: MCBinNames) {
            std::string& MCBinName = MCBinNamesElement.second;
            fullName = jetPropertyName + "_" + selectionRegionName + "_selectedCaloJets_closeToTruePhoton_MC_" + MCBinName;
            initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
            fullName = jetPropertyName + "_" + selectionRegionName + "_selectedCaloJets_awayFromTruePhotons_MC_" + MCBinName;
            initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
          }
        } // MC plots
        for (auto&& jetCriterionNamesElement: jetCriterionNames) {
          std::string& jetCriterionName = jetCriterionNamesElement.second;
          fullName = jetPropertyName + "_" + selectionRegionName + "_marginallyUnselectedCaloJets_" + jetCriterionName;
          initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
          if (isMC) {
            for (auto&& MCBinNamesElement: MCBinNames) {
              std::string& MCBinName = MCBinNamesElement.second;
              fullName = jetPropertyName + "_" + selectionRegionName + "_marginallyUnselectedCaloJets_" + jetCriterionName + "_closeToTruePhoton_MC_" + MCBinName;
              initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
              fullName = jetPropertyName + "_" + selectionRegionName + "_marginallyUnselectedCaloJets_" + jetCriterionName + "_awayFromTruePhotons_MC_" + MCBinName;
              initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
            }
          }
        }
      } // jet plots
    }
  }

  void fillStatisticsHistogramByName(const std::string& histogramName, const float& value, const float& weight) {
    if (stats.find(histogramName) == stats.end()) {
      std::cout << "ERROR: tried to fill statistics histogram with name \"" << histogramName << "\"; a histogram with this name was not initialized!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    else {
      (stats[histogramName]).Fill(value, weight);
    }
  }

  void fillStatisticsHistogramByName(const std::string& histogramName, const float& value) {
    fillStatisticsHistogramByName(histogramName, value, 1.0);
  }

  std::string getStatisticsHistogramName(const eventProperty& event_property, const selectionRegion& region) {
    return ((eventPropertyAttributes[event_property]).name + "_" + selectionRegionNames[region] + "_selectedEvents");
  }

  std::string getStatisticsHistogramName(const eventProperty& event_property, const eventSelectionCriterion& criterion) {
    return ((eventPropertyAttributes[event_property]).name + "_marginallyUnselectedEvents_" + eventSelectionCriterionNames[criterion]);
  }
  
  std::string getStatisticsHistogramName(const truthPhotonProperty& truthPhotonProperty, const selectionRegion& region, const int& MCBinIndex) {
    return ((truthPhotonPropertyAttributes[truthPhotonProperty]).name + "_" + selectionRegionNames[region] + "_truePhotons_MC_" + MCBinNames[MCBinIndex]);
  }
  
  std::string getStatisticsHistogramName(const photonProperty& property, const selectionRegion& region, const bool& trueMedium_falseFake) {
    std::string name = ((photonPropertyAttributes[property]).name + "_" + selectionRegionNames[region] + "_selected");
    if (trueMedium_falseFake) name += "MediumCaloPhotons";
    else name += "FakeCaloPhotons";
    return name;
  }

  std::string getStatisticsHistogramName(const photonProperty& property, const selectionRegion& region, const bool& trueMedium_falseFake, const bool& trueClose_falseAway, const int& MCBinIndex) {
    std::string name = ((photonPropertyAttributes[property]).name + "_" + selectionRegionNames[region] + "_selected");
    if (trueMedium_falseFake) name += "MediumCaloPhotons_";
    else name += "FakeCaloPhotons_";
    if (trueClose_falseAway) name += "closeToTruePhoton_MC_";
    else name += "awayFromTruePhotons_MC_";
    name += MCBinNames[MCBinIndex];
    return name;
  }

  std::string getStatisticsHistogramName(const photonProperty& property, const selectionRegion& region, const mediumPhotonCriterion& mediumCriterion) {
    return ((photonPropertyAttributes[property]).name + "_" + selectionRegionNames[region] + "_marginallyUnselectedMediumCaloPhotons_" + mediumPhotonCriterionNames[mediumCriterion]);
  }

  std::string getStatisticsHistogramName(const photonProperty& property, const selectionRegion& region, const mediumPhotonCriterion& mediumCriterion, const bool& trueClose_falseAway, const int& MCBinIndex) {
    std::string name = (photonPropertyAttributes[property]).name + "_" + selectionRegionNames[region] + "_marginallyUnselectedMediumCaloPhotons_" + mediumPhotonCriterionNames[mediumCriterion] + "_";
    if (trueClose_falseAway) name += "closeToTruePhoton_MC_";
    else name += "awayFromTruePhotons_MC_";
    name += MCBinNames[MCBinIndex];
    return name;
  }

  std::string getStatisticsHistogramName(const photonProperty& property, const selectionRegion& region, const fakePhotonCriterion& fakeCriterion) {
    return ((photonPropertyAttributes[property]).name + "_" + selectionRegionNames[region] + "_marginallyUnselectedFakeCaloPhotons_" + fakePhotonCriterionNames[fakeCriterion]);
  }

  std::string getStatisticsHistogramName(const photonProperty& property, const selectionRegion& region, const fakePhotonCriterion& fakeCriterion, const bool& trueClose_falseAway, const int& MCBinIndex) {
    std::string name = (photonPropertyAttributes[property]).name + "_" + selectionRegionNames[region] + "_marginallyUnselectedFakeCaloPhotons_" + fakePhotonCriterionNames[fakeCriterion] + "_";
    if (trueClose_falseAway) name += "closeToTruePhoton_MC_";
    else name += "awayFromTruePhotons_MC_";
    name += MCBinNames[MCBinIndex];
    return name;
  }

  std::string getStatisticsHistogramName(const jetProperty& property, const selectionRegion& region) {
    return ((jetPropertyAttributes[property]).name + "_" + selectionRegionNames[region] + "_selectedCaloJets");
  }

  std::string getStatisticsHistogramName(const jetProperty& property, const selectionRegion& region, const bool& trueClose_falseAway, const int& MCBinIndex) {
    std::string name = (jetPropertyAttributes[property]).name + "_" + selectionRegionNames[region] + "_selectedCaloJets_";
    if (trueClose_falseAway) name += "closeToTruePhoton_MC_";
    else name += "awayFromTruePhotons_MC_";
    name += MCBinNames[MCBinIndex];
    return name;
  }

  std::string getStatisticsHistogramName(const jetProperty& property, const selectionRegion& region, const jetCriterion& criterion) {
    return ((jetPropertyAttributes[property]).name + "_" + selectionRegionNames[region] + "_marginallyUnselectedCaloJets_" + jetCriterionNames[criterion]);
  }

  std::string getStatisticsHistogramName(const jetProperty& property, const selectionRegion& region, const jetCriterion& criterion, const bool& trueClose_falseAway, const int& MCBinIndex) {
    std::string name = (jetPropertyAttributes[property]).name + "_" + selectionRegionNames[region] + "_marginallyUnselectedCaloJets_" + jetCriterionNames[criterion] + "_";
    if (trueClose_falseAway) name += "closeToTruePhoton_MC_";
    else name += "awayFromTruePhotons_MC_";
    name += MCBinNames[MCBinIndex];
    return name;
  }

  void fillStatisticsHistograms(eventProperties& selectedEventPropertiesMap,
                                const bool& isMarginallyUnselectedEvent,
                                unselectedEventProperties& marginallyUnselectedEventPropertiesPair,
                                truthPhotonPropertiesCollection& selectedTruePhotonProperties,
                                photonPropertiesCollection& selectedMediumPhotonProperties,
                                photonPropertiesCollection& selectedMediumPhotonProperties_closeToTruePhoton,
                                photonPropertiesCollection& selectedMediumPhotonProperties_awayFromTruePhoton,
                                unselectedMediumPhotonPropertiesCollection& marginallyUnselectedMediumPhotonProperties,
                                unselectedMediumPhotonPropertiesCollection& marginallyUnselectedMediumPhotonProperties_closeToTruePhoton,
                                unselectedMediumPhotonPropertiesCollection& marginallyUnselectedMediumPhotonProperties_awayFromTruePhoton,
                                photonPropertiesCollection& selectedFakePhotonProperties,
                                photonPropertiesCollection& selectedFakePhotonProperties_closeToTruePhoton,
                                photonPropertiesCollection& selectedFakePhotonProperties_awayFromTruePhoton,
                                unselectedFakePhotonPropertiesCollection& marginallyUnselectedFakePhotonProperties,
                                unselectedFakePhotonPropertiesCollection& marginallyUnselectedFakePhotonProperties_closeToTruePhoton,
                                unselectedFakePhotonPropertiesCollection& marginallyUnselectedFakePhotonProperties_awayFromTruePhoton,
                                jetPropertiesCollection& selectedJetProperties,
                                jetPropertiesCollection& selectedJetProperties_closeToTruePhoton,
                                jetPropertiesCollection& selectedJetProperties_awayFromTruePhoton,
                                unselectedJetPropertiesCollection& marginallyUnselectedJetProperties,
                                unselectedJetPropertiesCollection& marginallyUnselectedJetProperties_closeToTruePhoton,
                                unselectedJetPropertiesCollection& marginallyUnselectedJetProperties_awayFromTruePhoton,
                                selectionRegion& region,
                                const bool& isMC,
                                const int& MCBinIndex) {
    if (region == selectionRegion::nSelectionRegions) return;
    if (!(isMarginallyUnselectedEvent)) {
      for (auto&& selectedEventPropertiesMapElement: selectedEventPropertiesMap) {
        auto& event_property = selectedEventPropertiesMapElement.first;
        auto& value = selectedEventPropertiesMapElement.second;
        fillStatisticsHistogramByName(getStatisticsHistogramName(event_property, region), value);
      }
    }
    else {
      auto& criterion = marginallyUnselectedEventPropertiesPair.first;
      auto& propertiesMap = marginallyUnselectedEventPropertiesPair.second;
      for (auto&& propertiesMapElement: propertiesMap) {
        auto& event_property = propertiesMapElement.first;
        auto& value = propertiesMapElement.second;
        fillStatisticsHistogramByName(getStatisticsHistogramName(event_property, criterion), value);
      }
    }

    for (auto&& selectedMediumPhotonPropertiesMap: selectedMediumPhotonProperties) {
      for (auto&& element: selectedMediumPhotonPropertiesMap) {
        auto& property = element.first;
        auto& value = element.second;
        fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, true), value);
      }
    }

    for (auto&& marginallyUnselectedMediumPhotonPropertiesPair: marginallyUnselectedMediumPhotonProperties) {
      auto& criterion = marginallyUnselectedMediumPhotonPropertiesPair.first;
      auto& propertiesMap = marginallyUnselectedMediumPhotonPropertiesPair.second;
      for (auto&& propertiesMapElement: propertiesMap) {
        auto& property = propertiesMapElement.first;
        auto& value = propertiesMapElement.second;
        fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, criterion), value);
      }
    }

    for (auto&& selectedFakePhotonPropertiesMap: selectedFakePhotonProperties) {
      for (auto&& element: selectedFakePhotonPropertiesMap) {
        auto& property = element.first;
        auto& value = element.second;
        fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, false), value);
      }
    }

    for (auto&& marginallyUnselectedFakePhotonPropertiesPair: marginallyUnselectedFakePhotonProperties) {
      auto& criterion = marginallyUnselectedFakePhotonPropertiesPair.first;
      auto& propertiesMap = marginallyUnselectedFakePhotonPropertiesPair.second;
      for (auto&& propertiesMapElement: propertiesMap) {
        auto& property = propertiesMapElement.first;
        auto& value = propertiesMapElement.second;
        fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, criterion), value);
      }
    }

    for (auto&& selectedJetPropertiesMap: selectedJetProperties) {
      for (auto&& element: selectedJetPropertiesMap) {
        auto& property = element.first;
        auto& value = element.second;
        fillStatisticsHistogramByName(getStatisticsHistogramName(property, region), value);
      }
    }

    for (auto&& marginallyUnselectedJetPropertiesPair: marginallyUnselectedJetProperties) {
      auto& criterion = marginallyUnselectedJetPropertiesPair.first;
      auto& propertiesMap = marginallyUnselectedJetPropertiesPair.second;
      for (auto&& element: propertiesMap) {
        auto& property = element.first;
        auto& value = element.second;
        fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, criterion), value);
      }
    }

    if (isMC && (MCBinIndex > 0)) {
      for (auto&& selectedTruePhotonPropertiesMap: selectedTruePhotonProperties) {
        for (auto&& element: selectedTruePhotonPropertiesMap) {
          auto& property = element.first;
          auto& value = element.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, MCBinIndex), value);
        }
      }

      for (auto&& selectedMediumPhotonProperties_closeToTruePhotonMap: selectedMediumPhotonProperties_closeToTruePhoton) {
        for (auto&& element: selectedMediumPhotonProperties_closeToTruePhotonMap) {
          auto& property = element.first;
          auto& value = element.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, true, true, MCBinIndex), value);
        }
      }

      for (auto&& selectedMediumPhotonProperties_awayFromTruePhotonMap: selectedMediumPhotonProperties_awayFromTruePhoton) {
        for (auto&& element: selectedMediumPhotonProperties_awayFromTruePhotonMap) {
          auto& property = element.first;
          auto& value = element.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, true, true, MCBinIndex), value);
        }
      }

      for (auto&& marginallyUnselectedMediumPhotonProperties_closeToTruePhotonPair: marginallyUnselectedMediumPhotonProperties_closeToTruePhoton) {
        auto& criterion = marginallyUnselectedMediumPhotonProperties_closeToTruePhotonPair.first;
        auto& propertiesMap = marginallyUnselectedMediumPhotonProperties_closeToTruePhotonPair.second;
        for (auto&& propertiesMapElement: propertiesMap) {
          auto& property = propertiesMapElement.first;
          auto& value = propertiesMapElement.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, criterion, true, MCBinIndex), value);
        }
      }

      for (auto&& marginallyUnselectedMediumPhotonProperties_awayFromTruePhotonPair: marginallyUnselectedMediumPhotonProperties_awayFromTruePhoton) {
        mediumPhotonCriterion& criterion = marginallyUnselectedMediumPhotonProperties_awayFromTruePhotonPair.first;
        auto& propertiesMap = marginallyUnselectedMediumPhotonProperties_awayFromTruePhotonPair.second;
        for (auto&& propertiesMapElement: propertiesMap) {
          auto& property = propertiesMapElement.first;
          auto& value = propertiesMapElement.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, criterion, false, MCBinIndex), value);
        }
      }

      for (auto&& selectedFakePhotonProperties_closeToTruePhotonMap: selectedFakePhotonProperties_closeToTruePhoton) {
        for (auto&& element: selectedFakePhotonProperties_closeToTruePhotonMap) {
          auto& property = element.first;
          auto& value = element.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, false, true, MCBinIndex), value);
        }
      }

      for (auto&& selectedFakePhotonProperties_awayFromTruePhotonMap: selectedFakePhotonProperties_awayFromTruePhoton) {
        for (auto&& element: selectedFakePhotonProperties_awayFromTruePhotonMap) {
          auto& property = element.first;
          auto& value = element.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, false, false, MCBinIndex), value);  
        }
      }

      for (auto&& marginallyUnselectedFakePhotonProperties_closeToTruePhotonPair: marginallyUnselectedFakePhotonProperties_closeToTruePhoton) {
        auto& criterion = marginallyUnselectedFakePhotonProperties_closeToTruePhotonPair.first;
        auto& propertiesMap = marginallyUnselectedFakePhotonProperties_closeToTruePhotonPair.second;
        for (auto&& propertiesMapElement: propertiesMap) {
          auto& property = propertiesMapElement.first;
          auto& value = propertiesMapElement.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, criterion, true, MCBinIndex), value);
        }
      }

      for (auto&& marginallyUnselectedFakePhotonProperties_awayFromTruePhotonPair: marginallyUnselectedFakePhotonProperties_awayFromTruePhoton) {
        auto& criterion = marginallyUnselectedFakePhotonProperties_awayFromTruePhotonPair.first;
        auto& propertiesMap = marginallyUnselectedFakePhotonProperties_awayFromTruePhotonPair.second;
        for (auto&& propertiesMapElement: propertiesMap) {
          auto& property = propertiesMapElement.first;
          auto& value = propertiesMapElement.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, criterion, false, MCBinIndex), value);
        }
      }

      for (auto&& selectedJetProperties_closeToTruePhotonMap: selectedJetProperties_closeToTruePhoton) {
        for (auto&& element: selectedJetProperties_closeToTruePhotonMap) {
          auto& property = element.first;
          auto& value = element.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, true, MCBinIndex), value);
        }
      }

      for (auto&& selectedJetProperties_awayFromTruePhotonMap: selectedJetProperties_awayFromTruePhoton) {
        for (auto&& element: selectedJetProperties_awayFromTruePhotonMap) {
          auto& property = element.first;
          auto& value = element.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, false, MCBinIndex), value);
        }
      }

      for (auto&& marginallyUnselectedJetProperties_closeToTruePhotonPair: marginallyUnselectedJetProperties_closeToTruePhoton) {
        auto& criterion = marginallyUnselectedJetProperties_closeToTruePhotonPair.first;
        auto& propertiesMap = marginallyUnselectedJetProperties_closeToTruePhotonPair.second;
        for (auto&& element: propertiesMap) {
          auto& property = element.first;
          auto& value = element.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, criterion, true, MCBinIndex), value);
        }
      }

      for (auto&& marginallyUnselectedJetProperties_awayFromTruePhotonPair: marginallyUnselectedJetProperties_awayFromTruePhoton) {
        auto& criterion = marginallyUnselectedJetProperties_awayFromTruePhotonPair.first;
        auto& propertiesMap = marginallyUnselectedJetProperties_awayFromTruePhotonPair.second;
        for (auto&& element: propertiesMap) {
          auto& property = element.first;
          auto& value = element.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, criterion, false, MCBinIndex), value);
        }
      }
    }
  }

  void writeToFile(const std::string& outputFileRelativePath) {
    TFile *outputFile = TFile::Open(outputFileRelativePath.c_str(), "RECREATE");
    if (!(outputFile->IsOpen()) || outputFile->IsZombie()) {
      std::cout << "ERROR: Unable to open output file to write. File path: " << outputFileRelativePath << std::endl;
    }
    for (auto&& statsElement: stats) {
      outputFile->WriteTObject(&(statsElement.second));
    }
    outputFile->Close();
  }
};

/* List of histograms: */
/* all: "<eventProperty>_<selectionRegion>_selectedEvents" */
/* all: "<eventProperty>_marginallyUnselectedEvents_<eventCriterion>" */
/*  MC: "<truthPhotonProperty>_<selectionRegion>_truePhotons_MCBin<bin>" */
/* *MC: "<truthJetProperty>_<selectionRegion>_trueJets_MCBin<bin>" */

/* all: "<photonProperty>_<selectionRegion>_selectedMediumCaloPhotons" */
/* all: "<photonProperty>_<selectionRegion>_marginallyUnselectedMediumCaloPhotons_<mediumPhotonCriterion>" */
/* all: "<photonProperty>_<selectionRegion>_selectedFakeCaloPhotons" */
/* all: "<photonProperty>_<selectionRegion>_marginallyUnselectedFakeCaloPhotons_<fakePhotonCriterion>" */
/*  MC: "<photonProperty>_<selectionRegion>_selectedMediumCaloPhotons_<close/away_truePhoton>_MCBin<bin>" */
/*  MC: "<photonProperty>_<selectionRegion>_marginallyUnselectedMediumCaloPhotons_<mediumPhotonCriterion>_<close/away_truePhoton>_MCBin<bin>" */
/*  MC: "<photonProperty>_<selectionRegion>_selectedFakeCaloPhotons_<close/away_truePhoton>_MCBin<bin>" */
/*  MC: "<photonProperty>_<selectionRegion>_marginallyUnselectedFakeCaloPhotons_<fakePhotonCriterion>_<close/away_truePhoton>_MCBin<bin>" */
/* *MC: "<photonProperty>_<selectionRegion>_selectedMediumCaloPhotons_<close/away_trueJet>_MCBin<bin>" */
/* *MC: "<photonProperty>_<selectionRegion>_marginallyUnselectedMediumCaloPhotons_<mediumPhotonCriterion>_<close/away_trueJet>_MCBin<bin>" */
/* *MC: "<photonProperty>_<selectionRegion>_selectedFakeCaloPhotons_<close/away_trueJet>_MCBin<bin>" */
/* *MC: "<photonProperty>_<selectionRegion>_marginallyUnselectedFakeCaloPhotons_<fakePhotonCriterion>_<close/away_trueJet>_MCBin<bin>" */

/* all: "<jetProperty>_<selectionRegion>_selectedCaloJets" */
/* all: "<jetProperty>_<selectionRegion>_marginallyUnselectedCaloJets_<jetCriterion>" */
/*  MC: "<jetProperty>_<selectionRegion>_selectedCaloJets_<close/away_truePhoton>_MCBin<bin>" */
/* *MC: "<jetProperty>_<selectionRegion>_selectedCaloJets_<close/away_trueJet>_MCBin<bin>" */
/*  MC: "<jetProperty>_<selectionRegion>_marginallyUnselectedCaloJets_<jetCriterion>_<close/away_truePhoton>_MCBin<bin>" */
/* *MC: "<jetProperty>_<selectionRegion>_marginallyUnselectedCaloJets_<jetCriterion>_<close/away_trueJet>_MCBin<bin>" */

/* improvements: */

/* add true jets */
/* add automatic lines for cuts */
/* add ranges */
/* event weights */

#endif
