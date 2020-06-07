#ifndef H_STATISTICSHISTOGRAMS
#define H_STATISTICSHISTOGRAMS

#include <iostream>
#include <cstdlib>
#include <map>

#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TEfficiency.h"

#include "objectProperties.h"
#include "selectionCriteria.h"
#include "MCRegions.h"

class statisticsHistograms {
 public:
  bool isMC, fillMarginallyUnselectedPhotonJetHistograms;
  std::map<std::string, TH1F*> stats;
  std::map<std::string, TEfficiency*> stats_HLTEfficiency;
  std::map<std::string, TEfficiency*> stats_IDEfficiency;
  float IDEfficiency_STMin = -1.;
  float IDEfficiency_STMax = -1.;

  std::string getStatisticsHistogramName(const eventProperty& event_property, const selectionRegion& region) {
    return ((eventPropertyAttributes.at(event_property)).name + "_" + selectionRegionNames.at(region) + "_selectedEvents");
  }

  std::string getStatisticsHistogramName(const eventProperty& event_property, const selectionRegion& region, const int& MCRegionIndex) {
    return ((eventPropertyAttributes.at(event_property)).name + "_" + selectionRegionNames.at(region) + "_selectedEvents_MC_" + MCRegions::regionNames.at(MCRegionIndex));
  }

  std::string getStatisticsHistogramName(const eventProperty& event_property, const eventSelectionCriterion& criterion) {
    return ((eventPropertyAttributes.at(event_property)).name + "_marginallyUnselectedEvents_" + eventSelectionCriterionNames.at(criterion));
  }

  std::string getStatisticsHistogramName(const eventProperty& event_property, const eventSelectionCriterion& criterion, const int& MCRegionIndex) {
    return ((eventPropertyAttributes.at(event_property)).name + "_marginallyUnselectedEvents_" + eventSelectionCriterionNames.at(criterion) + "_MC_" + MCRegions::regionNames.at(MCRegionIndex));
  }

  std::string getStatisticsHistogramName(const truthPhotonProperty& truth_photon_property, const selectionRegion& region, const int& MCRegionIndex) {
    return ((truthPhotonPropertyAttributes.at(truth_photon_property)).name + "_" + selectionRegionNames.at(region) + "_truePhotons_MC_" + MCRegions::regionNames.at(MCRegionIndex));
  }

  std::string getStatisticsHistogramName(const truthJetCandidateProperty& truth_jetCandidate_property, const selectionRegion& region, const int& MCRegionIndex) {
    return ((truthJetCandidatePropertyAttributes.at(truth_jetCandidate_property)).name + "_" + selectionRegionNames.at(region) + "_trueJetCandidates_all_MC_" + MCRegions::regionNames.at(MCRegionIndex));
  }

  std::string getStatisticsHistogramName(const truthJetCandidateProperty& truth_jetCandidate_property, const selectionRegion& region, const bool& trueFromEventProgenitor_falseFromSinglet, const int& MCRegionIndex) {
    std::string name = (truthJetCandidatePropertyAttributes.at(truth_jetCandidate_property)).name + "_" + selectionRegionNames.at(region) + "_trueJetCandidates_from";
    if (trueFromEventProgenitor_falseFromSinglet) name += "EventProgenitor";
    else name += "Singlet";
    name += "_MC_" + MCRegions::regionNames.at(MCRegionIndex);
    return name;
  }

  std::string getStatisticsHistogramName(const photonProperty& photon_property, const selectionRegion& region, const photonType& photon_type) {
    std::string name = ((photonPropertyAttributes.at(photon_property)).name + "_" + selectionRegionNames.at(region) + "_selected_" + photonTypeNames.at(photon_type) + "_caloPhotons");
    return name;
  }

  std::string getStatisticsHistogramName(const photonProperty& photon_property, const selectionRegion& region, const photonType& photon_type, const bool& trueClose_falseAway, const int& MCRegionIndex) {
    std::string name = ((photonPropertyAttributes.at(photon_property)).name + "_" + selectionRegionNames.at(region) + "_selected_" + photonTypeNames.at(photon_type) + "_caloPhotons_");
    if (trueClose_falseAway) name += "closeToTruePhoton_MC_";
    else name += "awayFromTruePhotons_MC_";
    name += MCRegions::regionNames.at(MCRegionIndex);
    return name;
  }

  std::string getStatisticsHistogramName(const photonProperty& photon_property, const selectionRegion& region, const mediumPhotonCriterion& medium_photon_criterion) {
    return ((photonPropertyAttributes.at(photon_property)).name + "_" + selectionRegionNames.at(region) + "_marginallyUnselectedMediumCaloPhotons_" + mediumPhotonCriterionNames.at(medium_photon_criterion));
  }

  std::string getStatisticsHistogramName(const photonProperty& photon_property, const selectionRegion& region, const mediumPhotonCriterion& medium_photon_criterion, const bool& trueClose_falseAway, const int& MCRegionIndex) {
    std::string name = (photonPropertyAttributes.at(photon_property)).name + "_" + selectionRegionNames.at(region) + "_marginallyUnselectedMediumCaloPhotons_" + mediumPhotonCriterionNames.at(medium_photon_criterion) + "_";
    if (trueClose_falseAway) name += "closeToTruePhoton_MC_";
    else name += "awayFromTruePhotons_MC_";
    name += MCRegions::regionNames.at(MCRegionIndex);
    return name;
  }

  std::string getStatisticsHistogramName(const photonProperty& photon_property, const selectionRegion& region, const vetoedPhotonCriterion& vetoed_photon_criterion) {
    return ((photonPropertyAttributes.at(photon_property)).name + "_" + selectionRegionNames.at(region) + "_marginallyUnselectedVetoedCaloPhotons_" + vetoedPhotonCriterionNames.at(vetoed_photon_criterion));
  }

  std::string getStatisticsHistogramName(const photonProperty& photon_property, const selectionRegion& region, const vetoedPhotonCriterion& vetoed_photon_criterion, const bool& trueClose_falseAway, const int& MCRegionIndex) {
    std::string name = (photonPropertyAttributes.at(photon_property)).name + "_" + selectionRegionNames.at(region) + "_marginallyUnselectedVetoedCaloPhotons_" + vetoedPhotonCriterionNames.at(vetoed_photon_criterion) + "_";
    if (trueClose_falseAway) name += "closeToTruePhoton_MC_";
    else name += "awayFromTruePhotons_MC_";
    name += MCRegions::regionNames.at(MCRegionIndex);
    return name;
  }

  std::string getStatisticsHistogramName(const photonProperty& photon_property, const selectionRegion& region, const fakePhotonCriterion& fake_photon_criterion) {
    return ((photonPropertyAttributes.at(photon_property)).name + "_" + selectionRegionNames.at(region) + "_marginallyUnselectedFakeCaloPhotons_" + fakePhotonCriterionNames.at(fake_photon_criterion));
  }

  std::string getStatisticsHistogramName(const photonProperty& photon_property, const selectionRegion& region, const fakePhotonCriterion& fake_photon_criterion, const bool& trueClose_falseAway, const int& MCRegionIndex) {
    std::string name = (photonPropertyAttributes.at(photon_property)).name + "_" + selectionRegionNames.at(region) + "_marginallyUnselectedFakeCaloPhotons_" + fakePhotonCriterionNames.at(fake_photon_criterion) + "_";
    if (trueClose_falseAway) name += "closeToTruePhoton_MC_";
    else name += "awayFromTruePhotons_MC_";
    name += MCRegions::regionNames.at(MCRegionIndex);
    return name;
  }

  std::string getStatisticsHistogramName(const jetProperty& jet_property, const selectionRegion& region) {
    return ((jetPropertyAttributes.at(jet_property)).name + "_" + selectionRegionNames.at(region) + "_selectedCaloJets");
  }

  std::string getStatisticsHistogramName(const jetProperty& jet_property, const selectionRegion& region, const bool& trueClose_falseAway, const int& MCRegionIndex) {
    std::string name = (jetPropertyAttributes.at(jet_property)).name + "_" + selectionRegionNames.at(region) + "_selectedCaloJets_";
    if (trueClose_falseAway) name += "closeToTruePhoton_MC_";
    else name += "awayFromTruePhotons_MC_";
    name += MCRegions::regionNames.at(MCRegionIndex);
    return name;
  }

  std::string getStatisticsHistogramName(const jetProperty& jet_property, const selectionRegion& region, const jetCriterion& criterion) {
    return ((jetPropertyAttributes.at(jet_property)).name + "_" + selectionRegionNames.at(region) + "_marginallyUnselectedCaloJets_" + jetCriterionNames.at(criterion));
  }

  std::string getStatisticsHistogramName(const jetProperty& jet_property, const selectionRegion& region, const jetCriterion& criterion, const bool& trueClose_falseAway, const int& MCRegionIndex) {
    std::string name = (jetPropertyAttributes.at(jet_property)).name + "_" + selectionRegionNames.at(region) + "_marginallyUnselectedCaloJets_" + jetCriterionNames.at(criterion) + "_";
    if (trueClose_falseAway) name += "closeToTruePhoton_MC_";
    else name += "awayFromTruePhotons_MC_";
    name += MCRegions::regionNames.at(MCRegionIndex);
    return name;
  }

  std::string getStatisticsHistogramName(const genJetProperty& gen_jet_property, const selectionRegion& region, const int& MCRegionIndex) {
    return ((genJetPropertyAttributes.at(gen_jet_property)).name + "_" + selectionRegionNames.at(region) + "_all_genJets_MC_" + MCRegions::regionNames.at(MCRegionIndex));
  }

  std::string getStatisticsHistogramName(const genJetProperty& gen_jet_property, const bool& mom_trueEventProgenitor_falseSinglet, const selectionRegion& region, const int& MCRegionIndex) {
    std::string name = (genJetPropertyAttributes.at(gen_jet_property)).name + "_" + selectionRegionNames.at(region) + "_";
    if (mom_trueEventProgenitor_falseSinglet) name += "eventProgenitorMom_";
    else name += "singletMom_";
    name += "genJets_MC_" + MCRegions::regionNames.at(MCRegionIndex);
    return name;
  }

  void initializeWithCheck(const std::string& name, const int& nBins, const float& xmin, const float& xmax) {
    if (stats.find(name) == stats.end()) {
      stats[name] = new TH1F(name.c_str(), name.c_str(), nBins, xmin, xmax);
      // stats[name].SetCanExtend(TH1::kAllAxes);
    }
    else {
      std::cout << "ERROR: tried to create new 1D statistics histogram with name \"" << name << "\", but it already exists!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  void initializeHLTEfficienciesWithCheck(const std::string& name, const std::vector<double>& xEdges, const std::string& xTitle) {
    if (stats_HLTEfficiency.find(name) == stats_HLTEfficiency.end()) {
      stats_HLTEfficiency[name] = new TEfficiency(name.c_str(), (name + ";" + xTitle + ";").c_str(), (xEdges.size()-1), &(xEdges.at(0)));
    }
    else {
      std::cout << "ERROR: tried to create new 1D HLT efficiency object with name \"" << name << "\", but it already exists!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  void initializeHLTEfficienciesWithCheck(const std::string& name, const std::vector<double>& xEdges, const std::string& xTitle, const std::vector<double>& yEdges, const std::string& yTitle) {
    if (stats_HLTEfficiency.find(name) == stats_HLTEfficiency.end()) {
      stats_HLTEfficiency[name] = new TEfficiency(name.c_str(), (name + ";" + xTitle + ";" + yTitle + ";").c_str(), (xEdges.size()-1), &(xEdges.at(0)), (yEdges.size()-1), &(yEdges.at(0)));
    }
    else {
      std::cout << "ERROR: tried to create new 2D HLT efficiency object with name \"" << name << "\", but it already exists!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  void initializeIDEfficienciesWithCheck(const std::string& name, const std::vector<double>& STBoundaries) {
    if (stats_IDEfficiency.find(name) == stats_IDEfficiency.end()) {
      stats_IDEfficiency[name] = new TEfficiency(name.c_str(), (name + ";ST;efficiency").c_str(), (STBoundaries.size()-1), &(STBoundaries.at(0)));
    }
    else {
      std::cout << "ERROR: tried to create new 1D ID efficiency object with name \"" << name << "\", but it already exists!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  statisticsHistograms(const bool& is_MC, const bool& fill_MarginallyUnselectedPhotonJetHistograms, const std::vector<double>& etaBinEdges, const std::vector<double>& pTBinEdges, const std::vector<double>& STBoundaries) {
    isMC = is_MC;
    fillMarginallyUnselectedPhotonJetHistograms = fill_MarginallyUnselectedPhotonJetHistograms;
    std::vector<double> STBoundariesModified;
    assert(STBoundaries.size() > 1);
    for (int STBoundaryIndex = 0; STBoundaryIndex < (static_cast<int>(STBoundaries.size()) - 1); ++STBoundaryIndex) STBoundariesModified.push_back(STBoundaries.at(STBoundaryIndex)); // Create a copy of STRegionBoundaries with all but the last bin, which is set at some unphysical high value.
    STBoundariesModified.push_back(3500.); // Set the last bin edge to 3500 GeV; this way the plots are readable.
    IDEfficiency_STMin = STBoundariesModified.at(0);
    IDEfficiency_STMax = STBoundariesModified.at(STBoundariesModified.size()-1);

    std::string fullName;
    for (auto&& eventPropertyAttributesElement: eventPropertyAttributes) {
      auto& event_property = (eventPropertyAttributesElement.first);
      auto& attributesElement = (eventPropertyAttributesElement.second);
      for (auto&& selectionRegionNamesElement: selectionRegionNames) {
        auto& region = selectionRegionNamesElement.first;
        fullName = getStatisticsHistogramName(event_property, region);
        initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
        if (isMC) {
          for (auto&& MCRegionNamesElement: MCRegions::regionNames) {
            auto& MCRegionIndex = MCRegionNamesElement.first;
            fullName = getStatisticsHistogramName(event_property, region, MCRegionIndex);
            initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
          }
        }
      }
      for (auto&& eventSelectionCriterionNamesElement: eventSelectionCriterionNames) {
        auto& event_selection_criterion = eventSelectionCriterionNamesElement.first;
        fullName = getStatisticsHistogramName(event_property, event_selection_criterion);
        initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
        if (isMC) {
          for (auto&& MCRegionNamesElement: MCRegions::regionNames) {
            auto& MCRegionIndex = MCRegionNamesElement.first;
            fullName = getStatisticsHistogramName(event_property, event_selection_criterion, MCRegionIndex);
            initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
          }
        }
      }
    }

    if (isMC) {
      for (auto&& truthPhotonPropertyAttributesElement: truthPhotonPropertyAttributes) {
        auto& truth_photon_property = (truthPhotonPropertyAttributesElement.first);
        auto& attributesElement = (truthPhotonPropertyAttributesElement.second);
        for (auto&& selectionRegionNamesElement: selectionRegionNames) {
          auto& region = selectionRegionNamesElement.first;
          for (auto&& MCRegionNamesElement: MCRegions::regionNames) {
            auto& MCRegionIndex = MCRegionNamesElement.first;
            fullName = getStatisticsHistogramName(truth_photon_property, region, MCRegionIndex);
            initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
          }
        }
      }
      for (auto&& jetCandidatePropertyAttributesElement: truthJetCandidatePropertyAttributes) {
        auto& truth_jetCandidate_property = (jetCandidatePropertyAttributesElement.first);
        auto& attributesElement = (jetCandidatePropertyAttributesElement.second);
        for (auto&& selectionRegionNamesElement: selectionRegionNames) {
          auto& region = selectionRegionNamesElement.first;
          for (auto&& MCRegionNamesElement: MCRegions::regionNames) {
            auto& MCRegionIndex = MCRegionNamesElement.first;
            fullName = getStatisticsHistogramName(truth_jetCandidate_property, region, MCRegionIndex);
            initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
            fullName = getStatisticsHistogramName(truth_jetCandidate_property, region, true, MCRegionIndex);
            initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
            fullName = getStatisticsHistogramName(truth_jetCandidate_property, region, false, MCRegionIndex);
            initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
          }
        }
      }
      for (auto&& genJetPropertyAttributesElement: genJetPropertyAttributes) {
        auto& gen_jet_property = (genJetPropertyAttributesElement.first);
        auto& attributesElement = (genJetPropertyAttributesElement.second);
        for (auto&& selectionRegionNamesElement: selectionRegionNames) {
          auto& region = selectionRegionNamesElement.first;
          for (auto&& MCRegionNamesElement: MCRegions::regionNames) {
            auto& MCRegionIndex = MCRegionNamesElement.first;
            fullName = getStatisticsHistogramName(gen_jet_property, region, MCRegionIndex);
            initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
            fullName = getStatisticsHistogramName(gen_jet_property, true, region, MCRegionIndex);
            initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
            fullName = getStatisticsHistogramName(gen_jet_property, false, region, MCRegionIndex);
            initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
          }
        }
      }
    } // end of condition isMC

    for (auto&& selectionRegionNamesElement: selectionRegionNames) {
      auto& region = selectionRegionNamesElement.first;
      for (auto&& photonPropertyAttributesElement: photonPropertyAttributes) {
        auto& photon_property = (photonPropertyAttributesElement.first);
        auto& attributesElement = (photonPropertyAttributesElement.second);
	for (int int_photonTypeIndex = photonTypeFirst; int_photonTypeIndex != static_cast<int>(photonType::nPhotonTypes); ++int_photonTypeIndex) {
	  photonType photonTypeIndex = static_cast<photonType>(int_photonTypeIndex);
	  fullName = getStatisticsHistogramName(photon_property, region, photonTypeIndex);
	  initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
	}
        if (isMC) {
          for (auto&& MCRegionNamesElement: MCRegions::regionNames) {
            auto& MCRegionIndex = MCRegionNamesElement.first;
            for (int int_photonTypeIndex = photonTypeFirst; int_photonTypeIndex != static_cast<int>(photonType::nPhotonTypes); ++int_photonTypeIndex) {
	      photonType photonTypeIndex = static_cast<photonType>(int_photonTypeIndex);
	      fullName = getStatisticsHistogramName(photon_property, region, photonTypeIndex, true, MCRegionIndex);
	      initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
	      fullName = getStatisticsHistogramName(photon_property, region, photonTypeIndex, false, MCRegionIndex);
	      initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
	    }
          }
        } // MC plots
	if (fillMarginallyUnselectedPhotonJetHistograms) {
	  for (auto&& mediumPhotonCriterionNamesElement: mediumPhotonCriterionNames) {
	    auto& medium_photon_criterion = mediumPhotonCriterionNamesElement.first;
	    fullName = getStatisticsHistogramName(photon_property, region, medium_photon_criterion);
	    initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
	    if (isMC) {
	      for (auto&& MCRegionNamesElement: MCRegions::regionNames) {
		auto& MCRegionIndex = MCRegionNamesElement.first;
		fullName = getStatisticsHistogramName(photon_property, region, medium_photon_criterion, true, MCRegionIndex);
		initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
		fullName = getStatisticsHistogramName(photon_property, region, medium_photon_criterion, false, MCRegionIndex);
		initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
	      }
	    } // MC plots
	  } // medium photon criteria
	  for (auto&& vetoedPhotonCriterionNamesElement: vetoedPhotonCriterionNames) {
	    auto& vetoed_photon_criterion = vetoedPhotonCriterionNamesElement.first;
	    fullName = getStatisticsHistogramName(photon_property, region, vetoed_photon_criterion);
	    initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
	    if (isMC) {
	      for (auto&& MCRegionNamesElement: MCRegions::regionNames) {
		auto& MCRegionIndex = MCRegionNamesElement.first;
		fullName = getStatisticsHistogramName(photon_property, region, vetoed_photon_criterion, true, MCRegionIndex);
		initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
		fullName = getStatisticsHistogramName(photon_property, region, vetoed_photon_criterion, false, MCRegionIndex);
		initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
	      }
	    } // MC plots
	  } // vetoed photon criteria
	  for (auto&& fakePhotonCriterionNamesElement: fakePhotonCriterionNames) {
	    auto& fake_photon_criterion = fakePhotonCriterionNamesElement.first;
	    fullName = getStatisticsHistogramName(photon_property, region, fake_photon_criterion);
	    initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
	    if (isMC) {
	      for (auto&& MCRegionNamesElement: MCRegions::regionNames) {
		auto& MCRegionIndex = MCRegionNamesElement.first;
		fullName = getStatisticsHistogramName(photon_property, region, fake_photon_criterion, true, MCRegionIndex);
		initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
		fullName = getStatisticsHistogramName(photon_property, region, fake_photon_criterion, false, MCRegionIndex);
		initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
	      }
	    } // MC plots
	  } // fake photon criteria
	}
      } // photon plots
      for (auto&& jetPropertyAttributesElement: jetPropertyAttributes) {
        auto& jet_property = (jetPropertyAttributesElement.first);
        auto& attributesElement = (jetPropertyAttributesElement.second);
        fullName = getStatisticsHistogramName(jet_property, region);
        initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
        if (isMC) {
          for (auto&& MCRegionNamesElement: MCRegions::regionNames) {
            auto& MCRegionIndex = MCRegionNamesElement.first;
            fullName = getStatisticsHistogramName(jet_property, region, true, MCRegionIndex);
            initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
            fullName = getStatisticsHistogramName(jet_property, region, false, MCRegionIndex);
            initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
          }
        } // MC plots
	if (fillMarginallyUnselectedPhotonJetHistograms) {
	  for (auto&& jetCriterionNamesElement: jetCriterionNames) {
	    auto& jet_criterion = jetCriterionNamesElement.first;
	    fullName = getStatisticsHistogramName(jet_property, region, jet_criterion);
	    initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
	    if (isMC) {
	      for (auto&& MCRegionNamesElement: MCRegions::regionNames) {
		auto& MCRegionIndex = MCRegionNamesElement.first;
		fullName = getStatisticsHistogramName(jet_property, region, jet_criterion, true, MCRegionIndex);
		initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
		fullName = getStatisticsHistogramName(jet_property, region, jet_criterion, false, MCRegionIndex);
		initializeWithCheck(fullName, attributesElement.plot_nBins, attributesElement.plot_minRange, attributesElement.plot_maxRange);
	      }
	    }
	  }
	}
      } // jet plots

      // HLT efficiencies
      fullName = std::string("hltEfficiency_leadingPhoton_" + selectionRegionNames.at(region));
      initializeHLTEfficienciesWithCheck(fullName, etaBinEdges, "eta", pTBinEdges, "pT");
      fullName = std::string("hltEfficiency_subLeadingPhoton_" + selectionRegionNames.at(region));
      initializeHLTEfficienciesWithCheck(fullName, etaBinEdges, "eta", pTBinEdges, "pT");
      fullName = std::string("hltEfficiency_pTBinned_" + selectionRegionNames.at(region));
      initializeHLTEfficienciesWithCheck(fullName, pTBinEdges, "pT_leading", pTBinEdges, "pT_subLeading");
      fullName = std::string("hltEfficiency1D_leadingPhoton_" + selectionRegionNames.at(region));
      initializeHLTEfficienciesWithCheck(fullName, pTBinEdges, "pT_leading");

      // ID efficiencies
      fullName = std::string("IDEfficiency_" + selectionRegionNames.at(region)); // unbinned in nJets
      initializeIDEfficienciesWithCheck(fullName, STBoundariesModified);
      for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) { // hardcoded for now
	fullName = std::string("IDEfficiency_" + std::to_string(nJetsBin) + "Jets_" + selectionRegionNames.at(region));
	initializeIDEfficienciesWithCheck(fullName, STBoundariesModified);
	if (isMC) {
	  for (auto&& MCRegionNamesElement: MCRegions::regionNames) {
	    fullName = std::string("IDEfficiency_" + std::to_string(nJetsBin) + "Jets_" + selectionRegionNames.at(region) + "_MC_" + MCRegionNamesElement.second);
	    initializeIDEfficienciesWithCheck(fullName, STBoundariesModified);
	  }
	}
      }
    } // ends loop over selection regions
  }

  void fillStatisticsHistogramByName(const std::string& histogramName, const float& value, const float& weight) {
    if (stats.find(histogramName) == stats.end()) {
      std::cout << "ERROR: tried to fill statistics histogram with name \"" << histogramName << "\"; a histogram with this name was not initialized!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    else {
      (stats[histogramName])->Fill(value, weight);
    }
  }

  void fillStatisticsHistogramByName(const std::string& histogramName, const float& value) {
    fillStatisticsHistogramByName(histogramName, value, 1.0);
  }

  void fillHLTEfficiencyByName(const std::string& efficiencyName, const bool& passesDiphotonHLT, const float& eta, const float& pT) {
    if (stats_HLTEfficiency.find(efficiencyName) == stats_HLTEfficiency.end()) {
      std::cout << "ERROR: tried to fill HLT efficiency statistics histogram with name \"" << efficiencyName << "\"; a histogram with this name was not initialized!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    else {
      (stats_HLTEfficiency[efficiencyName])->Fill(passesDiphotonHLT, eta, pT);
    }
  }

  void fillPTBinnedHLTEfficiencyByName(const std::string& efficiencyName, const bool& passesDiphotonHLT, const float& pT_leadingPhoton, const float& pT_subLeadingPhoton) {
    if (stats_HLTEfficiency.find(efficiencyName) == stats_HLTEfficiency.end()) {
      std::cout << "ERROR: tried to fill HLT efficiency statistics histogram with name \"" << efficiencyName << "\"; a histogram with this name was not initialized!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    else {
      (stats_HLTEfficiency[efficiencyName])->Fill(passesDiphotonHLT, pT_leadingPhoton, pT_subLeadingPhoton);
    }
  }

  void fillPTBinned1DHLTEfficiencyByName(const std::string& efficiencyName, const bool& passesDiphotonHLT, const float& pT) {
    if (stats_HLTEfficiency.find(efficiencyName) == stats_HLTEfficiency.end()) {
      std::cout << "ERROR: tried to fill HLT efficiency statistics histogram with name \"" << efficiencyName << "\"; a histogram with this name was not initialized!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    else {
      (stats_HLTEfficiency[efficiencyName])->Fill(passesDiphotonHLT, pT);
    }
  }

  void fillIDEfficiencyByName(const std::string& efficiencyName, const bool& passesID, const float& ST) {
    if (stats_IDEfficiency.find(efficiencyName) == stats_IDEfficiency.end()) {
      std::cout << "ERROR: tried to fill ID efficiency statistics histogram with name \"" << efficiencyName << "\"; a histogram with this name was not initialized!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    else {
      (stats_IDEfficiency[efficiencyName])->Fill(passesID, ST);
    }
  }

  void fill1DStatisticsHistograms(eventProperties& selectedEventPropertiesMap,
                                  const bool& isMarginallyUnselectedEvent,
                                  unselectedEventProperties& marginallyUnselectedEventPropertiesPair,
                                  truthPhotonPropertiesCollection& selectedTruePhotonProperties,
                                  truthJetCandidatePropertiesCollection& selectedTrueJetCandidateProperties_all,
                                  truthJetCandidatePropertiesCollection& selectedTrueJetCandidateProperties_fromEventProgenitor,
                                  truthJetCandidatePropertiesCollection& selectedTrueJetCandidateProperties_fromSinglet,
                                  photonPropertiesCollection& selectedMediumPhotonProperties,
                                  photonPropertiesCollection& selectedMediumPhotonProperties_closeToTruePhoton,
                                  photonPropertiesCollection& selectedMediumPhotonProperties_awayFromTruePhoton,
                                  unselectedMediumPhotonPropertiesCollection& marginallyUnselectedMediumPhotonProperties,
                                  unselectedMediumPhotonPropertiesCollection& marginallyUnselectedMediumPhotonProperties_closeToTruePhoton,
                                  unselectedMediumPhotonPropertiesCollection& marginallyUnselectedMediumPhotonProperties_awayFromTruePhoton,
				  photonPropertiesCollection& selectedVetoedPhotonProperties,
                                  photonPropertiesCollection& selectedVetoedPhotonProperties_closeToTruePhoton,
                                  photonPropertiesCollection& selectedVetoedPhotonProperties_awayFromTruePhoton,
                                  unselectedVetoedPhotonPropertiesCollection& marginallyUnselectedVetoedPhotonProperties,
                                  unselectedVetoedPhotonPropertiesCollection& marginallyUnselectedVetoedPhotonProperties_closeToTruePhoton,
                                  unselectedVetoedPhotonPropertiesCollection& marginallyUnselectedVetoedPhotonProperties_awayFromTruePhoton,
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
                                  genJetPropertiesCollection& gen_jet_properties_collection,
                                  genJetPropertiesCollection& eventProgenitor_mom_gen_jet_properties_collection,
                                  genJetPropertiesCollection& singlet_mom_gen_jet_properties_collection,
                                  selectionRegion& region,
                                  const int& MCRegionIndex) {
    if (region == selectionRegion::nSelectionRegions) return;

    if (!(isMarginallyUnselectedEvent)) {
      for (auto&& selectedEventPropertiesMapElement: selectedEventPropertiesMap) {
        auto& event_property = selectedEventPropertiesMapElement.first;
        float& value = selectedEventPropertiesMapElement.second;
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
        fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, photonType::medium), value);
      }
    }

    if (fillMarginallyUnselectedPhotonJetHistograms) {
      for (auto&& marginallyUnselectedMediumPhotonPropertiesPair: marginallyUnselectedMediumPhotonProperties) {
	mediumPhotonCriterion& criterion = marginallyUnselectedMediumPhotonPropertiesPair.first;
	photonProperties& propertiesMap = marginallyUnselectedMediumPhotonPropertiesPair.second;
	for (auto&& propertiesMapElement: propertiesMap) {
	  auto& property = propertiesMapElement.first;
	  auto& value = propertiesMapElement.second;
	  fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, criterion), value);
	}
      }
    }

    for (auto&& selectedVetoedPhotonPropertiesMap: selectedVetoedPhotonProperties) {
      for (auto&& element: selectedVetoedPhotonPropertiesMap) {
        auto& property = element.first;
        auto& value = element.second;
        fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, photonType::vetoed), value);
      }
    }

    if (fillMarginallyUnselectedPhotonJetHistograms) {
      for (auto&& marginallyUnselectedVetoedPhotonPropertiesPair: marginallyUnselectedVetoedPhotonProperties) {
	vetoedPhotonCriterion& criterion = marginallyUnselectedVetoedPhotonPropertiesPair.first;
	photonProperties& propertiesMap = marginallyUnselectedVetoedPhotonPropertiesPair.second;
	for (auto&& propertiesMapElement: propertiesMap) {
	  auto& property = propertiesMapElement.first;
	  auto& value = propertiesMapElement.second;
	  fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, criterion), value);
	}
      }
    }

    for (auto&& selectedFakePhotonPropertiesMap: selectedFakePhotonProperties) {
      for (auto&& element: selectedFakePhotonPropertiesMap) {
        auto& property = element.first;
        auto& value = element.second;
        fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, photonType::fake), value);
      }
    }

    if (fillMarginallyUnselectedPhotonJetHistograms) {
      for (auto&& marginallyUnselectedFakePhotonPropertiesPair: marginallyUnselectedFakePhotonProperties) {
	auto& criterion = marginallyUnselectedFakePhotonPropertiesPair.first;
	auto& propertiesMap = marginallyUnselectedFakePhotonPropertiesPair.second;
	for (auto&& propertiesMapElement: propertiesMap) {
	  auto& property = propertiesMapElement.first;
	  auto& value = propertiesMapElement.second;
	  fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, criterion), value);
	}
      }
    }

    for (auto&& selectedJetPropertiesMap: selectedJetProperties) {
      for (auto&& element: selectedJetPropertiesMap) {
        auto& property = element.first;
        auto& value = element.second;
        fillStatisticsHistogramByName(getStatisticsHistogramName(property, region), value);
      }
    }

    if (fillMarginallyUnselectedPhotonJetHistograms) {
      for (auto&& marginallyUnselectedJetPropertiesPair: marginallyUnselectedJetProperties) {
	auto& criterion = marginallyUnselectedJetPropertiesPair.first;
	auto& propertiesMap = marginallyUnselectedJetPropertiesPair.second;
	for (auto&& element: propertiesMap) {
	  auto& property = element.first;
	  auto& value = element.second;
	  fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, criterion), value);
	}
      }
    }

    if (isMC && (MCRegionIndex > 0)) {
      if (!(isMarginallyUnselectedEvent)) {
        for (auto&& selectedEventPropertiesMapElement: selectedEventPropertiesMap) {
          auto& event_property = selectedEventPropertiesMapElement.first;
          float& value = selectedEventPropertiesMapElement.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(event_property, region, MCRegionIndex), value);
        }
      }
      else {
        auto& criterion = marginallyUnselectedEventPropertiesPair.first;
        auto& propertiesMap = marginallyUnselectedEventPropertiesPair.second;
        for (auto&& propertiesMapElement: propertiesMap) {
          auto& event_property = propertiesMapElement.first;
          auto& value = propertiesMapElement.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(event_property, criterion, MCRegionIndex), value);
        }
      }

      for (auto&& selectedTruePhotonPropertiesMap: selectedTruePhotonProperties) {
        for (auto&& element: selectedTruePhotonPropertiesMap) {
          auto& property = element.first;
          auto& value = element.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, MCRegionIndex), value);
        }
      }

      for (auto&& selectedTrueJetCandidatePropertiesMap: selectedTrueJetCandidateProperties_all) {
        for (auto&& element: selectedTrueJetCandidatePropertiesMap) {
          auto& property = element.first;
          auto& value = element.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, MCRegionIndex), value);
        }
      }

      for (auto&& selectedTrueJetCandidatePropertiesMap: selectedTrueJetCandidateProperties_fromEventProgenitor) {
        for (auto&& element: selectedTrueJetCandidatePropertiesMap) {
          auto& property = element.first;
          auto& value = element.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, true, MCRegionIndex), value);
        }
      }

      for (auto&& selectedTrueJetCandidatePropertiesMap: selectedTrueJetCandidateProperties_fromSinglet) {
        for (auto&& element: selectedTrueJetCandidatePropertiesMap) {
          auto& property = element.first;
          auto& value = element.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, false, MCRegionIndex), value);
        }
      }

      for (auto&& selectedMediumPhotonProperties_closeToTruePhotonMap: selectedMediumPhotonProperties_closeToTruePhoton) {
        for (auto&& element: selectedMediumPhotonProperties_closeToTruePhotonMap) {
          auto& property = element.first;
          auto& value = element.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, photonType::medium, true, MCRegionIndex), value);
        }
      }

      for (auto&& selectedMediumPhotonProperties_awayFromTruePhotonMap: selectedMediumPhotonProperties_awayFromTruePhoton) {
        for (auto&& element: selectedMediumPhotonProperties_awayFromTruePhotonMap) {
          auto& property = element.first;
          auto& value = element.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, photonType::medium, false, MCRegionIndex), value);
        }
      }

      if (fillMarginallyUnselectedPhotonJetHistograms) {
	for (auto&& marginallyUnselectedMediumPhotonProperties_closeToTruePhotonPair: marginallyUnselectedMediumPhotonProperties_closeToTruePhoton) {
	  auto& criterion = marginallyUnselectedMediumPhotonProperties_closeToTruePhotonPair.first;
	  auto& propertiesMap = marginallyUnselectedMediumPhotonProperties_closeToTruePhotonPair.second;
	  for (auto&& propertiesMapElement: propertiesMap) {
	    auto& property = propertiesMapElement.first;
	    auto& value = propertiesMapElement.second;
	    fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, criterion, true, MCRegionIndex), value);
	  }
	}

	for (auto&& marginallyUnselectedMediumPhotonProperties_awayFromTruePhotonPair: marginallyUnselectedMediumPhotonProperties_awayFromTruePhoton) {
	  auto& criterion = marginallyUnselectedMediumPhotonProperties_awayFromTruePhotonPair.first;
	  auto& propertiesMap = marginallyUnselectedMediumPhotonProperties_awayFromTruePhotonPair.second;
	  for (auto&& propertiesMapElement: propertiesMap) {
	    auto& property = propertiesMapElement.first;
	    auto& value = propertiesMapElement.second;
	    fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, criterion, false, MCRegionIndex), value);
	  }
	}
      }

      for (auto&& selectedVetoedPhotonProperties_closeToTruePhotonMap: selectedVetoedPhotonProperties_closeToTruePhoton) {
        for (auto&& element: selectedVetoedPhotonProperties_closeToTruePhotonMap) {
          auto& property = element.first;
          auto& value = element.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, photonType::vetoed, true, MCRegionIndex), value);
        }
      }

      for (auto&& selectedVetoedPhotonProperties_awayFromTruePhotonMap: selectedVetoedPhotonProperties_awayFromTruePhoton) {
        for (auto&& element: selectedVetoedPhotonProperties_awayFromTruePhotonMap) {
          auto& property = element.first;
          auto& value = element.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, photonType::vetoed, false, MCRegionIndex), value);
        }
      }

      if (fillMarginallyUnselectedPhotonJetHistograms) {
	for (auto&& marginallyUnselectedVetoedPhotonProperties_closeToTruePhotonPair: marginallyUnselectedVetoedPhotonProperties_closeToTruePhoton) {
	  auto& criterion = marginallyUnselectedVetoedPhotonProperties_closeToTruePhotonPair.first;
	  auto& propertiesMap = marginallyUnselectedVetoedPhotonProperties_closeToTruePhotonPair.second;
	  for (auto&& propertiesMapElement: propertiesMap) {
	    auto& property = propertiesMapElement.first;
	    auto& value = propertiesMapElement.second;
	    fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, criterion, true, MCRegionIndex), value);
	  }
	}

	for (auto&& marginallyUnselectedVetoedPhotonProperties_awayFromTruePhotonPair: marginallyUnselectedVetoedPhotonProperties_awayFromTruePhoton) {
	  auto& criterion = marginallyUnselectedVetoedPhotonProperties_awayFromTruePhotonPair.first;
	  auto& propertiesMap = marginallyUnselectedVetoedPhotonProperties_awayFromTruePhotonPair.second;
	  for (auto&& propertiesMapElement: propertiesMap) {
	    auto& property = propertiesMapElement.first;
	    auto& value = propertiesMapElement.second;
	    fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, criterion, false, MCRegionIndex), value);
	  }
	}
      }

      for (auto&& selectedFakePhotonProperties_closeToTruePhotonMap: selectedFakePhotonProperties_closeToTruePhoton) {
        for (auto&& element: selectedFakePhotonProperties_closeToTruePhotonMap) {
          auto& property = element.first;
          auto& value = element.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, photonType::fake, true, MCRegionIndex), value);
        }
      }

      for (auto&& selectedFakePhotonProperties_awayFromTruePhotonMap: selectedFakePhotonProperties_awayFromTruePhoton) {
        for (auto&& element: selectedFakePhotonProperties_awayFromTruePhotonMap) {
          auto& property = element.first;
          auto& value = element.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, photonType::fake, false, MCRegionIndex), value);  
        }
      }

      if (fillMarginallyUnselectedPhotonJetHistograms) {
	for (auto&& marginallyUnselectedFakePhotonProperties_closeToTruePhotonPair: marginallyUnselectedFakePhotonProperties_closeToTruePhoton) {
	  auto& criterion = marginallyUnselectedFakePhotonProperties_closeToTruePhotonPair.first;
	  auto& propertiesMap = marginallyUnselectedFakePhotonProperties_closeToTruePhotonPair.second;
	  for (auto&& propertiesMapElement: propertiesMap) {
	    auto& property = propertiesMapElement.first;
	    auto& value = propertiesMapElement.second;
	    fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, criterion, true, MCRegionIndex), value);
	  }
	}

	for (auto&& marginallyUnselectedFakePhotonProperties_awayFromTruePhotonPair: marginallyUnselectedFakePhotonProperties_awayFromTruePhoton) {
	  auto& criterion = marginallyUnselectedFakePhotonProperties_awayFromTruePhotonPair.first;
	  auto& propertiesMap = marginallyUnselectedFakePhotonProperties_awayFromTruePhotonPair.second;
	  for (auto&& propertiesMapElement: propertiesMap) {
	    auto& property = propertiesMapElement.first;
	    auto& value = propertiesMapElement.second;
	    fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, criterion, false, MCRegionIndex), value);
	  }
	}
      }

      for (auto&& selectedJetProperties_closeToTruePhotonMap: selectedJetProperties_closeToTruePhoton) {
        for (auto&& element: selectedJetProperties_closeToTruePhotonMap) {
          auto& property = element.first;
          auto& value = element.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, true, MCRegionIndex), value);
        }
      }

      for (auto&& selectedJetProperties_awayFromTruePhotonMap: selectedJetProperties_awayFromTruePhoton) {
        for (auto&& element: selectedJetProperties_awayFromTruePhotonMap) {
          auto& property = element.first;
          auto& value = element.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, false, MCRegionIndex), value);
        }
      }

      if (fillMarginallyUnselectedPhotonJetHistograms) {
	for (auto&& marginallyUnselectedJetProperties_closeToTruePhotonPair: marginallyUnselectedJetProperties_closeToTruePhoton) {
	  auto& criterion = marginallyUnselectedJetProperties_closeToTruePhotonPair.first;
	  auto& propertiesMap = marginallyUnselectedJetProperties_closeToTruePhotonPair.second;
	  for (auto&& element: propertiesMap) {
	    auto& property = element.first;
	    auto& value = element.second;
	    fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, criterion, true, MCRegionIndex), value);
	  }
	}

	for (auto&& marginallyUnselectedJetProperties_awayFromTruePhotonPair: marginallyUnselectedJetProperties_awayFromTruePhoton) {
	  auto& criterion = marginallyUnselectedJetProperties_awayFromTruePhotonPair.first;
	  auto& propertiesMap = marginallyUnselectedJetProperties_awayFromTruePhotonPair.second;
	  for (auto&& element: propertiesMap) {
	    auto& property = element.first;
	    auto& value = element.second;
	    fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, criterion, false, MCRegionIndex), value);
	  }
	}
      }

      for (auto& gen_jet_properties_map: gen_jet_properties_collection) {
        for (auto&& element: gen_jet_properties_map) {
          auto& property = element.first;
          auto& value = element.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, region, MCRegionIndex), value);
        }
      }

      for (auto& eventProgenitor_mom_gen_jet_properties_map: eventProgenitor_mom_gen_jet_properties_collection) {
        for (auto&& element: eventProgenitor_mom_gen_jet_properties_map) {
          auto& property = element.first;
          auto& value = element.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, true, region, MCRegionIndex), value);
        }
      }

      for (auto& singlet_mom_gen_jet_properties_map: singlet_mom_gen_jet_properties_collection) {
        for (auto&& element: singlet_mom_gen_jet_properties_map) {
          auto& property = element.first;
          auto& value = element.second;
          fillStatisticsHistogramByName(getStatisticsHistogramName(property, false, region, MCRegionIndex), value);
        }
      }
    }
  }

  void fillHLTEfficiencyStatisticsHistograms(const float& eta_leadingPhoton,
                                            const float& pT_leadingPhoton,
                                            const float& eta_subLeadingPhoton,
                                            const float& pT_subLeadingPhoton,
                                            const bool& passesDiphotonHLT,
                                            const selectionRegion& region) {
    // HLT efficiency statistics
    fillHLTEfficiencyByName(std::string("hltEfficiency_leadingPhoton_" + selectionRegionNames.at(region)), passesDiphotonHLT, eta_leadingPhoton, pT_leadingPhoton);
    fillHLTEfficiencyByName(std::string("hltEfficiency_subLeadingPhoton_" + selectionRegionNames.at(region)), passesDiphotonHLT, eta_subLeadingPhoton, pT_subLeadingPhoton);
    fillPTBinnedHLTEfficiencyByName(std::string("hltEfficiency_pTBinned_" + selectionRegionNames.at(region)), passesDiphotonHLT, pT_leadingPhoton, pT_subLeadingPhoton);
    fillPTBinned1DHLTEfficiencyByName(std::string("hltEfficiency1D_leadingPhoton_" + selectionRegionNames.at(region)), passesDiphotonHLT, pT_leadingPhoton);
  }

  void fillIDEfficiencyStatisticsHistograms(const float& eventST, const int& nJetsDR, const bool& passesEventSelection, const selectionRegion& eventRegion, const int& MCRegionIndex) {
    if (eventST <= IDEfficiency_STMin) return;
    if (eventST >= IDEfficiency_STMax) return;
    if (nJetsDR < 2) return;

    int nJetsBin = nJetsDR;
    if (nJetsBin > 6) nJetsBin = 6;
    int eventRegionInt = static_cast<int>(eventRegion);
    for (int regionIndex = selectionRegionFirst; regionIndex < static_cast<int>(selectionRegion::nSelectionRegions); ++regionIndex) {
      selectionRegion region = static_cast<selectionRegion>(regionIndex);
      fillIDEfficiencyByName(std::string("IDEfficiency_" + selectionRegionNames.at(region)), ((regionIndex == eventRegionInt) && passesEventSelection), eventST);
      fillIDEfficiencyByName(std::string("IDEfficiency_" + std::to_string(nJetsBin) + "Jets_" + selectionRegionNames.at(region)), ((regionIndex == eventRegionInt) && passesEventSelection), eventST);
      if (isMC) {
	if (MCRegionIndex == 0) continue;
	fillIDEfficiencyByName(std::string("IDEfficiency_" + std::to_string(nJetsBin) + "Jets_" + selectionRegionNames.at(region) + "_MC_" + MCRegions::regionNames.at(MCRegionIndex)), ((regionIndex == eventRegionInt) && passesEventSelection), eventST);
      }
    }
  }

  void writeToFile(const std::string& outputFileRelativePath) {
    TFile *outputFile = TFile::Open(outputFileRelativePath.c_str(), "RECREATE");
    if (!(outputFile->IsOpen()) || outputFile->IsZombie()) {
      std::cout << "ERROR: Unable to open output file to write. File path: " << outputFileRelativePath << std::endl;
    }
    for (auto&& statsElement: stats) {
      outputFile->WriteTObject(statsElement.second);
    }
    for (auto&& statsElement: stats_HLTEfficiency) {
      outputFile->WriteTObject(statsElement.second);
    }
    for (auto&& IDElement: stats_IDEfficiency) {
      outputFile->WriteTObject(IDElement.second);
    }
    outputFile->Close();
  }
};

#endif
