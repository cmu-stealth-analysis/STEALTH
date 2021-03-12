#include "../include/getMCUncertainties.h"

outputHistogramsStruct* initializeOutputHistograms(optionsStruct& options, MCTemplateReader& templateReader, const STRegionsStruct& STRegions) {
  outputHistogramsStruct* outputHistograms = new outputHistogramsStruct();
  for (int STRegionIndex = 1; STRegionIndex <= (1+STRegions.nSTSignalBins); ++STRegionIndex) {
    for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
      if ((nJetsBin <= 3 || STRegionIndex == 1) || options.getSignalContaminationOutsideSidebands) { // Signal contamination is to be calculated only in the low nJets sideband or at all nJets in the normalization bin
        outputHistograms->h_signalContamination[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("signalContamination", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("signalContamination", STRegionIndex, nJetsBin, STRegions, "").c_str(), templateReader.nEventProgenitorMassBins, templateReader.minEventProgenitorMass, templateReader.maxEventProgenitorMass, templateReader.nNeutralinoMassBins, templateReader.minNeutralinoMass, templateReader.maxNeutralinoMass);
      }
      std::vector<std::string> UpDownShifts = {"Down", "Up"};
      for (const std::string& UpDownShift: UpDownShifts) {
        if(nJetsBin >= 4) { // the rest of the plots are only useful in the signal bins
          outputHistograms->h_MCStatisticsFractionalError[STRegionIndex][nJetsBin][UpDownShift] = new TH2F(("h_" + getHistogramName("MCStatisticsFractionalError" + UpDownShift, STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("MCStatisticsFractionalError", STRegionIndex, nJetsBin, STRegions, UpDownShift).c_str(), templateReader.nEventProgenitorMassBins, templateReader.minEventProgenitorMass, templateReader.maxEventProgenitorMass, templateReader.nNeutralinoMassBins, templateReader.minNeutralinoMass, templateReader.maxNeutralinoMass);
          outputHistograms->h_JECUncertainty[STRegionIndex][nJetsBin][UpDownShift] = new TH2F(("h_" + getHistogramName("JECUncertainty" + UpDownShift, STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("JECUncertainty", STRegionIndex, nJetsBin, STRegions, UpDownShift).c_str(), templateReader.nEventProgenitorMassBins, templateReader.minEventProgenitorMass, templateReader.maxEventProgenitorMass, templateReader.nNeutralinoMassBins, templateReader.minNeutralinoMass, templateReader.maxNeutralinoMass);
          outputHistograms->h_UnclusteredMETUncertainty[STRegionIndex][nJetsBin][UpDownShift] = new TH2F(("h_" + getHistogramName("UnclusteredMETUncertainty" + UpDownShift, STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("UnclusteredMETUncertainty", STRegionIndex, nJetsBin, STRegions, UpDownShift).c_str(), templateReader.nEventProgenitorMassBins, templateReader.minEventProgenitorMass, templateReader.maxEventProgenitorMass, templateReader.nNeutralinoMassBins, templateReader.minNeutralinoMass, templateReader.maxNeutralinoMass);
          outputHistograms->h_JERMETUncertainty[STRegionIndex][nJetsBin][UpDownShift] = new TH2F(("h_" + getHistogramName("JERMETUncertainty" + UpDownShift, STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("JERMETUncertainty", STRegionIndex, nJetsBin, STRegions, UpDownShift).c_str(), templateReader.nEventProgenitorMassBins, templateReader.minEventProgenitorMass, templateReader.maxEventProgenitorMass, templateReader.nNeutralinoMassBins, templateReader.minNeutralinoMass, templateReader.maxNeutralinoMass);
	  outputHistograms->h_missingHEMUncertainty[STRegionIndex][nJetsBin][UpDownShift] = new TH2F(("h_" + getHistogramName("missingHEMUncertainty" + UpDownShift, STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("missingHEMUncertainty", STRegionIndex, nJetsBin, STRegions, UpDownShift).c_str(), templateReader.nEventProgenitorMassBins, templateReader.minEventProgenitorMass, templateReader.maxEventProgenitorMass, templateReader.nNeutralinoMassBins, templateReader.minNeutralinoMass, templateReader.maxNeutralinoMass);
          outputHistograms->h_prefiringWeightsUncertainty[STRegionIndex][nJetsBin][UpDownShift] = new TH2F(("h_" + getHistogramName("prefiringWeightsUncertainty" + UpDownShift, STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("prefiringWeightsUncertainty", STRegionIndex, nJetsBin, STRegions, UpDownShift).c_str(), templateReader.nEventProgenitorMassBins, templateReader.minEventProgenitorMass, templateReader.maxEventProgenitorMass, templateReader.nNeutralinoMassBins, templateReader.minNeutralinoMass, templateReader.maxNeutralinoMass);
          outputHistograms->h_HLTUncertainty[STRegionIndex][nJetsBin][UpDownShift] = new TH2F(("h_" + getHistogramName("HLTUncertainty" + UpDownShift, STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("HLTUncertainty", STRegionIndex, nJetsBin, STRegions, UpDownShift).c_str(), templateReader.nEventProgenitorMassBins, templateReader.minEventProgenitorMass, templateReader.maxEventProgenitorMass, templateReader.nNeutralinoMassBins, templateReader.minNeutralinoMass, templateReader.maxNeutralinoMass);
          outputHistograms->h_photonMCScaleFactorUncertainty[STRegionIndex][nJetsBin][UpDownShift] = new TH2F(("h_" + getHistogramName("photonMCScaleFactorUncertainty" + UpDownShift, STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("photonMCScaleFactorUncertainty", STRegionIndex, nJetsBin, STRegions, UpDownShift).c_str(), templateReader.nEventProgenitorMassBins, templateReader.minEventProgenitorMass, templateReader.maxEventProgenitorMass, templateReader.nNeutralinoMassBins, templateReader.minNeutralinoMass, templateReader.maxNeutralinoMass);
        }
      } // end loop over up or down shifts
    } // end loop over nJetsBin
  } // end loop over STRegionIndex
  return outputHistograms;
}

inputHistogramsStruct* readInputHistograms(TFile *inputFile, const STRegionsStruct& STRegions) {
  inputHistogramsStruct *inputHistograms = new inputHistogramsStruct();
  for (int STRegionIndex = 1; STRegionIndex <= (1+STRegions.nSTSignalBins); ++STRegionIndex) {
    for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
      std::stringstream nameStreamTotalNEvents;
      nameStreamTotalNEvents << "h_totalNEvents_STRegion" << STRegionIndex << "_" << nJetsBin << "Jets";
      inputHistograms->h_totalNEvents[STRegionIndex][nJetsBin] = (TH2I*)(inputFile->Get(nameStreamTotalNEvents.str().c_str()));
      if (inputHistograms->h_totalNEvents[STRegionIndex][nJetsBin] == nullptr) {
        std::cout << "Unable to find histogram with name " << nameStreamTotalNEvents.str() << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if(inputHistograms->h_totalNEvents[STRegionIndex][nJetsBin]->GetBinErrorOption() != TH1::EBinErrorOpt::kPoisson) {
        std::cout << "ERROR: errors on histogram of total number of events not Poisson." << std::endl;
        std::exit(EXIT_FAILURE);
      }

      for (int shiftTypeIndex = shiftTypeFirst; shiftTypeIndex != static_cast<int>(shiftType::nShiftTypes); ++shiftTypeIndex) {
        shiftType typeIndex = static_cast<shiftType>(shiftTypeIndex);
        std::stringstream nameStreamTotalNEventsShifted;
        nameStreamTotalNEventsShifted << "h_totalNEvents_shifted_" << shiftTypeNames[typeIndex] << "_STRegion" << STRegionIndex << "_" << nJetsBin << "Jets";
        inputHistograms->h_totalNEvents_shifted[typeIndex][STRegionIndex][nJetsBin] = (TH2I*)(inputFile->Get(nameStreamTotalNEventsShifted.str().c_str()));
        if (inputHistograms->h_totalNEvents_shifted[typeIndex][STRegionIndex][nJetsBin] == nullptr) {
          std::cout << "Unable to find histogram with name " << nameStreamTotalNEventsShifted.str() << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }

      std::string commonPrefix = "h_lumiBasedYearWeightedNEvents";
      std::stringstream commonSuffixStringStream;
      commonSuffixStringStream << "_STRegion" << STRegionIndex << "_" << nJetsBin << "Jets";
      std::string commonSuffix = commonSuffixStringStream.str();

      inputHistograms->h_lumiBasedYearWeightedNEvents[STRegionIndex][nJetsBin] = (TH2F*)(inputFile->Get((commonPrefix + commonSuffix).c_str()));
      if (inputHistograms->h_lumiBasedYearWeightedNEvents[STRegionIndex][nJetsBin] == nullptr) {
        std::cout << "Unable to find histogram with name " << (commonPrefix + commonSuffix) << std::endl;
        std::exit(EXIT_FAILURE);
      }

      inputHistograms->h_lumiBasedYearWeightedNEvents_prefiringDown[STRegionIndex][nJetsBin] = (TH2F*)(inputFile->Get((commonPrefix + "_prefiringDown" + commonSuffix).c_str()));
      if (inputHistograms->h_lumiBasedYearWeightedNEvents_prefiringDown[STRegionIndex][nJetsBin] == nullptr) {
        std::cout << "Unable to find histogram with name " << (commonPrefix + "_prefiringDown" + commonSuffix) << std::endl;
        std::exit(EXIT_FAILURE);
      }
      inputHistograms->h_lumiBasedYearWeightedNEvents_prefiringUp[STRegionIndex][nJetsBin] = (TH2F*)(inputFile->Get((commonPrefix + "_prefiringUp" + commonSuffix).c_str()));
      if (inputHistograms->h_lumiBasedYearWeightedNEvents_prefiringUp[STRegionIndex][nJetsBin] == nullptr) {
        std::cout << "Unable to find histogram with name " << (commonPrefix + "_prefiringUp" + commonSuffix) << std::endl;
        std::exit(EXIT_FAILURE);
      }
      inputHistograms->h_lumiBasedYearWeightedNEvents_HLTDown[STRegionIndex][nJetsBin] = (TH2F*)(inputFile->Get((commonPrefix + "_HLTDown" + commonSuffix).c_str()));
      if (inputHistograms->h_lumiBasedYearWeightedNEvents_HLTDown[STRegionIndex][nJetsBin] == nullptr) {
        std::cout << "Unable to find histogram with name " << (commonPrefix + "_HLTDown" + commonSuffix) << std::endl;
        std::exit(EXIT_FAILURE);
      }
      inputHistograms->h_lumiBasedYearWeightedNEvents_HLTUp[STRegionIndex][nJetsBin] = (TH2F*)(inputFile->Get((commonPrefix + "_HLTUp" + commonSuffix).c_str()));
      if (inputHistograms->h_lumiBasedYearWeightedNEvents_HLTUp[STRegionIndex][nJetsBin] == nullptr) {
        std::cout << "Unable to find histogram with name " << (commonPrefix + "_HLTUp" + commonSuffix) << std::endl;
        std::exit(EXIT_FAILURE);
      }
      inputHistograms->h_lumiBasedYearWeightedNEvents_photonScaleFactorDown[STRegionIndex][nJetsBin] = (TH2F*)(inputFile->Get((commonPrefix + "_photonScaleFactorDown" + commonSuffix).c_str()));
      if (inputHistograms->h_lumiBasedYearWeightedNEvents_photonScaleFactorDown[STRegionIndex][nJetsBin] == nullptr) {
        std::cout << "Unable to find histogram with name " << (commonPrefix + "_photonScaleFactorDown" + commonSuffix) << std::endl;
        std::exit(EXIT_FAILURE);
      }
      inputHistograms->h_lumiBasedYearWeightedNEvents_photonScaleFactorUp[STRegionIndex][nJetsBin] = (TH2F*)(inputFile->Get((commonPrefix + "_photonScaleFactorUp" + commonSuffix).c_str()));
      if (inputHistograms->h_lumiBasedYearWeightedNEvents_photonScaleFactorUp[STRegionIndex][nJetsBin] == nullptr) {
        std::cout << "Unable to find histogram with name " << (commonPrefix + "_photonScaleFactorUp" + commonSuffix) << std::endl;
        std::exit(EXIT_FAILURE);
      }
    } // ends loop over NJetsBin
  } // ends loop over ST regions
  return inputHistograms;
}

void fillHistogramsWithAsymmetricErrorsFromWeightedNEvents(std::map<std::string, TH2F*> &histogramsToFill, const float& weightedNEventsDown, const float& weightedNEventsUp, const float& weightedNEventsNominal, const float& eventProgenitorMass, const float& neutralinoMass) {
  float fractionalErrorDown = (weightedNEventsDown-weightedNEventsNominal)/weightedNEventsNominal;
  histogramsToFill["Down"]->SetBinContent(histogramsToFill["Down"]->FindFixBin(eventProgenitorMass, neutralinoMass), fractionalErrorDown);
  float fractionalErrorUp = (weightedNEventsUp-weightedNEventsNominal)/weightedNEventsNominal;
  histogramsToFill["Up"]->SetBinContent(histogramsToFill["Up"]->FindFixBin(eventProgenitorMass, neutralinoMass), fractionalErrorUp);
}

void fillHistogramsWithAsymmetricErrorsFromNEvents(std::map<std::string, TH2F*> &histogramsToFill, const int& nEventsDown, const int& nEventsUp, const int& nEventsNominal, const float& eventProgenitorMass, const float& neutralinoMass) {
  fillHistogramsWithAsymmetricErrorsFromWeightedNEvents(histogramsToFill, static_cast<float>(nEventsDown), static_cast<float>(nEventsUp), static_cast<float>(nEventsNominal), eventProgenitorMass, neutralinoMass);
}

void fillSystematicsHistograms(outputHistogramsStruct *outputHistograms, optionsStruct& options, MCTemplateReader& templateReader, const STRegionsStruct& STRegions, inputNEventsStruct& inputNEvents// , inputDataUncertaintiesStruct& inputDataUncertainties, inputDataSTScalingUncertaintiesStruct& inputDataSTScalingUncertainties
                               ) {
  TFile *inputFile = TFile::Open(options.inputPath.c_str(), "READ");
  if (inputFile->IsZombie() || !(inputFile->IsOpen())) {
    std::cout << "Error in opening file " << options.inputPath << std::endl;
    std::exit(EXIT_FAILURE);
  }
  inputHistogramsStruct *inputHistograms = readInputHistograms(inputFile, STRegions);

  std::cout << "Getting systematics..." << std::endl;

  // Fill TGraphs and TH2s with the JEC fractional uncertainty and estimated error
  int nProblematicBins = 0;
  for (int STRegionIndex = 1; STRegionIndex <= (1+STRegions.nSTSignalBins); ++STRegionIndex) {
    for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
      for (int eventProgenitorBinIndex = 1; eventProgenitorBinIndex <= templateReader.nEventProgenitorMassBins; ++eventProgenitorBinIndex) {
        float eventProgenitorMass = (templateReader.eventProgenitorMasses).at(eventProgenitorBinIndex);
        for (int neutralinoBinIndex = 1; neutralinoBinIndex <= templateReader.nNeutralinoMassBins; ++neutralinoBinIndex) {
          if (!(templateReader.isValidBin(eventProgenitorBinIndex, neutralinoBinIndex))) continue;
          float neutralinoMass = (templateReader.neutralinoMasses).at(neutralinoBinIndex);
          float weightedNEvents_nominal = inputHistograms->h_lumiBasedYearWeightedNEvents[STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_lumiBasedYearWeightedNEvents[STRegionIndex][nJetsBin]->FindFixBin(eventProgenitorMass, neutralinoMass));
          // First the signal contamination plots
          if (((nJetsBin == 2) || (STRegionIndex == 1)) || options.getSignalContaminationOutsideSidebands) {
            int nBackgroundEvts = (inputNEvents.data).at(std::string("observedNEvents_STRegion" + std::to_string(STRegionIndex) + "_" + std::to_string(nJetsBin) + "Jets"));
            if (nBackgroundEvts > 0) {
              outputHistograms->h_signalContamination[STRegionIndex][nJetsBin]->SetBinContent(outputHistograms->h_signalContamination[STRegionIndex][nJetsBin]->FindFixBin(eventProgenitorMass, neutralinoMass), weightedNEvents_nominal/nBackgroundEvts);
            }
            else {
              outputHistograms->h_signalContamination[STRegionIndex][nJetsBin]->SetBinContent(outputHistograms->h_signalContamination[STRegionIndex][nJetsBin]->FindFixBin(eventProgenitorMass, neutralinoMass), weightedNEvents_nominal/0.693);
            }
          }
          // Systematics plots
          if(nJetsBin >= 4) {
            bool zeroMCEventsRecorded = false;
            int totalNEvents_nominal = inputHistograms->h_totalNEvents[STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_totalNEvents[STRegionIndex][nJetsBin]->FindFixBin(eventProgenitorMass, neutralinoMass));

            if (totalNEvents_nominal == 0) {
              std::cout << "WARNING: zero events recorded at eventProgenitor mass: " << eventProgenitorMass << ", neutralino mass: " << neutralinoMass << ", for STRegionIndex: " << STRegionIndex << ", nJets: " << nJetsBin << std::endl;
              zeroMCEventsRecorded = true;
            }

            if ((weightedNEvents_nominal == 0) && !(zeroMCEventsRecorded)) { // sanity check
              std::cout << "WARNING: total number of recorded events is nonzero but weighted number of recorded events is 0 at eventProgenitor mass: " << eventProgenitorMass << ", neutralino mass: " << neutralinoMass << ", for STRegionIndex: " << STRegionIndex << ", nJets: " << nJetsBin << std::endl;
              ++nProblematicBins;
            }

            if (zeroMCEventsRecorded || (weightedNEvents_nominal == 0)) {
              std::vector<std::string> UpDownShifts = {"Down", "Up"};
              std::map<std::string, float> UpDownShiftMultiplier = {
                {"Down", -1.0},
                {"Up", 1.0}
              };

	      // "Pretend" that one event is observed -- should make no difference to any signal prediction
	      outputHistograms->h_MCStatisticsFractionalError[STRegionIndex][nJetsBin]["Down"]->SetBinContent(outputHistograms->h_MCStatisticsFractionalError[STRegionIndex][nJetsBin]["Down"]->FindFixBin(eventProgenitorMass, neutralinoMass), -0.827);
	      outputHistograms->h_MCStatisticsFractionalError[STRegionIndex][nJetsBin]["Up"]->SetBinContent(outputHistograms->h_MCStatisticsFractionalError[STRegionIndex][nJetsBin]["Up"]->FindFixBin(eventProgenitorMass, neutralinoMass), 2.3);
	      // Fill the rest with some numbers -- DEFAULT_FRACTIONAL_ERROR is actually 0 in the latest iteration
              for (const std::string& UpDownShift: UpDownShifts) {
                outputHistograms->h_JECUncertainty[STRegionIndex][nJetsBin][UpDownShift]->SetBinContent(outputHistograms->h_JECUncertainty[STRegionIndex][nJetsBin][UpDownShift]->FindFixBin(eventProgenitorMass, neutralinoMass), UpDownShiftMultiplier.at(UpDownShift)*DEFAULT_FRACTIONAL_ERROR);
                outputHistograms->h_UnclusteredMETUncertainty[STRegionIndex][nJetsBin][UpDownShift]->SetBinContent(outputHistograms->h_UnclusteredMETUncertainty[STRegionIndex][nJetsBin][UpDownShift]->FindFixBin(eventProgenitorMass, neutralinoMass), UpDownShiftMultiplier.at(UpDownShift)*DEFAULT_FRACTIONAL_ERROR);
                outputHistograms->h_JERMETUncertainty[STRegionIndex][nJetsBin][UpDownShift]->SetBinContent(outputHistograms->h_JERMETUncertainty[STRegionIndex][nJetsBin][UpDownShift]->FindFixBin(eventProgenitorMass, neutralinoMass), UpDownShiftMultiplier.at(UpDownShift)*DEFAULT_FRACTIONAL_ERROR);
		outputHistograms->h_missingHEMUncertainty[STRegionIndex][nJetsBin][UpDownShift]->SetBinContent(outputHistograms->h_missingHEMUncertainty[STRegionIndex][nJetsBin][UpDownShift]->FindFixBin(eventProgenitorMass, neutralinoMass), UpDownShiftMultiplier.at(UpDownShift)*DEFAULT_FRACTIONAL_ERROR);
                outputHistograms->h_prefiringWeightsUncertainty[STRegionIndex][nJetsBin][UpDownShift]->SetBinContent(outputHistograms->h_prefiringWeightsUncertainty[STRegionIndex][nJetsBin][UpDownShift]->FindFixBin(eventProgenitorMass, neutralinoMass), UpDownShiftMultiplier.at(UpDownShift)*DEFAULT_FRACTIONAL_ERROR);
                outputHistograms->h_HLTUncertainty[STRegionIndex][nJetsBin][UpDownShift]->SetBinContent(outputHistograms->h_HLTUncertainty[STRegionIndex][nJetsBin][UpDownShift]->FindFixBin(eventProgenitorMass, neutralinoMass), UpDownShiftMultiplier.at(UpDownShift)*DEFAULT_FRACTIONAL_ERROR);
                outputHistograms->h_photonMCScaleFactorUncertainty[STRegionIndex][nJetsBin][UpDownShift]->SetBinContent(outputHistograms->h_photonMCScaleFactorUncertainty[STRegionIndex][nJetsBin][UpDownShift]->FindFixBin(eventProgenitorMass, neutralinoMass), UpDownShiftMultiplier.at(UpDownShift)*DEFAULT_FRACTIONAL_ERROR);
              }
            }
            else {
              double totalNEventsErrorDown = (-1.0)*inputHistograms->h_totalNEvents[STRegionIndex][nJetsBin]->GetBinErrorLow(inputHistograms->h_totalNEvents[STRegionIndex][nJetsBin]->FindFixBin(eventProgenitorMass, neutralinoMass)); // Note inversion in sign
              double totalNEventsErrorUp = inputHistograms->h_totalNEvents[STRegionIndex][nJetsBin]->GetBinErrorUp(inputHistograms->h_totalNEvents[STRegionIndex][nJetsBin]->FindFixBin(eventProgenitorMass, neutralinoMass));
              outputHistograms->h_MCStatisticsFractionalError[STRegionIndex][nJetsBin]["Down"]->SetBinContent(outputHistograms->h_MCStatisticsFractionalError[STRegionIndex][nJetsBin]["Down"]->FindFixBin(eventProgenitorMass, neutralinoMass), totalNEventsErrorDown/(static_cast<float>(totalNEvents_nominal)));
              outputHistograms->h_MCStatisticsFractionalError[STRegionIndex][nJetsBin]["Up"]->SetBinContent(outputHistograms->h_MCStatisticsFractionalError[STRegionIndex][nJetsBin]["Up"]->FindFixBin(eventProgenitorMass, neutralinoMass), totalNEventsErrorUp/(static_cast<float>(totalNEvents_nominal)));

              int totalNEvents_JECDown = inputHistograms->h_totalNEvents_shifted[shiftType::JECDown][STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_totalNEvents_shifted[shiftType::JECDown][STRegionIndex][nJetsBin]->FindFixBin(eventProgenitorMass, neutralinoMass));
              int totalNEvents_JECUp = inputHistograms->h_totalNEvents_shifted[shiftType::JECUp][STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_totalNEvents_shifted[shiftType::JECUp][STRegionIndex][nJetsBin]->FindFixBin(eventProgenitorMass, neutralinoMass));
              fillHistogramsWithAsymmetricErrorsFromNEvents(outputHistograms->h_JECUncertainty[STRegionIndex][nJetsBin], totalNEvents_JECDown, totalNEvents_JECUp, totalNEvents_nominal, eventProgenitorMass, neutralinoMass);

              int totalNEvents_UnclusteredMETDown = inputHistograms->h_totalNEvents_shifted[shiftType::UnclusteredMETDown][STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_totalNEvents_shifted[shiftType::UnclusteredMETDown][STRegionIndex][nJetsBin]->FindFixBin(eventProgenitorMass, neutralinoMass));
              int totalNEvents_UnclusteredMETUp = inputHistograms->h_totalNEvents_shifted[shiftType::UnclusteredMETUp][STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_totalNEvents_shifted[shiftType::UnclusteredMETUp][STRegionIndex][nJetsBin]->FindFixBin(eventProgenitorMass, neutralinoMass));
              fillHistogramsWithAsymmetricErrorsFromNEvents(outputHistograms->h_UnclusteredMETUncertainty[STRegionIndex][nJetsBin], totalNEvents_UnclusteredMETDown, totalNEvents_UnclusteredMETUp, totalNEvents_nominal, eventProgenitorMass, neutralinoMass);

              int totalNEvents_JERMETDown = inputHistograms->h_totalNEvents_shifted[shiftType::JERMETDown][STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_totalNEvents_shifted[shiftType::JERMETDown][STRegionIndex][nJetsBin]->FindFixBin(eventProgenitorMass, neutralinoMass));
              int totalNEvents_JERMETUp = inputHistograms->h_totalNEvents_shifted[shiftType::JERMETUp][STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_totalNEvents_shifted[shiftType::JERMETUp][STRegionIndex][nJetsBin]->FindFixBin(eventProgenitorMass, neutralinoMass));
              fillHistogramsWithAsymmetricErrorsFromNEvents(outputHistograms->h_JERMETUncertainty[STRegionIndex][nJetsBin], totalNEvents_JERMETDown, totalNEvents_JERMETUp, totalNEvents_nominal, eventProgenitorMass, neutralinoMass);

	      int totalNEvents_missingHEMDown = inputHistograms->h_totalNEvents_shifted[shiftType::missingHEMDown][STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_totalNEvents_shifted[shiftType::missingHEMDown][STRegionIndex][nJetsBin]->FindFixBin(eventProgenitorMass, neutralinoMass));
              int totalNEvents_missingHEMUp = inputHistograms->h_totalNEvents_shifted[shiftType::missingHEMUp][STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_totalNEvents_shifted[shiftType::missingHEMUp][STRegionIndex][nJetsBin]->FindFixBin(eventProgenitorMass, neutralinoMass));
	      fillHistogramsWithAsymmetricErrorsFromNEvents(outputHistograms->h_missingHEMUncertainty[STRegionIndex][nJetsBin], totalNEvents_missingHEMDown, totalNEvents_missingHEMUp, totalNEvents_nominal, eventProgenitorMass, neutralinoMass);

              float weightedNEvents_prefiringDown = inputHistograms->h_lumiBasedYearWeightedNEvents_prefiringDown[STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_lumiBasedYearWeightedNEvents_prefiringDown[STRegionIndex][nJetsBin]->FindFixBin(eventProgenitorMass, neutralinoMass));
              float weightedNEvents_prefiringUp = inputHistograms->h_lumiBasedYearWeightedNEvents_prefiringUp[STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_lumiBasedYearWeightedNEvents_prefiringUp[STRegionIndex][nJetsBin]->FindFixBin(eventProgenitorMass, neutralinoMass));
              fillHistogramsWithAsymmetricErrorsFromWeightedNEvents(outputHistograms->h_prefiringWeightsUncertainty[STRegionIndex][nJetsBin], weightedNEvents_prefiringDown, weightedNEvents_prefiringUp, weightedNEvents_nominal, eventProgenitorMass, neutralinoMass);

              float weightedNEvents_HLTDown = inputHistograms->h_lumiBasedYearWeightedNEvents_HLTDown[STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_lumiBasedYearWeightedNEvents_HLTDown[STRegionIndex][nJetsBin]->FindFixBin(eventProgenitorMass, neutralinoMass));
              float weightedNEvents_HLTUp = inputHistograms->h_lumiBasedYearWeightedNEvents_HLTUp[STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_lumiBasedYearWeightedNEvents_HLTUp[STRegionIndex][nJetsBin]->FindFixBin(eventProgenitorMass, neutralinoMass));
              fillHistogramsWithAsymmetricErrorsFromWeightedNEvents(outputHistograms->h_HLTUncertainty[STRegionIndex][nJetsBin], weightedNEvents_HLTDown, weightedNEvents_HLTUp, weightedNEvents_nominal, eventProgenitorMass, neutralinoMass);

              float weightedNEvents_photonScaleFactorDown = inputHistograms->h_lumiBasedYearWeightedNEvents_photonScaleFactorDown[STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_lumiBasedYearWeightedNEvents_photonScaleFactorDown[STRegionIndex][nJetsBin]->FindFixBin(eventProgenitorMass, neutralinoMass));
              float weightedNEvents_photonScaleFactorUp = inputHistograms->h_lumiBasedYearWeightedNEvents_photonScaleFactorUp[STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_lumiBasedYearWeightedNEvents_photonScaleFactorUp[STRegionIndex][nJetsBin]->FindFixBin(eventProgenitorMass, neutralinoMass));
              fillHistogramsWithAsymmetricErrorsFromWeightedNEvents(outputHistograms->h_photonMCScaleFactorUncertainty[STRegionIndex][nJetsBin], weightedNEvents_photonScaleFactorDown, weightedNEvents_photonScaleFactorUp, weightedNEvents_nominal, eventProgenitorMass, neutralinoMass);
            }
          } // end condition that nJets should be >= 4
        } // end loop over neutralino mass
      } // end loop over eventProgenitor mass
    } // end loop over nJetsBin
  } // end loop over STRegionIndex

  if (nProblematicBins > 0) std::cout << "WARNING: total number of recorded events is nonzero but weighted number of recorded events is 0 in a few bins." << std::endl;

  inputFile->Close();
}

void savePlots(outputHistogramsStruct *outputHistograms, optionsStruct &options, const STRegionsStruct& STRegions) {
  std::cout << "Saving output plots..." << std::endl;
  TFile *outputFile = TFile::Open((options.outputDirectory + "/" + options.outputPrefix + "_MCUncertainties_savedObjects.root").c_str(), "RECREATE");
  for (int STRegionIndex = 1; STRegionIndex <= (1+STRegions.nSTSignalBins); ++STRegionIndex) {
    for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
      if ((nJetsBin <= 3 || STRegionIndex == 1) || options.getSignalContaminationOutsideSidebands) {
        std::string histogramName_signalContamination = getHistogramName("signalContamination", STRegionIndex, nJetsBin);
        tmROOTSaverUtils::saveSingleObject(outputHistograms->h_signalContamination[STRegionIndex][nJetsBin], "c_h_" + histogramName_signalContamination, outputFile, options.outputDirectory_signalContamination + "/" + options.outputPrefix + "_" + histogramName_signalContamination + ".pdf", 1024, 768, 0, ".1e", "COLZ TEXT25", false, false, true, 0, 0, 0, 0, 0, 0);
      }
      if(nJetsBin >= 4) {
        std::vector<std::string> UpDownShifts = {"Down", "Up"};
        for (const std::string& UpDownShift: UpDownShifts) {
          std::string histogramName_MCStatisticsFractionalError = getHistogramName("MCStatisticsFractionalError" + UpDownShift, STRegionIndex, nJetsBin);
          tmROOTSaverUtils::saveSingleObject(outputHistograms->h_MCStatisticsFractionalError[STRegionIndex][nJetsBin][UpDownShift], "c_h_" + histogramName_MCStatisticsFractionalError, outputFile, options.outputDirectory + "/" + options.outputPrefix + "_" + histogramName_MCStatisticsFractionalError + ".pdf", 1024, 768, 0, ".0e", "COLZ TEXT25", false, false, true, 0, 0, 0, 0, 0, 0);

          std::string histogramName_JECUncertainty = getHistogramName("JECUncertainty" + UpDownShift, STRegionIndex, nJetsBin);
          tmROOTSaverUtils::saveSingleObject(outputHistograms->h_JECUncertainty[STRegionIndex][nJetsBin][UpDownShift], "c_h_" + histogramName_JECUncertainty, outputFile, options.outputDirectory + "/" + options.outputPrefix + "_" + histogramName_JECUncertainty + ".pdf", 1024, 768, 0, ".0e", "COLZ TEXT25", false, false, true, 0, 0, 0, 0, 0, 0);

          std::string histogramName_UnclusteredMETUncertainty = getHistogramName("UnclusteredMETUncertainty" + UpDownShift, STRegionIndex, nJetsBin);
          tmROOTSaverUtils::saveSingleObject(outputHistograms->h_UnclusteredMETUncertainty[STRegionIndex][nJetsBin][UpDownShift], "c_h_" + histogramName_UnclusteredMETUncertainty, outputFile, options.outputDirectory + "/" + options.outputPrefix + "_" + histogramName_UnclusteredMETUncertainty + ".pdf", 1024, 768, 0, ".0e", "COLZ TEXT25", false, false, true, 0, 0, 0, 0, 0, 0);

          std::string histogramName_JERMETUncertainty = getHistogramName("JERMETUncertainty" + UpDownShift, STRegionIndex, nJetsBin);
          tmROOTSaverUtils::saveSingleObject(outputHistograms->h_JERMETUncertainty[STRegionIndex][nJetsBin][UpDownShift], "c_h_" + histogramName_JERMETUncertainty, outputFile, options.outputDirectory + "/" + options.outputPrefix + "_" + histogramName_JERMETUncertainty + ".pdf", 1024, 768, 0, ".0e", "COLZ TEXT25", false, false, true, 0, 0, 0, 0, 0, 0);

	  std::string histogramName_missingHEMUncertainty = getHistogramName("missingHEMUncertainty" + UpDownShift, STRegionIndex, nJetsBin);
          tmROOTSaverUtils::saveSingleObject(outputHistograms->h_missingHEMUncertainty[STRegionIndex][nJetsBin][UpDownShift], "c_h_" + histogramName_missingHEMUncertainty, outputFile, options.outputDirectory + "/" + options.outputPrefix + "_" + histogramName_missingHEMUncertainty + ".pdf", 1024, 768, 0, ".0e", "COLZ TEXT25", false, false, true, 0, 0, 0, 0, 0, 0);

          std::string histogramName_prefiringWeightsUncertainty = getHistogramName("prefiringWeightsUncertainty" + UpDownShift, STRegionIndex, nJetsBin);
          tmROOTSaverUtils::saveSingleObject(outputHistograms->h_prefiringWeightsUncertainty[STRegionIndex][nJetsBin][UpDownShift], "c_h_" + histogramName_prefiringWeightsUncertainty, outputFile, options.outputDirectory + "/" + options.outputPrefix + "_" + histogramName_prefiringWeightsUncertainty + ".pdf", 1024, 768, 0, ".0e", "COLZ TEXT25", false, false, true, 0, 0, 0, 0, 0, 0);

          std::string histogramName_HLTUncertainty = getHistogramName("HLTUncertainty" + UpDownShift, STRegionIndex, nJetsBin);
          tmROOTSaverUtils::saveSingleObject(outputHistograms->h_HLTUncertainty[STRegionIndex][nJetsBin][UpDownShift], "c_h_" + histogramName_HLTUncertainty, outputFile, options.outputDirectory + "/" + options.outputPrefix + "_" + histogramName_HLTUncertainty + ".pdf", 1024, 768, 0, ".0e", "COLZ TEXT25", false, false, true, 0, 0, 0, 0, 0, 0);

          std::string histogramName_photonMCScaleFactorUncertainty = getHistogramName("photonMCScaleFactorUncertainty" + UpDownShift, STRegionIndex, nJetsBin);
          tmROOTSaverUtils::saveSingleObject(outputHistograms->h_photonMCScaleFactorUncertainty[STRegionIndex][nJetsBin][UpDownShift], "c_h_" + histogramName_photonMCScaleFactorUncertainty, outputFile, options.outputDirectory + "/" + options.outputPrefix + "_" + histogramName_photonMCScaleFactorUncertainty + ".pdf", 1024, 768, 0, ".0e", "COLZ TEXT25", false, false, true, 0, 0, 0, 0, 0, 0);
        }
      }
    }
  }
  outputFile->Close();
}

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  std::cout << "Starting program to get all MC systematics..." << std::endl;

  tmArgumentParser argumentParser = tmArgumentParser("Read in event histograms and calculate various kinds of MC systematics.");
  argumentParser.addArgument("inputPath", "analysis/MCEventHistograms/MC2018_savedObjects.root", true, "Path to ROOT file containing event histograms.");
  argumentParser.addArgument("MCTemplatePath", "", true, "Path to root file that contains a TH2F with bins containing points with generated masses set to 1 and all other bins set to 0.");
  argumentParser.addArgument("inputFile_STRegionBoundaries", "STRegionBoundaries.dat", false, "Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.");
  argumentParser.addArgument("inputNEventsFile", "analysis/dataSystematics/signal_observedEventCounters.dat", false, "Path to file with observations of the nEvents. Used for the signal contamination estimates.");
  argumentParser.addArgument("inputDataUncertaintiesFile", "analysis/dataSystematics/signal_dataSystematics.dat", false, "Path to file with the fractional errors on the background prediction. Used for the signal contamination estimates.");
  // argumentParser.addArgument("inputDataSTScalingUncertaintiesFile", "analysis/dataSystematics/control_dataSystematics_sTScaling.dat", false, "Path to file with the fractional errors on ST scaling. Used for the signal contamination estimates.");
  argumentParser.addArgument("outputDirectory", "analysis/MCSystematics/", false, "Output directory.");
  argumentParser.addArgument("outputDirectory_signalContamination", "analysis/signalContamination/", false, "Output directory for signal contamination plots.");
  argumentParser.addArgument("outputPrefix", "", true, "Prefix to output files.");
  argumentParser.addArgument("getSignalContaminationOutsideSidebands", "false", false, "If set to the string \"true\", then signal contamination is evaluated for all bins in nJets and ST.");
  argumentParser.setPassedStringValues(argc, argv);
  optionsStruct options = getOptionsFromParser(argumentParser);

  MCTemplateReader templateReader(options.MCTemplatePath);
  STRegionsStruct STRegions(options.inputFile_STRegionBoundaries, 20000.0);
  inputNEventsStruct inputNEvents(options.inputNEventsFile);
  // inputDataUncertaintiesStruct inputDataUncertainties(options.inputDataUncertaintiesFile);
  // inputDataSTScalingUncertaintiesStruct inputDataSTScalingUncertainties(options.inputDataSTScalingUncertaintiesFile);
  outputHistogramsStruct* outputHistograms = initializeOutputHistograms(options, templateReader, STRegions);
  fillSystematicsHistograms(outputHistograms, options, templateReader, STRegions, inputNEvents// , inputDataUncertainties, inputDataSTScalingUncertainties
                            );
  savePlots(outputHistograms, options, STRegions);
  return 0;
}
