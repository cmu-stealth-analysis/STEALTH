#include "../include/getMCUncertainties.h"

outputHistogramsStruct* initializeOutputHistograms(optionsStruct& options, const STRegionsStruct& STRegions) {
  outputHistogramsStruct* outputHistograms = new outputHistogramsStruct();
  for (int STRegionIndex = 1; STRegionIndex <= (1+STRegions.nSTSignalBins); ++STRegionIndex) {
    for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
      if ((nJetsBin <= 3 || STRegionIndex == 1) || options.unrestrictedSignalContamination) { // Signal contamination is to be calculated only in the low nJets sideband or at all nJets in the normalization bin
        outputHistograms->h_signalContamination[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("signalContamination", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("signalContamination", STRegionIndex, nJetsBin, STRegions).c_str(), options.nGluinoMassBins, options.minGluinoMass, options.maxGluinoMass, options.nNeutralinoMassBins, options.minNeutralinoMass, options.maxNeutralinoMass);
      }
      if(nJetsBin >= 4) { // the rest of the plots are only useful in the signal bins
        outputHistograms->h_MCStatisticsFractionalError[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("MCStatisticsFractionalError", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("MCStatisticsFractionalError", STRegionIndex, nJetsBin, STRegions).c_str(), options.nGluinoMassBins, options.minGluinoMass, options.maxGluinoMass, options.nNeutralinoMassBins, options.minNeutralinoMass, options.maxNeutralinoMass);
        outputHistograms->h_JECUncertainty[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("JECUncertainty", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("JECUncertainty", STRegionIndex, nJetsBin, STRegions).c_str(), options.nGluinoMassBins, options.minGluinoMass, options.maxGluinoMass, options.nNeutralinoMassBins, options.minNeutralinoMass, options.maxNeutralinoMass);
        outputHistograms->h_UnclusteredMETUncertainty[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("UnclusteredMETUncertainty", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("UnclusteredMETUncertainty", STRegionIndex, nJetsBin, STRegions).c_str(), options.nGluinoMassBins, options.minGluinoMass, options.maxGluinoMass, options.nNeutralinoMassBins, options.minNeutralinoMass, options.maxNeutralinoMass);
        outputHistograms->h_JERMETUncertainty[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("JERMETUncertainty", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("JERMETUncertainty", STRegionIndex, nJetsBin, STRegions).c_str(), options.nGluinoMassBins, options.minGluinoMass, options.maxGluinoMass, options.nNeutralinoMassBins, options.minNeutralinoMass, options.maxNeutralinoMass);
        outputHistograms->h_prefiringWeightsUncertainty[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("prefiringWeightsUncertainty", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("prefiringWeightsUncertainty", STRegionIndex, nJetsBin, STRegions).c_str(), options.nGluinoMassBins, options.minGluinoMass, options.maxGluinoMass, options.nNeutralinoMassBins, options.minNeutralinoMass, options.maxNeutralinoMass);
        outputHistograms->h_photonMCScaleFactorUncertainty[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("photonMCScaleFactorUncertainty", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("photonMCScaleFactorUncertainty", STRegionIndex, nJetsBin, STRegions).c_str(), options.nGluinoMassBins, options.minGluinoMass, options.maxGluinoMass, options.nNeutralinoMassBins, options.minNeutralinoMass, options.maxNeutralinoMass);
      }
    } // end loop over nJetsBin
  } // end loop over STRegionIndex
  return outputHistograms;
}

inputHistogramsStruct* readInputHistograms(TFile *inputFile, const STRegionsStruct& STRegions) {
  inputHistogramsStruct *inputHistograms = new inputHistogramsStruct();
  for (int STRegionIndex = 1; STRegionIndex <= (1+STRegions.nSTSignalBins); ++STRegionIndex) {
    for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
      std::stringstream nameStreamTotalNEvents;
      nameStreamTotalNEvents << "h_totalNEvents_" << nJetsBin << "Jets_STRegion" << STRegionIndex;
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
        nameStreamTotalNEventsShifted << "h_totalNEvents_shifted_" << shiftTypeNames[typeIndex] << "_regionIndex_" << STRegionIndex << "_" << nJetsBin << "Jets";
        inputHistograms->h_totalNEvents_shifted[typeIndex][STRegionIndex][nJetsBin] = (TH2I*)(inputFile->Get(nameStreamTotalNEventsShifted.str().c_str()));
        if (inputHistograms->h_totalNEvents_shifted[typeIndex][STRegionIndex][nJetsBin] == nullptr) {
          std::cout << "Unable to find histogram with name " << nameStreamTotalNEventsShifted.str() << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }

      std::string commonPrefix = "h_lumiBasedYearWeightedNEvents";
      std::stringstream commonSuffixStringStream;
      commonSuffixStringStream << "_" << nJetsBin << "Jets_STRegion" << STRegionIndex;
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

TH2F* readMCTemplate(TFile *inputFile) {
  TH2F *MCTemplateTH2 = (TH2F*)(inputFile->Get("h_susyMasses_template"));
  if (MCTemplateTH2 == nullptr) {
    std::cout << "ERROR: MC template not found" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return MCTemplateTH2;
}

float getError(float nominal, float variation1, float variation2) {
  float fractionalError1 = std::fabs(variation1/nominal - 1.0f);
  float fractionalError2 = std::fabs(variation2/nominal - 1.0f);
  float averageFractionalError = 0.5*(fractionalError1 + fractionalError2);
  return averageFractionalError;
}

float getErrorInt(int nominal, int variation1, int variation2) {
  return getError(static_cast<float>(nominal), static_cast<float>(variation1), static_cast<float>(variation2));
}

void fillSystematicsHistograms(outputHistogramsStruct *outputHistograms, optionsStruct& options, const STRegionsStruct& STRegions, inputNEventsStruct& inputNEvents) {
  TFile *inputFile = TFile::Open(options.inputPath.c_str(), "READ");
  if (inputFile->IsZombie() || !(inputFile->IsOpen())) {
    std::cout << "Error in opening file " << options.inputPath << std::endl;
    std::exit(EXIT_FAILURE);
  }
  inputHistogramsStruct *inputHistograms = readInputHistograms(inputFile, STRegions);

  TFile *MCTemplateFile = TFile::Open(options.MCTemplate.c_str(), "READ");
  if (MCTemplateFile->IsZombie() || !(MCTemplateFile->IsOpen())) {
    std::cout << "Error in opening file " << options.MCTemplate << std::endl;
    std::exit(EXIT_FAILURE);
  }
  TH2F *MCTemplateTH2 = readMCTemplate(MCTemplateFile);

  std::cout << "Getting systematics..." << std::endl;

  // Fill TGraphs and TH2s with the JEC fractional uncertainty and estimated error
  for (int STRegionIndex = 1; STRegionIndex <= (1+STRegions.nSTSignalBins); ++STRegionIndex) {
    for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
      // Loop over only those bins that show a "1" in the MC template file
      for (int templateGluinoMassIndex = 1; templateGluinoMassIndex <= MCTemplateTH2->GetXaxis()->GetNbins(); ++templateGluinoMassIndex) {
        double gluinoMass = MCTemplateTH2->GetXaxis()->GetBinCenter(templateGluinoMassIndex);
        for (int templateNeutralinoMassIndex = 1; templateNeutralinoMassIndex <= MCTemplateTH2->GetYaxis()->GetNbins(); ++templateNeutralinoMassIndex) {
          if (!(1 == static_cast<int>(0.5 + MCTemplateTH2->GetBinContent(templateGluinoMassIndex, templateNeutralinoMassIndex)))) continue;
          double neutralinoMass = MCTemplateTH2->GetYaxis()->GetBinCenter(templateNeutralinoMassIndex);
          float weightedNEvents_nominal = inputHistograms->h_lumiBasedYearWeightedNEvents[STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_lumiBasedYearWeightedNEvents[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass));
          if ((nJetsBin <= 3 || STRegionIndex == 1) || options.unrestrictedSignalContamination) {
            std::stringstream inputNEventsStringStream;
            inputNEventsStringStream << "observedNEvents_STRegion" << STRegionIndex << "_" << nJetsBin << "Jets";
            int nBackgroundEvts = ((inputNEvents.data)[inputNEventsStringStream.str()]);
            if (nBackgroundEvts > 0) {
              outputHistograms->h_signalContamination[STRegionIndex][nJetsBin]->SetBinContent(outputHistograms->h_signalContamination[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass), weightedNEvents_nominal/nBackgroundEvts);
            }
            else {
              outputHistograms->h_signalContamination[STRegionIndex][nJetsBin]->SetBinContent(outputHistograms->h_signalContamination[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass), weightedNEvents_nominal);
            }
          }
          if(nJetsBin >= 4) {
            bool zeroMCEventsRecorded = false;
            int totalNEvents_nominal = inputHistograms->h_totalNEvents[STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_totalNEvents[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass));
            double totalNEventsError_nominal = 0.5*((inputHistograms->h_totalNEvents[STRegionIndex][nJetsBin]->GetBinErrorUp(inputHistograms->h_totalNEvents[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass))) + (inputHistograms->h_totalNEvents[STRegionIndex][nJetsBin]->GetBinErrorLow(inputHistograms->h_totalNEvents[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass))));

            if (totalNEvents_nominal == 0) {
              std::cout << "WARNING: zero events recorded at gluino mass: " << gluinoMass << ", neutralino mass: " << neutralinoMass << ", for STRegionIndex: " << STRegionIndex << ", nJets: " << nJetsBin << std::endl;
              zeroMCEventsRecorded = true;
            }

            if ((weightedNEvents_nominal == 0) && !(zeroMCEventsRecorded)) { // sanity check
              std::cout << "ERROR: total number of recorded events is nonzero but weighted number of recorded events is 0 at gluino mass: " << gluinoMass << ", neutralino mass: " << neutralinoMass << ", for STRegionIndex: " << STRegionIndex << ", nJets: " << nJetsBin << std::endl;
              std::exit(EXIT_FAILURE);
            }

            if (zeroMCEventsRecorded) {
              outputHistograms->h_MCStatisticsFractionalError[STRegionIndex][nJetsBin]->SetBinContent(outputHistograms->h_MCStatisticsFractionalError[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass), DEFAULT_FRACTIONAL_ERROR);
              outputHistograms->h_JECUncertainty[STRegionIndex][nJetsBin]->SetBinContent(outputHistograms->h_JECUncertainty[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass), DEFAULT_FRACTIONAL_ERROR);
              outputHistograms->h_UnclusteredMETUncertainty[STRegionIndex][nJetsBin]->SetBinContent(outputHistograms->h_UnclusteredMETUncertainty[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass), DEFAULT_FRACTIONAL_ERROR);
              outputHistograms->h_JERMETUncertainty[STRegionIndex][nJetsBin]->SetBinContent(outputHistograms->h_JERMETUncertainty[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass), DEFAULT_FRACTIONAL_ERROR);
              outputHistograms->h_prefiringWeightsUncertainty[STRegionIndex][nJetsBin]->SetBinContent(outputHistograms->h_prefiringWeightsUncertainty[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass), DEFAULT_FRACTIONAL_ERROR);
              outputHistograms->h_photonMCScaleFactorUncertainty[STRegionIndex][nJetsBin]->SetBinContent(outputHistograms->h_photonMCScaleFactorUncertainty[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass), DEFAULT_FRACTIONAL_ERROR);
            }
            else {
              float fractionalMCStatisticsUncertainty = totalNEventsError_nominal/(static_cast<float>(totalNEvents_nominal));
              outputHistograms->h_MCStatisticsFractionalError[STRegionIndex][nJetsBin]->SetBinContent(outputHistograms->h_MCStatisticsFractionalError[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass), fractionalMCStatisticsUncertainty);
            
              int totalNEvents_JECDown = inputHistograms->h_totalNEvents_shifted[shiftType::JECDown][STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_totalNEvents_shifted[shiftType::JECDown][STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass));
              int totalNEvents_JECUp = inputHistograms->h_totalNEvents_shifted[shiftType::JECUp][STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_totalNEvents_shifted[shiftType::JECUp][STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass));
              outputHistograms->h_JECUncertainty[STRegionIndex][nJetsBin]->SetBinContent(outputHistograms->h_JECUncertainty[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass), getErrorInt(totalNEvents_nominal, totalNEvents_JECDown, totalNEvents_JECUp));

              int totalNEvents_UnclusteredMETDown = inputHistograms->h_totalNEvents_shifted[shiftType::UnclusteredMETDown][STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_totalNEvents_shifted[shiftType::UnclusteredMETDown][STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass));
              int totalNEvents_UnclusteredMETUp = inputHistograms->h_totalNEvents_shifted[shiftType::UnclusteredMETUp][STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_totalNEvents_shifted[shiftType::UnclusteredMETUp][STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass));
              outputHistograms->h_UnclusteredMETUncertainty[STRegionIndex][nJetsBin]->SetBinContent(outputHistograms->h_UnclusteredMETUncertainty[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass), getErrorInt(totalNEvents_nominal, totalNEvents_UnclusteredMETDown, totalNEvents_UnclusteredMETUp));

              int totalNEvents_JERMETDown = inputHistograms->h_totalNEvents_shifted[shiftType::JERMETDown][STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_totalNEvents_shifted[shiftType::JERMETDown][STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass));
              int totalNEvents_JERMETUp = inputHistograms->h_totalNEvents_shifted[shiftType::JERMETUp][STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_totalNEvents_shifted[shiftType::JERMETUp][STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass));
              outputHistograms->h_JERMETUncertainty[STRegionIndex][nJetsBin]->SetBinContent(outputHistograms->h_JERMETUncertainty[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass), getErrorInt(totalNEvents_nominal, totalNEvents_JERMETDown, totalNEvents_JERMETUp));
            
              float weightedNEvents_prefiringDown = inputHistograms->h_lumiBasedYearWeightedNEvents_prefiringDown[STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_lumiBasedYearWeightedNEvents_prefiringDown[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass));
              float weightedNEvents_prefiringUp = inputHistograms->h_lumiBasedYearWeightedNEvents_prefiringUp[STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_lumiBasedYearWeightedNEvents_prefiringUp[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass));
              outputHistograms->h_prefiringWeightsUncertainty[STRegionIndex][nJetsBin]->SetBinContent(outputHistograms->h_prefiringWeightsUncertainty[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass), getError(weightedNEvents_nominal, weightedNEvents_prefiringDown, weightedNEvents_prefiringUp));

              float weightedNEvents_photonScaleFactorDown = inputHistograms->h_lumiBasedYearWeightedNEvents_photonScaleFactorDown[STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_lumiBasedYearWeightedNEvents_photonScaleFactorDown[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass));
              float weightedNEvents_photonScaleFactorUp = inputHistograms->h_lumiBasedYearWeightedNEvents_photonScaleFactorUp[STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_lumiBasedYearWeightedNEvents_photonScaleFactorUp[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass));
              outputHistograms->h_photonMCScaleFactorUncertainty[STRegionIndex][nJetsBin]->SetBinContent(outputHistograms->h_photonMCScaleFactorUncertainty[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass), getError(weightedNEvents_nominal, weightedNEvents_photonScaleFactorDown, weightedNEvents_photonScaleFactorUp));
            }
          } // end condition that nJets should be >= 4
        } // end loop over neutralino mass
      } // end loop over gluino mass
    } // end loop over nJetsBin
  } // end loop over STRegionIndex

  MCTemplateFile->Close();
  inputFile->Close();
}

void savePlots(outputHistogramsStruct *outputHistograms, optionsStruct &options, const STRegionsStruct& STRegions) {
  std::cout << "Saving output plots..." << std::endl;
  TFile *outputFile = TFile::Open((options.outputDirectory + "/" + options.outputPrefix + "_MCUncertainties_savedObjects.root").c_str(), "RECREATE");
  for (int STRegionIndex = 1; STRegionIndex <= (1+STRegions.nSTSignalBins); ++STRegionIndex) {
    for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
      if ((nJetsBin <= 3 || STRegionIndex == 1) || options.unrestrictedSignalContamination) {
        std::string histogramName_signalContamination = getHistogramName("signalContamination", STRegionIndex, nJetsBin);
        tmROOTSaverUtils::saveSingleObject(outputHistograms->h_signalContamination[STRegionIndex][nJetsBin], "c_h_" + histogramName_signalContamination, outputFile, options.outputDirectory_signalContamination + "/" + options.outputPrefix + "_" + histogramName_signalContamination + ".png", 1024, 768, 0, ".0e", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);
      }
      if(nJetsBin >= 4) {
        std::string histogramName_MCStatisticsFractionalError = getHistogramName("MCStatisticsFractionalError", STRegionIndex, nJetsBin);
        tmROOTSaverUtils::saveSingleObject(outputHistograms->h_MCStatisticsFractionalError[STRegionIndex][nJetsBin], "c_h_" + histogramName_MCStatisticsFractionalError, outputFile, options.outputDirectory + "/" + options.outputPrefix + "_" + histogramName_MCStatisticsFractionalError + ".png", 1024, 768, 0, ".0e", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);

        std::string histogramName_JECUncertainty = getHistogramName("JECUncertainty", STRegionIndex, nJetsBin);
        tmROOTSaverUtils::saveSingleObject(outputHistograms->h_JECUncertainty[STRegionIndex][nJetsBin], "c_h_" + histogramName_JECUncertainty, outputFile, options.outputDirectory + "/" + options.outputPrefix + "_" + histogramName_JECUncertainty + ".png", 1024, 768, 0, ".0e", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);

        std::string histogramName_UnclusteredMETUncertainty = getHistogramName("UnclusteredMETUncertainty", STRegionIndex, nJetsBin);
        tmROOTSaverUtils::saveSingleObject(outputHistograms->h_UnclusteredMETUncertainty[STRegionIndex][nJetsBin], "c_h_" + histogramName_UnclusteredMETUncertainty, outputFile, options.outputDirectory + "/" + options.outputPrefix + "_" + histogramName_UnclusteredMETUncertainty + ".png", 1024, 768, 0, ".0e", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);

        std::string histogramName_JERMETUncertainty = getHistogramName("JERMETUncertainty", STRegionIndex, nJetsBin);
        tmROOTSaverUtils::saveSingleObject(outputHistograms->h_JERMETUncertainty[STRegionIndex][nJetsBin], "c_h_" + histogramName_JERMETUncertainty, outputFile, options.outputDirectory + "/" + options.outputPrefix + "_" + histogramName_JERMETUncertainty + ".png", 1024, 768, 0, ".0e", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);

        std::string histogramName_prefiringWeightsUncertainty = getHistogramName("prefiringWeightsUncertainty", STRegionIndex, nJetsBin);
        tmROOTSaverUtils::saveSingleObject(outputHistograms->h_prefiringWeightsUncertainty[STRegionIndex][nJetsBin], "c_h_" + histogramName_prefiringWeightsUncertainty, outputFile, options.outputDirectory + "/" + options.outputPrefix + "_" + histogramName_prefiringWeightsUncertainty + ".png", 1024, 768, 0, ".0e", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);

        std::string histogramName_photonMCScaleFactorUncertainty = getHistogramName("photonMCScaleFactorUncertainty", STRegionIndex, nJetsBin);
        tmROOTSaverUtils::saveSingleObject(outputHistograms->h_photonMCScaleFactorUncertainty[STRegionIndex][nJetsBin], "c_h_" + histogramName_photonMCScaleFactorUncertainty, outputFile, options.outputDirectory + "/" + options.outputPrefix + "_" + histogramName_photonMCScaleFactorUncertainty + ".png", 1024, 768, 0, ".0e", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);
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
  argumentParser.addArgument("MCTemplate", "plot_susyMasses_template.root", false, "Path to root file that contains a TH2F with bins containing points with generated masses set to 1 and all other bins set to 0.");
  argumentParser.addArgument("inputFile_STRegionBoundaries", "STRegionBoundaries.dat", false, "Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.");
  argumentParser.addArgument("inputNEventsFile", "analysis/dataSystematics/signal_observedEventCounters.dat", false, "Path to file with observations of the nEvents. Used for the signal contamination estimates.");
  argumentParser.addArgument("outputDirectory", "analysis/MCSystematics/", false, "Output directory.");
  argumentParser.addArgument("outputDirectory_signalContamination", "analysis/signalContamination/", false, "Output directory for signal contamination plots.");
  argumentParser.addArgument("outputPrefix", "", true, "Prefix to output files.");
  argumentParser.addArgument("nGluinoMassBins", "20", false, "nBins on the gluino mass axis."); // (800 - 25) GeV --> (1750 + 25) GeV in steps of 50 GeV
  argumentParser.addArgument("minGluinoMass", "775.0", false, "Min gluino mass for the 2D plots.");
  argumentParser.addArgument("maxGluinoMass", "1775.0", false, "Max gluino mass for the 2D plots.");
  argumentParser.addArgument("nNeutralinoMassBins", "133", false, "nBins on the neutralino mass axis.");
  argumentParser.addArgument("minNeutralinoMass", "93.75", false, "Min neutralino mass for the 2D plots.");
  argumentParser.addArgument("maxNeutralinoMass", "1756.25", false, "Max neutralino mass for the 2D plots."); // (100 - 6.25) GeV --> (1750 + 6.25) GeV in steps of 12.5 GeV
  argumentParser.addArgument("unrestrictedSignalContamination", "false", false, "If set to the string \"true\", then signal contamination is evaluated for all bins in nJets and ST.");
  argumentParser.setPassedStringValues(argc, argv);
  optionsStruct options = getOptionsFromParser(argumentParser);

  STRegionsStruct STRegions(options.inputFile_STRegionBoundaries);
  inputNEventsStruct inputNEvents(options.inputNEventsFile);
  outputHistogramsStruct* outputHistograms = initializeOutputHistograms(options, STRegions);
  fillSystematicsHistograms(outputHistograms, options, STRegions, inputNEvents);
  savePlots(outputHistograms, options, STRegions);
  return 0;
}
