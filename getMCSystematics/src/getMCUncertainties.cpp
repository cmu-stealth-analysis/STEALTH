#include "../include/getMCUncertainties.h"

optionsStruct getOptionsFromParser(tmArgumentParser& argumentParser) {
  optionsStruct options = optionsStruct();
  options.inputPath = argumentParser.getArgumentString("inputPath");
  options.MCTemplate = argumentParser.getArgumentString("MCTemplate");
  options.inputFile_STRegionBoundaries = argumentParser.getArgumentString("inputFile_STRegionBoundaries");
  options.inputNEventsFile = argumentParser.getArgumentString("inputNEventsFile");
  options.outputDirectory = argumentParser.getArgumentString("outputDirectory");
  options.outputPrefix = argumentParser.getArgumentString("outputPrefix");
  options.nGluinoMassBins = std::stoi(argumentParser.getArgumentString("nGluinoMassBins"));
  options.minGluinoMass = std::stod(argumentParser.getArgumentString("minGluinoMass"));
  options.maxGluinoMass = std::stod(argumentParser.getArgumentString("maxGluinoMass"));
  options.nNeutralinoMassBins = std::stoi(argumentParser.getArgumentString("nNeutralinoMassBins"));
  options.minNeutralinoMass = std::stod(argumentParser.getArgumentString("minNeutralinoMass"));
  options.maxNeutralinoMass = std::stod(argumentParser.getArgumentString("maxNeutralinoMass"));
  return options;
}

std::string getHistogramName(const std::string& histogramType, const int& STRegionIndex, const int& nJetsBin) {
  std::stringstream nameStream;
  nameStream << histogramType << "_" << nJetsBin << "Jets_STRegion" << STRegionIndex;
  return nameStream.str();
}

std::string getHistogramTitle(const std::string& histogramType, const int& STRegionIndex, const int& nJetsBin, const STRegionsStruct& STRegions) {
  std::string histogramTypeString;
  if (histogramType == "JECUncertainty") histogramTypeString = "JEC Uncertainty";
  else if (histogramType == "JECRatioUpToNominal") histogramTypeString = "nEvents passing: (JEC up)/(no JEC shift)";
  else if (histogramType == "JECRatioDownToNominal") histogramTypeString = "nEvents passing: (JEC down)/(no JEC shift)";
  else if (histogramType == "MCStatisticsFractionalError") histogramTypeString = "Fractional error due to MC statistics";
  else if (histogramType == "signalContamination") histogramTypeString = "Signal Contamination";
  else {
    std::cout << "ERROR: Unrecognized histogram type: " << histogramType << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::stringstream nJetsStringStream;
  if (nJetsBin >= 2 && nJetsBin < 6) nJetsStringStream << nJetsBin << " Jets";
  else if (nJetsBin == 6) nJetsStringStream << "#geq 6 Jets";
  else {
    std::cout << "ERROR: Inappropriate nJets bin: " << nJetsBin << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string nJetsString = nJetsStringStream.str();

  std::stringstream sTRangeStringStream;
  if (STRegionIndex >= 1 && STRegionIndex < (1+STRegions.nSTSignalBins)){ // all except the last bin
    sTRangeStringStream << std::fixed << std::setprecision(0) << (STRegions.STAxis).GetBinLowEdge(STRegionIndex) << " < #it{S}_{T} < " << (STRegions.STAxis).GetBinUpEdge(STRegionIndex);
  }
  else if (STRegionIndex == (1+STRegions.nSTSignalBins)) {
    sTRangeStringStream << "#it{S}_{T} > " << std::fixed << std::setprecision(0) << (STRegions.STAxis).GetBinLowEdge(STRegionIndex);
  }
  else {
    std::cout << "ERROR: unexpected region index: " << STRegionIndex << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string sTRangeString = sTRangeStringStream.str();

  std::string axesLabelsString = ";m_{#tilde{#it{g}}};m_{#tilde{#it{#chi_{1}^{0}}}}";

  std::stringstream titleStream;
  titleStream << histogramTypeString << ", " << nJetsString << ", " << sTRangeString << axesLabelsString;
  return titleStream.str();
}

outputHistogramsStruct* initializeOutputHistograms(optionsStruct& options, const STRegionsStruct& STRegions) {
  outputHistogramsStruct* outputHistograms = new outputHistogramsStruct();
  for (int STRegionIndex = 1; STRegionIndex <= (1+STRegions.nSTSignalBins); ++STRegionIndex) {
    for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
      if (nJetsBin <= 3 || STRegionIndex == 1) { // Signal contamination is to be calculated only in the low nJets sideband or at all nJets in the normalization bin
        outputHistograms->h_signalContamination[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("signalContamination", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("signalContamination", STRegionIndex, nJetsBin, STRegions).c_str(), options.nGluinoMassBins, options.minGluinoMass, options.maxGluinoMass, options.nNeutralinoMassBins, options.minNeutralinoMass, options.maxNeutralinoMass);
      }
      if(nJetsBin >= 4) { // the rest of the plots are only useful in the signal bins
        outputHistograms->h_JECUncertainty[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("JECUncertainty", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("JECUncertainty", STRegionIndex, nJetsBin, STRegions).c_str(), options.nGluinoMassBins, options.minGluinoMass, options.maxGluinoMass, options.nNeutralinoMassBins, options.minNeutralinoMass, options.maxNeutralinoMass);
        outputHistograms->h_ratios_JECUpToNominal[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("JECRatioUpToNominal", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("JECRatioUpToNominal", STRegionIndex, nJetsBin, STRegions).c_str(), options.nGluinoMassBins, options.minGluinoMass, options.maxGluinoMass, options.nNeutralinoMassBins, options.minNeutralinoMass, options.maxNeutralinoMass);
        outputHistograms->h_ratios_JECDownToNominal[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("JECRatioDownToNominal", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("JECRatioDownToNominal", STRegionIndex, nJetsBin, STRegions).c_str(), options.nGluinoMassBins, options.minGluinoMass, options.maxGluinoMass, options.nNeutralinoMassBins, options.minNeutralinoMass, options.maxNeutralinoMass);
        outputHistograms->h_MCStatisticsFractionalError[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("MCStatisticsFractionalError", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("MCStatisticsFractionalError", STRegionIndex, nJetsBin, STRegions).c_str(), options.nGluinoMassBins, options.minGluinoMass, options.maxGluinoMass, options.nNeutralinoMassBins, options.minNeutralinoMass, options.maxNeutralinoMass);
      }
    } // end loop over nJetsBin
  } // end loop over STRegionIndex
  return outputHistograms;
}

inputHistogramsStruct* readInputHistograms(TFile *inputFile, const std::vector<std::string>& allowedJECs, const STRegionsStruct& STRegions) {
  inputHistogramsStruct *inputHistograms = new inputHistogramsStruct();
  for (const auto& jec: allowedJECs) {
    for (int STRegionIndex = 1; STRegionIndex <= (1+STRegions.nSTSignalBins); ++STRegionIndex) {
      for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
        std::stringstream nameStreamTotalNEvents;
        nameStreamTotalNEvents << "h_total_nMCEvents_" << jec << "_" << nJetsBin << "Jets_STRegion" << STRegionIndex;
        inputHistograms->h_totalNEvents[jec][STRegionIndex][nJetsBin] = (TH2F*)(inputFile->Get(nameStreamTotalNEvents.str().c_str()));
        if (inputHistograms->h_totalNEvents[jec][STRegionIndex][nJetsBin] == nullptr) {
          std::cout << "Unable to find histogram with name " << nameStreamTotalNEvents.str() << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if(inputHistograms->h_totalNEvents[jec][STRegionIndex][nJetsBin]->GetBinErrorOption() != TH1::EBinErrorOpt::kPoisson) {
          std::cout << "ERROR: errors on histogram of total number of events not Poisson." << std::endl;
          std::exit(EXIT_FAILURE);
        }

        std::stringstream nameStreamWeightedNEvents;
        nameStreamWeightedNEvents << "h_weighted_nMCEvents_" << jec << "_" << nJetsBin << "Jets_STRegion" << STRegionIndex;
        inputHistograms->h_weightedNEvents[jec][STRegionIndex][nJetsBin] = (TH2F*)(inputFile->Get(nameStreamWeightedNEvents.str().c_str()));
        if (inputHistograms->h_weightedNEvents[jec][STRegionIndex][nJetsBin] == nullptr) {
          std::cout << "Unable to find histogram with name " << nameStreamWeightedNEvents.str() << std::endl;
          std::exit(EXIT_FAILURE);
        }
      } // ends loop over NJetsBin
    } // ends loop over ST regions
  } // ends loop over JECs
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

void fillSystematicsHistograms(outputHistogramsStruct *outputHistograms, optionsStruct& options, const std::vector<std::string>& allowedJECs, const STRegionsStruct& STRegions, inputNEventsStruct& inputNEvents) {
  TFile *inputFile = TFile::Open(options.inputPath.c_str(), "READ");
  if (inputFile->IsZombie() || !(inputFile->IsOpen())) {
    std::cout << "Error in opening file " << options.inputPath << std::endl;
    std::exit(EXIT_FAILURE);
  }
  inputHistogramsStruct *inputHistograms = readInputHistograms(inputFile, allowedJECs, STRegions);

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
          double weightedNEvents_nominal = inputHistograms->h_weightedNEvents["JECNominal"][STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_weightedNEvents["JECNominal"][STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass));
          if (nJetsBin <= 3 || STRegionIndex == 1) {
            std::stringstream inputNEventsStringStream;
            inputNEventsStringStream << "observedNEvents_STRegion" << STRegionIndex << "_" << nJetsBin << "Jets";
            outputHistograms->h_signalContamination[STRegionIndex][nJetsBin]->SetBinContent(outputHistograms->h_signalContamination[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass), weightedNEvents_nominal/((inputNEvents.data)[inputNEventsStringStream.str()]));
          }
          if(nJetsBin >= 4) {
            bool zeroMCEventsRecorded = false;
            double weightedNEvents_jecUp = inputHistograms->h_weightedNEvents["JECUp"][STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_weightedNEvents["JECUp"][STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass));
            double weightedNEvents_jecDown = inputHistograms->h_weightedNEvents["JECDown"][STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_weightedNEvents["JECDown"][STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass));
            double totalNEvents_nominal = inputHistograms->h_totalNEvents["JECNominal"][STRegionIndex][nJetsBin]->GetBinContent(inputHistograms->h_totalNEvents["JECNominal"][STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass));
            if (totalNEvents_nominal == 0) {
              std::cout << "WARNING: zero events recorded at gluino mass: " << gluinoMass << ", neutralino mass: " << neutralinoMass << ", for STRegionIndex: " << STRegionIndex << ", nJets: " << nJetsBin << std::endl;
              zeroMCEventsRecorded = true;
            }
            if ((weightedNEvents_nominal == 0) && !(zeroMCEventsRecorded)) {
              std::cout << "ERROR: total number of recorded events is nonzero but weighted number of recorded events is 0 at gluino mass: " << gluinoMass << ", neutralino mass: " << neutralinoMass << ", for STRegionIndex: " << STRegionIndex << ", nJets: " << nJetsBin << std::endl;
              std::exit(EXIT_FAILURE);
            }
            double totalNEventsError_nominal = 0.5*(inputHistograms->h_totalNEvents["JECNominal"][STRegionIndex][nJetsBin]->GetBinErrorUp(inputHistograms->h_totalNEvents["JECNominal"][STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass)) + inputHistograms->h_totalNEvents["JECNominal"][STRegionIndex][nJetsBin]->GetBinErrorLow(inputHistograms->h_totalNEvents["JECNominal"][STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass)));

            double fractionalMCStatisticsUncertainty = 0.;
            if (zeroMCEventsRecorded) fractionalMCStatisticsUncertainty = totalNEventsError_nominal; // fractional error on 0 is ill-defined...
            else fractionalMCStatisticsUncertainty = totalNEventsError_nominal/totalNEvents_nominal;
            outputHistograms->h_MCStatisticsFractionalError[STRegionIndex][nJetsBin]->SetBinContent(outputHistograms->h_MCStatisticsFractionalError[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass), fractionalMCStatisticsUncertainty);
            if (zeroMCEventsRecorded) continue;

            double ratioUpToNominal = weightedNEvents_jecUp / weightedNEvents_nominal;
            double deviationUp = std::fabs(ratioUpToNominal - 1.0);
            double ratioDownToNominal = weightedNEvents_jecDown / weightedNEvents_nominal;
            double deviationDown = std::fabs(ratioDownToNominal - 1.0);
            double averageDeviation = 0.5*(deviationUp + deviationDown);
            outputHistograms->h_JECUncertainty[STRegionIndex][nJetsBin]->SetBinContent(outputHistograms->h_JECUncertainty[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass), averageDeviation);
            outputHistograms->h_ratios_JECUpToNominal[STRegionIndex][nJetsBin]->SetBinContent(outputHistograms->h_ratios_JECUpToNominal[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass), ratioUpToNominal);
            outputHistograms->h_ratios_JECDownToNominal[STRegionIndex][nJetsBin]->SetBinContent(outputHistograms->h_ratios_JECDownToNominal[STRegionIndex][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass), ratioDownToNominal);
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
      if (nJetsBin <= 3 || STRegionIndex == 1) {
        std::string histogramName_signalContamination = getHistogramName("signalContamination", STRegionIndex, nJetsBin);
        tmROOTSaverUtils::saveSingleObject(outputHistograms->h_signalContamination[STRegionIndex][nJetsBin], "c_h_" + histogramName_signalContamination, outputFile, "analysis/signalContamination/" + options.outputPrefix + "_" + histogramName_signalContamination + "_hist.png", 1024, 768, 0, ".0e", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);
      }
      if(nJetsBin >= 4) {
        std::string histogramName_JECUncertainty = getHistogramName("JECUncertainty", STRegionIndex, nJetsBin);
        tmROOTSaverUtils::saveSingleObject(outputHistograms->h_JECUncertainty[STRegionIndex][nJetsBin], "c_h_" + histogramName_JECUncertainty, outputFile, options.outputDirectory + "/" + options.outputPrefix + "_" + histogramName_JECUncertainty + "_hist.png", 1024, 768, 0, ".0e", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);
        std::string histogramName_ratios_JECUpToNominal = getHistogramName("ratios_JECUpToNominal", STRegionIndex, nJetsBin);
        outputHistograms->h_ratios_JECUpToNominal[STRegionIndex][nJetsBin]->SetMarkerSize(0.7);
        tmROOTSaverUtils::saveSingleObject(outputHistograms->h_ratios_JECUpToNominal[STRegionIndex][nJetsBin], "c_h_" + histogramName_ratios_JECUpToNominal, outputFile, options.outputDirectory + "/" + options.outputPrefix + "_" + histogramName_ratios_JECUpToNominal + "_hist.png", 1024, 768, 0, ".1e", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0.97, 1.03);
        std::string histogramName_ratios_JECDownToNominal = getHistogramName("ratios_JECDownToNominal", STRegionIndex, nJetsBin);
        outputHistograms->h_ratios_JECDownToNominal[STRegionIndex][nJetsBin]->SetMarkerSize(0.7);
        tmROOTSaverUtils::saveSingleObject(outputHistograms->h_ratios_JECDownToNominal[STRegionIndex][nJetsBin], "c_h_" + histogramName_ratios_JECDownToNominal, outputFile, options.outputDirectory + "/" + options.outputPrefix + "_" + histogramName_ratios_JECDownToNominal + "_hist.png", 1024, 768, 0, ".1e", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0.97, 1.03);
        std::string histogramName_MCStatisticsFractionalError = getHistogramName("MCStatisticsFractionalError", STRegionIndex, nJetsBin);
        tmROOTSaverUtils::saveSingleObject(outputHistograms->h_MCStatisticsFractionalError[STRegionIndex][nJetsBin], "c_h_" + histogramName_MCStatisticsFractionalError, outputFile, options.outputDirectory + "/" + options.outputPrefix + "_" + histogramName_MCStatisticsFractionalError + "_hist.png", 1024, 768, 0, ".0e", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);
      }
    }
  }
  outputFile->Close();
}

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  std::cout << "Starting program to get systematics due to uncertainty on JECs..." << std::endl;

  std::vector<std::string> allowedJECs{"JECDown", "JECNominal", "JECUp"};

  tmArgumentParser argumentParser = tmArgumentParser("Read in event histograms and calculate systematics due to uncertainty on jet energy corrections.");
  argumentParser.addArgument("inputPath", "analysis/MCEventHistograms/MC2018_savedObjects.root", true, "Path to ROOT file containing event histograms.");
  argumentParser.addArgument("MCTemplate", "plot_susyMasses_template.root", false, "Path to root file that contains a TH2F with bins containing points with generated masses set to 1 and all other bins set to 0.");
  argumentParser.addArgument("inputFile_STRegionBoundaries", "STRegionBoundaries.dat", false, "Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.");
  argumentParser.addArgument("inputNEventsFile", "analysis/dataSystematics/signal_observedEventCounters.dat", false, "Path to file with observations of the nEvents. Used for the signal contamination estimates.");
  argumentParser.addArgument("outputDirectory", "analysis/MCSystematics/", false, "Prefix to output files.");
  argumentParser.addArgument("outputPrefix", "", true, "Prefix to output files.");
  argumentParser.addArgument("nGluinoMassBins", "20", false, "nBins on the gluino mass axis."); // (800 - 25) GeV --> (1750 + 25) GeV in steps of 50 GeV
  argumentParser.addArgument("minGluinoMass", "775.0", false, "Min gluino mass for the 2D plots.");
  argumentParser.addArgument("maxGluinoMass", "1775.0", false, "Max gluino mass for the 2D plots.");
  argumentParser.addArgument("nNeutralinoMassBins", "133", false, "nBins on the neutralino mass axis.");
  argumentParser.addArgument("minNeutralinoMass", "93.75", false, "Min neutralino mass for the 2D plots.");
  argumentParser.addArgument("maxNeutralinoMass", "1756.25", false, "Max neutralino mass for the 2D plots."); // (100 - 6.25) GeV --> (1750 + 6.25) GeV in steps of 12.5 GeV
  argumentParser.setPassedStringValues(argc, argv);
  optionsStruct options = getOptionsFromParser(argumentParser);

  STRegionsStruct STRegions(options.inputFile_STRegionBoundaries);
  inputNEventsStruct inputNEvents(options.inputNEventsFile);
  outputHistogramsStruct* outputHistograms = initializeOutputHistograms(options, STRegions);
  fillSystematicsHistograms(outputHistograms, options, allowedJECs, STRegions, inputNEvents);
  savePlots(outputHistograms, options, STRegions);
  return 0;
}
