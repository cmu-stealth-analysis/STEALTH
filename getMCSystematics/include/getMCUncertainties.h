#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <map>
#include <algorithm>
#include <cassert>

#include "tmArgumentParser.h"
#include "tmROOTSaverUtils.h"
#include "tmMiscellaneous.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH2I.h"
#include "TH2F.h"
#include "TGraph2D.h"

#include "../../eventSelection/include/STRegionsStruct.h"
#include "../../eventSelection/include/shiftedObservablesStruct.h"
#include "../../eventSelection/include/MCTemplateReader.h"

#define DEFAULT_FRACTIONAL_ERROR 0.0004

struct optionsStruct{
  std::string inputPath, MCTemplatePath, inputFile_STRegionBoundaries, inputNEventsFile, inputDataUncertaintiesFile, inputDataSTScalingUncertaintiesFile, outputDirectory, outputDirectory_signalContamination, outputPrefix;
  bool getSignalContaminationOutsideSidebands;
};

struct outputHistogramsStruct{
  // syntax for signal contamination: outputHistogram[STRegionIndex][nJetsBin] where regionIndex ranges from 1 to (1 + number of ST signal bins), where regionIndex 1 corresponds to the normalization bin
  // syntax for all others: same as for signal contamination but with an additional index that can be either "Up" or "Down" (for asymmetric systematics)
  std::map<int, std::map<int, TH2F* > > h_signalContamination;
  std::map<int, std::map<int, std::map<std::string, TH2F* > > > h_MCStatisticsFractionalError;
  std::map<int, std::map<int, std::map<std::string, TH2F* > > > h_JECUncertainty;
  std::map<int, std::map<int, std::map<std::string, TH2F* > > > h_UnclusteredMETUncertainty;
  std::map<int, std::map<int, std::map<std::string, TH2F* > > > h_JERMETUncertainty;
  std::map<int, std::map<int, std::map<std::string, TH2F* > > > h_missingHEMUncertainty;
  std::map<int, std::map<int, std::map<std::string, TH2F* > > > h_prefiringWeightsUncertainty;
  std::map<int, std::map<int, std::map<std::string, TH2F* > > > h_photonMCScaleFactorUncertainty;
};

struct inputHistogramsStruct{
  // histograms for nominal number of events
  // syntax: histograms[regionIndex][nJetsBin] where regionIndex ranges from 1 to (1 + number of ST signal bins), where regionIndex 1 corresponds to the normalization bin
  std::map<int, std::map<int, TH2I* > > h_totalNEvents;
  std::map<int, std::map<int, TH2F* > > h_lumiBasedYearWeightedNEvents; // nominal
  std::map<int, std::map<int, TH2F* > > h_lumiBasedYearWeightedNEvents_prefiringDown;
  std::map<int, std::map<int, TH2F* > > h_lumiBasedYearWeightedNEvents_prefiringUp;
  std::map<int, std::map<int, TH2F* > > h_lumiBasedYearWeightedNEvents_photonScaleFactorDown;
  std::map<int, std::map<int, TH2F* > > h_lumiBasedYearWeightedNEvents_photonScaleFactorUp;

  // shifted distributions
  // syntax: histograms[shiftType][regionIndex][nJetsBin] where shiftType is a predefined enum used in the event selection
  std::map< shiftType, std::map<int, std::map<int, TH2I* > > > h_totalNEvents_shifted;
  // no need for weighted histograms for shifted distributions
};

struct inputNEventsStruct{
  std::map<std::string, int> data;

  inputNEventsStruct(std::string inputNEventsFileName) {
    std::string inputLine;
    std::ifstream inputNEventsFileObject(inputNEventsFileName.c_str());
    if (inputNEventsFileObject.is_open()) {
      while (std::getline(inputNEventsFileObject, inputLine)) {
        if (inputLine.length() >= 1) {
          std::vector<std::string> splitOnSpace = tmMiscUtils::getSplitString(inputLine, " ");
          if (!(splitOnSpace.size() == 2)) {
            std::cout << "ERROR: More than one space present in input line: " << inputLine << std::endl;
            std::exit(EXIT_FAILURE);
          }
          std::vector<std::string> splitOnEqualSign = tmMiscUtils::getSplitString(splitOnSpace[1], "=");
          if (!(splitOnEqualSign.size() == 2)) {
            std::cout << "ERROR: More than one equal to sign present in assignment: " << splitOnSpace[1] << std::endl;
            std::exit(EXIT_FAILURE);
          }
          data[splitOnEqualSign[0]] = std::stoi(splitOnEqualSign[1]);
          std::cout << "Assigned: data[\"" << splitOnEqualSign[0] << "\"] = " << data[splitOnEqualSign[0]] << std::endl;
        }
      }
      inputNEventsFileObject.close();
    }
    else {
      std::cout << "ERROR: Unable to open file with name = " << inputNEventsFileName << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
};

struct inputDataUncertaintiesStruct{
  std::map<std::string, float> data;

  inputDataUncertaintiesStruct(std::string inputDataUncertaintiesFileName) {
    std::string inputLine;
    std::ifstream inputDataUncertaintiesFileObject(inputDataUncertaintiesFileName.c_str());
    if (inputDataUncertaintiesFileObject.is_open()) {
      while (std::getline(inputDataUncertaintiesFileObject, inputLine)) {
        if (inputLine.length() >= 1) {
          std::vector<std::string> splitOnSpace = tmMiscUtils::getSplitString(inputLine, " ");
          if (!(splitOnSpace.size() == 2)) {
            std::cout << "ERROR: More than one space present in input line: " << inputLine << std::endl;
            std::exit(EXIT_FAILURE);
          }
          std::vector<std::string> splitOnEqualSign = tmMiscUtils::getSplitString(splitOnSpace[1], "=");
          if (!(splitOnEqualSign.size() == 2)) {
            std::cout << "ERROR: More than one equal to sign present in assignment: " << splitOnSpace[1] << std::endl;
            std::exit(EXIT_FAILURE);
          }
          data[splitOnEqualSign[0]] = std::stof(splitOnEqualSign[1]);
          std::cout << "Assigned: data[\"" << splitOnEqualSign[0] << "\"] = " << data[splitOnEqualSign[0]] << std::endl;
        }
      }
      inputDataUncertaintiesFileObject.close();
    }
    else {
      std::cout << "ERROR: Unable to open file with name = " << inputDataUncertaintiesFileName << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
};

struct inputDataSTScalingUncertaintiesStruct{
  std::map<std::string, float> data;

  inputDataSTScalingUncertaintiesStruct(std::string inputDataSTScalingUncertaintiesFileName) {
    std::string inputLine;
    std::ifstream inputDataSTScalingUncertaintiesFileObject(inputDataSTScalingUncertaintiesFileName.c_str());
    if (inputDataSTScalingUncertaintiesFileObject.is_open()) {
      while (std::getline(inputDataSTScalingUncertaintiesFileObject, inputLine)) {
        if (inputLine.length() >= 1) {
          std::vector<std::string> splitOnSpace = tmMiscUtils::getSplitString(inputLine, " ");
          if (!(splitOnSpace.size() == 2)) {
            std::cout << "ERROR: More than one space present in input line: " << inputLine << std::endl;
            std::exit(EXIT_FAILURE);
          }
          std::vector<std::string> splitOnEqualSign = tmMiscUtils::getSplitString(splitOnSpace[1], "=");
          if (!(splitOnEqualSign.size() == 2)) {
            std::cout << "ERROR: More than one equal to sign present in assignment: " << splitOnSpace[1] << std::endl;
            std::exit(EXIT_FAILURE);
          }
          data[splitOnEqualSign[0]] = std::stof(splitOnEqualSign[1]);
          std::cout << "Assigned: data[\"" << splitOnEqualSign[0] << "\"] = " << data[splitOnEqualSign[0]] << std::endl;
        }
      }
      inputDataSTScalingUncertaintiesFileObject.close();
    }
    else {
      std::cout << "ERROR: Unable to open file with name = " << inputDataSTScalingUncertaintiesFileName << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
};

optionsStruct getOptionsFromParser(tmArgumentParser& argumentParser) {
  optionsStruct options = optionsStruct();
  options.inputPath = argumentParser.getArgumentString("inputPath");
  options.MCTemplatePath = argumentParser.getArgumentString("MCTemplatePath");
  options.inputFile_STRegionBoundaries = argumentParser.getArgumentString("inputFile_STRegionBoundaries");
  options.inputNEventsFile = argumentParser.getArgumentString("inputNEventsFile");
  options.inputDataUncertaintiesFile = argumentParser.getArgumentString("inputDataUncertaintiesFile");
  options.inputDataSTScalingUncertaintiesFile = argumentParser.getArgumentString("inputDataSTScalingUncertaintiesFile");
  options.outputDirectory = argumentParser.getArgumentString("outputDirectory");
  options.outputDirectory_signalContamination = argumentParser.getArgumentString("outputDirectory_signalContamination");
  options.outputPrefix = argumentParser.getArgumentString("outputPrefix");
  if (argumentParser.getArgumentString("getSignalContaminationOutsideSidebands") == "true") options.getSignalContaminationOutsideSidebands = true;
  else options.getSignalContaminationOutsideSidebands = false;
  return options;
}

std::string getHistogramName(const std::string& histogramType, const int& STRegionIndex, const int& nJetsBin) {
  std::stringstream nameStream;
  nameStream << histogramType << "_STRegion" << STRegionIndex << "_" << nJetsBin << "Jets";
  return nameStream.str();
}

std::string getHistogramTitle(const std::string& histogramType, const int& STRegionIndex, const int& nJetsBin, const STRegionsStruct& STRegions, const std::string& UpDownShift) {
  std::string histogramTypeString;
  if (histogramType == "signalContamination") histogramTypeString = "Signal Contamination";
  else if (histogramType == "MCStatisticsFractionalError") histogramTypeString = "Fractional error due to MC statistics";
  else if (histogramType == "JECUncertainty") histogramTypeString = "Fractional error due to JEC uncertainty";
  else if (histogramType == "UnclusteredMETUncertainty") histogramTypeString = "Fractional error due to unclustered energy uncertainty";
  else if (histogramType == "JERMETUncertainty") histogramTypeString = "Fractional error due to jet energy resolution";
  else if (histogramType == "missingHEMUncertainty") histogramTypeString = "Fractional error due to missing HEM";
  else if (histogramType == "prefiringWeightsUncertainty") histogramTypeString = "Fractional error due to uncertainty on prefiring weights";
  else if (histogramType == "photonMCScaleFactorUncertainty") histogramTypeString = "Fractional error due to uncertainty on photon MC scale factor";
  else {
    std::cout << "ERROR: Unrecognized histogram type: " << histogramType << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (!(UpDownShift == "")) histogramTypeString += (", shift " + UpDownShift);

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
