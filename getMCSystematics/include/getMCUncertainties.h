#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

#include "tmArgumentParser.h"
#include "tmROOTSaverUtils.h"
#include "tmMiscellaneous.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH2F.h"
#include "TGraph2D.h"

#include "STRegionsStruct.h"
#include "../../eventSelection/include/shiftedObservablesStruct.h"

struct optionsStruct{
  std::string inputPath, MCTemplate, inputFile_STRegionBoundaries, inputNEventsFile, outputDirectory, outputPrefix;
  int nGluinoMassBins, nNeutralinoMassBins;
  double minGluinoMass, maxGluinoMass, minNeutralinoMass, maxNeutralinoMass;
};

struct outputHistogramsStruct{
  // syntax: outputHistogram[STRegionIndex][nJetsBin] where regionIndex ranges from 1 to (1 + number of ST signal bins), where regionIndex 1 corresponds to the normalization bin
  std::map<int, std::map< int, TH2F* > > h_JECUncertainty;
  std::map<int, std::map< int, TH2F* > > h_ratios_JECUpToNominal;
  std::map<int, std::map< int, TH2F* > > h_ratios_JECDownToNominal;
  std::map<int, std::map< int, TH2F* > > h_MCStatisticsFractionalError;
  std::map<int, std::map< int, TH2F* > > h_signalContamination;
};

struct inputHistogramsStruct{
  // syntax: inputHistograms[JEC][STRegionIndex][nJetsBin] where JEC belongs to allowedJECs and regionIndex ranges from 1 to (1 + number of ST signal bins), where regionIndex 1 corresponds to the normalization bin
  std::map< std::string, std::map<int, std::map< int, TH2F* > > > h_totalNEvents;
  std::map< std::string, std::map<int, std::map< int, TH2F* > > > h_weightedNEvents;
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
