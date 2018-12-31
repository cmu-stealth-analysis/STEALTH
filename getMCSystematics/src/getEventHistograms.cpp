#include "../include/getEventHistograms.h"

std::string getHistogramName(std::string histogramType, std::string jec, int zoneID, int nJetsBin) {
  std::stringstream nameStream;
  if (histogramType == "sTDistribution") {
    nameStream << histogramType << "_" << jec << "_specialZone" << zoneID << "_" << nJetsBin << "Jets";
    return nameStream.str();
  }
  if (histogramType == "lumiBasedYearWeighted") nameStream << histogramType << "_nMCEvents_" << nJetsBin << "Jets_STRegion" << zoneID;
  else if ((histogramType == "averagePrescaleWeights") || (histogramType == "prescaleWeights1D")) nameStream << histogramType << "_" << nJetsBin << "Jets_STRegion" << zoneID;
  else nameStream << histogramType << "_nMCEvents_" << jec << "_" << nJetsBin << "Jets_STRegion" << zoneID;
  return nameStream.str();
}

std::string getHistogramTitle(std::string histogramType, std::string jec, int zoneID, int nJetsBin, argumentsStruct& arguments, const STRegionsStruct& STRegions) {
  if (histogramType == "sTDistribution") {
    std::stringstream histogramTitle;
    histogramTitle << "sT Distributions, " << std::fixed << std::setprecision(0)
                   << arguments.specialZonesFor_sTDistributions[zoneID].minGluinoMass
                   << " < m_{#tilde{#it{g}}} < "
                   << arguments.specialZonesFor_sTDistributions[zoneID].maxGluinoMass << ", "
                   << arguments.specialZonesFor_sTDistributions[zoneID].minNeutralinoMass
                   << " < m_{#tilde{#it{#chi_{1}^{0}}}} < "
                   << arguments.specialZonesFor_sTDistributions[zoneID].maxNeutralinoMass
                   << ", " << nJetsBin << " Jets;#it{S}_{T}(GeV);Weighted nEvents/("
                   << static_cast<int>(0.5 + ((arguments.sTMax_toPlot - STRegions.STNormRangeMin)/arguments.n_sTBinsToPlot))
                   << " GeV)";
    return histogramTitle.str();
  }
  std::string histogramTypeString;
  if (histogramType == "total") histogramTypeString = "Total MC Events";
  else if (histogramType == "weighted") histogramTypeString = "Weighted MC Events";
  else if (histogramType == "lumiBasedYearWeighted") histogramTypeString = "Lumi-based weighted MC Events";
  else if (histogramType == "averagePrescaleWeights") histogramTypeString = "Prescale weights profile";
  else if (histogramType == "prescaleWeights1D") histogramTypeString = "Prescale weights histogram";
  else {
    std::cout << "ERROR: Unrecognized histogram type: " << histogramType << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::string jecString;
  if (jec == "JECDown") jecString = "shift down by JECUnc";
  else if (jec == "JECNominal") jecString = "no JEC shift";
  else if (jec == "JECUp") jecString = "shift up by JECUnc";
  else {
    std::cout << "ERROR: Unrecognized jec: " << jec << std::endl;
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
  if (zoneID >= 1 && zoneID < (1+STRegions.nSTSignalBins)){ // all except the last bin
    sTRangeStringStream << std::fixed << std::setprecision(0) << (STRegions.STAxis).GetBinLowEdge(zoneID) << " < #it{S}_{T} < " << (STRegions.STAxis).GetBinUpEdge(zoneID);
  }
  else if (zoneID == (1+STRegions.nSTSignalBins)) {
    sTRangeStringStream << "#it{S}_{T} > " << std::fixed << std::setprecision(0) << (STRegions.STAxis).GetBinLowEdge(zoneID);
  }
  else {
    std::cout << "ERROR: unexpected zone ID: " << zoneID << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string sTRangeString = sTRangeStringStream.str();

  std::string axesLabelsString = ";m_{#tilde{#it{g}}};m_{#tilde{#it{#chi_{1}^{0}}}}";

  std::stringstream titleStream;
  if ((histogramType == "lumiBasedYearWeighted") || (histogramType == "averagePrescaleWeights")) titleStream << histogramTypeString << ", " << nJetsString << ", " << sTRangeString << axesLabelsString;
  else if (histogramType == "prescaleWeights1D") titleStream << histogramTypeString << ", " << nJetsString << ", " << sTRangeString << ";prescale weight;nEvents";
  else titleStream << histogramTypeString << ", " << jecString << ", " << nJetsString << ", " << sTRangeString << axesLabelsString;
  return titleStream.str();
}

argumentsStruct getArgumentsFromParser(tmArgumentParser& argumentParser) {
  argumentsStruct arguments = argumentsStruct();
  arguments.inputMCPathMain = argumentParser.getArgumentString("inputMCPathMain");
  arguments.integratedLuminosityMain = std::stod(argumentParser.getArgumentString("integratedLuminosityMain"));
  arguments.inputMCPathAux = argumentParser.getArgumentString("inputMCPathAux");
  arguments.integratedLuminosityAux = std::stod(argumentParser.getArgumentString("integratedLuminosityAux"));
  arguments.inputMCPath_JECUp = argumentParser.getArgumentString("inputMCPath_JECUp");
  arguments.inputMCPath_JECDown = argumentParser.getArgumentString("inputMCPath_JECDown");
  arguments.maxMCEvents = std::stol(argumentParser.getArgumentString("maxMCEvents"));
  arguments.crossSectionsFilePath = argumentParser.getArgumentString("crossSectionsFilePath");
  arguments.inputFile_STRegionBoundaries = argumentParser.getArgumentString("inputFile_STRegionBoundaries");
  arguments.sTMax_toPlot = std::stod(argumentParser.getArgumentString("sTMax_toPlot"));
  arguments.n_sTBinsToPlot = std::stoi(argumentParser.getArgumentString("n_sTBinsToPlot"));
  arguments.outputDirectory = argumentParser.getArgumentString("outputDirectory");
  arguments.outputPrefix = argumentParser.getArgumentString("outputPrefix");
  arguments.nGeneratedEventsPerBin = std::stoi(argumentParser.getArgumentString("nGeneratedEventsPerBin"));
  arguments.nGluinoMassBins = std::stoi(argumentParser.getArgumentString("nGluinoMassBins"));
  arguments.minGluinoMass = std::stod(argumentParser.getArgumentString("minGluinoMass"));
  arguments.maxGluinoMass = std::stod(argumentParser.getArgumentString("maxGluinoMass"));
  arguments.nNeutralinoMassBins = std::stoi(argumentParser.getArgumentString("nNeutralinoMassBins"));
  arguments.minNeutralinoMass = std::stod(argumentParser.getArgumentString("minNeutralinoMass"));
  arguments.maxNeutralinoMass = std::stod(argumentParser.getArgumentString("maxNeutralinoMass"));
  std::vector<std::string> regionArguments = tmMiscUtils::getSplitString(argumentParser.getArgumentString("regionsIn_sTHistograms"), "|");
  int specialZoneIndex = 1;
  for (const auto& regionArgument : regionArguments) {
    std::vector<std::string> massBoundaries = tmMiscUtils::getSplitString(regionArgument, ":");
    if (!(massBoundaries.size() == 4)) {
      std::cout << "ERROR: passed mass boundaries must have exactly four values in the format minGluinoMass:maxGluinoMass:minNeutralinoMas:maxNeutralinoMass. The problematic boundaries set passed is: " << regionArgument << std::endl;
      std::exit(EXIT_FAILURE);
    }
    arguments.specialZonesFor_sTDistributions[specialZoneIndex] = parameterSpaceRegion();
    arguments.specialZonesFor_sTDistributions[specialZoneIndex].setParameters(std::stod(massBoundaries[0]), std::stod(massBoundaries[1]), std::stod(massBoundaries[2]), std::stod(massBoundaries[3]));
    ++specialZoneIndex;
  }
  return arguments;
}

outputHistogramsStruct* initializeOutputHistograms(argumentsStruct& arguments, const std::vector<std::string>& allowedJECs, const STRegionsStruct& STRegions) {
  outputHistogramsStruct* outputHistograms = new outputHistogramsStruct();
  // 2D histograms for nEvents
  for (int STRegionIndex = 1; STRegionIndex <= (1+STRegions.nSTSignalBins); ++STRegionIndex) {
    for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
      for (const auto& jec: allowedJECs) {
        outputHistograms->h_totalNEvents[jec][STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("total", jec, STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("total", jec, STRegionIndex, nJetsBin, arguments, STRegions).c_str(), arguments.nGluinoMassBins, arguments.minGluinoMass, arguments.maxGluinoMass, arguments.nNeutralinoMassBins, arguments.minNeutralinoMass, arguments.maxNeutralinoMass);
        outputHistograms->h_totalNEvents[jec][STRegionIndex][nJetsBin]->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
        outputHistograms->h_weightedNEvents[jec][STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("weighted", jec, STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("weighted", jec, STRegionIndex, nJetsBin, arguments, STRegions).c_str(), arguments.nGluinoMassBins, arguments.minGluinoMass, arguments.maxGluinoMass, arguments.nNeutralinoMassBins, arguments.minNeutralinoMass, arguments.maxNeutralinoMass);
      }
      outputHistograms->h_lumiBasedYearWeightedNEvents[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("lumiBasedYearWeighted", "JECNominal", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("lumiBasedYearWeighted", "JECNominal", STRegionIndex, nJetsBin, arguments, STRegions).c_str(), arguments.nGluinoMassBins, arguments.minGluinoMass, arguments.maxGluinoMass, arguments.nNeutralinoMassBins, arguments.minNeutralinoMass, arguments.maxNeutralinoMass);
      outputHistograms->h_averagePrescaleWeights[STRegionIndex][nJetsBin] = new TProfile2D(("h_" + getHistogramName("averagePrescaleWeights", "JECNominal", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("averagePrescaleWeights", "JECNominal", STRegionIndex, nJetsBin, arguments, STRegions).c_str(), arguments.nGluinoMassBins, arguments.minGluinoMass, arguments.maxGluinoMass, arguments.nNeutralinoMassBins, arguments.minNeutralinoMass, arguments.maxNeutralinoMass);
      outputHistograms->h_prescaleWeights1D[STRegionIndex][nJetsBin] = new TH1F(("h_" + getHistogramName("prescaleWeights1D", "JECNominal", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("prescaleWeights1D", "JECNominal", STRegionIndex, nJetsBin, arguments, STRegions).c_str(), 510, -0.01, 1.01);
    }
  }

  // 1D sT distributions
  for (const auto& jec: allowedJECs) {
    for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
      for (const auto& keyValuePair : arguments.specialZonesFor_sTDistributions) {
        int specialZoneIndex = keyValuePair.first;
        outputHistograms->h_sTDistributions[specialZoneIndex][jec][nJetsBin] = new TH1F(("h_" + getHistogramName("sTDistribution", jec, specialZoneIndex, nJetsBin)).c_str(), getHistogramTitle("sTDistribution", jec, specialZoneIndex, nJetsBin, arguments, STRegions).c_str(), arguments.n_sTBinsToPlot, STRegions.STNormRangeMin, arguments.sTMax_toPlot);
      }
    }
  }
  return outputHistograms;
}

std::map<int, double> getCrossSectionsFromFile(std::string crossSectionsFilePath, bool printCrossSections = false) {
  std::map<int, double> crossSections;
  std::ifstream crossSectionsFileStream(crossSectionsFilePath.c_str());
  if (!(crossSectionsFileStream.is_open())) {
    std::cout << "ERROR: Unable to open file: " << crossSectionsFilePath << std::endl;
    std::exit(EXIT_FAILURE);
  }
  int gluinoMass;
  double crossSection, crossSectionPercentUncertainty;
  while (crossSectionsFileStream >> gluinoMass >> crossSection >> crossSectionPercentUncertainty) {
    crossSections[gluinoMass] = crossSection;
    // crossSectionsFractionalUncertainty[gluinoMass] = 0.01*crossSectionPercentUncertainty;
    if (printCrossSections) {
      std::cout << "crossSections[" << gluinoMass << "] = " << crossSection << ", "
                << "crossSectionsFractionalUncertainty[" << gluinoMass << "] = " << 0.01*crossSectionPercentUncertainty << std::endl;
    }
  }
  return crossSections;
}

void fillOutputHistogramsForJEC(outputHistogramsStruct *outputHistograms, argumentsStruct& arguments, const std::string& jec, std::map<int, double>& crossSections, const STRegionsStruct& STRegions ) {
  std::cout << "Filling events map for JEC type: " << jec << std::endl;
  std::string inputPath;
  if (jec == "JECDown") inputPath = arguments.inputMCPath_JECDown;
  else if (jec == "JECNominal") inputPath = arguments.inputMCPathMain;
  else if (jec == "JECUp") inputPath = arguments.inputMCPath_JECUp;
  else {
    std::cout << "ERROR: unknown JEC type: \"" << jec << "\"" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  TFile *inputFile = TFile::Open(inputPath.c_str(), "READ");
  if (!inputFile || inputFile->IsZombie()) {
    std::cout << "Error in opening file with path: " << inputPath << std::endl;
  }
  TChain *inputChain = new TChain("ggNtuplizer/EventTree");
  inputChain->Add(inputPath.c_str());
  long nEntriesMain = inputChain->GetEntries();
  std::cout << "Number of entries in main MC chain: " << nEntriesMain << std::endl;
  if (jec == "JECNominal") inputChain->Add((arguments.inputMCPathAux).c_str());
  long nEntriesSource = inputChain->GetEntries();
  // long nEntriesSource = inputTree->GetEntries();
  std::cout << "Total number of entries found: " << nEntriesSource << std::endl;
  // TTree *inputTree = (TTree*) inputFile->Get("ggNtuplizer/EventTree");
  // TTreeReader inputTreeReader(inputTree);
  TTreeReader inputTreeReader(inputChain);
  TTreeReaderValue<int> evt_nJets(inputTreeReader, "b_nJets");
  TTreeReaderValue<float> evt_ST(inputTreeReader, "b_evtST");
  TTreeReaderValue<float> scaleFactor(inputTreeReader, "b_evtScaleFactor");
  TTreeReaderValue<int> nMC(inputTreeReader, "nMC");
  TTreeReaderArray<int> mcPIDs(inputTreeReader, "mcPID");
  TTreeReaderArray<float> mcMasses(inputTreeReader, "mcMass");
  TTreeReaderArray<int> mcMomPIDs(inputTreeReader, "mcMomPID");
  TTreeReaderArray<float> mcMomMasses(inputTreeReader, "mcMomMass");
  long nEntriesToRead = (arguments.maxMCEvents > 0) ? (arguments.maxMCEvents < nEntriesSource ? arguments.maxMCEvents : nEntriesSource) : nEntriesSource;
  std::cout << "Looping over " << nEntriesToRead << " entries." << std::endl;
  tmProgressBar *progressBar = new tmProgressBar(nEntriesToRead);
  int tmp = static_cast<int>(0.5 + 1.0*nEntriesToRead/1000);
  int progressBarUpdatePeriod = tmp > 1 ? tmp : 1;
  int entryIndex = 0;
  progressBar->initialize();
  while (inputTreeReader.Next()) {
    if (entryIndex >= nEntriesToRead) {
      std::cout << "ERROR: entry index seems to be greater than available number of events" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // sanity checks begin
    if (*nMC != static_cast<int>(mcPIDs.GetSize()) || *nMC != static_cast<int>(mcMomPIDs.GetSize())) {
      std::cout << "ERROR: Something has gone terribly wrong, number of MC particles is not the same as the size of the vectors containing MCPIDs" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    if (*evt_nJets < 2) {
      std::cout << "ERROR: Fewer than 2 jets in the event, should not have passed selection..." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // sanity checks end
    
    if (entryIndex % progressBarUpdatePeriod == 0) progressBar->updateBar(1.0*entryIndex/nEntriesToRead, entryIndex);
    // std::cout << "At entry index = " << entryIndex << ", nJets = " << *evt_nJets << ", ST = " << *evt_ST << ", nMC = " << *nMC << ", mcPID size: " << mcPIDs.GetSize() << ", mcMomPIDs size: " << mcMomPIDs.GetSize() << std::endl;

    int STRegionIndex = (STRegions.STAxis).FindFixBin(*evt_ST);
    if (STRegionIndex == 0) continue;
    if (STRegionIndex > (1+STRegions.nSTSignalBins)) {
      std::cout << "ERROR: unexpected region index: " << STRegionIndex << std::endl;
      std::exit(EXIT_FAILURE);
    }

    int nJetsBin = *evt_nJets;
    if (*evt_nJets > 6) nJetsBin = 6;

    // get generated gluino, neutralino mass
    float generated_gluinoMass = 0;
    bool gluinoMassIsSet = false;
    float generated_neutralinoMass = 0;
    bool neutralinoMassIsSet = false;
    for (int generatedParticleIndex = 0; generatedParticleIndex < *nMC; ++generatedParticleIndex) {
      int particle_mcPID = mcPIDs[generatedParticleIndex];
      if (!(gluinoMassIsSet) && particle_mcPID == MCPID_GLUINO) {
        generated_gluinoMass = mcMasses[generatedParticleIndex];
        gluinoMassIsSet = true;
      }
      int particle_mcMomPID = mcMomPIDs[generatedParticleIndex];
      if (!(neutralinoMassIsSet) && particle_mcMomPID == MCPID_NEUTRALINO) {
        generated_neutralinoMass = mcMomMasses[generatedParticleIndex];
        neutralinoMassIsSet = true;
      }
    }
    if (!(gluinoMassIsSet && neutralinoMassIsSet)) {
      std::cout << "ERROR: Unable to find gluino or neutralino mass in an event." << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // get event weight
    int gMassInt = static_cast<int>(0.5 + generated_gluinoMass);
    double eventWeight = (*scaleFactor)*crossSections[gMassInt]*(arguments.integratedLuminosityMain + arguments.integratedLuminosityAux)/arguments.nGeneratedEventsPerBin;
    double yearWeight;
    bool isAux = false;
    if ((entryIndex >= nEntriesMain)) {
      isAux = true;
      yearWeight = arguments.integratedLuminosityAux/(arguments.integratedLuminosityMain + arguments.integratedLuminosityAux);
    }
    else {
      yearWeight = arguments.integratedLuminosityMain/(arguments.integratedLuminosityMain + arguments.integratedLuminosityAux);
    }

    if (jec == "JECNominal") {// Fill lumi-based weighted nEvent 2D histograms and average prescale weights from both MC samples
      outputHistograms->h_lumiBasedYearWeightedNEvents[STRegionIndex][nJetsBin]->Fill(static_cast<double>(generated_gluinoMass), static_cast<double>(generated_neutralinoMass), eventWeight*yearWeight);
      outputHistograms->h_averagePrescaleWeights[STRegionIndex][nJetsBin]->Fill(static_cast<double>(generated_gluinoMass), static_cast<double>(generated_neutralinoMass), (*scaleFactor));
      outputHistograms->h_prescaleWeights1D[STRegionIndex][nJetsBin]->Fill((*scaleFactor));
    }

    if (!(isAux)) {// Fill weighted and total nEvent 2D histograms and ST distributions only from main MC sample
      outputHistograms->h_totalNEvents[jec][STRegionIndex][nJetsBin]->Fill(static_cast<double>(generated_gluinoMass), static_cast<double>(generated_neutralinoMass), 1.0);
      outputHistograms->h_weightedNEvents[jec][STRegionIndex][nJetsBin]->Fill(static_cast<double>(generated_gluinoMass), static_cast<double>(generated_neutralinoMass), eventWeight);
      for (const auto& keyValuePair : arguments.specialZonesFor_sTDistributions) {
        int specialZoneIndex = keyValuePair.first;
        if (arguments.specialZonesFor_sTDistributions[specialZoneIndex].contains(generated_gluinoMass, generated_neutralinoMass)) outputHistograms->h_sTDistributions[specialZoneIndex][jec][nJetsBin]->Fill(*evt_ST, eventWeight);
      }
    }
    ++entryIndex;
  }
  delete progressBar;
  inputFile->Close();
}

void fillOutputHistograms(outputHistogramsStruct *outputHistograms, argumentsStruct& arguments, const std::vector<std::string>& allowedJECs, std::map<int, double>& crossSections, const STRegionsStruct& STRegions) {
  for (const auto& jec: allowedJECs) {
    fillOutputHistogramsForJEC(outputHistograms, arguments, jec, crossSections, STRegions);
  }
}

void saveHistograms(outputHistogramsStruct *outputHistograms, argumentsStruct& arguments, const std::vector<std::string>& allowedJECs, const STRegionsStruct& STRegions) {
  std::cout << "Saving histograms..." << std::endl;

  // First the 2D event histograms
  TFile *outputFile = TFile::Open((arguments.outputDirectory + "/" + arguments.outputPrefix + "_savedObjects.root").c_str(), "RECREATE");
  for (const auto& jec: allowedJECs) {
    for (int STRegionIndex = 1; STRegionIndex <= (1+STRegions.nSTSignalBins); ++STRegionIndex) {
      for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
        std::string histogramName_total = getHistogramName("total", jec, STRegionIndex, nJetsBin);
        tmROOTSaverUtils::saveSingleObject(outputHistograms->h_totalNEvents[jec][STRegionIndex][nJetsBin], "c_" + histogramName_total, outputFile, arguments.outputDirectory + "/" + arguments.outputPrefix + "_" + histogramName_total + ".png", 1024, 768, 0, ".0f", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);
        std::string histogramName_weighted = getHistogramName("weighted", jec, STRegionIndex, nJetsBin);
        tmROOTSaverUtils::saveSingleObject(outputHistograms->h_weightedNEvents[jec][STRegionIndex][nJetsBin], "c_" + histogramName_weighted, outputFile, arguments.outputDirectory + "/" + arguments.outputPrefix + "_" + histogramName_weighted + ".png", 1024, 768, 0, ".0f", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);
        if (jec == "JECNominal") {
          std::string histogramName_lumiBasedYearWeighted = getHistogramName("lumiBasedYearWeighted", jec, STRegionIndex, nJetsBin);
          tmROOTSaverUtils::saveSingleObject(outputHistograms->h_lumiBasedYearWeightedNEvents[STRegionIndex][nJetsBin], "c_" + histogramName_lumiBasedYearWeighted, outputFile, arguments.outputDirectory + "/" + arguments.outputPrefix + "_" + histogramName_lumiBasedYearWeighted + ".png", 1024, 768, 0, ".0f", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);
          std::string histogramName_averagePrescaleWeights = getHistogramName("averagePrescaleWeights", jec, STRegionIndex, nJetsBin);
          tmROOTSaverUtils::saveSingleObject(outputHistograms->h_averagePrescaleWeights[STRegionIndex][nJetsBin], "c_" + histogramName_averagePrescaleWeights, outputFile, arguments.outputDirectory + "/" + arguments.outputPrefix + "_" + histogramName_averagePrescaleWeights + ".png", 1024, 768, 0, ".2f", "TEXTCOLZ", false, false, false, 0, 0, 0, 0, 0, 0);
          std::string histogramName_prescaleWeights1D = getHistogramName("prescaleWeights1D", jec, STRegionIndex, nJetsBin);
          TCanvas *prescaleWeightsCanvas = tmROOTSaverUtils::saveSingleObject(outputHistograms->h_prescaleWeights1D[STRegionIndex][nJetsBin], "c_" + histogramName_prescaleWeights1D, outputFile, arguments.outputDirectory + "/" + arguments.outputPrefix + "_" + histogramName_prescaleWeights1D + ".png", 1024, 768, 112211, "", "", false, true, false, 0, 0, 0, 0, 0, 0); // Ugly hack. This call ensures that Draw() is called...
          TPaveStats *statsBox = (TPaveStats*)(prescaleWeightsCanvas->GetPrimitive("stats"));// Fetching the statistics box object...
          statsBox->SetX1NDC(0.1); // ...setting xlow...
          statsBox->SetX2NDC(0.4); // ...setting xhigh...
          prescaleWeightsCanvas->SaveAs((arguments.outputDirectory + "/" + arguments.outputPrefix + "_" + histogramName_prescaleWeights1D + ".png").c_str()); // and finally replacing the output file with the version with the correct position for the statistics box.
        }
      }
    }
  }

  // next the 1D sT distributions
  for (const auto& keyValuePair : arguments.specialZonesFor_sTDistributions) {
    int specialZoneIndex = keyValuePair.first;
    for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
      std::stringstream histogramNameStream;
      histogramNameStream << "sTDistribution_specialZone" << specialZoneIndex << "_" << nJetsBin << "Jets";

      TObjArray *sTDistributionsArray = new TObjArray(4);
      TLegend *legend = new TLegend(0.1, 0.7, 0.4, 0.9);
      for (const auto& jec: allowedJECs) {
        int indexToAddAt = (jec == "JECNominal"? 0 : (jec == "JECUp" ? 1 : (jec == "JECDown" ? 2 : -1)));
        if (indexToAddAt < 0) {
          std::cout << "ERROR: bad indexToAddAt = " << indexToAddAt << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (jec == "JECNominal") outputHistograms->h_sTDistributions[specialZoneIndex][jec][nJetsBin]->GetYaxis()->SetTitleOffset(1.4);
        sTDistributionsArray->AddAt(outputHistograms->h_sTDistributions[specialZoneIndex][jec][nJetsBin], indexToAddAt);
        outputHistograms->h_sTDistributions[specialZoneIndex][jec][nJetsBin]->SetLineColor(jec == "JECNominal"? kBlack : (jec == "JECUp" ? kRed : (jec == "JECDown" ? kBlue : kGreen)));
        TLegendEntry *legendEntry = legend->AddEntry(outputHistograms->h_sTDistributions[specialZoneIndex][jec][nJetsBin], jec == "JECNominal" ? "No Jet {p}_{T} shift" : (jec == "JECUp" ? "Jet {p}_{T} shifted up by JEC uncertainty" : (jec == "JECDown" ? "Jet {p}_{T} shifted down by JEC uncertainty" : "Unknown JEC")));
        legendEntry->SetTextColor(jec == "JECNominal"? kBlack : (jec == "JECUp" ? kRed : (jec == "JECDown" ? kBlue : kGreen)));
      }
      legend->SetFillStyle(0);
      sTDistributionsArray->AddAt(legend, 3);
      tmROOTSaverUtils::saveObjects(sTDistributionsArray, "c_" + histogramNameStream.str(), outputFile, arguments.outputDirectory + "/" + arguments.outputPrefix + "_" + histogramNameStream.str() + ".png", 1024, 768, 0, "", "", false, false, false, 0, 0, 0, 0, 0, 0);
    }
  }
  outputFile->Close();
}

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  std::cout << "Starting code to calculate systematics due to uncertainty on jet energy corrections..." << std::endl;
  std::cout << "Current working directory: " << tmMiscUtils::getCWD() << std::endl;

  std::vector<std::string> allowedJECs{"JECDown", "JECNominal", "JECUp"};
  // std::vector<std::string> allowedZones{"norm", "sub", "main"};

  tmArgumentParser argumentParser = tmArgumentParser("Calculate systematics due to uncertainty on jet energy corrections.");
  argumentParser.addArgument("inputMCPathMain", "root://cmseos.fnal.gov//store/user/lpcsusystealth/selections/combinedSignal/MC_2018Production_DoubleMedium_optimized2017.root", true, "Path to 2017-optimized MC with no change to jet energies.");
  argumentParser.addArgument("integratedLuminosityMain", "41900.0", false, "Integrated luminosity in main MC reference.");
  argumentParser.addArgument("inputMCPathAux", "root://cmseos.fnal.gov//store/user/lpcsusystealth/selections/combinedSignal/MC_2018Production_DoubleMedium_optimized2016.root", true, "Path to 2016-optimized MC with no change to jet energies.");
  argumentParser.addArgument("integratedLuminosityAux", "35920.0", false, "Integrated luminosity in auxiliary MC reference.");
  argumentParser.addArgument("inputMCPath_JECUp", "root://cmseos.fnal.gov//store/user/lpcsusystealth/selections/combined/MC_2018Production_DoubleMedium_optimized2017_JECUp.root", true, "Path to MC with all energies shifted up by the JEC uncertainty.");
  argumentParser.addArgument("inputMCPath_JECDown", "root://cmseos.fnal.gov//store/user/lpcsusystealth/selections/combined/MC_2018Production_DoubleMedium_optimized2017_JECDown.root", true, "Path to MC with all energies shifted up by the JEC uncertainty.");
  argumentParser.addArgument("maxMCEvents", "0", false, "Set a custom maximum number of MC events.");
  argumentParser.addArgument("crossSectionsFilePath", "SusyCrossSections13TevGluGlu.txt", false, "Path to dat file that contains cross-sections as a function of gluino mass, to use while weighting events.");
  argumentParser.addArgument("inputFile_STRegionBoundaries", "STRegionBoundaries.dat", false, "Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.");
  argumentParser.addArgument("sTMax_toPlot", "4000.0", false, "Max value of sT to plot.");
  argumentParser.addArgument("n_sTBinsToPlot", "29", false, "Number of sT bins to plot."); // default: 23 bins from 1200 to 3500 GeV in steps of 100 GeV
  argumentParser.addArgument("outputDirectory", "analysis/MCEventHistograms/", false, "Prefix to output files.");
  argumentParser.addArgument("outputPrefix", "", true, "Prefix to output files.");
  argumentParser.addArgument("nGeneratedEventsPerBin", "150000", false, "Number of generated events per bin, to use while calculating event weights.");
  argumentParser.addArgument("nGluinoMassBins", "20", false, "nBins on the gluino mass axis."); // (800 - 25) GeV --> (1750 + 25) GeV in steps of 50 GeV
  argumentParser.addArgument("minGluinoMass", "775.0", false, "Min gluino mass for the 2D plots.");
  argumentParser.addArgument("maxGluinoMass", "1775.0", false, "Max gluino mass for the 2D plots.");
  argumentParser.addArgument("nNeutralinoMassBins", "133", false, "nBins on the neutralino mass axis.");
  argumentParser.addArgument("minNeutralinoMass", "93.75", false, "Min neutralino mass for the 2D plots.");
  argumentParser.addArgument("maxNeutralinoMass", "1756.25", false, "Max neutralino mass for the 2D plots."); // (100 - 6.25) GeV --> (1750 + 6.25) GeV in steps of 12.5 GeV
  argumentParser.addArgument("regionsIn_sTHistograms", "1725.0:1775.0:650.0:950.0|1025.0:1075.0:975.0:1075.0", false, "List of the regions in which to fill and save the sT histograms. Each element of the list is in format minGluinoMass:maxGluinoMass:minNeutralinoMas:maxNeutralinoMass, and each element is separated from the next by the character \"|\". See also the default value for this argument.");
  argumentParser.setPassedStringValues(argc, argv);
  argumentsStruct arguments = getArgumentsFromParser(argumentParser);

  STRegionsStruct STRegions(arguments.inputFile_STRegionBoundaries);
  outputHistogramsStruct* outputHistograms = initializeOutputHistograms(arguments, allowedJECs, STRegions);
  std::map<int, double> crossSections = getCrossSectionsFromFile(arguments.crossSectionsFilePath);
  fillOutputHistograms(outputHistograms, arguments, allowedJECs, crossSections, STRegions);
  saveHistograms(outputHistograms, arguments, allowedJECs, STRegions);
  return 0;
}
