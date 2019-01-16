#include "../include/getEventHistograms.h"

outputHistogramsStruct* initializeOutputHistograms(argumentsStruct& arguments, const STRegionsStruct& STRegions) {
  outputHistogramsStruct* outputHistograms = new outputHistogramsStruct();
  // 2D histograms for nEvents
  for (int STRegionIndex = 1; STRegionIndex <= (1+STRegions.nSTSignalBins); ++STRegionIndex) {
    for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
      outputHistograms->h_totalNEvents[STRegionIndex][nJetsBin] = new TH2I(("h_" + getHistogramName("totalNEvents", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("totalNEvents", STRegionIndex, nJetsBin, arguments, STRegions).c_str(), arguments.nGluinoMassBins, arguments.minGluinoMass, arguments.maxGluinoMass, arguments.nNeutralinoMassBins, arguments.minNeutralinoMass, arguments.maxNeutralinoMass);
      outputHistograms->h_totalNEvents[STRegionIndex][nJetsBin]->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson); // this is the main plot that will be used to estimate statistical MC uncertainties
      outputHistograms->h_lumiBasedYearWeightedNEvents[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("lumiBasedYearWeightedNEvents", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("lumiBasedYearWeightedNEvents", STRegionIndex, nJetsBin, arguments, STRegions).c_str(), arguments.nGluinoMassBins, arguments.minGluinoMass, arguments.maxGluinoMass, arguments.nNeutralinoMassBins, arguments.minNeutralinoMass, arguments.maxNeutralinoMass);
      for (int shiftTypeIndex = shiftTypeFirst; shiftTypeIndex != static_cast<int>(shiftType::nShiftTypes); ++shiftTypeIndex) {
        shiftType typeIndex = static_cast<shiftType>(shiftTypeIndex);
        outputHistograms->h_totalNEvents_shifted[typeIndex][STRegionIndex][nJetsBin] = new TH2I(("h_" + getHistogramName(typeIndex, "totalNEvents_shifted", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle(typeIndex, "totalNEvents_shifted", STRegionIndex, nJetsBin, arguments, STRegions).c_str(), arguments.nGluinoMassBins, arguments.minGluinoMass, arguments.maxGluinoMass, arguments.nNeutralinoMassBins, arguments.minNeutralinoMass, arguments.maxNeutralinoMass);
      }
    }
  }

  // 1D sT distributions
  for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
    for (const auto& keyValuePair : arguments.specialZonesFor_sTDistributions) {
      int specialZoneIndex = keyValuePair.first;
      outputHistograms->h_sTDistributions[specialZoneIndex][nJetsBin] = new TH1F(("h_" + getHistogramName("sTDistribution", specialZoneIndex, nJetsBin)).c_str(), getHistogramTitle("sTDistribution", specialZoneIndex, nJetsBin, arguments, STRegions).c_str(), arguments.n_sTBinsToPlot, STRegions.STNormRangeMin, arguments.sTMax_toPlot);
      for (int shiftTypeIndex = shiftTypeFirst; shiftTypeIndex != static_cast<int>(shiftType::nShiftTypes); ++shiftTypeIndex) {
        shiftType typeIndex = static_cast<shiftType>(shiftTypeIndex);
        outputHistograms->h_sTDistributions_shifted[specialZoneIndex][typeIndex][nJetsBin] = new TH1F(("h_" + getHistogramName(typeIndex, "sTDistribution", specialZoneIndex, nJetsBin)).c_str(), getHistogramTitle(typeIndex, "sTDistribution", specialZoneIndex, nJetsBin, arguments, STRegions).c_str(), arguments.n_sTBinsToPlot, STRegions.STNormRangeMin, arguments.sTMax_toPlot);
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

void fillOutputHistograms(outputHistogramsStruct *outputHistograms, argumentsStruct& arguments, std::map<int, double>& crossSections, const STRegionsStruct& STRegions ) {
  std::cout << "Filling events output histograms..." << std::endl;

  TFile *inputFile = TFile::Open((arguments.inputMCPathMain).c_str(), "READ");
  if (!inputFile || inputFile->IsZombie()) {
    std::cout << "Error in opening file with path: " << (arguments.inputMCPathMain) << std::endl;
  }
  TChain *inputChain = new TChain("ggNtuplizer/EventTree");
  inputChain->Add((arguments.inputMCPathMain).c_str());
  long nEntriesMain = inputChain->GetEntries();
  std::cout << "Number of entries in main MC chain: " << nEntriesMain << std::endl;
  inputChain->Add((arguments.inputMCPathAux).c_str());
  long nEntriesSource = inputChain->GetEntries();
  // long nEntriesSource = inputTree->GetEntries();
  std::cout << "Total number of entries found: " << nEntriesSource << std::endl;
  // TTree *inputTree = (TTree*) inputFile->Get("ggNtuplizer/EventTree");
  // TTreeReader inputTreeReader(inputTree);
  TTreeReader inputTreeReader(inputChain);
  TTreeReaderValue<int> evt_nJets(inputTreeReader, "b_nJets");
  TTreeReaderValue<int> evt_nJets_shifted_JECDown(inputTreeReader, getShiftedVariableBranchName(shiftType::JECDown, "nJets").c_str());
  TTreeReaderValue<int> evt_nJets_shifted_JECUp(inputTreeReader, getShiftedVariableBranchName(shiftType::JECUp, "nJets").c_str());
  TTreeReaderValue<float> evt_ST(inputTreeReader, "b_evtST");
  TTreeReaderValue<float> evt_ST_shifted_JECDown(inputTreeReader, getShiftedVariableBranchName(shiftType::JECDown, "evtST").c_str());
  TTreeReaderValue<float> evt_ST_shifted_JECUp(inputTreeReader, getShiftedVariableBranchName(shiftType::JECUp, "evtST").c_str());
  TTreeReaderValue<float> evt_ST_shifted_UnclusteredMETDown(inputTreeReader, getShiftedVariableBranchName(shiftType::UnclusteredMETDown, "evtST").c_str());
  TTreeReaderValue<float> evt_ST_shifted_UnclusteredMETUp(inputTreeReader, getShiftedVariableBranchName(shiftType::UnclusteredMETUp, "evtST").c_str());
  TTreeReaderValue<float> evt_ST_shifted_JERMETDown(inputTreeReader, getShiftedVariableBranchName(shiftType::JERMETDown, "evtST").c_str());
  TTreeReaderValue<float> evt_ST_shifted_JERMETUp(inputTreeReader, getShiftedVariableBranchName(shiftType::JERMETUp, "evtST").c_str());
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

    // if (*evt_nJets < 2) {
    //   std::cout << "ERROR: Fewer than 2 jets in the event, should not have passed selection..." << std::endl;
    //   std::exit(EXIT_FAILURE);
    // } // commenting out because nJets distribution could be slightly different with JEC corrections
    // sanity checks end
    
    if (entryIndex % progressBarUpdatePeriod == 0) progressBar->updateBar(1.0*entryIndex/nEntriesToRead, entryIndex);
    // std::cout << "At entry index = " << entryIndex << ", nJets = " << *evt_nJets << ", ST = " << *evt_ST << ", nMC = " << *nMC << ", mcPID size: " << mcPIDs.GetSize() << ", mcMomPIDs size: " << mcMomPIDs.GetSize() << std::endl;

    int STRegionIndex = (STRegions.STAxis).FindFixBin(*evt_ST);
    int STRegionIndex_shifted_JECDown = (STRegions.STAxis).FindFixBin(*evt_ST_shifted_JECDown);
    int STRegionIndex_shifted_JECUp = (STRegions.STAxis).FindFixBin(*evt_ST_shifted_JECUp);
    int STRegionIndex_shifted_UnclusteredMETDown = (STRegions.STAxis).FindFixBin(*evt_ST_shifted_UnclusteredMETDown);
    int STRegionIndex_shifted_UnclusteredMETUp = (STRegions.STAxis).FindFixBin(*evt_ST_shifted_UnclusteredMETUp);
    int STRegionIndex_shifted_JERMETDown = (STRegions.STAxis).FindFixBin(*evt_ST_shifted_JERMETDown);
    int STRegionIndex_shifted_JERMETUp = (STRegions.STAxis).FindFixBin(*evt_ST_shifted_JERMETUp);
    if (STRegionIndex == 0) continue;
    if (STRegionIndex > (1+STRegions.nSTSignalBins)) {
      std::cout << "ERROR: unexpected region index: " << STRegionIndex << std::endl;
      std::exit(EXIT_FAILURE);
    }

    int nJetsBin = *evt_nJets;
    if (*evt_nJets > 6) nJetsBin = 6;
    int nJetsBin_JECDown = *evt_nJets_shifted_JECDown;
    if (*evt_nJets_shifted_JECDown > 6) nJetsBin_JECDown = 6;
    int nJetsBin_JECUp = *evt_nJets_shifted_JECUp;
    if (*evt_nJets_shifted_JECUp > 6) nJetsBin_JECUp = 6;

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
    double eventWeight = (*scaleFactor)*crossSections[gMassInt]*(arguments.integratedLuminosityMain + arguments.integratedLuminosityAux)/arguments.nGeneratedEventsPerBin;// factor of 4 to account for the fact that the MC production assumes a 50% branching ratio of neutralino to photon, while the analysis assumes 100% <-- 07 Jan 2019: factor of 4 removed until better understood...
    double yearWeight;
    bool isAux = false;
    if ((entryIndex >= nEntriesMain)) {
      isAux = true;
      yearWeight = arguments.integratedLuminosityAux/(arguments.integratedLuminosityMain + arguments.integratedLuminosityAux);
    }
    else {
      yearWeight = arguments.integratedLuminosityMain/(arguments.integratedLuminosityMain + arguments.integratedLuminosityAux);
    }

    double tmp_gM = static_cast<double>(generated_gluinoMass);
    double tmp_nM = static_cast<double>(generated_neutralinoMass);

    if (nJetsBin >= 2) outputHistograms->h_lumiBasedYearWeightedNEvents[STRegionIndex][nJetsBin]->Fill(tmp_gM, tmp_nM, eventWeight*yearWeight);

    if (!(isAux)) {// Fill total nEvent 2D histograms and ST distributions only from main MC sample
      if (nJetsBin >= 2) outputHistograms->h_totalNEvents[STRegionIndex][nJetsBin]->Fill(tmp_gM, tmp_nM);
      if (nJetsBin_JECDown >= 2) outputHistograms->h_totalNEvents_shifted[shiftType::JECDown][STRegionIndex_shifted_JECDown][nJetsBin_JECDown]->Fill(tmp_gM, tmp_nM);
      if (nJetsBin_JECUp >= 2) outputHistograms->h_totalNEvents_shifted[shiftType::JECUp][STRegionIndex_shifted_JECUp][nJetsBin_JECUp]->Fill(tmp_gM, tmp_nM);
      if (nJetsBin >= 2) {
        outputHistograms->h_totalNEvents_shifted[shiftType::UnclusteredMETDown][STRegionIndex_shifted_UnclusteredMETDown][nJetsBin]->Fill(tmp_gM, tmp_nM);
        outputHistograms->h_totalNEvents_shifted[shiftType::UnclusteredMETUp][STRegionIndex_shifted_UnclusteredMETUp][nJetsBin]->Fill(tmp_gM, tmp_nM);
        outputHistograms->h_totalNEvents_shifted[shiftType::JERMETDown][STRegionIndex_shifted_JERMETDown][nJetsBin]->Fill(tmp_gM, tmp_nM);
        outputHistograms->h_totalNEvents_shifted[shiftType::JERMETUp][STRegionIndex_shifted_JERMETUp][nJetsBin]->Fill(tmp_gM, tmp_nM);
      }

      for (const auto& keyValuePair : arguments.specialZonesFor_sTDistributions) {
        int specialZoneIndex = keyValuePair.first;
        if (arguments.specialZonesFor_sTDistributions[specialZoneIndex].contains(generated_gluinoMass, generated_neutralinoMass)) {
          if (nJetsBin >= 2) outputHistograms->h_sTDistributions[specialZoneIndex][nJetsBin]->Fill(*evt_ST, eventWeight);
          if (nJetsBin_JECDown >= 2) outputHistograms->h_sTDistributions_shifted[specialZoneIndex][shiftType::JECDown][nJetsBin_JECDown]->Fill(*evt_ST_shifted_JECDown, eventWeight);
          if (nJetsBin_JECUp >= 2) outputHistograms->h_sTDistributions_shifted[specialZoneIndex][shiftType::JECUp][nJetsBin_JECUp]->Fill(*evt_ST_shifted_JECUp, eventWeight);
          if (nJetsBin >= 2) {
            outputHistograms->h_sTDistributions_shifted[specialZoneIndex][shiftType::UnclusteredMETDown][nJetsBin]->Fill(*evt_ST_shifted_UnclusteredMETDown, eventWeight);
            outputHistograms->h_sTDistributions_shifted[specialZoneIndex][shiftType::UnclusteredMETUp][nJetsBin]->Fill(*evt_ST_shifted_UnclusteredMETUp, eventWeight);
            outputHistograms->h_sTDistributions_shifted[specialZoneIndex][shiftType::JERMETDown][nJetsBin]->Fill(*evt_ST_shifted_JERMETDown, eventWeight);
            outputHistograms->h_sTDistributions_shifted[specialZoneIndex][shiftType::JERMETUp][nJetsBin]->Fill(*evt_ST_shifted_JERMETUp, eventWeight);
          }
        }
      }
    }
    ++entryIndex;
  }
  delete progressBar;
  inputFile->Close();
}

void saveHistograms(outputHistogramsStruct *outputHistograms, argumentsStruct& arguments, const STRegionsStruct& STRegions) {
  std::cout << "Saving histograms..." << std::endl;

  // First the 2D event histograms
  TFile *outputFile = TFile::Open((arguments.outputDirectory + "/" + arguments.outputPrefix + "_savedObjects.root").c_str(), "RECREATE");
  for (int STRegionIndex = 1; STRegionIndex <= (1+STRegions.nSTSignalBins); ++STRegionIndex) {
    for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
      std::string histogramName_total = getHistogramName("totalNEvents", STRegionIndex, nJetsBin);
      tmROOTSaverUtils::saveSingleObject(outputHistograms->h_totalNEvents[STRegionIndex][nJetsBin], "c_" + histogramName_total, outputFile, arguments.outputDirectory + "/" + arguments.outputPrefix + "_" + histogramName_total + ".png", 1024, 768, 0, ".0f", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);
      std::string histogramName_weighted = getHistogramName("lumiBasedYearWeightedNEvents", STRegionIndex, nJetsBin);
      tmROOTSaverUtils::saveSingleObject(outputHistograms->h_lumiBasedYearWeightedNEvents[STRegionIndex][nJetsBin], "c_" + histogramName_weighted, outputFile, arguments.outputDirectory + "/" + arguments.outputPrefix + "_" + histogramName_weighted + ".png", 1024, 768, 0, ".0f", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);
      for (int shiftTypeIndex = shiftTypeFirst; shiftTypeIndex != static_cast<int>(shiftType::nShiftTypes); ++shiftTypeIndex) {
        shiftType typeIndex = static_cast<shiftType>(shiftTypeIndex);
        std::string histogramName_shifted = getHistogramName(typeIndex, "totalNEvents_shifted", STRegionIndex, nJetsBin);
        tmROOTSaverUtils::saveSingleObject(outputHistograms->h_totalNEvents_shifted[typeIndex][STRegionIndex][nJetsBin], "c_" + histogramName_shifted, outputFile, arguments.outputDirectory + "/" + arguments.outputPrefix + "_" + histogramName_shifted + ".png", 1024, 768, 0, ".0f", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);
      }
    }
  }

  // next the 1D sT distributions
  for (const auto& keyValuePair : arguments.specialZonesFor_sTDistributions) {
    int specialZoneIndex = keyValuePair.first;
    for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
      std::stringstream histogramNameStream;
      histogramNameStream << "sTDistribution_specialZone" << specialZoneIndex << "_" << nJetsBin << "Jets";

      TObjArray *sTDistributionsArray = new TObjArray(2+static_cast<int>(shiftType::nShiftTypes)); // one nominal ST distribution + nShiftTypes different distributions with shifts + one legend = 2+nTypes objects
      TLegend *legend = new TLegend(0.1, 0.7, 0.4, 0.9);
      outputHistograms->h_sTDistributions[specialZoneIndex][nJetsBin]->GetYaxis()->SetTitleOffset(1.4);
      sTDistributionsArray->AddAt(outputHistograms->h_sTDistributions[specialZoneIndex][nJetsBin], 0);
      TLegendEntry *legendEntry = legend->AddEntry(outputHistograms->h_sTDistributions[specialZoneIndex][nJetsBin], "Nominal");
      legendEntry->SetTextColor(kBlack);
      for (int shiftTypeIndex = shiftTypeFirst; shiftTypeIndex != static_cast<int>(shiftType::nShiftTypes); ++shiftTypeIndex) {
        shiftType typeIndex = static_cast<shiftType>(shiftTypeIndex);
        sTDistributionsArray->AddAt(outputHistograms->h_sTDistributions_shifted[specialZoneIndex][typeIndex][nJetsBin], 1+shiftTypeIndex);
        outputHistograms->h_sTDistributions_shifted[specialZoneIndex][typeIndex][nJetsBin]->SetLineColor(shiftTypeColors[typeIndex]);
        TLegendEntry *legendEntry = legend->AddEntry(outputHistograms->h_sTDistributions_shifted[specialZoneIndex][typeIndex][nJetsBin], ("Shifted: " + shiftTypeNames[typeIndex]).c_str());
        legendEntry->SetTextColor(shiftTypeColors[typeIndex]);
      }
      legend->SetFillStyle(0);
      sTDistributionsArray->AddAt(legend, 1+static_cast<int>(shiftType::nShiftTypes));
      tmROOTSaverUtils::saveObjects(sTDistributionsArray, "c_" + histogramNameStream.str(), outputFile, arguments.outputDirectory + "/" + arguments.outputPrefix + "_" + histogramNameStream.str() + ".png", 1024, 768, 0, "", "", false, false, false, 0, 0, 0, 0, 0, 0);
    }
  }
  outputFile->Close();
}

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  std::cout << "Starting code to calculate systematics due to uncertainty on jet energy corrections..." << std::endl;
  std::cout << "Current working directory: " << tmMiscUtils::getCWD() << std::endl;

  tmArgumentParser argumentParser = tmArgumentParser("Calculate systematics due to uncertainty on jet energy corrections.");
  argumentParser.addArgument("inputMCPathMain", "root://cmseos.fnal.gov//store/user/lpcsusystealth/selections/combinedSignal/MC_2018Production_DoubleMedium_optimized2017.root", true, "Path to 2017-optimized MC with no change to jet energies.");
  argumentParser.addArgument("integratedLuminosityMain", "41900.0", false, "Integrated luminosity in main MC reference.");
  argumentParser.addArgument("inputMCPathAux", "root://cmseos.fnal.gov//store/user/lpcsusystealth/selections/combinedSignal/MC_2018Production_DoubleMedium_optimized2016.root", true, "Path to 2016-optimized MC with no change to jet energies.");
  argumentParser.addArgument("integratedLuminosityAux", "35920.0", false, "Integrated luminosity in auxiliary MC reference.");
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
  outputHistogramsStruct* outputHistograms = initializeOutputHistograms(arguments, STRegions);
  std::map<int, double> crossSections = getCrossSectionsFromFile(arguments.crossSectionsFilePath);
  fillOutputHistograms(outputHistograms, arguments, crossSections, STRegions);
  saveHistograms(outputHistograms, arguments, STRegions);
  return 0;
}
