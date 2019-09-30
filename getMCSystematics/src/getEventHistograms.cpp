#include "../include/getEventHistograms.h"

outputHistogramsStruct* initializeOutputHistograms(argumentsStruct& arguments, const STRegionsStruct& STRegions) {
  outputHistogramsStruct* outputHistograms = new outputHistogramsStruct();
  // 2D histograms for nEvents
  for (int STRegionIndex = 1; STRegionIndex <= (1+STRegions.nSTSignalBins); ++STRegionIndex) {
    for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
      outputHistograms->h_totalNEvents[STRegionIndex][nJetsBin] = new TH2I(("h_" + getHistogramName("totalNEvents", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("totalNEvents", STRegionIndex, nJetsBin, arguments, STRegions).c_str(), arguments.nGluinoMassBins, arguments.minGluinoMass, arguments.maxGluinoMass, arguments.nNeutralinoMassBins, arguments.minNeutralinoMass, arguments.maxNeutralinoMass);
      outputHistograms->h_totalNEvents[STRegionIndex][nJetsBin]->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson); // this is the main plot that will be used to estimate statistical MC uncertainties
      outputHistograms->h_lumiBasedYearWeightedNEvents[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("lumiBasedYearWeightedNEvents", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("lumiBasedYearWeightedNEvents", STRegionIndex, nJetsBin, arguments, STRegions).c_str(), arguments.nGluinoMassBins, arguments.minGluinoMass, arguments.maxGluinoMass, arguments.nNeutralinoMassBins, arguments.minNeutralinoMass, arguments.maxNeutralinoMass);
      outputHistograms->h_lumiBasedYearWeightedNEvents_prefiringDown[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("lumiBasedYearWeightedNEvents_prefiringDown", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("lumiBasedYearWeightedNEvents_prefiringDown", STRegionIndex, nJetsBin, arguments, STRegions).c_str(), arguments.nGluinoMassBins, arguments.minGluinoMass, arguments.maxGluinoMass, arguments.nNeutralinoMassBins, arguments.minNeutralinoMass, arguments.maxNeutralinoMass);
      outputHistograms->h_lumiBasedYearWeightedNEvents_prefiringUp[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("lumiBasedYearWeightedNEvents_prefiringUp", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("lumiBasedYearWeightedNEvents_prefiringUp", STRegionIndex, nJetsBin, arguments, STRegions).c_str(), arguments.nGluinoMassBins, arguments.minGluinoMass, arguments.maxGluinoMass, arguments.nNeutralinoMassBins, arguments.minNeutralinoMass, arguments.maxNeutralinoMass);
      outputHistograms->h_lumiBasedYearWeightedNEvents_photonScaleFactorDown[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("lumiBasedYearWeightedNEvents_photonScaleFactorDown", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("lumiBasedYearWeightedNEvents_photonScaleFactorDown", STRegionIndex, nJetsBin, arguments, STRegions).c_str(), arguments.nGluinoMassBins, arguments.minGluinoMass, arguments.maxGluinoMass, arguments.nNeutralinoMassBins, arguments.minNeutralinoMass, arguments.maxNeutralinoMass);
      outputHistograms->h_lumiBasedYearWeightedNEvents_photonScaleFactorUp[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("lumiBasedYearWeightedNEvents_photonScaleFactorUp", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("lumiBasedYearWeightedNEvents_photonScaleFactorUp", STRegionIndex, nJetsBin, arguments, STRegions).c_str(), arguments.nGluinoMassBins, arguments.minGluinoMass, arguments.maxGluinoMass, arguments.nNeutralinoMassBins, arguments.minNeutralinoMass, arguments.maxNeutralinoMass);
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

void fillOutputHistogramsFromFile(const std::string& fileName, outputHistogramsStruct *outputHistograms, argumentsStruct& arguments, std::map<int, double>& crossSections, const STRegionsStruct& STRegions, const double& relativeMCWeight, const double& integratedLuminosityTotal, const bool& isMain, TEfficiency* inputEfficiency_leading, TEfficiency* inputEfficiency_subLeading) {
  TChain *inputChain = new TChain("ggNtuplizer/EventTree");
  inputChain->Add(fileName.c_str());
  long nEntries = inputChain->GetEntries();
  std::cout << "Number of entries in " << fileName << ": " << nEntries << std::endl;
  TTreeReader inputTreeReader(inputChain);
  TTreeReaderValue<int> evt_nJets(inputTreeReader, "b_nJets");
  TTreeReaderValue<int> evt_nJets_shifted_JECDown(inputTreeReader, getShiftedVariableBranchName(shiftType::JECDown, "nJets").c_str());
  TTreeReaderValue<int> evt_nJets_shifted_JECUp(inputTreeReader, getShiftedVariableBranchName(shiftType::JECUp, "nJets").c_str());
  TTreeReaderValue<float> evt_ST(inputTreeReader, "b_evtST");
  TTreeReaderValue<float> evt_photonPT_leading(inputTreeReader, "b_photonPT_leading");
  TTreeReaderValue<float> evt_photonPT_subLeading(inputTreeReader, "b_photonPT_subLeading");
  TTreeReaderValue<float> evt_photonEta_leading(inputTreeReader, "b_photonEta_leading");
  TTreeReaderValue<float> evt_photonEta_subLeading(inputTreeReader, "b_photonEta_subLeading");
  TTreeReaderValue<float> evt_ST_shifted_JECDown(inputTreeReader, getShiftedVariableBranchName(shiftType::JECDown, "evtST").c_str());
  TTreeReaderValue<float> evt_ST_shifted_JECUp(inputTreeReader, getShiftedVariableBranchName(shiftType::JECUp, "evtST").c_str());
  TTreeReaderValue<float> evt_ST_shifted_UnclusteredMETDown(inputTreeReader, getShiftedVariableBranchName(shiftType::UnclusteredMETDown, "evtST").c_str());
  TTreeReaderValue<float> evt_ST_shifted_UnclusteredMETUp(inputTreeReader, getShiftedVariableBranchName(shiftType::UnclusteredMETUp, "evtST").c_str());
  TTreeReaderValue<float> evt_ST_shifted_JERMETDown(inputTreeReader, getShiftedVariableBranchName(shiftType::JERMETDown, "evtST").c_str());
  TTreeReaderValue<float> evt_ST_shifted_JERMETUp(inputTreeReader, getShiftedVariableBranchName(shiftType::JERMETUp, "evtST").c_str());
  TTreeReaderValue<float> prefiringWeight(inputTreeReader, "b_evtPrefiringWeight");
  TTreeReaderValue<float> prefiringWeightDown(inputTreeReader, "b_evtPrefiringWeightDown");
  TTreeReaderValue<float> prefiringWeightUp(inputTreeReader, "b_evtPrefiringWeightUp");
  TTreeReaderValue<float> photonMCScaleFactor(inputTreeReader, "b_evtphotonMCScaleFactor");
  TTreeReaderValue<float> photonMCScaleFactorDown(inputTreeReader, "b_evtphotonMCScaleFactorDown");
  TTreeReaderValue<float> photonMCScaleFactorUp(inputTreeReader, "b_evtphotonMCScaleFactorUp");
  TTreeReaderValue<int> nMC(inputTreeReader, "nMC");
  TTreeReaderArray<int> mcPIDs(inputTreeReader, "mcPID");
  TTreeReaderArray<float> mcMasses(inputTreeReader, "mcMass");
  TTreeReaderArray<int> mcMomPIDs(inputTreeReader, "mcMomPID");
  TTreeReaderArray<float> mcMomMasses(inputTreeReader, "mcMomMass");
  tmProgressBar *progressBar = new tmProgressBar(nEntries);
  int tmp = static_cast<int>(0.5 + 1.0*nEntries/1000);
  int progressBarUpdatePeriod = tmp > 1 ? tmp : 1;
  int entryIndex = 0;
  progressBar->initialize();
  while (inputTreeReader.Next()) {
    if (entryIndex >= nEntries) {
      std::cout << "ERROR: entry index seems to be greater than available number of events" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // sanity checks begin
    if (*nMC != static_cast<int>(mcPIDs.GetSize()) || *nMC != static_cast<int>(mcMomPIDs.GetSize())) {
      std::cout << "ERROR: Something has gone terribly wrong, number of MC particles is not the same as the size of the vectors containing MCPIDs" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    int maxNJetsEvt = std::max({*evt_nJets, *evt_nJets_shifted_JECUp, *evt_nJets_shifted_JECDown});
    if (maxNJetsEvt < 2) {
      std::cout << "ERROR: Fewer than 2 jets with all JECs in event, should not have passed selection..." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // sanity checks end
    if (entryIndex % progressBarUpdatePeriod == 0) progressBar->updateBar(1.0*entryIndex/nEntries, entryIndex);

    int STRegionIndex = (STRegions.STAxis).FindFixBin(*evt_ST);
    int STRegionIndex_shifted_JECDown = (STRegions.STAxis).FindFixBin(*evt_ST_shifted_JECDown);
    int STRegionIndex_shifted_JECUp = (STRegions.STAxis).FindFixBin(*evt_ST_shifted_JECUp);
    int STRegionIndex_shifted_UnclusteredMETDown = (STRegions.STAxis).FindFixBin(*evt_ST_shifted_UnclusteredMETDown);
    int STRegionIndex_shifted_UnclusteredMETUp = (STRegions.STAxis).FindFixBin(*evt_ST_shifted_UnclusteredMETUp);
    int STRegionIndex_shifted_JERMETDown = (STRegions.STAxis).FindFixBin(*evt_ST_shifted_JERMETDown);
    int STRegionIndex_shifted_JERMETUp = (STRegions.STAxis).FindFixBin(*evt_ST_shifted_JERMETUp);
    int maxRegionIndex = std::max({STRegionIndex, STRegionIndex_shifted_JECDown, STRegionIndex_shifted_JECUp, STRegionIndex_shifted_UnclusteredMETDown, STRegionIndex_shifted_UnclusteredMETUp, STRegionIndex_shifted_JERMETDown, STRegionIndex_shifted_JERMETUp});
    if (maxRegionIndex == 0) continue;
    if (maxRegionIndex > (1+STRegions.nSTSignalBins)) {
      std::cout << "WARNING: unexpected event ST: " << *evt_ST << std::endl;
      continue;
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

    // get HLT efficiencies
    float HLTEfficiency_leading = inputEfficiency_leading->GetEfficiency(inputEfficiency_leading->FindFixBin(*evt_photonEta_leading, *evt_photonPT_leading));
    float HLTEfficiency_subLeading = inputEfficiency_subLeading->GetEfficiency(inputEfficiency_subLeading->FindFixBin(*evt_photonEta_subLeading, *evt_photonPT_subLeading));
    float netHLTEfficiency = HLTEfficiency_leading*HLTEfficiency_subLeading;

    // get event weight
    int gMassInt = static_cast<int>(0.5 + generated_gluinoMass);
    double unscaledWeight = crossSections.at(gMassInt)*(integratedLuminosityTotal)/arguments.nGeneratedEventsPerBin;// factor of 4 to account for the fact that the MC production assumes a 50% branching ratio of neutralino to photon, while the analysis assumes 100% <-- 07 Jan 2019: factor of 4 removed until better understood...
    double nominalWeight = unscaledWeight*netHLTEfficiency*(*prefiringWeight)*(*photonMCScaleFactor);

    double tmp_gM = static_cast<double>(generated_gluinoMass);
    double tmp_nM = static_cast<double>(generated_neutralinoMass);

    if (isMain) { // To be filled only for MAIN sample
      if ((nJetsBin >= 2) && (STRegionIndex > 0)) outputHistograms->h_totalNEvents[STRegionIndex][nJetsBin]->Fill(tmp_gM, tmp_nM);
      if ((nJetsBin_JECDown >= 2) && (STRegionIndex_shifted_JECDown > 0)) outputHistograms->h_totalNEvents_shifted[shiftType::JECDown][STRegionIndex_shifted_JECDown][nJetsBin_JECDown]->Fill(tmp_gM, tmp_nM);
      if ((nJetsBin_JECUp >= 2) && (STRegionIndex_shifted_JECUp > 0)) outputHistograms->h_totalNEvents_shifted[shiftType::JECUp][STRegionIndex_shifted_JECUp][nJetsBin_JECUp]->Fill(tmp_gM, tmp_nM);
      if (nJetsBin >= 2) {
	if (STRegionIndex_shifted_UnclusteredMETDown > 0) outputHistograms->h_totalNEvents_shifted[shiftType::UnclusteredMETDown][STRegionIndex_shifted_UnclusteredMETDown][nJetsBin]->Fill(tmp_gM, tmp_nM);
	if (STRegionIndex_shifted_UnclusteredMETUp > 0) outputHistograms->h_totalNEvents_shifted[shiftType::UnclusteredMETUp][STRegionIndex_shifted_UnclusteredMETUp][nJetsBin]->Fill(tmp_gM, tmp_nM);
	if (STRegionIndex_shifted_JERMETDown > 0) outputHistograms->h_totalNEvents_shifted[shiftType::JERMETDown][STRegionIndex_shifted_JERMETDown][nJetsBin]->Fill(tmp_gM, tmp_nM);
	if (STRegionIndex_shifted_JERMETUp > 0) outputHistograms->h_totalNEvents_shifted[shiftType::JERMETUp][STRegionIndex_shifted_JERMETUp][nJetsBin]->Fill(tmp_gM, tmp_nM);
      }

      for (const auto& keyValuePair : arguments.specialZonesFor_sTDistributions) {
	int specialZoneIndex = keyValuePair.first;
	if (arguments.specialZonesFor_sTDistributions[specialZoneIndex].contains(generated_gluinoMass, generated_neutralinoMass)) {
	  if (nJetsBin >= 2) outputHistograms->h_sTDistributions[specialZoneIndex][nJetsBin]->Fill(*evt_ST, nominalWeight);
	  if (nJetsBin_JECDown >= 2) outputHistograms->h_sTDistributions_shifted[specialZoneIndex][shiftType::JECDown][nJetsBin_JECDown]->Fill(*evt_ST_shifted_JECDown, nominalWeight);
	  if (nJetsBin_JECUp >= 2) outputHistograms->h_sTDistributions_shifted[specialZoneIndex][shiftType::JECUp][nJetsBin_JECUp]->Fill(*evt_ST_shifted_JECUp, nominalWeight);
	  if (nJetsBin >= 2) {
	    outputHistograms->h_sTDistributions_shifted[specialZoneIndex][shiftType::UnclusteredMETDown][nJetsBin]->Fill(*evt_ST_shifted_UnclusteredMETDown, nominalWeight);
	    outputHistograms->h_sTDistributions_shifted[specialZoneIndex][shiftType::UnclusteredMETUp][nJetsBin]->Fill(*evt_ST_shifted_UnclusteredMETUp, nominalWeight);
	    outputHistograms->h_sTDistributions_shifted[specialZoneIndex][shiftType::JERMETDown][nJetsBin]->Fill(*evt_ST_shifted_JERMETDown, nominalWeight);
	    outputHistograms->h_sTDistributions_shifted[specialZoneIndex][shiftType::JERMETUp][nJetsBin]->Fill(*evt_ST_shifted_JERMETUp, nominalWeight);
	  }
	}
      }
    } // end to be filled only for MAIN sample

    if ((nJetsBin >= 2) && (STRegionIndex > 0)) { // to be filled in regardless of whether or not sample is the main sample
      outputHistograms->h_lumiBasedYearWeightedNEvents[STRegionIndex][nJetsBin]->Fill(tmp_gM, tmp_nM, unscaledWeight*netHLTEfficiency*relativeMCWeight*(*prefiringWeight)*(*photonMCScaleFactor));
      outputHistograms->h_lumiBasedYearWeightedNEvents_prefiringDown[STRegionIndex][nJetsBin]->Fill(tmp_gM, tmp_nM, unscaledWeight*netHLTEfficiency*relativeMCWeight*(*prefiringWeightDown)*(*photonMCScaleFactor));
      outputHistograms->h_lumiBasedYearWeightedNEvents_prefiringUp[STRegionIndex][nJetsBin]->Fill(tmp_gM, tmp_nM, unscaledWeight*netHLTEfficiency*relativeMCWeight*(*prefiringWeightUp)*(*photonMCScaleFactor));
      outputHistograms->h_lumiBasedYearWeightedNEvents_photonScaleFactorDown[STRegionIndex][nJetsBin]->Fill(tmp_gM, tmp_nM, unscaledWeight*netHLTEfficiency*relativeMCWeight*(*prefiringWeight)*(*photonMCScaleFactorDown));
      outputHistograms->h_lumiBasedYearWeightedNEvents_photonScaleFactorUp[STRegionIndex][nJetsBin]->Fill(tmp_gM, tmp_nM, unscaledWeight*netHLTEfficiency*relativeMCWeight*(*prefiringWeight)*(*photonMCScaleFactorUp));
    }
    ++entryIndex;
  }
  delete progressBar;
}

void fillOutputHistograms(outputHistogramsStruct *outputHistograms, argumentsStruct& arguments, std::map<int, double>& crossSections, const STRegionsStruct& STRegions, const double& integratedLuminosityTotal) {
  std::cout << "Filling events output histograms..." << std::endl;
  std::cout << "Fetching efficiency histograms..." << std::endl;
  std::vector<std::string> HLTEfficiencySourcesSplit = tmMiscUtils::getSplitString(arguments.HLTEfficiencySources, ":");
  if (!(HLTEfficiencySourcesSplit.size() == 3)) {
    std::cout << "ERROR: argument HLTEfficiencySources in unexpected format. Current value: " << arguments.HLTEfficiencySources << std::endl;
    std::exit(EXIT_FAILURE);
  }
  TFile *HLTEfficiencies_inputFile = TFile::Open((HLTEfficiencySourcesSplit.at(0)).c_str(),"READ");
  TEfficiency *inputEfficiency_leading;
  HLTEfficiencies_inputFile->GetObject((HLTEfficiencySourcesSplit.at(1)).c_str(), inputEfficiency_leading);
  if (inputEfficiency_leading) {
    std::cout << "Opened HLT efficiency: " << HLTEfficiencySourcesSplit.at(1) << std::endl;
  }
  else {
    std::cout << "ERROR: Unable to open: " << HLTEfficiencySourcesSplit.at(1) << std::endl;
    std::exit(EXIT_FAILURE);
  }
  TEfficiency *inputEfficiency_subLeading;
  HLTEfficiencies_inputFile->GetObject((HLTEfficiencySourcesSplit.at(2)).c_str(), inputEfficiency_subLeading);
  if (inputEfficiency_leading) {
    std::cout << "Opened HLT efficiency: " << HLTEfficiencySourcesSplit.at(2) << std::endl;
  }
  else {
    std::cout << "ERROR: Unable to open: " << HLTEfficiencySourcesSplit.at(2) << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::cout << "First for the primary MC sample..." << std::endl;
  fillOutputHistogramsFromFile(arguments.inputMCPathMain, outputHistograms, arguments, crossSections, STRegions, (arguments.integratedLuminosityMain)/integratedLuminosityTotal, integratedLuminosityTotal, true, inputEfficiency_leading, inputEfficiency_subLeading);
  if (!(arguments.inputMCPathsAux.empty())) {
    std::cout << "Now the aux samples..." << std::endl;
    for (const std::string& inputMCPathAux: arguments.inputMCPathsAux) {
      fillOutputHistogramsFromFile(inputMCPathAux, outputHistograms, arguments, crossSections, STRegions, ((arguments.integratedLuminositiesAux).at(inputMCPathAux))/integratedLuminosityTotal, integratedLuminosityTotal, false, inputEfficiency_leading, inputEfficiency_subLeading);
    }
  }
  HLTEfficiencies_inputFile->Close();
}

void saveHistograms(outputHistogramsStruct *outputHistograms, argumentsStruct& arguments, const STRegionsStruct& STRegions) {
  std::cout << "Saving histograms..." << std::endl;

  // First the 2D event histograms
  TFile *outputFile = TFile::Open((arguments.outputDirectory + "/" + arguments.outputPrefix + "_savedObjects.root").c_str(), "RECREATE");
  for (int STRegionIndex = 1; STRegionIndex <= (1+STRegions.nSTSignalBins); ++STRegionIndex) {
    for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
      std::string histogramName_total = getHistogramName("totalNEvents", STRegionIndex, nJetsBin);
      tmROOTSaverUtils::saveSingleObject(outputHistograms->h_totalNEvents[STRegionIndex][nJetsBin], "c_" + histogramName_total, outputFile, "", 1024, 768, 0, ".0f", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);
      std::string histogramName_weighted = getHistogramName("lumiBasedYearWeightedNEvents", STRegionIndex, nJetsBin);
      tmROOTSaverUtils::saveSingleObject(outputHistograms->h_lumiBasedYearWeightedNEvents[STRegionIndex][nJetsBin], "c_" + histogramName_weighted, outputFile, arguments.outputDirectory + "/" + arguments.outputPrefix + "_" + histogramName_weighted + ".png", 1024, 768, 0, ".1f", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);
      tmROOTSaverUtils::saveSingleObject(outputHistograms->h_lumiBasedYearWeightedNEvents_prefiringDown[STRegionIndex][nJetsBin], "c_" + getHistogramName("lumiBasedYearWeightedNEvents_prefiringDown", STRegionIndex, nJetsBin), outputFile, "", 1024, 768, 0, ".1f", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);
      tmROOTSaverUtils::saveSingleObject(outputHistograms->h_lumiBasedYearWeightedNEvents_prefiringUp[STRegionIndex][nJetsBin], "c_" + getHistogramName("lumiBasedYearWeightedNEvents_prefiringUp", STRegionIndex, nJetsBin), outputFile, "", 1024, 768, 0, ".1f", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);
      tmROOTSaverUtils::saveSingleObject(outputHistograms->h_lumiBasedYearWeightedNEvents_photonScaleFactorDown[STRegionIndex][nJetsBin], "c_" + getHistogramName("lumiBasedYearWeightedNEvents_photonScaleFactorDown", STRegionIndex, nJetsBin), outputFile, "", 1024, 768, 0, ".1f", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);
      tmROOTSaverUtils::saveSingleObject(outputHistograms->h_lumiBasedYearWeightedNEvents_photonScaleFactorUp[STRegionIndex][nJetsBin], "c_" + getHistogramName("lumiBasedYearWeightedNEvents_photonScaleFactorUp", STRegionIndex, nJetsBin), outputFile, "", 1024, 768, 0, ".1f", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);
      for (int shiftTypeIndex = shiftTypeFirst; shiftTypeIndex != static_cast<int>(shiftType::nShiftTypes); ++shiftTypeIndex) {
        shiftType typeIndex = static_cast<shiftType>(shiftTypeIndex);
        std::string histogramName_shifted = getHistogramName(typeIndex, "totalNEvents_shifted", STRegionIndex, nJetsBin);
        tmROOTSaverUtils::saveSingleObject(outputHistograms->h_totalNEvents_shifted[typeIndex][STRegionIndex][nJetsBin], "c_" + histogramName_shifted, outputFile, "", 1024, 768, 0, ".0f", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);
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
  argumentParser.addArgument("inputMCPathsAux", "", false, "Semicolon-separated list of paths to other MCs.");
  argumentParser.addArgument("integratedLuminositiesAux", "", false, "Semicolon-separated list of integrated luminosities in auxiliary MC reference.");
  // argumentParser.addArgument("maxMCEvents", "0", false, "Set a custom maximum number of MC events.");
  argumentParser.addArgument("crossSectionsFilePath", "SusyCrossSections13TevGluGlu.txt", false, "Path to dat file that contains cross-sections as a function of gluino mass, to use while weighting events.");
  argumentParser.addArgument("inputFile_STRegionBoundaries", "STRegionBoundaries.dat", false, "Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.");
  argumentParser.addArgument("sTMax_toPlot", "4000.0", false, "Max value of sT to plot.");
  argumentParser.addArgument("n_sTBinsToPlot", "29", false, "Number of sT bins to plot."); // default: 23 bins from 1200 to 3500 GeV in steps of 100 GeV
  argumentParser.addArgument("outputDirectory", "analysis/MCEventHistograms/", false, "Prefix to output files.");
  argumentParser.addArgument("outputPrefix", "", true, "Prefix to output files.");
  argumentParser.addArgument("nGeneratedEventsPerBin", "150000", false, "Number of generated events per bin, to use while calculating event weights.");
  argumentParser.addArgument("nGluinoMassBins", "16", false, "nBins on the gluino mass axis."); // (800 - 25) GeV --> (1750 + 25) GeV in steps of 50 GeV
  argumentParser.addArgument("minGluinoMass", "975.0", false, "Min gluino mass for the 2D plots.");
  argumentParser.addArgument("maxGluinoMass", "1775.0", false, "Max gluino mass for the 2D plots.");
  argumentParser.addArgument("nNeutralinoMassBins", "133", false, "nBins on the neutralino mass axis.");
  argumentParser.addArgument("minNeutralinoMass", "93.75", false, "Min neutralino mass for the 2D plots.");
  argumentParser.addArgument("maxNeutralinoMass", "1756.25", false, "Max neutralino mass for the 2D plots."); // (100 - 6.25) GeV --> (1750 + 6.25) GeV in steps of 12.5 GeV
  argumentParser.addArgument("regionsIn_sTHistograms", "1725.0:1775.0:650.0:950.0|1025.0:1075.0:975.0:1075.0", false, "List of the regions in which to fill and save the sT histograms. Each element of the list is in format minGluinoMass:maxGluinoMass:minNeutralinoMas:maxNeutralinoMass, and each element is separated from the next by the character \"|\". See also the default value for this argument.");
  argumentParser.addArgument("HLTEfficiencySources", "", true, "Location of the TEfficiency objects from which to read in the HLT efficiency values for the leading and subleading photons. Format: \"A:B:C\", where A = path of ROOT file containing histograms, B = leading photon HLT efficiency, C = subleading photon HLT efficiency.");
  argumentParser.setPassedStringValues(argc, argv);
  argumentsStruct arguments = getArgumentsFromParser(argumentParser);
  double integratedLuminosityTotal = arguments.integratedLuminosityMain;
  if (!(arguments.inputMCPathsAux.empty())) {
    for (const auto& MCPathAux: arguments.inputMCPathsAux) {
      integratedLuminosityTotal += (arguments.integratedLuminositiesAux).at(MCPathAux);
    }
  }

  STRegionsStruct STRegions(arguments.inputFile_STRegionBoundaries);
  outputHistogramsStruct* outputHistograms = initializeOutputHistograms(arguments, STRegions);
  std::map<int, double> crossSections = getCrossSectionsFromFile(arguments.crossSectionsFilePath);
  fillOutputHistograms(outputHistograms, arguments, crossSections, STRegions, integratedLuminosityTotal);
  saveHistograms(outputHistograms, arguments, STRegions);
  return 0;
}
