#include "../include/getEventHistograms.h"

outputHistogramsStruct* initializeOutputHistograms(argumentsStruct& arguments, MCTemplateReader& templateReader, const STRegionsStruct& STRegions) {
  outputHistogramsStruct* outputHistograms = new outputHistogramsStruct();

  // 1D histograms for event weights
  outputHistograms->h_overallWeights = new TH1F("h_overallWeights", "overall weights;weight;events/0.01", 201, -0.005, 2.005);
  outputHistograms->h_hltEfficiency = new TH1F("h_hltEfficiency", "HLT efficiency;weight;events/0.01", 201, -0.005, 2.005);
  outputHistograms->h_pileupWeights = new TH1F("h_pileupWeights", "pileup weights;weight;events/0.1", 101, -0.05, 10.05);
  outputHistograms->h_prefiringWeights = new TH1F("h_prefiringWeights", "prefiring weights;weight;events/0.01", 201, -0.005, 2.005);
  outputHistograms->h_photonMCScaleFactorWeights = new TH1F("h_photonMCScaleFactorWeights", "photon MC scale factor weights;weight;events/0.01", 201, -0.005, 2.005);

  // 2D histograms for nEvents
  for (int STRegionIndex = 1; STRegionIndex <= (1+STRegions.nSTSignalBins); ++STRegionIndex) {
    for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
      outputHistograms->h_totalNEvents[STRegionIndex][nJetsBin] = new TH2I(("h_" + getHistogramName("totalNEvents", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("totalNEvents", STRegionIndex, nJetsBin, arguments, STRegions).c_str(), templateReader.nEventProgenitorMassBins, templateReader.minEventProgenitorMass, templateReader.maxEventProgenitorMass, templateReader.nNeutralinoMassBins, templateReader.minNeutralinoMass, templateReader.maxNeutralinoMass);
      outputHistograms->h_totalNEvents[STRegionIndex][nJetsBin]->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson); // this is the main plot that will be used to estimate statistical MC uncertainties
      outputHistograms->h_lumiBasedYearWeightedNEvents[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("lumiBasedYearWeightedNEvents", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("lumiBasedYearWeightedNEvents", STRegionIndex, nJetsBin, arguments, STRegions).c_str(), templateReader.nEventProgenitorMassBins, templateReader.minEventProgenitorMass, templateReader.maxEventProgenitorMass, templateReader.nNeutralinoMassBins, templateReader.minNeutralinoMass, templateReader.maxNeutralinoMass);
      outputHistograms->h_lumiBasedYearWeightedNEvents_prefiringDown[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("lumiBasedYearWeightedNEvents_prefiringDown", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("lumiBasedYearWeightedNEvents_prefiringDown", STRegionIndex, nJetsBin, arguments, STRegions).c_str(), templateReader.nEventProgenitorMassBins, templateReader.minEventProgenitorMass, templateReader.maxEventProgenitorMass, templateReader.nNeutralinoMassBins, templateReader.minNeutralinoMass, templateReader.maxNeutralinoMass);
      outputHistograms->h_lumiBasedYearWeightedNEvents_prefiringUp[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("lumiBasedYearWeightedNEvents_prefiringUp", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("lumiBasedYearWeightedNEvents_prefiringUp", STRegionIndex, nJetsBin, arguments, STRegions).c_str(), templateReader.nEventProgenitorMassBins, templateReader.minEventProgenitorMass, templateReader.maxEventProgenitorMass, templateReader.nNeutralinoMassBins, templateReader.minNeutralinoMass, templateReader.maxNeutralinoMass);
      outputHistograms->h_lumiBasedYearWeightedNEvents_HLTDown[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("lumiBasedYearWeightedNEvents_HLTDown", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("lumiBasedYearWeightedNEvents_HLTDown", STRegionIndex, nJetsBin, arguments, STRegions).c_str(), templateReader.nEventProgenitorMassBins, templateReader.minEventProgenitorMass, templateReader.maxEventProgenitorMass, templateReader.nNeutralinoMassBins, templateReader.minNeutralinoMass, templateReader.maxNeutralinoMass);
      outputHistograms->h_lumiBasedYearWeightedNEvents_HLTUp[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("lumiBasedYearWeightedNEvents_HLTUp", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("lumiBasedYearWeightedNEvents_HLTUp", STRegionIndex, nJetsBin, arguments, STRegions).c_str(), templateReader.nEventProgenitorMassBins, templateReader.minEventProgenitorMass, templateReader.maxEventProgenitorMass, templateReader.nNeutralinoMassBins, templateReader.minNeutralinoMass, templateReader.maxNeutralinoMass);
      outputHistograms->h_lumiBasedYearWeightedNEvents_photonScaleFactorDown[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("lumiBasedYearWeightedNEvents_photonScaleFactorDown", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("lumiBasedYearWeightedNEvents_photonScaleFactorDown", STRegionIndex, nJetsBin, arguments, STRegions).c_str(), templateReader.nEventProgenitorMassBins, templateReader.minEventProgenitorMass, templateReader.maxEventProgenitorMass, templateReader.nNeutralinoMassBins, templateReader.minNeutralinoMass, templateReader.maxNeutralinoMass);
      outputHistograms->h_lumiBasedYearWeightedNEvents_photonScaleFactorUp[STRegionIndex][nJetsBin] = new TH2F(("h_" + getHistogramName("lumiBasedYearWeightedNEvents_photonScaleFactorUp", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle("lumiBasedYearWeightedNEvents_photonScaleFactorUp", STRegionIndex, nJetsBin, arguments, STRegions).c_str(), templateReader.nEventProgenitorMassBins, templateReader.minEventProgenitorMass, templateReader.maxEventProgenitorMass, templateReader.nNeutralinoMassBins, templateReader.minNeutralinoMass, templateReader.maxNeutralinoMass);
      for (int shiftTypeIndex = shiftTypeFirst; shiftTypeIndex != static_cast<int>(shiftType::nShiftTypes); ++shiftTypeIndex) {
        shiftType typeIndex = static_cast<shiftType>(shiftTypeIndex);
        outputHistograms->h_totalNEvents_shifted[typeIndex][STRegionIndex][nJetsBin] = new TH2I(("h_" + getHistogramName(typeIndex, "totalNEvents_shifted", STRegionIndex, nJetsBin)).c_str(), getHistogramTitle(typeIndex, "totalNEvents_shifted", STRegionIndex, nJetsBin, arguments, STRegions).c_str(), templateReader.nEventProgenitorMassBins, templateReader.minEventProgenitorMass, templateReader.maxEventProgenitorMass, templateReader.nNeutralinoMassBins, templateReader.minNeutralinoMass, templateReader.maxNeutralinoMass);
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
  int eventProgenitorMass;
  double crossSection, crossSectionPercentUncertainty;
  while (crossSectionsFileStream >> eventProgenitorMass >> crossSection >> crossSectionPercentUncertainty) {
    crossSections[eventProgenitorMass] = crossSection;
    // crossSectionsFractionalUncertainty[eventProgenitorMass] = 0.01*crossSectionPercentUncertainty;
    if (printCrossSections) {
      std::cout << "crossSections[" << eventProgenitorMass << "] = " << crossSection << ", "
                << "crossSectionsFractionalUncertainty[" << eventProgenitorMass << "] = " << 0.01*crossSectionPercentUncertainty << std::endl;
    }
  }
  return crossSections;
}

void fillOutputHistogramsFromFile(const std::string& fileName, outputHistogramsStruct *outputHistograms, argumentsStruct& arguments, std::map<int, double>& crossSections, MCTemplateReader& templateReader, const STRegionsStruct& STRegions, const double& relativeMCWeight, const double& integratedLuminosityTotal, const std::string& reweightingWeightsSourceFilePath, const std::string& HLTEfficienciesSourcePath, const bool& isMain) {
  TFile* reweightingWeightsSourceFile = TFile::Open(reweightingWeightsSourceFilePath.c_str(), "READ");
  TH1D* PUReweightingHistogram;
  reweightingWeightsSourceFile->GetObject("pileupWeights", PUReweightingHistogram);
  assert(PUReweightingHistogram != nullptr);
  TFile* HLTEfficienciesSourceFile = TFile::Open(HLTEfficienciesSourcePath.c_str(), "READ");
  TEfficiency* HLTEfficiency_clean;
  HLTEfficienciesSourceFile->GetObject("hltEfficiency_clean", HLTEfficiency_clean);
  assert(HLTEfficiency_clean != nullptr);
  // TEfficiency* HLTEfficiency_unclean;
  // HLTEfficienciesSourceFile->GetObject("hltEfficiency_unclean", HLTEfficiency_unclean);
  // assert(HLTEfficiency_unclean != nullptr);
  TChain *inputChain = new TChain("ggNtuplizer/EventTree");
  inputChain->Add(fileName.c_str());
  long nEntries = inputChain->GetEntries();
  TAxis eventProgenitorAxis(templateReader.nEventProgenitorMassBins, templateReader.minEventProgenitorMass, templateReader.maxEventProgenitorMass);
  TAxis neutralinoAxis(templateReader.nNeutralinoMassBins, templateReader.minNeutralinoMass, templateReader.maxNeutralinoMass);
  std::cout << "Number of entries in " << fileName << ": " << nEntries << std::endl;
  TTreeReader inputTreeReader(inputChain);
  TTreeReaderArray<int> evt_BX_for_PU(inputTreeReader, "puBX");
  TTreeReaderArray<float> evt_PU(inputTreeReader, "puTrue");
  TTreeReaderValue<int> evt_nJets(inputTreeReader, "b_nJetsDR");
  TTreeReaderValue<int> evt_nJets_shifted_JECDown(inputTreeReader, getShiftedVariableBranchName(shiftType::JECDown, "nJetsDR").c_str());
  TTreeReaderValue<int> evt_nJets_shifted_JECUp(inputTreeReader, getShiftedVariableBranchName(shiftType::JECUp, "nJetsDR").c_str());
  TTreeReaderValue<int> evt_nJets_shifted_missingHEMDown(inputTreeReader, getShiftedVariableBranchName(shiftType::missingHEMDown, "nJetsDR").c_str());
  TTreeReaderValue<int> evt_nJets_shifted_missingHEMUp(inputTreeReader, getShiftedVariableBranchName(shiftType::missingHEMUp, "nJetsDR").c_str());
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
  TTreeReaderValue<float> evt_ST_shifted_missingHEMDown(inputTreeReader, getShiftedVariableBranchName(shiftType::missingHEMDown, "evtST").c_str());
  TTreeReaderValue<float> evt_ST_shifted_missingHEMUp(inputTreeReader, getShiftedVariableBranchName(shiftType::missingHEMUp, "evtST").c_str());
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
  int tmp = static_cast<int>(0.5 + 1.0*nEntries/50);
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

    int maxNJetsEvt = std::max({*evt_nJets, *evt_nJets_shifted_JECUp, *evt_nJets_shifted_JECDown, *evt_nJets_shifted_missingHEMUp, *evt_nJets_shifted_missingHEMDown});
    if (maxNJetsEvt < 2) {
      std::cout << "ERROR: Fewer than 2 jets with all corrections in event, should not have passed selection..." << std::endl;
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
    int STRegionIndex_shifted_missingHEMDown = (STRegions.STAxis).FindFixBin(*evt_ST_shifted_missingHEMDown);
    int STRegionIndex_shifted_missingHEMUp = (STRegions.STAxis).FindFixBin(*evt_ST_shifted_missingHEMUp);
    int maxRegionIndex = std::max({STRegionIndex, STRegionIndex_shifted_JECDown, STRegionIndex_shifted_JECUp, STRegionIndex_shifted_UnclusteredMETDown, STRegionIndex_shifted_UnclusteredMETUp, STRegionIndex_shifted_JERMETDown, STRegionIndex_shifted_JERMETUp, STRegionIndex_shifted_missingHEMDown, STRegionIndex_shifted_missingHEMUp});
    if (maxRegionIndex == 0) {
      ++entryIndex;
      continue;
    }
    if (maxRegionIndex > (1+STRegions.nSTSignalBins)) {
      std::cout << "WARNING: unexpected event ST: " << *evt_ST << std::endl;
      ++entryIndex;
      continue;
    }

    int nJetsBin = *evt_nJets;
    if (*evt_nJets > 6) nJetsBin = 6;
    int nJetsBin_JECDown = *evt_nJets_shifted_JECDown;
    if (*evt_nJets_shifted_JECDown > 6) nJetsBin_JECDown = 6;
    int nJetsBin_JECUp = *evt_nJets_shifted_JECUp;
    if (*evt_nJets_shifted_JECUp > 6) nJetsBin_JECUp = 6;
    int nJetsBin_missingHEMDown = *evt_nJets_shifted_missingHEMDown;
    if (*evt_nJets_shifted_missingHEMDown > 6) nJetsBin_missingHEMDown = 6;
    int nJetsBin_missingHEMUp = *evt_nJets_shifted_missingHEMUp;
    if (*evt_nJets_shifted_missingHEMUp > 6) nJetsBin_missingHEMUp = 6;

    // get HLT efficiencies
    float adjusted_leading_photonPT = *evt_photonPT_leading;
    if (adjusted_leading_photonPT > (HLTEmulation::pTBinEdges).at((HLTEmulation::pTBinEdges).size()-1)) adjusted_leading_photonPT = (HLTEmulation::pTBinEdges).at((HLTEmulation::pTBinEdges).size()-1);
    else if (adjusted_leading_photonPT < (HLTEmulation::pTBinEdges).at(0)) adjusted_leading_photonPT = (HLTEmulation::pTBinEdges).at(0);
    float eventHLTEfficiency_clean = HLTEfficiency_clean->GetEfficiency(HLTEfficiency_clean->FindFixBin(adjusted_leading_photonPT));
    assert(eventHLTEfficiency_clean >= 0.);
    assert(eventHLTEfficiency_clean <= 1.);
    float eventHLTEfficiencyErrorDown = HLTEfficiency_clean->GetEfficiencyErrorLow(HLTEfficiency_clean->FindFixBin(adjusted_leading_photonPT));
    float eventHLTEfficiency_down = eventHLTEfficiency_clean - eventHLTEfficiencyErrorDown;
    if (eventHLTEfficiency_down < 0.) eventHLTEfficiency_down = 0.;
    float eventHLTEfficiencyErrorUp = HLTEfficiency_clean->GetEfficiencyErrorUp(HLTEfficiency_clean->FindFixBin(adjusted_leading_photonPT));
    float eventHLTEfficiency_up = eventHLTEfficiency_clean + eventHLTEfficiencyErrorUp;
    if (eventHLTEfficiency_up > 1.) eventHLTEfficiency_up = 1.;
    // float eventHLTEfficiency_unclean = HLTEfficiency_unclean->GetEfficiency(HLTEfficiency_unclean->FindFixBin(adjusted_leading_photonPT));
    // (void)eventHLTEfficiency_unclean;

    // get generated eventProgenitor, neutralino mass
    float generated_eventProgenitorMass = 0;
    bool eventProgenitorMassIsSet = false;
    float generated_neutralinoMass = 0;
    bool neutralinoMassIsSet = false;
    for (int generatedParticleIndex = 0; generatedParticleIndex < *nMC; ++generatedParticleIndex) {
      int particle_mcPID = mcPIDs[generatedParticleIndex];
      bool isCorrectProgenitor = false;
      if (arguments.eventProgenitor == "squark") isCorrectProgenitor = PIDUtils::isSquarkPID(particle_mcPID);
      else if (arguments.eventProgenitor == "gluino") isCorrectProgenitor = PIDUtils::isGluinoPID(particle_mcPID);
      if (!(eventProgenitorMassIsSet) && isCorrectProgenitor) {
        generated_eventProgenitorMass = mcMasses[generatedParticleIndex];
        eventProgenitorMassIsSet = true;
      }
      int particle_mcMomPID = mcMomPIDs[generatedParticleIndex];
      if (!(neutralinoMassIsSet) && PIDUtils::isNeutralinoPID(particle_mcMomPID)) {
        generated_neutralinoMass = mcMomMasses[generatedParticleIndex];
        neutralinoMassIsSet = true;
      }
    }
    if (!(eventProgenitorMassIsSet && neutralinoMassIsSet)) {
      std::cout << "ERROR: Unable to find eventProgenitor or neutralino mass in an event." << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // get event weight
    int eventProgenitorMassBin = eventProgenitorAxis.FindFixBin(generated_eventProgenitorMass);
    int neutralinoMassBin = neutralinoAxis.FindFixBin(generated_neutralinoMass);
    if (!(templateReader.isValidBin(eventProgenitorMassBin, neutralinoMassBin))) {
      ++entryIndex;
      continue;
    }
    int gMassInt = static_cast<int>(0.5 + generated_eventProgenitorMass);
    float eventPU = -1.;
    for (unsigned int BXCounter = 0; BXCounter < evt_BX_for_PU.GetSize(); ++BXCounter) {
      int bx = evt_BX_for_PU[BXCounter];
      if (bx == 0) {
	eventPU = evt_PU[BXCounter];
	break;
      }
    }
    assert(eventPU > 0.);
    double pileupWeight = PUReweightingHistogram->GetBinContent(PUReweightingHistogram->GetXaxis()->FindFixBin(eventPU));
    double unscaledWeight = crossSections.at(gMassInt)*(integratedLuminosityTotal)/(templateReader.getTotalNEvents(eventProgenitorMassBin, neutralinoMassBin));// factor of 4 to account for the fact that the MC production assumes a 50% branching ratio of neutralino to photon, while the analysis assumes 100% <-- 07 Jan 2019: factor of 4 removed until better understood...
    double nominalWeight = unscaledWeight*eventHLTEfficiency_clean*pileupWeight*(*prefiringWeight)*(*photonMCScaleFactor);

    double tmp_gM = static_cast<double>(generated_eventProgenitorMass);
    double tmp_nM = static_cast<double>(generated_neutralinoMass);

    if (isMain) { // To be filled only for MAIN sample
      if ((nJetsBin >= 2) && (STRegionIndex > 0)) outputHistograms->h_totalNEvents[STRegionIndex][nJetsBin]->Fill(tmp_gM, tmp_nM);
      if ((nJetsBin_JECDown >= 2) && (STRegionIndex_shifted_JECDown > 0)) outputHistograms->h_totalNEvents_shifted[shiftType::JECDown][STRegionIndex_shifted_JECDown][nJetsBin_JECDown]->Fill(tmp_gM, tmp_nM);
      if ((nJetsBin_JECUp >= 2) && (STRegionIndex_shifted_JECUp > 0)) outputHistograms->h_totalNEvents_shifted[shiftType::JECUp][STRegionIndex_shifted_JECUp][nJetsBin_JECUp]->Fill(tmp_gM, tmp_nM);
      if ((nJetsBin_missingHEMDown >= 2) && (STRegionIndex_shifted_missingHEMDown > 0)) outputHistograms->h_totalNEvents_shifted[shiftType::missingHEMDown][STRegionIndex_shifted_missingHEMDown][nJetsBin_missingHEMDown]->Fill(tmp_gM, tmp_nM);
      if ((nJetsBin_missingHEMUp >= 2) && (STRegionIndex_shifted_missingHEMUp > 0)) outputHistograms->h_totalNEvents_shifted[shiftType::missingHEMUp][STRegionIndex_shifted_missingHEMUp][nJetsBin_missingHEMUp]->Fill(tmp_gM, tmp_nM);
      if (nJetsBin >= 2) {
	if (STRegionIndex_shifted_UnclusteredMETDown > 0) outputHistograms->h_totalNEvents_shifted[shiftType::UnclusteredMETDown][STRegionIndex_shifted_UnclusteredMETDown][nJetsBin]->Fill(tmp_gM, tmp_nM);
	if (STRegionIndex_shifted_UnclusteredMETUp > 0) outputHistograms->h_totalNEvents_shifted[shiftType::UnclusteredMETUp][STRegionIndex_shifted_UnclusteredMETUp][nJetsBin]->Fill(tmp_gM, tmp_nM);
	if (STRegionIndex_shifted_JERMETDown > 0) outputHistograms->h_totalNEvents_shifted[shiftType::JERMETDown][STRegionIndex_shifted_JERMETDown][nJetsBin]->Fill(tmp_gM, tmp_nM);
	if (STRegionIndex_shifted_JERMETUp > 0) outputHistograms->h_totalNEvents_shifted[shiftType::JERMETUp][STRegionIndex_shifted_JERMETUp][nJetsBin]->Fill(tmp_gM, tmp_nM);
      }

      for (const auto& keyValuePair : arguments.specialZonesFor_sTDistributions) {
	int specialZoneIndex = keyValuePair.first;
	if (arguments.specialZonesFor_sTDistributions[specialZoneIndex].contains(generated_eventProgenitorMass, generated_neutralinoMass)) {
	  if (nJetsBin >= 2) outputHistograms->h_sTDistributions[specialZoneIndex][nJetsBin]->Fill(*evt_ST, nominalWeight);
	  if (nJetsBin_JECDown >= 2) outputHistograms->h_sTDistributions_shifted[specialZoneIndex][shiftType::JECDown][nJetsBin_JECDown]->Fill(*evt_ST_shifted_JECDown, nominalWeight);
	  if (nJetsBin_JECUp >= 2) outputHistograms->h_sTDistributions_shifted[specialZoneIndex][shiftType::JECUp][nJetsBin_JECUp]->Fill(*evt_ST_shifted_JECUp, nominalWeight);
	  if (nJetsBin_missingHEMDown >= 2) outputHistograms->h_sTDistributions_shifted[specialZoneIndex][shiftType::missingHEMDown][nJetsBin_missingHEMDown]->Fill(*evt_ST_shifted_missingHEMDown, nominalWeight);
	  if (nJetsBin_missingHEMUp >= 2) outputHistograms->h_sTDistributions_shifted[specialZoneIndex][shiftType::missingHEMUp][nJetsBin_missingHEMUp]->Fill(*evt_ST_shifted_missingHEMUp, nominalWeight);
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
      outputHistograms->h_overallWeights->Fill(eventHLTEfficiency_clean*pileupWeight*(*prefiringWeight)*(*photonMCScaleFactor), relativeMCWeight);
      outputHistograms->h_hltEfficiency->Fill(eventHLTEfficiency_clean, relativeMCWeight);
      outputHistograms->h_pileupWeights->Fill(pileupWeight, relativeMCWeight);
      outputHistograms->h_prefiringWeights->Fill((*prefiringWeight), relativeMCWeight);
      outputHistograms->h_photonMCScaleFactorWeights->Fill((*photonMCScaleFactor), relativeMCWeight);

      outputHistograms->h_lumiBasedYearWeightedNEvents[STRegionIndex][nJetsBin]->Fill(tmp_gM, tmp_nM, unscaledWeight*eventHLTEfficiency_clean*pileupWeight*relativeMCWeight*(*prefiringWeight)*(*photonMCScaleFactor));
      outputHistograms->h_lumiBasedYearWeightedNEvents_prefiringDown[STRegionIndex][nJetsBin]->Fill(tmp_gM, tmp_nM, unscaledWeight*eventHLTEfficiency_clean*pileupWeight*relativeMCWeight*(*prefiringWeightDown)*(*photonMCScaleFactor));
      outputHistograms->h_lumiBasedYearWeightedNEvents_prefiringUp[STRegionIndex][nJetsBin]->Fill(tmp_gM, tmp_nM, unscaledWeight*eventHLTEfficiency_clean*pileupWeight*relativeMCWeight*(*prefiringWeightUp)*(*photonMCScaleFactor));
      outputHistograms->h_lumiBasedYearWeightedNEvents_HLTDown[STRegionIndex][nJetsBin]->Fill(tmp_gM, tmp_nM, unscaledWeight*eventHLTEfficiency_down*pileupWeight*relativeMCWeight*(*prefiringWeight)*(*photonMCScaleFactor));
      outputHistograms->h_lumiBasedYearWeightedNEvents_HLTUp[STRegionIndex][nJetsBin]->Fill(tmp_gM, tmp_nM, unscaledWeight*eventHLTEfficiency_up*pileupWeight*relativeMCWeight*(*prefiringWeight)*(*photonMCScaleFactor));
      outputHistograms->h_lumiBasedYearWeightedNEvents_photonScaleFactorDown[STRegionIndex][nJetsBin]->Fill(tmp_gM, tmp_nM, unscaledWeight*eventHLTEfficiency_clean*pileupWeight*relativeMCWeight*(*prefiringWeight)*(*photonMCScaleFactorDown));
      outputHistograms->h_lumiBasedYearWeightedNEvents_photonScaleFactorUp[STRegionIndex][nJetsBin]->Fill(tmp_gM, tmp_nM, unscaledWeight*eventHLTEfficiency_clean*pileupWeight*relativeMCWeight*(*prefiringWeight)*(*photonMCScaleFactorUp));
    }
    ++entryIndex;
  }
  delete progressBar;
  reweightingWeightsSourceFile->Close();
}

void fillOutputHistograms(outputHistogramsStruct *outputHistograms, argumentsStruct& arguments, std::map<int, double>& crossSections, MCTemplateReader& templateReader, const STRegionsStruct& STRegions, const double& integratedLuminosityTotal) {
  std::cout << "Filling events output histograms..." << std::endl;
  std::cout << "First for the primary MC sample..." << std::endl;
  fillOutputHistogramsFromFile(arguments.inputMCPathMain, outputHistograms, arguments, crossSections, templateReader, STRegions, (arguments.integratedLuminosityMain)/integratedLuminosityTotal, integratedLuminosityTotal, arguments.inputPUWeightsPathMain, arguments.inputHLTEfficienciesPathMain, true);

  if (!(arguments.inputMCPathsAux.empty())) {
    std::cout << "Now the aux samples..." << std::endl;
    for (unsigned int auxIndex = 0; auxIndex < static_cast<unsigned int>((arguments.inputMCPathsAux).size()); ++auxIndex) {
      std::string inputMCPathAux = (arguments.inputMCPathsAux).at(auxIndex);
      fillOutputHistogramsFromFile(inputMCPathAux, outputHistograms, arguments, crossSections, templateReader, STRegions, ((arguments.integratedLuminositiesAux).at(auxIndex))/integratedLuminosityTotal, integratedLuminosityTotal, (arguments.inputPUWeightsPathsAux).at(auxIndex), (arguments.inputHLTEfficienciesPathsAux).at(auxIndex), false);
    }
  }
}

void saveHistograms(outputHistogramsStruct *outputHistograms, argumentsStruct& arguments, const STRegionsStruct& STRegions) {
  std::cout << "Saving histograms..." << std::endl;
  TFile *outputFile = TFile::Open((arguments.outputDirectory + "/" + arguments.outputPrefix + "_savedObjects.root").c_str(), "RECREATE");

  // First the weights histograms
  tmROOTSaverUtils::saveSingleObject(outputHistograms->h_overallWeights, "c_overallWeights", outputFile, "", 1024, 768, 110110, "", "", false, false, false, 0, 0, 0, 0, 0, 0);
  tmROOTSaverUtils::saveSingleObject(outputHistograms->h_hltEfficiency, "c_hltEfficiency", outputFile, "", 1024, 768, 110110, "", "", false, false, false, 0, 0, 0, 0, 0, 0);
  tmROOTSaverUtils::saveSingleObject(outputHistograms->h_pileupWeights, "c_pileupWeights", outputFile, "", 1024, 768, 110110, "", "", false, false, false, 0, 0, 0, 0, 0, 0);
  tmROOTSaverUtils::saveSingleObject(outputHistograms->h_prefiringWeights, "c_prefiringWeights", outputFile, "", 1024, 768, 110110, "", "", false, false, false, 0, 0, 0, 0, 0, 0);
  tmROOTSaverUtils::saveSingleObject(outputHistograms->h_photonMCScaleFactorWeights, "c_photonMCScaleFactorWeights", outputFile, "", 1024, 768, 110110, "", "", false, false, false, 0, 0, 0, 0, 0, 0);

  // Next the 2D event histograms
  for (int STRegionIndex = 1; STRegionIndex <= (1+STRegions.nSTSignalBins); ++STRegionIndex) {
    for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
      std::string histogramName_total = getHistogramName("totalNEvents", STRegionIndex, nJetsBin);
      tmROOTSaverUtils::saveSingleObject(outputHistograms->h_totalNEvents[STRegionIndex][nJetsBin], "c_" + histogramName_total, outputFile, "", 1024, 768, 0, ".0f", "COLZ TEXT25", false, false, true, 0, 0, 0, 0, 0, 0);
      std::string histogramName_weighted = getHistogramName("lumiBasedYearWeightedNEvents", STRegionIndex, nJetsBin);
      tmROOTSaverUtils::saveSingleObject(outputHistograms->h_lumiBasedYearWeightedNEvents[STRegionIndex][nJetsBin], "c_" + histogramName_weighted, outputFile, arguments.outputDirectory + "/" + arguments.outputPrefix + "_" + histogramName_weighted + ".pdf", 1024, 768, 0, ".1f", "COLZ TEXT25", false, false, true, 0, 0, 0, 0, 0, 0);
      tmROOTSaverUtils::saveSingleObject(outputHistograms->h_lumiBasedYearWeightedNEvents_prefiringDown[STRegionIndex][nJetsBin], "c_" + getHistogramName("lumiBasedYearWeightedNEvents_prefiringDown", STRegionIndex, nJetsBin), outputFile, "", 1024, 768, 0, ".1f", "COLZ TEXT25", false, false, true, 0, 0, 0, 0, 0, 0);
      tmROOTSaverUtils::saveSingleObject(outputHistograms->h_lumiBasedYearWeightedNEvents_prefiringUp[STRegionIndex][nJetsBin], "c_" + getHistogramName("lumiBasedYearWeightedNEvents_prefiringUp", STRegionIndex, nJetsBin), outputFile, "", 1024, 768, 0, ".1f", "COLZ TEXT25", false, false, true, 0, 0, 0, 0, 0, 0);
      tmROOTSaverUtils::saveSingleObject(outputHistograms->h_lumiBasedYearWeightedNEvents_HLTDown[STRegionIndex][nJetsBin], "c_" + getHistogramName("lumiBasedYearWeightedNEvents_HLTDown", STRegionIndex, nJetsBin), outputFile, "", 1024, 768, 0, ".1f", "COLZ TEXT25", false, false, true, 0, 0, 0, 0, 0, 0);
      tmROOTSaverUtils::saveSingleObject(outputHistograms->h_lumiBasedYearWeightedNEvents_HLTUp[STRegionIndex][nJetsBin], "c_" + getHistogramName("lumiBasedYearWeightedNEvents_HLTUp", STRegionIndex, nJetsBin), outputFile, "", 1024, 768, 0, ".1f", "COLZ TEXT25", false, false, true, 0, 0, 0, 0, 0, 0);
      tmROOTSaverUtils::saveSingleObject(outputHistograms->h_lumiBasedYearWeightedNEvents_photonScaleFactorDown[STRegionIndex][nJetsBin], "c_" + getHistogramName("lumiBasedYearWeightedNEvents_photonScaleFactorDown", STRegionIndex, nJetsBin), outputFile, "", 1024, 768, 0, ".1f", "COLZ TEXT25", false, false, true, 0, 0, 0, 0, 0, 0);
      tmROOTSaverUtils::saveSingleObject(outputHistograms->h_lumiBasedYearWeightedNEvents_photonScaleFactorUp[STRegionIndex][nJetsBin], "c_" + getHistogramName("lumiBasedYearWeightedNEvents_photonScaleFactorUp", STRegionIndex, nJetsBin), outputFile, "", 1024, 768, 0, ".1f", "COLZ TEXT25", false, false, true, 0, 0, 0, 0, 0, 0);
      for (int shiftTypeIndex = shiftTypeFirst; shiftTypeIndex != static_cast<int>(shiftType::nShiftTypes); ++shiftTypeIndex) {
        shiftType typeIndex = static_cast<shiftType>(shiftTypeIndex);
        std::string histogramName_shifted = getHistogramName(typeIndex, "totalNEvents_shifted", STRegionIndex, nJetsBin);
        tmROOTSaverUtils::saveSingleObject(outputHistograms->h_totalNEvents_shifted[typeIndex][STRegionIndex][nJetsBin], "c_" + histogramName_shifted, outputFile, "", 1024, 768, 0, ".0f", "COLZ TEXT25", false, false, true, 0, 0, 0, 0, 0, 0);
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
      tmROOTSaverUtils::saveObjects(sTDistributionsArray, "c_" + histogramNameStream.str(), outputFile, arguments.outputDirectory + "/" + arguments.outputPrefix + "_" + histogramNameStream.str() + ".pdf", 1024, 768, 0, "", "", false, false, false, 0, 0, 0, 0, 0, 0);
    }
  }
  outputFile->Close();
}

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  std::cout << "Starting code to calculate systematics due to uncertainty on jet energy corrections..." << std::endl;
  std::cout << "Current working directory: " << tmMiscUtils::getCWD() << std::endl;

  tmArgumentParser argumentParser = tmArgumentParser("Calculate systematics due to uncertainty on jet energy corrections.");
  argumentParser.addArgument("inputMCPathMain", "", true, "Path to 2017-optimized MC with no change to jet energies.");
  argumentParser.addArgument("inputPUWeightsPathMain", "", true, "Path to PU reweighting histograms corresponding to inputMCPathMain.");
  argumentParser.addArgument("inputHLTEfficienciesPathMain", "", true, "Path to HLT efficiencies (e.g. higgs->gg) corresponding to inputMCPathMain.") ;
  argumentParser.addArgument("integratedLuminosityMain", "", true, "Integrated luminosity in main MC reference.");
  argumentParser.addArgument("inputMCPathsAux", "", false, "Semicolon-separated list of paths to other MCs.");
  argumentParser.addArgument("inputPUWeightsPathsAux", "", false, "Semicolon-separated list of paths to PU reweighting histograms corresponding to inputMCPathsAux.");
  argumentParser.addArgument("inputHLTEfficienciesPathsAux", "", false, "Semicolon-separated list of paths to HLT efficiencies corresponding to inputMCPathsAux.");
  argumentParser.addArgument("integratedLuminositiesAux", "", false, "Semicolon-separated list of integrated luminosities in auxiliary MC reference.");
  argumentParser.addArgument("crossSectionsFilePath", "", true, "Path to dat file that contains cross-sections as a function of eventProgenitor mass, to use while weighting events.");
  argumentParser.addArgument("eventProgenitor", "", true, "Type of stealth sample. Two possible values: \"squark\" or \"gluino\".");
  argumentParser.addArgument("inputFile_STRegionBoundaries", "STRegionBoundaries.dat", false, "Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity.");
  argumentParser.addArgument("sTMax_toPlot", "4000.0", false, "Max value of sT to plot.");
  argumentParser.addArgument("n_sTBinsToPlot", "28", false, "Number of sT bins to plot."); // default: 28 bins from 1200 to 4000 GeV in steps of 100 GeV
  argumentParser.addArgument("outputDirectory", "analysis/MCEventHistograms/", false, "Prefix to output files.");
  argumentParser.addArgument("outputPrefix", "", true, "Prefix to output files.");
  argumentParser.addArgument("regionsIn_sTHistograms", "1725.0:1775.0:650.0:950.0|1025.0:1075.0:975.0:1075.0", false, "List of the regions in which to fill and save the sT histograms. Each element of the list is in format minEventProgenitorMass:maxEventProgenitorMass:minNeutralinoMass:maxNeutralinoMass, and each element is separated from the next by the character \"|\". See also the default value for this argument.");
  argumentParser.addArgument("MCTemplatePath", "", true, "Path to root file that contains a 2d histogram of the generated masses.");
  argumentParser.setPassedStringValues(argc, argv);
  argumentsStruct arguments = getArgumentsFromParser(argumentParser);
  double integratedLuminosityTotal = arguments.integratedLuminosityMain;
  if (!(arguments.inputMCPathsAux.empty())) {
    for (unsigned int auxIndex = 0; auxIndex < static_cast<unsigned int>((arguments.integratedLuminositiesAux).size()); ++auxIndex) {
      integratedLuminosityTotal += (arguments.integratedLuminositiesAux).at(auxIndex);
    }
  }

  MCTemplateReader templateReader(arguments.MCTemplatePath);
  STRegionsStruct STRegions(arguments.inputFile_STRegionBoundaries, 20000.0);
  outputHistogramsStruct* outputHistograms = initializeOutputHistograms(arguments, templateReader, STRegions);
  std::map<int, double> crossSections = getCrossSectionsFromFile(arguments.crossSectionsFilePath);
  fillOutputHistograms(outputHistograms, arguments, crossSections, templateReader, STRegions, integratedLuminosityTotal);
  saveHistograms(outputHistograms, arguments, STRegions);
  return 0;
}
