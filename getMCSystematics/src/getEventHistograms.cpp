#include "../include/getEventHistograms.h"

std::string getHistogramName(std::string histogramType, std::string jec, std::string zone, int nJetsBin) {
  std::stringstream nameStream;
  if (histogramType == "sTDistribution") {
    nameStream << histogramType << "_" << jec << "_" << nJetsBin << "Jets";
    return nameStream.str();
  }
  nameStream << histogramType << "_nMCEvents_" << jec << "_" << nJetsBin << "Jets_" << zone;
  return nameStream.str();
}

std::string getHistogramTitle(std::string histogramType, std::string jec, std::string zone, int nJetsBin) {
  if (histogramType == "sTDistribution") {
    std::stringstream histogramTitle;
    histogramTitle << "sT Distributions, " << nJetsBin << "Jets;#it{S}_{T}(GeV);nEvents/(" << static_cast<int>(0.5 + ((sTMax_toPlot - sTMin_normWindow)/n_sTBinsToPlot)) << " GeV)";
    return histogramTitle.str();
  }
  std::string histogramTypeString;
  if (histogramType == "total") histogramTypeString = "Total MC Events";
  else if (histogramType == "weighted") histogramTypeString = "Weighted MC Events";
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
  if (zone == "norm") sTRangeStringStream << std::fixed << std::setprecision(0) << sTMin_normWindow << " < #it{S}_{T} < " << sTMax_normWindow;
  else if (zone == "sub") sTRangeStringStream << std::fixed << std::setprecision(0) << sTMax_normWindow << " < #it{S}_{T} < " << sTStartMainRegion;
  else if (zone == "main") sTRangeStringStream << "#it{S}_{T} > " << std::fixed << std::setprecision(0) << sTStartMainRegion;
  else {
    std::cout << "ERROR: Unknown zone: " << zone << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string sTRangeString = sTRangeStringStream.str();

  std::string axesLabelsString = ";m_{#tilde{#it{g}}};m_{#tilde{#it{#chi_{1}^{0}}}}";

  std::stringstream titleStream;
  titleStream << histogramTypeString << ", " << jecString << ", " << nJetsString << ", " << sTRangeString << axesLabelsString;
  return titleStream.str();
}

void setGlobalVariables(tmArgumentParser& argumentParser) {
  inputMCPath = argumentParser.getArgumentString("inputMCPath");
  inputMCPath_JECUp = argumentParser.getArgumentString("inputMCPath_JECUp");
  inputMCPath_JECDown = argumentParser.getArgumentString("inputMCPath_JECDown");
  maxMCEvents = std::stol(argumentParser.getArgumentString("maxMCEvents"));
  crossSectionsFilePath = argumentParser.getArgumentString("crossSectionsFilePath");
  MCTemplate = argumentParser.getArgumentString("MCTemplate");
  sTMin_normWindow = std::stod(argumentParser.getArgumentString("sTMin_normWindow"));
  sTMax_normWindow = std::stod(argumentParser.getArgumentString("sTMax_normWindow"));
  sTStartMainRegion = std::stod(argumentParser.getArgumentString("sTStartMainRegion"));
  sTMax_toPlot = std::stod(argumentParser.getArgumentString("sTMax_toPlot"));
  n_sTBinsToPlot = std::stoi(argumentParser.getArgumentString("n_sTBinsToPlot"));
  outputDirectory = argumentParser.getArgumentString("outputDirectory");
  outputPrefix = argumentParser.getArgumentString("outputPrefix");
  integratedLuminosity = std::stod(argumentParser.getArgumentString("integratedLuminosity"));
  nGeneratedEventsPerBin = std::stoi(argumentParser.getArgumentString("nGeneratedEventsPerBin"));
  nGluinoMassBins = std::stoi(argumentParser.getArgumentString("nGluinoMassBins"));
  minGluinoMass = std::stod(argumentParser.getArgumentString("minGluinoMass"));
  maxGluinoMass = std::stod(argumentParser.getArgumentString("maxGluinoMass"));
  nNeutralinoMassBins = std::stoi(argumentParser.getArgumentString("nNeutralinoMassBins"));
  minNeutralinoMass = std::stod(argumentParser.getArgumentString("minNeutralinoMass"));
  maxNeutralinoMass = std::stod(argumentParser.getArgumentString("maxNeutralinoMass"));
}

void initializeHistograms() {
  // 2D histograms for nEvents
  for (auto jec: allowedJECs) {
    std::map< std::string, std::map< int, TH2F*> > mapForJEC_total;
    std::map< std::string, std::map< int, TH2F*> > mapForJEC_weighted;
    for (auto zone: allowedZones) {
      std::map<int, TH2F*> mapForZone_total;
      std::map<int, TH2F*> mapForZone_weighted;
      for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
        mapForZone_total[nJetsBin] = new TH2F(("h_" + getHistogramName("total", jec, zone, nJetsBin)).c_str(), getHistogramTitle("total", jec, zone, nJetsBin).c_str(), nGluinoMassBins, minGluinoMass, maxGluinoMass, nNeutralinoMassBins, minNeutralinoMass, maxNeutralinoMass);
        mapForZone_weighted[nJetsBin] = new TH2F(("h_" + getHistogramName("weighted", jec, zone, nJetsBin)).c_str(), getHistogramTitle("weighted", jec, zone, nJetsBin).c_str(), nGluinoMassBins, minGluinoMass, maxGluinoMass, nNeutralinoMassBins, minNeutralinoMass, maxNeutralinoMass);
      }
      mapForJEC_total[zone] = mapForZone_total;
      mapForJEC_weighted[zone] = mapForZone_weighted;
    }
    h_totalNEvents[jec] = mapForJEC_total;
    h_weightedNEvents[jec] = mapForJEC_weighted;
  }

  // 1D sT distributions
  for (auto jec: allowedJECs) {
    for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
      h_sTDistributions[jec][nJetsBin] = new TH1F(("h_" + getHistogramName("sTDistribution", jec, "", nJetsBin)).c_str(), getHistogramTitle("sTDistribution", jec, "", nJetsBin).c_str(), n_sTBinsToPlot, sTMin_normWindow, sTMax_toPlot);
    }
  }
}

void setCrossSectionsFromFile(bool printCrossSections = false) {
  std::ifstream crossSectionsFileStream(crossSectionsFilePath.c_str());
  if (!(crossSectionsFileStream.is_open())) {
    std::cout << "ERROR: Unable to open file: " << crossSectionsFilePath << std::endl;
    std::exit(EXIT_FAILURE);
  }
  int gluinoMass;
  double crossSection, crossSectionPercentUncertainty;
  while (crossSectionsFileStream >> gluinoMass >> crossSection >> crossSectionPercentUncertainty) {
    crossSections[gluinoMass] = crossSection;
    crossSectionsFractionalUncertainty[gluinoMass] = 0.01*crossSectionPercentUncertainty;
    if (printCrossSections) {
      std::cout << "crossSections[" << gluinoMass << "] = " << crossSection << ", "
                << "crossSectionsFractionalUncertainty[" << gluinoMass << "] = " << 0.01*crossSectionPercentUncertainty << std::endl;
    }
  }
}

void fillHistogramsForJEC(std::string jec) {
  std::cout << "Filling events map for JEC type: " << jec << std::endl;
  std::string inputPath;
  if (jec == "JECDown") inputPath = inputMCPath_JECDown;
  else if (jec == "JECNominal") inputPath = inputMCPath;
  else if (jec == "JECUp") inputPath = inputMCPath_JECUp;
  else {
    std::cout << "ERROR: unknown JEC type: \"" << jec << "\"" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  TFile *inputFile = TFile::Open(inputPath.c_str(), "READ");
  if (!inputFile || inputFile->IsZombie()) {
    std::cout << "Error in opening file with path: " << inputPath << std::endl;
  }
  TTree *inputTree = (TTree*) inputFile->Get("ggNtuplizer/EventTree");
  TTreeReader inputTreeReader(inputTree);
  TTreeReaderValue<int> evt_nJets(inputTreeReader, "b_nJets");
  TTreeReaderValue<double> evt_ST(inputTreeReader, "b_evtST");
  TTreeReaderValue<int> nMC(inputTreeReader, "nMC");
  TTreeReaderArray<int> mcPIDs(inputTreeReader, "mcPID");
  TTreeReaderArray<float> mcMasses(inputTreeReader, "mcMass");
  TTreeReaderArray<int> mcMomPIDs(inputTreeReader, "mcMomPID");
  TTreeReaderArray<float> mcMomMasses(inputTreeReader, "mcMomMass");
  long nEntriesSource = inputTree->GetEntries();
  std::cout << "Number of entries found: " << nEntriesSource << std::endl;
  long nEntriesToRead = (maxMCEvents > 0) ? (maxMCEvents < nEntriesSource ? maxMCEvents : nEntriesSource) : nEntriesSource;
  std::cout << "Looping over " << nEntriesToRead << " entries." << std::endl;
  tmProgressBar *progressBar = new tmProgressBar(nEntriesToRead);
  int tmp = static_cast<int>(0.5 + 1.0*nEntriesToRead/1000);
  int progressBarUpdatePeriod = tmp > 0? tmp: 1;
  int entryIndex = 0;
  progressBar->initialize();
  while (inputTreeReader.Next()) {
    if (entryIndex >= nEntriesToRead) {
      if (maxMCEvents == 0) {
        std::cout << "ERROR: entry index seems to be greater than available number of events" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      break;
    }
    if (entryIndex % progressBarUpdatePeriod == 0) progressBar->updateBar(1.0*entryIndex/nEntriesToRead, entryIndex);
    // std::cout << "At entry index = " << entryIndex << ", nJets = " << *evt_nJets << ", ST = " << *evt_ST << ", nMC = " << *nMC << ", mcPID size: " << mcPIDs.GetSize() << ", mcMomPIDs size: " << mcMomPIDs.GetSize() << std::endl;

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
    
    if (*evt_ST < sTMin_normWindow) continue;

    std::string zone;
    if (*evt_ST >= sTMin_normWindow && *evt_ST < sTMax_normWindow) zone = "norm";
    else if (*evt_ST >= sTMax_normWindow && *evt_ST < sTStartMainRegion) zone = "sub";
    else if (*evt_ST >= sTStartMainRegion) zone = "main";
    else {
      std::cout << "ERROR: logic error" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    int nJetsBin = *evt_nJets;
    if (*evt_nJets > 6) nJetsBin = 6;
    h_sTDistributions[jec][nJetsBin]->Fill(*evt_ST);

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
    double eventWeight = crossSections[gMassInt]*integratedLuminosity/nGeneratedEventsPerBin;

    h_totalNEvents[jec][zone][nJetsBin]->Fill(static_cast<double>(generated_gluinoMass), static_cast<double>(generated_neutralinoMass), 1.0);
    h_weightedNEvents[jec][zone][nJetsBin]->Fill(static_cast<double>(generated_gluinoMass), static_cast<double>(generated_neutralinoMass), eventWeight);
    
    ++entryIndex;
  }
  delete progressBar;
  inputFile->Close();
}

void saveHistograms() {
  std::cout << "Saving histograms..." << std::endl;

  // First the 2D event histograms
  TFile *outputFile = TFile::Open((outputDirectory + "/" + outputPrefix + "_savedObjects.root").c_str(), "RECREATE");
  for (auto jec: allowedJECs) {
    for (auto zone: allowedZones) {
      for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
        std::string histogramName_total = getHistogramName("total", jec, zone, nJetsBin);
        tmROOTSaverUtils::saveSingleObject(h_totalNEvents[jec][zone][nJetsBin], "c_" + histogramName_total, outputFile, outputDirectory + "/" + outputPrefix + "_" + histogramName_total + ".png", 1024, 768, 0, ".0f", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);
        std::string histogramName_weighted = getHistogramName("weighted", jec, zone, nJetsBin);
        tmROOTSaverUtils::saveSingleObject(h_weightedNEvents[jec][zone][nJetsBin], "c_" + histogramName_weighted, outputFile, outputDirectory + "/" + outputPrefix + "_" + histogramName_weighted + ".png", 1024, 768, 0, ".0f", "TEXTCOLZ", false, false, true, 0, 0, 0, 0, 0, 0);
      }
    }
  }

  // next the 1D sT distributions
  for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
    std::stringstream histogramNameStream;
    histogramNameStream << "sTDistributions_" << nJetsBin << "Jets";
    TObjArray *sTDistributionsArray = new TObjArray(4);
    TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    std::stringstream legendTitleStream;
    legendTitleStream << "JEC uncertainty shifts, nJets = " << nJetsBin;
    legend->SetHeader(legendTitleStream.str().c_str());
    for (auto jec: allowedJECs) {
      int indexToAddAt = (jec == "JECNominal"? 0 : (jec == "JECUp" ? 1 : (jec == "JECDown" ? 2 : -1)));
      if (indexToAddAt < 0) {
        std::cout << "ERROR: bad indexToAddAt = " << indexToAddAt << std::endl;
        std::exit(EXIT_FAILURE);
      }
      sTDistributionsArray->AddAt(h_sTDistributions[jec][nJetsBin], indexToAddAt);
      h_sTDistributions[jec][nJetsBin]->SetLineColor(jec == "JECNominal"? kBlack : (jec == "JECUp" ? kRed : (jec == "JECDown" ? kBlue : kGreen)));
      TLegendEntry *legendEntry = legend->AddEntry(h_sTDistributions[jec][nJetsBin], jec == "JECNominal" ? "No JEC uncertainty corrections" : (jec == "JECUp" ? "Shifted up by JEC uncertainty" : (jec == "JECDown" ? "Shifted down by JEC uncertainty" : "Unknown JEC")));
      legendEntry->SetTextColor(jec == "JECNominal"? kBlack : (jec == "JECUp" ? kRed : (jec == "JECDown" ? kBlue : kGreen)));
    }
    sTDistributionsArray->AddAt(legend, 3);
    tmROOTSaverUtils::saveObjects(sTDistributionsArray, "c_" + histogramNameStream.str(), outputFile, outputDirectory + "/" + outputPrefix + "_" + histogramNameStream.str() + ".png", 1024, 768, 0, "", "", false, false, false, 0, 0, 0, 0, 0, 0);
  }
  outputFile->Close();
}

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  std::cout << "Starting code to calculate systematics due to uncertainty on jet energy corrections..." << std::endl;
  std::cout << "Current working directory: " << tmMiscUtils::getCWD() << std::endl;
  tmArgumentParser argumentParser = tmArgumentParser("Calculate systematics due to uncertainty on jet energy corrections.");
  argumentParser.addArgument("inputMCPath", "root://cmseos.fnal.gov//store/user/lpcsusystealth/selections/combined/MC_2018Production_DoubleMedium.root", true, "Path to MC with no change to jet energies.");
  argumentParser.addArgument("inputMCPath_JECUp", "root://cmseos.fnal.gov//store/user/lpcsusystealth/selections/combined/MC_JECUp_2018Production_DoubleMedium.root", true, "Path to MC with all energies shifted up by the JEC uncertainty.");
  argumentParser.addArgument("inputMCPath_JECDown", "root://cmseos.fnal.gov//store/user/lpcsusystealth/selections/combined/MC_JECDown_2018Production_DoubleMedium.root", true, "Path to MC with all energies shifted up by the JEC uncertainty.");
  argumentParser.addArgument("maxMCEvents", "0", false, "Set a custom maximum number of MC events.");
  argumentParser.addArgument("crossSectionsFilePath", "SusyCrossSections13TevGluGlu.txt", false, "Path to dat file that contains cross-sections as a function of gluino mass, to use while weighting events.");
  argumentParser.addArgument("MCTemplate", "plot_susyMasses_template.root", false, "Path to root file that contains a TH2F with bins containing points with generated masses set to 1 and all other bins set to 0.");
  argumentParser.addArgument("sTMin_normWindow", "1200.0", false, "Lower sT boundary of normalization window.");
  argumentParser.addArgument("sTMax_normWindow", "1300.0", false, "Upper sT boundary of normalization window.");
  argumentParser.addArgument("sTStartMainRegion", "2500.0", false, "Lowest value of sT in main observation bin.");
  argumentParser.addArgument("sTMax_toPlot", "3500.0", false, "Max value of sT to plot.");
  argumentParser.addArgument("n_sTBinsToPlot", "23", false, "Number of sT bins to plot."); // default: 23 bins from 1200 to 3500 GeV in steps of 100 GeV
  argumentParser.addArgument("outputDirectory", "analysis/JECUncertainties/", false, "Prefix to output files.");
  argumentParser.addArgument("outputPrefix", "", true, "Prefix to output files.");
  argumentParser.addArgument("integratedLuminosity", "83780.0", false, "Lowest value of sT in main observation bin.");
  argumentParser.addArgument("nGeneratedEventsPerBin", "150000", false, "Number of generated events per bin, to use while calculating event weights.");
  argumentParser.addArgument("nGluinoMassBins", "20", false, "nBins on the gluino mass axis."); // (800 - 25) GeV --> (1750 + 25) GeV in steps of 50 GeV
  argumentParser.addArgument("minGluinoMass", "775.0", false, "Min gluino mass for the 2D plots.");
  argumentParser.addArgument("maxGluinoMass", "1775.0", false, "Max gluino mass for the 2D plots.");
  argumentParser.addArgument("nNeutralinoMassBins", "133", false, "nBins on the neutralino mass axis.");
  argumentParser.addArgument("minNeutralinoMass", "93.75", false, "Min neutralino mass for the 2D plots.");
  argumentParser.addArgument("maxNeutralinoMass", "1756.25", false, "Max neutralino mass for the 2D plots."); // (100 - 6.25) GeV --> (1750 + 6.25) GeV in steps of 12.5 GeV
  argumentParser.setPassedStringValues(argc, argv);
  setGlobalVariables(argumentParser);
  initializeHistograms();

  setCrossSectionsFromFile();
  for (auto jec: allowedJECs) {
    fillHistogramsForJEC(jec);
  }
  saveHistograms();
  return 0;
}
