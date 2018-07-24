#include "../include/getJECUncertaintySystematics.h"

void setGlobalVariables(tmArgumentParser& argumentParser) {
  inputPath_ = argumentParser.getArgumentString("inputPath");
  MCTemplate_ = argumentParser.getArgumentString("MCTemplate");
  outputDirectory_ = argumentParser.getArgumentString("outputDirectory");
  outputPrefix_ = argumentParser.getArgumentString("outputPrefix");
  nGluinoMassBins_ = std::stoi(argumentParser.getArgumentString("nGluinoMassBins"));
  minGluinoMass_ = std::stod(argumentParser.getArgumentString("minGluinoMass"));
  maxGluinoMass_ = std::stod(argumentParser.getArgumentString("maxGluinoMass"));
  nNeutralinoMassBins_ = std::stoi(argumentParser.getArgumentString("nNeutralinoMassBins"));
  minNeutralinoMass_ = std::stod(argumentParser.getArgumentString("minNeutralinoMass"));
  maxNeutralinoMass_ = std::stod(argumentParser.getArgumentString("maxNeutralinoMass"));
  sTMin_normWindow_ = std::stod(argumentParser.getArgumentString("sTMin_normWindow"));
  sTMax_normWindow_ = std::stod(argumentParser.getArgumentString("sTMax_normWindow"));
  sTStartMainRegion_ = std::stod(argumentParser.getArgumentString("sTStartMainRegion"));
  minGluinoMassToPlot_ = std::stod(argumentParser.getArgumentString("minGluinoMassToPlot"));
}

std::string getHistogramName(std::string histogramType, std::string zone, int nJetsBin) {
  std::stringstream nameStream;
  nameStream << histogramType << "_" << nJetsBin << "Jets_" << zone;
  return nameStream.str();
}

std::string getHistogramTitle(std::string histogramType, std::string zone, int nJetsBin) {
  std::string histogramTypeString;
  if (histogramType == "JECUncertainty") histogramTypeString = "JEC Uncertainty";
  else if (histogramType == "JECUncertainty_fractionalError") histogramTypeString = "Estimated fractional error on JEC Uncertainty";
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
  if (zone == "sub") sTRangeStringStream << std::fixed << std::setprecision(0) << sTMax_normWindow_ << " < #it{S}_{T} < " << sTStartMainRegion_;
  else if (zone == "main") sTRangeStringStream << "#it{S}_{T} > " << std::fixed << std::setprecision(0) << sTStartMainRegion_;
  else {
    std::cout << "ERROR: Unknown zone: " << zone << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string sTRangeString = sTRangeStringStream.str();

  std::string axesLabelsString = ";m_{#tilde{#it{g}}};m_{#tilde{#it{#chi_{1}^{0}}}}";

  std::stringstream titleStream;
  titleStream << histogramTypeString << ", " << nJetsString << ", " << sTRangeString << axesLabelsString;
  return titleStream.str();
}

void initializePlots() {
  for (auto zone: allowedZones) {
    for (int nJetsBin = 4; nJetsBin <= 6; ++nJetsBin) {
      g_JECUncertainty[zone][nJetsBin] = new TGraph2D();
      g_JECUncertainty[zone][nJetsBin]->SetName(("g_" + getHistogramName("JECUncertainty", zone, nJetsBin)).c_str());
      g_JECUncertainty[zone][nJetsBin]->SetTitle(getHistogramTitle("JECUncertainty", zone, nJetsBin).c_str());
      g_JECUncertainty[zone][nJetsBin]->SetNpx(16);
      g_JECUncertainty[zone][nJetsBin]->SetNpy(133);
      h_JECUncertainty[zone][nJetsBin] = new TH2F(("h_" + getHistogramName("JECUncertainty", zone, nJetsBin)).c_str(), getHistogramTitle("JECUncertainty", zone, nJetsBin).c_str(), nGluinoMassBins_, minGluinoMass_, maxGluinoMass_, nNeutralinoMassBins_, minNeutralinoMass_, maxNeutralinoMass_);
      
      g_JECUncertainty_fractionalError[zone][nJetsBin] = new TGraph2D();
      g_JECUncertainty_fractionalError[zone][nJetsBin]->SetName(("g_" + getHistogramName("JECUncertainty_fractionalError", zone, nJetsBin)).c_str());
      g_JECUncertainty_fractionalError[zone][nJetsBin]->SetTitle(getHistogramTitle("JECUncertainty_fractionalError", zone, nJetsBin).c_str());
      g_JECUncertainty_fractionalError[zone][nJetsBin]->SetNpx(16);
      g_JECUncertainty_fractionalError[zone][nJetsBin]->SetNpy(133);
      h_JECUncertainty_fractionalError[zone][nJetsBin] = new TH2F(("h_" + getHistogramName("JECUncertainty_fractionalError", zone, nJetsBin)).c_str(), getHistogramTitle("JECUncertainty_fractionalError", zone, nJetsBin).c_str(), nGluinoMassBins_, minGluinoMass_, maxGluinoMass_, nNeutralinoMassBins_, minNeutralinoMass_, maxNeutralinoMass_);
    } // end loop over nJetsBin
  } // end loop over zone
}

void getSystematics() {
  std::cout << "Getting systematics..." << std::endl;

  // first fetch the histograms
  TFile *inputFile = TFile::Open(inputPath_.c_str(), "READ");
  for (auto jec: allowedJECs) {
    for (auto zone: allowedZones) {
      for (int nJetsBin = 4; nJetsBin <= 6; ++nJetsBin) {
        std::stringstream nameStreamTotalNEvents;
        nameStreamTotalNEvents << "h_total_nMCEvents_" << jec << "_" << nJetsBin << "Jets_" << zone;
        h_totalNEvents[jec][zone][nJetsBin] = (TH2F*)(inputFile->Get(nameStreamTotalNEvents.str().c_str()));
        if (h_totalNEvents[jec][zone][nJetsBin] == nullptr) {
          std::cout << "Unable to find histogram with name " << nameStreamTotalNEvents.str() << std::endl;
          std::exit(EXIT_FAILURE);
        }
        std::stringstream nameStreamWeightedNEvents;
        nameStreamWeightedNEvents << "h_weighted_nMCEvents_" << jec << "_" << nJetsBin << "Jets_" << zone;
        h_weightedNEvents[jec][zone][nJetsBin] = (TH2F*)(inputFile->Get(nameStreamWeightedNEvents.str().c_str()));
        if (h_weightedNEvents[jec][zone][nJetsBin] == nullptr) {
          std::cout << "Unable to find histogram with name " << nameStreamWeightedNEvents.str() << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
    }
  }

  // get MC template
  TFile *MCTemplateFile = TFile::Open(MCTemplate_.c_str(), "READ");
  TH2F *MCTemplateTH2 = (TH2F*)(MCTemplateFile->Get("h_susyMasses_template"));
  if (MCTemplateTH2 == nullptr) {
    std::cout << "ERROR: MC template not found" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  
  // Now fill TGraphs and TH2s with the JEC fractional uncertainty and estimated error
  for (auto zone: allowedZones) {
    for (int nJetsBin = 4; nJetsBin <= 6; ++nJetsBin) {
      // Loop over only those bins that show a "1" in the MC template file
      for (int templateGluinoMassIndex = 1; templateGluinoMassIndex <= MCTemplateTH2->GetXaxis()->GetNbins(); ++templateGluinoMassIndex) {
        double gluinoMass = MCTemplateTH2->GetXaxis()->GetBinCenter(templateGluinoMassIndex);
        // if (gluinoMass < minGluinoMassToPlot_) continue;
        for (int templateNeutralinoMassIndex = 1; templateNeutralinoMassIndex <= MCTemplateTH2->GetYaxis()->GetNbins(); ++templateNeutralinoMassIndex) {
          if (!(1 == static_cast<int>(0.5 + MCTemplateTH2->GetBinContent(templateGluinoMassIndex, templateNeutralinoMassIndex)))) continue;
          double neutralinoMass = MCTemplateTH2->GetYaxis()->GetBinCenter(templateNeutralinoMassIndex);
          double weightedNEvents_jecUp = h_weightedNEvents["JECUp"][zone][nJetsBin]->GetBinContent(h_weightedNEvents["JECUp"][zone][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass));
          double weightedNEvents_nominal = h_weightedNEvents["JECNominal"][zone][nJetsBin]->GetBinContent(h_weightedNEvents["JECNominal"][zone][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass));
          double weightedNEvents_jecDown = h_weightedNEvents["JECDown"][zone][nJetsBin]->GetBinContent(h_weightedNEvents["JECDown"][zone][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass));
          double totalNEvents_nominal = h_totalNEvents["JECNominal"][zone][nJetsBin]->GetBinContent(h_totalNEvents["JECNominal"][zone][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass));
          double totalNEventsError_nominal = h_totalNEvents["JECNominal"][zone][nJetsBin]->GetBinError(h_totalNEvents["JECNominal"][zone][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass));
          if (totalNEvents_nominal == 0) {
            std::cout << "WARNING: zero events recorded at gluino mass: " << gluinoMass << ", neutralino mass: " << neutralinoMass << ", for zone: " << zone << ", nJets: " << nJetsBin << std::endl;
            continue;
          }
          if (weightedNEvents_nominal == 0) {
            std::cout << "ERROR: total number of recorded events is nonzero but weighted number of recorded events is 0 at gluino mass: " << gluinoMass << ", neutralino mass: " << neutralinoMass << ", for zone: " << zone << ", nJets: " << nJetsBin << std::endl;
            std::exit(EXIT_FAILURE);
          }
          double ratioUpToNominal = weightedNEvents_jecUp / weightedNEvents_nominal;
          double deviationUp = std::fabs(ratioUpToNominal - 1.0);
          double ratioDownToNominal = weightedNEvents_jecDown / weightedNEvents_nominal;
          double deviationDown = std::fabs(ratioDownToNominal - 1.0);
          double averageDeviation = 0.5*(deviationUp + deviationDown);
          if (g_JECUncertainty[zone][nJetsBin] == nullptr) {
            g_JECUncertainty[zone][nJetsBin]->SetPoint(1, gluinoMass, neutralinoMass, averageDeviation);
          }
          else {
            g_JECUncertainty[zone][nJetsBin]->SetPoint(1 + g_JECUncertainty[zone][nJetsBin]->GetN(), gluinoMass, neutralinoMass, averageDeviation);
          }
          h_JECUncertainty[zone][nJetsBin]->SetBinContent(h_JECUncertainty[zone][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass), averageDeviation);
          double fractionalUncertaintyOnDeviation = std::sqrt(2)*totalNEventsError_nominal/totalNEvents_nominal;
          g_JECUncertainty_fractionalError[zone][nJetsBin]->SetPoint(g_JECUncertainty_fractionalError[zone][nJetsBin]->GetN(), gluinoMass, neutralinoMass, fractionalUncertaintyOnDeviation);
          h_JECUncertainty_fractionalError[zone][nJetsBin]->SetBinContent(h_JECUncertainty_fractionalError[zone][nJetsBin]->FindFixBin(gluinoMass, neutralinoMass), fractionalUncertaintyOnDeviation);
        } // end loop over neutralino mass
      } // end loop over gluino mass
    } // end loop over nJetsBin
  } // end loop over zone
  MCTemplateFile->Close();
  inputFile->Close();
}

void savePlots() {
  std::cout << "Saving output plots..." << std::endl;
  TFile *outputFile = TFile::Open((outputDirectory_ + "/" + outputPrefix_ + "_savedObjects.root").c_str(), "RECREATE");
  for (auto zone: allowedZones) {
    for (int nJetsBin = 4; nJetsBin <= 6; ++nJetsBin) {
      std::string histogramName_JECUncertainty = getHistogramName("JECUncertainty", zone, nJetsBin);
      tmROOTSaverUtils::saveSingleObject(g_JECUncertainty[zone][nJetsBin], "c_g_" + histogramName_JECUncertainty, outputFile, outputDirectory_ + "/" + outputPrefix_ + "_" + histogramName_JECUncertainty + "_graph.png", 1024, 768, 0, "", "COLZ", false, false, true, minGluinoMassToPlot_, maxGluinoMass_, 0, 0, 0, 0);
      tmROOTSaverUtils::saveSingleObject(h_JECUncertainty[zone][nJetsBin], "c_h_" + histogramName_JECUncertainty, outputFile, outputDirectory_ + "/" + outputPrefix_ + "_" + histogramName_JECUncertainty + "_hist.png", 1024, 768, 0, ".0e", "TEXTCOLZ", false, false, true, minGluinoMassToPlot_, maxGluinoMass_, 0, 0, 0, 0);
      std::string histogramName_JECUncertainty_fractionalError = getHistogramName("JECUncertainty_fractionalError", zone, nJetsBin);
      tmROOTSaverUtils::saveSingleObject(g_JECUncertainty_fractionalError[zone][nJetsBin], "c_g_" + histogramName_JECUncertainty_fractionalError, outputFile, outputDirectory_ + "/" + outputPrefix_ + "_" + histogramName_JECUncertainty_fractionalError + "_graph.png", 1024, 768, 0, "", "COLZ", false, false, true, minGluinoMassToPlot_, maxGluinoMass_, 0, 0, 0, 0);
      tmROOTSaverUtils::saveSingleObject(h_JECUncertainty_fractionalError[zone][nJetsBin], "c_h_" + histogramName_JECUncertainty_fractionalError, outputFile, outputDirectory_ + "/" + outputPrefix_ + "_" + histogramName_JECUncertainty_fractionalError + "_hist.png", 1024, 768, 0, ".0e", "TEXTCOLZ", false, false, true, minGluinoMassToPlot_, maxGluinoMass_, 0, 0, 0, 0);
    }
  }
  outputFile->Close();
}

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  std::cout << "Starting program to get systematics due to uncertainty on JECs..." << std::endl;
  tmArgumentParser argumentParser = tmArgumentParser("Read in event histograms and calculate systematics due to uncertainty on jet energy corrections.");
  argumentParser.addArgument("inputPath", "analysis/JECUncertainties/MC2018_withJECUncertainties_savedObjects.root", true, "Path to ROOT file containing event histograms.");
  argumentParser.addArgument("MCTemplate", "plot_susyMasses_template.root", false, "Path to root file that contains a TH2F with bins containing points with generated masses set to 1 and all other bins set to 0.");
  argumentParser.addArgument("outputDirectory", "analysis/MCSystematics/", false, "Prefix to output files.");
  argumentParser.addArgument("outputPrefix", "", true, "Prefix to output files.");
  argumentParser.addArgument("nGluinoMassBins", "20", false, "nBins on the gluino mass axis."); // (800 - 25) GeV --> (1750 + 25) GeV in steps of 50 GeV
  argumentParser.addArgument("minGluinoMass", "775.0", false, "Min gluino mass for the 2D plots.");
  argumentParser.addArgument("maxGluinoMass", "1775.0", false, "Max gluino mass for the 2D plots.");
  argumentParser.addArgument("nNeutralinoMassBins", "133", false, "nBins on the neutralino mass axis.");
  argumentParser.addArgument("minNeutralinoMass", "93.75", false, "Min neutralino mass for the 2D plots.");
  argumentParser.addArgument("maxNeutralinoMass", "1756.25", false, "Max neutralino mass for the 2D plots."); // (100 - 6.25) GeV --> (1750 + 6.25) GeV in steps of 12.5 GeV
  argumentParser.addArgument("sTMin_normWindow", "1200.0", false, "Lower sT boundary of normalization window.");
  argumentParser.addArgument("sTMax_normWindow", "1300.0", false, "Upper sT boundary of normalization window.");
  argumentParser.addArgument("sTStartMainRegion", "2500.0", false, "Lowest value of sT in main observation bin.");
  argumentParser.addArgument("minGluinoMassToPlot", "975.0", false, "Min gluino mass to use in the 2D plots.");
  
  argumentParser.setPassedStringValues(argc, argv);
  setGlobalVariables(argumentParser);
  initializePlots();
  getSystematics();
  savePlots();
  return 0;
}
