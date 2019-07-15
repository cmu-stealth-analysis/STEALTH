#include "../include/getHLTEmulationEfficiency.h"

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  std::cout << "Obtaining HLT emulation efficiency histograms..." << std::endl;
  tmArgumentParser argumentParser = tmArgumentParser("Calculate systematics due to uncertainty on jet energy corrections.");
  argumentParser.addArgument("inputFileName", "", true, "Path to file containing HLT emulation histograms.");
  argumentParser.addArgument("outputFolder", "mergedStatistics", false, "Output folder.");
  argumentParser.addArgument("outputFileName", "HLTEfficiencies.root", false, "Name of output file.");
  argumentParser.setPassedStringValues(argc, argv);
  optionsStruct options = getOptionsFromParser(argumentParser);

  std::map<selectionRegion, TH2F*> outputHistograms_leadingPhoton;
  std::map<selectionRegion, TH2F*> outputHistograms_subLeadingPhoton;

  TFile *inputFile = TFile::Open((options.inputFileName).c_str(), "READ");
  if (!(inputFile->IsOpen()) || inputFile->IsZombie()) {
    std::cout << "ERROR: Unable to open input file to read from. Attempted to read from file with path: " << options.inputFileName << std::endl;
  }
  for (auto&& selectionRegionNamesElement: selectionRegionNames) {
    auto& region = selectionRegionNamesElement.first;
    std::string histogramName;
    histogramName = "hltEfficiency_leadingPhoton_" + selectionRegionNames[region];
    outputHistograms_leadingPhoton[region] = new TH2F(histogramName.c_str(), histogramName.c_str(), HLTEmulation::nEtaBins, HLTEmulation::etaMin, HLTEmulation::etaMax, HLTEmulation::nPTBins, HLTEmulation::PTMin, HLTEmulation::PTMax);
    histogramName = "hltEfficiency_subLeadingPhoton_" + selectionRegionNames[region];
    outputHistograms_subLeadingPhoton[region] = new TH2F(histogramName.c_str(), histogramName.c_str(), HLTEmulation::nEtaBins, HLTEmulation::etaMin, HLTEmulation::etaMax, HLTEmulation::nPTBins, HLTEmulation::PTMin, HLTEmulation::PTMax);
    TH2I *inputHistogram_denominator_leadingPhoton;
    TH2I *inputHistogram_denominator_subLeadingPhoton;
    TH2I *inputHistogram_numerator_leadingPhoton;
    TH2I *inputHistogram_numerator_subLeadingPhoton;
    histogramName = std::string("hltEfficiency_leadingPhoton_totalEvents_" + selectionRegionNames[region]);
    inputFile->GetObject(histogramName.c_str(), inputHistogram_denominator_leadingPhoton);
    histogramName = std::string("hltEfficiency_subLeadingPhoton_totalEvents_" + selectionRegionNames[region]);
    inputFile->GetObject(histogramName.c_str(), inputHistogram_denominator_subLeadingPhoton);
    histogramName = std::string("hltEfficiency_leadingPhoton_passingEmulation_" + selectionRegionNames[region]);
    inputFile->GetObject(histogramName.c_str(), inputHistogram_numerator_leadingPhoton);
    histogramName = std::string("hltEfficiency_subLeadingPhoton_passingEmulation_" + selectionRegionNames[region]);
    inputFile->GetObject(histogramName.c_str(), inputHistogram_numerator_subLeadingPhoton);

    if (inputHistogram_denominator_leadingPhoton && inputHistogram_denominator_subLeadingPhoton && inputHistogram_numerator_leadingPhoton && inputHistogram_numerator_subLeadingPhoton) {
      std::cout << "Opened input histograms." << std::endl;
    }
    else {
      std::cout << "ERROR: unable to open one or more input histogram." << std::endl;
      std::exit(EXIT_FAILURE);
    }

    for (int xbincounter = 0; xbincounter <= 1+((outputHistograms_leadingPhoton[region])->GetXaxis()->GetNbins()); ++xbincounter) {
      for (int ybincounter = 0; ybincounter <= 1+((outputHistograms_leadingPhoton[region])->GetYaxis()->GetNbins()); ++ybincounter) {
        int leading_numerator = inputHistogram_numerator_leadingPhoton->GetBinContent(xbincounter, ybincounter);
        int leading_denominator = inputHistogram_denominator_leadingPhoton->GetBinContent(xbincounter, ybincounter);
        float leading_ratio = 0.;
        if (leading_denominator > 0.) leading_ratio = (1.0*leading_numerator)/(1.0*leading_denominator);
        outputHistograms_leadingPhoton[region]->SetBinContent(xbincounter, ybincounter, leading_ratio);

        int subLeading_numerator = inputHistogram_numerator_subLeadingPhoton->GetBinContent(xbincounter, ybincounter);
        int subLeading_denominator = inputHistogram_denominator_subLeadingPhoton->GetBinContent(xbincounter, ybincounter);
        float subLeading_ratio = 0.;
        if (subLeading_denominator > 0.) subLeading_ratio = (1.0*subLeading_numerator)/(1.0*subLeading_denominator);
        outputHistograms_subLeadingPhoton[region]->SetBinContent(xbincounter, ybincounter, subLeading_ratio);
      }
    }
  }
  // inputFile->Close();

  TFile *outputFile = TFile::Open((options.outputFolder + "/" + options.outputFileName).c_str(), "RECREATE");
  if (!(outputFile->IsOpen()) || outputFile->IsZombie()) {
    std::cout << "ERROR: Unable to open output file to write to. Attempted to write to file with path: " << (options.outputFolder + "/" + options.outputFileName) << std::endl;
  }
  for (auto&& selectionRegionNamesElement: selectionRegionNames) {
    auto& region = selectionRegionNamesElement.first;
    outputFile->WriteTObject(outputHistograms_leadingPhoton[region]);
    outputFile->WriteTObject(outputHistograms_subLeadingPhoton[region]);
  }
  outputFile->Close();
  std::cout << "Finished writing HLT emulation efficiency histograms to output." << std::endl;
  inputFile->Close();
  return 0;
}
