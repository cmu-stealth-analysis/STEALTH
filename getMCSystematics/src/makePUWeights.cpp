#include "../include/makePUWeights.h"

void fillMCHistogramBinsThatAreNonzeroInDataFromFile(TH1D* inputHistogram_MC, TH1D* inputHistogram_data, const std::string& inputMCPath) {
  TChain *inputChain = new TChain("ggNtuplizer/EventTree");
  inputChain->Add(inputMCPath.c_str());
  long nEntries = inputChain->GetEntries();
  assert((nEntries > 0));
  std::cout << "Number of entries in " << inputMCPath << ": " << nEntries << std::endl;
  TTreeReader inputTreeReader(inputChain);
  TTreeReaderArray<int> evt_BX_for_PU(inputTreeReader, "puBX");
  TTreeReaderArray<float> evt_PU(inputTreeReader, "puTrue");
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
    if (entryIndex % progressBarUpdatePeriod == 0) progressBar->updateBar(1.0*entryIndex/nEntries, entryIndex);

    float eventPU = -1.;
    for (unsigned int BXCounter = 0; BXCounter < evt_BX_for_PU.GetSize(); ++BXCounter) {
      int bx = evt_BX_for_PU[BXCounter];
      if (bx == 0) {
	eventPU = evt_PU[BXCounter];
	break;
      }
    }
    assert(eventPU > 0.);

    if (inputHistogram_data->GetBinContent(inputHistogram_data->FindFixBin(eventPU)) > 0.) {
      inputHistogram_MC->Fill(eventPU);
    }
    ++entryIndex;
  }
  std::cout << std::endl;
}

double getWeightedTotalNEntries(TH1D* inputHistogram) {
  double sumEntries = 0.;
  for (int binCounter = 0; binCounter <= (1+inputHistogram->GetXaxis()->GetNbins()); ++binCounter) {
    sumEntries += inputHistogram->GetBinContent(binCounter);
  }
  return sumEntries;
}

void normalizeHistogram(TH1D* inputHistogram) {
  double totalNEntries = getWeightedTotalNEntries(inputHistogram);
  std::cout << "For inputHistogram with name: " << inputHistogram->GetName() << ", weighted total nEntries: " << totalNEntries << std::endl;
  for (int binCounter = 0; binCounter <= (1+inputHistogram->GetXaxis()->GetNbins()); ++binCounter) {
    double originalContent = inputHistogram->GetBinContent(binCounter);
    double originalError = inputHistogram->GetBinError(binCounter);
    inputHistogram->SetBinContent(binCounter, originalContent/totalNEntries);
    if (originalContent > 0.) {
      inputHistogram->SetBinError(binCounter, originalError/totalNEntries);
    }
    else {
      inputHistogram->SetBinError(binCounter, 0.);
    }
  }
}

// void printHistogram(TH1D* inputHistogram) {
//   for (int binCounter = 0; binCounter <= (1+inputHistogram->GetXaxis()->GetNbins()); ++binCounter) {
//     std::cout << "At bin index = " << binCounter << " (" << inputHistogram->GetXaxis()->GetBinLowEdge(binCounter) << " < x < " << inputHistogram->GetXaxis()->GetBinUpEdge(binCounter) << "), content: " << inputHistogram->GetBinContent(binCounter) << std::endl;
//   }
// }

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  std::cout << "Starting code to make histograms with pileup weights..." << std::endl;
  std::cout << "Current working directory: " << tmMiscUtils::getCWD() << std::endl;

  tmArgumentParser argumentParser = tmArgumentParser("Find histogram of pileup weights to apply.");
  argumentParser.addArgument("inputDataPath", "", true, "Path to ROOT file containing a histogram of the pileup in data.");
  argumentParser.addArgument("inputMCPath", "", true, "Path to file containing MC ntuples.");
  argumentParser.addArgument("outputFolder", "root://cmseos.fnal.gov//store/user/lpcsusystealth/analysisEOSAreas/analysis", false, "Output folder.");
  argumentParser.addArgument("outputFileName", "", true, "Name of output file.");
  argumentParser.setPassedStringValues(argc, argv);
  argumentsStruct arguments = getArgumentsFromParser(argumentParser);

  TH1D* inputHistogram_data;
  TFile* inputDataFile = TFile::Open((arguments.inputDataPath).c_str(), "READ");
  assert((inputDataFile->IsOpen() && !(inputDataFile->IsZombie())));
  inputDataFile->GetObject("pileup", inputHistogram_data);
  // sanity checks
  assert((inputHistogram_data != nullptr));
  assert((inputHistogram_data->GetXaxis()->GetNbins() == 100));
  assert((std::fabs(inputHistogram_data->GetXaxis()->GetXmin()) < 0.1));
  assert((std::fabs(inputHistogram_data->GetXaxis()->GetXmax() - 100.0) < 0.1));
  normalizeHistogram(inputHistogram_data);

  TH1D* inputHistogram_MC = new TH1D("pileup_MC", "pileup_MC", 100, 0., 100.);
  fillMCHistogramBinsThatAreNonzeroInDataFromFile(inputHistogram_MC, inputHistogram_data, arguments.inputMCPath);
  normalizeHistogram(inputHistogram_MC);

  TH1D* puWeightsHistogram = new TH1D("pileupWeights", "pileup weights", 100, 0., 100.);
  for (int binCounter = 0; binCounter <= (1+puWeightsHistogram->GetXaxis()->GetNbins()); ++binCounter) {
    double denominator = inputHistogram_MC->GetBinContent(binCounter);
    if (denominator > 0.) {
      double numerator = inputHistogram_data->GetBinContent(binCounter);
      double numeratorError = inputHistogram_data->GetBinError(binCounter);
      double denominatorError = inputHistogram_MC->GetBinError(binCounter);
      double weight = numerator/denominator;
      puWeightsHistogram->SetBinContent(binCounter, weight);
      if (numerator > 0.) {
	puWeightsHistogram->SetBinError(binCounter, weight*std::sqrt(std::pow(numeratorError/numerator, 2) + std::pow(denominatorError/denominator, 2)));
      }
      else {
	puWeightsHistogram->SetBinError(binCounter, 0.);
      }
    }
    else {
      puWeightsHistogram->SetBinContent(binCounter, 0.);
      puWeightsHistogram->SetBinError(binCounter, 0.);
    }
  }
  TFile *outputFile = TFile::Open(("~/cmslpc_scratch/puWeights/" + arguments.outputFileName).c_str(), "RECREATE");
  outputFile->WriteTObject(puWeightsHistogram);
  outputFile->Close();
  inputDataFile->Close();

  int xrdcp_return_status = system(("set -x && xrdcp --verbose --force --path --streams 15 ~/cmslpc_scratch/puWeights/" + arguments.outputFileName + " " + arguments.outputFolder + "/" + arguments.outputFileName + " && rm -f ~/cmslpc_scratch/puWeights/" + arguments.outputFileName + " && set +x").c_str());
  if (xrdcp_return_status != 0) {
    std::cout << "ERROR: xrdcp likely failed with status "<< xrdcp_return_status << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::cout << "PU weights histogram saved." << std::endl;
  return 0;
}
