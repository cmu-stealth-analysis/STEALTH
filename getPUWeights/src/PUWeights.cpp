#include "../include/PUWeights.h"

void fillMCHistogramBinsThatAreNonzeroInDataFromFile(TH1D* inputHistogram_MC, TH1D* inputHistogram_data, const std::vector<std::string> & inputMCPaths, const bool& fetchMCWeights, const bool& addMCXSecWeight) {
  TChain *inputChain = new TChain("ggNtuplizer/EventTree");
  for (const std::string & inputMCPath : inputMCPaths) {
    std::cout << "Adding events from path: " << inputMCPath << std::endl;
    inputChain->Add(inputMCPath.c_str());
  }
  long nEntries = inputChain->GetEntries();
  assert((nEntries > 0));
  std::cout << "Number of available entries: " << nEntries << std::endl;

  inputChain->SetBranchStatus("*", 0); // so that only the needed branches, explicitly activated below, are read in per event

  std::vector<int> * evt_BX_for_PU = nullptr;
  inputChain->SetBranchStatus("puBX", 1);
  inputChain->SetBranchAddress("puBX", &(evt_BX_for_PU));

  std::vector<float> * evt_PU = nullptr;
  inputChain->SetBranchStatus("puTrue", 1);
  inputChain->SetBranchAddress("puTrue", &(evt_PU));

  double MCXSecWeight = -1.;
  float MCGenWeight = -1.;
  if (addMCXSecWeight) {
    inputChain->SetBranchStatus("b_MCXSecWeight", 1);
    inputChain->SetBranchAddress("b_MCXSecWeight", &(MCXSecWeight));
    inputChain->SetBranchStatus("genWeight", 1);
    inputChain->SetBranchAddress("genWeight", &(MCGenWeight));
  }

  float MCPrefiringWeight = -1.;
  float MCScaleFactorWeight = -1.;
  if (fetchMCWeights) {
    inputChain->SetBranchStatus("b_evtPrefiringWeight", 1);
    inputChain->SetBranchAddress("b_evtPrefiringWeight", &(MCPrefiringWeight));
    inputChain->SetBranchStatus("b_evtphotonMCScaleFactor", 1);
    inputChain->SetBranchAddress("b_evtphotonMCScaleFactor", &(MCScaleFactorWeight));
  }
    
  tmProgressBar *progressBar = new tmProgressBar(nEntries);
  int tmp = static_cast<int>(0.5 + 1.0*nEntries/50);
  int progressBarUpdatePeriod = tmp > 1 ? tmp : 1;
  progressBar->initialize();

  for (Long64_t entryIndex = 0; entryIndex < nEntries; ++entryIndex) {
    Long64_t loadStatus = inputChain->LoadTree(entryIndex);
    assert(loadStatus >= 0);
    int nBytesRead = inputChain->GetEntry(entryIndex, 0); // Get only the required branches
    assert(nBytesRead > 0);
    if ((entryIndex > 0) && (((entryIndex % static_cast<Long64_t>(progressBarUpdatePeriod)) == 0) || (entryIndex == (nEntries-1)))) progressBar->updateBar(static_cast<double>(1.0*entryIndex/nEntries), entryIndex);

    float eventPU = -1.;
    for (unsigned int BXCounter = 0; BXCounter < static_cast<unsigned int>((*evt_BX_for_PU).size()); ++BXCounter) {
      int bx = (*evt_BX_for_PU).at(BXCounter);
      if (bx == 0) {
	eventPU = (*evt_PU).at(BXCounter);
	break;
      }
    }
    assert(eventPU > 0.);

    double eventWeight = 1.0;
    if (fetchMCWeights) eventWeight *= (MCPrefiringWeight*MCScaleFactorWeight);
    if (addMCXSecWeight) {
      assert(MCXSecWeight > 0.);
      eventWeight *= (MCXSecWeight*MCGenWeight);
    }

    if (inputHistogram_data->GetBinContent(inputHistogram_data->FindFixBin(eventPU)) > 0.) {
      inputHistogram_MC->Fill(eventPU, eventWeight);
    }
  }
  std::cout << std::endl;
}

double getWeightedTotalNEntries(TH1D* inputHistogram, const bool& countUnderflow, const bool& countOverflow) {
  double sumEntries = 0.;
  if (countUnderflow) sumEntries += inputHistogram->GetBinContent(0);
  for (int binCounter = 1; binCounter <= (inputHistogram->GetXaxis()->GetNbins()); ++binCounter) {
    sumEntries += inputHistogram->GetBinContent(binCounter);
  }
  if (countOverflow) sumEntries += inputHistogram->GetBinContent(1+(inputHistogram->GetXaxis()->GetNbins()));
  return sumEntries;
}

void normalizeHistogram(TH1D* inputHistogram) {
  double totalNEntries = getWeightedTotalNEntries(inputHistogram, false, false);
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
  argumentParser.addArgument("inputMCPaths", "", true, "Path to files containing MC ntuples, separated by the character \",\".");
  argumentParser.addArgument("outputFolder", "root://cmseos.fnal.gov//store/user/lpcsusystealth/analysisEOSAreas/analysis", false, "Output folder.");
  argumentParser.addArgument("outputFileName", "", true, "Name of output file.");
  argumentParser.addArgument("addMCXSecWeight", "false", true, "If this argument is set, then relative weights are read in from an additional branch, used for GJet MC samples.");
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
  inputHistogram_MC->Sumw2();
  fillMCHistogramBinsThatAreNonzeroInDataFromFile(inputHistogram_MC, inputHistogram_data, arguments.inputMCPaths, true, arguments.addMCXSecWeight);
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

  int xrdcp_return_status = system(("set -x && xrdcp --nopbar --silent --force --path --streams 15 ~/cmslpc_scratch/puWeights/" + arguments.outputFileName + " " + arguments.outputFolder + "/" + arguments.outputFileName + " && rm -f ~/cmslpc_scratch/puWeights/" + arguments.outputFileName + " && set +x").c_str());
  if (xrdcp_return_status != 0) {
    std::cout << "ERROR: xrdcp likely failed with status "<< xrdcp_return_status << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::cout << "PU weights histogram saved." << std::endl;
  return 0;
}
