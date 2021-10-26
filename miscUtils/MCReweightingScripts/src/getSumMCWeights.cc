#include "../include/getSumMCWeights.h"

TChain * get_chain_from_input_paths_files(const std::vector<std::string> & inputPathsFiles) {
  TChain * inputChain = new TChain("ggNtuplizer/EventTree");
  inputChain->SetMaxTreeSize(10000000000000LL); // 10 TB
  for (const std::string & inputPathsFile : inputPathsFiles) {
    std::cout << "Adding paths from file: " << inputPathsFile << std::endl;
    std::ifstream inputPathsFileStream;
    inputPathsFileStream.open(inputPathsFile.c_str());
    assert(inputPathsFileStream.is_open());
    while (!(inputPathsFileStream.eof())) {
      std::string inputPath;
      inputPathsFileStream >> inputPath;
      if (!(inputPath.empty())) {
	if (!(file_has_zero_events(inputPath))) inputChain->Add(inputPath.c_str());
      }
    }
    inputPathsFileStream.close();
  }
  return inputChain;
}

void setup_chain(TChain * inputChain, eventInfoStruct & event_info) {
  inputChain->SetBranchStatus("*", 0);
  inputChain->SetBranchStatus("genWeight", 1);
  inputChain->SetBranchAddress("genWeight", &(event_info.MCWeight));
  inputChain->SetBranchStatus("puBX", 1);
  inputChain->SetBranchAddress("puBX", &(event_info.BX_for_PU));
  inputChain->SetBranchStatus("puTrue", 1);
  inputChain->SetBranchAddress("puTrue", &(event_info.PUTrue));
}

void loop_over_chain_events(TChain * inputChain, eventInfoStruct & event_info, outputInfoStruct & output_info) {
  int n_events_without_pu_info = 0;
  long nEntries = inputChain->GetEntries();
  tmProgressBar progressBar(nEntries);
  int tmp = static_cast<int>(0.5 + 1.0*nEntries/20);
  int progressBarUpdatePeriod = tmp > 1 ? tmp : 1;
  progressBar.initialize();
  for (Long64_t entryIndex = 0; entryIndex < nEntries; ++entryIndex) {
    Long64_t loadStatus = inputChain->LoadTree(entryIndex);
    assert(loadStatus >= 0);
    int nBytesRead = inputChain->GetEntry(entryIndex, 0); // Get only the required branches
    assert(nBytesRead > 0);
    if ((entryIndex == 0) ||
        (entryIndex == (nEntries-1)) ||
        ((entryIndex % static_cast<Long64_t>(progressBarUpdatePeriod)) == 0)) progressBar.updateBar(static_cast<double>(1.0*entryIndex/nEntries), entryIndex);
    ++(output_info.totalNEvts);

    float eventPU = -1.0;
    for (size_t BXCounter = 0; BXCounter < (*(event_info.BX_for_PU)).size(); ++BXCounter) {
      int bx = (*(event_info.BX_for_PU)).at(BXCounter);
      if (bx == 0) {
	eventPU = (*(event_info.PUTrue)).at(BXCounter);
	break;
      }
    }
    bool PUFound = ((eventPU >= PU_MINVAL) && (eventPU <= PU_MAXVAL));
    if (!(PUFound)) {
      ++n_events_without_pu_info;
      if (static_cast<double>(1.0*entryIndex) > static_cast<double>(1.0/MAX_FRAC_EVENTS_WITHOUT_PU_INFO)) {
	assert(static_cast<double>(1.0*n_events_without_pu_info/entryIndex) < static_cast<double>(MAX_FRAC_EVENTS_WITHOUT_PU_INFO));
      }
    }
    if (!(PUFound)) continue;
    ++(output_info.totalNEvtsWithPUInfo);
    (output_info.sumWeights) += (event_info.MCWeight);
    (output_info.pu_MC).Fill(eventPU, event_info.MCWeight);
  }
  if (n_events_without_pu_info > 0) {
    std::cout << "WARNING: " << n_events_without_pu_info << "/" << nEntries << " events (" << ((100.0*n_events_without_pu_info)/(1.0*nEntries)) << "%) without PU info." << std::endl;
  }
  progressBar.terminate();
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
  double weightedTotalNEntries = getWeightedTotalNEntries(inputHistogram, false, false);
  std::cout << "For inputHistogram with name: " << inputHistogram->GetName() << ", weighted total nEntries: " << weightedTotalNEntries << std::endl;
  assert(weightedTotalNEntries > 0.);
  for (int binCounter = 0; binCounter <= (1+inputHistogram->GetXaxis()->GetNbins()); ++binCounter) {
    double originalContent = inputHistogram->GetBinContent(binCounter);
    double originalError = inputHistogram->GetBinError(binCounter);
    inputHistogram->SetBinContent(binCounter, originalContent/weightedTotalNEntries);
    if (originalContent > 0.) {
      inputHistogram->SetBinError(binCounter, originalError/weightedTotalNEntries);
    }
    else {
      inputHistogram->SetBinError(binCounter, 0.);
    }
  }
}

void save_PU_weights_into_histogram(TH1D * puWeightsHistogram, TH1D * PU_MC, TH1D * PU_data) {
  normalizeHistogram(PU_MC);
  normalizeHistogram(PU_data);
  for (int binCounter = 0; binCounter <= (1+puWeightsHistogram->GetXaxis()->GetNbins()); ++binCounter) {
    double denominator = PU_MC->GetBinContent(binCounter);
    if (denominator > 0.) {
      double numerator = PU_data->GetBinContent(binCounter);
      double numeratorError = PU_data->GetBinError(binCounter);
      double denominatorError = PU_MC->GetBinError(binCounter);
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
}

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  TH1::AddDirectory(kFALSE);
  std::cout << "Saving sum of weights and PU info into outputs..." << std::endl;
  argumentsStruct arguments = get_command_line_arguments(argc, argv);
  TChain * inputChain = get_chain_from_input_paths_files(arguments.inputPathsFiles);
  eventInfoStruct event_info;
  setup_chain(inputChain, event_info);
  outputInfoStruct output_info;
  (output_info.pu_MC).Sumw2();
  loop_over_chain_events(inputChain, event_info, output_info);
  write_weight_outputs_to_json_file(arguments.outputFileNameWeights, output_info);

  TH1D * PU_sourceData = new TH1D();
  TFile* inputPUSourceDataFileHandle = TFile::Open((arguments.dataPUSourceWithXRDPrefix).c_str(), "READ");
  assert((inputPUSourceDataFileHandle->IsOpen() && !(inputPUSourceDataFileHandle->IsZombie())));
  inputPUSourceDataFileHandle->GetObject("pileup", PU_sourceData);
  // sanity checks
  assert((PU_sourceData != nullptr));
  assert((PU_sourceData->GetXaxis()->GetNbins() == 100));
  assert((std::fabs(PU_sourceData->GetXaxis()->GetXmin()) < 0.1));
  assert((std::fabs(PU_sourceData->GetXaxis()->GetXmax() - 100.0) < 0.1));
  // sanity checks end

  TH1D * puWeightsHistogram = new TH1D("pileupWeights", "pileup weights", 100, 0., 100.);
  puWeightsHistogram->Sumw2();
  save_PU_weights_into_histogram(puWeightsHistogram, &(output_info.pu_MC), PU_sourceData);
  write_histogram_outputs_to_root_file(arguments.outputFileNamePU, &(output_info.pu_MC), puWeightsHistogram);
  inputPUSourceDataFileHandle->Close();

  std::cout << "All done." << std::endl;
  return EXIT_SUCCESS;
}
