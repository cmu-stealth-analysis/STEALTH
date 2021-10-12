#include "../include/diphoton.h"

// std::string get_output_th1_name() {
//   // return (std::string("pT_leadingPhoton_") + std::to_string(nJetsBin) + std::string("JetsBin"));
//   return std::string("nJets_in_normST");
// }

void initialize_output_th1s_map(std::map<std::string, TH1D> & output_th1s) {
  std::string hname = std::string(DIPH_TH1_NAME);
  assert(output_th1s.find(hname) == output_th1s.end());
  output_th1s[hname] = TH1D(hname.c_str(), ("nJets bin, " + std::to_string(DIPH_STNORM_MIN) + " GeV < ST < " + std::to_string(DIPH_STNORM_MAX) + " GeV;nJets bin;nEvts").c_str(), 5, 1.5, 6.5);
  output_th1s[hname].Sumw2();
}

void setup_chain(TChain * inputChain, eventDataStruct & event_data, const bool & addMCWeights) {
  inputChain->SetBranchStatus("*", 0);
  inputChain->SetBranchStatus("b_evtST", 1);
  inputChain->SetBranchAddress("b_evtST", &(event_data.evtST));
  inputChain->SetBranchStatus("b_nJetsDR", 1);
  inputChain->SetBranchAddress("b_nJetsDR", &(event_data.nJetsDR));
  // inputChain->SetBranchStatus("phoIDbit", 1);
  // inputChain->SetBranchAddress("phoIDbit", &(event_data.phoID));
  // inputChain->SetBranchStatus("b_photonIndex_leading", 1);
  // inputChain->SetBranchAddress("b_photonIndex_leading", &(event_data.photonIndex_leading));
  // inputChain->SetBranchStatus("b_photonIndex_subLeading", 1);
  // inputChain->SetBranchAddress("b_photonIndex_subLeading", &(event_data.photonIndex_subLeading));
  if (addMCWeights) {
    inputChain->SetBranchStatus("b_MCXSecWeight", 1);
    inputChain->SetBranchAddress("b_MCXSecWeight", &(event_data.MCXSecWeight));
    inputChain->SetBranchStatus("genWeight", 1);
    inputChain->SetBranchAddress("genWeight", &(event_data.MCGenWeight));
    inputChain->SetBranchStatus("b_evtPrefiringWeight", 1);
    inputChain->SetBranchAddress("b_evtPrefiringWeight", &(event_data.prefiringWeight));
    inputChain->SetBranchStatus("b_evtphotonMCScaleFactor", 1);
    inputChain->SetBranchAddress("b_evtphotonMCScaleFactor", &(event_data.photonMCScaleFactor));
  }
}

bool passes_selection(eventDataStruct & event_data) {
  // if (event_data.nJetsDR != 2) return false;
  // bool leading_is_tight = static_cast<bool>((((event_data.phoID)->at(event_data.photonIndex_leading))>>2)&1);
  // bool subLeading_is_tight = static_cast<bool>((((event_data.phoID)->at(event_data.photonIndex_subLeading))>>2)&1);
  // return (leading_is_tight && subLeading_is_tight);
  return ((event_data.nJetsDR >= 2) &&
          (event_data.evtST >= DIPH_STNORM_MIN) &&
          (event_data.evtST <= DIPH_STNORM_MAX));
}

void fill_histograms(eventDataStruct & event_data, std::map<std::string, TH1D> & output_th1s, const bool & addMCWeights) {
  int nJetsBin = ((event_data.nJetsDR) <= 6) ? (event_data.nJetsDR) : 6;
  std::string hname = std::string(DIPH_TH1_NAME);
  double bin_width = (output_th1s.at(hname)).GetXaxis()->GetBinWidth((output_th1s.at(hname)).GetXaxis()->FindFixBin(event_data.evtST));
  assert(std::fabs(bin_width - 1.0) < DIPH_TOLERANCE_SANITY_CHECK);
  double weight = 1.0/bin_width;
  if (addMCWeights) {
    weight *= ((event_data.MCXSecWeight)*(event_data.MCGenWeight)*(event_data.prefiringWeight)*(event_data.photonMCScaleFactor));
  }
  (output_th1s.at(hname)).Fill(nJetsBin, weight);
}

void loop_over_chain_events(TChain * inputChain, eventDataStruct & event_data, std::map<std::string, TH1D> & output_th1s, const bool & addMCWeights) {
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
    if (passes_selection(event_data)) fill_histograms(event_data, output_th1s, addMCWeights);
  }
  progressBar.terminate();
}

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  TH1::AddDirectory(kFALSE);
  std::cout << "Saving distributions of event ST..." << std::endl;
  common::argumentsStruct arguments = common::get_command_line_arguments(argc, argv);
  TChain * inputChain = common::get_chain_from_input_paths_files(arguments.inputPathsFiles);
  std::map<std::string, TH1D> output_th1s;
  initialize_output_th1s_map(output_th1s);
  eventDataStruct event_data;
  setup_chain(inputChain, event_data, arguments.addMCWeights);
  loop_over_chain_events(inputChain, event_data, output_th1s, arguments.addMCWeights);
  common::write_output_th1s_to_file(std::string("~/cmslpc_scratch/MCNormsTemp/") + arguments.outputFileName, output_th1s);
  common::move_via_xrdcp("~/cmslpc_scratch/MCNormsTemp/" + arguments.outputFileName, arguments.outputFolder + "/" + arguments.outputFileName);
  std::cout << "All done." << std::endl;
  return EXIT_SUCCESS;
}
