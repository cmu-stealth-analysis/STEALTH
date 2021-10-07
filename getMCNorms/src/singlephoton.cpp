#include "../include/singlephoton.h"

std::string get_output_th1_name(const int & nJetsBin) {
  return (std::string("pT_leadingPhoton_") + std::to_string(nJetsBin) + std::string("JetsBin"));
}

void initialize_output_th1s_map(std::map<std::string, TH1D> & output_th1s) {
  for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
    std::string hname = get_output_th1_name(nJetsBin);
    assert(output_th1s.find(hname) == output_th1s.end());
    output_th1s[hname] = TH1D(hname.c_str(), ("pT of leading photon, " + std::to_string(nJetsBin) + " Jets Bin" ";pT;nEvts/bin").c_str(), 50, 200., 700.);
  }
}

void setup_chain(TChain * inputChain, eventDataStruct & event_data, const bool & addMCWeights) {
  inputChain->SetBranchStatus("*", 0);
  inputChain->SetBranchStatus("b_photonPT_leading", 1);
  inputChain->SetBranchAddress("b_photonPT_leading", &(event_data.pT_leadingPhoton));
  inputChain->SetBranchStatus("b_nJetsDR", 1);
  inputChain->SetBranchAddress("b_nJetsDR", &(event_data.nJetsDR));
  if (addMCWeights) {
    inputChain->SetBranchStatus("b_MCXSecWeight", 1);
    inputChain->SetBranchAddress("b_MCXSecWeight", &(event_data.MCXSecWeight));
    inputChain->SetBranchStatus("b_evtPrefiringWeight", 1);
    inputChain->SetBranchAddress("b_evtPrefiringWeight", &(event_data.prefiringWeight));
    inputChain->SetBranchStatus("b_evtphotonMCScaleFactor", 1);
    inputChain->SetBranchAddress("b_evtphotonMCScaleFactor", &(event_data.photonMCScaleFactor));
  }
}

bool passes_selection(eventDataStruct & event_data) {
  return ((event_data.nJetsDR >= 2) &&
	  (event_data.pT_leadingPhoton >= 200.));
}

void fill_histograms(eventDataStruct & event_data, std::map<std::string, TH1D> & output_th1s, const bool & addMCWeights) {
  int nJetsBin = ((event_data.nJetsDR) <= 6) ? (event_data.nJetsDR) : 6;
  std::string hname = get_output_th1_name(nJetsBin);
  double bin_width = (output_th1s.at(hname)).GetXaxis()->GetBinWidth((output_th1s.at(hname)).GetXaxis()->FindFixBin(event_data.pT_leadingPhoton));
  double weight = 1.0/bin_width;
  if (addMCWeights) weight *= ((event_data.MCXSecWeight)*(event_data.prefiringWeight)*(event_data.photonMCScaleFactor));
  (output_th1s.at(hname)).Fill(event_data.pT_leadingPhoton, weight);
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
  std::cout << "Saving distributions of leading photon pT..." << std::endl;
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