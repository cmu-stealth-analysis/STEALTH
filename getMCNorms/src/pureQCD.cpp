#include "../include/pureQCD.h"

std::string get_output_th1_name(const int & nJetsBin, const bool & switch_ST, const bool & STSwitch_fineBinned) {
  if (switch_ST) {
    std::string fineBinned = "";
    if (STSwitch_fineBinned) fineBinned = "fineBinned_";
    return (std::string("ST_") + fineBinned + std::to_string(nJetsBin) + std::string("JetsBin"));
  }
  return (std::string("pT_leadingJet_") + std::to_string(nJetsBin) + std::string("JetsBin"));
}

void initialize_output_th1s_map(std::map<std::string, TH1D> & output_th1s, const STRegionsStruct & STRegions, const STRegionsStruct & STRegionsFineBinned) {
  for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
    std::string hname;
    hname = get_output_th1_name(nJetsBin, false, false);
    assert(output_th1s.find(hname) == output_th1s.end());
    output_th1s[hname] = TH1D(hname.c_str(), ("pT of leading jet, " + std::to_string(nJetsBin) + " Jets Bin;pT;nEvts/GeV").c_str(), 80, 200., 1000.);
    output_th1s[hname].Sumw2();
    hname = get_output_th1_name(nJetsBin, true, false);
    assert(output_th1s.find(hname) == output_th1s.end());
    output_th1s[hname] = TH1D(hname.c_str(), ("ST, " + std::to_string(nJetsBin) + " Jets Bin;ST;nEvts/GeV").c_str(), (STRegions.STBoundaries.size()-1), &(STRegions.STBoundaries.at(0)));
    output_th1s[hname].Sumw2();
    hname = get_output_th1_name(nJetsBin, true, true);
    assert(output_th1s.find(hname) == output_th1s.end());
    output_th1s[hname] = TH1D(hname.c_str(), ("ST, fine-binned," + std::to_string(nJetsBin) + " Jets Bin;ST;nEvts/GeV").c_str(), (STRegionsFineBinned.STBoundaries.size()-1), &(STRegionsFineBinned.STBoundaries.at(0)));
    output_th1s[hname].Sumw2();
  }
}

void setup_chain(TChain * inputChain, eventDataStruct & event_data, const bool & addMCWeights) {
  inputChain->SetBranchStatus("*", 0);
  inputChain->SetBranchStatus("b_evtST", 1);
  inputChain->SetBranchAddress("b_evtST", &(event_data.evtST));
  inputChain->SetBranchStatus("b_jetPT_leading", 1);
  inputChain->SetBranchAddress("b_jetPT_leading", &(event_data.pT_leadingJet));
  inputChain->SetBranchStatus("b_nJetsDR", 1);
  inputChain->SetBranchAddress("b_nJetsDR", &(event_data.nJetsDR));
  inputChain->SetBranchStatus("b_nJetsAll", 1);
  inputChain->SetBranchAddress("b_nJetsAll", &(event_data.nJetsAll));
  inputChain->SetBranchStatus("phoIDbit", 1);
  inputChain->SetBranchAddress("phoIDbit", &(event_data.phoID));
  inputChain->SetBranchStatus("b_photonIndex_leading", 1);
  inputChain->SetBranchAddress("b_photonIndex_leading", &(event_data.photonIndex_leading));
  inputChain->SetBranchStatus("b_photonIndex_subLeading", 1);
  inputChain->SetBranchAddress("b_photonIndex_subLeading", &(event_data.photonIndex_subLeading));
  if (addMCWeights) {
    inputChain->SetBranchStatus("b_MCXSecWeight", 1);
    inputChain->SetBranchAddress("b_MCXSecWeight", &(event_data.MCXSecWeight));
    inputChain->SetBranchStatus("genWeight", 1);
    inputChain->SetBranchAddress("genWeight", &(event_data.MCGenWeight));
    inputChain->SetBranchStatus("b_evtPrefiringWeight", 1);
    inputChain->SetBranchAddress("b_evtPrefiringWeight", &(event_data.prefiringWeight));
    inputChain->SetBranchStatus("b_evtphotonMCScaleFactor", 1);
    inputChain->SetBranchAddress("b_evtphotonMCScaleFactor", &(event_data.photonMCScaleFactor));
    inputChain->SetBranchStatus("b_PUWeightNoSelection", 1);
    inputChain->SetBranchAddress("b_PUWeightNoSelection", &(event_data.MCPUWeight));
  }
}

bool passes_selection(eventDataStruct & event_data) {
  bool leading_is_medium = false;
  if (event_data.photonIndex_leading >= 0) leading_is_medium = static_cast<bool>((((event_data.phoID)->at(event_data.photonIndex_leading))>>1)&1);
  bool subLeading_is_medium = false;
  if (event_data.photonIndex_subLeading >= 0) subLeading_is_medium = static_cast<bool>((((event_data.phoID)->at(event_data.photonIndex_subLeading))>>1)&1);
  if (leading_is_medium && subLeading_is_medium) return false;
  return (event_data.nJetsDR >= 2);
}

void fill_histograms(eventDataStruct & event_data, std::map<std::string, TH1D> & output_th1s, const bool & addMCWeights, const double & STNormRangeMin, const double & STFineBinnedNormRangeMin) {
  int nJetsBin = ((event_data.nJetsDR) <= 6) ? (event_data.nJetsDR) : 6;
  std::string hname;
  double bin_width;
  double weight;

  if (event_data.pT_leadingJet >= 200.) {
    hname = get_output_th1_name(nJetsBin, false, false);
    bin_width = (output_th1s.at(hname)).GetXaxis()->GetBinWidth((output_th1s.at(hname)).GetXaxis()->FindFixBin(event_data.pT_leadingJet));
    weight = 1.0/bin_width;
    if (addMCWeights) {
      weight *= ((event_data.MCXSecWeight)*(event_data.MCGenWeight)*(event_data.prefiringWeight)*(event_data.photonMCScaleFactor)*(event_data.MCPUWeight));
    }
    (output_th1s.at(hname)).Fill(event_data.pT_leadingJet, weight);
  }

  if (static_cast<double>(event_data.evtST) >= STNormRangeMin) {
    hname = get_output_th1_name(nJetsBin, true, false);
    bin_width = (output_th1s.at(hname)).GetXaxis()->GetBinWidth((output_th1s.at(hname)).GetXaxis()->FindFixBin(event_data.evtST));
    weight = 1.0/bin_width;
    if (addMCWeights) {
      weight *= ((event_data.MCXSecWeight)*(event_data.MCGenWeight)*(event_data.prefiringWeight)*(event_data.photonMCScaleFactor)*(event_data.MCPUWeight));
    }
    (output_th1s.at(hname)).Fill(event_data.evtST, weight);
  }
  if (static_cast<double>(event_data.evtST) >= STFineBinnedNormRangeMin) {
    hname = get_output_th1_name(nJetsBin, true, true);
    bin_width = (output_th1s.at(hname)).GetXaxis()->GetBinWidth((output_th1s.at(hname)).GetXaxis()->FindFixBin(event_data.evtST));
    weight = 1.0/bin_width;
    if (addMCWeights) {
      weight *= ((event_data.MCXSecWeight)*(event_data.MCGenWeight)*(event_data.prefiringWeight)*(event_data.photonMCScaleFactor)*(event_data.MCPUWeight));
    }
    (output_th1s.at(hname)).Fill(event_data.evtST, weight);
  }
}

void loop_over_chain_events(TChain * inputChain, eventDataStruct & event_data, std::map<std::string, TH1D> & output_th1s, const bool & addMCWeights, const STRegionsStruct & STRegions, const STRegionsStruct & STRegionsFineBinned) {
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
    if (passes_selection(event_data)) fill_histograms(event_data, output_th1s, addMCWeights, STRegions.STNormRangeMin, STRegionsFineBinned.STNormRangeMin);
  }
  progressBar.terminate();
}

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  TH1::AddDirectory(kFALSE);
  std::cout << "Saving distributions of leading jet pT..." << std::endl;
  common::argumentsStruct arguments = common::get_command_line_arguments(argc, argv);
  TChain * inputChain = common::get_chain_from_input_paths_files(arguments.inputPathsFiles);
  std::map<std::string, TH1D> output_th1s;
  initialize_output_th1s_map(output_th1s, arguments.STRegions, arguments.STRegionsFineBinned);
  eventDataStruct event_data;
  setup_chain(inputChain, event_data, arguments.addMCWeights);
  loop_over_chain_events(inputChain, event_data, output_th1s, arguments.addMCWeights, arguments.STRegions, arguments.STRegionsFineBinned);
  common::write_output_th1s_to_file(std::string("~/cmslpc_scratch/MCNormsTemp/") + arguments.outputFileName, output_th1s);
  common::move_via_xrdcp("~/cmslpc_scratch/MCNormsTemp/" + arguments.outputFileName, arguments.outputFolder + "/" + arguments.outputFileName);
  std::cout << "All done." << std::endl;
  return EXIT_SUCCESS;
}
