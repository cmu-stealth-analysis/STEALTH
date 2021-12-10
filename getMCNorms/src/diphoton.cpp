#include "../include/diphoton.h"

std::string get_output_st_distribution_name(const int & nJetsBin, const bool & STSwitch_fineBinned) {
  std::string fineBinned = "";
  if (STSwitch_fineBinned) fineBinned = "fineBinned_";
  return (std::string("ST_") + fineBinned + std::to_string(nJetsBin) + std::string("JetsBin"));
}

void initialize_output_th1s_map(std::map<std::string, TH1D> & output_th1s, const STRegionsStruct & STRegions, const STRegionsStruct & STRegionsFineBinned) {
  std::string hname;
  hname = std::string(DIPH_NJETS_TH1_NAME);
  assert(output_th1s.find(hname) == output_th1s.end());
  output_th1s[hname] = TH1D(hname.c_str(), ("nJets bin, " + std::to_string(DIPH_STNORM_MIN) + " GeV < ST < " + std::to_string(DIPH_STNORM_MAX) + " GeV;nJets bin;nEvts").c_str(), 5, 1.5, 6.5);
  output_th1s[hname].Sumw2();
  hname = std::string(DIPH_INVMASS_TH1_NAME);
  assert(output_th1s.find(hname) == output_th1s.end());
  output_th1s[hname] = TH1D(hname.c_str(), "invariant mass;m;nEvts/GeV", DIPH_INVMASS_NBINS, DIPH_INVMASS_MIN, DIPH_INVMASS_MAX);
  output_th1s[hname].Sumw2();
  for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
    hname = get_output_st_distribution_name(nJetsBin, false);
    assert(output_th1s.find(hname) == output_th1s.end());
    output_th1s[hname] = TH1D(hname.c_str(), ("ST, " + std::to_string(nJetsBin) + " Jets").c_str(), (STRegions.STBoundaries.size()-1), &(STRegions.STBoundaries.at(0)));
    output_th1s[hname].Sumw2();
    hname = get_output_st_distribution_name(nJetsBin, true);
    assert(output_th1s.find(hname) == output_th1s.end());
    output_th1s[hname] = TH1D(hname.c_str(), ("ST, fine-binned, " + std::to_string(nJetsBin) + " Jets").c_str(), (STRegionsFineBinned.STBoundaries.size()-1), &(STRegionsFineBinned.STBoundaries.at(0)));
    output_th1s[hname].Sumw2();
  }
}

void setup_chain(TChain * inputChain, eventDataStruct & event_data, const bool & addMCWeights) {
  inputChain->SetBranchStatus("*", 0);
  inputChain->SetBranchStatus("b_evtST", 1);
  inputChain->SetBranchAddress("b_evtST", &(event_data.evtST));
  inputChain->SetBranchStatus("b_invMass", 1);
  inputChain->SetBranchAddress("b_invMass", &(event_data.invMass));
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

bool passes_selection1(eventDataStruct & event_data) {
  return ((event_data.nJetsDR >= 2) &&
          (event_data.evtST >= DIPH_STNORM_MIN) &&
          (event_data.evtST <= DIPH_STNORM_MAX));
}

bool passes_selection2(eventDataStruct & event_data) {
  if (event_data.nJetsDR > 0) return false;
  bool leading_is_tight = static_cast<bool>((((event_data.phoID)->at(event_data.photonIndex_leading))>>2)&1);
  bool subLeading_is_tight = static_cast<bool>((((event_data.phoID)->at(event_data.photonIndex_subLeading))>>2)&1);
  return ((leading_is_tight && subLeading_is_tight) && ((event_data.invMass > DIPH_INVMASS_MIN) && (event_data.invMass < DIPH_INVMASS_MAX)));
}

bool passes_selection3(eventDataStruct & event_data, const double & STNormRangeMin, const double & STFineBinnedNormRangeMin) {
  if (event_data.nJetsDR < 2) return false;
  if ((event_data.evtST < STNormRangeMin) && (event_data.evtST < STFineBinnedNormRangeMin)) return false;
  return true;
}

void fill_histograms1(eventDataStruct & event_data, std::map<std::string, TH1D> & output_th1s, const bool & addMCWeights) {
  int nJetsBin = ((event_data.nJetsDR) <= 6) ? (event_data.nJetsDR) : 6;
  std::string hname = std::string(DIPH_NJETS_TH1_NAME);
  double bin_width = (output_th1s.at(hname)).GetXaxis()->GetBinWidth((output_th1s.at(hname)).GetXaxis()->FindFixBin(nJetsBin));
  assert(std::fabs(bin_width - 1.0) < DIPH_TOLERANCE_SANITY_CHECK);
  double weight = 1.0/bin_width;
  if (addMCWeights) {
    weight *= ((event_data.MCXSecWeight)*(event_data.MCGenWeight)*(event_data.prefiringWeight)*(event_data.photonMCScaleFactor)*(event_data.MCPUWeight));
  }
  (output_th1s.at(hname)).Fill(nJetsBin, weight);
}

void fill_histograms2(eventDataStruct & event_data, std::map<std::string, TH1D> & output_th1s, const bool & addMCWeights) {
  assert(event_data.nJetsDR == 0);
  std::string hname = std::string(DIPH_INVMASS_TH1_NAME);
  double bin_width = (output_th1s.at(hname)).GetXaxis()->GetBinWidth((output_th1s.at(hname)).GetXaxis()->FindFixBin(event_data.invMass));
  double weight = 1.0/bin_width;
  if (addMCWeights) {
    weight *= ((event_data.MCXSecWeight)*(event_data.MCGenWeight)*(event_data.prefiringWeight)*(event_data.photonMCScaleFactor)*(event_data.MCPUWeight));
  }
  (output_th1s.at(hname)).Fill(event_data.invMass, weight);
}

void fill_histograms3(eventDataStruct & event_data, std::map<std::string, TH1D> & output_th1s, const bool & addMCWeights, const double & STNormRangeMin, const double & STFineBinnedNormRangeMin) {
  assert(event_data.nJetsDR >= 2);
  int nJetsBin = ((event_data.nJetsDR) <= 6) ? (event_data.nJetsDR) : 6;
  std::string hname;
  double bin_width;
  double weight;

  if (event_data.evtST >= STNormRangeMin) {
    hname = get_output_st_distribution_name(nJetsBin, false);
    bin_width = (output_th1s.at(hname)).GetXaxis()->GetBinWidth((output_th1s.at(hname)).GetXaxis()->FindFixBin(event_data.evtST));
    weight = 1.0/bin_width;
    if (addMCWeights) {
      weight *= ((event_data.MCXSecWeight)*(event_data.MCGenWeight)*(event_data.prefiringWeight)*(event_data.photonMCScaleFactor)*(event_data.MCPUWeight));
    }
    (output_th1s.at(hname)).Fill(event_data.evtST, weight);
  }

  if (event_data.evtST >= STFineBinnedNormRangeMin) {
    hname = get_output_st_distribution_name(nJetsBin, true);
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
    if (passes_selection1(event_data)) fill_histograms1(event_data, output_th1s, addMCWeights);
    if (passes_selection2(event_data)) fill_histograms2(event_data, output_th1s, addMCWeights);
    if (passes_selection3(event_data, STRegions.STNormRangeMin, STRegionsFineBinned.STNormRangeMin)) fill_histograms3(event_data, output_th1s, addMCWeights, STRegions.STNormRangeMin, STRegionsFineBinned.STNormRangeMin);
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
  initialize_output_th1s_map(output_th1s, arguments.STRegions, arguments.STRegionsFineBinned);
  eventDataStruct event_data;
  setup_chain(inputChain, event_data, arguments.addMCWeights);
  loop_over_chain_events(inputChain, event_data, output_th1s, arguments.addMCWeights, arguments.STRegions, arguments.STRegionsFineBinned);
  common::write_output_th1s_to_file(std::string("~/cmslpc_scratch/MCNormsTemp/") + arguments.outputFileName, output_th1s);
  common::move_via_xrdcp("~/cmslpc_scratch/MCNormsTemp/" + arguments.outputFileName, arguments.outputFolder + "/" + arguments.outputFileName);
  std::cout << "All done." << std::endl;
  return EXIT_SUCCESS;
}
