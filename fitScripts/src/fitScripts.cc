#include "../include/fitScripts.h"

template<typename T>
std::vector<T> getVectorCopyStartingFromIndex(const int& startIndex, const std::vector<T>& vin) {
  assert (startIndex <= static_cast<int>(vin.size()));
  std::vector<T> vout;
  for (int index=startIndex; index < static_cast<int>(vin.size()); ++index) {
    vout.push_back(vin.at(index));
  }
  return vout;
}

void printSeparator(const int& nCharactersToPrint=240) {
  std::cout << std::endl;
  for (int characterCounter = 0; characterCounter < nCharactersToPrint; ++characterCounter) std::cout << std::string("_");
  std::cout << std::endl;
  std::cout << std::endl;
}

double parseLineForFloatWithCheck(const std::string& inputLine, const std::string& targetName, const bool& print_debug=false) {
  if (print_debug) std::cout << "parseLineForFloatWithCheck called with inputLine: " << inputLine << ", targetName: " << targetName << std::endl;
  std::vector<std::string> components;
  std::stringstream runningComponent;
  const char spaceCharacter = ' ';
  const char equalSignCharacter = '=';
  for (unsigned int stringIndex = 0; stringIndex < static_cast<unsigned int>(inputLine.size()); ++stringIndex) {
    const char &character = inputLine.at(stringIndex);
    if ((character == spaceCharacter) || (character == equalSignCharacter)) {
      components.push_back(runningComponent.str());
      runningComponent.str(std::string());
      runningComponent.clear();
    }
    else {
      runningComponent << character;
    }
  }
  components.push_back(runningComponent.str());
  if (print_debug) {
    std::cout << "components: " << std::endl;
    printVector(components);
  }
  assert(static_cast<int>(components.size()) == 3);
  assert(components.at(0) == "float");
  assert(components.at(1) == targetName);
  return std::stod(components.at(2));
}

double get_fTest_prob(const double& chi2_1, const double& chi2_2, const int& ndf_1, const int& ndf_2, const std::string& id_string) {
  std::string outputTmpFileName = "get_fTest_prob_tmp_" + id_string + ".txt";

  // Step 1: Get the environment variable EOSTMPAREA, a folder where temp files can be stored
  std::string tmpFolder = std::string(getenv("EOSTMPAREA"));
  if (tmpFolder == "") {
    std::cout << "ERROR: env variable \"EOSTMPAREA\" does not appear to be set." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  
  // Step 2: Generate python command that gets this value
  std::string commandToRun = "python -c \"import tmStatsUtils; print(tmStatsUtils.get_fTest_prob(chi2_1=" + std::to_string(chi2_1) + ", chi2_2=" + std::to_string(chi2_2) + ", ndf_1=" + std::to_string(ndf_1) + ", ndf_2=" + std::to_string(ndf_2) + "))\" > " + tmpFolder + "/" + outputTmpFileName + " 2>&1";
  int ftest_command_return_status = system(commandToRun.c_str());
  if (ftest_command_return_status != 0) {
    std::cout << "ERROR in running command: " << commandToRun << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Step 3: Open output file and read in the value
  std::vector<double> valuesFromFile;
  double valueFromFile;
  std::ifstream inputFileObject((tmpFolder + "/" + outputTmpFileName).c_str());
  if (inputFileObject.is_open()) {
    while (inputFileObject >> valueFromFile) {
      valuesFromFile.push_back(valueFromFile);
    }
    inputFileObject.close();
  }
  else {
    std::cout << "ERROR: Unable to open file: " << (tmpFolder + "/" + outputTmpFileName) << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Step 4: Check that there is indeed just one value
  if (!(valuesFromFile.size() == 1)) {
    std::cout << "ERROR: this tmp file is in an unexpected format: " << (tmpFolder + "/" + outputTmpFileName) << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Step 5: Delete the temp file
  std::string file_remove_command = "rm " + (tmpFolder + "/" + outputTmpFileName);
  int rm_return_status = system(file_remove_command.c_str());
  if (rm_return_status != 0) {
    std::cout << "ERROR in running command: " << file_remove_command << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // finally return the value read in from step 3
  return (valuesFromFile.at(0));
}

double get_nll_prob(const double& nll_null, const double& nll_alt, const int& n_extra_parameters_in_alternative) {
  // Step 1: Get the environment variable EOSTMPAREA, a folder where temp files can be stored
  std::string tmpFolder = std::string(getenv("EOSTMPAREA"));
  if (tmpFolder == "") {
    std::cout << "ERROR: env variable \"EOSTMPAREA\" does not appear to be set." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Step 2: Generate python command that gets this value
  std::string commandToRun = "python -c \"import tmStatsUtils; print(tmStatsUtils.get_pVal_of_alternative_with_wilks(nll_null=" + std::to_string(nll_null) + ", nll_alternative=" + std::to_string(nll_alt) + ", n_extra_parameters_in_alternative=" + std::to_string(n_extra_parameters_in_alternative) + "))\" > " + tmpFolder + "/get_nll_prob_tmp.txt 2>&1";
  int ftest_command_return_status = system(commandToRun.c_str());
  if (ftest_command_return_status != 0) {
    std::cout << "ERROR in running command: " << commandToRun << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Step 3: Open output file and read in the value
  std::vector<double> valuesFromFile;
  double valueFromFile;
  std::ifstream inputFileObject((tmpFolder + "/get_nll_prob_tmp.txt").c_str());
  if (inputFileObject.is_open()) {
    while (inputFileObject >> valueFromFile) {
      valuesFromFile.push_back(valueFromFile);
    }
    inputFileObject.close();
  }
  else {
    std::cout << "ERROR: Unable to open file: " << (tmpFolder + "/get_nll_prob_tmp.txt") << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Step 4: Check that there is indeed just one value
  if (!(valuesFromFile.size() == 1)) {
    std::cout << "ERROR: this tmp file is in an unexpected format: " << (tmpFolder + "/get_nll_prob_tmp.txt") << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Step 5: Delete the temp file
  std::string file_remove_command = "rm " + (tmpFolder + "/get_nll_prob_tmp.txt");
  int rm_return_status = system(file_remove_command.c_str());
  if (rm_return_status != 0) {
    std::cout << "ERROR in running command: " << file_remove_command << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // finally return the value read in from step 3
  return (valuesFromFile.at(0));
}

TGraph get_TF1_as_TGraph(TF1* inputTF1, int nGraphPoints, double xMin, double xMax, double additionalScale=1.0) {
  TGraph outputGraph = TGraph(nGraphPoints);
  for (int xCounter = 0; xCounter <= nGraphPoints; ++xCounter) {
    double x = xMin + (1.0*xCounter/nGraphPoints)*(xMax-xMin);
    outputGraph.SetPoint(xCounter, x, additionalScale*(inputTF1->Eval(x)));
  }
  return outputGraph;
}

std::map<int, double> getNominalAdjustmentsFromBinIntegralMaps(const STRegionsStruct &regions, const std::map<int, double> &bin_integrals_divided_by_bin_widths_nominal, const std::map<int, double> &bin_integrals_divided_by_bin_widths_from_low_njets_kernel) {
  std::map<int, double> adjustments;
  for (int regionIndex = 1; regionIndex <= regions.STAxis.GetNbins(); ++regionIndex) {
    adjustments[regionIndex] = (bin_integrals_divided_by_bin_widths_nominal.at(regionIndex))/(bin_integrals_divided_by_bin_widths_from_low_njets_kernel.at(regionIndex));
  }
  return adjustments;
}

std::map<int, double> getAdjustmentFractionalCorrections(const STRegionsStruct &regions, const std::map<int, double> &bin_integrals_divided_by_bin_widths_nominal, const std::map<int, double> &bin_integrals_divided_by_bin_widths_shifted) {
  std::map<int, double> fractionalCorrections;
  for (int regionIndex = 1; regionIndex <= regions.STAxis.GetNbins(); ++regionIndex) {
    fractionalCorrections[regionIndex] = (((bin_integrals_divided_by_bin_widths_shifted.at(regionIndex))/(bin_integrals_divided_by_bin_widths_nominal.at(regionIndex))) - 1.0);
  }
  return fractionalCorrections;
}

std::map<int, double> generate_eigenmode_fluctuations_map(const int& nEigenmodes, TRandom3 *random_generator, const bool &print_debug=false) {
  if (print_debug) std::cout << "generate_eigenmode_fluctuations_map called with nEigenmodes = " << nEigenmodes << std::endl;
  std::map<int, double> fluctuations_map;
  for (int eigen_index = 0; eigen_index < nEigenmodes; ++eigen_index) {
    fluctuations_map[eigen_index] = random_generator->Gaus(0., 1.);
    if (print_debug) {
      std::cout << "output_map[" << eigen_index << "] = " << fluctuations_map.at(eigen_index) << std::endl;
    }
  }
  return fluctuations_map;
}

std::map<int, double> get_adjustments_from_second_order_pol(const double & slope, const double & quad, const STRegionsStruct & regions, TF1 * unadjusted_tf1, const bool & print_debug=false) {
  // adjustment in bin i would be easy to calculate if we assumed it to be equal to a linear perturbation at the bin center:
  // adjustment(bin i) = 1.0 + slope*(ST_midpoint/ST_norm - 1)
  // but this isn't quite relevant for our use case... we have to consider that the "bin barycenter" doesn't
  // lie at the midpoint given our distribution
  // another way to think about it: suppose our linear correction as a value of ST is as follows:
  // correction(ST) = 1.0 + m*(ST/ST_norm - 1)
  // then, what we really care about is the difference in the prediction of the final nEvents:
  // adjustment(bin i) = (integral(pdf(ST)*correction(ST)) from lo to hi) / (integral(pdf(ST)) from lo to hi)
  //                   = (integral(pdf(ST)*(1+slope*(ST/ST_norm - 1))) from lo to hi) / (integral(pdf(ST)) from lo to hi)
  // slopes are small in our case, so corrections to the normalization need not be accounted for
  // might seem like an overly complex way to do things, and indeed here we probably don't need it
  // but this is more flexible and can be modified easily to implement corrections more complex than a linear correction
  // update: turns out we did need the more complex case after all :-)

  if (print_debug) std::cout << "get_adjustment_from_slope called with slope: " << slope << ", quad: " << quad << std::endl;
  double STNormTargetTmp = 0.5*(regions.STNormRangeMin + regions.STNormRangeMax); // different from options.STNormTarget because this is a different STRegionsStruct
  if (print_debug) std::cout << "STNormTargetTmp = " << STNormTargetTmp << std::endl;
  TF1 *adjusted_tf1 = new TF1((std::string("slope_quad_adjusted_") + unadjusted_tf1->GetName()).c_str(), [&](double *x, double *p){ (void)p; return ((1.0 + slope*(((x[0])/STNormTargetTmp) - 1.0) + quad*(std::pow((x[0])/STNormTargetTmp, 2) - 1.0))*(unadjusted_tf1->Eval(x[0]))); }, regions.STNormRangeMin, ST_MAX_RANGE, 0);
  // in the lambda expression above, the (void)p is just to avoid a compilation error with "gcc -Werror=unused-parameter"
  std::map<int, double> adjustments_from_second_order_pol;
  for (int regionIndex = 1; regionIndex <= regions.STAxis.GetNbins(); ++regionIndex) {
    double bin_low_edge = regions.STAxis.GetBinLowEdge(regionIndex);
    double bin_up_edge = regions.STAxis.GetBinUpEdge(regionIndex);
    double integral_adjusted = adjusted_tf1->Integral(bin_low_edge, bin_up_edge, TF1_INTEGRAL_REL_TOLERANCE);
    double integral_unadjusted = unadjusted_tf1->Integral(bin_low_edge, bin_up_edge, TF1_INTEGRAL_REL_TOLERANCE);
    adjustments_from_second_order_pol[regionIndex] = integral_adjusted/integral_unadjusted;
    if (print_debug) std::cout << "At regionIndex: " << regionIndex << ", bin_low_edge: " << bin_low_edge << ", bin_up_edge: " << bin_up_edge << ", integral_adjusted: " << integral_adjusted << ", integral_unadjusted: " << integral_unadjusted << ", adjustment: " << adjustments_from_second_order_pol.at(regionIndex) << std::endl;
  }
  delete adjusted_tf1;
  return adjustments_from_second_order_pol;
}

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  TH1::AddDirectory(kFALSE);
  RooMsgService::instance().setGlobalKillBelow(MsgLevel::WARNING);
  ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(TF1_INTEGRAL_REL_TOLERANCE);
  do_sanity_checks_customizationTypes();

  tmArgumentParser argumentParser = tmArgumentParser("Run script that prints useful info about the normalization.");
  argumentParser.addArgument("sourceData", "", true, "Comma-separated list of input specifications. An input specification can be either: one single file path, which is used as the source for event data; or two strings separated by the \"!\" symbol, in which case the first string is taken as the source file path as in the preceding case, and the second string is read as a list of five nJets-dependent corrections to the weights (at nJets=2,3,4,5, and >=6 respectively), separated by the \"#\" character; or three strings separated by the \"!\" symbol with two arguments as in the preceding case, and a floating point number used as an extra weight applied to all events from the source.");
  argumentParser.addArgument("outputFolder", "", true, "Output folder.");
  argumentParser.addArgument("selection", "", true, "Name of selection: \"singlemedium\", \"signal_loose\", etc.");
  // argumentParser.addArgument("fetchMCWeights", "false", false, "If this argument is set, then MC weights are read in from the input file.");
  argumentParser.addArgument("getJECShiftedDistributions", "false", false, "If this argument is set, then JEC-shifted distributions are also saved.");
  argumentParser.addArgument("identifier", "", true, "Identifier: \"MC_GJet17\", \"MC_GJet\", etc.");
  argumentParser.addArgument("nJetsNorm", "2", false, "nJets bin to use for normalization.");
  argumentParser.addArgument("STBoundariesSourceFile", "STRegionBoundaries_normOptimization.dat", false, "Source file for reading in ST region boundaries.");
  argumentParser.addArgument("yearString", "all", false, "String with year: can take values \"2016\", \"2017\", \"2018\", or \"all\".");
  argumentParser.addArgument("PDF_nSTBins", "25", false, "Number of bins for plotting datasets.");
  argumentParser.addArgument("rhoNominal", "", true, "Value of the AGK parameter rho to use for the low njets shape.");
  argumentParser.addArgument("preNormalizationBuffer", "200.0", false, "Buffer in ST to use before normalization bin for the kernel.");
  argumentParser.addArgument("adjustmentPlots_min", "0.5", false, "Min of y-range for bottom adjustments plot.");
  argumentParser.addArgument("adjustmentPlots_max", "5.5", false, "Max of y-range for bottom adjustments plot.");
  argumentParser.addArgument("minAllowedEMST", "-1.0", false, "Minimum allowable value of the electromagnetic component of ST. Useful for single photon selections.");
  argumentParser.addArgument("readParametersFromFiles", "/dev/null,/dev/null", false, "If this argument is set, then no fits are performed; instead, the fit parameters is read in from the file locations given as the value of this argument. This should be a list of precisely two files separated by a comma: in order, the binned parameters, and a file containing ST region boundaries to use for saving the (observed/best-fit) ratio adjustments.");
  argumentParser.addArgument("plotConcise", "false", false, "If this argument is set, then only the linear fit and associated errors are plotted.");
  argumentParser.addArgument("disableStrictChecks", "false", false, "If this argument is set, then some strict checks are disabled. Used primarily for plots that don't contribute paramaters to the analysis, to allow some leeway with low stats.");
  argumentParser.setPassedStringValues(argc, argv);
  optionsStruct options = getOptionsFromParser(argumentParser);
  std::cout << "Options passed:" << std::endl << options << std::endl;

  int mkdir_return_status = system(("set -x && mkdir -p " + options.outputFolder + " && set +x").c_str());
  if (mkdir_return_status != 0) {
    std::cout << "ERROR in creating output folder with path: " << options.outputFolder << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // For some strange reason, the executable sometimes crashes at the end of this function
  // with something like:
  // *** Error in `./fitScripts/bin/runFits': free(): invalid pointer: 0x0000000008c3dbb8 ***
  // Workaround: don't check exit status passed to bash, just create an execution status file, set its contents to the number "1"
  // when the script starts, and rewrite contents to "0" at the end... read this status file from any calling script
  std::ofstream executionStatusOutputFileStart((options.outputFolder + "/execution_status_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".txt").c_str());
  assert(executionStatusOutputFileStart.is_open());
  executionStatusOutputFileStart << 1 << std::endl;
  executionStatusOutputFileStart.close();

  TFile* outputFile = TFile::Open((options.outputFolder + "/dataSetStorage_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".root").c_str(), "RECREATE");
  if ((!(outputFile->IsOpen())) or (outputFile->IsZombie())) {
    std::cout << "ERROR: unable to open file with name: " << (options.outputFolder + "/dataSetStorage_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".root") << std::endl;
    std::exit(EXIT_FAILURE);
  }
  RooAbsData::setDefaultStorageType(RooAbsData::Tree);
  RooRealVar rooVar_ST("roo_ST", "roo_ST", (options.STRegions.STNormRangeMin - options.preNormalizationBuffer), ST_MAX_RANGE, "GeV");
  RooRealVar rooVar_weight("roo_weight", "roo_weight", 1., 0., 100000.);

  // RooDataSet STDataSet_2Jets = RooDataSet("STDataSet_2JetsBin", "STDataSet_2JetsBin", RooArgSet(rooVar_ST, rooVar_weight), WeightVar(rooVar_weight));
  // // RooDataSet STDataSet_3Jets = RooDataSet("STDataSet_3JetsBin", "STDataSet_3JetsBin", RooArgSet(rooVar_ST, rooVar_weight), WeightVar(rooVar_weight));
  // // RooDataSet STDataSet_4Jets = RooDataSet("STDataSet_4JetsBin", "STDataSet_4JetsBin", RooArgSet(rooVar_ST, rooVar_weight), WeightVar(rooVar_weight));
  // // RooDataSet STDataSet_5Jets = RooDataSet("STDataSet_5JetsBin", "STDataSet_5JetsBin", RooArgSet(rooVar_ST, rooVar_weight), WeightVar(rooVar_weight));
  // // RooDataSet STDataSet_6Jets = RooDataSet("STDataSet_6JetsBin", "STDataSet_6JetsBin", RooArgSet(rooVar_ST, rooVar_weight), WeightVar(rooVar_weight));
  // std::map<int, RooDataSet*> STDataSets = {
  //   {2, &(STDataSet_2Jets)},
  //   // {3, &(STDataSet_3Jets)},
  //   // {4, &(STDataSet_4Jets)},
  //   // {5, &(STDataSet_5Jets)},
  //   // {6, &(STDataSet_6Jets)}
  // };
  // // idiotic, I know, but the following compiles and results in a segfault that I've been unable to debug:
  // // std::map<int, RooDataSet> STDataSets;
  // // for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
  // //   STDataSets[nJetsBin] = RooDataSet(("STDataSet_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("STDataSet_" + std::to_string(nJetsBin) + "JetsBin").c_str(), RooArgSet(rooVar_ST, rooVar_weight), WeightVar(rooVar_weight));
  // // }
  // All the above is commented out because it isn't flexible enough... the following works just as well
  RooDataSet STDataSet_low_njets = RooDataSet("STDataSet_lowNJetsBin", "STDataSet_lowNJetsBin", RooArgSet(rooVar_ST, rooVar_weight), WeightVar(rooVar_weight));

  std::map<int, TH1D> STHistograms;
  std::map<int, TH1D> STHistograms_JECDown;
  std::map<int, TH1D> STHistograms_JECUp;
  for (int nJetsBin = options.nJetsNorm; nJetsBin <= 6; ++nJetsBin) {
    STHistograms[nJetsBin] = TH1D(("STHistogram_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("ST distribution, " + std::to_string(nJetsBin) + " Jets bin;ST(GeV);weighted events/bin").c_str(), (options.STRegions.STBoundaries.size()-1), &(options.STRegions.STBoundaries.at(0)));
    (STHistograms.at(nJetsBin)).Sumw2();
    STHistograms_JECDown[nJetsBin] = TH1D(("STHistogram_" + std::to_string(nJetsBin) + "JetsBin_JECDown").c_str(), ("ST distribution, " + std::to_string(nJetsBin) + " Jets bin, JEC Down;ST(GeV);weighted events/bin").c_str(), (options.STRegions.STBoundaries.size()-1), &(options.STRegions.STBoundaries.at(0)));
    (STHistograms_JECDown.at(nJetsBin)).Sumw2();
    STHistograms_JECUp[nJetsBin] = TH1D(("STHistogram_" + std::to_string(nJetsBin) + "JetsBin_JECUp").c_str(), ("ST distribution, " + std::to_string(nJetsBin) + " Jets bin, JEC Up;ST(GeV);weighted events/bin").c_str(), (options.STRegions.STBoundaries.size()-1), &(options.STRegions.STBoundaries.at(0)));
    (STHistograms_JECUp.at(nJetsBin)).Sumw2();
  }

  for (int source_data_index = 0; source_data_index < static_cast<int>((options.sourceData).size()); ++source_data_index) {
    // TH1D * pileup_weights = nullptr;
    // TFile * pu_reweighting_source_file = nullptr;
    // if (((options.sourceData).at(source_data_index)).PUReweightingNeeded) {
    //   pu_reweighting_source_file = TFile::Open((((options.sourceData).at(source_data_index)).PUWeightsPath).c_str(), "READ");
    //   assert((pu_reweighting_source_file->IsOpen() && !(pu_reweighting_source_file->IsZombie())));
    //   pu_reweighting_source_file->GetObject("pileupWeights", pileup_weights);
    //   assert(pileup_weights != nullptr);
    // }

    if (((options.sourceData).at(source_data_index)).custom_weight_overall == 0) continue;

    TChain * inputChain = new TChain("ggNtuplizer/EventTree");
    std::cout << "Adding events from file: " << ((options.sourceData).at(source_data_index)).sourceFilePath << std::endl;
    inputChain->Add((((options.sourceData).at(source_data_index)).sourceFilePath).c_str());
    inputChain->SetBranchStatus("*", 0); // so that only the needed branches, explicitly activated below, are read in per event
    float evt_ST = -1.;
    inputChain->SetBranchStatus("b_evtST", 1);
    inputChain->SetBranchAddress("b_evtST", &(evt_ST));
    float evt_ST_JECDown = -1.;
    float evt_ST_JECUp = -1.;
    if (options.getJECShiftedDistributions) {
      inputChain->SetBranchStatus("b_evtST_shifted_JECDown", 1);
      inputChain->SetBranchAddress("b_evtST_shifted_JECDown", &(evt_ST_JECDown));
      inputChain->SetBranchStatus("b_evtST_shifted_JECUp", 1);
      inputChain->SetBranchAddress("b_evtST_shifted_JECUp", &(evt_ST_JECUp));
    }
    float evt_ST_EM = -1.;
    inputChain->SetBranchStatus("b_evtST_electromagnetic", 1);
    inputChain->SetBranchAddress("b_evtST_electromagnetic", &(evt_ST_EM));
    int evt_nJets = -1;
    inputChain->SetBranchStatus("b_nJetsDR", 1);
    inputChain->SetBranchAddress("b_nJetsDR", &(evt_nJets));
    int evt_nJets_JECDown = -1;
    int evt_nJets_JECUp = -1;
    if (options.getJECShiftedDistributions) {
      inputChain->SetBranchStatus("b_nJetsDR_shifted_JECDown", 1);
      inputChain->SetBranchAddress("b_nJetsDR_shifted_JECDown", &(evt_nJets_JECDown));
      inputChain->SetBranchStatus("b_nJetsDR_shifted_JECUp", 1);
      inputChain->SetBranchAddress("b_nJetsDR_shifted_JECUp", &(evt_nJets_JECUp));
    }
    double MCXSecWeight = -1.;
    float MCGenWeight = -1.;
    float MCPrefiringWeight = -1.;
    float MCScaleFactorWeight = -1.;
    double MCPUWeight = -1.;
    // if (options.fetchMCWeights) {
    if (((options.sourceData).at(source_data_index)).fetchMCWeights) {
      inputChain->SetBranchStatus("b_MCXSecWeight", 1);
      inputChain->SetBranchAddress("b_MCXSecWeight", &(MCXSecWeight));
      inputChain->SetBranchStatus("genWeight", 1);
      inputChain->SetBranchAddress("genWeight", &(MCGenWeight));
      inputChain->SetBranchStatus("b_evtPrefiringWeight", 1);
      inputChain->SetBranchAddress("b_evtPrefiringWeight", &(MCPrefiringWeight));
      inputChain->SetBranchStatus("b_evtphotonMCScaleFactor", 1);
      inputChain->SetBranchAddress("b_evtphotonMCScaleFactor", &(MCScaleFactorWeight));
      inputChain->SetBranchStatus("b_PUWeightNoSelection", 1);
      inputChain->SetBranchAddress("b_PUWeightNoSelection", &(MCPUWeight));
    }

    // std::vector<int> * evt_BX_for_PU = nullptr;
    // std::vector<float> * evt_PU = nullptr;
    // if (((options.sourceData).at(source_data_index)).PUReweightingNeeded) {
    //   inputChain->SetBranchStatus("puBX", 1);
    //   inputChain->SetBranchAddress("puBX", &(evt_BX_for_PU));
    //   inputChain->SetBranchStatus("puTrue", 1);
    //   inputChain->SetBranchAddress("puTrue", &(evt_PU));
    // }

    long nEntries = inputChain->GetEntries();

    tmProgressBar *progressBar = new tmProgressBar(nEntries);
    int tmp = static_cast<int>(0.5 + 1.0*nEntries/20);
    int progressBarUpdatePeriod = tmp > 1 ? tmp : 1;
    progressBar->initialize();
    for (Long64_t entryIndex = 0; entryIndex < nEntries; ++entryIndex) {
      Long64_t loadStatus = inputChain->LoadTree(entryIndex);
      assert(loadStatus >= 0);
      int nBytesRead = inputChain->GetEntry(entryIndex, 0); // Get only the required branches
      assert(nBytesRead > 0);
      if ((entryIndex > 0) && (((entryIndex % static_cast<Long64_t>(progressBarUpdatePeriod)) == 0) || (entryIndex == (nEntries-1)))) progressBar->updateBar(static_cast<double>(1.0*entryIndex/nEntries), entryIndex);

      if ((evt_ST < (options.STRegions.STNormRangeMin - options.preNormalizationBuffer)) || (evt_ST > ST_MAX_RANGE)) continue;

      int nJetsBin = (evt_nJets <= 6) ? evt_nJets : 6;
      // if (nJetsBin < options.nJetsNorm) continue;
      int nJetsBin_JECDown = -1;
      int nJetsBin_JECUp = -1;
      if (options.getJECShiftedDistributions) {
        nJetsBin_JECDown = (evt_nJets_JECDown <= 6) ? evt_nJets_JECDown : 6;
        nJetsBin_JECUp = (evt_nJets_JECUp <= 6) ? evt_nJets_JECUp : 6;
      }

      if ((options.minAllowedEMST > 0.) && (evt_ST_EM <= options.minAllowedEMST)) continue;

      double eventWeight = 1.0;
      // if (options.fetchMCWeights) {
      if (((options.sourceData).at(source_data_index)).fetchMCWeights) {
        eventWeight *= (MCXSecWeight*MCGenWeight*MCPrefiringWeight*MCScaleFactorWeight*MCPUWeight);
      }
      // if (((options.sourceData).at(source_data_index)).PUReweightingNeeded) {
      //   float eventPU = -1.;
      //   for (unsigned int BXCounter = 0; BXCounter < static_cast<unsigned int>((*evt_BX_for_PU).size()); ++BXCounter) {
      //     int bx = (*evt_BX_for_PU).at(BXCounter);
      //     if (bx == 0) {
      //       eventPU = (*evt_PU).at(BXCounter);
      //       break;
      //     }
      //   }
      //   assert(eventPU > 0.);
      //   eventWeight *= (pileup_weights->GetBinContent(pileup_weights->GetXaxis()->FindFixBin(eventPU)));
      // }
      if ((((options.sourceData).at(source_data_index)).nJetsReweightingNeeded) && (nJetsBin >= options.nJetsNorm)) {
	eventWeight *= (((options.sourceData).at(source_data_index)).custom_nJets_weights).at(nJetsBin);
      }
      if (((options.sourceData).at(source_data_index)).customWeightingNeeded) {
	eventWeight *= ((options.sourceData).at(source_data_index)).custom_weight_overall;
      }
      double eventWeight_histograms = -1.;
      if (nJetsBin >= options.nJetsNorm) eventWeight_histograms = eventWeight/((STHistograms.at(nJetsBin)).GetXaxis()->GetBinWidth((STHistograms.at(nJetsBin)).FindFixBin(evt_ST)));

      if (nJetsBin == options.nJetsNorm) {
        rooVar_ST.setVal(evt_ST);
        // (STDataSets.at(nJetsBin))->add(RooArgSet(rooVar_ST), eventWeight);
        STDataSet_low_njets.add(RooArgSet(rooVar_ST), eventWeight);
      }

      if ((evt_ST >= options.STRegions.STNormRangeMin) && nJetsBin >= options.nJetsNorm) {
        (STHistograms.at(nJetsBin)).Fill(evt_ST, eventWeight_histograms);
      } // no "pre-norm buffer" needed for histograms

      if (options.getJECShiftedDistributions) {
        if ((evt_ST_JECDown >= options.STRegions.STNormRangeMin) && nJetsBin_JECDown >= options.nJetsNorm) {
          (STHistograms_JECDown.at(nJetsBin_JECDown)).Fill(evt_ST_JECDown, eventWeight_histograms);
        } // no "pre-norm buffer" needed for histograms

        if ((evt_ST_JECUp >= options.STRegions.STNormRangeMin) && nJetsBin_JECUp >= options.nJetsNorm) {
          (STHistograms_JECUp.at(nJetsBin_JECUp)).Fill(evt_ST_JECUp, eventWeight_histograms);
        } // no "pre-norm buffer" needed for histograms
      }
    }
    progressBar->terminate();
    // if (((options.sourceData).at(source_data_index)).PUReweightingNeeded) {
    //   pu_reweighting_source_file->Close();
    // }
  }

  // for (int nJetsBin = options.nJetsNorm; nJetsBin <= 6; ++nJetsBin) {
  //   (STDataSets.at(nJetsBin))->Print();
  // }
  // (STDataSets.at(options.nJetsNorm))->Print();
  STDataSet_low_njets.Print();

  // A few useful initializations
  std::map<std::string, double> fitParametersBinned;
  std::map<std::string, double> fitParameterErrorsBinned;
  std::map<customizationType, std::map<int, double> > fit_pvalues;
  std::map<std::string, std::map<int, double> > ftest_pValues;
  std::vector<std::string> adjustments_slope_sqrt_fit_forOutputFile;
  std::vector<std::string> ratio_adjustments_forOutputFile;
  customizationType customization_type_for_adjustments_output = customizationType::Slope;
  customizationType customization_type_denominator_for_ratios = customizationType::ScaleOnly;
  double scale_minVal = 0.0;
  double scale_maxVal = 15.0;
  double slope_minVal = -1.0/(((ST_MAX_RANGE)/(options.STNormTarget)) - 1.0);
  double slope_maxVal = 5.0;
  double sqrt_minVal = -1.0/(std::sqrt((ST_MAX_RANGE)/(options.STNormTarget)) - 1.0);
  double sqrt_maxVal = 25.0;
  double quad_minVal = -1.0/((std::pow((ST_MAX_RANGE)/(options.STNormTarget), 2)) - 1.0);
  double quad_maxVal = 3.0;
  TRandom3 *random_generator = new TRandom3(1234); // for repeatability
  rooVar_ST.setRange("normRange", options.STRegions.STNormRangeMin, options.STRegions.STNormRangeMax);
  rooVar_ST.setRange("fitRange", options.STRegions.STNormRangeMin, ST_MAX_RANGE);
  rooVar_ST.setRange("plotRange", (options.STRegions.STNormRangeMin - options.preNormalizationBuffer), ST_MAX_RANGE);

  // If fits are to be read in from another file, do it now
  if (options.readParametersFromFiles) {
    std::string lineFromFile;

    std::ifstream inputFileObject_binned(options.inputBinnedParametersFileName.c_str());
    assert(inputFileObject_binned.is_open());
    std::string fit_parameter_name;
    for (int nJetsBin = (1+options.nJetsNorm); nJetsBin <= 6; ++nJetsBin) {
      for (int customization_type_index = customizationTypeFirst; customization_type_index < static_cast<int>(customizationType::nCustomizationTypes); ++customization_type_index) {
        customizationType customization_type = static_cast<customizationType>(customization_type_index);
        for (int parameter_index = 0; parameter_index < customizationTypeNPars.at(customization_type); ++parameter_index) {
          fit_parameter_name = get_parameter_name(customization_type, parameter_index, nJetsBin);
          assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned[fit_parameter_name] = parseLineForFloatWithCheck(lineFromFile, fit_parameter_name);
        }
        for (int eigen_index = 0; eigen_index < customizationTypeNPars.at(customization_type); ++eigen_index) {
          for (int parameter_index = 0; parameter_index < customizationTypeNPars.at(customization_type); ++parameter_index) {
            fit_parameter_name = get_eigencoefficient_name(customization_type, eigen_index, parameter_index, nJetsBin);
            assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned[fit_parameter_name] = parseLineForFloatWithCheck(lineFromFile, fit_parameter_name);
          }
          fit_parameter_name = get_eigenerror_name(customization_type, eigen_index, nJetsBin);
          assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned[fit_parameter_name] = parseLineForFloatWithCheck(lineFromFile, fit_parameter_name);
        }
      }
    }
    assert(!(getline(inputFileObject_binned, lineFromFile))); // makes sure there's nothing else in the input file
    inputFileObject_binned.close();
  }

  // First get the low njets RooKeysPdf
  RooKeysPdf pdf_low_njets = RooKeysPdf("pdf_low_njets", "pdf_low_njets", rooVar_ST, STDataSet_low_njets, RooKeysPdf::MirrorLeft, options.rhoNominal);

  // Plot low njets shape and dataset
  TCanvas binned_pdfCanvas_low_nJets = TCanvas("c_dataSetAndPdf_binned_low_nJets", "c_dataSetAndPdf_binned_low_nJets", 2560, 1440);
  gStyle->SetOptStat(0);
  TLegend legend_dataSetsAndPdf_low_nJets_binned = TLegend(0.5, 0.6, 0.9, 0.9);
  (STHistograms.at(options.nJetsNorm)).SetLineColor(static_cast<EColor>(kBlack)); (STHistograms.at(options.nJetsNorm)).Draw();
  (STHistograms.at(options.nJetsNorm)).GetYaxis()->SetRange(((STHistograms.at(options.nJetsNorm)).GetMaximum())/10000., ((STHistograms.at(options.nJetsNorm)).GetMaximum()));
  binned_pdfCanvas_low_nJets.Update();
  TLegendEntry *legendEntry_dataset_low_nJets = legend_dataSetsAndPdf_low_nJets_binned.AddEntry(&(STHistograms.at(options.nJetsNorm)), ("data, " + std::to_string(options.nJetsNorm) + " Jets").c_str());
  legendEntry_dataset_low_nJets->SetMarkerColor(static_cast<EColor>(kBlack)); legendEntry_dataset_low_nJets->SetLineColor(static_cast<EColor>(kBlack)); legendEntry_dataset_low_nJets->SetTextColor(static_cast<EColor>(kBlack));

  customizedPDF pdf_low_njets_customized(&pdf_low_njets, &rooVar_ST, options.STNormTarget, customizationType::ScaleOnly);
  pdf_low_njets_customized.setNominalScale("fitRange", ((STHistograms.at(options.nJetsNorm)).Integral(1, (STHistograms.at(options.nJetsNorm)).GetXaxis()->GetNbins(), "width")));
  TF1 pdf_low_njets_customized_TF1 = TF1("pdf_low_njets_customized_TF1", pdf_low_njets_customized, options.STRegions.STNormRangeMin, ST_MAX_RANGE, 1);
  pdf_low_njets_customized_TF1.SetParameter(0, 1.0);
  pdf_low_njets_customized_TF1.SetLineColor(static_cast<EColor>(kBlue));
  pdf_low_njets_customized_TF1.SetLineWidth(2);
  pdf_low_njets_customized_TF1.Draw("CSAME");
  TLegendEntry *legendEntry_low_nJets_kernel = legend_dataSetsAndPdf_low_nJets_binned.AddEntry(&(pdf_low_njets_customized_TF1), ("kernel, " + std::to_string(options.nJetsNorm) + " Jets").c_str());
  legendEntry_low_nJets_kernel->SetMarkerColor(static_cast<EColor>(kBlue)); legendEntry_low_nJets_kernel->SetLineColor(static_cast<EColor>(kBlue)); legendEntry_low_nJets_kernel->SetTextColor(static_cast<EColor>(kBlue));
  binned_pdfCanvas_low_nJets.Update();
  gPad->SetLogy();
  binned_pdfCanvas_low_nJets.Update();
  legend_dataSetsAndPdf_low_nJets_binned.Draw(); binned_pdfCanvas_low_nJets.Update();
  binned_pdfCanvas_low_nJets.SaveAs((options.outputFolder + "/binned_pdfAndData_" + std::to_string(options.nJetsNorm) + "JetsBin_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".pdf").c_str());

  for (int nJetsBin = (1+options.nJetsNorm); nJetsBin <= 6; ++nJetsBin) {
    printSeparator();
    std::cout << "Finding fits at nJetsBin = " << nJetsBin << std::endl;

    // some useful initializations
    std::map<customizationType, std::map<int, parameter_initialization_struct> > parameter_initializations = {
      {customizationType::ScaleOnly, {
          {0, parameter_initialization_struct(get_parameter_name(customizationType::ScaleOnly, 0, nJetsBin), 1.0, scale_minVal, scale_maxVal)}
        }
      },
      {customizationType::Slope, {
          {0, parameter_initialization_struct(get_parameter_name(customizationType::Slope, 0, nJetsBin), 1.0, scale_minVal, scale_maxVal)},
          {1, parameter_initialization_struct(get_parameter_name(customizationType::Slope, 1, nJetsBin), 0., slope_minVal, slope_maxVal)}
        }
      },
      {customizationType::Sqrt, {
          {0, parameter_initialization_struct(get_parameter_name(customizationType::Sqrt, 0, nJetsBin), 1.0, scale_minVal, scale_maxVal)},
          {1, parameter_initialization_struct(get_parameter_name(customizationType::Sqrt, 1, nJetsBin), 0., sqrt_minVal, sqrt_maxVal)}
        }
      },
      {customizationType::SlopeSqrt, {
          {0, parameter_initialization_struct(get_parameter_name(customizationType::SlopeSqrt, 0, nJetsBin), 1.0, scale_minVal, scale_maxVal)},
          {1, parameter_initialization_struct(get_parameter_name(customizationType::SlopeSqrt, 1, nJetsBin), 0., slope_minVal, slope_maxVal)},
          {2, parameter_initialization_struct(get_parameter_name(customizationType::SlopeSqrt, 2, nJetsBin), 0., sqrt_minVal, sqrt_maxVal)}
        }
      },
      {customizationType::SlopeSqrtQuad, {
          {0, parameter_initialization_struct(get_parameter_name(customizationType::SlopeSqrtQuad, 0, nJetsBin), 1.0, scale_minVal, scale_maxVal)},
          {1, parameter_initialization_struct(get_parameter_name(customizationType::SlopeSqrtQuad, 1, nJetsBin), 0., slope_minVal, slope_maxVal)},
          {2, parameter_initialization_struct(get_parameter_name(customizationType::SlopeSqrtQuad, 2, nJetsBin), 0., sqrt_minVal, sqrt_maxVal)},
          {3, parameter_initialization_struct(get_parameter_name(customizationType::SlopeSqrtQuad, 3, nJetsBin), 0., quad_minVal, quad_maxVal)}
        }
      }
    };
    assert(static_cast<int>(parameter_initializations.size()) == static_cast<int>(customizationType::nCustomizationTypes));
    for (int customization_type_index = customizationTypeFirst; customization_type_index < static_cast<int>(customizationType::nCustomizationTypes); ++customization_type_index) {
      customizationType customization_type = static_cast<customizationType>(customization_type_index);
      assert(static_cast<int>((parameter_initializations.at(customization_type)).size()) == customizationTypeNPars.at(customization_type));
    }

    std::map<customizationType, customizedPDF*> customized_pdfs;
    std::map<customizationType, customizedTF1*> customized_tf1s;
    for (int customization_type_index = customizationTypeFirst; customization_type_index < static_cast<int>(customizationType::nCustomizationTypes); ++customization_type_index) {
      customizationType customization_type = static_cast<customizationType>(customization_type_index);
      customizedPDF *customized_pdf = new customizedPDF(&pdf_low_njets, &rooVar_ST, options.STNormTarget, customization_type);
      customized_pdf->setNominalScale("normRange", ((STHistograms.at(nJetsBin)).GetBinContent(1))*((STHistograms.at(nJetsBin)).GetBinWidth(1)));
      customizedTF1 *customized_tf1 = new customizedTF1("pdf_customized_" + std::to_string(nJetsBin) + "JetsBin", customized_pdf, options.STRegions.STNormRangeMin, ST_MAX_RANGE, customization_type);
      customized_tf1->initializeParameters(parameter_initializations.at(customization_type));
      if (options.readParametersFromFiles) {
        customized_tf1->setFitResultsFromSource(fitParametersBinned, nJetsBin);
      }
      else {
        customized_tf1->fitToTH1(STHistograms.at(nJetsBin));
        fit_pvalues[customization_type][nJetsBin] = (customized_tf1->fit_result).pvalue;
        for (int parameter_index = 0; parameter_index < customizationTypeNPars.at(customization_type); ++parameter_index) {
          fitParametersBinned[get_parameter_name(customization_type, parameter_index, nJetsBin)] = ((customized_tf1->fit_result).best_fit_values).at(parameter_index);
          fitParameterErrorsBinned[get_parameter_name(customization_type, parameter_index, nJetsBin)] = ((customized_tf1->fit_result).best_fit_errors).at(parameter_index);
        }
        for (int eigen_index = 0; eigen_index < customizationTypeNPars.at(customization_type); ++eigen_index) {
          for (int parameter_index = 0; parameter_index < customizationTypeNPars.at(customization_type); ++parameter_index) {
            fitParametersBinned[get_eigencoefficient_name(customization_type, eigen_index, parameter_index, nJetsBin)] = ((((customized_tf1->fit_result).eigenmodes).at(eigen_index)).eigenvector).at(parameter_index);
          }
          fitParametersBinned[get_eigenerror_name(customization_type, eigen_index, nJetsBin)] = std::sqrt((((customized_tf1->fit_result).eigenmodes).at(eigen_index)).eigenvalue);
        }
      }
      customized_pdfs[customization_type] = customized_pdf;
      customized_tf1s[customization_type] = customized_tf1;
    }

    // data distributions
    TGraphErrors ratioGraph_binned_nJetsDistribution_to_unadjusted = TGraphErrors();
    ratioGraph_binned_nJetsDistribution_to_unadjusted.SetName(("ratioGraph_binned_nJetsDistribution_to_unadjusted_at" + std::to_string(nJetsBin) + "Jets").c_str());
    ratioGraph_binned_nJetsDistribution_to_unadjusted.SetTitle(("ST distribution at " + std::to_string(nJetsBin) + " Jets / unadjusted").c_str());
    TGraphErrors ratioGraph_binned_nJetsDistribution_to_unadjusted_JECDown = TGraphErrors();
    ratioGraph_binned_nJetsDistribution_to_unadjusted_JECDown.SetName(("ratioGraph_binned_nJetsDistribution_to_unadjusted_at" + std::to_string(nJetsBin) + "Jets_JECDown").c_str());
    ratioGraph_binned_nJetsDistribution_to_unadjusted_JECDown.SetTitle(("ST distribution at " + std::to_string(nJetsBin) + " Jets / unadjusted, JEC Down").c_str());
    TGraphErrors ratioGraph_binned_nJetsDistribution_to_unadjusted_JECUp = TGraphErrors();
    ratioGraph_binned_nJetsDistribution_to_unadjusted_JECUp.SetName(("ratioGraph_binned_nJetsDistribution_to_unadjusted_at" + std::to_string(nJetsBin) + "Jets_JECUp").c_str());
    ratioGraph_binned_nJetsDistribution_to_unadjusted_JECUp.SetTitle(("ST distribution at " + std::to_string(nJetsBin) + " Jets / unadjusted, JEC Up").c_str());

    std::map<int, double> bin_integrals_divided_by_bin_widths_from_low_njets_kernel;
    if (!(options.readParametersFromFiles)) {
      (customized_tf1s.at(customization_type_denominator_for_ratios))->set_TF_parameters_to_nominal();
      bin_integrals_divided_by_bin_widths_from_low_njets_kernel = (customized_tf1s.at(customization_type_denominator_for_ratios))->getBinIntegralsDividedByBinWidthFromTF1(options.STRegions);
    }
    for (int binCounter = 1; binCounter <= (STHistograms.at(nJetsBin)).GetXaxis()->GetNbins(); ++binCounter) {
      double STMidpoint = (STHistograms.at(nJetsBin)).GetXaxis()->GetBinCenter(binCounter);
      double binWidth = (STHistograms.at(nJetsBin)).GetXaxis()->GetBinWidth(binCounter);
      double numerator = (STHistograms.at(nJetsBin)).GetBinContent(binCounter);
      double numerator_JECDown = -1.;
      double numerator_JECUp = -1.;
      if (options.getJECShiftedDistributions) {
        numerator_JECDown = (STHistograms_JECDown.at(nJetsBin)).GetBinContent(binCounter);
        numerator_JECUp = (STHistograms_JECUp.at(nJetsBin)).GetBinContent(binCounter);
      }
      (customized_tf1s.at(customization_type_denominator_for_ratios))->set_TF_parameters_to_nominal();
      double denominator = ((customized_tf1s.at(customization_type_denominator_for_ratios))->getTFIntegral((STHistograms.at(nJetsBin)).GetXaxis()->GetBinLowEdge(binCounter), (STHistograms.at(nJetsBin)).GetXaxis()->GetBinUpEdge(binCounter)))/((STHistograms.at(nJetsBin)).GetXaxis()->GetBinWidth(binCounter));
      assert(denominator > 0.);
      double ratio = numerator/denominator;
      double ratio_JECDown = -1.;
      double ratio_JECUp = -1.;
      if (options.getJECShiftedDistributions) {
        ratio_JECDown = numerator_JECDown/denominator;
        ratio_JECUp = numerator_JECUp/denominator;
      }
      double numeratorError = (STHistograms.at(nJetsBin)).GetBinError(binCounter);
      double denominatorError = 0.; // might change later
      double ratioError = ratio*std::sqrt(std::pow(numeratorError/numerator, 2) + std::pow(denominatorError/denominator, 2));
      int graph_currentPointIndex = ratioGraph_binned_nJetsDistribution_to_unadjusted.GetN();
      ratioGraph_binned_nJetsDistribution_to_unadjusted.SetPoint(graph_currentPointIndex, STMidpoint, ratio);
      ratioGraph_binned_nJetsDistribution_to_unadjusted.SetPointError(graph_currentPointIndex, binWidth/(std::sqrt(12)), ratioError);
      if (options.getJECShiftedDistributions) {
        ratioGraph_binned_nJetsDistribution_to_unadjusted_JECDown.SetPoint(graph_currentPointIndex, STMidpoint, ratio_JECDown);
        ratioGraph_binned_nJetsDistribution_to_unadjusted_JECDown.SetPointError(graph_currentPointIndex, binWidth/(std::sqrt(12)), 0.);
        ratioGraph_binned_nJetsDistribution_to_unadjusted_JECUp.SetPoint(graph_currentPointIndex, STMidpoint, ratio_JECUp);
        ratioGraph_binned_nJetsDistribution_to_unadjusted_JECUp.SetPointError(graph_currentPointIndex, binWidth/(std::sqrt(12)), 0.);
      }
    }

    // initialize some variables useful for plots
    double fractionalError_normBin = ((STHistograms.at(nJetsBin)).GetBinError(1))/((STHistograms.at(nJetsBin)).GetBinContent(1));
    // std::cout << "bin error: " << ((STHistograms.at(nJetsBin)).GetBinError(1)) << ", bin content: " << ((STHistograms.at(nJetsBin)).GetBinContent(1)) << ", fractionalError_normBin: " << fractionalError_normBin << std::endl;
    if (options.disableStrictChecks) { // overly strict checks are mainly disabled for plots that don't contribute to analysis parameters
      if (fractionalError_normBin > 1.0) std::cout << "WARNING: fractionalError_normBin = " << fractionalError_normBin << std::endl;
    }
    else {
      assert (fractionalError_normBin <= 1.0); // sanity check, to make sure weights aren't affecting the errors in weird ways...
    }

    // plot the raw shapes
    TCanvas pdfCanvas_binned = TCanvas(("c_dataSetAndPdf_binned_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("c_dataSetAndPdf_binned_" + std::to_string(nJetsBin) + "JetsBin").c_str(), 1600, 1280);
    TLegend legend_dataSetsAndPdf_binned = TLegend(0.5, 0.6, 0.9, 0.9);
    gStyle->SetOptStat(0);

    (STHistograms.at(nJetsBin)).SetLineColor(static_cast<EColor>(kBlack)); (STHistograms.at(nJetsBin)).Draw(); pdfCanvas_binned.Update();
    TLegendEntry *legendEntry_binned_nJetsDistribution = legend_dataSetsAndPdf_binned.AddEntry(&(STHistograms.at(nJetsBin)), (std::to_string(nJetsBin) + " jets distribution").c_str());
    legendEntry_binned_nJetsDistribution->SetMarkerColor(static_cast<EColor>(kBlack)); legendEntry_binned_nJetsDistribution->SetLineColor(static_cast<EColor>(kBlack)); legendEntry_binned_nJetsDistribution->SetTextColor(static_cast<EColor>(kBlack));
    (STHistograms.at(nJetsBin)).GetYaxis()->SetRange(((STHistograms.at(nJetsBin)).GetMaximum())/10000., ((STHistograms.at(nJetsBin)).GetMaximum())); pdfCanvas_binned.Update();
    // draw data with JEC shifts
    if (options.getJECShiftedDistributions) {
      (STHistograms_JECDown.at(nJetsBin)).SetLineColor(static_cast<EColor>(kOrange-3)); (STHistograms_JECDown.at(nJetsBin)).Draw("HIST SAME"); pdfCanvas_binned.Update();
      TLegendEntry *legendEntry_binned_nJetsDistribution_JECShifted = legend_dataSetsAndPdf_binned.AddEntry(&(STHistograms_JECDown.at(nJetsBin)), (std::to_string(nJetsBin) + " jets distribution, JEC up/down").c_str());
      legendEntry_binned_nJetsDistribution_JECShifted->SetMarkerColor(static_cast<EColor>(kOrange-3)); legendEntry_binned_nJetsDistribution_JECShifted->SetLineColor(static_cast<EColor>(kOrange-3)); legendEntry_binned_nJetsDistribution_JECShifted->SetTextColor(static_cast<EColor>(kOrange-3));
      (STHistograms_JECUp.at(nJetsBin)).SetLineColor(static_cast<EColor>(kOrange-3)); (STHistograms_JECUp.at(nJetsBin)).Draw("HIST SAME"); pdfCanvas_binned.Update();
    }

    std::map<customizationType, TGraph> function_graphs;
    std::map<customizationType, TGraph> function_graphs_fluctuationUp;
    std::map<customizationType, TGraph> function_graphs_fluctuationDown;
    std::map<customizationType, std::vector<TGraph> > function_graphs_randomFluctuations;
    std::map<customizationType, std::vector<std::map<int, double> > > generated_random_eigenfluctuations;
    for (int customization_type_index = customizationTypeFirst; customization_type_index < static_cast<int>(customizationType::nCustomizationTypes); ++customization_type_index) {
      customizationType customization_type = static_cast<customizationType>(customization_type_index);
      function_graphs[customization_type] = (customized_tf1s.at(customization_type))->get_nominal_fit_as_TGraph(1000, options.STRegions.STNormRangeMin, ST_MAX_RANGE);
      format_TGraph_as_nominal_fit(function_graphs.at(customization_type), customization_type);
      if (customizationTypeNPars.at(customization_type) >= 1) {
        function_graphs_fluctuationUp[customization_type] = (customized_tf1s.at(customization_type))->get_eigenmode_fluctuation_as_TGraph(0, 1.0, 1000, options.STRegions.STNormRangeMin, ST_MAX_RANGE); // dominant eigenmode only
        format_TGraph_as_fluctuation(function_graphs_fluctuationUp.at(customization_type), customization_type);
        function_graphs_fluctuationDown[customization_type] = (customized_tf1s.at(customization_type))->get_eigenmode_fluctuation_as_TGraph(0, -1.0, 1000, options.STRegions.STNormRangeMin, ST_MAX_RANGE); // dominant eigenmode only
        format_TGraph_as_fluctuation(function_graphs_fluctuationDown.at(customization_type), customization_type);
        generated_random_eigenfluctuations[customization_type] = std::vector<std::map<int, double> >();
        function_graphs_randomFluctuations[customization_type] = std::vector<TGraph>();
        for (int random_fluctuation_counter = 0; random_fluctuation_counter < N_FLUCTUATIONS_TO_PLOT; ++random_fluctuation_counter) {
          std::map<int, double> random_eigenfluctuation = generate_eigenmode_fluctuations_map(customizationTypeNPars.at(customization_type), random_generator);
          (generated_random_eigenfluctuations.at(customization_type)).push_back(random_eigenfluctuation);
          TGraph graph_eigenfluctuation = (customized_tf1s.at(customization_type))->get_eigenmode_fluctuation_as_TGraph(random_eigenfluctuation, 1000, options.STRegions.STNormRangeMin, ST_MAX_RANGE, random_fluctuation_counter);
          format_TGraph_as_random_eigenfluctuation(graph_eigenfluctuation, customization_type);
          (function_graphs_randomFluctuations.at(customization_type)).push_back(graph_eigenfluctuation);
        }
      }
      if ((!(options.plotConcise)) || (options.plotConcise && customizationTypeActiveInConciseWorkflow.at(customization_type))) {
        if ((customizationTypeNPars.at(customization_type) >= 1) && (customizationTypePlotEigenfluctuations.at(customization_type))) {
          for (int random_fluctuation_counter = 0; random_fluctuation_counter < N_FLUCTUATIONS_TO_PLOT; ++random_fluctuation_counter) {
            (((function_graphs_randomFluctuations).at(customization_type)).at(random_fluctuation_counter)).Draw("CSAME"); pdfCanvas_binned.Update();
          }
        }
      }
      if ((customization_type == customization_type_for_adjustments_output) && (!(options.readParametersFromFiles))) {
        (customized_tf1s.at(customization_type))->set_TF_parameters_to_nominal();
        std::map<int, double> bin_integrals_divided_by_bin_widths_nominal = (customized_tf1s.at(customization_type))->getBinIntegralsDividedByBinWidthFromTF1(options.STRegions);
        std::map<int, double> adjustments_nominal = getNominalAdjustmentsFromBinIntegralMaps(options.STRegions, bin_integrals_divided_by_bin_widths_nominal, bin_integrals_divided_by_bin_widths_from_low_njets_kernel);
        std::map<int, std::map<int, double> > adjustmentFractionalCorrections_oneSigmaUp; // first index: eigenmode index
        std::map<int, std::map<int, double> > adjustmentFractionalCorrections_oneSigmaDown; // first index: eigenmode index
        for (int eigen_index = 0; eigen_index < customizationTypeNPars.at(customization_type); ++eigen_index) {
          (customized_tf1s.at(customization_type))->set_TF_parameters_to_eigenmode_fluctuation(eigen_index, 1.0);
          std::map<int, double> bin_integrals_divided_by_bin_widths_withEigenmodeFluctuationsUp = (customized_tf1s.at(customization_type))->getBinIntegralsDividedByBinWidthFromTF1(options.STRegions);
          (customized_tf1s.at(customization_type))->set_TF_parameters_to_eigenmode_fluctuation(eigen_index, -1.0);
          std::map<int, double> bin_integrals_divided_by_bin_widths_withEigenmodeFluctuationsDown = (customized_tf1s.at(customization_type))->getBinIntegralsDividedByBinWidthFromTF1(options.STRegions);
          adjustmentFractionalCorrections_oneSigmaUp[eigen_index] = getAdjustmentFractionalCorrections(options.STRegions, bin_integrals_divided_by_bin_widths_nominal, bin_integrals_divided_by_bin_widths_withEigenmodeFluctuationsUp);
          adjustmentFractionalCorrections_oneSigmaDown[eigen_index] = getAdjustmentFractionalCorrections(options.STRegions, bin_integrals_divided_by_bin_widths_nominal, bin_integrals_divided_by_bin_widths_withEigenmodeFluctuationsDown);
        }
        for (int regionIndex = 2; regionIndex <= options.STRegions.STAxis.GetNbins(); ++regionIndex) {
          adjustments_slope_sqrt_fit_forOutputFile.push_back("float nominalAdjustment_STRegion" + std::to_string(regionIndex) + "_" + std::to_string(nJetsBin) + "Jets=" + std::to_string(adjustments_nominal.at(regionIndex)));
          for (int eigen_index = 0; eigen_index < customizationTypeNPars.at(customization_type); ++eigen_index) {
            adjustments_slope_sqrt_fit_forOutputFile.push_back("float fractionalUncertaintyUp_mode" + std::to_string(eigen_index) + "_STRegion" + std::to_string(regionIndex) + "_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((adjustmentFractionalCorrections_oneSigmaUp.at(eigen_index)).at(regionIndex)));
            adjustments_slope_sqrt_fit_forOutputFile.push_back("float fractionalUncertaintyDown_mode" + std::to_string(eigen_index) + "_STRegion" + std::to_string(regionIndex) + "_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((adjustmentFractionalCorrections_oneSigmaDown.at(eigen_index)).at(regionIndex)));
          }
        }
      }
    }

    // now draw the more important plots
    for (int customization_type_index = customizationTypeFirst; customization_type_index < static_cast<int>(customizationType::nCustomizationTypes); ++customization_type_index) {
      customizationType customization_type = static_cast<customizationType>(customization_type_index);
      if ((!(options.plotConcise)) || (options.plotConcise && customizationTypeActiveInConciseWorkflow.at(customization_type))) {
        function_graphs.at(customization_type).Draw("CSAME"); pdfCanvas_binned.Update();
        TLegendEntry *legendEntry = legend_dataSetsAndPdf_binned.AddEntry(&(function_graphs.at(customization_type)), (customizationTypeLegendLabels.at(customization_type)).c_str());
        set_legend_entry_color(legendEntry, customization_type);
        if ((customizationTypeNPars.at(customization_type) >= 1) && (customizationTypePlotEigenfluctuations.at(customization_type))) {
          (function_graphs_fluctuationUp.at(customization_type)).Draw("CSAME"); pdfCanvas_binned.Update();
          (function_graphs_fluctuationDown.at(customization_type)).Draw("CSAME"); pdfCanvas_binned.Update();
        }
      }
    }

    // draw the data again so the datapoints aren't obscured
    if (options.getJECShiftedDistributions) {
      (STHistograms_JECDown.at(nJetsBin)).Draw("HIST SAME"); (STHistograms_JECUp.at(nJetsBin)).Draw("HIST SAME");
    }
    (STHistograms.at(nJetsBin)).Draw("SAME"); pdfCanvas_binned.Update();
    gPad->SetLogy(); pdfCanvas_binned.Update();
    legend_dataSetsAndPdf_binned.SetFillStyle(0); legend_dataSetsAndPdf_binned.Draw(); pdfCanvas_binned.Update();
    pdfCanvas_binned.SaveAs((options.outputFolder + "/binned_pdfAndData_" + std::to_string(nJetsBin) + "JetsBin_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".pdf").c_str());

    // calculate shape ratios
    std::map<customizationType, TGraph> ratioGraphs_customized_to_unadjusted;
    std::map<customizationType, TGraph> ratioGraphs_customized_to_unadjusted_fluctuationUp;
    std::map<customizationType, TGraph> ratioGraphs_customized_to_unadjusted_fluctuationDown;
    std::map<customizationType, std::vector<TGraph> > ratioGraphs_customized_to_unadjusted_randomFluctuations;
    for (int customization_type_index = customizationTypeFirst; customization_type_index < static_cast<int>(customizationType::nCustomizationTypes); ++customization_type_index) {
      customizationType customization_type = static_cast<customizationType>(customization_type_index);
      (customized_tf1s.at(customization_type))->set_TF_parameters_to_nominal();
      if (customization_type == customization_type_denominator_for_ratios) continue; // this is the denominator wrt which all other adjustments are calculated
      ratioGraphs_customized_to_unadjusted[customization_type] = TGraph();
      (ratioGraphs_customized_to_unadjusted.at(customization_type)).SetName(("ratioGraph_binned_" + customizationTypeNames.at(customization_type) + "_to_" + customizationTypeNames.at(customization_type_denominator_for_ratios) + "_at_" + std::to_string(nJetsBin) + "Jets").c_str());
      ratioGraphs_customized_to_unadjusted_fluctuationUp[customization_type] = TGraph();
      (ratioGraphs_customized_to_unadjusted_fluctuationUp.at(customization_type)).SetName(("ratioGraph_binned_fluctuationUp_" + customizationTypeNames.at(customization_type) + "_to_" + customizationTypeNames.at(customization_type_denominator_for_ratios) + "_at_" + std::to_string(nJetsBin) + "Jets").c_str());
      ratioGraphs_customized_to_unadjusted_fluctuationDown[customization_type] = TGraph();
      (ratioGraphs_customized_to_unadjusted_fluctuationDown.at(customization_type)).SetName(("ratioGraph_binned_fluctuationDown_" + customizationTypeNames.at(customization_type) + "_to_" + customizationTypeNames.at(customization_type_denominator_for_ratios) + "_at_" + std::to_string(nJetsBin) + "Jets").c_str());
      ratioGraphs_customized_to_unadjusted_randomFluctuations[customization_type] = std::vector<TGraph>();
      for (int random_fluctuation_counter = 0; random_fluctuation_counter < N_FLUCTUATIONS_TO_PLOT; ++random_fluctuation_counter) {
        TGraph ratioGraphs_customized_to_unadjusted_randomFluctuation = TGraph();
        ratioGraphs_customized_to_unadjusted_randomFluctuation.SetName(("ratioGraph_binned_randomFluctuation_" + std::to_string(random_fluctuation_counter) + "_" + customizationTypeNames.at(customization_type) + "_to_" + customizationTypeNames.at(customization_type_denominator_for_ratios) + "_at_" + std::to_string(nJetsBin) + "Jets").c_str());
        (ratioGraphs_customized_to_unadjusted_randomFluctuations.at(customization_type)).push_back(ratioGraphs_customized_to_unadjusted_randomFluctuation);
      }
    }

    for (int STCounter = 0; STCounter <= 1000; ++STCounter) {
      double STVal = options.STRegions.STNormRangeMin + (1.0*STCounter/1000)*(ST_MAX_RANGE - options.STRegions.STNormRangeMin);
      double common_denominator = (customized_tf1s.at(customization_type_denominator_for_ratios))->evaluate_TF_at(STVal);

      for (int customization_type_index = customizationTypeFirst; customization_type_index < static_cast<int>(customizationType::nCustomizationTypes); ++customization_type_index) {
        customizationType customization_type = static_cast<customizationType>(customization_type_index);
        if (customization_type == customization_type_denominator_for_ratios) continue; // this is the denominator wrt which all other adjustments are calculated
        (customized_tf1s.at(customization_type))->set_TF_parameters_to_nominal();
        double pdf_nominal = (customized_tf1s.at(customization_type))->evaluate_TF_at(STVal);
        (customized_tf1s.at(customization_type))->set_TF_parameters_to_eigenmode_fluctuation(0, 1.0); // only dominant eigenmode
        double pdf_fluctuation_up = (customized_tf1s.at(customization_type))->evaluate_TF_at(STVal);
        (customized_tf1s.at(customization_type))->set_TF_parameters_to_eigenmode_fluctuation(0, -1.0); // only dominant eigenmode
        double pdf_fluctuation_down = (customized_tf1s.at(customization_type))->evaluate_TF_at(STVal);
        double ratio = pdf_nominal/common_denominator;
        double ratio_fluctuation_up = pdf_fluctuation_up/common_denominator;
        double ratio_fluctuation_down = pdf_fluctuation_down/common_denominator;
        (ratioGraphs_customized_to_unadjusted.at(customization_type)).SetPoint(STCounter, STVal, std::max(0., ratio));
        (ratioGraphs_customized_to_unadjusted_fluctuationUp.at(customization_type)).SetPoint(STCounter, STVal, std::max(0., ratio_fluctuation_up));
        (ratioGraphs_customized_to_unadjusted_fluctuationDown.at(customization_type)).SetPoint(STCounter, STVal, std::max(0., ratio_fluctuation_down));
        for (int random_fluctuation_counter = 0; random_fluctuation_counter < N_FLUCTUATIONS_TO_PLOT; ++random_fluctuation_counter) {
          (customized_tf1s.at(customization_type))->set_TF_parameters_to_eigenmode_fluctuation((generated_random_eigenfluctuations.at(customization_type)).at(random_fluctuation_counter));
          double pdf_fluctuation = (customized_tf1s.at(customization_type))->evaluate_TF_at(STVal);
          double ratio_fluctuation = pdf_fluctuation/common_denominator;
          ((ratioGraphs_customized_to_unadjusted_randomFluctuations.at(customization_type)).at(random_fluctuation_counter)).SetPoint(STCounter, STVal, std::max(0., ratio_fluctuation));
        }
      }
    }

    // plot shape ratios
    TMultiGraph binned_shape_ratios_multigraph = TMultiGraph(("binned_shape_ratios_multigraph_at" + std::to_string(nJetsBin) + "Jets").c_str(), ("Shape ratios (binned), " + std::to_string(nJetsBin) + " Jets bin").c_str());
    TCanvas binned_shape_ratios_canvas = TCanvas(("c_binnedShapeRatios_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("c_binnedShapeRatios_" + std::to_string(nJetsBin) + "JetsBin").c_str(), 1600, 1280);
    TLegend legend_binned_shape_ratios_multigraph = TLegend(0.1, 0.6, 0.5, 0.9);

    for (int customization_type_index = customizationTypeFirst; customization_type_index < static_cast<int>(customizationType::nCustomizationTypes); ++customization_type_index) {
      customizationType customization_type = static_cast<customizationType>(customization_type_index);
      if (customization_type == customization_type_denominator_for_ratios) continue; // this is the denominator wrt which all other adjustments are calculated
      // first add the random eigenfluctuations
      if ((!(options.plotConcise)) || (options.plotConcise && customizationTypeActiveInConciseWorkflow.at(customization_type))) {
        if ((customizationTypeNPars.at(customization_type) >= 1) && (customizationTypePlotEigenfluctuations.at(customization_type))) {
          for (int random_fluctuation_counter = 0; random_fluctuation_counter < N_FLUCTUATIONS_TO_PLOT; ++random_fluctuation_counter) {
            format_ratio_TGraph_as_random_fluctuation_and_add_to_multigraph((ratioGraphs_customized_to_unadjusted_randomFluctuations.at(customization_type)).at(random_fluctuation_counter), binned_shape_ratios_multigraph, customization_type);
          }
        }
      }
    }
    // then add the fits, so they are plotted on top of the random eigenfluctuations
    for (int customization_type_index = customizationTypeFirst; customization_type_index < static_cast<int>(customizationType::nCustomizationTypes); ++customization_type_index) {
      customizationType customization_type = static_cast<customizationType>(customization_type_index);
      if (customization_type == customization_type_denominator_for_ratios) continue; // this is the denominator wrt which all other adjustments are calculated
      if ((!(options.plotConcise)) || (options.plotConcise && customizationTypeActiveInConciseWorkflow.at(customization_type))) {
        format_ratio_TGraph_as_nominal_and_add_to_multigraph(ratioGraphs_customized_to_unadjusted.at(customization_type), binned_shape_ratios_multigraph, legend_binned_shape_ratios_multigraph, customization_type);
        if ((customizationTypeNPars.at(customization_type) >= 1) && (customizationTypePlotEigenfluctuations.at(customization_type))) {
          format_ratio_TGraph_as_fluctuation_and_add_to_multigraph(ratioGraphs_customized_to_unadjusted_fluctuationUp.at(customization_type), binned_shape_ratios_multigraph, customization_type);
          format_ratio_TGraph_as_fluctuation_and_add_to_multigraph(ratioGraphs_customized_to_unadjusted_fluctuationDown.at(customization_type), binned_shape_ratios_multigraph, customization_type);
        }
      }
    }

    if (options.getJECShiftedDistributions) {
      ratioGraph_binned_nJetsDistribution_to_unadjusted_JECDown.SetLineColor(static_cast<EColor>(kOrange-3)); ratioGraph_binned_nJetsDistribution_to_unadjusted_JECDown.SetDrawOption("P"); binned_shape_ratios_multigraph.Add(&ratioGraph_binned_nJetsDistribution_to_unadjusted_JECDown);
      TLegendEntry *legendEntry_binned_nJetsDistribution_to_unadjusted_JECShifted = legend_binned_shape_ratios_multigraph.AddEntry(&ratioGraph_binned_nJetsDistribution_to_unadjusted_JECDown, (std::to_string(nJetsBin) + " jets distribution, JEC up/down / " + std::to_string(options.nJetsNorm) + " jets kernel").c_str());
      legendEntry_binned_nJetsDistribution_to_unadjusted_JECShifted->SetMarkerColor(static_cast<EColor>(kOrange-3)); legendEntry_binned_nJetsDistribution_to_unadjusted_JECShifted->SetLineColor(static_cast<EColor>(kOrange-3)); legendEntry_binned_nJetsDistribution_to_unadjusted_JECShifted->SetTextColor(static_cast<EColor>(kOrange-3));
      ratioGraph_binned_nJetsDistribution_to_unadjusted_JECUp.SetLineColor(static_cast<EColor>(kOrange-3)); ratioGraph_binned_nJetsDistribution_to_unadjusted_JECUp.SetDrawOption("P"); binned_shape_ratios_multigraph.Add(&ratioGraph_binned_nJetsDistribution_to_unadjusted_JECUp);
    }
    ratioGraph_binned_nJetsDistribution_to_unadjusted.SetLineColor(static_cast<EColor>(kBlack)); ratioGraph_binned_nJetsDistribution_to_unadjusted.SetDrawOption("P"); binned_shape_ratios_multigraph.Add(&ratioGraph_binned_nJetsDistribution_to_unadjusted);
    TLegendEntry *legendEntry_binned_nJetsDistribution_to_unadjusted = legend_binned_shape_ratios_multigraph.AddEntry(&ratioGraph_binned_nJetsDistribution_to_unadjusted, (std::to_string(nJetsBin) + " jets distribution / " + std::to_string(options.nJetsNorm) + " jets kernel").c_str());
    legendEntry_binned_nJetsDistribution_to_unadjusted->SetMarkerColor(static_cast<EColor>(kBlack)); legendEntry_binned_nJetsDistribution_to_unadjusted->SetLineColor(static_cast<EColor>(kBlack)); legendEntry_binned_nJetsDistribution_to_unadjusted->SetTextColor(static_cast<EColor>(kBlack));

    binned_shape_ratios_multigraph.Draw("A");
    legend_binned_shape_ratios_multigraph.SetFillStyle(0);
    legend_binned_shape_ratios_multigraph.Draw();
    binned_shape_ratios_multigraph.GetXaxis()->SetTitle("ST (GeV)");
    binned_shape_ratios_multigraph.GetXaxis()->SetRangeUser((STHistograms.at(nJetsBin)).GetXaxis()->GetBinLowEdge(1), (STHistograms.at(nJetsBin)).GetXaxis()->GetBinUpEdge((STHistograms.at(nJetsBin)).GetXaxis()->GetNbins()));
    binned_shape_ratios_multigraph.GetYaxis()->SetTitle("ratio");
    binned_shape_ratios_multigraph.SetMinimum(options.adjustmentPlots_min);
    binned_shape_ratios_multigraph.SetMaximum(options.adjustmentPlots_max);
    binned_shape_ratios_canvas.Update();
    binned_shape_ratios_canvas.SaveAs((options.outputFolder + "/binned_shapeRatios_" + std::to_string(nJetsBin) + "JetsBin_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".pdf").c_str());

    // if doing a comparison (e.g. reading parameters from MC and applying to data), calculate and save the ratios of the data wrt the nominal adjustment and fits to a straight line
    if (options.readParametersFromFiles) {
      TGraphErrors ratios_wrt_chosen_adjustment = TGraphErrors();
      ratios_wrt_chosen_adjustment.SetName(("ratios_wrt_chosen_adjustment_" + std::to_string(nJetsBin) + "JetsBin").c_str());
      ratios_wrt_chosen_adjustment.SetTitle(";ST [GeV];ratio");
      TGraphErrors ratios_wrt_chosen_adjustment_eigenfluctuation_errors = TGraphErrors();
      ratios_wrt_chosen_adjustment_eigenfluctuation_errors.SetName(("ratios_wrt_chosen_adjustment_eigenfluctuation_errors_"+ std::to_string(nJetsBin) + "JetsBin").c_str());
      ratios_wrt_chosen_adjustment_eigenfluctuation_errors.SetTitle(";ST [GeV];ratio");
      ratios_wrt_chosen_adjustment_eigenfluctuation_errors.SetFillColorAlpha(kGreen+3, 0.5);
      ratios_wrt_chosen_adjustment_eigenfluctuation_errors.SetFillStyle(1001);
      for (int binCounter = 1; binCounter <= (STHistograms.at(nJetsBin)).GetXaxis()->GetNbins(); ++binCounter) {
        double STMidpoint = (STHistograms.at(nJetsBin)).GetXaxis()->GetBinCenter(binCounter);
        double binWidth = (STHistograms.at(nJetsBin)).GetXaxis()->GetBinWidth(binCounter);
        double numerator = (STHistograms.at(nJetsBin)).GetBinContent(binCounter);
        (customized_tf1s.at(customization_type_for_adjustments_output))->set_TF_parameters_to_nominal();
        double denominator = ((customized_tf1s.at(customization_type_for_adjustments_output))->getTFIntegral((STHistograms.at(nJetsBin)).GetXaxis()->GetBinLowEdge(binCounter), (STHistograms.at(nJetsBin)).GetXaxis()->GetBinUpEdge(binCounter)))/((STHistograms.at(nJetsBin)).GetXaxis()->GetBinWidth(binCounter));
	double sum_squares_fractionalErrors = 0.;
	for (int eigen_index = 0; eigen_index < customizationTypeNPars.at(customization_type_for_adjustments_output); ++eigen_index) {
	  (customized_tf1s.at(customization_type_for_adjustments_output))->set_TF_parameters_to_eigenmode_fluctuation(eigen_index, 1.0);
	  double denominator_with_eigenfluctuation_up   = ((customized_tf1s.at(customization_type_for_adjustments_output))->getTFIntegral((STHistograms.at(nJetsBin)).GetXaxis()->GetBinLowEdge(binCounter), (STHistograms.at(nJetsBin)).GetXaxis()->GetBinUpEdge(binCounter)))/((STHistograms.at(nJetsBin)).GetXaxis()->GetBinWidth(binCounter));
	  double denominator_fractional_error_eigenfluctuation_up = std::fabs(denominator_with_eigenfluctuation_up/denominator - 1.0);
	  (customized_tf1s.at(customization_type_for_adjustments_output))->set_TF_parameters_to_eigenmode_fluctuation(eigen_index, -1.0);
	  double denominator_with_eigenfluctuation_down = ((customized_tf1s.at(customization_type_for_adjustments_output))->getTFIntegral((STHistograms.at(nJetsBin)).GetXaxis()->GetBinLowEdge(binCounter), (STHistograms.at(nJetsBin)).GetXaxis()->GetBinUpEdge(binCounter)))/((STHistograms.at(nJetsBin)).GetXaxis()->GetBinWidth(binCounter));
	  double denominator_fractional_error_eigenfluctuation_down = std::fabs(denominator_with_eigenfluctuation_down/denominator - 1.0);
	  double avg_fractional_error = 0.5*(denominator_fractional_error_eigenfluctuation_up + denominator_fractional_error_eigenfluctuation_down);
	  // std::cout << "For eigenindex = " << eigen_index << ", Fractional error: " << avg_fractional_error << std::endl;
	  sum_squares_fractionalErrors += avg_fractional_error*avg_fractional_error;
	}
	double net_fractional_error_denominator = std::sqrt(sum_squares_fractionalErrors);
	// std::cout << "Net fractional error in denominator: " << net_fractional_error_denominator << std::endl;
	(customized_tf1s.at(customization_type_for_adjustments_output))->set_TF_parameters_to_nominal();
        assert(denominator > 0.);
        double ratio = numerator/denominator;
        double numeratorError = (STHistograms.at(nJetsBin)).GetBinError(binCounter);
        double denominatorError = 0.; // might change later
        double ratioError = ratio*std::sqrt(std::pow(numeratorError/numerator, 2) + std::pow(denominatorError/denominator, 2));
        int graph_currentPointIndex = ratios_wrt_chosen_adjustment.GetN();
        ratios_wrt_chosen_adjustment.SetPoint(graph_currentPointIndex, STMidpoint, ratio);
        ratios_wrt_chosen_adjustment.SetPointError(graph_currentPointIndex, binWidth/(std::sqrt(12)), ratioError);
	int graph_errors_currentPointIndex = ratios_wrt_chosen_adjustment_eigenfluctuation_errors.GetN();
	ratios_wrt_chosen_adjustment_eigenfluctuation_errors.SetPoint(graph_errors_currentPointIndex, STMidpoint, 1.0);
        ratios_wrt_chosen_adjustment_eigenfluctuation_errors.SetPointError(graph_errors_currentPointIndex, 0.5*binWidth, net_fractional_error_denominator);
      }
      std::string functional_form_for_TF1 = "";
      if (constants::fit_type_ratios_wrt_chosen_adjustment == fitType_ratios_wrt_chosen_adjustment::Linear) functional_form_for_TF1 = "[0] + [1]*((x/" + std::to_string(options.STNormTarget) + ") - 1.0)";
      else if (constants::fit_type_ratios_wrt_chosen_adjustment == fitType_ratios_wrt_chosen_adjustment::Quad) functional_form_for_TF1 = "[0] + [1]*((x/" + std::to_string(options.STNormTarget) + ") - 1.0) + [2]*(((x/" + std::to_string(options.STNormTarget) + ")**2)-1.0)";
      TF1 fitFunction_ratios_wrt_chosen_adjustment(("ratios_wrt_chosen_adjustment_at" + std::to_string(nJetsBin) + "Jets").c_str(), functional_form_for_TF1.c_str(), options.STRegions.STNormRangeMin, ST_MAX_RANGE);
      fitFunction_ratios_wrt_chosen_adjustment.SetParName(0, ("ratios_wrt_chosen_adjustment_const_" + std::to_string(nJetsBin) + "JetsBin").c_str());
      fitFunction_ratios_wrt_chosen_adjustment.SetParameter(0, 1.0);
      fitFunction_ratios_wrt_chosen_adjustment.SetParLimits(0, -9.5, 10.5);
      fitFunction_ratios_wrt_chosen_adjustment.SetParName(1, ("ratios_wrt_chosen_adjustment_slope_" + std::to_string(nJetsBin) + "JetsBin").c_str());
      fitFunction_ratios_wrt_chosen_adjustment.SetParameter(1, 0.);
      fitFunction_ratios_wrt_chosen_adjustment.SetParLimits(1, -10.0, 10.0);
      if (constants::fit_type_ratios_wrt_chosen_adjustment == fitType_ratios_wrt_chosen_adjustment::Quad) {
	fitFunction_ratios_wrt_chosen_adjustment.SetParName(2, ("ratios_wrt_chosen_adjustment_quad_" + std::to_string(nJetsBin) + "JetsBin").c_str());
	fitFunction_ratios_wrt_chosen_adjustment.SetParameter(2, 0.);
	fitFunction_ratios_wrt_chosen_adjustment.SetParLimits(2, -5.0, 5.0);
      }
      TFitResultPtr ratios_wrt_chosen_adjustment_fit_result;
      ratios_wrt_chosen_adjustment_fit_result = ratios_wrt_chosen_adjustment.Fit(&fitFunction_ratios_wrt_chosen_adjustment, (constants::binnedFitOptions_ratios_wrt_chosen_adjustment).c_str());
      if (ratios_wrt_chosen_adjustment_fit_result->Status() != 0) {
	ratios_wrt_chosen_adjustment_fit_result = ratios_wrt_chosen_adjustment.Fit(&fitFunction_ratios_wrt_chosen_adjustment, (constants::binnedFitOptions_ratios_wrt_chosen_adjustment_backup).c_str());
	assert(ratios_wrt_chosen_adjustment_fit_result->Status() == 0);
      }
      double best_fit_const = ratios_wrt_chosen_adjustment_fit_result->Value(0);
      double best_fit_const_error = ratios_wrt_chosen_adjustment_fit_result->ParError(0);
      double best_fit_slope = ratios_wrt_chosen_adjustment_fit_result->Value(1);
      double best_fit_slope_error = ratios_wrt_chosen_adjustment_fit_result->ParError(1);
      double best_fit_quad = 0.;
      double best_fit_quad_error = 0.;
      if (constants::fit_type_ratios_wrt_chosen_adjustment == fitType_ratios_wrt_chosen_adjustment::Quad) {
	best_fit_quad = ratios_wrt_chosen_adjustment_fit_result->Value(2);
	best_fit_quad_error = ratios_wrt_chosen_adjustment_fit_result->ParError(2);
      }
      TCanvas canvas_ratios_wrt_chosen_adjustment = TCanvas(("c_ratios_wrt_chosen_adjustment_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("c_ratios_wrt_chosen_adjustment_" + std::to_string(nJetsBin) + "JetsBin").c_str(), 1600, 1080);
      gStyle->SetOptStat(0);
      TLegend legend_ratios_wrt_chosen_adjustment = TLegend(0.1, 0.7, 0.9, 0.9);
      legend_ratios_wrt_chosen_adjustment.SetFillStyle(0);
      ratios_wrt_chosen_adjustment_eigenfluctuation_errors.Draw("A20"); canvas_ratios_wrt_chosen_adjustment.Update();
      ratios_wrt_chosen_adjustment_eigenfluctuation_errors.GetXaxis()->SetRangeUser((STHistograms.at(nJetsBin)).GetXaxis()->GetBinLowEdge(1), (STHistograms.at(nJetsBin)).GetXaxis()->GetBinUpEdge((STHistograms.at(nJetsBin)).GetXaxis()->GetNbins())); canvas_ratios_wrt_chosen_adjustment.Update();
      ratios_wrt_chosen_adjustment_eigenfluctuation_errors.GetYaxis()->SetRangeUser(-0.5, 3.5); canvas_ratios_wrt_chosen_adjustment.Update();
      ratios_wrt_chosen_adjustment.Draw("P0"); canvas_ratios_wrt_chosen_adjustment.Update();
      ratios_wrt_chosen_adjustment.SetLineColor(static_cast<EColor>(kBlack)); ratios_wrt_chosen_adjustment.SetLineWidth(2);
      TLegendEntry *legendEntry_graph_ratios_wrt_chosen_adjustment = legend_ratios_wrt_chosen_adjustment.AddEntry(&ratios_wrt_chosen_adjustment, (std::to_string(nJetsBin) + " jets distribution / " + customizationTypeLegendLabels.at(customization_type_for_adjustments_output)).c_str());
      legendEntry_graph_ratios_wrt_chosen_adjustment->SetMarkerColor(static_cast<EColor>(kBlack)); legendEntry_graph_ratios_wrt_chosen_adjustment->SetLineColor(static_cast<EColor>(kBlack)); legendEntry_graph_ratios_wrt_chosen_adjustment->SetTextColor(static_cast<EColor>(kBlack));
      fitFunction_ratios_wrt_chosen_adjustment.SetLineColor(static_cast<EColor>(kBlue)); fitFunction_ratios_wrt_chosen_adjustment.SetLineWidth(2);
      fitFunction_ratios_wrt_chosen_adjustment.Draw("C SAME"); canvas_ratios_wrt_chosen_adjustment.Update();

      TF1 *fitFunction_ratios_wrt_chosen_adjustment_up = (TF1*)(fitFunction_ratios_wrt_chosen_adjustment.Clone());
      fitFunction_ratios_wrt_chosen_adjustment_up->SetName((std::string(fitFunction_ratios_wrt_chosen_adjustment.GetName()) + std::string("_up")).c_str());
      fitFunction_ratios_wrt_chosen_adjustment_up->SetLineColor(static_cast<EColor>(kBlue)); fitFunction_ratios_wrt_chosen_adjustment_up->SetLineWidth(2); fitFunction_ratios_wrt_chosen_adjustment_up->SetLineStyle(kDashed);
      fitFunction_ratios_wrt_chosen_adjustment_up->SetParameter(0, best_fit_const + best_fit_const_error);
      fitFunction_ratios_wrt_chosen_adjustment_up->SetParameter(1, best_fit_slope + best_fit_slope_error);
      if (constants::fit_type_ratios_wrt_chosen_adjustment == fitType_ratios_wrt_chosen_adjustment::Quad) fitFunction_ratios_wrt_chosen_adjustment_up->SetParameter(2, best_fit_quad + best_fit_quad_error);
      fitFunction_ratios_wrt_chosen_adjustment_up->Draw("C SAME"); canvas_ratios_wrt_chosen_adjustment.Update();

      TF1 *fitFunction_ratios_wrt_chosen_adjustment_down = (TF1*)(fitFunction_ratios_wrt_chosen_adjustment.Clone());
      fitFunction_ratios_wrt_chosen_adjustment_down->SetName((std::string(fitFunction_ratios_wrt_chosen_adjustment.GetName()) + std::string("_down")).c_str());
      fitFunction_ratios_wrt_chosen_adjustment_down->SetLineColor(static_cast<EColor>(kBlue)); fitFunction_ratios_wrt_chosen_adjustment_down->SetLineWidth(2); fitFunction_ratios_wrt_chosen_adjustment_down->SetLineStyle(kDashed);
      fitFunction_ratios_wrt_chosen_adjustment_down->SetParameter(0, best_fit_const - best_fit_const_error);
      fitFunction_ratios_wrt_chosen_adjustment_down->SetParameter(1, best_fit_slope - best_fit_slope_error);
      if (constants::fit_type_ratios_wrt_chosen_adjustment == fitType_ratios_wrt_chosen_adjustment::Quad) fitFunction_ratios_wrt_chosen_adjustment_down->SetParameter(2, best_fit_quad - best_fit_quad_error);
      fitFunction_ratios_wrt_chosen_adjustment_down->Draw("C SAME"); canvas_ratios_wrt_chosen_adjustment.Update();

      std::string text_fit_description = ("nominal fit: (" + get_string_precision_n(3, best_fit_const) + " #pm " + get_string_precision_n(3, best_fit_const_error) + ") + (" + get_string_precision_n(3, best_fit_slope) + " #pm " + get_string_precision_n(3, best_fit_slope_error) + ")*(ST/" + get_string_precision_n(5, options.STNormTarget) + " - 1.0)");
      if (constants::fit_type_ratios_wrt_chosen_adjustment == fitType_ratios_wrt_chosen_adjustment::Quad) {
	text_fit_description += (" + (" + get_string_precision_n(3, best_fit_quad) + " #pm " + get_string_precision_n(3, best_fit_quad_error) + ")*((ST/" + get_string_precision_n(5, options.STNormTarget) + ")^2 - 1.0)");
      }
      TLegendEntry *legendEntry_nominal_fit = legend_ratios_wrt_chosen_adjustment.AddEntry(&ratios_wrt_chosen_adjustment, text_fit_description.c_str());
      legendEntry_nominal_fit->SetMarkerColor(static_cast<EColor>(kBlue)); legendEntry_nominal_fit->SetLineColor(static_cast<EColor>(kBlue)); legendEntry_nominal_fit->SetTextColor(static_cast<EColor>(kBlue));
      legend_ratios_wrt_chosen_adjustment.Draw();
      canvas_ratios_wrt_chosen_adjustment.Update();
      canvas_ratios_wrt_chosen_adjustment.SaveAs((options.outputFolder + "/binned_ratios_wrt_chosen_adjustment_" + std::to_string(nJetsBin) + "JetsBin_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".pdf").c_str());
      std::map<int, double> adjustments_from_second_order_pol = get_adjustments_from_second_order_pol(best_fit_slope, best_fit_quad, options.STRegions_for_ratio_wrt_chosen_adjustment, (customized_tf1s.at(customization_type_for_adjustments_output))->raw_TF1);
      for (int regionIndex = 1; regionIndex <= options.STRegions_for_ratio_wrt_chosen_adjustment.STAxis.GetNbins(); ++regionIndex) {
        ratio_adjustments_forOutputFile.push_back("float ratio_adjustment_STRegion" + std::to_string(regionIndex) + "_" + std::to_string(nJetsBin) + "Jets=" + std::to_string(adjustments_from_second_order_pol.at(regionIndex)));
      }
    }

    if (!(options.readParametersFromFiles)) {
      std::cout << "Getting p-values using chi2 values from binned fits..." << std::endl;
      // sanity checks
      assert((((customized_tf1s.at(customizationType::ScaleOnly))->fit_result).ndf) == (1 + (((customized_tf1s.at(customizationType::Slope))->fit_result).ndf)));
      assert((((customized_tf1s.at(customizationType::ScaleOnly))->fit_result).ndf) == (1 + (((customized_tf1s.at(customizationType::Sqrt))->fit_result).ndf)));
      assert((((customized_tf1s.at(customizationType::ScaleOnly))->fit_result).ndf) == (2 + (((customized_tf1s.at(customizationType::SlopeSqrt))->fit_result).ndf)));
      assert((((customized_tf1s.at(customizationType::ScaleOnly))->fit_result).ndf) == (3 + (((customized_tf1s.at(customizationType::SlopeSqrtQuad))->fit_result).ndf)));

      ftest_pValues["unadjusted_vs_slope"][nJetsBin] = get_fTest_prob(((customized_tf1s.at(customizationType::ScaleOnly))->fit_result).chi_sq, ((customized_tf1s.at(customizationType::Slope))->fit_result).chi_sq, ((customized_tf1s.at(customizationType::ScaleOnly))->fit_result).ndf, ((customized_tf1s.at(customizationType::Slope))->fit_result).ndf, options.yearString + "_" + options.identifier + "_" + options.selection);
      ftest_pValues["slope_vs_slope_sqrt"][nJetsBin] = get_fTest_prob(((customized_tf1s.at(customizationType::Slope))->fit_result).chi_sq, ((customized_tf1s.at(customizationType::SlopeSqrt))->fit_result).chi_sq, ((customized_tf1s.at(customizationType::Slope))->fit_result).ndf, ((customized_tf1s.at(customizationType::SlopeSqrt))->fit_result).ndf, options.yearString + "_" + options.identifier + "_" + options.selection);
      ftest_pValues["unadjusted_vs_sqrt"][nJetsBin] = get_fTest_prob(((customized_tf1s.at(customizationType::ScaleOnly))->fit_result).chi_sq, ((customized_tf1s.at(customizationType::Sqrt))->fit_result).chi_sq, ((customized_tf1s.at(customizationType::ScaleOnly))->fit_result).ndf, ((customized_tf1s.at(customizationType::Sqrt))->fit_result).ndf, options.yearString + "_" + options.identifier + "_" + options.selection);
      ftest_pValues["sqrt_vs_slope_sqrt"][nJetsBin] = get_fTest_prob(((customized_tf1s.at(customizationType::Sqrt))->fit_result).chi_sq, ((customized_tf1s.at(customizationType::SlopeSqrt))->fit_result).chi_sq, ((customized_tf1s.at(customizationType::Sqrt))->fit_result).ndf, ((customized_tf1s.at(customizationType::SlopeSqrt))->fit_result).ndf, options.yearString + "_" + options.identifier + "_" + options.selection);
      ftest_pValues["slope_sqrt_vs_slope_sqrt_quad"][nJetsBin] = get_fTest_prob(((customized_tf1s.at(customizationType::SlopeSqrt))->fit_result).chi_sq, ((customized_tf1s.at(customizationType::SlopeSqrtQuad))->fit_result).chi_sq, ((customized_tf1s.at(customizationType::SlopeSqrt))->fit_result).ndf, ((customized_tf1s.at(customizationType::SlopeSqrtQuad))->fit_result).ndf, options.yearString + "_" + options.identifier + "_" + options.selection);
    }

    std::cout << "end1 for nJetsBin: " << nJetsBin << std::endl;
    for (int customization_type_index = customizationTypeFirst; customization_type_index < static_cast<int>(customizationType::nCustomizationTypes); ++customization_type_index) {
      customizationType customization_type = static_cast<customizationType>(customization_type_index);
      delete customized_tf1s.at(customization_type);
      delete customized_pdfs.at(customization_type);
    }
    std::cout << "end2 for nJetsBin: " << nJetsBin << std::endl;

    printSeparator();
  } // ends loop over nJetsBin
  std::cout << "HERE1" << std::endl;
  delete random_generator;

  std::cout << "HERE2" << std::endl;
  if (options.readParametersFromFiles) {
    // write out adjustment values
    std::cout << "HERE3" << std::endl;
    std::ofstream ratio_adjustment_outputFile((options.outputFolder + "/ratio_adjustment_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".dat").c_str());
    assert(ratio_adjustment_outputFile.is_open());
    for (int ratio_adjustment_index = 0; ratio_adjustment_index < static_cast<int>(ratio_adjustments_forOutputFile.size()); ++ratio_adjustment_index) {
      ratio_adjustment_outputFile << ratio_adjustments_forOutputFile.at(ratio_adjustment_index) << std::endl;
    }
    ratio_adjustment_outputFile.close();
    std::cout << "HERE4" << std::endl;
    std::cout << "Ratio adjustments written to file: " << (options.outputFolder + "/ratio_adjustment_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".dat") << std::endl;
    std::cout << "HERE5" << std::endl;
  }
  else {
    // print ftest pvalues from binned chi2 fits in a LaTeX-formatted table
    std::cout << std::endl;
    std::cout << "p-values for binned fit comparisons using f-statistic:" << std::endl;
    std::cout << "\\begin{tabular}{|p{0.14\\textwidth}|p{0.14\\textwidth}|p{0.14\\textwidth}|p{0.14\\textwidth}|p{0.14\\textwidth}|p{0.14\\textwidth}|}" << std::endl;
    std::cout << "  \\hline" << std::endl;
    std::cout << "  p-values \\newline (f-statistic) & unadjusted \\newline vs \\newline linear & linear \\newline vs \\newline (linear+sqrt) & unadjusted \\newline vs \\newline sqrt & sqrt \\newline vs \\newline (linear+sqrt) & (linear+sqrt) \\newline vs \\newline (linear+sqrt \\newline +quad) \\\\ \\hline" << std::endl;
    if (options.nJetsNorm < 3) std::cout << std::fixed << std::setprecision(3) << "  nJets = 3 & " << (ftest_pValues.at("unadjusted_vs_slope")).at(3) << " & " << (ftest_pValues.at("slope_vs_slope_sqrt")).at(3) << " & " << (ftest_pValues.at("unadjusted_vs_sqrt")).at(3) << " & " << (ftest_pValues.at("sqrt_vs_slope_sqrt")).at(3) << " & " << (ftest_pValues.at("slope_sqrt_vs_slope_sqrt_quad")).at(3) << " \\\\ \\hline" << std::endl;
    std::cout << std::fixed << std::setprecision(3) << "  nJets = 4 & " << (ftest_pValues.at("unadjusted_vs_slope")).at(4) << " & " << (ftest_pValues.at("slope_vs_slope_sqrt")).at(4) << " & " << (ftest_pValues.at("unadjusted_vs_sqrt")).at(4) << " & " << (ftest_pValues.at("sqrt_vs_slope_sqrt")).at(4) << " & " << (ftest_pValues.at("slope_sqrt_vs_slope_sqrt_quad")).at(4) << " \\\\ \\hline" << std::endl;
    std::cout << std::fixed << std::setprecision(3) << "  nJets = 5 & " << (ftest_pValues.at("unadjusted_vs_slope")).at(5) << " & " << (ftest_pValues.at("slope_vs_slope_sqrt")).at(5) << " & " << (ftest_pValues.at("unadjusted_vs_sqrt")).at(5) << " & " << (ftest_pValues.at("sqrt_vs_slope_sqrt")).at(5) << " & " << (ftest_pValues.at("slope_sqrt_vs_slope_sqrt_quad")).at(5) << " \\\\ \\hline" << std::endl;
    std::cout << std::fixed << std::setprecision(3) << "  nJets $\\geq$ 6 & " << (ftest_pValues.at("unadjusted_vs_slope")).at(6) << " & " << (ftest_pValues.at("slope_vs_slope_sqrt")).at(6) << " & " << (ftest_pValues.at("unadjusted_vs_sqrt")).at(6) << " & " << (ftest_pValues.at("sqrt_vs_slope_sqrt")).at(6) << " & " << (ftest_pValues.at("slope_sqrt_vs_slope_sqrt_quad")).at(6) << " \\\\ \\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << std::endl;

    // print p-values for best fits in a LaTeX-formatted table
    std::cout << std::endl;
    std::cout << "p-values for binned fits:" << std::endl;
    // create tabular environment
    int n_columns = 1+static_cast<int>(customizationType::nCustomizationTypes); // leftmost column for labels + one for each fit function
    double total_textwidth_fraction_for_pvalue_columns = 0.7;
    double total_textwidth_fraction_per_pvalue_column = total_textwidth_fraction_for_pvalue_columns/n_columns;
    std::cout << "\\begin{tabular}{";
    for (int column_counter = 0; column_counter < n_columns; ++column_counter) std::cout << std::setprecision(3) << "|p{" << total_textwidth_fraction_per_pvalue_column << "\\textwidth}" << std::fixed;
    std::cout << "|}" << std::endl;
    std::cout << "  \\hline" << std::endl;
    // print column header
    std::cout << "  p-values \\newline (fits)";
    for (int customization_type_index = customizationTypeFirst; customization_type_index < static_cast<int>(customizationType::nCustomizationTypes); ++customization_type_index) {
      customizationType customization_type = static_cast<customizationType>(customization_type_index);
      std::cout << " & " << customizationTypeHumanReadableNames.at(customization_type);
    }
    std::cout << " \\\\ \\hline" << std::endl;
    // print p-values
    for (int nJetsBin = (1+options.nJetsNorm); nJetsBin <= 6; ++nJetsBin) {
      if (nJetsBin == 6) std::cout << "  nJets $\\geq$ 6";
      else std::cout << "  nJets = " << nJetsBin;
      for (int customization_type_index = customizationTypeFirst; customization_type_index < static_cast<int>(customizationType::nCustomizationTypes); ++customization_type_index) {
        customizationType customization_type = static_cast<customizationType>(customization_type_index);
        std::cout << std::setprecision(3) << " & " << (fit_pvalues.at(customization_type)).at(nJetsBin) << std::fixed;
      }
      std::cout << " \\\\ \\hline" << std::endl;
    }
    // end tabular environment
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << std::endl;

    // // print best fit values for sqrt fit in a LaTeX-formatted table
    // std::cout << std::endl;
    // std::cout << "Best fit values for sqrt fit:" << std::endl;
    // std::cout << "\\begin{tabular}{|p{0.2\\textwidth}|p{0.2\\textwidth}|p{0.2\\textwidth}|}" << std::endl;
    // std::cout << "  \\hline" << std::endl;
    // std::cout << "  Best fit values & $A$ & $p$ \\\\ \\hline" << std::endl;
    // for (int nJetsBin = (1+options.nJetsNorm); nJetsBin <= 6; ++nJetsBin) {
    //   if (nJetsBin == 6) std::cout << "  nJets $\\geq$ 6";
    //   else std::cout << "  nJets = " << nJetsBin;
    //   std::cout << std::setprecision(3) << " & " << fitParametersBinned.at(get_parameter_name(customizationType::Sqrt, 0, nJetsBin)) << " $\\pm$ " << fitParameterErrorsBinned.at(get_parameter_name(customizationType::Sqrt, 0, nJetsBin)) << " & " << fitParametersBinned.at(get_parameter_name(customizationType::Sqrt, 1, nJetsBin)) << " $\\pm$ " << fitParameterErrorsBinned.at(get_parameter_name(customizationType::Sqrt, 1, nJetsBin)) << std::fixed;
    //   std::cout << " \\\\ \\hline" << std::endl;
    // }
    // // end tabular environment
    // std::cout << "\\end{tabular}" << std::endl;
    // std::cout << std::endl;

    // print best fit values for linear fit in a LaTeX-formatted table
    std::cout << std::endl;
    std::cout << "Best fit values for linear fit:" << std::endl;
    std::cout << "\\begin{tabular}{|p{0.2\\textwidth}p{0.2\\textwidth}p{0.2\\textwidth}|}" << std::endl;
    std::cout << "  \\hline" << std::endl;
    std::cout << "  Best fit values & $A$ & $m$ \\\\ \\hline" << std::endl;
    for (int nJetsBin = 4; nJetsBin <= 6; ++nJetsBin) {
      if (nJetsBin == 6) std::cout << "  \\nJets $\\geq$ 6";
      else std::cout << "  \\nJets = " << nJetsBin;
      std::cout << std::setprecision(2) << " & " << fitParametersBinned.at(get_parameter_name(customizationType::Slope, 0, nJetsBin)) << " $\\pm$ " << fitParameterErrorsBinned.at(get_parameter_name(customizationType::Slope, 0, nJetsBin)) << " & " << fitParametersBinned.at(get_parameter_name(customizationType::Slope, 1, nJetsBin)) << " $\\pm$ " << fitParameterErrorsBinned.at(get_parameter_name(customizationType::Slope, 1, nJetsBin)) << std::fixed;
      std::cout << " \\\\ \\hline" << std::endl;
    }
    // end tabular environment
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << std::endl;

    // write parameters for binned fit
    std::ofstream fitParametersBinnedFile((options.outputFolder + "/binned_fitParameters_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".dat").c_str());
    assert(fitParametersBinnedFile.is_open());
    std::string fit_parameter_name;
    for (int nJetsBin = (1+options.nJetsNorm); nJetsBin <= 6; ++nJetsBin) {
      for (int customization_type_index = customizationTypeFirst; customization_type_index < static_cast<int>(customizationType::nCustomizationTypes); ++customization_type_index) {
        customizationType customization_type = static_cast<customizationType>(customization_type_index);
        for (int parameter_index = 0; parameter_index < customizationTypeNPars.at(customization_type); ++parameter_index) {
          fit_parameter_name = get_parameter_name(customization_type, parameter_index, nJetsBin);
          fitParametersBinnedFile << "float " << fit_parameter_name << "=" << fitParametersBinned.at(fit_parameter_name) << std::endl;
        }
        for (int eigen_index = 0; eigen_index < customizationTypeNPars.at(customization_type); ++eigen_index) {
          for (int parameter_index = 0; parameter_index < customizationTypeNPars.at(customization_type); ++parameter_index) {
            fit_parameter_name = get_eigencoefficient_name(customization_type, eigen_index, parameter_index, nJetsBin);
            fitParametersBinnedFile << "float " << fit_parameter_name << "=" << fitParametersBinned.at(fit_parameter_name) << std::endl;
          }
          fit_parameter_name = get_eigenerror_name(customization_type, eigen_index, nJetsBin);
          fitParametersBinnedFile << "float " << fit_parameter_name << "=" << fitParametersBinned.at(fit_parameter_name) << std::endl;
        }
      }
    }
    fitParametersBinnedFile.close();
    std::cout << "Binned fit parameters written to file: " << (options.outputFolder + "/binned_fitParameters_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".dat") << std::endl;

    // write out adjustment values
    std::ofstream adjustmentsOutputFile((options.outputFolder + "/adjustments_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".dat").c_str());
    assert(adjustmentsOutputFile.is_open());
    for (int adjustments_slope_sqrt_fit_forOutputFile_index = 0; adjustments_slope_sqrt_fit_forOutputFile_index < static_cast<int>(adjustments_slope_sqrt_fit_forOutputFile.size()); ++adjustments_slope_sqrt_fit_forOutputFile_index) {
      adjustmentsOutputFile << adjustments_slope_sqrt_fit_forOutputFile.at(adjustments_slope_sqrt_fit_forOutputFile_index) << std::endl;
    }
    adjustmentsOutputFile.close();
    std::cout << "Scaling adjustments written to file: " << (options.outputFolder + "/adjustments_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".dat") << std::endl;
  }

  std::cout << "All done!" << std::endl;

  // For some strange reason, the executable sometimes crashes at the end of this function
  // with something like:
  // *** Error in `./fitScripts/bin/runFits': free(): invalid pointer: 0x0000000008c3dbb8 ***
  // Workaround: don't check exit status passed to bash, just create an execution status file, set its contents to the number "1"
  // when the script starts, and rewrite contents to "0" at the end... read this status file from any calling script
  std::ofstream executionStatusOutputFileEnd((options.outputFolder + "/execution_status_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".txt").c_str());
  assert(executionStatusOutputFileEnd.is_open());
  executionStatusOutputFileEnd << 0 << std::endl;
  executionStatusOutputFileEnd.close();

  return EXIT_SUCCESS;
}
