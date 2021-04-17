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

double get_fTest_prob(const double& chi2_1, const double& chi2_2, const int& ndf_1, const int& ndf_2) {
  // Step 1: Get the environment variable EOSTMPAREA, a folder where temp files can be stored
  std::string tmpFolder = std::string(getenv("EOSTMPAREA"));
  if (tmpFolder == "") {
    std::cout << "ERROR: env variable \"EOSTMPAREA\" does not appear to be set." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  
  // Step 2: Generate python command that gets this value
  std::string commandToRun = "python -c \"import tmStatsUtils; print(tmStatsUtils.get_fTest_prob(chi2_1=" + std::to_string(chi2_1) + ", chi2_2=" + std::to_string(chi2_2) + ", ndf_1=" + std::to_string(ndf_1) + ", ndf_2=" + std::to_string(ndf_2) + "))\" > " + tmpFolder + "/get_fTest_prob_tmp.txt 2>&1";
  int ftest_command_return_status = system(commandToRun.c_str());
  if (ftest_command_return_status != 0) {
    std::cout << "ERROR in running command: " << commandToRun << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Step 3: Open output file and read in the value
  std::vector<double> valuesFromFile;
  double valueFromFile;
  std::ifstream inputFileObject((tmpFolder + "/get_fTest_prob_tmp.txt").c_str());
  if (inputFileObject.is_open()) {
    while (inputFileObject >> valueFromFile) {
      valuesFromFile.push_back(valueFromFile);
    }
    inputFileObject.close();
  }
  else {
    std::cout << "ERROR: Unable to open file: " << (tmpFolder + "/get_fTest_prob_tmp.txt") << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Step 4: Check that there is indeed just one value
  if (!(valuesFromFile.size() == 1)) {
    std::cout << "ERROR: this tmp file is in an unexpected format: " << (tmpFolder + "/get_fTest_prob_tmp.txt") << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Step 5: Delete the temp file
  std::string file_remove_command = "rm " + (tmpFolder + "/get_fTest_prob_tmp.txt");
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

std::map<int, double> getNominalAdjustmentsFromBinIntegralMaps(const STRegionsStruct &regions, const std::map<int, double> &bin_integrals_divided_by_bin_widths_nominal, const std::map<int, double> &bin_integrals_divided_by_bin_widths_from_2_jets_kernel) {
  std::map<int, double> adjustments;
  for (int regionIndex = 1; regionIndex <= regions.STAxis.GetNbins(); ++regionIndex) {
    adjustments[regionIndex] = (bin_integrals_divided_by_bin_widths_nominal.at(regionIndex))/(bin_integrals_divided_by_bin_widths_from_2_jets_kernel.at(regionIndex));
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

std::map<int, double> get_adjustments_from_slope(const double & slope, const STRegionsStruct & regions, TF1 * unadjusted_tf1, const bool & print_debug=false) {
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

  if (print_debug) std::cout << "get_adjustment_from_slope called with slope=" << slope << std::endl;
  double STNormTargetTmp = 0.5*(regions.STNormRangeMin + regions.STNormRangeMax); // different from options.STNormTarget because this is a different STRegionsStruct
  if (print_debug) std::cout << "STNormTargetTmp = " << STNormTargetTmp << std::endl;
  TF1 *adjusted_tf1 = new TF1((std::string("slope_adjusted_") + unadjusted_tf1->GetName()).c_str(), [&](double *x, double *p){ (void)p; return ((1.0 + slope*(((x[0])/STNormTargetTmp) - 1.0))*(unadjusted_tf1->Eval(x[0]))); }, regions.STNormRangeMin, ST_MAX_RANGE, 0);
  // in the lambda expression above, the (void)p is just to avoid a compilation error with "gcc -Werror=unused-parameter"
  std::map<int, double> adjustments_from_slope;
  for (int regionIndex = 1; regionIndex <= regions.STAxis.GetNbins(); ++regionIndex) {
    double bin_low_edge = regions.STAxis.GetBinLowEdge(regionIndex);
    double bin_up_edge = regions.STAxis.GetBinUpEdge(regionIndex);
    double integral_adjusted = adjusted_tf1->Integral(bin_low_edge, bin_up_edge, TF1_INTEGRAL_REL_TOLERANCE);
    double integral_unadjusted = unadjusted_tf1->Integral(bin_low_edge, bin_up_edge, TF1_INTEGRAL_REL_TOLERANCE);
    adjustments_from_slope[regionIndex] = integral_adjusted/integral_unadjusted;
    if (print_debug) std::cout << "At regionIndex: " << regionIndex << ", bin_low_edge: " << bin_low_edge << ", bin_up_edge: " << bin_up_edge << ", integral_adjusted: " << integral_adjusted << ", integral_unadjusted: " << integral_unadjusted << ", adjustment: " << adjustments_from_slope.at(regionIndex) << std::endl;
  }
  delete adjusted_tf1;
  return adjustments_from_slope;
}

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  TH1::AddDirectory(kFALSE);
  RooMsgService::instance().setGlobalKillBelow(MsgLevel::WARNING);
  ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(TF1_INTEGRAL_REL_TOLERANCE);
  do_sanity_checks_customizationTypes();

  tmArgumentParser argumentParser = tmArgumentParser("Run script that prints useful info about the normalization.");
  argumentParser.addArgument("sourceFilePath", "", true, "Path to file containing list of paths with n-tuplized events.");
  argumentParser.addArgument("outputFolder", "", true, "Output folder.");
  argumentParser.addArgument("selection", "", true, "Name of selection: \"singlemedium\", \"signal_loose\", etc.");
  argumentParser.addArgument("identifier", "", true, "Identifier: \"MC_GJet17\", \"MC_GJet\", etc.");
  argumentParser.addArgument("STBoundariesSourceFile", "STRegionBoundaries_normOptimization.dat", false, "Source file for reading in ST region boundaries.");
  argumentParser.addArgument("yearString", "all", false, "String with year: can take values \"2016\", \"2017\", \"2018\", or \"all\".");
  argumentParser.addArgument("PDF_nSTBins", "25", false, "Number of bins for plotting datasets.");
  argumentParser.addArgument("preNormalizationBuffer", "200.0", false, "Buffer in ST to use before normalization bin for the kernel.");
  argumentParser.addArgument("readParametersFromFiles", "/dev/null,/dev/null,/dev/null", false, "If this argument is set, then no fits are performed; instead, the fit parameters is read in from the file locations given as the value of this argument. This should be a list of precisely three files separated by a comma: in order, the unbinned parameters, the binned parameters, and a file containing ST region boundaries to use for saving the (observed/best-fit) ratio adjustments.");
  argumentParser.addArgument("plotConcise", "false", false, "If this argument is set, then only the (linear+sqrt) fit and associated errors are plotted.");
  argumentParser.setPassedStringValues(argc, argv);
  optionsStruct options = getOptionsFromParser(argumentParser);
  std::cout << "Options passed:" << std::endl << options << std::endl;

  int mkdir_return_status = system(("set -x && mkdir -p " + options.outputFolder + " && set +x").c_str());
  if (mkdir_return_status != 0) {
    std::cout << "ERROR in creating output folder with path: " << options.outputFolder << std::endl;
    std::exit(EXIT_FAILURE);
  }

  TFile* outputFile = TFile::Open((options.outputFolder + "/dataSetStorage_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".root").c_str(), "RECREATE");
  if ((!(outputFile->IsOpen())) or (outputFile->IsZombie())) {
    std::cout << "ERROR: unable to open file with name: " << (options.outputFolder + "/dataSetStorage_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".root") << std::endl;
    std::exit(EXIT_FAILURE);
  }
  RooAbsData::setDefaultStorageType(RooAbsData::Tree);
  RooRealVar rooVar_ST("roo_ST", "roo_ST", (options.STRegions.STNormRangeMin - options.preNormalizationBuffer), ST_MAX_RANGE, "GeV");
  RooRealVar rooVar_weight("roo_weight", "roo_weight", 1., 0., 100000.);

  RooDataSet STDataSet_2Jets = RooDataSet("STDataSet_2JetsBin", "STDataSet_2JetsBin", RooArgSet(rooVar_ST, rooVar_weight), WeightVar(rooVar_weight));
  RooDataSet STDataSet_3Jets = RooDataSet("STDataSet_3JetsBin", "STDataSet_3JetsBin", RooArgSet(rooVar_ST, rooVar_weight), WeightVar(rooVar_weight));
  RooDataSet STDataSet_4Jets = RooDataSet("STDataSet_4JetsBin", "STDataSet_4JetsBin", RooArgSet(rooVar_ST, rooVar_weight), WeightVar(rooVar_weight));
  RooDataSet STDataSet_5Jets = RooDataSet("STDataSet_5JetsBin", "STDataSet_5JetsBin", RooArgSet(rooVar_ST, rooVar_weight), WeightVar(rooVar_weight));
  RooDataSet STDataSet_6Jets = RooDataSet("STDataSet_6JetsBin", "STDataSet_6JetsBin", RooArgSet(rooVar_ST, rooVar_weight), WeightVar(rooVar_weight));
  std::map<int, RooDataSet*> STDataSets = {
    {2, &(STDataSet_2Jets)},
    {3, &(STDataSet_3Jets)},
    {4, &(STDataSet_4Jets)},
    {5, &(STDataSet_5Jets)},
    {6, &(STDataSet_6Jets)}
  };
  // idiotic, I know, but the following compiles and results in a segfault that I've been unable to debug:
  // std::map<int, RooDataSet> STDataSets;
  // for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
  //   STDataSets[nJetsBin] = RooDataSet(("STDataSet_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("STDataSet_" + std::to_string(nJetsBin) + "JetsBin").c_str(), RooArgSet(rooVar_ST, rooVar_weight), WeightVar(rooVar_weight));
  // }

  std::map<int, TH1D> STHistograms;
  for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
    STHistograms[nJetsBin] = TH1D(("STHistogram_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("ST distribution, " + std::to_string(nJetsBin) + " Jets bin;ST(GeV);weighted events/bin").c_str(), (options.STRegions.STBoundaries.size()-1), &(options.STRegions.STBoundaries.at(0)));
    (STHistograms.at(nJetsBin)).Sumw2();
  }

  TFile *sourceFile = TFile::Open((options.sourceFilePath).c_str(), "READ");
  if ((!(sourceFile->IsOpen())) or (sourceFile->IsZombie())) {
    std::cout << "ERROR: unable to open file with name: " << options.sourceFilePath << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::cout << "Fetching ST trees from file: " << options.sourceFilePath << std::endl;
  for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
    TTree* STTree = new TTree();
    std::string inputTreeName = ("STTree_" + std::to_string(nJetsBin) + "JetsBin");
    std::cout << "Getting tree with name: " << inputTreeName << std::endl;
    sourceFile->GetObject(inputTreeName.c_str(), STTree);
    assert(STTree != nullptr);
    STTree->SetName(inputTreeName.c_str());
    TTreeReader inputTreeReader(STTree);
    TTreeReaderValue<double> ST(inputTreeReader, "ST");
    TTreeReaderValue<double> weight(inputTreeReader, "weight");
    while (inputTreeReader.Next()) {
      rooVar_ST.setVal(*ST);
      double eventWeight = *weight;
      if (((*ST) < (options.STRegions.STNormRangeMin - options.preNormalizationBuffer)) || ((*ST) > ST_MAX_RANGE)) continue;
      (STDataSets.at(nJetsBin))->add(RooArgSet(rooVar_ST), eventWeight);
      if ((*ST) < options.STRegions.STNormRangeMin) continue; // no "pre-norm buffer" needed for histograms
      double binWidth = (STHistograms.at(nJetsBin)).GetXaxis()->GetBinWidth((STHistograms.at(nJetsBin)).GetXaxis()->FindFixBin(*ST));
      (STHistograms.at(nJetsBin)).Fill(*ST, eventWeight/binWidth);
    }
    (STDataSets.at(nJetsBin))->Print();
  }

  // A few useful initializations
  std::map<std::string, std::map<int, double> > fitParametersUnbinned;
  std::vector<std::string> fitParametersUnbinnedList;
  // std::map<std::string, std::map<int, double> > fitParametersBinned;
  // std::vector<std::string> fitParametersBinnedList;
  // std::map<std::string, std::map<int, goodnessOfFitStruct> > fit_qualities_binned;
  std::map<std::string, double> fitParametersBinned;
  std::map<customizationType, std::map<int, double> > fit_pvalues;
  std::map<std::string, std::map<int, double> > ftest_pValues;
  std::vector<std::string> adjustments_slope_sqrt_fit_forOutputFile;
  std::vector<std::string> ratio_adjustments_forOutputFile;
  customizationType customization_type_for_adjustments_output = customizationType::Sqrt;
  customizationType customization_type_denominator_for_ratios = customizationType::ScaleOnly;
  double scale_minVal = 0.0;
  double scale_maxVal = 5.0;
  double slope_minVal = -1.0/(((ST_MAX_RANGE)/(options.STNormTarget)) - 1.0);
  double slope_maxVal = 5.0;
  double sqrt_minVal = -1.0/(std::sqrt((ST_MAX_RANGE)/(options.STNormTarget)) - 1.0);
  double sqrt_maxVal = 25.0;
  double quad_minVal = -1.0/((std::pow((ST_MAX_RANGE)/(options.STNormTarget), 2)) - 1.0);
  double quad_maxVal = 3.0;
  TRandom3 *random_generator = new TRandom3(1234); // for repeatability
  // for (int boundaryIndex = 0; boundaryIndex < static_cast<int>((options.STRegions.STBoundaries.size() - 1)); ++boundaryIndex) {
  //   double thisBoundary = options.STRegions.STBoundaries.at(boundaryIndex);
  //   double nextBoundary = options.STRegions.STBoundaries.at(1+boundaryIndex);
  //   rooVar_ST.setRange(("range_STBinIndex_" + std::to_string(1+boundaryIndex)).c_str(), thisBoundary, nextBoundary);
  // }
  rooVar_ST.setRange("normRange", options.STRegions.STNormRangeMin, options.STRegions.STNormRangeMax);
  rooVar_ST.setRange("fitRange", options.STRegions.STNormRangeMin, ST_MAX_RANGE);
  rooVar_ST.setRange("plotRange", (options.STRegions.STNormRangeMin - options.preNormalizationBuffer), ST_MAX_RANGE);

  // If fits are to be read in from another file, do it now
  if (options.readParametersFromFiles) {
    std::string lineFromFile;

    std::ifstream inputFileObject_unbinned(options.inputUnbinnedParametersFileName.c_str());
    assert(inputFileObject_unbinned.is_open());
    for (int nJetsBin = 3; nJetsBin <= 6; ++nJetsBin) {
      assert(getline(inputFileObject_unbinned, lineFromFile)); fitParametersUnbinned["slope_fit_slope"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_fit_slope_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_unbinned, lineFromFile)); fitParametersUnbinned["slope_fit_slopeError"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_fit_slopeError_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_unbinned, lineFromFile)); fitParametersUnbinned["sqrt_fit_sqrt"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "sqrt_fit_sqrt_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_unbinned, lineFromFile)); fitParametersUnbinned["sqrt_fit_sqrtError"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "sqrt_fit_sqrtError_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_unbinned, lineFromFile)); fitParametersUnbinned["slope_sqrt_fit_slope"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_fit_slope_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_unbinned, lineFromFile)); fitParametersUnbinned["slope_sqrt_fit_sqrt"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_fit_sqrt_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_unbinned, lineFromFile)); fitParametersUnbinned["slope_sqrt_fit_mode1_slopeCoefficient"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_fit_mode1_slopeCoefficient_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_unbinned, lineFromFile)); fitParametersUnbinned["slope_sqrt_fit_mode1_sqrtCoefficient"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_fit_mode1_sqrtCoefficient_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_unbinned, lineFromFile)); fitParametersUnbinned["slope_sqrt_fit_mode1_error"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_fit_mode1_error_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_unbinned, lineFromFile)); fitParametersUnbinned["slope_sqrt_fit_mode2_slopeCoefficient"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_fit_mode2_slopeCoefficient_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_unbinned, lineFromFile)); fitParametersUnbinned["slope_sqrt_fit_mode2_sqrtCoefficient"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_fit_mode2_sqrtCoefficient_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_unbinned, lineFromFile)); fitParametersUnbinned["slope_sqrt_fit_mode2_error"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_fit_mode2_error_" + std::to_string(nJetsBin) + "Jets");
    }
    assert(!(getline(inputFileObject_unbinned, lineFromFile))); // makes sure there's nothing else in the input file
    inputFileObject_unbinned.close();

    std::ifstream inputFileObject_binned(options.inputBinnedParametersFileName.c_str());
    assert(inputFileObject_binned.is_open());
    std::string fit_parameter_name;
    for (int nJetsBin = 3; nJetsBin <= 6; ++nJetsBin) {
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

  // First get the 2-jets RooKeysPdf
  RooKeysPdf pdf_2Jets = RooKeysPdf("pdf_2Jets", "pdf_2Jets", rooVar_ST, *(STDataSets.at(2)), RooKeysPdf::MirrorLeft, 1.5);

  // Plot 2-jets shape and dataset
  TCanvas unbinned_pdfCanvas_2Jets = TCanvas("c_dataSetAndPdf_unbinned_2Jets", "c_dataSetAndPdf_unbinned_2Jets", 2560, 1440);
  RooPlot* rooFrame = rooVar_ST.frame();
  (STDataSets.at(2))->plotOn(rooFrame, Binning(options.PDF_nSTBins, (options.STRegions.STNormRangeMin - options.preNormalizationBuffer), ST_MAX_RANGE));
  pdf_2Jets.plotOn(rooFrame, NormRange("normRange"), LineColor(static_cast<EColor>(kBlue)), LineWidth(1));
  rooFrame->Draw();
  unbinned_pdfCanvas_2Jets.Update();
  rooFrame->SetMinimum((rooFrame->GetMaximum())/10000.);
  gPad->SetLogy();
  unbinned_pdfCanvas_2Jets.Update();
  unbinned_pdfCanvas_2Jets.SaveAs((options.outputFolder + "/unbinned_pdfAndData_2JetsBin_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".pdf").c_str());

  TCanvas binned_pdfCanvas_2Jets = TCanvas("c_dataSetAndPdf_binned_2Jets", "c_dataSetAndPdf_binned_2Jets", 2560, 1440);
  gStyle->SetOptStat(0);
  (STHistograms.at(2)).Draw();
  (STHistograms.at(2)).GetYaxis()->SetRange(((STHistograms.at(2)).GetMaximum())/10000., ((STHistograms.at(2)).GetMaximum()));
  binned_pdfCanvas_2Jets.Update();
  customizedPDF pdf_2Jets_customized(&pdf_2Jets, &rooVar_ST, options.STNormTarget, customizationType::ScaleOnly);
  pdf_2Jets_customized.setNominalScale("fitRange", ((STHistograms.at(2)).Integral(1, (STHistograms.at(2)).GetXaxis()->GetNbins(), "width")));
  TF1 pdf_2Jets_customized_TF1 = TF1("pdf_2Jets_customized_TF1", pdf_2Jets_customized, options.STRegions.STNormRangeMin, ST_MAX_RANGE, 1);
  pdf_2Jets_customized_TF1.SetParameter(0, 1.0);
  pdf_2Jets_customized_TF1.SetLineColor(static_cast<EColor>(kBlue));
  pdf_2Jets_customized_TF1.SetLineWidth(2);
  pdf_2Jets_customized_TF1.Draw("CSAME");
  binned_pdfCanvas_2Jets.Update();
  gPad->SetLogy();
  binned_pdfCanvas_2Jets.Update();
  binned_pdfCanvas_2Jets.SaveAs((options.outputFolder + "/binned_pdfAndData_2JetsBin_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".pdf").c_str());

  for (int nJetsBin = 3; nJetsBin <= 6; ++nJetsBin) {
    printSeparator();
    std::cout << "Finding fits at nJetsBin = " << nJetsBin << std::endl;
    std::cout << "First the unbinned analysis:" << std::endl;

    // unadjusted
    TLegend legend_dataSetsAndPdf_unbinned = TLegend(0.6, 0.7, 0.9, 0.9);
    RooPlot* rooFrame = rooVar_ST.frame();
    (STDataSets.at(nJetsBin))->plotOn(rooFrame, Binning(options.PDF_nSTBins, (options.STRegions.STNormRangeMin - options.preNormalizationBuffer), ST_MAX_RANGE));
    pdf_2Jets.plotOn(rooFrame, NormRange("normRange"), LineColor(static_cast<EColor>(kBlue)), LineWidth(1));
    // RooKeysPdf *pdf_2Jets_copyForNLL = (RooKeysPdf*)(pdf_2Jets.Clone());
    rooVar_ST.setRange((options.STRegions.STNormRangeMin - options.preNormalizationBuffer), ST_MAX_RANGE);
    TLegendEntry *legendEntry_unbinned_unadjusted_slope = legend_dataSetsAndPdf_unbinned.AddEntry(&pdf_2Jets, "2 jets kernel, unadjusted");
    legendEntry_unbinned_unadjusted_slope->SetMarkerColor(static_cast<EColor>(kBlue));
    legendEntry_unbinned_unadjusted_slope->SetLineColor(static_cast<EColor>(kBlue));
    legendEntry_unbinned_unadjusted_slope->SetTextColor(static_cast<EColor>(kBlue));
    rooVar_ST.setVal(options.STNormTarget);
    double pdf_2Jets_at_STNorm = pdf_2Jets.getVal(rooVar_ST);

    RooKeysPdf pdf_nJets_kernel = RooKeysPdf(("pdf_nJets_kernel_at_" + std::to_string(nJetsBin) + "Jets").c_str(), ("pdf_nJets_kernel_at_" + std::to_string(nJetsBin) + "Jets").c_str(), rooVar_ST, *(STDataSets.at(nJetsBin)), RooKeysPdf::MirrorLeft, 1.5);
    TGraph ratioGraph_unbinned_nJetsKernelToUnadjusted = TGraph();
    ratioGraph_unbinned_nJetsKernelToUnadjusted.SetName(("ratioGraph_unbinned_nJetsKernelToUnadjusted_at" + std::to_string(nJetsBin) + "Jets").c_str());
    ratioGraph_unbinned_nJetsKernelToUnadjusted.SetTitle(("kernel from data at " + std::to_string(nJetsBin) + " Jets / unadjusted").c_str());
    rooVar_ST.setVal(options.STNormTarget);
    double pdf_nJets_kernel_at_STNorm = pdf_nJets_kernel.getVal(rooVar_ST);
    for (int STCounter = 0; STCounter <= 1000; ++STCounter) {
      double STVal = options.STRegions.STNormRangeMin + (1.0*STCounter/1000)*(ST_MAX_RANGE - options.STRegions.STNormRangeMin);
      rooVar_ST.setVal(STVal);
      double ratio_nJetsKernelToUnadjusted = (pdf_nJets_kernel.getVal(rooVar_ST)/(pdf_nJets_kernel_at_STNorm))/((pdf_2Jets.getVal(rooVar_ST))/(pdf_2Jets_at_STNorm));
      ratioGraph_unbinned_nJetsKernelToUnadjusted.SetPoint(STCounter, STVal, ratio_nJetsKernelToUnadjusted);
    }

    // slope correction only
    std::string slopeVar_name = "rooVar_slope_" + std::to_string(nJetsBin) + "JetsBin";
    RooRealVar rooVar_slope(slopeVar_name.c_str(), slopeVar_name.c_str(), 0.0, slope_minVal, slope_maxVal);
    std::string function_slopeAdjustment = "1.0 + (" + slopeVar_name + "*((roo_ST/" + std::to_string(options.STNormTarget) + ") - 1.0))";
    RooGenericPdf pdf_slopeAdjustment("slopeAdjustment", "slopeAdjustment", function_slopeAdjustment.c_str(), RooArgSet(rooVar_ST, rooVar_slope));
    RooProdPdf pdf_nJets_adjusted_slopeOnly(("pdf_2Jets_slopeOnlyAdjustment_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("pdf_2Jets_slopeOnlyAdjustment_" + std::to_string(nJetsBin) + "JetsBin").c_str(), RooArgSet(pdf_2Jets, pdf_slopeAdjustment));
    // double nll_slopeOnly;
    if (options.readParametersFromFiles) {
      rooVar_slope.setVal((fitParametersUnbinned.at("slope_fit_slope")).at(nJetsBin));
      // rooVar_slope.setError((fitParametersUnbinned.at("slope_fit_slopeError")).at(nJetsBin));
      // nll_slopeOnly = (pdf_nJets_adjusted_slopeOnly.createNLL(*(STDataSets.at(nJetsBin)), Verbose(kFALSE), Range(options.STRegions.STNormRangeMin, ST_MAX_RANGE)))->getVal();
    }
    else {
      RooFitResult *fitResult = pdf_nJets_adjusted_slopeOnly.fitTo(*(STDataSets.at(nJetsBin)), Range(options.STRegions.STNormRangeMin, ST_MAX_RANGE), Optimize(kFALSE), PrintLevel(-1), Verbose(kFALSE), PrintEvalErrors(2), SumW2Error(kTRUE), Save(kTRUE));
      assert(fitResult->status() == 0);
      // nll_slopeOnly = fitResult->minNll();
      fitParametersUnbinned["slope_fit_slope"][nJetsBin] = rooVar_slope.getValV();
      fitParametersUnbinnedList.push_back(std::string("float slope_fit_slope_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("slope_fit_slope")).at(nJetsBin))));
      fitParametersUnbinned["slope_fit_slopeError"][nJetsBin] = rooVar_slope.getError();
      fitParametersUnbinnedList.push_back(std::string("float slope_fit_slopeError_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("slope_fit_slopeError")).at(nJetsBin))));
    }
    rooVar_slope.Print();
    // std::map<int, double> integralsOverBins_slopeOnly;
    // for (int dataHistBinIndex = 1; dataHistBinIndex <= dataHist.GetXaxis()->GetNbins(); ++dataHistBinIndex) {
    //   rooVar_ST.setRange((options.STRegions.STNormRangeMin - options.preNormalizationBuffer), ST_MAX_RANGE);
    //   integralsOverBins_slopeOnly[dataHistBinIndex] = (pdf_nJets_adjusted_slopeOnly.createIntegral(rooVar_ST, NormSet(rooVar_ST), Range(("range_STBinIndex_" + std::to_string(dataHistBinIndex)).c_str())))->getVal();
    // }
    // totalIntegral = (pdf_nJets_adjusted_slopeOnly.createIntegral(rooVar_ST, NormSet(rooVar_ST), Range("plotRange")))->getVal();
    // assert(std::fabs(totalIntegral - 1.0) < CHECK_TOLERANCE);
    pdf_nJets_adjusted_slopeOnly.plotOn(rooFrame, NormRange("normRange"), LineColor(static_cast<EColor>(kRed+1)), LineWidth(1));
    // plus and minus one-sigma plotted with dashed linestyle
    rooVar_slope.setVal((fitParametersUnbinned.at("slope_fit_slope")).at(nJetsBin) + (fitParametersUnbinned.at("slope_fit_slopeError")).at(nJetsBin));
    pdf_nJets_adjusted_slopeOnly.plotOn(rooFrame, NormRange("normRange"), LineColor(static_cast<EColor>(kRed+1)), LineStyle(kDashed), LineWidth(1));
    rooVar_slope.setVal((fitParametersUnbinned.at("slope_fit_slope")).at(nJetsBin) - (fitParametersUnbinned.at("slope_fit_slopeError")).at(nJetsBin));
    pdf_nJets_adjusted_slopeOnly.plotOn(rooFrame, NormRange("normRange"), LineColor(static_cast<EColor>(kRed+1)), LineStyle(kDashed), LineWidth(1));
    rooVar_slope.setVal((fitParametersUnbinned.at("slope_fit_slope")).at(nJetsBin));
    TLegendEntry *legendEntry_unbinned_adjusted_slopeOnly = legend_dataSetsAndPdf_unbinned.AddEntry(&pdf_nJets_adjusted_slopeOnly, "2 jets kernel + slope adjustment");
    legendEntry_unbinned_adjusted_slopeOnly->SetMarkerColor(static_cast<EColor>(kRed+1));
    legendEntry_unbinned_adjusted_slopeOnly->SetLineColor(static_cast<EColor>(kRed+1));
    legendEntry_unbinned_adjusted_slopeOnly->SetTextColor(static_cast<EColor>(kRed+1));
    TGraph ratioGraph_unbinned_slopeOnlyToUnadjusted = TGraph();
    TGraph ratioGraph_unbinned_slopeOnlyToUnadjusted_high_estimate = TGraph();
    TGraph ratioGraph_unbinned_slopeOnlyToUnadjusted_low_estimate = TGraph();
    ratioGraph_unbinned_slopeOnlyToUnadjusted.SetName(("ratioGraph_unbinned_slopeOnlyToUnadjusted_at" + std::to_string(nJetsBin) + "Jets").c_str());
    ratioGraph_unbinned_slopeOnlyToUnadjusted.SetTitle(("slope-only fit at " + std::to_string(nJetsBin) + " Jets / unadjusted").c_str());
    rooVar_ST.setVal(options.STNormTarget);
    double pdf_nJets_adjusted_slopeOnly_at_STNorm = pdf_nJets_adjusted_slopeOnly.getVal(rooVar_ST);
    for (int STCounter = 0; STCounter <= 1000; ++STCounter) {
      double STVal = options.STRegions.STNormRangeMin + (1.0*STCounter/1000)*(ST_MAX_RANGE - options.STRegions.STNormRangeMin);
      rooVar_ST.setVal(STVal);
      rooVar_slope.setVal((fitParametersUnbinned.at("slope_fit_slope")).at(nJetsBin));
      double pdf_nominal = (pdf_nJets_adjusted_slopeOnly.getVal(rooVar_ST)/(pdf_nJets_adjusted_slopeOnly_at_STNorm));
      rooVar_slope.setVal((fitParametersUnbinned.at("slope_fit_slope")).at(nJetsBin) + (fitParametersUnbinned.at("slope_fit_slopeError")).at(nJetsBin));
      double pdf_plus_one_sigma = (pdf_nJets_adjusted_slopeOnly.getVal(rooVar_ST)/(pdf_nJets_adjusted_slopeOnly_at_STNorm));
      rooVar_slope.setVal((fitParametersUnbinned.at("slope_fit_slope")).at(nJetsBin) - (fitParametersUnbinned.at("slope_fit_slopeError")).at(nJetsBin));
      double pdf_minus_one_sigma = (pdf_nJets_adjusted_slopeOnly.getVal(rooVar_ST)/(pdf_nJets_adjusted_slopeOnly_at_STNorm));
      double pdf_higher = std::max({pdf_plus_one_sigma, pdf_nominal, pdf_minus_one_sigma});
      double pdf_lower = std::min({pdf_plus_one_sigma, pdf_nominal, pdf_minus_one_sigma});
      double commonDenominator = ((pdf_2Jets.getVal(rooVar_ST))/(pdf_2Jets_at_STNorm));
      double ratio_slopeOnlyToUnadjusted = pdf_nominal/commonDenominator;
      double ratio_higher = pdf_higher/commonDenominator;
      double ratio_lower = pdf_lower/commonDenominator;
      ratioGraph_unbinned_slopeOnlyToUnadjusted.SetPoint(STCounter, STVal, ratio_slopeOnlyToUnadjusted);
      ratioGraph_unbinned_slopeOnlyToUnadjusted_high_estimate.SetPoint(STCounter, STVal, ratio_higher);
      ratioGraph_unbinned_slopeOnlyToUnadjusted_low_estimate.SetPoint(STCounter, STVal, ratio_lower);
    }

    // sqrt correction only
    std::string sqrtVar_name = "rooVar_sqrt_" + std::to_string(nJetsBin) + "JetsBin";
    RooRealVar rooVar_sqrt(sqrtVar_name.c_str(), sqrtVar_name.c_str(), 0.0, sqrt_minVal, sqrt_maxVal);
    std::string function_sqrtAdjustment = "1.0 + (" + sqrtVar_name + "*(sqrt(roo_ST/" + std::to_string(options.STNormTarget) + ") - 1.0))";
    RooGenericPdf pdf_sqrtAdjustment("sqrtAdjustment", "sqrtAdjustment", function_sqrtAdjustment.c_str(), RooArgSet(rooVar_ST, rooVar_sqrt));
    RooProdPdf pdf_nJets_adjusted_sqrtOnly(("pdf_2Jets_sqrtOnlyAdjustment_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("pdf_2Jets_sqrtOnlyAdjustment" + std::to_string(nJetsBin) + "JetsBin").c_str(), RooArgSet(pdf_2Jets, pdf_sqrtAdjustment));
    // double nll_sqrtOnly;
    if (options.readParametersFromFiles) {
      rooVar_sqrt.setVal((fitParametersUnbinned.at("sqrt_fit_sqrt")).at(nJetsBin));
      // rooVar_sqrt.setError((fitParametersUnbinned.at("sqrt_fit_sqrtError")).at(nJetsBin));
      // nll_sqrtOnly = (pdf_nJets_adjusted_sqrtOnly.createNLL(*(STDataSets.at(nJetsBin)), Verbose(kFALSE), Range(options.STRegions.STNormRangeMin, ST_MAX_RANGE)))->getVal();
    }
    else {
      RooFitResult *fitResult = pdf_nJets_adjusted_sqrtOnly.fitTo(*(STDataSets.at(nJetsBin)), Range(options.STRegions.STNormRangeMin, ST_MAX_RANGE), Optimize(kFALSE), PrintLevel(-1), Verbose(kFALSE), PrintEvalErrors(2), SumW2Error(kTRUE), Save(kTRUE));
      assert(fitResult->status() == 0);
      // nll_sqrtOnly = fitResult->minNll();
      fitParametersUnbinned["sqrt_fit_sqrt"][nJetsBin] = rooVar_sqrt.getValV();
      fitParametersUnbinnedList.push_back(std::string("float sqrt_fit_sqrt_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("sqrt_fit_sqrt")).at(nJetsBin))));
      fitParametersUnbinned["sqrt_fit_sqrtError"][nJetsBin] = rooVar_sqrt.getError();
      fitParametersUnbinnedList.push_back(std::string("float sqrt_fit_sqrtError_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("sqrt_fit_sqrtError")).at(nJetsBin))));
    }
    rooVar_sqrt.Print();
    // std::map<int, double> integralsOverBins_sqrtOnly;
    // for (int dataHistBinIndex = 1; dataHistBinIndex <= dataHist.GetXaxis()->GetNbins(); ++dataHistBinIndex) {
    //   rooVar_ST.setRange((options.STRegions.STNormRangeMin - options.preNormalizationBuffer), ST_MAX_RANGE);
    //   integralsOverBins_sqrtOnly[dataHistBinIndex] = (pdf_nJets_adjusted_sqrtOnly.createIntegral(rooVar_ST, NormSet(rooVar_ST), Range(("range_STBinIndex_" + std::to_string(dataHistBinIndex)).c_str())))->getVal();
    // }
    // totalIntegral = (pdf_nJets_adjusted_sqrtOnly.createIntegral(rooVar_ST, NormSet(rooVar_ST), Range("plotRange")))->getVal();
    // assert(std::fabs(totalIntegral - 1.0) < CHECK_TOLERANCE);
    pdf_nJets_adjusted_sqrtOnly.plotOn(rooFrame, NormRange("normRange"), LineColor(static_cast<EColor>(kGreen+3)), LineWidth(1));
    // plus and minus one-sigma plotted with dashed linestyle
    rooVar_sqrt.setVal((fitParametersUnbinned.at("sqrt_fit_sqrt")).at(nJetsBin) + (fitParametersUnbinned.at("sqrt_fit_sqrtError")).at(nJetsBin));
    pdf_nJets_adjusted_sqrtOnly.plotOn(rooFrame, NormRange("normRange"), LineColor(static_cast<EColor>(kGreen+3)), LineStyle(kDashed), LineWidth(1));
    rooVar_sqrt.setVal((fitParametersUnbinned.at("sqrt_fit_sqrt")).at(nJetsBin) - (fitParametersUnbinned.at("sqrt_fit_sqrtError")).at(nJetsBin));
    pdf_nJets_adjusted_sqrtOnly.plotOn(rooFrame, NormRange("normRange"), LineColor(static_cast<EColor>(kGreen+3)), LineStyle(kDashed), LineWidth(1));
    rooVar_sqrt.setVal((fitParametersUnbinned.at("sqrt_fit_sqrt")).at(nJetsBin));
    TLegendEntry *legendEntry_unbinned_adjusted_sqrtOnly = legend_dataSetsAndPdf_unbinned.AddEntry(&pdf_nJets_adjusted_sqrtOnly, "2 jets kernel + sqrt adjustment");
    legendEntry_unbinned_adjusted_sqrtOnly->SetMarkerColor(static_cast<EColor>(kGreen+3));
    legendEntry_unbinned_adjusted_sqrtOnly->SetLineColor(static_cast<EColor>(kGreen+3));
    legendEntry_unbinned_adjusted_sqrtOnly->SetTextColor(static_cast<EColor>(kGreen+3));
    TGraph ratioGraph_unbinned_sqrtOnlyToUnadjusted = TGraph();
    TGraph ratioGraph_unbinned_sqrtOnlyToUnadjusted_high_estimate = TGraph();
    TGraph ratioGraph_unbinned_sqrtOnlyToUnadjusted_low_estimate = TGraph();
    ratioGraph_unbinned_sqrtOnlyToUnadjusted.SetName(("ratioGraph_unbinned_sqrtOnlyToUnadjusted_at" + std::to_string(nJetsBin) + "Jets").c_str());
    ratioGraph_unbinned_sqrtOnlyToUnadjusted.SetTitle(("sqrt-only fit at " + std::to_string(nJetsBin) + " Jets / unadjusted").c_str());
    rooVar_ST.setVal(options.STNormTarget);
    double pdf_nJets_adjusted_sqrtOnly_at_STNorm = pdf_nJets_adjusted_sqrtOnly.getVal(rooVar_ST);
    for (int STCounter = 0; STCounter <= 1000; ++STCounter) {
      double STVal = options.STRegions.STNormRangeMin + (1.0*STCounter/1000)*(ST_MAX_RANGE - options.STRegions.STNormRangeMin);
      rooVar_ST.setVal(STVal);
      rooVar_sqrt.setVal((fitParametersUnbinned.at("sqrt_fit_sqrt")).at(nJetsBin));
      double pdf_nominal = (pdf_nJets_adjusted_sqrtOnly.getVal(rooVar_ST)/(pdf_nJets_adjusted_sqrtOnly_at_STNorm));
      rooVar_sqrt.setVal((fitParametersUnbinned.at("sqrt_fit_sqrt")).at(nJetsBin) + (fitParametersUnbinned.at("sqrt_fit_sqrtError")).at(nJetsBin));
      double pdf_plus_one_sigma = (pdf_nJets_adjusted_sqrtOnly.getVal(rooVar_ST)/(pdf_nJets_adjusted_sqrtOnly_at_STNorm));
      rooVar_sqrt.setVal((fitParametersUnbinned.at("sqrt_fit_sqrt")).at(nJetsBin) - (fitParametersUnbinned.at("sqrt_fit_sqrtError")).at(nJetsBin));
      double pdf_minus_one_sigma = (pdf_nJets_adjusted_sqrtOnly.getVal(rooVar_ST)/(pdf_nJets_adjusted_sqrtOnly_at_STNorm));
      double pdf_higher = std::max({pdf_plus_one_sigma, pdf_nominal, pdf_minus_one_sigma});
      double pdf_lower = std::min({pdf_plus_one_sigma, pdf_nominal, pdf_minus_one_sigma});
      double commonDenominator = ((pdf_2Jets.getVal(rooVar_ST))/(pdf_2Jets_at_STNorm));
      double ratio_sqrtOnlyToUnadjusted = pdf_nominal/commonDenominator;
      double ratio_higher = pdf_higher/commonDenominator;
      double ratio_lower = pdf_lower/commonDenominator;
      ratioGraph_unbinned_sqrtOnlyToUnadjusted.SetPoint(STCounter, STVal, ratio_sqrtOnlyToUnadjusted);
      ratioGraph_unbinned_sqrtOnlyToUnadjusted_high_estimate.SetPoint(STCounter, STVal, ratio_higher);
      ratioGraph_unbinned_sqrtOnlyToUnadjusted_low_estimate.SetPoint(STCounter, STVal, ratio_lower);
    }

    // slope correction + sqrt correction
    std::string slopeVar_combinedFit_name = "rooVar_slope_combinedFit_" + std::to_string(nJetsBin) + "JetsBin";
    RooRealVar rooVar_slope_combinedFit(slopeVar_combinedFit_name.c_str(), slopeVar_combinedFit_name.c_str(), 0., slope_minVal, slope_maxVal);
    function_slopeAdjustment = "1.0 + (" + slopeVar_combinedFit_name + "*((roo_ST/" + std::to_string(options.STNormTarget) + ") - 1.0))";
    RooGenericPdf pdf_slopeAdjustment_combinedFit("slopeAdjustment_combinedFit", "slopeAdjustment_combinedFit", function_slopeAdjustment.c_str(), RooArgSet(rooVar_ST, rooVar_slope_combinedFit));
    std::string sqrtVar_combinedFit_name = "rooVar_sqrt_combinedFit_" + std::to_string(nJetsBin) + "JetsBin";
    RooRealVar rooVar_sqrt_combinedFit(sqrtVar_combinedFit_name.c_str(), sqrtVar_combinedFit_name.c_str(), (fitParametersUnbinned.at("sqrt_fit_sqrt")).at(nJetsBin), sqrt_minVal, sqrt_maxVal);
    function_sqrtAdjustment = "1.0 + (" + sqrtVar_combinedFit_name + "*(sqrt(roo_ST/" + std::to_string(options.STNormTarget) + ") - 1.0))";
    RooGenericPdf pdf_sqrtAdjustment_combinedFit("sqrtAdjustment_combinedFit", "sqrtAdjustment_combinedFit", function_sqrtAdjustment.c_str(), RooArgSet(rooVar_ST, rooVar_sqrt_combinedFit));
    RooProdPdf pdf_nJets_adjusted_combined(("pdf_2Jets_slopeAndSqrtAdjustment_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("pdf_2Jets_slopeAndSqrtAdjustment_" + std::to_string(nJetsBin) + "JetsBin").c_str(), RooArgSet(pdf_2Jets, pdf_slopeAdjustment_combinedFit, pdf_sqrtAdjustment_combinedFit));
    // double nll_combined;
    if (options.readParametersFromFiles) {
      rooVar_slope_combinedFit.setVal((fitParametersUnbinned.at("slope_sqrt_fit_slope")).at(nJetsBin));
      // rooVar_slope_combinedFit.setError((fitParametersUnbinned.at("slopeError_combinedFit")).at(nJetsBin));
      rooVar_sqrt_combinedFit.setVal((fitParametersUnbinned.at("slope_sqrt_fit_sqrt")).at(nJetsBin));
      // rooVar_sqrt_combinedFit.setError((fitParametersUnbinned.at("sqrtError_combinedFit")).at(nJetsBin));
      // nll_combined = (pdf_nJets_adjusted_combined.createNLL(*(STDataSets.at(nJetsBin)), Verbose(kFALSE), Range(options.STRegions.STNormRangeMin, ST_MAX_RANGE)))->getVal();
    }
    else {
      RooFitResult *fitResult = pdf_nJets_adjusted_combined.fitTo(*(STDataSets.at(nJetsBin)), Range(options.STRegions.STNormRangeMin, ST_MAX_RANGE), Optimize(kFALSE), PrintLevel(-1), Verbose(kFALSE), PrintEvalErrors(2), SumW2Error(kTRUE), Save(kTRUE));
      // check that the parameters are indexed as assumed
      RooArgList fit_parameters = fitResult->floatParsFinal();
      assert(std::string(fit_parameters.at(0)->GetName()) == slopeVar_combinedFit_name);
      assert(std::string(fit_parameters.at(1)->GetName()) == sqrtVar_combinedFit_name);
      assert(fit_parameters.getSize() == static_cast<int>(2));
      assert(fitResult->status() == 0);
      // nll_combined = fitResult->minNll();
      // step 1: get covariance matrix
      TMatrixDSym covarianceMatrix = fitResult->covarianceMatrix();
      std::cout << "For combined slope + sqrt fit, covarianceMatrix: ";
      printSquareMatrix(covarianceMatrix, fit_parameters.getSize());
      // step 2: get eigendecomposition
      TMatrixDSymEigen eigendecomposition_setup = TMatrixDSymEigen(covarianceMatrix);
      TVectorD eigenvalues = eigendecomposition_setup.GetEigenValues();
      std::cout << "eigenvalues: ";
      printTVector(eigenvalues);
      TMatrixD eigenvectors = eigendecomposition_setup.GetEigenVectors();
      std::cout << "eigenvectors: ";
      printSquareMatrix(eigenvectors, fit_parameters.getSize());
      std::vector<eigenmode_struct> eigenmodes;
      for (int eigen_index = 0; eigen_index < fit_parameters.getSize(); ++eigen_index) {
        double eigenvalue = eigenvalues(eigen_index);
        std::vector<double> eigenvector = getColumnFromTMatrixD(eigenvectors, eigen_index, fit_parameters.getSize());
        eigenmode_struct current_pair = eigenmode_struct(eigenvalue, eigenvector);
        check_eigendecomposition(current_pair, covarianceMatrix);
        eigenmodes.push_back(current_pair);
      }
      fitParametersUnbinned["slope_sqrt_fit_slope"][nJetsBin] = rooVar_slope_combinedFit.getValV();
      fitParametersUnbinnedList.push_back(std::string("float slope_sqrt_fit_slope_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("slope_sqrt_fit_slope")).at(nJetsBin))));
      // fitParametersUnbinned["slopeError_combinedFit"][nJetsBin] = rooVar_slope_combinedFit.getError();
      // fitParametersUnbinnedList.push_back(std::string("float slopeError_combinedFit_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("slopeError_combinedFit")).at(nJetsBin))));
      // std::cout << "best-fit slope: " << (fitParametersUnbinned.at("slope_sqrt_fit_slope")).at(nJetsBin) << " +/- " << (fitParametersUnbinned.at("slopeError_combinedFit")).at(nJetsBin) << std::endl;
      fitParametersUnbinned["slope_sqrt_fit_sqrt"][nJetsBin] = rooVar_sqrt_combinedFit.getValV();
      fitParametersUnbinnedList.push_back(std::string("float slope_sqrt_fit_sqrt_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("slope_sqrt_fit_sqrt")).at(nJetsBin))));
      // fitParametersUnbinned["sqrtError_combinedFit"][nJetsBin] = rooVar_sqrt_combinedFit.getError();
      // fitParametersUnbinnedList.push_back(std::string("float sqrtError_combinedFit_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("sqrtError_combinedFit")).at(nJetsBin))));
      // std::cout << "best-fit sqrt: " << (fitParametersUnbinned.at("slope_sqrt_fit_sqrt")).at(nJetsBin) << " +/- " << (fitParametersUnbinned.at("sqrtError_combinedFit")).at(nJetsBin) << std::endl;
      fitParametersUnbinned["slope_sqrt_fit_mode1_slopeCoefficient"][nJetsBin] = (eigenmodes.at(0)).eigenvector.at(0);
      fitParametersUnbinnedList.push_back(std::string("float slope_sqrt_fit_mode1_slopeCoefficient_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("slope_sqrt_fit_mode1_slopeCoefficient")).at(nJetsBin))));
      fitParametersUnbinned["slope_sqrt_fit_mode1_sqrtCoefficient"][nJetsBin] = (eigenmodes.at(0)).eigenvector.at(1);
      fitParametersUnbinnedList.push_back(std::string("float slope_sqrt_fit_mode1_sqrtCoefficient_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("slope_sqrt_fit_mode1_sqrtCoefficient")).at(nJetsBin))));
      fitParametersUnbinned["slope_sqrt_fit_mode1_error"][nJetsBin] = std::sqrt((eigenmodes.at(0)).eigenvalue);
      fitParametersUnbinnedList.push_back(std::string("float slope_sqrt_fit_mode1_error_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("slope_sqrt_fit_mode1_error")).at(nJetsBin))));
      fitParametersUnbinned["slope_sqrt_fit_mode2_slopeCoefficient"][nJetsBin] = (eigenmodes.at(1)).eigenvector.at(0);
      fitParametersUnbinnedList.push_back(std::string("float slope_sqrt_fit_mode2_slopeCoefficient_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("slope_sqrt_fit_mode2_slopeCoefficient")).at(nJetsBin))));
      fitParametersUnbinned["slope_sqrt_fit_mode2_sqrtCoefficient"][nJetsBin] = (eigenmodes.at(1)).eigenvector.at(1);
      fitParametersUnbinnedList.push_back(std::string("float slope_sqrt_fit_mode2_sqrtCoefficient_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("slope_sqrt_fit_mode2_sqrtCoefficient")).at(nJetsBin))));
      fitParametersUnbinned["slope_sqrt_fit_mode2_error"][nJetsBin] = std::sqrt((eigenmodes.at(1)).eigenvalue);
      fitParametersUnbinnedList.push_back(std::string("float slope_sqrt_fit_mode2_error_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("slope_sqrt_fit_mode2_error")).at(nJetsBin))));
    }
    rooVar_slope_combinedFit.Print();
    rooVar_sqrt_combinedFit.Print();
    // std::map<int, double> integralsOverBins_combined;
    // for (int dataHistBinIndex = 1; dataHistBinIndex <= dataHist.GetXaxis()->GetNbins(); ++dataHistBinIndex) {
    //   rooVar_ST.setRange((options.STRegions.STNormRangeMin - options.preNormalizationBuffer), ST_MAX_RANGE);
    //   integralsOverBins_combined[dataHistBinIndex] = (pdf_nJets_adjusted_combined.createIntegral(rooVar_ST, NormSet(rooVar_ST), Range(("range_STBinIndex_" + std::to_string(dataHistBinIndex)).c_str())))->getVal();
    // }
    // totalIntegral = (pdf_nJets_adjusted_combined.createIntegral(rooVar_ST, NormSet(rooVar_ST), Range("plotRange")))->getVal();
    // assert(std::fabs(totalIntegral - 1.0) < CHECK_TOLERANCE);
    pdf_nJets_adjusted_combined.plotOn(rooFrame, NormRange("normRange"), LineColor(static_cast<EColor>(kViolet)), LineWidth(1));
    // plus and minus one-sigma plotted with dashed linestyle
    rooVar_slope_combinedFit.setVal((fitParametersUnbinned.at("slope_sqrt_fit_slope")).at(nJetsBin) + ((fitParametersUnbinned.at("slope_sqrt_fit_mode1_error")).at(nJetsBin))*((fitParametersUnbinned.at("slope_sqrt_fit_mode1_slopeCoefficient")).at(nJetsBin)));
    rooVar_sqrt_combinedFit.setVal((fitParametersUnbinned.at("slope_sqrt_fit_sqrt")).at(nJetsBin) + ((fitParametersUnbinned.at("slope_sqrt_fit_mode1_error")).at(nJetsBin))*((fitParametersUnbinned.at("slope_sqrt_fit_mode1_sqrtCoefficient")).at(nJetsBin)));
    pdf_nJets_adjusted_combined.plotOn(rooFrame, NormRange("normRange"), LineColor(static_cast<EColor>(kViolet)), LineStyle(kDashed), LineWidth(1));
    rooVar_slope_combinedFit.setVal((fitParametersUnbinned.at("slope_sqrt_fit_slope")).at(nJetsBin) - ((fitParametersUnbinned.at("slope_sqrt_fit_mode1_error")).at(nJetsBin))*((fitParametersUnbinned.at("slope_sqrt_fit_mode1_slopeCoefficient")).at(nJetsBin)));
    rooVar_sqrt_combinedFit.setVal((fitParametersUnbinned.at("slope_sqrt_fit_sqrt")).at(nJetsBin) - ((fitParametersUnbinned.at("slope_sqrt_fit_mode1_error")).at(nJetsBin))*((fitParametersUnbinned.at("slope_sqrt_fit_mode1_sqrtCoefficient")).at(nJetsBin)));
    pdf_nJets_adjusted_combined.plotOn(rooFrame, NormRange("normRange"), LineColor(static_cast<EColor>(kViolet)), LineStyle(kDashed), LineWidth(1));
    rooVar_slope_combinedFit.setVal((fitParametersUnbinned.at("slope_sqrt_fit_slope")).at(nJetsBin));
    rooVar_sqrt_combinedFit.setVal((fitParametersUnbinned.at("slope_sqrt_fit_sqrt")).at(nJetsBin));
    TLegendEntry *legendEntry_unbinned_adjusted_combined = legend_dataSetsAndPdf_unbinned.AddEntry(&pdf_nJets_adjusted_combined, "2 jets kernel + slope adjustment + sqrt adjustment");
    legendEntry_unbinned_adjusted_combined->SetMarkerColor(static_cast<EColor>(kViolet));
    legendEntry_unbinned_adjusted_combined->SetLineColor(static_cast<EColor>(kViolet));
    legendEntry_unbinned_adjusted_combined->SetTextColor(static_cast<EColor>(kViolet));
    TGraph ratioGraph_unbinned_combinedToUnadjusted = TGraph();
    TGraph ratioGraph_unbinned_combinedToUnadjusted_high_estimate = TGraph();
    TGraph ratioGraph_unbinned_combinedToUnadjusted_low_estimate = TGraph();
    ratioGraph_unbinned_combinedToUnadjusted.SetName(("ratioGraph_unbinned_combinedToUnadjusted_at" + std::to_string(nJetsBin) + "Jets").c_str());
    ratioGraph_unbinned_combinedToUnadjusted.SetTitle(("combined fit at " + std::to_string(nJetsBin) + " Jets / unadjusted").c_str());
    rooVar_ST.setVal(options.STNormTarget);
    double pdf_nJets_adjusted_combined_at_STNorm = pdf_nJets_adjusted_combined.getVal(rooVar_ST);
    for (int STCounter = 0; STCounter <= 1000; ++STCounter) {
      double STVal = options.STRegions.STNormRangeMin + (1.0*STCounter/1000)*(ST_MAX_RANGE - options.STRegions.STNormRangeMin);
      rooVar_ST.setVal(STVal);
      rooVar_slope_combinedFit.setVal((fitParametersUnbinned.at("slope_sqrt_fit_slope")).at(nJetsBin));
      rooVar_sqrt_combinedFit.setVal((fitParametersUnbinned.at("slope_sqrt_fit_sqrt")).at(nJetsBin));
      double pdf_nominal = (pdf_nJets_adjusted_combined.getVal(rooVar_ST)/(pdf_nJets_adjusted_combined_at_STNorm));
      rooVar_slope_combinedFit.setVal((fitParametersUnbinned.at("slope_sqrt_fit_slope")).at(nJetsBin) + ((fitParametersUnbinned.at("slope_sqrt_fit_mode1_error")).at(nJetsBin))*((fitParametersUnbinned.at("slope_sqrt_fit_mode1_slopeCoefficient")).at(nJetsBin)));
      rooVar_sqrt_combinedFit.setVal((fitParametersUnbinned.at("slope_sqrt_fit_sqrt")).at(nJetsBin) + ((fitParametersUnbinned.at("slope_sqrt_fit_mode1_error")).at(nJetsBin))*((fitParametersUnbinned.at("slope_sqrt_fit_mode1_sqrtCoefficient")).at(nJetsBin)));
      double pdf_plus_one_sigma = (pdf_nJets_adjusted_combined.getVal(rooVar_ST)/(pdf_nJets_adjusted_combined_at_STNorm));
      rooVar_slope_combinedFit.setVal((fitParametersUnbinned.at("slope_sqrt_fit_slope")).at(nJetsBin) - ((fitParametersUnbinned.at("slope_sqrt_fit_mode1_error")).at(nJetsBin))*((fitParametersUnbinned.at("slope_sqrt_fit_mode1_slopeCoefficient")).at(nJetsBin)));
      rooVar_sqrt_combinedFit.setVal((fitParametersUnbinned.at("slope_sqrt_fit_sqrt")).at(nJetsBin) - ((fitParametersUnbinned.at("slope_sqrt_fit_mode1_error")).at(nJetsBin))*((fitParametersUnbinned.at("slope_sqrt_fit_mode1_sqrtCoefficient")).at(nJetsBin)));
      double pdf_minus_one_sigma = (pdf_nJets_adjusted_combined.getVal(rooVar_ST)/(pdf_nJets_adjusted_combined_at_STNorm));
      double pdf_higher = std::max({pdf_plus_one_sigma, pdf_nominal, pdf_minus_one_sigma});
      double pdf_lower = std::min({pdf_plus_one_sigma, pdf_nominal, pdf_minus_one_sigma});
      double commonDenominator = ((pdf_2Jets.getVal(rooVar_ST))/(pdf_2Jets_at_STNorm));
      double ratio_combinedToUnadjusted = pdf_nominal/commonDenominator;
      double ratio_higher = pdf_higher/commonDenominator;
      double ratio_lower = pdf_lower/commonDenominator;
      ratioGraph_unbinned_combinedToUnadjusted.SetPoint(STCounter, STVal, ratio_combinedToUnadjusted);
      ratioGraph_unbinned_combinedToUnadjusted_high_estimate.SetPoint(STCounter, STVal, ratio_higher);
      ratioGraph_unbinned_combinedToUnadjusted_low_estimate.SetPoint(STCounter, STVal, ratio_lower);
    }

    TCanvas pdfCanvas_unbinned = TCanvas(("c_dataSetAndPdf_unbinned_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("c_dataSetAndPdf_unbinned_" + std::to_string(nJetsBin) + "JetsBin").c_str(), 1600, 1280);
    rooFrame->Draw();
    pdfCanvas_unbinned.Update();
    rooFrame->SetMinimum((rooFrame->GetMaximum())/10000.);
    gPad->SetLogy();
    pdfCanvas_unbinned.Update();
    legend_dataSetsAndPdf_unbinned.SetFillStyle(0);
    legend_dataSetsAndPdf_unbinned.Draw();
    pdfCanvas_unbinned.Update();
    pdfCanvas_unbinned.SaveAs((options.outputFolder + "/unbinned_pdfAndData_" + std::to_string(nJetsBin) + "JetsBin_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".pdf").c_str());

    TMultiGraph unbinned_shape_ratios_multigraph = TMultiGraph(("unbinned_shape_ratios_multigraph_at" + std::to_string(nJetsBin) + "Jets").c_str(), ("Shape ratios, " + std::to_string(nJetsBin) + " Jets bin").c_str());
    TCanvas unbinned_shape_ratios_canvas = TCanvas(("c_unbinnedShapeRatios_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("c_unbinnedShapeRatios_" + std::to_string(nJetsBin) + "JetsBin").c_str(), 1600, 1280);
    TLegend legend_unbinned_shape_ratios_multigraph = TLegend(0.1, 0.6, 0.4, 0.9);
    ratioGraph_unbinned_slopeOnlyToUnadjusted.SetLineColor(static_cast<EColor>(kRed+1)); unbinned_shape_ratios_multigraph.Add(&ratioGraph_unbinned_slopeOnlyToUnadjusted); ratioGraph_unbinned_slopeOnlyToUnadjusted.SetDrawOption("C");
    ratioGraph_unbinned_slopeOnlyToUnadjusted_high_estimate.SetLineColor(static_cast<EColor>(kRed+1)); ratioGraph_unbinned_slopeOnlyToUnadjusted_high_estimate.SetLineStyle(kDashed); unbinned_shape_ratios_multigraph.Add(&ratioGraph_unbinned_slopeOnlyToUnadjusted_high_estimate); ratioGraph_unbinned_slopeOnlyToUnadjusted_high_estimate.SetDrawOption("C");
    ratioGraph_unbinned_slopeOnlyToUnadjusted_low_estimate.SetLineColor(static_cast<EColor>(kRed+1)); ratioGraph_unbinned_slopeOnlyToUnadjusted_low_estimate.SetLineStyle(kDashed); unbinned_shape_ratios_multigraph.Add(&ratioGraph_unbinned_slopeOnlyToUnadjusted_low_estimate); ratioGraph_unbinned_slopeOnlyToUnadjusted_low_estimate.SetDrawOption("C");
    TLegendEntry *legendEntry_unbinned_slopeOnlyToUnadjusted = legend_unbinned_shape_ratios_multigraph.AddEntry(&ratioGraph_unbinned_slopeOnlyToUnadjusted, "slope-only adjustment / 2 jets kernel");
    legendEntry_unbinned_slopeOnlyToUnadjusted->SetMarkerColor(static_cast<EColor>(kRed+1)); legendEntry_unbinned_slopeOnlyToUnadjusted->SetLineColor(static_cast<EColor>(kRed+1)); legendEntry_unbinned_slopeOnlyToUnadjusted->SetTextColor(static_cast<EColor>(kRed+1));
    ratioGraph_unbinned_sqrtOnlyToUnadjusted.SetLineColor(static_cast<EColor>(kGreen+3)); unbinned_shape_ratios_multigraph.Add(&ratioGraph_unbinned_sqrtOnlyToUnadjusted); ratioGraph_unbinned_sqrtOnlyToUnadjusted.SetDrawOption("C");
    ratioGraph_unbinned_sqrtOnlyToUnadjusted_high_estimate.SetLineColor(static_cast<EColor>(kGreen+3)); ratioGraph_unbinned_sqrtOnlyToUnadjusted_high_estimate.SetLineStyle(kDashed); unbinned_shape_ratios_multigraph.Add(&ratioGraph_unbinned_sqrtOnlyToUnadjusted_high_estimate); ratioGraph_unbinned_sqrtOnlyToUnadjusted_high_estimate.SetDrawOption("C");
    ratioGraph_unbinned_sqrtOnlyToUnadjusted_low_estimate.SetLineColor(static_cast<EColor>(kGreen+3)); ratioGraph_unbinned_sqrtOnlyToUnadjusted_low_estimate.SetLineStyle(kDashed); unbinned_shape_ratios_multigraph.Add(&ratioGraph_unbinned_sqrtOnlyToUnadjusted_low_estimate); ratioGraph_unbinned_sqrtOnlyToUnadjusted_low_estimate.SetDrawOption("C");
    TLegendEntry *legendEntry_unbinned_sqrtOnlyToUnadjusted = legend_unbinned_shape_ratios_multigraph.AddEntry(&ratioGraph_unbinned_sqrtOnlyToUnadjusted, "sqrt-only adjustment / 2 jets kernel");
    legendEntry_unbinned_sqrtOnlyToUnadjusted->SetMarkerColor(static_cast<EColor>(kGreen+3)); legendEntry_unbinned_sqrtOnlyToUnadjusted->SetLineColor(static_cast<EColor>(kGreen+3)); legendEntry_unbinned_sqrtOnlyToUnadjusted->SetTextColor(static_cast<EColor>(kGreen+3));
    ratioGraph_unbinned_combinedToUnadjusted.SetLineColor(static_cast<EColor>(kViolet)); unbinned_shape_ratios_multigraph.Add(&ratioGraph_unbinned_combinedToUnadjusted); ratioGraph_unbinned_combinedToUnadjusted.SetDrawOption("C");
    ratioGraph_unbinned_combinedToUnadjusted_high_estimate.SetLineColor(static_cast<EColor>(kViolet)); ratioGraph_unbinned_combinedToUnadjusted_high_estimate.SetLineStyle(kDashed); unbinned_shape_ratios_multigraph.Add(&ratioGraph_unbinned_combinedToUnadjusted_high_estimate); ratioGraph_unbinned_combinedToUnadjusted_high_estimate.SetDrawOption("C");
    ratioGraph_unbinned_combinedToUnadjusted_low_estimate.SetLineColor(static_cast<EColor>(kViolet)); ratioGraph_unbinned_combinedToUnadjusted_low_estimate.SetLineStyle(kDashed); unbinned_shape_ratios_multigraph.Add(&ratioGraph_unbinned_combinedToUnadjusted_low_estimate); ratioGraph_unbinned_combinedToUnadjusted_low_estimate.SetDrawOption("C");
    TLegendEntry *legendEntry_unbinned_combinedToUnadjusted = legend_unbinned_shape_ratios_multigraph.AddEntry(&ratioGraph_unbinned_combinedToUnadjusted, "combined adjustment / 2 jets kernel");
    legendEntry_unbinned_combinedToUnadjusted->SetMarkerColor(static_cast<EColor>(kViolet)); legendEntry_unbinned_combinedToUnadjusted->SetLineColor(static_cast<EColor>(kViolet)); legendEntry_unbinned_combinedToUnadjusted->SetTextColor(static_cast<EColor>(kViolet));
    ratioGraph_unbinned_nJetsKernelToUnadjusted.SetLineColor(static_cast<EColor>(kBlack)); unbinned_shape_ratios_multigraph.Add(&ratioGraph_unbinned_nJetsKernelToUnadjusted); ratioGraph_unbinned_nJetsKernelToUnadjusted.SetDrawOption("C");
    TLegendEntry *legendEntry_unbinned_nJetsKernelToUnadjusted = legend_unbinned_shape_ratios_multigraph.AddEntry(&ratioGraph_unbinned_nJetsKernelToUnadjusted, (std::to_string(nJetsBin) + " jets kernel / 2 jets kernel").c_str());
    legendEntry_unbinned_nJetsKernelToUnadjusted->SetMarkerColor(static_cast<EColor>(kBlack)); legendEntry_unbinned_nJetsKernelToUnadjusted->SetLineColor(static_cast<EColor>(kBlack)); legendEntry_unbinned_nJetsKernelToUnadjusted->SetTextColor(static_cast<EColor>(kBlack));
    unbinned_shape_ratios_multigraph.Draw("A");
    legend_unbinned_shape_ratios_multigraph.SetFillStyle(0);
    legend_unbinned_shape_ratios_multigraph.Draw();
    unbinned_shape_ratios_multigraph.GetXaxis()->SetTitle("ST (GeV)");
    unbinned_shape_ratios_multigraph.GetYaxis()->SetTitle("ratio");
    unbinned_shape_ratios_canvas.SaveAs((options.outputFolder + "/unbinned_shapeRatios_" + std::to_string(nJetsBin) + "JetsBin_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".pdf").c_str());

    std::cout << "Now the binned analysis: " << std::endl;

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
      customizedPDF *customized_pdf = new customizedPDF(&pdf_2Jets, &rooVar_ST, options.STNormTarget, customization_type);
      customized_pdf->setNominalScale("normRange", ((STHistograms.at(nJetsBin)).GetBinContent(1))*((STHistograms.at(nJetsBin)).GetBinWidth(1)));
      customizedTF1 *customized_tf1 = new customizedTF1(std::string("pdf_2Jets_"), customized_pdf, options.STRegions.STNormRangeMin, ST_MAX_RANGE, customization_type);
      customized_tf1->initializeParameters(parameter_initializations.at(customization_type));
      if (options.readParametersFromFiles) {
        customized_tf1->setFitResultsFromSource(fitParametersBinned, nJetsBin);
      }
      else {
        customized_tf1->fitToTH1(STHistograms.at(nJetsBin));
        fit_pvalues[customization_type][nJetsBin] = (customized_tf1->fit_result).pvalue;
        for (int parameter_index = 0; parameter_index < customizationTypeNPars.at(customization_type); ++parameter_index) {
          fitParametersBinned[get_parameter_name(customization_type, parameter_index, nJetsBin)] = ((customized_tf1->fit_result).best_fit_values).at(parameter_index);
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

    // data distribution
    TGraphErrors ratioGraph_binned_nJetsDistribution_to_unadjusted = TGraphErrors();
    ratioGraph_binned_nJetsDistribution_to_unadjusted.SetName(("ratioGraph_binned_nJetsDistribution_to_unadjusted_at" + std::to_string(nJetsBin) + "Jets").c_str());
    ratioGraph_binned_nJetsDistribution_to_unadjusted.SetTitle(("ST distribution at " + std::to_string(nJetsBin) + " Jets / unadjusted").c_str());
    std::map<int, double> bin_integrals_divided_by_bin_widths_from_2_jets_kernel;
    if (!(options.readParametersFromFiles)) {
      (customized_tf1s.at(customization_type_denominator_for_ratios))->set_TF_parameters_to_nominal();
      bin_integrals_divided_by_bin_widths_from_2_jets_kernel = (customized_tf1s.at(customization_type_denominator_for_ratios))->getBinIntegralsDividedByBinWidthFromTF1(options.STRegions);
    }
    for (int binCounter = 1; binCounter <= (STHistograms.at(nJetsBin)).GetXaxis()->GetNbins(); ++binCounter) {
      double STMidpoint = (STHistograms.at(nJetsBin)).GetXaxis()->GetBinCenter(binCounter);
      double binWidth = (STHistograms.at(nJetsBin)).GetXaxis()->GetBinWidth(binCounter);
      double numerator = (STHistograms.at(nJetsBin)).GetBinContent(binCounter);
      (customized_tf1s.at(customization_type_denominator_for_ratios))->set_TF_parameters_to_nominal();
      double denominator = ((customized_tf1s.at(customization_type_denominator_for_ratios))->getTFIntegral((STHistograms.at(nJetsBin)).GetXaxis()->GetBinLowEdge(binCounter), (STHistograms.at(nJetsBin)).GetXaxis()->GetBinUpEdge(binCounter)))/((STHistograms.at(nJetsBin)).GetXaxis()->GetBinWidth(binCounter));
      assert(denominator > 0.);
      double ratio = numerator/denominator;
      double numeratorError = (STHistograms.at(nJetsBin)).GetBinError(binCounter);
      double denominatorError = 0.; // might change later
      double ratioError = ratio*std::sqrt(std::pow(numeratorError/numerator, 2) + std::pow(denominatorError/denominator, 2));
      int graph_currentPointIndex = ratioGraph_binned_nJetsDistribution_to_unadjusted.GetN();
      ratioGraph_binned_nJetsDistribution_to_unadjusted.SetPoint(graph_currentPointIndex, STMidpoint, ratio);
      ratioGraph_binned_nJetsDistribution_to_unadjusted.SetPointError(graph_currentPointIndex, binWidth/(std::sqrt(12)), ratioError);
    }

    // initialize some variables useful for plots
    double fractionalError_normBin = ((STHistograms.at(nJetsBin)).GetBinError(1))/((STHistograms.at(nJetsBin)).GetBinContent(1));
    assert (fractionalError_normBin < 1.0); // sanity check, to make sure weights aren't affecting the errors in weird ways...

    // plot the raw shapes
    TCanvas pdfCanvas_binned = TCanvas(("c_dataSetAndPdf_binned_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("c_dataSetAndPdf_binned_" + std::to_string(nJetsBin) + "JetsBin").c_str(), 1600, 1280);
    TLegend legend_dataSetsAndPdf_binned = TLegend(0.5, 0.6, 0.9, 0.9);
    gStyle->SetOptStat(0);

    (STHistograms.at(nJetsBin)).SetLineColor(static_cast<EColor>(kBlack)); (STHistograms.at(nJetsBin)).Draw(); pdfCanvas_binned.Update();
    TLegendEntry *legendEntry_binned_nJetsDistribution = legend_dataSetsAndPdf_binned.AddEntry(&(STHistograms.at(nJetsBin)), (std::to_string(nJetsBin) + " jets distribution").c_str());
    legendEntry_binned_nJetsDistribution->SetMarkerColor(static_cast<EColor>(kBlack)); legendEntry_binned_nJetsDistribution->SetLineColor(static_cast<EColor>(kBlack)); legendEntry_binned_nJetsDistribution->SetTextColor(static_cast<EColor>(kBlack));
    (STHistograms.at(nJetsBin)).GetYaxis()->SetRange(((STHistograms.at(nJetsBin)).GetMaximum())/10000., ((STHistograms.at(nJetsBin)).GetMaximum())); pdfCanvas_binned.Update();

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
        std::map<int, double> adjustments_nominal = getNominalAdjustmentsFromBinIntegralMaps(options.STRegions, bin_integrals_divided_by_bin_widths_nominal, bin_integrals_divided_by_bin_widths_from_2_jets_kernel);
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
    (STHistograms.at(nJetsBin)).Draw("SAME"); pdfCanvas_binned.Update(); // draw the data again so the datapoints aren't obscured by later plots
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

    ratioGraph_binned_nJetsDistribution_to_unadjusted.SetLineColor(static_cast<EColor>(kBlack)); ratioGraph_binned_nJetsDistribution_to_unadjusted.SetDrawOption("P"); binned_shape_ratios_multigraph.Add(&ratioGraph_binned_nJetsDistribution_to_unadjusted);
    TLegendEntry *legendEntry_binned_nJetsDistribution_to_unadjusted = legend_binned_shape_ratios_multigraph.AddEntry(&ratioGraph_binned_nJetsDistribution_to_unadjusted, (std::to_string(nJetsBin) + " jets distribution / 2 jets kernel").c_str());
    legendEntry_binned_nJetsDistribution_to_unadjusted->SetMarkerColor(static_cast<EColor>(kBlack)); legendEntry_binned_nJetsDistribution_to_unadjusted->SetLineColor(static_cast<EColor>(kBlack)); legendEntry_binned_nJetsDistribution_to_unadjusted->SetTextColor(static_cast<EColor>(kBlack));

    binned_shape_ratios_multigraph.Draw("A");
    legend_binned_shape_ratios_multigraph.SetFillStyle(0);
    legend_binned_shape_ratios_multigraph.Draw();
    binned_shape_ratios_multigraph.GetXaxis()->SetTitle("ST (GeV)");
    binned_shape_ratios_multigraph.GetYaxis()->SetTitle("ratio");
    binned_shape_ratios_multigraph.SetMinimum(-0.5);
    binned_shape_ratios_multigraph.SetMaximum(5.5);
    binned_shape_ratios_canvas.Update();
    binned_shape_ratios_canvas.SaveAs((options.outputFolder + "/binned_shapeRatios_" + std::to_string(nJetsBin) + "JetsBin_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".pdf").c_str());

    // if doing a comparison (e.g. reading parameters from MC and applying to data), calculate and save the ratios of the data wrt the nominal adjustment and fits to a straight line
    if (options.readParametersFromFiles) {
      TGraphErrors ratios_wrt_chosen_adjustment = TGraphErrors();
      for (int binCounter = 1; binCounter <= (STHistograms.at(nJetsBin)).GetXaxis()->GetNbins(); ++binCounter) {
        double STMidpoint = (STHistograms.at(nJetsBin)).GetXaxis()->GetBinCenter(binCounter);
        double binWidth = (STHistograms.at(nJetsBin)).GetXaxis()->GetBinWidth(binCounter);
        double numerator = (STHistograms.at(nJetsBin)).GetBinContent(binCounter);
        (customized_tf1s.at(customization_type_for_adjustments_output))->set_TF_parameters_to_nominal();
        double denominator = ((customized_tf1s.at(customization_type_for_adjustments_output))->getTFIntegral((STHistograms.at(nJetsBin)).GetXaxis()->GetBinLowEdge(binCounter), (STHistograms.at(nJetsBin)).GetXaxis()->GetBinUpEdge(binCounter)))/((STHistograms.at(nJetsBin)).GetXaxis()->GetBinWidth(binCounter));
        assert(denominator > 0.);
        double ratio = numerator/denominator;
        double numeratorError = (STHistograms.at(nJetsBin)).GetBinError(binCounter);
        double denominatorError = 0.; // might change later
        double ratioError = ratio*std::sqrt(std::pow(numeratorError/numerator, 2) + std::pow(denominatorError/denominator, 2));
        int graph_currentPointIndex = ratios_wrt_chosen_adjustment.GetN();
        ratios_wrt_chosen_adjustment.SetPoint(graph_currentPointIndex, STMidpoint, ratio);
        ratios_wrt_chosen_adjustment.SetPointError(graph_currentPointIndex, binWidth/(std::sqrt(12)), ratioError);
      }
      std::string straight_line_functional_form_for_TF1 = "[0] + [1]*((x/" + std::to_string(options.STNormTarget) + ") - 1.0)";
      TF1 fitFunction_ratios_wrt_chosen_adjustment(("ratios_wrt_chosen_adjustment_at" + std::to_string(nJetsBin) + "Jets").c_str(), straight_line_functional_form_for_TF1.c_str(), options.STRegions.STNormRangeMin, ST_MAX_RANGE);
      fitFunction_ratios_wrt_chosen_adjustment.SetParName(0, ("ratios_wrt_chosen_adjustment_const_" + std::to_string(nJetsBin) + "JetsBin").c_str());
      fitFunction_ratios_wrt_chosen_adjustment.SetParameter(0, 1.0);
      fitFunction_ratios_wrt_chosen_adjustment.SetParLimits(0, scale_minVal, scale_maxVal);
      fitFunction_ratios_wrt_chosen_adjustment.SetParName(1, ("ratios_wrt_chosen_adjustment_slope_" + std::to_string(nJetsBin) + "JetsBin").c_str());
      fitFunction_ratios_wrt_chosen_adjustment.SetParameter(1, 0.);
      fitFunction_ratios_wrt_chosen_adjustment.SetParLimits(1, slope_minVal, slope_maxVal);
      TFitResultPtr ratios_wrt_chosen_adjustment_fit_result = ratios_wrt_chosen_adjustment.Fit(&fitFunction_ratios_wrt_chosen_adjustment, (constants::binnedFitOptions_ratios_wrt_chosen_adjustment).c_str());
      assert(ratios_wrt_chosen_adjustment_fit_result->Status() == 0);
      double best_fit_const = ratios_wrt_chosen_adjustment_fit_result->Value(0);
      double best_fit_const_error = ratios_wrt_chosen_adjustment_fit_result->ParError(0);
      double best_fit_slope = ratios_wrt_chosen_adjustment_fit_result->Value(1);
      double best_fit_slope_error = ratios_wrt_chosen_adjustment_fit_result->ParError(1);
      TCanvas canvas_ratios_wrt_chosen_adjustment = TCanvas(("c_ratios_wrt_chosen_adjustment_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("c_ratios_wrt_chosen_adjustment_" + std::to_string(nJetsBin) + "JetsBin").c_str(), 1600, 1280);
      gStyle->SetOptStat(0);
      TLegend legend_ratios_wrt_chosen_adjustment = TLegend(0.1, 0.7, 0.9, 0.9);
      legend_ratios_wrt_chosen_adjustment.SetFillStyle(0);
      ratios_wrt_chosen_adjustment.Draw("AP0"); canvas_ratios_wrt_chosen_adjustment.Update();
      ratios_wrt_chosen_adjustment.GetYaxis()->SetRangeUser(-0.5, 3.5);
      ratios_wrt_chosen_adjustment.SetLineColor(static_cast<EColor>(kBlack)); ratios_wrt_chosen_adjustment.SetLineWidth(2);
      TLegendEntry *legendEntry_graph_ratios_wrt_chosen_adjustment = legend_ratios_wrt_chosen_adjustment.AddEntry(&ratios_wrt_chosen_adjustment, (std::to_string(nJetsBin) + " jets distribution / " + customizationTypeLegendLabels.at(customization_type_for_adjustments_output)).c_str());
      legendEntry_graph_ratios_wrt_chosen_adjustment->SetMarkerColor(static_cast<EColor>(kBlack)); legendEntry_graph_ratios_wrt_chosen_adjustment->SetLineColor(static_cast<EColor>(kBlack)); legendEntry_graph_ratios_wrt_chosen_adjustment->SetTextColor(static_cast<EColor>(kBlack));
      fitFunction_ratios_wrt_chosen_adjustment.Draw("C SAME"); canvas_ratios_wrt_chosen_adjustment.Update();
      fitFunction_ratios_wrt_chosen_adjustment.SetLineColor(static_cast<EColor>(kBlue)); fitFunction_ratios_wrt_chosen_adjustment.SetLineWidth(2);
      TLegendEntry *legendEntry_nominal_fit = legend_ratios_wrt_chosen_adjustment.AddEntry(&ratios_wrt_chosen_adjustment, ("nominal fit: (" + get_string_precision_n(4, best_fit_const) + " #pm " + get_string_precision_n(4, best_fit_const_error) + ") + (" + get_string_precision_n(4, best_fit_slope) + " #pm " + get_string_precision_n(4, best_fit_slope_error) + ")*(ST/" + get_string_precision_n(5, options.STNormTarget) + " - 1.0)").c_str());
      legendEntry_nominal_fit->SetMarkerColor(static_cast<EColor>(kBlue)); legendEntry_nominal_fit->SetLineColor(static_cast<EColor>(kBlue)); legendEntry_nominal_fit->SetTextColor(static_cast<EColor>(kBlue));
      legend_ratios_wrt_chosen_adjustment.Draw();
      canvas_ratios_wrt_chosen_adjustment.Update();
      canvas_ratios_wrt_chosen_adjustment.SaveAs((options.outputFolder + "/binned_ratios_wrt_chosen_adjustment_" + std::to_string(nJetsBin) + "JetsBin_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".pdf").c_str());
      std::map<int, double> adjustments_from_slope = get_adjustments_from_slope(best_fit_slope, options.STRegions_for_ratio_wrt_chosen_adjustment, (customized_tf1s.at(customization_type_for_adjustments_output))->raw_TF1);
      for (int regionIndex = 1; regionIndex <= options.STRegions_for_ratio_wrt_chosen_adjustment.STAxis.GetNbins(); ++regionIndex) {
        ratio_adjustments_forOutputFile.push_back("float ratio_adjustment_STRegion" + std::to_string(regionIndex) + "_" + std::to_string(nJetsBin) + "Jets=" + std::to_string(adjustments_from_slope.at(regionIndex)));
      }
    }

    if (!(options.readParametersFromFiles)) {
      std::cout << "Getting p-values using chi2 values from binned fits..." << std::endl;
      // sanity checks
      assert((((customized_tf1s.at(customizationType::ScaleOnly))->fit_result).ndf) == (1 + (((customized_tf1s.at(customizationType::Slope))->fit_result).ndf)));
      assert((((customized_tf1s.at(customizationType::ScaleOnly))->fit_result).ndf) == (1 + (((customized_tf1s.at(customizationType::Sqrt))->fit_result).ndf)));
      assert((((customized_tf1s.at(customizationType::ScaleOnly))->fit_result).ndf) == (2 + (((customized_tf1s.at(customizationType::SlopeSqrt))->fit_result).ndf)));
      assert((((customized_tf1s.at(customizationType::ScaleOnly))->fit_result).ndf) == (3 + (((customized_tf1s.at(customizationType::SlopeSqrtQuad))->fit_result).ndf)));

      ftest_pValues["unadjusted_vs_slope"][nJetsBin] = get_fTest_prob(((customized_tf1s.at(customizationType::ScaleOnly))->fit_result).chi_sq, ((customized_tf1s.at(customizationType::Slope))->fit_result).chi_sq, ((customized_tf1s.at(customizationType::ScaleOnly))->fit_result).ndf, ((customized_tf1s.at(customizationType::Slope))->fit_result).ndf);
      ftest_pValues["slope_vs_slope_sqrt"][nJetsBin] = get_fTest_prob(((customized_tf1s.at(customizationType::Slope))->fit_result).chi_sq, ((customized_tf1s.at(customizationType::SlopeSqrt))->fit_result).chi_sq, ((customized_tf1s.at(customizationType::Slope))->fit_result).ndf, ((customized_tf1s.at(customizationType::SlopeSqrt))->fit_result).ndf);
      ftest_pValues["unadjusted_vs_sqrt"][nJetsBin] = get_fTest_prob(((customized_tf1s.at(customizationType::ScaleOnly))->fit_result).chi_sq, ((customized_tf1s.at(customizationType::Sqrt))->fit_result).chi_sq, ((customized_tf1s.at(customizationType::ScaleOnly))->fit_result).ndf, ((customized_tf1s.at(customizationType::Sqrt))->fit_result).ndf);
      ftest_pValues["sqrt_vs_slope_sqrt"][nJetsBin] = get_fTest_prob(((customized_tf1s.at(customizationType::Sqrt))->fit_result).chi_sq, ((customized_tf1s.at(customizationType::SlopeSqrt))->fit_result).chi_sq, ((customized_tf1s.at(customizationType::Sqrt))->fit_result).ndf, ((customized_tf1s.at(customizationType::SlopeSqrt))->fit_result).ndf);
      ftest_pValues["slope_sqrt_vs_slope_sqrt_quad"][nJetsBin] = get_fTest_prob(((customized_tf1s.at(customizationType::SlopeSqrt))->fit_result).chi_sq, ((customized_tf1s.at(customizationType::SlopeSqrtQuad))->fit_result).chi_sq, ((customized_tf1s.at(customizationType::SlopeSqrt))->fit_result).ndf, ((customized_tf1s.at(customizationType::SlopeSqrtQuad))->fit_result).ndf);
    }

    for (int customization_type_index = customizationTypeFirst; customization_type_index < static_cast<int>(customizationType::nCustomizationTypes); ++customization_type_index) {
      customizationType customization_type = static_cast<customizationType>(customization_type_index);
      delete customized_tf1s.at(customization_type);
      delete customized_pdfs.at(customization_type);
    }

    printSeparator();
  }
  delete random_generator;

  if (options.readParametersFromFiles) {
    // write out adjustment values
    std::ofstream ratio_adjustment_outputFile((options.outputFolder + "/ratio_adjustment_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".dat").c_str());
    assert(ratio_adjustment_outputFile.is_open());
    for (int ratio_adjustment_index = 0; ratio_adjustment_index < static_cast<int>(ratio_adjustments_forOutputFile.size()); ++ratio_adjustment_index) {
      ratio_adjustment_outputFile << ratio_adjustments_forOutputFile.at(ratio_adjustment_index) << std::endl;
    }
    ratio_adjustment_outputFile.close();
    std::cout << "Ratio adjustments written to file: " << (options.outputFolder + "/ratio_adjustment_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".dat") << std::endl;
  }
  else {
    // print ftest pvalues from binned chi2 fits in a LaTeX-formatted table
    std::cout << "p-values for binned fit comparisons using f-statistic:" << std::endl;
    std::cout << "\\begin{tabular}{|p{0.14\\textwidth}|p{0.14\\textwidth}|p{0.14\\textwidth}|p{0.14\\textwidth}|p{0.14\\textwidth}|p{0.14\\textwidth}|}" << std::endl;
    std::cout << "  \\hline" << std::endl;
    std::cout << "  p-values \\newline (f-statistic) & unadjusted \\newline vs \\newline linear & linear \\newline vs \\newline (linear+sqrt) & unadjusted \\newline vs \\newline sqrt & sqrt \\newline vs \\newline (linear+sqrt) & (linear+sqrt) \\newline vs \\newline (linear+sqrt \\newline +quad) \\\\ \\hline" << std::endl;
    std::cout << std::fixed << std::setprecision(3) << "  nJets = 3 & " << (ftest_pValues.at("unadjusted_vs_slope")).at(3) << " & " << (ftest_pValues.at("slope_vs_slope_sqrt")).at(3) << " & " << (ftest_pValues.at("unadjusted_vs_sqrt")).at(3) << " & " << (ftest_pValues.at("sqrt_vs_slope_sqrt")).at(3) << " & " << (ftest_pValues.at("slope_sqrt_vs_slope_sqrt_quad")).at(3) << " \\\\ \\hline" << std::endl;
    std::cout << std::fixed << std::setprecision(3) << "  nJets = 4 & " << (ftest_pValues.at("unadjusted_vs_slope")).at(4) << " & " << (ftest_pValues.at("slope_vs_slope_sqrt")).at(4) << " & " << (ftest_pValues.at("unadjusted_vs_sqrt")).at(4) << " & " << (ftest_pValues.at("sqrt_vs_slope_sqrt")).at(4) << " & " << (ftest_pValues.at("slope_sqrt_vs_slope_sqrt_quad")).at(4) << " \\\\ \\hline" << std::endl;
    std::cout << std::fixed << std::setprecision(3) << "  nJets = 5 & " << (ftest_pValues.at("unadjusted_vs_slope")).at(5) << " & " << (ftest_pValues.at("slope_vs_slope_sqrt")).at(5) << " & " << (ftest_pValues.at("unadjusted_vs_sqrt")).at(5) << " & " << (ftest_pValues.at("sqrt_vs_slope_sqrt")).at(5) << " & " << (ftest_pValues.at("slope_sqrt_vs_slope_sqrt_quad")).at(5) << " \\\\ \\hline" << std::endl;
    std::cout << std::fixed << std::setprecision(3) << "  nJets $\\geq$ 6 & " << (ftest_pValues.at("unadjusted_vs_slope")).at(6) << " & " << (ftest_pValues.at("slope_vs_slope_sqrt")).at(6) << " & " << (ftest_pValues.at("unadjusted_vs_sqrt")).at(6) << " & " << (ftest_pValues.at("sqrt_vs_slope_sqrt")).at(6) << " & " << (ftest_pValues.at("slope_sqrt_vs_slope_sqrt_quad")).at(6) << " \\\\ \\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;

    // print p-values for best fits in a LaTeX-formatted table
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
    for (int nJetsBin = 3; nJetsBin <= 6; ++nJetsBin) {
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

    // following commented out because we no longer use the unbinned fits
    // // Print best-fit values for combined fit in a LaTeX-formatted table
    // // example:
    // std::cout << "Best fit values for unbinned combined fit:" << std::endl;
    // std::cout << "\\begin{tabular}{|p{0.14\\textwidth}|p{0.1\\textwidth}|p{0.1\\textwidth}|p{0.25\\textwidth}|p{0.25\\textwidth}|}" << std::endl;
    // std::cout << "  \\hline" << std::endl;
    // std::cout << "  best-fits & $m$ & $p$ & $\\sqrt{\\lambda_1}$; eigenmode 1 & $\\sqrt{\\lambda_2}$; eigenmode 2 \\\\ \\hline" << std::endl;
    // std::cout << std::setprecision(3) << "  nJets = 3 & " << (fitParametersUnbinned.at("slope_sqrt_fit_slope")).at(3) << " & " << (fitParametersUnbinned.at("slope_sqrt_fit_sqrt")).at(3) << " & " << (fitParametersUnbinned.at("slope_sqrt_fit_mode1_error")).at(3) << "; (" << std::setprecision(2) << (fitParametersUnbinned.at("slope_sqrt_fit_mode1_slopeCoefficient")).at(3) << ", " << (fitParametersUnbinned.at("slope_sqrt_fit_mode1_sqrtCoefficient")).at(3) << std::setprecision(3) << ") & " << (fitParametersUnbinned.at("slope_sqrt_fit_mode2_error")).at(3) << std::setprecision(2) << "; (" << (fitParametersUnbinned.at("slope_sqrt_fit_mode2_slopeCoefficient")).at(3) << ", " << (fitParametersUnbinned.at("slope_sqrt_fit_mode2_sqrtCoefficient")).at(3) << ") \\\\ \\hline" << std::endl;
    // std::cout << std::setprecision(3) << "  nJets = 4 & " << (fitParametersUnbinned.at("slope_sqrt_fit_slope")).at(4) << " & " << (fitParametersUnbinned.at("slope_sqrt_fit_sqrt")).at(4) << " & " << (fitParametersUnbinned.at("slope_sqrt_fit_mode1_error")).at(4) << "; (" << std::setprecision(2) << (fitParametersUnbinned.at("slope_sqrt_fit_mode1_slopeCoefficient")).at(4) << ", " << (fitParametersUnbinned.at("slope_sqrt_fit_mode1_sqrtCoefficient")).at(4) << std::setprecision(3) << ") & " << (fitParametersUnbinned.at("slope_sqrt_fit_mode2_error")).at(4) << std::setprecision(2) << "; (" << (fitParametersUnbinned.at("slope_sqrt_fit_mode2_slopeCoefficient")).at(4) << ", " << (fitParametersUnbinned.at("slope_sqrt_fit_mode2_sqrtCoefficient")).at(4) << ") \\\\ \\hline" << std::endl;
    // std::cout << std::setprecision(3) << "  nJets = 5 & " << (fitParametersUnbinned.at("slope_sqrt_fit_slope")).at(5) << " & " << (fitParametersUnbinned.at("slope_sqrt_fit_sqrt")).at(5) << " & " << (fitParametersUnbinned.at("slope_sqrt_fit_mode1_error")).at(5) << "; (" << std::setprecision(2) << (fitParametersUnbinned.at("slope_sqrt_fit_mode1_slopeCoefficient")).at(5) << ", " << (fitParametersUnbinned.at("slope_sqrt_fit_mode1_sqrtCoefficient")).at(5) << std::setprecision(3) << ") & " << (fitParametersUnbinned.at("slope_sqrt_fit_mode2_error")).at(5) << std::setprecision(2) << "; (" << (fitParametersUnbinned.at("slope_sqrt_fit_mode2_slopeCoefficient")).at(5) << ", " << (fitParametersUnbinned.at("slope_sqrt_fit_mode2_sqrtCoefficient")).at(5) << ") \\\\ \\hline" << std::endl;
    // std::cout << std::setprecision(3) << "  nJets $\\geq$ 6 & " << (fitParametersUnbinned.at("slope_sqrt_fit_slope")).at(6) << " & " << (fitParametersUnbinned.at("slope_sqrt_fit_sqrt")).at(6) << " & " << (fitParametersUnbinned.at("slope_sqrt_fit_mode1_error")).at(6) << std::setprecision(2) << "; (" << (fitParametersUnbinned.at("slope_sqrt_fit_mode1_slopeCoefficient")).at(6) << ", " << (fitParametersUnbinned.at("slope_sqrt_fit_mode1_sqrtCoefficient")).at(6) << std::setprecision(3) << ") & " << (fitParametersUnbinned.at("slope_sqrt_fit_mode2_error")).at(6) << std::setprecision(2) << "; (" << (fitParametersUnbinned.at("slope_sqrt_fit_mode2_slopeCoefficient")).at(6) << ", " << (fitParametersUnbinned.at("slope_sqrt_fit_mode2_sqrtCoefficient")).at(6) << ") \\\\ \\hline" << std::endl;
    // std::cout << std::fixed << "\\end{tabular}" << std::endl;

    // following commented out because it was suitable for an older version of this script in which the slope+sqrt fit only had two free parameters
    // std::cout << "Best fit values for binned combined fit:" << std::endl;
    // std::cout << "\\begin{tabular}{|p{0.14\\textwidth}|p{0.1\\textwidth}|p{0.1\\textwidth}|p{0.25\\textwidth}|p{0.25\\textwidth}|}" << std::endl;
    // std::cout << "  \\hline" << std::endl;
    // std::cout << "  best-fits & $m$ & $p$ & $\\sqrt{\\lambda_1}$; eigenmode 1 & $\\sqrt{\\lambda_2}$; eigenmode 2 \\\\ \\hline" << std::endl;

    // std::cout << std::setprecision(3) << "  nJets = 3 & " << fitParametersBinned.at(get_parameter_name(customizationType::SlopeSqrt, 0, 3)) << " & " << fitParametersBinned.at(get_parameter_name(customizationType::SlopeSqrt, 1, 3)) << " & " << fitParametersBinned.at(get_eigenerror_name(customizationType::SlopeSqrt, 0, 3)) << "; (" << std::setprecision(2) << fitParametersBinned.at(get_eigencoefficient_name(customizationType::SlopeSqrt, 0, 0, 3)) << ", " << fitParametersBinned.at(get_eigencoefficient_name(customizationType::SlopeSqrt, 0, 1, 3)) << std::setprecision(3) << ") & " << fitParametersBinned.at(get_eigenerror_name(customizationType::SlopeSqrt, 1, 3)) << std::setprecision(2) << "; (" << fitParametersBinned.at(get_eigencoefficient_name(customizationType::SlopeSqrt, 1, 0, 3)) << ", " << fitParametersBinned.at(get_eigencoefficient_name(customizationType::SlopeSqrt, 1, 1, 3)) << ") \\\\ \\hline" << std::endl;
    // std::cout << std::setprecision(3) << "  nJets = 4 & " << fitParametersBinned.at(get_parameter_name(customizationType::SlopeSqrt, 0, 4)) << " & " << fitParametersBinned.at(get_parameter_name(customizationType::SlopeSqrt, 1, 4)) << " & " << fitParametersBinned.at(get_eigenerror_name(customizationType::SlopeSqrt, 0, 4)) << "; (" << std::setprecision(2) << fitParametersBinned.at(get_eigencoefficient_name(customizationType::SlopeSqrt, 0, 0, 4)) << ", " << fitParametersBinned.at(get_eigencoefficient_name(customizationType::SlopeSqrt, 0, 1, 4)) << std::setprecision(3) << ") & " << fitParametersBinned.at(get_eigenerror_name(customizationType::SlopeSqrt, 1, 4)) << std::setprecision(2) << "; (" << fitParametersBinned.at(get_eigencoefficient_name(customizationType::SlopeSqrt, 1, 0, 4)) << ", " << fitParametersBinned.at(get_eigencoefficient_name(customizationType::SlopeSqrt, 1, 1, 4)) << ") \\\\ \\hline" << std::endl;
    // std::cout << std::setprecision(3) << "  nJets = 5 & " << fitParametersBinned.at(get_parameter_name(customizationType::SlopeSqrt, 0, 5)) << " & " << fitParametersBinned.at(get_parameter_name(customizationType::SlopeSqrt, 1, 5)) << " & " << fitParametersBinned.at(get_eigenerror_name(customizationType::SlopeSqrt, 0, 5)) << "; (" << std::setprecision(2) << fitParametersBinned.at(get_eigencoefficient_name(customizationType::SlopeSqrt, 0, 0, 5)) << ", " << fitParametersBinned.at(get_eigencoefficient_name(customizationType::SlopeSqrt, 0, 1, 5)) << std::setprecision(3) << ") & " << fitParametersBinned.at(get_eigenerror_name(customizationType::SlopeSqrt, 1, 5)) << std::setprecision(2) << "; (" << fitParametersBinned.at(get_eigencoefficient_name(customizationType::SlopeSqrt, 1, 0, 5)) << ", " << fitParametersBinned.at(get_eigencoefficient_name(customizationType::SlopeSqrt, 1, 1, 5)) << ") \\\\ \\hline" << std::endl;
    // std::cout << std::setprecision(3) << "  nJets $\\geq$ 6 & " << fitParametersBinned.at(get_parameter_name(customizationType::SlopeSqrt, 0, 6)) << " & " << fitParametersBinned.at(get_parameter_name(customizationType::SlopeSqrt, 1, 6)) << " & " << fitParametersBinned.at(get_eigenerror_name(customizationType::SlopeSqrt, 0, 6)) << "; (" << std::setprecision(2) << fitParametersBinned.at(get_eigencoefficient_name(customizationType::SlopeSqrt, 0, 0, 6)) << ", " << fitParametersBinned.at(get_eigencoefficient_name(customizationType::SlopeSqrt, 0, 1, 6)) << std::setprecision(3) << ") & " << fitParametersBinned.at(get_eigenerror_name(customizationType::SlopeSqrt, 1, 6)) << std::setprecision(2) << "; (" << fitParametersBinned.at(get_eigencoefficient_name(customizationType::SlopeSqrt, 1, 0, 6)) << ", " << fitParametersBinned.at(get_eigencoefficient_name(customizationType::SlopeSqrt, 1, 1, 6)) << ") \\\\ \\hline" << std::endl;

    // write parameters for unbinned fit
    std::ofstream fitParametersUnbinnedFile((options.outputFolder + "/unbinned_fitParameters_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".dat").c_str());
    assert(fitParametersUnbinnedFile.is_open());
    for (int fitParametersUnbinnedListIndex = 0; fitParametersUnbinnedListIndex < static_cast<int>(fitParametersUnbinnedList.size()); ++fitParametersUnbinnedListIndex) {
      fitParametersUnbinnedFile << fitParametersUnbinnedList.at(fitParametersUnbinnedListIndex) << std::endl;
    }
    fitParametersUnbinnedFile.close();
    std::cout << "Unbinned fit parameters written to file: " << (options.outputFolder + "/unbinned_fitParameters_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".dat") << std::endl;

    // write parameters for binned fit
    std::ofstream fitParametersBinnedFile((options.outputFolder + "/binned_fitParameters_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".dat").c_str());
    assert(fitParametersBinnedFile.is_open());
    std::string fit_parameter_name;
    for (int nJetsBin = 3; nJetsBin <= 6; ++nJetsBin) {
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
  return EXIT_SUCCESS;
}
