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

template<typename T>
void printSquareMatrix(const T& matrixToPrint, const int& size) {
  std::cout << std::setprecision(3);
  std::cout << std::endl;
  for (int row_index = 0; row_index < size; ++row_index) {
    if (row_index == 0) std::cout << "(   ";
    else std::cout << "    ";
    for (int column_index = 0; column_index < size; ++column_index) {
      std::cout << matrixToPrint(row_index, column_index) << "    ";
    }
    if (row_index == (size - 1)) std::cout << ")";
    std::cout << std::endl;
  }
  std::cout << std::fixed;
}

template<typename T>
void printTVector(const T& vectorToPrint) {
  std::cout << std::setprecision(3);
  std::cout << std::endl << "(";
  for (int index = 0; index < vectorToPrint.GetNoElements(); ++index) std::cout << vectorToPrint(index) << "; ";
  std::cout << ")" << std::endl << std::fixed;
}

template<typename T>
void printVector(const std::vector<T> &vectorToPrint) {
  std::cout << std::setprecision(3);
  std::cout << std::endl << "(";
  for (int index = 0; index < static_cast<int>(vectorToPrint.size()); ++index) std::cout << vectorToPrint.at(index) << "; ";
  std::cout << ")" << std::endl << std::fixed;
}

std::vector<double> getColumnFromTMatrixD(const TMatrixD &source_matrix, const int &column_index, const int& size) {
  std::vector<double> column;
  for (int row_index = 0; row_index < size; ++row_index) column.push_back(source_matrix(row_index, column_index));
  return column;
}

template<typename T>
void check_eigendecomposition(const eigenvalue_eigenvector_pair_struct& pair, const T& matrix_to_check, const bool& print_debug=false) {
  std::vector<double> matrix_times_eigenvector;
  std::vector<double> eigenvalue_times_eigenvector;
  if (print_debug) {
    std::cout << "Checking eigendecomposition..." << std::endl;
    std::cout << "Testing eigenvalue: " << pair.eigenvalue << std::endl;
    std::cout << "Testing eigenvector: ";
    printVector(pair.eigenvector);
    std::cout << "Testing matrix: ";
    printSquareMatrix(matrix_to_check, static_cast<int>(pair.eigenvector.size()));
  }
  for (int row_index = 0; row_index < static_cast<int>(pair.eigenvector.size()); ++row_index) {
    double sum = 0.;
    for (int column_index = 0; column_index < static_cast<int>(pair.eigenvector.size()); ++column_index) {
      sum += matrix_to_check(row_index, column_index)*((pair.eigenvector).at(column_index));
    }
    if (print_debug) {
      std::cout << "At i: " << row_index << ", i'th component of matrix times eigenvector: " << sum << ", while i'th component of eigenvalue times eigenvector: " << pair.eigenvalue*((pair.eigenvector).at(row_index)) << std::endl;
      matrix_times_eigenvector.push_back(sum);
      eigenvalue_times_eigenvector.push_back(pair.eigenvalue*((pair.eigenvector).at(row_index)));
    }
    assert(std::fabs(sum - (pair.eigenvalue)*((pair.eigenvector).at(row_index))) < CHECK_TOLERANCE);
  }
  if (print_debug) {
    std::cout << "matrix_times_eigenvector: ";
    printVector(matrix_times_eigenvector);
    std::cout << "eigenvalue_times_eigenvector: ";
    printVector(eigenvalue_times_eigenvector);
  }
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

// double getDensityHistogramChiSquareWRTFunction(TH1D& inputDensityHistogram, std::map<int, double>& functionIntegralsOverBins) {
//   // first get number of entries and total integral for normalization
//   double weightedNEntries = 0.;
//   double sumFunctionIntegrals = 0.;
//   for (int binIndex = 1; binIndex <= inputDensityHistogram.GetXaxis()->GetNbins(); ++binIndex) {
//     weightedNEntries += (inputDensityHistogram.GetBinContent(binIndex))*(inputDensityHistogram.GetXaxis()->GetBinWidth(binIndex));
//     sumFunctionIntegrals += functionIntegralsOverBins.at(binIndex);
//   }
//   double functionNormalization = weightedNEntries/sumFunctionIntegrals;

//   double chiSquare = 0.;
//   for (int binIndex = 1; binIndex <= inputDensityHistogram.GetXaxis()->GetNbins(); ++binIndex) {
//     double observationFromDensityHistogram = (inputDensityHistogram.GetBinContent(binIndex))*(inputDensityHistogram.GetXaxis()->GetBinWidth(binIndex));
//     double error = (inputDensityHistogram.GetBinError(binIndex))*(inputDensityHistogram.GetXaxis()->GetBinWidth(binIndex));
//     double expectation_fromFunction = functionNormalization*functionIntegralsOverBins.at(binIndex);
//     // std::cout << "At binIndex = " << binIndex << ", bin center: " << inputDensityHistogram.GetXaxis()->GetBinCenter(binIndex) << ", observationFromDensityHistogram: " << observationFromDensityHistogram << ", functionIntegralsOverBins: " << functionIntegralsOverBins.at(binIndex) << ", expectation_fromFunction: " << expectation_fromFunction << ", error: " << error << std::endl;
//     if ((error > 0) && (observationFromDensityHistogram > 0)) {
//       chiSquare += (std::pow(observationFromDensityHistogram - expectation_fromFunction, 2))/error;
//     }
//   }
//   // std::cout << "chiSquare = " << chiSquare << std::endl;
//   return chiSquare;
// }

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

int getNNonEmptyBins(TH1D* inputHistogram) {
  int n_nonempty_bins = 0;
  for (int binCounter = 1; binCounter <= static_cast<int>(inputHistogram->GetXaxis()->GetNbins()); ++binCounter) {
    if ((inputHistogram->GetBinContent(binCounter)) > 0.) ++n_nonempty_bins;
  }
  return n_nonempty_bins;
}

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  TH1::AddDirectory(kFALSE);
  RooMsgService::instance().setGlobalKillBelow(MsgLevel::WARNING);
  ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.e-3);

  tmArgumentParser argumentParser = tmArgumentParser("Run script that prints useful info about the normalization.");
  argumentParser.addArgument("sourceFilePath", "", true, "Path to file containing list of paths with n-tuplized events.");
  argumentParser.addArgument("outputFolder", "", true, "Output folder.");
  argumentParser.addArgument("selection", "", true, "Name of selection: \"singlemedium\", \"signal_loose\", etc.");
  argumentParser.addArgument("identifier", "", true, "Identifier: \"MC_GJet17\", \"MC_GJet\", etc.");
  argumentParser.addArgument("STBoundariesSourceFile", "STRegionBoundaries_normOptimization.dat", false, "Identifier: \"MC_GJet17\", \"MC_GJet\", etc.");
  argumentParser.addArgument("yearString", "all", false, "String with year: can take values \"2016\", \"2017\", \"2018\", or \"all\".");
  argumentParser.addArgument("PDF_nSTBins", "25", false, "Number of bins for plotting datasets.");
  argumentParser.addArgument("preNormalizationBuffer", "200.0", false, "Buffer in ST to use before normalization bin for the kernel.");
  argumentParser.addArgument("readParametersFromFiles", "/dev/null,/dev/null", false, "If this argument is set, then no fits are performed; instead, the fit parameters is read in from the file locations given as the value of this argument. There should be precisely two files separated by a comma.");
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
  std::map<std::string, std::map<int, double> > fitParametersBinned;
  std::vector<std::string> fitParametersBinnedList;
  std::map<std::string, std::map<int, goodnessOfFitStruct> > fit_qualities_binned;
  std::map<std::string, std::map<int, double> > ftest_pValues;
  double slope_minVal = -1.0/(((ST_MAX_RANGE)/(options.STNormTarget)) - 1.0);
  double slope_maxVal = 5.0;
  double sqrt_minVal = -1.0/(std::sqrt((ST_MAX_RANGE)/(options.STNormTarget)) - 1.0);
  double sqrt_maxVal = 25.0;
  double quad_minVal = -1.0/((std::pow((ST_MAX_RANGE)/(options.STNormTarget), 2)) - 1.0);
  double quad_maxVal = 3.0;
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
    for (int nJetsBin = 3; nJetsBin <= 6; ++nJetsBin) {
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_fit_slope"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_fit_slope_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_fit_slopeError"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_fit_slopeError_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["sqrt_fit_sqrt"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "sqrt_fit_sqrt_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["sqrt_fit_sqrtError"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "sqrt_fit_sqrtError_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_sqrt_fit_slope"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_fit_slope_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_sqrt_fit_sqrt"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_fit_sqrt_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_sqrt_fit_mode1_slopeCoefficient"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_fit_mode1_slopeCoefficient_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_sqrt_fit_mode1_sqrtCoefficient"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_fit_mode1_sqrtCoefficient_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_sqrt_fit_mode1_error"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_fit_mode1_error_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_sqrt_fit_mode2_slopeCoefficient"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_fit_mode2_slopeCoefficient_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_sqrt_fit_mode2_sqrtCoefficient"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_fit_mode2_sqrtCoefficient_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_sqrt_fit_mode2_error"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_fit_mode2_error_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_sqrt_quad_fit_slope"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_quad_fit_slope_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_sqrt_quad_fit_sqrt"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_quad_fit_sqrt_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_sqrt_quad_fit_quad"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_quad_fit_quad_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_sqrt_quad_fit_mode1_slopeCoefficient"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_quad_fit_mode1_slopeCoefficient_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_sqrt_quad_fit_mode1_sqrtCoefficient"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_quad_fit_mode1_sqrtCoefficient_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_sqrt_quad_fit_mode1_quadCoefficient"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_quad_fit_mode1_quadCoefficient_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_sqrt_quad_fit_mode1_error"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_quad_fit_mode1_error_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_sqrt_quad_fit_mode2_slopeCoefficient"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_quad_fit_mode2_slopeCoefficient_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_sqrt_quad_fit_mode2_sqrtCoefficient"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_quad_fit_mode2_sqrtCoefficient_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_sqrt_quad_fit_mode2_quadCoefficient"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_quad_fit_mode2_quadCoefficient_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_sqrt_quad_fit_mode2_error"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_quad_fit_mode2_error_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_sqrt_quad_fit_mode3_slopeCoefficient"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_quad_fit_mode3_slopeCoefficient_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_sqrt_quad_fit_mode3_sqrtCoefficient"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_quad_fit_mode3_sqrtCoefficient_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_sqrt_quad_fit_mode3_quadCoefficient"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_quad_fit_mode3_quadCoefficient_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject_binned, lineFromFile)); fitParametersBinned["slope_sqrt_quad_fit_mode3_error"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_sqrt_quad_fit_mode3_error_" + std::to_string(nJetsBin) + "Jets");
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
  customizedPDF pdf_2Jets_customized(&pdf_2Jets, &rooVar_ST, options.STNormTarget, customizedPDF::customizationType::ScaleOnly);
  pdf_2Jets_customized.setScale("fitRange", ((STHistograms.at(2)).Integral(1, (STHistograms.at(2)).GetXaxis()->GetNbins(), "width")));
  TF1 pdf_2Jets_customized_TF1 = TF1("pdf_2Jets_customized_TF1", pdf_2Jets_customized, options.STRegions.STNormRangeMin, ST_MAX_RANGE, 0);
  pdf_2Jets_customized_TF1.SetLineColor(static_cast<EColor>(kBlue));
  pdf_2Jets_customized_TF1.SetLineWidth(1);
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
      std::vector<eigenvalue_eigenvector_pair_struct> eigenvalues_and_eigenvectors;
      for (int eigen_index = 0; eigen_index < fit_parameters.getSize(); ++eigen_index) {
        double eigenvalue = eigenvalues(eigen_index);
        std::vector<double> eigenvector = getColumnFromTMatrixD(eigenvectors, eigen_index, fit_parameters.getSize());
        eigenvalue_eigenvector_pair_struct current_pair = eigenvalue_eigenvector_pair_struct(eigenvalue, eigenvector);
        check_eigendecomposition(current_pair, covarianceMatrix);
        eigenvalues_and_eigenvectors.push_back(current_pair);
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
      fitParametersUnbinned["slope_sqrt_fit_mode1_slopeCoefficient"][nJetsBin] = (eigenvalues_and_eigenvectors.at(0)).eigenvector.at(0);
      fitParametersUnbinnedList.push_back(std::string("float slope_sqrt_fit_mode1_slopeCoefficient_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("slope_sqrt_fit_mode1_slopeCoefficient")).at(nJetsBin))));
      fitParametersUnbinned["slope_sqrt_fit_mode1_sqrtCoefficient"][nJetsBin] = (eigenvalues_and_eigenvectors.at(0)).eigenvector.at(1);
      fitParametersUnbinnedList.push_back(std::string("float slope_sqrt_fit_mode1_sqrtCoefficient_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("slope_sqrt_fit_mode1_sqrtCoefficient")).at(nJetsBin))));
      fitParametersUnbinned["slope_sqrt_fit_mode1_error"][nJetsBin] = std::sqrt((eigenvalues_and_eigenvectors.at(0)).eigenvalue);
      fitParametersUnbinnedList.push_back(std::string("float slope_sqrt_fit_mode1_error_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("slope_sqrt_fit_mode1_error")).at(nJetsBin))));
      fitParametersUnbinned["slope_sqrt_fit_mode2_slopeCoefficient"][nJetsBin] = (eigenvalues_and_eigenvectors.at(1)).eigenvector.at(0);
      fitParametersUnbinnedList.push_back(std::string("float slope_sqrt_fit_mode2_slopeCoefficient_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("slope_sqrt_fit_mode2_slopeCoefficient")).at(nJetsBin))));
      fitParametersUnbinned["slope_sqrt_fit_mode2_sqrtCoefficient"][nJetsBin] = (eigenvalues_and_eigenvectors.at(1)).eigenvector.at(1);
      fitParametersUnbinnedList.push_back(std::string("float slope_sqrt_fit_mode2_sqrtCoefficient_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("slope_sqrt_fit_mode2_sqrtCoefficient")).at(nJetsBin))));
      fitParametersUnbinned["slope_sqrt_fit_mode2_error"][nJetsBin] = std::sqrt((eigenvalues_and_eigenvectors.at(1)).eigenvalue);
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

    TCanvas pdfCanvas_unbinned = TCanvas(("c_dataSetAndPdf_unbinned_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("c_dataSetAndPdf_unbinned_" + std::to_string(nJetsBin) + "JetsBin").c_str(), 2560, 1440);
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
    TCanvas unbinned_shape_ratios_canvas = TCanvas(("c_unbinnedShapeRatios_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("c_unbinnedShapeRatios_" + std::to_string(nJetsBin) + "JetsBin").c_str(), 2560, 1440);
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
    std::string binnedFitOptions = "QSI+";

    // unadjusted
    TGraphErrors ratioGraph_binned_nJetsDistribution_to_unadjusted = TGraphErrors();
    ratioGraph_binned_nJetsDistribution_to_unadjusted.SetName(("ratioGraph_binned_nJetsDistribution_to_unadjusted_at" + std::to_string(nJetsBin) + "Jets").c_str());
    ratioGraph_binned_nJetsDistribution_to_unadjusted.SetTitle(("ST distribution at " + std::to_string(nJetsBin) + " Jets / unadjusted").c_str());
    customizedPDF pdf_2Jets_scaled(&pdf_2Jets, &rooVar_ST, options.STNormTarget, customizedPDF::customizationType::ScaleOnly);
    pdf_2Jets_scaled.setScale("normRange", ((STHistograms.at(nJetsBin)).GetBinContent(1))*((STHistograms.at(nJetsBin)).GetBinWidth(1)));
    TF1 pdf_2Jets_scaled_TF1 = TF1("pdf_2Jets_scaled_TF1", pdf_2Jets_scaled, options.STRegions.STNormRangeMin, ST_MAX_RANGE, 0);
    if (!(options.readParametersFromFiles)) {
      double chisq_value = (STHistograms.at(nJetsBin)).Chisquare(&pdf_2Jets_scaled_TF1, "R");
      fit_qualities_binned["unadjusted"][nJetsBin] = goodnessOfFitStruct(chisq_value, getNNonEmptyBins(&(STHistograms.at(nJetsBin))));
    }
    for (int binCounter = 1; binCounter <= (STHistograms.at(nJetsBin)).GetXaxis()->GetNbins(); ++binCounter) {
      double STMidpoint = (STHistograms.at(nJetsBin)).GetXaxis()->GetBinCenter(binCounter);
      double binWidth = (STHistograms.at(nJetsBin)).GetXaxis()->GetBinWidth(binCounter);
      double numerator = (STHistograms.at(nJetsBin)).GetBinContent(binCounter);
      double denominator = (pdf_2Jets_scaled_TF1.Integral((STHistograms.at(nJetsBin)).GetXaxis()->GetBinLowEdge(binCounter), (STHistograms.at(nJetsBin)).GetXaxis()->GetBinUpEdge(binCounter), 1.e-6))/((STHistograms.at(nJetsBin)).GetXaxis()->GetBinWidth(binCounter));
      assert(denominator > 0.);
      double ratio = numerator/denominator;
      double numeratorError = (STHistograms.at(nJetsBin)).GetBinError(binCounter);
      double denominatorError = 0.; // might change later
      double ratioError = ratio*std::sqrt(std::pow(numeratorError/numerator, 2) + std::pow(denominatorError/denominator, 2));
      int graph_currentPointIndex = ratioGraph_binned_nJetsDistribution_to_unadjusted.GetN();
      ratioGraph_binned_nJetsDistribution_to_unadjusted.SetPoint(graph_currentPointIndex, STMidpoint, ratio);
      ratioGraph_binned_nJetsDistribution_to_unadjusted.SetPointError(graph_currentPointIndex, binWidth/(std::sqrt(12)), ratioError);
    }

    // slope correction
    customizedPDF pdf_2Jets_scaled_slope(&pdf_2Jets, &rooVar_ST, options.STNormTarget, customizedPDF::customizationType::Slope);
    pdf_2Jets_scaled_slope.setScale("normRange", ((STHistograms.at(nJetsBin)).GetBinContent(1))*((STHistograms.at(nJetsBin)).GetBinWidth(1)));
    TF1 pdf_2Jets_scaled_slope_TF1 = TF1("pdf_2Jets_scaled_slope_TF1", pdf_2Jets_scaled_slope, options.STRegions.STNormRangeMin, ST_MAX_RANGE, 1);
    pdf_2Jets_scaled_slope_TF1.SetParName(0, ("slope_fit_slope_" + std::to_string(nJetsBin) + "JetsBin").c_str());
    pdf_2Jets_scaled_slope_TF1.SetParameter(0, (fitParametersUnbinned.at("slope_fit_slope")).at(nJetsBin)); // bootstrap fit
    pdf_2Jets_scaled_slope_TF1.SetParLimits(0, slope_minVal, slope_maxVal);
    if (options.readParametersFromFiles) {
      pdf_2Jets_scaled_slope_TF1.SetParameter(0, (fitParametersBinned.at("slope_fit_slope")).at(nJetsBin));
    }
    else {
      TFitResultPtr binned_fit_result = (STHistograms.at(nJetsBin)).Fit(&pdf_2Jets_scaled_slope_TF1, binnedFitOptions.c_str());
      assert(binned_fit_result->Status() == 0);
      assert(binned_fit_result->NTotalParameters() == 1);
      fit_qualities_binned["slope"][nJetsBin] = goodnessOfFitStruct(binned_fit_result->Chi2(), binned_fit_result->Ndf());
      // binned_fit_result->Print();
      fitParametersBinned["slope_fit_slope"][nJetsBin] = binned_fit_result->Parameter(0);
      fitParametersBinned["slope_fit_slopeError"][nJetsBin] = binned_fit_result->ParError(0);
      fitParametersBinnedList.push_back(std::string("float slope_fit_slope_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_fit_slope")).at(nJetsBin))));
      fitParametersBinnedList.push_back(std::string("float slope_fit_slopeError_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_fit_slopeError")).at(nJetsBin))));
    }

    // sqrt correction
    customizedPDF pdf_2Jets_scaled_sqrt(&pdf_2Jets, &rooVar_ST, options.STNormTarget, customizedPDF::customizationType::Sqrt);
    pdf_2Jets_scaled_sqrt.setScale("normRange", ((STHistograms.at(nJetsBin)).GetBinContent(1))*((STHistograms.at(nJetsBin)).GetBinWidth(1)));
    TF1 pdf_2Jets_scaled_sqrt_TF1 = TF1("pdf_2Jets_scaled_sqrt_TF1", pdf_2Jets_scaled_sqrt, options.STRegions.STNormRangeMin, ST_MAX_RANGE, 1);
    pdf_2Jets_scaled_sqrt_TF1.SetParName(0, ("sqrt_fit_sqrt_" + std::to_string(nJetsBin) + "JetsBin").c_str());
    pdf_2Jets_scaled_sqrt_TF1.SetParameter(0, (fitParametersUnbinned.at("sqrt_fit_sqrt")).at(nJetsBin)); // bootstrap fit
    pdf_2Jets_scaled_sqrt_TF1.SetParLimits(0, sqrt_minVal, sqrt_maxVal);
    if (options.readParametersFromFiles) {
      pdf_2Jets_scaled_sqrt_TF1.SetParameter(0, (fitParametersBinned.at("sqrt_fit_sqrt")).at(nJetsBin));
    }
    else {
      TFitResultPtr binned_fit_result = (STHistograms.at(nJetsBin)).Fit(&pdf_2Jets_scaled_sqrt_TF1, binnedFitOptions.c_str());
      assert(binned_fit_result->Status() == 0);
      assert(binned_fit_result->NTotalParameters() == 1);
      fit_qualities_binned["sqrt"][nJetsBin] = goodnessOfFitStruct(binned_fit_result->Chi2(), binned_fit_result->Ndf());
      // binned_fit_result->Print();
      fitParametersBinned["sqrt_fit_sqrt"][nJetsBin] = binned_fit_result->Parameter(0);
      fitParametersBinned["sqrt_fit_sqrtError"][nJetsBin] = binned_fit_result->ParError(0);
      fitParametersBinnedList.push_back(std::string("float sqrt_fit_sqrt_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("sqrt_fit_sqrt")).at(nJetsBin))));
      fitParametersBinnedList.push_back(std::string("float sqrt_fit_sqrtError_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("sqrt_fit_sqrtError")).at(nJetsBin))));
    }

    // slope correction + sqrt correction
    customizedPDF pdf_2Jets_scaled_slope_sqrt(&pdf_2Jets, &rooVar_ST, options.STNormTarget, customizedPDF::customizationType::SlopeSqrt);
    pdf_2Jets_scaled_slope_sqrt.setScale("normRange", ((STHistograms.at(nJetsBin)).GetBinContent(1))*((STHistograms.at(nJetsBin)).GetBinWidth(1)));
    TF1 pdf_2Jets_scaled_slope_sqrt_TF1 = TF1("pdf_2Jets_scaled_slope_sqrt_TF1", pdf_2Jets_scaled_slope_sqrt, options.STRegions.STNormRangeMin, ST_MAX_RANGE, 2);
    pdf_2Jets_scaled_slope_sqrt_TF1.SetParName(0, ("slope_sqrt_fit_slope_" + std::to_string(nJetsBin) + "JetsBin").c_str());
    pdf_2Jets_scaled_slope_sqrt_TF1.SetParameter(0, (fitParametersUnbinned.at("slope_sqrt_fit_slope")).at(nJetsBin)); // bootstrap fit
    pdf_2Jets_scaled_slope_sqrt_TF1.SetParLimits(0, slope_minVal, slope_maxVal);
    pdf_2Jets_scaled_slope_sqrt_TF1.SetParName(1, ("slope_sqrt_fit_sqrt_" + std::to_string(nJetsBin) + "JetsBin").c_str());
    pdf_2Jets_scaled_slope_sqrt_TF1.SetParameter(1, (fitParametersUnbinned.at("slope_sqrt_fit_sqrt")).at(nJetsBin)); // bootstrap fit
    pdf_2Jets_scaled_slope_sqrt_TF1.SetParLimits(1, sqrt_minVal, sqrt_maxVal);
    if (options.readParametersFromFiles) {
      pdf_2Jets_scaled_slope_sqrt_TF1.SetParameter(0, (fitParametersBinned.at("slope_sqrt_fit_slope")).at(nJetsBin));
      pdf_2Jets_scaled_slope_sqrt_TF1.SetParameter(1, (fitParametersBinned.at("slope_sqrt_fit_sqrt")).at(nJetsBin));
    }
    else {
      TFitResultPtr binned_fit_result = (STHistograms.at(nJetsBin)).Fit(&pdf_2Jets_scaled_slope_sqrt_TF1, binnedFitOptions.c_str());
      assert(binned_fit_result->Status() == 0);
      assert(binned_fit_result->NTotalParameters() == 2);
      fit_qualities_binned["slope_sqrt"][nJetsBin] = goodnessOfFitStruct(binned_fit_result->Chi2(), binned_fit_result->Ndf());
      // step 1: get covariance matrix
      TMatrixDSym covarianceMatrix = binned_fit_result->GetCovarianceMatrix();
      std::cout << "For binned combined slope + sqrt fit, covarianceMatrix: ";
      printSquareMatrix(covarianceMatrix, binned_fit_result->NTotalParameters());
      // step 2: get eigendecomposition
      TMatrixDSymEigen eigendecomposition_setup = TMatrixDSymEigen(covarianceMatrix);
      TVectorD eigenvalues = eigendecomposition_setup.GetEigenValues();
      std::cout << "eigenvalues: ";
      printTVector(eigenvalues);
      TMatrixD eigenvectors = eigendecomposition_setup.GetEigenVectors();
      std::cout << "eigenvectors: ";
      printSquareMatrix(eigenvectors, binned_fit_result->NTotalParameters());
      std::vector<eigenvalue_eigenvector_pair_struct> eigenvalues_and_eigenvectors;
      for (int eigen_index = 0; eigen_index < static_cast<int>(binned_fit_result->NTotalParameters()); ++eigen_index) {
        double eigenvalue = eigenvalues(eigen_index);
        std::vector<double> eigenvector = getColumnFromTMatrixD(eigenvectors, eigen_index, binned_fit_result->NTotalParameters());
        eigenvalue_eigenvector_pair_struct current_pair = eigenvalue_eigenvector_pair_struct(eigenvalue, eigenvector);
        check_eigendecomposition(current_pair, covarianceMatrix);
        eigenvalues_and_eigenvectors.push_back(current_pair);
      }
      fitParametersBinned["slope_sqrt_fit_slope"][nJetsBin] = binned_fit_result->Parameter(0);
      fitParametersBinned["slope_sqrt_fit_sqrt"][nJetsBin] = binned_fit_result->Parameter(1);
      fitParametersBinned["slope_sqrt_fit_mode1_slopeCoefficient"][nJetsBin] = (eigenvalues_and_eigenvectors.at(0)).eigenvector.at(0);
      fitParametersBinned["slope_sqrt_fit_mode1_sqrtCoefficient"][nJetsBin] = (eigenvalues_and_eigenvectors.at(0)).eigenvector.at(1);
      fitParametersBinned["slope_sqrt_fit_mode1_error"][nJetsBin] = std::sqrt((eigenvalues_and_eigenvectors.at(0)).eigenvalue);
      fitParametersBinned["slope_sqrt_fit_mode2_slopeCoefficient"][nJetsBin] = (eigenvalues_and_eigenvectors.at(1)).eigenvector.at(0);
      fitParametersBinned["slope_sqrt_fit_mode2_sqrtCoefficient"][nJetsBin] = (eigenvalues_and_eigenvectors.at(1)).eigenvector.at(1);
      fitParametersBinned["slope_sqrt_fit_mode2_error"][nJetsBin] = std::sqrt((eigenvalues_and_eigenvectors.at(1)).eigenvalue);
      fitParametersBinnedList.push_back(std::string("float slope_sqrt_fit_slope_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_sqrt_fit_slope")).at(nJetsBin))));
      fitParametersBinnedList.push_back(std::string("float slope_sqrt_fit_sqrt_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_sqrt_fit_sqrt")).at(nJetsBin))));
      fitParametersBinnedList.push_back(std::string("float slope_sqrt_fit_mode1_slopeCoefficient_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_sqrt_fit_mode1_slopeCoefficient")).at(nJetsBin))));
      fitParametersBinnedList.push_back(std::string("float slope_sqrt_fit_mode1_sqrtCoefficient_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_sqrt_fit_mode1_sqrtCoefficient")).at(nJetsBin))));
      fitParametersBinnedList.push_back(std::string("float slope_sqrt_fit_mode1_error_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_sqrt_fit_mode1_error")).at(nJetsBin))));
      fitParametersBinnedList.push_back(std::string("float slope_sqrt_fit_mode2_slopeCoefficient_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_sqrt_fit_mode2_slopeCoefficient")).at(nJetsBin))));
      fitParametersBinnedList.push_back(std::string("float slope_sqrt_fit_mode2_sqrtCoefficient_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_sqrt_fit_mode2_sqrtCoefficient")).at(nJetsBin))));
      fitParametersBinnedList.push_back(std::string("float slope_sqrt_fit_mode2_error_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_sqrt_fit_mode2_error")).at(nJetsBin))));
    }

    // slope + sqrt + quad correction
    customizedPDF pdf_2Jets_scaled_slope_sqrt_quad(&pdf_2Jets, &rooVar_ST, options.STNormTarget, customizedPDF::customizationType::SlopeSqrtQuad);
    pdf_2Jets_scaled_slope_sqrt_quad.setScale("normRange", ((STHistograms.at(nJetsBin)).GetBinContent(1))*((STHistograms.at(nJetsBin)).GetBinWidth(1)));
    TF1 pdf_2Jets_scaled_slope_sqrt_quad_TF1 = TF1("pdf_2Jets_scaled_slope_sqrt_quad_TF1", pdf_2Jets_scaled_slope_sqrt_quad, options.STRegions.STNormRangeMin, ST_MAX_RANGE, 3);
    pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParName(0, ("slope_sqrt_quad_fit_slope_" + std::to_string(nJetsBin) + "JetsBin").c_str());
    pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(0, (fitParametersBinned.at("slope_sqrt_fit_slope")).at(nJetsBin)); // bootstrap fit
    pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParLimits(0, slope_minVal, slope_maxVal);
    pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParName(1, ("slope_sqrt_quad_fit_sqrt_" + std::to_string(nJetsBin) + "JetsBin").c_str());
    pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(1, (fitParametersBinned.at("slope_sqrt_fit_sqrt")).at(nJetsBin)); // bootstrap fit
    pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParLimits(1, sqrt_minVal, sqrt_maxVal);
    pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParName(2, ("slope_sqrt_quad_fit_quad_" + std::to_string(nJetsBin) + "JetsBin").c_str());
    pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(2, 0.);
    pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParLimits(2, quad_minVal, quad_maxVal);
    if (options.readParametersFromFiles) {
      pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(0, (fitParametersBinned.at("slope_sqrt_quad_fit_slope")).at(nJetsBin));
      pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(1, (fitParametersBinned.at("slope_sqrt_quad_fit_sqrt")).at(nJetsBin));
      pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(2, (fitParametersBinned.at("slope_sqrt_quad_fit_quad")).at(nJetsBin));
    }
    else {
      TFitResultPtr binned_fit_result = (STHistograms.at(nJetsBin)).Fit(&pdf_2Jets_scaled_slope_sqrt_quad_TF1, binnedFitOptions.c_str());
      assert(binned_fit_result->Status() == 0);
      assert(binned_fit_result->NTotalParameters() == 3);
      fit_qualities_binned["slope_sqrt_quad"][nJetsBin] = goodnessOfFitStruct(binned_fit_result->Chi2(), binned_fit_result->Ndf());
      // step 1: get covariance matrix
      TMatrixDSym covarianceMatrix = binned_fit_result->GetCovarianceMatrix();
      std::cout << "For binned combined slope + sqrt fit, covarianceMatrix: ";
      printSquareMatrix(covarianceMatrix, binned_fit_result->NTotalParameters());
      // step 2: get eigendecomposition
      TMatrixDSymEigen eigendecomposition_setup = TMatrixDSymEigen(covarianceMatrix);
      TVectorD eigenvalues = eigendecomposition_setup.GetEigenValues();
      std::cout << "eigenvalues: ";
      printTVector(eigenvalues);
      TMatrixD eigenvectors = eigendecomposition_setup.GetEigenVectors();
      std::cout << "eigenvectors: ";
      printSquareMatrix(eigenvectors, binned_fit_result->NTotalParameters());
      std::vector<eigenvalue_eigenvector_pair_struct> eigenvalues_and_eigenvectors;
      for (int eigen_index = 0; eigen_index < static_cast<int>(binned_fit_result->NTotalParameters()); ++eigen_index) {
        double eigenvalue = eigenvalues(eigen_index);
        std::vector<double> eigenvector = getColumnFromTMatrixD(eigenvectors, eigen_index, binned_fit_result->NTotalParameters());
        eigenvalue_eigenvector_pair_struct current_pair = eigenvalue_eigenvector_pair_struct(eigenvalue, eigenvector);
        check_eigendecomposition(current_pair, covarianceMatrix);
        eigenvalues_and_eigenvectors.push_back(current_pair);
      }
      fitParametersBinned["slope_sqrt_quad_fit_slope"][nJetsBin] = binned_fit_result->Parameter(0);
      fitParametersBinned["slope_sqrt_quad_fit_sqrt"][nJetsBin] = binned_fit_result->Parameter(1);
      fitParametersBinned["slope_sqrt_quad_fit_quad"][nJetsBin] = binned_fit_result->Parameter(2);
      fitParametersBinned["slope_sqrt_quad_fit_mode1_slopeCoefficient"][nJetsBin] = (eigenvalues_and_eigenvectors.at(0)).eigenvector.at(0);
      fitParametersBinned["slope_sqrt_quad_fit_mode1_sqrtCoefficient"][nJetsBin] = (eigenvalues_and_eigenvectors.at(0)).eigenvector.at(1);
      fitParametersBinned["slope_sqrt_quad_fit_mode1_quadCoefficient"][nJetsBin] = (eigenvalues_and_eigenvectors.at(0)).eigenvector.at(2);
      fitParametersBinned["slope_sqrt_quad_fit_mode1_error"][nJetsBin] = std::sqrt((eigenvalues_and_eigenvectors.at(0)).eigenvalue);
      fitParametersBinned["slope_sqrt_quad_fit_mode2_slopeCoefficient"][nJetsBin] = (eigenvalues_and_eigenvectors.at(1)).eigenvector.at(0);
      fitParametersBinned["slope_sqrt_quad_fit_mode2_sqrtCoefficient"][nJetsBin] = (eigenvalues_and_eigenvectors.at(1)).eigenvector.at(1);
      fitParametersBinned["slope_sqrt_quad_fit_mode2_quadCoefficient"][nJetsBin] = (eigenvalues_and_eigenvectors.at(1)).eigenvector.at(2);
      fitParametersBinned["slope_sqrt_quad_fit_mode2_error"][nJetsBin] = std::sqrt((eigenvalues_and_eigenvectors.at(1)).eigenvalue);
      fitParametersBinned["slope_sqrt_quad_fit_mode3_slopeCoefficient"][nJetsBin] = (eigenvalues_and_eigenvectors.at(2)).eigenvector.at(0);
      fitParametersBinned["slope_sqrt_quad_fit_mode3_sqrtCoefficient"][nJetsBin] = (eigenvalues_and_eigenvectors.at(2)).eigenvector.at(1);
      fitParametersBinned["slope_sqrt_quad_fit_mode3_quadCoefficient"][nJetsBin] = (eigenvalues_and_eigenvectors.at(2)).eigenvector.at(2);
      fitParametersBinned["slope_sqrt_quad_fit_mode3_error"][nJetsBin] = std::sqrt((eigenvalues_and_eigenvectors.at(2)).eigenvalue);
      fitParametersBinnedList.push_back(std::string("float slope_sqrt_quad_fit_slope_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_sqrt_quad_fit_slope")).at(nJetsBin))));
      fitParametersBinnedList.push_back(std::string("float slope_sqrt_quad_fit_sqrt_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_sqrt_quad_fit_sqrt")).at(nJetsBin))));
      fitParametersBinnedList.push_back(std::string("float slope_sqrt_quad_fit_quad_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_sqrt_quad_fit_quad")).at(nJetsBin))));
      fitParametersBinnedList.push_back(std::string("float slope_sqrt_quad_fit_mode1_slopeCoefficient_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_slopeCoefficient")).at(nJetsBin))));
      fitParametersBinnedList.push_back(std::string("float slope_sqrt_quad_fit_mode1_sqrtCoefficient_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_sqrtCoefficient")).at(nJetsBin))));
      fitParametersBinnedList.push_back(std::string("float slope_sqrt_quad_fit_mode1_quadCoefficient_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_quadCoefficient")).at(nJetsBin))));
      fitParametersBinnedList.push_back(std::string("float slope_sqrt_quad_fit_mode1_error_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_error")).at(nJetsBin))));
      fitParametersBinnedList.push_back(std::string("float slope_sqrt_quad_fit_mode2_slopeCoefficient_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_sqrt_quad_fit_mode2_slopeCoefficient")).at(nJetsBin))));
      fitParametersBinnedList.push_back(std::string("float slope_sqrt_quad_fit_mode2_sqrtCoefficient_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_sqrt_quad_fit_mode2_sqrtCoefficient")).at(nJetsBin))));
      fitParametersBinnedList.push_back(std::string("float slope_sqrt_quad_fit_mode2_quadCoefficient_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_sqrt_quad_fit_mode2_quadCoefficient")).at(nJetsBin))));
      fitParametersBinnedList.push_back(std::string("float slope_sqrt_quad_fit_mode2_error_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_sqrt_quad_fit_mode2_error")).at(nJetsBin))));
      fitParametersBinnedList.push_back(std::string("float slope_sqrt_quad_fit_mode3_slopeCoefficient_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_sqrt_quad_fit_mode3_slopeCoefficient")).at(nJetsBin))));
      fitParametersBinnedList.push_back(std::string("float slope_sqrt_quad_fit_mode3_sqrtCoefficient_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_sqrt_quad_fit_mode3_sqrtCoefficient")).at(nJetsBin))));
      fitParametersBinnedList.push_back(std::string("float slope_sqrt_quad_fit_mode3_quadCoefficient_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_sqrt_quad_fit_mode3_quadCoefficient")).at(nJetsBin))));
      fitParametersBinnedList.push_back(std::string("float slope_sqrt_quad_fit_mode3_error_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersBinned.at("slope_sqrt_quad_fit_mode3_error")).at(nJetsBin))));
    }

    // initialize some variables useful for plots
    double fractionalError_normBin = ((STHistograms.at(nJetsBin)).GetBinError(1))/((STHistograms.at(nJetsBin)).GetBinContent(1));
    assert (fractionalError_normBin < 1.0); // sanity check, to make sure weights aren't affecting the errors in weird ways...
    double normBinCorrectionUp = 1.0 + fractionalError_normBin;
    double normBinCorrectionDown = 1.0 - fractionalError_normBin;

    // plot the raw shapes
    TCanvas pdfCanvas_binned = TCanvas(("c_dataSetAndPdf_binned_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("c_dataSetAndPdf_binned_" + std::to_string(nJetsBin) + "JetsBin").c_str(), 2560, 1440);
    TLegend legend_dataSetsAndPdf_binned = TLegend(0.6, 0.7, 0.9, 0.9);
    gStyle->SetOptStat(0);

    (STHistograms.at(nJetsBin)).SetLineColor(static_cast<EColor>(kBlack)); (STHistograms.at(nJetsBin)).Draw(); pdfCanvas_binned.Update();
    TLegendEntry *legendEntry_binned_nJetsDistribution = legend_dataSetsAndPdf_binned.AddEntry(&(STHistograms.at(nJetsBin)), (std::to_string(nJetsBin) + " jets distribution").c_str());
    legendEntry_binned_nJetsDistribution->SetMarkerColor(static_cast<EColor>(kBlack)); legendEntry_binned_nJetsDistribution->SetLineColor(static_cast<EColor>(kBlack)); legendEntry_binned_nJetsDistribution->SetTextColor(static_cast<EColor>(kBlack));
    (STHistograms.at(nJetsBin)).GetYaxis()->SetRange(((STHistograms.at(nJetsBin)).GetMaximum())/10000., ((STHistograms.at(nJetsBin)).GetMaximum())); pdfCanvas_binned.Update();

    TGraph pdf_2Jets_scaled_TF1_asTGraph = get_TF1_as_TGraph(&pdf_2Jets_scaled_TF1, 1000, options.STRegions.STNormRangeMin, ST_MAX_RANGE);
    pdf_2Jets_scaled_TF1_asTGraph.SetLineColor(static_cast<EColor>(kBlue)); pdf_2Jets_scaled_TF1_asTGraph.SetLineWidth(1);
    pdf_2Jets_scaled_TF1_asTGraph.Draw("CSAME"); pdfCanvas_binned.Update();
    TLegendEntry *legendEntry_binned_2JetsKernel = legend_dataSetsAndPdf_binned.AddEntry(&pdf_2Jets_scaled_TF1_asTGraph, "2 jets kernel, normalized");
    legendEntry_binned_2JetsKernel->SetMarkerColor(static_cast<EColor>(kBlue)); legendEntry_binned_2JetsKernel->SetLineColor(static_cast<EColor>(kBlue)); legendEntry_binned_2JetsKernel->SetTextColor(static_cast<EColor>(kBlue));

    pdf_2Jets_scaled_slope_TF1.SetParameter(0, (fitParametersBinned.at("slope_fit_slope")).at(nJetsBin));
    TGraph pdf_2Jets_scaled_slope_TF1_as_TGraph = get_TF1_as_TGraph(&pdf_2Jets_scaled_slope_TF1, 1000, options.STRegions.STNormRangeMin, ST_MAX_RANGE);
    pdf_2Jets_scaled_slope_TF1_as_TGraph.SetLineColor(static_cast<EColor>(kRed+1)); pdf_2Jets_scaled_slope_TF1_as_TGraph.SetLineWidth(1);
    pdf_2Jets_scaled_slope_TF1_as_TGraph.Draw("C"); pdfCanvas_binned.Update();
    pdf_2Jets_scaled_slope_TF1.SetParameter(0, ((fitParametersBinned.at("slope_fit_slope")).at(nJetsBin) + (fitParametersBinned.at("slope_fit_slopeError")).at(nJetsBin)));
    TGraph pdf_2Jets_scaled_slope_TF1_as_TGraph_up = get_TF1_as_TGraph(&pdf_2Jets_scaled_slope_TF1, 1000, options.STRegions.STNormRangeMin, ST_MAX_RANGE, normBinCorrectionUp);
    pdf_2Jets_scaled_slope_TF1_as_TGraph_up.SetLineStyle(kDashed); pdf_2Jets_scaled_slope_TF1_as_TGraph_up.SetLineColor(static_cast<EColor>(kRed+1)); pdf_2Jets_scaled_slope_TF1_as_TGraph_up.SetLineWidth(1);
    pdf_2Jets_scaled_slope_TF1_as_TGraph_up.Draw("C"); pdfCanvas_binned.Update();
    pdf_2Jets_scaled_slope_TF1.SetParameter(0, ((fitParametersBinned.at("slope_fit_slope")).at(nJetsBin) - (fitParametersBinned.at("slope_fit_slopeError")).at(nJetsBin)));
    TGraph pdf_2Jets_scaled_slope_TF1_as_TGraph_down = get_TF1_as_TGraph(&pdf_2Jets_scaled_slope_TF1, 1000, options.STRegions.STNormRangeMin, ST_MAX_RANGE, normBinCorrectionDown);
    pdf_2Jets_scaled_slope_TF1_as_TGraph_down.SetLineStyle(kDashed); pdf_2Jets_scaled_slope_TF1_as_TGraph_down.SetLineColor(static_cast<EColor>(kRed+1)); pdf_2Jets_scaled_slope_TF1_as_TGraph_down.SetLineWidth(1);
    pdf_2Jets_scaled_slope_TF1_as_TGraph_down.Draw("C"); pdfCanvas_binned.Update();
    TLegendEntry *legendEntry_binned_slope = legend_dataSetsAndPdf_binned.AddEntry(&pdf_2Jets_scaled_slope_TF1_as_TGraph, "2 jets kernel + linear adjustment");
    legendEntry_binned_slope->SetMarkerColor(static_cast<EColor>(kRed+1)); legendEntry_binned_slope->SetLineColor(static_cast<EColor>(kRed+1)); legendEntry_binned_slope->SetTextColor(static_cast<EColor>(kRed+1));

    pdf_2Jets_scaled_sqrt_TF1.SetParameter(0, (fitParametersBinned.at("sqrt_fit_sqrt")).at(nJetsBin));
    TGraph pdf_2Jets_scaled_sqrt_TF1_as_TGraph = get_TF1_as_TGraph(&pdf_2Jets_scaled_sqrt_TF1, 1000, options.STRegions.STNormRangeMin, ST_MAX_RANGE);
    pdf_2Jets_scaled_sqrt_TF1_as_TGraph.SetLineColor(static_cast<EColor>(kGreen+3)); pdf_2Jets_scaled_sqrt_TF1_as_TGraph.SetLineWidth(1);
    pdf_2Jets_scaled_sqrt_TF1_as_TGraph.Draw("C"); pdfCanvas_binned.Update();
    pdf_2Jets_scaled_sqrt_TF1.SetParameter(0, ((fitParametersBinned.at("sqrt_fit_sqrt")).at(nJetsBin) + (fitParametersBinned.at("sqrt_fit_sqrtError")).at(nJetsBin)));
    TGraph pdf_2Jets_scaled_sqrt_TF1_as_TGraph_up = get_TF1_as_TGraph(&pdf_2Jets_scaled_sqrt_TF1, 1000, options.STRegions.STNormRangeMin, ST_MAX_RANGE, normBinCorrectionUp);
    pdf_2Jets_scaled_sqrt_TF1_as_TGraph_up.SetLineStyle(kDashed); pdf_2Jets_scaled_sqrt_TF1_as_TGraph_up.SetLineColor(static_cast<EColor>(kGreen+3)); pdf_2Jets_scaled_sqrt_TF1_as_TGraph_up.SetLineWidth(1);
    pdf_2Jets_scaled_sqrt_TF1_as_TGraph_up.Draw("C"); pdfCanvas_binned.Update();
    pdf_2Jets_scaled_sqrt_TF1.SetParameter(0, ((fitParametersBinned.at("sqrt_fit_sqrt")).at(nJetsBin) - (fitParametersBinned.at("sqrt_fit_sqrtError")).at(nJetsBin)));
    TGraph pdf_2Jets_scaled_sqrt_TF1_as_TGraph_down = get_TF1_as_TGraph(&pdf_2Jets_scaled_sqrt_TF1, 1000, options.STRegions.STNormRangeMin, ST_MAX_RANGE, normBinCorrectionDown);
    pdf_2Jets_scaled_sqrt_TF1_as_TGraph_down.SetLineStyle(kDashed); pdf_2Jets_scaled_sqrt_TF1_as_TGraph_down.SetLineColor(static_cast<EColor>(kGreen+3)); pdf_2Jets_scaled_sqrt_TF1_as_TGraph_down.SetLineWidth(1);
    pdf_2Jets_scaled_sqrt_TF1_as_TGraph_down.Draw("C"); pdfCanvas_binned.Update();
    TLegendEntry *legendEntry_binned_sqrt = legend_dataSetsAndPdf_binned.AddEntry(&pdf_2Jets_scaled_sqrt_TF1_as_TGraph, "2 jets kernel + sqrt adjustment");
    legendEntry_binned_sqrt->SetMarkerColor(static_cast<EColor>(kGreen+3)); legendEntry_binned_sqrt->SetLineColor(static_cast<EColor>(kGreen+3)); legendEntry_binned_sqrt->SetTextColor(static_cast<EColor>(kGreen+3));

    pdf_2Jets_scaled_slope_sqrt_TF1.SetParameter(0, (fitParametersBinned.at("slope_sqrt_fit_slope")).at(nJetsBin));
    pdf_2Jets_scaled_slope_sqrt_TF1.SetParameter(1, (fitParametersBinned.at("slope_sqrt_fit_sqrt")).at(nJetsBin));
    TGraph pdf_2Jets_scaled_slope_sqrt_TF1_as_TGraph = get_TF1_as_TGraph(&pdf_2Jets_scaled_slope_sqrt_TF1, 1000, options.STRegions.STNormRangeMin, ST_MAX_RANGE);
    pdf_2Jets_scaled_slope_sqrt_TF1_as_TGraph.SetLineColor(static_cast<EColor>(kViolet)); pdf_2Jets_scaled_slope_sqrt_TF1_as_TGraph.SetLineWidth(1);
    pdf_2Jets_scaled_slope_sqrt_TF1_as_TGraph.Draw("C"); pdfCanvas_binned.Update();
    pdf_2Jets_scaled_slope_sqrt_TF1.SetParameter(0, (fitParametersBinned.at("slope_sqrt_fit_slope")).at(nJetsBin) + ((fitParametersBinned.at("slope_sqrt_fit_mode1_error")).at(nJetsBin))*((fitParametersBinned.at("slope_sqrt_fit_mode1_slopeCoefficient")).at(nJetsBin)));
    pdf_2Jets_scaled_slope_sqrt_TF1.SetParameter(1, (fitParametersBinned.at("slope_sqrt_fit_sqrt")).at(nJetsBin) + ((fitParametersBinned.at("slope_sqrt_fit_mode1_error")).at(nJetsBin))*((fitParametersBinned.at("slope_sqrt_fit_mode1_sqrtCoefficient")).at(nJetsBin)));
    TGraph pdf_2Jets_scaled_slope_sqrt_TF1_as_TGraph_up = get_TF1_as_TGraph(&pdf_2Jets_scaled_slope_sqrt_TF1, 1000, options.STRegions.STNormRangeMin, ST_MAX_RANGE, normBinCorrectionUp);
    pdf_2Jets_scaled_slope_sqrt_TF1_as_TGraph_up.SetLineStyle(kDashed); pdf_2Jets_scaled_slope_sqrt_TF1_as_TGraph_up.SetLineColor(static_cast<EColor>(kViolet)); pdf_2Jets_scaled_slope_sqrt_TF1_as_TGraph_up.SetLineWidth(1);
    pdf_2Jets_scaled_slope_sqrt_TF1_as_TGraph_up.Draw("C"); pdfCanvas_binned.Update();
    pdf_2Jets_scaled_slope_sqrt_TF1.SetParameter(0, (fitParametersBinned.at("slope_sqrt_fit_slope")).at(nJetsBin) - ((fitParametersBinned.at("slope_sqrt_fit_mode1_error")).at(nJetsBin))*((fitParametersBinned.at("slope_sqrt_fit_mode1_slopeCoefficient")).at(nJetsBin)));
    pdf_2Jets_scaled_slope_sqrt_TF1.SetParameter(1, (fitParametersBinned.at("slope_sqrt_fit_sqrt")).at(nJetsBin) - ((fitParametersBinned.at("slope_sqrt_fit_mode1_error")).at(nJetsBin))*((fitParametersBinned.at("slope_sqrt_fit_mode1_sqrtCoefficient")).at(nJetsBin)));
    TGraph pdf_2Jets_scaled_slope_sqrt_TF1_as_TGraph_down = get_TF1_as_TGraph(&pdf_2Jets_scaled_slope_sqrt_TF1, 1000, options.STRegions.STNormRangeMin, ST_MAX_RANGE, normBinCorrectionDown);
    pdf_2Jets_scaled_slope_sqrt_TF1_as_TGraph_down.SetLineStyle(kDashed); pdf_2Jets_scaled_slope_sqrt_TF1_as_TGraph_down.SetLineColor(static_cast<EColor>(kViolet)); pdf_2Jets_scaled_slope_sqrt_TF1_as_TGraph_down.SetLineWidth(1);
    pdf_2Jets_scaled_slope_sqrt_TF1_as_TGraph_down.Draw("C"); pdfCanvas_binned.Update();
    TLegendEntry *legendEntry_binned_slope_sqrt = legend_dataSetsAndPdf_binned.AddEntry(&pdf_2Jets_scaled_slope_sqrt_TF1_as_TGraph, "2 jets kernel + (linear+sqrt) adjustment");
    legendEntry_binned_slope_sqrt->SetMarkerColor(static_cast<EColor>(kViolet)); legendEntry_binned_slope_sqrt->SetLineColor(static_cast<EColor>(kViolet)); legendEntry_binned_slope_sqrt->SetTextColor(static_cast<EColor>(kViolet));

    pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(0, (fitParametersBinned.at("slope_sqrt_quad_fit_slope")).at(nJetsBin));
    pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(1, (fitParametersBinned.at("slope_sqrt_quad_fit_sqrt")).at(nJetsBin));
    pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(2, (fitParametersBinned.at("slope_sqrt_quad_fit_quad")).at(nJetsBin));
    TGraph pdf_2Jets_scaled_slope_sqrt_quad_TF1_as_TGraph = get_TF1_as_TGraph(&pdf_2Jets_scaled_slope_sqrt_quad_TF1, 1000, options.STRegions.STNormRangeMin, ST_MAX_RANGE);
    pdf_2Jets_scaled_slope_sqrt_quad_TF1_as_TGraph.SetLineColor(static_cast<EColor>(kYellow+2)); pdf_2Jets_scaled_slope_sqrt_quad_TF1_as_TGraph.SetLineWidth(1);
    pdf_2Jets_scaled_slope_sqrt_quad_TF1_as_TGraph.Draw("C"); pdfCanvas_binned.Update();
    pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(0, (fitParametersBinned.at("slope_sqrt_quad_fit_slope")).at(nJetsBin) + ((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_error")).at(nJetsBin))*((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_slopeCoefficient")).at(nJetsBin)));
    pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(1, (fitParametersBinned.at("slope_sqrt_quad_fit_sqrt")).at(nJetsBin) + ((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_error")).at(nJetsBin))*((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_sqrtCoefficient")).at(nJetsBin)));
    pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(2, (fitParametersBinned.at("slope_sqrt_quad_fit_quad")).at(nJetsBin) + ((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_error")).at(nJetsBin))*((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_quadCoefficient")).at(nJetsBin)));
    TGraph pdf_2Jets_scaled_slope_sqrt_quad_TF1_as_TGraph_up = get_TF1_as_TGraph(&pdf_2Jets_scaled_slope_sqrt_quad_TF1, 1000, options.STRegions.STNormRangeMin, ST_MAX_RANGE, normBinCorrectionUp);
    pdf_2Jets_scaled_slope_sqrt_quad_TF1_as_TGraph_up.SetLineStyle(kDashed); pdf_2Jets_scaled_slope_sqrt_quad_TF1_as_TGraph_up.SetLineColor(static_cast<EColor>(kYellow+2)); pdf_2Jets_scaled_slope_sqrt_quad_TF1_as_TGraph_up.SetLineWidth(1);
    pdf_2Jets_scaled_slope_sqrt_quad_TF1_as_TGraph_up.Draw("C"); pdfCanvas_binned.Update();
    pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(0, (fitParametersBinned.at("slope_sqrt_quad_fit_slope")).at(nJetsBin) - ((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_error")).at(nJetsBin))*((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_slopeCoefficient")).at(nJetsBin)));
    pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(1, (fitParametersBinned.at("slope_sqrt_quad_fit_sqrt")).at(nJetsBin) - ((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_error")).at(nJetsBin))*((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_sqrtCoefficient")).at(nJetsBin)));
    pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(2, (fitParametersBinned.at("slope_sqrt_quad_fit_quad")).at(nJetsBin) - ((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_error")).at(nJetsBin))*((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_quadCoefficient")).at(nJetsBin)));
    TGraph pdf_2Jets_scaled_slope_sqrt_quad_TF1_as_TGraph_down = get_TF1_as_TGraph(&pdf_2Jets_scaled_slope_sqrt_quad_TF1, 1000, options.STRegions.STNormRangeMin, ST_MAX_RANGE, normBinCorrectionDown);
    pdf_2Jets_scaled_slope_sqrt_quad_TF1_as_TGraph_down.SetLineStyle(kDashed); pdf_2Jets_scaled_slope_sqrt_quad_TF1_as_TGraph_down.SetLineColor(static_cast<EColor>(kYellow+2)); pdf_2Jets_scaled_slope_sqrt_quad_TF1_as_TGraph_down.SetLineWidth(1);
    pdf_2Jets_scaled_slope_sqrt_quad_TF1_as_TGraph_down.Draw("C"); pdfCanvas_binned.Update();
    TLegendEntry *legendEntry_binned_slope_sqrt_quad = legend_dataSetsAndPdf_binned.AddEntry(&pdf_2Jets_scaled_slope_sqrt_quad_TF1_as_TGraph, "2 jets kernel + (linear+sqrt+quad) adjustment");
    legendEntry_binned_slope_sqrt_quad->SetMarkerColor(static_cast<EColor>(kYellow+2)); legendEntry_binned_slope_sqrt_quad->SetLineColor(static_cast<EColor>(kYellow+2)); legendEntry_binned_slope_sqrt_quad->SetTextColor(static_cast<EColor>(kYellow+2));

    gPad->SetLogy(); pdfCanvas_binned.Update();
    legend_dataSetsAndPdf_binned.SetFillStyle(0); legend_dataSetsAndPdf_binned.Draw(); pdfCanvas_binned.Update();
    pdfCanvas_binned.SaveAs((options.outputFolder + "/binned_pdfAndData_" + std::to_string(nJetsBin) + "JetsBin_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".pdf").c_str());

    // calculate shape ratios
    double pdf_2Jets_scaled_TF1_atNorm = pdf_2Jets_scaled_TF1.Eval(options.STNormTarget);

    TGraph ratioGraph_binned_slope_to_unadjusted = TGraph();
    ratioGraph_binned_slope_to_unadjusted.SetName(("ratioGraph_binned_slope_to_unadjusted_at" + std::to_string(nJetsBin) + "Jets").c_str());
    ratioGraph_binned_slope_to_unadjusted.SetTitle((std::to_string(nJetsBin) + " Jets linear fit / 2 Jets kernel").c_str());
    TGraph ratioGraph_binned_slope_to_unadjusted_high_estimate = TGraph();
    ratioGraph_binned_slope_to_unadjusted_high_estimate.SetName(("ratioGraph_binned_slope_to_unadjusted_high_estimate_at" + std::to_string(nJetsBin) + "Jets").c_str());
    TGraph ratioGraph_binned_slope_to_unadjusted_low_estimate = TGraph();
    ratioGraph_binned_slope_to_unadjusted_low_estimate.SetName(("ratioGraph_binned_slope_to_unadjusted_low_estimate_at" + std::to_string(nJetsBin) + "Jets").c_str());
    pdf_2Jets_scaled_slope_TF1.SetParameter(0, ((fitParametersBinned.at("slope_fit_slope")).at(nJetsBin)));
    double pdf_2Jets_scaled_slope_TF1_atNorm = pdf_2Jets_scaled_slope_TF1.Eval(options.STNormTarget);

    TGraph ratioGraph_binned_sqrt_to_unadjusted = TGraph();
    ratioGraph_binned_sqrt_to_unadjusted.SetName(("ratioGraph_binned_sqrt_to_unadjusted_at" + std::to_string(nJetsBin) + "Jets").c_str());
    ratioGraph_binned_sqrt_to_unadjusted.SetTitle((std::to_string(nJetsBin) + " Jets sqrt fit / 2 Jets kernel").c_str());
    TGraph ratioGraph_binned_sqrt_to_unadjusted_high_estimate = TGraph();
    ratioGraph_binned_sqrt_to_unadjusted_high_estimate.SetName(("ratioGraph_binned_sqrt_to_unadjusted_high_estimate_at" + std::to_string(nJetsBin) + "Jets").c_str());
    TGraph ratioGraph_binned_sqrt_to_unadjusted_low_estimate = TGraph();
    ratioGraph_binned_sqrt_to_unadjusted_low_estimate.SetName(("ratioGraph_binned_sqrt_to_unadjusted_low_estimate_at" + std::to_string(nJetsBin) + "Jets").c_str());
    pdf_2Jets_scaled_sqrt_TF1.SetParameter(0, ((fitParametersBinned.at("sqrt_fit_sqrt")).at(nJetsBin)));
    double pdf_2Jets_scaled_sqrt_TF1_atNorm = pdf_2Jets_scaled_sqrt_TF1.Eval(options.STNormTarget);

    TGraph ratioGraph_binned_slope_sqrt_to_unadjusted = TGraph();
    ratioGraph_binned_slope_sqrt_to_unadjusted.SetName(("ratioGraph_binned_slope_sqrt_to_unadjusted_at" + std::to_string(nJetsBin) + "Jets").c_str());
    ratioGraph_binned_slope_sqrt_to_unadjusted.SetTitle((std::to_string(nJetsBin) + " Jets slope+sqrt fit / 2 Jets kernel").c_str());
    TGraph ratioGraph_binned_slope_sqrt_to_unadjusted_high_estimate = TGraph();
    ratioGraph_binned_slope_sqrt_to_unadjusted_high_estimate.SetName(("ratioGraph_binned_slope_sqrt_to_unadjusted_high_estimate_at" + std::to_string(nJetsBin) + "Jets").c_str());
    TGraph ratioGraph_binned_slope_sqrt_to_unadjusted_low_estimate = TGraph();
    ratioGraph_binned_slope_sqrt_to_unadjusted_low_estimate.SetName(("ratioGraph_binned_slope_sqrt_to_unadjusted_low_estimate_at" + std::to_string(nJetsBin) + "Jets").c_str());
    pdf_2Jets_scaled_slope_sqrt_TF1.SetParameter(0, (fitParametersBinned.at("slope_sqrt_fit_slope")).at(nJetsBin));
    pdf_2Jets_scaled_slope_sqrt_TF1.SetParameter(1, (fitParametersBinned.at("slope_sqrt_fit_sqrt")).at(nJetsBin));
    double pdf_2Jets_scaled_slope_sqrt_TF1_atNorm = pdf_2Jets_scaled_slope_sqrt_TF1.Eval(options.STNormTarget);

    TGraph ratioGraph_binned_slope_sqrt_quad_to_unadjusted = TGraph();
    ratioGraph_binned_slope_sqrt_quad_to_unadjusted.SetName(("ratioGraph_binned_slope_sqrt_quad_to_unadjusted_at" + std::to_string(nJetsBin) + "Jets").c_str());
    ratioGraph_binned_slope_sqrt_quad_to_unadjusted.SetTitle((std::to_string(nJetsBin) + " Jets slope+sqrt+quad fit / 2 Jets kernel").c_str());
    TGraph ratioGraph_binned_slope_sqrt_quad_to_unadjusted_high_estimate = TGraph();
    ratioGraph_binned_slope_sqrt_quad_to_unadjusted_high_estimate.SetName(("ratioGraph_binned_slope_sqrt_quad_to_unadjusted_high_estimate_at" + std::to_string(nJetsBin) + "Jets").c_str());
    TGraph ratioGraph_binned_slope_sqrt_quad_to_unadjusted_low_estimate = TGraph();
    ratioGraph_binned_slope_sqrt_quad_to_unadjusted_low_estimate.SetName(("ratioGraph_binned_slope_sqrt_quad_to_unadjusted_low_estimate_at" + std::to_string(nJetsBin) + "Jets").c_str());
    pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(0, (fitParametersBinned.at("slope_sqrt_quad_fit_slope")).at(nJetsBin));
    pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(1, (fitParametersBinned.at("slope_sqrt_quad_fit_sqrt")).at(nJetsBin));
    pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(2, (fitParametersBinned.at("slope_sqrt_quad_fit_quad")).at(nJetsBin));
    double pdf_2Jets_scaled_slope_sqrt_quad_TF1_atNorm = pdf_2Jets_scaled_slope_sqrt_quad_TF1.Eval(options.STNormTarget);

    for (int STCounter = 0; STCounter <= 1000; ++STCounter) {
      double STVal = options.STRegions.STNormRangeMin + (1.0*STCounter/1000)*(ST_MAX_RANGE - options.STRegions.STNormRangeMin);
      double common_denominator = pdf_2Jets_scaled_TF1.Eval(STVal)/pdf_2Jets_scaled_TF1_atNorm;

      pdf_2Jets_scaled_slope_TF1.SetParameter(0, ((fitParametersBinned.at("slope_fit_slope")).at(nJetsBin)));
      double pdf_slope_nominal = pdf_2Jets_scaled_slope_TF1.Eval(STVal)/pdf_2Jets_scaled_slope_TF1_atNorm;
      pdf_2Jets_scaled_slope_TF1.SetParameter(0, ((fitParametersBinned.at("slope_fit_slope")).at(nJetsBin) + (fitParametersBinned.at("slope_fit_slopeError")).at(nJetsBin)));
      double pdf_slope_plus_one_sigma = (pdf_2Jets_scaled_slope_TF1.Eval(STVal))/pdf_2Jets_scaled_slope_TF1_atNorm;
      pdf_2Jets_scaled_slope_TF1.SetParameter(0, ((fitParametersBinned.at("slope_fit_slope")).at(nJetsBin) - (fitParametersBinned.at("slope_fit_slopeError")).at(nJetsBin)));
      double pdf_slope_minus_one_sigma = (pdf_2Jets_scaled_slope_TF1.Eval(STVal))/pdf_2Jets_scaled_slope_TF1_atNorm;
      double pdf_slope_higher = std::max({pdf_slope_plus_one_sigma, pdf_slope_nominal, pdf_slope_minus_one_sigma});
      double pdf_slope_lower = std::min({pdf_slope_plus_one_sigma, pdf_slope_nominal, pdf_slope_minus_one_sigma});
      double ratio_slope_to_unadjusted = pdf_slope_nominal/common_denominator;
      double ratio_slope_to_unadjusted_higher = pdf_slope_higher/common_denominator;
      double ratio_slope_to_unadjusted_lower = pdf_slope_lower/common_denominator;
      ratioGraph_binned_slope_to_unadjusted.SetPoint(STCounter, STVal, ratio_slope_to_unadjusted);
      ratioGraph_binned_slope_to_unadjusted_high_estimate.SetPoint(STCounter, STVal, normBinCorrectionUp*ratio_slope_to_unadjusted_higher);
      ratioGraph_binned_slope_to_unadjusted_low_estimate.SetPoint(STCounter, STVal, normBinCorrectionDown*ratio_slope_to_unadjusted_lower);

      pdf_2Jets_scaled_sqrt_TF1.SetParameter(0, ((fitParametersBinned.at("sqrt_fit_sqrt")).at(nJetsBin)));
      double pdf_sqrt_nominal = pdf_2Jets_scaled_sqrt_TF1.Eval(STVal)/pdf_2Jets_scaled_sqrt_TF1_atNorm;
      pdf_2Jets_scaled_sqrt_TF1.SetParameter(0, ((fitParametersBinned.at("sqrt_fit_sqrt")).at(nJetsBin) + (fitParametersBinned.at("sqrt_fit_sqrtError")).at(nJetsBin)));
      double pdf_sqrt_plus_one_sigma = (pdf_2Jets_scaled_sqrt_TF1.Eval(STVal))/pdf_2Jets_scaled_sqrt_TF1_atNorm;
      pdf_2Jets_scaled_sqrt_TF1.SetParameter(0, ((fitParametersBinned.at("sqrt_fit_sqrt")).at(nJetsBin) - (fitParametersBinned.at("sqrt_fit_sqrtError")).at(nJetsBin)));
      double pdf_sqrt_minus_one_sigma = (pdf_2Jets_scaled_sqrt_TF1.Eval(STVal))/pdf_2Jets_scaled_sqrt_TF1_atNorm;
      double pdf_sqrt_higher = std::max({pdf_sqrt_plus_one_sigma, pdf_sqrt_nominal, pdf_sqrt_minus_one_sigma});
      double pdf_sqrt_lower = std::min({pdf_sqrt_plus_one_sigma, pdf_sqrt_nominal, pdf_sqrt_minus_one_sigma});
      double ratio_sqrt_to_unadjusted = pdf_sqrt_nominal/common_denominator;
      double ratio_sqrt_to_unadjusted_higher = pdf_sqrt_higher/common_denominator;
      double ratio_sqrt_to_unadjusted_lower = pdf_sqrt_lower/common_denominator;
      ratioGraph_binned_sqrt_to_unadjusted.SetPoint(STCounter, STVal, ratio_sqrt_to_unadjusted);
      ratioGraph_binned_sqrt_to_unadjusted_high_estimate.SetPoint(STCounter, STVal, normBinCorrectionUp*ratio_sqrt_to_unadjusted_higher);
      ratioGraph_binned_sqrt_to_unadjusted_low_estimate.SetPoint(STCounter, STVal, normBinCorrectionDown*ratio_sqrt_to_unadjusted_lower);

      pdf_2Jets_scaled_slope_sqrt_TF1.SetParameter(0, (fitParametersBinned.at("slope_sqrt_fit_slope")).at(nJetsBin));
      pdf_2Jets_scaled_slope_sqrt_TF1.SetParameter(1, (fitParametersBinned.at("slope_sqrt_fit_sqrt")).at(nJetsBin));
      double pdf_slope_sqrt_nominal = pdf_2Jets_scaled_slope_sqrt_TF1.Eval(STVal)/pdf_2Jets_scaled_slope_sqrt_TF1_atNorm;
      pdf_2Jets_scaled_slope_sqrt_TF1.SetParameter(0, (fitParametersBinned.at("slope_sqrt_fit_slope")).at(nJetsBin) + ((fitParametersBinned.at("slope_sqrt_fit_mode1_error")).at(nJetsBin))*((fitParametersBinned.at("slope_sqrt_fit_mode1_slopeCoefficient")).at(nJetsBin)));
      pdf_2Jets_scaled_slope_sqrt_TF1.SetParameter(1, (fitParametersBinned.at("slope_sqrt_fit_sqrt")).at(nJetsBin) + ((fitParametersBinned.at("slope_sqrt_fit_mode1_error")).at(nJetsBin))*((fitParametersBinned.at("slope_sqrt_fit_mode1_sqrtCoefficient")).at(nJetsBin)));
      double pdf_slope_sqrt_plus_one_sigma = (pdf_2Jets_scaled_slope_sqrt_TF1.Eval(STVal))/pdf_2Jets_scaled_slope_sqrt_TF1_atNorm;
      pdf_2Jets_scaled_slope_sqrt_TF1.SetParameter(0, (fitParametersBinned.at("slope_sqrt_fit_slope")).at(nJetsBin) - ((fitParametersBinned.at("slope_sqrt_fit_mode1_error")).at(nJetsBin))*((fitParametersBinned.at("slope_sqrt_fit_mode1_slopeCoefficient")).at(nJetsBin)));
      pdf_2Jets_scaled_slope_sqrt_TF1.SetParameter(1, (fitParametersBinned.at("slope_sqrt_fit_sqrt")).at(nJetsBin) - ((fitParametersBinned.at("slope_sqrt_fit_mode1_error")).at(nJetsBin))*((fitParametersBinned.at("slope_sqrt_fit_mode1_sqrtCoefficient")).at(nJetsBin)));
      double pdf_slope_sqrt_minus_one_sigma = (pdf_2Jets_scaled_slope_sqrt_TF1.Eval(STVal))/pdf_2Jets_scaled_slope_sqrt_TF1_atNorm;
      double pdf_slope_sqrt_higher = std::max({pdf_slope_sqrt_plus_one_sigma, pdf_slope_sqrt_nominal, pdf_slope_sqrt_minus_one_sigma});
      double pdf_slope_sqrt_lower = std::min({pdf_slope_sqrt_plus_one_sigma, pdf_slope_sqrt_nominal, pdf_slope_sqrt_minus_one_sigma});
      double ratio_slope_sqrt_to_unadjusted = pdf_slope_sqrt_nominal/common_denominator;
      double ratio_slope_sqrt_to_unadjusted_higher = pdf_slope_sqrt_higher/common_denominator;
      double ratio_slope_sqrt_to_unadjusted_lower = pdf_slope_sqrt_lower/common_denominator;
      ratioGraph_binned_slope_sqrt_to_unadjusted.SetPoint(STCounter, STVal, ratio_slope_sqrt_to_unadjusted);
      ratioGraph_binned_slope_sqrt_to_unadjusted_high_estimate.SetPoint(STCounter, STVal, normBinCorrectionUp*ratio_slope_sqrt_to_unadjusted_higher);
      ratioGraph_binned_slope_sqrt_to_unadjusted_low_estimate.SetPoint(STCounter, STVal, normBinCorrectionDown*ratio_slope_sqrt_to_unadjusted_lower);

      pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(0, (fitParametersBinned.at("slope_sqrt_quad_fit_slope")).at(nJetsBin));
      pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(1, (fitParametersBinned.at("slope_sqrt_quad_fit_sqrt")).at(nJetsBin));
      pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(2, (fitParametersBinned.at("slope_sqrt_quad_fit_quad")).at(nJetsBin));
      double pdf_slope_sqrt_quad_nominal = pdf_2Jets_scaled_slope_sqrt_quad_TF1.Eval(STVal)/pdf_2Jets_scaled_slope_sqrt_quad_TF1_atNorm;
      pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(0, (fitParametersBinned.at("slope_sqrt_quad_fit_slope")).at(nJetsBin) + ((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_error")).at(nJetsBin))*((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_slopeCoefficient")).at(nJetsBin)));
      pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(1, (fitParametersBinned.at("slope_sqrt_quad_fit_sqrt")).at(nJetsBin) + ((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_error")).at(nJetsBin))*((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_sqrtCoefficient")).at(nJetsBin)));
      pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(2, (fitParametersBinned.at("slope_sqrt_quad_fit_quad")).at(nJetsBin) + ((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_error")).at(nJetsBin))*((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_quadCoefficient")).at(nJetsBin)));
      double pdf_slope_sqrt_quad_plus_one_sigma = (pdf_2Jets_scaled_slope_sqrt_quad_TF1.Eval(STVal))/pdf_2Jets_scaled_slope_sqrt_quad_TF1_atNorm;
      pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(0, (fitParametersBinned.at("slope_sqrt_quad_fit_slope")).at(nJetsBin) - ((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_error")).at(nJetsBin))*((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_slopeCoefficient")).at(nJetsBin)));
      pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(1, (fitParametersBinned.at("slope_sqrt_quad_fit_sqrt")).at(nJetsBin) - ((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_error")).at(nJetsBin))*((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_sqrtCoefficient")).at(nJetsBin)));
      pdf_2Jets_scaled_slope_sqrt_quad_TF1.SetParameter(2, (fitParametersBinned.at("slope_sqrt_quad_fit_quad")).at(nJetsBin) - ((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_error")).at(nJetsBin))*((fitParametersBinned.at("slope_sqrt_quad_fit_mode1_quadCoefficient")).at(nJetsBin)));
      double pdf_slope_sqrt_quad_minus_one_sigma = (pdf_2Jets_scaled_slope_sqrt_quad_TF1.Eval(STVal))/pdf_2Jets_scaled_slope_sqrt_quad_TF1_atNorm;
      double pdf_slope_sqrt_quad_higher = std::max({pdf_slope_sqrt_quad_plus_one_sigma, pdf_slope_sqrt_quad_nominal, pdf_slope_sqrt_quad_minus_one_sigma});
      double pdf_slope_sqrt_quad_lower = std::min({pdf_slope_sqrt_quad_plus_one_sigma, pdf_slope_sqrt_quad_nominal, pdf_slope_sqrt_quad_minus_one_sigma});
      double ratio_slope_sqrt_quad_to_unadjusted = pdf_slope_sqrt_quad_nominal/common_denominator;
      double ratio_slope_sqrt_quad_to_unadjusted_higher = pdf_slope_sqrt_quad_higher/common_denominator;
      double ratio_slope_sqrt_quad_to_unadjusted_lower = pdf_slope_sqrt_quad_lower/common_denominator;
      ratioGraph_binned_slope_sqrt_quad_to_unadjusted.SetPoint(STCounter, STVal, ratio_slope_sqrt_quad_to_unadjusted);
      ratioGraph_binned_slope_sqrt_quad_to_unadjusted_high_estimate.SetPoint(STCounter, STVal, normBinCorrectionUp*ratio_slope_sqrt_quad_to_unadjusted_higher);
      ratioGraph_binned_slope_sqrt_quad_to_unadjusted_low_estimate.SetPoint(STCounter, STVal, normBinCorrectionDown*ratio_slope_sqrt_quad_to_unadjusted_lower);
    }

    // plot shape ratios
    TMultiGraph binned_shape_ratios_multigraph = TMultiGraph(("binned_shape_ratios_multigraph_at" + std::to_string(nJetsBin) + "Jets").c_str(), ("Shape ratios (binned), " + std::to_string(nJetsBin) + " Jets bin").c_str());
    TCanvas binned_shape_ratios_canvas = TCanvas(("c_binnedShapeRatios_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("c_binnedShapeRatios_" + std::to_string(nJetsBin) + "JetsBin").c_str(), 2560, 1440);
    TLegend legend_binned_shape_ratios_multigraph = TLegend(0.1, 0.6, 0.4, 0.9);

    ratioGraph_binned_nJetsDistribution_to_unadjusted.SetLineColor(static_cast<EColor>(kBlack)); ratioGraph_binned_nJetsDistribution_to_unadjusted.SetDrawOption("P"); binned_shape_ratios_multigraph.Add(&ratioGraph_binned_nJetsDistribution_to_unadjusted);
    TLegendEntry *legendEntry_binned_nJetsDistribution_to_unadjusted = legend_binned_shape_ratios_multigraph.AddEntry(&ratioGraph_binned_nJetsDistribution_to_unadjusted, (std::to_string(nJetsBin) + " jets distribution / 2 jets kernel").c_str());
    legendEntry_binned_nJetsDistribution_to_unadjusted->SetMarkerColor(static_cast<EColor>(kBlack)); legendEntry_binned_nJetsDistribution_to_unadjusted->SetLineColor(static_cast<EColor>(kBlack)); legendEntry_binned_nJetsDistribution_to_unadjusted->SetTextColor(static_cast<EColor>(kBlack));

    ratioGraph_binned_slope_to_unadjusted.SetLineColor(static_cast<EColor>(kRed+1)); ratioGraph_binned_slope_to_unadjusted.SetDrawOption("C"); binned_shape_ratios_multigraph.Add(&ratioGraph_binned_slope_to_unadjusted);
    TLegendEntry *legendEntry_binned_slope_to_unadjusted = legend_binned_shape_ratios_multigraph.AddEntry(&ratioGraph_binned_slope_to_unadjusted, "linear adjustment");
    legendEntry_binned_slope_to_unadjusted->SetMarkerColor(static_cast<EColor>(kRed+1)); legendEntry_binned_slope_to_unadjusted->SetLineColor(static_cast<EColor>(kRed+1)); legendEntry_binned_slope_to_unadjusted->SetTextColor(static_cast<EColor>(kRed+1));
    ratioGraph_binned_slope_to_unadjusted_high_estimate.SetLineColor(static_cast<EColor>(kRed+1)); ratioGraph_binned_slope_to_unadjusted_high_estimate.SetLineStyle(kDashed); ratioGraph_binned_slope_to_unadjusted_high_estimate.SetDrawOption("C"); binned_shape_ratios_multigraph.Add(&ratioGraph_binned_slope_to_unadjusted_high_estimate);
    ratioGraph_binned_slope_to_unadjusted_low_estimate.SetLineColor(static_cast<EColor>(kRed+1)); ratioGraph_binned_slope_to_unadjusted_low_estimate.SetLineStyle(kDashed); ratioGraph_binned_slope_to_unadjusted_low_estimate.SetDrawOption("C"); binned_shape_ratios_multigraph.Add(&ratioGraph_binned_slope_to_unadjusted_low_estimate);

    ratioGraph_binned_sqrt_to_unadjusted.SetLineColor(static_cast<EColor>(kGreen+3)); ratioGraph_binned_sqrt_to_unadjusted.SetDrawOption("C"); binned_shape_ratios_multigraph.Add(&ratioGraph_binned_sqrt_to_unadjusted);
    TLegendEntry *legendEntry_binned_sqrt_to_unadjusted = legend_binned_shape_ratios_multigraph.AddEntry(&ratioGraph_binned_sqrt_to_unadjusted, "sqrt adjustment");
    legendEntry_binned_sqrt_to_unadjusted->SetMarkerColor(static_cast<EColor>(kGreen+3)); legendEntry_binned_sqrt_to_unadjusted->SetLineColor(static_cast<EColor>(kGreen+3)); legendEntry_binned_sqrt_to_unadjusted->SetTextColor(static_cast<EColor>(kGreen+3));
    ratioGraph_binned_sqrt_to_unadjusted_high_estimate.SetLineColor(static_cast<EColor>(kGreen+3)); ratioGraph_binned_sqrt_to_unadjusted_high_estimate.SetLineStyle(kDashed); ratioGraph_binned_sqrt_to_unadjusted_high_estimate.SetDrawOption("C"); binned_shape_ratios_multigraph.Add(&ratioGraph_binned_sqrt_to_unadjusted_high_estimate);
    ratioGraph_binned_sqrt_to_unadjusted_low_estimate.SetLineColor(static_cast<EColor>(kGreen+3)); ratioGraph_binned_sqrt_to_unadjusted_low_estimate.SetLineStyle(kDashed); ratioGraph_binned_sqrt_to_unadjusted_low_estimate.SetDrawOption("C"); binned_shape_ratios_multigraph.Add(&ratioGraph_binned_sqrt_to_unadjusted_low_estimate);

    ratioGraph_binned_slope_sqrt_to_unadjusted.SetLineColor(static_cast<EColor>(kViolet)); ratioGraph_binned_slope_sqrt_to_unadjusted.SetDrawOption("C"); binned_shape_ratios_multigraph.Add(&ratioGraph_binned_slope_sqrt_to_unadjusted);
    TLegendEntry *legendEntry_binned_slope_sqrt_to_unadjusted = legend_binned_shape_ratios_multigraph.AddEntry(&ratioGraph_binned_slope_sqrt_to_unadjusted, "(linear+sqrt) adjustment");
    legendEntry_binned_slope_sqrt_to_unadjusted->SetMarkerColor(static_cast<EColor>(kViolet)); legendEntry_binned_slope_sqrt_to_unadjusted->SetLineColor(static_cast<EColor>(kViolet)); legendEntry_binned_slope_sqrt_to_unadjusted->SetTextColor(static_cast<EColor>(kViolet));
    ratioGraph_binned_slope_sqrt_to_unadjusted_high_estimate.SetLineColor(static_cast<EColor>(kViolet)); ratioGraph_binned_slope_sqrt_to_unadjusted_high_estimate.SetLineStyle(kDashed); ratioGraph_binned_slope_sqrt_to_unadjusted_high_estimate.SetDrawOption("C"); binned_shape_ratios_multigraph.Add(&ratioGraph_binned_slope_sqrt_to_unadjusted_high_estimate);
    ratioGraph_binned_slope_sqrt_to_unadjusted_low_estimate.SetLineColor(static_cast<EColor>(kViolet)); ratioGraph_binned_slope_sqrt_to_unadjusted_low_estimate.SetLineStyle(kDashed); ratioGraph_binned_slope_sqrt_to_unadjusted_low_estimate.SetDrawOption("C"); binned_shape_ratios_multigraph.Add(&ratioGraph_binned_slope_sqrt_to_unadjusted_low_estimate);

    ratioGraph_binned_slope_sqrt_quad_to_unadjusted.SetLineColor(static_cast<EColor>(kYellow+2)); ratioGraph_binned_slope_sqrt_quad_to_unadjusted.SetDrawOption("C"); binned_shape_ratios_multigraph.Add(&ratioGraph_binned_slope_sqrt_quad_to_unadjusted);
    TLegendEntry *legendEntry_binned_slope_sqrt_quad_to_unadjusted = legend_binned_shape_ratios_multigraph.AddEntry(&ratioGraph_binned_slope_sqrt_quad_to_unadjusted, "(linear+sqrt+quad) adjustment");
    legendEntry_binned_slope_sqrt_quad_to_unadjusted->SetMarkerColor(static_cast<EColor>(kYellow+2)); legendEntry_binned_slope_sqrt_quad_to_unadjusted->SetLineColor(static_cast<EColor>(kYellow+2)); legendEntry_binned_slope_sqrt_quad_to_unadjusted->SetTextColor(static_cast<EColor>(kYellow+2));
    ratioGraph_binned_slope_sqrt_quad_to_unadjusted_high_estimate.SetLineColor(static_cast<EColor>(kYellow+2)); ratioGraph_binned_slope_sqrt_quad_to_unadjusted_high_estimate.SetLineStyle(kDashed); ratioGraph_binned_slope_sqrt_quad_to_unadjusted_high_estimate.SetDrawOption("C"); binned_shape_ratios_multigraph.Add(&ratioGraph_binned_slope_sqrt_quad_to_unadjusted_high_estimate);
    ratioGraph_binned_slope_sqrt_quad_to_unadjusted_low_estimate.SetLineColor(static_cast<EColor>(kYellow+2)); ratioGraph_binned_slope_sqrt_quad_to_unadjusted_low_estimate.SetLineStyle(kDashed); ratioGraph_binned_slope_sqrt_quad_to_unadjusted_low_estimate.SetDrawOption("C"); binned_shape_ratios_multigraph.Add(&ratioGraph_binned_slope_sqrt_quad_to_unadjusted_low_estimate);

    binned_shape_ratios_multigraph.Draw("A");
    legend_binned_shape_ratios_multigraph.SetFillStyle(0);
    legend_binned_shape_ratios_multigraph.Draw();
    binned_shape_ratios_multigraph.GetXaxis()->SetTitle("ST (GeV)");
    binned_shape_ratios_multigraph.GetYaxis()->SetTitle("ratio");
    binned_shape_ratios_canvas.Update();
    binned_shape_ratios_canvas.SaveAs((options.outputFolder + "/binned_shapeRatios_" + std::to_string(nJetsBin) + "JetsBin_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".pdf").c_str());

    // double chi2_unadjusted = getDensityHistogramChiSquareWRTFunction(dataHist, integralsOverBins_unadjusted);
    // double chi2_slopeOnly = getDensityHistogramChiSquareWRTFunction(dataHist, integralsOverBins_slopeOnly);
    // double chi2_sqrtOnly = getDensityHistogramChiSquareWRTFunction(dataHist, integralsOverBins_sqrtOnly);
    // double chi2_combined = getDensityHistogramChiSquareWRTFunction(dataHist, integralsOverBins_combined);

    if (!(options.readParametersFromFiles)) {
      std::cout << "Getting p-values using chi2 values from binned fits..." << std::endl;
      // sanity checks
      assert((((fit_qualities_binned.at("unadjusted")).at(nJetsBin)).ndf) == (1 + ((fit_qualities_binned.at("slope")).at(nJetsBin)).ndf));
      assert((((fit_qualities_binned.at("unadjusted")).at(nJetsBin)).ndf) == (1 + ((fit_qualities_binned.at("sqrt")).at(nJetsBin)).ndf));
      assert((((fit_qualities_binned.at("unadjusted")).at(nJetsBin)).ndf) == (2 + ((fit_qualities_binned.at("slope_sqrt")).at(nJetsBin)).ndf));
      assert((((fit_qualities_binned.at("unadjusted")).at(nJetsBin)).ndf) == (3 + ((fit_qualities_binned.at("slope_sqrt_quad")).at(nJetsBin)).ndf));

      ftest_pValues["unadjusted_vs_slope"][nJetsBin] = get_fTest_prob(((fit_qualities_binned.at("unadjusted")).at(nJetsBin)).chi2, ((fit_qualities_binned.at("slope")).at(nJetsBin)).chi2, ((fit_qualities_binned.at("unadjusted")).at(nJetsBin)).ndf, ((fit_qualities_binned.at("slope")).at(nJetsBin)).ndf);
      ftest_pValues["slope_vs_slope_sqrt"][nJetsBin] = get_fTest_prob(((fit_qualities_binned.at("slope")).at(nJetsBin)).chi2, ((fit_qualities_binned.at("slope_sqrt")).at(nJetsBin)).chi2, ((fit_qualities_binned.at("slope")).at(nJetsBin)).ndf, ((fit_qualities_binned.at("slope_sqrt")).at(nJetsBin)).ndf);
      ftest_pValues["unadjusted_vs_sqrt"][nJetsBin] = get_fTest_prob(((fit_qualities_binned.at("unadjusted")).at(nJetsBin)).chi2, ((fit_qualities_binned.at("sqrt")).at(nJetsBin)).chi2, ((fit_qualities_binned.at("unadjusted")).at(nJetsBin)).ndf, ((fit_qualities_binned.at("sqrt")).at(nJetsBin)).ndf);
      ftest_pValues["sqrt_vs_slope_sqrt"][nJetsBin] = get_fTest_prob(((fit_qualities_binned.at("sqrt")).at(nJetsBin)).chi2, ((fit_qualities_binned.at("slope_sqrt")).at(nJetsBin)).chi2, ((fit_qualities_binned.at("sqrt")).at(nJetsBin)).ndf, ((fit_qualities_binned.at("slope_sqrt")).at(nJetsBin)).ndf);
      ftest_pValues["slope_sqrt_vs_slope_sqrt_quad"][nJetsBin] = get_fTest_prob(((fit_qualities_binned.at("slope_sqrt")).at(nJetsBin)).chi2, ((fit_qualities_binned.at("slope_sqrt_quad")).at(nJetsBin)).chi2, ((fit_qualities_binned.at("slope_sqrt")).at(nJetsBin)).ndf, ((fit_qualities_binned.at("slope_sqrt_quad")).at(nJetsBin)).ndf);
    }

    printSeparator();
  }

  // // Print f-test prob values in a LaTeX-formatted table
  // std::cout << "f-test prob values (formatted):" << std::endl;
  // std::cout << "\\begin{tabular}{|l|c|c|c|}" << std::endl;
  // std::cout << "  \\hline" << std::endl;
  // std::cout << "  f-prob & const\\_vs\\_lin & lin\\_vs\\_quad & fixed\\_lin\\_vs\\_lin \\\\ \\hline" << std::endl;
  // std::cout << std::fixed << std::setprecision(3) << "  nJets = 3 & " << (fTestProbValues.at("const_vs_lin")).at(3) << " & " << (fTestProbValues.at("lin_vs_quad")).at(3) << " & " << (fTestProbValues.at("constrained_lin_vs_lin")).at(3) << " \\\\ \\hline" << std::endl;
  // std::cout << std::fixed << std::setprecision(3) << "  nJets = 4 & " << (fTestProbValues.at("const_vs_lin")).at(4) << " & " << (fTestProbValues.at("lin_vs_quad")).at(4) << " & " << (fTestProbValues.at("constrained_lin_vs_lin")).at(4) << " \\\\ \\hline" << std::endl;
  // std::cout << std::fixed << std::setprecision(3) << "  nJets = 5 & " << (fTestProbValues.at("const_vs_lin")).at(5) << " & " << (fTestProbValues.at("lin_vs_quad")).at(5) << " & " << (fTestProbValues.at("constrained_lin_vs_lin")).at(5) << " \\\\ \\hline" << std::endl;
  // std::cout << std::fixed << std::setprecision(3) << "  nJets $\\geq$ 6 & " << (fTestProbValues.at("const_vs_lin")).at(6) << " & " << (fTestProbValues.at("lin_vs_quad")).at(6) << " & " << (fTestProbValues.at("constrained_lin_vs_lin")).at(6) << " \\\\ \\hline" << std::endl;
  // std::cout << "\\end{tabular}" << std::endl;

  // // Print best fit slopes in LaTeX-formatted table
  // std::cout << "slope values (formatted):" << std::endl;
  // std::cout << "\\begin{tabular}{|l|c|}" << std::endl;
  // std::cout << "  \\hline" << std::endl;
  // std::cout << "  nJets Bin & Best-fit slope \\\\ \\hline" << std::endl;
  // std::cout << std::fixed << std::setprecision(3) << "  nJets = 3 & " << (fitParametersBinned.at("slope")).at(3) << " $\\pm$ " << (fitParametersBinned.at("slopeError")).at(3) << " \\\\ \\hline" << std::endl;
  // std::cout << std::fixed << std::setprecision(3) << "  nJets = 4 & " << (fitParametersBinned.at("slope")).at(4) << " $\\pm$ " << (fitParametersBinned.at("slopeError")).at(4) << " \\\\ \\hline" << std::endl;
  // std::cout << std::fixed << std::setprecision(3) << "  nJets = 5 & " << (fitParametersBinned.at("slope")).at(5) << " $\\pm$ " << (fitParametersBinned.at("slopeError")).at(5) << " \\\\ \\hline" << std::endl;
  // std::cout << std::fixed << std::setprecision(3) << "  nJets $\\geq$ 6 & " << (fitParametersBinned.at("slope")).at(6) << " $\\pm$ " << (fitParametersBinned.at("slopeError")).at(6) << " \\\\ \\hline" << std::endl;
  // std::cout << "\\end{tabular}" << std::endl;

  // Print p-values from nlls in a LaTeX-formatted table
  // std::cout << "p-values from nlls (formatted):" << std::endl;
  // std::cout << "\\begin{tabular}{|l|p{0.18\\textwidth}|p{0.18\\textwidth}|p{0.18\\textwidth}|p{0.18\\textwidth}|}" << std::endl;
  // std::cout << "  \\hline" << std::endl;
  // std::cout << "  p-values & unadjusted vs slope-only & slope-only vs combined & unadjusted vs sqrt-only & sqrt-only vs combined \\\\ \\hline" << std::endl;
  // std::cout << std::fixed << std::setprecision(3) << "  nJets = 3 & " << (nll_pValues.at("unadjusted_vs_slopeOnly")).at(3) << " & " << (nll_pValues.at("slopeOnly_vs_combined")).at(3) << " & " << (nll_pValues.at("unadjusted_vs_sqrtOnly")).at(3) << " & " << (nll_pValues.at("sqrtOnly_vs_combined")).at(3) << " \\\\ \\hline" << std::endl;
  // std::cout << std::fixed << std::setprecision(3) << "  nJets = 4 & " << (nll_pValues.at("unadjusted_vs_slopeOnly")).at(4) << " & " << (nll_pValues.at("slopeOnly_vs_combined")).at(4) << " & " << (nll_pValues.at("unadjusted_vs_sqrtOnly")).at(4) << " & " << (nll_pValues.at("sqrtOnly_vs_combined")).at(4) << " \\\\ \\hline" << std::endl;
  // std::cout << std::fixed << std::setprecision(3) << "  nJets = 5 & " << (nll_pValues.at("unadjusted_vs_slopeOnly")).at(5) << " & " << (nll_pValues.at("slopeOnly_vs_combined")).at(5) << " & " << (nll_pValues.at("unadjusted_vs_sqrtOnly")).at(5) << " & " << (nll_pValues.at("sqrtOnly_vs_combined")).at(5) << " \\\\ \\hline" << std::endl;
  // std::cout << std::fixed << std::setprecision(3) << "  nJets $\\geq$ 6 & " << (nll_pValues.at("unadjusted_vs_slopeOnly")).at(6) << " & " << (nll_pValues.at("slopeOnly_vs_combined")).at(6) << " & " << (nll_pValues.at("unadjusted_vs_sqrtOnly")).at(6) << " & " << (nll_pValues.at("sqrtOnly_vs_combined")).at(6) << " \\\\ \\hline" << std::endl;
  // std::cout << "\\end{tabular}" << std::endl;

  // Print ftest pvalues from binned chi2 fits in a LaTeX-formatted table
  if (!(options.readParametersFromFiles)) {
    std::cout << "p-values from binned chi2 fits:" << std::endl;
    std::cout << "\\begin{tabular}{|p{0.14\\textwidth}|p{0.14\\textwidth}|p{0.14\\textwidth}|p{0.14\\textwidth}|p{0.14\\textwidth}|p{0.14\\textwidth}|}" << std::endl;
    std::cout << "  \\hline" << std::endl;
    std::cout << "  p-values & unadjusted \\newline vs \\newline linear & linear \\newline vs \\newline (linear+sqrt) & unadjusted \\newline vs \\newline sqrt & sqrt \\newline vs \\newline (linear+sqrt) & (linear+sqrt) \\newline vs \\newline (linear+sqrt \\newline +quad) \\\\ \\hline" << std::endl;
    std::cout << std::fixed << std::setprecision(3) << "  nJets = 3 & " << (ftest_pValues.at("unadjusted_vs_slope")).at(3) << " & " << (ftest_pValues.at("slope_vs_slope_sqrt")).at(3) << " & " << (ftest_pValues.at("unadjusted_vs_sqrt")).at(3) << " & " << (ftest_pValues.at("sqrt_vs_slope_sqrt")).at(3) << " & " << (ftest_pValues.at("slope_sqrt_vs_slope_sqrt_quad")).at(3) << " \\\\ \\hline" << std::endl;
    std::cout << std::fixed << std::setprecision(3) << "  nJets = 4 & " << (ftest_pValues.at("unadjusted_vs_slope")).at(4) << " & " << (ftest_pValues.at("slope_vs_slope_sqrt")).at(4) << " & " << (ftest_pValues.at("unadjusted_vs_sqrt")).at(4) << " & " << (ftest_pValues.at("sqrt_vs_slope_sqrt")).at(4) << " & " << (ftest_pValues.at("slope_sqrt_vs_slope_sqrt_quad")).at(4) << " \\\\ \\hline" << std::endl;
    std::cout << std::fixed << std::setprecision(3) << "  nJets = 5 & " << (ftest_pValues.at("unadjusted_vs_slope")).at(5) << " & " << (ftest_pValues.at("slope_vs_slope_sqrt")).at(5) << " & " << (ftest_pValues.at("unadjusted_vs_sqrt")).at(5) << " & " << (ftest_pValues.at("sqrt_vs_slope_sqrt")).at(5) << " & " << (ftest_pValues.at("slope_sqrt_vs_slope_sqrt_quad")).at(5) << " \\\\ \\hline" << std::endl;
    std::cout << std::fixed << std::setprecision(3) << "  nJets $\\geq$ 6 & " << (ftest_pValues.at("unadjusted_vs_slope")).at(6) << " & " << (ftest_pValues.at("slope_vs_slope_sqrt")).at(6) << " & " << (ftest_pValues.at("unadjusted_vs_sqrt")).at(6) << " & " << (ftest_pValues.at("sqrt_vs_slope_sqrt")).at(6) << " & " << (ftest_pValues.at("slope_sqrt_vs_slope_sqrt_quad")).at(6) << " \\\\ \\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
  }

  // Print best fit slopes from unbinned fit in LaTeX-formatted table
  // std::cout << "slope values from unbinned fit (formatted):" << std::endl;
  // std::cout << "\\begin{tabular}{|l|c|}" << std::endl;
  // std::cout << "  \\hline" << std::endl;
  // std::cout << "  nJets Bin & Best-fit slope \\\\ \\hline" << std::endl;
  // std::cout << std::fixed << std::setprecision(3) << "  nJets = 3 & " << (fitParametersUnbinned.at("slope_fit_slope")).at(3) << " $\\pm$ " << (fitParametersUnbinned.at("slope_fit_slopeError")).at(3) << " \\\\ \\hline" << std::endl;
  // std::cout << std::fixed << std::setprecision(3) << "  nJets = 4 & " << (fitParametersUnbinned.at("slope_fit_slope")).at(4) << " $\\pm$ " << (fitParametersUnbinned.at("slope_fit_slopeError")).at(4) << " \\\\ \\hline" << std::endl;
  // std::cout << std::fixed << std::setprecision(3) << "  nJets = 5 & " << (fitParametersUnbinned.at("slope_fit_slope")).at(5) << " $\\pm$ " << (fitParametersUnbinned.at("slope_fit_slopeError")).at(5) << " \\\\ \\hline" << std::endl;
  // std::cout << std::fixed << std::setprecision(3) << "  nJets $\\geq$ 6 & " << (fitParametersUnbinned.at("slope_fit_slope")).at(6) << " $\\pm$ " << (fitParametersUnbinned.at("slope_fit_slopeError")).at(6) << " \\\\ \\hline" << std::endl;
  // std::cout << "\\end{tabular}" << std::endl;

  // // write parameters for binned fit
  // std::ofstream fitParametersBinnedFile((options.outputFolder + "/fitParameters_binnedFit_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".dat").c_str());
  // assert(fitParametersBinnedFile.is_open());
  // for (int fitParametersBinnedListIndex = 0; fitParametersBinnedListIndex < static_cast<int>(fitParametersBinnedList.size()); ++fitParametersBinnedListIndex) {
  //   fitParametersBinnedFile << fitParametersBinnedList.at(fitParametersBinnedListIndex) << std::endl;
  // }
  // fitParametersBinnedFile.close();
  // std::cout << "Binned fit parameters written to file: " << (options.outputFolder + "/fitParameters_binnedFit_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".dat") << std::endl;

  // // write parameters for unbinned fit
  if (!(options.readParametersFromFiles)) {
    std::ofstream fitParametersUnbinnedFile((options.outputFolder + "/unbinned_fitParameters_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".dat").c_str());
    assert(fitParametersUnbinnedFile.is_open());
    for (int fitParametersUnbinnedListIndex = 0; fitParametersUnbinnedListIndex < static_cast<int>(fitParametersUnbinnedList.size()); ++fitParametersUnbinnedListIndex) {
      fitParametersUnbinnedFile << fitParametersUnbinnedList.at(fitParametersUnbinnedListIndex) << std::endl;
    }
    fitParametersUnbinnedFile.close();
    std::cout << "Unbinned fit parameters written to file: " << (options.outputFolder + "/unbinned_fitParameters_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".dat") << std::endl;
    std::ofstream fitParametersBinnedFile((options.outputFolder + "/binned_fitParameters_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".dat").c_str());
    assert(fitParametersBinnedFile.is_open());
    for (int fitParametersBinnedListIndex = 0; fitParametersBinnedListIndex < static_cast<int>(fitParametersBinnedList.size()); ++fitParametersBinnedListIndex) {
      fitParametersBinnedFile << fitParametersBinnedList.at(fitParametersBinnedListIndex) << std::endl;
    }
    fitParametersBinnedFile.close();
    std::cout << "Binned fit parameters written to file: " << (options.outputFolder + "/binned_fitParameters_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".dat") << std::endl;
  }

  std::cout << "All done!" << std::endl;
  return EXIT_SUCCESS;
}
