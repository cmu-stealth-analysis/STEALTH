#include "../include/optimizeNormBin.h"

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

double parseLineForFloatWithCheck(const std::string& inputLine, const std::string& targetName) {
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
  assert(static_cast<int>(components.size()) == 3);
  assert(components.at(0) == "float");
  assert(components.at(1) == targetName);
  return std::stod(components.at(2));
}

double getDensityHistogramChiSquareWRTFunction(TH1F& inputDensityHistogram, std::map<int, double>& functionIntegralsOverBins) {
  // first get number of entries and total integral for normalization
  double weightedNEntries = 0.;
  double sumFunctionIntegrals = 0.;
  for (int binIndex = 1; binIndex <= inputDensityHistogram.GetXaxis()->GetNbins(); ++binIndex) {
    weightedNEntries += (inputDensityHistogram.GetBinContent(binIndex))*(inputDensityHistogram.GetXaxis()->GetBinWidth(binIndex));
    sumFunctionIntegrals += functionIntegralsOverBins.at(binIndex);
  }
  double functionNormalization = weightedNEntries/sumFunctionIntegrals;

  double chiSquare = 0.;
  for (int binIndex = 1; binIndex <= inputDensityHistogram.GetXaxis()->GetNbins(); ++binIndex) {
    double observationFromDensityHistogram = (inputDensityHistogram.GetBinContent(binIndex))*(inputDensityHistogram.GetXaxis()->GetBinWidth(binIndex));
    double error = (inputDensityHistogram.GetBinError(binIndex))*(inputDensityHistogram.GetXaxis()->GetBinWidth(binIndex));
    double expectation_fromFunction = functionNormalization*functionIntegralsOverBins.at(binIndex);
    // std::cout << "At binIndex = " << binIndex << ", bin center: " << inputDensityHistogram.GetXaxis()->GetBinCenter(binIndex) << ", observationFromDensityHistogram: " << observationFromDensityHistogram << ", functionIntegralsOverBins: " << functionIntegralsOverBins.at(binIndex) << ", expectation_fromFunction: " << expectation_fromFunction << ", error: " << error << std::endl;
    if ((error > 0) && (observationFromDensityHistogram > 0)) {
      chiSquare += (std::pow(observationFromDensityHistogram - expectation_fromFunction, 2))/error;
    }
  }
  // std::cout << "chiSquare = " << chiSquare << std::endl;
  return chiSquare;
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

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  TH1::AddDirectory(kFALSE);
  do_sanity_checks_fitTypes();
  RooMsgService::instance().setGlobalKillBelow(MsgLevel::WARNING);

  tmArgumentParser argumentParser = tmArgumentParser("Run script that prints useful info about the normalization.");
  argumentParser.addArgument("sourceFilePath", "", true, "Path to file containing list of paths with n-tuplized events.");
  argumentParser.addArgument("outputFolder", "", true, "Output folder.");
  argumentParser.addArgument("selection", "", true, "Name of selection: \"singlemedium\", \"signal_loose\", etc.");
  argumentParser.addArgument("identifier", "", true, "Identifier: \"MC_GJet17\", \"MC_GJet\", etc.");
  argumentParser.addArgument("STBoundariesSourceFile", "STRegionBoundaries_normOptimization.dat", false, "Identifier: \"MC_GJet17\", \"MC_GJet\", etc.");
  argumentParser.addArgument("yearString", "all", false, "String with year: can take values \"2016\", \"2017\", \"2018\", or \"all\".");
  argumentParser.addArgument("STNormTarget", "1250.0", false, "Target value to use for ST norm.");
  argumentParser.addArgument("STNormMax", "1450.0", false, "Bins with a central value larger than this value are not considered.");
  argumentParser.addArgument("PDF_nSTBins", "25", false, "Number of bins for plotting datasets.");
  argumentParser.addArgument("PDF_STMin", "1000.0", false, "Min of range for plotting datasets.");
  argumentParser.addArgument("PDF_STMax", "3500.0", false, "Max of range for plotting datasets.");
  argumentParser.addArgument("readParametersFromFile", "/dev/null", false, "If this argument is set, then no unbinned fits are performed; instead, the data is read in from the file location given as the value of this argument.");
  argumentParser.setPassedStringValues(argc, argv);
  optionsStruct options = getOptionsFromParser(argumentParser);
  std::cout << "Options passed:" << std::endl << options << std::endl;

  std::string constrainedLinFunction = "1.0+[0]*((x/" + std::to_string(options.STNormTarget) + ") - 1.0)";
  std::map<fitType, std::string> fitFunctions = {
    {fitType::fitConst, "pol0"},
    {fitType::fitLin, "pol1"},
    {fitType::fitQuad, "pol2"},
    {fitType::fitConstrainedLin, constrainedLinFunction}
  };
  assert(static_cast<int>(fitFunctions.size()) == static_cast<int>(fitType::nFitTypes));

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
  RooRealVar rooVar_ST("roo_ST", "roo_ST", options.PDF_STMin, options.PDF_STMax, "GeV");
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
  // idiotic, I know, but the following compiles but results in a segfault that I've been unable to debug:
  // std::map<int, RooDataSet> STDataSets;
  // for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
  //   STDataSets[nJetsBin] = RooDataSet(("STDataSet_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("STDataSet_" + std::to_string(nJetsBin) + "JetsBin").c_str(), RooArgSet(rooVar_ST, rooVar_weight), WeightVar(rooVar_weight));
  // }

  TFile *sourceFile = TFile::Open((options.sourceFilePath).c_str(), "READ");
  if ((!(sourceFile->IsOpen())) or (sourceFile->IsZombie())) {
    std::cout << "ERROR: unable to open file with name: " << options.sourceFilePath << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::cout << "Fetching ST datasets and histograms from file: " << options.sourceFilePath << std::endl;
  std::map<int, TH1F*> STDistributions;
  for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
    STDistributions[nJetsBin] = new TH1F();
    std::string inputHistogramName = "h_ST_" + std::to_string(nJetsBin) + "JetsBin";
    std::cout << "Getting histogram with name: " << inputHistogramName << std::endl;
    sourceFile->GetObject(inputHistogramName.c_str(), STDistributions[nJetsBin]);
    assert(STDistributions[nJetsBin] != nullptr);
    STDistributions[nJetsBin]->SetName(inputHistogramName.c_str());

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
      if (((*ST) < options.PDF_STMin) || ((*ST) > options.PDF_STMax)) continue;
      (STDataSets.at(nJetsBin))->add(RooArgSet(rooVar_ST), eventWeight);
    }
    (STDataSets.at(nJetsBin))->Print();
  }

  // Following isn't necessary?
  // std::map<int, TH1F*> ratioHistograms;
  // std::map<int, std::map<int, bool> > isZeroBin;
  // for (int nJetsBin = 3; nJetsBin <= 6; ++nJetsBin) {
  //   ratioHistograms[nJetsBin] = new TH1F(("h_ST_ratio_" + std::to_string(nJetsBin) + "JetsTo2Jets").c_str(), ("ST distribution: " + std::to_string(nJetsBin) + " Jets/2 Jets;ST").c_str(), (options.STRegions.STBoundaries.size()-1), &(options.STRegions.STBoundaries.at(0)));
  //   for (int binCounter = 1; binCounter <= (ratioHistograms.at(nJetsBin))->GetXaxis()->GetNbins(); ++binCounter) {
  //     double numerator = (STDistributions.at(nJetsBin))->GetBinContent(binCounter);
  //     double denominator = (STDistributions.at(2))->GetBinContent(binCounter);
  //     double ratio = 1.0;
  //     double ratioError = 1.0;
  //     if ((denominator > 0) && (numerator > 0)) {
  //       double numeratorError = (STDistributions.at(nJetsBin))->GetBinError(binCounter);
  //       double denominatorError = (STDistributions.at(2))->GetBinError(binCounter);
  //       ratio = numerator/denominator;
  //       ratioError = ratio*std::sqrt(std::pow(numeratorError/numerator, 2) + std::pow(denominatorError/denominator, 2));
  //       isZeroBin[nJetsBin][binCounter] = false;
  //     }
  //     else {
  //       isZeroBin[nJetsBin][binCounter] = true;
  //     }
  //     ratioHistograms.at(nJetsBin)->SetBinContent(binCounter, ratio);
  //     ratioHistograms.at(nJetsBin)->SetBinError(binCounter, ratioError);
  //   }
  // }

  int normBinIndex = 5;
  int STNormTargetBin = STDistributions.at(2)->GetXaxis()->FindFixBin(options.STNormTarget);

  std::map<fitType, std::map<int, TGraph> > chiSqPerNDFGraphs;
  for (int fitTypeIndex = fitTypeFirst; fitTypeIndex < static_cast<int>(fitType::nFitTypes); ++fitTypeIndex) {
    fitType fit_type = static_cast<fitType>(fitTypeIndex);
    for (int nJetsBin = 3; nJetsBin <= 6; ++nJetsBin) {
      chiSqPerNDFGraphs[fit_type][nJetsBin] = TGraph();
      chiSqPerNDFGraphs[fit_type][nJetsBin].SetName(("chiSqPerNDFs_fitType" + fitTypeNames.at(fit_type) + "_" + std::to_string(nJetsBin) + "JetsTo2Jets").c_str());
    }
  }

  std::vector<std::string> fitParametersBinnedList;
  std::map<std::string, std::map<int, double> > fitParametersBinned;
  std::map<std::string, std::map<int, float> > fTestProbValues;

  while (true) {
    // TODO fix normBinIndex, off by 1?
    
    double STNorm = STDistributions.at(2)->GetXaxis()->GetBinCenter(normBinIndex);
    bool saveRatios = (STDistributions.at(2)->GetXaxis()->FindFixBin(STNorm) == STNormTargetBin);
    if (STNorm > options.STNormMax) break;
    std::cout << "Trying ST norm: " << std::to_string(STNorm) << std::endl;
    std::cout << "saveRatios: " << (saveRatios? "true" : "false") << std::endl;

    std::vector<double> st_boundaries_copy = getVectorCopyStartingFromIndex(normBinIndex, options.STRegions.STBoundaries);
    TH1F histogramCopy_2Jets(("h_ST_2JetsBin_normBin" + std::to_string(normBinIndex)).c_str(), "ST distribution: 2 Jets;ST", (st_boundaries_copy.size()-1), &(st_boundaries_copy.at(0)));
    histogramCopy_2Jets.Sumw2();

    for (int binCounter=1; binCounter <= histogramCopy_2Jets.GetXaxis()->GetNbins(); ++binCounter) {
      double binCenter = histogramCopy_2Jets.GetXaxis()->GetBinCenter(binCounter);
      int original_binIndex = STDistributions.at(2)->FindFixBin(binCenter);
      histogramCopy_2Jets.SetBinContent(binCounter, STDistributions.at(2)->GetBinContent(original_binIndex));
      histogramCopy_2Jets.SetBinError(binCounter, STDistributions.at(2)->GetBinError(original_binIndex));
    }

    for (int nJetsBin = 3; nJetsBin <= 6; ++nJetsBin) {
      TH1F histogramCopy(("h_ST_" + std::to_string(nJetsBin) + "JetsBin_normBin" + std::to_string(normBinIndex)).c_str(), ("ST distribution: " + std::to_string(nJetsBin) + " Jets;ST").c_str(), (st_boundaries_copy.size()-1), &(st_boundaries_copy.at(0)));
      histogramCopy.Sumw2();
      TGraphErrors ratioGraph = TGraphErrors();
      ratioGraph.SetName(("ratioGraph_" + std::to_string(nJetsBin) + "JetsTo2Jets_normBin" + std::to_string(normBinIndex)).c_str());
      std::vector<TGraphErrorsPointStruct> ratioGraphPoints;
      for (int binCounter=1; binCounter <= histogramCopy.GetXaxis()->GetNbins(); ++binCounter) {
        double binCenter = histogramCopy.GetXaxis()->GetBinCenter(binCounter);
        int original_binIndex = STDistributions.at(nJetsBin)->FindFixBin(binCenter);
        histogramCopy.SetBinContent(binCounter, STDistributions.at(nJetsBin)->GetBinContent(original_binIndex));
        histogramCopy.SetBinError(binCounter, STDistributions.at(nJetsBin)->GetBinError(original_binIndex));
        double numerator = (STDistributions.at(nJetsBin)->GetBinContent(original_binIndex))/(STDistributions.at(nJetsBin)->GetBinContent(normBinIndex));
        double denominator = (STDistributions.at(2)->GetBinContent(original_binIndex))/(STDistributions.at(2)->GetBinContent(normBinIndex));
        if ((numerator > 0) && (denominator > 0)) {
          double fracError_numerator_thisBin = (STDistributions.at(nJetsBin)->GetBinError(original_binIndex))/(STDistributions.at(nJetsBin)->GetBinContent(original_binIndex));
          double fracError_numerator_normBin = (STDistributions.at(nJetsBin)->GetBinError(normBinIndex))/(STDistributions.at(nJetsBin)->GetBinContent(normBinIndex));

          double fracError_denominator_thisBin = (STDistributions.at(2)->GetBinError(original_binIndex))/(STDistributions.at(2)->GetBinContent(original_binIndex));
          double fracError_denominator_normBin = (STDistributions.at(2)->GetBinError(normBinIndex))/(STDistributions.at(2)->GetBinContent(normBinIndex));

          double ratio = numerator/denominator;
          double ratioError = ratio*std::sqrt(std::pow(fracError_numerator_thisBin, 2) + std::pow(fracError_numerator_normBin, 2) + std::pow(fracError_denominator_thisBin, 2) + std::pow(fracError_denominator_normBin, 2));
          int ratioGraphIndex = ratioGraph.GetN();
          ratioGraph.SetPoint(ratioGraphIndex, binCenter, ratio);
          ratioGraph.SetPointError(ratioGraphIndex, 0.5*(histogramCopy.GetXaxis()->GetBinUpEdge(binCounter) - histogramCopy.GetXaxis()->GetBinLowEdge(binCounter)), ratioError);
          ratioGraphPoints.push_back(TGraphErrorsPointStruct(binCenter, ratio, 0.5*(histogramCopy.GetXaxis()->GetBinUpEdge(binCounter) - histogramCopy.GetXaxis()->GetBinLowEdge(binCounter)), ratioError));
        }
      }
      std::map<fitType, TF1*> fits;
      std::map<fitType, TFitResultPtr> fitResults;
      for (int fit_type_index = fitTypeFirst; fit_type_index < static_cast<int>(fitType::nFitTypes); ++fit_type_index) {
        fitType fit_type = static_cast<fitType>(fit_type_index);
        fits[fit_type] = new TF1(("fit_type" + fitTypeNames.at(fit_type) + "_" + std::to_string(nJetsBin) + "JetsBin_normBin" + std::to_string(normBinIndex)).c_str(), (fitFunctions.at(fit_type)).c_str(), histogramCopy.GetXaxis()->GetBinLowEdge(2), histogramCopy.GetXaxis()->GetBinUpEdge(histogramCopy.GetXaxis()->GetNbins()));

        if (fit_type == fitType::fitConstrainedLin) {
          (fits.at(fit_type))->SetParName(0, "slope");
          (fits.at(fit_type))->SetParameter(0, 0.);
          (fits.at(fit_type))->SetParLimits(0, -0.5, 5.);
          (fits.at(fit_type))->SetLineColor(colors.at(nJetsBin));
          (fits.at(fit_type))->SetLineStyle(kDashed);
        }

        fitResults[fit_type] = ratioGraph.Fit(fits.at(fit_type), "EX0QREMS+");
        double fitChi2 = fitResults.at(fit_type)->Chi2();
        double fitNDF = fitResults.at(fit_type)->Ndf();
        if (fitNDF > 0.) {
          double chiSqPerNDF = fitChi2/fitNDF;
          ((chiSqPerNDFGraphs.at(fit_type)).at(nJetsBin)).SetPoint(((chiSqPerNDFGraphs.at(fit_type)).at(nJetsBin)).GetN(), STNorm, chiSqPerNDF);
        }
        else {
          std::cout << "WARNING: Zero division error for nJetsBin: " << nJetsBin << ", fitType: " << fitTypeNames.at(fit_type) << ", STNorm: " << STNorm << std::endl;
        }
        if (saveRatios && (fit_type == fitType::fitConstrainedLin)) {
          double bestFitSlope = fitResults.at(fit_type)->Parameter(0);
          fitParametersBinnedList.push_back(std::string("float bestFitSlope_" + std::to_string(nJetsBin) + "Jets=" + std::to_string(bestFitSlope)));
          fitParametersBinned["slope"][nJetsBin] = bestFitSlope;
          double bestFitSlopeError = fitResults.at(fit_type)->ParError(0);
          fitParametersBinnedList.push_back(std::string("float bestFitSlopeError_" + std::to_string(nJetsBin) + "Jets=" + std::to_string(bestFitSlopeError)));
          fitParametersBinned["slopeError"][nJetsBin] = bestFitSlopeError;
        }
      }
      if (saveRatios) {
        for (int fit_type_index = fitTypeFirst; fit_type_index < static_cast<int>(fitType::nFitTypes); ++fit_type_index) {
          fitType fit_type = static_cast<fitType>(fit_type_index);
          TCanvas outputCanvas = TCanvas(("c_residuals_" + fitTypeNames.at(fit_type) + "_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("c_residuals_" + fitTypeNames.at(fit_type) + "_" + std::to_string(nJetsBin) + "JetsBin").c_str(), 1024, 256);
          TGraphErrors residualsGraph = TGraphErrors();
          residualsGraph.SetName(("residuals_" + fitTypeNames.at(fit_type) + "_" + std::to_string(nJetsBin) + "JetsBin").c_str());
          for (int index=0; index < ratioGraph.GetN(); ++index) {
            TGraphErrorsPointStruct& point = ratioGraphPoints.at(index);
            double fitFunctionValue = 1.;
            if (index > 0) fitFunctionValue = (fits.at(fit_type))->Eval(point.x_val);
            double residual = point.y_val - fitFunctionValue;
            residualsGraph.SetPoint(index, point.x_val, residual);
            residualsGraph.SetPointError(index, point.x_err, point.y_err);
          }
          residualsGraph.Draw("AP");
          TLine lineAt0 = TLine(histogramCopy.GetXaxis()->GetBinLowEdge(2), 0., histogramCopy.GetXaxis()->GetBinUpEdge(histogramCopy.GetXaxis()->GetNbins()), 0.);
          lineAt0.SetLineStyle(kSolid);
          lineAt0.Draw();
          residualsGraph.SetTitle(("Residuals for fitType: " + fitTypeNames.at(fit_type)).c_str());
          gPad->Update();
          outputCanvas.SaveAs((options.outputFolder + "/residuals_targetNorm_fitType_" + fitTypeNames.at(fit_type) + "_" + std::to_string(nJetsBin) + "JetsBin_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".pdf").c_str());
        }

        TCanvas outputCanvas = TCanvas(("c_ratios_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("c_ratios_" + std::to_string(nJetsBin) + "JetsBin").c_str(), 2560, 1440);
        ratioGraph.SetLineColor(colors.at(nJetsBin));
        ratioGraph.Draw("AP");
        outputCanvas.SaveAs((options.outputFolder + "/ratios_targetNorm_" + std::to_string(nJetsBin) + "JetsBin_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".pdf").c_str());
        std::cout << "Getting f-test values for nJetsBin = " << nJetsBin << std::endl;
        fTestProbValues["const_vs_lin"][nJetsBin] = get_fTest_prob((fitResults.at(fitType::fitConst))->Chi2(), (fitResults.at(fitType::fitLin))->Chi2(), (fitResults.at(fitType::fitConst))->Ndf(), (fitResults.at(fitType::fitLin))->Ndf());
        fTestProbValues["lin_vs_quad"][nJetsBin] = get_fTest_prob((fitResults.at(fitType::fitLin))->Chi2(), (fitResults.at(fitType::fitQuad))->Chi2(), (fitResults.at(fitType::fitLin))->Ndf(), (fitResults.at(fitType::fitQuad))->Ndf());
        fTestProbValues["constrained_lin_vs_lin"][nJetsBin] = get_fTest_prob((fitResults.at(fitType::fitConstrainedLin))->Chi2(), (fitResults.at(fitType::fitLin))->Chi2(), (fitResults.at(fitType::fitConstrainedLin))->Ndf(), (fitResults.at(fitType::fitLin))->Ndf());
      }
    }
    ++normBinIndex;
  }

  for (int fit_type_index = fitTypeFirst; fit_type_index < static_cast<int>(fitType::nFitTypes); ++fit_type_index) {
    fitType fit_type = static_cast<fitType>(fit_type_index);
    TCanvas outputCanvas = TCanvas(("c_chiSqPerNDFs_fitType" + fitTypeNames.at(fit_type)).c_str(), ("c_chiSqPerNDFs_fitType" + fitTypeNames.at(fit_type)).c_str(), 2560, 1440);
    TMultiGraph chiSqPerNDFsMultigraph = TMultiGraph();
    for (int nJetsBin = 3; nJetsBin <= 6; ++nJetsBin) {
      ((chiSqPerNDFGraphs.at(fit_type)).at(nJetsBin)).SetLineColor(colors.at(nJetsBin));
      ((chiSqPerNDFGraphs.at(fit_type)).at(nJetsBin)).SetMarkerStyle(kCircle);
      ((chiSqPerNDFGraphs.at(fit_type)).at(nJetsBin)).SetMarkerSize(0.75);
      chiSqPerNDFsMultigraph.Add(&((chiSqPerNDFGraphs.at(fit_type)).at(nJetsBin)));
      chiSqPerNDFsMultigraph.Draw("APL");
      chiSqPerNDFsMultigraph.GetXaxis()->SetTitle("ST norm");
      chiSqPerNDFsMultigraph.GetYaxis()->SetTitle("#chi^{2}/NDF");
      chiSqPerNDFsMultigraph.SetTitle(chiSqPerNDFGraphTitles.at(fit_type).c_str());
      outputCanvas.Update();
      outputCanvas.SaveAs((options.outputFolder + "/chiSqPerNDFs_fitType_" + fitTypeNames.at(fit_type) + "_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".pdf").c_str());
    }
  }
    
  TCanvas outputCanvas = TCanvas("c_STDistributions", "c_STDistributions", 2560, 1440);
  TLegend legend = TLegend(0.7, 0.6, 0.9, 0.9);
  STDistributions.at(2)->SetLineColor(colors.at(2));
  STDistributions.at(2)->Draw();
  TLegendEntry *legendEntry = legend.AddEntry(STDistributions.at(2), "nJets = 2");
  legendEntry->SetMarkerColor(colors.at(2));
  legendEntry->SetLineColor(colors.at(2));
  legendEntry->SetTextColor(colors.at(2));
  for (int nJetsBin = 3; nJetsBin <= 6; ++nJetsBin) {
    (STDistributions.at(nJetsBin))->Scale((STDistributions.at(2))->GetBinContent((STDistributions.at(2))->FindFixBin(options.STNormTarget))/(STDistributions.at(nJetsBin))->GetBinContent((STDistributions.at(nJetsBin))->FindFixBin(options.STNormTarget)));
    (STDistributions.at(nJetsBin))->SetLineColor(colors.at(nJetsBin));
    (STDistributions.at(nJetsBin))->Draw("SAME");
    std::string legendText = "nJets = " + std::to_string(nJetsBin);
    if (nJetsBin == 6) legendText = "nJets #geq 6";
    TLegendEntry *legendEntry_nJetsBin = legend.AddEntry((STDistributions.at(nJetsBin)), legendText.c_str());
    legendEntry_nJetsBin->SetMarkerColor(colors.at(nJetsBin));
    legendEntry_nJetsBin->SetLineColor(colors.at(nJetsBin));
    legendEntry_nJetsBin->SetTextColor(colors.at(nJetsBin));
  }

  legend.Draw();
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  outputCanvas.Update();
  outputCanvas.SaveAs((options.outputFolder + "/STDistributions_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".pdf").c_str());

  // A few useful initializations
  std::map<std::string, std::map<int, double> > fitParametersUnbinned;
  std::vector<std::string> fitParametersUnbinnedList;
  double rooVar_slope_minVal = -1.0/(((options.PDF_STMax)/(options.STNormTarget)) - 1.0);
  double rooVar_sqrt_minVal = -1.0/(std::sqrt((options.PDF_STMax)/(options.STNormTarget)) - 1.0);
  std::map<std::string, std::map<int, double> > nll_pValues;
  std::map<std::string, std::map<int, double> > ftest_pValues_unbinnedFit;
  // double effective_STRange = options.PDF_STMax - options.STRegions.STAxis.GetBinLowEdge(options.STRegions.STAxis.FindFixBin(options.STNormTarget));
  // double effective_STWidth = (options.PDF_STMax-options.PDF_STMin)/(options.PDF_nSTBins);
  // int effective_nbins_for_histogram = static_cast<int>(0.5 + effective_STRange/effective_STWidth); // -1 because of the overall normalization;
  std::vector<double> st_boundaries_for_unbinned_chi2;
  int st_boundaries_for_unbinned_chi2_runningIndex = 1;
  for (int boundaryIndex = 0; boundaryIndex < static_cast<int>((options.STRegions.STBoundaries.size() - 1)); ++boundaryIndex) {
    double thisBoundary = options.STRegions.STBoundaries.at(boundaryIndex);
    double nextBoundary = options.STRegions.STBoundaries.at(1+boundaryIndex);
    if (nextBoundary >= options.STNormTarget) {
      st_boundaries_for_unbinned_chi2.push_back(thisBoundary);
      rooVar_ST.setRange(("range_STBinIndex_" + std::to_string(st_boundaries_for_unbinned_chi2_runningIndex)).c_str(), thisBoundary, nextBoundary);
      ++st_boundaries_for_unbinned_chi2_runningIndex;
    }
  }
  st_boundaries_for_unbinned_chi2.push_back(options.STRegions.STBoundaries.at(options.STRegions.STBoundaries.size() - 1));
  rooVar_ST.setRange("normRange", options.STRegions.STAxis.GetBinLowEdge(options.STRegions.STAxis.FindFixBin(options.STNormTarget)), options.STRegions.STAxis.GetBinUpEdge(options.STRegions.STAxis.FindFixBin(options.STNormTarget)));
  rooVar_ST.setRange("fitRange", options.STRegions.STAxis.GetBinLowEdge(options.STRegions.STAxis.FindFixBin(options.STNormTarget)), options.PDF_STMax);
  rooVar_ST.setRange("plotRange", options.PDF_STMin, options.PDF_STMax);

  // If fits are to be read in from another file, do it now
  if (options.readParametersFromFile) {
    std::string lineFromFile;
    std::ifstream inputFileObject(options.inputParametersFileName.c_str());
    assert(inputFileObject.is_open());
    for (int nJetsBin = 3; nJetsBin <= 6; ++nJetsBin) {
      assert(getline(inputFileObject, lineFromFile)); fitParametersUnbinned["slope_slopeOnlyFit"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_slopeOnlyFit_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject, lineFromFile)); fitParametersUnbinned["slopeError_slopeOnlyFit"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slopeError_slopeOnlyFit_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject, lineFromFile)); fitParametersUnbinned["sqrt_sqrtOnlyFit"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "sqrt_sqrtOnlyFit_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject, lineFromFile)); fitParametersUnbinned["sqrtError_sqrtOnlyFit"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "sqrtError_sqrtOnlyFit_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject, lineFromFile)); fitParametersUnbinned["slope_combinedFit"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slope_combinedFit_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject, lineFromFile)); fitParametersUnbinned["slopeError_combinedFit"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "slopeError_combinedFit_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject, lineFromFile)); fitParametersUnbinned["sqrt_combinedFit"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "sqrt_combinedFit_" + std::to_string(nJetsBin) + "Jets");
      assert(getline(inputFileObject, lineFromFile)); fitParametersUnbinned["sqrtError_combinedFit"][nJetsBin] = parseLineForFloatWithCheck(lineFromFile, "sqrtError_combinedFit_" + std::to_string(nJetsBin) + "Jets");
    }
    assert(!(getline(inputFileObject, lineFromFile))); // makes sure there's nothing else in the input file
    inputFileObject.close();
  }

  // First get the 2-jets RooKeysPdf
  RooKeysPdf pdf_2Jets = RooKeysPdf("pdf_2Jets", "pdf_2Jets", rooVar_ST, *(STDataSets.at(2)), RooKeysPdf::MirrorLeft, 1.5);

  // Plot 2-jets shape and dataset
  TCanvas pdfCanvas = TCanvas("c_dataSetAndPdf_2Jets", "c_dataSetAndPdf_2Jets", 2560, 1440);
  RooPlot* rooFrame = rooVar_ST.frame();
  (STDataSets.at(2))->plotOn(rooFrame, Binning(options.PDF_nSTBins, options.PDF_STMin, options.PDF_STMax));
  pdf_2Jets.plotOn(rooFrame, NormRange("normRange"), LineColor(static_cast<EColor>(kBlue)), LineWidth(1));
  rooFrame->Draw();
  pdfCanvas.Update();
  rooFrame->SetMinimum((rooFrame->GetMaximum())/10000.);
  gPad->SetLogy();
  pdfCanvas.Update();
  pdfCanvas.SaveAs((options.outputFolder + "/pdfAndData_2JetsBin_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".pdf").c_str());

  // Next plot the 2-jets pdf on top of other nJets bins
  for (int nJetsBin = 3; nJetsBin <= 6; ++nJetsBin) {
    printSeparator();
    std::cout << "Doing unbinned analysis at nJetsBin = " << nJetsBin << std::endl;

    // initialize now to be able to use binning
    TH1F dataHist(("dataHistCopy_forChi2_" + std::to_string(nJetsBin) + "Jets").c_str(), ("dataHistCopy_forChi2_" + std::to_string(nJetsBin) + "Jets").c_str(), (st_boundaries_for_unbinned_chi2.size()-1), &(st_boundaries_for_unbinned_chi2.at(0)));
    double totalIntegral;

    // unadjusted
    TCanvas pdfCanvas = TCanvas(("c_dataSetAndPdf_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("c_dataSetAndPdf_" + std::to_string(nJetsBin) + "JetsBin").c_str(), 2560, 1440);
    TLegend legend_dataSetsAndPdf = TLegend(0.6, 0.8, 0.9, 0.9);
    RooPlot* rooFrame = rooVar_ST.frame();
    (STDataSets.at(nJetsBin))->plotOn(rooFrame, Binning(options.PDF_nSTBins, options.PDF_STMin, options.PDF_STMax));
    pdf_2Jets.plotOn(rooFrame, NormRange("normRange"), LineColor(static_cast<EColor>(kBlue)), LineWidth(1));
    RooKeysPdf pdf_2Jets_copyForNLL = RooKeysPdf(("pdf_2Jets_copyForNLL_at" + std::to_string(nJetsBin) + "Jets").c_str(), ("pdf_2Jets_copyForNLL_at" + std::to_string(nJetsBin) + "Jets").c_str(), rooVar_ST, *(STDataSets.at(2)), RooKeysPdf::MirrorLeft, 1.5);
    double nll_unadjusted = (pdf_2Jets_copyForNLL.createNLL(*(STDataSets.at(nJetsBin)), Verbose(kFALSE), Range(options.STRegions.STAxis.GetBinLowEdge(options.STRegions.STAxis.FindFixBin(options.STNormTarget)), options.PDF_STMax)))->getVal();
    std::map<int, double> integralsOverBins_unadjusted;
    for (int dataHistBinIndex = 1; dataHistBinIndex <= dataHist.GetXaxis()->GetNbins(); ++dataHistBinIndex) {
      rooVar_ST.setRange(options.PDF_STMin, options.PDF_STMax);
      integralsOverBins_unadjusted[dataHistBinIndex] = (pdf_2Jets_copyForNLL.createIntegral(rooVar_ST, NormSet(rooVar_ST), Range(("range_STBinIndex_" + std::to_string(dataHistBinIndex)).c_str())))->getVal();
    }
    totalIntegral = (pdf_2Jets_copyForNLL.createIntegral(rooVar_ST, NormSet(rooVar_ST), Range("plotRange")))->getVal();
    assert(totalIntegral > 0.);
    TLegendEntry *legendEntry_unadjusted_slope = legend_dataSetsAndPdf.AddEntry(&pdf_2Jets, "2 jets kernel, unadjusted");
    legendEntry_unadjusted_slope->SetMarkerColor(static_cast<EColor>(kBlue));
    legendEntry_unadjusted_slope->SetLineColor(static_cast<EColor>(kBlue));
    legendEntry_unadjusted_slope->SetTextColor(static_cast<EColor>(kBlue));

    // slope correction only
    std::string slopeVar_name = "rooVar_slope_" + std::to_string(nJetsBin) + "JetsBin";
    RooRealVar rooVar_slope(slopeVar_name.c_str(), slopeVar_name.c_str(), 0.0, rooVar_slope_minVal, 5.0);
    std::string function_slopeAdjustment = "1.0 + (" + slopeVar_name + "*((roo_ST/" + std::to_string(options.STNormTarget) + ") - 1.0))";
    std::cout << "function_slopeAdjustment: " << function_slopeAdjustment << std::endl;
    RooGenericPdf pdf_slopeAdjustment("slopeAdjustment", "slopeAdjustment", function_slopeAdjustment.c_str(), RooArgSet(rooVar_ST, rooVar_slope));
    RooProdPdf pdf_nJets_adjusted_slopeOnly(("pdf_2Jets_slopeOnlyAdjustment_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("pdf_2Jets_slopeOnlyAdjustment_" + std::to_string(nJetsBin) + "JetsBin").c_str(), RooArgSet(pdf_2Jets, pdf_slopeAdjustment));
    double nll_slopeOnly;
    if (options.readParametersFromFile) {
      rooVar_slope.setVal((fitParametersUnbinned.at("slope_slopeOnlyFit")).at(nJetsBin));
      rooVar_slope.setError((fitParametersUnbinned.at("slopeError_slopeOnlyFit")).at(nJetsBin));
      nll_slopeOnly = (pdf_nJets_adjusted_slopeOnly.createNLL(*(STDataSets.at(nJetsBin)), Verbose(kFALSE), Range(options.STRegions.STAxis.GetBinLowEdge(options.STRegions.STAxis.FindFixBin(options.STNormTarget)), options.PDF_STMax)))->getVal();
    }
    else {
      RooFitResult *fitResult = pdf_nJets_adjusted_slopeOnly.fitTo(*(STDataSets.at(nJetsBin)), Range(options.STRegions.STAxis.GetBinLowEdge(options.STRegions.STAxis.FindFixBin(options.STNormTarget)), options.PDF_STMax), Optimize(kFALSE), PrintLevel(-1), Verbose(kFALSE), PrintEvalErrors(2), SumW2Error(kTRUE), Save(kTRUE));
      nll_slopeOnly = fitResult->minNll();
      fitParametersUnbinned["slope_slopeOnlyFit"][nJetsBin] = rooVar_slope.getValV();
      fitParametersUnbinnedList.push_back(std::string("float slope_slopeOnlyFit_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("slope_slopeOnlyFit")).at(nJetsBin))));
      fitParametersUnbinned["slopeError_slopeOnlyFit"][nJetsBin] = rooVar_slope.getError();
      fitParametersUnbinnedList.push_back(std::string("float slopeError_slopeOnlyFit_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("slopeError_slopeOnlyFit")).at(nJetsBin))));
    }
    rooVar_slope.Print();
    std::map<int, double> integralsOverBins_slopeOnly;
    for (int dataHistBinIndex = 1; dataHistBinIndex <= dataHist.GetXaxis()->GetNbins(); ++dataHistBinIndex) {
      rooVar_ST.setRange(options.PDF_STMin, options.PDF_STMax);
      integralsOverBins_slopeOnly[dataHistBinIndex] = (pdf_nJets_adjusted_slopeOnly.createIntegral(rooVar_ST, NormSet(rooVar_ST), Range(("range_STBinIndex_" + std::to_string(dataHistBinIndex)).c_str())))->getVal();
    }
    totalIntegral = (pdf_nJets_adjusted_slopeOnly.createIntegral(rooVar_ST, NormSet(rooVar_ST), Range("plotRange")))->getVal();
    assert(totalIntegral > 0.);
    pdf_nJets_adjusted_slopeOnly.plotOn(rooFrame, NormRange("normRange"), LineColor(static_cast<EColor>(kRed+1)), LineWidth(1));
    // plus and minus one-sigma plotted with dashed linestyle
    rooVar_slope.setVal((fitParametersUnbinned.at("slope_slopeOnlyFit")).at(nJetsBin) + (fitParametersUnbinned.at("slopeError_slopeOnlyFit")).at(nJetsBin));
    pdf_nJets_adjusted_slopeOnly.plotOn(rooFrame, NormRange("normRange"), LineColor(static_cast<EColor>(kRed+1)), LineStyle(kDashed), LineWidth(1));
    rooVar_slope.setVal((fitParametersUnbinned.at("slope_slopeOnlyFit")).at(nJetsBin) - (fitParametersUnbinned.at("slopeError_slopeOnlyFit")).at(nJetsBin));
    pdf_nJets_adjusted_slopeOnly.plotOn(rooFrame, NormRange("normRange"), LineColor(static_cast<EColor>(kRed+1)), LineStyle(kDashed), LineWidth(1));
    rooVar_slope.setVal((fitParametersUnbinned.at("slope_slopeOnlyFit")).at(nJetsBin));
    TLegendEntry *legendEntry_adjusted_slopeOnly = legend_dataSetsAndPdf.AddEntry(&pdf_nJets_adjusted_slopeOnly, "2 jets kernel + slope adjustment");
    legendEntry_adjusted_slopeOnly->SetMarkerColor(static_cast<EColor>(kRed+1));
    legendEntry_adjusted_slopeOnly->SetLineColor(static_cast<EColor>(kRed+1));
    legendEntry_adjusted_slopeOnly->SetTextColor(static_cast<EColor>(kRed+1));

    // sqrt correction only
    std::string sqrtVar_name = "rooVar_sqrt_" + std::to_string(nJetsBin) + "JetsBin";
    RooRealVar rooVar_sqrt(sqrtVar_name.c_str(), sqrtVar_name.c_str(), 0.0, rooVar_sqrt_minVal, 25.0);
    std::string function_sqrtAdjustment = "1.0 + (" + sqrtVar_name + "*(sqrt(roo_ST/" + std::to_string(options.STNormTarget) + ") - 1.0))";
    std::cout << "function_sqrtAdjustment: " << function_sqrtAdjustment << std::endl;
    RooGenericPdf pdf_sqrtAdjustment("sqrtAdjustment", "sqrtAdjustment", function_sqrtAdjustment.c_str(), RooArgSet(rooVar_ST, rooVar_sqrt));
    RooProdPdf pdf_nJets_adjusted_sqrtOnly(("pdf_2Jets_sqrtOnlyAdjustment_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("pdf_2Jets_sqrtOnlyAdjustment" + std::to_string(nJetsBin) + "JetsBin").c_str(), RooArgSet(pdf_2Jets, pdf_sqrtAdjustment));
    double nll_sqrtOnly;
    if (options.readParametersFromFile) {
      rooVar_sqrt.setVal((fitParametersUnbinned.at("sqrt_sqrtOnlyFit")).at(nJetsBin));
      rooVar_sqrt.setError((fitParametersUnbinned.at("sqrtError_sqrtOnlyFit")).at(nJetsBin));
      nll_sqrtOnly = (pdf_nJets_adjusted_sqrtOnly.createNLL(*(STDataSets.at(nJetsBin)), Verbose(kFALSE), Range(options.STRegions.STAxis.GetBinLowEdge(options.STRegions.STAxis.FindFixBin(options.STNormTarget)), options.PDF_STMax)))->getVal();
    }
    else {
      RooFitResult *fitResult = pdf_nJets_adjusted_sqrtOnly.fitTo(*(STDataSets.at(nJetsBin)), Range(options.STRegions.STAxis.GetBinLowEdge(options.STRegions.STAxis.FindFixBin(options.STNormTarget)), options.PDF_STMax), Optimize(kFALSE), PrintLevel(-1), Verbose(kFALSE), PrintEvalErrors(2), SumW2Error(kTRUE), Save(kTRUE));
      nll_sqrtOnly = fitResult->minNll();
      fitParametersUnbinned["sqrt_sqrtOnlyFit"][nJetsBin] = rooVar_sqrt.getValV();
      fitParametersUnbinnedList.push_back(std::string("float sqrt_sqrtOnlyFit_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("sqrt_sqrtOnlyFit")).at(nJetsBin))));
      fitParametersUnbinned["sqrtError_sqrtOnlyFit"][nJetsBin] = rooVar_sqrt.getError();
      fitParametersUnbinnedList.push_back(std::string("float sqrtError_sqrtOnlyFit_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("sqrtError_sqrtOnlyFit")).at(nJetsBin))));
    }
    rooVar_sqrt.Print();
    std::map<int, double> integralsOverBins_sqrtOnly;
    for (int dataHistBinIndex = 1; dataHistBinIndex <= dataHist.GetXaxis()->GetNbins(); ++dataHistBinIndex) {
      rooVar_ST.setRange(options.PDF_STMin, options.PDF_STMax);
      integralsOverBins_sqrtOnly[dataHistBinIndex] = (pdf_nJets_adjusted_sqrtOnly.createIntegral(rooVar_ST, NormSet(rooVar_ST), Range(("range_STBinIndex_" + std::to_string(dataHistBinIndex)).c_str())))->getVal();
    }
    totalIntegral = (pdf_nJets_adjusted_sqrtOnly.createIntegral(rooVar_ST, NormSet(rooVar_ST), Range("plotRange")))->getVal();
    assert(totalIntegral > 0.);
    pdf_nJets_adjusted_sqrtOnly.plotOn(rooFrame, NormRange("normRange"), LineColor(static_cast<EColor>(kGreen+3)), LineWidth(1));
    // plus and minus one-sigma plotted with dashed linestyle
    rooVar_sqrt.setVal((fitParametersUnbinned.at("sqrt_sqrtOnlyFit")).at(nJetsBin) + (fitParametersUnbinned.at("sqrtError_sqrtOnlyFit")).at(nJetsBin));
    pdf_nJets_adjusted_sqrtOnly.plotOn(rooFrame, NormRange("normRange"), LineColor(static_cast<EColor>(kGreen+3)), LineStyle(kDashed), LineWidth(1));
    rooVar_sqrt.setVal((fitParametersUnbinned.at("sqrt_sqrtOnlyFit")).at(nJetsBin) - (fitParametersUnbinned.at("sqrtError_sqrtOnlyFit")).at(nJetsBin));
    pdf_nJets_adjusted_sqrtOnly.plotOn(rooFrame, NormRange("normRange"), LineColor(static_cast<EColor>(kGreen+3)), LineStyle(kDashed), LineWidth(1));
    rooVar_sqrt.setVal((fitParametersUnbinned.at("sqrt_sqrtOnlyFit")).at(nJetsBin));
    TLegendEntry *legendEntry_adjusted_sqrtOnly = legend_dataSetsAndPdf.AddEntry(&pdf_nJets_adjusted_sqrtOnly, "2 jets kernel + sqrt adjustment");
    legendEntry_adjusted_sqrtOnly->SetMarkerColor(static_cast<EColor>(kGreen+3));
    legendEntry_adjusted_sqrtOnly->SetLineColor(static_cast<EColor>(kGreen+3));
    legendEntry_adjusted_sqrtOnly->SetTextColor(static_cast<EColor>(kGreen+3));

    // slope correction + sqrt correction
    std::string slopeVar_combinedFit_name = "rooVar_slope_combinedFit_" + std::to_string(nJetsBin) + "JetsBin";
    RooRealVar rooVar_slope_combinedFit(slopeVar_combinedFit_name.c_str(), slopeVar_combinedFit_name.c_str(), 0.5*(fitParametersUnbinned.at("slope_slopeOnlyFit")).at(nJetsBin), rooVar_slope_minVal, 5.0);
    function_slopeAdjustment = "1.0 + (" + slopeVar_combinedFit_name + "*((roo_ST/" + std::to_string(options.STNormTarget) + ") - 1.0))";
    std::cout << "For combined fit, function_slopeAdjustment: " << function_slopeAdjustment << std::endl;
    RooGenericPdf pdf_slopeAdjustment_combinedFit("slopeAdjustment_combinedFit", "slopeAdjustment_combinedFit", function_slopeAdjustment.c_str(), RooArgSet(rooVar_ST, rooVar_slope_combinedFit));
    std::string sqrtVar_combinedFit_name = "rooVar_sqrt_combinedFit_" + std::to_string(nJetsBin) + "JetsBin";
    RooRealVar rooVar_sqrt_combinedFit(sqrtVar_combinedFit_name.c_str(), sqrtVar_combinedFit_name.c_str(), 0.5*(fitParametersUnbinned.at("sqrt_sqrtOnlyFit")).at(nJetsBin), rooVar_sqrt_minVal, 25.0);
    function_sqrtAdjustment = "1.0 + (" + sqrtVar_combinedFit_name + "*(sqrt(roo_ST/" + std::to_string(options.STNormTarget) + ") - 1.0))";
    std::cout << "For combined fit, function_sqrtAdjustment: " << function_sqrtAdjustment << std::endl;
    RooGenericPdf pdf_sqrtAdjustment_combinedFit("sqrtAdjustment_combinedFit", "sqrtAdjustment_combinedFit", function_sqrtAdjustment.c_str(), RooArgSet(rooVar_ST, rooVar_sqrt_combinedFit));
    RooProdPdf pdf_nJets_adjusted_combined(("pdf_2Jets_slopeAndSqrtAdjustment_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("pdf_2Jets_slopeAndSqrtAdjustment_" + std::to_string(nJetsBin) + "JetsBin").c_str(), RooArgSet(pdf_2Jets, pdf_slopeAdjustment_combinedFit, pdf_sqrtAdjustment_combinedFit));
    double nll_combined;
    if (options.readParametersFromFile) {
      rooVar_slope_combinedFit.setVal((fitParametersUnbinned.at("slope_combinedFit")).at(nJetsBin));
      rooVar_slope_combinedFit.setError((fitParametersUnbinned.at("slopeError_combinedFit")).at(nJetsBin));
      rooVar_sqrt_combinedFit.setVal((fitParametersUnbinned.at("sqrt_combinedFit")).at(nJetsBin));
      rooVar_sqrt_combinedFit.setError((fitParametersUnbinned.at("sqrtError_combinedFit")).at(nJetsBin));
      nll_combined = (pdf_nJets_adjusted_combined.createNLL(*(STDataSets.at(nJetsBin)), Verbose(kFALSE), Range(options.STRegions.STAxis.GetBinLowEdge(options.STRegions.STAxis.FindFixBin(options.STNormTarget)), options.PDF_STMax)))->getVal();
    }
    else {
      RooFitResult *fitResult = pdf_nJets_adjusted_combined.fitTo(*(STDataSets.at(nJetsBin)), Range(options.STRegions.STAxis.GetBinLowEdge(options.STRegions.STAxis.FindFixBin(options.STNormTarget)), options.PDF_STMax), Optimize(kFALSE), PrintLevel(-1), Verbose(kFALSE), PrintEvalErrors(2), SumW2Error(kTRUE), Save(kTRUE));
      nll_combined = fitResult->minNll();
      fitParametersUnbinned["slope_combinedFit"][nJetsBin] = rooVar_slope_combinedFit.getValV();
      fitParametersUnbinnedList.push_back(std::string("float slope_combinedFit_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("slope_combinedFit")).at(nJetsBin))));
      fitParametersUnbinned["slopeError_combinedFit"][nJetsBin] = rooVar_slope_combinedFit.getError();
      fitParametersUnbinnedList.push_back(std::string("float slopeError_combinedFit_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("slopeError_combinedFit")).at(nJetsBin))));
      fitParametersUnbinned["sqrt_combinedFit"][nJetsBin] = rooVar_sqrt_combinedFit.getValV();
      fitParametersUnbinnedList.push_back(std::string("float sqrt_combinedFit_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("sqrt_combinedFit")).at(nJetsBin))));
      fitParametersUnbinned["sqrtError_combinedFit"][nJetsBin] = rooVar_sqrt_combinedFit.getError();
      fitParametersUnbinnedList.push_back(std::string("float sqrtError_combinedFit_" + std::to_string(nJetsBin) + "Jets=" + std::to_string((fitParametersUnbinned.at("sqrtError_combinedFit")).at(nJetsBin))));
    }
    rooVar_slope_combinedFit.Print();
    rooVar_sqrt_combinedFit.Print();
    std::map<int, double> integralsOverBins_combined;
    for (int dataHistBinIndex = 1; dataHistBinIndex <= dataHist.GetXaxis()->GetNbins(); ++dataHistBinIndex) {
      rooVar_ST.setRange(options.PDF_STMin, options.PDF_STMax);
      integralsOverBins_combined[dataHistBinIndex] = (pdf_nJets_adjusted_combined.createIntegral(rooVar_ST, NormSet(rooVar_ST), Range(("range_STBinIndex_" + std::to_string(dataHistBinIndex)).c_str())))->getVal();
    }
    totalIntegral = (pdf_nJets_adjusted_combined.createIntegral(rooVar_ST, NormSet(rooVar_ST), Range("plotRange")))->getVal();
    assert(totalIntegral > 0.);
    pdf_nJets_adjusted_combined.plotOn(rooFrame, NormRange("normRange"), LineColor(static_cast<EColor>(kViolet)), LineWidth(1));
    // plus and minus one-sigma plotted with dashed linestyle
    rooVar_slope_combinedFit.setVal((fitParametersUnbinned.at("slope_combinedFit")).at(nJetsBin) + std::sqrt(0.5)*(fitParametersUnbinned.at("slopeError_combinedFit")).at(nJetsBin));
    rooVar_sqrt_combinedFit.setVal((fitParametersUnbinned.at("sqrt_combinedFit")).at(nJetsBin) + std::sqrt(0.5)*(fitParametersUnbinned.at("sqrtError_combinedFit")).at(nJetsBin));
    pdf_nJets_adjusted_combined.plotOn(rooFrame, NormRange("normRange"), LineColor(static_cast<EColor>(kViolet)), LineStyle(kDashed), LineWidth(1));
    rooVar_slope_combinedFit.setVal((fitParametersUnbinned.at("slope_combinedFit")).at(nJetsBin) - std::sqrt(0.5)*(fitParametersUnbinned.at("slopeError_combinedFit")).at(nJetsBin));
    rooVar_sqrt_combinedFit.setVal((fitParametersUnbinned.at("sqrt_combinedFit")).at(nJetsBin) - std::sqrt(0.5)*(fitParametersUnbinned.at("sqrtError_combinedFit")).at(nJetsBin));
    pdf_nJets_adjusted_combined.plotOn(rooFrame, NormRange("normRange"), LineColor(static_cast<EColor>(kViolet)), LineStyle(kDashed), LineWidth(1));
    rooVar_slope_combinedFit.setVal((fitParametersUnbinned.at("slope_combinedFit")).at(nJetsBin));
    rooVar_sqrt_combinedFit.setVal((fitParametersUnbinned.at("sqrt_combinedFit")).at(nJetsBin));
    TLegendEntry *legendEntry_adjusted_combined = legend_dataSetsAndPdf.AddEntry(&pdf_nJets_adjusted_combined, "2 jets kernel + slope adjustment + sqrt adjustment");
    legendEntry_adjusted_combined->SetMarkerColor(static_cast<EColor>(kViolet));
    legendEntry_adjusted_combined->SetLineColor(static_cast<EColor>(kViolet));
    legendEntry_adjusted_combined->SetTextColor(static_cast<EColor>(kViolet));

    rooFrame->Draw();
    pdfCanvas.Update();
    rooFrame->SetMinimum((rooFrame->GetMaximum())/10000.);
    gPad->SetLogy();
    pdfCanvas.Update();
    legend_dataSetsAndPdf.SetFillStyle(0);
    legend_dataSetsAndPdf.Draw();
    pdfCanvas.Update();
    pdfCanvas.SaveAs((options.outputFolder + "/pdfAndData_" + std::to_string(nJetsBin) + "JetsBin_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".pdf").c_str());

    std::cout << "Getting p-values from computed nlls..." << std::endl;
    nll_pValues["unadjusted_vs_slopeOnly"][nJetsBin] = get_nll_prob(nll_unadjusted, nll_slopeOnly, 1);
    nll_pValues["slopeOnly_vs_combined"][nJetsBin] = get_nll_prob(nll_slopeOnly, nll_combined, 1);
    nll_pValues["unadjusted_vs_sqrtOnly"][nJetsBin] = get_nll_prob(nll_unadjusted, nll_sqrtOnly, 1);
    nll_pValues["sqrtOnly_vs_combined"][nJetsBin] = get_nll_prob(nll_sqrtOnly, nll_combined, 1);

    // Copy STDistributions into dataHist for the chi2 analysis
    dataHist.Sumw2();
    int ndf_unadjusted_forChi2 = -1; // start from (-1) because one degree of freedom is eaten up by the normalization
    for (int dataHistBinIndex = 1; dataHistBinIndex <= dataHist.GetXaxis()->GetNbins(); ++dataHistBinIndex) {
      double binCenter = dataHist.GetXaxis()->GetBinCenter(dataHistBinIndex);
      double binContent = (STDistributions.at(nJetsBin))->GetBinContent((STDistributions.at(nJetsBin))->FindFixBin(binCenter));
      double binError = (STDistributions.at(nJetsBin))->GetBinError((STDistributions.at(nJetsBin))->FindFixBin(binCenter));
      dataHist.SetBinContent(dataHistBinIndex, binContent);
      dataHist.SetBinError(dataHistBinIndex, binError);
      if (binContent > 0.) {
        assert(binError > 0.);
        ++ndf_unadjusted_forChi2;
      }
    }
    // getDensityHistogramChiSquareWRTFunction(TH1F& inputDensityHistogram, std::map<int, double> functionIntegralsOverBins)
    double chi2_unadjusted = getDensityHistogramChiSquareWRTFunction(dataHist, integralsOverBins_unadjusted);
    double chi2_slopeOnly = getDensityHistogramChiSquareWRTFunction(dataHist, integralsOverBins_slopeOnly);
    double chi2_sqrtOnly = getDensityHistogramChiSquareWRTFunction(dataHist, integralsOverBins_sqrtOnly);
    double chi2_combined = getDensityHistogramChiSquareWRTFunction(dataHist, integralsOverBins_combined);

    std::cout << "Getting p-values from computed chi2 values..." << std::endl;
    ftest_pValues_unbinnedFit["unadjusted_vs_slopeOnly"][nJetsBin] = get_fTest_prob(chi2_unadjusted, chi2_slopeOnly, ndf_unadjusted_forChi2, ndf_unadjusted_forChi2-1);
    ftest_pValues_unbinnedFit["slopeOnly_vs_combined"][nJetsBin] = get_fTest_prob(chi2_slopeOnly, chi2_combined, ndf_unadjusted_forChi2-1, ndf_unadjusted_forChi2-2);
    ftest_pValues_unbinnedFit["unadjusted_vs_sqrtOnly"][nJetsBin] = get_fTest_prob(chi2_unadjusted, chi2_sqrtOnly, ndf_unadjusted_forChi2, ndf_unadjusted_forChi2-1);
    ftest_pValues_unbinnedFit["sqrtOnly_vs_combined"][nJetsBin] = get_fTest_prob(chi2_sqrtOnly, chi2_combined, ndf_unadjusted_forChi2-1, ndf_unadjusted_forChi2-2);
    printSeparator();
  }

  // Print f-test prob values in a LaTeX-formatted table
  std::cout << "f-test prob values (formatted):" << std::endl;
  std::cout << "\\begin{tabular}{|l|c|c|c|}" << std::endl;
  std::cout << "  \\hline" << std::endl;
  std::cout << "  f-prob & const\\_vs\\_lin & lin\\_vs\\_quad & fixed\\_lin\\_vs\\_lin \\\\ \\hline" << std::endl;
  std::cout << std::fixed << std::setprecision(3) << "  nJets = 3 & " << (fTestProbValues.at("const_vs_lin")).at(3) << " & " << (fTestProbValues.at("lin_vs_quad")).at(3) << " & " << (fTestProbValues.at("constrained_lin_vs_lin")).at(3) << " \\\\ \\hline" << std::endl;
  std::cout << std::fixed << std::setprecision(3) << "  nJets = 4 & " << (fTestProbValues.at("const_vs_lin")).at(4) << " & " << (fTestProbValues.at("lin_vs_quad")).at(4) << " & " << (fTestProbValues.at("constrained_lin_vs_lin")).at(4) << " \\\\ \\hline" << std::endl;
  std::cout << std::fixed << std::setprecision(3) << "  nJets = 5 & " << (fTestProbValues.at("const_vs_lin")).at(5) << " & " << (fTestProbValues.at("lin_vs_quad")).at(5) << " & " << (fTestProbValues.at("constrained_lin_vs_lin")).at(5) << " \\\\ \\hline" << std::endl;
  std::cout << std::fixed << std::setprecision(3) << "  nJets $\\geq$ 6 & " << (fTestProbValues.at("const_vs_lin")).at(6) << " & " << (fTestProbValues.at("lin_vs_quad")).at(6) << " & " << (fTestProbValues.at("constrained_lin_vs_lin")).at(6) << " \\\\ \\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;

  // Print best fit slopes in LaTeX-formatted table
  std::cout << "slope values (formatted):" << std::endl;
  std::cout << "\\begin{tabular}{|l|c|}" << std::endl;
  std::cout << "  \\hline" << std::endl;
  std::cout << "  nJets Bin & Best-fit slope \\\\ \\hline" << std::endl;
  std::cout << std::fixed << std::setprecision(3) << "  nJets = 3 & " << (fitParametersBinned.at("slope")).at(3) << " $\\pm$ " << (fitParametersBinned.at("slopeError")).at(3) << " \\\\ \\hline" << std::endl;
  std::cout << std::fixed << std::setprecision(3) << "  nJets = 4 & " << (fitParametersBinned.at("slope")).at(4) << " $\\pm$ " << (fitParametersBinned.at("slopeError")).at(4) << " \\\\ \\hline" << std::endl;
  std::cout << std::fixed << std::setprecision(3) << "  nJets = 5 & " << (fitParametersBinned.at("slope")).at(5) << " $\\pm$ " << (fitParametersBinned.at("slopeError")).at(5) << " \\\\ \\hline" << std::endl;
  std::cout << std::fixed << std::setprecision(3) << "  nJets $\\geq$ 6 & " << (fitParametersBinned.at("slope")).at(6) << " $\\pm$ " << (fitParametersBinned.at("slopeError")).at(6) << " \\\\ \\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;

  // Print p-values from nlls in a LaTeX-formatted table
  std::cout << "p-values from nlls (formatted):" << std::endl;
  std::cout << "\\begin{tabular}{|l|p{0.18\\textwidth}|p{0.18\\textwidth}|p{0.18\\textwidth}|p{0.18\\textwidth}|}" << std::endl;
  std::cout << "  \\hline" << std::endl;
  std::cout << "  p-values & unadjusted vs slope-only & slope-only vs combined & unadjusted vs sqrt-only & sqrt-only vs combined \\\\ \\hline" << std::endl;
  std::cout << std::fixed << std::setprecision(3) << "  nJets = 3 & " << (nll_pValues.at("unadjusted_vs_slopeOnly")).at(3) << " & " << (nll_pValues.at("slopeOnly_vs_combined")).at(3) << " & " << (nll_pValues.at("unadjusted_vs_sqrtOnly")).at(3) << " & " << (nll_pValues.at("sqrtOnly_vs_combined")).at(3) << " \\\\ \\hline" << std::endl;
  std::cout << std::fixed << std::setprecision(3) << "  nJets = 4 & " << (nll_pValues.at("unadjusted_vs_slopeOnly")).at(4) << " & " << (nll_pValues.at("slopeOnly_vs_combined")).at(4) << " & " << (nll_pValues.at("unadjusted_vs_sqrtOnly")).at(4) << " & " << (nll_pValues.at("sqrtOnly_vs_combined")).at(4) << " \\\\ \\hline" << std::endl;
  std::cout << std::fixed << std::setprecision(3) << "  nJets = 5 & " << (nll_pValues.at("unadjusted_vs_slopeOnly")).at(5) << " & " << (nll_pValues.at("slopeOnly_vs_combined")).at(5) << " & " << (nll_pValues.at("unadjusted_vs_sqrtOnly")).at(5) << " & " << (nll_pValues.at("sqrtOnly_vs_combined")).at(5) << " \\\\ \\hline" << std::endl;
  std::cout << std::fixed << std::setprecision(3) << "  nJets $\\geq$ 6 & " << (nll_pValues.at("unadjusted_vs_slopeOnly")).at(6) << " & " << (nll_pValues.at("slopeOnly_vs_combined")).at(6) << " & " << (nll_pValues.at("unadjusted_vs_sqrtOnly")).at(6) << " & " << (nll_pValues.at("sqrtOnly_vs_combined")).at(6) << " \\\\ \\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;

  // Print ftest pvalues from estimated chi2 in a LaTeX-formatted table
  std::cout << "p-values from estimated chi2 (formatted):" << std::endl;
  std::cout << "\\begin{tabular}{|l|p{0.18\\textwidth}|p{0.18\\textwidth}|p{0.18\\textwidth}|p{0.18\\textwidth}|}" << std::endl;
  std::cout << "  \\hline" << std::endl;
  std::cout << "  p-values & unadjusted vs slope-only & slope-only vs combined & unadjusted vs sqrt-only & sqrt-only vs combined \\\\ \\hline" << std::endl;
  std::cout << std::fixed << std::setprecision(3) << "  nJets = 3 & " << (ftest_pValues_unbinnedFit.at("unadjusted_vs_slopeOnly")).at(3) << " & " << (ftest_pValues_unbinnedFit.at("slopeOnly_vs_combined")).at(3) << " & " << (ftest_pValues_unbinnedFit.at("unadjusted_vs_sqrtOnly")).at(3) << " & " << (ftest_pValues_unbinnedFit.at("sqrtOnly_vs_combined")).at(3) << " \\\\ \\hline" << std::endl;
  std::cout << std::fixed << std::setprecision(3) << "  nJets = 4 & " << (ftest_pValues_unbinnedFit.at("unadjusted_vs_slopeOnly")).at(4) << " & " << (ftest_pValues_unbinnedFit.at("slopeOnly_vs_combined")).at(4) << " & " << (ftest_pValues_unbinnedFit.at("unadjusted_vs_sqrtOnly")).at(4) << " & " << (ftest_pValues_unbinnedFit.at("sqrtOnly_vs_combined")).at(4) << " \\\\ \\hline" << std::endl;
  std::cout << std::fixed << std::setprecision(3) << "  nJets = 5 & " << (ftest_pValues_unbinnedFit.at("unadjusted_vs_slopeOnly")).at(5) << " & " << (ftest_pValues_unbinnedFit.at("slopeOnly_vs_combined")).at(5) << " & " << (ftest_pValues_unbinnedFit.at("unadjusted_vs_sqrtOnly")).at(5) << " & " << (ftest_pValues_unbinnedFit.at("sqrtOnly_vs_combined")).at(5) << " \\\\ \\hline" << std::endl;
  std::cout << std::fixed << std::setprecision(3) << "  nJets $\\geq$ 6 & " << (ftest_pValues_unbinnedFit.at("unadjusted_vs_slopeOnly")).at(6) << " & " << (ftest_pValues_unbinnedFit.at("slopeOnly_vs_combined")).at(6) << " & " << (ftest_pValues_unbinnedFit.at("unadjusted_vs_sqrtOnly")).at(6) << " & " << (ftest_pValues_unbinnedFit.at("sqrtOnly_vs_combined")).at(6) << " \\\\ \\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;

  // Print best fit slopes from unbinned fit in LaTeX-formatted table
  std::cout << "slope values from unbinned fit (formatted):" << std::endl;
  std::cout << "\\begin{tabular}{|l|c|}" << std::endl;
  std::cout << "  \\hline" << std::endl;
  std::cout << "  nJets Bin & Best-fit slope \\\\ \\hline" << std::endl;
  std::cout << std::fixed << std::setprecision(3) << "  nJets = 3 & " << (fitParametersUnbinned.at("slope_slopeOnlyFit")).at(3) << " $\\pm$ " << (fitParametersUnbinned.at("slopeError_slopeOnlyFit")).at(3) << " \\\\ \\hline" << std::endl;
  std::cout << std::fixed << std::setprecision(3) << "  nJets = 4 & " << (fitParametersUnbinned.at("slope_slopeOnlyFit")).at(4) << " $\\pm$ " << (fitParametersUnbinned.at("slopeError_slopeOnlyFit")).at(4) << " \\\\ \\hline" << std::endl;
  std::cout << std::fixed << std::setprecision(3) << "  nJets = 5 & " << (fitParametersUnbinned.at("slope_slopeOnlyFit")).at(5) << " $\\pm$ " << (fitParametersUnbinned.at("slopeError_slopeOnlyFit")).at(5) << " \\\\ \\hline" << std::endl;
  std::cout << std::fixed << std::setprecision(3) << "  nJets $\\geq$ 6 & " << (fitParametersUnbinned.at("slope_slopeOnlyFit")).at(6) << " $\\pm$ " << (fitParametersUnbinned.at("slopeError_slopeOnlyFit")).at(6) << " \\\\ \\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;

  // write parameters for binned fit
  std::ofstream fitParametersBinnedFile((options.outputFolder + "/fitParameters_binnedFit_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".dat").c_str());
  assert(fitParametersBinnedFile.is_open());
  for (int fitParametersBinnedListIndex = 0; fitParametersBinnedListIndex < static_cast<int>(fitParametersBinnedList.size()); ++fitParametersBinnedListIndex) {
    fitParametersBinnedFile << fitParametersBinnedList.at(fitParametersBinnedListIndex) << std::endl;
  }
  fitParametersBinnedFile.close();
  std::cout << "Binned fit parameters written to file: " << (options.outputFolder + "/fitParameters_binnedFit_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".dat") << std::endl;

  // write parameters for unbinned fit
  if (!(options.readParametersFromFile)) {
    std::ofstream fitParametersUnbinnedFile((options.outputFolder + "/fitParameters_unbinnedFit_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".dat").c_str());
    assert(fitParametersUnbinnedFile.is_open());
    for (int fitParametersUnbinnedListIndex = 0; fitParametersUnbinnedListIndex < static_cast<int>(fitParametersUnbinnedList.size()); ++fitParametersUnbinnedListIndex) {
      fitParametersUnbinnedFile << fitParametersUnbinnedList.at(fitParametersUnbinnedListIndex) << std::endl;
    }
    fitParametersUnbinnedFile.close();
    std::cout << "Unbinned fit parameters written to file: " << (options.outputFolder + "/fitParameters_unbinnedFit_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".dat") << std::endl;    
  }

  std::cout << "All done!" << std::endl;
  return EXIT_SUCCESS;
}
