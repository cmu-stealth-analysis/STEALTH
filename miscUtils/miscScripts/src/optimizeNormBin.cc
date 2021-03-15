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

double get_fTest_prob(const double& chi2_1, const double& chi2_2, const int& ndf_1, const int& ndf_2) {
  // Step 1: Get the environment variable EOSTMPAREA, a folder where temp files can be stored
  std::string tmpFolder = std::string(getenv("EOSTMPAREA"));
  if (tmpFolder == "") {
    std::cout << "ERROR: env variable \"EOSTMPAREA\" does not appear to be set." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  
  // Step 2: Generate python command that gets this value
  std::string commandToRun = "python -c \"import tmStatsUtils; print(tmStatsUtils.get_fTest_prob(chi2_1=" + std::to_string(chi2_1) + ", chi2_2=" + std::to_string(chi2_2) + ", ndf_1=" + std::to_string(ndf_1) + ", ndf_2=" + std::to_string(ndf_2) + "))\" > " + tmpFolder + "/optimizeNormBin_tmp.txt 2>&1";
  int ftest_command_return_status = system(commandToRun.c_str());
  if (ftest_command_return_status != 0) {
    std::cout << "ERROR in running command: " << commandToRun << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Step 3: Open output file and read in the value
  std::vector<double> valuesFromFile;
  double valueFromFile;
  std::ifstream inputFileObject((tmpFolder + "/optimizeNormBin_tmp.txt").c_str());
  if (inputFileObject.is_open()) {
    while (inputFileObject >> valueFromFile) {
      valuesFromFile.push_back(valueFromFile);
    }
    inputFileObject.close();
  }
  else {
    std::cout << "ERROR: Unable to open file: " << (tmpFolder + "/optimizeNormBin_tmp.txt") << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Step 4: Check that there is indeed just one value
  if (!(valuesFromFile.size() == 1)) {
    std::cout << "ERROR: this tmp file is in an unexpected format: " << (tmpFolder + "/optimizeNormBin_tmp.txt") << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Step 5: Delete the temp file
  std::string file_remove_command = "rm " + (tmpFolder + "/optimizeNormBin_tmp.txt");
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

  tmArgumentParser argumentParser = tmArgumentParser("Run script that prints useful info about the normalization.");
  argumentParser.addArgument("sourceFilePath", "", true, "Path to file containing list of paths with n-tuplized events.");
  argumentParser.addArgument("outputFolder", "", true, "Output folder.");
  argumentParser.addArgument("selection", "", true, "Name of selection: \"singlemedium\", \"signal_loose\", etc.");
  argumentParser.addArgument("identifier", "", true, "Identifier: \"MC_GJet17\", \"MC_GJet\", etc.");
  argumentParser.addArgument("STBoundariesSourceFile", "STRegionBoundaries_normOptimization.dat", false, "Identifier: \"MC_GJet17\", \"MC_GJet\", etc.");
  argumentParser.addArgument("yearString", "all", false, "String with year: can take values \"2016\", \"2017\", \"2018\", or \"all\".");
  argumentParser.addArgument("STNormTarget", "1250.0", false, "Target value to use for ST norm.");
  argumentParser.addArgument("STNormMax", "1450.0", false, "Bins with a central value larger than this value are not considered.");
  argumentParser.setPassedStringValues(argc, argv);
  optionsStruct options = getOptionsFromParser(argumentParser);
  std::cout << "Options passed:" << std::endl << options << std::endl;

  std::string constrainedLinFunction = "1.0+[0]*(x-" + std::to_string(options.STNormTarget) + ")/1000.0";
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

  TFile *sourceFile = TFile::Open((options.sourceFilePath).c_str(), "READ");
  if ((!(sourceFile->IsOpen())) or (sourceFile->IsZombie())) {
    std::cout << "ERROR: unable to open file with name: " << options.sourceFilePath << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::cout << "Fetching ST datasets and histograms from file: " << options.sourceFilePath << std::endl;
  std::map<int, TH1F*> STDistributions;
  // dataSets = {}
  // STDistributions = {}
  for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
    STDistributions[nJetsBin] = new TH1F();
    std::string inputHistogramName = "h_ST_" + std::to_string(nJetsBin) + "JetsBin";
    std::cout << "Getting histogram with name: " << inputHistogramName << std::endl;
    sourceFile->GetObject(inputHistogramName.c_str(), STDistributions[nJetsBin]);
    STDistributions[nJetsBin]->SetName(inputHistogramName.c_str());
    assert(STDistributions[nJetsBin] != nullptr);
  }
  // dataSets[nJetsBin] = ROOT.RooDataSet()
  //   sourceFile.GetObject("dataSet_{n}Jets".format(n=nJetsBin), dataSets[nJetsBin])
  //   dataSets[nJetsBin].SetName("dataSet_{n}Jets".format(n=nJetsBin))


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

  std::vector<std::string> fitParametersList;
  std::map<std::string, std::map<int, double> > fitParameters;
  std::map<std::string, std::map<int, float> > fTestProbValues;
  // fTestProbValues = {"const_vs_lin": {}, "lin_vs_quad": {}, "constrained_lin_vs_lin": {}};

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
          (fits.at(fit_type))->SetParName(0, "slopeDividedBy1000");
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
          double bestFitSlopeDividedBy1000 = fitResults.at(fit_type)->Parameter(0);
          fitParametersList.push_back(std::string("float bestFitSlopeDividedBy1000_" + std::to_string(nJetsBin) + "Jets=" + std::to_string(bestFitSlopeDividedBy1000)));
          fitParameters["slopeDividedBy1000"][nJetsBin] = bestFitSlopeDividedBy1000;
          double bestFitSlopeErrorDividedBy1000 = fitResults.at(fit_type)->ParError(0);
          fitParametersList.push_back(std::string("float bestFitSlopeErrorDividedBy1000_" + std::to_string(nJetsBin) + "Jets=" + std::to_string(bestFitSlopeErrorDividedBy1000)));
          fitParameters["slopeErrorDividedBy1000"][nJetsBin] = bestFitSlopeErrorDividedBy1000;
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
            // (binCenter, ratioValue, binWidth, ratioError) = ratioGraphPoints[index]
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

        TCanvas outputCanvas = TCanvas(("c_ratios_" + std::to_string(nJetsBin) + "JetsBin").c_str(), ("c_ratios_" + std::to_string(nJetsBin) + "JetsBin").c_str(), 1024, 768);
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

  // write fit parameters
  std::ofstream fitParametersFile((options.outputFolder + "/fitParameters_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".dat").c_str());
  assert(fitParametersFile.is_open());
  for (int fitParametersListIndex = 0; fitParametersListIndex < static_cast<int>(fitParametersList.size()); ++fitParametersListIndex) {
    fitParametersFile << fitParametersList.at(fitParametersListIndex) << std::endl;
  }
  fitParametersFile.close();
  std::cout << "Configuration parameters written to file: " << (options.outputFolder + "/fitParameters_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".dat") << std::endl;

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
  std::cout << std::fixed << std::setprecision(3) << "  nJets = 3 & " << (fitParameters.at("slopeDividedBy1000")).at(3) << " $\\pm$ " << (fitParameters.at("slopeErrorDividedBy1000")).at(3) << " \\\\ \\hline" << std::endl;
  std::cout << std::fixed << std::setprecision(3) << "  nJets = 4 & " << (fitParameters.at("slopeDividedBy1000")).at(4) << " $\\pm$ " << (fitParameters.at("slopeErrorDividedBy1000")).at(4) << " \\\\ \\hline" << std::endl;
  std::cout << std::fixed << std::setprecision(3) << "  nJets = 5 & " << (fitParameters.at("slopeDividedBy1000")).at(5) << " $\\pm$ " << (fitParameters.at("slopeErrorDividedBy1000")).at(5) << " \\\\ \\hline" << std::endl;
  std::cout << std::fixed << std::setprecision(3) << "  nJets $\\geq$ 6 & " << (fitParameters.at("slopeDividedBy1000")).at(6) << " $\\pm$ " << (fitParameters.at("slopeErrorDividedBy1000")).at(6) << " \\\\ \\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;

  for (int fit_type_index = fitTypeFirst; fit_type_index < static_cast<int>(fitType::nFitTypes); ++fit_type_index) {
    fitType fit_type = static_cast<fitType>(fit_type_index);
    TCanvas outputCanvas = TCanvas(("c_chiSqPerNDFs_fitType" + fitTypeNames.at(fit_type)).c_str(), ("c_chiSqPerNDFs_fitType" + fitTypeNames.at(fit_type)).c_str(), 1024, 768);
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
    
  TCanvas outputCanvas = TCanvas("c_STDistributions", "c_STDistributions", 1024, 768);
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
    TLegendEntry *legendEntry = legend.AddEntry((STDistributions.at(nJetsBin)), legendText.c_str());
    legendEntry->SetMarkerColor(colors.at(nJetsBin));
    legendEntry->SetLineColor(colors.at(nJetsBin));
    legendEntry->SetTextColor(colors.at(nJetsBin));
  }

  legend.Draw();
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  outputCanvas.Update();
  outputCanvas.SaveAs((options.outputFolder + "/STDistributions_" + options.yearString + "_" + options.identifier + "_" + options.selection + ".pdf").c_str());

  std::cout << "All done!" << std::endl;

  return EXIT_SUCCESS;
}
