#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <map>
#include <cassert>

#include "TROOT.h"
#include "TDirectory.h"
#include "Rtypes.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TMatrixT.h"
#include "TMatrixDfwd.h"
#include "TMatrixTSym.h"
#include "TMatrixDSymfwd.h"
#include "TMatrixDSymEigen.h"
#include "TVectorT.h"
#include "TVectorDfwd.h"
#include "RooMsgService.h"
#include "RooGlobalFunc.h"
#include "RooCmdArg.h"
#include "RooBinning.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooArgSet.h"
#include "RooAbsArg.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooKeysPdf.h"
#include "RooGenericPdf.h"
#include "RooProdPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooArgList.h"

#include "tmArgumentParser.h"
#include "tmProgressBar.h"

#include "../../eventSelection/include/STRegionsStruct.h"

using namespace RooFit;

#define CHECK_TOLERANCE 0.001

#define ST_MAX_RANGE 3500.0

struct optionsStruct {
  std::string sourceFilePath, outputFolder, selection, identifier, yearString, inputParametersFileName;
  double preNormalizationBuffer;
  STRegionsStruct STRegions;
  double STNormTarget; // found implicitly from STRegions
  int PDF_nSTBins;
  bool readParametersFromFile;

  friend std::ostream& operator<< (std::ostream& out, const optionsStruct& options) {
    out << "sourceFilePath: " << options.sourceFilePath << std::endl
	<< "outputFolder: " << options.outputFolder << std::endl
        << "selection: " << options.selection << std::endl
        << "identifier: " << options.identifier << std::endl
        << "yearString: " << options.yearString << std::endl
        << "STRegions: " << options.STRegions << std::endl
        << "STNormTarget: " << options.STNormTarget << std::endl
        << "PDF_nSTBins: " << options.PDF_nSTBins << std::endl
        << "preNormalizationBuffer: " << options.preNormalizationBuffer << std::endl
        << "readParametersFromFile: " << (options.readParametersFromFile? "true": "false") << std::endl
        << "inputParametersFileName: " << options.inputParametersFileName << std::endl;
    return out;
  }
};

optionsStruct getOptionsFromParser(tmArgumentParser& argumentParser) {
  optionsStruct options = optionsStruct();
  options.sourceFilePath = argumentParser.getArgumentString("sourceFilePath");
  options.outputFolder = argumentParser.getArgumentString("outputFolder");
  options.selection = argumentParser.getArgumentString("selection");
  options.identifier = argumentParser.getArgumentString("identifier");
  options.yearString = argumentParser.getArgumentString("yearString");
  std::string STBoundariesSourceFile = argumentParser.getArgumentString("STBoundariesSourceFile");
  options.PDF_nSTBins = std::stoi(argumentParser.getArgumentString("PDF_nSTBins"));
  options.STRegions = STRegionsStruct(STBoundariesSourceFile, ST_MAX_RANGE);
  options.STNormTarget = 0.5*(options.STRegions.STNormRangeMin + options.STRegions.STNormRangeMax);
  options.preNormalizationBuffer = std::stod(argumentParser.getArgumentString("preNormalizationBuffer"));
  options.inputParametersFileName = argumentParser.getArgumentString("readParametersFromFile");
  options.readParametersFromFile = (options.inputParametersFileName != "/dev/null");
  return options;
}

struct eigenvalue_eigenvector_pair_struct {
  double eigenvalue;
  std::vector<double> eigenvector;

  eigenvalue_eigenvector_pair_struct() {
    eigenvalue = 0.;
    assert(static_cast<int>(eigenvector.size()) == static_cast<int>(0));
  }

  eigenvalue_eigenvector_pair_struct(const double& eigenvalue_, const std::vector<double>& eigenvector_) {
    eigenvalue = eigenvalue_;
    eigenvector = eigenvector_;
  }
};
