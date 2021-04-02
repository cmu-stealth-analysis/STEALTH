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
#include "TMath.h"
#include "Math/IntegratorOptions.h"

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
#define TF1_INTEGRAL_REL_TOLERANCE 1.e-4
#define ST_MAX_RANGE 3500.0

struct optionsStruct {
  std::string sourceFilePath, outputFolder, selection, identifier, yearString, inputUnbinnedParametersFileName, inputBinnedParametersFileName;
  double preNormalizationBuffer;
  STRegionsStruct STRegions;
  double STNormTarget; // found implicitly from STRegions
  int PDF_nSTBins;
  bool readParametersFromFiles, plotConcise;

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
        << "readParametersFromFiles: " << (options.readParametersFromFiles? "true": "false") << std::endl
        << "inputUnbinnedParametersFileName: " << options.inputUnbinnedParametersFileName << std::endl
        << "inputBinnedParametersFileName: " << options.inputBinnedParametersFileName << std::endl
        << "plotConcise: " << (options.plotConcise? "true": "false") << std::endl;
    return out;
  }
};

std::vector<std::string> getComponentsOfCommaSeparatedString(const std::string &inputString) {
  std::vector<std::string> components;
  std::stringstream runningComponent;
  const char commaCharacter = ',';
  for (unsigned int stringIndex = 0; stringIndex < static_cast<unsigned int>(inputString.size()); ++stringIndex) {
    const char &character = inputString.at(stringIndex);
    if (character == commaCharacter) {
      components.push_back(runningComponent.str());
      runningComponent.str(std::string());
      runningComponent.clear();
    }
    else {
      runningComponent << character;
    }
  }
  components.push_back(runningComponent.str());
  return components;
}

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
  std::string readParametersFromFilesRaw = argumentParser.getArgumentString("readParametersFromFiles");
  options.readParametersFromFiles = (readParametersFromFilesRaw != "/dev/null,/dev/null");
  std::vector<std::string> readParametersFromFiles_components = getComponentsOfCommaSeparatedString(readParametersFromFilesRaw);
  assert(static_cast<int>(readParametersFromFiles_components.size()) == 2);
  options.inputUnbinnedParametersFileName = readParametersFromFiles_components.at(0);
  options.inputBinnedParametersFileName = readParametersFromFiles_components.at(1);;
  std::string plotConciseRaw = argumentParser.getArgumentString("plotConcise");
  if (plotConciseRaw == "true") options.plotConcise = true;
  else if (plotConciseRaw == "false") options.plotConcise = false;
  else {
    std::cout << "ERROR: unrecognized value for argument plotConcise, needs to be \"true\" or \"false\". Currently, value: " << plotConciseRaw << std::endl;
    std::exit(EXIT_FAILURE);
  }
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

enum class customizationType{ScaleOnly=0, Slope, Sqrt, SlopeSqrt, SlopeSqrtQuad, nCustomizationTypes};
int customizationTypeFirst = static_cast<int>(customizationType::ScaleOnly);
std::map<customizationType, std::string> customizationTypeNames = {
  {customizationType::ScaleOnly, "scaled"},
  {customizationType::Slope, "linear"},
  {customizationType::Sqrt, "sqrt"},
  {customizationType::SlopeSqrt, "slope_sqrt"},
  {customizationType::SlopeSqrtQuad, "linear_sqrt_quad"}
};

void do_sanity_checks_customizationTypes() {
  assert(static_cast<int>(customizationTypeNames.size()) == static_cast<int>(customizationType::nCustomizationTypes));
}

class customizedPDF {
 public:
  RooAbsPdf* pdf;
  RooRealVar* var;
  double scale;
  double normTarget;
  customizationType customization_type;
  double (customizedPDF::*getPDFTimesAdjustmentsAt)(double, double *);

  void setScale(std::string targetRangeName, double targetIntegralValue) {
    assert(targetIntegralValue > 0.);
    double pdfIntegralOverTargetRange = (pdf->createIntegral(*var, NormSet(*var), Range(targetRangeName.c_str())))->getVal();
    scale = targetIntegralValue/pdfIntegralOverTargetRange;
  }

  double evaluatePDFAt(double x) {
    var->setVal(x);
    double fvalue = pdf->getVal(*var);
    return (fvalue);
  }

  double getSlopeAdjustmentAt(double x, double slope) {
    return (1.0 + slope*((x/normTarget) - 1.0));
  }

  double getSqrtAdjustmentAt(double x, double sqrtTerm) {
    return (1.0 + sqrtTerm*((std::sqrt(x/normTarget)) - 1.0));
  }

  double getQuadAdjustmentAt(double x, double quadTerm) {
    return (1.0 + quadTerm*((std::pow(x/normTarget, 2)) - 1.0));
  }

  double PDFTimesAdjustment_ScaleOnly(double x, double *p) {
    (void)p;
    return scale*evaluatePDFAt(x);
  }

  double PDFTimesAdjustment_Slope(double x, double *p) {
    return scale*evaluatePDFAt(x)*getSlopeAdjustmentAt(x, p[0]); /* p[0] is interpreted as the slope */
  }

  double PDFTimesAdjustment_Sqrt(double x, double *p) {
    return scale*evaluatePDFAt(x)*getSqrtAdjustmentAt(x, p[0]); /* p[0] is interpreted as the sqrt term */
  }

  double PDFTimesAdjustment_SlopeSqrt(double x, double *p) {
    return scale*evaluatePDFAt(x)*getSlopeAdjustmentAt(x, p[0])*getSqrtAdjustmentAt(x, p[1]); /* p[0] is interpreted as the slope, p[1] as the sqrt term */
  }

  double PDFTimesAdjustment_SlopeSqrtQuad(double x, double *p) {
    return scale*evaluatePDFAt(x)*getSlopeAdjustmentAt(x, p[0])*getSqrtAdjustmentAt(x, p[1])*getQuadAdjustmentAt(x, p[2]); /* p[0] is interpreted as the slope, p[1] as the sqrt term, p[2] as the quad term */
  }

  customizedPDF(RooAbsPdf* pdf_, RooRealVar* var_, double normTarget_, customizationType customization_type_) {
    pdf = pdf_;
    var = var_;
    scale = 1.0;
    normTarget = normTarget_;
    customization_type = customization_type_;
    switch(customization_type) {
    case customizationType::ScaleOnly:
      getPDFTimesAdjustmentsAt = &customizedPDF::PDFTimesAdjustment_ScaleOnly;
      break;
    case customizationType::Slope:
      getPDFTimesAdjustmentsAt = &customizedPDF::PDFTimesAdjustment_Slope;
      break;
    case customizationType::Sqrt:
      getPDFTimesAdjustmentsAt = &customizedPDF::PDFTimesAdjustment_Sqrt;
      break;
    case customizationType::SlopeSqrt:
      getPDFTimesAdjustmentsAt = &customizedPDF::PDFTimesAdjustment_SlopeSqrt;
      break;
    case customizationType::SlopeSqrtQuad:
      getPDFTimesAdjustmentsAt = &customizedPDF::PDFTimesAdjustment_SlopeSqrtQuad;
      break;
    default:
      std::cout << "ERROR: unexpected customization type" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  double operator()(double *x, double *p) {
    return (this->*getPDFTimesAdjustmentsAt)(x[0], p);
  }
};

struct goodnessOfFitStruct {
  double chi2;
  int ndf;

  goodnessOfFitStruct() {
    chi2 = 0.; ndf = 0;
  }

  goodnessOfFitStruct(double chi2_, int ndf_) {
    chi2 = chi2_; ndf = ndf_;
  }
};
