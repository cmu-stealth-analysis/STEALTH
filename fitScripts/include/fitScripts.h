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
#include "TRandom3.h"

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
#define N_FLUCTUATIONS_TO_PLOT 100
#define FLUCTUATIONS_TRANSPARENCY 0.15

std::vector<std::string> splitStringByCharacter(const std::string & inputString, const char & split_character) {
  std::vector<std::string> components;
  std::stringstream runningComponent;
  // const char commaCharacter = ',';
  for (unsigned int stringIndex = 0; stringIndex < static_cast<unsigned int>(inputString.size()); ++stringIndex) {
    const char &character = inputString.at(stringIndex);
    if (character == split_character) {
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


struct sourceDataStruct {
  std::string sourceFilePath;
  bool PUReweightingNeeded;
  std::string PUWeightsPath;

  sourceDataStruct(const std::string & init_string) {
    const char exclamation_mark_character = '!';
    std::vector<std::string> init_string_split = splitStringByCharacter(init_string, exclamation_mark_character);
    if (static_cast<int>(init_string_split.size()) == 1) {
      sourceFilePath = init_string_split.at(0);
      PUReweightingNeeded = false;
      PUWeightsPath = "";
    }
    else if (static_cast<int>(init_string_split.size()) == 2) {
      sourceFilePath = init_string_split.at(0);
      PUReweightingNeeded = true;
      PUWeightsPath = init_string_split.at(1);;
    }
    else {
      std::cout << "ERROR: Tried to initialize sourceDataStruct in unrecognized format: " << init_string << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  friend std::ostream& operator<< (std::ostream& out, const sourceDataStruct & source_data) {
    out << "sourceFilePath: " << source_data.sourceFilePath << std::endl;
    out << "PUReweightingNeeded: " << (source_data.PUReweightingNeeded ? "true" : "false") << std::endl;
    if (source_data.PUReweightingNeeded) out << "PUWeightsPath: " << source_data.PUWeightsPath << std::endl;
    return out;
  }
};

struct optionsStruct {
  std::vector<sourceDataStruct> sourceData;
  std::string outputFolder, selection, identifier, yearString, /* inputUnbinnedParametersFileName,  */inputBinnedParametersFileName;
  double rhoNominal, preNormalizationBuffer, minAllowedEMST;
  STRegionsStruct STRegions, STRegions_for_ratio_wrt_chosen_adjustment;
  double STNormTarget; // found implicitly from STRegions
  int nJetsNorm, PDF_nSTBins;
  bool fetchMCWeights, getJECShiftedDistributions, readParametersFromFiles, plotConcise;

  friend std::ostream& operator<< (std::ostream& out, const optionsStruct& options) {
    out << "sourceData: " << std::endl;
    for (int source_data_index = 0; source_data_index < static_cast<int>((options.sourceData).size()); ++source_data_index) {
      out << "At index: " << source_data_index << std::endl;
      out << (options.sourceData).at(source_data_index);
    }
    out << "outputFolder: " << options.outputFolder << std::endl
        << "selection: " << options.selection << std::endl
        << "identifier: " << options.identifier << std::endl
        << "nJetsNorm: " << options.nJetsNorm << std::endl
        << "yearString: " << options.yearString << std::endl
        << "STRegions: " << options.STRegions << std::endl
        << "STNormTarget: " << options.STNormTarget << std::endl
        << "PDF_nSTBins: " << options.PDF_nSTBins << std::endl
        << "rhoNominal: " << options.rhoNominal << std::endl
        << "preNormalizationBuffer: " << options.preNormalizationBuffer << std::endl
        << "minAllowedEMST: " << options.minAllowedEMST << std::endl
        << "fetchMCWeights: " << (options.fetchMCWeights? "true": "false") << std::endl
        << "getJECShiftedDistributions: " << (options.getJECShiftedDistributions? "true": "false") << std::endl
        << "readParametersFromFiles: " << (options.readParametersFromFiles? "true": "false") << std::endl;
    if (options.readParametersFromFiles) {
      out /* << "inputUnbinnedParametersFileName: " << options.inputUnbinnedParametersFileName << std::endl */
          << "inputBinnedParametersFileName: " << options.inputBinnedParametersFileName << std::endl
          << "STRegions_for_ratio_wrt_chosen_adjustment: " << options.STRegions_for_ratio_wrt_chosen_adjustment << std::endl;
    }
    out << "plotConcise: " << (options.plotConcise? "true": "false") << std::endl; 
    return out;
  }
};

optionsStruct getOptionsFromParser(tmArgumentParser& argumentParser) {
  optionsStruct options = optionsStruct();
  std::string sourceDataRaw = argumentParser.getArgumentString("sourceData");
  const char comma_character = ',';
  std::vector<std::string> sourceDataRawSplit = splitStringByCharacter(sourceDataRaw, comma_character);
  for (int source_data_index = 0; source_data_index < static_cast<int>(sourceDataRawSplit.size()); ++source_data_index) {
    (options.sourceData).push_back(sourceDataStruct(sourceDataRawSplit.at(source_data_index)));
  }
  assert((options.sourceData).size() >= 1);
  options.outputFolder = argumentParser.getArgumentString("outputFolder");
  options.selection = argumentParser.getArgumentString("selection");
  std::string fetchMCWeightsRaw = argumentParser.getArgumentString("fetchMCWeights");
  if (fetchMCWeightsRaw == "true") options.fetchMCWeights = true;
  else if (fetchMCWeightsRaw == "false") options.fetchMCWeights = false;
  else {
    std::cout << "ERROR: unrecognized value for argument fetchMCWeights, needs to be \"true\" or \"false\". Currently, value: " << fetchMCWeightsRaw << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string getJECShiftedDistributionsRaw = argumentParser.getArgumentString("getJECShiftedDistributions");
  if (getJECShiftedDistributionsRaw == "true") options.getJECShiftedDistributions = true;
  else if (getJECShiftedDistributionsRaw == "false") options.getJECShiftedDistributions = false;
  else {
    std::cout << "ERROR: unrecognized value for argument getJECShiftedDistributions, needs to be \"true\" or \"false\". Currently, value: " << getJECShiftedDistributionsRaw << std::endl;
    std::exit(EXIT_FAILURE);
  }
  options.identifier = argumentParser.getArgumentString("identifier");
  options.nJetsNorm = std::stoi(argumentParser.getArgumentString("nJetsNorm"));
  options.yearString = argumentParser.getArgumentString("yearString");
  std::string STBoundariesSourceFile = argumentParser.getArgumentString("STBoundariesSourceFile");
  options.PDF_nSTBins = std::stoi(argumentParser.getArgumentString("PDF_nSTBins"));
  options.STRegions = STRegionsStruct(STBoundariesSourceFile, ST_MAX_RANGE);
  options.STNormTarget = 0.5*(options.STRegions.STNormRangeMin + options.STRegions.STNormRangeMax);
  options.rhoNominal = std::stod(argumentParser.getArgumentString("rhoNominal"));
  options.preNormalizationBuffer = std::stod(argumentParser.getArgumentString("preNormalizationBuffer"));
  options.minAllowedEMST = std::stod(argumentParser.getArgumentString("minAllowedEMST"));
  std::string readParametersFromFilesRaw = argumentParser.getArgumentString("readParametersFromFiles");
  options.readParametersFromFiles = (readParametersFromFilesRaw != "/dev/null,/dev/null");
  if (options.readParametersFromFiles) {
    std::vector<std::string> readParametersFromFiles_components = splitStringByCharacter(readParametersFromFilesRaw, comma_character);
    assert(static_cast<int>(readParametersFromFiles_components.size()) == 2);
    /* options.inputUnbinnedParametersFileName = readParametersFromFiles_components.at(0); */
    options.inputBinnedParametersFileName = readParametersFromFiles_components.at(0);
    options.STRegions_for_ratio_wrt_chosen_adjustment = STRegionsStruct(readParametersFromFiles_components.at(1), ST_MAX_RANGE);
  }
  std::string plotConciseRaw = argumentParser.getArgumentString("plotConcise");
  if (plotConciseRaw == "true") options.plotConcise = true;
  else if (plotConciseRaw == "false") options.plotConcise = false;
  else {
    std::cout << "ERROR: unrecognized value for argument plotConcise, needs to be \"true\" or \"false\". Currently, value: " << plotConciseRaw << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return options;
}

enum class fitType_ratios_wrt_chosen_adjustment{Linear=0, Quad, nFitTypes};

namespace constants {
  std::string binnedFitOptions = "QSI0+";
  std::string binnedFitOptions_ratios_wrt_chosen_adjustment = "QS0+";
  fitType_ratios_wrt_chosen_adjustment fit_type_ratios_wrt_chosen_adjustment = fitType_ratios_wrt_chosen_adjustment::Linear;
}

struct eigenmode_struct {
  double eigenvalue;
  std::vector<double> eigenvector;

  eigenmode_struct() {
    eigenvalue = 0.;
    assert(static_cast<int>(eigenvector.size()) == static_cast<int>(0));
  }

  eigenmode_struct(const double& eigenvalue_, const std::vector<double>& eigenvector_) {
    eigenvalue = eigenvalue_;
    eigenvector = eigenvector_;
  }
};

enum class customizationType{ScaleOnly=0, Slope, Sqrt, SlopeSqrt, SlopeSqrtQuad, nCustomizationTypes};
int customizationTypeFirst = static_cast<int>(customizationType::ScaleOnly);
std::map<customizationType, std::string> customizationTypeNames = {
  {customizationType::ScaleOnly, "scaled"},
  {customizationType::Slope, "slope"},
  {customizationType::Sqrt, "sqrt"},
  {customizationType::SlopeSqrt, "slope_sqrt"},
  {customizationType::SlopeSqrtQuad, "slope_sqrt_quad"}
};
std::map<customizationType, std::string> customizationTypeHumanReadableNames = {
  {customizationType::ScaleOnly, "scale-only"},
  {customizationType::Slope, "scale + slope"},
  {customizationType::Sqrt, "scale + sqrt"},
  {customizationType::SlopeSqrt, "scale + slope + sqrt"},
  {customizationType::SlopeSqrtQuad, "scale + slope + sqrt + quad"}
};
std::map<customizationType, int> customizationTypeNPars = {
  {customizationType::ScaleOnly, 1},
  {customizationType::Slope, 2},
  {customizationType::Sqrt, 2},
  {customizationType::SlopeSqrt, 3},
  {customizationType::SlopeSqrtQuad, 4}
};
std::map<customizationType, std::map<int, std::string> > customizationTypeParameterLabels = {
  {customizationType::ScaleOnly, {{0, "scale"}}},
  {customizationType::Slope, {{0, "scale"}, {1, "slope"}}},
  {customizationType::Sqrt, {{0, "scale"}, {1, "sqrt"}}},
  {customizationType::SlopeSqrt, {{0, "scale"}, {1, "slope"}, {2, "sqrt"}}},
  {customizationType::SlopeSqrtQuad, {{0, "scale"}, {1, "slope"}, {2, "sqrt"}, {3, "quad"}}}
};
std::map<customizationType, bool> customizationTypeActiveInConciseWorkflow = {
  {customizationType::ScaleOnly, true},
  {customizationType::Slope, false},
  {customizationType::Sqrt, true},
  {customizationType::SlopeSqrt, false},
  {customizationType::SlopeSqrtQuad, false}
};
std::map<customizationType, bool> customizationTypePlotEigenfluctuations = {
  {customizationType::ScaleOnly, false},
  {customizationType::Slope, false},
  {customizationType::Sqrt, true},
  {customizationType::SlopeSqrt, false},
  {customizationType::SlopeSqrtQuad, false}
};
std::map<customizationType, EColor> customizationTypeColors = {
  {customizationType::ScaleOnly, static_cast<EColor>(kBlue)},
  {customizationType::Slope, static_cast<EColor>(kRed+1)},
  {customizationType::Sqrt, static_cast<EColor>(kGreen+3)},
  {customizationType::SlopeSqrt, static_cast<EColor>(kViolet)},
  {customizationType::SlopeSqrtQuad, static_cast<EColor>(kYellow+2)}
};
std::map<customizationType, std::string> customizationTypeLegendLabels = {
  {customizationType::ScaleOnly, "low nJets template, normalized"},
  {customizationType::Slope, "low nJets template + linear adjustment"},
  {customizationType::Sqrt, "low nJets template + sqrt adjustment"},
  {customizationType::SlopeSqrt, "low nJets template + (linear+sqrt) adjustment"},
  {customizationType::SlopeSqrtQuad, "low nJets template + (linear+sqrt+quad) adjustment"}
};
std::map<customizationType, std::string> customizationTypeRatioLegendLabels = {
  {customizationType::ScaleOnly, "dummy string, not required"},
  {customizationType::Slope, "linear adjustment"},
  {customizationType::Sqrt, "sqrt adjustment"},
  {customizationType::SlopeSqrt, "(linear+sqrt) adjustment"},
  {customizationType::SlopeSqrtQuad, "(linear+sqrt+quad) adjustment"}
};

void do_sanity_checks_customizationTypes() {
  int n_customization_types = static_cast<int>(customizationType::nCustomizationTypes);
  assert(static_cast<customizationType>(n_customization_types) == customizationType::nCustomizationTypes);
  assert(static_cast<int>(customizationTypeNames.size()) == n_customization_types);
  assert(static_cast<int>(customizationTypeHumanReadableNames.size()) == n_customization_types);
  assert(static_cast<int>(customizationTypeNPars.size()) == n_customization_types);
  assert(static_cast<int>(customizationTypeParameterLabels.size()) == n_customization_types);
  for (int customization_type_index = customizationTypeFirst; customization_type_index < static_cast<int>(customizationType::nCustomizationTypes); ++customization_type_index) {
    customizationType customization_type = static_cast<customizationType>(customization_type_index);
    assert(static_cast<int>((customizationTypeParameterLabels.at(customization_type)).size()) == customizationTypeNPars.at(customization_type));
  }
  assert(static_cast<int>(customizationTypeActiveInConciseWorkflow.size()) == n_customization_types);
  assert(static_cast<int>(customizationTypePlotEigenfluctuations.size()) == n_customization_types);
  assert(static_cast<int>(customizationTypeColors.size()) == n_customization_types);
  assert(static_cast<int>(customizationTypeLegendLabels.size()) == n_customization_types);
  assert(static_cast<int>(customizationTypeRatioLegendLabels.size()) == n_customization_types);
}

struct parameter_initialization_struct{
  std::string name;
  double initial_value, range_min, range_max;

  parameter_initialization_struct() {
    name = "default";
    initial_value = 0.;
    range_min = 0.;
    range_max = 0.;
  }

  parameter_initialization_struct(std::string name_, double initial_value_, double range_min_, double range_max_) {
    name = name_;
    initial_value = initial_value_;
    range_min = range_min_;
    range_max = range_max_;
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

class customizedPDF {
 public:
  RooAbsPdf* pdf;
  RooRealVar* var;
  double nominal_scale;
  double norm_target;
  customizationType customization_type;
  double (customizedPDF::*getPDFTimesAdjustmentsAt)(double, double *);

  void setNominalScale(std::string targetRangeName, double targetIntegralValue) {
    assert(targetIntegralValue > 0.);
    double pdfIntegralOverTargetRange = (pdf->createIntegral(*var, NormSet(*var), Range(targetRangeName.c_str())))->getVal();
    nominal_scale = targetIntegralValue/pdfIntegralOverTargetRange;
  }

  double evaluatePDFAt(double x) {
    var->setVal(x);
    double fvalue = pdf->getVal(*var);
    return (fvalue);
  }

  double getSlopeAdjustmentAt(double x, double slope) {
    return (slope*((x/norm_target) - 1.0));
  }

  double getSqrtAdjustmentAt(double x, double sqrtTerm) {
    return (sqrtTerm*((std::sqrt(x/norm_target)) - 1.0));
  }

  double getQuadAdjustmentAt(double x, double quadTerm) {
    return (quadTerm*((std::pow(x/norm_target, 2)) - 1.0));
  }

  double PDFTimesAdjustment_ScaleOnly(double x, double *p) {
    (void)p;
    return nominal_scale*evaluatePDFAt(x)*p[0]; // p[0] is the overall scale
  }

  double PDFTimesAdjustment_Slope(double x, double *p) {
    return nominal_scale*evaluatePDFAt(x)*(p[0] + getSlopeAdjustmentAt(x, p[1])); /* p[0] is the overall scale, p[1] is the slope */
  }

  double PDFTimesAdjustment_Sqrt(double x, double *p) {
    return nominal_scale*evaluatePDFAt(x)*(p[0] + getSqrtAdjustmentAt(x, p[1])); /* p[0] is the overall scale, p[1] is the sqrt term */
  }

  double PDFTimesAdjustment_SlopeSqrt(double x, double *p) {
    return nominal_scale*evaluatePDFAt(x)*(p[0] + getSlopeAdjustmentAt(x, p[1]) + getSqrtAdjustmentAt(x, p[2])); /* p[0] is the overall scale, p[1] is the slope, p[2] is the sqrt term */
  }

  double PDFTimesAdjustment_SlopeSqrtQuad(double x, double *p) {
    return nominal_scale*evaluatePDFAt(x)*(p[0] + getSlopeAdjustmentAt(x, p[1]) + getSqrtAdjustmentAt(x, p[2]) + getQuadAdjustmentAt(x, p[3])); /* p[0] is the overall scale, p[1] is the slope, p[2] is the sqrt term, p[3] is the quad term */
  }

  customizedPDF(RooAbsPdf* pdf_, RooRealVar* var_, double norm_target_, customizationType customization_type_) {
    pdf = pdf_;
    var = var_;
    nominal_scale = 1.0;
    norm_target = norm_target_;
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

  customizedPDF() {
    pdf = nullptr;
    var = nullptr;
    nominal_scale = -1.0;
    norm_target = -1.0;
    customization_type = customizationType::ScaleOnly;
    getPDFTimesAdjustmentsAt = &customizedPDF::PDFTimesAdjustment_ScaleOnly;
  }

  double operator()(double *x, double *p) {
    return (this->*getPDFTimesAdjustmentsAt)(x[0], p);
  }
};

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
void check_eigendecomposition(const eigenmode_struct& pair, const T& matrix_to_check, const bool& print_debug=false) {
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

int getNNonEmptyBins(TH1D& inputHistogram) {
  int n_nonempty_bins = 0;
  for (int binCounter = 1; binCounter <= static_cast<int>(inputHistogram.GetXaxis()->GetNbins()); ++binCounter) {
    if ((inputHistogram.GetBinContent(binCounter)) > 0.) ++n_nonempty_bins;
  }
  return n_nonempty_bins;
}

std::string get_parameter_name(const customizationType &customization_type, const int& parameter_index, const int &n_jets_bin) {
  return std::string(customizationTypeNames.at(customization_type) + "_fit_" + (customizationTypeParameterLabels.at(customization_type)).at(parameter_index) + "_" + std::to_string(n_jets_bin) + "JetsBin");
}

std::string get_eigencoefficient_name(const customizationType &customization_type, const int& eigen_index, const int& parameter_index, const int &n_jets_bin) {
  return std::string(customizationTypeNames.at(customization_type) + "_fit_eigenmode_" + std::to_string(eigen_index) + "_coefficient_" + (customizationTypeParameterLabels.at(customization_type)).at(parameter_index) + "_" + std::to_string(n_jets_bin) + "JetsBin");
}

std::string get_eigenerror_name(const customizationType &customization_type, const int& eigen_index, const int &n_jets_bin) {
  return std::string(customizationTypeNames.at(customization_type) + "_fit_eigenmode_" + std::to_string(eigen_index) + "_error_" + std::to_string(n_jets_bin) + "JetsBin");
}

struct fit_result_struct {
  double chi_sq;
  int ndf;
  double pvalue;
  std::map<int, double> best_fit_values;
  std::map<int, double> best_fit_errors;
  std::vector<eigenmode_struct> eigenmodes;

  fit_result_struct() {
    chi_sq = -1.0;
    ndf = -1;
    pvalue = -1.;
  }

  fit_result_struct(double chi_sq_, int ndf_, double pvalue_, std::map<int, double> best_fit_values_, std::map<int, double> best_fit_errors_, std::vector<eigenmode_struct> eigenmodes_) {
    chi_sq = chi_sq_;
    ndf = ndf_;
    pvalue = pvalue_;
    best_fit_values = best_fit_values_;
    best_fit_errors = best_fit_errors_;
    eigenmodes = eigenmodes_;
  }
};

class customizedTF1 {
 private:
  customizationType customization_type;

  double getChisquareWRTHistogram(TH1D& inputHistogram) {
    return inputHistogram.Chisquare(raw_TF1, "R");
  }

 public:
  TF1 *raw_TF1;
  fit_result_struct fit_result;

  customizedTF1(std::string prefix_, customizedPDF* basePDF_, double rangeMin_, double rangeMax_, customizationType customization_type_) {
    raw_TF1 = new TF1((prefix_ + "_" + customizationTypeNames.at(customization_type_) + "_TF1").c_str(), *basePDF_, rangeMin_, rangeMax_, customizationTypeNPars.at(customization_type_));
    customization_type = customization_type_;
  }

  customizedTF1() {
    raw_TF1 = nullptr;
    customization_type = customizationType::nCustomizationTypes;
  }

  ~customizedTF1() {
    delete raw_TF1;
  }

  double getTFIntegral(const double &min, const double &max) {
    return raw_TF1->Integral(min, max, TF1_INTEGRAL_REL_TOLERANCE);
  }

  void initializeParameters(const std::map<int, parameter_initialization_struct> &parameter_initialization_map) {
    for (int parameter_index = 0; parameter_index < customizationTypeNPars.at(customization_type); ++parameter_index) {
      const parameter_initialization_struct& parameter_initialization = parameter_initialization_map.at(parameter_index);
      raw_TF1->SetParName(parameter_index, (parameter_initialization.name).c_str());
      raw_TF1->SetParameter(parameter_index, parameter_initialization.initial_value);
      raw_TF1->SetParLimits(parameter_index, parameter_initialization.range_min, parameter_initialization.range_max);
    }
  }

  void setFitResultsFromSource(const std::map<std::string, double> &fitParametersBinned, const int& n_jets_bin) {
    double chi_sq = -1.0;
    int ndf = -1;
    double pvalue = -1.0;
    std::map<int, double> best_fit_values;
    std::map<int, double> best_fit_errors;
    std::vector<eigenmode_struct> eigenmodes;
    for (int parameter_index = 0; parameter_index < customizationTypeNPars.at(customization_type); ++parameter_index) {
      best_fit_values[parameter_index] = fitParametersBinned.at(get_parameter_name(customization_type, parameter_index, n_jets_bin));
      best_fit_errors[parameter_index] = -1.; // we don't really need these errors
    }
    for (int eigen_index = 0; eigen_index < customizationTypeNPars.at(customization_type); ++eigen_index) {
      double eigenvalue;
      std::vector<double> eigenvector;
      for (int parameter_index = 0; parameter_index < customizationTypeNPars.at(customization_type); ++parameter_index) {
        eigenvector.push_back(fitParametersBinned.at(get_eigencoefficient_name(customization_type, eigen_index, parameter_index, n_jets_bin)));
      }
      eigenvalue = std::pow(fitParametersBinned.at(get_eigenerror_name(customization_type, eigen_index, n_jets_bin)), 2);
      eigenmode_struct eigenmode = eigenmode_struct(eigenvalue, eigenvector);
      eigenmodes.push_back(eigenmode);
    }
    fit_result = fit_result_struct(chi_sq, ndf, pvalue, best_fit_values, best_fit_errors, eigenmodes);
  }

  void fitToTH1(TH1D& inputHistogram, bool print_verbose=true) {
    int n_parameters = customizationTypeNPars.at(customization_type);
    double chisquare;
    int ndf;
    double pvalue;
    std::map<int, double> best_fit_values;
    std::map<int, double> best_fit_errors;
    std::vector<eigenmode_struct> eigenmodes;
    if (n_parameters == 0) {
      // there's nothing to fit, just calculate the chisquare and be done with it
      chisquare = getChisquareWRTHistogram(inputHistogram);
      ndf = getNNonEmptyBins(inputHistogram);
      pvalue = 1.0;
      fit_result = fit_result_struct(chisquare, ndf, pvalue, best_fit_values, best_fit_errors, eigenmodes);
      return;
    }
    TFitResultPtr root_fit_result_ptr = inputHistogram.Fit(raw_TF1, (constants::binnedFitOptions).c_str());
    assert(root_fit_result_ptr->Status() == 0);
    assert(static_cast<int>(root_fit_result_ptr->NTotalParameters()) == n_parameters);

    // step 1: get covariance matrix
    TMatrixDSym covarianceMatrix = root_fit_result_ptr->GetCovarianceMatrix();
    if (print_verbose) {
      std::cout << "For customization type: " << customizationTypeNames.at(customization_type) << ", covarianceMatrix: ";
      printSquareMatrix(covarianceMatrix, n_parameters);
    }

    // step 2: get eigendecomposition
    TMatrixDSymEigen eigendecomposition_setup = TMatrixDSymEigen(covarianceMatrix);
    TVectorD eigenvalues = eigendecomposition_setup.GetEigenValues();
    if (print_verbose) {
      std::cout << "eigenvalues: ";
      printTVector(eigenvalues);
    }
    TMatrixD eigenvectors = eigendecomposition_setup.GetEigenVectors();
    if (print_verbose) {
      std::cout << "eigenvectors: ";
      printSquareMatrix(eigenvectors, n_parameters);
    }
    for (int parameter_index = 0; parameter_index < n_parameters; ++parameter_index) {
      best_fit_values[parameter_index] = root_fit_result_ptr->Parameter(parameter_index);
      best_fit_errors[parameter_index] = root_fit_result_ptr->ParError(parameter_index);
    }
    for (int eigen_index = 0; eigen_index < n_parameters; ++eigen_index) {
      double eigenvalue = eigenvalues(eigen_index);
      std::vector<double> eigenvector = getColumnFromTMatrixD(eigenvectors, eigen_index, n_parameters);
      eigenmode_struct current_pair = eigenmode_struct(eigenvalue, eigenvector);
      check_eigendecomposition(current_pair, covarianceMatrix);
      eigenmodes.push_back(current_pair);
    }
    chisquare = root_fit_result_ptr->Chi2();
    ndf = root_fit_result_ptr->Ndf();
    pvalue = root_fit_result_ptr->Prob();
    if (print_verbose) {
      std::cout << "Fit chi^2: " << chisquare << ", ndf: " << ndf << ", pvalue: " << pvalue << std::endl;
    }
    fit_result = fit_result_struct(chisquare, ndf, pvalue, best_fit_values, best_fit_errors, eigenmodes);
  }

  void set_TF_parameters_to_nominal() {
    for (int parameter_index = 0; parameter_index < customizationTypeNPars.at(customization_type); ++parameter_index) {
      raw_TF1->SetParameter(parameter_index, fit_result.best_fit_values.at(parameter_index));
    }
  }

  void set_TF_parameters_to_eigenmode_fluctuation(const int &eigenmode_index, const double &fluctuation_nsigmas) {
    eigenmode_struct& fit_eigenmode = (fit_result.eigenmodes).at(eigenmode_index);
    for (int parameter_index = 0; parameter_index < customizationTypeNPars.at(customization_type); ++parameter_index) {
      double parameter_value = fit_result.best_fit_values.at(parameter_index) + fluctuation_nsigmas*std::sqrt(fit_eigenmode.eigenvalue)*((fit_eigenmode.eigenvector).at(parameter_index));
      raw_TF1->SetParameter(parameter_index, parameter_value);
    }
  }

  void set_TF_parameters_to_eigenmode_fluctuation(const std::map<int, double> &fluctuation_nsigmas_map) {
    for (int parameter_index = 0; parameter_index < customizationTypeNPars.at(customization_type); ++parameter_index) {
      double parameter_value = fit_result.best_fit_values.at(parameter_index);
      for (int eigen_index = 0; eigen_index < customizationTypeNPars.at(customization_type); ++eigen_index) {
        eigenmode_struct& fit_eigenmode = (fit_result.eigenmodes).at(eigen_index);
        parameter_value += (fluctuation_nsigmas_map.at(eigen_index))*std::sqrt(fit_eigenmode.eigenvalue)*((fit_eigenmode.eigenvector).at(parameter_index));
      }
      raw_TF1->SetParameter(parameter_index, parameter_value);
    }
  }

  TGraph get_nominal_fit_as_TGraph(const int &nGraphPoints, const double &xMin, const double &xMax) {
    set_TF_parameters_to_nominal();
    TGraph outputGraph = TGraph(nGraphPoints);
    outputGraph.SetName(("graph_nominal_fit_" + customizationTypeNames.at(customization_type)).c_str());
    for (int xCounter = 0; xCounter <= nGraphPoints; ++xCounter) {
      double x = xMin + (1.0*xCounter/nGraphPoints)*(xMax-xMin);
      outputGraph.SetPoint(xCounter, x, std::max(0., raw_TF1->Eval(x)));
    }
    return outputGraph;
  }

  TGraph get_eigenmode_fluctuation_as_TGraph(const int &eigenmode_index, const double &fluctuation_nsigmas, const int &nGraphPoints, const double &xMin, const double &xMax) {
    assert(fluctuation_nsigmas != 0);
    assert(eigenmode_index < customizationTypeNPars.at(customization_type));
    set_TF_parameters_to_eigenmode_fluctuation(eigenmode_index, fluctuation_nsigmas);
    TGraph outputGraph = TGraph(nGraphPoints);
    std::string fluctuationType = (fluctuation_nsigmas < 0? "down" : "up");
    outputGraph.SetName(("graph_eigenmode_fluctuation_mode_" + std::to_string(eigenmode_index) + "_fluctuation_" + fluctuationType + "_" + std::to_string(fluctuation_nsigmas) + "_sigmas_" + customizationTypeNames.at(customization_type)).c_str());
    for (int xCounter = 0; xCounter <= nGraphPoints; ++xCounter) {
      double x = xMin + (1.0*xCounter/nGraphPoints)*(xMax-xMin);
      outputGraph.SetPoint(xCounter, x, std::max(0., raw_TF1->Eval(x)));
    }
    return outputGraph;
  }

  TGraph get_eigenmode_fluctuation_as_TGraph(const std::map<int, double> &fluctuation_nsigmas_map, const int &nGraphPoints, const double &xMin, const double &xMax, const int& fluctuation_index) {
    set_TF_parameters_to_eigenmode_fluctuation(fluctuation_nsigmas_map);
    TGraph outputGraph = TGraph(nGraphPoints);
    outputGraph.SetName(("graph_random_fluctuation_" + std::to_string(fluctuation_index) + "_" + customizationTypeNames.at(customization_type)).c_str());
    for (int xCounter = 0; xCounter <= nGraphPoints; ++xCounter) {
      double x = xMin + (1.0*xCounter/nGraphPoints)*(xMax-xMin);
      outputGraph.SetPoint(xCounter, x, std::max(0., raw_TF1->Eval(x)));
    }
    return outputGraph;
  }

  std::map<int, double> getBinIntegralsDividedByBinWidthFromTF1(const STRegionsStruct &regions) {
    std::map<int, double> bin_integrals_divided_by_bin_widths;
    for (int regionIndex = 1; regionIndex <= regions.STAxis.GetNbins(); ++regionIndex) {
      double bin_low_edge = regions.STAxis.GetBinLowEdge(regionIndex);
      double bin_up_edge = regions.STAxis.GetBinUpEdge(regionIndex);
      double integral = (raw_TF1->Integral(bin_low_edge, bin_up_edge, TF1_INTEGRAL_REL_TOLERANCE));
      double bin_width = regions.STAxis.GetBinWidth(regionIndex);
      bin_integrals_divided_by_bin_widths[regionIndex] = integral/bin_width;
    }
    return bin_integrals_divided_by_bin_widths;
  }

  double evaluate_TF_at(const double& x) {
    return raw_TF1->Eval(x);
  }
};

void format_TGraph_as_nominal_fit(TGraph &inputGraph, const customizationType &customization_type) {
  inputGraph.SetLineColor(customizationTypeColors.at(customization_type)); inputGraph.SetLineWidth(2);
}

void format_TGraph_as_fluctuation(TGraph &inputGraph, const customizationType &customization_type) {
  inputGraph.SetLineStyle(kDashed); inputGraph.SetLineColor(customizationTypeColors.at(customization_type)); inputGraph.SetLineWidth(1);
}

void format_TGraph_as_random_eigenfluctuation(TGraph &inputGraph, const customizationType &customization_type) {
  inputGraph.SetLineColorAlpha(customizationTypeColors.at(customization_type), FLUCTUATIONS_TRANSPARENCY); inputGraph.SetLineWidth(3);
}

void set_legend_entry_color(TLegendEntry *legendEntry, const customizationType &customization_type) {
  legendEntry->SetMarkerColor(customizationTypeColors.at(customization_type)); legendEntry->SetLineColor(customizationTypeColors.at(customization_type)); legendEntry->SetTextColor(customizationTypeColors.at(customization_type));
}

void format_ratio_TGraph_as_nominal_and_add_to_multigraph(TGraph &inputGraph, TMultiGraph &inputMultigraph, TLegend &inputLegend, const customizationType &customization_type) {
  inputGraph.SetLineColor(customizationTypeColors.at(customization_type)); inputGraph.SetLineWidth(2); inputGraph.SetDrawOption("C"); inputMultigraph.Add(&inputGraph);
  TLegendEntry *legendEntry = inputLegend.AddEntry(&inputGraph, (customizationTypeRatioLegendLabels.at(customization_type)).c_str());
  set_legend_entry_color(legendEntry, customization_type);
}

void format_ratio_TGraph_as_fluctuation_and_add_to_multigraph(TGraph &inputGraph, TMultiGraph &inputMultigraph, const customizationType &customization_type) {
  inputGraph.SetLineColor(customizationTypeColors.at(customization_type)); inputGraph.SetLineStyle(kDashed); inputGraph.SetLineWidth(1); inputGraph.SetDrawOption("C"); inputMultigraph.Add(&inputGraph);
}

void format_ratio_TGraph_as_random_fluctuation_and_add_to_multigraph(TGraph &inputGraph, TMultiGraph &inputMultigraph, const customizationType &customization_type) {
  inputGraph.SetLineColorAlpha(customizationTypeColors.at(customization_type), FLUCTUATIONS_TRANSPARENCY); inputGraph.SetLineWidth(3); inputGraph.SetDrawOption("C"); inputMultigraph.Add(&inputGraph);
}

template<typename T>
std::string get_string_precision_n(const int & precision, const T & streamable_source) {
  std::stringstream ss;
  ss << std::setprecision(precision) << streamable_source << std::fixed;
  return (ss.str());
}
