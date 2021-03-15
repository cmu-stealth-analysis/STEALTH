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
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include "tmArgumentParser.h"
#include "tmProgressBar.h"

#include "../../../eventSelection/include/STRegionsStruct.h"

struct optionsStruct {
  std::string sourceFilePath, outputFolder, selection, identifier, yearString;
  double STNormTarget, STNormMax;
  STRegionsStruct STRegions;

  friend std::ostream& operator<< (std::ostream& out, const optionsStruct& options) {
    out << "sourceFilePath: " << options.sourceFilePath << std::endl
	<< "outputFolder: " << options.outputFolder << std::endl
        << "selection: " << options.selection << std::endl
        << "identifier: " << options.identifier << std::endl
        << "yearString: " << options.yearString << std::endl
        << "STRegions: " << options.STRegions << std::endl
        << "STNormTarget: " << options.STNormTarget << std::endl
        << "STNormMax: " << options.STNormMax << std::endl;
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
  options.STRegions = STRegionsStruct(STBoundariesSourceFile, 3500.0);
  options.STNormTarget = std::stod(argumentParser.getArgumentString("STNormTarget"));
  options.STNormMax = std::stod(argumentParser.getArgumentString("STNormMax"));
  return options;
}

std::map<int, EColor> colors = {{2, static_cast<EColor>(kBlack)},
                                {3, static_cast<EColor>(kBlue+2)},
                                {4, static_cast<EColor>(kRed+1)},
                                {5, static_cast<EColor>(kGreen+3)},
                                {6, static_cast<EColor>(kViolet)}};


enum class fitType{fitConst=0, fitLin, fitQuad, fitConstrainedLin, nFitTypes};
int fitTypeFirst = static_cast<int>(fitType::fitConst);
std::map<fitType, std::string> fitTypeNames = {
  {fitType::fitConst, "const"},
  {fitType::fitLin, "lin"},
  {fitType::fitQuad, "quad"},
  {fitType::fitConstrainedLin, "constrained_lin"}
};
// fitFunctions has to be set at runtime because norm bin is set at runtime
std::map<fitType, std::string> chiSqPerNDFGraphTitles = {
  {fitType::fitConst, "#chi^{2}/NDF, constant fit"},
  {fitType::fitLin, "#chi^{2}/NDF, linear fit"},
  {fitType::fitQuad, "#chi^{2}/NDF, quadratic fit"},
  {fitType::fitConstrainedLin, "#chi^{2}/NDF, constrained linear fit"}
};

void do_sanity_checks_fitTypes() {
  assert(static_cast<int>(fitTypeNames.size()) == static_cast<int>(fitType::nFitTypes));
  assert(static_cast<int>(chiSqPerNDFGraphTitles.size()) == static_cast<int>(fitType::nFitTypes));
}

struct TGraphErrorsPointStruct {
  double x_val, y_val, x_err, y_err;

  TGraphErrorsPointStruct(double x_val_, double y_val_, double x_err_, double y_err_) {
    x_val = x_val_;
    y_val = y_val_;
    x_err = x_err_;
    y_err = y_err_;
  }

  TGraphErrorsPointStruct() {
    x_val = 0.;
    y_val = 0.;
    x_err = 0.;
    y_err = 0.;
  }
};
