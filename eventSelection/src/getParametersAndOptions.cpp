#include "../include/getParametersAndOptions.h"

optionsStruct getOptionsFromParser(tmArgumentParser& argumentParser) {
  optionsStruct options = optionsStruct();
  options.inputFilePath = argumentParser.getArgumentString("inputFilePath");
  options.optionsStringPath = argumentParser.getArgumentString("optionsStringPath");
  options.parametersStringPath = argumentParser.getArgumentString("parametersStringPath");
  return options;
}

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  tmArgumentParser argumentParser = tmArgumentParser("Calculate systematics due to uncertainty on jet energy corrections.");
  argumentParser.addArgument("inputFilePath", "", true, "Path to input file.");
  argumentParser.addArgument("optionsStringPath", "optionsString", false, "Name of TObject in which options are stored.");
  argumentParser.addArgument("parametersStringPath", "parametersString", false, "Name of TObject in which parameters are stored.");
  argumentParser.setPassedStringValues(argc, argv);
  optionsStruct options = getOptionsFromParser(argumentParser);
  TFile *inputFile = TFile::Open(options.inputFilePath.c_str(), "READ");
  TNamed *optionsStringObject = (TNamed*)(inputFile->Get(options.optionsStringPath.c_str()));
  std::string optionsString = std::string(optionsStringObject->GetTitle());
  std::cout << getNDashes(100) << std::endl
            << "Options:" << std::endl
            << optionsString << std::endl
            << getNDashes(100) << std::endl;
  TNamed *parametersStringObject = (TNamed*)(inputFile->Get(options.parametersStringPath.c_str()));
  std::string parametersString = std::string(parametersStringObject->GetTitle());
  std::cout << getNDashes(100) << std::endl
            << "Parameters:" << std::endl
            << parametersString << std::endl
            << getNDashes(100) << std::endl;
  inputFile->Close();
  return 0;
}
