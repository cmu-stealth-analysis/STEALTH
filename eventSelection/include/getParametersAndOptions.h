#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

#include "tmArgumentParser.h"
#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"
#include "TNamed.h"

std::string getNDashes(const int& n) {
  std::stringstream dashes;
  for (int counter = 0; counter < n; ++counter) dashes << "-";
  return dashes.str();
}

struct optionsStruct {
  std::string inputFilePath, optionsStringPath, parametersStringPath;
};
