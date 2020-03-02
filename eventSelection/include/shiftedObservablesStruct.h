#ifndef H_SHIFTEDOBSERVABLESSTRUCT
#define H_SHIFTEDOBSERVABLESSTRUCT

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

enum class shiftType{JECDown=0, JECUp, UnclusteredMETDown, UnclusteredMETUp, JERMETDown, JERMETUp, missingHEMDown, missingHEMUp, nShiftTypes};
int shiftTypeFirst = static_cast<int>(shiftType::JECDown);
std::map<shiftType, std::string> shiftTypeNames = {
  {shiftType::JECDown, "JECDown"},
  {shiftType::JECUp, "JECUp"},
  {shiftType::UnclusteredMETDown, "UnclusteredMETDown"},
  {shiftType::UnclusteredMETUp, "UnclusteredMETUp"},
  {shiftType::JERMETDown, "JERMETDown"},
  {shiftType::JERMETUp, "JERMETUp"},
  {shiftType::missingHEMDown, "missingHEMDown"},
  {shiftType::missingHEMUp, "missingHEMUp"}
};

std::string getShiftedVariableBranchName(shiftType typeIndex, std::string variable) {
  std::string prefix = "b_";
  return (prefix + variable + "_shifted_" + shiftTypeNames[typeIndex]);
}

template<typename printableType>
void printShiftedVariablesMap(const std::map<shiftType, printableType>& mapToPrint, const std::string& mapName) {
  for (int shiftTypeIndex = shiftTypeFirst; shiftTypeIndex != static_cast<int>(shiftType::nShiftTypes); ++shiftTypeIndex) {
    shiftType typeIndex = static_cast<shiftType>(shiftTypeIndex);
    std::cout << mapName << "[" << shiftTypeNames[typeIndex] << "]: " << mapToPrint.at(typeIndex) << ", ";
  }
  std::cout << "end." << std::endl;
}

#endif
