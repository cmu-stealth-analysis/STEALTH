#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

enum class shiftType{JECDown=0, JECUp, UnclusteredMETDown, UnclusteredMETUp, JERMETDown, JERMETUp, nShiftTypes};
int shiftTypeFirst = static_cast<int>(shiftType::JECDown);
std::map<shiftType, std::string> shiftTypeNames = {
  {shiftType::JECDown, "JECDown"},
  {shiftType::JECUp, "JECUp"},
  {shiftType::UnclusteredMETDown, "UnclusteredMETDown"},
  {shiftType::UnclusteredMETUp, "UnclusteredMETUp"},
  {shiftType::JERMETDown, "JERMETDown"},
  {shiftType::JERMETUp, "JERMETUp"}
};

std::string getShiftedVariableBranchName(shiftType typeIndex, std::string variable) {
  std::string prefix = "b_";
  return (prefix + variable + "_shifted_" + shiftTypeNames[typeIndex]);
}
