#include "../include/JEC.h"

int main(int argc, char* argv[]) {
  std::cout << "Hello world!" << std::endl;
  tmArgumentParser argumentParser = tmArgumentParser("Template parser.");
  argumentParser.addArgument("arg", "defaultValue", false, "argument example.");
  argumentParser.setPassedStringValues(argc, argv);

  std::string argValue = argumentParser.getArgumentString("arg");
  std::cout << "Argument \"arg\" is set to \"" << argValue << "\"" << std::endl;
  return 0;
}
