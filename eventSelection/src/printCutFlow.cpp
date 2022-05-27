#include "../include/printCutFlow.h"

void fillCutflowCountersFromFiles(const optionsStruct & options, std::map<iyear, cutflowCountersStruct> & counters_for_year) {
  for (iyear y = 2016; y <= 2018; ++y) {
    TFile *input_file_handle = TFile::Open(options.inputFilePaths.at(y).c_str(), "READ");
    assert((input_file_handle->IsOpen() && !(input_file_handle->IsZombie())));
    TH1D * counts = nullptr;
    input_file_handle->GetObject("counters", counts);
    assert(counts != nullptr);
    counters_for_year[y] = cutflowCountersStruct(counts);
    input_file_handle->Close();
    std::cout << "Read in cut-flow counters:" << std::endl;
    std::cout << counters_for_year.at(y) << std::endl;
  }
}

void combineCounters(const std::map<iyear, cutflowCountersStruct> & counters_for_year, cutflowCountersStruct & combined_counters) {
  for (iyear y = 2016; y <= 2018; ++y) {
    combined_counters.N_analyzed += (counters_for_year.at(y)).N_analyzed;
    combined_counters.N_selected += (counters_for_year.at(y)).N_selected;
    for (int criterionIndex = eventSelectionCriterionFirst; criterionIndex != static_cast<int>(eventSelectionCriterion::nEventSelectionCriteria); ++criterionIndex) {
      eventSelectionCriterion criterion = static_cast<eventSelectionCriterion>(criterionIndex);
      (combined_counters.N_passing_cut)[criterion] += ((counters_for_year.at(y)).N_passing_cut).at(criterion);
      (combined_counters.N_passing_all_cuts_upto)[criterion] += ((counters_for_year.at(y)).N_passing_all_cuts_upto).at(criterion);
      (combined_counters.N_passing_all_cuts_besides)[criterion] += ((counters_for_year.at(y)).N_passing_all_cuts_besides).at(criterion);
    }
  }
}

void printCutFlowToFile(const std::string & output_file_path, const cutflowCountersStruct & combined_counters) {
  std::ofstream output_stream(output_file_path.c_str());
  assert(output_stream.is_open());
  output_stream << combined_counters.N_analyzed << std::endl;
  output_stream << combined_counters.N_selected << std::endl;
  for (int criterionIndex = eventSelectionCriterionFirst; criterionIndex != static_cast<int>(eventSelectionCriterion::nEventSelectionCriteria); ++criterionIndex) {
    eventSelectionCriterion criterion = static_cast<eventSelectionCriterion>(criterionIndex);
    output_stream << eventSelectionCriterionNames.at(criterion) << "    " << (combined_counters.N_passing_cut).at(criterion) << "    " << (combined_counters.N_passing_all_cuts_upto).at(criterion) << "    " << (combined_counters.N_passing_all_cuts_besides).at(criterion) << std::endl;
  }
  output_stream.close();
}

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  TH1::AddDirectory(kFALSE);
  std::cout << "Saving cut flow..." << std::endl;
  optionsStruct options(argc, argv);
  std::cout << "Options: " << options << std::endl;
  std::map<iyear, cutflowCountersStruct> counters_for_year;
  fillCutflowCountersFromFiles(options, counters_for_year);
  cutflowCountersStruct combined_counters;
  combineCounters(counters_for_year, combined_counters);
  printCutFlowToFile(options.outputFilePath, combined_counters);
  std::cout << "Cut flow script finished successfully." << std::endl;
  return EXIT_SUCCESS;
}
