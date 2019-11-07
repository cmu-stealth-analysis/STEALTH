#ifndef H_MCTEMPLATEREADER
#define H_MCTEMPLATEREADER

#include <vector>
#include <map>
#include <cassert>

#include "TFile.h"
#include "TH2F.h"
#include "TAxis.h"

class MCTemplateReader {
 public:
  const float nEventsFractionThreshold = 0.1;
  int nEventProgenitorMassBins, nNeutralinoMassBins;
  float minEventProgenitorMass, maxEventProgenitorMass, minNeutralinoMass, maxNeutralinoMass;
  float maxNEvents; // TH2F stores floats, to allow for weighted events
  std::map<int, float> eventProgenitorMasses; // map from eventProgenitor mass bin index to eventProgenitor mass
  std::map<int, float> neutralinoMasses; // map from neutralino mass bin index to neutralino mass
  std::map<int, std::map<int, float> > generated_nEvents; // first index: eventProgenitor mass bin, second index: neutralino mass bin; stores the number of events in a particular bin, useful for events weights

  MCTemplateReader(const std::string& templateSourceFilePath) {
    TFile* templateFile;
    templateFile = TFile::Open(templateSourceFilePath.c_str(), "READ");
    assert(templateFile != nullptr);
    TH2F* h_template;
    templateFile->GetObject("h_masses", h_template);
    assert(h_template != nullptr);
    nEventProgenitorMassBins = h_template->GetXaxis()->GetNbins();
    minEventProgenitorMass = h_template->GetXaxis()->GetXmin();
    maxEventProgenitorMass = h_template->GetXaxis()->GetXmax();
    nNeutralinoMassBins = h_template->GetYaxis()->GetNbins();
    minNeutralinoMass = h_template->GetYaxis()->GetXmin();
    maxNeutralinoMass = h_template->GetYaxis()->GetXmax();

    maxNEvents = -1;
    for (int eventProgenitorBinIndex = 1; eventProgenitorBinIndex <= nEventProgenitorMassBins; ++eventProgenitorBinIndex) {
      float xBinCenter = h_template->GetXaxis()->GetBinCenter(eventProgenitorBinIndex);
      eventProgenitorMasses[eventProgenitorBinIndex] = xBinCenter;
      for (int neutralinoBinIndex = 1; neutralinoBinIndex <= nNeutralinoMassBins; ++neutralinoBinIndex) {
	float yBinCenter = h_template->GetYaxis()->GetBinCenter(neutralinoBinIndex);
	neutralinoMasses[neutralinoBinIndex] = yBinCenter;
	float binContent = h_template->GetBinContent(eventProgenitorBinIndex, neutralinoBinIndex);
	// std::cout << "At (eventProgenitorMass, neutralinoMass) = (" << xBinCenter << ", " << yBinCenter << "), templateContents: " << binContent << std::endl;
	generated_nEvents[eventProgenitorBinIndex][neutralinoBinIndex] = binContent;
	if ((maxNEvents < 0) || (maxNEvents < binContent)) maxNEvents = binContent;
      }
    }
    templateFile->Close();
  }

  float getTotalNEvents(const int& eventProgenitorBinIndex, const int& neutralinoBinIndex) {
    return ((generated_nEvents.at(eventProgenitorBinIndex)).at(neutralinoBinIndex));
  }

  bool isValidBin(const int& eventProgenitorBinIndex, const int& neutralinoBinIndex) {
    float nEvents = getTotalNEvents(eventProgenitorBinIndex, neutralinoBinIndex);
    return (nEvents > nEventsFractionThreshold*maxNEvents);
  }

  void test() {
    for (int eventProgenitorBinIndex = 1; eventProgenitorBinIndex <= nEventProgenitorMassBins; ++eventProgenitorBinIndex) {
      for (int neutralinoBinIndex = 1; neutralinoBinIndex <= nNeutralinoMassBins; ++neutralinoBinIndex) {
	if (isValidBin(eventProgenitorBinIndex, neutralinoBinIndex)) {
	  std::cout << "Found valid bin at (eventProgenitorMass, neutralinoMass): (" << eventProgenitorMasses.at(eventProgenitorBinIndex) << ", " << neutralinoMasses.at(neutralinoBinIndex) << "); number of events = " << getTotalNEvents(eventProgenitorBinIndex, neutralinoBinIndex) << std::endl;
	}
      }
    }
  }
};

#endif
