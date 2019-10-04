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
  int nGluinoMassBins, nNeutralinoMassBins;
  float minGluinoMass, maxGluinoMass, minNeutralinoMass, maxNeutralinoMass;
  float maxNEvents; // TH2F stores floats, to allow for weighted events
  std::map<int, float> gluinoMasses; // map from gluino mass bin index to gluino mass
  std::map<int, float> neutralinoMasses; // map from neutralino mass bin index to neutralino mass
  std::map<int, std::map<int, float> > generated_nEvents; // first index: gluino mass bin, second index: neutralino mass bin; stores the number of events in a particular bin, useful for events weights

  MCTemplateReader(const std::string& templateSourceFilePath) {
    TFile* templateFile;
    templateFile = TFile::Open(templateSourceFilePath.c_str(), "READ");
    assert(templateFile != nullptr);
    TH2F* h_template;
    templateFile->GetObject("h_masses", h_template);
    assert(h_template != nullptr);
    nGluinoMassBins = h_template->GetXaxis()->GetNbins();
    minGluinoMass = h_template->GetXaxis()->GetXmin();
    maxGluinoMass = h_template->GetXaxis()->GetXmax();
    nNeutralinoMassBins = h_template->GetYaxis()->GetNbins();
    minNeutralinoMass = h_template->GetYaxis()->GetXmin();
    maxNeutralinoMass = h_template->GetYaxis()->GetXmax();

    maxNEvents = -1;
    for (int gluinoBinIndex = 1; gluinoBinIndex <= nGluinoMassBins; ++gluinoBinIndex) {
      float xBinCenter = h_template->GetXaxis()->GetBinCenter(gluinoBinIndex);
      gluinoMasses[gluinoBinIndex] = xBinCenter;
      for (int neutralinoBinIndex = 1; neutralinoBinIndex <= nNeutralinoMassBins; ++neutralinoBinIndex) {
	float yBinCenter = h_template->GetYaxis()->GetBinCenter(neutralinoBinIndex);
	neutralinoMasses[neutralinoBinIndex] = yBinCenter;
	float binContent = h_template->GetBinContent(gluinoBinIndex, neutralinoBinIndex);
	// std::cout << "At (gluinoMass, neutralinoMass) = (" << xBinCenter << ", " << yBinCenter << "), templateContents: " << binContent << std::endl;
	generated_nEvents[gluinoBinIndex][neutralinoBinIndex] = binContent;
	if ((maxNEvents < 0) || (maxNEvents < binContent)) maxNEvents = binContent;
      }
    }
    templateFile->Close();
  }

  float getTotalNEvents(const int& gluinoBinIndex, const int& neutralinoBinIndex) {
    return ((generated_nEvents.at(gluinoBinIndex)).at(neutralinoBinIndex));
  }

  bool isValidBin(const int& gluinoBinIndex, const int& neutralinoBinIndex) {
    float nEvents = getTotalNEvents(gluinoBinIndex, neutralinoBinIndex);
    return (nEvents > nEventsFractionThreshold*maxNEvents);
  }

  void test() {
    for (int gluinoBinIndex = 1; gluinoBinIndex <= nGluinoMassBins; ++gluinoBinIndex) {
      for (int neutralinoBinIndex = 1; neutralinoBinIndex <= nNeutralinoMassBins; ++neutralinoBinIndex) {
	if (isValidBin(gluinoBinIndex, neutralinoBinIndex)) {
	  std::cout << "Found valid bin at (gluinoMass, neutralinoMass): (" << gluinoMasses.at(gluinoBinIndex) << ", " << neutralinoMasses.at(neutralinoBinIndex) << "); number of events = " << getTotalNEvents(gluinoBinIndex, neutralinoBinIndex) << std::endl;
	}
      }
    }
  }
};

#endif
