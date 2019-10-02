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
  const float nEventsFractionThreshold = 0.5;
  int nGluinoMassBins, nNeutralinoMassBins;
  float minGluinoMass, maxGluinoMass, minNeutralinoMass, maxNeutralinoMass;
  float maxNEvents; // TH2F stores floats, to allow for weighted events
  TFile* templateFile;
  TH2F* h_template;
  std::map<int, int> gluinoMasses; // map from gluino mass bin index to gluino mass
  std::map<int, int> neutralinoMasses; // map from neutralino mass bin index to neutralino mass

  MCTemplateReader(const std::string& templateSourceFilePath, const std::string& templateName) {
    templateFile = TFile::Open(templateSourceFilePath.c_str(), "READ");
    assert(templateFile != nullptr);
    templateFile->GetObject(templateName.c_str(), h_template);
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
      int gluinoMass = static_cast<int>(0.5+xBinCenter);
      gluinoMasses[gluinoBinIndex] = gluinoMass;
      for (int neutralinoBinIndex = 1; neutralinoBinIndex <= nNeutralinoMassBins; ++neutralinoBinIndex) {
	float yBinCenter = h_template->GetYaxis()->GetBinCenter(neutralinoBinIndex);
	int neutralinoMass = static_cast<int>(0.5+yBinCenter);
	neutralinoMasses[neutralinoBinIndex] = neutralinoMass;
	float binContent = h_template->GetBinContent(gluinoBinIndex, neutralinoBinIndex);
	std::cout << "At (gluinoMass, neutralinoMass) = (" << xBinCenter << ", " << yBinCenter << "), templateContents: " << binContent << std::endl;
	if ((maxNEvents < 0) || (maxNEvents < binContent)) maxNEvents = binContent;
      }
    }
  }

  bool getTotalNEvents(const int& gluinoBinIndex, const int& neutralinoBinIndex) {
    return (h_template->GetBinContent(gluinoBinIndex, neutralinoBinIndex));
  }

  bool isValidBin(const int& gluinoBinIndex, const int& neutralinoBinIndex) {
    float nEvents = getTotalNEvents(gluinoBinIndex, neutralinoBinIndex);
    return (nEvents > nEventsFractionThreshold*maxNEvents);
  }

  ~MCTemplateReader() {
    templateFile->Close();
  }
};

#endif
