#ifndef H_HLTEMULATION
#define H_HLTEMULATION

#include "objectProperties.h"

namespace hltEmulation{
  bool passesHLTEmulation(const int& year, photonProperties& properties_leadingPhoton, photonProperties& properties_subLeadingPhoton, const int& triggerBit) {
    bool passesEmulation = false;
    switch(year) {
    case 2016:
      {
	if (triggerBit == 16) {
	  bool leadingPassesEmulation = false;
	  float& leading_photon_pT = properties_leadingPhoton[photonProperty::pT];
	  float& leading_photon_R9 = properties_leadingPhoton[photonProperty::R9];
	  float& leading_photon_HOverE = properties_leadingPhoton[photonProperty::hOverE];
	  float& leading_photon_sigmaIEtaIEta = properties_leadingPhoton[photonProperty::sigmaIEtaIEta];
	  float& leading_photon_ecalClusIso = properties_leadingPhoton[photonProperty::ecalClusIso];

	  if (leading_photon_R9 >= 0.85) {
	    leadingPassesEmulation = ((leading_photon_pT > 30.) &&
				      (leading_photon_HOverE <= 0.1));
	  }
	  else if ((leading_photon_R9 >= 0.5) && (leading_photon_R9 < 0.85)) {
	    leadingPassesEmulation = ((leading_photon_pT > 30.) &&
				      (leading_photon_HOverE <= 0.1) &&
				      (leading_photon_sigmaIEtaIEta <= 0.015) &&
				      (leading_photon_ecalClusIso <= (6.0 + 0.012*leading_photon_pT)));
	  }

	  bool subLeadingPassesEmulation = false;
	  float& subLeading_photon_pT = properties_subLeadingPhoton[photonProperty::pT];
	  float& subLeading_photon_R9 = properties_subLeadingPhoton[photonProperty::R9];
	  float& subLeading_photon_HOverE = properties_subLeadingPhoton[photonProperty::hOverE];
	  float& subLeading_photon_sigmaIEtaIEta = properties_subLeadingPhoton[photonProperty::sigmaIEtaIEta];
	  float& subLeading_photon_ecalClusIso = properties_subLeadingPhoton[photonProperty::ecalClusIso];
	  float& subLeading_photon_trkIso = properties_subLeadingPhoton[photonProperty::trkIso];

	  if (subLeading_photon_R9 >= 0.85) {
	    subLeadingPassesEmulation = ((subLeading_photon_pT > 18.) &&
					 (subLeading_photon_HOverE <= 0.1));
	  }
	  else if ((subLeading_photon_R9 >= 0.5) && (subLeading_photon_R9 < 0.85)) {
	    subLeadingPassesEmulation = ((subLeading_photon_pT > 18.) &&
					 (subLeading_photon_HOverE <= 0.1) &&
					 (subLeading_photon_sigmaIEtaIEta <= 0.015) &&
					 (subLeading_photon_ecalClusIso <= (6.0 + 0.012*subLeading_photon_pT)) &&
					 (subLeading_photon_trkIso <= (6.0 + 0.002*subLeading_photon_pT)));
	  }

	  passesEmulation = leadingPassesEmulation && subLeadingPassesEmulation;
	}
	else if (triggerBit == 22) {
	  bool leadingPassesEmulation = false;
	  float& leading_photon_pT = properties_leadingPhoton[photonProperty::pT];
	  float& leading_photon_HOverE = properties_leadingPhoton[photonProperty::hOverE];

	  leadingPassesEmulation = ((leading_photon_pT > 60.) &&
				    (leading_photon_HOverE <= 0.15));

	  bool subLeadingPassesEmulation = false;
	  float& subLeading_photon_pT = properties_subLeadingPhoton[photonProperty::pT];
	  float& subLeading_photon_HOverE = properties_subLeadingPhoton[photonProperty::hOverE];

	  subLeadingPassesEmulation = ((subLeading_photon_pT > 60.) &&
				       (subLeading_photon_HOverE) <= 0.15);

	  passesEmulation = leadingPassesEmulation && subLeadingPassesEmulation;
	}
	else {
	  std::cout << "ERROR: unsupported trigger bit: " << triggerBit << std::endl;
	  std::exit(EXIT_FAILURE);
	}
	break;
      }
    case 2017:
    case 2018:
      {
	if (triggerBit == 37) {
	  bool leadingPassesEmulation = false;
	  float& leading_photon_pT = properties_leadingPhoton[photonProperty::pT];
	  float& leading_photon_R9 = properties_leadingPhoton[photonProperty::R9];
	  float& leading_photon_HOverE = properties_leadingPhoton[photonProperty::hOverE];
	  float& leading_photon_sigmaIEtaIEta = properties_leadingPhoton[photonProperty::sigmaIEtaIEta];
	  float& leading_photon_ecalClusIso = properties_leadingPhoton[photonProperty::ecalClusIso];

	  if (leading_photon_R9 >= 0.85) {
	    leadingPassesEmulation = ((leading_photon_pT > 30.) &&
				      (leading_photon_HOverE <= 0.1));
	  }
	  else if ((leading_photon_R9 >= 0.5) && (leading_photon_R9 < 0.85)) {
	    leadingPassesEmulation = ((leading_photon_pT > 30.) &&
				      (leading_photon_HOverE <= 0.1) &&
				      (leading_photon_sigmaIEtaIEta <= 0.015) &&
				      (leading_photon_ecalClusIso <= (6.0 + 0.012*leading_photon_pT)));
	  }

	  bool subLeadingPassesEmulation = false;
	  float& subLeading_photon_pT = properties_subLeadingPhoton[photonProperty::pT];
	  float& subLeading_photon_R9 = properties_subLeadingPhoton[photonProperty::R9];
	  float& subLeading_photon_HOverE = properties_subLeadingPhoton[photonProperty::hOverE];
	  float& subLeading_photon_sigmaIEtaIEta = properties_subLeadingPhoton[photonProperty::sigmaIEtaIEta];
	  float& subLeading_photon_ecalClusIso = properties_subLeadingPhoton[photonProperty::ecalClusIso];
	  float& subLeading_photon_trkIso = properties_subLeadingPhoton[photonProperty::trkIso];

	  if (subLeading_photon_R9 >= 0.85) {
	    subLeadingPassesEmulation = ((subLeading_photon_pT > 18.) &&
					 (subLeading_photon_HOverE <= 0.1));
	  }
	  else if ((subLeading_photon_R9 >= 0.5) && (subLeading_photon_R9 < 0.85)) {
	    subLeadingPassesEmulation = ((subLeading_photon_pT > 18.) &&
					 (subLeading_photon_HOverE <= 0.1) &&
					 (subLeading_photon_sigmaIEtaIEta <= 0.015) &&
					 (subLeading_photon_ecalClusIso <= (6.0 + 0.012*subLeading_photon_pT)) &&
					 (subLeading_photon_trkIso <= (6.0 + 0.002*subLeading_photon_pT)));
	  }

	  passesEmulation = leadingPassesEmulation && subLeadingPassesEmulation;
	}
	else if (triggerBit == 22) {
	  bool leadingPassesEmulation = false;
	  float& leading_photon_pT = properties_leadingPhoton[photonProperty::pT];
	  float& leading_photon_HOverE = properties_leadingPhoton[photonProperty::hOverE];

	  leadingPassesEmulation = ((leading_photon_pT > 70.) &&
				    (leading_photon_HOverE <= 0.15));

	  bool subLeadingPassesEmulation = false;
	  float& subLeading_photon_pT = properties_subLeadingPhoton[photonProperty::pT];
	  float& subLeading_photon_HOverE = properties_subLeadingPhoton[photonProperty::hOverE];

	  subLeadingPassesEmulation = ((subLeading_photon_pT > 70.) &&
				       (subLeading_photon_HOverE) <= 0.15);

	  passesEmulation = leadingPassesEmulation && subLeadingPassesEmulation;
	}
	else {
	  std::cout << "ERROR: unsupported trigger bit: " << triggerBit << std::endl;
	  std::exit(EXIT_FAILURE);
	}
	break;
      }
    default:
      {
	std::cout << "ERROR: unsupported year: " << year << std::endl;
	std::exit(EXIT_FAILURE);
      }
    }
    return passesEmulation;
  }
}

#endif
