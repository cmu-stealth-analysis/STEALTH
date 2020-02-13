#ifndef H_HLTEMULATION
#define H_HLTEMULATION

#include <iostream>

#include "miscDataStructures.h"
#include "objectProperties.h"

namespace hltEmulation{
  bool passesHLTEmulation(const int& year, const triggerType& trigger_type, photonProperties& properties_leadingPhoton, photonProperties& properties_subLeadingPhoton, const float& eventHT, const int& triggerBit) {
    /* std::cout << "Called for year: " << year << ", trigger_type: " << triggerTypeNames.at(trigger_type) << std::endl; */
    /* std::cout << "properties_leadingPhoton: " << std::endl; */
    /* print_photon_properties(properties_leadingPhoton); */
    /* std::cout << "properties_subLeadingPhoton: " << std::endl; */
    /* print_photon_properties(properties_subLeadingPhoton); */
    /* std::cout << "eventHT: " << eventHT << ", triggerBit: " << triggerBit << std::endl; */
    bool passesEmulation = false;
    switch(year) {
    case 2016:
      {
	if (trigger_type == triggerType::photon) {
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
	  else if (triggerBit == 7) {
	    float& leading_photon_pT = properties_leadingPhoton[photonProperty::pT];
	    float& leading_photon_HOverE = properties_leadingPhoton[photonProperty::hOverE];
	    passesEmulation = ((leading_photon_pT > 175.) &&
			       (leading_photon_HOverE <= 0.15));
	  }
	  else if (triggerBit == 1) {
	    float& leading_photon_pT = properties_leadingPhoton[photonProperty::pT];
	    float& leading_photon_HOverE = properties_leadingPhoton[photonProperty::hOverE];
	    passesEmulation = ((leading_photon_pT > 30.) &&
			       (leading_photon_HOverE <= 0.15));
	  }
	  else if (triggerBit == 21) {
	    passesEmulation = true;
	  }
	  else {
	    std::cout << "ERROR: unsupported trigger bit: " << triggerBit << std::endl;
	    std::exit(EXIT_FAILURE);
	  }
	  break;
	}
	else if (trigger_type == triggerType::jet) {
	  if (triggerBit == 33) {
	    passesEmulation = (eventHT > 900.);
	  }
	  else {
	    std::cout << "ERROR: unsupported trigger bit: " << triggerBit << std::endl;
	    std::exit(EXIT_FAILURE);
	  }
	  break;
	}
	else {
	  std::cout << "ERROR: unsupported trigger type!" << std::endl;
	  std::exit(EXIT_FAILURE);
	}
      }
    case 2017:
    case 2018:
      {
	if (trigger_type == triggerType::photon) {
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
	  else if (triggerBit == 2) {
	    float& leading_photon_pT = properties_leadingPhoton[photonProperty::pT];
	    float& leading_photon_HOverE = properties_leadingPhoton[photonProperty::hOverE];
	    passesEmulation = ((leading_photon_pT > 30.) &&
			       (leading_photon_HOverE <= 0.15));
	  }
	  else if (triggerBit == 10) {
	    float& leading_photon_pT = properties_leadingPhoton[photonProperty::pT];
	    float& leading_photon_HOverE = properties_leadingPhoton[photonProperty::hOverE];
	    passesEmulation = ((leading_photon_pT > 200.) &&
			       (leading_photon_HOverE <= 0.15));
	  }
	  else {
	    std::cout << "ERROR: unsupported trigger bit: " << triggerBit << std::endl;
	    std::exit(EXIT_FAILURE);
	  }
	  break;
	}
	else if (trigger_type == triggerType::jet) {
	  if (triggerBit == 37) {
	    passesEmulation = (eventHT > 1050.);
	  }
	  else {
	    std::cout << "ERROR: unsupported trigger bit: " << triggerBit << std::endl;
	    std::exit(EXIT_FAILURE);
	  }
	  break;
	}
	else {
	  std::cout << "ERROR: unsupported trigger type!" << std::endl;
	  std::exit(EXIT_FAILURE);
	}
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
