#ifndef H_HLTEMULATION
#define H_HLTEMULATION

#include <iostream>

#include "miscDataStructures.h"
#include "objectProperties.h"

namespace hltEmulation{
  bool passesLeading_2016_photon_16(const photonProperties & properties) {
    bool passes = false;
    const float & leading_photon_pT = properties.at(photonProperty::pT);
    const float & leading_photon_R9 = properties.at(photonProperty::R9);
    const float & leading_photon_HOverE = properties.at(photonProperty::hOverE);
    const float & leading_photon_sigmaIEtaIEta = properties.at(photonProperty::sigmaIEtaIEta);
    const float & leading_photon_ecalClusIso = properties.at(photonProperty::ecalClusIso);

    if (leading_photon_R9 >= 0.85) {
      passes = ((leading_photon_pT > 30.) &&
                (leading_photon_HOverE <= 0.1));
    }
    else if ((leading_photon_R9 >= 0.5) && (leading_photon_R9 < 0.85)) {
      passes = ((leading_photon_pT > 30.) &&
                (leading_photon_HOverE <= 0.1) &&
                (leading_photon_sigmaIEtaIEta <= 0.015) &&
                (leading_photon_ecalClusIso <= (6.0 + 0.012*leading_photon_pT)));
    }
    return passes;
  }

  bool passesSubLeading_2016_photon_16(const photonProperties & properties) {
    bool passes = false;
    const float & subLeading_photon_pT = properties.at(photonProperty::pT);
    const float & subLeading_photon_R9 = properties.at(photonProperty::R9);
    const float & subLeading_photon_HOverE = properties.at(photonProperty::hOverE);
    const float & subLeading_photon_sigmaIEtaIEta = properties.at(photonProperty::sigmaIEtaIEta);
    const float & subLeading_photon_ecalClusIso = properties.at(photonProperty::ecalClusIso);
    const float & subLeading_photon_trkIso = properties.at(photonProperty::trkIso);

    if (subLeading_photon_R9 >= 0.85) {
      passes = ((subLeading_photon_pT > 18.) &&
                                   (subLeading_photon_HOverE <= 0.1));
    }
    else if ((subLeading_photon_R9 >= 0.5) && (subLeading_photon_R9 < 0.85)) {
      passes = ((subLeading_photon_pT > 18.) &&
                                   (subLeading_photon_HOverE <= 0.1) &&
                                   (subLeading_photon_sigmaIEtaIEta <= 0.015) &&
                                   (subLeading_photon_ecalClusIso <= (6.0 + 0.012*subLeading_photon_pT)) &&
                                   (subLeading_photon_trkIso <= (6.0 + 0.002*subLeading_photon_pT)));
    }
    return passes;
  }

  bool passesLeading_2017_2018_photon_37(const photonProperties & properties) {
    bool passes = false;
    const float & leading_photon_pT = properties.at(photonProperty::pT);
    const float & leading_photon_R9 = properties.at(photonProperty::R9);
    const float & leading_photon_HOverE = properties.at(photonProperty::hOverE);
    const float & leading_photon_sigmaIEtaIEta = properties.at(photonProperty::sigmaIEtaIEta);
    const float & leading_photon_ecalClusIso = properties.at(photonProperty::ecalClusIso);

    if (leading_photon_R9 >= 0.85) {
      passes = ((leading_photon_pT > 30.) &&
                (leading_photon_HOverE <= 0.1));
    }
    else if ((leading_photon_R9 >= 0.5) && (leading_photon_R9 < 0.85)) {
      passes = ((leading_photon_pT > 30.) &&
                (leading_photon_HOverE <= 0.1) &&
                (leading_photon_sigmaIEtaIEta <= 0.015) &&
                (leading_photon_ecalClusIso <= (6.0 + 0.012*leading_photon_pT)));
    }
    return passes;
  }

  bool passesSubLeading_2017_2018_photon_37(const photonProperties & properties) {
    bool passes = false;
    const float & subLeading_photon_pT = properties.at(photonProperty::pT);
    const float & subLeading_photon_R9 = properties.at(photonProperty::R9);
    const float & subLeading_photon_HOverE = properties.at(photonProperty::hOverE);
    const float & subLeading_photon_sigmaIEtaIEta = properties.at(photonProperty::sigmaIEtaIEta);
    const float & subLeading_photon_ecalClusIso = properties.at(photonProperty::ecalClusIso);
    const float & subLeading_photon_trkIso = properties.at(photonProperty::trkIso);

    if (subLeading_photon_R9 >= 0.85) {
      passes = ((subLeading_photon_pT > 18.) &&
                (subLeading_photon_HOverE <= 0.1));
    }
    else if ((subLeading_photon_R9 >= 0.5) && (subLeading_photon_R9 < 0.85)) {
      passes = ((subLeading_photon_pT > 18.) &&
                (subLeading_photon_HOverE <= 0.1) &&
                (subLeading_photon_sigmaIEtaIEta <= 0.015) &&
                (subLeading_photon_ecalClusIso <= (6.0 + 0.012*subLeading_photon_pT)) &&
                (subLeading_photon_trkIso <= (6.0 + 0.002*subLeading_photon_pT)));
    }
    return passes;
  }

  bool passesHLTEmulation(const int& year, const triggerType& trigger_type, const photonPropertiesCollection & preselectedPhotonProperties, // const float& eventHT,
                          const int& triggerBit) {
    /* std::cout << "Called for year: " << year << ", trigger_type: " << triggerTypeNames.at(trigger_type) << std::endl; */
    /* std::cout << "properties_leadingPhoton: " << std::endl; */
    /* print_photon_properties(properties_leadingPhoton); */
    /* std::cout << "properties_subLeadingPhoton: " << std::endl; */
    /* print_photon_properties(properties_subLeadingPhoton); */
    /* std::cout << "eventHT: " << eventHT << ", triggerBit: " << triggerBit << std::endl; */
    // bool passesEmulation = false;
    if (year == 2016) {
      if (trigger_type == triggerType::photon) {
        if (triggerBit == 16) {
          int n_preselected_photons = static_cast<int>(preselectedPhotonProperties.size());
          if (n_preselected_photons < 2) return false;
          std::map<int, bool> passesLeading;
          std::map<int, bool> passesSubLeading;
          for (int i = 0; i < n_preselected_photons; ++i) {
            passesLeading[i] = passesLeading_2016_photon_16(preselectedPhotonProperties.at(i));
            passesSubLeading[i] = passesSubLeading_2016_photon_16(preselectedPhotonProperties.at(i));
          }
          for (int i1 = 0; i1 < n_preselected_photons; ++i1) {
            if (passesLeading.at(i1)) {
              for (int i2 = 0; i2 < n_preselected_photons; ++i2) {
                if (i1 == i2) continue;
                if (passesSubLeading.at(i2)) return true;
              }
            }
          }
          return false;
        }
        // else if (triggerBit == 22) {
        //   bool leadingPassesEmulation = false;
        //   float& leading_photon_pT = properties_leadingPhoton[photonProperty::pT];
        //   float& leading_photon_HOverE = properties_leadingPhoton[photonProperty::hOverE];

        //   leadingPassesEmulation = ((leading_photon_pT > 60.) &&
        //                             (leading_photon_HOverE <= 0.15));

        //   bool subLeadingPassesEmulation = false;
        //   float& subLeading_photon_pT = properties_subLeadingPhoton[photonProperty::pT];
        //   float& subLeading_photon_HOverE = properties_subLeadingPhoton[photonProperty::hOverE];

        //   subLeadingPassesEmulation = ((subLeading_photon_pT > 60.) &&
        //                                (subLeading_photon_HOverE) <= 0.15);

        //   passesEmulation = leadingPassesEmulation && subLeadingPassesEmulation;
        // }
        // else if (triggerBit == 7) {
        //   float& leading_photon_pT = properties_leadingPhoton[photonProperty::pT];
        //   float& leading_photon_HOverE = properties_leadingPhoton[photonProperty::hOverE];
        //   passesEmulation = ((leading_photon_pT > 175.) &&
        //                      (leading_photon_HOverE <= 0.15));
        // }
        // else if (triggerBit == 1) {
        //   float& leading_photon_pT = properties_leadingPhoton[photonProperty::pT];
        //   float& leading_photon_HOverE = properties_leadingPhoton[photonProperty::hOverE];
        //   passesEmulation = ((leading_photon_pT > 30.) &&
        //                      (leading_photon_HOverE <= 0.15));
        // }
        // else if (triggerBit == 21) {
        //   passesEmulation = true;
        // }
        // else {
        //   std::cout << "ERROR: unsupported trigger bit: " << triggerBit << std::endl;
        //   std::exit(EXIT_FAILURE);
        // }
        // break;
      }
      // else if (trigger_type == triggerType::jet) {
      //   if (triggerBit == 33) {
      //     passesEmulation = (eventHT > 900.);
      //   }
      //   else {
      //     std::cout << "ERROR: unsupported trigger bit: " << triggerBit << std::endl;
      //     std::exit(EXIT_FAILURE);
      //   }
      //   break;
      // }
      else {
        std::cout << "ERROR: unsupported trigger type!" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    else if ((year == 2017) || (year == 2018)) {
      if (trigger_type == triggerType::photon) {
        if (triggerBit == 37) {
          int n_preselected_photons = static_cast<int>(preselectedPhotonProperties.size());
          if (n_preselected_photons < 2) return false;
          std::map<int, bool> passesLeading;
          std::map<int, bool> passesSubLeading;
          for (int i = 0; i < n_preselected_photons; ++i) {
            passesLeading[i] = passesLeading_2017_2018_photon_37(preselectedPhotonProperties.at(i));
            passesSubLeading[i] = passesSubLeading_2017_2018_photon_37(preselectedPhotonProperties.at(i));
          }
          for (int i1 = 0; i1 < n_preselected_photons; ++i1) {
            if (passesLeading.at(i1)) {
              for (int i2 = 0; i2 < n_preselected_photons; ++i2) {
                if (i1 == i2) continue;
                if (passesSubLeading.at(i2)) return true;
              }
            }
          }
          return false;
        }
        //   else if (triggerBit == 22) {
        //     bool leadingPassesEmulation = false;
        //     float& leading_photon_pT = properties_leadingPhoton[photonProperty::pT];
        //     float& leading_photon_HOverE = properties_leadingPhoton[photonProperty::hOverE];

        //     leadingPassesEmulation = ((leading_photon_pT > 70.) &&
        //                               (leading_photon_HOverE <= 0.15));

        //     bool subLeadingPassesEmulation = false;
        //     float& subLeading_photon_pT = properties_subLeadingPhoton[photonProperty::pT];
        //     float& subLeading_photon_HOverE = properties_subLeadingPhoton[photonProperty::hOverE];

        //     subLeadingPassesEmulation = ((subLeading_photon_pT > 70.) &&
        //                                  (subLeading_photon_HOverE) <= 0.15);

        //     passesEmulation = leadingPassesEmulation && subLeadingPassesEmulation;
        //   }
        //   else if (triggerBit == 2) {
        //     float& leading_photon_pT = properties_leadingPhoton[photonProperty::pT];
        //     float& leading_photon_HOverE = properties_leadingPhoton[photonProperty::hOverE];
        //     passesEmulation = ((leading_photon_pT > 30.) &&
        //                        (leading_photon_HOverE <= 0.15));
        //   }
        //   else if (triggerBit == 10) {
        //     float& leading_photon_pT = properties_leadingPhoton[photonProperty::pT];
        //     float& leading_photon_HOverE = properties_leadingPhoton[photonProperty::hOverE];
        //     passesEmulation = ((leading_photon_pT > 200.) &&
        //                        (leading_photon_HOverE <= 0.15));
        //   }
        //   else {
        //     std::cout << "ERROR: unsupported trigger bit: " << triggerBit << std::endl;
        //     std::exit(EXIT_FAILURE);
        //   }
        //   break;
      }
      // else if (trigger_type == triggerType::jet) {
      //   if (triggerBit == 37) {
      //     passesEmulation = (eventHT > 1050.);
      //   }
      //   else {
      //     std::cout << "ERROR: unsupported trigger bit: " << triggerBit << std::endl;
      //     std::exit(EXIT_FAILURE);
      //   }
      //   break;
      // }
      else {
        std::cout << "ERROR: unsupported trigger type!" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    else {
      std::cout << "ERROR: unsupported year: " << year << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::cout << "ERROR: Control should not reach here!" << std::endl;
    return false;
  }
}

#endif
