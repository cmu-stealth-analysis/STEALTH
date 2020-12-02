#ifndef H_SELECTIONREGIONUTILS
#define H_SELECTIONREGIONUTILS

#include "selectionCriteria.h"

struct selectionRegionDetailsStruct{
  int indexLeadingPhoton = -1;
  int indexSubLeadingPhoton = -1;
  selectionRegion selection_region = selectionRegion::nSelectionRegions;
};

namespace selectionRegionUtils{
  selectionRegionDetailsStruct getSelectionRegion(const bool& doSinglePhotonSelection, const int& n_mediumPhotons, const int& n_mediumPhotonsPassingLeadingPTCut, const std::map<int, int>& selectedMediumPhotonIndices, const int& n_vetoedPhotons, const int& n_vetoedPhotonsPassingLeadingPTCut, const std::map<int, int>& selectedVetoedPhotonIndices, const int& n_fakePhotons, const int& n_fakePhotonsPassingLeadingPTCut, const std::map<int, int>& selectedFakePhotonIndices, const std::map<int, float>& selectedPhotonPTs) {
    selectionRegionDetailsStruct selection_region_details;
    selection_region_details.selection_region = selectionRegion::nSelectionRegions;

    /* Check if there is only one photon (medium, loose, or fake in that order), only if "singlephoton" selection is requested*/
    if (doSinglePhotonSelection) {
      if (n_mediumPhotons == 1) {
        if ((n_mediumPhotonsPassingLeadingPTCut >= 1)) {
          selection_region_details.indexLeadingPhoton = selectedMediumPhotonIndices.at(0);
          selection_region_details.indexSubLeadingPhoton = -1;
          selection_region_details.selection_region = selectionRegion::control_singlemedium;
          return selection_region_details;
        }
      }
      if (n_vetoedPhotons == 1) {
        if ((n_vetoedPhotonsPassingLeadingPTCut >= 1)) {
          selection_region_details.indexLeadingPhoton = selectedVetoedPhotonIndices.at(0);
          selection_region_details.indexSubLeadingPhoton = -1;
          selection_region_details.selection_region = selectionRegion::control_singleloose;
          return selection_region_details;
        }
      }
      if (n_fakePhotons == 1) {
        if ((n_fakePhotonsPassingLeadingPTCut >= 1)) {
          selection_region_details.indexLeadingPhoton = selectedFakePhotonIndices.at(0);
          selection_region_details.indexSubLeadingPhoton = -1;
          selection_region_details.selection_region = selectionRegion::control_singlefake;
          return selection_region_details;
        }
      }
      selection_region_details.indexLeadingPhoton = -1;
      selection_region_details.indexSubLeadingPhoton = -1;
      selection_region_details.selection_region = selectionRegion::nSelectionRegions;
      return selection_region_details;
    }

    /* A selected photon can be medium(1), vetoed(2), or fake(3). */
    /* Combinations making up the selection regions: */
    /* 1 + 1: signal */
    /* 1 + 2, 2 + 1, 2 + 2: signal_loose */
    /* 3 + 3: control */
    if (n_mediumPhotons == 2) {
      if (n_mediumPhotonsPassingLeadingPTCut >= 1) { /* 1 + 1 */
        selection_region_details.indexLeadingPhoton = selectedMediumPhotonIndices.at(0);
        selection_region_details.indexSubLeadingPhoton = selectedMediumPhotonIndices.at(1);
        selection_region_details.selection_region = selectionRegion::signal;
        return selection_region_details;
      }
      /* n_mediumPhotons = 2 but both are subleading, check if there is a vetoed or fake candidate for leading photon */
      if (n_vetoedPhotonsPassingLeadingPTCut >= 1) { /* 2 + 1 */
        selection_region_details.indexLeadingPhoton = selectedVetoedPhotonIndices.at(0);
        selection_region_details.indexSubLeadingPhoton = selectedMediumPhotonIndices.at(0);
        selection_region_details.selection_region = selectionRegion::signal_loose;
        return selection_region_details;
      }
    }
    if (n_mediumPhotons == 1) {
      if (n_mediumPhotonsPassingLeadingPTCut >= 1) {
        float pT_medium = selectedPhotonPTs.at(selectedMediumPhotonIndices.at(0));
        /* We already have a leading medium photon candidate, now search for the other candidate */
        if (n_vetoedPhotons >= 1) {
          float pT_vetoed = selectedPhotonPTs.at(selectedVetoedPhotonIndices.at(0));
          if (pT_medium >= pT_vetoed) { /* 1 + 2 */
            selection_region_details.indexLeadingPhoton = selectedMediumPhotonIndices.at(0);
            selection_region_details.indexSubLeadingPhoton = selectedVetoedPhotonIndices.at(0);
            selection_region_details.selection_region = selectionRegion::signal_loose;
            return selection_region_details;
          }
          else { /* 2 + 1 */
            selection_region_details.indexLeadingPhoton = selectedVetoedPhotonIndices.at(0);
            selection_region_details.indexSubLeadingPhoton = selectedMediumPhotonIndices.at(0);
            selection_region_details.selection_region = selectionRegion::signal_loose;
            return selection_region_details;
          }
        }
      }
      else {
        /* We have a medium photon but it is subleading, the other one must be a vetoed leading candidate */
        if (n_vetoedPhotonsPassingLeadingPTCut >= 1) { /* 2 + 1 */
          selection_region_details.indexLeadingPhoton = selectedVetoedPhotonIndices.at(0);
          selection_region_details.indexSubLeadingPhoton = selectedMediumPhotonIndices.at(0);
          selection_region_details.selection_region = selectionRegion::signal_loose;
          return selection_region_details;
        }
      }
    }

    /* Remaining case for the signal_loose region: 2 + 2 */
    if ((n_vetoedPhotons >= 2) && (n_vetoedPhotonsPassingLeadingPTCut >= 1)) { /* 2 + 2 */
      selection_region_details.indexLeadingPhoton = selectedVetoedPhotonIndices.at(0);
      selection_region_details.indexSubLeadingPhoton = selectedVetoedPhotonIndices.at(1);
      selection_region_details.selection_region = selectionRegion::signal_loose;
      return selection_region_details;
    }

    /* Check if the event belongs to the double fake control region*/
    if ((n_fakePhotons >= 2) && (n_fakePhotonsPassingLeadingPTCut >= 1)) { /* 3 + 3 */
      selection_region_details.indexLeadingPhoton = selectedFakePhotonIndices.at(0);
      selection_region_details.indexSubLeadingPhoton = selectedFakePhotonIndices.at(1);
      selection_region_details.selection_region = selectionRegion::control_fakefake;
      return selection_region_details;
    }

    selection_region_details.indexLeadingPhoton = -1;
    selection_region_details.indexSubLeadingPhoton = -1;
    selection_region_details.selection_region = selectionRegion::nSelectionRegions;
    return selection_region_details;
  }
}

#endif
