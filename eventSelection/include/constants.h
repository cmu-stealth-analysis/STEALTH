#ifndef H_CONSTANTS
#define H_CONSTANTS

#include <cstdlib>
#include <iostream>
#include "TMath.h"

namespace constants{ // for readability
  const float PI = static_cast<float>(TMath::Pi());
  const float TWOPI = static_cast<float>(2.0*TMath::Pi());
}

namespace PIDUtils {
  const int quark_d = 1;
  const int quark_u = 2;
  const int quark_s = 3;
  const int quark_c = 4;
  const int quark_b = 5;
  const int quark_t = 6;
  const int gluon = 21;
  const int photon = 22;
  const int higgs = 25;
  const int electron = 11;
  const int electron_neutrino = 12;
  const int muon = 13;
  const int muon_neutrino = 14;
  const int tau = 15;
  const int tau_neutrino = 16;
  const int WBoson = 24;
  const int mesonRangeMin = 100;
  const int mesonRangeMax = 999;
  const int baryonRangeMin = 1000;
  const int baryonRangeMax = 9999;
  const int gluino = 1000021;
  const int neutralino = 1000022;
  const int chargino = 1000024;
  const int gravitino = 1000039;
  const int singlino = 3000001;
  const int singlet = 3000002;
  const int unusual1Min = 10000;
  const int unusual1Max = 999999;
  const int unusual2Min = 4000000;
  const int unusual2Max = 9999999;

  bool isJetCandidatePID(const int& candidate_id) {
    int abs_candidate_id = std::abs(candidate_id);
    return ((abs_candidate_id == quark_d) ||
            (abs_candidate_id == quark_u) ||
            (abs_candidate_id == quark_s) ||
            (abs_candidate_id == quark_c) ||
            (abs_candidate_id == quark_b) ||
            (abs_candidate_id == quark_t) ||
            (abs_candidate_id == gluon));
  }
  bool isPhotonPID(const int& candidate_id) {
    return (candidate_id == photon);
  }
  bool isHiggsPID(const int& candidate_id) {
    return (candidate_id == higgs);
  }
  bool isLeptonPID(const int& candidate_id) {
    int abs_candidate_id = std::abs(candidate_id);
    return ((abs_candidate_id == electron) ||
            (abs_candidate_id == muon) ||
            (abs_candidate_id == tau));
  }
  bool isLeptonNeutrinoPID(const int& candidate_id) {
    return ((candidate_id == electron_neutrino) ||
            (candidate_id == muon_neutrino) ||
            (candidate_id == tau_neutrino));
  }
  bool isWBosonPID(const int& candidate_id) {
    int abs_candidate_id = std::abs(candidate_id);
    return (abs_candidate_id == WBoson);
  }
  bool isMesonPID(const int& candidate_id) {
    int abs_candidate_id = std::abs(candidate_id);
    return ((abs_candidate_id >= mesonRangeMin) &&
            (abs_candidate_id <= mesonRangeMax));
  }
  bool isBaryonPID(const int& candidate_id) {
    int abs_candidate_id = std::abs(candidate_id);
    return ((abs_candidate_id >= baryonRangeMin) &&
            (abs_candidate_id <= baryonRangeMax));
  }
  bool isGluinoPID(const int& candidate_id) {
    return (candidate_id == gluino);
  }
  bool isNeutralinoPID(const int& candidate_id) {
    return (candidate_id == neutralino);
  }
  bool isCharginoPID(const int& candidate_id) {
    int abs_candidate_id = std::abs(candidate_id);
    return (abs_candidate_id == chargino);
  }
  bool isGravitinoPID(const int& candidate_id) {
    return (candidate_id == gravitino);
  }
  bool isSinglinoPID(const int& candidate_id) {
    return (candidate_id == singlino);
  }
  bool isSingletPID(const int& candidate_id) {
    return (candidate_id == singlet);
  }
  bool isUnusualPID(const int& candidate_id) {
    int abs_candidate_id = std::abs(candidate_id);
    return (((abs_candidate_id >= unusual1Min) &&
             (abs_candidate_id <= unusual1Max)) ||
            ((abs_candidate_id >= unusual2Min) &&
             (abs_candidate_id <= unusual2Max)));
  }
  bool isGGNtuplizerKnownDefault(const int& candidate_id) {
    return ((candidate_id == -99) || (candidate_id == -999));
  }

  /* The default MCPIDs are not convenient to plot on a histogram.
     Instead we plot the ID defined by
     the following translation.*/
  int getCustomParticleID(const int& candidate_id) {
    if (isJetCandidatePID(candidate_id)) return 1;
    else if (isPhotonPID(candidate_id)) return 2;
    else if (isHiggsPID(candidate_id)) return 3;
    else if (isLeptonPID(candidate_id)) return 4;
    else if (isLeptonNeutrinoPID(candidate_id)) return 5;
    else if (isWBosonPID(candidate_id)) return 6;
    else if (isMesonPID(candidate_id)) return 7;
    else if (isBaryonPID(candidate_id)) return 8;
    else if (isGluinoPID(candidate_id)) return 9;
    else if (isNeutralinoPID(candidate_id)) return 10;
    else if (isCharginoPID(candidate_id)) return 11;
    else if (isGravitinoPID(candidate_id)) return 12;
    else if (isSinglinoPID(candidate_id)) return 13;
    else if (isSingletPID(candidate_id)) return 14;
    else if (isUnusualPID(candidate_id)) return 15;
    else if (isGGNtuplizerKnownDefault(candidate_id)) return 16;
    return 17;
  }
}

namespace HLTEmulation{
  /* const int nEtaBins = 30; */
  /* const float etaMin = -1.5; */
  /* const float etaMax = 1.5; */
  /* const int nPTBins = 20; */
  /* const float PTMin = 15.; */
  /* const float PTMax = 115.; */
  std::vector<double> etaBinEdges = {-1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5};
  std::vector<double> pTBinEdges = {10., 30., 50., 70., 90., 110., 130., 150., 170., 190., 210., 230., 250., 300., 350., 400., 450., 500., 600., 700., 800., 1000.};
}

#endif
