#ifndef H_CONSTANTS
#define H_CONSTANTS

#include <cstdlib>
#include "TMath.h"

namespace constants{ // for readability
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
  const int gluino = 1000021;
  const int neutralino = 1000022;
  const int gravitino = 1000039;
  const int singlino = 3000001;
  const int singlet = 3000002;

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
  bool isGluinoPID(const int& candidate_id) {
    return (candidate_id == gluino);
  }
  bool isNeutralinoPID(const int& candidate_id) {
    return (candidate_id == neutralino);
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
}

#endif
