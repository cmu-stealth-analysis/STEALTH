#ifndef H_CONSTANTS
#define H_CONSTANTS

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

  /* TODO: set these */
  /* const int singlet; */
  /* const int singlino; */
  bool isJetCandidatePID(const int& candidate_id) {
    return ((candidate_id == quark_d) ||
            (candidate_id == quark_u) ||
            (candidate_id == quark_s) ||
            (candidate_id == quark_c) ||
            (candidate_id == quark_b) ||
            (candidate_id == quark_t) ||
            (candidate_id == gluon));
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

  /* TODO: set these */
  /* bool isSingletPID(const int& candidate_id) { */
  /*   return (candidate_id == singlet); */
  /* } */
  /* bool isSinglinoPID(const int& candidate_id) { */
  /*   return (candidate_id == singlino); */
  /* } */
  /* bool isStealthParticlePID(const int& candidate_id) { */
  /*   return (isGluinoPID(candidate_id) || */
  /*           isNeutralinoPID(candidate_id) || */
  /*           isGravitinoPID(candidate_id) || */
  /*           isSingletPID(candidate_id) || */
  /*           isSinglinoPID(candidate_id)); */
  /* } */
}

#endif
