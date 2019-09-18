#ifndef H_MISCDATASTRUCTURES
#define H_MISCDATASTRUCTURES

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

#include "constants.h"

struct quadraticPolynomialStruct{
  float constCoefficient, linearCoefficient, squareCoefficient;
  quadraticPolynomialStruct () : constCoefficient(0.),
    linearCoefficient(0.),
    squareCoefficient(0.) {}

  quadraticPolynomialStruct (float constCoefficient_, float linearCoefficient_, float squareCoefficient_) : constCoefficient(constCoefficient_),
    linearCoefficient(linearCoefficient_),
    squareCoefficient(squareCoefficient_) {}
  float getPolynomialValue(float pT) {
    return (constCoefficient + linearCoefficient*pT + squareCoefficient*pT*pT);
  }
  friend std::ostream& operator<< (std::ostream& out, const quadraticPolynomialStruct& polynomial) {
    out << "const --> " << polynomial.constCoefficient << ", "
        << "linear--> " << polynomial.linearCoefficient << ", "
        << "square --> " << polynomial.squareCoefficient;
    return out;
  }
};

struct angularVariablesStruct{
  float eta, phi;
  angularVariablesStruct (float eta_, float phi_) : eta(eta_), phi(phi_) {}

  float get_deltaR(const angularVariablesStruct& angularVariables) {
    float deltaEta = angularVariables.eta - this->eta;
    float phi1 = angularVariables.phi;
    float phi2 = this->phi;
    if (phi2 > phi1) { // make sure phi1 > phi2
      phi1 = this->phi;
      phi2 = angularVariables.phi;
    }
    float deltaPhi_direction1 = phi1 - phi2;
    float deltaPhi_direction2 = constants::TWOPI - deltaPhi_direction1;
    float deltaPhi = (deltaPhi_direction1 < deltaPhi_direction2) ? deltaPhi_direction1 : deltaPhi_direction2;
    return std::sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
  }

  float getMinDeltaR(std::vector<angularVariablesStruct>& otherAngularVariables) {
    float min_dR = -0.005;
    for (auto&& angularVariablesToCompare: otherAngularVariables) {
      float dR = this->get_deltaR(angularVariablesToCompare);
      if ((min_dR < 0) || (dR < min_dR)) min_dR = dR;
    }
    return min_dR;
  }
};

struct photonQualityCutsStruct{
  float towerHOverE, sigmaIEtaIEta, sigmaIEtaIEtaLoose, chargedIsolation, chargedIsolationLoose;
  quadraticPolynomialStruct neutralIsolation, photonIsolation, photonIsolationLoose;

  photonQualityCutsStruct () : towerHOverE(0.),
    sigmaIEtaIEta(0.),
    sigmaIEtaIEtaLoose(0.),
    chargedIsolation(0.),
    chargedIsolationLoose(0.) {
    neutralIsolation = quadraticPolynomialStruct(0., 0., 0.);
    photonIsolation = quadraticPolynomialStruct(0., 0., 0.);
    photonIsolationLoose = quadraticPolynomialStruct(0., 0., 0.);
  }

  photonQualityCutsStruct (float towerHOverE_, float sigmaIEtaIEta_, float sigmaIEtaIEtaLoose_, float chargedIsolation_, float chargedIsolationLoose_, float neutralIsolationConst_, float neutralIsolationLinear_, float neutralIsolationSquare_, float photonIsolationConst_, float photonIsolationLinear_, float photonIsolationConstLoose_, float photonIsolationLinearLoose_) : towerHOverE(towerHOverE_),
    sigmaIEtaIEta(sigmaIEtaIEta_),
    sigmaIEtaIEtaLoose(sigmaIEtaIEtaLoose_),
    chargedIsolation(chargedIsolation_),
    chargedIsolationLoose(chargedIsolationLoose_)
  {
    neutralIsolation = quadraticPolynomialStruct(neutralIsolationConst_, neutralIsolationLinear_, neutralIsolationSquare_);
    photonIsolation = quadraticPolynomialStruct(photonIsolationConst_, photonIsolationLinear_, 0.);
    photonIsolationLoose = quadraticPolynomialStruct(photonIsolationConstLoose_, photonIsolationLinearLoose_, 0.);
  }

  friend std::ostream& operator<< (std::ostream& out, const photonQualityCutsStruct& cuts) {
    out << "towerHOverE: " << cuts.towerHOverE << ", "
        << "sigmaIEtaIEta: " << cuts.sigmaIEtaIEta << ", "
        << "sigmaIEtaIEtaLoose: " << cuts.sigmaIEtaIEtaLoose << ", "
        << "chargedIsolation: " << cuts.chargedIsolation << ", "
        << "chargedIsolationLoose: " << cuts.chargedIsolationLoose << ", "
        << "neutral isolation coefficients: " << cuts.neutralIsolation << ", "
        << "photon isolation coefficients: " << cuts.photonIsolation << ", "
	<< "photon isolation coefficients (loose): " << cuts.photonIsolationLoose;
    return out;
  }
};

enum class PFTypesForEA{chargedHadron=0, neutralHadron, photon};
std::map<PFTypesForEA, std::string> PFTypesForEANames = {
  {PFTypesForEA::chargedHadron, "chargedHadron"},
  {PFTypesForEA::neutralHadron, "neutralHadron"},
  {PFTypesForEA::photon, "photon"}
};

struct EAValuesStruct{
  float regionUpperBound, chargedHadronsEA, neutralHadronsEA, photonsEA;

  EAValuesStruct () : regionUpperBound(0.),
    chargedHadronsEA(0.),
    neutralHadronsEA(0.),
    photonsEA(0.) {}
  
  EAValuesStruct (float regionUpperBound_, float chargedHadronsEA_, float neutralHadronsEA_, float photonsEA_) : regionUpperBound(regionUpperBound_),
    chargedHadronsEA(chargedHadronsEA_),
    neutralHadronsEA(neutralHadronsEA_),
    photonsEA(photonsEA_) {}

  float getEffectiveArea(const PFTypesForEA& PFType) const{
    float effectiveArea = 0.0;
    switch(PFType) {
    case (PFTypesForEA::chargedHadron) :
      effectiveArea = chargedHadronsEA;
      break;
    case (PFTypesForEA::neutralHadron) :
      effectiveArea = neutralHadronsEA;
      break;
    case (PFTypesForEA::photon) :
      effectiveArea = photonsEA;
      break;
    default :
      std::cout << "ERROR: Unknown PF type for EA!"<< std::endl;
      std::exit(EXIT_FAILURE);
    }
    return effectiveArea;
  }
  
  friend std::ostream& operator<< (std::ostream& out, const EAValuesStruct& EAValues) {
    out << "region upper eta bound --> " << EAValues.regionUpperBound << ", effective areas: "
        << "charged hadrons --> " << EAValues.chargedHadronsEA << ", "
        << "neutral hadrons --> " << EAValues.neutralHadronsEA << ", "
        << "photons --> " << EAValues.photonsEA;
    return out;
  }
};

struct eventWeightsStruct{
  float nominal, down, up;
  eventWeightsStruct () : nominal(1.0f), down(1.0f), up(1.0f) {} // empty constructor
  eventWeightsStruct (float nominal_, float down_, float up_) : nominal(nominal_), down(down_), up(up_) {}
};

#endif
