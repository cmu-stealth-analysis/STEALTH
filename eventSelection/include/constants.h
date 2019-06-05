#ifndef H_CONSTANTS
#define H_CONSTANTS

namespace constants{ // for readability
  const float TWOPI = static_cast<float>(2.0*TMath::Pi());
}

struct PIDsStruct {
  const int photon = 22;
  const int gluino = 1000021;
  const int neutralino = 1000022;
  friend std::ostream& operator<< (std::ostream& out, const PIDsStruct& PIDs) {
    out << "photon --> " << PIDs.photon << ", "
        << "gluino --> " << PIDs.gluino << ", "
        << "neutralino --> " << PIDs.neutralino;
    return out;
  }
};

#endif
