#ifndef H_MISCUTILS
#define H_MISCUTILS

#include "Rtypes.h"

namespace miscUtils{
  bool passesBitMask(const UShort_t& bitCollection, const UShort_t& bitMask) {
    return ((bitCollection&bitMask) == bitMask);
  }

  bool checkHLTBit(const ULong64_t& inputHLTBits, const int& indexOfBitToCheck) {
    return (((inputHLTBits>>indexOfBitToCheck)&1) == 1);
  }

  template<typename criterion>
  int getNFalseBits(std::map<criterion, bool>& bits) {
    int nFalseBits = 0;
    for (auto&& bitsElement: bits) {
      if (!(bitsElement.second)) ++nFalseBits;
    }
    return nFalseBits;
  }

  template<typename criterion>
  criterion getFirstFalseCriterion(std::map<criterion, bool>& bits) {
    for (auto&& bitsElement: bits) {
      if (!(bitsElement.second)) return (bitsElement.first);
    }
    // Control shouldn't reach here
    std::cout << "ERROR: getFirstFalseCriterion called with a collection of bits of which none is false." << std::endl;
    std::exit(EXIT_FAILURE);

    // Formality just to get code to compile
    return (bits.begin())->first;
  }
}

#endif
