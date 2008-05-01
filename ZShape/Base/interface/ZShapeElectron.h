#ifndef ZShapeElectron_h_included
#define ZShapeElectron_h_included 1

#include "DataFormats/Math/interface/LorentzVector.h"
#include <ostream>
#include <string>
#include <map>

class ZShapeElectron {
public:
  bool passed(const std::string& cutName) const; // false if not present!
  void cutResult(const std::string& cutName, bool didpass);
  // four-vector of electron (standardized)
  math::PtEtaPhiMLorentzVector p4_;
  void dump(std::ostream& s) const;
  void clear();
private:
  // pass/fail for each possible cut
  std::map<std::string, bool> cutMap_;
  
};


#endif // ZShapeElectron_h_included