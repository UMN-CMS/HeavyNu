#ifndef ZSHAPE_BASE_ZSHAPEZDEF_H
#define ZSHAPE_BASE_ZSHAPEZDEF_H 1

#include <vector>
#include <string>
class ZShapeEvent;

/** \class ZShapeZDef
  *
  * Basically, a Z-definition is a list of cuts (by name) which must pass.
  * These cuts apply to the first electron, the second electron, and then
  * the di-electron. 
  * 
  * $Date: 2009/05/20 22:39:05 $
  * $Revision: 1.3 $
  * \author J. Mans - Minnesota
  */
class ZShapeZDef {
public:
  enum CriterionType { crit_Z=0, crit_E1=1, crit_E2=2 };
  /// add a criterion
  void addCriterion(CriterionType ctype, const std::string& cit);
  /// get the number of criteria
  int criteriaCount(CriterionType ctype) const;
  /// get the n-th criterion
  const std::string& criteria(CriterionType ctype, int n) const;
  /// does this event pass all criteria of this Z Definition?
  bool pass(const ZShapeEvent& evt, bool* swapPairing=0) const;
  /// does this event pass all criteria of this Z Definition up-to-and-including the n-th for each item?
  bool pass(const ZShapeEvent& evt, int ne1, int ne2=1000, int nz=1000, bool* swapPairing=0) const;
private:
  bool pass(const ZShapeElectron& e1, const ZShapeElectron& e2, int ne1, int ne2) const;
  bool passZ(const ZShapeEvent& evt, int nz) const;
  std::vector<std::string> m_requirements[3];
  
  enum ZReqType { zr_Undef, zr_Mass, zr_Rapidity, zr_NGsf20 };
  struct ZReqStruct {
    ZReqType reqType;
    double lower, upper;
  };
  std::vector<ZReqStruct> m_zCuts;

};

#endif
