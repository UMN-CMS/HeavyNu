#include "ZShape/Base/interface/ZShapeEvent.h"
#include "ZShape/Base/interface/ZShapeZDef.h"

void ZShapeZDef::addCriterion(ZShapeZDef::CriterionType ct, const std::string& crit) {
  m_requirements[ct].push_back(crit);
  if (ct==crit_Z) {
    char pref[10];
    float a,b;
    ZReqStruct zrs;
    zrs.reqType=zr_Undef;

    int cvt=sscanf(crit.c_str(),"%[^(](%f,%f)",pref,&a,&b);
    zrs.lower=a;
    zrs.upper=b;

    if (cvt==3) {
      if (!strcasecmp(pref,"m")) {
	zrs.reqType=zr_Mass;
      }
      if (!strcasecmp(pref,"Y")) {
	zrs.reqType=zr_Rapidity;
      }
      if (!strcasecmp(pref,"G")) {
	zrs.reqType=zr_NGsf20;
      }
    }
    m_zCuts.push_back(zrs);
  }
}

int ZShapeZDef::criteriaCount(ZShapeZDef::CriterionType ct) const {
  return (int)(m_requirements[ct].size());
}

static const std::string BAD("Overflow Error");
const std::string& ZShapeZDef::criteria(ZShapeZDef::CriterionType ct, int n) const {
  return (n<0 || n>=criteriaCount(ct))?(BAD):(m_requirements[ct][n]);
}

bool ZShapeZDef::pass(const ZShapeEvent& evt, bool* swapPairing) const {
  return pass(evt,m_requirements[0].size(),m_requirements[1].size(),m_requirements[2].size(),swapPairing);
}

bool ZShapeZDef::pass(const ZShapeEvent& evt, int ne1, int ne2, int nz,bool* swapPairing) const {
  // try straight-up assignment first
  bool p=pass(evt.elec(0),evt.elec(1),ne1,ne2);
  if (!p) { // try swapping electrons!
    p=pass(evt.elec(1),evt.elec(0),ne1,ne2);
    if (swapPairing!=0) *swapPairing=true;
  } else if (swapPairing!=0) *swapPairing=false;

  if (p) {
    p=passZ(evt, nz);
  }

  return p;
}

bool ZShapeZDef::pass(const ZShapeElectron& e1, const ZShapeElectron& e2, int ne1, int ne2) const {
  bool ok=true;
  // concrete assignment
  for (unsigned int i=0; ok && int(i)<ne1 && i<m_requirements[crit_E1].size(); ++i) {
    if (!e1.passed(m_requirements[crit_E1][i])) ok=false;
  }
  for (unsigned int i=0; ok && int(i)<ne2 && i<m_requirements[crit_E2].size(); ++i) {
    if (!e2.passed(m_requirements[crit_E2][i])) ok=false;
  }
  
  return ok;
}

bool ZShapeZDef::passZ(const ZShapeEvent& evt, int nz) const {
  bool ok=true;

  for (unsigned int i=0; ok && int(i)<nz && i<m_requirements[crit_Z].size(); ++i) {
    bool rule=true;
    if (m_zCuts[i].reqType==zr_Mass) {
      rule=evt.m()>=m_zCuts[i].lower && evt.m()<m_zCuts[i].upper;
    } else if (m_zCuts[i].reqType==zr_Rapidity) {
      rule=evt.Y()>=m_zCuts[i].lower && evt.Y()<m_zCuts[i].upper;
    } else if (m_zCuts[i].reqType==zr_NGsf20) {
      rule=evt.n_gsf20 >=m_zCuts[i].lower && evt.n_gsf20 <=m_zCuts[i].upper;
    }
    ok=ok && rule;
  }
  return ok;
}
