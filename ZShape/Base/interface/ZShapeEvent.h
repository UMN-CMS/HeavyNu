#ifndef ZSHAPE_BASE_ZSHAPEEVENT_H
#define ZSHAPE_BASE_ZSHAPEEVENT_H 1

#include "DataFormats/Math/interface/Vector3D.h"
#include "ZShape/Base/interface/ZShapeElectron.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <map>
#include <string>
#include <ostream>

/** \class ZShapeEvent
  *  
  * $Date: 2009/05/20 22:39:05 $
  * $Revision: 1.5 $
  * \author J. Mans - Minnesota
  */
class ZShapeEvent {
public:
  ZShapeEvent();
  // primary vertex
  ::math::XYZPoint vtx_;
  int n_elec;
  int n_TLelec;
  int n_gsf20;

  double m() const { return mass_; }
  double Y() const { return rap_; }

  ZShapeElectron& elec(int i) { return ((i==0)?(*e1_):(*e2_)); }
  const ZShapeElectron& elec(int i) const { return ((i==0)?(*e1_):(*e2_)); }
  reco::GenParticle& elecTreeLevel(int i) { return ((i==0)?(*eTL1_):(*eTL2_)); }
  const reco::GenParticle& elecTreeLevel(int i) const { return ((i==0)?(*eTL1_):(*eTL2_)); }

  virtual ~ZShapeEvent();
  virtual void clear();
  virtual void dump(std::ostream& s) const;
  virtual void dump() const;
  void afterLoad();
private:
  virtual ZShapeElectron* allocateElectron();
  virtual reco::GenParticle* allocateTreeLevelElectron();
  ZShapeElectron* e1_;
  ZShapeElectron* e2_;
  reco::GenParticle * eTL1_;
  reco::GenParticle * eTL2_;
  double mass_, rap_;
};

#endif
