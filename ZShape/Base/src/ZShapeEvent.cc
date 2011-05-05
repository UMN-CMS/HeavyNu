#include "ZShape/Base/interface/ZShapeEvent.h"
#include <DataFormats/Math/interface/LorentzVector.h>
#include <iostream>

ZShapeElectron* ZShapeEvent::allocateElectron() {
  return new ZShapeElectron(); // vanilla
}

reco::GenParticle* ZShapeEvent::allocateTreeLevelElectron() {
  return new reco::GenParticle(); // vanilla(?)
}

void ZShapeEvent::dump() const {
  dump(std::cout);
}
void ZShapeEvent::dump(std::ostream& s) const {
  s << "Vertex: " << vtx_.x() << ", " << vtx_.y() << ", " << vtx_.z() << std::endl;
  for (int i=0; i<2; i++) {
    s << "Electron " << i+1;
    elec(i).dump(s);
  }
}

void ZShapeEvent::clear() {
  if (e1_==0) e1_=allocateElectron();
  if (e2_==0) e2_=allocateElectron();
  elec(0).clear();
  elec(1).clear();
  n_elec=0;
  n_gsf20=0;

  if(eTL1_==0) eTL1_=allocateTreeLevelElectron();
  if(eTL2_==0) eTL2_=allocateTreeLevelElectron();
  n_TLelec=0;

  // initializing
  eTL1_->setP4(reco::Particle::LorentzVector(0,0,0,0));
  eTL1_->setCharge(0);
  eTL1_->setVertex(math::XYZPoint(0,0,0));
  eTL1_->setPdgId(0);
  eTL1_->setStatus(0);
  eTL2_->setP4(reco::Particle::LorentzVector(0,0,0,0));
  eTL2_->setCharge(0);
  eTL2_->setVertex(math::XYZPoint(0,0,0));
  eTL2_->setPdgId(0);
  eTL2_->setStatus(0);

}

ZShapeEvent::ZShapeEvent() : e1_(0),e2_(0),eTL1_(0),eTL2_(0) {
}

ZShapeEvent::~ZShapeEvent() {
  if (e1_!=0) delete e1_;
  if (e2_!=0) delete e2_;
  if (eTL1_!=0) delete eTL1_;
  if (eTL2_!=0) delete eTL2_;
}

void ZShapeEvent::afterLoad() {
  math::PtEtaPhiMLorentzVector ZP4 = elec(0).p4_ + elec(1).p4_;
  mass_=ZP4.M();
  rap_=ZP4.Rapidity();
}

