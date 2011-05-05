#include "ZShape/Base/interface/ZShapeStandardCuts.h"

void ZShapeStandardCuts::acceptanceCuts(ZShapeElectron& elec) const {
  elec.cutResult("ACC(EB)",theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_EB));
  elec.cutResult("ACC(EE)",theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_EE));
  elec.cutResult("ACC(EE+TRK)",theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_EE_tracker));
  elec.cutResult("ACC(ECAL)",theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_EE));
  elec.cutResult("ACC(ECAL+TRK)",theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_ECAL) && (!theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_EE_notracker)));
  elec.cutResult("ACC(ECAL-TRK)",theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_ECAL) && (theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_EE_notracker)));
  elec.cutResult("ACC(HF)",theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_HF));
  elec.cutResult("ACC(ANY)",theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_ANY));
  elec.cutResult("ACC(ANY-EEnotrk)",theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_ANY) && (!theEtaSelector_.isInAcceptance(elec,EtaAcceptance::zone_EE_notracker)));
}

void ZShapeStandardCuts::ptCuts(ZShapeElectron& elec) const {
  double pt=elec.p4_.Pt();
  elec.cutResult("PT5",pt>=5);
  elec.cutResult("PT10",pt>=10);
  elec.cutResult("PT15",pt>=15);
  elec.cutResult("PT20",pt>=20);
  elec.cutResult("PT25",pt>=25);
}

void ZShapeStandardCuts::dummyCuts(ZShapeElectron& elec) const {
  elec.cutResult("DummyPass",true);
  elec.cutResult("DummyFail",false);
}
