/** -- C++ -- **/
#ifndef HeavyNuTrigger_h_included
#define HeavyNuTrigger_h_included

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include <string>
#include <map>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TRandom.h"

class HeavyNuTrigger {
 public:

  struct trigHistos_t {
    TH2D *trigMatchPtCorrel;
    TH1D *trigMatchDR2;
    TH2D *trigMatchDRDPt;
    TH2D *trigMatchDetaPhi;
    TH1D *trigUnmatchedPt;
    TH1D *trigAllCandMuPt;
    TH2D *trigUnmatchedEtaPhi;
    TH2D *trigAllCandMuEtaPhi;
  };

  explicit HeavyNuTrigger(const edm::ParameterSet & iConfig);
  void book(const TFileDirectory& tdir, trigHistos_t *thist);

  inline bool matchingEnabled() { return matchingEnabled_; }

  bool isTriggerMatched(const pat::MuonRef & muon,
			const edm::Event   & iEvent,
			trigHistos_t *thist = NULL);

  bool simulateForMC(double pt,double eta,int signOfError2apply=0);
  bool simulateForMC(double pt,int signOfError2apply=0) { return simulateForMC(pt,0.,signOfError2apply); }  

  int trigEra() { return year_ ; } 
  void endJob();

 private:

  bool          matchingEnabled_;
  TRandom      *triggerRandom_;
  edm::InputTag trigEventTag_;
  std::string   muonMatch_;
  int           johnnyApple_;
  int           year_;
};

#endif
