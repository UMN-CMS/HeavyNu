/** -- C++ -- **/
#ifndef HeavyNuID_h_included
#define HeavyNuID_h_included

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include <string>
#include <map>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

class HeavyNuID {
 public:

  explicit HeavyNuID(const edm::ParameterSet & iConfig);
  double   weightForMC(double pt,int signOfError2apply=0);
  double   weightForMCbyEta(double eta,int signOfError2apply);
  double   weightElectronsForMC(double eta,int signOfError2apply=0);
  int      idEra() { return idEra_ ; } 
  void endJob();

 private:
  int idEra_ ; 

};

#endif
