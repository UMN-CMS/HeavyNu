// -*- C++ -*-
//
// Package:    EleWithMet
// Class:      EleWithMet
// 
/**\class EleWithMet EleWithMet.cc ZShape/HFZeeVBTF/src/EleWithMet.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  G. Franzoni
//         Created:  Wed Aug 04 00:21:26 CDT 2010
// $Id: EleWithMet.cc,v 1.5 2010/08/13 15:09:14 franzoni Exp $
//
//


// system include files
#include <memory>

#include <algorithm>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/FileBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "DataFormats/EgammaReco/interface/HFEMClusterShape.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShapeAssociation.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShape.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "TH1.h"
#include "TH2.h"
#include <TMath.h>
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "ZShape/HFZeeVBTF/interface/SelectElectron.h"

//////////////////////////////////////////////////////////////////
// generically maximum
template <class T> const T& max ( const T& a, const T& b ) {
  return (b<a)?a:b;     // or: return comp(b,a)?a:b; for the comp version
}

class EleWithMet : public edm::EDAnalyzer {
public:
  explicit EleWithMet(const edm::ParameterSet&);
  ~EleWithMet();


private:
  virtual void respondToOpenInputFile(edm::FileBlock const& fb) {
    currentFile_=fb.fileName();
  }

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  std::string ecalID_, currentFile_;
  bool dolog_;
  std::string myName_;
  std::vector<int> acceptedElectronIDs_;
  selectElectron w;
  float ETCut_;
  float METCut_;
  edm::InputTag metCollectionTag_;


  // ----------member data ---------------------------

  // gf set of histograms per each == W definition ==

  struct HistPerDef {
    //book histogram set w/ common suffix inside the provided TFileDirectory
    void book(TFileDirectory td,const std::string&);
    // fill all histos of the set with the two electron candidates
    void fill(pat::ElectronCollection::const_iterator ecalE);
    //    void fill(pat::ElectronCollection::const_iterator ecalE,  const reco::RecoEcalCandidate& hfE, const reco::HFEMClusterShape& hfshape, const edm::ParameterSet& myPs, const float hf_2d_cut);
    
    TH1* eleEta, *elePhi, *eleEt;
    
  };

  // gf set of histo for all Z definitios in a stack
  struct HistStruct {
    TH1 *nelec, *met, *metEta, *metPhi;
    HistPerDef ele90metLoose, ele80metLoose, ele70metLoose, ele60metLoose,
      ele90metTight, ele80metTight, ele70metTight,    ele60metTight;
  } hists;

};


void EleWithMet::HistPerDef::book(TFileDirectory td, const std::string& post) {
  std::string title;

  title=std::string("eta_{el} ")+post+std::string(";#eta_{el}");
  eleEta=td.make<TH1D>("eta_{el}",title.c_str(),30,-3,3);
  title=std::string("phi_{el} ")+post+std::string(";#phi_{el}");
  elePhi=td.make<TH1D>("phi_{el}",title.c_str(),30,-1*TMath::Pi(),TMath::Pi());
  title=std::string("et_{el} ")+post+std::string(";Et_{el}");
  eleEt=td.make<TH1D>("et_{el}",title.c_str(),50,0,100);

}

void EleWithMet::HistPerDef::fill(pat::ElectronCollection::const_iterator ecalE//,  
				  //const reco::RecoEcalCandidate& hfE, const reco::HFEMClusterShape& hfshape, 
				  //				  const edm::ParameterSet& myPs, const float hf_2d_cut
				  ) {
  
  //pat::ElectronCollection::const_iterator ecalE,  
  eleEta -> Fill(ecalE->eta());
  elePhi -> Fill(ecalE->phi());
  eleEt  -> Fill(ecalE->pt());

}// end of fill()

// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
EleWithMet::EleWithMet(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  ecalID_=iConfig.getParameter<std::string>("ECALid");
  dolog_=iConfig.getParameter<bool>("DoLog");
  acceptedElectronIDs_ = iConfig.getParameter< std::vector<int> >("acceptedElectronIDs");
  myName_=iConfig.getParameter<std::string>("myName");
  myName_+=std::string("    ");

  metCollectionTag_ = iConfig.getParameter<edm::InputTag>("metCollectionTag");

  ETCut_  = iConfig.getParameter<double>("ETCut");
  METCut_ = iConfig.getParameter<double>("METCut");

  std::vector<int>::const_iterator it;
  for(it=acceptedElectronIDs_.begin(); it!=acceptedElectronIDs_.end(); it++)
    {    w.add(*it);//      std::cout << "adding: " << (*it) << std::endl;
    }

  edm::Service<TFileService> fs;

  // these are general histograms filled regardless of the eleId/Met selection
  hists.nelec=fs ->make<TH1D>("nelec","N_Elec; num electron",10,-0.5,9.5);
  hists.met=fs   ->make<TH1D>("met","met;mE_{T} [GeV]",100,0,100);
  hists.metEta=fs->make<TH1D>("met eta","metEta;#eta_{mE_{T}}",3,-5,5);
  hists.metPhi=fs->make<TH1D>("met phi","metPhi;#phi_{mE_{T}}",36,-1*TMath::Pi(),TMath::Pi());

  hists.ele90metLoose.book(fs->mkdir("ele90metLoose"),"(ele90,metLoose)");
  hists.ele80metLoose.book(fs->mkdir("ele80metLoose"),"(ele80,metLoose)");
  hists.ele70metLoose.book(fs->mkdir("ele70metLoose"),"(ele70,metLoose)");  
  hists.ele60metLoose.book(fs->mkdir("ele60metLoose"),"(ele60,metLoose)");  

  hists.ele90metTight.book(fs->mkdir("ele90metTight"),"(ele90,metTight)");  
  hists.ele80metTight.book(fs->mkdir("ele80metTight"),"(ele80,metTight)");  
  hists.ele70metTight.book(fs->mkdir("ele70metTight"),"(ele70,metTight)");  
  hists.ele60metTight.book(fs->mkdir("ele60metTight"),"(ele60,metTight)");  

}
  
EleWithMet::~EleWithMet()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
EleWithMet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::Handle<reco::RecoEcalCandidateCollection> HFElectrons;
  edm::Handle<reco::SuperClusterCollection> SuperClusters;
  edm::Handle<reco::HFEMClusterShapeAssociationCollection> ClusterAssociation;
  
  iEvent.getByLabel("hfEMClusters",SuperClusters);
  iEvent.getByLabel("hfEMClusters",ClusterAssociation);
  iEvent.getByLabel("hfRecoEcalCandidate",HFElectrons);
  
  edm::Handle<pat::ElectronCollection> patElectron;
  iEvent.getByLabel("patElectrons", patElectron);
  if ( ! patElectron.isValid()) {
    std::cout << "No electrons found in this event with tag " 
	      << "patElectrons" << std::endl;
    return ; // RETURN if no elecs in the event
  }

  const pat::ElectronCollection& pElecs = *(patElectron.product());
  pat::ElectronCollection::const_iterator i, ecalE=pElecs.end();
  
  for (i=pElecs.begin(); i!=pElecs.end(); i++) {
    reco::SuperClusterRef scr=i->superCluster();
    
    if (scr.isNull()) {
      std::cout << "I'm NULL!!!\n";
      continue;
    }
    
    // gf: ECAL eta acceptance cut on supercluster. What about phi acceptance and central crack?
    // explanation of electronID value: https://twiki.cern.ch/twiki/bin/view/CMS/SimpleCutBasedEleID
    double eta_det=scr.get()->eta();
    
    // ECAL acceptance and minEnergy cut
    if (
	ecalE==pElecs.end()   // 3: ECAL acceptance cut
	&& (fabs(eta_det)<1.4442 || fabs(eta_det)>1.560)
	&&  i->pt() > ETCut_
	) {
      ecalE=i;
    }// end if
    
  }// end loop over electrons
  

  // reqiure presence of at least one GSF electron within acceptance and min energy
  if (ecalE==pElecs.end() ) return;

  
  hists.nelec->Fill(pElecs.size());
  
  // Met used to select electrons from W decay
  edm::Handle<pat::METCollection> patMET;
  iEvent.getByLabel(metCollectionTag_,  patMET);
  const pat::METCollection *pMet = patMET.product();
  const pat::METCollection::const_iterator met = pMet->begin();
  const pat::MET theMET = *met;

  hists.metEta->Fill( met->eta() ); //(eta is pointless for met - jsut a consistency check that it's always 0)
  hists.met   ->Fill( met->et() );  // add a plot which has Met phi above a threshold!
  hists.metPhi->Fill( met->phi() );
  
  if( met->et()  > METCut_){
    if(w.doesElePass( ecalE->electronID("simpleEleId90relIso") ))  hists.ele90metLoose.fill(ecalE);
    if (w.doesElePass( ecalE->electronID("simpleEleId80relIso") )) hists.ele80metLoose.fill(ecalE);
    if (w.doesElePass( ecalE->electronID("simpleEleId70relIso") )) hists.ele70metLoose.fill(ecalE);
    if (w.doesElePass( ecalE->electronID("simpleEleId60relIso") )) hists.ele60metLoose.fill(ecalE);
  }
  if( met->et()  > (METCut_+10) ){
    if(w.doesElePass( ecalE->electronID("simpleEleId90relIso") ))  hists.ele90metTight.fill(ecalE);
    if (w.doesElePass( ecalE->electronID("simpleEleId80relIso") )) hists.ele80metTight.fill(ecalE);
    if (w.doesElePass( ecalE->electronID("simpleEleId70relIso") )) hists.ele70metTight.fill(ecalE);
    if (w.doesElePass( ecalE->electronID("simpleEleId60relIso") )) hists.ele60metTight.fill(ecalE);
  }

}


// reminder re electronID(""):
//0: fails
//1: passes electron ID only
//2: passes electron Isolation only
//3: passes electron ID and Isolation only   <<< 
//4: passes conversion rejection
//5: passes conversion rejection and ID
//6: passes conversion rejection and Isolation
//7: passes the whole selection
//https://twiki.cern.ch/twiki/bin/view/CMS/SimpleCutBasedEleID 

// ------------ method called once each job just before starting event loop  ------------
void 
EleWithMet::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EleWithMet::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(EleWithMet);
