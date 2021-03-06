// -*- C++ -*-
//
// Package:    MuJetBackground
// Class:      MuJetBackground
// 
/**\class MuJetBackground MuJetBackground.cc HeavyNu/AnalyzerModules/src/MuJetBackground.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jeremy M Mans
//         Created:  Mon May 31 07:00:26 CDT 2010
// $Id: MuJetBackground.cc,v 1.29 2013/02/21 22:17:39 bdahmes Exp $
//
//

// system include files
#include <memory>
#include <iostream>
#include <algorithm>
#include <vector>

// According to
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
// this must be included before Frameworkfwd.h
//
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/FileBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile2D.h"
#include "TVector3.h"

#include "HeavyNu/AnalysisModules/src/HeavyNuEvent.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuID.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuTrigger.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuCommon.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//////////////////////////////////////////////////////////////////
// generic maximum/minimum
template <class T> const T& max ( const T& a, const T& b ) {
  return (b<a)?a:b;
}
template <class T> const T& min ( const T& a, const T& b ) {
  return (b<a)?b:a;
}

template <class T>
inline std::string int2str(T i) {
  std::ostringstream ss;
  ss << i;
  return ss.str();
}

//============================================================

class compare {
public:
  template <class T> bool operator() (const T& a, const T& b) { return a.pt() > b.pt() ; } 
};

//============================================================

class MuJetBackground : public edm::EDFilter {
public:
  explicit MuJetBackground(const edm::ParameterSet&);
  ~MuJetBackground();


private:
  virtual void respondToOpenInputFile(edm::FileBlock const& fb) {
    currentFile_=fb.fileName();
  }
  
  virtual void beginJob          ();
  virtual bool filter            ( edm::Event&, const edm::EventSetup& );
  virtual void endJob            ();

  void initializeHNE(HeavyNuEvent& hne, edm::Handle< std::vector<PileupSummaryInfo> >& pPU, 
		     edm::Handle<reco::VertexCollection>& pvHandle, bool isMC, bool isPF, double weight) ;

//   virtual bool selectJets        ( std::vector< std::pair<pat::Jet,float> >& jets,
// 				   HeavyNuEvent& hne, int minNjets=2 );
  virtual void selectMuonsInJets ( std::vector<pat::Muon>& muons,
				   std::vector< std::pair<pat::Jet,float> >& jets,
				   HeavyNuEvent& hne );
  virtual void selectMuonPairs   ( std::vector<pat::Muon>& muons,
                                   HeavyNuEvent& hne );
  std::vector<pat::Muon> getNonIsolatedMuons ( std::vector<pat::Muon>& muons ) ; 
  bool findQCDmuon               ( const std::vector<pat::Muon>& muons, 
				   edm::Handle<pat::MuonCollection>& pMuons,
				   HeavyNuEvent& hne );
  bool secondQualityMuon         (const pat::Muon& mu1, const pat::Muon& mu2);
  virtual bool findQCDjet        ( const std::vector< std::pair<pat::Jet,float> >& jets,
				   edm::Handle<pat::MuonCollection>& pMuons,
				   HeavyNuEvent& hne );
  virtual bool isDijetCandidate  ( HeavyNuEvent& hne,pat::MET& theMET ); 
  virtual TH1 *bookRunHisto      ( uint32_t runNumber );
  
  edm::InputTag muonTag_;
  edm::InputTag jetTag_, metTag_;
  edm::InputTag jptTag_;

  std::string currentFile_;
  bool dolog_;
  bool firstEvent_;
  // int theMETtype ; 
  std::vector<double> rwLowPtbin, rwHighPtbin ; 
  std::vector<double> rw2011A, rw2011B, rwQCD ; 

  // double trigEtaLimit_ ; 
  bool calcSurvival_ ; 
  // bool doClosure_, doQuadJet_ ; 

  JetCorrectionUncertainty *jecuObj_;

  bool   applyMuIDCorrections_ ;
    
  // HeavyNu_NNIF *nnif_;
  HeavyNuTrigger *trig_;
  HeavyNuID *muid_ ; 

  int    pileupEra_;
  bool isPFJets_; 
  edm::LumiReWeighting MCweightByVertex_;
  std::map<uint32_t,TH1 *> m_runHistos_;

  // ----------member data ---------------------------

  struct HistPerDef {
    //book histogram set w/ common suffix inside the provided TFileDirectory
    void book(TFileDirectory *, const std::string&, const std::vector<hNuMassHypothesis>&) ;
    // fill all histos of the set with the mu-jet candidate
    void fill(pat::Muon& theMuon, pat::Jet& theJet, pat::MET& theMET, bool isMC, double trkIso) ;
    // void fill(reco::SuperCluster theSC) ; 
    // void fill(pat::MuonCollection muons, pat::JetCollection jets,bool isMC) ;
    // fill all histos of the set with the two muon candidates
    void fill(const HeavyNuEvent& hne, double w1=1.0, double w2=1.0) ;

    TH1 *ptMu1, *ptMu2, *ptJet1, *ptJet2;
    TH1 *etaMu1, *etaMu2, *etaJet1, *etaJet2 ;
    TH1 *phiMu1, *phiMu2, *phiJet1, *phiJet2 ;
    TH1 *ptMET, *phiMET ; 
    TH1 *ptMu1trkIso, *etaMu1trkIso, *phiMu1trkIso ; 

    TH3 *ptMu_etaMu_trkIso ; 
      
    TH1 *ecalIso_muPt100, *hcalIso_muPt100, *hoeIso_muPt100 ; 
    TH1 *ecalIso_muPt200, *hcalIso_muPt200, *hoeIso_muPt200 ;
      
    TH1 *mWR, *mNuR1, *mNuR2, *mMuMu, *mMuMuZoom, *mJJ ; 
    TH2 *mNuR2D, *mNuR2D_raw ; 

    TH1 *mWR_raw, *mNuR1_raw, *mNuR2_raw, *mMuMu_raw, *mMuMuZoom_raw, *mJJ_raw ; 

    TFileDirectory *mydir;
    // TFileDirectory *nndir;

    // HeavyNuTrigger::trigHistos_t trigHistos;
  };

  bool init_;

  // gf set of histo for all Z definitions in a stack
  struct HistStruct {
    TH1 *nmu, *njet ;

    TFileDirectory *rundir;

    HistPerDef NoCuts;
    HistPerDef dPhiCuts;
    HistPerDef dPhiJetCuts;
    HistPerDef dPhiCoreCuts;
    HistPerDef dPhi10pctCuts;
    HistPerDef dPhiTrigCuts;
    HistPerDef dPhiJetTrigCuts;
    HistPerDef dPhiCoreTrigCuts;
    HistPerDef dPhi10pctTrigCuts;
    HistPerDef dPhiTightCuts;
    HistPerDef dPhiJetTightCuts;
    HistPerDef dPhiCoreTightCuts;
    HistPerDef dPhi10pctTightCuts;
    HistPerDef dPhiTightTrigCuts;
    HistPerDef dPhiJetTightTrigCuts;
    HistPerDef dPhiCoreTightTrigCuts;
    HistPerDef dPhi10pctTightTrigCuts;

//     HistPerDef LLJJpTCuts;
//     HistPerDef TrigMatches;
//     HistPerDef VertexCuts;
//     HistPerDef Mu1HighPtCut;
//     HistPerDef diLmassCut;

//     // Closure tests using two muons --> run automatically with quadjet
//     HistPerDef mmJClosure; 
//     HistPerDef mMuJClosure;
//     HistPerDef mmJClosureMETlt20; 
//     HistPerDef mMuJClosureMETlt20;

//     HistPerDef ssmmjjClosure; 
//     HistPerDef ssMuMujjClosure; 

//     HistPerDef MuMuJClosure;

//     // Closure test using one muon
//     HistPerDef LJJJClosure; 
//     HistPerDef L2JClosure;
  } hists;

  struct CutsStruct {
    double minimum_mu1_pt;
    double minimum_mu2_pt;
    double minimum_jet_pt;
    double maximum_mu_abseta;
    double maximum_jet_abseta;
    // double minimum_mumu_mass;
    double minimum_muon_jet_dR;
    double muon_trackiso_limit;
    // double maxVertexZsep;
    // double maxJetVZsepCM;
    double minimum_dijet_dPhi; 
    double minimum_dijet_pt; 
    // double minimum_SCEt; 
    double minimum_extraJet_dR; 
  } cuts;
  
};

//======================================================================

void MuJetBackground::initializeHNE(HeavyNuEvent& hne, 
				    edm::Handle< std::vector<PileupSummaryInfo> >& pPU, 
				    edm::Handle<reco::VertexCollection>& pvHandle, 
				    bool isMC, bool isPF, double weight) {
  hne.isMC   = isMC ; 
  hne.pfJets = isPF ; 

  if (hne.isMC) { 
    std::pair<float,double> pileup = hnu::pileupReweighting(pPU,MCweightByVertex_) ; 
    hne.n_pue     = int(pileup.first) ; 
    hne.eventWgt *= pileup.second ; 
  }
  hne.n_primary_vertex = hnu::numberOfPrimaryVertices(pvHandle) ;
}


double qcdScaleFactor(double pt, 
		   std::vector<double> ptLow, 
		   std::vector<double> ptHigh, 
		   std::vector<double> corr) {
  double correction = 1.0 ; 
  for (unsigned int i=0; i<ptLow.size(); i++) 
    if ( (pt >= ptLow.at(i)) && (pt < ptHigh.at(i)) ) return corr.at(i) ; 
  return correction ; 
}

//======================================================================

void
MuJetBackground::HistPerDef::book(TFileDirectory *td,
				  const std::string& post,
				  const std::vector<hNuMassHypothesis>& v_masspts )
{
  std::string t; // histogram title string;
  
  TH1::SetDefaultSumw2();

  mydir = td;

  // ----------  Muon histograms  ----------

  t="p_{T}(#mu_{1}) "+post;   ptMu1=td->make<TH1D>("ptMu1",t.c_str(),200,0.,1000.);
  t="p_{T}(#mu_{2}) "+post;   ptMu2=td->make<TH1D>("ptMu2",t.c_str(),200,0.,1000.);
  t="#eta(#mu_{1}) " +post;  etaMu1=td->make<TH1D>("etaMu1",t.c_str(),40,-2.5,2.5);
  t="#eta(#mu_{2}) " +post;  etaMu2=td->make<TH1D>("etaMu2",t.c_str(),40,-2.5,2.5);
  t="#phi(#mu_{1}) " +post;  phiMu1=td->make<TH1D>("phiMu1",t.c_str(),30,-3.14159,3.14159);
  t="#phi(#mu_{2}) " +post;  phiMu2=td->make<TH1D>("phiMu2",t.c_str(),30,-3.14159,3.14159);

  t="p_{T}(#mu_{1}) "+post;   ptMu1trkIso=td->make<TH1D>("ptMu1trkIso",t.c_str(),200,0.,1000.);
  t="#eta(#mu_{1}) " +post;  etaMu1trkIso=td->make<TH1D>("etaMu1trkIso",t.c_str(),40,-2.5,2.5);
  t="#phi(#mu_{1}) " +post;  phiMu1trkIso=td->make<TH1D>("phiMu1trkIso",t.c_str(),30,-3.14159,3.14159);

  t="p_{T}(#mu) vs. #eta(#mu) vs. relative tracker isolation "+post;
  ptMu_etaMu_trkIso = td->make<TH3D>("ptMu_etaMu_trkIso",t.c_str(),100,0.,1000.,50,-2.5,2.5,100.,0.,10.);
  
  // ----------  Jet histograms ----------

  t="p_{T}(j_{1}) "            +post;     ptJet1=td->make<TH1D>("ptJet1",  t.c_str(),50,0.,500.);
  t="p_{T}(j_{2}) "            +post;     ptJet2=td->make<TH1D>("ptJet2",  t.c_str(),50,0.,500.);
  t= "#eta(j_{1}) "            +post;    etaJet1=td->make<TH1D>("etaJet1", t.c_str(),40,-5,5);
  t= "#eta(j_{2}) "            +post;    etaJet2=td->make<TH1D>("etaJet2", t.c_str(),40,-5,5);
  t= "#phi(j_{1}) "            +post;    phiJet1=td->make<TH1D>("phiJet1", t.c_str(),30,-3.14159,3.14159);
  t= "#phi(j_{2}) "            +post;    phiJet2=td->make<TH1D>("phiJet2", t.c_str(),30,-3.14159,3.14159);

  t="Missing E_{T} "+post;        ptMET=td->make<TH1D>("ptMET",t.c_str(),50,0.,100.);
  t="#phi (Missing E_{T}) "+post; phiMET=td->make<TH1D>("phiMET",t.c_str(),30,-3.14159,3.14159);
  
  // ----------  Composite histograms  ----------
  t="M(W_{R}) "                    +post;        mWR=td->make<TH1D>("mWR",   t.c_str(),70,0,2800);
  t="M(N_{R}) with #mu_{1} "       +post;      mNuR1=td->make<TH1D>("mNuR1", t.c_str(),70,0,2800);
  t="M(N_{R}) with #mu_{2} "       +post;      mNuR2=td->make<TH1D>("mNuR2", t.c_str(),70,0,1400);
  t="M(W_{R}) "                    +post;    mWR_raw=td->make<TH1D>("mWR_raw",   t.c_str(),70,0,2800);
  t="M(N_{R}) with #mu_{1} "       +post;  mNuR1_raw=td->make<TH1D>("mNuR1_raw", t.c_str(),70,0,2800);
  t="M(N_{R}) with #mu_{2} "       +post;  mNuR2_raw=td->make<TH1D>("mNuR2_raw", t.c_str(),70,0,1400);
  t="M(N_{R}) #mu_{1} vs. #mu_{2} "+post;     mNuR2D=td->make<TH2D>("mNuR2D",t.c_str(),70,0,2800,70,0,1400);
  t="M(N_{R}) #mu_{1} vs. #mu_{2} "+post; mNuR2D_raw=td->make<TH2D>("mNuR2D_raw",t.c_str(),70,0,2800,70,0,1400);

  mNuR2D->Sumw2() ; 

  t="M(#mu #mu)"                   +post;         mMuMu=td->make<TH1D>("mMuMu",    t.c_str(),50,0,2000);
  t="M(#mu #mu)"                   +post;     mMuMuZoom=td->make<TH1D>("mMuMuZoom",t.c_str(),50,0,200);
  t="M(jj)"                        +post;           mJJ=td->make<TH1D>("mJJ",      t.c_str(),50,0,2000);
  t="M(#mu #mu)"                   +post;     mMuMu_raw=td->make<TH1D>("mMuMu_raw",    t.c_str(),50,0,2000);
  t="M(#mu #mu)"                   +post; mMuMuZoom_raw=td->make<TH1D>("mMuMuZoom_raw",t.c_str(),50,0,200);
  t="M(jj)"                        +post;       mJJ_raw=td->make<TH1D>("mJJ_raw",      t.c_str(),50,0,2000);

  t="ECAL iso vs. #mu p_{T} > 100 GeV "      +post; ecalIso_muPt100=td->make<TH2D>("ecalIso_muPt100",t.c_str(),100,0.,500.,100,100.,1000.) ; 
  t="HCAL iso vs. #mu p_{T} > 100 GeV "      +post; hcalIso_muPt100=td->make<TH2D>("hcalIso_muPt100",t.c_str(),100,0.,500.,100,100.,1000.) ; 
  t="H/E iso vs. #mu p_{T} > 100 GeV "       +post; hoeIso_muPt100=td->make<TH2D>("hoeIso_muPt100",t.c_str(),100,0.,10.,100,100.,1000.) ; 
  ecalIso_muPt100->Sumw2() ; hcalIso_muPt100->Sumw2() ; hoeIso_muPt100->Sumw2() ; 
  t="ECAL iso vs. #mu p_{T} > 200 GeV "      +post; ecalIso_muPt200=td->make<TH2D>("ecalIso_muPt200",t.c_str(),100,0.,500.,100,100.,1000.) ; 
  t="HCAL iso vs. #mu p_{T} > 200 GeV "      +post; hcalIso_muPt200=td->make<TH2D>("hcalIso_muPt200",t.c_str(),100,0.,500.,100,100.,1000.) ; 
  t="H/E iso vs. #mu p_{T} > 200 GeV "       +post; hoeIso_muPt200=td->make<TH2D>("hoeIso_muPt200",t.c_str(),100,0.,10.,100,100.,1000.) ; 
  ecalIso_muPt200->Sumw2() ; hcalIso_muPt200->Sumw2() ; hoeIso_muPt200->Sumw2() ; 

}// end of book()

//======================================================================

void MuJetBackground::HistPerDef::fill(pat::Muon& theMuon,
				       pat::Jet& theJet, 
				       pat::MET& theMET, 
				       bool isMC,
				       double trkIsoLimit)
{  
  // Muons 
  ptMu1->Fill(theMuon.pt()) ; 
  etaMu1->Fill(theMuon.eta()) ; 
  phiMu1->Fill(theMuon.phi()) ; 

  // Jets 
  ptJet1->Fill(theJet.pt()) ; 
  etaJet1->Fill(theJet.eta()) ; 
  phiJet1->Fill(theJet.phi()) ; 

  // MET 
  ptMET->Fill(theMET.pt()) ; 
  phiMET->Fill(theMET.phi()) ; 

  // Muon relative track isolation as a function of muon pT
  double relTrkIso = theMuon.trackIso()/theMuon.pt() ;
  ptMu_etaMu_trkIso->Fill(theMuon.pt(),theMuon.eta(),relTrkIso) ; 
  
  if ( (theMuon.trackIso()/theMuon.pt()) < trkIsoLimit ) { // Relative track isolation
    ptMu1trkIso->Fill(theMuon.pt()) ; 
    etaMu1trkIso->Fill(theMuon.eta()) ; 
    phiMu1trkIso->Fill(theMuon.phi()) ;

    // Examine some of the high pT muon events
    if ( theMuon.pt() > 100.0 ) {
        ecalIso_muPt100->Fill( theMuon.ecalIso(),theMuon.pt() ) ; 
        hcalIso_muPt100->Fill( theMuon.hcalIso(),theMuon.pt() ) ; 
        hoeIso_muPt100->Fill( (theMuon.hcalIso()/theMuon.ecalIso()),theMuon.pt() ) ; 
        if ( theMuon.pt() > 200.0 ) {
            ecalIso_muPt200->Fill( theMuon.ecalIso(),theMuon.pt() ) ; 
            hcalIso_muPt200->Fill( theMuon.hcalIso(),theMuon.pt() ) ; 
            hoeIso_muPt200->Fill( (theMuon.hcalIso()/theMuon.ecalIso()),theMuon.pt() ) ; 
        }
    }
  }
}// end of fill()

//======================================================================

void
MuJetBackground::HistPerDef::fill(const HeavyNuEvent& hne, double w1, double w2)
{
  double weight = w1 * w2 ; 

  // Muons 
  if (hne.nMuons > 0) {
    ptMu1->Fill(hne.mu1.pt(),hne.eventWgt) ;     
    etaMu1->Fill(hne.mu1.eta(),hne.eventWgt) ; 
    phiMu1->Fill(hne.mu1.phi(),hne.eventWgt) ; 

    if (hne.nMuons > 1 ) {  
      ptMu2->Fill(hne.mu2.pt(),hne.eventWgt) ; 
      etaMu2->Fill(hne.mu2.eta(),hne.eventWgt) ; 
      phiMu2->Fill(hne.mu2.phi(),hne.eventWgt) ; 
    }
  }
  
  // Jets 
  if (hne.nJets > 0) {
    ptJet1->Fill(hne.j1.pt(),hne.eventWgt) ; 
    etaJet1->Fill(hne.j1.eta(),hne.eventWgt) ; 
    phiJet1->Fill(hne.j1.phi(),hne.eventWgt) ; 

    if (hne.nJets > 1) {
      ptJet2->Fill(hne.j2.pt(),hne.eventWgt) ; 
      etaJet2->Fill(hne.j2.eta(),hne.eventWgt) ; 
      phiJet2->Fill(hne.j2.phi(),hne.eventWgt) ; 

      mNuR1->Fill ( hne.mNuR1,w1*hne.eventWgt ) ; 
      mJJ->Fill   ( hne.mJJ,weight*hne.eventWgt );

      mNuR1_raw->Fill ( hne.mNuR1,hne.eventWgt ) ; 
      mJJ_raw->Fill   ( hne.mJJ,hne.eventWgt   );

      if (hne.nMuons > 1) {
        mWR->Fill   ( hne.mWR,weight*hne.eventWgt   ) ; 
        mNuR2->Fill ( hne.mNuR2,w2*hne.eventWgt ) ; 
        mNuR2D->Fill( hne.mNuR1, hne.mNuR2,weight*hne.eventWgt );

        mWR_raw->Fill   ( hne.mWR,hne.eventWgt   ) ; 
        mNuR2_raw->Fill ( hne.mNuR2,hne.eventWgt ) ; 
        mNuR2D_raw->Fill( hne.mNuR1, hne.mNuR2,hne.eventWgt );
      }
    }

  }

  if ( hne.nMuons > 1 ) { 
    mMuMu->Fill( hne.mLL,weight*hne.eventWgt );
    mMuMuZoom->Fill( hne.mLL,weight*hne.eventWgt );

    mMuMu_raw->Fill( hne.mLL,hne.eventWgt );
    mMuMuZoom_raw->Fill( hne.mLL,hne.eventWgt );
  }
}// end of fill()

//======================================================================

//
// constants, enums and typedefs
//
const std::vector<hNuMassHypothesis> v_null;

//
// static data member definitions
//

//======================================================================

//
// constructors and destructor
//
MuJetBackground::MuJetBackground(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  dolog_=iConfig.getParameter<bool>("DoLog");

  muonTag_ = iConfig.getParameter< edm::InputTag >( "muonTag" );
  jetTag_  = iConfig.getParameter< edm::InputTag >( "jetTag"  );
  metTag_  = iConfig.getParameter< edm::InputTag >( "metTag"  );
  jptTag_  = iConfig.getParameter< edm::InputTag >( "jptTag"  );

  trig_ = new HeavyNuTrigger(iConfig.getParameter<edm::ParameterSet>("trigMatchPset"));
  muid_ = new HeavyNuID(iConfig.getParameter<edm::ParameterSet>("muIDPset"));

  edm::Service<TFileService> fs;
  hists.nmu      = fs->make<TH1D>("nmu",   "N(#mu^{#pm})",10,-0.5,9.5);
  hists.njet     = fs->make<TH1D>("njet",  "N(Jet)",50,-0.5,49.5);

  // trigEtaLimit_ = iConfig.getParameter<double>("trigEtaLimit") ; 
  // Histos per cut:
  calcSurvival_ = iConfig.getParameter<bool>("getSurvivalRate") ;
//   doClosure_    = iConfig.getParameter<bool>("doClosureTest") ;
//   doQuadJet_    = iConfig.getParameter<bool>("doQuadJetTest") ;

  if ( calcSurvival_ ) { 
    hists.NoCuts.book        ( new TFileDirectory(fs->mkdir("NoCuts")), 
                               "(mu+jet)", v_null );
    hists.dPhiCuts.book      ( new TFileDirectory(fs->mkdir("dPhiCuts")), 
                               "(Basic)", v_null );
    hists.dPhiJetCuts.book   ( new TFileDirectory(fs->mkdir("dPhiJetCuts")), 
                               "(Jet overlap)", v_null );
    hists.dPhiCoreCuts.book  ( new TFileDirectory(fs->mkdir("dPhiCoreCuts")), 
                               "(Jet core overlap)", v_null );
    hists.dPhi10pctCuts.book ( new TFileDirectory(fs->mkdir("dPhi10pctCuts")), 
                               "(10+% calorimeter energy)", v_null );
    hists.dPhiTrigCuts.book      ( new TFileDirectory(fs->mkdir("dPhiTrigCuts")), 
                                   "(Basic, Trigger)", v_null );
    hists.dPhiJetTrigCuts.book   ( new TFileDirectory(fs->mkdir("dPhiJetTrigCuts")), 
                                   "(Jet overlap, Trigger)", v_null );
    hists.dPhiCoreTrigCuts.book  ( new TFileDirectory(fs->mkdir("dPhiCoreTrigCuts")), 
                                   "(Jet core overlap, Trigger)", v_null );
    hists.dPhi10pctTrigCuts.book ( new TFileDirectory(fs->mkdir("dPhi10pctTrigCuts")), 
                                   "(10+% calorimeter energy, Trigger)", v_null );
    hists.dPhiTightCuts.book      ( new TFileDirectory(fs->mkdir("dPhiTightCuts")), 
                                    "(Tight)", v_null );
    hists.dPhiJetTightCuts.book   ( new TFileDirectory(fs->mkdir("dPhiJetTightCuts")), 
                                    "(Jet overlap, Tight)", v_null );
    hists.dPhiCoreTightCuts.book  ( new TFileDirectory(fs->mkdir("dPhiCoreTightCuts")), 
                                    "(Jet core overlap, Tight)", v_null );
    hists.dPhi10pctTightCuts.book ( new TFileDirectory(fs->mkdir("dPhi10pctTightCuts")), 
                                    "(10+% calorimeter energy, Tight)", v_null );
    hists.dPhiTightTrigCuts.book      ( new TFileDirectory(fs->mkdir("dPhiTightTrigCuts")), 
                                        "(Tight, Trigger)", v_null );
    hists.dPhiJetTightTrigCuts.book   ( new TFileDirectory(fs->mkdir("dPhiJetTightTrigCuts")), 
                                        "(Jet overlap, Tight, Trigger)", v_null );
    hists.dPhiCoreTightTrigCuts.book  ( new TFileDirectory(fs->mkdir("dPhiCoreTightTrigCuts")), 
                                        "(Jet core overlap, Tight, Trigger)", v_null );
    hists.dPhi10pctTightTrigCuts.book ( new TFileDirectory(fs->mkdir("dPhi10pctTightTrigCuts")), 
                                        "(10+% calorimeter energy, Tight, Trigger)", v_null );
  } 
//   if ( doQuadJet_ ) { 
//     hists.LLJJpTCuts.book   ( new TFileDirectory(fs->mkdir("LLJJpTCuts")), 
//                               "(two muons, two jets)", v_null );
//     hists.TrigMatches.book  ( new TFileDirectory(fs->mkdir("TrigMatches")), 
//                               "(at least one trigger matched muon)", v_null );
//     hists.VertexCuts.book   ( new TFileDirectory(fs->mkdir("VertexCuts")), 
//                               "(all objects share common vtx)", v_null );
//     hists.Mu1HighPtCut.book ( new TFileDirectory(fs->mkdir("Mu1HighPtCut")), 
//                               "(passes mu1 pT rqmt)", v_null );
//     hists.diLmassCut.book   ( new TFileDirectory(fs->mkdir("diLmassCuts")), 
//                               "(passes mumu mass cut)", v_null );
//     hists.mmJClosure.book   ( new TFileDirectory(fs->mkdir("mmJClosure")), 
//                               "(1+ jets, 2 fake muons)", v_null );
//     hists.mmJClosureMETlt20.book   ( new TFileDirectory(fs->mkdir("mmJClosureMETlt20")), 
//                               "(1+ jets, 2 fake muons, MET < 20 GeV)", v_null );
//     hists.mMuJClosure.book  ( new TFileDirectory(fs->mkdir("mMuJClosure")), 
//                               "(1+ jets, one fake muon, one muon)", v_null );
//     hists.mMuJClosureMETlt20.book  ( new TFileDirectory(fs->mkdir("mMuJClosureMETlt20")), 
//                               "(1+ jets, one fake muon, one muon, MET < 20)", v_null );
//   }
//   if ( doClosure_ ) { 
//     hists.LJJJClosure.book ( new TFileDirectory(fs->mkdir("LJJJClosure")), 
//                              "(3 jets, 1 muon in jet)", v_null );
//     hists.L2JClosure.book  ( new TFileDirectory(fs->mkdir("L2JClosure")), 
//                              "(2 jets, tight muon)", v_null );
//   }

  hists.rundir = new TFileDirectory(fs->mkdir("RunDir"));

  init_=false;

  cuts.minimum_mu1_pt       = iConfig.getParameter<double>("minMu1pt");
  cuts.minimum_mu2_pt       = iConfig.getParameter<double>("minMu2pt");
  cuts.minimum_jet_pt       = iConfig.getParameter<double>("minJetPt");
  cuts.maximum_mu_abseta    = iConfig.getParameter<double>("maxMuAbsEta");
  cuts.maximum_jet_abseta   = iConfig.getParameter<double>("maxJetAbsEta");
  // cuts.minimum_mumu_mass    = iConfig.getParameter<double>("minMuMuMass");
  cuts.minimum_muon_jet_dR  = iConfig.getParameter<double>("minMuonJetdR");
  cuts.muon_trackiso_limit  = iConfig.getParameter<double>("muonTrackRelIsoLimit");
  // cuts.maxVertexZsep        = iConfig.getParameter<double>("dimuonMaxVertexZsepCM");
  // cuts.maxJetVZsepCM        = iConfig.getParameter<double>("maxJetVZsepCM");
  cuts.minimum_dijet_dPhi   = iConfig.getParameter<double>("minimumMuJetdPhi");
  cuts.minimum_dijet_pt     = iConfig.getParameter<double>("minimumJetPtForDijets");
  cuts.minimum_extraJet_dR  = iConfig.getParameter<double>("minimumDeltaRforExtraJets");
  // cuts.minimum_SCEt         = iConfig.getParameter<double>("minimumSuperClusterEt");

//   theMETtype  = iConfig.getParameter<int>("METvariety") ; 
//   rwLowPtbin  = iConfig.getParameter< std::vector<double> >("reweightPtLow") ; 
//   rwHighPtbin = iConfig.getParameter< std::vector<double> >("reweightPtHigh") ; 
//   rw2011A     = iConfig.getParameter< std::vector<double> >("reweight2011A") ; 
//   rw2011B     = iConfig.getParameter< std::vector<double> >("reweight2011B") ; 

//   // Special check: Make sure all vectors are of the same size
//   unsigned int vecsize = rwLowPtbin.size() ; 
//   if ( ( doClosure_ || doQuadJet_ ) && 
//        ( rwHighPtbin.size() != vecsize ||
//          rw2011A.size() != vecsize ||
//  	 rw2011B.size() != vecsize ) )
//     throw cms::Exception( "Please ensure that all QCD reweighting vectors are equal size");
  
  pileupEra_ = iConfig.getParameter<int>("pileupEra");
//   if ( pileupEra_ == 20111 ) { // 2011 A
//       rwQCD = rw2011A ;
//   } else if ( pileupEra_ == 20112 ) { // 2011 B
//       rwQCD = rw2011B ;
//   } else {
//       std::cout << "WARNING: Unknown era for QCD corrections requested.  Assigning 2011 A" << std::endl ;
//       rwQCD = rw2011A ;
//   }      
  
  isPFJets_ = iConfig.getParameter<bool>("isPFJets") ; 
  MCweightByVertex_ = edm::LumiReWeighting(hnu::get_standard_pileup_mc(pileupEra_),hnu::get_standard_pileup_data(pileupEra_));

  // For the record...
  std::cout << "Configurable cut values applied:" << std::endl;
  std::cout << "muonTag          = " << muonTag_                 << std::endl;
  // std::cout << "trigEtaLimit     = " << trigEtaLimit_            << std::endl;
  std::cout << "jetTag           = " << jetTag_                  << std::endl;
  std::cout << "metTag           = " << metTag_                  << std::endl;
  std::cout << "jptTag           = " << jptTag_                  << std::endl;
  std::cout << "minMu1pt         = " << cuts.minimum_mu1_pt      << " GeV" << std::endl;
  std::cout << "minMu2pt         = " << cuts.minimum_mu2_pt      << " GeV" << std::endl;
  std::cout << "minJetPt         = " << cuts.minimum_jet_pt      << " GeV" << std::endl;
  std::cout << "maxMuAbsEta      = " << cuts.maximum_mu_abseta   << std::endl;
  std::cout << "maxJetAbsEta     = " << cuts.maximum_jet_abseta  << std::endl;
  std::cout << "minMuonJetdR     = " << cuts.minimum_muon_jet_dR << std::endl;
  std::cout << "muonTrackIso     = " << 100 * cuts.muon_trackiso_limit << "%" << std::endl;
  std::cout << "minimumMuJetdPhi = " << cuts.minimum_dijet_dPhi  << std::endl; 
  std::cout << "minimumQCDjetPt  = " << cuts.minimum_dijet_pt   << " GeV " << std::endl; 
  std::cout << "minExtraJetdR    = " << cuts.minimum_extraJet_dR << std::endl;
  // std::cout << "dimuonMaxVertexZsepCM = " << cuts.maxVertexZsep << std::endl ;
  // std::cout << "maxJetVZsepCM         = " << cuts.maxJetVZsepCM << std::endl ;
  std::cout << "pileup era       = " << pileupEra_ << std::endl;

}
  
MuJetBackground::~MuJetBackground()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

//======================================================================

//
// member functions
//
//======================================================================

//======================================================================

TH1 *
MuJetBackground::bookRunHisto(uint32_t runNumber)
{
  std::string runstr = int2str<uint32_t>(runNumber);
  return hists.rundir->make <TH1I> (runstr.c_str(), runstr.c_str(),1,1,2);
}

//======================================================================

// bool 
// MuJetBackground::selectJets(std::vector< std::pair<pat::Jet,float> >& jets,
// 			    HeavyNuEvent& hne, int minNjets) {
//   for (unsigned int i=0; i<jets.size(); i++) { 
//     if ( hne.nJets == 2 ) return true ; 
//     pat::Jet iJ = jets.at(i).first ;

//     // Must have jets of sufficient transverse momentum
//     if ( iJ.pt() < cuts.minimum_jet_pt ) continue ; 

//     // Jets must be separated from candidate muon(s)
//     double dRm1j = deltaR(iJ.eta(), iJ.phi(), hne.mu1.eta(), hne.mu1.phi()) ; 
//     double dRm2j = ( hne.nMuons > 1 ) ? 
//       ( deltaR(iJ.eta(), iJ.phi(), hne.mu2.eta(), hne.mu2.phi()) ) : ( 100.0 ) ; 
//     if (std::min(dRm1j,dRm2j) > cuts.minimum_muon_jet_dR) {
//       hne.nJets++ ; 
//       if      ( hne.nJets == 1 ) {
//           hne.j1 = iJ ;
//           hne.j1scale = 1.0 ;
//       }
//       else if ( hne.nJets == 2 ) {
//           hne.j2 = iJ ;
//           hne.j2scale = 1.0 ;
//       }
//       else    std::cout << "WARNING: Expected empty jet position" << std::endl ; 
//     }
//   }
//   return ( hne.nJets >= minNjets ) ; 
// }


bool
MuJetBackground::findQCDjet(const std::vector< std::pair<pat::Jet,float> >& jets,
			    edm::Handle<pat::MuonCollection>& pMuons,
			    HeavyNuEvent& hne) {
  pat::Jet iJqcd ; 
  bool foundJet = false ; unsigned int jetLoc = jets.size() ;  
  float maxdPhi_jet_mu = cuts.minimum_dijet_dPhi ; 
  for (unsigned int i=0; i<jets.size(); i++) {
    pat::Jet iJ = jets.at(i).first ; 
    float dPhi = fabs( deltaPhi(iJ.phi(),hne.mu1.phi()) ) ; 
    if ( dPhi > maxdPhi_jet_mu ) { 
      foundJet = true ; 
      maxdPhi_jet_mu = dPhi ; 
      iJqcd = iJ ; 
      jetLoc = i ; 
    } 
  }
  if ( !foundJet ) return false ; 

  // Reject the event if extra high-energy jets are present outside mu/jet axis
  for (unsigned int i=0; i<jets.size(); i++) {
    if ( i == jetLoc ) continue ; 
    pat::Jet iJ = jets.at(i).first ; 
    if ( iJ.pt() > cuts.minimum_jet_pt ) {  
      float dRm = deltaR(iJ.eta(),iJ.phi(),hne.mu1.eta(),hne.mu1.phi()) ; 
      float dRj = deltaR(iJ.eta(),iJ.phi(),iJqcd.eta(),iJqcd.phi()) ; 
      if ( std::min(dRm,dRj) > cuts.minimum_extraJet_dR ) return false ; 
    }
  }

  // Reject event if muon found with most of jet energy 
  for (unsigned int i=0; i<pMuons->size(); i++) {
    pat::MuonRef iM=pat::MuonRef( pMuons,i ) ;
    float dR  = deltaR(iM->eta(),iM->phi(),hne.mu1.eta(),hne.mu1.phi()) ; 
    if ( dR < 0.001 ) continue ; // Matched to original muon
    dR = deltaR(iM->eta(),iM->phi(),iJqcd.eta(),iJqcd.phi()) ; 
    if ( dR < 0.5 ) { // muon is in the jet 
      double pTratio = iM->pt() / iJqcd.pt() ; 
      if ( pTratio > 0.75 ) return false ; 
    }
  }

  hne.j1 = iJqcd ; 
  return true ; 
} //MuJetBackground::findQCDjet

//======================================================================

void MuJetBackground::selectMuonPairs(std::vector<pat::Muon>& muons,
                                      HeavyNuEvent& hne) {

    hne.mu1 = muons.at(0) ;
    hne.nMuons = 1 ; 
    for (unsigned int i=1; i<muons.size(); i++) {
        double dR = deltaR(hne.mu1.eta(), hne.mu1.phi(), muons.at(i).eta(), muons.at(i).phi()) ;
        std::cout << "Delta R between muon 0 and muon " << i << " is " << dR << std::endl ; 
        if ( dR > 0.3 ) { // Require some separation between muons
            hne.mu2 = muons.at(i) ;
            hne.nMuons++ ; 
            return ; 
        }
    }
}

void
MuJetBackground::selectMuonsInJets(std::vector<pat::Muon>& muons,
				   std::vector< std::pair<pat::Jet,float> >& jets,
				   HeavyNuEvent& hne) {
  double mu1wgt = 1.0 ; 
  double mu2wgt = 1.0 ; 
  for (unsigned int i=0; i<muons.size(); i++) {
    pat::Muon iM = muons.at(i) ; 
    bool muInJet = false ; 
    for (unsigned int j=0; j<jets.size(); j++) {
      if ( muInJet ) break ; 
      pat::Jet iJ = jets.at(j).first;
      if ( iJ.pt() < cuts.minimum_jet_pt ) continue ;
      if ( hnu::jetID(iJ) < 1 ) continue ; 
      double dR = deltaR(iJ.eta(),iJ.phi(),iM.eta(),iM.phi()) ;
      if ( dR < cuts.minimum_muon_jet_dR ) muInJet = true ; 
    }

    if ( muInJet ) { 
      hne.nMuons++ ; 
      if ( hne.nMuons == 1 ) { 
	hne.mu1 = iM ; 
	if ( hne.isMC && applyMuIDCorrections_ ) mu1wgt = muid_->weightForMC( iM.pt(),0 ) ; 
      } 
      if ( hne.nMuons == 2 ) { 
	hne.mu2 = iM ; 
	if ( hne.isMC && applyMuIDCorrections_ ) mu2wgt = muid_->weightForMC( iM.pt(),0 ) ; 
      } 
    }
  }

  hne.eventWgt *= mu1wgt * mu2wgt ;
}

std::vector<pat::Muon> MuJetBackground::getNonIsolatedMuons(std::vector<pat::Muon>& muons) {

    std::vector<pat::Muon> dirtyMuons ;
    for (unsigned int i=0; i<muons.size(); i++) {
        pat::Muon iM = muons.at(i) ; 
        bool isDirty = false ;
        if ( iM.pt() < 100.0 ) isDirty = ( (iM.ecalIso()+iM.hcalIso()) > 10. ) ; 
        else                   isDirty = ( ((iM.ecalIso()+iM.hcalIso())/iM.pt()) > 0.10 ) ;
        if ( isDirty ) dirtyMuons.push_back( iM ) ;
    }
    if ( dirtyMuons.size() > 1 ) 
        std::sort( dirtyMuons.begin(),dirtyMuons.end(),compare() ) ; 

    return dirtyMuons ;
}

bool MuJetBackground::secondQualityMuon(const pat::Muon& mu1, const pat::Muon& mu2) { 

  if ( !hnu::isLooseMuonNoPF(mu2) ) return false ; 

  float relIso = ( mu2.trackIso() + mu2.hcalIso() + mu2.ecalIso() ) / mu2.pt() ;
  if ( relIso < 0.15 ) return true ;

  reco::Particle::LorentzVector vMuMu = mu1.p4() + mu2.p4() ;
  if ( vMuMu.M() > 70. && vMuMu.M() < 110. ) return true ; 

  return false ; 
}

bool
MuJetBackground::findQCDmuon(const std::vector<pat::Muon>& muons, 
			     edm::Handle<pat::MuonCollection>& pMuons,
			     HeavyNuEvent& hne) {

  pat::Muon qcdCand = muons.at(0) ; 
  for (unsigned int i=1; i<muons.size(); i++) { // Skip first candidate
    if ( secondQualityMuon(qcdCand,muons.at(i)) ) return false ;  
  }
  for (size_t iMuon=0; iMuon<pMuons->size(); iMuon++) {
    pat::MuonRef iM=pat::MuonRef(pMuons,iMuon);
    if ( deltaR(qcdCand.eta(),qcdCand.phi(),iM->eta(),iM->phi()) < 0.001 ) continue ; // same muon 
    if ( secondQualityMuon(qcdCand,(*iM)) ) return false ;  
  }
  hne.mu1 = qcdCand ; 
  return true ; 
} // MuJetBackground::findQCDmuon

bool MuJetBackground::isDijetCandidate(HeavyNuEvent& hne,pat::MET& theMET) {

  if ( theMET.et() > 20. ) return false ; // Absolute MET requirement
  if ( hnu::jetID(hne.j1) < 1 ) return false ; 
  return true ; 

} // MuJetBackground::isDijetCandidate

//======================================================================

// ------------ method called to for each event  ------------
bool
MuJetBackground::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  bool keepThisEvent = false ; 

  applyMuIDCorrections_ = !iEvent.isRealData();

  if (iEvent.isRealData())
  {
    uint32_t runn = iEvent.id().run();
    std::map<uint32_t,TH1 *>::const_iterator it = m_runHistos_.find(runn);
    TH1 *runh;
    if (it == m_runHistos_.end()) {
      runh = bookRunHisto(runn);
      m_runHistos_[runn] = runh;
    } else
      runh = it->second;
    runh->Fill(1);
  }

  edm::Handle<pat::MuonCollection> pMuons ; 
  iEvent.getByLabel(muonTag_,pMuons) ; 

  edm::Handle<pat::JetCollection> pJets ;
  iEvent.getByLabel(jetTag_, pJets) ;

  edm::Handle<reco::JPTJetCollection> jptJets;
  iEvent.getByLabel(jptTag_, jptJets); 
    
  // edm::Handle<reco::SuperClusterCollection> hybridClusters ; 
  // iEvent.getByLabel(hybridSClabel_, hybridClusters) ; 

  // edm::Handle<reco::SuperClusterCollection> multi5x5Clusters ; 
  // iEvent.getByLabel(multiSClabel_, multi5x5Clusters) ; 

  edm::Handle<std::vector<PileupSummaryInfo> > pPU;
  iEvent.getByLabel("addPileupInfo", pPU);    

  edm::Handle<reco::VertexCollection> pvHandle;
  iEvent.getByLabel("offlinePrimaryVertices", pvHandle);

  if ( !pMuons.isValid() || 
       !pJets.isValid() ) {
    std::cout << "Exiting as valid PAT objects not found: " 
	      << pMuons.isValid() << " " 
	      << pJets.isValid() << " " 
	      << std::endl ;
    return false; 
  }

//   edm::Handle<pat::METCollection> patMetCollection ; 
//   iEvent.getByLabel("patMETs", patMetCollection) ; 

  edm::Handle<reco::PFMETCollection> pfMetCollection ; 
  iEvent.getByLabel(metTag_, pfMetCollection) ;

  //Shirpa reweighting info
  double genweight = 1.0;
  edm::Handle<GenEventInfoProduct> geneventinfo;
  iEvent.getByLabel("generator", geneventinfo);
  if(!iEvent.isRealData()) genweight = geneventinfo->weight();

  if ( !pfMetCollection.isValid() ) { 
    std::cout << "Exiting as valid MET collection not found" << std::endl ; 
    return false ; 
  }
  std::auto_ptr<pat::METCollection> pMET(new pat::METCollection) ; 
  reco::MET thePFMET(pfMetCollection.product()->at(0).p4(),pfMetCollection.product()->at(0).vertex()) ; 
  pMET->push_back(pat::MET(thePFMET)) ;

  if (firstEvent_) {
    // handle the jet corrector parameters collection,
    // get the jet corrector parameters collection from the global tag
    //
    edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
    if ( isPFJets_ ) iSetup.get<JetCorrectionsRecord > ().get("AK5PFchs", JetCorParColl) ; 
    else             iSetup.get<JetCorrectionsRecord > ().get("AK5Calo", JetCorParColl) ;
    
    // get the uncertainty parameters from the collection,
    // instantiate the jec uncertainty object
    //
//     if ( !iEvent.isRealData() ) {
//         int pileupYear = pileupEra_ ;
//         int idYear     = muid_->idEra() ;
        
//         bool allErasMatch = ( pileupYear == idYear ) ;
//         if ( !allErasMatch ) {
//             std::cout << "WARNING: You do not appear to have consistent corrections applied!" << std::endl ;
//             std::cout << "         pileup year is " << pileupEra_ << ", year for mu ID is " << idYear
// 		      << std::endl ; 
//         } else {
//             std::cout << "Looking at corrections, I assume you are running with the " << pileupYear << " year settings" << std::endl ; 
//         }
//         std::cout << "==================================" << std::endl ; 
//     }  
    firstEvent_ = false;
  }

  hists.nmu  ->Fill(pMuons->size()) ;
  hists.njet ->Fill(pJets->size()) ;

  int muonID = ( (muid_->idEra() != 0) ? (muid_->idEra()/abs(muid_->idEra())) : 0 ) ; 
  
  std::vector<pat::Muon> muCands = 
    hnu::getMuonList(pMuons,pvHandle,muonID,cuts.minimum_mu2_pt,cuts.maximum_mu_abseta,1.0) ;
  
  std::vector< std::pair<pat::Jet,float> > jetCands = 
    hnu::getJetList(pJets,jecuObj_,cuts.minimum_dijet_pt,cuts.maximum_jet_abseta,0) ; 
  if ( muCands.size() < 1 || jetCands.size() < 1 ) return false ;
  
  if ( calcSurvival_ ) { 
    HeavyNuEvent hnuDijet(HeavyNuEvent::QCD) ;  
    initializeHNE(hnuDijet,pPU,pvHandle,!iEvent.isRealData(),isPFJets_, genweight) ;
    hists.NoCuts.fill( muCands.at(0),jetCands.at(0).first,pMET->at(0),!iEvent.isRealData(),
                       cuts.muon_trackiso_limit ) ; 

    if ( findQCDmuon( muCands, pMuons, hnuDijet ) ) { 
      if ( findQCDjet( jetCands, pMuons, hnuDijet ) ) { 
	if ( isDijetCandidate( hnuDijet,pMET->at(0) ) ) { 
	  bool mu1trig = false ; 
	  if ( trig_->matchingEnabled() && iEvent.isRealData() ) {
            mu1trig = trig_->isTriggerMatched( hnuDijet.mu1, iEvent ) ; 
	  } else if ( !iEvent.isRealData() ) {
            if ( fabs(hnuDijet.mu1.eta()) > 2.1 && fabs(hnuDijet.mu1.eta()) < 2.4 ) mu1trig = true ; 
	    else mu1trig = trig_->simulateForMC( hnuDijet.mu1.pt(),hnuDijet.mu1.eta(),0 );
	  }

          reco::TrackRef cktTrack = (muon::tevOptimized(hnuDijet.mu1, 200, 40., 17., 0.25)).first;
          bool isTightMuon = ( (!cktTrack.isNull()) ? hnu::isTightHighPtMuon(hnuDijet.mu1, pvHandle) : false ) ; 
          
          hists.dPhiCuts.fill( hnuDijet.mu1,hnuDijet.j1,pMET->at(0),!iEvent.isRealData(),
                               cuts.muon_trackiso_limit ) ;
          if ( isTightMuon ) hists.dPhiTightCuts.fill( hnuDijet.mu1,hnuDijet.j1,pMET->at(0),!iEvent.isRealData(),
                                                       cuts.muon_trackiso_limit ) ;
          if ( mu1trig ) {
            hists.dPhiTrigCuts.fill( hnuDijet.mu1,hnuDijet.j1,pMET->at(0),!iEvent.isRealData(),
                                     cuts.muon_trackiso_limit ) ; 
            if ( isTightMuon ) hists.dPhiTightTrigCuts.fill( hnuDijet.mu1,hnuDijet.j1,pMET->at(0),!iEvent.isRealData(),
                                                             cuts.muon_trackiso_limit ) ; 
          }

          // Check: Does the muon overlap with a jet?
          bool jetOverlap = false, jetOverlapCore = false ; 
          for (unsigned int i=0; i<jetCands.size(); i++) {
            pat::Jet iJ = jetCands.at(i).first ; 
            if ( iJ.pt() < cuts.minimum_jet_pt ) continue ; 
            float dR = deltaR(iJ.eta(),iJ.phi(),hnuDijet.mu1.eta(),hnuDijet.mu1.phi()) ;
            if ( dR < 0.5 ) {
              jetOverlap = true ; 
              if ( dR < 0.3 ) {
                jetOverlapCore = true ; 
                break ;
              }
            }
          }
          
          if ( jetOverlap ) {
            hists.dPhiJetCuts.fill( hnuDijet.mu1,hnuDijet.j1,pMET->at(0),!iEvent.isRealData(),
                                   cuts.muon_trackiso_limit ) ; 
            if ( isTightMuon ) hists.dPhiJetTightCuts.fill( hnuDijet.mu1,hnuDijet.j1,pMET->at(0),!iEvent.isRealData(),
                                                           cuts.muon_trackiso_limit ) ; 
            if ( mu1trig ) {
              hists.dPhiJetTrigCuts.fill( hnuDijet.mu1,hnuDijet.j1,pMET->at(0),!iEvent.isRealData(),
                                         cuts.muon_trackiso_limit ) ; 
              if ( isTightMuon ) hists.dPhiJetTightTrigCuts.fill( hnuDijet.mu1,hnuDijet.j1,pMET->at(0),!iEvent.isRealData(),
                                                                  cuts.muon_trackiso_limit ) ; 
            }
            if ( jetOverlapCore ) { 
              hists.dPhiCoreCuts.fill( hnuDijet.mu1,hnuDijet.j1,pMET->at(0),!iEvent.isRealData(),
                                       cuts.muon_trackiso_limit ) ; 
              if ( isTightMuon ) hists.dPhiCoreTightCuts.fill( hnuDijet.mu1,hnuDijet.j1,pMET->at(0),!iEvent.isRealData(),
                                                               cuts.muon_trackiso_limit ) ; 
              if ( mu1trig ) {
                hists.dPhiCoreTrigCuts.fill( hnuDijet.mu1,hnuDijet.j1,pMET->at(0),!iEvent.isRealData(),
                                             cuts.muon_trackiso_limit ) ; 
                if ( isTightMuon ) hists.dPhiCoreTightTrigCuts.fill( hnuDijet.mu1,hnuDijet.j1,pMET->at(0),!iEvent.isRealData(),
                                                                     cuts.muon_trackiso_limit ) ; 
              }
            }
          }
          
          if ( ((hnuDijet.mu1.ecalIso()+hnuDijet.mu1.hcalIso())/hnuDijet.mu1.pt()) > 0.10 ) {
            hists.dPhi10pctCuts.fill( hnuDijet.mu1,hnuDijet.j1,pMET->at(0),!iEvent.isRealData(),
                                      cuts.muon_trackiso_limit ) ; 
            if ( isTightMuon ) hists.dPhi10pctTightCuts.fill( hnuDijet.mu1,hnuDijet.j1,pMET->at(0),!iEvent.isRealData(),
                                                              cuts.muon_trackiso_limit ) ; 
            if ( mu1trig ) {
              hists.dPhi10pctTrigCuts.fill( hnuDijet.mu1,hnuDijet.j1,pMET->at(0),!iEvent.isRealData(),
                                            cuts.muon_trackiso_limit ) ; 
              if ( isTightMuon ) hists.dPhi10pctTightTrigCuts.fill( hnuDijet.mu1,hnuDijet.j1,pMET->at(0),!iEvent.isRealData(),
                                                                    cuts.muon_trackiso_limit ) ; 
            }
          }
          
          if ( hnuDijet.mu1.pt() > 100.0 ) {
            keepThisEvent = true ;   
            std::cout << "Dijet event: " << iEvent.id() << std::endl;
            std::cout << "\tMuon: ";
            std::cout << "pt=" << hnuDijet.mu1.pt() << " GeV, eta=" << hnuDijet.mu1.eta() << ", phi=" << hnuDijet.mu1.phi();
            std::cout << ", H+E: " << (hnuDijet.mu1.ecalIso()+hnuDijet.mu1.hcalIso())
                      << ", (H+E)/pT: " << ((hnuDijet.mu1.ecalIso()+hnuDijet.mu1.hcalIso())/hnuDijet.mu1.pt())
		      << ", Trig: " << mu1trig  
                      << std::endl;
          }
	}
      }
    }
  }

  // Deprecated 3 Apr 2013 (BMD)
  // Other functionality, if useful, incorporated into HeavyNu
  /* 
  double MET = pMET->at(0).pt() ; 
  std::vector<pat::Muon> dirtyMuons = getNonIsolatedMuons( muCands ) ;
  if ( doQuadJet_ ) { 
      // if ( muCands.size() >= 2 && jetCands.size() >= 4 ) { 
      if ( dirtyMuons.size() >= 2 && jetCands.size() >= 2 ) {
          HeavyNuEvent hnuQuadjet(HeavyNuEvent::QCD) ;  
          initializeHNE(hnuQuadjet,pPU,pvHandle,!iEvent.isRealData(),isPFJets_,genweight) ;
          selectMuonPairs( dirtyMuons,hnuQuadjet ) ; 
//               hnuQuadjet.mu1 = dirtyMuons.at(0) ; 
//               hnuQuadjet.mu2 = dirtyMuons.at(1) ;
//               hnuQuadjet.nMuons = dirtyMuons.size() ; 
          if ( hnuQuadjet.isMC ) 
              hnuQuadjet.eventWgt *= muid_->weightForMC( hnuQuadjet.mu1.pt(),0 ) * muid_->weightForMC( hnuQuadjet.mu2.pt(),0 ) ; 
          if ( hnuQuadjet.nMuons > 1 && selectJets( jetCands,hnuQuadjet ) ) {
              if ( hnu::jetID(hnuQuadjet.j1) > 0 && hnu::jetID(hnuQuadjet.j2) > 0 ) { 
                  
                  hnuQuadjet.tjV1 = hnu::caloJetVertex(hnuQuadjet.j1, *jptJets);
                  hnuQuadjet.tjV2 = hnu::caloJetVertex(hnuQuadjet.j2, *jptJets);
                  
                  hnuQuadjet.scaleMuE() ; 
                  hnuQuadjet.regularize();
                  hnuQuadjet.calculate() ;
                  
                  double mu1scaleFactor, mu2scaleFactor ; 
                  mu1scaleFactor = qcdScaleFactor(hnuQuadjet.mu1.pt(),
                                                  rwLowPtbin,rwHighPtbin,rwQCD) ; 
                  mu2scaleFactor = qcdScaleFactor(hnuQuadjet.mu2.pt(),
                                                  rwLowPtbin,rwHighPtbin,rwQCD) ; 
                  
                  hists.LLJJpTCuts.fill( hnuQuadjet,mu1scaleFactor,mu2scaleFactor ) ; 
                  
                  bool mu1trig = false ; bool mu2trig = false ;
                  if ( trig_->matchingEnabled() && iEvent.isRealData() ) {
                      mu1trig = trig_->isTriggerMatched( hnuQuadjet.mu1, iEvent) ; 
                      mu2trig = trig_->isTriggerMatched( hnuQuadjet.mu2, iEvent) ; 
                  } else if ( !iEvent.isRealData() ) {
                      mu1trig = trig_->simulateForMC( hnuQuadjet.mu1.pt(),hnuQuadjet.mu1.eta(),0 );
                      mu2trig = trig_->simulateForMC( hnuQuadjet.mu2.pt(),hnuQuadjet.mu2.eta(),0 );
                  }
                  
                  if ( mu1trig || mu2trig ) { 
                      hists.TrigMatches.fill( hnuQuadjet,mu1scaleFactor,mu2scaleFactor ) ; 
                      
                      //--- Impose vertex requirement here ---//
                      float deltaVzJ1J2 = fabs(hnuQuadjet.tjV1-hnuQuadjet.tjV2);
                      float deltaVzJ1M1 = fabs(hnuQuadjet.tjV1-hnuQuadjet.mu1.vertex().Z());
                      float deltaVzJ2M2 = fabs(hnuQuadjet.tjV2-hnuQuadjet.mu2.vertex().Z());
                      float deltaVzJ1M2 = fabs(hnuQuadjet.tjV1-hnuQuadjet.mu2.vertex().Z());
                      float deltaVzJ2M1 = fabs(hnuQuadjet.tjV2-hnuQuadjet.mu1.vertex().Z());
                      float deltaVzM1M2 = fabs(hnuQuadjet.mu1.vertex().Z()-hnuQuadjet.mu2.vertex().Z());
                      if ( (cuts.maxJetVZsepCM < 0 || cuts.maxVertexZsep < 0) || 
                           ((deltaVzJ1J2 < cuts.maxJetVZsepCM) && (deltaVzJ1M1 < cuts.maxJetVZsepCM) &&
                            (deltaVzJ2M2 < cuts.maxJetVZsepCM) && (deltaVzJ1M2 < cuts.maxJetVZsepCM) &&
                            (deltaVzJ2M1 < cuts.maxJetVZsepCM) && (deltaVzM1M2 < cuts.maxVertexZsep)) ) {
                          
                          hists.VertexCuts.fill( hnuQuadjet,mu1scaleFactor,mu2scaleFactor ) ; 
                          
                          keepThisEvent = true ;
                          if(iEvent.isRealData()) {
                              std::cout << "\tQUADJET: " << iEvent.id() << std::endl;
                              std::cout << "\tM(W_R)  = " << hnuQuadjet.mWR << " GeV";
                              std::cout << ", M(NuR1) = " << hnuQuadjet.mNuR1 << " GeV";
                              std::cout << ", M(NuR2) = " << hnuQuadjet.mNuR2 << " GeV" << std::endl;
                              std::cout << "\tM(mumu) = " << hnuQuadjet.mLL << " GeV";
                              std::cout << ", M(JJ) = " << hnuQuadjet.mJJ << " GeV" << std::endl;
                              std::cout << "\tJets:   j1 ";
                              std::cout << "pt=" << hnuQuadjet.j1.pt() << " GeV, eta=" << hnuQuadjet.j1.eta() << ", phi=" << hnuQuadjet.j1.phi();
                              std::cout << ", j2 ";
                              std::cout << "pt=" << hnuQuadjet.j2.pt() << " GeV, eta=" << hnuQuadjet.j2.eta() << ", phi=" << hnuQuadjet.j2.phi();
                              std::cout << std::endl;
                              std::cout << "\tMuons: mu1 ";
                              std::cout << "pt=" << hnuQuadjet.mu1.pt() << " GeV, eta=" << hnuQuadjet.mu1.eta() << ", phi=" << hnuQuadjet.mu1.phi();
                              std::cout << ", mu2 ";
                              std::cout << "pt=" << hnuQuadjet.mu2.pt() << " GeV, eta=" << hnuQuadjet.mu2.eta() << ", phi=" << hnuQuadjet.mu2.phi();
                              std::cout << std::endl;
                          }
                          
                          if ( hnuQuadjet.mu1.pt() > cuts.minimum_mu1_pt ) {
                              hists.Mu1HighPtCut.fill( hnuQuadjet,mu1scaleFactor,mu2scaleFactor ) ;
                              if ( hnuQuadjet.mLL > cuts.minimum_mumu_mass ) { // dimuon mass rqmt
                                  hists.diLmassCut.fill( hnuQuadjet,mu1scaleFactor,mu2scaleFactor ) ; 
                              }
                          }
                      }
                  }
              }
          }
      }
  }

  // Closure test using two muons, not one
  if ( dirtyMuons.size() >= 2 ) { 
      if ( jetCands.size() >= 1 ) {
          HeavyNuEvent mmjet(HeavyNuEvent::QCD) ;
          initializeHNE(mmjet,pPU,pvHandle,!iEvent.isRealData(),isPFJets_,genweight) ;
          selectMuonPairs( dirtyMuons,mmjet ) ; 
//           mmjet.mu1 = dirtyMuons.at(0) ; 
//           mmjet.mu2 = dirtyMuons.at(1) ;
//           mmjet.nMuons = dirtyMuons.size() ; 
          if ( mmjet.isMC ) 
              mmjet.eventWgt *= muid_->weightForMC( mmjet.mu1.pt(),0 ) * muid_->weightForMC( mmjet.mu2.pt(),0 ) ; 
          // selectMuonsInJets( muCands,jetCands,mmjet ) ; 
          // if ( mmjet.nMuons >= 2 ) {  
          if ( mmjet.nMuons >= 2 && selectJets( jetCands,mmjet,1 ) ) {
              if ( mmjet.nJets < 2 ) mmjet.j2 = mmjet.j1 ; 
              if ( hnu::jetID(mmjet.j1) > 0 ) { 
                  mmjet.tjV1 = hnu::caloJetVertex(mmjet.j1, *jptJets);
                  
                  mmjet.scaleMuE() ; 
                  mmjet.regularize();
                  mmjet.calculate() ;
                  
                  double mu1scaleFactor, mu2scaleFactor ; 
                  mu1scaleFactor = qcdScaleFactor(mmjet.mu1.pt(),
                                                  rwLowPtbin,rwHighPtbin,rwQCD) ; 
                  mu2scaleFactor = qcdScaleFactor(mmjet.mu2.pt(),
                                                  rwLowPtbin,rwHighPtbin,rwQCD) ; 
                  
                  bool mu1trig = false ; bool mu2trig = false ;
                  if ( trig_->matchingEnabled() && iEvent.isRealData() ) {
                      mu1trig = trig_->isTriggerMatched( mmjet.mu1, iEvent) ; 
                      mu2trig = trig_->isTriggerMatched( mmjet.mu2, iEvent) ; 
                  } else if ( !iEvent.isRealData() ) {
                      mu1trig = trig_->simulateForMC( mmjet.mu1.pt(),mmjet.mu1.eta(),0 );
                      mu2trig = trig_->simulateForMC( mmjet.mu2.pt(),mmjet.mu2.eta(),0 );
                  }
                  if ( mu1trig || mu2trig ) {
                      //--- Impose vertex requirement here ---//
                      float deltaVzJ1M1 = fabs(mmjet.tjV1-mmjet.mu1.vertex().Z());
                      float deltaVzJ1M2 = fabs(mmjet.tjV1-mmjet.mu2.vertex().Z());
                      float deltaVzM1M2 = fabs(mmjet.mu1.vertex().Z()-mmjet.mu2.vertex().Z());
                      if ( (cuts.maxJetVZsepCM < 0 || cuts.maxVertexZsep < 0) || 
                           ((deltaVzJ1M1 < cuts.maxJetVZsepCM) && (deltaVzJ1M2 < cuts.maxJetVZsepCM) &&
                            (deltaVzM1M2 < cuts.maxVertexZsep)) ) {
                          if (iEvent.isRealData()) {
                              std::cout << "\tTWOMUINJET: " << iEvent.id() << std::endl;
                              std::cout << "\tM(mumu) = " << mmjet.mLL << " GeV";
                              std::cout << "\tJets:   j1 ";
                              std::cout << "pt=" << mmjet.j1.pt() << " GeV, eta=" << mmjet.j1.eta() << ", phi=" << mmjet.j1.phi();
                              std::cout << std::endl;
                              std::cout << "\tMET: " << MET << std::endl ; 
                              std::cout << "\tMuons: mu1 ";
                              std::cout << "pt=" << mmjet.mu1.pt() << " GeV, eta=" << mmjet.mu1.eta() << ", phi=" << mmjet.mu1.phi();
                              std::cout << ", mu2 ";
                              std::cout << "pt=" << mmjet.mu2.pt() << " GeV, eta=" << mmjet.mu2.eta() << ", phi=" << mmjet.mu2.phi();
                              std::cout << std::endl;
                          }
                          if ( mmjet.mu1.pt() > cuts.minimum_mu1_pt ) {
                              hists.mmJClosure.fill( mmjet,mu1scaleFactor,mu2scaleFactor ) ;
                              if ( MET < 20.0 ) hists.mmJClosureMETlt20.fill( mmjet,mu1scaleFactor,mu2scaleFactor ) ;
                          }
                      }
                  }
              }
          }
      }
  }
  if ( dirtyMuons.size() >= 1 ) {
      if ( jetCands.size() >= 1 ) {
          HeavyNuEvent mMjet(HeavyNuEvent::QCD) ;
          initializeHNE(mMjet,pPU,pvHandle,!iEvent.isRealData(),isPFJets_,genweight) ;
          mMjet.mu1 = dirtyMuons.at(0) ;
	  std::cout << "mu1 has pT " << mMjet.mu1.pt() << std::endl ; 
	  int dirtyMuonPosition = 1 ; 
          if ( mMjet.isMC ) 
              mMjet.eventWgt *= muid_->weightForMC( mMjet.mu1.pt(),0 ) ; 
          if ( selectJets( jetCands,mMjet,1 ) ) {
              if ( mMjet.nJets < 2 ) mMjet.j2 = mMjet.j1 ;
              int nStandardMuons = 0 ;
              double mcWeight = 1.0 ; 
              for (unsigned int i=0; i<muCands.size(); i++) {
                  pat::Muon isoMuon = muCands.at(i) ;
		  std::cout << "Investigating potential iso muon with pT " << isoMuon.pt() << std::endl ; 
                  bool overlapDirty = false ; 
                  for (unsigned int j=0; j<dirtyMuons.size(); j++) { // Same dR = 0.3 standard for dirty/clean muons
                      if ( deltaR(isoMuon.eta(), isoMuon.phi(), dirtyMuons.at(j).eta(), dirtyMuons.at(j).phi()) < 0.3 )
                          overlapDirty = true ;
                  }
                  if ( overlapDirty ) continue ;
                  if ( hnu::muIsolation(isoMuon,1.0) < cuts.muon_trackiso_limit ) {
		    std::cout << "Found an isolated muon with pT " << isoMuon.pt() << std::endl ; 
                      nStandardMuons++ ;
                      if ( applyMuIDCorrections_ )
                          mcWeight = muid_->weightForMC(isoMuon.pt(),0) ;
                      if ( isoMuon.pt() > mMjet.mu1.pt() ) {
                          mMjet.mu2 = mMjet.mu1 ;
                          mMjet.mu1 = isoMuon ;
			  dirtyMuonPosition = 2 ; 
                      } else {
                          mMjet.mu2 = isoMuon ;
                      }
		      std::cout << "Muons: mu1 with pT " << mMjet.mu1.pt() << " and mu2 with pT " << mMjet.mu2.pt() << std::endl ; 
                      break ;
                  }
              } // HERE!!!
              if ( nStandardMuons > 1 ) std::cout << "WARNING!  Found more than one standard muon" << std::endl ;
              mMjet.nMuons = dirtyMuons.size() + nStandardMuons ; 
              if ( applyMuIDCorrections_ ) mMjet.eventWgt *= mcWeight ;
              if ( nStandardMuons > 0 && hnu::jetID(mMjet.j1) > 0 ) { 
                  mMjet.tjV1 = hnu::caloJetVertex(mMjet.j1, *jptJets);
                
                  mMjet.scaleMuE() ; 
                  mMjet.regularize();
                  mMjet.calculate() ;
                  
                  double muScaleFactor ; 
                  muScaleFactor = qcdScaleFactor(((dirtyMuonPosition == 1)?(mMjet.mu1.pt()):(mMjet.mu2.pt())),
						 rwLowPtbin,rwHighPtbin,rwQCD) ; 
                
                  bool mu1trig = false ; bool mu2trig = false ;
                  if ( trig_->matchingEnabled() && iEvent.isRealData() ) {
                      mu1trig = trig_->isTriggerMatched( mMjet.mu1, iEvent) ; 
                      mu2trig = trig_->isTriggerMatched( mMjet.mu2, iEvent) ; 
                  } else if ( !iEvent.isRealData() ) {
                      mu1trig = trig_->simulateForMC( mMjet.mu1.pt(),mMjet.mu1.eta(),0 );
                      mu2trig = trig_->simulateForMC( mMjet.mu2.pt(),mMjet.mu2.eta(),0 );
                  }
                  if ( mu1trig || mu2trig ) {
                    //--- Impose vertex requirement here ---//
                    float deltaVzJ1M1 = fabs(mMjet.tjV1-mMjet.mu1.vertex().Z());
                    float deltaVzJ1M2 = fabs(mMjet.tjV1-mMjet.mu2.vertex().Z());
                    float deltaVzM1M2 = fabs(mMjet.mu1.vertex().Z()-mMjet.mu2.vertex().Z());
                    if ( (cuts.maxJetVZsepCM < 0 || cuts.maxVertexZsep < 0) || 
                         ((deltaVzJ1M1 < cuts.maxJetVZsepCM) && (deltaVzJ1M2 < cuts.maxJetVZsepCM) &&
                          (deltaVzM1M2 < cuts.maxVertexZsep)) ) {

                        if (iEvent.isRealData()) {
                            std::cout << "\tONEMUINJET: " << iEvent.id() << std::endl;
                            std::cout << "\tM(mumu) = " << mMjet.mLL << " GeV";
                            std::cout << "\tJets:   j1 ";
                            std::cout << "pt=" << mMjet.j1.pt() << " GeV, eta=" << mMjet.j1.eta() << ", phi=" << mMjet.j1.phi();
                            std::cout << std::endl;
			    std::cout << "\tMET: " << MET << std::endl ; 
                            std::cout << "\tMuons: ";
                            std::cout << "pt=" << mMjet.mu1.pt() << " GeV, eta=" << mMjet.mu1.eta() << ", phi=" << mMjet.mu1.phi();
                            std::cout << ", pt=" << mMjet.mu2.pt() << " GeV, eta=" << mMjet.mu2.eta() << ", phi=" << mMjet.mu2.phi();
			    std::cout << " dirty muon: " << dirtyMuonPosition ; 
                            std::cout << std::endl;
                        }

                        if ( mMjet.mu1.pt() > cuts.minimum_mu1_pt ) { 
			  hists.mMuJClosure.fill( mMjet,muScaleFactor ) ; 
			  if ( MET < 20.0 ) hists.mMuJClosureMETlt20.fill( mMjet,muScaleFactor ) ; 
                        }
                    }
                  }
              }
          }
        }
      }

  // 
  // Closure test
  // 
  if ( doClosure_ ) { 
      if ( pMuons->size() >= 1 && pJets->size() >= 3
	   && pMET->at(0).pt() < 20. ) { 

	HeavyNuEvent hnuClosure3jet1muon(HeavyNuEvent::CLO) ;  
	initializeHNE(hnuClosure3jet1muon,pPU,pvHandle,!iEvent.isRealData(),isPFJets_,genweight);
      
	selectMuonsInJets( muCands,jetCands,hnuClosure3jet1muon ) ; 
	if ( hnuClosure3jet1muon.nMuons == 1 ) {  
	  if ( selectJets( jetCands,hnuClosure3jet1muon ) ) {
            if ( hnu::jetID(hnuClosure3jet1muon.j1) > 0 && hnu::jetID(hnuClosure3jet1muon.j2) > 0 ) { 

              hnuClosure3jet1muon.tjV1 = hnu::caloJetVertex(hnuClosure3jet1muon.j1, *jptJets);
              hnuClosure3jet1muon.tjV2 = hnu::caloJetVertex(hnuClosure3jet1muon.j2, *jptJets);

              hnuClosure3jet1muon.regularize(); 
              hnuClosure3jet1muon.scaleMuE();
              hnuClosure3jet1muon.calculate() ; 

              bool mu1trig = false ; 
              if ( trig_->matchingEnabled() && iEvent.isRealData() ) {
	        mu1trig = trig_->isTriggerMatched( hnuClosure3jet1muon.mu1, iEvent) ; 
              } else if (!iEvent.isRealData()) {
	        mu1trig = trig_->simulateForMC( hnuClosure3jet1muon.mu1.pt(),
                                                hnuClosure3jet1muon.mu1.eta(),0 );
              }

              if ( mu1trig ) { 
	        double mu1scaleFactor = qcdScaleFactor(hnuClosure3jet1muon.mu1.pt(),
                                                       rwLowPtbin,rwHighPtbin,rwQCD) ; 

                //--- Impose vertex requirement here ---//
                float deltaVzJ1J2 = fabs(hnuClosure3jet1muon.tjV1-hnuClosure3jet1muon.tjV2);
                float deltaVzJ1M1 = fabs(hnuClosure3jet1muon.tjV1-hnuClosure3jet1muon.mu1.vertex().Z());
                float deltaVzJ2M1 = fabs(hnuClosure3jet1muon.tjV2-hnuClosure3jet1muon.mu1.vertex().Z());
                if ( (cuts.maxJetVZsepCM < 0) ||
                     ((deltaVzJ1J2 < cuts.maxJetVZsepCM) && (deltaVzJ1M1 < cuts.maxJetVZsepCM) &&
                      (deltaVzJ2M1 < cuts.maxJetVZsepCM)) ) { 
		  hists.LJJJClosure.fill( hnuClosure3jet1muon, mu1scaleFactor ) ; 
                  if(iEvent.isRealData()) {
                    std::cout << "\t" << iEvent.id() << std::endl;
                    std::cout << "3 jet: M(NuR1) = " << hnuClosure3jet1muon.mNuR1 << " GeV" << std::endl;
                    std::cout << "3 jet: M(JJ)   = " << hnuClosure3jet1muon.mJJ << " GeV" << std::endl;
                    std::cout << "\tJets:   j1 ";
                    std::cout << "pt=" << hnuClosure3jet1muon.j1.pt() << " GeV, eta=" << hnuClosure3jet1muon.j1.eta() << ", phi=" << hnuClosure3jet1muon.j1.phi();
                    std::cout << ", j2 ";
                    std::cout << "pt=" << hnuClosure3jet1muon.j2.pt() << " GeV, eta=" << hnuClosure3jet1muon.j2.eta() << ", phi=" << hnuClosure3jet1muon.j2.phi();
                    std::cout << std::endl;
                    std::cout << "\tMuons: mu1 ";
                    std::cout << "pt=" << hnuClosure3jet1muon.mu1.pt() << " GeV, eta=" << hnuClosure3jet1muon.mu1.eta() << ", phi=" << hnuClosure3jet1muon.mu1.phi();
                    std::cout << std::endl;
                  }
                }
              }
	    }
	  }
	}
      }
    
      if ( muCands.size() >= 1 && jetCands.size() >= 2 
	   && pMET->at(0).pt() < 20. ) { // MET cut to remove W+2 jets

	HeavyNuEvent hnuClosure2jet1muon(HeavyNuEvent::CLO) ;
	initializeHNE(hnuClosure2jet1muon,pPU,pvHandle,!iEvent.isRealData(),isPFJets_,genweight);

        for (unsigned int i=0; i<jetCands.size(); i++) {
          pat::Jet jet = jetCands.at(i).first ;
          if ( jet.pt() > cuts.minimum_jet_pt ) {
            hnuClosure2jet1muon.nJets++ ;
            if ( hnuClosure2jet1muon.nJets == 1 ) { 
              hnuClosure2jet1muon.j1 = jet ; 
              hnuClosure2jet1muon.j1scale = 1.0 ;
              hnuClosure2jet1muon.tjV1 = hnu::caloJetVertex(hnuClosure2jet1muon.j1, *jptJets);
            } else if ( hnuClosure2jet1muon.nJets == 2 ) { 
              hnuClosure2jet1muon.j2 = jet ; 
              hnuClosure2jet1muon.j2scale = 1.0 ;
              hnuClosure2jet1muon.tjV2 = hnu::caloJetVertex(hnuClosure2jet1muon.j2, *jptJets);
            } else break ; 
          }
        }

        if ( hnuClosure2jet1muon.nJets == 2 &&
             hnu::jetID(hnuClosure2jet1muon.j1) > 0 && hnu::jetID(hnuClosure2jet1muon.j2) > 0 ) { 
          for (unsigned int i=0; i<muCands.size(); i++) { 
	    if ( hnuClosure2jet1muon.nMuons == 1 ) break ; 
            pat::Muon iM = muCands.at(i) ; 
            if ( hnu::muIsolation(iM,1.0) < cuts.muon_trackiso_limit ) {
	      double dRj1 = deltaR(iM.eta(), iM.phi(), hnuClosure2jet1muon.j1.eta(), hnuClosure2jet1muon.j1.phi()) ; 
              double dRj2 = deltaR(iM.eta(), iM.phi(), hnuClosure2jet1muon.j2.eta(), hnuClosure2jet1muon.j2.phi()) ; 
              if (dRj1 > cuts.minimum_muon_jet_dR && dRj2 > cuts.minimum_muon_jet_dR) { 
	        hnuClosure2jet1muon.nMuons++ ; 
                if ( hnuClosure2jet1muon.nMuons == 1 ) {
                  hnuClosure2jet1muon.mu1 = iM ;
                  if ( applyMuIDCorrections_ )
                      hnuClosure2jet1muon.eventWgt *= muid_->weightForMC((hnuClosure2jet1muon.mu1.pt()),0) ;
                }
                else std::cout << "WARNING: Expected empty muon position" << std::endl ; 
              }
            }
          }
          if ( hnuClosure2jet1muon.nMuons > 0 ) {

            hnuClosure2jet1muon.regularize(); 
            hnuClosure2jet1muon.scaleMuE();
            hnuClosure2jet1muon.calculate() ; 

            bool mu1trig = false ; 
            if ( trig_->matchingEnabled() && iEvent.isRealData() ) {
              mu1trig = trig_->isTriggerMatched( hnuClosure2jet1muon.mu1, iEvent) ; 
            } else if (!iEvent.isRealData()) {
	      mu1trig = trig_->simulateForMC( hnuClosure2jet1muon.mu1.pt(),hnuClosure2jet1muon.mu1.eta(),0 ) ; 
            }

            if ( mu1trig ) { 
	      //--- Impose vertex requirement here ---//
              float deltaVzJ1J2 = fabs(hnuClosure2jet1muon.tjV1-hnuClosure2jet1muon.tjV2);
              float deltaVzJ1M1 = fabs(hnuClosure2jet1muon.tjV1-hnuClosure2jet1muon.mu1.vertex().Z());
              float deltaVzJ2M1 = fabs(hnuClosure2jet1muon.tjV2-hnuClosure2jet1muon.mu1.vertex().Z());

              if ( (cuts.maxJetVZsepCM < 0) ||
                   ((deltaVzJ1J2 < cuts.maxJetVZsepCM) && (deltaVzJ1M1 < cuts.maxJetVZsepCM) &&
                    (deltaVzJ2M1 < cuts.maxJetVZsepCM)) ) { 
                hists.L2JClosure.fill( hnuClosure2jet1muon ) ; 
                if(iEvent.isRealData()) {
                  std::cout << "\t" << iEvent.id() << std::endl;
                  std::cout << "2 jet: M(NuR1) = " << hnuClosure2jet1muon.mNuR1 << " GeV" << std::endl;
                  std::cout << "2 jet: M(JJ)   = " << hnuClosure2jet1muon.mJJ << " GeV" << std::endl;
                  std::cout << "\tJets:   j1 ";
                  std::cout << "pt=" << hnuClosure2jet1muon.j1.pt() << " GeV, eta=" << hnuClosure2jet1muon.j1.eta() << ", phi=" << hnuClosure2jet1muon.j1.phi();
                  std::cout << ", j2 ";
                  std::cout << "pt=" << hnuClosure2jet1muon.j2.pt() << " GeV, eta=" << hnuClosure2jet1muon.j2.eta() << ", phi=" << hnuClosure2jet1muon.j2.phi();
                  std::cout << std::endl;
                  std::cout << "\tMuons: mu1 ";
                  std::cout << "pt=" << hnuClosure2jet1muon.mu1.pt() << " GeV, eta=" << hnuClosure2jet1muon.mu1.eta() << ", phi=" << hnuClosure2jet1muon.mu1.phi();
                  std::cout << std::endl;
                }
              }
            }
	  }
	}
      }
  }
  */
  
  return keepThisEvent ; 
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuJetBackground::beginJob() {
  firstEvent_ = true;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuJetBackground::endJob() {
  // nnif_->endJob();
  trig_->endJob();
  muid_->endJob();
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuJetBackground);

