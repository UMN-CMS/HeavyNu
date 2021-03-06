#include "HeavyNu/AnalysisModules/src/HeavyNuTrigger.h"

// #include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "HeavyNuCommon.h"

//======================================================================

HeavyNuTrigger::HeavyNuTrigger(const edm::ParameterSet & iConfig) :
  trigEventTag_     ( iConfig.getParameter< edm::InputTag > ( "trigEventTag" ) ),
  muonTriggers_     ( iConfig.getParameter< std::vector<std::string> > ( "muonTriggers" ) ),
  electronTriggers_ ( iConfig.getParameter< std::vector<std::string> > ( "electronTriggers" ) ),
  beginRun_         ( iConfig.getParameter< std::vector<int> > ( "firstRun" ) ),
  endRun_           ( iConfig.getParameter< std::vector<int> > ( "lastRun" ) ),
  // muonMatch_        ( iConfig.getParameter< std::string >  ( "muonMatch"     ) ),
  electronFilters_  ( iConfig.getParameter< std::vector<std::string > >  ( "electronFilters" ) ),
  triggerPt_        ( iConfig.getParameter< double >       ( "triggerPt" ) ),
  muTriggerEta_     ( iConfig.getParameter< double >       ( "muTriggerEta" ) ),
  // trigEra_          ( iConfig.getParameter< int >          ( "trigEra"   ) ),
  johnnyApple_      ( iConfig.getParameter< int >          ( "randomSeed" ) )
{
  matchingEnabled_ = false;
  if ( trigEventTag_.label().size() &&
       ( muonTriggers_.size() || (electronTriggers_.size() && electronFilters_.size()) ) )
    matchingEnabled_ = true;

  if (!matchingEnabled_) {
    std::cout << "Trigger matching is === DISABLED ===" << std::endl;
    std::cout << "   (Trigger sim. random seed = "<<johnnyApple_<<")"<<std::endl;
    triggerRandom_ = new TRandom(johnnyApple_);
    // std::cout << "   (Trigger Era              = "<<trigEra_<<")"<<std::endl;
  } else {
    std::cout << "Trigger matching is === ENABLED ===" << std::endl;
  }
  std::cout << "   (Trigger pT               = " << triggerPt_    << ")" << std::endl;
  std::cout << "   (Trigger eta (muons)      = " << muTriggerEta_ << ")" << std::endl;

  // Safety: first/last run vectors are not the same size as the trigger list
  if ( muonTriggers_.size() != beginRun_.size() || 
       muonTriggers_.size() != endRun_.size() || 
       beginRun_.size() != endRun_.size() ) { 
    std::cout << "INFO: Size of trigger initialization vectors are not equal." << std::endl ; 
    std::cout << "      Resetting all vectors to size of trigger list and dropping run restriction " << std::endl ; 

    unsigned int maxsize = std::max( muonTriggers_.size(),
				     std::max( beginRun_.size(),endRun_.size() ) ) ; 
    
    beginRun_.clear() ; 
    endRun_.clear() ; 
    for (unsigned int i=0; i<maxsize; i++) { 
      beginRun_.push_back(0) ; 
      endRun_.push_back(999999) ; 
    }
  }
  
}                                      // HeavyNuTrigger::HeavyNuTrigger

//======================================================================

void
HeavyNuTrigger::book(const TFileDirectory& tdir, trigHistos_t *thist)
{
  if (!matchingEnabled_) return;

  assert(thist);

  thist->trigMatchPtCorrel = tdir.make< TH2D >( "h2d_trigMatchPtCorrel","", 60,0.,300.,60,0.,300. );
  thist->trigMatchPtCorrel->SetTitle ( "#mu_{reco} vs. #mu_{trig} p_{T} (GeV)" );
  thist->trigMatchPtCorrel->SetXTitle( "p_{T}(#mu_{reco}) (GeV)" );
  thist->trigMatchPtCorrel->SetYTitle( "p_{T}(#mu_{trig}) (GeV)" );

  thist->trigMatchDR2      = tdir.make< TH1D >( "h1d_trigMatchDR2", "",50,0.,0.25 );
  thist->trigMatchDR2      ->SetTitle ( "(#Delta R(#mu_{reco}, #mu_{trig}))^{2}" );
  thist->trigMatchDR2      ->SetXTitle( "(#Delta R(#mu_{reco}, #mu_{trig}))^{2}" );

  thist->trigMatchDRDPt    = tdir.make< TH2D >( "h2d_trigMatchDRDPt","",50,0.,0.25,50,-0.5,0.5 );
  thist->trigMatchDRDPt    ->SetTitle ( "(#Delta R(#mu_{reco}, #mu_{trig}))^{2} vs #Delta p_{T}_{rel}(#mu_{reco}, #mu_{trig})" );
  thist->trigMatchDRDPt    ->SetXTitle( "(#Delta R(#mu_{reco}, #mu_{trig}))^{2}" );
  thist->trigMatchDRDPt    ->SetYTitle( "#Delta p_{T}_{rel}(#mu_{reco}, #mu_{trig})" );

  thist->trigMatchDetaPhi  = tdir.make< TH2D >( "h2d_trigMatchDetaPhi","",50,-0.25,0.25,50,-0.25,0.25 );
  thist->trigMatchDetaPhi  ->SetTitle ( "#Delta #phi(#mu_{reco}, #mu_{trig}) vs #Delta #eta(#mu_{reco}, #mu_{trig})" );
  thist->trigMatchDetaPhi  ->SetXTitle( "#Delta #eta(#mu_{reco}, #mu_{trig})" );
  thist->trigMatchDetaPhi  ->SetYTitle( "#Delta #phi(#mu_{reco}, #mu_{trig})" );

  thist->trigUnmatchedPt   = tdir.make< TH1D >( "h1d_trigUnmatchedPt","", 200,0,2000 );
  thist->trigUnmatchedPt   ->SetTitle( "Unmatched muon p_{T}; #mu p_{T} (GeV)" );

  thist->trigAllCandMuPt   = tdir.make< TH1D >( "h1d_trigAllCandMuPt","", 200,0,2000 );
  thist->trigAllCandMuPt   ->SetTitle( "All trigger-match muon candidates p_{T}; #mu p_{T} (GeV)" );

  thist->trigUnmatchedEtaPhi = tdir.make< TH2D >( "h2d_trigUnmatchEtaPhi","", 50,-2.5,2.5,63,0.,6.3 );
  thist->trigUnmatchedEtaPhi ->SetTitle( "Unmatched muon #eta/#phi; #mu #eta; #mu #phi" );

  thist->trigAllCandMuEtaPhi = tdir.make< TH2D >( "h2d_trigAllCandMuEtaPhi","", 50,-2.5,2.5,63,0.,6.3 );
  thist->trigAllCandMuEtaPhi ->SetTitle( "All trigger-match muon candidates; #mu #eta; #mu #phi" );

}                                                // HeavyNuTrigger::book

//======================================================================

bool
HeavyNuTrigger::isTriggerMatched(const pat::Muon&  m,
                                 const edm::Event& iEvent,
                                 trigHistos_t *thist)
{
    bool matched = false;
    // bool passTrig=false;

    // Only one trigger can be used for matching in a given run
    int run = iEvent.run() ;
    std::vector<std::string> validHLTpaths ;

    for (unsigned int i = 0; i < muonTriggers_.size(); i++)
        if (run >= beginRun_.at(i) && run <= endRun_.at(i)) validHLTpaths.push_back(muonTriggers_.at(i)) ;

       //std::cout << "Looking for trigger match for muon with pT " << m.pt() 
       // 	    << " and eta " << m.eta() << std::endl ; 

    // muon trigger matching is only allowed within |eta|<2.1
    if ( matchingEnabled_ && (fabs(m.eta()) < muTriggerEta_) )
    {
        // PAT trigger information
        //     edm::Handle< pat::TriggerEvent > triggerEvent;

        //     iEvent.getByLabel( trigEventTag_, triggerEvent );
        //     if ( !triggerEvent.isValid() ) {
        //       std::cerr << "triggerEvent not found " << std::endl;
        //       return false;
        //     }

        const pat::TriggerObjectStandAloneCollection muonMatchCollection = m.triggerObjectMatches();
        //std::cout << "Trigger object matches size: " << muonMatchCollection.size() << std::endl ; 

        for (unsigned int i = 0; i < muonMatchCollection.size(); i++)
        {
            if ( matched ) break ; // Quit as soon as we find a match
            // std::cout << "Trigger object " << i+1 << " of " << muonMatchCollection.size() << std::endl ; 
            pat::TriggerObject muonTrigger = muonMatchCollection.at(i) ;
            pat::TriggerObjectStandAlone muonTriggerInfo = muonMatchCollection.at(i) ;
            // Look for a match with one of our paths
            std::vector<std::string> hltPaths = muonTriggerInfo.pathNames(true, false) ;
            bool hltPathMatch = false ;
            for (unsigned int j = 0; j < hltPaths.size(); j++)
            {
                if (hltPathMatch) break ;
                for (unsigned int k = 0; k < validHLTpaths.size(); k++)
                {
                    if (hltPaths.at(j) == validHLTpaths.at(k))
                    {
                        // std::cout << "Found a match to HLT path: " << muonTriggers_.at(k) << std::endl ; 
                        hltPathMatch = true ;
                        break ;
                    }
                }
            }
            // Finding a trigger object is not enough.  Need to impose the last filter (pT) 
            // Requirements to see if the trigger would have accepted the event based on this muon
            if ( hltPathMatch && muonTrigger.pt() > triggerPt_ )
            {
                double dr2  = reco::deltaR2 <pat::Muon, pat::TriggerObject > ( m, muonTrigger );
                //double dpt  = 1. - (muonTrigger.pt() / m.pt());
                //double dphi = reco::deltaPhi( m.phi(), muonTrigger.phi() );
                //double deta = m.eta() - muonTrigger.eta();

                // One more requirement: make sure that the muon is nearby the trigger object
                if ( sqrt(dr2) < 0.1 ) matched = true ;

                //         std::cout << "pT is " << muonTrigger.pt() << " with dpt = " << dpt << std::endl ; 
                //         std::cout << "dR is " << sqrt(dr2) << std::endl ; 

                /*if (thist)
                {
                    thist->trigMatchPtCorrel->Fill( m.pt(), muonTrigger.pt() );
                    thist->trigMatchDR2     ->Fill( dr2 );
                    thist->trigMatchDRDPt   ->Fill( dr2, dpt );
                    thist->trigMatchDetaPhi ->Fill( deta, dphi );
                }*/
            }
            /*if ( !matched )
            {
                if ( thist )
                {
                    thist->trigUnmatchedPt->Fill( m.pt() );
                    thist->trigUnmatchedEtaPhi->Fill( m.eta(), m.phi() );
                }
            }*/
        }

        /*if ( thist )
        {
            thist->trigAllCandMuPt->Fill( m.pt() );
            thist->trigAllCandMuEtaPhi->Fill( m.eta(), m.phi() );
        }*/
    }
    return ( matched );
}                                    // HeavyNuTrigger::isTriggerMatched

bool HeavyNuTrigger::isTriggerMatched(const pat::Electron& e1,
                                      const pat::Electron& e2,
                                      const edm::Event& iEvent,
                                      int nMatchesNeeded,
                                      trigHistos_t *thist)
{
    if ( !matchingEnabled_ ) return false ;

    // This code is adopted from HeavyNuEleTrigger.cc

    edm::InputTag hltTrigInfoTag("hltTriggerSummaryAOD","","HLT");
    edm::Handle<trigger::TriggerEvent> trigEvent;
    iEvent.getByLabel(hltTrigInfoTag, trigEvent);
    if ( !trigEvent.isValid() )
    {
        std::cout << "hltTrigInfoTag not found. Bailing out. " << std::endl;
        assert(0);
    }

    bool matched_e1 = false, matched_e2 = false;
    for(std::vector<std::string>::const_iterator ieFilters = electronFilters_.begin(); ieFilters != electronFilters_.end(); ++ieFilters)
    {
        std::vector<math::XYZTLorentzVector> trigObjects;
        trigtools::getP4sOfObsPassingFilter(trigObjects, *trigEvent, *ieFilters, hltTrigInfoTag.process());
        for(std::vector<math::XYZTLorentzVector>::const_iterator iTP = trigObjects.begin(); iTP != trigObjects.end(); ++iTP)
        {
            double dr2_e1 = reco::deltaR2<pat::Electron, math::XYZTLorentzVector > (e1, *iTP);
            double dr2_e2 = reco::deltaR2<pat::Electron, math::XYZTLorentzVector > (e2, *iTP);
            if ( sqrt(dr2_e1) < 0.1 && (dr2_e1 <= dr2_e2) )
            {
                // std::cout << "e1 Matched to this object" << std::endl ;
                matched_e1 = true ;
            }
            if ( sqrt(dr2_e2) < 0.1 && (dr2_e2 <= dr2_e1) )
            {
                // std::cout << "e2 Matched to this object" << std::endl ;
                matched_e2 = true ;
            }
        }
    }

    if(nMatchesNeeded <= 1) 
    {
        //std::cout << "One Trig Match: " << matched_e1 << "\t" << matched_e2 << std::endl;
        return matched_e1 || matched_e2;
    }
    else
    {
        //std::cout << "Two Trig Match: " << matched_e1 << "\t" << matched_e2 << std::endl;
        return matched_e1 && matched_e2;
    }

}                                    // HeavyNuTrigger::isTriggerMatched


//======================================================================

bool
HeavyNuTrigger::simulateForMC(double pt,double eta,int signOfError2apply)
{
  if (matchingEnabled_)
    throw cms::Exception("invalid trigger configuration");

  // Triggers outside |eta| < 2.1 not allowed
  if ( fabs(eta) >= 2.1 ) return false ; 
  // Cannot trigger if you do not meet the minimum pT threshold
  if ( pt < triggerPt_ ) return false ;

  // Run 2012 results use Mu40_eta2p1 only

  // For 2012, the trigger corrections are made as a function of eta
  //const double effslo2012[]   = {0.776442,0.844443,0.931687,0.941174,0.846101,0.807243};
  //const double effsnom2012[]  = {0.786102,0.856497,0.936184,0.945360,0.858314,0.816294};
  //const double effshi2012[]   = {0.795555,0.868014,0.940492,0.949355,0.869968,0.825120};
  //const double upedge2012[]   = {    -1.2,    -0.9,     0.0,     0.9,     1.2,     2.1};
  
  
  // Trigger Efficiencies for 2012A, B, and C
  // https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs#2012_data
  //const double effslo2012ABC[]   = {0.000000,0.000000,0.000000};  //UNCERTAINTIES +/- 0.2%
  // https://twiki.cern.ch/twiki/bin/view/CMS/MuonTagAndProbe
  //const double effsnom2012ABC[]  = {0.940100,0.843685,0.830823};
  //const double effshi2012ABC[]   = {0.000000,0.000000,0.000000};
  // Jan22rereco reference - Mu40 HighPTID
  // https://indico.cern.ch/getFile.py/access?contribId=2&resId=0&materialId=slides&confId=257000
  const double effsnom2012ABC[]  = {0.928,0.8302,0.8018};
  const double upedge2012ABC[]   = {  0.9,   1.2,   2.1};

  // 2011 A is the default
  const double *effs = effsnom2012ABC; 
  double scale = 1.0;
  if ( signOfError2apply )
  {
      if(pt < 45) // add extra uncertainty for possible trigger efficiency turn-on below 45 GeV
      {
          scale = (signOfError2apply > 0) ? (1 + 0.002) : (1 - 0.10);
      }
      else scale = (signOfError2apply > 0) ? (1 + 0.002) : (1 - 0.002);
  }

  int i;
  const double * upedge = upedge2012ABC ;
  for (i = 0; abs(eta) < 2.1 && upedge[i] < eta; i++);
  double eff=effs[i];
    
  return (triggerRandom_->Uniform()<(eff*scale));
}

//======================================================================

bool
HeavyNuTrigger::simulateForMCElePt(double pt,double eta,int signOfError2apply)
{
  // what follows is somewhat of a dummy method, which returns true basucally always
  // it's put in place in case we will want to switch over from using dielectron mass to single electron pt
  // to parameterize trigger efficiencies

  if (matchingEnabled_)
    throw cms::Exception("invalid trigger configuration");

  // Triggers outside |eta| < 2.5 not allowed
  if ( fabs(eta) >= 2.5 ) return false ; 
  // Cannot trigger if you do not meet the minimum pT threshold
  if ( pt < triggerPt_ ) return false ;

  // Trigger studies updated 28 May 2012;
  // so far we don't have measurements for 2011 => set them the same as 2012 and throw a warning 
  // if you want to make them sensible, modify these arrays BUT ALSO the code below

  //const double effElelo2011a[]  = {1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000};
  //const double effElenom2011a[] = {1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000};
  //const double effElehi2011a[]  = {1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000};
  //const double upedgeEle2011a[]  = {      50,      60,      80,     100,     200,    3500,      -1};
  //
  //const double effElelo2011b[]  = {1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000};
  //const double effElenom2011b[] = {1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000};
  //const double effElehi2011b[]  = {1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000};
  //const double upedgeEle2011b[]  = {      50,      60,      80,     100,     200,    3500,      -1};

  const double effElelo2012[]   = {1.,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000};
  const double effElenom2012[]  = {1.,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000};
  const double effElehi2012[]   = {1.,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000};
  const double upedgeEle2012[]   = {0,      50,      60,      80,     100,     200,    3500,      -1};

  // 2012 is the default
  // actually, at the moment we only have 2012, so use twose in all cases (there was a warning about trigEra_ )
  const double *         effs = effElenom2012 ; 
//   if (trigEra_ == 20111) effs = effElenom2012 ;
//   if (trigEra_ == 20112) effs = effElenom2012 ;
  if ( signOfError2apply ) {
//     if ( trigEra_ == 20111 ) effs = (signOfError2apply > 0) ? effElehi2012 : effElelo2012;
//     if ( trigEra_ == 20112 ) effs = (signOfError2apply > 0) ? effElehi2012 : effElelo2012;
    effs = (signOfError2apply > 0) ? effElehi2012 : effElelo2012;
  }

  int i;
  const double *         upedge = upedgeEle2012 ; 
//   if (trigEra_ == 20111) upedge = upedgeEle2012 ; 
//   if (trigEra_ == 20112) upedge = upedgeEle2012 ; 
  for (i=0; upedge[i]>0 && upedge[i]<pt; i++);
  double eff=effs[i];
    
  return (triggerRandom_->Uniform()<eff);
}

//======================================================================


bool
HeavyNuTrigger::simulateForMCdiEleMass(double m,double eta,int signOfError2apply)
{
  // what follows incorporates the trigger efficiency measurements
  // as shown here:  http://homepages.spa.umn.edu/~franzoni/reps/12/may25/ 

  if (matchingEnabled_)
    throw cms::Exception("invalid trigger configuration");

  // Triggers outside |eta| < 2.5 not allowed
  if ( fabs(eta) >= 2.5 ) return false ; 
  // eta may end up being used to separate EB,EE combinations?

  // Trigger studies updated 28 May 2012;
  // so far we don't have measurements for 2011 => set them the same as 2012 and throw a warning 
  // if you want to make them sensible, modify these arrays BUT ALSO the code below
  //const double effEleMasslo2011a[]  = {1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000};
  //const double effEleMassnom2011a[] = {1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000};
  //const double effEleMasshi2011a[]  = {1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000};
  //const double upedgeEleMass2011a[]  = {      50,      60,      80,     100,     200,    3500,      -1};
  //
  //const double effEleMasslo2011b[]  = {1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000};
  //const double effEleMassnom2011b[] = {1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000};
  //const double effEleMasshi2011b[]  = {1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000};
  //const double upedgeEleMass2011b[]  = {      50,      60,      80,     100,     200,    3500,      -1};
  
  //const double upedgeEleMass2012[]   = {70,80,90,100,110,120,130,140,160,180,200,300,-1};
  //const double effEleMassnom2012[]  = {1,1,0.979488,0.962445,0.918403,0.81,1,0.81,1,1,1,1,1};
  //const double effEleMasshi2012[]   = {1,1,0.992667,0.975271,0.984532,0.925816,1,0.959204,1,1,1,1,1};
  //const double effEleMasslo2012[]   = {0.383681,0.735776,0.95288,0.944351,0.74699,0.599341,0.753433,0.460385,0.69163,0.58867,0.261936,0.663539,0.471301};

  // central value is _two legs_ electron efficiency from the 2.4fb-1 0jet sample shown here:
  // http://homepages.spa.umn.edu/~franzoni/reps/12/jun11/trEffplots-superBins/efficiency-vs-mass-all-cases.png
  // errors of 0jet measurement added in quadrature to 0.94% taken as uncertainty on the approximation 0jetEff == 2jetEff
  //const double upedgeEleMass2012[]  = {0.     ,70      ,120,      200,      -1};
  //const double effEleMassnom2012[]  = {1       ,0.981803, 0.984283, 0.988131,  0.988131};
  //const double effEleMasshi2012[]   = {1.      ,0.991383, 0.997178, 1., 1.};
  //const double effEleMasslo2012[]   = {0.81523,0.972204,  0.969239, 0.955896, 0.955896};

  // update on Jun 21, with 2.9 fb-1 and V02-01-11
  //  ==> m.p.v. for single leg efficiency: hNuEtriggerEff is: 99.5345 %
  //  'transport' systematic between 0jets and >=2jets: 0.3%
  //  ==> arrays below are for DOUBLE leg efficiency: 
  //const double upedgeEleMass2012[]  = {0.,       70,        120,      200       ,-1};
  //const double effEleMassnom2012[]  = {0.965815, 0.990676,  0.993142, 0.987408  ,0.987408};
  //const double effEleMasshi2012[]   = {1.,       0.99389,   0.999,    1.        ,1.};
  //const double effEleMasslo2012[]   = {0.81685,  0.987434,  0.985218, 0.963095  ,0.963095};
  
  //we have decided to use flat 1 for the trigger for electrons and allow the scale normalization to deal with it
  const double upedgeEleMass2012[]  = {0.,       70,        120,      200       ,-1};
  const double effEleMassnom2012[]  = {1.00000, 1.00000, 1.00000, 1.00000, 1.00000};
  const double effEleMasshi2012[]   = {1.00000, 1.00000, 1.00000, 1.00000, 1.00000};
  const double effEleMasslo2012[]   = {1.00000, 1.00000, 1.00000, 1.00000, 1.00000};

  // 2012 is the default
  // actually, at the moment we only have 2012, so use twose in all cases (there was a warning about trigEra_ )
   const double *effs = effEleMassnom2012 ;
  if ( signOfError2apply ) {
    effs = (signOfError2apply > 0) ? effEleMasshi2012 : effEleMasslo2012;
  }

  int i;
  const double *upedge = upedgeEleMass2012 ; 
  for (i=0; upedge[i]>=0 && upedge[i]<m; i++);
  double eff=effs[i];

  return (triggerRandom_->Uniform()<eff);
}

//======================================================================


void HeavyNuTrigger::endJob()
{
}
