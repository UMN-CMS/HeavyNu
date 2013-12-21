// -*- C++ -*-
//
// Package:    HeavyNuPreFilter
// Class:      HeavyNuPreFilter
// 
/**\class HeavyNuPreFilter HeavyNuPreFilter.cc HeavyNu/HeavyNuPreFilter/src/HeavyNuPreFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
 */
//
// Original Author:  Nathaniel Pastika
//         Created:  Fri May 17 09:43:10 CDT 2013
// $Id: HeavyNuPreFilter.cc,v 1.1 2013/05/17 22:20:15 pastika Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"

#include "HeavyNu/AnalysisModules/src/HeavyNuCommon.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuEvent.h"

class HeavyNuPreFilter : public edm::EDFilter
{
public:
    explicit HeavyNuPreFilter(const edm::ParameterSet&);
    ~HeavyNuPreFilter();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    // ----------member data ---------------------------
    edm::InputTag muonTag_;
    edm::InputTag elecTag_;
    edm::InputTag elecRhoTag_;

    double elecRho_, ptL1_, ptL2_;
} ;

HeavyNuPreFilter::HeavyNuPreFilter(const edm::ParameterSet& iConfig)
{
    muonTag_ = iConfig.getParameter< edm::InputTag > ("muonTag");
    elecTag_ = iConfig.getParameter< edm::InputTag > ("electronTag");
    elecRhoTag_ = iConfig.getParameter< edm::InputTag > ("electronRho");
    ptL1_ = iConfig.getUntrackedParameter< double > ("ptL1");
    ptL2_ = iConfig.getUntrackedParameter< double > ("ptL2");
}

HeavyNuPreFilter::~HeavyNuPreFilter(){

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}

// ------------ method called on each new Event  ------------

bool
HeavyNuPreFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    edm::Handle<pat::MuonCollection> pMuons;
    iEvent.getByLabel(muonTag_, pMuons);

    edm::Handle<pat::ElectronCollection> pElecs;
    iEvent.getByLabel(elecTag_, pElecs);

    edm::Handle<reco::VertexCollection> pvHandle;
    iEvent.getByLabel("offlinePrimaryVertices", pvHandle);

    edm::Handle<double> electronRhoHandle ;
    iEvent.getByLabel(elecRhoTag_, electronRhoHandle) ;
    elecRho_ = ((electronRhoHandle.isValid()) ? (*(electronRhoHandle.product())) : 0.);

    // Look for valid muons
    std::vector<pat::Muon> muCands = hnu::getMuonList(pMuons, pvHandle, 20122, 40, 2.4, 1.0, 0, false);

    // Look for valid electrons
    std::vector< std::pair<pat::Electron, float> > eCands = hnu::getElectronList(pElecs, 2.5, 40, 40, 41+100, pvHandle, elecRho_);
    
    int nMu = muCands.size();
    int nE  = eCands.size();

    return (nMu >= 2) || (nE >= 2) || (nMu >= 1 && nE >= 1);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

void
HeavyNuPreFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(HeavyNuPreFilter);
