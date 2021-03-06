import FWCore.ParameterSet.Config as cms

process = cms.Process("ElectronsSKIMHeavyNu")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

from PhysicsTools.PatAlgos.tools.metTools import *
from PhysicsTools.PatAlgos.tools.tauTools import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.pfTools import *

from PhysicsTools.PatAlgos.selectionLayer1.leptonCountFilter_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.photonCountFilter_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.electronCountFilter_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.jetCountFilter_cfi import *

from HLTrigger.HLTfilters.triggerResultsFilter_cfi import *

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# source
process.source = cms.Source("PoolSource",      
                            fileNames=cms.untracked.vstring('file:/data/franzoni/data/data-Run2012A-DoubleElectron/Run2012A-DoubleElectron-AOD-PromptReco-v1-000-191-833-5ADFF104-9C8C-E111-A37B-003048D3C982.root')
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Load additional processes
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')

#process.GlobalTag.globaltag = 'GR_R_50_V13::All'
#process.GlobalTag.globaltag = 'FT_53_V6_AN1::All'
process.GlobalTag.globaltag = 'GR_R_53_V11::All'


### Output needs to be created before working with PAT objects ###
process.out = cms.OutputModule( "PoolOutputModule",
   fileName = cms.untracked.string("poolout.1.root"),
   maxSize = cms.untracked.int32(3000000),
   SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring('pA'),
   ),
   outputCommands = cms.untracked.vstring("keep *")
)



###########################################################################
###  P r e p a r a t i o n      o f    t h e    P A T    O b j e c t s  ###
###########################################################################
#------------------
#Load PAT sequences
process.load("PhysicsTools.PatAlgos.patSequences_cff")



# cleaning: scraping filter
process.scrapingFilter = cms.EDFilter("FilterOutScraping",
                                      applyfilter = cms.untracked.bool(True),
                                      debugOn = cms.untracked.bool(False),
                                      numtrack = cms.untracked.uint32(10),
                                      thresh = cms.untracked.double(0.25)
                                      )


# Get a list of good primary vertices, in 42x, these are DAF vertices
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(3.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
    )

#
# this is to check event content
process.dumpEvContent = cms.EDAnalyzer("EventContentAnalyzer")


#
# pack all pat-related modules here
process.patCandidateSummary.candidates = cms.VInputTag( cms.InputTag("patMuons") , cms.InputTag("patElectrons"))
process.patCandidates      = cms.Sequence(
    #process.makePatMuons +
    process.makePatElectrons +
    #process.makePatJets +
    process.patCandidateSummary
    )
process.patDefaultSequence = cms.Sequence( process.patCandidates )

process.HLTTagAndProbe = cms.EDFilter("TriggerResultsFilter",
                                      hltResults              = cms.InputTag('TriggerResults','','HLT'),   # HLT results   - set to empty to ignore HLT
                                      l1tResults              = cms.InputTag(''),                 # L1 GT results - set to empty to ignore L1
                                      l1tIgnoreMask           = cms.bool(False),                  # use L1 mask
                                      l1techIgnorePrescales   = cms.bool(False),                  # read L1 technical bits from PSB#9, bypassing the prescales
                                      daqPartitions           = cms.uint32(0x01),                 # used by the definition of the L1 mask
                                      throw                   = cms.bool(False),                  # if HLT path not in the table, crash/ignore according to true/false
                                      triggerConditions       = cms.vstring(
                                                                             'HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v3',
                                                                             'HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v4',
                                                                             'HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v5'
                                                                             )
                                      )


# this does not filter event, only reduces the collection
process.ElectronsAbove35 = cms.EDFilter("CandViewSelector",
                                        src = cms.InputTag("patElectrons"),
                                        cut = cms.string('(ecalEnergy*sin(superClusterPosition.theta)) >' + str(35) )
                                        )
#here's the real filter
process.oneElectronAbove35 = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("ElectronsAbove35"),
                                       minNumber = cms.uint32(1)
                                       )


# class to filter on a combination of N electrons and M (distinct) superclusters - mass cut available as well  
import HeavyNu.MuFilter.EleFilter_cfi
process.elePLusSCFilter35     = HeavyNu.MuFilter.EleFilter_cfi.eleFilter.clone()
process.elePLusSCFilter35.minElePt  = cms.double( 35. )
process.elePLusSCFilter35.minSCEt   = cms.double( 35. )




process.pA = cms.Path( 

    process.scrapingFilter *
    process.goodOfflinePrimaryVertices   *

    process.HLTTagAndProbe *

    process.patDefaultSequence *
    
    # this does not filter event, only reduces the collection
    process.ElectronsAbove35  *
    # filter on presence of 1electron>35
    process.oneElectronAbove35  * 

    # require presence of an electron>35 and a separate SC>35; no mass requirements, at this point 
    process.elePLusSCFilter35  

    # * process.dumpEvContent
    
    )

process.outpath = cms.EndPath(process.out)

from PhysicsTools.PatAlgos.tools.coreTools import *
runOnData(process, ['All'])
