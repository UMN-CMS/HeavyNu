import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")


process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)


import FWCore.Modules.printContent_cfi
process.dumpEv = FWCore.Modules.printContent_cfi.printContent.clone()

# source
process.source = cms.Source("PoolSource",      
    fileNames=cms.untracked.vstring(     )
                  
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.options = cms.untracked.PSet(
                            fileMode = cms.untracked.string('NOMERGE')
)



## Load additional processes
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
## global tags:
process.GlobalTag.globaltag = cms.string('GR_R_36X_V11A::All') ## GF changed here
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')


################################################################################################
###    P r e p a r a t i o n      o f    t h e    P A T    O b j e c t s   f r o m    A O D  ###
################################################################################################

## pat sequences to be loaded:
#process.load("PhysicsTools.PFCandProducer.PF2PAT_cff")
process.load("PhysicsTools.PatAlgos.patSequences_cff")

#
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## MET creation     <=== WARNING: YOU MAY WANT TO MODIFY THIS PART OF THE CODE       %%%%%%%%%%%%%
##                                specify the names of the MET collections that you need here %%%%
##                                                                                             #%%
## if you don't specify anything the default MET is the raw Calo MET                           #%%
process.caloMET = process.patMETs.clone(                                                       #%%
        metSource = cms.InputTag("met","","RECO"),
            addTrigMatch = cms.bool(False),
            addMuonCorrections = cms.bool(False),
            addGenMET = cms.bool(False),
        )
process.tcMET = process.patMETs.clone(                                                         #%%
        metSource = cms.InputTag("tcMet","","RECO"),
            addTrigMatch = cms.bool(False),
            addMuonCorrections = cms.bool(False),
            addGenMET = cms.bool(False),
        )
process.pfMET = process.patMETs.clone(                                                         #%%
        metSource = cms.InputTag("pfMet","","RECO"),
            addTrigMatch = cms.bool(False),
            addMuonCorrections = cms.bool(False),
            addGenMET = cms.bool(False),
        )
## specify here what you want to have on the plots! <===== MET THAT YOU WANT ON THE PLOTS  %%%%%%%
myMetCollection   = 'caloMET'
myPfMetCollection =   'pfMET'
myTcMetCollection =   'tcMET'

## specify here what you want to have on the plots! <===== MET THAT YOU WANT ON THE PLOTS  %%%%%%%
## myDesiredMetCollection = 'layer1RawCaloMETs'
## modify the sequence of the MET creation:                                                    #%%
process.makePatMETs = cms.Sequence(process.caloMET*process.tcMET*process.pfMET)

 ## GF changed here: more than one pat done here
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## modify the final pat sequence: keep only electrons + METS (muons are needed for met corrections)
process.load("RecoEgamma.EgammaIsolationAlgos.egammaIsolationSequence_cff")
#process.patElectronIsolation = cms.Sequence(process.egammaIsolationSequence)

process.patElectrons.isoDeposits = cms.PSet()
process.patElectrons.userIsolation = cms.PSet()
process.patElectrons.addElectronID = cms.bool(True)
process.patElectrons.electronIDSources = cms.PSet(
    simpleEleId95relIso= cms.InputTag("simpleEleId95relIso"),
    simpleEleId90relIso= cms.InputTag("simpleEleId90relIso"),
    simpleEleId85relIso= cms.InputTag("simpleEleId85relIso"),
    simpleEleId80relIso= cms.InputTag("simpleEleId80relIso"),
    simpleEleId70relIso= cms.InputTag("simpleEleId70relIso"),
    simpleEleId60relIso= cms.InputTag("simpleEleId60relIso"),
    simpleEleId95cIso= cms.InputTag("simpleEleId95cIso"),
    simpleEleId90cIso= cms.InputTag("simpleEleId90cIso"),
    simpleEleId85cIso= cms.InputTag("simpleEleId85cIso"),
    simpleEleId80cIso= cms.InputTag("simpleEleId80cIso"),
    simpleEleId70cIso= cms.InputTag("simpleEleId70cIso"),
    simpleEleId60cIso= cms.InputTag("simpleEleId60cIso"),    
    )
##
process.patElectrons.addGenMatch = cms.bool(False)
process.patElectrons.embedGenMatch = cms.bool(False)
process.patElectrons.usePV = cms.bool(False)
##
process.load("ElectroWeakAnalysis.WENu.simpleEleIdSequence_cff")
 ## GF changed here: 
# you have to tell the ID that it is data
process.simpleEleId95relIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId90relIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId85relIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId80relIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId70relIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId60relIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId95cIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId90cIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId85cIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId80cIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId70cIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId60cIso.dataMagneticFieldSetUp = cms.bool(True)
#

process.patElectronIDs = cms.Sequence(process.simpleEleIdSequence)
process.makePatElectrons = cms.Sequence(process.patElectronIDs*process.patElectrons)
# process.makePatMuons may be needed depending on how you calculate the MET
process.makePatCandidates = cms.Sequence(process.makePatElectrons+process.makePatMETs)
process.patDefaultSequence = cms.Sequence(process.makePatCandidates)
##
##  ################################################################################
##
##  the filter to select the candidates from the data samples: one electron
##
##
## WARNING: you may want to modify this item:
HLT_process_name = "HLT"   # 
# trigger path selection


# this path was introduced at run 140058
HLT_path_name     = "HLT_Ele15_SW_L1R"
# trigger filter name
HLT_filter_name  =  "hltL1NonIsoHLTNonIsoSingleElectronEt15PixelMatchFilter"

# keep this for runs <140058 when HLT_Ele15_SW_L1R was not present. HLT_Ele15_LW_L1R was prescaled starting from run 141956
HLT_path_name_extra0   = "HLT_Ele15_LW_L1R"
HLT_filter_name_extra0 = cms.untracked.InputTag("hltL1NonIsoHLTNonIsoSingleElectronLWEt15PixelMatchFilter","",HLT_process_name)

HLT_path_name_extra1   = "HLT_Ele15_SiStrip_L1R"
HLT_filter_name_extra1 = cms.untracked.InputTag("hltL1NonIsoHLTNonIsoSingleElectronSiStripEt15PixelMatchFilter","",HLT_process_name)


process.z1LegFilter =cms.EDFilter('WenuCandidateFilter',
                                  ### the input collections needed:
                                  electronCollectionTag = cms.untracked.InputTag("patElectrons","","PAT"),
                                  metCollectionTag = cms.untracked.InputTag(myMetCollection,"","PAT"),
                                  pfMetCollectionTag = cms.untracked.InputTag(myPfMetCollection,"","PAT"),
                                  tcMetCollectionTag = cms.untracked.InputTag(myTcMetCollection,"","PAT"),
                                  triggerCollectionTag = cms.untracked.InputTag("TriggerResults","",HLT_process_name),
                                  triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLT_process_name),
                                  hltpath = cms.untracked.string(HLT_path_name),
                                  hltpathFilter = cms.untracked.InputTag(HLT_filter_name,"",HLT_process_name),
                                  ebRecHits = cms.untracked.InputTag("reducedEcalRecHitsEB"),
                                  eeRecHits = cms.untracked.InputTag("reducedEcalRecHitsEE"),
                                  PrimaryVerticesCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
                                  ### here the preselection is applied
                                  # fiducial cuts:
                                  BarrelMaxEta = cms.untracked.double(1.4442),
                                  EndCapMinEta = cms.untracked.double(1.566),
                                  EndCapMaxEta = cms.untracked.double(2.5),
                                  # demand ecal driven electron:
                                  useEcalDrivenElectrons = cms.untracked.bool(True),
                                  # demand offline spike cleaning with the Swiss Cross criterion:
                                  useSpikeRejection = cms.untracked.bool(True),
                                  spikeCleaningSwissCrossCut = cms.untracked.double(0.95),
                                  # demand geometrically matched to an HLT object with ET>15GeV
                                  useTriggerInfo = cms.untracked.bool(False), # GF >>>>>>>>>>>> changed here on Thu Sep  9 09:25:44 CDT 2010
                                  electronMatched2HLT = cms.untracked.bool(False), # GF changed Thu Sep  9 09:25:44 CDT 2010
                                  electronMatched2HLT_DR = cms.untracked.double(0.2), # GF 
                                  useHLTObjectETCut = cms.untracked.bool(False), # GF 
                                  hltObjectETCut = cms.untracked.double(15.),
                                  useExtraTrigger = cms.untracked.bool(False), # GF  >>>>>>>>>> Thu Sep  9 09:25:44 CDT 2010
                                  vHltpathExtra = cms.untracked.vstring(HLT_path_name_extra0,HLT_path_name_extra1),
                                  vHltpathFilterExtra = cms.untracked.VInputTag(HLT_filter_name_extra0, HLT_filter_name_extra1),
                                  # ET Cut in the SC
                                  ETCut = cms.untracked.double(18.), # GF : changed here Thu Sep  9 09:25:44 CDT 2010
                                  METCut = cms.untracked.double(0.),
                                  # reject events with a 2nd electron with ET > 20 that passes the WP95%
                                  vetoSecondElectronEvents = cms.untracked.bool(False),
                                  storeSecondElectron = cms.untracked.bool(False),
                                  ETCut2ndEle = cms.untracked.double(25.),
                                  vetoSecondElectronIDType = cms.untracked.string("simpleEleId80relIso"), # GF
                                  vetoSecondElectronIDSign = cms.untracked.string("="),
                                  vetoSecondElectronIDValue = cms.untracked.double(7.),
                                  # Other parameters of the code - leave them as they are
                                  useValidFirstPXBHit = cms.untracked.bool(False),
                                  useConversionRejection = cms.untracked.bool(False),
                                  useExpectedMissingHits = cms.untracked.bool(False),
                                  maxNumberOfExpectedMissingHits = cms.untracked.int32(1),
                                  # calculate some new cuts
                                  calculateValidFirstPXBHit = cms.untracked.bool(True),
                                  calculateConversionRejection = cms.untracked.bool(True),
                                  calculateExpectedMissingHits = cms.untracked.bool(True),
                                  # we are dealing with DATA
                                  dataMagneticFieldSetUp = cms.untracked.bool(True),
                                  dcsTag = cms.untracked.InputTag("scalersRawToDigi"),
                                  )


##
##  ################################################################################
##
##  the filter to select the candidates from the data samples: one TIGHT electron
##
##
HLT_path_name_tight   = "HLT_Ele25_SW_L1R"
HLT_filter_name_tight = "hltL1NonIsoHLTNonIsoSingleElectronEt25PixelMatchFilter"

process.z1TightLegFilter =cms.EDFilter('WenuCandidateFilter',
                                       ### the input collections needed:
                                       electronCollectionTag = cms.untracked.InputTag("patElectrons","","PAT"),
                                       metCollectionTag = cms.untracked.InputTag(myMetCollection,"","PAT"),
                                       pfMetCollectionTag = cms.untracked.InputTag(myPfMetCollection,"","PAT"),
                                       tcMetCollectionTag = cms.untracked.InputTag(myTcMetCollection,"","PAT"),
                                       triggerCollectionTag = cms.untracked.InputTag("TriggerResults","",HLT_process_name),
                                       triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLT_process_name),
                                       hltpath = cms.untracked.string(HLT_path_name_tight),
                                       hltpathFilter = cms.untracked.InputTag(HLT_filter_name_tight,"",HLT_process_name),
                                       ebRecHits = cms.untracked.InputTag("reducedEcalRecHitsEB"),
                                       eeRecHits = cms.untracked.InputTag("reducedEcalRecHitsEE"),
                                       PrimaryVerticesCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
                                       ### here the preselection is applied
                                       # fiducial cuts:
                                       BarrelMaxEta = cms.untracked.double(1.4442),
                                       EndCapMinEta = cms.untracked.double(1.566),
                                       EndCapMaxEta = cms.untracked.double(2.5),
                                       # demand ecal driven electron:
                                       useEcalDrivenElectrons = cms.untracked.bool(True),
                                       # demand offline spike cleaning with the Swiss Cross criterion:
                                       useSpikeRejection = cms.untracked.bool(True),
                                       spikeCleaningSwissCrossCut = cms.untracked.double(0.95),
                                       # demand geometrically matched to an HLT object with ET>15GeV
                                       useTriggerInfo = cms.untracked.bool(False), # GF >>>>>>>>>>>> changed here on Thu Sep  9 09:25:44 CDT 2010
                                       electronMatched2HLT = cms.untracked.bool(False), # GF changed Thu Sep  9 09:25:44 CDT 2010
                                       electronMatched2HLT_DR = cms.untracked.double(0.2), # GF 
                                       useHLTObjectETCut = cms.untracked.bool(False), # GF 
                                       hltObjectETCut = cms.untracked.double(15.),
                                       useExtraTrigger = cms.untracked.bool(False), # GF  >>>>>>>>>> Wed Aug  4 11:52:45 CDT 2010
                                       vHltpathExtra = cms.untracked.vstring(HLT_path_name_extra0,HLT_path_name_extra1),
                                       vHltpathFilterExtra = cms.untracked.VInputTag(HLT_filter_name_extra0, HLT_filter_name_extra1),
                                       # ET Cut in the SC
                                       ETCut = cms.untracked.double(25.),
                                       METCut = cms.untracked.double(0.),
                                       # reject events with a 2nd electron with ET > 20 that passes the WP95%
                                       vetoSecondElectronEvents = cms.untracked.bool(False),
                                       storeSecondElectron = cms.untracked.bool(False),
                                       ETCut2ndEle = cms.untracked.double(20.),
                                       vetoSecondElectronIDType = cms.untracked.string("simpleEleId80relIso"), # GF
                                       vetoSecondElectronIDSign = cms.untracked.string("="),
                                       vetoSecondElectronIDValue = cms.untracked.double(7.),
                                       # Other parameters of the code - leave them as they are
                                       useValidFirstPXBHit = cms.untracked.bool(False),
                                       useConversionRejection = cms.untracked.bool(False),
                                       useExpectedMissingHits = cms.untracked.bool(False),
                                       maxNumberOfExpectedMissingHits = cms.untracked.int32(1),
                                       # calculate some new cuts
                                       calculateValidFirstPXBHit = cms.untracked.bool(True),
                                       calculateConversionRejection = cms.untracked.bool(True),
                                       calculateExpectedMissingHits = cms.untracked.bool(True),
                                       # we are dealing with DATA
                                       dataMagneticFieldSetUp = cms.untracked.bool(True),
                                       dcsTag = cms.untracked.InputTag("scalersRawToDigi"),
                                       )

####################################################################################
# turn superclusters into candidates and require that there's a second "Z loose leg", in HF or in ECAL"


# turn superclusters into candidates
process.hybridCandidateSuperClusters = cms.EDProducer("ConcreteEcalCandidateProducer",
                                                      src = cms.InputTag("correctedHybridSuperClusters"),
                                                      particleType = cms.string("gamma")
                                                      )
process.multi5x5CandidateSuperClusters = cms.EDProducer("ConcreteEcalCandidateProducer",
                                                        src = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"),
                                                        particleType = cms.string("gamma")
                                                        )
# select separately EB and EE candidates
process.ebSC = cms.EDFilter("CandViewSelector",
                            src = cms.InputTag("hybridCandidateSuperClusters"),
                            cut = cms.string('abs( eta ) < 1.4442')
                            )
process.eeSC = cms.EDFilter("CandViewSelector",
                            src = cms.InputTag("multi5x5CandidateSuperClusters"),
                            cut = cms.string('abs( eta ) > 1.560 & abs( eta ) < 3')
                            )
# join them in ECAL calidates (no-duplication guaranteed)
process.ecalSC = cms.EDFilter("CandViewMerger",
                              src = cms.VInputTag(cms.InputTag("ebSC"), cms.InputTag("eeSC"))
                              )                                                  
# select ECAL sucperclusters above a certain Et
process.ecalCut = cms.EDFilter("EtaPtMinCandViewSelector",
                               src = cms.InputTag("ecalSC"),
                               etaMin = cms.double(-3),
                               etaMax = cms.double(3),
                               ptMin = cms.double(15.0)
                               )
# this is the FILTER which requires at least N=2 ecal SC above a certain Et 
process.ecalSel = cms.EDFilter("CandViewCountFilter",
                               src = cms.InputTag("ecalCut"),
                               minNumber = cms.uint32(2)
                               )



# turn HF superclusters into candidates
process.hfCandidateSuperClusters = cms.EDProducer("ConcreteEcalCandidateProducer",
                                                  src = cms.InputTag("hfEMClusters"),
                                                  particleType = cms.string("gamma")
                                                  )
process.hfSC = cms.EDFilter("CandViewSelector",
                            src = cms.InputTag("hfCandidateSuperClusters"),
                            cut = cms.string('abs( eta ) > 3 & abs( eta ) < 5')
                            )
# select HF sucperclusters above a certain Et
process.hfCut = cms.EDFilter("EtaPtMinCandViewSelector",
                             src = cms.InputTag("hfSC"),
                             etaMin = cms.double(-5),
                             etaMax = cms.double(5),
                             ptMin = cms.double(13.0)
                             )
# this is the FILTER which requires at least N=1 HF SC above a certain Et 
process.hfSel = cms.EDFilter("CandViewCountFilter",
                               src = cms.InputTag("hfCut"),
                               minNumber = cms.uint32(1)
                               )

####################################################################################



process.z2EcalLegsPath = cms.Path(process.patDefaultSequence
                                  * process.hybridCandidateSuperClusters
                                  * process.multi5x5CandidateSuperClusters
                                  * process.ebSC
                                  * process.eeSC
                                  * process.ecalSC
                                  * process.ecalCut
                                  * process.ecalSel      # this requires 2 ECAL SC > 13 GeV
                                  * process.z1LegFilter  # this requires one electron > 15 with HLT matching
                                  # * process.dumpEv
                                  )

####################################################################################
process.z1Ecal1HfLegPath = cms.Path(process.patDefaultSequence
                                    * process.hfCandidateSuperClusters
                                    * process.hfSC
                                    * process.hfCut
                                    * process.hfSel        # this requires 1 HF SC > 10 GeV
                                    * process.z1LegFilter  # this requires one electron > 15 with HLT matching
                                    #  * process.dumpEv
                                    )

####################################################################################
process.z1EcalLegTightPath = cms.Path(process.patDefaultSequence
                                   * process.z1TightLegFilter  # this requires one electron > 25 with HLT matching
                                   #  * process.dumpEv
                                   )


process.z1lOutputModule = cms.OutputModule( "PoolOutputModule",
                                            fileName = cms.untracked.string(""),
                                            maxSize = cms.untracked.int32(500000),
                                            SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('z2EcalLegsPath','z1EcalLegTightPath','z1Ecal1HfLegPath'),
                                                                              )
                                            )


process.outpath = cms.EndPath(process.z1lOutputModule)

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
