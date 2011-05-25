import FWCore.ParameterSet.Config as cms

import os

#import sys
#isMC=sys.modules['__main__'].isMC
#isMCsignal=sys.modules['__main__'].isMCsignal
#process = sys.modules['__main__'].process

isMC=False
isMCsignal=False
Training=False
isRun2010LoLumi=False
isRun2011=False

isData=not isMC

## Low and high lumi data selection is controlled by the JSON-derived cfi's imported
## below. For run 2010, the low lumi data is that for which the HLT_Mu9 trigger path
## was active and unprescaled, (uncertified) run range 133446 - 147116. Certification
## restricts this run range further.
##
isRun2010HiLumi=not isRun2010LoLumi

process = cms.Process("PAT");

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)
# source
process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring('file:/hdfs/cms/skim/mu/39X/Dec22ReReco/Run2010B/HiLumi/pool_1_1_4Cv.root')
#    fileNames=cms.untracked.vstring('file:/local/cms/phedex/store/mc/Fall10/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/START38_V12-v3/0001/02D6EA60-9D02-E011-A091-90E6BA442F11.root')
)
#process.load('HeavyNu.AnalysisModules.in_cff')

if isData:
    if isRun2010LoLumi:
        print "===========> Flag is SET for LOW luminosity data <============"
    else:
        print "===========> Flag is SET for HIGH luminosity data <============"
    
    from HeavyNu.AnalysisModules.goodRunList_cfi import lumisToProcess
    process.source.lumisToProcess = lumisToProcess

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Load additional processes
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
## global tags:
#
# Use 38X for MC, 39X for data
#
if (isMC):
    print "=================> MC flag is SET <===================="
    process.GlobalTag.globaltag = cms.string('START38_V14::All')
else:
    process.GlobalTag.globaltag = cms.string('GR_R_39X_V6::All')

process.load("Configuration.StandardSequences.MagneticField_cff")

################################################################################################
###    P r e p a r a t i o n      o f    t h e    P A T    O b j e c t s   f r o m    A O D  ###
################################################################################################

## pat sequences to be loaded:
process.load("PhysicsTools.PatAlgos.patSequences_cff")

########################################
# Output module - has to be defined before PAT python tools will work
########################################

if isData:
    from PhysicsTools.PatAlgos.tools.coreTools import *
    removeMCMatching(process, ['All'], outputInProcess = False)

########################################
# PAT Jet Energy Corrections - MC vs Data
########################################
# 

from PhysicsTools.PatAlgos.tools.jetTools import switchJetCollection
if isMC:
    switchJetCollection( process,
                         jetCollection=cms.InputTag('ak5CaloJets'),
                         jetCorrLabel=('AK5Calo', ['L2Relative','L3Absolute']))
else:
    switchJetCollection( process,
                         jetCollection=cms.InputTag('ak5CaloJets'),
                         jetCorrLabel=('AK5Calo', ['L2Relative','L3Absolute','L2L3Residual']))

########################################
# PAT Trigger matching
########################################
# imported directly from PhysicsTools/PatExamples/test/analyzePatTrigger_onTheFly_cfg.py
#
process.load("HeavyNu.AnalysisModules.hnutrigmatch_cfi")

### ============
### Python tools
### ============
### Attention: order matters!

## --
## Switch to selected PAT objects in the main work flow
## --
from PhysicsTools.PatAlgos.tools.coreTools import removeCleaning
removeCleaning( process, False )

## --
## Switch on PAT trigger - but only for data!
## --
from PhysicsTools.PatAlgos.tools.trigTools import *
if isData:
    switchOnTriggerMatching( process, triggerMatchers = [ 'muonTriggerMatchHLTMuons' ], outputModule = '' )
    removeCleaningFromTriggerMatching( process, outputModule = '' )
    if isRun2010LoLumi: process.muonTriggerMatchHLTMuons.pathNames = cms.vstring('HLT_Mu9')
    else:               process.muonTriggerMatchHLTMuons.pathNames = cms.vstring('HLT_Mu15_v1')

##########################################
## Add analysis
##########################################

process.TFileService = cms.Service("TFileService",
       fileName = cms.string("anal.root"),
)

process.load("HeavyNu.AnalysisModules.heavynuanalysis_cfi")
process.hNu.minMu2pt         = cms.double(30.)
process.hNu.studyScaleFactor = cms.bool(True)
if isMC:
    process.hNu.applyMuIDEffcorr = cms.bool(True)
else:
    process.hNu.applyMuIDEffcorr = cms.bool(False)

if isData:
    # turn on trigger match requirement
    process.hNu.trigMatchPset.trigEventTag=cms.InputTag("patTriggerEvent")
    process.hNu.trigMatchPset.muonMatch=cms.string('muonTriggerMatchHLTMuons')
else:
    # turn on MC trigger simulation
    process.hNu.trigMatchPset.randomSeed=cms.int32(os.getpid())

if Training:
    process.hNu.trainingFileName=cms.untracked.string("changeme_nntraining.txt")
    
## ---
## Define the paths
## ---
if isMC:
    process.hNu2011Z20 = process.hNu.clone(minMu2pt = cms.double(20.))
    process.hNu2011Z30 = process.hNu.clone(minMu2pt = cms.double(30.))
    process.hNu2010Z20 = process.hNu.clone(minMu2pt = cms.double(20.))
    process.hNu2010Z30 = process.hNu.clone(minMu2pt = cms.double(30.))
    process.hNu2010Z20.applyMuIDEffcorr     = cms.bool(False)
    process.hNu2010Z20.trigMatchPset.year   = cms.int32(2010)
    process.hNu2010Z30.applyMuIDEffcorr     = cms.bool(False)
    process.hNu2010Z30.trigMatchPset.year   = cms.int32(2010)
    process.hNu2011Z20.highestPtTriggerOnly = cms.bool(True)
else:
    process.hNu20 = process.hNu.clone(minMu2pt = cms.double(20.))
    if isRun2011:
        process.hNu20.highestPtTriggerOnly = cms.bool(True)
    process.hNu30 = process.hNu.clone(minMu2pt = cms.double(30.))

if isMC:
    process.p2010Z20 = cms.Path(process.patDefaultSequence+process.hNu2010Z20)
    process.p2010Z30 = cms.Path(process.patDefaultSequence+process.hNu2010Z30)
    process.p2011Z20 = cms.Path(process.patDefaultSequence+process.hNu2011Z20)
    process.p2011Z30 = cms.Path(process.patDefaultSequence+process.hNu2011Z30)
    process.s = cms.Schedule(process.p2010Z20,
                             process.p2010Z30,
                             process.p2011Z20,
                             process.p2011Z30)
else:
    process.pZ20 = cms.Path(process.patDefaultSequence+process.hNu20)
    process.pZ30 = cms.Path(process.patDefaultSequence+process.hNu30)
    process.s = cms.Schedule(process.pZ20,
                             process.pZ30)