import FWCore.ParameterSet.Config as cms

import os

#import sys
#isMC=sys.modules['__main__'].isMC
#isMCsignal=sys.modules['__main__'].isMCsignal
#process = sys.modules['__main__'].process

isMC=True
isMCsignal=False
Training=False
isRun2010LoLumi=False
isRun2011=True
isPileupMC=True

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
    fileNames=cms.untracked.vstring('file:/hdfs/cms/skim/mu/39X/Dec22ReReco/Run2010B/Mu_Run2010B_Dec22ReReco_v1_AOD_068-muSkim-pool.root')
#This file is AT FNAL:
#    fileNames=cms.untracked.vstring('/store/mc/Fall10/Z1Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola/GEN-SIM-RECO/START38_V12-v2/0018/027050CA-D50B-E011-91DD-00261894391C.root')
#    fileNames=cms.untracked.vstring( "file:/hdfs/cms/user/heavynu/HeavyNuRecoFromHLT/38X/WR1000_nuRmu100/HeavyNuGenHLT_WR1000_nuRmu100_1-reco-pool.root" )
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
    if (isPileupMC):
        process.GlobalTag.globaltag=cms.string('START311_V1G1::All')
        print "=============> isPileupMC flag is SET <================"
    else:
        process.GlobalTag.globaltag = cms.string('START38_V14::All')
else:
    process.GlobalTag.globaltag = cms.string('GR_R_39X_V6::All')

process.load("Configuration.StandardSequences.MagneticField_cff")

################################################################################################
###    P r e p a r a t i o n      o f    t h e    P A T    O b j e c t s   f r o m    A O D  ###
################################################################################################

## pat sequences to be loaded:
process.load("PhysicsTools.PatAlgos.patSequences_cff")

## ---
## Define the path
## ---
process.p = cms.Path(
  process.patDefaultSequence
)

########################################
# Output module - has to be defined before PAT python tools will work
########################################

if isData:
    process.out = cms.OutputModule("PoolOutputModule",
                                   fileName = cms.untracked.string('heavynu_candevents.root'),
    # save only events passing the full path
                                   SelectEvents=cms.untracked.PSet(SelectEvents=cms.vstring('p')),
                                   outputCommands = cms.untracked.vstring("keep *")
                                   )
    process.outpath  = cms.EndPath(process.out)
    from PhysicsTools.PatAlgos.tools.coreTools import *
    removeMCMatching(process, ['All'])

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
removeCleaning( process, isData )

## Special change for saving good products in data
if isData:
    process.out.outputCommands = cms.untracked.vstring("keep *")

## --
## Switch on PAT trigger - but only for data!
## --
from PhysicsTools.PatAlgos.tools.trigTools import *
if isData:
    switchOnTriggerMatching( process, triggerMatchers = [ 'muonTriggerMatchHLTMuons' ] )
    removeCleaningFromTriggerMatching( process )
    if isRun2010LoLumi: process.muonTriggerMatchHLTMuons.pathNames = cms.vstring('HLT_Mu9')
    else:               process.muonTriggerMatchHLTMuons.pathNames = cms.vstring('HLT_Mu15_v1')

##########################################
## Add analysis
##########################################

process.TFileService = cms.Service("TFileService",
       fileName = cms.string("anal.root"),
)

process.load("HeavyNu.AnalysisModules.heavynuanalysis_cfi")
if isRun2011:
    process.hNu.minMu2pt = cms.double(30.)
    process.hNu.jecEra   = cms.int32(3)
    if isMC:
        process.hNu.applyMuIDEffcorr = cms.bool(True)
else:
    process.hNu.minMu2pt = cms.double(20.)
    process.hNu.jecEra   = cms.int32(0)

if isData:
    # turn on trigger match requirement
    process.hNu.trigMatchPset.trigEventTag=cms.InputTag("patTriggerEvent")
    process.hNu.trigMatchPset.muonMatch=cms.string('muonTriggerMatchHLTMuons')
else:
    # turn on MC trigger simulation
    process.hNu.trigMatchPset.randomSeed=cms.int32(os.getpid())

if Training:
    process.hNu.trainingFileName=cms.untracked.string("changeme_nntraining.txt")
    
if isMCsignal:
    process.hNu.isSignal = cms.bool(True)


process.p += process.hNu
