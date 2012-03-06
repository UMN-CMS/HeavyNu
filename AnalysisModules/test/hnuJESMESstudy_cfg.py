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
    fileNames=cms.untracked.vstring('file:input.root')
)
# process.load('HeavyNu.AnalysisModules.in_cff')

if isData:
    if isRun2010LoLumi:
        print "===========> Flag is SET for LOW luminosity data <============"
        from HeavyNu.AnalysisModules.run2010loLumiRunList_cfi import lumisToProcess
    else:
        print "===========> Flag is SET for HIGH luminosity data <============"
        from HeavyNu.AnalysisModules.run2010hiLumiRunList_cfi import lumisToProcess
    
    process.source.lumisToProcess = lumisToProcess

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Load additional processes
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
## global tags:
if (isMC):
    print "=================> MC flag is SET <===================="
    process.GlobalTag.globaltag = cms.string('START38_V14::All')
else:
    process.GlobalTag.globaltag = cms.string('GR_R_38X_V15::All')

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

##########################################

# JEC uncertainty application
# applyJECUsign=0 means don't apply it, which is the default
#
process.hNuJEShi = process.hNu.clone(applyJECUsign = cms.int32(1))
process.hNuJESlo = process.hNu.clone(applyJECUsign = cms.int32(-1))

process.pJEShi = cms.Path(process.hNuJEShi)
process.pJESlo = cms.Path(process.hNuJESlo)

##########################################

# MES uncertainty application
# applyMESfactor=1.0 is the obvious default
#
process.hNuMEShi = process.hNu.clone(applyMESfactor = cms.double(1.01))
process.hNuMESlo = process.hNu.clone(applyMESfactor = cms.double(0.99))

process.pMEShi = cms.Path(process.hNuMEShi)
process.pMESlo = cms.Path(process.hNuMESlo)

##########################################

# MuID uncertainty application
# applyMESfactor=1.0 is the obvious default
#
process.hNuMuIDhi = process.hNu.clone(applyMuIDEffsign = cms.int32(1))
process.hNuMuIDlo = process.hNu.clone(applyMuIDEffsign = cms.int32(-1))

process.pMuIDhi = cms.Path(process.hNuMuIDhi)
process.pMuIDlo = cms.Path(process.hNuMuIDlo)

##########################################

process.s = cms.Schedule(process.p,
                         process.pJEShi,process.pJESlo,
                         process.pMEShi,process.pMESlo,
                         process.pMuIDhi,process.pMuIDlo)