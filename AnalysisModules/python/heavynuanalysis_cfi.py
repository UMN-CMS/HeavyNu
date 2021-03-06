import FWCore.ParameterSet.Config as cms

hNu = cms.EDFilter(
    "HeavyNu",
    trigMatchPset = cms.PSet(
        trigEventTag = cms.InputTag( "" ),
        muonTriggers = cms.vstring( '' ),
        electronTriggers = cms.vstring( '' ),
        # trigEra      = cms.int32( 20111 ),
        firstRun     = cms.vint32( 0 ),
        lastRun      = cms.vint32( 999999 ),
        muonMatch    = cms.string( '' ),
        electronMatch = cms.string( '' ),
        triggerPt    = cms.double( 40. ),
        muTriggerEta = cms.double( 2.1 ),
        randomSeed   = cms.int32( 0 ),  # for MC
        electronFilters = cms.vstring('')
    ),
    muIDPset = cms.PSet(
        eraForId     = cms.int32( 1 )
    ),
    DoLog        = cms.bool( False ),
    muonTag      = cms.InputTag( 'selectedPatMuons' ),
    jetTag       = cms.InputTag( 'selectedPatJets' ),
    metTag       = cms.InputTag( 'patMETs' ),
    electronTag  = cms.InputTag( 'selectedPatElectrons' ),
    trackTag     = cms.InputTag( 'patTracksPt10' ),
    electronRho  = cms.InputTag( 'kt6PFJetsForIsolation','rho' ),
    BtagName     = cms.string('combinedSecondaryVertexBJetTags'),
    minBtagDiscr = cms.double(0.679), # yields 2.5% fake rate
    minMu1pt     = cms.double(60.),
    minMu2pt     = cms.double(40.), # 2012 defaults
    minJetPt     = cms.double(40),
    maxMuAbsEta  = cms.double(2.4),
    maxElecAbsEta = cms.double(2.5),
    maxJetAbsEta = cms.double(2.5),
    minMuMuMass  = cms.double(200),
    min4objMass  = cms.double(600),
    minMuonJetdR = cms.double(0.5),
    qcdMuJetOverlap = cms.double(0.3),
    
    muonTrackRelIsoLimit  = cms.double(0.1), # 10.0),
    maxVertexZsepCM       = cms.double(0.03),
    maxJetVZsepCM         = cms.double(0.1),
    
    ZmassWinMinGeV= cms.double(86.),
    ZmassWinMaxGeV= cms.double(96.),

    pileupEra         = cms.int32(20112),
    systPileupShift   = cms.int32(0),
    DisableTriggerCorrection = cms.bool(False),

    jecEra            = cms.int32(0),
    applyJECUsign     = cms.int32(0),
    applyJERsign      = cms.int32(0),
    applyTrigEffsign  = cms.int32(0),
    applyMESfactor    = cms.double(1.0),
    checkMERUnc       = cms.bool(False),
    correctEscale     = cms.bool(False),
    applyMuIDEffcorr  = cms.bool(False),
    applyMuIDEffsign  = cms.int32(0),
    applyEERUnc       = cms.bool(False),

    studyMuSelectEff      = cms.bool(True),
    oneTPcand             = cms.bool(False),
    tpRandomSeed          = cms.int32(122674),
    studyScaleFactor      = cms.bool(False),
    studyRatePerRun       = cms.bool(False),
    alternativeSelections = cms.bool(False),
    highestPtTriggerOnly  = cms.bool(False),
    isSignal              = cms.bool(False),
    mNuRnormalization     = cms.double(1000.0),

    doPDFReweight = cms.untracked.bool(False),
#    pdfReweightBaseName = cms.untracked.string('cteq6ll.LHpdf'),
    pdfReweightBaseName = cms.untracked.string('CT10.LHgrid'),
    pdfReweightBaseId = cms.untracked.int32(0),
    pdfReweightTargetName = cms.untracked.string('CT10'),
    pdfReweightTargetId = cms.untracked.int32(0),    

    isPFJets = cms.bool(False),

    #options are HNUMU, HNUE, TOP, QCD, CLO
    analysisMode = cms.untracked.string('HNUMU'),

    nFakeLeptons = cms.untracked.int32(0), 

    heepVersion = cms.untracked.int32(40),
    addSlopeTree = cms.untracked.bool(False),
    noElecIso = cms.untracked.bool(False)
    )
