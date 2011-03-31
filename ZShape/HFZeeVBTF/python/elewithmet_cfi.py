import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('EleWithMet',
# source
  acceptedElectronIDs = cms.vint32( 7 )
)
