import FWCore.ParameterSet.Config as cms
readFiles = cms.untracked.vstring()
source = cms.Source("PoolSource",fileNames=readFiles)
readFiles.extend( [
'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1300_nuRmu400/HeavyNuGenHLT_WR1300_nuRmu400_1-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1300_nuRmu400/HeavyNuGenHLT_WR1300_nuRmu400_2-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1300_nuRmu400/HeavyNuGenHLT_WR1300_nuRmu400_3-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1300_nuRmu400/HeavyNuGenHLT_WR1300_nuRmu400_4-reco-pool.root'
, 'file:/local/cms/user/dudero/HeavyNuRecoFromHLT/WR1300_nuRmu400/HeavyNuGenHLT_WR1300_nuRmu400_5-reco-pool.root'
])
