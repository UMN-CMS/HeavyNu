import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 2000

process.maxEvents = cms.untracked.PSet( 
	 input = cms.untracked.int32(100000)
)

process.options = cms.untracked.PSet(
                            fileMode = cms.untracked.string('NOMERGE')
)

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames=cms.untracked.vstring(     )
)


# # in this file I collect the lumi sections intervals I am interested in
# from ZShape.HFZeeVBTF.myGoodlumi_cfi import *
# # this is the equivalent of the JSON file 
# process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(
#     #                                                                  propmtRecoJsonAug4
#     # + propmtRecoJsonAug4
#     fullPromptRecoUpTo144114
#     )


# see explanations in : https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideGoodLumiSectionsJSONFile#cmsRun
# in 36x need tags: V01-00-02 FWCore/PythonUtilities     and      V00-04-00 PhysicsTools/PythonAnalysis
import PhysicsTools.PythonAnalysis.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
myLumis = LumiList.LumiList(filename = '/home/hep/franzoni/cmssw/CMSSW_3_6_1_patch3_Z_HF_update7/src/ZShape/HFZeeVBTF/fileLists/Cert_132440-144114_7TeV_StreamExpress_Collisions10_JSON_v3.txt').getCMSSWString().split(',')
process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
process.source.lumisToProcess.extend(myLumis)



process.load("RecoEgamma.EgammaHFProducers.hfRecoEcalCandidate_cfi")
process.hfRecoEcalCandidate.intercept2DCut=0.3


process.TFileService = cms.Service("TFileService",
       fileName = cms.string(''),
)

# to access values of EldId cuts
import ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi

process.IdIsoRej = cms.EDAnalyzer('HFZeeVBTF',
    ECALid = cms.string('simpleEleId90relIso'),
    myName = cms.string('IdIsoRej'),
    DoLog = cms.bool(True),
#   this is instance of the analysis code which matters; keep only enectrons that pass the full selection
    acceptedElectronIDs = cms.vint32( 7 ),
                              
#    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
                              minimal = cms.PSet(
    barrel =  cms.vdouble(9999., 1.0e-02, 9999., 9999., -1, -1, 9999., 9999., 9999., 9999., 9999., 9999., 
                          9999., 1.2e-01, 9999., 9999., 9999., 9999., 0.0, -9999., 9999., 9999., 1, -1, 0.0, 0.0, ),
    endcap =  cms.vdouble(9999., 3.0e-02, 9999., 9999., -1, -1, 9999., 9999., 9999., 9999., 9999., 9999., 
                          9999., 5.0e-02, 9999., 9999., 9999., 9999., 0.0, -9999., 9999., 9999., 1, -1, 0.0, 0.0, ),
    ),
    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
    robust90relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90relIsoEleIDCutsV04.clone()),
    robust85relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85relIsoEleIDCutsV04.clone()),
    robust80relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80relIsoEleIDCutsV04.clone()),
    robust70relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70relIsoEleIDCutsV04.clone()),
    robust60relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60relIsoEleIDCutsV04.clone()),
    robust95cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95cIsoEleIDCutsV04.clone()),
    robust90cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90cIsoEleIDCutsV04.clone()),
    robust85cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85cIsoEleIDCutsV04.clone()),
    robust80cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80cIsoEleIDCutsV04.clone()),
    robust70cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70cIsoEleIDCutsV04.clone()),
    robust60cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60cIsoEleIDCutsV04.clone())
    
                              
)


# EWK analyzer: different Z definitions with one ECAL electron and one HF electron
# ---> this is the instance to run AFTER the Wenu EWK filter
process.IdIso = cms.EDAnalyzer('HFZeeVBTF',
    ECALid = cms.string('simpleEleId90cIso'),
    myName = cms.string('IdIso'),
    DoLog = cms.bool(True),
#   this is instance of the analysis code which I carry along to allow the keeping also of electrons that have not passed conversion rejection
    acceptedElectronIDs = cms.vint32( 3, 7 ),
                              
#    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
    robust90relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90relIsoEleIDCutsV04.clone()),
    robust85relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85relIsoEleIDCutsV04.clone()),
    robust80relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80relIsoEleIDCutsV04.clone()),
    robust70relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70relIsoEleIDCutsV04.clone()),
    robust60relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60relIsoEleIDCutsV04.clone()),
    robust95cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95cIsoEleIDCutsV04.clone()),
    robust90cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90cIsoEleIDCutsV04.clone()),
    robust85cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85cIsoEleIDCutsV04.clone()),
    robust80cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80cIsoEleIDCutsV04.clone()),
    robust70cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70cIsoEleIDCutsV04.clone()),
    robust60cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60cIsoEleIDCutsV04.clone())


)

process.IdRej = cms.EDAnalyzer('HFZeeVBTF',
    ECALid = cms.string('simpleEleId90relIso'),
    myName = cms.string('IdRej'),
    DoLog = cms.bool(True),
#   this is instance of the analysis code which matters; keep only enectrons that pass the full selection
    acceptedElectronIDs = cms.vint32( 5 , 7),
                              
#    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
                              minimal = cms.PSet(
    barrel =  cms.vdouble(9999., 1.0e-02, 9999., 9999., -1, -1, 9999., 9999., 9999., 9999., 9999., 9999., 
                          9999., 1.2e-01, 9999., 9999., 9999., 9999., 0.0, -9999., 9999., 9999., 1, -1, 0.0, 0.0, ),
    endcap =  cms.vdouble(9999., 3.0e-02, 9999., 9999., -1, -1, 9999., 9999., 9999., 9999., 9999., 9999., 
                          9999., 5.0e-02, 9999., 9999., 9999., 9999., 0.0, -9999., 9999., 9999., 1, -1, 0.0, 0.0, ),
    ),
    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
    robust90relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90relIsoEleIDCutsV04.clone()),
    robust85relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85relIsoEleIDCutsV04.clone()),
    robust80relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80relIsoEleIDCutsV04.clone()),
    robust70relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70relIsoEleIDCutsV04.clone()),
    robust60relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60relIsoEleIDCutsV04.clone()),
    robust95cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95cIsoEleIDCutsV04.clone()),
    robust90cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90cIsoEleIDCutsV04.clone()),
    robust85cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85cIsoEleIDCutsV04.clone()),
    robust80cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80cIsoEleIDCutsV04.clone()),
    robust70cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70cIsoEleIDCutsV04.clone()),
    robust60cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60cIsoEleIDCutsV04.clone())
    
                              
)

process.IsoRej = cms.EDAnalyzer('HFZeeVBTF',
    ECALid = cms.string('simpleEleId90relIso'),
    myName = cms.string('IsoRej'),
    DoLog = cms.bool(True),
#   this is instance of the analysis code which matters; keep only enectrons that pass the full selection
    acceptedElectronIDs = cms.vint32( 6 , 7),
                              
#    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
                              minimal = cms.PSet(
    barrel =  cms.vdouble(9999., 1.0e-02, 9999., 9999., -1, -1, 9999., 9999., 9999., 9999., 9999., 9999., 
                          9999., 1.2e-01, 9999., 9999., 9999., 9999., 0.0, -9999., 9999., 9999., 1, -1, 0.0, 0.0, ),
    endcap =  cms.vdouble(9999., 3.0e-02, 9999., 9999., -1, -1, 9999., 9999., 9999., 9999., 9999., 9999., 
                          9999., 5.0e-02, 9999., 9999., 9999., 9999., 0.0, -9999., 9999., 9999., 1, -1, 0.0, 0.0, ),
    ),
    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
    robust90relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90relIsoEleIDCutsV04.clone()),
    robust85relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85relIsoEleIDCutsV04.clone()),
    robust80relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80relIsoEleIDCutsV04.clone()),
    robust70relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70relIsoEleIDCutsV04.clone()),
    robust60relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60relIsoEleIDCutsV04.clone()),
    robust95cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95cIsoEleIDCutsV04.clone()),
    robust90cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90cIsoEleIDCutsV04.clone()),
    robust85cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85cIsoEleIDCutsV04.clone()),
    robust80cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80cIsoEleIDCutsV04.clone()),
    robust70cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70cIsoEleIDCutsV04.clone()),
    robust60cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60cIsoEleIDCutsV04.clone())
    
                              
)






process.Id = cms.EDAnalyzer('HFZeeVBTF',
    ECALid = cms.string('simpleEleId90relIso'),
    myName = cms.string('Id'),
    DoLog = cms.bool(True),
#   this is instance of the analysis code which matters; keep only enectrons that pass the full selection
    acceptedElectronIDs = cms.vint32( 1 , 3, 5, 7),
                              
#    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
                              minimal = cms.PSet(
    barrel =  cms.vdouble(9999., 1.0e-02, 9999., 9999., -1, -1, 9999., 9999., 9999., 9999., 9999., 9999., 
                          9999., 1.2e-01, 9999., 9999., 9999., 9999., 0.0, -9999., 9999., 9999., 1, -1, 0.0, 0.0, ),
    endcap =  cms.vdouble(9999., 3.0e-02, 9999., 9999., -1, -1, 9999., 9999., 9999., 9999., 9999., 9999., 
                          9999., 5.0e-02, 9999., 9999., 9999., 9999., 0.0, -9999., 9999., 9999., 1, -1, 0.0, 0.0, ),
    ),
    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
    robust90relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90relIsoEleIDCutsV04.clone()),
    robust85relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85relIsoEleIDCutsV04.clone()),
    robust80relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80relIsoEleIDCutsV04.clone()),
    robust70relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70relIsoEleIDCutsV04.clone()),
    robust60relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60relIsoEleIDCutsV04.clone()),
    robust95cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95cIsoEleIDCutsV04.clone()),
    robust90cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90cIsoEleIDCutsV04.clone()),
    robust85cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85cIsoEleIDCutsV04.clone()),
    robust80cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80cIsoEleIDCutsV04.clone()),
    robust70cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70cIsoEleIDCutsV04.clone()),
    robust60cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60cIsoEleIDCutsV04.clone())
    
                              
)


process.Iso = cms.EDAnalyzer('HFZeeVBTF',
    ECALid = cms.string('simpleEleId90relIso'),
    myName = cms.string('Iso'),
    DoLog = cms.bool(True),
#   this is instance of the analysis code which matters; keep only enectrons that pass the full selection
    acceptedElectronIDs = cms.vint32( 2, 3, 6, 7 ),
                              
#    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
                              minimal = cms.PSet(
    barrel =  cms.vdouble(9999., 1.0e-02, 9999., 9999., -1, -1, 9999., 9999., 9999., 9999., 9999., 9999., 
                          9999., 1.2e-01, 9999., 9999., 9999., 9999., 0.0, -9999., 9999., 9999., 1, -1, 0.0, 0.0, ),
    endcap =  cms.vdouble(9999., 3.0e-02, 9999., 9999., -1, -1, 9999., 9999., 9999., 9999., 9999., 9999., 
                          9999., 5.0e-02, 9999., 9999., 9999., 9999., 0.0, -9999., 9999., 9999., 1, -1, 0.0, 0.0, ),
    ),
    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
    robust90relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90relIsoEleIDCutsV04.clone()),
    robust85relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85relIsoEleIDCutsV04.clone()),
    robust80relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80relIsoEleIDCutsV04.clone()),
    robust70relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70relIsoEleIDCutsV04.clone()),
    robust60relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60relIsoEleIDCutsV04.clone()),
    robust95cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95cIsoEleIDCutsV04.clone()),
    robust90cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90cIsoEleIDCutsV04.clone()),
    robust85cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85cIsoEleIDCutsV04.clone()),
    robust80cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80cIsoEleIDCutsV04.clone()),
    robust70cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70cIsoEleIDCutsV04.clone()),
    robust60cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60cIsoEleIDCutsV04.clone())
    
                              
)


process.Rej = cms.EDAnalyzer('HFZeeVBTF',
    ECALid = cms.string('simpleEleId90relIso'),
    myName = cms.string('Rej'),
    DoLog = cms.bool(True),
#   this is instance of the analysis code which matters; keep only enectrons that pass the full selection
    acceptedElectronIDs = cms.vint32( 4, 5, 6 , 7 ),
                              
#    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
                              minimal = cms.PSet(
    barrel =  cms.vdouble(9999., 1.0e-02, 9999., 9999., -1, -1, 9999., 9999., 9999., 9999., 9999., 9999., 
                          9999., 1.2e-01, 9999., 9999., 9999., 9999., 0.0, -9999., 9999., 9999., 1, -1, 0.0, 0.0, ),
    endcap =  cms.vdouble(9999., 3.0e-02, 9999., 9999., -1, -1, 9999., 9999., 9999., 9999., 9999., 9999., 
                          9999., 5.0e-02, 9999., 9999., 9999., 9999., 0.0, -9999., 9999., 9999., 1, -1, 0.0, 0.0, ),
    ),
    robust95relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95relIsoEleIDCutsV04.clone()),
    robust90relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90relIsoEleIDCutsV04.clone()),
    robust85relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85relIsoEleIDCutsV04.clone()),
    robust80relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80relIsoEleIDCutsV04.clone()),
    robust70relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70relIsoEleIDCutsV04.clone()),
    robust60relIsoEleIDCutsV04 = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60relIsoEleIDCutsV04.clone()),
    robust95cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust95cIsoEleIDCutsV04.clone()),
    robust90cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust90cIsoEleIDCutsV04.clone()),
    robust85cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust85cIsoEleIDCutsV04.clone()),
    robust80cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust80cIsoEleIDCutsV04.clone()),
    robust70cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust70cIsoEleIDCutsV04.clone()),
    robust60cIsoEleIDCutsV04   = cms.PSet(ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi.simpleCutBasedElectronID.robust60cIsoEleIDCutsV04.clone())
    
                              
)


from ZShape.HFZeeVBTF.elewithmet_cfi import *
#myDesiredMetCollection = "caloMET"
myDesiredMetCollection = "pfMET"
#myDesiredMetCollection = "tcMET"

process.metEleIdIsoRej = cms.EDAnalyzer('EleWithMet',
                                 # source
                                 acceptedElectronIDs = cms.vint32( 7 ),
                                 ECALid = cms.string('simpleEleId90relIso'),
                                 myName = cms.string('metEleIdIsoRej'),
                                 DoLog = cms.bool(True),
                                 ETCut = cms.double(20.),
                                 METCut = cms.double(25.),
                                 metCollectionTag = cms.InputTag(myDesiredMetCollection,"","PAT"),
                                 )


process.metEleIdIso = cms.EDAnalyzer('EleWithMet',
                                 # source
                                 acceptedElectronIDs = cms.vint32( 3, 7 ),
                                 ECALid = cms.string('simpleEleId90relIso'),
                                 myName = cms.string('metEleIdIso'),
                                 DoLog = cms.bool(True),
                                 ETCut = cms.double(20.),
                                 METCut = cms.double(25.),
                                 metCollectionTag = cms.InputTag(myDesiredMetCollection,"","PAT"),
                                 )


process.metEleIdRej = cms.EDAnalyzer('EleWithMet',
                                 # source
                                 acceptedElectronIDs = cms.vint32( 5, 7 ),
                                 ECALid = cms.string('simpleEleId90relIso'),
                                 myName = cms.string('metEleIdRej'),
                                 DoLog = cms.bool(True),
                                 ETCut = cms.double(20.),
                                 METCut = cms.double(25.),
                                 metCollectionTag = cms.InputTag(myDesiredMetCollection,"","PAT"),
                                 )


process.metEleIsoRej = cms.EDAnalyzer('EleWithMet',
                                 # source
                                 acceptedElectronIDs = cms.vint32( 6, 7 ),
                                 ECALid = cms.string('simpleEleId90relIso'),
                                 myName = cms.string('metEleIsoRej'),
                                 DoLog = cms.bool(True),
                                 ETCut = cms.double(20.),
                                 METCut = cms.double(25.),
                                 metCollectionTag = cms.InputTag(myDesiredMetCollection,"","PAT"),
                                 )


process.metEleId = cms.EDAnalyzer('EleWithMet',
                                 # source
                                 acceptedElectronIDs = cms.vint32( 1, 3, 5 ,7 ),
                                 ECALid = cms.string('simpleEleId90relIso'),
                                 myName = cms.string('metEleId'),
                                 DoLog = cms.bool(True),
                                 ETCut = cms.double(20.),
                                 METCut = cms.double(25.),
                                 metCollectionTag = cms.InputTag(myDesiredMetCollection,"","PAT"),
                                 )


process.metEleIso = cms.EDAnalyzer('EleWithMet',
                                 # source
                                 acceptedElectronIDs = cms.vint32( 2, 3, 6, 7 ),
                                 ECALid = cms.string('simpleEleId90relIso'),
                                 myName = cms.string('metEleIso'),
                                 DoLog = cms.bool(True),
                                 ETCut = cms.double(20.),
                                 METCut = cms.double(25.),
                                 metCollectionTag = cms.InputTag(myDesiredMetCollection,"","PAT"),
                                 )


process.metEleRej = cms.EDAnalyzer('EleWithMet',
                                 # source
                                 acceptedElectronIDs = cms.vint32( 4, 5, 6, 7 ),
                                 ECALid = cms.string('simpleEleId90relIso'),
                                 myName = cms.string('metEleRej'),
                                 DoLog = cms.bool(True),
                                 ETCut = cms.double(20.),
                                 METCut = cms.double(25.),
                                 metCollectionTag = cms.InputTag(myDesiredMetCollection,"","PAT"),
                                 )


import FWCore.Modules.printContent_cfi
process.dumpEv = FWCore.Modules.printContent_cfi.printContent.clone()

process.p = cms.Path(
    # process.EG_1e28 *
    process.hfRecoEcalCandidate
    *process.IdIso
    *process.IdRej
    *process.IsoRej
    *process.Id
    *process.Iso
    *process.Rej
    *process.IdIsoRej
    * process.metEleIdIsoRej
    * process.metEleIdIso
    * process.metEleIdRej
    * process.metEleIsoRej
    * process.metEleId
    * process.metEleIso
    * process.metEleRej
    # *process.dumpEv
    )
