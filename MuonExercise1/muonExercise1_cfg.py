import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonExercise1")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(0)
)

process.maxEvents = cms.untracked.PSet( 
    #input = cms.untracked.int32(-1) 
    input = cms.untracked.int32(20000)
)

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False),
SkipEvent = cms.untracked.vstring('ProductNotFound'))

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring('root://cmseos.fnal.gov///store/user/cmsdas/2018/short_exercises/Muons/Samples/dymm.root')
)

process.demo = cms.EDAnalyzer('MuonExercise1')


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('histos3.root')
)

process.p = cms.Path(process.demo)
