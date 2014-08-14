import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:/uscms_data/d3/cgrud85/new/CMSSW_5_3_18/src/MyAnalysis/Vertex/AOD/res/thisone_81_1_rBO.root'
#		"/store/user/cgrud/HTo2LongLivedTo4L/HiggsToFourMu_10GeV_17mm_RECO_VertObjects/09f68ee19c332f4c7af327b8c7882771/thisone_100_1_yGu.root"
#		'file:/uscms_data/d3/cgrud85/new/CMSSW_5_3_18/src/MyAnalysis/Vertex/this_other_one.root'
#		'file:/uscms_data/d3/cgrud85/new/CMSSW_5_3_18/src/MyAnalysis/Vertex/Mu_2_GeV_1_mm_with_2_new_cut.root'
		'file:/eos/uscms/store/user/cgrud85/ROOTFiles/Mu_2_GeV_10_cm_with_2_new_cut.root' 
    )
)

process.TFileService = cms.Service("TFileService",
      fileName = cms.string("thatplot_2_GeV_10_cm_with_2_new_cut.root")
)

process.demo = cms.EDAnalyzer('SecondaryVertices'
)


process.p = cms.Path(process.demo)
