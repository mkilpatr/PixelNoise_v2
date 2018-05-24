import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process("RECO",eras.Run2_2017)

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#process.source = cms.Source("PoolSource",
#    # replace 'myfile.root' with the source file you want to use
#    fileNames = cms.untracked.vstring(
#        'file:/store/data/Run2017B/ZeroBias/ALCARECO/LumiPixelsMinBias-PromptReco-v1/000/297/723/00000/10BFE3CA-DF5E-E711-ABD5-02163E01A6F7.root'
#    )
#)
my_output_file_name = cms.string( 'PixelDigiRateRecoHits_cff.root' )

my_input_file_names = cms.untracked.vstring()
my_input_file_names.extend( [    
	#'file:/store/data/Run2017B/ZeroBias/ALCARECO/LumiPixelsMinBias-PromptReco-v1/000/297/723/00000/10BFE3CA-DF5E-E711-ABD5-02163E01A6F7.root'
        '/store/data/Run2017B/ZeroBias1/RAW/v1/000/299/178/00000/8A927A10-286A-E711-A4F7-02163E019BBE.root'
	
] );

#Import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
#process.load('Configuration.Geometry.GeometryExtended2017Plan1Reco_cff')
#process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometryExtended2017_cff')
#process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")


process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.11 $'),
    annotation = cms.untracked.string('RE nevts:1'),
    name = cms.untracked.string('Applications')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 10 ) )

#process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool( True ) )
process.options = cms.untracked.PSet( )

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource", fileNames = my_input_file_names )

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    fileName = cms.untracked.string('file:valid_reco.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW')
    )
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_dataRun2_Express_v2', '')

process.digiHits = cms.EDAnalyzer('George',
				  Verbosity = cms.untracked.bool(False),
                                  src1 = cms.InputTag("simSiPixelDigis"),
                                  src2 = cms.InputTag("siPixelRecHits"),
                                  OutputFile = my_output_file_name,
                                  associatePixel = cms.bool(True),
                                  associateStrip = cms.bool(False),
                                  associateRecoTracks = cms.bool(False),
                                  ROUList = cms.vstring('g4SimHitsTrackerHitsPixelBarrelLowTof',
                                                        'g4SimHitsTrackerHitsPixelBarrelHighTof',
                                                        'g4SimHitsTrackerHitsPixelEndcapLowTof',
                                                        'g4SimHitsTrackerHitsPixelEndcapHighTof')
)

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.user_step = cms.Path(process.digiHits)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.user_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.customise_mixing
from SLHCUpgradeSimulations.Configuration.customise_mixing import customise_pixelMixing
#process = customise_pixelMixing(process)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.phase1TkCustoms
from SLHCUpgradeSimulations.Configuration.phase1TkCustoms import customise
#process = customise(process)


#process.p = cms.Path(process.demo)
