# Auto generated configuration file
# using: 
# Revision: 1.11 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: RE -s DIGI,L1,DIGI2RAW,RAW2DIGI,L1Reco,RECO --eventcontent FEVTDEBUGHLT --datatier GEN-SIM-DIGI-RAW --conditions DESIGN61_V10::All --geometry ExtendedPhaseIPixel,ExtendedPhaseIPixelReco --customise SLHCUpgradeSimulations/Configuration/customise_mixing.customise_NoCrossing,SLHCUpgradeSimulations/Configuration/phase1TkCustoms.customise --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

my_output_file_name = cms.string( 'PixelDigiRateRecoHits_cff.root' )

my_input_file_names = cms.untracked.vstring()
my_input_file_names.extend( [
    '/store/t0streamer/Data/PhysicsZeroBias2/000/297/562/run297562_ls0075_streamPhysicsZeroBias2_StorageManager.dat'
    #'/store/user/jzabel/Pixel_Upgrade/CMSSW_6_1_1/TuneZ2Star_Generated_MC/MinBias_TuneZ2Star_14TeV_pythia6_R30F12_GEN_SIM_13_02_19/MinBias_TuneZ2Star_14TeV_pythia6_R30F12_GEN_SIM_13_02_19_9_1_E4Q.root'
] );

process = cms.Process('RECO',eras.Run2_2017)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2017Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

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
process.GlobalTag = GlobalTag(process.GlobalTag, 'DESIGN61_V10::All', '')

process.digiHits = cms.EDAnalyzer("PixelDigiRateRecoHits",
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
process = customise_pixelMixing(process)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.phase1TkCustoms
from SLHCUpgradeSimulations.Configuration.phase1TkCustoms import customise 
process = customise(process)
