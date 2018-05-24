import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("SiPixelTemplateDBUpload")
#process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.load("CondCore.CondDB.CondDB_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")

#process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")
process.load("Configuration.Geometry.GeometryExtended2017Reco_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '75X_upgrade2017_design_v4', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '80X_upgrade2017_design_v10', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_design', '')

files_to_upload = cms.vstring(
        "../../DB/MC/template_summary_zp0030.out",
        "../../DB/MC/template_summary_zp0031.out"
)
theDetIds      = cms.vuint32( 1, 2)
theTemplateIds = cms.vuint32(30,31)

process.source = cms.Source("EmptyIOVSource",
                            timetype = cms.string('runnumber'),
                            firstValue = cms.uint64(1),
                            lastValue = cms.uint64(1),
                            interval = cms.uint64(1)
                            )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
    )

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
   DBParameters = cms.PSet(
    messageLevel = cms.untracked.int32(0),
    authenticationPath = cms.untracked.string('.')
    ),
    timetype = cms.untracked.string('runnumber'),
    connect = cms.string('sqlite_file:SiPixelTemplateDBObject_phase1_38T_mc_v1.db'),
    toPut = cms.VPSet(cms.PSet(
     record = cms.string('SiPixelTemplateDBObjectRcd'),
     tag = cms.string('SiPixelTemplateDBObject_phase1_38T_mc_v1')
    ))
   )

#process.uploader = cms.EDAnalyzer("SiPixelTemplateDBObjectUploader", # CondTools 
process.uploader = cms.EDAnalyzer("SiPixelTemplateDBLoader", # DPGAnalysis 
                                  siPixelTemplateCalibrations = files_to_upload,
                                  #theTemplateBaseString = cms.string(template_base),
                                  Version = cms.double("1.0"),
                                  #MagField = cms.double(MagFieldValue),
                                  detIds = theDetIds,
                                  templateIds = theTemplateIds
)

process.p = cms.Path(process.uploader)
