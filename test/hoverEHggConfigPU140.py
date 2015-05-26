import FWCore.ParameterSet.Config as cms

process = cms.Process("HGCSimHitsAnalysis")

#process.load('Configuration.StandardSequences.Services_cff')
#process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')    
process.load('FWCore.MessageService.MessageLogger_cfi')
#v6 geometry
#process.load('Configuration.Geometry.GeometryExtended2023HGCalV6MuonReco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023HGCalV6Muon_cff')
#v5 geometry
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')
#v4 geometry
#process.load('Configuration.Geometry.GeometryExtended2023HGCalV4MuonReco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023HGCalV4Muon_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("RecoParticleFlow.PFClusterProducer.particleFlowRecHitHGCEE_cfi")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')


process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False)
                                        #SkipEvent = cms.untracked.vstring('ProductNotFound')
                                        ) 

# configure from command line
# cmsRun test/runHGCHitsAnalyzer_cfg.py tag
# where tag can be any sub-directory under /store/cmst3/group/hgcal/CMSSW
#           or any upgrade relval sample (may need tweaking for new releases...)
ffile=0
step=-1
preFix='Single13_CMSSW_6_2_0_SLHC18'
doFullAnalysis=True
import os,sys
#if(len(sys.argv)<3):
#    print '\ncmsRun runHGCHitsAnalyzer_cfg.py doFullAnalysis tag first_file step\n'
#    print '\ttag - process tag'
#    print '\tfirst_file - first file to process'
#    print '\tstep - number of files to process\n'
#    sys.exit()

#preFix=sys.argv[2]
#if(len(sys.argv)>3):
#    if(sys.argv[3].isdigit()) : ffile=int(sys.argv[3])
#if(len(sys.argv)>4):
#    if(sys.argv[4].isdigit()) : step=int(sys.argv[4])
#print '[runHGCHitsAnalyzer] processing %d files of %s, starting from %d'%(step,preFix,ffile)

#configure the source (list all files in directory within range [ffile,ffile+step[
from UserCode.HGCanalysis.storeTools_cff import fillFromStore
#process.source = cms.Source("PoolSource",                            
#                            fileNames=cms.untracked.vstring()
#                            )
#if preFix.find('/store')>=0 :
#    process.source.fileNames=fillFromStore(preFix,ffile,step)
#else :
#process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("file:/afs/cern.ch/user/l/lcorpe/work/public/HGCAL/SingleElectronPt35_PU0_RECO_1.root",
#"file:/afs/cern.ch/user/l/lcorpe/work/public/HGCAL/SingleElectronPt35_PU0_RECO_2.root",
#"file:/afs/cern.ch/user/l/lcorpe/work/public/HGCAL/SingleElectronPt35_PU0_RECO_3.root",
#"file:/afs/cern.ch/user/l/lcorpe/work/public/HGCAL/SingleElectronPt35_PU0_RECO_4.root",
#"file:/afs/cern.ch/user/l/lcorpe/work/public/HGCAL/SingleElectronPt35_PU0_RECO_5.root",
#"file:/afs/cern.ch/user/l/lcorpe/work/public/HGCAL/SingleElectronPt35_PU0_RECO_6.root",
#"file:/afs/cern.ch/user/l/lcorpe/work/public/HGCAL/SingleElectronPt35_PU0_RECO_7.root",
#"file:/afs/cern.ch/user/l/lcorpe/work/public/HGCAL/SingleElectronPt35_PU0_RECO_8.root",
#"file:/afs/cern.ch/user/l/lcorpe/work/public/HGCAL/SingleElectronPt35_PU0_RECO_9.root"))

#fileNames = open("LCFilenames.txt","r")
#fileNames = open("hovereFile.txt","r")
fileNames = open("hovereFileCombinedPU140.txt","r")
#fileNames = open("hovereFileQCD.txt","r")
#fileNames = open("relval140PU.txt","r")

process_ =0

import os,sys
if(len(sys.argv)>2):
	#print sys.argv[2]
	process_= int(sys.argv[2])
	print 'index %d'%(process_)
if (len(sys.argv) ==0):
	print 'no index! default is 0'

#process.GlobalTag.globaltag = 'auto:upgradePLS3'


process.source = cms.Source("PoolSource",
                            #fileNames=cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/relval/CMSSW_6_2_0_SLHC22/RelValH130GGgluonfusion_14TeV/GEN-SIM-RECO/PH2_1K_FB_V6_UPGHGCalV5-v1/00000/1CC2630B-6A8F-E411-95D3-0025905A48BA.root"),
														#fileNames=cms.untracked.vstring("/store/relval/CMSSW_6_2_0_SLHC23_patch1/RelValH130GGgluonfusion_14TeV/GEN-SIM-RECO/PU_PH2_1K_FB_V6_HGCalV5PU140-v5/00000/0E088A59-4CA1-E411-98F6-003048FFD744.root"),

                            fileNames=cms.untracked.vstring(fileNames),
                            #fileNames=cms.untracked.vstring("file:HggRelval.root"),
                            #fileNames=cms.untracked.vstring("file:/afs/cern.ch/user/l/lcorpe/work/private/HGCALreco3/CMSSW_6_2_0_SLHC22/src/Hgg0PU-1kEvents_1.root"),
                            skipEvents=cms.untracked.uint32(process_*500+250))

#process.source.fileNames=fillFromStore('/store/cmst3/group/hgcal/CMSSW/%s'%preFix,ffile,step)
#process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(150) )

#load the analyzer
import getpass
whoami=getpass.getuser()
outputTag=preFix.replace('/','_')
#process.TFileService = cms.Service("TFileService", fileName = cms.string('/tmp/%s/%s_Hits_%d.root'%(whoami,outputTag,ffile)))
#process.TFileService = cms.Service("TFileService", fileName = cms.string('HoverEHggQCD.root'))
process.TFileService = cms.Service("TFileService", fileName = cms.string('HoverEHggPU140_%d_b_v2b.root'%(process_)))
#process.load('UserCode.HGCanalysis.hgcHitsAnalyzer_cfi')

process.hgg = cms.EDAnalyzer("HoverEAnalyzer",
                        #geometrySource   = cms.untracked.vstring('HGCalEESensitive','HGCalHESiliconSensitive',  'HGCalHEScintillatorSensitive')
												geometrySource = cms.untracked.string('HGCalEESensitive'),
												endcapRecHitCollection = cms.untracked.InputTag("HGCalRecHit:HGCEERecHits"),
												endcapSuperClusterCollection = cms.untracked.InputTag("particleFlowSuperClusterHGCEE"),
												endcapClusterCollection = cms.untracked.InputTag("particleFlowClusterHGCEE"),
												genParticlesTag =  cms.untracked.InputTag("genParticles"),
												hcalTowers = cms.InputTag("towerMaker"),
											#	eeRecHitCollection = cms.untracked.InputTag("particleFlowRecHitHGCEELC"),
												hOverEPtMin = cms.double(2.),
												hOverEMethodEndcap = cms.int32(3),
												hOverEConeSize = cms.double(0.15),
												endcapHCALClusters= cms.InputTag("particleFlowClusterHGCHEF"),
												PU = cms.int32(140),
												process_ = cms.int32(process_)
                          )


#run it
process.p = cms.Path(#process.analysis
                    # process.particleFlowRecHitHGCEELC*
										 process.hgg
)

