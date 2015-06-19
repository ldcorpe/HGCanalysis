from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'phoID_270515_HGC'
config.General.workArea = 'crab_projects'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'hoverEHggConfig_MVA.py'
#config.JobType.pyCfgParams = ['output=alcareco.root', 'skim=none', 'type=ALCARECO','doTree=0', 'jsonFile=json.txt', 'secondaryOutput=ntuple.root', 'isCrab=1']
config.JobType.allowUndistributedCMSSW = True
config.Data.inputDataset = '/GJet_Pt-15to3000_Tune4C_14TeV_pythia8/TP2023HGCALDR-HGCALMar26_PU140BX25_PH2_1K_FB_V6-v1/GEN-SIM-RECO'
#config.JobType.inputFiles = ['PH1_A0_PU50_PhotonID_BDTG.weights.xml']
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '' #you need to specific your output folder.
config.Data.publication = False
config.Site.storageSite = "T2_CH_CERN"
