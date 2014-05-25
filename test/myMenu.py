'''

Creates L1ExtraNtuples (L1 Style) using a UCT->GT jump

Authors: L. Dodd, N. Woods, T. Perry, A. Levine,, S. Dasu, M. Cepeda, E. Friis (UW Madison)

'''

import FWCore.ParameterSet.Config as cms
import os

from FWCore.ParameterSet.VarParsing import VarParsing
process = cms.Process("ReRunningL1")

process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/C2D9064D-FB7F-E311-BFCC-003048679228.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/389085BE-CF7F-E311-A4C7-0026189438EF.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/98E05035-0A80-E311-A784-003048D42D92.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/328D17F8-F87F-E311-ADA9-0025905A6122.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/08B46ECD-CA7F-E311-847A-002618943946.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/20003/18F85924-3D7C-E311-AF32-0026189438D7.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/02B79593-F47F-E311-8FF6-003048FFD796.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/CAD02E5D-A67F-E311-843F-003048FFD756.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/20AE3CAC-CF7F-E311-B9DE-003048678E24.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/3CDB2CF8-0880-E311-8878-002618943935.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/0098522E-D77F-E311-89CC-0025905A610A.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/20002/C2F9C279-117A-E311-A4EB-00261894389A.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/6AF2C1E2-DF7F-E311-B452-003048679162.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/9643B5C2-0280-E311-A66E-00261894394D.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/A0B9E0A6-C87F-E311-92E4-0026189438EF.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/4AEC3425-0580-E311-B5A2-0025905A611E.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/20002/9673477E-D57A-E311-ACD1-002618943983.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/20003/341F0B18-157B-E311-A26C-0025905A6088.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/EAC88847-C27F-E311-8824-0026189438F5.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/0A53FC3E-0780-E311-9D84-0025905A60AA.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/B08574C9-E57F-E311-99A5-003048D15DDA.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/4A4B2971-C07F-E311-843E-0025905A6094.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/1EFB29C6-BE7F-E311-83C7-003048FF9AA6.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/6AD89D7E-A57F-E311-B867-0025905A611E.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/AA99D298-B47F-E311-AB06-002618943901.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/9EF8E9D8-AA7F-E311-8589-0025905A60CA.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/20003/BAB14DA6-CF7C-E311-BDBA-0026189438A7.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/20002/5068EC78-B47A-E311-8CEB-003048678FAE.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/20001/306835C8-A679-E311-9A07-00261894387B.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/0CB7AE7E-A57F-E311-A1FF-00261894382A.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/BA95FC23-A87F-E311-BA73-001A928116D2.root",
"/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/683C7D19-C97F-E311-AF6E-003048FF86CA.root",

                             )
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Tested on Monte Carlo, for a test with data edit ahead
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'POSTLS161_V2::All'

# Load emulation and RECO sequences
process.load("L1Trigger.UCT2015.emulationMC_cfi") 
#process.load("L1Trigger.UCT2015.emulation_cfi") # For running on data 
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("L1Trigger.UCT2015.uctl1extraparticles_cfi")
process.uctGctDigis.jetSource =cms.InputTag("UCT2015Producer","CorrJetUnpacked")


process.l1ExtraTreeProducerUCT = cms.EDAnalyzer("L1ExtraTreeProducer",
   nonIsoEmLabel = cms.untracked.InputTag("l1extraParticlesUCT:NonIsolated"),
   isoEmLabel = cms.untracked.InputTag("l1extraParticlesUCT:Isolated"),
   tauJetLabel = cms.untracked.InputTag("l1extraParticlesUCT:Tau"),
   cenJetLabel = cms.untracked.InputTag("l1extraParticlesUCT:Central"),
   fwdJetLabel = cms.untracked.InputTag("l1extraParticlesUCT:Forward"),
   muonLabel = cms.untracked.InputTag("l1extraParticlesUCT"),
   metLabel = cms.untracked.InputTag("l1extraParticlesUCT:MET"),
   mhtLabel = cms.untracked.InputTag("l1extraParticlesUCT:MHT"),
   hfRingsLabel = cms.untracked.InputTag("l1extraParticlesUCT"),
   maxL1Extra = cms.uint32(20)
)

process.l1ExtraTreeProducer = cms.EDAnalyzer("L1ExtraTreeProducer",
   nonIsoEmLabel = cms.untracked.InputTag("l1extraParticles:NonIsolated"),
   isoEmLabel = cms.untracked.InputTag("l1extraParticles:Isolated"),
   tauJetLabel = cms.untracked.InputTag("l1extraParticles:Tau"),
   cenJetLabel = cms.untracked.InputTag("l1extraParticles:Central"),
   fwdJetLabel = cms.untracked.InputTag("l1extraParticles:Forward"),
   muonLabel = cms.untracked.InputTag("l1extraParticles"),
   metLabel = cms.untracked.InputTag("l1extraParticles:MET"),
   mhtLabel = cms.untracked.InputTag("l1extraParticles:MHT"),
   hfRingsLabel = cms.untracked.InputTag("l1extraParticles"),
   maxL1Extra = cms.uint32(20)
)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('genRate2.root')
)

process.createGenParticlesEle =cms.EDProducer("FilterGenParticles",
        MinPtThreshold=cms.untracked.double(5),
        GenLevelSelect=cms.untracked.int32(11),
        MaxIsolation=cms.untracked.double(100)
)


process.createGenParticlesEleIso =cms.EDProducer("FilterGenParticles",
        MinPtThreshold=cms.untracked.double(5),
        GenLevelSelect=cms.untracked.int32(11),
        MaxIsolation=cms.untracked.double(1)
)

process.createGenParticlesTau =cms.EDProducer("FilterGenParticles",
        MinPtThreshold=cms.untracked.double(10),
        GenLevelSelect=cms.untracked.int32(15),
        GenLevelStatus=cms.untracked.int32(3),
        MaxIsolation=cms.untracked.double(100)
)


common_ntuple_branches = cms.PSet(
    # Run, lumi, event number
    run = cms.string("id.run"),
    lumi = cms.string("id.luminosityBlock"),
    evt = cms.string("id.event"),

    recoPt = cms.string("reco.pt"),
    recoEta = cms.string("reco.eta"),
    recoPhi = cms.string("reco.phi"),

    # Whether there exists a L1/UCT object matched to reco
    l1Match = cms.string("l1Match"),
    l1gMatch = cms.string("l1gMatch"),

    l1Pt = cms.string("? l1Match ? l1.pt : 0"),
    l1Eta = cms.string("? l1Match ? l1.eta : 0"),
    l1Phi = cms.string("? l1Match ? l1.phi : 0"),
    l1Type = cms.string("? l1Match ? l1.type() : -1"),
    # TODO add L1extra eta/phi indices

    l1DPhi = cms.string("? l1Match ? deltaPhi(l1.phi, reco.phi) : -1"),
    l1DR = cms.string("? l1Match ? deltaR(l1.eta, l1.phi, reco.eta, reco.phi) : -1"),

    l1gPt = cms.string("? l1gMatch ? l1g.pt : 0"),
    l1gEta = cms.string("? l1gMatch ? l1g.eta : 0"),
    l1gPhi = cms.string("? l1gMatch ? l1g.phi : 0"),

    l1gDPhi = cms.string("? l1gMatch ? deltaPhi(l1g.phi, reco.phi) : -1"),
    l1gDEta = cms.string("? l1gMatch ? l1g.eta - reco.eta : -10"),
    l1gDR = cms.string("? l1gMatch ? deltaR(l1g.eta, l1g.phi, reco.eta, reco.phi) : -1"),

)

# Specific to EG tau objects
egtau_branches = cms.PSet(
    l1gSecondRegionEt = cms.string("? l1gMatch ? l1g.getFloat('associatedSecondRegionEt', -4) : -2"),
    l1gThirdRegionEt = cms.string("? l1gMatch ? l1g.getFloat('associatedThirdRegionEt', -4) : -2"),
    l1gJetPt = cms.string("? l1gMatch ? l1g.getFloat('associatedJetPt', -4) : -2"),
    l1gEllIso = cms.string("? l1gMatch ? l1g.getInt('ellIsolation', -4) : -2"),
    l1gTauVeto = cms.string("? l1gMatch ? l1g.getInt('tauVeto', -4) : -2"),
    l1gMIP = cms.string("? l1gMatch ? l1g.getInt('mipBit', -4) : -2"),
    l1gIsEle = cms.string("? l1gMatch ? l1g.getInt('isEle', -4) : -2"),
)


process.rlxTauEfficiency = cms.EDAnalyzer(
    "EfficiencyGenTree",
    recoSrc = cms.VInputTag("createGenParticlesTau"),
    l1Src = cms.VInputTag(
        # These two collections
        cms.InputTag("l1extraParticles", "Tau"),
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("l1extraParticlesUCT","Tau")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevel"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches#,egtau_branches
    )
)


process.rlxTauL1TauPlusJetsEfficiency = cms.EDAnalyzer(
    "EfficiencyGenTree",
    recoSrc = cms.VInputTag("createGenParticlesTau"),
    l1Src = cms.VInputTag(
        # These two collections
        cms.InputTag("l1extraParticles", "Tau"),
        cms.InputTag("l1extraParticles", "Central"),
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("l1extraParticlesUCT","Tau")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevel"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches#,egtau_branches
    )
)



process.compareEfficiency = cms.EDAnalyzer(
    "EfficiencyGenTree",
    recoSrc = cms.VInputTag("createGenParticlesEle"),
    l1Src = cms.VInputTag(
        # These two collections
        cms.InputTag("l1extraParticlesUCT","NonIsolated")
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("UCT2015Producer", "RelaxedEGUnpacked")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevel"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches,
    )
)


process.rlxEGEfficiency = cms.EDAnalyzer(
    "EfficiencyGenTree",
    recoSrc = cms.VInputTag("createGenParticlesEle"),
    l1Src = cms.VInputTag(
        # These two collections
        cms.InputTag("l1extraParticles", "NonIsolated"),
        cms.InputTag("l1extraParticles", "Isolated"),
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("l1extraParticlesUCT","NonIsolated")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevel"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches#,egtau_branches
    )
)

process.isoEGEfficiency = cms.EDAnalyzer(
    "EfficiencyGenTree",
    recoSrc = cms.VInputTag("createGenParticlesEleIso"),
    l1Src = cms.VInputTag(
        # These two collections
        cms.InputTag("l1extraParticles", "Isolated"),
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("l1extraParticlesUCT","Isolated")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevel"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches#,egtau_branches
    )
)

process.isoEGEfficiencyGENRLX = cms.EDAnalyzer(
    "EfficiencyGenTree",
    recoSrc = cms.VInputTag("createGenParticlesEle"),
    l1Src = cms.VInputTag(
        # These two collections
        cms.InputTag("l1extraParticles", "Isolated"),
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("l1extraParticlesUCT","Isolated")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevel"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches#,egtau_branches
    )
)

process.jetEfficiency = cms.EDAnalyzer(
    "EfficiencyGenTree",
    recoSrc = cms.VInputTag("cleanGenJets"),
    l1Src = cms.VInputTag(
        # Combine central jets + tau + forward jets
        cms.InputTag("l1extraParticles", "Central"),
        cms.InputTag("l1extraParticles", "Forward"),
    ),
    l1GSrc = cms.VInputTag(cms.InputTag("l1extraParticlesUCT", "Central"), cms.InputTag("l1extraParticlesUCT","Forward")),
    l1GPUSrc = cms.InputTag("UCT2015Producer", "PULevel"),
    # Max DR for RECO-trigger matching
    maxDR = cms.double(0.5),
    # Ntuple configuration
    ntuple = cms.PSet(
        common_ntuple_branches,
    )
)

process.cleanGenJets = cms.EDProducer("GenJetClean")

process.l1SumsEfficiency = cms.EDAnalyzer(
    "SumsEfficiencyTree",
    tree2015 =cms.bool(False),
    l1MHTSrc = cms.InputTag("l1extraParticles", "MHT"),
    l1METSrc = cms.InputTag("l1extraParticles", "MET"),
    # Evan said change l1METSigSrc to match recoMETSigSrc
    l1METSigSrc = cms.InputTag("UCT2015Producer", "METSIGUnpacked"),
    #l1METSigSrc = cms.InputTag("metsignificance"),
    # fixme
    l1SHTSrc = cms.InputTag("l1extraParticles", "MHT"),
    l1SETSrc = cms.InputTag("l1extraParticles", "MET"),
    recoMHTSrc = cms.InputTag("genMetCalo"),
    recoMETSrc = cms.InputTag("genMetCalo"), # calomet
    recoMETSigSrc  = cms.InputTag("genMetCalo"),
    recoSHTSrc = cms.InputTag("genMetCalo"),
    recoSETSrc = cms.InputTag("genMetCalo"),
    recoPFMETSrc = cms.InputTag("genMetTrue"), # pfmet
    recoFile  =cms.bool(False)
)
process.uctSumsEfficiency = cms.EDAnalyzer(
    "SumsEfficiencyTree",
    tree2015 =cms.bool(False),
    l1MHTSrc = cms.InputTag("l1extraParticlesUCT", "MHT"),
    l1METSrc = cms.InputTag("l1extraParticlesUCT", "MET"),
    l1METSigSrc = cms.InputTag("UCT2015Producer", "METSIGUnpacked"),
    l1SHTSrc = cms.InputTag("l1extraParticlesUCT", "MHT"),
    l1SETSrc = cms.InputTag("l1extraParticlesUCT", "MET"),
    recoMHTSrc = cms.InputTag("genMetCalo"),
    recoMETSrc = cms.InputTag("genMetCalo"), # calomet
    recoMETSigSrc  = cms.InputTag("genMetCalo"),
    recoSHTSrc = cms.InputTag("genMetCalo"),
    recoSETSrc = cms.InputTag("genMetCalo"),
    recoPFMETSrc = cms.InputTag("genMetTrue"), # pfmet
    recoFile  =cms.bool(False)
)

process.FullMenuUCT=cms.EDAnalyzer(
    "Menu",
    isUCT = cms.bool(True),
    srcEG = cms.VInputTag(cms.InputTag("UCT2015Producer", "RelaxedEGUnpacked")),
    srcIsoEG = cms.VInputTag(cms.InputTag("UCT2015Producer", "IsolatedEGUnpacked")),
    srcJET = cms.VInputTag(cms.InputTag("UCT2015Producer", "CorrJetUnpacked")),
    srcTAU = cms.VInputTag(cms.InputTag("UCT2015Producer", "RelaxedTauUnpacked")),
    srcMET = cms.InputTag("UCT2015Producer", "METUnpacked"),
    srcMHT = cms.InputTag("UCT2015Producer", "MHTUnpacked"),
    srcSET = cms.InputTag("UCT2015Producer", "SETUnpacked"),
    srcSHT = cms.InputTag("UCT2015Producer", "SHTUnpacked"),
)

process.FullMenuUCTL1Extra=cms.EDAnalyzer(
    "Menu",
    isUCT = cms.bool(False),
    srcEG = cms.VInputTag(cms.InputTag("l1extraParticlesUCT", "NonIsolated"),cms.InputTag("l1extraParticlesUCT", "Isolated")),
    srcIsoEG = cms.VInputTag(cms.InputTag("l1extraParticlesUCT", "Isolated")),
    srcJET = cms.VInputTag(cms.InputTag("l1extraParticlesUCT", "Central"),cms.InputTag("l1extraParticlesUCT", "Forward")),
    srcTAU = cms.VInputTag(cms.InputTag("l1extraParticlesUCT", "Tau")),
    srcMET = cms.InputTag("l1extraParticlesUCT", "MET"),
    srcMHT = cms.InputTag("l1extraParticlesUCT", "MHT"),
    srcSET = cms.InputTag("l1extraParticlesUCT", "MET"),
    srcSHT = cms.InputTag("l1extraParticlesUCT", "MHT"),
)

process.p1 = cms.Path(
    process.emulationSequence *
    process.uct2015L1Extra 
    *process.FullMenuUCTL1Extra
    *process.FullMenuUCT    
)

# Make the framework shut up.
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('out.root'),
    outputCommands = cms.untracked.vstring('drop *',
          'keep *_*_*_ReRunningL1',
          'keep *_l1extraParticles*_*_*') 
)

#process.out = cms.EndPath(process.output)


