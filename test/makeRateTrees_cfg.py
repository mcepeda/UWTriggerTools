#flake8: noqa
'''

Make ntuples with information about the rates.

Usage:

    ./makeRateTrees_cfg.py

Optional arguments:

    inputFiles=myFile.root outputFile=outputFile.root maxEvents=-1

Authors: L. Dodd, N. Woods, I. Ojalvo, S. Dasu, M. Cepeda, E. Friis (UW Madison)

'''

import FWCore.ParameterSet.Config as cms
import os

# Get command line options
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
# Set useful defaults
#options.inputFiles = '/store/user/tapas/ETauSkim/skim_12_1_erV.root'
#options.inputFiles = '/store/user/tapas/2012-08-01-CRAB_ZEESkim/skim_10_1_wd2.root'
options.outputFile = "uct_rate_tree.root"
options.register(
    'eicIsolationThreshold',
    3,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "EIC Isolation threshold")
options.register(
    'hActivityCut',
    0.5,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "HCAL activity threshold")
options.register(
    'ecalCalib',
    'CALIB_V4',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    'Can be CALIB_V1, CALIB_V3, or CALIB_V4')
options.register(
    'eicCardHcalOnly',
    0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    'If 1, turn off the ECAL for the stage1 EGTau path.')
options.register(
    'isMC',
    0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    'Set to 1 for simulated samples - updates GT, emulates HCAL TPGs.')

options.parseArguments()

process = cms.Process("L1UCTRates")

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
# Load the correct global tag, based on the release
if 'CMSSW_6' in os.environ['CMSSW_VERSION']:
    process.GlobalTag.globaltag = 'POSTLS161_V2::All'
    print "Using global tag for upgrade MC: %s" % process.GlobalTag.globaltag
    if not options.isMC:
        raise ValueError("There is no data in CMSSW 6, you must mean isMC=1")
else:
    if not options.isMC:
        # CMSSW 5 data
        process.GlobalTag.globaltag = 'GR_H_V28::All'
    else:
        # CMSSW 5 MC
        process.GlobalTag.globaltag = 'START53_V7B::All'
    process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
    process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
    print "Using global tag for 52X data: %s" % process.GlobalTag.globaltag

# UNCOMMENT THIS LINE TO RUN ON SETTINGS FROM THE DATABASE
# process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource', 'GlobalTag')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(options.outputFile)
)

# Load emulation and RECO sequences
if not options.isMC:
    process.load("L1Trigger.UCT2015.emulation_cfi")
else:
    process.load("L1Trigger.UCT2015.emulationMC_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")

process.load("L1Trigger.L1ExtraFromDigis.l1extraParticles_cfi")

# Read inst. lumi. info from the scalers
process.load("EventFilter.ScalersRawToDigi.ScalersRawToDigi_cfi")
process.scalersRawToDigi.scalersInputTag = 'rawDataCollector'

# Define the tree producers
process.tauL1Rate = cms.EDAnalyzer(
    "RateTree",
    src = cms.VInputTag(cms.InputTag("l1extraParticles", "Tau")),
    isUCT = cms.bool(False),
)
process.rlxTauUCTRate = cms.EDAnalyzer(
    "RateTree",
    src = cms.VInputTag(cms.InputTag("UCT2015Producer", "RelaxedTauUnpacked")),
    isUCT = cms.bool(True),
)
process.isoTauUCTRate = cms.EDAnalyzer(
    "RateTree",
    src = cms.VInputTag(cms.InputTag("UCT2015Producer", "IsolatedTauUnpacked")),
    isUCT = cms.bool(True),

)
process.isoEGL1Rate = cms.EDAnalyzer(
    "RateTree",
    src = cms.VInputTag(cms.InputTag("l1extraParticles", "Isolated")),
    isUCT = cms.bool(False),
)
process.isoEGUCTRate = cms.EDAnalyzer(
    "RateTree",
    src = cms.VInputTag(cms.InputTag("UCT2015Producer", "IsolatedEGUnpacked")),
    isUCT = cms.bool(True),
)
process.rlxEGL1Rate = cms.EDAnalyzer(
    "RateTree",
    src = cms.VInputTag(
        cms.InputTag("l1extraParticles", "NonIsolated"),
        cms.InputTag("l1extraParticles", "Isolated")
    ),
    isUCT = cms.bool(False),
)
process.rlxEGUCTRate = cms.EDAnalyzer(
    "RateTree",
    src = cms.VInputTag(cms.InputTag("UCT2015Producer", "RelaxedEGUnpacked")),
    isUCT = cms.bool(True),
)

process.uctLeptonRates = cms.Sequence(
    process.rlxTauUCTRate *
    process.isoTauUCTRate *
    process.rlxEGUCTRate *
    process.isoEGUCTRate
)
process.jetL1Rate = cms.EDAnalyzer(
    "RateTree",
    src = cms.VInputTag(
        # Combine central and forward
        cms.InputTag("l1extraParticles", "Central"),
        cms.InputTag("l1extraParticles", "Tau"),
        cms.InputTag("l1extraParticles", "Forward"),
    ),
    isUCT = cms.bool(False),
)
process.jetUCTRate = cms.EDAnalyzer(
    "RateTree",
    src = cms.VInputTag(cms.InputTag("UCT2015Producer", "JetUnpacked")),
    isUCT = cms.bool(True),
)

process.corrjetUCTRate = cms.EDAnalyzer(
    "RateTree",
    src = cms.VInputTag(cms.InputTag("UCT2015Producer", "CorrJetUnpacked")),
    isUCT = cms.bool(True),
)

process.sumsL1Rates = cms.EDAnalyzer(
    "SumsRateTree",
    l1MHTSrc = cms.InputTag("l1extraParticles", "MHT"),
    l1METSrc = cms.InputTag("l1extraParticles", "MET"),
    # Placeholder
    l1METSigSrc = cms.InputTag("UCT2015Producer", "METSIGUnpacked"),
    # fixme
    #l1METSigSrc = cms.InputTag("UCT2015Producer", "METSIGUnpacked"),

    l1SHTSrc = cms.InputTag("l1extraParticles", "MHT"),
    l1SETSrc = cms.InputTag("l1extraParticles", "MET"),
)

process.sumsUCTRates = cms.EDAnalyzer(
    "SumsRateTree",
    l1MHTSrc = cms.InputTag("UCT2015Producer", "MHTUnpacked"),
    l1METSrc = cms.InputTag("UCT2015Producer", "METUnpacked"),
    l1METSigSrc = cms.InputTag("UCT2015Producer", "METSIGUnpacked"),
    # fixme
    #l1METSigSrc = cms.InputTag("UCT2015Producer", "METSIGUnpacked"),

    l1SHTSrc = cms.InputTag("UCT2015Producer", "MHTUnpacked"),
    l1SETSrc = cms.InputTag("UCT2015Producer", "METUnpacked"),
)

process.uctHadronicRates = cms.Sequence(
    process.corrjetUCTRate*
    process.jetUCTRate *
    process.sumsUCTRates*
    process.jetL1Rate
    #process.sumsL1Rates  this breaks - why?
)


process.p1 = cms.Path(
    process.emulationSequence *
    process.scalersRawToDigi *
    process.l1extraParticles*
    process.tauL1Rate *
    process.isoEGL1Rate *
    process.rlxEGL1Rate 
)

print "Building stage1 trees"
process.p1 += process.uctLeptonRates
process.p1 += process.uctHadronicRates

# Make the framework shut up.
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Spit out filter efficiency at the end.
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

#eic = options.eicIsolationThreshold
#print "Setting EIC threshold to %i" % eic
#process.RCTConfigProducers.eicIsolationThreshold = eic
#hActivity = options.hActivityCut
#print "Setting hActivity threshold to %f" % hActivity
#process.RCTConfigProducers.hActivityCut = hActivity
