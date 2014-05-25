#!/bin/bash

# Submit rate ntuple jobs on the ZeroBias3 dataset (2012C)
EXPECTED_ARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: $0 JOB_NAME"
  exit 1
fi

farmoutAnalysisJobs $1-RATES \
  --infer-cmssw-path \
  --input-file-list=NeutrinoGun13_Pu4025ns.txt  \
  --input-files-per-job=5 \
  ../myMenu.py \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName' 

farmoutAnalysisJobs $1-EFFIGEN-DY \
  --infer-cmssw-path \
  --input-file-list=DY25ns40PUGENRAW.txt \
  ../makeGenEffi_cfg.py  \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName'
                                
farmoutAnalysisJobs $1-EFFIGEN-H2TAU \
  --infer-cmssw-path \
  --input-file-list=H2TauSkim13TeV.txt \
  ../makeGenEffi_cfg.py  \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName'

farmoutAnalysisJobs $1-EFFIGEN-TT \
  --infer-cmssw-path \
  --input-file-list=TTbar25ns.txt \
  ../makeGenEffi_cfg.py  \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName'

                                        
