# masking correction
initializePuCorr -t data_8TeV_17.2_VtxLumi_201351 >& log_initializePuCorr_201351_04112014_v3.txt &
initializePuCorr -t data_8TeV_17.2_VtxLumi_207216 >& log_initializePuCorr_207216_04112014_v3.txt &
initializePuCorr -t data_8TeV_17.2_VtxLumi_207219 >& log_initializePuCorr_207219_04112014_v3.txt &
initializePuCorr -t data_8TeV_17.2_VtxLumi_214984 >& log_initializePuCorr_214984_04112014_v3.txt &
initializePuCorr -t data_8TeV_17.2_VtxLumi_215021 >& log_initializePuCorr_215021_04112014_v3.txt &
initializePuCorr -t mc_8TeV_17.2_VtxLumi_BothSamples >& log_initializePuCorr_MC_BothSamples_04112014_v3.txt &

#fake correction
initializeFakeCorrection -t mc_8TeV_17.2_VtxLumi_BothSamples >& log_initializeFakeCorrection_MC_BothSamples_04112014_v3.txt &

#Masking correction plots (the one from fake correction gets done in the cxx file)
python MaskingCorrectionPlots.py

#Make correction plots (masking and fake together so mu_rec vs mu_vis)
python PileupCorrectionPlots.py