# masking correction
initializePuCorr -t data_8TeV_17.2_VtxLumi_201351 >& log_initializePuCorr_201351_05112014.txt &
initializePuCorr -t data_8TeV_17.2_VtxLumi_207216 >& log_initializePuCorr_207216_05112014.txt &
initializePuCorr -t data_8TeV_17.2_VtxLumi_207219 >& log_initializePuCorr_207219_05112014.txt &
initializePuCorr -t data_8TeV_17.2_VtxLumi_214984 >& log_initializePuCorr_214984_05112014.txt &
initializePuCorr -t data_8TeV_17.2_VtxLumi_215021 >& log_initializePuCorr_215021_05112014.txt &
initializePuCorr -t mc_8TeV_17.2_VtxLumi_BothSamples >& log_initializePuCorr_MC_BothSamples_05112014.txt &

#fake correction
initializeFakeCorrection -t mc_8TeV_17.2_VtxLumi_BothSamples >& log_initializeFakeCorrection_MC_BothSamples_05112014.txt &

#Masking correction plots (the one from fake correction gets done in the cxx file)
python MaskingCorrectionPlots.py

#Make correction plots (masking and fake together so mu_rec vs mu_vis)
python OverallPileupCorrectionPlots.py

#Run VdM analysis
#Careful with the Hack for April (skip last points to fit)
runVdM -r 201351 -s 17.2-VtxLumi -v >& log_runVdM_201351_05112014.txt &
runVdM -r 207216 -s 17.2-VtxLumi -v >& log_runVdM_207216_05112014.txt &
runVdM -r 207219 -s 17.2-VtxLumi -v >& log_runVdM_207219_05112014.txt &
runVdM -r 214984 -s 17.2-VtxLumi -v >& log_runVdM_214984_05112014.txt &
runVdM -r 215021 -s 17.2-VtxLumi -v >& log_runVdM_215021_05112014.txt &

#tables
tables -r 201351 -s 17.2-VtxLumi -m NVtx
tables -r 207216 -s 17.2-VtxLumi -m NVtx
tables -r 207219 -s 17.2-VtxLumi -m NVtx
tables -r 214984 -s 17.2-VtxLumi -m NVtx
tables -r 215021 -s 17.2-VtxLumi -m NVtx

#check consistency
python PrintSigmaVisAndLumiSpVertexCounting.py
python VertexCountingParametersComparisonForNote.py 
python SigmaVisAndLumiSpComparisonForNote.py

#make fits plots
python vdmscans_fits_vtxcounting_perBCID.py