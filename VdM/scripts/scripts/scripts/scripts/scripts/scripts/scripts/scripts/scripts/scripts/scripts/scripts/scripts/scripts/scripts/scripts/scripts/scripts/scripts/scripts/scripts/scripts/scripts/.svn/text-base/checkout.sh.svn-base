# Setup external packages for VtxLumi studies
# Setup test area first

#if [ -z "${TestArea}" ]; then
#    echo "Need to set \$TestArea first"
#    exit
#fi
#
#cd ${TestArea}

#Check-out general packages
svn co $SVNOFF/PhysicsAnalysis/D3PDTools/RootCore/trunk RootCore
svn co $SVNOFF/InnerDetector/InDetValidation/InDetTruthVertexValidation/trunk InDetTruthVertexValidation
svn co $SVNOFF/Trigger/TrigAnalysis/TrigRootAnalysis/tags/TrigRootAnalysis-00-00-08 TrigRootAnalysis
svn co $SVNOFF/DataQuality/GoodRunsLists/tags/GoodRunsLists-00-01-03 GoodRunsLists

#Check-out our packages
#svn co svn+ssh://svn.cern.ch/reps/atlasinst/Institutes/LBNL/Luminosity/VdM/trunk VdM #This package
svn co $SVNINST/Institutes/LBNL/Luminosity/atlasstyle/trunk atlasstyle
svn co $SVNINST/Institutes/LBNL/Luminosity/GlobalSettings/trunk GlobalSettings
svn co $SVNINST/Institutes/LBNL/Luminosity/PileupCorrections/trunk PileupCorrections
svn co $SVNINST/Institutes/LBNL/Luminosity/LumiVtx/trunk LumiVtx
svn co $SVNINST/Institutes/LBNL/Luminosity/D3PDData/trunk D3PDData
svn co $SVNINST/Institutes/LBNL/Luminosity/D3PDMC/trunk D3PDMC
#svn co svn+ssh://svn.cern.ch/reps/atlasinst/Institutes/LBNL/Luminosity/VtxLumiUnfold/trunk VtxLumiUnfold

#Deprecated packages
#svn co $SVNINST/Institutes/LBNL/Luminosity/FakeCorrection/trunk FakeCorrection
#svn co $SVNINST/Institutes/LBNL/Luminosity/PileupMaskingCorrection/trunk PileupMaskingCorrection
#svn co $SVNINST/Institutes/LBNL/Luminosity/VtxCalibration/trunk VtxCalibration



#cd RootCore
#./configure
#cd ../
#source RootCore/scripts/setup.sh
#
#${ROOTCOREBIN}/scripts/find_packages.sh
#${ROOTCOREBIN}/scripts/compile.sh
#
#echo "All Done."
#
