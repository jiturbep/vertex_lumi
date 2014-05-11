#!/bin/sh

#Put somthing here for test-subfolders to be created, otherwise leave ""

VERSION="v2"; #EXTRAOPT="-v -e 2000000" #" -e 1000000" 
EXTRAOPT=""

EXEPATH=./

if [ -z "$1" ]; then
    echo "usage: $0 tag option"
    exit
else 
    TAG=$1
fi


SUBMIT_BATCH="false"

if [ "$2" == "test" ]; then
    echo "Test run"
    TESTSUBDIR="test"    
    EXTRAOPT="-p -v -e 10000 ${EXTRAOPT}"
elif [ "$2" == "batch" ]; then
    SUBMIT_BATCH="true"
elif [ "$2" == "local" ]; then
    TESTSUBDIR="${VERSION}"
else
    echo "usage: $0 tag (test | local | batch) ( | chi2ndf)"
fi

MAX_CHI2NDF="1.5"
if [ "$3" == "chi2ndf" ]; then
    echo "Setting max chi2/ndf = ${MAX_CHI2NDF}"
    EXTRAOPT="-c "

RESULT_DIR="/eliza18/atlas/dryu/Luminosity/VertexCounts/MC"


case ${TAG} in        
    "CT_7TeV_ND")
        INPUT_FILE="/eliza2/atlas/atlasdata/atlaslocalgroupdisk/user/spagan/valid1/user.spagan.valid1.105001.pythia_minbias.recon.VTXD3PD.e574_s1194_r3376.v1.0.120330145127/*root*"
        RESULT_DIR="${RESULT_DIR}/mc_7TeV_17.2_default_pythia6_singleevent/ND";;
    "CT_7TeV_DD")
        INPUT_FILE="/eliza2/atlas/atlasdata/atlaslocalgroupdisk/user/spagan/valid1/user.spagan.valid1.105004.pythia_ddiff.recon.VTXD3PD.e574_s1194_r3376.v1.0.120330144648/*root*"
        RESULT_DIR="${RESULT_DIR}/mc_7TeV_17.2_default_pythia6_singleevent/DD";;
    "CT_7TeV_SD")
        INPUT_FILE="/eliza2/atlas/atlasdata/atlaslocalgroupdisk/user/spagan/valid1/user.spagan.valid1.105003.pythia_sdiff.recon.VTXD3PD.e574_s1194_r3376.v1.0.120330144904/*root*"
        RESULT_DIR="${RESULT_DIR}/mc_7TeV_17.2_default_pythia6_singleevent/SD";;
    "CT_8TeV_ND")
        INPUT_FILE="/eliza2/atlas/atlasdata/atlaslocalgroupdisk/user/spagan/valid1/user.spagan.valid1.105001.pythia_minbias.recon.VTXD3PD.e671_s1194_r3376.v1.0.120330145852/*root*"
        RESULT_DIR="${RESULT_DIR}/mc_8TeV_17.2_default_pythia6_singleevent/ND";;
    "CT_8TeV_DD")
        INPUT_FILE="/eliza2/atlas/atlasdata/atlaslocalgroupdisk/user/spagan/valid1/user.spagan.valid1.105004.pythia_ddiff.recon.VTXD3PD.e671_s1194_r3376.v1.0.120330145400/*root*"
        RESULT_DIR="${RESULT_DIR}/mc_8TeV_17.2_default_pythia6_singleevent/DD";;
    "CT_8TeV_SD")
        INPUT_FILE="/eliza2/atlas/atlasdata/atlaslocalgroupdisk/user/spagan/valid1/user.spagan.valid1.105003.pythia_sdiff.recon.VTXD3PD.e671_s1194_r3376.v1.0.120330145624/*root* /eliza2/atlas/atlasdata/atlaslocalgroupdisk/user/spagan/valid1/user.spagan.valid1.105003.pythia_sdiff.recon.VTXD3PD.e671_s1194_r3376.v1.0.120330145624_r1/*root*"
        RESULT_DIR="${RESULT_DIR}/mc_8TeV_17.2_default_pythia6_singleevent/SD";;
    "CT_7TeV_MB_pythia6")
        INPUT_FILE="/eliza18/atlas/data/data12_8TeV/VtxNtuples_17.2.1.3/user.beate.trkReco_default.mc11_7TeV.105000.pythia_minbias_inelastic.recon.ESD.e816_s1310_s1300_r3043.v5.120424111716_D3PD/*root*"
        RESULT_DIR="${RESULT_DIR}/mc_8TeV_17.2_default_pythia6_highpu/ND";;
    "CT_7TeV_MB_pythia8")
        INPUT_FILE="/eliza18/atlas/data/data12_8TeV/VtxNtuples_17.2.1.3/user.beate.trkReco_default.mc11_7TeV.119994.Pythia8_A2M_MSTW2008LO_minbias_Inelastic.recon.ESD.e1051_s1372_s1370_r3157.v5.120424111529_D3PD/*root*"
        RESULT_DIR="${RESULT_DIR}/mc_7TeV_17.2_default_pythia8_pu";;
    "CT_8TeV_ND_pythia8")
        INPUT_FILE="/eliza2/atlas/atlasdata/atlaslocalgroupdisk/user/spagan/mc12_8TeV/user.spagan.mc12_8TeV.119997.Pythia8_A2MSTW2008LO_minbias_ND.merge.VTXD3PD.e1119_s1468_s1470_r3498_r3493.v1.0.120408174013/*root*"
        RESULT_DIR="${RESULT_DIR}/mc_8TeV_17.2_default_pythia8_singleevent/ND";;
    "CT_8TeV_DD_pythia8")
        INPUT_FILE="/eliza2/atlas/atlasdata/atlaslocalgroupdisk/user/spagan/mc12_8TeV/user.spagan.mc12_8TeV.119999.Pythia8_A2MSTW2008LO_minbias_DD.merge.VTXD3PD.e1119_s1468_s1470_r3498_r3493.v1.0.120411140915/*root*"
        RESULT_DIR="${RESULT_DIR}/mc_8TeV_17.2_default_pythia8_singleevent/DD";;
    "CT_8TeV_SD_pythia8")
        INPUT_FILE="/eliza2/atlas/atlasdata/atlaslocalgroupdisk/user/spagan/mc12_8TeV/user.spagan.mc12_8TeV.119998.Pythia8_A2MSTW2008LO_minbias_SD.merge.VTXD3PD.e1119_s1468_s1470_r3498_r3493.v1.0.120411140832/*root*"
        RESULT_DIR="${RESULT_DIR}/mc_8TeV_17.2_default_pythia8_singleevent/SD";;
    "mc_7TeV_16.X_hybrid_pythia6_pu")
        #INPUT_FILE="/eliza2/atlas/atlasdata/atlaslocalgroupdisk/user/spagan/mc10_7TeV/user.spagan.mc10_7TeV.105000.pythia_minbias_inelastic.merge.VTXD3PD.e723_s932_s946_r2302_r2300.v4.0.120401064511/*root*"
        INPUT_FILE="/eliza18/atlas/dryu/Luminosity/mc10_7TeV/user.spagan.mc10_7TeV.105000.pythia_minbias_inelastic.merge.VTXD3PD.e723_s932_s946_r2302_r2300.v6.0/*root*"
        RESULT_DIR="${RESULT_DIR}/mc_7TeV_16.X_hybrid_pythia6_pu";;
    "mc_7TeV_16.X_default_pythia6_pu")
        #Fake correction for old 7 TeV analysis.
        INPUT_FILE="/eliza18/atlas/spagan/VtxStudies/VtxPU/VTX_MON/mc10b/user.spagan.mc10_7TeV.105000.pythia_minbias_inelastic.merge.VTX_MON.e723_s932_s946_r2302_r2300.v1.1.120217055812/*root*"
        RESULT_DIR="${RESULT_DIR}/${TAG}";;
    "mc_7TeV_17.2_default_pythia6_pu")
        INPUT_FILE="/eliza18/atlas/spagan/VtxStudies/VtxPU/VTXD3PD/mc_Beate/MC11_7TeV.107499.singlepart_empty.pileup_pythia6*root*"
        RESULT_DIR="${RESULT_DIR}/${TAG}"
        EXTRAOPT="-a ${EXTRAOPT}";;
    "mc_7TeV_17.2_VtxLumi_pythia6_pu")
        echo "You don't have a sample for 7 TeV, 17.2 VtxLumi, pythia 6 yet!"
        exit;;
    "mc_7TeV_17.2_default_pythia8_pu")
        INPUT_FILE="/eliza18/atlas/spagan/VtxStudies/VtxPU/VTXD3PD/mc_Beate/MC11_7TeV.107499.singlepart_empty.pileup_Pythia8*root*"
        RESULT_DIR="${RESULT_DIR}/${TAG}"
        EXTRAOPT="-a ${EXTRAOPT}";;
    "mc_7TeV_17.2_VtxLumi_pythia8_pu")
        echo "You don't have a sample for 7 TeV, 17.2 VtxLumi, pythia 8 yet!"
        exit;;
    "mc_8TeV_17.2_default_pythia8_pu")
        INPUT_FILE="/eliza2/atlas/atlasdata/atlasscratchdisk/user/spagan/valid1/user.spagan.valid1.107499.singlepart_empty.recon.VTXD3PD.e603_s1469_r3560.v1.0.120419110055/*root*"
        RESULT_DIR="${RESULT_DIR}/${TAG}"
        EXTRAOPT="-a ${EXTRAOPT}";;
    "mc_8TeV_17.2_VtxLumi_pythia8_pu")
        INPUT_FILE="/eliza18/atlas/data/data12_8TeV/VtxNtuples_17.2.1.3/user.beate.trkReco_default.valid1.107499.singlepart_empty.recon.ESD.e603_s1469_r3560.v11.120430123934_D3PD/*root*"
        RESULT_DIR="${RESULT_DIR}/${TAG}"
        EXTRAOPT="-a ${EXTRAOPT}";;
    *)
        echo "Unknown tag: ${TAG}"
        exit;;
esac


if ! [ -z "${TESTSUBDIR}" ]; then
    RESULT_DIR="${RESULT_DIR}/${TESTSUBDIR}"
elif ! [ -z "${VERSION}" ]; then
    RESULT_DIR="${RESULT_DIR}/${VERSION}"
fi

RESULT_FILE="${RESULT_DIR}/InDetTrackD3PD_results"

COMMAND="${EXEPATH}/TrackingAna.exe -b -q ${EXTRAOPT} -o ${RESULT_FILE} ${INPUT_FILE}"

if [ "${SUBMIT_BATCH}" == "true" ]; then
    #Make batch submission script
    touch tmp.sh
    echo "#!/bin/sh" >> tmp.sh
    echo "${COMMAND}" >> tmp.sh
    echo "Batch submission:"
    cat tmp.sh

    echo ">>> Press any key to run. 'q' to exit."
    read ans
    if [ "${ans}" == "q" ] || [ "${ans}" == "Q" ]; then
        echo "Aborting"
        rm tmp.sh
        exit
    fi

    echo "Submitting batch job!"
    mkdir -pv ${RESULT_DIR}
    qsub -V -l h_vmem=2G tmp.sh
    rm tmp.sh

else
    LOGFILE="${RESULT_DIR}/InDetTrackD3PD_TrackingAna.log"
    echo "COMMAND: ${COMMAND}"
    echo "LOG: ${LOGFILE}"
    echo ">>> Press any key to run. 'q' to exit."
    read ans
    if [ "${ans}" == "q" ] || [ "${ans}" == "Q" ]; then
        echo "Aborting"
        exit
    fi
    echo "Running!"
    mkdir -pv ${RESULT_DIR}
    ${COMMAND}  >& ${LOGFILE} &

    COMMAND="ps aux | grep TrackingAna.exe | grep -v 'grep' | awk '{print \$1,\$2,\$3,\$4,\$11,\$15}'"
    echo ${COMMAND}
fi
