#!/bin/sh

#Put somthing here for test-subfolders to be created, otherwise leave ""

VERSION="v1"; #EXTRAOPT="-v -e 2000000" #" -e 1000000" 
EXTRAOPT=""

EXEPATH=./

if [ -z "$1" ]; then
    echo "usage: $0 energy settings option"
    exit
else 
    ENERGY=$1
fi

if [ -z "$2" ]; then
    echo "usage: $0 energy settings option"
    exit
else
    SETTINGS=$2
fi

if [ "${ENERGY}" == "7" -a "${SETTINGS}" == "16.X-tight" ]; then
	INPUT_FILE="/eliza18/atlas/dryu/Luminosity/mc10_7TeV/user.spagan.mc10_7TeV.105000.pythia_minbias_inelastic.merge.VTXD3PD.e723_s932_s946_r2302_r2300.v6.0/*root*"
elif [ "${ENERGY}" == "7" -a "${SETTINGS}" == "17.2-normal" ]; then
    #INPUT_FILE="/eliza18/atlas/data/data12_8TeV/VtxNtuples_17.2.1.3/user.beate.trkReco_default.mc11_7TeV.119994.Pythia8_A2M_MSTW2008LO_minbias_Inelastic.recon.ESD.e1051_s1372_s1370_r3157.v4.120418172907_D3PD/*root*"
    INPUT_FILE="/eliza18/atlas/data/data12_8TeV/VtxNtuples_17.2.1.3/user.beate.trkReco_default.mc11_7TeV.105000.pythia_minbias_inelastic.recon.ESD.e816_s1310_s1300_r3043.v2.120418111618_D3PD/*root*"
elif [ "${ENERGY}" == "7" -a "${SETTINGS}" == "17.2-VtxLumi" ]; then
    #INPUT_FILE="/eliza18/atlas/data/data12_8TeV/VtxNtuples_17.2.1.3/user.beate.trkReco_VtxLumi.mc11_7TeV.119994.Pythia8_A2M_MSTW2008LO_minbias_Inelastic.recon.ESD.e1051_s1372_s1370_r3157.v4.120418173257_D3PD/*root*"
    #INPUT_FILE="/eliza18/atlas/data/data12_8TeV/VtxNtuples_17.2.1.3/user.beate.trkReco_VtxLumi.mc11_7TeV.105000.pythia_minbias_inelastic.recon.ESD.e816_s1310_s1300_r3043.v2.120418113020_D3PD/*root*"
    INPUT_FILE="/eliza18/atlas/data/data12_8TeV/VtxNtuples_17.2.1.3/user.beate.trkReco_VtxLumi.mc11_7TeV.119994.Pythia8_A2M_MSTW2008LO_minbias_Inelastic.recon.ESD.e1051_s1372_s1370_r3157.v5.120424093604_D3PD/user.beate.002426.D3PD._00034.root"
elif [ "${ENERGY}" == "8" -a "${SETTINGS}" == "17.2-normal" ]; then
    INPUT_FILE="/eliza2/atlas/atlasdata/atlasscratchdisk/user/spagan/valid1/user.spagan.valid1.107499.singlepart_empty.recon.VTXD3PD.e603_s1469_r3560.v1.0.120419110055/*root*"
    EXTRAOPT="-a ${EXTRAOPT}"
else
	echo "Sample does not exist for energy ${ENERGY} and settings ${SETTINGS}"
	exit
fi

SUBMIT_BATCH="false"

if [ "$3" == "test" ]; then
    echo "Test run"
    TESTSUBDIR="test"    
    EXTRAOPT="-p -v -e 10000 ${EXTRAOPT}"
elif [ "$3" == "batch" ]; then
    SUBMIT_BATCH="true"
else
    TESTSUBDIR="${VERSION}"
fi


RESULT_DIR="/eliza18/atlas/dryu/Luminosity/VertexCounts/mc_${ENERGY}TeV/${SETTINGS}"

if ! [ -z "${TESTSUBDIR=}" ]; then
    RESULT_DIR="${RESULT_DIR}/${TESTSUBDIR}"
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
