#!/bin/sh

#Put somthing here for test-subfolders to be created, otherwise leave ""

VERSION="v9"; #EXTRAOPT="-v -e 2000000" #" -e 1000000" 
EXTRAOPT=""

EXEPATH=./

if [ -z "$1" ]; then
    echo "usage: $0 run_number settings option"
    exit
else 
    RUN=$1
fi

if [ -z "$2" ]; then
    echo "usage: $0 run_number settings option"
    exit
else
    SETTINGS=$2
fi

if ! [ -d input/$RUN/$SETTINGS ]; then
    echo "Directory input/${RUN}/${SETTINGS} does not exist. Exiting."
    exit
fi

SUBMIT_BATCH="false"

if [ "$3" == "test" ]; then
    echo "Test run"
    TESTSUBDIR="test"    
    EXTRAOPT="-v -e 2000000"
elif [ "$3" == "batch" ]; then
    SUBMIT_BATCH="true"
else
    TESTSUBDIR="${VERSION}"
fi

RESULT_DIR="/eliza18/atlas/dryu/Luminosity/VertexCounts/${RUN}/${SETTINGS}"

INPUT_FILE="input/${RUN}/${SETTINGS}/*root*"


if ! [ -z "${TESTSUBDIR=}" ]; then
    RESULT_DIR="${RESULT_DIR}/${TESTSUBDIR}"
fi
RESULT_FILE="${RESULT_DIR}/InDetTrackD3PD_results"

COMMAND="${EXEPATH}/TrackingAna.exe -b -q ${EXTRAOPT} -t L1_RDO_FILLED -o ${RESULT_FILE} ${INPUT_FILE}"

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
    qsub -V tmp.sh
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
