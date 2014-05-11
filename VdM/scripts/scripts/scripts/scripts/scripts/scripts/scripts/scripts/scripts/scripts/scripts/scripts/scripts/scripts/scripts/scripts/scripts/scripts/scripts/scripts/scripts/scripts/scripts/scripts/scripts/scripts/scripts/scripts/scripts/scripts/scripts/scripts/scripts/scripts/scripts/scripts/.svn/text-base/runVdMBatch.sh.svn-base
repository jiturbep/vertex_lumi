#!/bin/bash

EXTRAOPT="-n -v"

if [ -z "$1" ]; then
    echo "usage: $0 run settings"
    exit
else 
    RUN=$1
fi

if [ -z "$2" ]; then
    echo "usage: $0 run settings"
    exit
else 
    SETTINGS=$2
fi

if [ "$3" == "NEvt" ]; then
	VTX_METHOD="NEvt"
	EXTRAOPT="-e ${EXTRAOPT}"
else
	VTX_METHOD="NVtx"
fi

COMMAND="runVdM -r ${RUN} -s ${SETTINGS} ${EXTRAOPT} >& /u/dryu/Luminosity/Data/VdMCalibration/${RUN}/${SETTINGS}/log_${VTX_METHOD}.txt"

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
qsub -V -l h_vmem=2G -l eliza18io=1 tmp.sh
rm tmp.sh
