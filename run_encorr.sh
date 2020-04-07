#!/bin/bash

INDIR=$1
OUTDIR=$2
DATASET=$3
N=$4
SAMPLING_RATE=$5
TETINFO="${INDIR}/${DATASET}_tet_info.mat"
CCG_PATH="${OUTDIR}/CCG"; mkdir -p ${CCG_PATH}
CCF_PATH="${OUTDIR}/CCF"; mkdir -p ${CCF_PATH}
CORR_STAT_PATH="${OUTDIR}/corr_stat"; mkdir -p ${CORR_STAT_PATH}
CONN_STAT_PATH="${OUTDIR}/conn_stat"; mkdir -p ${CONN_STAT_PATH}
CORR_PATH="${OUTDIR}/correlograms"; mkdir -p ${CORR_PATH}

correlate=1
if [ "$correlate" = 1 ] ; then
    for i in $(seq 1 ${N})
    do
        for j in $(seq ${i} ${N})
        do 
            python3 encorr.py correlate "${INDIR}/${DATASET}_tet${i}.mat" \
                                        "${INDIR}/${DATASET}_tet${j}.mat" \
                                        "${INDIR}/${DATASET}_params.mat" \
                                        ${i} ${j} \
                                        ${SAMPLING_RATE} \
                                        "${CCG_PATH}/tet_${i}_tet_${j}.ccg"
        done
    done
fi

call=1
if [ "$call" = 1 ] ; then
    for i in $(seq 1 ${N})
    do
        for j in $(seq ${i} ${N})
        do 
            python3 encorr.py call "${CCG_PATH}/tet_${i}_tet_${j}.ccg" "${CCF_PATH}/tet_${i}_tet_${j}.ccf"
        done
    done
fi

corr_stat=1
if [ "$corr_stat" = 1 ] ; then
    python3 encorr.py corr-stat ${CCF_PATH} ${CORR_STAT_PATH}
fi

correlogram=1
if [ "$correlogram" = 1 ] ; then
    for i in $(seq 1 ${N})
    do
        for j in $(seq ${i} ${N})
        do 
            python3 encorr.py correlogram "${CCG_PATH}/tet_${i}_tet_${j}.ccg" "${CCF_PATH}/tet_${i}_tet_${j}.ccf" ${CORR_PATH}
        done
    done
fi

conn_stat=1
if [ "$conn_stat" = 1 ] ; then
    python3 encorr.py conn-stat ${CCF_PATH} ${TETINFO} ${CONN_STAT_PATH}/${DATASET}
fi
