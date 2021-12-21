#!/usr/bin/bash

#Paths to execs
PREFETCH=prefetch
PFD=parallel-fastq-dump

NCORES=40

OUTDIR=/workdir/dwm269/ma_ipsc/data/fastqs/Silva_CSC
mkdir -p ${OUTDIR}
cd ${OUTDIR}


# ${PREFETCH} \
# --verify yes \
# --max-size 999999999999 \
# --output-directory ${OUTDIR} \
# ${OUTDIR}/SRR12075831.sra

# ${PFD} \
# --sra-id SRR12075831 \
# --threads ${NCORES} \
# --outdir ${OUTDIR} \
# --tmpdir ${OUTDIR} \
# --split-files \
# --gzip


${PREFETCH} \
--verify yes \
--max-size 999999999999 \
--output-directory ${OUTDIR} \
${OUTDIR}/SRR12075832.sra

${PFD} \
--sra-id SRR12075832 \
--threads ${NCORES} \
--outdir ${OUTDIR} \
--tmpdir ${OUTDIR} \
--split-files \
--gzip
