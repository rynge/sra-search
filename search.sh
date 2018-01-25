#!/bin/bash

set -e

export PATH=/usr/local/sratoolkit/current/bin:/usr/local/bowtie2/current:/usr/local/RAPSearch2/bin:/usr/local/bin:/usr/bin

REF_BASENAME=$1
SRA_ID=$2

DIR=fastq$$
mkdir -p $DIR

# ensure the wrangler filesystem is mounted
if [ ! -e /nas/wrangler/NCBI/SRA/Downloads/sra ]; then
    echo "ERROR: Wrangler mount is not available!" 1>&2
    exit 1
fi

{
   echo

    # check wrangler cache first
    WRANGLER_LOC=/nas/wrangler/NCBI/SRA/Downloads/sra/$SRA_ID.sra
    if [ -e $WRANGLER_LOC ]; then
        SRA_SOURCE="$WRANGLER_LOC"
        echo "Will read $SRA_ID from $WRANGLER_LOC"
    else
        # default is download from SRA
        SRA_SOURCE="$SRA_ID"
        echo "WARNING: $SRA_ID not found on Wrangler - skipping..."
        touch $SRA_ID.bam $SRA_ID.bam.bai
        exit 0 
    fi

    fastq-dump --outdir $DIR --skip-technical  --readids --read-filter pass --dumpbase --split-files --clip $SRA_SOURCE

    READS=$(ls $DIR/* | tr \\n \, | sed -e 's/,$//')
    bowtie2 -p 2 -q --no-unal -x $REF_BASENAME -U $READS | samtools view -bS - | samtools sort - $SRA_ID

    samtools index $SRA_ID.bam

    echo

} 2>&1

rm -rf $DIR
rm -f $HOME/ncbi/public/sra/$SRA_ID.sra*

