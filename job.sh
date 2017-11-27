#!/bin/bash

set -e

export PATH=/usr/local/sratoolkit/current/bin:/usr/local/bowtie2/current:/usr/local/RAPSearch2/bin:/usr/local/bin:/usr/bin

REF_BASENAME=$1
SRA_ID=$2

ODIR=bam.$REF_BASENAME
mkdir -p $ODIR

DIR=fastq$$
mkdir -p $DIR

echo
echo
pwd
ls -al
echo
echo

{
    # check wrangler cache first
    WRANGLER_LOC=DISABLED/nas/wrangler/NCBI/SRA/Downloads/sra/$SRA_ID.sra
    if [ -e $WRANGLER_LOC ]; then
        cp $WRANGLER_LOC $DIR/
    else
        # fall back - download directly
        fastq-dump --outdir $DIR --skip-technical  --readids --read-filter pass --dumpbase --split-files --clip $SRA_ID
    fi

    READS=$(ls $DIR/* | tr \\n \, | sed -e 's/,$//')
    bowtie2 -p 6 -q --no-unal -x $REF_BASENAME -U $READS | samtools view -bS - | samtools sort - $ODIR/$SRA_ID

    samtools index $ODIR/$SRA_ID.bam

} 2>&1

rm -rf $DIR
rm -f $HOME/ncbi/public/sra/$SRA_ID.sra*

