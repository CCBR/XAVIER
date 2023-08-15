#! /bin/zsh
# source https://www.biostars.org/p/141256/#142145

if (( $# == 2 )) ; then

    paste - - - - - - - - \
    | tee \
        >(cut -f 1-4 \
        | tr "\t" "\n" \
        > $1) \
    | cut -f 5-8 \
    | tr "\t" "\n" \
    > $2

else

    echo -e "\
${0##*/}
GY140920

Deinterleaves a FASTQ file of paired reads

Usage:

    cat <in.fastq> | deinterleave.sh <out1.fq> <out2.fq>

    zcat <in.fq.gz> | deinterleave.sh <out1.fq> <out2.fq>

    samtools bam2fq -ns singles.fq <in.bam> \\
    | deinterleave.sh <out1.fq> <out2.fq>"

fi
