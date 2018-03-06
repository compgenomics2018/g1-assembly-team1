#!/bin/bash

inp_fold=/projects/data/team1_genomeAssembly/trimming2/fastq/trimmed
outp_fold=/projects/data/team1_genomeAssembly/reference_based_assembly/mash_sketches_samples
ref_gen_path=/projects/data/team1_genomeAssembly/reference_based_assembly/reference_genomes

for f in `ls $inp_fold | grep 'SRR' | cut -d '_' -f 1 | uniq`
do
    echo $f
    cat ${inp_fold}/${f}* > ${outp_fold}/${f}.fastq
    mash sketch -m 2 ${outp_fold}/${f}.fastq
    rm ${outp_fold}/${f}.fastq
done

mash paste /projects/data/team1_genomeAssembly/reference_based_assembly/mash_sketches_samples/all ${outp_fold}/*

mash sketch ${ref_gen_path}/*.fna -o /projects/data/team1_genomeAssembly/reference_based_assembly/all_genomes

mash dist /projects/data/team1_genomeAssembly/reference_based_assembly/mash_sketches_samples/all.msh /projects/data/team1_genomeAssembly/reference_based_assembly/all_genomes.msh > /projects/data/team1_genomeAssembly/reference_based_assembly/ref_gen_vs_sample_50.dist

for f in `ls $inp_fold | grep 'SRR' | cut -d '_' -f 1 | uniq`
do
   # cat ref_gen_vs_sample.dist | awk '/SRR5666618/ {print}' | awk 'BEGIN{a=1000}{if ($3<0+a) a=$3} END{print $2}' | awk '{split($1,a,"/"); print a[3]}'
    g=`cat /projects/data/team1_genomeAssembly/reference_based_assembly/ref_gen_vs_sample_50.dist | awk '/'${f}'/ {print}' | sort -k3 -n | head -n1 | awk '{print $2}' | awk '{split($1,a,"/"); print a[7]}'`
    #g=`cat /projects/data/team1_genomeAssembly/reference_based_assembly/ref_gen_vs_sample_50.dist | awk '/'${f}'/ {print}' | sort -k3 -n | head -n1 | awk '{print $2}'`
    echo ${f} ${g} >> /projects/data/team1_genomeAssembly/reference_based_assembly/best_genomes_50.txt
    #echo $f "\t" cat ref_gen_vs_sample.dist | awk '/${f}/ {print}' | awk 'BEGIN{a=1000}{if ($3<0+a) a=$3} END{print $2}' | awk '{split($1,a,"/"); print a[3]}'
done

