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

mash dist /projects/data/team1_genomeAssembly/reference_based_assembly/mash_sketches_samples/all.msh /projects/data/team1_genomeAssembly/reference_based_assembly/all_genomes.msh > /projects/data/team1_genomeAssembly/reference_based_assembly/ref_gen_vs_sample_all.dist

mash dist /projects/data/team1_genomeAssembly/reference_based_assembly/mash_sketches_samples/all.msh /projects/data/team1_genomeAssembly/reference_based_assembly/mash_sketches_samples/all.msh > /projects/data/team1_genomeAssembly/reference_based_assembly/samples_vs_samples.dist

cat "" > /projects/data/team1_genomeAssembly/reference_based_assembly/best_genomes_all.txt
for f in `ls $inp_fold | grep 'SRR' | cut -d '_' -f 1 | uniq`
do
    g=`cat /projects/data/team1_genomeAssembly/reference_based_assembly/ref_gen_vs_sample_all.dist | awk '/'${f}'/ {print}' | sort -k3 -n | head -n1 | awk '{print $2}'`
    echo ${f} ${g} >> /projects/data/team1_genomeAssembly/reference_based_assembly/best_genomes_all.txt
done
