#!/bin/bash

inp_fold=/projects/data/team1_genomeAssembly/trimming2/fastq/trimmed
outp_fold=/projects/data/team1_genomeAssembly/reference_based_assembly/assembly_all

cat '' > /projects/data/team1_genomeAssembly/reference_based_assembly/N50_summary.txt
for id in `ls $inp_fold | grep 'SRR' | cut -d '_' -f 1 | uniq`
do
    echo ${id} `cat ${outp_fold}/${id}/ref_10_and_denovo.summaryfile.txt  | grep 'After scaffolding' -A 8 | grep 'N50' | awk '{print $3}'` >> /projects/data/team1_genomeAssembly/reference_based_assembly/N50_summary.txt
done

cat '' > /projects/data/team1_genomeAssembly/reference_based_assembly/number_of_scaffolds_summary.txt
for id in `ls $inp_fold | grep 'SRR' | cut -d '_' -f 1 | uniq`
do
    echo ${id} `cat ${outp_fold}/${id}/ref_10_and_denovo.summaryfile.txt  | grep 'After scaffolding' -A 8 | grep 'Total number of scaffolds' | awk '{print $6}'` >> /projects/data/team1_genomeAssembly/reference_based_assembly/number_of_scaffolds_summary.txt
done

cat '' > /projects/data/team1_genomeAssembly/reference_based_assembly/number_of_Ns_summary.txt
for id in `ls $inp_fold | grep 'SRR' | cut -d '_' -f 1 | uniq`
do
    echo ${id} `cat ${outp_fold}/${id}/ref_10_and_denovo.summaryfile.txt  | grep 'After scaffolding' -A 8 | grep 'Total number of N' | awk '{print $6}'` >> /projects/data/team1_genomeAssembly/reference_based_assembly/number_of_Ns_summary.txt
done

cat '' > /projects/data/team1_genomeAssembly/reference_based_assembly/max_scaffold_length_summary.txt
for id in `ls $inp_fold | grep 'SRR' | cut -d '_' -f 1 | uniq`
do
    echo ${id} `cat ${outp_fold}/${id}/ref_10_and_denovo.summaryfile.txt  | grep 'After scaffolding' -A 8 | grep 'Max scaffold size' | awk '{print $5}'` >> /projects/data/team1_genomeAssembly/reference_based_assembly/max_scaffold_length_summary.txt
done

