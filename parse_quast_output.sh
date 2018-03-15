#!/bin/bash

inp_fold="/projects/data/team1_genomeAssembly/reference_based_assembly/quast_scaffolds"
outp_fold="/projects/data/team1_genomeAssembly/reference_based_assembly"

cat "" > ${outp_fold}"/ref_N50.txt"
cat "" > ${outp_fold}"/ref_number_of_contigs.txt"
cat "" > ${outp_fold}"/ref_largest_contig.txt"
cat "" > ${outp_fold}"/ref_Ns.txt"
cat "" > ${outp_fold}"/ref_total_length.txt"

for f in `ls $inp_fold`
do
    echo $f `cat ${inp_fold}"/"${f}"/"report.txt | grep 'N50' | awk '{print $2}'` >> ${outp_fold}"/ref_N50.txt"
    echo $f `cat ${inp_fold}"/"${f}"/"report.txt | grep "# contigs (>= 0 bp)" | awk '!/All statistics/' | awk '{print $6}'` >> ${outp_fold}"/ref_number_of_contigs.txt"
    echo $f `cat ${inp_fold}"/"${f}"/"report.txt | grep 'Largest contig' | awk '{print $3}'` >> ${outp_fold}"/ref_largest_contig.txt"
    echo $f `cat ${inp_fold}"/"${f}"/"report.txt | awk '/s per 100 kbp/ {print}' | awk '{print $6}'` >> ${outp_fold}"/ref_Ns.txt"
    echo $f `cat ${inp_fold}"/"${f}"/"report.txt | awk '/Total length/ {print}' | head -3 | tail -1 | awk '{print $6}'` >> ${outp_fold}"/ref_total_length.txt"
done

