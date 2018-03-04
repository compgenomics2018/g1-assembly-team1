#!/bin/bash

#Reference+denovo+Re-scaffolded with unmapped reads at 10Ns

referencebasedassembly=/projects/data/team1_genomeAssembly/reference_based_assembly/assembly_all

for a in $(cat /projects/home/nshah377/IDs.txt)
do
        echo "$a"
        mash sketch /projects/data/team1_genomeAssembly/reference_based_assembly/assembly_all/$a/ref_10_and_denovo.final.scaffolds.fasta \
        -o /projects/data/team1_genomeAssembly/MASH/refdenovo/$a.refdenovo
done

#SKESA
referencebasedassembly=/projects/data/team1_genomeAssembly/reference_based_assembly/assembly_all

for a in $(cat /projects/home/nshah377/IDs.txt)
do
        echo "$a"
        mash sketch /projects/data/team1_genomeAssembly/denovo_skesa/sspaceOutput/$a.sspace.final.scaffolds.fasta \
        -o /projects/data/team1_genomeAssembly/MASH/skesa/$a.skesa
done

#SPADES
for a in $(cat /projects/home/nshah377/IDs.txt)
do
        for b in $(ls /projects/data/team1_genomeAssembly/transfer/$a/output/fasta)
        do
                echo "$a,$b"
                mash sketch /projects/data/team1_genomeAssembly/transfer/$a/output/fasta/$b \
                -o /projects/data/team1_genomeAssembly/MASH/spades/$a.$b
        done
done
