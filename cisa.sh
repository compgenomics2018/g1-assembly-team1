#!/bin/bash

ref_fold=/projects/data/team1_genomeAssembly/reference_based_assembly/assembly_all
denovo_fold=/projects/data/team1_genomeAssembly/denovo_skesa/sspaceOutput
CISA_fold=/projects/home/nshah377/merging/CISAconfigs

mkdir -p $CISA_fold
echo "CISA folder created"
#a=SRR5666627
for a in $(cat /projects/home/nshah377/IDs.txt)
do
        echo "####################################################################################################################################"
        echo "$a"
        echo "Performing merging"
#Making merge.config files - output will be merge.<filename>.config
        mkdir -p /projects/home/nshah377/merging/output/$a
        cd /projects/home/nshah377/merging/output/$a
        echo "merge.$a.config"
        echo "count=2" > /projects/home/nshah377/merging/CISAconfigs/merge.$a.config
        echo "data=/projects/data/team1_genomeAssembly/reference_based_assembly/assembly_all/$a/ref_10_and_denovo.final.scaffolds.fasta,title=refNscaffolded" >> /projects/home/nshah377/merging/CISAconfigs/merge.$a.config
        echo "data=/projects/data/team1_genomeAssembly/denovo_skesa/sspaceOutput/$a.sspace.final.scaffolds.fasta,title=skesa" >> /projects/home/nshah377/merging/CISAconfigs/merge.$a.config
        echo "Master_file=/projects/home/nshah377/merging/output/$a/$a.contigs_m.fa" >> /projects/home/nshah377/merging/CISAconfigs/merge.$a.config
        echo "min_length=100" >> /projects/home/nshah377/merging/CISAconfigs/merge.$a.config
        echo "Gap=11" >> /projects/home/nshah377/merging/CISAconfigs/merge.$a.config

#Making cisa.config files - output will be cisa.$a.config
        echo "$a"
        echo "Performing CISA configuration"

        echo "cisa.$a.config"
        echo "genome=5629435" > /projects/home/nshah377/merging/CISAconfigs/cisa.$a.config
        echo "infile=/projects/home/nshah377/merging/output/$a/$a.contigs_m.fa" >> /projects/home/nshah377/merging/CISAconfigs/cisa.$a.config
        echo "outfile=/projects/home/nshah377/merging/output/$a/$a.mergedcontigs.fa" >> /projects/home/nshah377/merging/CISAconfigs/cisa.$a.config
        echo "nucmer=/projects/data/team1_genomeAssembly/bin/MUMmer3.23/nucmer" >> /projects/home/nshah377/merging/CISAconfigs/cisa.$a.config
        echo "R2_Gap=0.95" >> /projects/home/nshah377/merging/CISAconfigs/cisa.$a.config
        echo "CISA=/projects/data/team1_genomeAssembly/bin/CISA1.3" >> /projects/home/nshah377/merging/CISAconfigs/cisa.$a.config
        echo "makeblastdb=/projects/data/team1_genomeAssembly/bin/makeblastdb" >> /projects/home/nshah377/merging/CISAconfigs/cisa.$a.config
        echo "blastn=/projects/data/team1_genomeAssembly/bin/blastn" >> /projects/home/nshah377/merging/CISAconfigs/cisa.$a.config

        Merge.py /projects/home/nshah377/merging/CISAconfigs/merge.$a.config
        time yes "y" |CISA.py /projects/home/nshah377/merging/CISAconfigs/cisa.$a.config
        echo "######################################################################################################################################"

done
