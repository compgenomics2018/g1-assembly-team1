#!/bin/bash

ref_fold=/projects/home/nshah377/refbasedassembly
CISA_fold=/projects/home/nshah377/merging/CISAconfigs

mkdir -p $CISA_fold
echo "CISA folder created"

for a in $(cat /projects/home/nshah377/IDs.txt)
do
#Making merge.config files - output will be merge.<filename>.config
        echo "merge.$a.config"
        echo "count=2" > /projects/home/nshah377/merging/CISAconfigs/merge.$a.config
        echo "data=/projects/data/team1_genomeAssembly/reference_based_assembly/assembly_all/$a/sspace_insert_500.final.scaffolds.fasta" >> /projects/home/nshah377/merging/CISAconfigs/merge.$a.config
        echo "data=/projects/data/team1_genomeAssembly/denovo_skesa/skesaoutput/$a.skesa.fa" >> /projects/home/nshah377/merging/CISAconfigs/merge.$a.config
        echo "Master_file=/projects/home/nshah377/merging/$a.contigs_m.fa" >> /projects/home/nshah377/merging/CISAconfigs/merge.$a.config
        echo "min_length=100" >> /projects/home/nshah377/merging/CISAconfigs/merge.$a.config
        echo "Gap=11" >> /projects/home/nshah377/merging/CISAconfigs/merge.$a.config

#Making cisa.config files - output will be cisa.$a.config
        echo "cisa.$a.config"
        echo "genome=5629435" > /projects/home/nshah377/merging/CISAconfigs/cisa.$a.config
        echo "infile=/projects/home/nshah377/merging/$a.contigs_m.fa" >> /projects/home/nshah377/merging/CISAconfigs/cisa.$a.config
        echo "outfile=/projects/home/nshah377/merging/$a.mergedcontigs.fa" >> /projects/home/nshah377/merging/CISAconfigs/cisa.$a.config
        echo "nucmer=/projects/data/team1_genomeAssembly/bin/MUMmer3.23/nucmer" >> /projects/home/nshah377/merging/CISAconfigs/cisa.$a.config
        echo "R2_Gap=0.95" >> /projects/home/nshah377/merging/CISAconfigs/cisa.$a.config
        echo "CISA=/projects/data/team1_genomeAssembly/bin/CISA1.3" >> /projects/home/nshah377/merging/CISAconfigs/cisa.$a.config
        echo "makeblastdb=/projects/data/team1_genomeAssembly/bin/makeblastdb" >> /projects/home/nshah377/merging/CISAconfigs/cisa.$a.config
        echo "blastn=/projects/data/team1_genomeAssembly/bin/blastn" >> /projects/home/nshah377/merging/CISAconfigs/cisa.$a.config

#       Merge.py /projects/home/nshah377/merging/CISAconfigs/merge.$a.config
#       CISA.py /projects/home/nshah377/merging/CISAconfigs/cisa.$a.config

done
