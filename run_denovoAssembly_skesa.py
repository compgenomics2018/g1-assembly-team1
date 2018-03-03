import os

sspaceFolder = "/projects/data/team1_genomeAssembly/denovo_skesa/skesaoutput"

for file in os.listdir(sspaceFolder):

        if file.endswith ('.skesa.fa'):
                quastCmd = "quast.py -t 6 -m 1000 /projects/data/team1_genomeAssembly/denovo_skesa/skesaoutput/%s -o /projects/data/team1_genomeAssembly/denovo_skesa/quastOutput2/%s" %(file,file)
                os.system(quastCmd)


