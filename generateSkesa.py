import os
import shutil


fileDir = "/projects/data/team1_genomeAssembly/trimming2/fastq/trimmed"
#libFile = "/projects/data/team1_genomeAssembly/denovo_skesa/library"
libFile = "/projects/data/team1_genomeAssembly/denovo_skesa/sspaceLibrary"


# preparing a list of all the IDs
def generateGeneList():
	geneList = []
	for file in os.listdir(fileDir):
		if file.endswith(".fq"):
			a,b,c = file.split('_')
			geneList.append(a)
        	

	unique_geneList = (list(set(geneList)))
	return unique_geneList


# Running skesa (denovo assembly) command on all the forward and reverse paired files
def runSkesa(geneList):
	for a in geneList:
		fFile = '%s_forward_paired.fq' %(a)
		rFile = '%s_reverse_paired.fq' %(a)	
		forwardFile = os.path.join(fileDir,fFile)
		reverseFile = os.path.join(fileDir,rFile)
		#print (forwardFile,reverseFile)
		skesaCmd = 'skesa --fastq %s --fastq %s --contigs_out /projects/data/team1_genomeAssembly/denovo_skesa/skesaoutput/%s.skesa.fa' %(forwardFile,reverseFile,a)
		os.system(skesaCmd)

# Building a library with forward and reverse paired files, insert size for paired reads, error rate, ForwardReverse order, to run SSPACE
	
def generateLibFiles(geneList):
	for gene in geneList:
	
		libFileName  = '%s/%s'%(libFile,gene)
	
		libText="%s_lib /projects/data/team1_genomeAssembly/trimming2/fastq/trimmed/%s_forward_paired.fq /projects/data/team1_genomeAssembly/trimming2/fastq/trimmed/%s_reverse_paired.fq 250 0.75 FR" %(gene,gene,gene)
		if not os.path.exists(libFileName):
			with open(libFileName,'w') as fh:
				fh.write(libText)
				fh.close()
# Running the SSPACE command for scaffolding using default parameters and contig extension (-x 1)

		sspaceCmd = "perl /projects/data/team1_genomeAssembly/SSPACE/sspace_basic/SSPACE_Basic.pl -l /projects/data/team1_genomeAssembly/denovo_skesa/sspaceLibrary/%s -s /projects/data/team1_genomeAssembly/denovo_skesa/skesaoutput/%s.skesa.fa -x 1 -T 8 -b %s.sspace -m 20 -o 15 -a 0.8 -n 12 -g 3 -p 1" %(gene,gene,gene)


        	os.system(sspaceCmd)
		print("Done scaffolding")

# Remove all the folders created during scaffolding to minimize use of space

		shutil.rmtree('/projects/data/team1_genomeAssembly/denovo_skesa/reads')
		shutil.rmtree('/projects/data/team1_genomeAssembly/denovo_skesa/bowtieoutput')
		shutil.rmtree('/projects/data/team1_genomeAssembly/denovo_skesa/pairinfo')
		shutil.rmtree('/projects/data/team1_genomeAssembly/denovo_skesa/intermediate_results')
		shutil.rmtree('/projects/data/team1_genomeAssembly/denovo_skesa/dotfiles')


geneList = generateGeneList()
runSkesa(geneList)
generateLibFiles(geneList)
