List of scripts and supplementary files:
IDs.txt - file with SRR IDs of samples
adapters.fa - file with adapters, that we used for trimming
download_ref_genomes.sh - script for upload of complete reference genomes of Klebsiella from ncbi
generate_mash_sketch.sh - script for calculation of distances with MASH; it compares samples vs. samples and samples vs. reference genomes; it also chooses the best reference genome for each sample
parse_quast_output.sh - script for extraction of QC statistics from Quast output files
parse_sspace_summary.sh - script for extraction of QC statistics from SSPACE output files
run_fastqc.sh - script to run FastQC aster trimming
run_gzip.sh - script to compress data
run_reference_based_assembly_final.sh - pipeline for mixed reference based assembly and Skesa, scaffolding and quality control
run_trimmomatic.sh - script for trimming

Direcroties with data:
Raw reads: /projects/data/team1_genomeAssembly/downloads
Trimmed data: /projects/data/team1_genomeAssembly/trimming2/fastq/trimmed
Mash sketches: /projects/data/team1_genomeAssembly/reference_based_assembly/mash_sketches_samples, /projects/data/team1_genomeAssembly/reference_based_assembly
Reference genomes: /projects/data/team1_genomeAssembly/reference_based_assembly/reference_genomes
Final reference based assemblies: /projects/data/team1_genomeAssembly/reference_based_assembly/assembly_all
Final SPADEs results: /projects/data/team1_genomeAssembly/denovo_spades_result
Final Skesa results: /projects/data/team1_genomeAssembly/denovo_skesa/sspaceOutput
