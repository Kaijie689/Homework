# Homework description
# Rice transcriptome analysis
# Download data form NCBI
prefetch -O {SRA_DIR} {SRR_ID}

# Convert SRA data to FASTA
fasterq-dump {input.sra} -O {FASTQ_DIR} -p -3

# Quality Control Check the quality of sequencing and select reads with a high-queality greater than 25 for alignment
mkdir -p {QC_DIR}
fastqc -o {QC_DIR} {input.fastq1}
fastqc -o {QC_DIR} {input.fastq2}
      

# Build index Used to construct indexes for genomes. Its role is to construct the sequence of the reference genome into an index file suitable for HISAT2 to quickly align it.

# sequence alignment, HISAT2 was used to compare the filtered fastq file with the reference genome to obtain the bam file required for downstream analysis 

# Gene Count, Count data of genes were extracted from the alignment results for further analysis of gene expression
