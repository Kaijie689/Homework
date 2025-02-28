# 定义链接和其他参数
sra_id = "SRR2089753"  # 我们要分析的 SRA 文件

# 目标目录
RAW_DATA_DIR = "raw_data"
CLEAN_DATA_DIR = "clean_data"
QC_DIR = "qc_reports"
ALIGN_DIR = "aligned"
GENE_COUNT_DIR = "gene_counts"
REFERENCE_DIR = "reference"

rule all:
    input:
        f"{GENE_COUNT_DIR}/{sra_id}_gene_counts.tsv"


# 转换 SRA 文件为 FASTQ 文件
rule fasterq_dump:
    input:
        sra = f"sra/{sra_id}.sra"
    output:
        fastq_1 = f"{RAW_DATA_DIR}/{sra_id}_1.fastq",
        fastq_2 = f"{RAW_DATA_DIR}/{sra_id}_2.fastq"
    shell:
        """
        mkdir -p {RAW_DATA_DIR}
        fasterq-dump -p -3 -O {RAW_DATA_DIR} {input.sra}
        """

# 下载参考基因组文件
rule download_reference_files:
    output:
        genome = f"{REFERENCE_DIR}/all.con",
        annotation = f"{REFERENCE_DIR}/all.gff3"
    shell:
        """
        wget -O {output.genome} https://rice.uga.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.con
        wget -O {output.annotation} https://rice.uga.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.gff3
        """

# 构建参考基因组索引
rule hisat_index:
    input:
        genome = f"{REFERENCE_DIR}/all.con"
    output:
        genome_index_dir = directory(f"{REFERENCE_DIR}/genome_indices/o_sativa_index")
    shell:
        """
        mkdir -p {output.genome_index_dir}  # 确保 o_sativa_index 子目录存在
        hisat2-build -p 24 -f {input.genome} {output.genome_index_dir}/o_sativa_index  # 输出到子目录
        """

# 进行比对：
rule hisat2_align:
    input:
        fastq_1 = f"{RAW_DATA_DIR}/{sra_id}_1.fastq",  # 修改为原始 FASTQ 文件
        fastq_2 = f"{RAW_DATA_DIR}/{sra_id}_2.fastq",  # 修改为原始 FASTQ 文件
        index = f"{REFERENCE_DIR}/genome_indices/o_sativa_index"
    output:
        sam = f"{ALIGN_DIR}/{sra_id}_mapped_reads.sam",
        report = f"{ALIGN_DIR}/{sra_id}_mapping_report.txt"
    log:
        f"{QC_DIR}/{sra_id}_mapping.log"
    threads:
        6
    shell:
        """
        hisat2 --dta --fr --no-mixed --no-discordant --time --new-summary --no-unal \
            -p {threads} -x {input.index} -1 {input.fastq_1} -2 {input.fastq_2} \
            -S {output.sam} --summary-file {output.report} 2>> {log}
        """


# 将 SAM 文件转换为 BAM 文件，并进行排序
rule sam_to_bam:
    input:
        sam = f"{ALIGN_DIR}/{sra_id}_mapped_reads.sam"
    output:
        bam = f"{ALIGN_DIR}/{sra_id}.bam",
        bam_sorted = f"{ALIGN_DIR}/{sra_id}_sorted.bam",
        index = f"{ALIGN_DIR}/{sra_id}_sorted.bam.bai",
    log:
        f"{QC_DIR}/{sra_id}_samtools.log"
    shell:
        '''
        samtools view {input.sam} -b -o {output.bam} 2>> {log}
        samtools sort {output.bam} -O bam -o {output.bam_sorted} 2>> {log}
        samtools index -b {output.bam_sorted} -o {output.index} 2>> {log}
        '''

# 将 GFF 文件转换为 GTF 格式
rule gff_to_gtf:
    input:
        gff = f"{REFERENCE_DIR}/all.gff3"
    output:
        gtf = f"{REFERENCE_DIR}/all.gtf"
    shell:
        """
        gffread {input.gff} -T -o {output.gtf}
        """

# 基因计数：使用 featureCounts 计算基因计数
rule gene_counts:
    input:
        bam_sorted = f"{ALIGN_DIR}/{sra_id}_sorted.bam",
        annotation = f"{REFERENCE_DIR}/all.gtf"
    output:
        gene_counts = f"{GENE_COUNT_DIR}/{sra_id}_gene_counts.tsv",
        gene_summary = f"{GENE_COUNT_DIR}/{sra_id}_gene_counts.summary"
    log:
        f"{QC_DIR}/{sra_id}_featurecounts.log"
    shell:
        """
        featureCounts -t exon -g gene_id -s 2 -p -B -C --largestOverlap --verbose -F GTF \
            -a {input.annotation} -o {output.gene_counts} {input.bam_sorted} &>> {log} && \
            mv {output.gene_counts}.summary {output.gene_summary}
        """
