import os

# wgs pipeline -----------------
class GenomePipeline:
    def __init__(self, reference: str, thread: int, flag: str, output: str, prefix: str):
        """
        Initialize the pipeline with necessary parameters.

        :param reference: Path to the reference genome.
        :param thread: Number of threads to use for parallel processing.
        :param flag: Identifier flag (e.g., sample information).
        :param output: Path to the output directory.
        :param prefix: Prefix for output files.
        """
        
        self.reference = reference
        self.thread = thread
        self.flag = flag
        self.output = output
        self.prefix = prefix

    def bwa_mem(self, input1: str, input2: str) -> str:
        """
        Generate the BWA MEM command for sequence alignment.

        :param input1: Path to the first input FASTQ file.
        :param input2: Path to the second input FASTQ file.
        :return: BWA MEM command as a string.
        """
        bwa_mem_cmd = (
            f"bwa mem -t {self.thread} {self.reference} {input1} {input2} -R '@RG\\tID:{self.flag}\\tSM:{self.prefix}' "
            f"> {self.output}/{self.prefix}.sam"
        )
        return bwa_mem_cmd

    def samtools_sort(self) -> str:
        """
        Generate the SAMtools sort command to sort SAM files and convert to BAM format.

        :return: SAMtools sort command as a string.
        """
        samtools_sort_cmd = (
            f"samtools sort -@ {self.thread} -o {self.output}/{self.prefix}.bam {self.output}/{self.prefix}.sam"
        )
        return samtools_sort_cmd

    def samtools_index(self) -> str:
        """
        Generate the SAMtools index command to index the BAM file.

        :return: SAMtools index command as a string.
        """
        samtools_index_cmd = f"samtools index -@ {self.thread} {self.output}/{self.prefix}.bam"
        return samtools_index_cmd

    def samtools_flagstat(self) -> str:
        """
        Generate the SAMtools flagstat command to compute statistics on the BAM file.

        :return: SAMtools flagstat command as a string.
        """
        samtools_flagstat_cmd = (
            f"samtools flagstat -@ {self.thread} {self.output}/{self.prefix}.bam > {self.output}/{self.prefix}.bam.flagstat"
        )
        return samtools_flagstat_cmd

    def rm_dup(self) -> str:
        """
        Generate the GATK MarkDuplicates command to remove duplicate reads.

        :return: GATK MarkDuplicates command as a string.
        """
        gatk_markdup_cmd = (
            f"gatk MarkDuplicates -I {self.output}/{self.prefix}.bam -O {self.output}/{self.prefix}.rmdup.bam "
            f"-M {self.output}/{self.prefix}.rmdup.metrics.txt"
        )
        return gatk_markdup_cmd
      
    def samtools_sort_index_c(self) -> str:
        samtools_sort_index_c = (
            f"samtools sort -t {self.thread} -o {self.output}/{self.prefix}_sorted_rmdup.bam {self.output}/{self.prefix}.rmdup.bam\n"
            f"samtools index -c {self.output}/{self.prefix}_sorted_rmdup.bam"
        )
        return samtools_sort_index_c
      
    def gatk_hapcaller(self) -> str:
        gatk_hapcaller = (
            f"gatk HaplotypeCaller -R {self.reference} -I {self.output}/{self.prefix}_sorted_rmdup.bam -O {self.output}/{self.prefix}.vcf"
        )
        return gatk_hapcaller

    def select_variants(self) -> str:
        """
        Generate the GATK SelectVariants command to extract SNPs and INDELs from the VCF file.

        :return: GATK SelectVariants commands as a string.
        """
        select_snps_cmd = (
            f"gatk SelectVariants -R {self.reference} -V {self.output}/{self.prefix}.vcf "
            f"--select-type-to-include SNP -O {self.output}/{self.prefix}_snps.vcf"
        )
        select_indels_cmd = (
            f"gatk SelectVariants -R {self.reference} -V {self.output}/{self.prefix}.vcf "
            f"--select-type-to-include INDEL -O {self.output}/{self.prefix}_indels.vcf"
        )
        return select_snps_cmd + "\n" + select_indels_cmd

    def filter_variants(self) -> str:
        """
        Generate the GATK VariantFiltration command to filter SNPs and INDELs based on quality metrics.

        :return: GATK VariantFiltration commands as a string.
        """
        filter_snps_cmd = (
            f"gatk VariantFiltration -R {self.reference} -V {self.output}/{self.prefix}_snps.vcf "
            f"--filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0' --filter-name 'FILTER' "
            f"-O {self.output}/{self.prefix}_filtered_snps.vcf"
        )
        filter_indels_cmd = (
            f"gatk VariantFiltration -R {self.reference} -V {self.output}/{self.prefix}_indels.vcf "
            f"--filter-expression 'QD < 2.0 || FS > 200.0' --filter-name 'FILTER' "
            f"-O {self.output}/{self.prefix}_filtered_indels.vcf"
        )
        return filter_snps_cmd + "\n" + filter_indels_cmd

    def generate_slurm_script(self, input1: str, input2: str, slurm_file: str):
        """
        Generate a SLURM script for the pipeline.

        :param input1: Path to the first input FASTQ file.
        :param input2: Path to the second input FASTQ file.
        :param slurm_file: Path to the SLURM script to generate.
        """
        try:
            with open(slurm_file, 'w') as f:
                f.write("#!/bin/bash\n")
                f.write(f"#SBATCH --job-name={self.prefix}_pipeline\n")
                f.write(f"#SBATCH --auks=yes\n")
                f.write(f"#SBATCH --output={self.output}/{self.prefix}_%j.out\n")
                f.write(f"#SBATCH --error={self.output}/{self.prefix}_%j.err\n")
                f.write(f"#SBATCH --ntasks=1\n")
                f.write(f"#SBATCH --cpus-per-task={self.thread}\n")
                f.write(f"#SBATCH --mem=200G\n")
                f.write(f"#SBATCH --mail-type=END\n")
                f.write(f"#SBATCH --mail-user=tan@ipk-gatersleben.de\n")
                f.write(f"#SBATCH --time=168:00:00\n")
                f.write(f"module load bwa\n")
                f.write(f"module load samtools\n")
                f.write(f"module load gatk\n")
                f.write("\n")
                # Write BWA MEM command
                f.write(f"{self.bwa_mem(input1, input2)}\n")
                
                # Write SAMtools Sort command
                f.write(f"{self.samtools_sort()}\n")
                
                # Write SAMtools Index command
                f.write(f"{self.samtools_index()}\n")
                
                # Write SAMtools Flagstat command
                f.write(f"{self.samtools_flagstat()}\n")
                
                # Write GATK MarkDuplicates command
                f.write(f"{self.rm_dup()}\n")
                f.write(f"{self.samtools_sort_index_c()}\n")
                f.write(f"{self.gatk_hapcaller()}\n")
                # Write GATK SelectVariants command
                f.write(f"{self.select_variants()}\n")

                # Write GATK VariantFiltration command
                f.write(f"{self.filter_variants()}\n")

            print(f"SLURM script generated: {slurm_file}")
        except Exception as e:
            print(f"Failed to generate SLURM script: {e}")
# 初始化参数
sample = "sample2"  # 可以更改为任何样本名称
reference = "/filer-5/agruppen/PBP/tan/indexDir/barley_index/220830_Bowman_pseudomolecules_and_unplaced_contigs_CPclean.fasta"
output_dir = "/filer-5/agruppen/PBP/tan/wgs/workspace/bowman/second_time"
prefix = f"{sample}_output"
slurm_file = f"U://workDir/script/{sample}_bowman_second_time.sh"

pipeline = GenomePipeline(
    reference=reference,
    thread=32,
    flag=sample,
    output=output_dir,
    prefix=prefix
)

# 生成 SLURM 脚本
input1 = f"/filer-5/agruppen/PBP/tan/wgs/workspace/{sample}_1.fq.gz"
input2 = f"/filer-5/agruppen/PBP/tan/wgs/workspace/{sample}_2.fq.gz"

pipeline.generate_slurm_script(
    input1=input1,
    input2=input2,
    slurm_file=slurm_file
)

)
# reference:
# /filer-5/agruppen/PBP/tan/indexDir/barley_index/220809_Foma_pseudomolecules_and_unplaced_contigs_CPclean.fasta
# /filer-5/agruppen/PBP/tan/indexDir/barley_index/220830_Bowman_pseudomolecules_and_unplaced_contigs_CPclean.fasta
# /filer-5/agruppen/PBP/tan/indexDir/barley_index/Barley_Morex_V2_pseudomolecules.fa



# rnaseq pipeline --------------
import os

class RNASeqPipeline:
    def __init__(self, raw_data_dir: str, processed_data_dir: str, index_dir: str, thread: int):
        """
        Initialize the RNA-Seq pipeline with necessary parameters.

        :param raw_data_dir: Path to the directory containing raw FASTQ files.
        :param processed_data_dir: Path to the directory for processed data output.
        :param index_dir: Path to the directory containing the Kallisto index.
        :param thread: Number of threads to use for parallel processing.
        """
        self.raw_data_dir = raw_data_dir
        self.processed_data_dir = processed_data_dir
        self.index_dir = index_dir
        self.thread = thread
        self.array = array

    def fastp_command(self, sample_folder: str, sample_name: str, pair1: str, pair2: str) -> str:
        """
        Generate the fastp command for quality control of FASTQ files.

        :param sample_folder: Path to the sample folder.
        :param sample_name: Name of the sample.
        :param pair1: Suffix for the first read pair.
        :param pair2: Suffix for the second read pair.
        :return: fastp command as a string.
        """
        fastp_cmd = (
            f"fastp -i {sample_folder}/{sample_name}_{pair1}.fq.gz "
            f"-I {sample_folder}/{sample_name}_{pair2}.fq.gz "
            f"-o {self.processed_data_dir}/{sample_name}_clean_{pair1}.fq.gz "
            f"-O {self.processed_data_dir}/{sample_name}_clean_{pair2}.fq.gz "
            f"--thread {self.thread}"
        )
        return fastp_cmd

    def kallisto_command(self, sample_name: str, pair1: str, pair2: str) -> str:
        """
        Generate the Kallisto command for transcript quantification.

        :param sample_name: Name of the sample.
        :param pair1: Suffix for the first read pair.
        :param pair2: Suffix for the second read pair.
        :return: Kallisto command as a string.
        """
        kallisto_cmd = (
            f"kallisto quant -i {self.index_dir}/kallisto_index.idx "
            f"-o {self.processed_data_dir}/{sample_name}_output "
            f"{self.processed_data_dir}/{sample_name}_clean_{pair1}.fq.gz "
            f"{self.processed_data_dir}/{sample_name}_clean_{pair2}.fq.gz "
            f"--paired -t {self.thread}"
        )
        return kallisto_cmd

    def generate_slurm_script(self, slurm_file: str, pair1: str = "1", pair2: str = "2"):
        """
        Generate a SLURM script for the RNA-Seq pipeline.

        :param slurm_file: Path to the SLURM script to generate.
        :param pair1: Suffix for the first read pair.
        :param pair2: Suffix for the second read pair.
        """
        try:
            with open(slurm_file, 'w') as f:
                f.write("#!/bin/bash\n")
                f.write(f"#SBATCH --nodes=1\n")
                f.write(f"#SBATCH --cpus-per-task={self.thread}\n")
                f.write(f"#SBATCH --output={self.processed_data_dir}/job_%j.out\n")
                f.write(f"#SBATCH --error={self.processed_data_dir}/job_%j.err\n")
                f.write(f"#SBATCH --mem=200G\n")
                f.write(f"#SBATCH --mail-type=END\n")
                f.write(f"#SBATCH --mail-user=tan@ipk-gatersleben.de\n")
                f.write(f"#SBATCH --time=168:00:00\n")
                f.write(f"#SBATCH --mail-type=END\n")
                f.write(f"#SBATCH --mail-user=your_email@example.com\n")
                f.write("\n")
                f.write("#!/bin/bash\n")
                f.write(f"#SBATCH --job-name={self.prefix}_pipeline\n")
                f.write(f"#SBATCH --auks=yes\n")
                f.write(f"#SBATCH --output={self.output}/{self.prefix}_%j.out\n")
                f.write(f"#SBATCH --error={self.output}/{self.prefix}_%j.err\n")
                f.write(f"#SBATCH --ntasks=1\n")
                f.write(f"#SBATCH --cpus-per-task={self.thread}\n")

                f.write(f"module load bwa\n")
                f.write(f"module load samtools\n")
                f.write(f"module load gatk\n")
                f.write("\n")
                f.write("source /etc/profile\n")
                f.write("module load fastp\n")
                f.write("module load kallisto/0.46.1\n")
                f.write("\n")

                f.write(f"raw_data_dir={self.raw_data_dir}\n")
                f.write(f"processed_data_dir={self.processed_data_dir}\n")
                f.write(f"index_dir={self.index_dir}\n")
                f.write(f"pair1={pair1}\n")
                f.write(f"pair2={pair2}\n")
                f.write("\n")

                f.write("mkdir -p $processed_data_dir/kallisto\n")
                f.write("\n")

                f.write("for sample_folder in $raw_data_dir/*; do\n")
                f.write("    if [ -d \"$sample_folder\" ]; then\n")
                f.write("        sample_name=$(basename $sample_folder)\n")
                f.write("        echo Processing $sample_name\n")
                f.write(f"        {self.fastp_command('$sample_folder', 'sample_name', '$pair1', '$pair2')}\n")
                f.write(f"        {self.kallisto_command('sample_name', '$pair1', '$pair2')}\n")
                f.write("    fi\n")
                f.write("done\n")

            print(f"SLURM script generated: {slurm_file}")
        except Exception as e:
            print(f"Failed to generate SLURM script: {e}")


