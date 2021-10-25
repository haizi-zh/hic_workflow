from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()


def get_ref_genome(ref_genome):
    ref_map = {
        "hg19": [os.path.normpath(v) for v in expand(f"data/ref_genome/hg19/human_g1k_v37.fa.gz{{ext}}", ext=["", ".amb", ".ann", ".bwt", ".fai", ".pac", ".sa"])],
        "hg38": [os.path.normpath(v) for v in expand(f"data/ref_genome/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz{{ext}}", ext=["", ".amb", ".ann", ".bwt", ".fai", ".pac", ".sa"])],
    }
    return ref_map[ref_genome]


def get_chrom_sizes(ref_genome):
    ref_map = {
        "hg19": "data/ref_genome/hg19/human_g1k_v37.chrom.sizes",
        "hg38": "data/ref_genome/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.chrom.sizes"
    }
    return ref_map[ref_genome]


# >>> Configuration >>>
FULL_CORES = config.get("FULL_CORES", 16)
PART_CORES = config.get("PART_CORES", 8)
MAX_CORES = config.get("MAX_CORES", 9999)
# PART_CORES_S = PART_CORES
# PART_CORES_M = config.get("PART_CORES_M", 8)
MEM_PER_CORE = config.get("MEM_PER_CORE", 1800)
WALL_TIME_MAX = config.get("WALL_TIME_MAX", 2880)
WALL_TIME_MIN = config.get("WALL_TIME_MIN", 720)

# # Default base directory for data files. Default: ./data
# DATA_DIR = config.get("DATA_DIR", os.path.abspath("data"))
# # Base directory for reference genome
# REF_GENOME_DIR = config.get("REF_GENOME_DIR", os.path.normpath(os.path.join(DATA_DIR, "ref_genome")))

# # Read group ID for bwa mem
# BWA_RG_ID = config.get("BWA_RG_ID", None)
# BWA_RG_ID_SUFFIX = config.get("BWA_RG_ID_SUFFIX", "")

# <<< Configuration <<<



# Generate FastQC reports
rule fastqc:
    input: expand("fastq/{{sample}}.R{rg}.fastq.gz", rg=[1, 2])
    output: expand("fastqc/{{sample}}.R{rg}_fastqc.{ext}", rg=[1, 2], ext=["html", "zip"])
    params:
        slurm_job_label=lambda wildcards: f"fastqc.{wildcards.sample}",
    threads: 2
    resources:
        time=WALL_TIME_MAX,
        time_min=WALL_TIME_MIN,
        attempt=lambda wildcards, threads, attempt: attempt
    conda: "conda/conda.fastqc.yaml"
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        fastqc -t {threads} -o $tmpdir -f fastq {input}

        cp -av $tmpdir/*_fastqc.html $tmpdir/*_fastqc.zip `dirname {output[1]}`
        """


# Generate FastQC reports
rule fastqc_trimmed:
    input: expand("trimmomatic/{{sample}}.R{rg}.fastq.gz", rg=[1, 2])
    output: expand("fastqc_trimmed/{{sample}}.R{rg}_fastqc.{ext}", rg=[1, 2], ext=["html", "zip"])
    params:
        slurm_job_label=lambda wildcards: f"fastqc_trimmed.{wildcards.sample}",
    threads: 2
    resources:
        time=WALL_TIME_MAX,
        time_min=WALL_TIME_MIN,
        attempt=lambda wildcards, threads, attempt: attempt
    conda: "conda/conda.fastqc.yaml"
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        fastqc -t {threads} -o $tmpdir -f fastq {input}

        cp -av $tmpdir/*_fastqc.html $tmpdir/*_fastqc.zip `dirname {output[1]}`
        """


rule trimmomatic:
    input: 
        fastq=expand("fastq/{{sample}}.{read_name}.fastq.gz", read_name=["R1", "R2"]),
        adapter=HTTP.remote("https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE-2.fa"),
    output:
        expand("trimmomatic/{{sample}}.R{rg}{pairing}fastq.gz", rg=[1, 2], pairing=[".", ".unpaired."]),
    log: "trimmomatic/{sample}.trimmomatic.summary.txt",
    params:
        slurm_job_label=lambda wildcards: f"trimmomatic.{wildcards.sample}",
        trimmer=lambda wildcards, input: config.get("TRIMMER", f"ILLUMINACLIP:{input.adapter}:2:30:10:2:keepBothReads LEADING:5 TRAILING:5 AVGQUAL:20 SLIDINGWINDOW:4:10 MINLEN:36"),
    threads: lambda wildcards: min(PART_CORES, MAX_CORES)
    resources:
        time=WALL_TIME_MAX,
        time_min=WALL_TIME_MIN,
        attempt=lambda wildcards, threads, attempt: attempt
    conda: "conda/conda.trimmomatic.yaml"
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        trimmomatic PE -threads {threads} \
        {input.fastq} \
        $tmpdir/r1.fastq.gz $tmpdir/r1.unpaired.fastq.gz \
        $tmpdir/r2.fastq.gz $tmpdir/r2.unpaired.fastq.gz \
        {params.trimmer} \
        2>&1 | tee {log}

        mv $tmpdir/r1.fastq.gz {output[0]}.tmp
        mv $tmpdir/r1.unpaired.fastq.gz {output[1]}.tmp
        mv $tmpdir/r2.fastq.gz {output[2]}.tmp
        mv $tmpdir/r2.unpaired.fastq.gz {output[3]}.tmp

        mv {output[0]}.tmp {output[0]}
        mv {output[1]}.tmp {output[1]}
        mv {output[2]}.tmp {output[2]}
        mv {output[3]}.tmp {output[3]}
        """


# Align raw reads using bwa-mem
rule bwa:
    input: 
        ref=lambda wildcards: get_ref_genome(wildcards.ref_genome),
        fastq=expand("trimmomatic/{{sample}}.R{rg}.fastq.gz", rg=[1, 2]),
    output: 
        "temp/{sample}.{ref_genome,[^\\.]+}.raw.bam",
    params:
        slurm_job_label=lambda wildcards: f"bwa_raw.{wildcards.sample}.{wildcards.ref_genome}",
    threads: lambda wildcards, attempt: min(int(FULL_CORES * (0.5 + 0.5 * attempt)), MAX_CORES)
    resources:
        time=WALL_TIME_MAX,
        time_min=WALL_TIME_MIN,
        attempt=lambda wildcards, threads, attempt: attempt
    conda: "conda/conda.bwa.yaml"
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        bwa mem -t {threads} -SP5M {input.ref[0]} {input.fastq} | samtools view -Shb - > $tmpdir/output.bam

        mv $tmpdir/output.bam {output}.tmp
        mv {output}.tmp {output}
        """


rule pairsam_parse_sort:
    input: 
        bam="temp/{sample}.{ref_genome}.raw.bam",
        chrom_sizes=lambda wildcards: get_chrom_sizes(wildcards.ref_genome)
    output:
        temp("temp/{sample}.{ref_genome,[^\\.]+}.sam.pairs.gz")
    params:
        slurm_job_label=lambda wildcards: f"pairsam_parse_sort.{wildcards.sample}.{wildcards.ref_genome}",
        min_mapq=config.get("PAIRSAM_MIN_MAPQ", 30),
    threads: lambda wildcards: min(config.get("PAIRSAM_PARSE_SORT_CORES", 16), MAX_CORES)
    resources:
        mem_mb=lambda wildcards, threads: int(threads * MEM_PER_CORE - 3200),
        time=WALL_TIME_MAX,
        time_min=WALL_TIME_MAX,
        attempt=lambda wildcards, threads, attempt: attempt
    conda: "conda.yaml"
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        echo "Total threads: {threads}"
        sort_threads=$(({threads}-2))
        sort_mem=$(($sort_threads*{MEM_PER_CORE}))
        echo "Sorting threads: $sort_threads"

        samtools view -h {input.bam} |
        pairtools parse -c {input.chrom_sizes} --add-columns mapq \
        --assembly {wildcards.ref_genome} \
        --min-mapq {params.min_mapq} |
        pairtools sort --nproc $sort_threads --memory "$sort_mem"M \
        --compress-program gzip \
        --tmpdir $tmpdir \
        --output $tmpdir/output.sam.pairs.gz

        mv $tmpdir/output.sam.pairs.gz {output}.tmp
        mv {output}.tmp {output}
        """


rule pairsam_markdup:
    input: 
        "temp/{sample}.{ref_genome}.sam.pairs.gz",
    output:
        pairsam=temp("temp/{sample}.{ref_genome,[^\\.]+}.marked.sam.pairs.gz"),
        pairsam_idx=temp("temp/{sample}.{ref_genome}.marked.sam.pairs.gz.px2"),
    params:
        slurm_job_label=lambda wildcards: f"pairsam_markdup.{wildcards.sample}.{wildcards.ref_genome}",
    threads: 4
    resources:
        time=WALL_TIME_MAX,
        time_min=WALL_TIME_MIN,
        attempt=lambda wildcards, threads, attempt: attempt
    conda: "conda.yaml"
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        pairtools dedup --mark-dups --output-dups - --output-unmapped - --output $tmpdir/marked.sam.pairs.gz {input}
        pairix $tmpdir/marked.sam.pairs.gz

        mv $tmpdir/marked.sam.pairs.gz {output.pairsam}.tmp
        mv $tmpdir/marked.sam.pairs.gz.px2 {output.pairsam_idx}.tmp
        mv {output.pairsam}.tmp {output.pairsam}
        mv {output.pairsam_idx}.tmp {output.pairsam_idx}
        """


rule pairsam_filter:
    input: 
        pairsam="temp/{sample}.{ref_genome}.marked.sam.pairs.gz",
        pairsam_idx="temp/{sample}.{ref_genome}.marked.sam.pairs.gz.px2",
        chrom_sizes=lambda wildcards: get_chrom_sizes(wildcards.ref_genome)
    output:
        bam="pairsam/{sample}.{ref_genome,[^\\.]+}.lossless.bam",
        dedup="pairsam/{sample}.{ref_genome}.dedup.sam.pairs.gz",
        dedup_idx="pairsam/{sample}.{ref_genome}.dedup.sam.pairs.gz.px2",
        unmapped="pairsam/{sample}.{ref_genome}.unmapped.sam.pairs.gz",
    params:
        slurm_job_label=lambda wildcards: f"pairsam_filter.{wildcards.sample}.{wildcards.ref_genome}",
    threads: 4
    resources:
        time=WALL_TIME_MAX,
        time_min=WALL_TIME_MIN,
        attempt=lambda wildcards, threads, attempt: attempt
    conda: "conda.yaml"
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        # Generate lossless bam
        pairtools split --output-sam $tmpdir/lossless.bam {input.pairsam}
        
        # Select UU, UR, RU reads
        pairtools select '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' \
        --output-rest $tmpdir/unmapped.pairs.gz \
        --output $tmpdir/output.gz \
        {input.pairsam}
        
        pairtools split --output-pairs $tmpdir/output1.gz $tmpdir/output.gz
        pairtools select 'True' --chrom-subset {input.chrom_sizes} -o $tmpdir/dedup.pairs.gz $tmpdir/output1.gz

        # sanity check & indexing    
        pairix $tmpdir/dedup.pairs.gz

        mv $tmpdir/lossless.bam {output.bam}.tmp
        mv $tmpdir/dedup.pairs.gz {output.dedup}.tmp
        mv $tmpdir/unmapped.pairs.gz {output.unmapped}.tmp

        mv {output.bam}.tmp {output.bam}
        mv {output.dedup}.tmp {output.dedup}
        mv {output.unmapped}.tmp {output.unmapped}
        mv $tmpdir/dedup.pairs.gz.px2 {output.dedup_idx}
        """


rule pairsam2hic:
    input: 
        pairsam="pairsam/{sample}.{ref_genome}.dedup.sam.pairs.gz",
        pairsam_idx="pairsam/{sample}.{ref_genome}.dedup.sam.pairs.gz.px2",
        chrom_sizes=lambda wildcards: get_chrom_sizes(wildcards.ref_genome),
        juicer=HTTP.remote("https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar"),
    output:
        temp("temp/{sample}.{ref_genome,[^\\.]+}.raw.hic"),
    params:
        slurm_job_label=lambda wildcards: f"pairsam2hic.{wildcards.sample}.{wildcards.ref_genome}",
        min_mapq=config.get("PAIRSAM_MIN_MAPQ", 30),
        custom_res=config.get("CUSTOM_RES", ""),
        min_res=config.get("MIN_RES", ""),
    threads: 4
    resources:
        mem_mb=lambda wildcards, threads: int(threads * MEM_PER_CORE - 1200),
        time=WALL_TIME_MAX,
        time_min=WALL_TIME_MIN,
        attempt=lambda wildcards, threads, attempt: attempt
    conda: "conda.yaml"
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        custom_res={params.custom_res}
        min_res={params.min_res}

        if [[ ! -z $custom_res && ! -z $min_res ]]
        then
            lowest_custom_res=${{custom_res//,*/}}
            if [[ $lowest_custom_res != $min_res ]]
            then
                echo "Lowest custom res (-u) does not match min_res (-r)."
                exit 1
            fi
        fi

        if [[ ! -z $custom_res ]]
        then
            custom_res_opt="-r $custom_res"
        else
            custom_res_opt=""
        fi

        java -Xmx{resources.mem_mb}M -Xms{resources.mem_mb}M -jar {input.juicer} \
        pre -n {input.pairsam} --threads {threads} -q {params.min_mapq} $custom_res_opt \
        $tmpdir/output.hic {input.chrom_sizes} \

        mv $tmpdir/output.hic {output}.tmp
        mv {output}.tmp {output}
        """


rule hic_stats:
    input: 
        pairsam="pairsam/{sample}.{ref_genome}.dedup.sam.pairs.gz",
        pairsam_idx="pairsam/{sample}.{ref_genome}.dedup.sam.pairs.gz.px2",
    output:
        "hic/{sample}.{ref_genome,[^\\.]+}.stats.txt"
    params:
        slurm_job_label=lambda wildcards: f"hic_stats.{wildcards.sample}.{wildcards.ref_genome}",
    threads: 1
    resources:
        time=WALL_TIME_MAX,
        time_min=WALL_TIME_MIN,
        attempt=lambda wildcards, threads, attempt: attempt
    conda: "conda.yaml"
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        pairtools stats {input.pairsam} > $tmpdir/output.stats.txt
        mv $tmpdir/output.stats.txt {output}
        """


rule juicer_addnorm:
    input: 
        hic="temp/{sample}.{ref_genome}.raw.hic",
        juicer=HTTP.remote("https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar"),        
    output:
        "hic/{sample}.{ref_genome,[^\\.]+}.hic",
    params:
        slurm_job_label=lambda wildcards: f"juicer_addnorm.{wildcards.sample}.{wildcards.ref_genome}",
        custom_res=config.get("CUSTOM_RES", ""),
        min_res=config.get("MIN_RES", ""),
    threads: lambda wildcards: min(config.get("JUICER_ADDNORM_CORES", 6), MAX_CORES)
    resources:
        mem_mb=lambda wildcards, threads: int(threads * MEM_PER_CORE - 1200),
        time=WALL_TIME_MAX,
        time_min=WALL_TIME_MIN,
        attempt=lambda wildcards, threads, attempt: attempt
    conda: "conda.yaml"
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        custom_res={params.custom_res}
        min_res={params.min_res}

        if [[ ! -z $custom_res && ! -z $min_res ]]
        then
            lowest_custom_res=${{custom_res//,*/}}
            if [[ $lowest_custom_res != $min_res ]]
            then
                echo "Lowest custom res (-u) does not match min_res (-r)."
                exit 1
            fi
        fi

        if [[ ! -z $min_res ]]
        then
            min_res_opt="-w $min_res"
        else
            min_res_opt=""
        fi

        cp -a {input.hic} $tmpdir/output.hic

        java -Xmx{resources.mem_mb}M -Xms{resources.mem_mb}M -jar {input.juicer} \
        addNorm $min_res_opt -d -F --threads {threads} $tmpdir/output.hic

        mv $tmpdir/output.hic {output}.tmp
        mv {output}.tmp {output}
        """
