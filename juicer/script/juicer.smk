import os
import math


localrules: all, selected, get_ver_info, merge_abnorm_unmapped

# load_bwa = "module load  GCCcore/7.3.0 BWA/0.7.17"
# load_java = "module load Java/1.8.0_162" 


# >>> Configuration >>>
# The ratio between total memory in MB to the number of CPU cores.
#
# This config item is important for some clusters such as PSC Bridges-2,
# where the total memory available to a job depends on the number of
# cores requested. In this case, you specify the memory requirement for a 
# job, and snakemake knows how many cores to request. Alternatively, you
# can also specify the number of cores, and snakemake knows how much memory is 
MEM_PER_CORE = int(config.get("MEM_PER_CORE", 2000))

WALL_TIME_MAX = int(config.get("WALL_TIME_MAX", 2880))

# Number of threads for bwa mem
BWA_THREADS = int(config.get("BWA_THREADS", 4))

# Memory allocation for sorting, in MB
SORT_MEM = int(config.get("SORT_MEM", 4000))

# Memory allocation for Juicer tools pre, in MB
JUICER_PRE_MEM = int(config.get("JUICER_PRE_MEM", 4000))
JUICER_HICCUPS_MEM = int(config.get("JUICER_HICCUPS_MEM", 4000))
JUICER_ARROWHEAD_MEM = int(config.get("JUICER_ARROWHEAD_MEM", 4000))

JUICER_PRE_RES = config.get("JUICER_PRE_RES", "2500000,1000000,500000,250000,100000,50000,25000,10000,5000")

SITE = config.get("SITE", "")
SITE_FILE = config.get("SITE_FILE", "")
LIGATION = config.get("LIGATION", "")
JUST_EXACT = config.get("JUST_EXACT", "")
FRAG = config.get("FRAG", "")

GENOME = config.get("GENOME")
assert GENOME in ["mm9", "mm10", "hg38", "hg19"]
GENOME_PATH = config.get("GENOME_PATH", "")
REF_GENOME = config.get("REF_GENOME", "")


SAMPLE_IDS = config.get("SAMPLE_IDS", "")
# <<< Configuration <<<


# Get the path to Juicer scripts, which is the same directory as the snakemake file
def get_snakemake_file_path():
    snakefile = globals()['workflow'].snakefile
    return os.path.split(snakefile)[0]


## Set ligation junction based on restriction enzyme
def get_ligation_junction():
    if LIGATION != "":
        return LIGATION

    assert SITE != ""

    ligation_map = {
        "HindIII": "AAGCTAGCTT",
        "DpnII": "GATCGATC",
        "MboI": "GATCGATC",
        "NcoI": "CCATGCATGG",
        "Arima": "'(GAATAATC|GAATACTC|GAATAGTC|GAATATTC|GAATGATC|GACTAATC|GACTACTC|GACTAGTC|GACTATTC|GACTGATC|GAGTAATC|GAGTACTC|GAGTAGTC|GAGTATTC|GAGTGATC|GATCAATC|GATCACTC|GATCAGTC|GATCATTC|GATCGATC|GATTAATC|GATTACTC|GATTAGTC|GATTATTC|GATTGATC)'",
        "none": "XXXX"
    }

    if SITE in ligation_map:
        return ligation_map[SITE]
    else:
        print(f"{SITE} not listed as recognized enzyme. Using {SITE_FILE} as site file")
        print("Ligation junction is undefined")
        return "XXXX"


# Return the path to reference genome files based on genome ID
def get_ref_genome(genome_id):
    if genome_id == "hg19":
        return "references/Homo_sapiens_assembly19.fasta"
    elif genome_id == "mm10":
        return "references/Mus_musculus_assembly10.fasta"
    else:
        sys.exit(f"Invalid genome ID: {genome_id}")


def get_sample_subid(sample):
    sample_subid, = glob_wildcards(f"fastq/{sample}/{{sample_subid}}_R1.fastq.gz")
    return sample_subid


def mem_for_attempt(mem_mb, attempt):
    return int(mem_mb * (0.5 + 0.5 * attempt))


def threads_for_mem(mem_mb, mem_per_core=MEM_PER_CORE):
    return int(math.ceil(mem_mb / mem_per_core) + 2)


juicer_version = "1.6"
juicer_home = get_snakemake_file_path()

if SITE == "none":
    FRAG = ""

if SITE_FILE == "":
    SITE_FILE = f"restriction_sites/{GENOME}_{SITE}.txt"

ligation = get_ligation_junction()


if GENOME_PATH == "":
    GENOME_PATH = GENOME


# Aggregation: all fastq files are processed
sample_ids, = glob_wildcards("fastq/{sample,[^/]+}")
sample_ids_selected = config.get("SAMPLE_IDS", "").split(",")

rule all:
    input:
        [expand(f"splits/{sample}/{{sample_subid}}.bam", sample_subid=get_sample_subid(sample)) for sample in sample_ids],
        # expand("splits/{sample}.fastq.gz_linecount.txt", sample=sample_ids),
        expand("aligned/{sample}.{inter}.hic", sample=sample_ids, inter=["inter_30"]),
        expand("aligned/{sample}.{inter}.hiccups.zip", sample=sample_ids, inter=["inter_30"]),
        expand("aligned/{sample}.{inter}.domain.bedpe", sample=sample_ids, inter=["inter_30"]),
        expand("aligned/{sample}.merged_nodups.txt.gz", sample=sample_ids),
        expand("aligned/{sample}.dups.txt.gz", sample=sample_ids),
        expand("aligned/{sample}.opt_dups.txt.gz", sample=sample_ids),
        expand("aligned/{sample}.collisions.txt", sample=sample_ids),
        expand("aligned/{sample}.fastq.gz_abnormal.sam.gz", sample=sample_ids),


rule selected:
    input:
        [expand(f"splits/{sample}/{{sample_subid}}.bam", sample_subid=get_sample_subid(sample)) for sample in sample_ids_selected],
        # expand("splits/{sample}.fastq.gz_linecount.txt", sample=sample_ids_selected),
        expand("aligned/{sample}.{inter}.hic", sample=sample_ids_selected, inter=["inter_30"]),
        expand("aligned/{sample}.{inter}.hiccups.zip", sample=sample_ids_selected, inter=["inter_30"]),
        expand("aligned/{sample}.{inter}.domain.bedpe", sample=sample_ids_selected, inter=["inter_30"]),
        expand("aligned/{sample}.dups.txt.gz", sample=sample_ids_selected),
        expand("aligned/{sample}.opt_dups.txt.gz", sample=sample_ids_selected),
        expand("aligned/{sample}.merged_nodups.txt.gz", sample=sample_ids_selected),
        expand("aligned/{sample}.collisions.txt", sample=sample_ids_selected),
        expand("aligned/{sample}.fastq.gz_abnormal.sam.gz", sample=sample_ids_selected),


rule get_ver_info:
    output: "aligned/header"
    params:
        threads=BWA_THREADS,
        juicer_version=juicer_version,
        juicer_home=juicer_home
    shell:
        """
        date > {output}

        # Get version numbers of all software
        echo -ne "Juicer version {params.juicer_version};" >> {output}

        set +o pipefail
        bwa 2>&1 | awk '$1=="Version:"{{printf(" BWA %s; ", $2)}}' >> {output}

        java -version 2>&1 | awk 'NR==1{{printf("%s; ", $0);}}' >> {output}
        java -jar {params.juicer_home}/juicer_tools.jar -V 2>&1 | awk '$1=="Juicer" && $2=="Tools"{{printf("%s; ", $0);}}' >> {output}
        """


rule count_ligations:
    input:
        r1="fastq/{sample}/{sample_subid}_R1.fastq.gz",
        r2="fastq/{sample}/{sample_subid}_R2.fastq.gz",
    output:
        norm=temp("temp/{sample}/{sample_subid}.fastq.gz_norm0.txt.res.txt"),
        linecount="splits/{sample}/{sample_subid}.fastq.gz_linecount.txt",
    params:
        ligation=ligation,
        label=lambda wildcards: f"count_ligation.{wildcards.sample}.{wildcards.sample_subid}",
    threads: lambda wildcards, attempt: int(1 * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: int(threads * MEM_PER_CORE),
        time=WALL_TIME_MAX,
        time_min=720,
        attempt=lambda wildcards, threads, attempt: attempt
    shell:
        """
        num1=$(paste <(gunzip -c {input.r1}) <(gunzip -c {input.r2}) | awk '!((NR+2)%4)' | grep -cE {params.ligation})
        num2=$(gunzip -c {input.r1} | wc -l | awk '{{print $1}}')
        echo -ne "$num1 " > {output.norm} 
        echo "$num2" > {output.linecount}
        """


rule alignment:
    input:
        r1="fastq/{sample}/{sample_subid}_R1.fastq.gz",
        r2="fastq/{sample}/{sample_subid}_R2.fastq.gz",
        ref=lambda wildcards: REF_GENOME if REF_GENOME != "" else get_ref_genome(GENOME)
    output:
        sam=temp("temp/{sample}/{sample_subid}.sam")
    threads: lambda wildcards, attempt: int(BWA_THREADS * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: int(threads * MEM_PER_CORE),
        time=WALL_TIME_MAX,
        time_min=1440,
        attempt=lambda wildcards, threads, attempt: attempt
    params:
        label=lambda wildcards: f"alignment.{wildcards.sample}.{wildcards.sample_subid}",
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u
        bwa mem -SP5M -t {threads} {input.ref} {input.r1} {input.r2} > $tmpdir/output.sam

        mv $tmpdir/output.sam {output.sam}.tmp
        mv {output.sam}.tmp {output.sam}
        rm -rf $tmpdir
        """


rule alignment_bam:
    input: "temp/{sample}/{sample_subid}.sam"
    output: "splits/{sample}/{sample_subid}.bam"
    threads: lambda wildcards, attempt: int(2 * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: int(threads * MEM_PER_CORE),
        time=WALL_TIME_MAX,
        time_min=1440,
        attempt=lambda wildcards, threads, attempt: attempt
    params:
        label=lambda wildcards: f"alignment_bam.{wildcards.sample}.{wildcards.sample_subid}",
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        samtools view -@ {threads} -b {input} -o $tmpdir/output.bam

        mv $tmpdir/output.bam {output}.tmp
        mv {output}.tmp {output}
        """


rule chimeric:
    input:
        norm_res="temp/{sample}/{sample_subid}.fastq.gz_norm0.txt.res.txt",
        sam="temp/{sample}/{sample_subid}.sam"
    output:
        frag=temp("temp/{sample}/{sample_subid}.fastq.gz.frag.txt"),
        norm_res="splits/{sample}/{sample_subid}.fastq.gz_norm.txt.res.txt",
        abnorm=temp("temp/{sample}/{sample_subid}.fastq.gz_abnormal.sam.gz"),
        unmapped=temp("temp/{sample}/{sample_subid}.fastq.gz_unmapped.sam.gz"),
    threads: lambda wildcards, attempt: int(2 * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: int(threads * MEM_PER_CORE),
        time=WALL_TIME_MAX,
        time_min=1440,
        attempt=lambda wildcards, threads, attempt: attempt
    params:
        label=lambda wildcards: f"chimeric.{wildcards.sample}.{wildcards.sample_subid}",
        juicer_home=juicer_home,
        site=SITE,
        site_file=SITE_FILE,
        frag=FRAG,
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        norm_reads=$tmpdir/norm.txt
        abnorm_sam=$tmpdir/abnormal.sam
        unmapped_sam=$tmpdir/unmapped.sam

        cp -av {input.norm_res} "$norm_reads".res.txt

        awk -v "fname1"=$norm_reads -v "fname2"=$abnorm_sam -v "fname3"=$unmapped_sam -f {params.juicer_home}/chimeric_blacklist.awk {input.sam}
        
        # if any normal reads were written, find what fragment they correspond to 
        # and store that
        if [ -e $norm_reads ] && [ "{params.site}" != "none" ] && [ -e "{params.site_file}" ]
        then
            {params.juicer_home}/fragment.pl $norm_reads $tmpdir/output.frag.txt {params.site_file}
        elif [ "{params.site}" == "none" ] || [ -z {params.frag} ] 
        then
            awk '{{printf("%s %s %s %d %s %s %s %d", $1, $2, $3, 0, $4, $5, $6, 1); for (i=7; i<=NF; i++) {{printf(" %s",$i);}}printf("\n");}}' $norm_reads > $tmpdir/output.frag.txt
        else                                                                    
            echo "***! No $norm_reads file created"
            exit 1
        fi

        mv $tmpdir/output.frag.txt {output.frag}.tmp
        pigz -p {threads} -c $abnorm_sam > {output.abnorm}.tmp
        pigz -p {threads} -c $unmapped_sam > {output.unmapped}.tmp

        cp -av "$norm_reads".res.txt {output.norm_res}
        mv {output.frag}.tmp {output.frag}
        mv {output.abnorm}.tmp {output.abnorm}
        mv {output.unmapped}.tmp {output.unmapped}

        rm -rf $tmpdir
        """


rule chimeric_sort:
    input: "temp/{sample}/{sample_subid}.fastq.gz.frag.txt",
    output: temp("temp/{sample}/{sample_subid}.sort.frag.txt"),
    threads: lambda wildcards, attempt: threads_for_mem(mem_for_attempt(SORT_MEM, attempt))
    resources:
        mem_mb=lambda wildcards, threads: int(threads * MEM_PER_CORE),
        time=WALL_TIME_MAX,
        time_min=1440,
        attempt=lambda wildcards, threads, attempt: attempt,
        sort_mem=lambda wildcards, threads, attempt: mem_for_attempt(SORT_MEM, attempt)
    params:
        label=lambda wildcards: f"chimeric_sort.{wildcards.sample}.{wildcards.sample_subid}",
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        sort -S {resources.sort_mem}M --parallel {threads} -T $tmpdir -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n {input} > $tmpdir/output.txt

        mv $tmpdir/output.txt {output}.tmp
        mv {output}.tmp {output}
        rm -rf $tmpdir
        """


rule frag_merge:
    input: lambda wildcards: expand("temp/{{sample}}/{sample_subid}.sort.frag.txt", sample_subid=get_sample_subid(wildcards.sample))
    output: temp("temp/{sample}.merged_sort.txt")
    threads: lambda wildcards, attempt: threads_for_mem(mem_for_attempt(SORT_MEM, attempt))
    resources:
        mem_mb=lambda wildcards, threads: int(threads * MEM_PER_CORE),
        time=WALL_TIME_MAX,
        time_min=1440,
        attempt=lambda wildcards, threads, attempt: attempt,
        sort_mem=lambda wildcards, threads, attempt: mem_for_attempt(SORT_MEM, attempt)
    params:
        label=lambda wildcards: f"frag_merge.{wildcards.sample}",
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        sort -S {resources.sort_mem}M --parallel {threads} -T $tmpdir -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n {input} > $tmpdir/merged_sort.txt

        mv $tmpdir/merged_sort.txt {output}.tmp
        mv {output}.tmp {output}
        rm -rf $tmpdir
        """


rule dedup:
    input: "temp/{sample}.merged_sort.txt"
    output:
        dups=temp("temp/{sample}.dups.txt"),
        optdups=temp("temp/{sample}.opt_dups.txt"),
        merged_nodups=temp("temp/{sample}.merged_nodups.txt"),
    threads: lambda wildcards, attempt: int(4 * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: int(threads * MEM_PER_CORE),
        time=WALL_TIME_MAX,
        time_min=1440,
        attempt=lambda wildcards, threads, attempt: attempt
    params:
        label=lambda wildcards: f"dedup.{wildcards.sample}",
        juicer_home=juicer_home,
        justexact=JUST_EXACT,
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        if [ -z {params.justexact} ]
        then
            awk -f {params.juicer_home}/dups.awk -v name=$tmpdir/ {input}
        else
            awk -f {params.juicer_home}/dups.awk -v name=$tmpdir/ -v nowobble=1 {input}
        fi

        # dups.awk may not generate these files. To make sure this rule doesn't fail, we need placeholders in this case.
        if [ -e $tmpdir/dups.txt ]; then mv $tmpdir/dups.txt {output.dups}.tmp; else touch {output.dups}.tmp; fi
        if [ -e $tmpdir/optdups.txt ]; then mv $tmpdir/optdups.txt {output.optdups}.tmp; else touch {output.optdups}.tmp; fi
        if [ -e $tmpdir/merged_nodups.txt ]; then mv $tmpdir/merged_nodups.txt {output.merged_nodups}.tmp; else touch {output.merged_nodups}.tmp; fi

        mv {output.dups}.tmp {output.dups}
        mv {output.optdups}.tmp {output.optdups}
        mv {output.merged_nodups}.tmp {output.merged_nodups}
        rm -rf $tmpdir
        """


rule stats_libcomplexity:
    input: 
        header="aligned/header",
        norm_res=lambda wildcards: expand("splits/{{sample}}/{sample_subid}.fastq.gz_norm.txt.res.txt", sample_subid=get_sample_subid(wildcards.sample)),
        dups="temp/{sample}.dups.txt",
        optdups="temp/{sample}.opt_dups.txt",
        merged_nodups="temp/{sample}.merged_nodups.txt",
    output: 
        inter=temp("temp/{sample}.inter.txt"),
    threads: lambda wildcards, attempt: int(3 * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: int(threads * MEM_PER_CORE),
        time=WALL_TIME_MAX,
        time_min=1440,
        attempt=lambda wildcards, threads, attempt: attempt
    params:
        label=lambda wildcards: f"stats_libcomplexity.{wildcards.sample}",
        java_mem=lambda wildcards, threads: int(threads * MEM_PER_CORE - 2000),
        juicer_home=juicer_home,
        site=SITE,
        site_file=SITE_FILE,
        ligation=ligation,
    shell:
        """
        export _JAVA_OPTIONS=-Xmx{params.java_mem}m

        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        ln -s $(readlink -f {input.merged_nodups}) $tmpdir/merged_nodups.txt
        ln -s $(readlink -f {input.dups}) $tmpdir/dups.txt
        ln -s $(readlink -f {input.optdups}) $tmpdir/opt_dups.txt

        tail -n1 {input.header} | awk '{{printf"%-1000s\\n", $0}}' > $tmpdir/inter.txt
        cat {input.norm_res} | awk -f {params.juicer_home}/stats_sub.awk >> $tmpdir/inter.txt

        java $_JAVA_OPTIONS -jar {params.juicer_home}/juicer_tools.jar LibraryComplexity $tmpdir inter.txt >> $tmpdir/inter.txt 

        mv $tmpdir/inter.txt {output.inter}

        rm -rf $tmpdir
        """


rule stats_inter:
    input: 
        header="aligned/header",
        merged_nodups="temp/{sample}.merged_nodups.txt",
        inter="temp/{sample}.inter.txt",
    output: 
        inter="aligned/{sample}.{inter_type,(inter|inter_30)}.txt",
        hist="aligned/{sample}.{inter_type,(inter|inter_30)}_hists.m",
    threads: lambda wildcards, attempt: int(3 * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: int(threads * MEM_PER_CORE),
        time=WALL_TIME_MAX,
        time_min=1440,
        attempt=lambda wildcards, threads, attempt: attempt
    params:
        label=lambda wildcards: f"stats_inter.{wildcards.sample}.{wildcards.inter_type}",
        juicer_home=juicer_home,
        site=SITE,
        site_file=SITE_FILE,
        ligation=ligation,
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        ln -s $(readlink -f {input.merged_nodups}) $tmpdir/merged_nodups.txt
        cp -av {input.inter} $tmpdir/inter.txt

        if [ "{wildcards.inter_type}" == "inter_30" ]
        then
            mapq=30
        else
            mapq=1
        fi

        {params.juicer_home}/statistics.pl -s {params.site_file} -l {params.ligation} -o $tmpdir/inter.txt -q $mapq $tmpdir/merged_nodups.txt

        mv $tmpdir/inter.txt {output.inter}
        mv $tmpdir/inter_hists.m {output.hist}

        rm -rf $tmpdir
        """


rule stats_collisions:
    input: 
        abnorm="aligned/{sample}.fastq.gz_abnormal.sam.gz",
    output: 
        collisions="aligned/{sample}.collisions.txt",
        collisions_dups="aligned/{sample}.collisions_dups.txt",
        collisions_nodups="aligned/{sample}.collisions_nodups.txt",
    threads: lambda wildcards, attempt: int(4 * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: int(threads * MEM_PER_CORE),
        time=WALL_TIME_MAX,
        time_min=1440,
        attempt=lambda wildcards, threads, attempt: attempt
    params:
        label=lambda wildcards: f"stats_collisions.{wildcards.sample}",
        java_mem=lambda wildcards, threads: int(threads * MEM_PER_CORE - 2000),
        sort_mem=lambda wildcards, threads: int(threads * MEM_PER_CORE - 2000),
        juicer_home=juicer_home,
        site=SITE,
        site_file=SITE_FILE,
        ligation=ligation,
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        pigz -p {threads} -dc {input.abnorm} | awk -f {params.juicer_home}/collisions.awk > $tmpdir/collisions.txt

        # Collisions dedupping: two pass algorithm, ideally would make one pass
        gawk -v fname=$tmpdir/collisions.txt -f {params.juicer_home}/collisions_dedup_rearrange_cols.awk $tmpdir/collisions.txt | \
        sort -S {params.sort_mem}M --parallel {threads} -T $tmpdir -k3,3n -k4,4n -k10,10n -k11,11n -k17,17n -k18,18n -k24,24n -k25,25n -k31,31n -k32,32n | \
        awk -v name=$tmpdir/ -f {params.juicer_home}/collisions_dups.awk
        
        ls -lah $tmpdir
        mv $tmpdir/collisions.txt {output.collisions}.tmp
        mv $tmpdir/collisions_dups.txt {output.collisions_dups}.tmp
        mv $tmpdir/collisions_nodups.txt {output.collisions_nodups}.tmp
        mv {output.collisions}.tmp {output.collisions}
        mv {output.collisions_dups}.tmp {output.collisions_dups}
        mv {output.collisions_nodups}.tmp {output.collisions_nodups}

        rm -rf $tmpdir
        """


rule merge_abnorm_unmapped:
    input: 
        abnorm=lambda wildcards: expand("temp/{{sample}}/{sample_subid}.fastq.gz_abnormal.sam.gz", sample_subid=get_sample_subid(wildcards.sample)),
        unmapped=lambda wildcards: expand("temp/{{sample}}/{sample_subid}.fastq.gz_unmapped.sam.gz", sample_subid=get_sample_subid(wildcards.sample)),
    output: 
        abnorm="aligned/{sample}.fastq.gz_abnormal.sam.gz",
        unmapped="aligned/{sample}.fastq.gz_unmapped.sam.gz",
    threads: lambda wildcards, attempt: int(1 * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: int(threads * MEM_PER_CORE),
        time=WALL_TIME_MAX,
        time_min=1440,
        attempt=lambda wildcards, threads, attempt: attempt
    params:
        label=lambda wildcards: f"merge_abnorm_unmapped.{wildcards.sample}",
    shell:
        """
        cat {input.abnorm} > {output.abnorm}.tmp
        cat {input.unmapped} > {output.unmapped}.tmp
        mv {output.abnorm}.tmp {output.abnorm}
        mv {output.unmapped}.tmp {output.unmapped}
        """


rule fragments_gz:
    input: 
        dups="temp/{sample}.dups.txt",
        optdups="temp/{sample}.opt_dups.txt",
        merged_nodups="temp/{sample}.merged_nodups.txt",
    output: 
        dups="aligned/{sample}.dups.txt.gz",
        optdups="aligned/{sample}.opt_dups.txt.gz",
        merged_nodups="aligned/{sample}.merged_nodups.txt.gz",
    threads: lambda wildcards, attempt: int(8 * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: int(threads * MEM_PER_CORE),
        time=WALL_TIME_MAX,
        time_min=1440,
        attempt=lambda wildcards, threads, attempt: attempt
    params:
        label=lambda wildcards: f"merged_nodups_gz.{wildcards.sample}",
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        pigz -p {threads} -c {input.dups} > $tmpdir/dups.txt.gz
        pigz -p {threads} -c {input.optdups} > $tmpdir/opt_dups.txt.gz
        pigz -p {threads} -c {input.merged_nodups} > $tmpdir/merged_nodups.txt.gz

        mv $tmpdir/dups.txt.gz {output.dups}.tmp
        mv $tmpdir/opt_dups.txt.gz {output.optdups}.tmp
        mv $tmpdir/merged_nodups.txt.gz {output.merged_nodups}.tmp
        mv {output.dups}.tmp {output.dups}
        mv {output.optdups}.tmp {output.optdups}
        mv {output.merged_nodups}.tmp {output.merged_nodups}

        rm -rf $tmpdir
        """


rule hic:
    input: 
        inter="aligned/{sample}.{inter_type}.txt",
        hist="aligned/{sample}.{inter_type}_hists.m",
        merged_nodups="temp/{sample}.merged_nodups.txt",
    output: "aligned/{sample}.{inter_type,(inter|inter_30)}.hic"
    threads: lambda wildcards, attempt: threads_for_mem(mem_for_attempt(JUICER_PRE_MEM, attempt))
    resources:
        mem_mb=lambda wildcards, threads: int(threads * MEM_PER_CORE),
        time=WALL_TIME_MAX,
        time_min=1440,
        attempt=lambda wildcards, threads, attempt: attempt,
        java_mem=lambda wildcards, threads, attempt: mem_for_attempt(JUICER_PRE_MEM, attempt),
        hic_threads=lambda wildcards, threads: min(8, threads)
    params:
        label=lambda wildcards: f"hic.{wildcards.sample}.{wildcards.inter_type}",
        frag=FRAG,
        juicer_home=juicer_home,
        site=SITE,
        site_file=SITE_FILE,
        genome_path=GENOME_PATH,
        res=config.get("JUICER_PRE_RES", "")
    shell:
        """
        export _JAVA_OPTIONS="-Xmx{resources.java_mem}m -Xms{resources.java_mem}m"

        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        if [ "{wildcards.inter_type}" == "inter_30" ]
        then
            mapq=30
        else
            mapq=1
        fi

        if [ -z "{params.res}" ]
        then
            res=""
        else
            res="-r {params.res}"
        fi

        if [ -z {params.frag} ]
        then
            java $_JAVA_OPTIONS -jar {params.juicer_home}/juicer_tools.jar pre -s {input.inter} -g {input.hist} -q $mapq -t $tmpdir --threads {resources.hic_threads} $res \
            {input.merged_nodups} $tmpdir/inter.hic {params.genome_path}
        else
            java $_JAVA_OPTIONS -jar {params.juicer_home}/juicer_tools.jar pre -f {params.site_file} -s {input.inter} -g {input.hist} -q $mapq -t $tmpdir --threads {resources.hic_threads} $res \
            {input.merged_nodups} $tmpdir/inter.hic {params.genome_path}
        fi

        mv $tmpdir/inter.hic {output}.tmp
        mv {output}.tmp {output}

        rm -rf $tmpdir
        """


rule arrowhead:
    input: 
        "aligned/{sample}.{inter_type}.hic"
    output:
        "aligned/{sample}.{inter_type,(inter|inter_30)}.domain.bedpe"
    threads: lambda wildcards, attempt: threads_for_mem(mem_for_attempt(JUICER_ARROWHEAD_MEM, attempt))
    resources:
        mem_mb=lambda wildcards, threads: int(threads * MEM_PER_CORE),
        time=WALL_TIME_MAX,
        time_min=1440,
        attempt=lambda wildcards, threads, attempt: attempt,
        java_mem=lambda wildcards, threads, attempt: mem_for_attempt(JUICER_ARROWHEAD_MEM, attempt),
    params:
        label=lambda wildcards: f"arrowhead.{wildcards.sample}.{wildcards.inter_type}",
        juicer_home=juicer_home,
    shell:
        """
        export _JAVA_OPTIONS="-Xmx{resources.java_mem}m -Xms{resources.java_mem}m"

        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        java $_JAVA_OPTIONS -jar {params.juicer_home}/juicer_tools.jar arrowhead -r 5000 -k KR --threads {threads} {input} $tmpdir
        mv $tmpdir/5000_blocks.bedpe {output}
        """


rule hiccups:
    input: 
        "aligned/{sample}.{inter_type}.hic"
    output:
        "aligned/{sample}.{inter_type,(inter|inter_30)}.hiccups.zip"
    threads: lambda wildcards, attempt: threads_for_mem(mem_for_attempt(JUICER_HICCUPS_MEM, attempt))
    resources:
        mem_mb=lambda wildcards, threads: int(threads * MEM_PER_CORE),
        time=WALL_TIME_MAX,
        time_min=1440,
        attempt=lambda wildcards, threads, attempt: attempt,
        java_mem=lambda wildcards, threads, attempt: mem_for_attempt(JUICER_HICCUPS_MEM, attempt),
    params:
        label=lambda wildcards: f"hiccups.{wildcards.sample}.{wildcards.inter_type}",
        juicer_home=juicer_home,
    shell:
        """
        export _JAVA_OPTIONS="-Xmx{resources.java_mem}m -Xms{resources.java_mem}m"

        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        java $_JAVA_OPTIONS -jar {params.juicer_home}/juicer_tools.jar hiccups --cpu --threads {threads} -m 512 \
        -r 5000,10000,25000 -k KR -f .1,.1,.1 -p 4,2,1 -i 7,5,3 -t 0.02,1.5,1.75,2 -d 20000,20000,50000 \
        {input} $tmpdir

        pushd $tmpdir
        zip -r hiccups.zip *
        popd

        mv $tmpdir/hiccups.zip {output}
        """