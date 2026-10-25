#!/usr/bin/env nextflow

/*
 * mappability
 *
 * Empirical paired-end mappability of a set of loci (by default every
 * protein-coding gene) with UniMap.jl: tile the loci with error-free read
 * pairs, align them genome-wide, and score how often the fragments covering
 * each base map back to where they came from.
 *
 *   GFF + FASTA
 *        |
 *    BUILD_LOCI            gene IDs -> loci.tsv (gene span + flanks)
 *        |
 *    SIMULATE_AND_BATCH    unimap simulate; per-region FASTQs packed into
 *        |                 n_batches balanced batches per read length
 *        |---- sim_meta/ ---------------------------+
 *        |                                          |
 *    ALIGN_BATCH x (n_batches x read lengths)       |
 *        |    bowtie2 -> unsorted shard BAM         |
 *        |                                          |
 *    CAT_AND_SCORE  <-------------------------------+
 *             samtools cat per L -> unimap score -> tar mappability.zarr
 *
 * Nearly all the compute is bowtie2, and sharding it is exact: each pair is
 * aligned independently and carries its origin in its read name, so the shards
 * concatenated in any order score identically to one monolithic alignment.
 * What must never happen is sorting or filtering the BAMs -- `unimap score`
 * needs every pair, mates adjacent. See modules/local/ for the details of each
 * step and the checks that guard it.
 */

nextflow.enable.dsl = 2

include { SAMTOOLS_FAIDX     } from './modules/local/samtools_faidx'
include { BOWTIE2_BUILD      } from './modules/local/bowtie2_build'
include { BUILD_LOCI         } from './modules/local/build_loci'
include { SIMULATE_AND_BATCH } from './modules/local/simulate_and_batch'
include { ALIGN_BATCH        } from './modules/local/align_batch'
include { CAT_AND_SCORE      } from './modules/local/cat_and_score'

// ---------------------------------------------------------------------------

def helpMessage() {
    log.info """
    mappability ${workflow.manifest.version}

      nextflow run . -profile local --unimap_project ../UniMap.jl \\
          --fasta ref.fna --gff ref.gff --bowtie2_index ref/bowtie2/ref
      nextflow run . -profile unity  --fasta ... --gff ... --bowtie2_index ...
      nextflow run . -profile test,local --unimap_project ../UniMap.jl

    Input
      --fasta             reference FASTA, uncompressed; must be the one the
                          bowtie2 index was built from. A .fai beside it is reused.
      --gff               GFF3 (or GTF with --genes), optionally .gz
      --bowtie2_index     bowtie2 index prefix, e.g. ref/bowtie2/GRCz12tu. Omit
                          to build one from --fasta (small genomes only).
      --outdir            results directory                     [${params.outdir}]

    Loci
      --gene_biotype      take every GFF3 gene with this gene_biotype
                                                                [${params.gene_biotype}]
      --genes             or a file of IDs, one per line        [${params.genes ?: 'none'}]
      --upstream          flank upstream of the TSS (bp)        [${params.upstream}]
      --downstream        flank downstream of the TES (bp)      [${params.downstream}]

    Simulation grid
      --read_lengths      [${params.read_lengths}]
      --fragment_lengths  [${params.fragment_lengths}]
                          bowtie2 -X is set to the largest of these
      --step              bp between fragment starts            [${params.step}]
      --orientation       plus, minus or both                   [${params.orientation}]

    Sharding and alignment
      --n_batches         alignment batches per read length     [${params.n_batches}]
      --bowtie2_args      extra bowtie2 arguments               [${params.bowtie2_args ?: 'none'}]

    Scoring
      --cutoffs           fragment MAPQ cutoffs                 [${params.cutoffs}]
      --frag_mapq         min or max of the two mates' MAPQ     [${params.frag_mapq}]
      --tolerance         bp allowed from the origin            [${params.tolerance}]

    Output
      --tar_zarr          publish mappability.zarr as one tar   [${params.tar_zarr}]
      --publish_bams      also publish the per-L BAMs (large)   [${params.publish_bams}]

    Software
      --unimap_project    UniMap.jl checkout, when running without containers
                                                                [${params.unimap_project ?: 'container'}]
      --julia_threads     julia threads for loci/simulate/score  [task cpus]
                          1 only to run the amd64 image under emulation

    Profiles: local, docker, apptainer, unity, aws, seqera, test
    """.stripIndent()
}

def parseInts(String name, value) {
    def items = value.toString().tokenize(',').collect { it.trim() }
    if (!items || !items.every { it ==~ /\d+/ }) {
        error "--${name} must be a comma-separated list of non-negative integers, got '${value}'"
    }
    def ints = items.collect { it as int }
    if (ints.unique(false).size() != ints.size()) {
        error "--${name} has repeated values: '${value}'"
    }
    return ints
}

/*
 * Arguments the pipeline sets itself, or that would make the BAMs unscorable.
 * -X is derived from the fragment grid: a fragment longer than -X cannot align
 * as a proper pair, so a mismatch silently scores the longest fragments as
 * unmappable.
 */
def checkBowtie2Args(String args) {
    def forbidden = ['--no-unal', '-p', '--threads', '-x', '-1', '-2', '-U', '-S', '-b',
                     '-X', '--maxins', '--un', '--al', '--un-conc', '--al-conc',
                     '--un-gz', '--al-gz', '--un-conc-gz', '--al-conc-gz']
    def bad = args.tokenize(' ').collect { it.tokenize('=')[0] }.findAll { it in forbidden }
    if (bad) {
        error "--bowtie2_args must not set ${bad.unique().join(', ')}: threads, index, inputs " +
              "and output are handled by the pipeline, -X is the largest fragment length, " +
              "and every pair must reach the BAM."
    }
}

/*
 * The six files of a bowtie2 index, small (.bt2) or large (.bt2l). Resolved
 * through file(), so a local path, a path relative to the launch directory and
 * an s3:// prefix all work.
 */
def bowtie2IndexFiles(String prefix) {
    def suffixes = ['1', '2', '3', '4', 'rev.1', 'rev.2']
    for (ext in ['bt2', 'bt2l']) {
        def files = suffixes.collect { file("${prefix}.${it}.${ext}") }
        if (files.every { it.exists() }) { return files }
    }
    error "--bowtie2_index: no complete bowtie2 index at prefix '${prefix}' " +
          "(expected ${prefix}.{1,2,3,4,rev.1,rev.2}.bt2 or .bt2l)"
}

def validateParams() {
    def missing = ['fasta', 'gff'].findAll { !params[it] }
    if (missing) {
        helpMessage()
        error "missing required parameter(s): ${missing.collect { '--' + it }.join(', ')}"
    }
    if (params.fasta.toString().endsWith('.gz')) {
        error "--fasta must be uncompressed: UniMap reads it through a samtools .fai index"
    }
    if (!params.genes && params.gff.toString() =~ /\.gtf(\.gz)?$/) {
        error "--gene_biotype selection reads GFF3 attributes; with a GTF give --genes"
    }
    def Ls = parseInts('read_lengths', params.read_lengths)
    def Fs = parseInts('fragment_lengths', params.fragment_lengths)
    if ((Ls + Fs).any { it < 1 }) {
        error "--read_lengths and --fragment_lengths must be positive"
    }
    parseInts('cutoffs', params.cutoffs)
    if (!(params.step.toString() ==~ /\d+/) || (params.step as int) < 1) {
        error "--step must be a positive integer"
    }
    if ((params.step as int) > Ls.min().intdiv(5)) {
        log.warn "--step ${params.step} is more than a fifth of the shortest read length " +
                 "(${Ls.min()}); per-base tracks get coarse at that spacing"
    }
    if (!(params.orientation in ['plus', 'minus', 'both'])) {
        error "--orientation must be plus, minus or both"
    }
    if (params.julia_threads && !(params.julia_threads.toString() ==~ /[1-9]\d*/)) {
        error "--julia_threads must be a positive integer"
    }
    if (!(params.frag_mapq in ['min', 'max'])) {
        error "--frag_mapq must be min or max"
    }
    if (!(params.n_batches.toString() ==~ /\d+/) || !((params.n_batches as int) in 1..999)) {
        error "--n_batches must be between 1 and 999"
    }
    checkBowtie2Args(params.bowtie2_args ?: '')
    if ('aws' in workflow.profile.tokenize(',') && !params.aws_queue) {
        error "-profile aws needs --aws_queue (under Seqera use -profile seqera instead)"
    }
    if (!workflow.containerEngine && !params.unimap_project) {
        error "no container engine is enabled, so UniMap must come from a local checkout: " +
              "pass --unimap_project <UniMap.jl dir> (or export UNIMAP_PROJECT), or use " +
              "-profile docker / apptainer"
    }
    if (params.unimap_project && !workflow.containerEngine &&
            !file("${params.unimap_project}/Project.toml").exists()) {
        error "--unimap_project ${params.unimap_project} is not a Julia project (no Project.toml)"
    }
    return Fs.max()
}

// ---------------------------------------------------------------------------

workflow {
    if (params.help) {
        helpMessage()
        return
    }

    def fmax = validateParams()

    def fasta = file(params.fasta, checkIfExists: true)
    def gff   = file(params.gff,   checkIfExists: true)
    def fai   = file("${params.fasta}.fai")
    def genes = params.genes ? file(params.genes, checkIfExists: true)
                             : file("${projectDir}/assets/NO_GENES")
    def index_files  = params.bowtie2_index ? bowtie2IndexFiles(params.bowtie2_index.toString()) : null
    def index_name   = index_files ? index_files[0].name.replaceFirst(/\.1\.bt2l?$/, '') : fasta.baseName
    def aligner_args = (["-X ${fmax}"] + (params.bowtie2_args ? [params.bowtie2_args.toString().trim()] : [])).join(' ')

    log.info "mappability ${workflow.manifest.version}: L=${params.read_lengths} " +
             "F=${params.fragment_lengths} step=${params.step} orientation=${params.orientation} " +
             "batches<=${params.n_batches} bowtie2 ${aligner_args}"
    if (!index_files) {
        log.warn "no --bowtie2_index: building one from ${fasta.name}. Fine for a test genome; " +
                 "for a large one build it once and pass --bowtie2_index."
    }

    fai_ch   = fai.exists() ? Channel.value(fai) : SAMTOOLS_FAIDX(fasta).fai
    index_ch = index_files ? Channel.value(index_files) : BOWTIE2_BUILD(fasta, index_name).index

    BUILD_LOCI(fasta, fai_ch, gff, genes)
    SIMULATE_AND_BATCH(fasta, fai_ch, BUILD_LOCI.out.loci)

    // batch_L{L}_{batch}_R{1,2}.fq.gz -> [[L, batch], [R1, R2]]
    reads = SIMULATE_AND_BATCH.out.fastq
        .flatten()
        .map { f ->
            def m = f.name =~ /^batch_L(\d+)_(\d+)_R[12]\.fq\.gz$/
            if (!m) { error "unexpected batch file name: ${f.name}" }
            [ [ m[0][1] as int, m[0][2] ], f ]
        }
        .groupTuple(size: 2)
        .map { key, fq -> [ key, fq.sort { it.name } ] }

    // The pairs each batch's FASTQ must hold, from the batch plan.
    expected = SIMULATE_AND_BATCH.out.manifest
        .splitCsv(sep: '\t', header: true)
        .map { row -> [ [ row.L as int, row.batch ], row.n_pairs as long ] }

    batches = reads
        .join(expected, failOnMismatch: true, failOnDuplicate: true)
        .map { key, fq, n ->
            if (!(fq[0].name.endsWith('_R1.fq.gz') && fq[1].name.endsWith('_R2.fq.gz'))) {
                error "batch ${key}: mates out of order: ${fq*.name}"
            }
            [ key[0], key[1], fq[0], fq[1], n ]
        }

    ALIGN_BATCH(batches, index_ch, index_name, aligner_args)

    CAT_AND_SCORE(
        SIMULATE_AND_BATCH.out.meta,
        SIMULATE_AND_BATCH.out.manifest,
        ALIGN_BATCH.out.bam.collect(),
        ALIGN_BATCH.out.version.collect(),
        fasta,
        fai_ch,
        index_name,
        aligner_args
    )
}
