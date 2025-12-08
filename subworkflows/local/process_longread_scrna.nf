//
// Performs alignment
//

// SUBWORKFLOWS
include { ALIGN_LONGREADS         } from '../../subworkflows/local/align_longreads'
include { QUANTIFY_SCRNA_ISOQUANT } from '../../subworkflows/local/quantify_scrna_isoquant'
include { QUANTIFY_SCRNA_OARFISH  } from '../../subworkflows/local/quantify_scrna_oarfish'
include { DEDUP_UMIS              } from '../../subworkflows/local/dedup_umis'

// MODULES
include { PICARD_MARKDUPLICATES                         } from '../../modules/nf-core/picard/markduplicates'
include { PICARD_COLLECTRNASEQMETRICS                   } from '../../modules/nf-core/picard/collectrnaseqmetrics'                                                                                           
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_TAGGED } from '../../modules/nf-core/samtools/flagstat'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_DEDUP  } from '../../modules/nf-core/samtools/flagstat'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TAGGED       } from '../../modules/nf-core/samtools/index'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUP        } from '../../modules/nf-core/samtools/index'
include { SAMTOOLS_VIEW as SAMTOOLS_FILTER_DEDUP        } from '../../modules/nf-core/samtools/view'

include { TAG_BARCODES                  } from '../../modules/local/tag_barcodes'
include { ADD_TAGS as ADD_TAGS_PARSE    } from '../../modules/local/add_tags/main'
include { ADD_TAGS as ADD_TAGS_ARGENTAG } from '../../modules/local/add_tags/main'
include { RRNA_FOR_COLLECTRNASEQMETRICS } from '../../Bionextflow3/modules/local/rrna_for_collectrnaseqmetrics'

workflow PROCESS_LONGREAD_SCRNA {
    take:
        fasta           // channel: [ val(meta), path(fasta) ]
        fai             // channel: [ val(meta), path(fai) ]
        gtf             // channel: [ val(meta), path(gtf) ]
        fastq           // channel: [ val(meta), path(fastq) ]
        picard_refflat  // channel: [ val(meta), path(picard_refflat) ]
        rseqc_bed       // channel: [ val(meta), path(rseqc_bed) ]
        read_bc_info    // channel: [ val(meta), path(read_barcode_info) ]
        quant_list      // list: List of quantifiers to use
        dedup_tool      // str: Name of deduplication tool to use
        genome_aligned  // bool: Whether the bam is aligned to the genome or not
        fasta_delimiter // str: Delimiter character used in sequence id in fasta

        platform // str: Technology used to generate the datasets

        skip_save_minimap2_index // bool: Skip saving the minimap2 index
        skip_qc                  // bool: Skip qc steps
        skip_rseqc               // bool: Skip RSeQC
        skip_bam_nanocomp        // bool: Skip Nanocomp
        skip_seurat              // bool: Skip seurat qc steps
        skip_dedup               // bool: Skip deduplication

    main:
        ch_versions = Channel.empty()

        //
        // SUBWORKFLOW: Align long Read Data
        //

        ALIGN_LONGREADS(
            fasta,
            fai,
            gtf,
            fastq,
            rseqc_bed,
            skip_save_minimap2_index,
            skip_qc,
            skip_rseqc,
            skip_bam_nanocomp
        )
        ch_versions = ch_versions.mix(ALIGN_LONGREADS.out.versions)

        //
        // MODULE: Tag Barcodes
        //
        if (platform == '10X') {
            TAG_BARCODES (
            ALIGN_LONGREADS.out.sorted_bam
                .join( ALIGN_LONGREADS.out.sorted_bai, by: 0 )
                .join( read_bc_info, by: 0)
            ).set {ch_tagged_bam}
            ch_versions = ch_versions.mix(TAG_BARCODES.out.versions)

        } else if (platform == 'Parse') {
            ADD_TAGS_PARSE(ALIGN_LONGREADS.out.sorted_bam).set {ch_tagged_bam}
            ch_versions = ch_versions.mix(ADD_TAGS_ARGENTAG.out.versions)

        } else if (platform == 'Argentag') {
            ADD_TAGS_ARGENTAG(ALIGN_LONGREADS.out.sorted_bam).set {ch_tagged_bam}
            ch_versions = ch_versions.mix(ADD_TAGS_ARGENTAG.out.versions)

        } else {
            exit 1, "Single cell platform not recognized. You can choose either 10X, Parse or Argentag.\n"
        }
        
        //
        // MODULE: Index Tagged Bam
        //
        SAMTOOLS_INDEX_TAGGED ( ch_tagged_bam.tagged_bam )
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_TAGGED.out.versions)

        //
        // MODULE: Flagstat Tagged Bam
        //
        SAMTOOLS_FLAGSTAT_TAGGED (
            ch_tagged_bam.tagged_bam
                .join( SAMTOOLS_INDEX_TAGGED.out.bai, by: [0])
        )
        ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT_TAGGED.out.versions)

        ch_bam = Channel.empty()
        ch_bai = Channel.empty()
        ch_flagstat = Channel.empty()
        ch_idxstats = Channel.empty()

        if (!skip_dedup) {
            DEDUP_UMIS (
                fasta,
                fai,
                gtf,
                ch_tagged_bam.tagged_bam,
                SAMTOOLS_INDEX_TAGGED.out.bai,
                true, // Used to split the bam
                genome_aligned,
                dedup_tool,
                fasta_delimiter
            )

            ch_bam = DEDUP_UMIS.out.dedup_bam
            ch_bai = DEDUP_UMIS.out.dedup_bai
            ch_flagstat = DEDUP_UMIS.out.dedup_flagstat
            ch_versions = DEDUP_UMIS.out.versions

            // QC: base proportion belonging to different kind of features:
            rrnas = RRNA_FOR_COLLECTRNASEQMETRICS(fai, gtf).rRRNAs.map{ it[1] }
            rnaseq_metrics = PICARD_COLLECTRNASEQMETRICS(ch_bam, picard_refflat, fasta.map{ it[1] }, rrnas)

        } else {

            ch_bam = ch_tagged_bam.tagged_bam
            ch_bai = SAMTOOLS_INDEX_TAGGED.out.bai
            ch_flagstat = SAMTOOLS_FLAGSTAT_TAGGED.out.flagstat
                .map{
                    meta, flagstat ->
                        id = ['id': meta.id]
                    [id, flagstat]
                }

        }
        //
        // SUBWORKFLOW: Quantify Features
        //

        ch_gene_qc_stats = Channel.empty()
        ch_transcript_qc_stats = Channel.empty()

        if (quant_list.contains("oarfish")) {
            QUANTIFY_SCRNA_OARFISH (
                ch_bam,
                ch_bai,
                ch_flagstat,
                fasta,
                skip_qc,
                skip_seurat
            )
            ch_versions = ch_versions.mix(QUANTIFY_SCRNA_OARFISH.out.versions)
            ch_transcript_qc_stats = QUANTIFY_SCRNA_OARFISH.out.transcript_qc_stats
        }

        if (quant_list.contains("isoquant")) {
            QUANTIFY_SCRNA_ISOQUANT (
                ch_bam,
                ch_bai,
                ch_flagstat,
                fasta,
                fai,
                gtf,
                skip_qc,
                skip_seurat
            )

            ch_versions = ch_versions.mix(QUANTIFY_SCRNA_ISOQUANT.out.versions)
            ch_gene_qc_stats = QUANTIFY_SCRNA_ISOQUANT.out.gene_qc_stats
            ch_transcript_qc_stats = QUANTIFY_SCRNA_ISOQUANT.out.transcript_qc_stats
        }
        
    emit:
        // Versions
        versions                 = ch_versions

        // Minimap results + qc's
        minimap_bam              = ALIGN_LONGREADS.out.sorted_bam
        minimap_bai              = ALIGN_LONGREADS.out.sorted_bai
        minimap_stats            = ALIGN_LONGREADS.out.stats
        minimap_flagstat         = ALIGN_LONGREADS.out.flagstat
        minimap_idxstats         = ALIGN_LONGREADS.out.idxstats
        minimap_rseqc_read_dist  = ALIGN_LONGREADS.out.rseqc_read_dist
        minimap_nanocomp_bam_txt = ALIGN_LONGREADS.out.nanocomp_bam_txt

        // Barcode tagging results + qc's
        //bc_tagged_bam            = TAG_BARCODES.out.tagged_bam
        bc_tagged_bam            = ch_tagged_bam.tagged_bam
        bc_tagged_bai            = SAMTOOLS_INDEX_TAGGED.out.bai
        bc_tagged_flagstat       = SAMTOOLS_FLAGSTAT_TAGGED.out.flagstat

        // Deduplication results
        dedup_bam                = ch_bam
        dedup_bai                = ch_bai
        dedup_flagstat           = ch_flagstat
        dedup_idxstats           = ch_idxstats

        // Picard QC Stats
        picard_metrics           = rnaseq_metrics.metrics

        // Seurat QC Stats
        gene_qc_stats            = ch_gene_qc_stats
        transcript_qc_stats      = ch_transcript_qc_stats
    
}
