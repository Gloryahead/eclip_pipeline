#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * scATAC-seq Pipeline — FASTQ → annotated cell types
 *
 * Tutorial : https://ngs101.com/how-to-analyze-single-cell-atac-seq-data
 * Dataset  : GSE282769 — COVID-19 PBMCs (2 severe, 2 healthy)
 * Reference: GRCh38-2024-A (10x Genomics Cell Ranger ARC)
 *
 * Stages
 *   1  SRA download → read-length detection → Cell Ranger ATAC → MACS3 peaks
 *   2  Per-sample Signac QC + AMULET doublet removal
 *   3  Consensus peak set → TF-IDF / LSI → Harmony → clustering
 *   4  Gene activity + label transfer → cell type annotation + motifs
 *
 * Usage (on HPC)
 *   nextflow run /home/u11/maarowosegbe/scATAC-seq_pipeline/main.nf \
 *       -c /home/u11/maarowosegbe/scATAC-seq_pipeline/nextflow.config \
 *       -profile slurm \
 *       -resume
 */

include { SRA_DOWNLOAD          } from './modules/local/sra_download'
include { RENAME_FASTQS         } from './modules/local/rename_fastqs'
include { CELLRANGER_ATAC       } from './modules/local/cellranger_atac'
include { MACS3_CALLPEAK        } from './modules/local/macs3'
include { TABIX_INDEX           } from './modules/local/tabix_index'
include { SIGNAC_QC             } from './modules/local/signac_qc'
include { CONSENSUS_PEAKS       } from './modules/local/consensus_peaks'
include { INTEGRATION           } from './modules/local/integration'
include { CELL_TYPE_ANNOTATION  } from './modules/local/cell_type_annotation'

workflow {

    // ── Load sample sheet (samples.csv → [sample_id, srr, condition]) ──────
    Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true, strip: true)
        .map { row -> tuple(row.sample_id, row.srr_accession, row.condition) }
        .set { samples_ch }

    // ── Stage 1: Preprocessing ─────────────────────────────────────────────
    SRA_DOWNLOAD(samples_ch)

    RENAME_FASTQS(SRA_DOWNLOAD.out.raw_fastqs)

    CELLRANGER_ATAC(
        RENAME_FASTQS.out.named_fastqs,
        params.cellranger_ref
    )

    MACS3_CALLPEAK(CELLRANGER_ATAC.out.fragments)

    TABIX_INDEX(CELLRANGER_ATAC.out.fragments)

    // ── Stage 2: Per-sample QC ──────────────────────────────────────────────
    // Join Cell Ranger outputs with tabix-indexed fragments on sample_id
    CELLRANGER_ATAC.out.cr_outs
        .join(TABIX_INDEX.out.indexed_frags, by: 0)
        .set { qc_input_ch }

    SIGNAC_QC(
        qc_input_ch,
        params.blacklist_bed,
        params.ensdb_cache
    )

    // ── Stage 3: Integration + clustering ──────────────────────────────────
    // Collect per-sample outputs into single lists for multi-sample steps
    CONSENSUS_PEAKS(
        SIGNAC_QC.out.qc_rds.collect { it[1] },
        MACS3_CALLPEAK.out.peak_bed.collect { it[1] },
        CELLRANGER_ATAC.out.fragments.collect { it[1] },
        TABIX_INDEX.out.indexed_frags.collect { it[1] },
        params.blacklist_bed
    )

    INTEGRATION(CONSENSUS_PEAKS.out.merged_rds)

    // ── Stage 4: Cell type annotation ──────────────────────────────────────
    CELL_TYPE_ANNOTATION(
        INTEGRATION.out.clustered_rds,
        params.rna_ref_rds
    )

    // Print final outputs
    CELL_TYPE_ANNOTATION.out.annotated_rds
        | view { rds -> "Pipeline complete. Annotated object: ${rds}" }
}
