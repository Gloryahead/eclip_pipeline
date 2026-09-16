# eCLIP-seq Analysis Pipeline
## Rhine et al. 2025, *Nature Neuroscience*
> "Neuronal aging causes mislocalization of splicing proteins and unchecked cellular stress"
> doi: [10.1038/s41593-025-01952-z](https://www.nature.com/articles/s41593-025-01952-z)

Reproduces the eCLIP-seq analysis from Rhine et al. 2025 using the **Skipper** pipeline
(Boyle et al. 2023, *Cell Genomics*). Designed for UArizona HPC (Puma cluster, SLURM),
but the core Skipper+Snakemake layer runs on any SLURM cluster with minor path edits.

---

## What this pipeline does

| Experiment | Targets | Samples | Genome |
|---|---|---|---|
| Human cell eCLIP | TDP-43, G3BP1, Caprin1 | Tdiff neurons vs iPSC-diff neurons | hg38 (GENCODE v47) |
| Human brain eCLIP | TDP-43, G3BP1 | Frontal cortex: mid-age (30–50 yr) vs old-age (80–90 yr) | hg38 (GENCODE v47) |
| Mouse eCLIP | TDP-43 | Cerebellum: 1.5, 6, 24 months | mm10 (GENCODE M35) |

**Pipeline steps (Skipper):**
1. Adapter trimming (cutadapt)
2. STAR alignment to hg38 or mm10
3. PCR deduplication using 10 nt UMIs (umi_tools) — critical for eCLIP, see note below
4. Statistical peak calling per genomic window vs size-matched input (Skipper model)
5. Reproducibility filtering across replicates
6. HOMER motif analysis for RNA-binding protein sequence preference
7. MultiQC quality report

> **Why UMI deduplication matters for eCLIP:** eCLIP reads naturally cluster at
> UV crosslink sites, so standard duplicate removal (flag-based) would discard real
> signal. UMI-based dedup (10 nt random UMIs in the RT primer) distinguishes true
> PCR duplicates from independent crosslinking events at the same position.

---

## Prerequisites

- UArizona HPC account with access to `/xdisk/haining/maarowosegbe`
- Internet access from a login or transfer node (for downloads)
- ~300 GB disk space (genomes, indices, FASTQ files, results)

---

## Quick start

### Step 1 — One-time environment setup
```bash
cd /xdisk/haining/maarowosegbe/eclip_pipeline
bash setup/00_install_micromamba.sh
```
This installs micromamba to `~/software/micromamba`, creates the `snakemake9` conda
environment, and clones the Skipper pipeline to `~/software/skipper`.

### Step 2 — Download reference genomes and build STAR indices
Submit as a SLURM job (needs 64 GB RAM for STAR indexing):
```bash
sbatch setup/01_download_references.sh
```
Takes ~3–4 hours. Downloads hg38 + mm10 genomes, GFF3 annotations, filters the GFF,
builds STAR indices, and downloads the ENCODE eCLIP blacklist.

### Step 3 — Get the raw FASTQ files
Find the GEO accession in the paper's Data Availability section, fill in the SRR IDs
in `scripts/download_geo_fastqs.sh`, then run:
```bash
bash scripts/download_geo_fastqs.sh
```

### Step 4 — Update FASTQ paths in the manifest
Edit `config/manifest_human.csv` and `config/manifest_mouse.csv` to replace the
placeholder paths with actual paths to your downloaded FASTQ files.

The downloaded files will be named `<LABEL>.fastq.gz` in
`/xdisk/haining/maarowosegbe/eclip_pipeline/data/fastq/` — just confirm the paths match.

### Step 5 — Run the analysis
```bash
# Human experiments (TDP-43, G3BP1, Caprin1 in neurons + brain tissue)
bash run_analysis.sh human

# Mouse age-series (TDP-43 in cerebellum)
bash run_analysis.sh mouse

# Dry run first to confirm job graph:
bash run_analysis.sh human --dry-run
```
Skipper submits each step as a SLURM job automatically. A full run takes ~6–24 hours
depending on queue wait times. Monitor with `squeue -u $USER`.

---

## File structure

```
eclip_pipeline/
├── README.md                        ← this file
├── run_analysis.sh                  ← main entry point
├── setup/
│   ├── 00_install_micromamba.sh     ← step 1: env + Skipper install
│   └── 01_download_references.sh   ← step 2: genome + STAR index
├── config/
│   ├── skipper_config_human.yaml    ← Skipper config, human (hg38)
│   ├── skipper_config_mouse.yaml    ← Skipper config, mouse (mm10)
│   ├── manifest_human.csv           ← sample manifest with FASTQ paths
│   ├── manifest_mouse.csv           ← mouse sample manifest
│   ├── adapters/
│   │   ├── eclip_ip_adapter.fasta   ← 3' adapter for IP samples
│   │   └── eclip_input_adapter.fasta← 3' adapter for SMInput samples
│   └── profiles/
│       └── slurm/
│           └── config.yaml          ← SLURM resource defaults
├── envs/
│   └── snakemake_env.yaml           ← conda env: snakemake + SLURM plugin
└── scripts/
    └── download_geo_fastqs.sh       ← SRA download helper
```

**Results** go to:
- Human: `/xdisk/haining/maarowosegbe/eclip_pipeline/results/human/`
- Mouse:  `/xdisk/haining/maarowosegbe/eclip_pipeline/results/mouse/`

Key output files (produced by Skipper):
```
results/human/
├── reproducible_enriched_windows/    ← final binding sites (TSV)
├── QC/multiqc/*/multiqc_report.html  ← QC report
├── homer/finemapped_results/*/       ← motif analysis
├── secondary_results/bams/           ← deduplicated BAM files
└── secondary_results/bigwigs/        ← coverage tracks (UCSC genome browser)
```

---

## Reusing this pipeline for new samples

To run new eCLIP experiments (same or different RBPs):

1. **Same genome (hg38):** Add rows to `config/manifest_human.csv` — one row per replicate.
   Each unique `Experiment` + `Sample` combination becomes its own analysis.

2. **Different genome:** Copy `config/skipper_config_human.yaml`, update `GFF`,
   `GENOME_FASTA`, `STAR_DIR`, `PARTITION`, `FEATURE_ANNOTATIONS` to the new genome's paths,
   and add a new branch to `run_analysis.sh`.

3. **Different UMI length:** Change `UMI_SIZE` in the config. Most ENCODE4 eCLIP uses 10 nt;
   ENCODE3 uses 5 nt.

---

## Software versions

| Software | Version | Role |
|---|---|---|
| Snakemake | 9.12.0 | Workflow manager |
| snakemake-executor-plugin-slurm | 1.4.0 | SLURM integration |
| Skipper | see `~/software/skipper_commit.txt` | eCLIP peak calling |
| STAR | 2.7.11b (Skipper-managed) | Alignment |
| UMI-tools | 1.1.5 (Skipper-managed) | PCR deduplication |
| cutadapt | 4.9 (Skipper-managed) | Adapter trimming |
| HOMER | 4.11 (Skipper-managed) | Motif analysis |
| GENCODE human | v47 | hg38 annotation |
| GENCODE mouse | M35 | mm10 annotation |

---

## Key paper parameters reproduced here

| Parameter | Value | Source |
|---|---|---|
| eCLIP protocol | ENCODE4 | Paper methods (adapter sequences) |
| UMI length | 10 nt | RT primer: 10 × N random bases |
| Genome (human) | GRCh38 (hg38), GENCODE v47 | Paper methods |
| Genome (mouse) | GRCm38 (mm10) | Paper methods |
| Peak caller | Skipper (Boyle et al. 2023) | Paper methods ref [73] |
| Background | Size-matched input (SMInput) | Standard ENCODE4 eCLIP |
| Replicates | 2 per condition (cells); 2–3 per cohort (tissue) | Paper |
| Reproducibility filter | Window enriched in ≥2 replicates | Skipper default |
| Motif analysis | HOMER | Paper Fig. 2e, 2l |

---

## Citation

If you use this pipeline, please cite:
- Rhine et al. (2025). Neuronal aging causes mislocalization of splicing proteins and
  unchecked cellular stress. *Nature Neuroscience*, 28, 1174–1184.
- Boyle et al. (2023). Skipper analysis of eCLIP datasets enables sensitive detection
  of constrained translation factor binding sites. *Cell Genomics*, 3, 100317.
