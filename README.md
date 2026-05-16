# Methylation Matters
### Mapping the Epigenome — from Raw Reads to Biological Clocks

*A bioinformatics deep-dive into DNA methylation, split across two complementary approaches: a full WGBS sequencing pipeline built on Galaxy, and a multi-clock epigenetic aging benchmark powered by the Bio-Learn Python library.*

---

## What This Project Is About

Methylation is quiet. It doesn't change the DNA sequence — it just adds a small chemical tag to cytosine bases, mostly at CpG dinucleotides. But those quiet tags control whether genes get read or silenced, and their patterns shift dramatically in cancer and aging. This project explores both.

**Two problems. Two toolkits. One molecule.**

| | WGBS Pipeline | Epigenetic Clocks |
|---|---|---|
| **Where it runs** | Galaxy (usegalaxy.org) | Google Colab / Python |
| **Input data** | FASTQ sequencing reads | EPIC/450K array β-values |
| **What we're doing** | End-to-end methylation calling + DMR detection | Benchmarking 8 aging clocks across 2 cohorts |
| **Core tools** | bwameth · MethylDackel · Metilene · deepTools | Bio-Learn · pandas · seaborn · scipy |
| **Biology angle** | Breast cancer epigenomics | Epigenetic aging & clock comparison |
| **Reference papers** | Lin et al. (2015) | Hannum et al. (2013) + clock papers |

---

## Contents

- [The Biology, Briefly](#the-biology-briefly)
- [Part 1 — Building a WGBS Pipeline on Galaxy](#part-1--building-a-wgbs-pipeline-on-galaxy)
  - [Why Breast Cancer?](#why-breast-cancer)
  - [Sample Set](#sample-set)
  - [Pipeline Walkthrough](#pipeline-walkthrough)
  - [What We Found](#what-we-found)
- [Part 2 — Benchmarking Epigenetic Clocks](#part-2--benchmarking-epigenetic-clocks)
  - [The Clock Concept](#the-clock-concept)
  - [Cohorts Used](#cohorts-used)
  - [The Eight Clocks](#the-eight-clocks)
  - [How the Analysis Works](#how-the-analysis-works)
  - [Results & Visuals](#results--visuals)
  - [What Changes Across Datasets](#what-changes-across-datasets)
- [Project Layout](#project-layout)
- [Setup & Requirements](#setup--requirements)
- [References](#references)

---

## The Biology, Briefly

DNA methylation happens when a methyl group attaches to cytosine — almost always at a CpG dinucleotide. In healthy cells, CpG islands (dense CpG clusters at ~70% of gene promoters) stay unmethylated, keeping those promoters open and active. Outside islands, background methylation sits high (~70–80%) across the genome.

In **cancer**, this balance breaks. Two things happen simultaneously: global methylation drops (destabilizing the genome), while island methylation rises at the promoters of tumour suppressor genes — silencing them epigenetically. In **aging**, methylation levels at hundreds of specific CpGs shift in predictable, clock-like ways that can predict how old someone is, or how fast they're aging.

---

## Part 1 — Building a WGBS Pipeline on Galaxy

### Why Breast Cancer?

This analysis follows **Lin et al. (2015)**, which used hierarchical clustering of breast cancer methylomes to uncover genes silenced by epigenetic reprogramming. The data is available via [Zenodo record 557099](https://zenodo.org/record/557099), and the tutorial framework comes from the [Galaxy Training Network](https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/methylation-seq/tutorial.html).

Breast cancer is a good model here — it's well-characterized epigenetically, has well-defined tumour suppressor loci, and available matched normal/tumour samples make DMR detection meaningful.

### Sample Set

Five tissue types representing the healthy-to-malignant spectrum:

| Sample | Tissue | What to expect |
|---|---|---|
| NB1, NB2 | Normal breast | Canonical patterns — unmethylated CGIs at active promoters |
| BT089 | Fibroadenoma (benign) | Minor early-stage epigenetic changes |
| BT126, BT198 | Invasive ductal carcinoma | Widespread reprogramming expected |
| MCF7 | Breast adenocarcinoma cell line | In vitro cancer reference |

> A small genomic subset (`subset_1.fastq`, `subset_2.fastq`) is used for compute feasibility. Pre-computed alignments are provided for downstream steps.

### Pipeline Walkthrough
