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
> Raw FASTQ reads
│
├─ Falco ──────────────── Quality control
│                         (base quality · GC content · bisulfite signature)
│
├─ bwameth ────────────── Bisulfite-aware alignment → hg38
│                         (handles OT / OB / CTOT / CTOB strand complexity)
│
├─ MethylDackel mbias ─── Positional bias check
│                         (SVG plots per strand; flag 5′/3′ artifacts)
│
├─ MethylDackel extract ─ CpG methylation fractions → bedGraph
│                         (mergeContext: per-dinucleotide values)
│
├─ BedGraph → bigWig ──── Format conversion for deepTools
│
├─ computeMatrix ──────── Methylation matrix centred on CGIs
│
├─ plotProfile ─────────── Average methylation profile across CGI window
│
└─ Metilene ────────────── DMR detection (Normal NB1/NB2 vs. Tumour BT198)
(binary segmentation + Mann-Whitney U)

**Why can't we use a standard aligner?** Bisulfite treatment converts unmethylated cytosines to uracil (→ thymine after PCR). The genome is no longer symmetric — standard aligners like HISAT2 or STAR assume Watson-Crick pairing and fail completely. `bwameth` converts both reads and reference in silico, correctly handling all four strand orientations.

**Step 1 — QC with Falco**
Falco (a faster FastQC reimplementation) checks read quality, GC content, adapter contamination, and base composition. In bisulfite libraries, the dramatic shift in T/C ratio — C plummets, T spikes — isn't a problem. It's the signal that bisulfite conversion worked.

**Step 2 — Alignment with bwameth**
Paired-end reads aligned to GRCh38 (hg38full). Output: sorted, indexed BAM.

**Step 3 — Bias assessment with MethylDackel mbias**
Generates SVG plots of methylation level across read positions per strand. End-repair artifacts show up as upticks at the 5′ end of OT strand reads. Singleton and discordant alignments retained (`--keepSingleton --keepDiscordant`) to maximize coverage on this small subset.

**Step 4 — Methylation extraction**
Per-CpG methylation fractions extracted to bedGraph. The `--mergeContext` flag pools both strands of each CpG into a single value — the standard for CpG-level analysis.

**Step 5 — Visualization around CpG islands**
bedGraph → bigWig → `computeMatrix` (reference-point mode, centred on CGI midpoints) → `plotProfile`. Normal samples show the expected methylation dip at CGI centres. In cancer, some of those dips flatten — epigenetic silencing in action.

> Chromosome naming note: pre-computed files use Ensembl-style names (`1`, `2`) rather than UCSC-style (`chr1`, `chr2`). The `Replace column` tool handles remapping before bigWig conversion.

**Step 6 — DMR detection with Metilene**
Metilene scans for stretches of consecutively different CpGs between groups using binary segmentation. Restricted to annotated CGI regions to focus on functionally relevant changes. Outputs: DMR BED file with methylation deltas + q-values, plus a PDF summary.

Groups:
- Normal → NB1_CpG.meth.bedGraph, NB2_CpG.meth.bedGraph
- Tumour → BT198_CpG.meth.bedGraph

### What We Found

- **Conversion confirmed:** T/C skew in QC validates successful bisulfite treatment.
- **5′ bias present:** OT strand reads show elevated methylation at the first ~5 bases — a common end-repair artifact. Corrected by trimming during extraction.
- **Normal tissue looks normal:** deepTools profiles show clean CGI hypomethylation dips in NB1 and NB2.
- **Cancer disrupts the pattern:** Tumour samples show altered CGI profiles with some islands gaining methylation — consistent with epigenetic silencing of tumour suppressors.
- **DMRs detected:** Metilene identified significant DMRs at CGIs distinguishing NB1/NB2 from BT198 — candidate loci for methylation-driven gene silencing.

---

## Part 2 — Benchmarking Epigenetic Clocks

### The Clock Concept

In 2013, Steve Horvath showed that methylation at a small set of CpGs could predict chronological age across virtually any tissue. That was surprising — and useful. Since then, the field has moved fast, producing clocks optimised for different outcomes: biological age, disease risk, pace of aging, tissue-specific prediction.

The problem is that most papers validate a new clock on the same data it was trained on (or adjacent data). This benchmark runs eight clocks on two independent cohorts and compares the results — testing what actually generalizes.

**Bio-Learn** provides a unified Python interface to published aging clocks and GEO datasets. It handles CpG imputation, normalization, and prediction through a consistent API, so swapping between clocks is just changing a string.

### Cohorts Used

**GSE40279 — Hannum et al. (2013)**

| Property | Value |
|---|---|
| Tissue | Whole blood |
| Platform | Illumina 450K |
| Full size | 656 samples × 473,034 CpGs |
| Used | First 200 samples |
| Notes | The Hannum clock was trained on this dataset — it's the "home ground" benchmark |

**GSE51057**

| Property | Value |
|---|---|
| Tissue | Whole blood |
| Platform | Illumina 450K / EPIC |
| Full size | 329 samples × 485,577 CpGs |
| Used | First 200 samples |
| Notes | Independent cohort — no clock here was trained on this data, so it's the real generalization test |

### The Eight Clocks

**Horvath v1 (2013)** — The original pan-tissue clock. Elastic net regression on 353 CpGs, trained across 51 tissue types. Still the canonical baseline. Works on almost any tissue, which was genuinely revolutionary.

**Horvath v2 / SkinBloodClock (2018)** — Skin and blood optimized update using 391 CpGs. Corrects the systematic underestimation of age in older individuals that v1 shows.

**Hannum (2013)** — Blood-specific lasso clock trained on 71 CpGs, built from this exact dataset (GSE40279). Expected to perform best on Dataset 1 — useful "in-sample" reference point.

**PhenoAge (2018)** — Two-stage design: first builds a "phenotypic age" score from nine clinical biomarkers via a mortality model, then trains methylation to predict that score. Captures disease risk and biological aging rather than calendar time.

**DunedinPACE (2022)** — Measures *speed* of aging, not absolute age. Output ~1.0 = average pace. Built from the longitudinal Dunedin Study where the same people were tracked for decades. Better for intervention studies than point-in-time clocks.

**YingCausAge (2024)** — Uses causal inference to identify CpGs whose changes are *caused by* aging rather than merely correlated with it. Designed to separate genuine aging signal from lifestyle confounders like smoking or BMI.

**Knight (2016)** — Trained on cord blood to estimate gestational age. Intentionally included as a control: when applied to adult whole blood, it produces compressed, biologically uninformative predictions. Illustrates why tissue specificity matters.

**Bocklandt (2011)** — One of the earliest epigenetic clocks, trained on saliva. Outputs fractional values (~0.3–0.4) in blood instead of years — another tissue mismatch control, and a useful historical anchor showing how much the field has advanced.

### How the Analysis Works

**Environment:** Google Colab. Install: `pip install biolearn matplotlib seaborn pandas numpy scipy`

**Data loading:** `DataLibrary().get(accession).load()` downloads and caches methylation matrices from GEO. Trimmed to 200 samples immediately after loading + explicit `gc.collect()` to stay within Colab's RAM limits (~6–8 GB peak with both matrices loaded).

**Inference loop:** A `run_clocks()` function instantiates each clock via `ModelGallery().get(clock)`, calls `model.predict(data)`, and collects results. Failed clocks (try/except) fill with NaN so downstream analysis continues uninterrupted.

**Compatibility check:** Clocks verified against 450K array CpG coverage before the main analysis. Knight, Bocklandt, YingCausAge, YingDamAge, YingAdaptAge, and DNAmTL all work. Vidalin, MEAT2, and Stubbs were dropped — missing CpG coverage on 450K.

### Results & Visuals

**Correlation matrix** (`corr_GSE40279.png`, `corr_GSE51057.png`)

Pairwise Pearson correlations across all clock predictions and chronological age, shown as annotated seaborn heatmaps. Clocks working from the same methylation-age signal in blood (Horvath v1/v2, Hannum, PhenoAge, YingCausAge) cluster together with high mutual correlations. DunedinPACE sits apart — its pace scale is conceptually different from absolute age. Knight and Bocklandt barely correlate with anything, as expected.

![Correlation matrix GSE40279](https://github.com/user-attachments/assets/9d197957-c313-43c9-a115-8d7cda3c7b47)
![Correlation matrix GSE51057](https://github.com/user-attachments/assets/aa43287a-f16d-4824-80ab-3de4789d45b2)

---

**Age deviation heatmap** (`deviation_GSE40279.png`, `deviation_GSE51057.png`)

Per-sample age deviation = `predicted − chronological`. Positive = biologically older than calendar age. Shown as a clustermap with hierarchical clustering on both axes. Highlights: samples consistently deviating high across clocks, clock clusters with similar bias patterns, and systematic over/underestimation that might reflect cohort-level biology or distributional shift from training data.

![Deviation heatmap](https://github.com/user-attachments/assets/bb941924-a30a-4a5e-bde1-5b085ee59cbc)

---

**Predicted vs. chronological age scatter** (`scatter_GSE40279.png`, `scatter_GSE51057.png`)

Eight scatter plots in a grid — predicted age on y, chronological on x, identity diagonal included. Pearson r and p-value annotated per subplot. A tight diagonal = accurate. A shifted diagonal = calibration bias. A flat scatter = DunedinPACE (expected — it's a pace, not an age). A random scatter = Knight/Bocklandt (also expected — tissue mismatch).

![Scatter plots](https://github.com/user-attachments/assets/f7fb55cd-0dda-4e01-b89e-1fdcc4be7aa2)

### What Changes Across Datasets

Comparing GSE40279 (Hannum's home ground) vs. GSE51057 (independent cohort) is where the real information lives:

- Clocks with consistent r across both datasets are genuinely robust.
- Clocks that degrade on GSE51057 may be overfit to the Hannum dataset's distribution.
- Deviation patterns that replicate across cohorts likely reflect true biological age differences, not dataset artifacts.
- Single-dataset evaluation is insufficient — external validation isn't optional.

---

## Project Layout
.
├── README.md                       ← you are here
├── biological_clocks.ipynb         ← Part 2 full Colab notebook
│
├── outputs/
│   ├── corr_GSE40279.png           ← Correlation matrix, Dataset 1
│   ├── corr_GSE51057.png           ← Correlation matrix, Dataset 2
│   ├── deviation_GSE40279.png      ← Deviation heatmap, Dataset 1
│   ├── deviation_GSE51057.png      ← Deviation heatmap, Dataset 2
│   ├── scatter_GSE40279.png        ← Scatter grid, Dataset 1
│   └── scatter_GSE51057.png        ← Scatter grid, Dataset 2
│
└── galaxy/
└── methylation-seq-tutorial.html   ← GTN tutorial reference (offline)

---

## Setup & Requirements

### Part 1 — Galaxy

All tools run at [usegalaxy.org](https://usegalaxy.org/) — no local install needed.

| Tool | Version | Purpose |
|---|---|---|
| Falco | 1.2.4+galaxy0 | Read QC |
| bwameth | 0.2.7+galaxy0 | Bisulfite-aware alignment to hg38 |
| MethylDackel | 0.5.2+galaxy0 | Bias assessment + methylation extraction |
| Replace column | 0.2 | Chromosome name remapping (Ensembl → UCSC) |
| Wig/BedGraph-to-bigWig | — | Format conversion |
| computeMatrix | 3.5.4+galaxy0 | CGI-centred methylation matrix |
| plotProfile | 3.5.4+galaxy0 | Profile visualization |
| Metilene | 0.2.6.1 | DMR detection |

Reference genome: GRCh38 / hg38 (full build for alignment; standard for extraction)

### Part 2 — Python

```bash
pip install biolearn matplotlib seaborn pandas numpy scipy
```

Open `biological_clocks.ipynb` in [Google Colab](https://colab.research.google.com/) and run cells top to bottom. Bio-Learn downloads and caches GEO datasets automatically on first run.

| Package | Version | Role |
|---|---|---|
| biolearn | 0.9.1 | Datasets + clock models |
| pandas | 2.2.2 | DataFrames + result storage |
| numpy | 2.0.2 | Numerics + NaN handling |
| matplotlib | 3.10.0 | Plot backend |
| seaborn | 0.13.2 | Heatmaps + scatter + correlations |
| scipy | 1.16.3 | Pearson r + statistical annotation |
| scikit-learn | 1.6.1 | biolearn dependency |
| cvxpy / ecos / osqp | various | Optimization deps of biolearn |

Estimated runtime: 10–25 min in Colab depending on download speed. Peak RAM: ~6–8 GB with both methylation matrices loaded.

---

## References

1. Lin, I.-H., Chen, D.-T., Chang, Y.-F., Lee, Y.-L., Su, C.-H. et al. (2015). Hierarchical Clustering of Breast Cancer Methylomes Revealed Differentially Methylated and Expressed Breast Cancer Genes. *PLOS ONE* 10(1): e0118453. https://doi.org/10.1371/journal.pone.0118453

2. Hannum, G., Guinney, J., Zhao, L., Zhang, L., Hughes, G. et al. (2013). Genome-wide Methylation Profiles Reveal Quantitative Views of Human Aging Rates. *Molecular Cell* 49(2): 359–367. https://doi.org/10.1016/j.molcel.2012.10.016

3. Horvath, S. (2013). DNA methylation age of human tissues and cell types. *Genome Biology* 14: R115. https://doi.org/10.1186/gb-2013-14-10

