# HLA Typing of 1000 Genomes Project Trios with SpecHLA

## Summary

This project performs full-resolution HLA typing on ~600 trios (~1,800 individuals) from the high-coverage 1000 Genomes Project (1kGP) cohort, using **SpecHLA** as the typing tool. The cohort is the expanded 3,202-sample resource from Byrska-Bishop et al. (Cell, 2022), sequenced at 30x on Illumina NovaSeq 6000 and aligned to GRCh38.

The pipeline derives the trio structure directly from the canonical 1kGP pedigree file, resolves each member's CRAM file location on the EBI 1kGP HTTPS server, remotely slices the MHC region from each CRAM, extracts HLA-related reads using SpecHLA's bundled extractor, and runs SpecHLA per sample to produce HLA typing calls. The entire pipeline runs inside a Docker container with pinned tool versions for reproducibility. The final deliverable is a master table of HLA alleles for every individual, joined to family/role/population metadata.

**Data flow at a glance:**

```
pedigree file  +  index files (2)
       │
       ▼
   Stage 1: Build families.json (one entry per trio, with CRAM/CRAI URLs)
       │
       ▼
   Stage 2: Remote slice MHC region from CRAMs → mini-CRAMs (~100-300 MB each)
       │
       ▼
   Stage 3: SpecHLA's extractor → paired FASTQs per sample
       │
       ▼
   Stage 4: Organize FASTQs into family-structured directory tree
       │
       ▼
   Stage 5: Run SpecHLA per sample → typing results
       │
       ▼
   Stage 6: Aggregate into master HLA table
```

**Scope decisions made:**

- **Pedigree-only.** Trio structure is derived directly from the canonical 1kGP pedigree file. No external scope CSV.
- **Per-sample typing only.** SpecHLA's pedigree-phasing pass-2 is **not** in scope, but data is organized so it could be added later.
- **Two-step extraction.** Slice region first, then run SpecHLA's extractor on the local mini-CRAM. Safer and resumable.
- **Region.** Extended MHC plus chr6 alt contigs and HLA contigs. Exact list of HLA-related reads comes from SpecHLA's bundled extractor — we do not hand-roll a BED file.
- **Containerized execution.** Pipeline runs inside Docker with pinned versions. Code in the image, data through mounted volumes.
- **Duos and quads:** Skip duos (incomplete trios). Quads (two children sharing parents) emit as two separate trio entries sharing parental data.

---

## 1. Background

### 1.1 The 1kGP high-coverage resource

- 3,202 samples total: 2,504 original unrelated samples + 698 additional related samples completing 602 trios.
- Sequenced to 30x mean coverage with Illumina NovaSeq 6000, paired 2x150bp.
- Aligned to GRCh38 (full analysis set with decoy and HLA contigs).
- LCL-derived DNA. ~5% of singletons are estimated to be cell-line artifacts (relevant as a caveat, not a blocker for HLA typing).
- All data is publicly available with no access restrictions.

### 1.2 Key data sources

| Resource | URL |
|---|---|
| 1kGP high-coverage data collection root | https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/ |
| Pedigree file (3,202 samples) | https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt |
| Index file: 2,504 unrelated samples | https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index |
| Index file: 698 related samples | https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_698_related_high_coverage.sequence.index |
| Alternative mirror (AWS S3) | `s3://1000genomes/1000G_2504_high_coverage/` |
| Alternative mirror (NCBI) | https://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/ |

### 1.3 SpecHLA

- GitHub: https://github.com/deepomicslab/SpecHLA
- Paper: Wang et al., Cell Reports Methods (2023)
- Genes typed: HLA-A, -B, -C, -DPA1, -DPB1, -DQA1, -DQB1, -DRB1
- Input format used here: paired-end FASTQ (R1/R2)
- Supports Linux and Windows WSL only (containerized Linux is the path here).

---

## 2. Environment

### 2.1 Host

- **Machine:** Remote Windows server, 80 GB RAM, GPU available but unused (none of these tools use GPU).
- **OS layer:** WSL2 on Windows, with Docker Desktop using the WSL2 backend.
- **Working directory on host:** Inside the WSL2 native filesystem (e.g., `~/hla_project/`). **Do not** put working data under `/mnt/c/...` — cross-filesystem I/O is severely slower.
- **Disk budget:** ~1 TB working space recommended (~300 GB mini-CRAMs, ~150-300 GB FASTQs, plus SpecHLA outputs and overhead).

### 2.2 Container

The pipeline runs inside a Docker container built from a Dockerfile pinned to specific versions of every tool. This is the unit of reproducibility.

**Image contents:**
- Pinned base image (e.g., `ubuntu:22.04`, not `latest`).
- Miniconda with the SpecHLA environment installed via `conda -c bioconda -c conda-forge`.
- Pinned versions of `samtools`, `bcftools`, Python, and any other tooling.
- The pipeline's Python scripts under `/app/`.
- Non-root user for sensible file ownership on mounted volumes.

**Not in the image:**
- The reference genome (~3 GB, changes infrequently — mount it).
- Any input data (mounted).
- Any output (written to mounted volume).

**Volume layout:**

```
host (WSL2 native fs)               container
~/hla_project/reference/      →    /data/reference/    (read-only)
~/hla_project/inputs/         →    /data/inputs/       (read-only)
~/hla_project/work/           →    /data/work/         (read-write)
```

Read-only mounts on `reference/` and `inputs/` catch bugs where the pipeline accidentally writes where it shouldn't.

**Image versioning:** Tag images explicitly (e.g., `hla-pipeline:v1.0`). When something downstream changes (newer IMGT/HLA database, updated SpecHLA, etc.), bump to `v1.1` and the diff is fully attributable to known changes.

### 2.3 Reference genome

- File: `GRCh38_full_analysis_set_plus_decoy_hla.fa` (the 1kGP-published reference).
- Indexed with `samtools faidx`.
- Mounted into the container at `/data/reference/`.
- CRAMs cannot be decoded without the matching reference. Wrong reference = silent corruption.

---

## 3. The pipeline

### Stage 0 — Setup (one-time)

1. Install Docker Desktop on the Windows host with WSL2 backend enabled.
2. Build the pipeline image: `docker build -t hla-pipeline:v1.0 .`
3. Download `GRCh38_full_analysis_set_plus_decoy_hla.fa` to `~/hla_project/reference/` and run `samtools faidx` on it (either inside the container or with a host-installed samtools, doesn't matter as long as the `.fai` is alongside the `.fa`).
4. Download the pedigree file and both `.sequence.index` files to `~/hla_project/inputs/`.
5. Run a smoke test: invoke the container against SpecHLA's bundled `example/` data. Do not proceed if this fails.

### Stage 1 — Build `families_1000G_trios.json`

**Inputs (downloaded to `families_fastq_collections/trios/`):**
- `20130606_g1k_3202_samples_ped_population.txt` — canonical pedigree (source of truth for sex, parentage, population).
- `1000G_2504_high_coverage.sequence.index` — CRAM URLs for unrelated samples.
- `1000G_698_related_high_coverage.sequence.index` — CRAM URLs for related samples.

Files are downloaded once on first run and cached locally; subsequent runs reuse the cached copies.

**Logic:**
1. Read the pedigree file. For each row, classify as founder (both parent IDs `0`), duo (exactly one parent ID `0` or one parent missing from cohort), or trio child (both parent IDs present in cohort). Founders are skipped silently; duos are counted and skipped; trio children are emitted as JSON entries.
2. Read both index files (skipping `##` metadata, capturing the `#`-prefixed header). Build a single `cram_index` keyed by `SAMPLE_NAME` (the HG/NA ID, *not* `SAMPLE_ID` which is an ENA accession). Halt on any duplicate or conflicting entry across files.
3. For each trio, look up child/father/mother in `cram_index` and enrich the member records with `cram_url`, `crai_url`, `cram_md5`, raw sub-population, and SpecHLA-friendly population.
4. Write `families_1000G_trios.json` to the same directory.
5. Log a summary: trios written, duos skipped, unique individuals.

**Trio key:** Composite `{FamilyID}_{FatherID}_{MotherID}_{ChildID}` to keep siblings distinct (one entry per child; siblings naturally share parent records by `sample_id`). The code raises immediately if a key collides — should be impossible given the over-determined key, but cheap insurance.

**Validations (halt on failure):**
- Every pedigree row has the expected number of fields.
- Index file headers are well-formed (`#`-prefixed line found before any data line).
- Every CRAM URL begins with `ftp://` (raises if EBI ever changes scheme).
- Sample ID appears as a substring of its own CRAM URL (catches column-shift bugs).
- No duplicate trio keys.
- No conflicting `cram_index` entries across the two index files.
- Every trio member resolves in `cram_index`.

**Sub-population → SpecHLA mapping** (per-member, derived from geographic regions):

| Region | 1kGP codes | SpecHLA `-p` |
|---|---|---|
| East Asian | CHB, JPT, CHS, CDX, KHV, CHD | Asian |
| South Asian | GIH, PJL, BEB, STU, ITU | Asian |
| European | CEU, TSI, GBR, FIN, IBS | Caucasian |
| African | YRI, LWK, GWD, MSL, ESN | Black |
| African-diaspora (Americas) | ASW, ACB | Black |
| Admixed American | MXL, PUR, CLM, PEL | Unknown |

This is finer-grained than mapping at the super-population level: ASW and ACB (which 1kGP groups under AMR) are mapped to `Black` rather than `Unknown`, since their HLA frequencies are dominated by African ancestry.

**Sex translation:** Pedigree uses `1`=M, `2`=F. Translated to `M`/`F` strings in JSON.

**HTTP vs FTP:** Index files give `ftp://ftp.sra.ebi.ac.uk/...`; we rewrite to `https://...` for stage 2 (htslib remote slicing prefers HTTPS). Path is otherwise identical.

**CRAI URL construction:** `cram_url + ".crai"`. Existence on EBI was verified once across all 3202 samples (HEAD requests); the check is left commented in the code for re-enablement on new cohorts.

**JSON schema:**

```json
{
  "metadata": {
    "generated_at": "ISO 8601 UTC timestamp",
    "n_trios": 1116,
    "n_unique_samples": 2418,
    "unique_samples": ["HG00403", "HG00404", "HG00405", "..."],
    "n_duos_skipped": 12,
    "source_files": {
      "pedigree":   "20130606_g1k_3202_samples_ped_population.txt",
      "index_2504": "1000G_2504_high_coverage.sequence.index",
      "index_698":  "1000G_698_related_high_coverage.sequence.index"
    }
  },
  "family_trios": {
    "SH089_HG00403_HG00404_HG00405": {
      "family_id": "SH089",
      "population": "CHS",
      "superpopulation": "EAS",
      "members": {
        "child":  { "sample_id": "HG00405", "sex": "F",
                    "cram_url": "https://...", "crai_url": "https://...",
                    "cram_md5": "...",
                    "member_sample_population": "CHS",
                    "member_spechla_pop": "Asian" },
        "father": { "sample_id": "HG00403", "sex": "M",
                    "cram_url": "https://...", "crai_url": "https://...",
                    "cram_md5": "...",
                    "member_sample_population": "CHS",
                    "member_spechla_pop": "Asian" },
        "mother": { "sample_id": "HG00404", "sex": "F",
                    "cram_url": "https://...", "crai_url": "https://...",
                    "cram_md5": "...",
                    "member_sample_population": "CHS",
                    "member_spechla_pop": "Asian" }
      }
    }
  }
}
```

**Logging:** All progress and summary information is mirrored to `families_fastq_collections/trios/collect_trios_py.log` (overwritten per run) and to stderr. Pipeline failures log the full traceback before re-raising.

### Stage 2 — Remote CRAM slicing

For each individual, pull only the MHC region from the remote CRAM.

**Region to slice:**
- chr6 primary MHC: `chr6:25,000,000-35,000,000` (generous extended MHC window).
- All `chr6_*_alt` contigs (full length).
- All `HLA-*` contigs (full length).

Build the region BED file once by parsing the `.fai` of the reference. Reuse for every sample.

**Per-sample command (conceptual, runs inside container):**
```bash
samtools view -C \
  -T /data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  -L /data/work/mhc_regions.bed \
  -o /data/work/mini_crams/<sample>.mhc.cram \
  <cram_url>
samtools index /data/work/mini_crams/<sample>.mhc.cram
```

**Operational requirements:**
- **Resumable:** write a `.done` marker per sample only after slice + index + sanity check pass. Driver script skips any sample with `.done`.
- **Sanity check:** `samtools quickcheck` plus a minimum read count threshold.
- **Concurrency:** 4-8 parallel jobs against EBI. More risks throttling.
- **Logging:** append-only log per sample with status, duration, file size, read count.
- **Run in `tmux`** on the host (or as a background container), expecting 1-3 days wall time for 1,800 samples.

**Failure modes to expect and handle:**
- Transient HTTPS errors (retry with backoff).
- Reference mismatch (manifests as decoder errors — caught by quickcheck).
- Truncated downloads (caught by quickcheck + read count threshold).
- Long-running connections to EBI from inside container (default Docker networking handles this, but retry/timeout logic must be robust).

### Stage 3 — HLA read extraction

For each mini-CRAM, run SpecHLA's bundled extractor. **Do not** hand-roll this with `samtools fastq` — let SpecHLA's extractor decide what counts as an HLA-related read.

```bash
spechla-extract-hla-reads \
  -s <sample_id> \
  -b /data/work/mini_crams/<sample>.mhc.cram \
  -r hg38 \
  -o /data/work/extracted_reads/<sample>/
```

Output: paired FASTQs ready for SpecHLA proper.

Same operational pattern as Stage 2: per-sample `.done` markers, resumable, parallel. CPU-bound, so concurrency = number of cores you want to dedicate (probably 8-16).

### Stage 4 — Family-organized directory tree

After extraction, organize FASTQs into a tree keyed by family. Symlinks save disk.

```
data/
  <superpop>/<family_id>/
    child/<child_id>_1.fq.gz
    child/<child_id>_2.fq.gz
    father/<father_id>_1.fq.gz
    father/<father_id>_2.fq.gz
    mother/<mother_id>_1.fq.gz
    mother/<mother_id>_2.fq.gz
    family.tsv                   # subset of families.json for this family
```

This makes per-sample SpecHLA calls and any later trio analysis trivially scriptable.

### Stage 5 — Per-sample SpecHLA typing

For each individual:

```bash
spechla \
  -n <sample_id> \
  -1 <r1.fq.gz> \
  -2 <r2.fq.gz> \
  -o /data/work/spechla_output/<sample_id>/ \
  -p <spechla_pop> \
  -j <threads>
```

Use the `spechla_pop` field from `families.json` for `-p`. Defaulting to `Unknown` works but population-aware annotation gives better calls.

Same resumability pattern. SpecHLA per sample is in the tens of minutes ballpark on a single core; with 8-16 parallel samples, ~1-2 days wall time for 1,800 samples.

**Outputs per sample** (in `spechla_output/<sample_id>/`):
- `hla.result.txt` — main typing result, one row per sample with all 8 genes × 2 alleles.
- `hla.result.details.txt` — best-matching alleles per gene.
- `hla.result.g.group.txt` — G-group resolution.
- `hla.allele.*.HLA_*.fasta` — reconstructed allele sequences.
- `HLA_*_freq.txt` — haplotype frequencies.
- `HLA_*.rephase.vcf.gz` — phased VCFs per gene.

### Stage 6 — Aggregate results

Walk the output tree, concatenate every `hla.result.txt` into one master table joined to family/role/population metadata from `families.json`.

**QC pass on the aggregated table:** for each trio, check Mendelian consistency — every allele in the child should appear in at least one parent for each gene. Mendelian errors are a free QC signal on the typing itself. Build this into Stage 6 even though we're not doing pedigree-phasing.

---

## 4. Running the pipeline

### Build the image

```bash
docker build -t hla-pipeline:v1.0 .
```

### Stage invocations

Each stage is one `docker run` against the same image with a stage-specific entrypoint or argument:

```bash
docker run --rm \
  -v ~/hla_project/reference:/data/reference:ro \
  -v ~/hla_project/inputs:/data/inputs:ro \
  -v ~/hla_project/work:/data/work \
  hla-pipeline:v1.0 \
  python /app/stage1_build_json.py
```

(Equivalent invocations for Stages 2-6, varying only the script name.)

For long-running stages (2, 3, 5), run with `-d` and `--name <stage>` so they can be detached, monitored, and inspected.

---

## 5. Critical risks and mitigations

| Risk | Mitigation |
|---|---|
| Reference mismatch when decoding CRAMs | Use the exact `GRCh38_full_analysis_set_plus_decoy_hla.fa` from the 1kGP. Run one full sample end-to-end before scaling up. |
| EBI throttling on parallel downloads | Cap at 4-8 concurrent connections. Retry with backoff. |
| WSL2 cross-filesystem I/O bottleneck | All host data lives in WSL2 native fs (`~/...`), never under `/mnt/c/`. |
| Missing samples discovered mid-pipeline | Validate everything in Stage 1. Skip log is the source of truth for what's actually being processed. |
| Interrupted runs losing progress | Per-sample `.done` markers in Stages 2, 3, 5. Driver scripts skip completed samples. |
| Hand-rolled HLA region BED missing reads on alts | Use SpecHLA's bundled extractor in Stage 3. Generous chr6 window + all alt + HLA contigs in Stage 2. |
| Tool version drift between runs | Pinned versions in Dockerfile. Image tagged explicitly. |
| Quads collapsed into single trios | Stage 1 detects shared `(FatherID, MotherID)` pairs and emits both as separate entries. |

---

## 6. Considerations not addressed

- **Novoalign license:** Optional. Falls back to bowtie2 if absent. Speed gain probably not worth the license cost.
- **Cloud migration:** If EBI throughput becomes painful, the same data is on AWS S3 (`s3://1000genomes/...`). Free egress within `us-east-1`. Only worth it if running on AWS already.
- **Trio pedigree phasing (SpecHLA pass-2):** Out of scope for this run. The directory layout and JSON support adding it later without restructuring.
- **HLA LOH detection:** SpecHLA supports it (`spechla-loh`). Not relevant here — these are LCL samples, not tumor.
- **IMGT/HLA database version:** SpecHLA ships with a default. If a newer version is desired, run `bin/renew_HLA_annotation_db.pl` once during image build (will change the image hash — bump version tag).

---

## 7. Quick reference

**Project layout (host side):**

```
hla_project/
  Dockerfile
  app/
    stage1_build_json.py
    stage2_slice_crams.py
    stage3_extract_reads.py
    stage4_organize.py
    stage5_run_spechla.py
    stage6_aggregate.py
  reference/
    GRCh38_full_analysis_set_plus_decoy_hla.fa(.fai)
  inputs/
    20130606_g1k_3202_samples_ped_population.txt
    1000G_2504_high_coverage.sequence.index
    1000G_698_related_high_coverage.sequence.index
  work/
    metadata/
      families.json
      skipped_families.log
    mhc_regions.bed
    mini_crams/
      <sample>.mhc.cram(.crai)
      <sample>.done
    extracted_reads/
      <sample>/<sample>_1.fq.gz
      <sample>/<sample>_2.fq.gz
      <sample>/.done
    data/
      <superpop>/<family_id>/...
    spechla_output/
      <sample_id>/...
    results/
      master_hla_table.tsv
      mendelian_qc.tsv
    logs/
      stage2.log
      stage3.log
      stage5.log
```

**Sample counts to expect:**
- 602 candidate trios in the pedigree → ~600 trios written (minus a small number skipped, plus any quad expansions).
- ~1,800 individuals processed end-to-end.
- One row per individual in the final master table.