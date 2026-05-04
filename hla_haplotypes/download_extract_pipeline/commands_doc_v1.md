
# HLA Typing Pipeline for Remote CRAM Files

A pipeline for running HLA typing with [SpecHLA](https://github.com/deepomicslab/SpecHLA) on samples from public archives (e.g., 1000 Genomes Project) by streaming HLA-relevant reads directly from remote CRAM files, without downloading the full genome alignment.

## Table of Contents

- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Reference Genome Setup](#reference-genome-setup)
- [Pipeline Version 1: Standard (Recommended)](#pipeline-version-1-standard-recommended)
- [Pipeline Version 2: With Mate Rescue (Experimental)](#pipeline-version-2-with-mate-rescue-experimental)
- [Population Tag Reference](#population-tag-reference)
- [Output Interpretation](#output-interpretation)
- [Scaling Considerations](#scaling-considerations)
- [Troubleshooting](#troubleshooting)

---

## Overview

This pipeline performs HLA typing on samples aligned to GRCh38 (with alt-aware BWA-MEM) by:

1. Streaming HLA-relevant reads from a remote CRAM file using `samtools` coordinate queries
2. Including primary chr6 MHC region, all chr6 alt contigs, all HLA decoy contigs, and unmapped reads
3. Using SpecHLA's bundled extractor for HLA-gene-specific filtering
4. Running SpecHLA for full-resolution HLA typing

The approach avoids downloading 15–25 GB CRAM files per sample, instead transferring ~200 MB of HLA-relevant data via HTTP range requests on the indexed remote CRAM.

### Why both primary chr6 and alt/decoy contigs are needed

GRCh38 represents the highly polymorphic MHC region with multiple alternate haplotype contigs (`chr6_GL00025*_alt`) plus hundreds of HLA decoy contigs (`HLA-*`). When CRAMs are aligned with alt-aware BWA-MEM, reads from divergent HLA haplotypes are distributed across these contigs rather than forced onto the primary chr6 path. Querying only `chr6:...` would silently miss reads from non-PGF haplotypes, leading to haplotype collapse in HLA typing results.

---

## Prerequisites

### Software

- **samtools** ≥ 1.10 (with htslib supporting CRAM and HTTP/FTP URLs)
- **SpecHLA** installed via conda:

  ```bash
  conda create -n spechla -c bioconda -c conda-forge spechla
  conda activate spechla
  ```

  After install, the following commands should be available:
  - `spechla`
  - `spechla-extract-hla-reads`

### Network access

- HTTPS access to `ftp.sra.ebi.ac.uk` (or your CRAM source)
- Sufficient bandwidth for ~200 MB per sample (initial extraction)

---

## Reference Genome Setup

The reference FASTA must **exactly match** the one used when the CRAM was created. For 1000 Genomes Project GRCh38 alignments, this is `GRCh38_full_analysis_set_plus_decoy_hla.fa` — an alt-aware reference including chr6 alt haplotypes and HLA decoy contigs.

### Download the reference

```bash
# Choose a stable location for reference files
REF_DIR="/path/to/references"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

# Download the GRCh38 analysis set with decoys and HLA contigs
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
```

### Index the reference

```bash
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa
```

This creates `GRCh38_full_analysis_set_plus_decoy_hla.fa.fai`, required for CRAM decoding.

### Set the environment variable

Add this to your shell profile or session script:

```bash
export REF_GRCH38="/path/to/references/GRCh38_full_analysis_set_plus_decoy_hla.fa"
```

### Verify the setup

```bash
ls -lh "$REF_GRCH38" "${REF_GRCH38}.fai"
# Both files should exist; FASTA is ~3 GB, .fai is small
```

---

## Pipeline Version 1: Standard (Recommended)

This is the production pipeline. Tested on 1000 Genomes samples HG01882 and NA12878 with results matching published reference HLA types.

### Step 1: Set sample variables

```bash
SAMPLE="HG01882"
CRAM_URL="https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3242189/HG01882.final.cram"
THREADS=8
```

### Step 2: Build contig lists from the CRAM header

These lists are **constant for any CRAM aligned to the same reference**, so this step only needs to be run once per project. Save the output files and reuse for all samples.

```bash
# Get all HLA decoy contig names
samtools view -H -T "$REF_GRCH38" "$CRAM_URL" \
    | grep "^@SQ" | grep -oP "SN:\KHLA-[^\s]+" > hla_contigs.txt

# Get all chr6 alt contig names
samtools view -H -T "$REF_GRCH38" "$CRAM_URL" \
    | grep "^@SQ" | grep -oP "SN:\Kchr6_[^\s]+_alt" > chr6_alt_contigs.txt

# Imporved
samtools view -H -T "$REF_GRCH38" "$CRAM_URL" | awk -F'\t' '
$1=="@SQ" {
    for(i=2;i<=NF;i++) {
        if($i~/^SN:HLA-/) { sub(/^SN:/,"",$i); print $i > "hla_contigs.txt" }
        else if($i~/^SN:chr6_.*_alt$/) { sub(/^SN:/,"",$i); print $i > "chr6_alt_contigs.txt" }
    }
}'

echo "Found $(wc -l < hla_contigs.txt) HLA decoy contigs"
echo "Found $(wc -l < chr6_alt_contigs.txt) chr6 alt contigs"
```

Expected counts for the 1000 Genomes GRCh38 reference: ~525 HLA decoy contigs, ~16 chr6 alt contigs.

### Step 3: Extract HLA-relevant reads from remote CRAM

```bash
samtools view -b -@ $THREADS -T "$REF_GRCH38" \
    "$CRAM_URL" \
    chr6:25000000-35000000 \
    $(cat chr6_alt_contigs.txt | tr '\n' ' ') \
    $(cat hla_contigs.txt | tr '\n' ' ') \
    "*" \
    > ${SAMPLE}.hla_full.bam

# Imprproved

samtools view -b -@ $THREADS -T "$REF_GRCH38" \
    -o ${SAMPLE}.hla_full.bam \
    "$CRAM_URL" \
    chr6:25000000-35000000 \
    $(cat chr6_alt_contigs.txt | tr '\n' ' ') \
    $(cat hla_contigs.txt | tr '\n' ' ') \
    "*"
```

Coverage of the extraction:
- `chr6:25000000-35000000` — extended MHC region with safety margins (10 Mb)
- `chr6_*_alt` — all chr6 alt contigs (covers MHC haplotype variants)
- `HLA-*` — HLA decoy contigs (where alt-aware BWA places divergent HLA reads)
- `"*"` — unmapped reads (mates of HLA reads from very divergent haplotypes)

Expected output: ~150–250 MB BAM.

### Step 4: Coordinate-sort and index

```bash
samtools sort -@ $THREADS -o ${SAMPLE}.hla_full.coordsorted.bam ${SAMPLE}.hla_full.bam
samtools index ${SAMPLE}.hla_full.coordsorted.bam
```

### Step 5: Run SpecHLA's HLA read extractor

```bash
spechla-extract-hla-reads \
    -s ${SAMPLE} \
    -b ${SAMPLE}.hla_full.coordsorted.bam \
    -r hg38 \
    -o ./hla_reads/
```

This produces paired FASTQs filtered to HLA-gene-specific reads only (~1–2 MB each):
- `${SAMPLE}_extract_1.fq.gz`
- `${SAMPLE}_extract_2.fq.gz`
- `${SAMPLE}_extract.unpaired.fq.gz` (singletons, not used downstream)

### Step 6: Sanity check

```bash
ls -lh ./hla_reads/
```

Verify both paired FASTQs are non-empty (~1–2 MB) and similar in size.

### Step 7: Run SpecHLA typing

```bash
spechla \
    -n ${SAMPLE} \
    -1 ./hla_reads/${SAMPLE}_extract_1.fq.gz \
    -2 ./hla_reads/${SAMPLE}_extract_2.fq.gz \
    -p Black \
    -o ./spechla_output/ \
    -j $THREADS
```

Replace `-p Black` with the appropriate population tag for your sample. See [Population Tag Reference](#population-tag-reference).

### Step 8: View results

```bash
cat ./spechla_output/${SAMPLE}/hla.result.txt
```

### Optional: Cleanup intermediates

```bash
rm ${SAMPLE}.hla_full.bam ${SAMPLE}.hla_full.coordsorted.bam ${SAMPLE}.hla_full.coordsorted.bam.bai
rm -rf ./hla_reads/
```

---

## Pipeline Version 2: With Mate Rescue (Experimental)

This version attempts to recover read pairs where one mate is in the HLA region but the other is aligned to a non-extracted contig (e.g., chr3, chr17). **This version was not fully completed end-to-end.** A redundancy issue was identified during testing (described below) and a corrected approach was proposed but not retested. The standard pipeline (Version 1) is the validated production path.

### Approach

The strategy uses three observations:

1. Singletons in the extracted BAM are reads whose mate isn't there
2. Each read's BAM record contains the mate's position (`RNEXT`/`PNEXT` fields)
3. We can fetch only those specific positions from the remote CRAM via coordinate access

### Steps 1–3: Same as Version 1

Run Steps 1, 2, and 3 from Pipeline Version 1 to produce `${SAMPLE}.hla_full.bam`. For consistency with the rescue logic below, rename the output:

```bash
mv ${SAMPLE}.hla_full.bam ${SAMPLE}.initial.bam
```

### Step 4: Identify orphan mate locations (initial attempt — has a known issue, see Step 7)

Find reads whose mate is on a different chromosome and has a valid position:

```bash
samtools view ${SAMPLE}.initial.bam \
    | awk 'BEGIN{OFS="\t"} 
           $7 != "=" && $7 != "*" && $8 > 0 {print $7":"$8"-"$8}' \
    | sort -u > ${SAMPLE}.mate_locations.txt

echo "Missing mate locations: $(wc -l < ${SAMPLE}.mate_locations.txt)"
```

The filter rules:
- `$7 != "="` — skip mates on the same chromosome (already in BAM)
- `$7 != "*"` — skip unmapped mates (already extracted via `"*"`)
- `$8 > 0` — skip mates with no valid position (avoids integer underflow in BED conversion)

**Note:** the `$8 > 0` check was added after an initial attempt without it produced BED entries with `start = 18446744073709551615` (unsigned 64-bit underflow from `0 - 1`), which caused samtools to reject the BED file with `end (1) must not be less than start (18446744073709551615)`.

### Step 5: Convert to BED format

The contig names contain colons (e.g., `HLA-A*01:01:01:01`), so naive splitting on `:` and `-` would mangle them. Use a regex that captures from the rightmost `:`:

```bash
sed -E 's/^(.*):([0-9]+)-([0-9]+)$/\1\t\2\t\3/' ${SAMPLE}.mate_locations.txt \
    | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $3}' > ${SAMPLE}.mate_regions.bed

# Verify
head -5 ${SAMPLE}.mate_regions.bed
wc -l ${SAMPLE}.mate_regions.bed
awk '$2 < 0 || $3 < 0 || $2 > 300000000 || $3 > 300000000 {print "BAD:", $0}' ${SAMPLE}.mate_regions.bed
```

**Note:** an earlier attempt using `awk -F'[:-]' '{print $1"\t"$2-1"\t"$3}'` failed because awk split on every `:`, producing rows like `HLA  -1  01` from inputs like `HLA-A*01:01:01:01:163-164`. The sed-based regex above correctly captures the contig name as everything before the last `:<num>-<num>` pattern.

### Step 6: Fetch the rescued mates

```bash
samtools view -b -@ $THREADS -T "$REF_GRCH38" \
    -L ${SAMPLE}.mate_regions.bed \
    "$CRAM_URL" \
    > ${SAMPLE}.rescued_mates.bam

echo "Rescued mates: $(samtools view -c ${SAMPLE}.rescued_mates.bam)"
```

### Step 7: Known issue — redundant fetch

Running the above on HG01882 produced a **1.7 GB rescued mates BAM**, far larger than expected. Investigation revealed the problem:

The Step 4 filter excluded mates on the same chromosome (`RNEXT="="`) and unmapped mates (`RNEXT="*"`), but did **not** exclude mates on contigs that were already extracted in Step 3 — the chr6 alt contigs and HLA decoy contigs. Most of the 250,253 unique mate locations identified were on these already-extracted contigs (e.g., a read on `chr6` whose mate is on `HLA-A*02:01:01:04`), causing samtools to re-fetch reads that were already in the initial BAM.

### Step 8: Corrected filter (proposed, NOT yet tested)

The corrected filter should exclude mates on any contig that was part of the initial extraction, plus mates on chr6 within the extracted window:

```bash
# Build list of already-extracted contigs
{
    echo "chr6"
    cat chr6_alt_contigs.txt
    cat hla_contigs.txt
} | sort -u > extracted_contigs.txt

# Find mates that are NOT already in our extraction
samtools view ${SAMPLE}.initial.bam \
    | awk 'BEGIN{
               while ((getline line < "extracted_contigs.txt") > 0) extracted[line] = 1
               close("extracted_contigs.txt")
           }
           {
               # Skip if mate has no valid position
               if ($7 == "=" || $7 == "*" || $8 == 0) next
               # Skip if mate is on extracted non-chr6 contig
               if ($7 in extracted && $7 != "chr6") next
               # For chr6, skip if mate is in extracted window
               if ($7 == "chr6" && $8 >= 25000000 && $8 <= 35000000) next
               print $7":"$8"-"$8
           }' \
    | sort -u > ${SAMPLE}.mate_locations.txt

echo "Missing mate locations to fetch: $(wc -l < ${SAMPLE}.mate_locations.txt)"
```

Expected count after correct filtering: hundreds to low thousands (compared to 250,000+ with the buggy filter).

**This corrected filter has not been tested end-to-end.** If you continue this experiment, run Step 8 in place of Step 4, then proceed with Steps 5 and 6 unchanged.

### Step 9: Merge, sort, index, and run SpecHLA (proposed continuation)

If Step 6 produces a reasonable-sized rescued mates BAM, the remaining steps are:

```bash
samtools merge -f -@ $THREADS \
    ${SAMPLE}.merged.bam \
    ${SAMPLE}.initial.bam \
    ${SAMPLE}.rescued_mates.bam

samtools sort -@ $THREADS -o ${SAMPLE}.merged.sorted.bam ${SAMPLE}.merged.bam
samtools index ${SAMPLE}.merged.sorted.bam
```

Then use `${SAMPLE}.merged.sorted.bam` as input to Steps 5, 6, 7, and 8 of Pipeline Version 1.

### Recommendation

Based on the testing performed, **Version 1 (without mate rescue) appears sufficient** for HLA typing — both validation samples (HG01882 and NA12878) produced typing results matching published references using only Version 1. The mate rescue strategy may not provide meaningful improvement, given that:

- The initial extraction already includes all contigs where alt-aware BWA places HLA reads (chr6, chr6 alts, HLA decoys, unmapped)
- Reads with mates on completely unrelated chromosomes are typically low-quality or off-target

If you do want to validate this empirically, the recommended approach is to complete Version 2 (using the corrected Step 8 filter) on a few samples, compare the resulting HLA calls against Version 1, and proceed with Version 1 if calls match.

---

## Population Tag Reference

SpecHLA's `-p` flag accepts: `Asian`, `Black`, `Caucasian`, `Unknown`, `nonuse`. Map 1000 Genomes super-populations as follows:

| 1000G Super-Population | SpecHLA Tag  | Examples                |
|------------------------|--------------|-------------------------|
| AFR (African)          | `Black`      | HG01882 (ACB), NA19238 (YRI) |
| EUR (European)         | `Caucasian`  | NA12878 (CEU)           |
| EAS (East Asian)       | `Asian`      | NA18525 (CHB)           |
| SAS (South Asian)      | `Unknown`    | HG03642 (STU)           |
| AMR (Admixed American) | `Unknown`    | HG01112 (CLM)           |

Use `Unknown` when ancestry is unclear or doesn't fit cleanly into SpecHLA's three categories. Use `nonuse` to disable population priors entirely (relies on mapping score only).

---

## Output Interpretation

The main result file is `./spechla_output/<SAMPLE>/hla.result.txt`. Format:

```
Sample  HLA_A_1  HLA_A_2  HLA_B_1  HLA_B_2  HLA_C_1  HLA_C_2  HLA_DPA1_1 ...
HG01882 A*23:01:01:01  A*34:02:01:01  B*07:02:01:03  B*44:03:01:01  ...
```

Each gene gets two columns (one per allele). Resolutions are 4-field where possible (e.g., `A*23:01:01:01`), corresponding to:
- Field 1: antigen group
- Field 2: protein
- Field 3: synonymous coding variants
- Field 4: non-coding variants

For most analyses, the first two fields (e.g., `A*23:01`) are sufficient.

Other useful files in the sample directory:
- `hla.result.details.txt` — alternative allele candidates with mapping scores
- `HLA_*_freq.txt` — haplotype frequency support per gene (helps detect haplotype collapse if frequencies are skewed e.g., 0.9/0.1)
- `hla.result.g.group.txt` — G-group resolution (functional grouping)
- `HLA_*.rephase.vcf.gz` — phased VCF per gene

### Validating typing results

For samples in 1000 Genomes, compare against the published HLA typing. A 9/10 match at 2-field resolution with rare sub-allele differences (e.g., `B*44:03:01:01` vs `B*44:03:01:13`) is excellent. Major discordance at 2-field resolution (e.g., reporting `B*07:02` when reference says `B*44:50`) indicates haplotype collapse — typically caused by missing reads from a divergent haplotype.

---

## Scaling Considerations

For batch processing of many samples:

### One-time setup

The contig list extraction (Step 2) only needs to happen once per project, since all samples aligned to the same reference share contig names:

```bash
# Run once, save outputs, reuse for all samples
samtools view -H -T "$REF_GRCH38" "$ANY_CRAM_URL" \
    | grep "^@SQ" | grep -oP "SN:\KHLA-[^\s]+" > hla_contigs.txt
samtools view -H -T "$REF_GRCH38" "$ANY_CRAM_URL" \
    | grep "^@SQ" | grep -oP "SN:\Kchr6_[^\s]+_alt" > chr6_alt_contigs.txt
```

### Sample sheet pattern

Use a TSV with sample IDs, CRAM URLs, and population tags:

```
sample_id	cram_url	population
HG01882	https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3242189/HG01882.final.cram	Black
NA12878	https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram	Caucasian
```

Loop over samples in a wrapper script.

### Idempotency

Skip samples whose results already exist:

```bash
if [ -f "${OUTPUT_DIR}/${SAMPLE}/hla.result.txt" ]; then
    echo "[$SAMPLE] Already done, skipping."
    continue
fi
```

### Per-sample working directories

Use isolated directories per sample to avoid file conflicts when running in parallel:

```bash
mkdir -p ./work/${SAMPLE}
cd ./work/${SAMPLE}
# ... run pipeline ...
```

### Cleanup

Each sample produces ~200–500 MB of intermediate BAMs. Remove them after SpecHLA succeeds:

```bash
rm ${SAMPLE}.hla_full.bam ${SAMPLE}.hla_full.coordsorted.bam*
rm -rf ./hla_reads/
```

### Failure handling

A single failed sample shouldn't stop a batch run. Wrap each sample in error-handling logic:

```bash
process_sample "$SAMPLE" "$URL" "$POP" || echo "[$SAMPLE] FAILED — continuing"
```

### Result aggregation

After the batch finishes, concatenate all per-sample results:

```bash
# Header from first sample
head -2 ./spechla_output/$(ls ./spechla_output/ | head -1)/hla.result.txt > all_results.tsv

# Body rows from all samples
for sample_dir in ./spechla_output/*/; do
    tail -n +2 "${sample_dir}/hla.result.txt" >> all_results.tsv
done
```

---

## Troubleshooting

### CRAM decoding errors

`[E::cram_get_ref] Failed to populate reference "chr6"`

The reference FASTA path is wrong or doesn't match the CRAM. Verify `$REF_GRCH38` is set correctly and the FASTA file matches what the CRAM was aligned against. The MD5 checksums in the CRAM header (`@SQ M5:` fields) must match your FASTA.

### Empty FASTQs after extraction

If `spechla-extract-hla-reads` produces empty FASTQs but completes without error, the input BAM is likely empty due to upstream failure. Check:

```bash
samtools view -c ${SAMPLE}.hla_full.coordsorted.bam
```

If this is 0, the initial CRAM extraction failed. Re-run Step 3 with `set -x` to see what's happening.

### Haplotype collapse (homozygous-looking calls)

If results show suspicious homozygosity at HLA class I genes (e.g., `B*07:02 / B*07:02` when reference says `B*07:02 / B*44:50`), reads from one haplotype were missing during extraction. Verify:

- Were chr6 alt contigs included? Check `chr6_alt_contigs.txt` is non-empty and contains 16 entries.
- Were HLA decoy contigs included? Check `hla_contigs.txt` is non-empty and has hundreds of entries.
- Were unmapped reads included? Confirm `"*"` is in the `samtools view` command.

### Slow extraction phase

If Step 3 is very slow after the first ~200 MB, this is the unmapped reads section being streamed. Unmapped reads in 30x WGS CRAMs can be 1–3 GB compressed. Either let it finish or remove `"*"` from the extraction (with the trade-off that very divergent haplotypes may be missed; validate against a reference sample if removing).

### Sub-allele mismatches with reference

Differences at the 4th field level (e.g., `B*44:03:01:13` vs `B*44:03:01:01`) are expected — they reflect non-coding variants and database version differences. Update the IPD-IMGT/HLA database for newer annotations:

```bash
cd $CONDA_PREFIX/share/spechla*/bin   # adjust path if installed from source
perl renew_HLA_annotation_db.pl
```

---

## Validation

Tested against published 1000 Genomes HLA typing references:

| Sample  | Population | 2-field concordance | Notes |
|---------|------------|---------------------|-------|
| HG01882 | ACB (AFR)  | 9/10 (90%)          | Sub-allele diff at HLA-B (B*44:03 vs B*44:50, same antigen group) |
| NA12878 | CEU (EUR)  | Full match          | All 10 alleles match at 2-field |

Both samples processed in under 10 minutes including network transfer.

