# genomic-array-nextflow

A Nextflow pipeline that processes raw genomic microarray data (ChIP-chip and array CGH) into genome-coordinate signal tracks and, for ChIP-chip, called peaks.

## Overview

This pipeline converts raw, probe-level microarray intensity data into genome-referenced results for two genomic array assay types used in VEuPathDB's genomic data workflows:

- **ChIP-chip** (`chipChip`) — transforms raw probe intensities to genome coordinates using a probe-to-genome BAM mapping, calls peaks and smoothed signal profiles with a custom Java peak finder, merges overlapping peak regions with BEDTools, and produces bigWig tracks plus a study-config file describing the peak results for loading into VEuPathDB's data warehouse.
- **Array CGH** (`cghArray`) — transforms raw probe intensities to genome coordinates, merges overlapping regions with BEDTools, and produces bigWig tracks of the merged signal.

Both assay types share the same raw-data-to-genome-coordinates and BEDTools-merge steps; ChIP-chip additionally runs peak calling and produces study-config output describing the results for downstream loading.

## Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2)
- A container engine: [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/)/[Apptainer](https://apptainer.org/)

Processes run from a mix of published containers (`veupathdb/bioperl`, `biocontainers/bedtools`, `biocontainers/tabix`, `biocontainers/ucsc-bedgraphtobigwig`) and a pipeline-specific container, `veupathdb/genomicarray:1.0.0`, built from the included `Dockerfile`. That image compiles `src/java/org/apidb/ggtools/array/ChIP_Chip_Peak_Finder.java` into `ChIPChipPeakFinder.jar`, which performs ChIP-chip peak finding and signal smoothing.

Execution profiles are provided under `conf/`:
- `conf/docker.config` — run with Docker (default, included via `nextflow.config`)
- `conf/singularity.config` — run with Singularity
- `conf/lsf.config` — run with Singularity on an LSF cluster

## Usage

The pipeline has a single, unnamed entry point (no `-entry` flag is needed); the assay is selected with `params.assayType`. Samples are read from a samplesheet (`params.samplesheetFileName`, a two-column CSV of sample id and raw data file name) located in `params.input`.

ChIP-chip:

```bash
nextflow run VEuPathDB/genomic-array-nextflow -r main -resume \
  --assayType chipChip \
  --input /path/to/input \
  --samplesheetFileName samplesheet.csv \
  --platformBamFile /path/to/probes.bam \
  --seqSizeFile /path/to/seqSizes.txt \
  --profileSetName "My ChIP-chip profile set" \
  --outDir /path/to/output \
  -C conf/docker.config
```

Array CGH:

```bash
nextflow run VEuPathDB/genomic-array-nextflow -r main -resume \
  --assayType cghArray \
  --input /path/to/input \
  --samplesheetFileName samplesheet.csv \
  --platformBamFile /path/to/probes.bam \
  --seqSizeFile /path/to/seqSizes.txt \
  --outDir /path/to/output \
  -C conf/docker.config
```

## Key Parameters

| Parameter              | Default              | Description |
| ----------------------- | --------------------- | ------------ |
| `assayType`              | `cghArray`             | Assay to run: `chipChip` or `cghArray`. |
| `input`                  | `input/`                | Directory containing the samplesheet and raw per-sample data files. |
| `samplesheetFileName`    | `samplesheet.csv`       | CSV listing samples: sample id in column 1, raw data file name (relative to `input`) in column 2. |
| `platformBamFile`        | `probes.bam`            | BAM of probe-to-genome alignments (e.g. produced by `array-probes-nextflow`) used to map raw probe data to genome coordinates. |
| `seqSizeFile`            | `input/seqSizes.txt`    | Chromosome/sequence size file used when writing bigWig output. |
| `outDir`                 | `output`                | Directory results are published to. |
| `peakFinderArgs`         | `NA`                    | Extra arguments passed to the ChIP-chip peak finder Java tool (`chipChip` only). |
| `profileSetName`         | `test profile set name` | Profile set name recorded in the study-config output for ChIP-chip peak results. |

## Output

- `outDir/*.bw` — bigWig tracks of merged, genome-coordinate signal for every sample (both assay types).
- For `chipChip`:
  - `outDir/*.peaks.bed.gz` and `.tbi` — bgzipped, tabix-indexed BED of called peaks per sample.
  - `outDir/insert_study_results` — collected study-config file describing the peak results (name, protocol, source ID type, profile set) for downstream loading.
