# Querying and Annotations

MEI-Fiber uses one query API for ONT and PacBio databases. The installed Python
package supplies the query functions; the HDF5 file supplies the samples,
chromosome-partitioned feature arrays, per-fiber slice indices, and optional shared
annotations. A neighboring `.index.pkl` file supplies the spatial index used by fast
coordinate queries. Queries still work without that sidecar, but fall back to an
array scan.

## Choose the query that matches the question

| Question | Interface | Result |
| --- | --- | --- |
| Which fibers overlap one coordinate interval? | `mei-fiber query` or `get_fibers_at()` | Fiber IDs |
| What is stored for one fiber? | `get_nucleosomes()`, `get_msp()`, `get_methylation()`, `get_fire_accessibility()` | NumPy arrays for that fiber |
| Which annotation classes are available? | `list_annotations()` | Annotation-class names |
| Where are all regions of one class? | `get_annotation_regions_with_ids()` | Coordinates and annotation IDs |
| What molecular features overlap one stored region? | `query_by_id()` | One row per overlapping fiber |
| What molecular features overlap every region in a class? | `query_annotation_fast()` | One row per `(region, fiber)` pair |
| How does one annotation class compare across samples? | `compare_samples()` | Concatenated per-fiber table with a `sample` column |
| Which fibers span most of a region? | `load_chromosome_data()` plus `get_spanning_fiber_ids()` | Integer fiber IDs passing a coverage threshold |

`query_annotation()` is a compatibility alias for `query_annotation_fast()`.

## Inspect a database first

```bash
mei-fiber info /path/to/fiber_database.h5
```

In Python:

```python
from mei_fiber import FiberDatabase

with FiberDatabase("/path/to/fiber_database.h5") as db:
    print("samples:", db.samples)
    print("annotations:", db.list_annotations())
```

Always provide a sample when a database contains more than one sample. Python query
methods default to the first sample when `sample` is omitted.

## Query one coordinate interval

The command-line query returns the IDs of fibers whose aligned spans overlap a
0-based, half-open genomic interval:

```bash
mei-fiber query \
  --db /path/to/fiber_database.h5 \
  --sample d0 \
  --region chr1:3000000-3050000 \
  --out d0_chr1_fiber_ids.txt
```

The equivalent Python query is:

```python
from mei_fiber import FiberDatabase

with FiberDatabase("/path/to/fiber_database.h5") as db:
    fibers = db.get_fibers_at("chr1", 3_000_000, 3_050_000, sample="d0")
    print(len(fibers), fibers[:5])
```

This operation finds fibers only. It does not summarize their molecular layers.
Use the per-fiber accessors or an annotation query for feature summaries.

## Inspect one fiber

```python
from mei_fiber import FiberDatabase

with FiberDatabase("/path/to/fiber_database.h5") as db:
    chrom = "chr1"
    sample = "d0"
    fibers = db.get_fibers_at(chrom, 3_000_000, 3_050_000, sample=sample)

    if fibers:
        fiber_id = fibers[0]
        print(db.get_nucleosomes(fiber_id, chrom, sample=sample))
        print(db.get_msp(fiber_id, chrom, sample=sample))
        print(db.get_methylation(fiber_id, chrom, "5mC", sample=sample))
        print(db.get_methylation(fiber_id, chrom, "6mA", sample=sample))
```

For a PacBio database built with `fire_accessibility`, also use:

```python
print(db.get_fire_accessibility(fiber_id, chrom, sample=sample))
```

PacBio 5hmC is not produced by the current extraction pathway. ONT and PacBio share
the same accessors for layers that are present in both databases.

## Query an annotation class

The following example finds fibers overlapping all stored CGIs and computes
per-fiber chromatin and modification summaries:

```python
from mei_fiber import FiberDatabase

with FiberDatabase("/path/to/fiber_database.h5") as db:
    result = db.query_annotation_fast(
        "CGI",
        sample="d0",
        feature_types=["nucleosomes", "msp", "5mC", "5hmC", "6mA"],
    )

print(result.head())
result.to_csv("d0_CGI_per_fiber.csv", index=False)
```

The output contains one row for every `(annotation region, overlapping fiber)` pair.
Depending on the requested layers, columns include:

- nucleosomes: `n_nucleosomes`, `mean_nuc_width`, `mean_linker`, and
  `mean_repeat_length`;
- 5mC: `n_cpg` and `pct_methylated`;
- 5hmC: `n_5hmC` and `pct_hydroxymethylated`;
- 6mA: `n_6mA`;
- MSPs: `n_msp` and `mean_msp_width`; and
- FIRE accessibility: `n_fire_accessibility`, `min_fire_score`,
  `mean_fire_score`, and `fire_overlap_bp`.

Use `max_regions=50` for a smoke test. Omit it for the full annotation class.

## Compare an annotation across samples

`compare_samples()` runs the same annotation query for every sample in the database:

```python
from mei_fiber import FiberDatabase

with FiberDatabase("/path/to/fiber_database.h5") as db:
    result = db.compare_samples(
        "CGI",
        feature_types=["nucleosomes", "msp", "5mC"],
    )

summary = result.groupby("sample").agg(
    fibers=("fiber_id", "nunique"),
    median_msp_count=("n_msp", "median"),
    median_nucleosome_count=("n_nucleosomes", "median"),
    median_pct_methylated=("pct_methylated", "median"),
)
print(summary)
```

For ONT, MSP count and width are the main packaged accessibility summaries. For
PacBio, MSPs are also available; databases built with the optional FIRE layer can
add `fire_accessibility` to `feature_types` and compare FIRE scores and overlap.

These summaries describe the fibers returned by each annotation query. Formal
differential testing, covariate adjustment, and biological interpretation remain
downstream analysis decisions.

## Query one stored annotation

First list annotation IDs:

```python
from mei_fiber import FiberDatabase

with FiberDatabase("/path/to/fiber_database.h5") as db:
    regions = db.get_annotation_regions_with_ids("Promoter")
    print(regions[:5])
```

Then query one ID:

```python
with FiberDatabase("/path/to/fiber_database.h5") as db:
    result = db.query_by_id(
        "Promoter_42",
        sample="d0",
        feature_types=["nucleosomes", "msp", "5mC"],
    )
```

## Query a specific gene

MEI-Fiber currently stores annotation coordinates and IDs, but not the optional gene
name column from the source BED. Gene-symbol lookup therefore begins in the original
annotation file. For the 9-column format described below, column 5 is the gene name:

```bash
MASTER=/path/to/master_annotations.uniqueID.bed

awk -F '\t' '$5 == "Sox2" {print $1, $2, $3, $4, $5, $6, $7}' \
  OFS='\t' "$MASTER"
```

Use the returned column-4 ID with `query_by_id()`. A promoter record asks about the
TSS-proximal interval; a gene-body record asks a different, gene-wide question.

For a coordinate-only annotation file, locate the gene in the source GTF/GFF/BED,
then use `get_fibers_at()` or the CLI coordinate query.

## Build an annotation file

### Required BED4 fields

The HDF5 builder reads the first four tab-delimited columns:

| Column | Meaning | Example |
| --- | --- | --- |
| 1 | Chromosome | `chr3` |
| 2 | Start, 0-based inclusive | `100000` |
| 3 | End, 0-based exclusive | `102500` |
| 4 | Unique annotation ID | `Promoter_17` |

Every ID should be unique and end in `_<integer>`. MEI-Fiber removes the final
numeric suffix to derive the annotation class. For example:

- `CGI_12` belongs to annotation class `CGI`;
- `Typical_Enhancer_51` belongs to `Typical_Enhancer`; and
- `Sox2_Promoter_1` belongs to `Sox2_Promoter`.

Example BED4:

```text
chr1	3531624	3531843	CGI_1
chr1	3669498	3671998	Promoter_1
chr3	100000	102500	GeneA_Promoter_1
```

Use coordinates from the same reference assembly as the aligned BAM. The `GeneA`
record above is illustrative rather than a reference annotation.

### Recommended BED9 fields

For strand-aware promoter/TSS heatmaps, retain nine columns:

| Column | Meaning |
| --- | --- |
| 1-4 | BED4 fields above |
| 5 | Gene name, or `.` |
| 6 | Strand encoded as `1`, `-1`, or `0` |
| 7 | TSS coordinate, or `-1` |
| 8 | Transcript ID, or `.` |
| 9 | External transcript ID, or `.` |

Example:

```text
chr3	100000	102500	Promoter_42	GeneA	1	100500	ENSMUST...	NM_...
```

Only columns 1-4 are embedded in the current HDF5 schema. Keep the original BED and
pass it to `examples/centered_heatmap.py --bed ...` when strand-aware orientation is
required.

### Validate the BED

Check the column count, interval validity, duplicate IDs, and sorting before build:

```bash
awk 'BEGIN{FS=OFS="\t"} NF < 4 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $2 >= $3 {
  print "invalid annotation row", NR, $0 > "/dev/stderr"; bad=1
} END{exit bad}' master_annotations.bed

cut -f4 master_annotations.bed | sort | uniq -d | head

LC_ALL=C sort -c -k1,1 -k2,2n master_annotations.bed
```

An empty duplicate-ID command and a successful `sort -c` indicate that those checks
passed. Sort the file if needed:

```bash
LC_ALL=C sort -k1,1 -k2,2n master_annotations.bed > master_annotations.sorted.bed
```

### Add annotations during build

Point the YAML configuration to the annotation file:

```yaml
annotations:
  master: /absolute/path/master_annotations.sorted.bed
```

Then build the database:

```bash
mei-fiber build --config /absolute/path/config.yaml
```

Annotations are incorporated during HDF5 construction. MEI-Fiber does not currently
modify an existing HDF5 in place, so adding or changing annotations requires a
rebuild. Confirm the result with:

```python
from mei_fiber import FiberDatabase

with FiberDatabase("/path/to/fiber_database.h5") as db:
    print(db.list_annotations())
    print(db.get_annotation_regions_with_ids("CGI")[:3])
```

## Platform compatibility

Coordinate, per-fiber, annotation, and cross-sample queries are platform-neutral.
They operate on whichever layers are present in the HDF5:

| Layer | ONT | PacBio |
| --- | --- | --- |
| Nucleosomes | yes | yes |
| MSPs | yes | yes |
| 5mC | yes | yes |
| 5hmC | yes when extracted | not produced by the current pathway |
| 6mA | yes | yes |
| FIRE accessibility | no | optional |

Requesting a layer that was not built cannot create that information. Inspect the
configuration and use a small test query before running a genome-wide annotation
analysis.
