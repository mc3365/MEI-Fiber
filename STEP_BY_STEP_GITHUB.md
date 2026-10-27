# Step-by-step MEI-Fiber workflow

This guide starts with an aligned ONT or PacBio Fiber-seq BAM and ends with a
queryable MEI-Fiber database and example figures. Extraction is platform-specific;
HDF5 construction, annotations, indexed queries, and shared visualizations use the
same workflow after extraction.

PacBio FIRE modeling remains an upstream process. MEI-Fiber can extract and package
an existing FIRE-annotated BAM and can run the downstream co-accessibility workflow.
The detailed FIRE preparation, ranking, and validation steps remain in
[`docs/tutorial.md`](docs/tutorial.md#fire-co-accessibility-preparation).

For a complete list of command-line and YAML options, see
[`docs/parameters.md`](docs/parameters.md).

## 1. Clone the repository

Choose a directory where you keep code projects, then clone MEI-Fiber:

```bash
git clone https://github.com/mc3365/MEI-Fiber.git
cd MEI-Fiber
```

To update an existing clone:

```bash
git checkout main
git pull --ff-only
```

## 2. Install MEI-Fiber and the long-read tools

The shared Conda environment includes the versions used during validation:

- samtools 1.22.1
- modkit 0.5.0
- fibertools-rs 0.8.0

`samtools` and `fibertools-rs` are used by both workflows. `modkit` is used by ONT
extraction and can remain installed in PacBio-only projects.

Create the environment and install MEI-Fiber:

```bash
conda env create -f environments/mei_fiber.yml
conda activate mei-fiber
pip install -e ".[viz]"
```

Confirm that the command-line tools are available:

```bash
mei-fiber --version
samtools --version | head -1
modkit --version
ft --version
```

On an HPC system, the Python package and long-read tools may live in separate
environments. In that case, use absolute executable paths in the YAML configuration:
`ft_executable` and `samtools_executable` for both platforms, plus
`modkit_executable` for ONT.

## 3. Check the input BAM

Every input BAM must:

- be aligned and coordinate sorted;
- have a BAM index (`.bai`);
- use chromosome names compatible with the reference FASTA.

ONT BAMs must contain paired MM/ML modification tags. PacBio BAMs must be compatible
with fibertools extraction; use a FIRE-annotated `*.fire.bam` when the optional
`fire_accessibility` layer and downstream co-accessibility analysis are required.

MEI-Fiber checks these requirements before extraction. Basic checks can also be run
directly:

```bash
samtools quickcheck sample.bam
samtools view -H sample.bam | grep '^@HD'
samtools view sample.bam | head -1
```

The BAM and reference FASTA must describe the same genome assembly. For the validated
mouse workflow, the reference path was an mm10 FASTA.

## 4. Prepare a platform-specific configuration

Copy the template matching the input platform:

```bash
# ONT
cp configs/ont_template.yaml configs/my_ont.yaml

# PacBio
cp configs/pacbio_template.yaml configs/my_pacbio.yaml
```

Edit these shared fields in the selected YAML:

- `output_dir`: directory for the HDF5 database;
- `output_file`: HDF5 filename;
- `reference`: reference FASTA used for alignment;
- `samples[].name`: short, unique sample name;
- `samples[].bam`: aligned, modification-tagged BAM;
- `samples[].layers`: paths for the platform's extraction outputs;
- `annotations.master`: optional shared BED annotations; and
- `extraction`: threads, overwrite behavior, and executable paths.

For ONT, the 5mC and 5hmC entries normally point to the same modkit TSV because both
call types are stored in that file. Keep `methylation_threshold: 0.5` to reproduce
the validated binarization rule. Raw probabilities are retained in the HDF5 database.

For PacBio, the extractor normalizes fibertools nucleosome, MSP, 6mA, and 5mC outputs
into the same builder schema. PacBio does not produce the current 5hmC layer. Add
`fire_accessibility` only when starting from a FIRE-annotated BAM or an existing
`ft fire --extract` BED.

For a second sample, add another item under `samples` with its own BAM and output
paths. Sample names become top-level groups in the database.

Common YAML settings:

| Field | Typical value | When to change it |
| --- | --- | --- |
| `parameters.methylation_threshold` | `0.5` | Binary 5mC/5hmC threshold; keep for the validated ONT behavior and cross-platform consistency. |
| `build.build_spatial_index` | `true` | Keep enabled for fast coordinate queries. |
| `extraction.threads` | `16` | Match available CPU cores. |
| `extraction.overwrite` | `false` | Set to `true` only when regenerating extraction products. |
| `extraction.validate_bam` | `true` | Keep enabled unless doing a very controlled rerun. |
| `extraction.write_manifest` | `true` | Keep enabled for provenance and troubleshooting. |
| `extraction.modkit_executable` | `modkit` or absolute path | ONT only. Use an absolute path on HPC if needed. |
| `extraction.ft_executable` | `ft` or absolute path | Required by both platform extractors. |
| `extraction.samtools_executable` | `samtools` or absolute path | Used for BAM validation on both platforms. |

## 5. Run a small extraction test

Before processing a full BAM, it is useful to test a small genomic interval. This
checks the environment, reference, BAM tags, and output paths at low cost:

```bash
samtools view -@ 4 -h -b sample.bam chr1:3000000-5000000 \
  | samtools sort -@ 4 -o sample_test.bam
samtools index sample_test.bam
samtools view -c sample_test.bam
```

Create a second YAML file that points to `sample_test.bam` and writes to a separate
test directory. Then run the matching extractor:

```bash
# ONT
mei-fiber extract \
  --platform ont \
  --config configs/ont_smoke_test.yaml \
  --samples sample_test

# PacBio
mei-fiber extract \
  --platform pacbio \
  --config configs/pacbio_smoke_test.yaml \
  --samples sample_test
```

Extraction command options:

| Option | Meaning |
| --- | --- |
| `--platform ont/pacbio` | Select the extractor matching the input BAM. |
| `--config` | YAML configuration file. |
| `--samples` | Optional sample filter. Repeat the flag to process multiple selected samples. |

Successful ONT extraction produces:

- a shared modkit TSV for 5mC and 5hmC;
- BED output for 6mA;
- BED output for MSPs;
- a flattened nucleosome CSV; and
- `MEI-Fiber_manifest_<sample>.json` with BAM QC, tool versions, commands, and output
  sizes.

Successful PacBio extraction produces normalized 5mC, 6mA, MSP, and nucleosome files
for the same HDF5 builder. When configured from a FIRE-annotated BAM, it also produces
the FIRE accessibility BED. Both platforms write the extraction manifest.

Inspect the manifest before continuing:

```bash
python -m json.tool /path/to/MEI-Fiber_manifest_sample_test.json | less
```

## 6. Extract the full samples

Run full extraction on a compute node or submit it as a Slurm job. To process every
sample in the YAML, choose one command:

```bash
# ONT
mei-fiber extract --platform ont --config configs/my_ont.yaml

# PacBio
mei-fiber extract --platform pacbio --config configs/my_pacbio.yaml
```

To process selected samples, repeat `--samples`. For example:

```bash
mei-fiber extract \
  --platform ont \
  --config configs/my_ont.yaml \
  --samples d0 \
  --samples d4
```

With `overwrite: false`, valid existing outputs are reused. Set `overwrite: true` only
when the extraction products should be regenerated.

## 7. Build the HDF5 database

Run the unit tests before a full build when developing from a Git checkout:

```bash
python -m pytest tests/test_utils.py tests/test_builder.py -q
```

Build every configured sample. The command is identical after either extraction
workflow:

```bash
CONFIG=configs/my_ont.yaml       # or configs/my_pacbio.yaml
mei-fiber build --config "$CONFIG"
```

Build command options:

| Option | Meaning |
| --- | --- |
| `--config` | YAML file describing extracted layer paths and database output. |
| `--samples` | Optional sample filter. Repeat to build multiple selected samples. |
| `--build-index` | Build the spatial-index sidecar after HDF5 creation. |
| `--no-build-index` | Skip spatial-index creation for this run. |

The builder creates `output_dir/output_file`. When `build.build_spatial_index: true`
is set in the YAML, it also creates a cached spatial-index sidecar named
`<output>.index.pkl` for fast region queries. Keep that file beside the HDF5 database.
The builder opens the HDF5 path for writing, so use a new filename or preserve the
previous HDF5 file before rebuilding.

Inspect the completed database:

```bash
mei-fiber info /path/to/output/fiber_database.h5
ls -lh /path/to/output/fiber_database.h5
ls -lh /path/to/output/fiber_database.index.pkl
```

Confirm that the expected samples and annotations appear and that the fiber counts are
consistent with the extraction inputs.

## 8. Query a region

Use the command line for a quick overlap check:

```bash
mei-fiber query \
  --db /path/to/output/fiber_database.h5 \
  --region chr1:3000000-5000000 \
  --sample d0
```

Query command options:

| Option | Meaning |
| --- | --- |
| `--db` | HDF5 database path. |
| `--region` | Region in `chrom:start-end` form. |
| `--sample` | Optional sample filter. |
| `--out` | Optional text file for overlapping fiber IDs. |

For repeated regional analysis, create the spatial index once:

```python
from mei_fiber.db import FiberDatabase

with FiberDatabase("/path/to/output/fiber_database.h5", build_index=True) as db:
    print(db.get_summary())
```

The index is saved beside the HDF5 database and loaded by later sessions.

## 9. Plot individual molecules

The example plot shows nucleosomes, MSPs, and one modification layer across fibers
overlapping a region:

```bash
python examples/region_plot.py \
  --db /path/to/output/fiber_database.h5 \
  --region chr1:3000000-3050000 \
  --sample d0 \
  --layer 5mC \
  --tracks full \
  --max-fibers 100 \
  --hide-unmodified \
  --out d0_chr1_region.png
```

Use a focused interval, usually 10-50 kb, so individual molecular patterns remain
readable. The same command accepts `5hmC` or `6mA` as the modification layer.

Region plot options:

| Option | Meaning |
| --- | --- |
| `--db` | HDF5 database path. |
| `--region` | Focused interval in `chrom:start-end` form. |
| `--sample` | Optional sample name. |
| `--layer` | Modification layer: `5mC`, `5hmC`, or `6mA`. |
| `--tracks` | `full` shows nucleosome/MSP plus the selected modification; `modification` shows only the selected modification; `chromatin` shows nucleosome/MSP only. |
| `--max-fibers` | Maximum number of fiber rows shown. |
| `--hide-unmodified` | Hide unmodified 5mC/5hmC calls for a cleaner modification-focused view. |
| `--out` | Output figure path. |

## 10. Plot summary figures

Global ECDF plots summarize per-fiber feature fractions:

```bash
python examples/feature_ecdf.py \
  --db /path/to/output/fiber_database.h5 \
  --samples d0 d4 \
  --outdir figures/ecdf
```

Useful options:

| Option | Meaning |
| --- | --- |
| `--samples` | One or more samples. Omit to include all samples. |
| `--max-fibers-per-chrom` | Optional downsampling cap for exploratory plots. |

For annotation-body heatmaps and metaplots:

```bash
python examples/annotation_heatmap.py \
  --db /path/to/output/fiber_database.h5 \
  --annotation CGI \
  --sample d0 \
  --max-regions 50 \
  --min-fibers 10 \
  --outdir figures/CGI_d0_test
```

Omit `--max-regions` for full annotation plots. Keep it for smoke tests.

For fixed-window center/TSS heatmaps and metaplots:

```bash
python examples/centered_heatmap.py \
  --db /path/to/output/fiber_database.h5 \
  --bed /path/to/master_annotations_v4.uniqueID.bed \
  --annotation CGI \
  --samples d0 d4 \
  --methylation-display binned \
  --min-fibers 10 \
  --outdir figures/CGI_center
```

Centered plot options:

| Option | Meaning |
| --- | --- |
| `--annotation` | Annotation class, such as `CGI`. |
| `--samples` | One or more samples to plot together in the metaplot. |
| `--bed` | Optional 9-column annotation BED for strand-aware TSS/gene orientation. |
| `--max-regions` | Optional smoke-test limit. Omit for full figures. |
| `--min-fibers` | Minimum spanning fibers required for a region. Default is `10`. |
| `--methylation-display` | `binned` for 50 bp bins or `smoothed` for smoothed base-resolution 5mC. |

For more detail on plotting outputs and defaults, see
[`docs/parameters.md`](docs/parameters.md#plotting).

## 11. Run benchmarks

Use the annotation benchmark for a biologically meaningful speed test:

```bash
python -m mei_fiber.benchmark.annotation \
  --db /path/to/output/fiber_database.h5 \
  --sample d0 \
  --annotation CGI \
  --layers nucleosomes 5mC msp \
  --max-regions 50 \
  --repeats 3 \
  --warmups 1 \
  --outdir benchmark/d0_CGI
```

Use the random-region benchmark for technical scaling across window sizes:

```bash
python -m mei_fiber.benchmark.database \
  --db /path/to/output/fiber_database.h5 \
  --config configs/my_platform.yaml \
  --sample d0 \
  --outdir benchmark/d0 \
  --region-sizes-kb 1 10 100 1000 \
  --n-regions 25 \
  --repeats 3 \
  --query-modes spatial_index \
  --build-index
```

For full production databases, avoid the slow array-scan baseline unless you
specifically need it. Use `--query-modes spatial_index` for routine benchmarking.
The database benchmark is platform-neutral. ONT annotation benchmarks can add
`5hmC`; PacBio benchmarks can add `6mA` or `fire_accessibility` when those layers
exist in the database.

Benchmark options and outputs are described in
[`docs/benchmark.md`](docs/benchmark.md) and
[`docs/parameters.md`](docs/parameters.md#random-region-benchmark).

## 12. Verify a development checkout

Contributors should install the development dependencies and run the complete checks:

```bash
pip install -e ".[dev,viz]"
python -m pytest -q
python -m ruff check src tests examples
```

GitHub Actions runs the automated tests after each push. The local checks should pass
before changes are committed.

## Troubleshooting

### `conda: command not found`

Initialize Conda for the current shell or use the full path to the environment's
executables. On batch systems, shell startup files may not be loaded by Slurm.

### An executable is not found

Activate the `mei-fiber` environment, add its `bin` directory to `PATH`, or set the
required executable paths explicitly in the YAML configuration.

### `gzip.BadGzipFile`

Some older extraction outputs have a `.gz` suffix but contain plain text. MEI-Fiber can
read both forms, but new modkit outputs should be checked with:

```bash
file sample_raw_mods.tsv.gz
gzip -t sample_raw_mods.tsv.gz
```

### Region queries report that no spatial index was found

The database is still usable, but queries fall back to an array scan. Open the database
once with `build_index=True` as shown in Step 8.

### GitHub authentication fails

Refresh the GitHub CLI login, then retry the push:

```bash
gh auth login -h github.com
git push
```

## Development history

MEI-Fiber began as a Python package skeleton and was filled in against the most recent
legacy V8 scripts. The database reader was first checked for output parity, followed by
the builder and then the ONT extraction workflow. Real d0 and d4 data were used to
validate the HDF5 build, regional queries, modification binarization, BAM QC, manifests,
and extraction outputs. PacBio normalization, the shared schema, FIRE accessibility,
and downstream co-accessibility were subsequently validated against legacy outputs.
The validated ONT milestone is tagged `v0.4.0-ont`.
