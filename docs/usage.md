# SparK usage guide

The [README](../README.md) covers installation and the shortest commands. All examples below assume the environment is activated and inputs use the same genome assembly.

## Files and output

`-cf` accepts control BigWig (`.bw`, `.bigWig`) or bedGraph files. `-tf` adds treatment files. No conversion from BigWig is needed. SparK plots supplied signal; prepare coverage and normalization upstream.

Every successful plot writes:

- `.svg`: editable vector master.
- `.jpg`: white-background image at 300 DPI and quality 95 by default.
- `.run.json`: code version/hash, settings, resolved coordinates, dependency versions, input identifiers, and output hashes.

Use `--formats pdf png` to add a vector PDF and a PNG. `--dpi 600` changes raster resolution; `--jpg-quality 90` changes JPEG compression. `--no-jpg` skips JPG; `--no-run-record` skips the run record. SVG is retained for all format selections. Outputs are staged before replacing previous files.

Run records identify inputs by path, size, modification/change times, and inode by default. `--hash-inputs` also computes full SHA-256 hashes. Hashes are cached against file metadata, but the first hash of a large input requires reading the whole file. Plotting and indexing never copy or modify genomic inputs.

## Saved settings

Add `--save-config settings.json` to a working command, then use:

```bash
python SparK.py --config settings.json
python SparK.py --config settings.json -pr chr1:20000-30000 -o another_figure
```

Explicit command-line options override saved values. Relative input and output paths in a JSON file are resolved against that file's folder. Unknown settings are errors.

For multiple named conditions, use a `tracks` list:

```json
{
  "region": "chr1:10000-20000",
  "show_plots": "averages",
  "output_name": "comparison",
  "tracks": [
    {"file": "control_1.bw", "condition": "Control", "group": "H3K27ac", "color": "#005CFF"},
    {"file": "control_2.bw", "condition": "Control", "group": "H3K27ac", "color": "#005CFF"},
    {"file": "drug_a.bw", "condition": "Drug A", "group": "H3K27ac", "color": "#DD5500"},
    {"file": "drug_b.bw", "condition": "Drug B", "group": "H3K27ac", "color": "#009E73"}
  ]
}
```

Files with the same `condition` and `group` are replicates. Each group is a separate panel. A track can also include `label` and a positive `scale_factor`; the factor is applied before replicate averaging. There is no hidden library-size normalization.

For a short command, repeat `--track CONDITION=FILE`:

```bash
python SparK.py --track Control=c1.bw --track Control=c2.bw \
  --track DrugA=a.bw --track DrugB=b.bw -ps averages -pr chr1:10000-20000
```

## Existing control/treatment commands

The original options continue to work:

```bash
python SparK.py -cf c1.bw c2.bw -tf t1.bw t2.bw \
  -cg 1 1 -tg 1 1 -ps averages -l Control Treated -pr chr1:10000-20000
```

Supply exactly one group number per file whenever using `-cg`/`-tg`. With standard plots and `-ps all`, inputs get separate panels by default. Means, STD, and sine plots default to one panel.

`-pt standard` draws individual signals or means. `-pt STD` draws mean ± population SD and requires at least two replicates per displayed condition. `-pt sine` mirrors two conditions on a shared scale and supports individual replicates or means. Use standard/STD for three or more conditions. Signed inputs retain their sign before control-side mirroring.

## Colors and scales

- `-f 333333`: one color for all signals.
- `-f 005CFF DD5500`: control first, treatment second.
- `--track_colors FF0000 00AA00 0000FF`: one color per input, controls first then treatments.
- More than two colors with `-f` also means per-input colors.
- Hex colors accept three or six digits. Quote colors starting with `#` in a shell.
- Presets: `blue/red`, `blue/grey`, `all_grey`, `blue/green`.

For named tracks, put colors in the JSON track entries or use a `condition_colors` mapping. Replicates combined into a mean/SD band must share a color. Named-condition legends use condition names; `condition_labels` can rename them. The legacy `-l` legend uses the first control/treatment colors; use `-gl` for individual panel labels.

Each panel scales independently. `-gs yes` shares a scale across panels; `-es 2` lets panel 2 keep its own scale; `-eg 2` excludes it. `-cs 10 D 100` sets the normalization magnitudes for panels 1 and 3 and leaves panel 2 automatic. Supply an entry for every original group number, including excluded groups. Standard signals occupy 90% of track height; axis labels include this margin. Sine ticks indicate the chosen positive/negative magnitude directly.

`--max_points 2000` controls spatial downsampling. `--point_stat average` uses overlap-weighted bin means; `--point_stat max` retains each bin's numeric maximum, useful for positive spikes. It does not change replicate averaging. `-sm 10` applies a trailing 10-bin mean, padding the start with the first value.

## Gene lookup and batch plotting

```bash
python SparK.py -cf sample.bw --gene GAPDH -gtf genes.gtf --flank 10000 -o GAPDH
python SparK.py -cf sample.bw --genes-file genes.txt -gtf genes.gtf -o gene_figures
python SparK.py -cf sample.bw --regions-file regions.bed --batch-shared-scale -o region_figures
```

`genes.txt` contains one exact gene name or gene ID per line. Ambiguous names require a unique ID or explicit coordinates. Gene queries cover the annotated gene body plus the requested flank. GTF coordinates are converted from one-based inclusive to zero-based half-open coordinates.

Batch regions accept BED columns `chromosome start end [name]`, or text lines `chr1:10000-20000 [name]`. Blank lines/comments are ignored. Batch output includes a browsable `index.html`, per-figure files/records, and `batch.run.json`. Duplicate names get numbered filenames. Per-region errors are reported while other jobs continue; the command returns nonzero if any fail.

`--batch-shared-scale` keeps corresponding panels comparable across the batch. Explicit `-cs` values take priority. JSON configurations can alternatively contain a `regions` list of coordinate strings or a `genes` list.

## Difference and ratio tracks

```bash
python SparK.py --config comparison.json --transform difference --compare "Drug A" Control
python SparK.py --config comparison.json --transform log2ratio --compare "Drug A" Control --pseudocount 1
```

Comparisons use condition means after applying explicit scale factors:

- `difference`: target mean − reference mean.
- `ratio`: `(target mean + pseudocount) / (reference mean + pseudocount)`.
- `log2ratio`: log₂ of that ratio.

The default pseudocount is 1 in the input signal's units. Choose it deliberately for your normalization. Undefined/non-finite ratios are errors. Difference supports signed values; ratio/log₂ inputs must satisfy their nonnegative/positive requirements after adding the pseudocount. Transforms use standard plots without sparks. With more than two conditions, `--compare TARGET REFERENCE` is required.

`-sp yes` adds the existing descriptive mean/SD highlights and requires at least two replicates of both comparison conditions. With multiple conditions, select the pair with `--compare`. This heuristic is not a significance test or p-value. `-sc` sets treatment-up/target-up then control-up/reference-up colors.

## Annotations and appearance

| Option | Purpose |
| --- | --- |
| `-bed file1.bed file2.bed` | BED interval tracks |
| `-bedcol`, `-bedlab` | BED colors and labels |
| `-gtf genes.gtf` | GTF annotations; gzip accepted |
| `-dg GAPDH IFFO1` | Filter annotation labels inside the selected region; use `--gene` to find coordinates |
| `-dt all` or `-dt ID1 ID2` | All or selected transcripts |
| `-wg yes` | Gene names for transcript labels |
| `-tss no` | Hide legacy translation start-codon arrows |
| `-genestart yes` | Show annotated gene/transcript-start arrows |
| `-dc top_left` | Region-label position; also top_right, bottom_left, bottom_right, no |
| `-w`, `--x_scale`, `--y_scale` | Figure width and scaling |
| `-scale no` | Hide the genomic scale bar |

GTF requires nine tab-separated columns, feature type in column 3 and attributes in column 9. Gene labels use `gene_name`, falling back to `gene_id`. Transcript plotting needs `transcript_id`. Recognized features are gene, transcript, exon, CDS and start_codon. Convert GFF3/mRNA annotations to GTF first.

Uncompressed bedGraph/GTF files get small temporary byte-range indexes. The first run builds them; subsequent runs reuse them. `SPARK_DISABLE_TEXT_INDEX=1` forces streaming. BigWig missing/non-finite values become zero; a missing chromosome is an error. bedGraph gaps are zero, overlapping rows use the last encountered value, and malformed relevant records are errors.

For compressed bedGraph, use HTSlib `bgzip` and `tabix -p bed file.bdg.gz`. A `.tbi` or `.csi` file requires an installed `tabix` executable; retrieval errors are reported. Without an index, gzip input is streamed.

## Python notebooks

```python
from SparK import plot, plot_many

result = plot(
    "chr1:10000-20000",
    tracks={"Control": ["c1.bw", "c2.bw"], "Drug": ["t1.bw", "t2.bw"]},
    show_plots="averages",
    output="comparison",
)
print(result["outputs"]["jpg"])

batch = plot_many(
    genes=["GAPDH", "IFFO1"],
    tracks={"Control": ["sample.bw"]},
    gtffile="genes.gtf",
    output="gene_figures",
)
```

The API returns output paths and raises normal Python exceptions for invalid single plots. Batch results include per-job errors. Calls have independent renderer state. The Python API does not start update checks or a browser. For an installed command, `python -m pip install .` provides `spark-plot`; running the standalone `SparK.py` continues to work.

## Local visual interface

`python SparK.py --gui --data-dir /path/to/data` opens a local browser interface. Browse existing files, choose a gene or coordinates, set conditions/colors, preview, and download exports/settings. Large genomic files are read in place. The "Try demo" button creates tiny synthetic inputs. The server listens only on loopback with an ephemeral access token; it uses no hosted service. Stop it with Ctrl+C. `--no-open` prints the URL without opening a browser.

## Update notices and tests

CLI runs check the latest published stable GitHub release in a separate background process. Cached successful checks are reused for 24 hours; offline/failed checks back off for 10 minutes. A newer version is announced on stderr after plotting finishes. The plotting process never waits for the worker, and piped streams are not held open. Updates are never installed automatically. Use `--no-update-check` or `SPARK_NO_UPDATE_CHECK=1` to disable checking. Only a public release request is sent; no input paths, regions or signal data are transmitted.

```bash
SPARK_NO_UPDATE_CHECK=1 python -m unittest discover -s tests -v
```

GitHub Actions tests Python 3.10/3.13 on Linux/macOS. The original 12 bedGraph and 12 BigWig commands are in [examples](../examples). Large reference inputs stay outside GitHub. See [CHANGELOG.md](../CHANGELOG.md) for version history.
