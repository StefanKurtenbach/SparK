# SparK 3.1.0

SparK plots genomic signal from **BigWig and bedGraph** files as editable SVG figures. It supports replicate overlays, replicate means, standard-deviation bands, mirrored sine plots, BED intervals, and GTF gene/transcript annotations.

Please cite the [SparK preprint](https://www.biorxiv.org/content/10.1101/845529v1.full). Open an [issue](https://github.com/StefanKurtenbach/SparK/issues) for questions, bugs, or feature requests.

<img src="Images/Summary.png" width="900" alt="Examples of SparK genomic signal plots">

## Install and run

Use Python 3.10 or newer. NumPy is required; pyBigWig is required only for `.bw` and `.bigWig` inputs. Installing the supplied requirements installs both.

```bash
git clone https://github.com/StefanKurtenbach/SparK.git
cd SparK
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python SparK.py -cf sample.bw -pr chr1:10000-20000 -o sample
python SparK.py --help
python SparK.py --version
```

On Windows, activate the environment with `.venv\Scripts\activate` instead. The example produces `sample.svg`. An output name that already ends in `.svg` is also accepted.

Regions use **zero-based, half-open coordinates**, matching bedGraph, BED, and BigWig: `chr1:0-100` contains 100 bases. Use a nonnegative start and an end greater than the start. GTF input coordinates are converted from their [one-based inclusive convention](https://www.ensembl.org/info/website/upload/gff.html). Use the same genome assembly for all inputs.

## Groups and replicates

With standard plots and `-ps all` (the default), each input gets its own track unless groups are supplied. `-ps averages` puts the inputs in one group by default. STD and sine plots also default to one shared group.

Supply **one group number for every file**, separately for controls and treatments. Missing or extra assignments are errors, so files cannot silently disappear. Group numbers must be positive; unused group numbers are skipped. Labels and custom scales are indexed by the original group numbers, including excluded groups.

```bash
python SparK.py \
  -cf control_H3K27ac_1.bw control_H3K27ac_2.bw control_H3K4me3_1.bw control_H3K4me3_2.bw \
  -tf treated_H3K27ac_1.bw treated_H3K27ac_2.bw treated_H3K4me3_1.bw treated_H3K4me3_2.bw \
  -cg 1 1 2 2 -tg 1 1 2 2 \
  -gl H3K27ac H3K4me3 -l Control Treated \
  -ps averages -gtf genes.gtf -dg GAPDH -pr chr12:6527512-6550512 -o comparison
```

| Plot | Option | Behavior |
| --- | --- | --- |
| Standard | `-pt standard` | Individual signals, or means with `-ps averages` |
| Standard deviation | `-pt STD` | Mean line and mean ± population SD; at least two replicates for each plotted condition |
| Sine | `-pt sine` | Controls below zero, treatments above zero, using one shared magnitude scale; requires both conditions in each group |

Sine plots honor `-ps averages`. Signed input values retain their sign; sine plots then mirror the control side. For a conventional mirror of two positive coverage signals, provide positive-valued inputs.

## Track colors

Hex colors accept three or six digits, with or without `#`. **Quote colors beginning with `#` in the shell**, otherwise the shell can treat them as comments.

```bash
# One color for every signal
python SparK.py -cf a.bw b.bw -pr chr1:0-1000 -f 333333

# Legacy two-color form: CONTROL first, TREATMENT second
python SparK.py -cf control.bw -tf treated.bw -pr chr1:0-1000 -f "#0066CC" "#DD5500"

# Four input files, four individual colors (issue #28)
python SparK.py -cf a.bw b.bw c.bw d.bw -pr chr1:0-1000 \
  -f 000000 FF0000 00AA00 0000FF

# Explicit per-input colors, including the otherwise ambiguous two-file case
python SparK.py -cf a.bw b.bw -pr chr1:0-1000 --track_colors 0066CC DD5500
```

More than two colors with `-f` means one color per input. `--track_colors` always means one color per input. In both forms, list **all control colors first, then all treatment colors**. Do not combine `-f` and `--track_colors`.

When replicates are combined into a mean or SD band, they must share a color within each condition/group. The `-l` two-condition legend uses the first control and treatment colors; use `-gl` to label separate tracks with individual colors.

Presets are `blue/red`, `blue/grey`, `all_grey`, and `blue/green`. Their historical colors and opacity are retained. Overlaid controls/treatments use transparency so both signals remain visible. The old README reversed the documented order of two custom colors; the code's control-then-treatment order is retained.

## Scaling, downsampling, and smoothing

Each group has its own y scale by default. `-gs yes` uses a shared scale across included groups, including negative-valued signals. `-es 2` lets group 2 keep its own scale. `-eg 2` removes group 2 from the plot and scaling calculations.

`-cs 10 D 100` sets the signal magnitude used to normalize groups 1 and 3, retaining automatic scaling for group 2. Supply one entry per group number. Values must be positive and finite. Signals occupy 90% of the track height; standard-axis tick labels account for this margin. Sine ticks are placed directly at the chosen positive/negative magnitude.

Long regions are reduced to `--max_points 2000` bins by default:

- `--point_stat average` computes the overlap-weighted average within each bin.
- `--point_stat max` preserves the maximum value overlapping each bin, useful for narrow positive peaks. For signed signals, this is the numeric maximum, not the largest absolute magnitude.
- `--point_stat` changes **spatial downsampling**, not the statistic used to combine replicates.
- `-sm 10` applies a trailing 10-bin mean after downsampling, padding the beginning with the first value. It does not add bins to short regions.

## Other options

| Option | Purpose |
| --- | --- |
| `-cf`, `-tf` | Control and treatment files, separated by spaces; controls and `-pr` are required |
| `-cg`, `-tg` | One group assignment per corresponding input file |
| `-gl` | One label per group number |
| `-l Control Treated` | Two-condition legend |
| `-dc` | `top_left`, `top_right`, `bottom_left`, `bottom_right`, or `no` for the region label |
| `-w 150` | Track width in SVG units |
| `--x_scale`, `--y_scale` | Positive horizontal/vertical scale factors; default 1 |
| `-scale no` | Hide the genomic scale bar |
| `-bed regions.bed` | One or more BED annotation tracks |
| `-bedcol`, `-bedlab` | One color or one per BED file; one label per BED file |
| `-gtf genes.gtf` | GTF annotations, optionally gzip-compressed |
| `-dg GAPDH IFFO1` | Gene names or fallback gene IDs to display |
| `-dt all` | Display individual transcripts rather than merged genes |
| `-dt ID1 ID2` | Display selected transcript IDs |
| `-wg yes` | Label transcripts with their gene names |
| `-tss no` | Hide translation start-codon arrows; this legacy option name does not mean transcription start |
| `-genestart yes` | Add arrows at annotated gene/transcript starts (strand-aware TSS) |
| `-sp yes` | Highlight regions where the difference of means exceeds the sum of SDs; requires at least two replicates per condition in each group |
| `-sc FF0000 00FF00` | Highlight colors: treatment-up first, control-up second |

“Sparks” are a descriptive mean/SD rule, **not a statistical significance test or p-value**. Choose appropriate biological replicates, normalization, and statistical methods for your analysis.

## Input formats and performance

Read BigWig files directly; conversion to bedGraph is unnecessary. For BAM input, prepare coverage with an appropriate upstream tool, such as [deepTools bamCoverage](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html), choosing normalization and bin size for the experiment. SparK plots supplied signal and does not normalize it.

bedGraph requires tab-separated chromosome, start, end, and finite numeric signal. Empty lines, comments, and standard track/browser headers are ignored. Overlapping records use the last value encountered. Missing bases are zero. For BigWig, missing/non-finite values are zero; a missing chromosome is an error. Chromosome prefixes `chr` and `Chr` are recognized. Input intervals are clipped before assigning signal, preventing values from wrapping around a selected region.

Uncompressed bedGraph and GTF files get small byte-range indexes in the system temporary directory. Data files are not copied or modified. The first run builds the index; repeat runs reuse it. File metadata changes invalidate it. Set `SPARK_DISABLE_TEXT_INDEX=1` to use streaming for troubleshooting.

For compressed bedGraph, install [HTSlib bgzip/tabix](https://www.htslib.org/doc/tabix.html), then:

```bash
bgzip -c sorted.bdg > sorted.bdg.gz
tabix -p bed sorted.bdg.gz
python SparK.py -cf sorted.bdg.gz -pr chr1:10000-20000
```

Gzip bedGraph without an index is streamed. An installed `tabix` executable is required when a `.tbi` index is used; retrieval failures are reported instead of silently producing an empty signal.

GTF requires nine tab-separated columns, with feature type in column 3 and attributes in column 9. Recognized features are `gene`, `transcript`, `exon`, `CDS`, and `start_codon`. Gene labels use `gene_name`, falling back to `gene_id`. Transcript plots require `transcript_id`. GFF3 attributes and `mRNA` feature records need conversion to GTF before use.

## Tests and example commands

```bash
python -m unittest discover -s tests -v
```

The automated suite uses tiny generated inputs and checks numeric downsampling, reader boundaries, input/group validation, individual colors, averaging, signed scales, annotations, complete SVG output, and highlight direction. GitHub Actions runs it for pushes and pull requests. It does not download large genomic datasets.

The [bedGraph example commands](examples/spark_test_command.sh) and [BigWig example commands](examples/spark_test_command_bw_only.sh) each generate 12 plots using the original ENCODE example filenames. Run them from a directory containing `SparK.py` and those inputs. They run sequentially and stop on errors. The sine commands now include all four replicates in each condition (`-cg 1 1 2 2 -tg 1 1 2 2`).

The original figures below illustrate the plot styles; corrected input boundaries, y labels, and annotations can differ from these historical examples.

| Individual tracks | Replicate means | Standard deviation |
| --- | --- | --- |
| ![Individual signals](Images/FigureA.png) | ![Means](Images/FigureB.png) | ![SD bands](Images/FigureF.png) |

See [CHANGELOG.md](CHANGELOG.md) for behavior changes and migration notes.
