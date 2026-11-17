# Changes

## 3.2.0

- Save JPG beside SVG automatically. Add optional PNG and vector PDF exports, with staged file replacement and configurable raster DPI/quality.
- Check the latest published stable GitHub release in a separate background process. Announce newer versions after plotting, without waiting, installing, or failing offline; add opt-out flags.
- Add reusable JSON settings, explicit command-line overrides, and multiple named conditions with replicates, colors, panel groups and per-input scale factors.
- Add gene-name/ID lookup, coordinate/BED/gene batches, shared batch scales, HTML previews, and per-job error reporting.
- Add mean-based difference, ratio and log2-ratio tracks with an explicit comparison pair and pseudocount. These are descriptive transformations, not significance tests.
- Save reproducible run records with code/output hashes, settings, coordinates and input metadata; optionally hash complete inputs.
- Isolate renderer state and expose `plot`/`plot_many` for Python notebooks. Add an optional local visual interface and tiny synthetic examples.
- Simplify the README and move advanced instructions into `docs/usage.md`. Add installation metadata and the `spark-plot` command.
- Reject figure/configuration destinations that would overwrite genomic inputs or the running script.
- Preserve all 24 established 3.1.0 reference SVGs byte-for-byte. Large genomic inputs remain in place.

### Upgrade

Install the updated dependencies with `python -m pip install -r requirements.txt`. Python 3.10 or newer is required. Existing plotting commands still work; they now also save `.jpg` and `.run.json`. Use `--no-jpg` or `--no-run-record` to skip those files. The update checker only reports published releases and never changes installed code.

## 3.1.0

- Add per-input colors for `-f` lists longer than two, and `--track_colors` for any input count. Validate hex colors and color counts. One/two-color `-f` behavior remains control/treatment based. Fixes #28.
- Include the per-bin `--point_stat max` and shared sine-axis changes from #27. Spatial reduction no longer changes replicate averaging.
- Fix weighted-average downsampling dropping the final fractional bin. Use integer bin boundaries and exact overlap weights.
- Clip bedGraph intervals before assigning values, preventing negative-index wraparound at the left boundary. Correct offset direction and support non-contiguous chromosome records.
- Preserve the actual first and last signal samples when closing SVG polygons.
- Honor `-ps averages` for sine plots, support treatment-only standard groups and control-only SD plots, and validate replicate requirements before plotting.
- Reject missing/extra group assignments instead of silently omitting files. Correct all eight group assignments in the sine examples.
- Skip excluded/empty groups safely; keep remaining tracks, labels, and annotation rows aligned. Handle shared scaling with negative signals and with all groups excluded from shared scaling.
- Fix short-region smoothing padding the result to 2,000 points. Match output length to the sampled signal.
- Correct standard y-axis labels for the 90% track-height margin. Include SD envelopes in automatic scales.
- Split highlights at nonsignificant gaps, use consistent control-up/treatment-up colors, and use the shared scale for sine highlights.
- Convert GTF one-based inclusive coordinates to the plot's zero-based half-open convention. Add gene-ID fallback, XML-escape labels, retain annotation row placement, and implement the existing `-genestart` option.
- Anchor the top chromosome label to the first track; expand the SVG canvas for wide tracks and long annotation stacks.
- Preserve buffered SVG writing and small cached text indexes from the Dropbox 3.0.1 version. Reuse loaded signal for shared-scale plots.
- Write output atomically. Invalid arguments and input errors return a nonzero exit status without replacing a previous figure.
- Update installation, colors, coordinates, plotting, test commands, and the distinction between descriptive highlights and statistical significance in the README. Add a generated-fixture test suite and GitHub Actions workflow.

### Migration notes

Previously saved figures are retained as historical references. Corrected reader boundaries, endpoint samples, binning, axes, GTF coordinates, and sine replicate assignments can intentionally change output. Review those differences before replacing an analysis baseline. Inputs are never copied into the repository.

Two custom `-f` colors remain **control first, treatment second**; older documentation incorrectly reversed this order. Use `--track_colors` for two independently colored files in the same condition. Colors beginning with `#` must be quoted in shell commands.
