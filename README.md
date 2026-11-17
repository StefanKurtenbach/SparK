# SparK 3.2.0

Plot BigWig and bedGraph signal files as editable genomic figures.

**Every plot saves an SVG, a JPG, and a small run record automatically.** PDF and PNG exports are also available.

## Install

Use Python 3.10 or newer. Linux and macOS are tested.

```bash
git clone https://github.com/StefanKurtenbach/SparK.git
cd SparK
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

## Make a plot

Replace `sample.bw` with your BigWig or bedGraph file:

```bash
python SparK.py -cf sample.bw -pr chr1:10000-20000 -o figure
```

This creates `figure.svg`, `figure.jpg`, and `figure.run.json`. Coordinates are zero-based: `chr1:0-100` covers 100 bases.

Try the included tiny example:

```bash
python SparK.py --config examples/demo.json
```

## Common tasks

```bash
# Open the visual interface. Files stay on your computer.
python SparK.py --gui --data-dir /path/to/your/data

# Find a gene and include 10 kb on each side.
python SparK.py -cf sample.bw --gene GAPDH -gtf genes.gtf --flank 10000 -o GAPDH

# Plot every region in a BED file.
python SparK.py -cf sample.bw --regions-file regions.bed -o figures

# Add PDF and PNG copies; SVG and JPG are still saved.
python SparK.py -cf sample.bw -pr chr1:10000-20000 --formats pdf png -o figure
```

Save any working command's settings with `--save-config settings.json`. Run it again with `python SparK.py --config settings.json`. Use `--no-jpg` to skip the JPG copy.

SparK checks for a newer published release in the background and prints a notice after plotting. It never installs updates or waits for the network. Disable this with `--no-update-check`.

For replicates, multiple conditions, colors, differences/ratios, batch scales, or Python notebooks, see the [usage guide](docs/usage.md). Run `python SparK.py --help` for all options.

Please cite the [SparK preprint](https://www.biorxiv.org/content/10.1101/845529v1.full). Report problems in [GitHub Issues](https://github.com/StefanKurtenbach/SparK/issues).
