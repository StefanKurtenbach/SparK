SparK_Version = "3.1.0"

import numpy as np
import copy
import os
import argparse
import math
import sys
import subprocess
import gzip
import atexit
import contextlib
import hashlib
import heapq
import json
import mmap
import re
import tempfile
from html import escape
try:
    import pyBigWig
except ImportError:
    pyBigWig = None

def get_gene_name():
    attributes = dict(re.findall(r'(?:^|;)\s*(\S+)\s+"([^"]*)"', line_split[8]))
    name = attributes.get('gene_name') or attributes.get('gene_id')
    if not name:
        raise ValueError('GTF annotation requires gene_name or gene_id in column 9')
    return name
def get_transcript_name():
    match = re.search(r'(?:^|;)\s*transcript_id\s+"([^"]+)"', line_split[8])
    if not match:
        raise ValueError('Transcript plotting requires transcript_id in GTF column 9')
    return match.group(1)
def draw_line(coordinates, thickness, color):
    output = '''<path d="'''
    for x, i in enumerate(coordinates):
        if x == 0:
            output += "M "
        else:
            output += " L "
        output += str(i[1]) + " " + str(y_start - i[0])
    output += '''" style="stroke:'''+ color + '''; stroke-width:'''+ str(thickness) +'''; fill:none "/>'''
    return(output)

def get_lines_tabix(filename, stretch):
    """Tabix query strings are one-based; input regions are zero-based half-open."""
    chromname = stretch[3] if len(stretch) > 3 else stretch[0]
    query = '{}:{}-{}'.format(chromname, stretch[1] + 1, stretch[2])
    result = subprocess.run(['tabix', filename, query], capture_output=True, text=True, check=True)
    yield from result.stdout.splitlines()

def get_lines_stream(filename):
    """
    helper function for make_raw_data_filled().
    If a data file is not tabix-indexed, go through it line by line.
    """
    if filename[-3:] == '.gz':
        f = gzip.open(filename, 'rt')
    else:
        f = open(filename, 'r')
    try:
        for line in f:
            yield line.rstrip()
    finally:
        f.close()

def _chromosome_offsets(mapped, chromosome):
    """Find complete chromosome fields in C, without decoding unrelated rows."""
    def positions(prefix):
        if mapped[:len(prefix)] == prefix:
            yield 0
        needle = b'\n' + prefix
        position = mapped.find(needle, 0)
        while position >= 0:
            yield position + 1
            position = mapped.find(needle, position + len(needle))
    variants = {chromosome, 'chr' + chromosome, 'Chr' + chromosome}
    streams = [positions(name.encode('utf-8') + b'\t') for name in sorted(variants)]
    return heapq.merge(*streams)


def _text_index(filename, chromosome, coordinate_columns, first_block_only):
    """Cache byte ranges, not data; stat changes invalidate the tiny index."""
    realpath = os.path.realpath(filename)
    stat = os.stat(realpath)
    signature = [realpath, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns,
                 stat.st_ino, chromosome, coordinate_columns, first_block_only, 1]
    key = hashlib.sha256(json.dumps(signature).encode('utf-8')).hexdigest()
    cache_dir = os.path.join(tempfile.gettempdir(),
                             'spark-text-index-v1-' + str(getattr(os, 'getuid', lambda: 0)()))
    cache_file = os.path.join(cache_dir, key + '.json')
    try:
        with open(cache_file) as handle:
            cached = json.load(handle)
        if cached['signature'] == signature:
            return cached['blocks']
    except (OSError, ValueError, KeyError, TypeError):
        pass
    blocks = []
    previous_end = None
    count = 0
    if stat.st_size:
        with open(realpath, 'rb') as source:
            with mmap.mmap(source.fileno(), 0, access=mmap.ACCESS_READ) as mapped:
                if re.search(rb'\r(?!\n)', mapped):
                    # Universal-newline text input also supports old-style
                    # CR records; retain that behavior through streaming.
                    raise ValueError('Standalone carriage-return record separators')
                for position in _chromosome_offsets(mapped, chromosome):
                    # The original bedGraph reader stops at the end of its
                    # first contiguous chromosome block. GTF may be unsorted.
                    if first_block_only and previous_end is not None and position != previous_end:
                        break
                    mapped.seek(position)
                    row = mapped.readline()
                    end = mapped.tell()
                    if not blocks or position != previous_end or (coordinate_columns and count >= 4096):
                        blocks.append([position, end, None, None, False])
                        count = 0
                    block = blocks[-1]
                    block[1] = end
                    if coordinate_columns:
                        try:
                            fields = row.split(b'\t')
                            start = int(fields[coordinate_columns[0]])
                            stop = int(fields[coordinate_columns[1]])
                            float(fields[3])
                            block[2] = start if block[2] is None else min(block[2], start)
                            block[3] = stop if block[3] is None else max(block[3], stop)
                        except (ValueError, IndexError):
                            # Never skip malformed records: let the existing
                            # parser report them as it did before indexing.
                            block[4] = True
                    count += 1
                    previous_end = end
    after = os.stat(realpath)
    if (after.st_size, after.st_mtime_ns, after.st_ctime_ns, after.st_ino) != (
            stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns, stat.st_ino):
        raise OSError('Input changed while its text index was being built')
    temporary = None
    try:
        os.makedirs(cache_dir, mode=0o700, exist_ok=True)
        with tempfile.NamedTemporaryFile(mode='w', dir=cache_dir, prefix=key + '.',
                                         suffix='.tmp', delete=False) as handle:
            temporary = handle.name
            json.dump({'signature': signature, 'blocks': blocks}, handle)
        os.replace(temporary, cache_file)
    except OSError:
        # Read-only/full temporary storage must not prevent a plot.
        if temporary is not None:
            try:
                os.unlink(temporary)
            except OSError:
                pass
    return blocks


def get_lines_indexed(filename, chromosome, lower=None, upper=None,
                      coordinate_columns=None, first_block_only=False):
    """Read the relevant plain-text ranges; preserve the streaming fallback."""
    if filename.endswith('.gz') or os.environ.get('SPARK_DISABLE_TEXT_INDEX') == '1':
        if filename.endswith('.gz'):
            with gzip.open(filename, 'rt') as source:
                yield from source
        else:
            with open(filename) as source:
                yield from source
        return
    try:
        # Lists rather than tuples keep the JSON signature stable on reload.
        columns = list(coordinate_columns) if coordinate_columns else None
        blocks = _text_index(filename, chromosome, columns, first_block_only)
    except (OSError, ValueError, TypeError, KeyError):
        with open(filename) as source:
            yield from source
        return
    with open(filename, 'rb') as source:
        for start, stop, min_start, max_stop, unsafe in blocks:
            if coordinate_columns and not unsafe and min_start is not None:
                if max_stop < lower or min_start > upper:
                    continue
            source.seek(start)
            while source.tell() < stop:
                # Match TextIOWrapper's universal-newline behavior.
                yield source.readline().decode('utf-8').replace('\r\n', '\n').replace('\r', '\n')


def make_raw_data_filled(stretch, files, offset):
    """Read zero-based half-open intervals, clipping before assigning signal."""
    lower, upper = stretch[1] + offset, stretch[2] + offset
    length = upper - lower
    raw_data_filled = []
    for filename in files:
        data = [0.0] * length
        if filename.lower().endswith(('.bw', '.bigwig')):
            if pyBigWig is None:
                raise ImportError('Install pyBigWig to read BigWig files: python -m pip install pyBigWig')
            with contextlib.closing(pyBigWig.open(filename)) as bw:
                chroms = bw.chroms()
                chromname = next((c for c in (stretch[0], 'chr' + stretch[0], 'Chr' + stretch[0]) if c in chroms), None)
                if chromname is None:
                    raise ValueError('Chromosome {} is absent from {}'.format(stretch[0], filename))
                start, stop = max(0, lower), min(upper, chroms[chromname])
                if start < stop:
                    values = bw.values(chromname, start, stop)
                    data[start-lower:stop-lower] = [float(v) if v is not None and math.isfinite(v) else 0.0 for v in values]
        else:
            query = [stretch[0], lower, upper] + stretch[3:]
            if os.path.isfile(filename + '.tbi') or os.path.isfile(filename + '.csi'):
                rows = get_lines_tabix(filename, query)
            else:
                rows = get_lines_indexed(filename, stretch[0], lower, upper, (1, 2))
            with contextlib.closing(rows):
                for line in rows:
                    if not line.strip() or line.startswith(('#', 'track ', 'browser ')):
                        continue
                    fields = line.rstrip().split('\t')
                    chrom = fields[0][3:] if fields[0].startswith(('chr', 'Chr')) else fields[0]
                    if chrom != stretch[0]:
                        continue
                    try:
                        start, stop, value = int(fields[1]), int(fields[2]), float(fields[3])
                    except (ValueError, IndexError) as exc:
                        raise ValueError('Invalid bedGraph record in {}: {}'.format(filename, line.rstrip())) from exc
                    if start < 0 or stop <= start or not math.isfinite(value):
                        raise ValueError('Invalid bedGraph interval or value in {}: {}'.format(filename, line.rstrip()))
                    start, stop = max(start, lower), min(stop, upper)
                    if start < stop:
                        data[start-lower:stop-lower] = [value] * (stop-start)
        if length > max_points:
            data = downsample(data, max_points, point_stat)
        if smoothen_tracks is not None:
            # A trailing window padded with the first value, with no extra bins.
            prefix = np.concatenate(([0.0], np.cumsum(data, dtype=float)))
            window = smoothen_tracks
            first = data[0]
            data = [(prefix[i+1] - prefix[max(0, i-window+1)] + max(0, window-i-1)*first) / window
                    for i in range(len(data))]
        raw_data_filled.append(data)
    return raw_data_filled


def downsample(data, bins, statistic='average'):
    """Exact integer bin boundaries avoid dropping fractional final bins."""
    size = len(data)
    if bins <= 0 or statistic not in ('average', 'max'):
        raise ValueError('Downsampling requires positive bins and average or max')
    if size <= bins:
        return list(data)
    output = []
    for position in range(bins):
        left, right = position * size, (position + 1) * size
        start, stop = left // bins, (right + bins - 1) // bins
        if statistic == 'max':
            output.append(max(data[start:stop]))
        else:
            last = stop - 1
            first_weight = (start + 1) * bins - left
            last_weight = right - last * bins
            value = (data[start] * first_weight + math.fsum(data[start+1:last]) * bins + data[last] * last_weight) / size
            output.append(value)
    return output
_svg_output = None
_svg_temporary = None

def write_to_file(row):
    global _svg_output, _svg_temporary
    if _svg_output is None:
        _svg_output = tempfile.NamedTemporaryFile(mode='w', encoding='utf-8',
            dir=os.path.dirname(os.path.abspath(output_filename)), prefix='.spark-', suffix='.svg', delete=False)
        _svg_temporary = _svg_output.name
        atexit.register(discard_incomplete_svg)
    _svg_output.write(row + "\n")


def discard_incomplete_svg():
    if _svg_output is not None and not _svg_output.closed:
        _svg_output.close()
    if _svg_temporary and os.path.exists(_svg_temporary):
        os.unlink(_svg_temporary)


def finish_svg(width, height):
    """Publish only a complete SVG, retaining the previous output on failure."""
    _svg_output.close()
    with open(_svg_temporary, encoding='utf-8') as handle:
        content = handle.read()
    content = re.sub(r'viewBox="[^"]+"', 'viewBox="0 0 {} {}"'.format(width, height), content, count=1)
    with open(_svg_temporary, 'w', encoding='utf-8') as handle:
        handle.write(content)
    os.replace(_svg_temporary, output_filename)
def get_max_value(datasets1, datasets2):
    maxima = [0.0]
    for datasets in (datasets1, datasets2):
        if not datasets:
            continue
        values = np.asarray(datasets, dtype=float)
        if plot_type == 'STD':
            mean, std = np.mean(values, axis=0), np.std(values, axis=0)
            maxima.extend([float(np.max(np.abs(mean + std))), float(np.max(np.abs(mean - std)))])
        elif show_plots == 'averages' and plot_type != 'sine':
            maxima.append(float(np.max(np.abs(np.mean(values, axis=0)))))
        else:
            maxima.append(float(np.max(np.abs(values))))
    return max(maxima)
def get_relative_hight(raw_value): # FIX make sure maxvalue can be 0 too
    if raw_value == 0:
        return(0)
    else:
        return((raw_value * hight * relative_track_hight_percentage) / max_value) # to not go up to the max
def get_relative_hight_custom(raw_value, denom):
    if raw_value == 0 or denom == 0:
        return 0
    return (raw_value * hight * relative_track_hight_percentage) / denom
def draw_rect(x_coord, y_0, color, width, hight1, opacity):
    return '''<rect x="''' + str(x_coord) + '''" opacity="''' + str(opacity) + '''" y="''' + str(y_0 - hight1) + '''" fill="''' + color + '''" width="''' + str(width) + '''" height="''' + str(hight1) + '''"/>'''
def draw_polygon(coordinates, opacity, color, stroke_width):
    for q, w in enumerate(coordinates):
        coordinates[q][0] = y_start - coordinates[q][0]
    string = '''<polygon points="'''
    for h, c in enumerate(coordinates):
        if h == 0:
            string += str(c[1]) + "," + str(c[0])
        else:
            string += " " + str(c[1]) + "," + str(c[0])
    string += '''" opacity="''' + str(opacity) + '''" fill="''' + color + '''"''' + ''' stroke="black" stroke-width="''' + str(stroke_width) + '''"/>'''
    return string
def draw_axis_for_group(y_start_val, max_value_val, has_negative_axis):
    axis_label = round(max_value_val / relative_track_hight_percentage, 1)
    if has_negative_axis:
        write_to_file('''<line x1="''' + str(x_start - 10) + '''" y1="''' + str(y_start_val + hight) + '''" x2="''' + str(x_start - 10) + '''" y2="''' + str(y_start_val - hight) + '''" stroke="black" stroke-width="1" />''')
        write_to_file('''<line x1="''' + str(x_start - 10.5) + '''" y1="''' + str(y_start_val) + '''" x2="''' + str(x_start - 6.5) + '''" y2="''' + str(y_start_val) + '''" stroke="black" stroke-width="1" />''')
        write_to_file('''<line x1="''' + str(x_start - 10.5) + '''" y1="''' + str(y_start_val - hight) + '''" x2="''' + str(x_start - 6.5) + '''" y2="''' + str(y_start_val - hight) + '''" stroke="black" stroke-width="1" />''')
        write_to_file('''<line x1="''' + str(x_start - 10.5) + '''" y1="''' + str(y_start_val + hight) + '''" x2="''' + str(x_start - 6.5) + '''" y2="''' + str(y_start_val + hight) + '''" stroke="black" stroke-width="1" />''')

        write_to_file('''<text text-anchor="end" font-family="Arial" x="''' + str(x_start - 14) + '''" y="''' + str(y_start_val + hight + 4) + '''" font-size="''' + str(font_size_axis_y) + '''" >''' + str(-axis_label) + '''</text>''')
        write_to_file('''<text text-anchor="end" font-family="Arial" x="''' + str(x_start - 14) + '''" y="''' + str(y_start_val + 4) + '''" font-size="''' + str(font_size_axis_y) + '''" >0</text>''')
        write_to_file('''<text text-anchor="end" font-family="Arial" x="''' + str(x_start - 14) + '''" y="''' + str(y_start_val - hight + 4) + '''" font-size="''' + str(font_size_axis_y) + '''" >''' + str(axis_label) + '''</text>''')
    else:
        write_to_file('''<line x1="''' + str(x_start - 10) + '''" y1="''' + str(y_start_val) + '''" x2="''' + str(x_start - 10) + '''" y2="''' + str(y_start_val - hight) + '''" stroke="black" stroke-width="1" />''')
        write_to_file('''<line x1="''' + str(x_start - 10.5) + '''" y1="''' + str(y_start_val) + '''" x2="''' + str(x_start - 6.5) + '''" y2="''' + str(y_start_val) + '''" stroke="black" stroke-width="1" />''')
        write_to_file('''<line x1="''' + str(x_start - 10.5) + '''" y1="''' + str(y_start_val - hight) + '''" x2="''' + str(x_start - 6.5) + '''" y2="''' + str(y_start_val - hight) + '''" stroke="black" stroke-width="1" />''')

        write_to_file('''<text text-anchor="end" font-family="Arial" x="''' + str(x_start - 14) + '''" y="''' + str(y_start_val + 4) + '''" font-size="''' + str(font_size_axis_y) + '''" >0</text>''')
        write_to_file('''<text text-anchor="end" font-family="Arial" x="''' + str(x_start - 14) + '''" y="''' + str(y_start_val - hight + 4) + '''" font-size="''' + str(font_size_axis_y) + '''" >''' + str(axis_label) + '''</text>''')
def draw_sine_axis_for_group(y_start_val, shared_max_val):
    # Shared axis with symmetric labels for positive/negative magnitudes
    if shared_max_val >= 0:
        axis_height = hight * relative_track_hight_percentage
        pos_tick_y = y_start_val - axis_height
        neg_tick_y = y_start_val + axis_height
        # Draw full shared axis
        write_to_file('''<line x1="''' + str(x_start - 10) + '''" y1="''' + str(neg_tick_y) + '''" x2="''' + str(x_start - 10) + '''" y2="''' + str(pos_tick_y) + '''" stroke="black" stroke-width="1" />''')
        # 0 tick
        write_to_file('''<line x1="''' + str(x_start - 10.5) + '''" y1="''' + str(y_start_val) + '''" x2="''' + str(x_start - 6.5) + '''" y2="''' + str(y_start_val) + '''" stroke="black" stroke-width="1" />''')
        # Positive tick and label
        write_to_file('''<line x1="''' + str(x_start - 10.5) + '''" y1="''' + str(pos_tick_y) + '''" x2="''' + str(x_start - 6.5) + '''" y2="''' + str(pos_tick_y) + '''" stroke="black" stroke-width="1" />''')
        write_to_file('''<text text-anchor="end" font-family="Arial" x="''' + str(x_start - 14) + '''" y="''' + str(pos_tick_y + 4) + '''" font-size="''' + str(font_size_axis_y) + '''" >''' + str(round(shared_max_val, 1)) + '''</text>''')
        # Negative tick and label
        write_to_file('''<line x1="''' + str(x_start - 10.5) + '''" y1="''' + str(neg_tick_y) + '''" x2="''' + str(x_start - 6.5) + '''" y2="''' + str(neg_tick_y) + '''" stroke="black" stroke-width="1" />''')
        write_to_file('''<text text-anchor="end" font-family="Arial" x="''' + str(x_start - 14) + '''" y="''' + str(neg_tick_y + 4) + '''" font-size="''' + str(font_size_axis_y) + '''" >-''' + str(round(shared_max_val, 1)) + '''</text>''')
def draw_standard_spark(sine=False, shared_max=None):
    if len(control_data) < 2 or len(treat_data) < 2:
        raise ValueError('Sparks require at least two control and two treatment files in every plotted group')
    coordinates, direction = [], None

    def emit():
        if not coordinates:
            return
        color = spark_color[1] if direction == 'control' else spark_color[0]
        points = list(coordinates)
        if sine:
            points = [[0, points[0][1]]] + points + [[points[-1][0], points[-1][1] + quantile], [0, points[-1][1] + quantile]]
        else:
            points.extend([[points[-1][0], points[-1][1] + quantile], [points[0][0], points[0][1] + quantile]])
        write_to_file(draw_polygon(points, spark_opacity, color, stroke_width_spark))

    for x, (ctrl, treat) in enumerate(zip(np.asarray(control_data).T, np.asarray(treat_data).T)):
        ctrl_mean, treat_mean = float(np.mean(ctrl)), float(np.mean(treat))
        ctrl_std, treat_std = float(np.std(ctrl)), float(np.std(treat))
        delta = ctrl_mean - treat_mean
        current = ('control' if delta > 0 else 'treatment') if abs(delta) > ctrl_std + treat_std else None
        if current != direction or current is None:
            emit()
            coordinates, direction = [], current
        if current is None:
            continue
        x_pos = x_start + x * quantile
        if sine:
            value = abs(delta) - ctrl_std - treat_std
            value = -value if current == 'control' else value
            coordinates.append([get_relative_hight_custom(value, shared_max), x_pos])
        else:
            low, high = ((treat_mean + treat_std, ctrl_mean - ctrl_std) if current == 'control'
                         else (ctrl_mean + ctrl_std, treat_mean - treat_std))
            coordinates.append([get_relative_hight(high), x_pos])
            coordinates.insert(0, [get_relative_hight(low), x_pos])
    emit()
def get_region_to_draw():
    region_to_draw = 0
    if line_split[4] > region[1] and line_split[3] < region[2]: # check if there is something to draw
        region_to_draw = [float(line_split[3]), float(line_split[4])]
        if line_split[3] < region[1]:
            region_to_draw[0] = float(region[1])
        if line_split[4] > region[2]:
            region_to_draw[1] = float(region[2])
    return(region_to_draw)


def hex_color(value):
    value = value.lstrip('#')
    if not re.fullmatch(r'(?:[0-9a-fA-F]{3}|[0-9a-fA-F]{6})', value):
        raise ValueError('Invalid hex color {!r}; use three or six hex digits, optionally prefixed by #'.format(value))
    return '#' + value


def parse_region(value):
    match = re.fullmatch(r'([^:\s]+):(\d+)-(\d+)', value)
    if not match:
        raise ValueError('Region must use chr1:10000-20000 format (zero-based, half-open)')
    chrom, start, stop = match.group(1), int(match.group(2)), int(match.group(3))
    if stop <= start:
        raise ValueError('Region end must be greater than start')
    normalized = chrom[3:] if chrom.startswith(('chr', 'Chr')) else chrom
    if not normalized:
        raise ValueError('Chromosome name is empty')
    return [normalized, start, stop, chrom]


def configure_colors(args):
    nctrl, ntreat = len(args['control_files']), len(args['treat_files'])
    values, explicit = args['fills'], args['track_colors']
    if values and explicit:
        raise ValueError('Use either -f or --track_colors, not both')
    presets = {'blue/red': (['#FF1800', '#005CFF'], 0.6),
               'blue/grey': (['#005CFF', '#A3A3A3'], 0.3),
               'all_grey': (['#848484', '#848484'], 0.6),
               'blue/green': (['#00FF12', '#005CFF'], 0.5)}
    individual = explicit is not None or (values is not None and len(values) > 2)
    if individual:
        colors = [hex_color(v) for v in (explicit or values)]
        if len(colors) != nctrl + ntreat:
            raise ValueError('Per-track colors require one color per input: controls first, then treatments')
        control, treatment = colors[:nctrl], colors[nctrl:]
        fills = [control[0], treatment[0] if treatment else control[0]]
        opacity = 1 if not ntreat else 0.6
    else:
        if values is None:
            fills, opacity = (['#FF1800', '#005CFF'], 0.6) if ntreat else (['#0000C1', '#0000C1'], 1)
        elif len(values) == 1 and values[0] in presets:
            fills, opacity = presets[values[0]]
        else:
            fills = [hex_color(v) for v in values]
            opacity = 1 if len(fills) == 1 else 0.6
            if len(fills) == 1:
                fills *= 2
        control, treatment = [fills[0]] * nctrl, [fills[1]] * ntreat
    return fills, opacity, control, treatment, individual


def validate_arguments(args):
    for key in ('max_points', 'track_width'):
        if args[key] <= 0:
            raise ValueError('{} must be positive'.format(key))
    if args['smoothen'] is not None and args['smoothen'] <= 0:
        raise ValueError('smoothen must be a positive integer')
    for key in ('x_scale', 'y_scale'):
        if not math.isfinite(args[key]) or args[key] <= 0:
            raise ValueError('{} must be finite and positive'.format(key))
    if args['track_width'] * args['x_scale'] < 1:
        raise ValueError('track_width multiplied by x_scale must be at least one SVG unit')
    region = parse_region(args['region'])
    controls, treatments = args['control_files'], args['treat_files']
    cg, tg = args['control_groups'], args['treat_groups']
    if cg or tg:
        if len(cg) != len(controls) or len(tg) != len(treatments):
            raise ValueError('Provide exactly one group assignment per file in -cg and -tg; no inputs may be omitted')
    elif args['show_plots'] == 'averages' or args['plot_type'] in ('STD', 'sine'):
        cg, tg = [1] * len(controls), [1] * len(treatments)
    else:
        cg = list(range(1, len(controls) + 1))
        tg = list(range(len(controls) + 1, len(controls) + len(treatments) + 1))
    if any(group <= 0 for group in cg + tg):
        raise ValueError('Group numbers must be positive integers')
    args['control_groups'], args['treat_groups'] = cg, tg
    groups = set(cg + tg)
    excluded = set(args['exclude_groups'] or [])
    excluded_scale = set(args['exclude_from_group_autoscale'] or [])
    if not excluded.issubset(groups) or not excluded_scale.issubset(groups):
        raise ValueError('Excluded group numbers must refer to an input group')
    active = sorted(groups - excluded)
    if not active:
        raise ValueError('At least one group must remain after -eg exclusions')
    nr = max(groups)
    for key in ('group_labels', 'custom_y_axis_scales'):
        if args[key] is not None and len(args[key]) != nr:
            raise ValueError('{} requires one entry for each group number from 1 to {}'.format(key, nr))
    if args['custom_y_axis_scales'] is not None:
        for value in args['custom_y_axis_scales']:
            if value != 'D' and (not math.isfinite(float(value)) or float(value) <= 0):
                raise ValueError('Custom scales must be finite and positive, or D for automatic')
    if args['labels'] is not None and len(args['labels']) != 2:
        raise ValueError('-l requires two labels: control then treatment')
    for group in active:
        nc, nt = cg.count(group), tg.count(group)
        if args['plot_type'] == 'sine' and (not nc or not nt):
            raise ValueError('Sine plots require control and treatment inputs in every plotted group')
        if args['plot_type'] == 'STD' and any(n == 1 for n in (nc, nt)):
            raise ValueError('STD plots require at least two replicates for each plotted condition')
        if args['spark'] == 'yes' and (nc < 2 or nt < 2):
            raise ValueError('Sparks require at least two control and two treatment files in every plotted group')
    bed_files = args['bed_files'] or []
    if args['bed_labels'] is not None and len(args['bed_labels']) != len(bed_files):
        raise ValueError('-bedlab requires one label per BED file')
    if args['bed_color'] is not None:
        if len(args['bed_color']) not in (1, len(bed_files)) or not bed_files:
            raise ValueError('-bedcol requires one color or one color per BED file')
        args['bed_color'] = [hex_color(v) for v in args['bed_color']]
        if len(args['bed_color']) == 1:
            args['bed_color'] *= len(bed_files)
    if args['spark_color'] is not None:
        if len(args['spark_color']) != 2:
            raise ValueError('-sc requires two colors: treatment-up then control-up')
        args['spark_color'] = [hex_color(v) for v in args['spark_color']]
    colors = configure_colors(args)
    if colors[-1] and (args['show_plots'] == 'averages' or args['plot_type'] == 'STD'):
        for ids, palette in ((cg, colors[2]), (tg, colors[3])):
            for group in active:
                if len({color for g, color in zip(ids, palette) if g == group}) > 1:
                    raise ValueError('Replicates averaged into one trace must share a color; use -f with one/two colors')
    paths = controls + treatments + bed_files + ([args['gtffile']] if args['gtffile'] else [])
    output = args['output_name'] or 'SparK_graph'
    output = output if output.lower().endswith('.svg') else output + '.svg'
    for path in paths:
        if not os.path.isfile(path):
            raise ValueError('Input file does not exist: {}'.format(path))
        if os.path.realpath(path) == os.path.realpath(output):
            raise ValueError('Output must not overwrite an input file')
    return region, active, colors, output


def draw_signal(data, color, denominator=None, sign=1):
    denominator = max_value if denominator is None else denominator
    points = [[get_relative_hight_custom(value * sign, denominator), x_start + x * quantile]
              for x, value in enumerate(data)]
    # Close at the baseline without erasing the actual first or last sample.
    points = [[0, x_start]] + points + [[points[-1][0], x_start + total_width], [0, x_start + total_width]]
    write_to_file(draw_polygon(points, opacity, color, stroke_width))


def draw_start_arrow(position, strand, row_y):
    global y_start
    direction = 1 if strand == '+' else -1
    if strand not in ('+', '-'):
        return
    x = x_start + (position - region[1]) * total_width / (region[2] - region[1])
    points = [[4, x], [8, x], [8, x + direction * 4.9248], [9, x + direction * 4.9248],
              [7.65, x + direction * 6.156], [6.3, x + direction * 4.9248],
              [7.3, x + direction * 4.9248], [7.3, x + direction * 0.7], [4, x + direction * 0.7]]
    previous = y_start
    try:
        y_start = row_y
        write_to_file(draw_polygon(points, 1, '#000000', 0))
    finally:
        y_start = previous


def build_parser():
    parser = argparse.ArgumentParser(description='SparK: genomic signal plots as SVG')
    parser.add_argument('-pt','--plot_type', help='choices: standard, STD, sine', choices=['standard', 'STD', 'sine'], required=False, type=str, default="standard")
    parser.add_argument('-ps','--show_plots', help='choices: all, averages', choices=['all', 'averages'], required=False, type=str, default="all")
    parser.add_argument('-pr','--region', help='zero-based half-open interval, e.g. chr1:10000-20000', required=True, type=str)
    parser.add_argument('-cf','--control_files', help='separate by space', required=True, nargs='+', type=str)
    parser.add_argument('-tf','--treat_files', help='separate by space', required=False, nargs='+', type=str, default=[])
    parser.add_argument('-cg','--control_groups', help='group numbers separate by space', required=False, nargs='+', type=int, default=[])
    parser.add_argument('-tg','--treat_groups', help='group numbers separate by space', required=False, nargs='+', type=int, default=[])
    parser.add_argument('-gl','--group_labels', help='set group labels', required=False, nargs='+', type=str)
    parser.add_argument('-l','--labels', help='set labels for controls and treatment', required=False, nargs='+', type=str)
    parser.add_argument('-gs','--group_autoscale', help='set to "yes" to autoscale all groups(plots), except the groups excluded in -eg', choices=['yes', 'no'], required=False, type=str)
    parser.add_argument('-es','--exclude_from_group_autoscale', help='group numbers of groups to be excluded from autoscale', required=False, nargs='+', type=int)
    parser.add_argument('-eg','--exclude_groups', help='Exclude groups from the analysis', required=False, nargs='+', type=int)
    parser.add_argument('-f','--fills', help='one color, control/treatment colors, a preset, or more than two per-input colors', required=False, nargs='+', type=str, default=None)
    parser.add_argument('-cs','--custom_y_axis_scales', help='Enter y-axis values for all groups(plots). All groups need a value. Enter "D" for each group no value should be assigned, e.g. to keep autoscaling functionality for some groups', required=False, nargs='+', type=str)
    parser.add_argument('-dc','--display_chrom_location', help='set to "no" if you do not want to plot the chromosomal coordinates', choices=['no', 'top_left', 'top_right', 'bottom_left', 'bottom_right'], required=False, type=str, default="top_left")
    parser.add_argument('-gtf', '--gtffile', help='link gtf file for drawing genes here', required=False, type=str)
    parser.add_argument('-tss', '--drawtss', help='set to no to hide translation start-codon arrows (legacy -tss behavior)', choices=['yes', 'no'], required=False, type=str, default="yes")
    parser.add_argument('-genestart', '--draw_genestart', help='set to yes to indicate gene/transcript starts (TSS)', choices=['yes', 'no'], required=False, type=str, default="no")
    parser.add_argument('-sp', '--spark', help='display significant change "yes"', choices=['yes', 'no'], required=False, type=str)
    parser.add_argument('-sc', '--spark_color', help='spark color', required=False, type=str, nargs='+')
    parser.add_argument('-sm', '--smoothen', help='smoothen tracks, int', required=False, type=int)
    parser.add_argument('-o','--output_name', help='output graph name, str', required=False, type=str)
    parser.add_argument('-bed','--bed_files', help='bed files to be plotted', required=False, type=str, nargs='+')
    parser.add_argument('-bedcol','--bed_color', help='colors of bed files in hex', required=False, type=str, nargs='+')
    parser.add_argument('-bedlab','--bed_labels', help='set labels for bed tracks', required=False, nargs='+', type=str)
    parser.add_argument('-w','--track_width', help='width of the track, default = 150, int', required=False, type=int, default=150)
    parser.add_argument('-dg','--display_genes', help='genes to display from the gtf file', nargs='+', required=False, type=str)
    parser.add_argument('--max_points', help='maximum datapoints per plot', required=False, type=int, default=2000)
    parser.add_argument('--point_stat', help='per-bin statistic for downsampling: average or max', choices=['average', 'max'], required=False, type=str, default='average')
    parser.add_argument('--y_scale', help='scale factor for y-axis', required=False, type=float, default=1.0)
    parser.add_argument('--x_scale', help='scale factor for x-axis', required=False, type=float, default=1.0)
    parser.add_argument('-dt','--display_transcripts', help='display custom transcripts. By default, all transcripts annotated in the gtf file will be merged and displayed as one gene. Alternatively all can be plotted seperatelly by setting this to "all". Further, Transcript IDs can be listed to plot only certain transcripts', nargs='+', required=False, type=str, default=["mergeall"])
    parser.add_argument('-wg','--write_genenames', help='write genename instead of transcript ID when transcripts are plotted. Set to "yes".', choices=['yes', 'no'], required=False, type=str, default="no")

    parser.add_argument('-scale','--display_scalebar', help='set to "no" to remove scalebar', choices=['yes', 'no'], required=False, type=str, default="yes")
    parser.add_argument('--track_colors', nargs='+', help='one hex color per input, controls first then treatments; unambiguous for one or two inputs')
    parser.add_argument('--version', action='version', version=SparK_Version)
    return parser


if __name__ == '__main__':
    parser = build_parser()
    try:
        args = vars(parser.parse_args())
        region, active_groups, color_config, output_filename = validate_arguments(args)
        print('SparK Version ' + SparK_Version + ' initiated')
        print('Plotting region: ' + args['region'])
        x_scale, y_scale = args['x_scale'], args['y_scale']
        hight = track_hight = 30 * y_scale
        font_size_axis_y, font_size_axis_x = 8, 7
        x_start = 50
        spark_opacity, stroke_width, stroke_width_spark = 1, 0, 0
        plot_all_TSS = False
        relative_track_hight_percentage = 0.9
        bed_files, bed_labels = args['bed_files'], args['bed_labels'] or args['bed_files']
        bed_color = args['bed_color'] or ['#0B34FF'] * len(bed_files or [])
        draw_TSS = args['drawtss']
        display_genestart = args['draw_genestart']
        display_transcripts = args['display_transcripts']
        plot_gene_name_instead_transcriptID = args['write_genenames']
        smoothen_tracks, max_points, point_stat = args['smoothen'], args['max_points'], args['point_stat']
        display_genes, gff_file = args['display_genes'], args['gtffile']
        total_width = int(args['track_width'] * x_scale)
        fills, opacity, all_control_colors, all_treat_colors, individual_colors = color_config
        plot_type, show_plots, spark = args['plot_type'], args['show_plots'], args['spark']
        display_scalebar = args['display_scalebar']
        spark_color = args['spark_color'] or ['#FF0000', '#00FF00']
        if args['spark_color'] is None:
            stroke_width_spark, spark_opacity = 0.05, 0.5
        all_control_files, all_treat_files = args['control_files'], args['treat_files']
        control_groups, treat_groups = args['control_groups'], args['treat_groups']
        nr_of_groups = max(control_groups + treat_groups)
        custom_scales, group_labels, labels = args['custom_y_axis_scales'], args['group_labels'], args['labels']
        group_autoscale = args['group_autoscale']
        group_autoscale_excluded = args['exclude_from_group_autoscale'] or []
        exclude_groups = args['exclude_groups'] or []

        # Read each input once, including when a shared scale is requested.
        loaded = {}
        for filename in all_control_files + all_treat_files:
            if filename not in loaded:
                loaded[filename] = make_raw_data_filled(region, [filename], 0)[0]
        group_datasets = {}
        for group_id in active_groups:
            ctrl = [loaded[f] for f, g in zip(all_control_files, control_groups) if g == group_id]
            treat = [loaded[f] for f, g in zip(all_treat_files, treat_groups) if g == group_id]
            group_datasets[group_id] = ctrl, treat
        autoscale_max = max([get_max_value(*group_datasets[g]) for g in active_groups
                             if g not in group_autoscale_excluded] or [0])

        initial_canvas_height = 150 + track_hight * 2 * len(active_groups)
        write_to_file('<svg style="background:#ffffff" viewBox="0 0 320 {}" xmlns="http://www.w3.org/2000/svg">'.format(initial_canvas_height))
        quantile = float(total_width) / min(region[2] - region[1], max_points)
        additional_hight = 0
        first_track_y_start = None
        group_layouts = []
        for plot_index, group_id in enumerate(active_groups):
            group = group_id - 1  # Custom scales and labels retain their original group numbers.
            y_start = 100 + plot_index * hight * 1.5 + additional_hight
            control_data, treat_data = group_datasets[group_id]
            control_colors = [color for color, g in zip(all_control_colors, control_groups) if g == group_id]
            treat_colors = [color for color, g in zip(all_treat_colors, treat_groups) if g == group_id]
            has_negative = any(any(v < 0 for v in data) for data in control_data + treat_data)
            if plot_type == 'STD':
                has_negative = any(np.any(np.mean(data, axis=0) - np.std(data, axis=0) < 0)
                                   for data in (control_data, treat_data) if data)
            if first_track_y_start is None:
                first_track_y_start = y_start
            max_value = autoscale_max if group_autoscale == 'yes' and group_id not in group_autoscale_excluded else get_max_value(control_data, treat_data)
            if custom_scales is not None and custom_scales[group] != 'D':
                max_value = float(custom_scales[group])

            if plot_type == 'standard':
                for datasets, colors in ((control_data, control_colors), (treat_data, treat_colors)):
                    if not datasets:
                        continue
                    if show_plots == 'averages':
                        draw_signal(np.mean(datasets, axis=0), colors[0])
                    else:
                        for data, color in zip(datasets, colors):
                            draw_signal(data, color)
                if spark == 'yes':
                    draw_standard_spark()
                draw_axis_for_group(y_start, max_value, has_negative)
            elif plot_type == 'STD':
                for datasets, colors in ((control_data, control_colors), (treat_data, treat_colors)):
                    if not datasets:
                        continue
                    values = np.asarray(datasets)
                    means, deviations = np.mean(values, axis=0), np.std(values, axis=0)
                    coords, mean_line = [], []
                    for x, (mean, deviation) in enumerate(zip(means, deviations)):
                        x_pos = x_start + x * quantile
                        coords.append([get_relative_hight(mean + deviation), x_pos])
                        coords.insert(0, [get_relative_hight(mean - deviation), x_pos])
                        mean_line.append([get_relative_hight(mean), x_pos])
                    coords.extend([[get_relative_hight(means[-1] + deviations[-1]), x_start + total_width],
                                   [get_relative_hight(means[-1] - deviations[-1]), x_start + total_width]])
                    mean_line.append([get_relative_hight(means[-1]), x_start + total_width])
                    write_to_file(draw_polygon(coords, 0.4, colors[0], stroke_width))
                    write_to_file(draw_line(mean_line, 0.3, colors[0]))
                if spark == 'yes':
                    draw_standard_spark()
                draw_axis_for_group(y_start, max_value, has_negative)
            else:  # sine: both directions share one magnitude scale.
                shared_max_value = max_value
                for datasets, colors, sign in ((control_data, control_colors, -1), (treat_data, treat_colors, 1)):
                    if show_plots == 'averages':
                        draw_signal(np.mean(datasets, axis=0), colors[0], shared_max_value, sign)
                    else:
                        for data, color in zip(datasets, colors):
                            draw_signal(data, color, shared_max_value, sign)
                if spark == 'yes':
                    draw_standard_spark(sine=True, shared_max=shared_max_value)
                draw_sine_axis_for_group(y_start, shared_max_value)
            lower_half = has_negative or plot_type == 'sine'
            group_layouts.append((group_id, y_start, lower_half))
            if lower_half:
                additional_hight += hight

        if display_scalebar == 'yes':
            span = region[2] - region[1]
            target = max(1, span * 0.1)
            order = 10 ** math.floor(math.log10(target))
            scale_length = max(1, max(v for v in (1, 2, 5, 10) if v * order <= target) * order)
            normalized_scalebar_width = scale_length * total_width / span
            units = ('bp', 'kb', 'Mb', 'Gb', 'Tb', 'Pb')
            unit_index = min(int(math.log10(scale_length) // 3), len(units) - 1)
            scalebar_display = '{:g} {}'.format(scale_length / 1000 ** unit_index, units[unit_index])
            write_to_file(draw_rect(total_width + x_start - normalized_scalebar_width, 100-hight+8, '#000000', normalized_scalebar_width, 1, 1))
            write_to_file('<text text-anchor="middle" font-family="Arial" x="{}" y="{}" font-size="{}" >{}</text>'.format(
                total_width + x_start - normalized_scalebar_width/2, 96-hight+8, font_size_axis_x, scalebar_display))

        if group_labels is not None:
            for group_id, track_y, _ in group_layouts:
                write_to_file('<text text-anchor="start" font-family="Arial" x="{}" y="{}" font-size="8" >{}</text>'.format(
                    x_start + total_width + 15, track_y - int(hight/2) + 2.401, escape(group_labels[group_id-1])))

        # Squares and labels for groups
        if labels is not None:
            write_to_file(draw_rect(x_start - 10.5, 34, fills[0], 10, 10, opacity))
            write_to_file('''<text text-anchor="start" font-family="Arial" x="''' + str(x_start + 3) + '''" y="''' + str(34 - 1.788) + '''" font-size="8" >''' + escape(str(labels[0])) + '''</text>''')
            write_to_file(draw_rect(x_start - 10.5, 47, fills[1], 10, 10, opacity))
            write_to_file('''<text text-anchor="start" font-family="Arial" x="''' + str(x_start + 3) + '''" y="''' + str(47 - 1.788) + '''" font-size="8" >''' + escape(str(labels[1])) + '''</text>''')

            if show_plots == "averages":
                write_to_file(draw_rect(x_start - 10.5, 60, fills[0], 10, 10, opacity))
                write_to_file('''<text text-anchor="start" font-family="Arial" x="''' + str(x_start + 3) + '''" y="''' + str(60 - 1.788) + '''" font-size="8" >''' + "Overlap" + '''</text>''')
                write_to_file(draw_rect(x_start - 10.5, 60, fills[1], 10, 10, opacity))

            if spark == "yes":
                write_to_file(draw_rect(x_start + 59.5, 34, spark_color[1], 10, 10, 0.75))
                write_to_file('''<text text-anchor="start" font-family="Arial" x="''' + str(x_start + 73) + '''" y="''' + str(34 - 1.788) + '''" font-size="8" >''' + escape(str(labels[0])) + ''' up</text>''')
                write_to_file(draw_rect(x_start + 59.5, 47, spark_color[0], 10, 10, 0.75))
                write_to_file('''<text text-anchor="start" font-family="Arial" x="''' + str(x_start + 73) + '''" y="''' + str(47 - 1.788) + '''" font-size="8" >''' + escape(str(labels[1])) + ''' up</text>''')

        # chromosome location
        if args['display_chrom_location'] != "no":
            if args['display_chrom_location'] == "bottom_left":
                y_value_chr_label = y_start + (hight if has_negative else 0) + 9
                x_value_chr_label = x_start - 10.5
                text_anchor = "start"
            elif args['display_chrom_location'] == "top_left":
                y_value_chr_label = (first_track_y_start if first_track_y_start is not None else 100) - hight - 4
                x_value_chr_label = x_start - 10.5
                text_anchor = "start"
            elif args['display_chrom_location'] == "bottom_right":
                y_value_chr_label = y_start + (hight if has_negative else 0) + 9
                x_value_chr_label = x_start + total_width
                text_anchor = "end"
            elif args['display_chrom_location'] == "top_right":
                y_value_chr_label = (first_track_y_start if first_track_y_start is not None else 100) - hight - 4
                x_value_chr_label = x_start + total_width
                text_anchor = "end"
            write_to_file('''<text text-anchor="''' + text_anchor + '''''''" font-family="Arial" x="''' + str(x_value_chr_label) + '''" y="''' + str(y_value_chr_label) + '''" font-size="''' + str(font_size_axis_x) + '''" >Chr''' + escape(str(region[0])) + ''': ''' + str(f"{region[1]:,}") + '''-''' + str(f"{region[2]:,}") + '''</text>''')

        # add bed files
        y_position_bed = 110 + (len(active_groups) - 1) * hight * 1.5 + additional_hight # dirty. y_position_bed is also start for gtf genes
        if bed_files is not None:
            for nr_bed, bed_file in enumerate(bed_files):
                with open(bed_file) as b:
                    for line in b:
                        if not line.strip() or line.startswith(("#", "track ", "browser ")):
                            continue
                        region_to_draw = [0, 0]
                        line_split = line.split("\t")
                        if len(line_split) < 3:
                            raise ValueError("BED requires at least three tab-separated columns")
                        if line_split[0] != "":
                            if line_split[0][:3] == "chr" or line_split[0][:3] == "Chr":
                                line_split[0] = line_split[0][3:]
                            else:
                                pass
                        if line_split[0] == region[0]:  # if same chromosome
                            line_split[1] = int(line_split[1])
                            line_split[2] = int(line_split[2])
                            if line_split[1] < 0 or line_split[2] <= line_split[1]:
                                raise ValueError("Invalid BED interval")
                            if line_split[2] > region[1] and line_split[1] < region[2]:  # check if there is something to draw
                                region_to_draw = [line_split[1], line_split[2]]
                                if line_split[1] < region[1]:
                                    region_to_draw[0] = region[1]
                                if line_split[2] > region[2]:
                                    region_to_draw[1] = region[2]
                        if region_to_draw != [0, 0]:
                            write_to_file(draw_rect(x_start + (((region_to_draw[0] - region[1]) * total_width) / float((region[2] - region[1]))), y_position_bed - 0.3 + (2 / 2), bed_color[nr_bed], ((region_to_draw[1] - region_to_draw[0]) * total_width) / float(region[2] - region[1]), 2, 1))
                write_to_file('''<text text-anchor="start" font-family="Arial" x="''' + str(x_start + total_width + 15) + '''" y="''' + str(y_position_bed + 3) + '''" font-size="8" >''' + escape(bed_labels[nr_bed]) + '''</text>''')
                y_position_bed += 8

        # add gene plots from gtf file
        if gff_file is not None:
            plot_transcripts = False

            with contextlib.closing(get_lines_indexed(gff_file, region[0])) as f:
                gene_or_transcript_names_to_plot = []
                gene_rows = {}
                gene_row_bottom = y_position_bed
                transcripts = []
                tss_plotted_genes = []
                y_genestart = y_position_bed
                exons = []
                UTRs = []
                for line in f:
                    if not line.strip() or line.startswith("#"):
                        continue
                    line_split = line.rstrip().split("\t")
                    if len(line_split) != 9:
                        raise ValueError("GTF requires nine tab-separated columns")
                    try:
                        line_split[3] = int(line_split[3]) - 1
                    except ValueError as exc:
                        raise ValueError("GTF start must be an integer") from exc
                    try:
                        line_split[4] = int(line_split[4])
                    except ValueError as exc:
                        raise ValueError("GTF end must be an integer") from exc
                    if line_split[3] < 0 or line_split[4] <= line_split[3]:
                        raise ValueError("Invalid GTF interval")
                    try:
                        line_split[len(line_split) - 1] = line_split[len(line_split) - 1].split("\n")[0]
                    except:
                        pass
                    try:
                        if line_split[0][:3] == "chr" or line_split[0][:3] == "Chr":
                            line_split[0] = line_split[0][3:]
                    except:
                        pass
                    if line_split[0] == region[0]:  # check if same chr
                        if display_genes is not None and get_gene_name() not in display_genes:
                            continue
                        to_draw = get_region_to_draw()
                        if to_draw != 0:

                            # the following gets either the gene name, or the transcript name, plot_transcripts indicates what
                            if line_split[2] in ["gene", "CDS", "exon", "transcript", "start_codon"]:
                                gene_or_transcript_name = ""
                                if display_transcripts[0] == "mergeall":
                                    if display_genes is None:
                                        gene_or_transcript_name = get_gene_name() # Get gene name if all merged genes are to be plotted
                                    else:
                                        if get_gene_name() in display_genes:
                                            gene_or_transcript_name = get_gene_name() # Get gene name only certain merged genes are to be plotted
                                else: # if individual transcripts should be plotted
                                    if "transcript_id " in line_split[8]:
                                        if display_transcripts[0] == "all":
                                            gene_or_transcript_name = get_transcript_name() # Get transcript ID if all individual transcripts instead of genes are to be plotted
                                            plot_transcripts = True # mark that instead of only genes, transcripts will be plotted
                                        else:
                                            if get_transcript_name() in display_transcripts:
                                                gene_or_transcript_name = get_transcript_name()
                                                plot_transcripts = True  # mark that instead of only genes, transcripts will be plotted

                                plot_gene = False
                                if gene_or_transcript_name != "":  # only if transcript name is not empty it will be plotted
                                    if gene_or_transcript_name not in gene_or_transcript_names_to_plot: # checks which genes have been plotted
                                        if display_transcripts[0] != "mergeall":
                                            if display_transcripts[0] == "all":
                                                y_genestart += 12
                                                gene_or_transcript_names_to_plot.append(gene_or_transcript_name)
                                            else:
                                                if gene_or_transcript_name in display_transcripts:
                                                    y_genestart += 12
                                                    gene_or_transcript_names_to_plot.append(gene_or_transcript_name)
                                        else:
                                            y_genestart += 12
                                            gene_or_transcript_names_to_plot.append(gene_or_transcript_name)


                                    if display_transcripts[0] != "mergeall": # if transcripts to plot were defined
                                        if gene_or_transcript_name in display_transcripts: # check if this one should be plotted
                                            plot_gene = True
                                        elif display_transcripts[0] == "all":
                                            plot_gene = True

                                    elif display_genes is not None:
                                        if gene_or_transcript_name in display_genes:
                                            plot_gene = True
                                    else:
                                        plot_gene = True

                                if plot_gene == True:
                                    if gene_or_transcript_name not in gene_rows:
                                        gene_row_bottom += 12
                                        gene_rows[gene_or_transcript_name] = gene_row_bottom
                                    y_genestart = gene_rows[gene_or_transcript_name]
                                    draw_gene_or_transcript = False
                                    if plot_transcripts == True:
                                        if line_split[2] == "transcript":
                                            draw_gene_or_transcript = True #plot only transcripts
                                    elif line_split[2] == "gene":
                                        draw_gene_or_transcript = True # plot only genes if transcripts are not defined

                                    if draw_gene_or_transcript == True:
                                        if plot_gene_name_instead_transcriptID == "yes":
                                            gene_label = get_gene_name()
                                        else:
                                            gene_label = gene_or_transcript_name
                                        write_to_file('''<text text-anchor="start" font-family="Arial" x="''' + str(x_start + total_width + 15) + '''" y="''' + str(y_genestart + 3) + '''" font-size="8" >''' + escape(gene_label) + '''</text>''')
                                        drawhight = 1.0
                                    elif line_split[2] == "CDS":
                                        drawhight = 6.0
                                    elif line_split[2] == "exon":
                                        drawhight = 3.0

                                    if draw_gene_or_transcript and display_genestart == "yes":
                                        gene_start = line_split[3] if line_split[6] == "+" else line_split[4] - 1
                                        if region[1] <= gene_start < region[2]:
                                            draw_start_arrow(gene_start, line_split[6], y_genestart)

                                    if draw_gene_or_transcript == True or line_split[2] in ["CDS", "exon"]:
                                        write_to_file(draw_rect(x_start + (((to_draw[0] - region[1]) * total_width) / float((region[2] - region[1]))), y_genestart - 0.3 + (drawhight / 2), "#0B34FF", ((to_draw[1] - to_draw[0]) * total_width) / float(region[2] - region[1]), drawhight, 1))

                                    if line_split[2] == "start_codon" and region[1] <= line_split[3] < region[2]:
                                        if draw_TSS == "yes":
                                            plot = True
                                            if plot_transcripts == False:  # makes sure there is only one TSS plotted if merged gene is selected
                                                plot = False
                                                plotted_TSS = True
                                                if gene_or_transcript_name not in tss_plotted_genes:
                                                    plotted_TSS = False
                                                    tss_plotted_genes.append(gene_or_transcript_name)
                                                plot = False
                                                if plotted_TSS == False or plot_all_TSS == True:
                                                    plot = True
                                                if display_genes is not None and gene_or_transcript_name not in display_genes:
                                                    plot = False
                                            if plot == True:
                                                codon_start = line_split[3] if line_split[6] == "+" else line_split[4] - 1
                                                draw_start_arrow(codon_start, line_split[6], y_genestart)
                                                plotted_TSS = True



        write_to_file("</svg>")
        last_bottom = max(y + (track_hight if negative else 0) for _, y, negative in group_layouts)
        annotation_bottom = max(y_position_bed, globals().get('gene_row_bottom', y_position_bed))
        canvas_height = max(initial_canvas_height, last_bottom + 30, annotation_bottom + 25)
        finish_svg(max(320, x_start + total_width + 120), canvas_height)
        print('Saved ' + output_filename)

    except (ValueError, OSError, ImportError, RuntimeError, subprocess.SubprocessError) as exc:
        parser.exit(2, 'SparK error: {}\n'.format(exc))
