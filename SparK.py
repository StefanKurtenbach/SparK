SparK_Version = "3.2.0"

import argparse
import atexit
import contextlib
import copy
import datetime
import hashlib
import heapq
import importlib.metadata
import io
import json
import math
import mmap
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time
import urllib.request
from html import escape
import numpy as np
try:
    import pyBigWig
except ImportError:
    pyBigWig = None
import gzip


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
            yield (position + 1)
            position = mapped.find(needle, position + len(needle))
    variants = {chromosome, 'chr' + chromosome, 'Chr' + chromosome}
    streams = [positions(name.encode('utf-8') + b'\t') for name in sorted(variants)]
    return heapq.merge(*streams)

def _text_index(filename, chromosome, coordinate_columns, first_block_only):
    """Cache byte ranges, not data; stat changes invalidate the tiny index."""
    realpath = os.path.realpath(filename)
    stat = os.stat(realpath)
    signature = [realpath, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns, stat.st_ino, chromosome, coordinate_columns, first_block_only, 1]
    key = hashlib.sha256(json.dumps(signature).encode('utf-8')).hexdigest()
    cache_dir = os.path.join(tempfile.gettempdir(), 'spark-text-index-v1-' + str(getattr(os, 'getuid', lambda: 0)()))
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
                if re.search(b'\\r(?!\\n)', mapped):
                    raise ValueError('Standalone carriage-return record separators')
                for position in _chromosome_offsets(mapped, chromosome):
                    if first_block_only and previous_end is not None and (position != previous_end):
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
                            block[4] = True
                    count += 1
                    previous_end = end
    after = os.stat(realpath)
    if (after.st_size, after.st_mtime_ns, after.st_ctime_ns, after.st_ino) != (stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns, stat.st_ino):
        raise OSError('Input changed while its text index was being built')
    temporary = None
    try:
        os.makedirs(cache_dir, mode=448, exist_ok=True)
        with tempfile.NamedTemporaryFile(mode='w', dir=cache_dir, prefix=key + '.', suffix='.tmp', delete=False) as handle:
            temporary = handle.name
            json.dump({'signature': signature, 'blocks': blocks}, handle)
        os.replace(temporary, cache_file)
    except OSError:
        if temporary is not None:
            try:
                os.unlink(temporary)
            except OSError:
                pass
    return blocks

def get_lines_indexed(filename, chromosome, lower=None, upper=None, coordinate_columns=None, first_block_only=False):
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
        columns = list(coordinate_columns) if coordinate_columns else None
        blocks = _text_index(filename, chromosome, columns, first_block_only)
    except (OSError, ValueError, TypeError, KeyError):
        with open(filename) as source:
            yield from source
        return
    with open(filename, 'rb') as source:
        for start, stop, min_start, max_stop, unsafe in blocks:
            if coordinate_columns and (not unsafe) and (min_start is not None):
                if max_stop < lower or min_start > upper:
                    continue
            source.seek(start)
            while source.tell() < stop:
                yield source.readline().decode('utf-8').replace('\r\n', '\n').replace('\r', '\n')

def read_signals(stretch, files, offset=0, max_points=2000, point_stat='average', smoothen_tracks=None):
    """Read zero-based half-open intervals, clipping before assigning signal."""
    lower, upper = (stretch[1] + offset, stretch[2] + offset)
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
                start, stop = (max(0, lower), min(upper, chroms[chromname]))
                if start < stop:
                    values = bw.values(chromname, start, stop)
                    data[start - lower:stop - lower] = [float(v) if v is not None and math.isfinite(v) else 0.0 for v in values]
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
                        start, stop, value = (int(fields[1]), int(fields[2]), float(fields[3]))
                    except (ValueError, IndexError) as exc:
                        raise ValueError('Invalid bedGraph record in {}: {}'.format(filename, line.rstrip())) from exc
                    if start < 0 or stop <= start or (not math.isfinite(value)):
                        raise ValueError('Invalid bedGraph interval or value in {}: {}'.format(filename, line.rstrip()))
                    start, stop = (max(start, lower), min(stop, upper))
                    if start < stop:
                        data[start - lower:stop - lower] = [value] * (stop - start)
        if length > max_points:
            data = downsample(data, max_points, point_stat)
        if smoothen_tracks is not None:
            prefix = np.concatenate(([0.0], np.cumsum(data, dtype=float)))
            window = smoothen_tracks
            first = data[0]
            data = [(prefix[i + 1] - prefix[max(0, i - window + 1)] + max(0, window - i - 1) * first) / window for i in range(len(data))]
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
        left, right = (position * size, (position + 1) * size)
        start, stop = (left // bins, (right + bins - 1) // bins)
        if statistic == 'max':
            output.append(max(data[start:stop]))
        else:
            last = stop - 1
            first_weight = (start + 1) * bins - left
            last_weight = right - last * bins
            value = (data[start] * first_weight + math.fsum(data[start + 1:last]) * bins + data[last] * last_weight) / size
            output.append(value)
    return output

def hex_color(value):
    value = value.lstrip('#')
    if not re.fullmatch('(?:[0-9a-fA-F]{3}|[0-9a-fA-F]{6})', value):
        raise ValueError('Invalid hex color {!r}; use three or six hex digits, optionally prefixed by #'.format(value))
    return '#' + value

def parse_region(value):
    match = re.fullmatch('([^:\\s]+):(\\d+)-(\\d+)', value)
    if not match:
        raise ValueError('Region must use chr1:10000-20000 format (zero-based, half-open)')
    chrom, start, stop = (match.group(1), int(match.group(2)), int(match.group(3)))
    if stop <= start:
        raise ValueError('Region end must be greater than start')
    normalized = chrom[3:] if chrom.startswith(('chr', 'Chr')) else chrom
    if not normalized:
        raise ValueError('Chromosome name is empty')
    return [normalized, start, stop, chrom]

def configure_colors(args):
    nctrl, ntreat = (len(args['control_files']), len(args['treat_files']))
    values, explicit = (args['fills'], args['track_colors'])
    if values and explicit:
        raise ValueError('Use either -f or --track_colors, not both')
    presets = {'blue/red': (['#FF1800', '#005CFF'], 0.6), 'blue/grey': (['#005CFF', '#A3A3A3'], 0.3), 'all_grey': (['#848484', '#848484'], 0.6), 'blue/green': (['#00FF12', '#005CFF'], 0.5)}
    individual = explicit is not None or (values is not None and len(values) > 2)
    if individual:
        colors = [hex_color(v) for v in explicit or values]
        if len(colors) != nctrl + ntreat:
            raise ValueError('Per-track colors require one color per input: controls first, then treatments')
        control, treatment = (colors[:nctrl], colors[nctrl:])
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
        control, treatment = ([fills[0]] * nctrl, [fills[1]] * ntreat)
    return (fills, opacity, control, treatment, individual)


def make_raw_data_filled(stretch, files, offset):
    """Compatibility helper; new callers should use read_signals with explicit options."""
    return read_signals(stretch, files, offset, globals().get('max_points', 2000),
                        globals().get('point_stat', 'average'), globals().get('smoothen_tracks'))


class PlotRenderer:
    """An independent SVG renderer; safe for separate notebook/server calls."""
    def __init__(self, quiet=True):
        self.quiet = quiet
        self._svg_output = None
        self._svg_temporary = None

    def log(self, *values, **kwargs):
        if not self.quiet:
            print(*values, **kwargs, flush=True)

    def read_data(self, stretch, files, offset):
        return read_signals(stretch, files, offset, self.max_points, self.point_stat, self.smoothen_tracks)

    def draw_condition_legend(self, args, groups):
        names = display_conditions(args, groups)
        for index, (name, color) in enumerate(names):
            y = 34 + index * 13
            self.write_to_file(self.draw_rect(self.x_start - 10.5, y, color, 10, 10, self.opacity))
            self.write_to_file('<text text-anchor="start" font-family="Arial" x="{}" y="{}" font-size="8" >{}</text>'.format(self.x_start+3, y-1.788, escape(name)))


    def get_gene_name(self):
        attributes = dict(re.findall('(?:^|;)\\s*(\\S+)\\s+"([^"]*)"', self.line_split[8]))
        name = attributes.get('gene_name') or attributes.get('gene_id')
        if not name:
            raise ValueError('GTF annotation requires gene_name or gene_id in column 9')
        return name

    def get_transcript_name(self):
        match = re.search('(?:^|;)\\s*transcript_id\\s+"([^"]+)"', self.line_split[8])
        if not match:
            raise ValueError('Transcript plotting requires transcript_id in GTF column 9')
        return match.group(1)

    def draw_line(self, coordinates, thickness, color):
        output = '<path d="'
        for x, i in enumerate(coordinates):
            if x == 0:
                output += 'M '
            else:
                output += ' L '
            output += str(i[1]) + ' ' + str(self.y_start - i[0])
        output += '" style="stroke:' + color + '; stroke-width:' + str(thickness) + '; fill:none "/>'
        return output

    def write_to_file(self, row):
        if self._svg_output is None:
            self._svg_output = tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=os.path.dirname(os.path.abspath(self.output_filename)), prefix='.spark-', suffix='.svg', delete=False)
            self._svg_temporary = self._svg_output.name
            atexit.register(self.discard_incomplete_svg)
        self._svg_output.write(row + '\n')

    def discard_incomplete_svg(self):
        if self._svg_output is not None and (not self._svg_output.closed):
            self._svg_output.close()
        if self._svg_temporary and os.path.exists(self._svg_temporary):
            os.unlink(self._svg_temporary)

    def finish_svg(self, width, height):
        """Publish only a complete SVG, retaining the previous output on failure."""
        self._svg_output.close()
        with open(self._svg_temporary, encoding='utf-8') as handle:
            content = handle.read()
        content = re.sub('viewBox="[^"]+"', 'viewBox="0 0 {} {}"'.format(width, height), content, count=1)
        with open(self._svg_temporary, 'w', encoding='utf-8') as handle:
            handle.write(content)
        os.replace(self._svg_temporary, self.output_filename)

    def get_max_value(self, datasets1=(), datasets2=(), *more):
        maxima = [0.0]
        for datasets in (datasets1, datasets2) + more:
            if not datasets:
                continue
            values = np.asarray(datasets, dtype=float)
            if self.plot_type == 'STD':
                mean, std = (np.mean(values, axis=0), np.std(values, axis=0))
                maxima.extend([float(np.max(np.abs(mean + std))), float(np.max(np.abs(mean - std)))])
            elif self.show_plots == 'averages' and self.plot_type != 'sine':
                maxima.append(float(np.max(np.abs(np.mean(values, axis=0)))))
            else:
                maxima.append(float(np.max(np.abs(values))))
        return max(maxima)

    def get_relative_hight(self, raw_value):
        if raw_value == 0:
            return 0
        else:
            return raw_value * self.hight * self.relative_track_hight_percentage / self.max_value

    def get_relative_hight_custom(self, raw_value, denom):
        if raw_value == 0 or denom == 0:
            return 0
        return raw_value * self.hight * self.relative_track_hight_percentage / denom

    def draw_rect(self, x_coord, y_0, color, width, hight1, opacity):
        return '<rect x="' + str(x_coord) + '" opacity="' + str(opacity) + '" y="' + str(y_0 - hight1) + '" fill="' + color + '" width="' + str(width) + '" height="' + str(hight1) + '"/>'

    def draw_polygon(self, coordinates, opacity, color, stroke_width):
        for q, w in enumerate(coordinates):
            coordinates[q][0] = self.y_start - coordinates[q][0]
        string = '<polygon points="'
        for h, c in enumerate(coordinates):
            if h == 0:
                string += str(c[1]) + ',' + str(c[0])
            else:
                string += ' ' + str(c[1]) + ',' + str(c[0])
        string += '" opacity="' + str(opacity) + '" fill="' + color + '"' + ' stroke="black" stroke-width="' + str(stroke_width) + '"/>'
        return string

    def draw_axis_for_group(self, y_start_val, max_value_val, has_negative_axis):
        axis_label = round(max_value_val / self.relative_track_hight_percentage, 1)
        if has_negative_axis:
            self.write_to_file('<line x1="' + str(self.x_start - 10) + '" y1="' + str(y_start_val + self.hight) + '" x2="' + str(self.x_start - 10) + '" y2="' + str(y_start_val - self.hight) + '" stroke="black" stroke-width="1" />')
            self.write_to_file('<line x1="' + str(self.x_start - 10.5) + '" y1="' + str(y_start_val) + '" x2="' + str(self.x_start - 6.5) + '" y2="' + str(y_start_val) + '" stroke="black" stroke-width="1" />')
            self.write_to_file('<line x1="' + str(self.x_start - 10.5) + '" y1="' + str(y_start_val - self.hight) + '" x2="' + str(self.x_start - 6.5) + '" y2="' + str(y_start_val - self.hight) + '" stroke="black" stroke-width="1" />')
            self.write_to_file('<line x1="' + str(self.x_start - 10.5) + '" y1="' + str(y_start_val + self.hight) + '" x2="' + str(self.x_start - 6.5) + '" y2="' + str(y_start_val + self.hight) + '" stroke="black" stroke-width="1" />')
            self.write_to_file('<text text-anchor="end" font-family="Arial" x="' + str(self.x_start - 14) + '" y="' + str(y_start_val + self.hight + 4) + '" font-size="' + str(self.font_size_axis_y) + '" >' + str(-axis_label) + '</text>')
            self.write_to_file('<text text-anchor="end" font-family="Arial" x="' + str(self.x_start - 14) + '" y="' + str(y_start_val + 4) + '" font-size="' + str(self.font_size_axis_y) + '" >0</text>')
            self.write_to_file('<text text-anchor="end" font-family="Arial" x="' + str(self.x_start - 14) + '" y="' + str(y_start_val - self.hight + 4) + '" font-size="' + str(self.font_size_axis_y) + '" >' + str(axis_label) + '</text>')
        else:
            self.write_to_file('<line x1="' + str(self.x_start - 10) + '" y1="' + str(y_start_val) + '" x2="' + str(self.x_start - 10) + '" y2="' + str(y_start_val - self.hight) + '" stroke="black" stroke-width="1" />')
            self.write_to_file('<line x1="' + str(self.x_start - 10.5) + '" y1="' + str(y_start_val) + '" x2="' + str(self.x_start - 6.5) + '" y2="' + str(y_start_val) + '" stroke="black" stroke-width="1" />')
            self.write_to_file('<line x1="' + str(self.x_start - 10.5) + '" y1="' + str(y_start_val - self.hight) + '" x2="' + str(self.x_start - 6.5) + '" y2="' + str(y_start_val - self.hight) + '" stroke="black" stroke-width="1" />')
            self.write_to_file('<text text-anchor="end" font-family="Arial" x="' + str(self.x_start - 14) + '" y="' + str(y_start_val + 4) + '" font-size="' + str(self.font_size_axis_y) + '" >0</text>')
            self.write_to_file('<text text-anchor="end" font-family="Arial" x="' + str(self.x_start - 14) + '" y="' + str(y_start_val - self.hight + 4) + '" font-size="' + str(self.font_size_axis_y) + '" >' + str(axis_label) + '</text>')

    def draw_sine_axis_for_group(self, y_start_val, shared_max_val):
        if shared_max_val >= 0:
            axis_height = self.hight * self.relative_track_hight_percentage
            pos_tick_y = y_start_val - axis_height
            neg_tick_y = y_start_val + axis_height
            self.write_to_file('<line x1="' + str(self.x_start - 10) + '" y1="' + str(neg_tick_y) + '" x2="' + str(self.x_start - 10) + '" y2="' + str(pos_tick_y) + '" stroke="black" stroke-width="1" />')
            self.write_to_file('<line x1="' + str(self.x_start - 10.5) + '" y1="' + str(y_start_val) + '" x2="' + str(self.x_start - 6.5) + '" y2="' + str(y_start_val) + '" stroke="black" stroke-width="1" />')
            self.write_to_file('<line x1="' + str(self.x_start - 10.5) + '" y1="' + str(pos_tick_y) + '" x2="' + str(self.x_start - 6.5) + '" y2="' + str(pos_tick_y) + '" stroke="black" stroke-width="1" />')
            self.write_to_file('<text text-anchor="end" font-family="Arial" x="' + str(self.x_start - 14) + '" y="' + str(pos_tick_y + 4) + '" font-size="' + str(self.font_size_axis_y) + '" >' + str(round(shared_max_val, 1)) + '</text>')
            self.write_to_file('<line x1="' + str(self.x_start - 10.5) + '" y1="' + str(neg_tick_y) + '" x2="' + str(self.x_start - 6.5) + '" y2="' + str(neg_tick_y) + '" stroke="black" stroke-width="1" />')
            self.write_to_file('<text text-anchor="end" font-family="Arial" x="' + str(self.x_start - 14) + '" y="' + str(neg_tick_y + 4) + '" font-size="' + str(self.font_size_axis_y) + '" >-' + str(round(shared_max_val, 1)) + '</text>')

    def draw_standard_spark(self, sine=False, shared_max=None):
        if len(self.control_data) < 2 or len(self.treat_data) < 2:
            raise ValueError('Sparks require at least two control and two treatment files in every plotted group')
        coordinates, direction = ([], None)

        def emit():
            if not coordinates:
                return
            color = self.spark_color[1] if direction == 'control' else self.spark_color[0]
            points = list(coordinates)
            if sine:
                points = [[0, points[0][1]]] + points + [[points[-1][0], points[-1][1] + self.quantile], [0, points[-1][1] + self.quantile]]
            else:
                points.extend([[points[-1][0], points[-1][1] + self.quantile], [points[0][0], points[0][1] + self.quantile]])
            self.write_to_file(self.draw_polygon(points, self.spark_opacity, color, self.stroke_width_spark))
        for x, (ctrl, treat) in enumerate(zip(np.asarray(self.control_data).T, np.asarray(self.treat_data).T)):
            ctrl_mean, treat_mean = (float(np.mean(ctrl)), float(np.mean(treat)))
            ctrl_std, treat_std = (float(np.std(ctrl)), float(np.std(treat)))
            delta = ctrl_mean - treat_mean
            current = ('control' if delta > 0 else 'treatment') if abs(delta) > ctrl_std + treat_std else None
            if current != direction or current is None:
                emit()
                coordinates, direction = ([], current)
            if current is None:
                continue
            x_pos = self.x_start + x * self.quantile
            if sine:
                value = abs(delta) - ctrl_std - treat_std
                value = -value if current == 'control' else value
                coordinates.append([self.get_relative_hight_custom(value, shared_max), x_pos])
            else:
                low, high = (treat_mean + treat_std, ctrl_mean - ctrl_std) if current == 'control' else (ctrl_mean + ctrl_std, treat_mean - treat_std)
                coordinates.append([self.get_relative_hight(high), x_pos])
                coordinates.insert(0, [self.get_relative_hight(low), x_pos])
        emit()

    def get_region_to_draw(self):
        region_to_draw = 0
        if self.line_split[4] > self.region[1] and self.line_split[3] < self.region[2]:
            region_to_draw = [float(self.line_split[3]), float(self.line_split[4])]
            if self.line_split[3] < self.region[1]:
                region_to_draw[0] = float(self.region[1])
            if self.line_split[4] > self.region[2]:
                region_to_draw[1] = float(self.region[2])
        return region_to_draw

    def draw_signal(self, data, color, denominator=None, sign=1):
        denominator = self.max_value if denominator is None else denominator
        points = [[self.get_relative_hight_custom(value * sign, denominator), self.x_start + x * self.quantile] for x, value in enumerate(data)]
        points = [[0, self.x_start]] + points + [[points[-1][0], self.x_start + self.total_width], [0, self.x_start + self.total_width]]
        self.write_to_file(self.draw_polygon(points, self.opacity, color, self.stroke_width))

    def draw_start_arrow(self, position, strand, row_y):
        direction = 1 if strand == '+' else -1
        if strand not in ('+', '-'):
            return
        x = self.x_start + (position - self.region[1]) * self.total_width / (self.region[2] - self.region[1])
        points = [[4, x], [8, x], [8, x + direction * 4.9248], [9, x + direction * 4.9248], [7.65, x + direction * 6.156], [6.3, x + direction * 4.9248], [7.3, x + direction * 4.9248], [7.3, x + direction * 0.7], [4, x + direction * 0.7]]
        previous = self.y_start
        try:
            self.y_start = row_y
            self.write_to_file(self.draw_polygon(points, 1, '#000000', 0))
        finally:
            self.y_start = previous

    def render(self, args):
        self.region, active_groups, color_config, self.output_filename = finalize_arguments(args)
        self.log('SparK Version ' + SparK_Version + ' initiated')
        self.log('Plotting region: ' + args['region'])
        x_scale, y_scale = (args['x_scale'], args['y_scale'])
        self.hight = track_hight = 30 * y_scale
        self.font_size_axis_y, font_size_axis_x = (8, 7)
        self.x_start = 50
        self.spark_opacity, self.stroke_width, self.stroke_width_spark = (1, 0, 0)
        plot_all_TSS = False
        self.relative_track_hight_percentage = 0.9
        bed_files, bed_labels = (args['bed_files'], args['bed_labels'] or args['bed_files'])
        bed_color = args['bed_color'] or ['#0B34FF'] * len(bed_files or [])
        draw_TSS = args['drawtss']
        display_genestart = args['draw_genestart']
        display_transcripts = args['display_transcripts']
        plot_gene_name_instead_transcriptID = args['write_genenames']
        self.smoothen_tracks, self.max_points, self.point_stat = (args['smoothen'], args['max_points'], args['point_stat'])
        display_genes, gff_file = (args['display_genes'], args['gtffile'])
        self.total_width = int(args['track_width'] * x_scale)
        fills, self.opacity, all_control_colors, all_treat_colors, individual_colors = color_config
        self.plot_type, self.show_plots, spark = (args['plot_type'], args['show_plots'], args['spark'])
        display_scalebar = args['display_scalebar']
        self.spark_color = args['spark_color'] or ['#FF0000', '#00FF00']
        if args['spark_color'] is None:
            self.stroke_width_spark, self.spark_opacity = (0.05, 0.5)
        all_control_files, all_treat_files = (args['control_files'], args['treat_files'])
        control_groups, treat_groups = (args['control_groups'], args['treat_groups'])
        nr_of_groups = max(control_groups + treat_groups)
        custom_scales, group_labels, labels = (args['custom_y_axis_scales'], args['group_labels'], args['labels'])
        group_autoscale = args['group_autoscale']
        group_autoscale_excluded = args['exclude_from_group_autoscale'] or []
        exclude_groups = args['exclude_groups'] or []
        loaded = {}
        for filename in all_control_files + all_treat_files:
            if filename not in loaded:
                loaded[filename] = self.read_data(self.region, [filename], 0)[0]
        series_by_group = build_group_series(args, loaded)
        autoscale_max = max([self.get_max_value(*[entry['data'] for entry in series_by_group[g]]) for g in active_groups if g not in group_autoscale_excluded] or [0])
        header_offset = max(0, (len(args.get('_conditions_order', [])) - 3) * 13) if args.get('_named_tracks') else 0
        initial_canvas_height = 150 + header_offset + track_hight * 2 * len(active_groups)
        self.write_to_file('<svg style="background:#ffffff" viewBox="0 0 320 {}" xmlns="http://www.w3.org/2000/svg">'.format(initial_canvas_height))
        self.quantile = float(self.total_width) / min(self.region[2] - self.region[1], self.max_points)
        additional_hight = 0
        first_track_y_start = None
        group_layouts = []
        for plot_index, group_id in enumerate(active_groups):
            group = group_id - 1
            self.y_start = 100 + header_offset + plot_index * self.hight * 1.5 + additional_hight
            series = series_by_group[group_id]
            self.control_data, self.treat_data = comparison_data(args, series)
            control_colors = next((entry['colors'] for entry in series if entry['name'] == args['_compare_reference']), [])
            treat_colors = next((entry['colors'] for entry in series if entry['name'] == args['_compare_target']), [])
            has_negative = any((any((v < 0 for v in data)) for data in [data for entry in series for data in entry['data']]))
            if self.plot_type == 'STD':
                has_negative = any((np.any(np.mean(data, axis=0) - np.std(data, axis=0) < 0) for data in [entry['data'] for entry in series] if data))
            if first_track_y_start is None:
                first_track_y_start = self.y_start
            self.max_value = autoscale_max if group_autoscale == 'yes' and group_id not in group_autoscale_excluded else self.get_max_value(*[entry['data'] for entry in series])
            if custom_scales is not None and custom_scales[group] != 'D':
                self.max_value = float(custom_scales[group])
            if self.plot_type == 'standard':
                for datasets, colors in ((entry['data'], entry['colors']) for entry in series):
                    if not datasets:
                        continue
                    if self.show_plots == 'averages':
                        self.draw_signal(np.mean(datasets, axis=0), colors[0])
                    else:
                        for data, color in zip(datasets, colors):
                            self.draw_signal(data, color)
                if spark == 'yes':
                    self.draw_standard_spark()
                self.draw_axis_for_group(self.y_start, self.max_value, has_negative)
            elif self.plot_type == 'STD':
                for datasets, colors in ((entry['data'], entry['colors']) for entry in series):
                    if not datasets:
                        continue
                    values = np.asarray(datasets)
                    means, deviations = (np.mean(values, axis=0), np.std(values, axis=0))
                    coords, mean_line = ([], [])
                    for x, (mean, deviation) in enumerate(zip(means, deviations)):
                        x_pos = self.x_start + x * self.quantile
                        coords.append([self.get_relative_hight(mean + deviation), x_pos])
                        coords.insert(0, [self.get_relative_hight(mean - deviation), x_pos])
                        mean_line.append([self.get_relative_hight(mean), x_pos])
                    coords.extend([[self.get_relative_hight(means[-1] + deviations[-1]), self.x_start + self.total_width], [self.get_relative_hight(means[-1] - deviations[-1]), self.x_start + self.total_width]])
                    mean_line.append([self.get_relative_hight(means[-1]), self.x_start + self.total_width])
                    self.write_to_file(self.draw_polygon(coords, 0.4, colors[0], self.stroke_width))
                    self.write_to_file(self.draw_line(mean_line, 0.3, colors[0]))
                if spark == 'yes':
                    self.draw_standard_spark()
                self.draw_axis_for_group(self.y_start, self.max_value, has_negative)
            else:
                shared_max_value = self.max_value
                for datasets, colors, sign in ((entry['data'], entry['colors'], -1 if entry['name'] == args['_compare_reference'] else 1) for entry in series):
                    if self.show_plots == 'averages':
                        self.draw_signal(np.mean(datasets, axis=0), colors[0], shared_max_value, sign)
                    else:
                        for data, color in zip(datasets, colors):
                            self.draw_signal(data, color, shared_max_value, sign)
                if spark == 'yes':
                    self.draw_standard_spark(sine=True, shared_max=shared_max_value)
                self.draw_sine_axis_for_group(self.y_start, shared_max_value)
            lower_half = has_negative or self.plot_type == 'sine'
            group_layouts.append((group_id, self.y_start, lower_half))
            if lower_half:
                additional_hight += self.hight
        if display_scalebar == 'yes':
            span = self.region[2] - self.region[1]
            target = max(1, span * 0.1)
            order = 10 ** math.floor(math.log10(target))
            scale_length = max(1, max((v for v in (1, 2, 5, 10) if v * order <= target)) * order)
            normalized_scalebar_width = scale_length * self.total_width / span
            units = ('bp', 'kb', 'Mb', 'Gb', 'Tb', 'Pb')
            unit_index = min(int(math.log10(scale_length) // 3), len(units) - 1)
            scalebar_display = '{:g} {}'.format(scale_length / 1000 ** unit_index, units[unit_index])
            self.write_to_file(self.draw_rect(self.total_width + self.x_start - normalized_scalebar_width, 100 + header_offset - self.hight + 8, '#000000', normalized_scalebar_width, 1, 1))
            self.write_to_file('<text text-anchor="middle" font-family="Arial" x="{}" y="{}" font-size="{}" >{}</text>'.format(self.total_width + self.x_start - normalized_scalebar_width / 2, 96 + header_offset - self.hight + 8, font_size_axis_x, scalebar_display))
        if group_labels is not None:
            for group_id, track_y, _ in group_layouts:
                self.write_to_file('<text text-anchor="start" font-family="Arial" x="{}" y="{}" font-size="8" >{}</text>'.format(self.x_start + self.total_width + 15, track_y - int(self.hight / 2) + 2.401, escape(group_labels[group_id - 1])))
        if args.get('_named_tracks') or args['transform'] != 'none':
            self.draw_condition_legend(args, series_by_group)
        elif labels is not None:
            self.write_to_file(self.draw_rect(self.x_start - 10.5, 34, fills[0], 10, 10, self.opacity))
            self.write_to_file('<text text-anchor="start" font-family="Arial" x="' + str(self.x_start + 3) + '" y="' + str(34 - 1.788) + '" font-size="8" >' + escape(str(labels[0])) + '</text>')
            self.write_to_file(self.draw_rect(self.x_start - 10.5, 47, fills[1], 10, 10, self.opacity))
            self.write_to_file('<text text-anchor="start" font-family="Arial" x="' + str(self.x_start + 3) + '" y="' + str(47 - 1.788) + '" font-size="8" >' + escape(str(labels[1])) + '</text>')
            if self.show_plots == 'averages':
                self.write_to_file(self.draw_rect(self.x_start - 10.5, 60, fills[0], 10, 10, self.opacity))
                self.write_to_file('<text text-anchor="start" font-family="Arial" x="' + str(self.x_start + 3) + '" y="' + str(60 - 1.788) + '" font-size="8" >' + 'Overlap' + '</text>')
                self.write_to_file(self.draw_rect(self.x_start - 10.5, 60, fills[1], 10, 10, self.opacity))
            if spark == 'yes':
                self.write_to_file(self.draw_rect(self.x_start + 59.5, 34, self.spark_color[1], 10, 10, 0.75))
                self.write_to_file('<text text-anchor="start" font-family="Arial" x="' + str(self.x_start + 73) + '" y="' + str(34 - 1.788) + '" font-size="8" >' + escape(str(labels[0])) + ' up</text>')
                self.write_to_file(self.draw_rect(self.x_start + 59.5, 47, self.spark_color[0], 10, 10, 0.75))
                self.write_to_file('<text text-anchor="start" font-family="Arial" x="' + str(self.x_start + 73) + '" y="' + str(47 - 1.788) + '" font-size="8" >' + escape(str(labels[1])) + ' up</text>')
        if args['display_chrom_location'] != 'no':
            if args['display_chrom_location'] == 'bottom_left':
                y_value_chr_label = self.y_start + (self.hight if has_negative else 0) + 9
                x_value_chr_label = self.x_start - 10.5
                text_anchor = 'start'
            elif args['display_chrom_location'] == 'top_left':
                y_value_chr_label = (first_track_y_start if first_track_y_start is not None else 100) - self.hight - 4
                x_value_chr_label = self.x_start - 10.5
                text_anchor = 'start'
            elif args['display_chrom_location'] == 'bottom_right':
                y_value_chr_label = self.y_start + (self.hight if has_negative else 0) + 9
                x_value_chr_label = self.x_start + self.total_width
                text_anchor = 'end'
            elif args['display_chrom_location'] == 'top_right':
                y_value_chr_label = (first_track_y_start if first_track_y_start is not None else 100) - self.hight - 4
                x_value_chr_label = self.x_start + self.total_width
                text_anchor = 'end'
            self.write_to_file('<text text-anchor="' + text_anchor + '" font-family="Arial" x="' + str(x_value_chr_label) + '" y="' + str(y_value_chr_label) + '" font-size="' + str(font_size_axis_x) + '" >Chr' + escape(str(self.region[0])) + ': ' + str(f'{self.region[1]:,}') + '-' + str(f'{self.region[2]:,}') + '</text>')
        y_position_bed = 110 + header_offset + (len(active_groups) - 1) * self.hight * 1.5 + additional_hight
        if bed_files is not None:
            for nr_bed, bed_file in enumerate(bed_files):
                with open(bed_file) as b:
                    for line in b:
                        if not line.strip() or line.startswith(('#', 'track ', 'browser ')):
                            continue
                        region_to_draw = [0, 0]
                        self.line_split = line.split('\t')
                        if len(self.line_split) < 3:
                            raise ValueError('BED requires at least three tab-separated columns')
                        if self.line_split[0] != '':
                            if self.line_split[0][:3] == 'chr' or self.line_split[0][:3] == 'Chr':
                                self.line_split[0] = self.line_split[0][3:]
                            else:
                                pass
                        if self.line_split[0] == self.region[0]:
                            self.line_split[1] = int(self.line_split[1])
                            self.line_split[2] = int(self.line_split[2])
                            if self.line_split[1] < 0 or self.line_split[2] <= self.line_split[1]:
                                raise ValueError('Invalid BED interval')
                            if self.line_split[2] > self.region[1] and self.line_split[1] < self.region[2]:
                                region_to_draw = [self.line_split[1], self.line_split[2]]
                                if self.line_split[1] < self.region[1]:
                                    region_to_draw[0] = self.region[1]
                                if self.line_split[2] > self.region[2]:
                                    region_to_draw[1] = self.region[2]
                        if region_to_draw != [0, 0]:
                            self.write_to_file(self.draw_rect(self.x_start + (region_to_draw[0] - self.region[1]) * self.total_width / float(self.region[2] - self.region[1]), y_position_bed - 0.3 + 2 / 2, bed_color[nr_bed], (region_to_draw[1] - region_to_draw[0]) * self.total_width / float(self.region[2] - self.region[1]), 2, 1))
                self.write_to_file('<text text-anchor="start" font-family="Arial" x="' + str(self.x_start + self.total_width + 15) + '" y="' + str(y_position_bed + 3) + '" font-size="8" >' + escape(bed_labels[nr_bed]) + '</text>')
                y_position_bed += 8
        if gff_file is not None:
            plot_transcripts = False
            with contextlib.closing(get_lines_indexed(gff_file, self.region[0])) as f:
                gene_or_transcript_names_to_plot = []
                gene_rows = {}
                gene_row_bottom = y_position_bed
                transcripts = []
                tss_plotted_genes = []
                y_genestart = y_position_bed
                exons = []
                UTRs = []
                for line in f:
                    if not line.strip() or line.startswith('#'):
                        continue
                    self.line_split = line.rstrip().split('\t')
                    if len(self.line_split) != 9:
                        raise ValueError('GTF requires nine tab-separated columns')
                    try:
                        self.line_split[3] = int(self.line_split[3]) - 1
                    except ValueError as exc:
                        raise ValueError('GTF start must be an integer') from exc
                    try:
                        self.line_split[4] = int(self.line_split[4])
                    except ValueError as exc:
                        raise ValueError('GTF end must be an integer') from exc
                    if self.line_split[3] < 0 or self.line_split[4] <= self.line_split[3]:
                        raise ValueError('Invalid GTF interval')
                    try:
                        self.line_split[len(self.line_split) - 1] = self.line_split[len(self.line_split) - 1].split('\n')[0]
                    except:
                        pass
                    try:
                        if self.line_split[0][:3] == 'chr' or self.line_split[0][:3] == 'Chr':
                            self.line_split[0] = self.line_split[0][3:]
                    except:
                        pass
                    if self.line_split[0] == self.region[0]:
                        if display_genes is not None and self.get_gene_name() not in display_genes:
                            continue
                        to_draw = self.get_region_to_draw()
                        if to_draw != 0:
                            if self.line_split[2] in ['gene', 'CDS', 'exon', 'transcript', 'start_codon']:
                                gene_or_transcript_name = ''
                                if display_transcripts[0] == 'mergeall':
                                    if display_genes is None:
                                        gene_or_transcript_name = self.get_gene_name()
                                    elif self.get_gene_name() in display_genes:
                                        gene_or_transcript_name = self.get_gene_name()
                                elif 'transcript_id ' in self.line_split[8]:
                                    if display_transcripts[0] == 'all':
                                        gene_or_transcript_name = self.get_transcript_name()
                                        plot_transcripts = True
                                    elif self.get_transcript_name() in display_transcripts:
                                        gene_or_transcript_name = self.get_transcript_name()
                                        plot_transcripts = True
                                plot_gene = False
                                if gene_or_transcript_name != '':
                                    if gene_or_transcript_name not in gene_or_transcript_names_to_plot:
                                        if display_transcripts[0] != 'mergeall':
                                            if display_transcripts[0] == 'all':
                                                y_genestart += 12
                                                gene_or_transcript_names_to_plot.append(gene_or_transcript_name)
                                            elif gene_or_transcript_name in display_transcripts:
                                                y_genestart += 12
                                                gene_or_transcript_names_to_plot.append(gene_or_transcript_name)
                                        else:
                                            y_genestart += 12
                                            gene_or_transcript_names_to_plot.append(gene_or_transcript_name)
                                    if display_transcripts[0] != 'mergeall':
                                        if gene_or_transcript_name in display_transcripts:
                                            plot_gene = True
                                        elif display_transcripts[0] == 'all':
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
                                        if self.line_split[2] == 'transcript':
                                            draw_gene_or_transcript = True
                                    elif self.line_split[2] == 'gene':
                                        draw_gene_or_transcript = True
                                    if draw_gene_or_transcript == True:
                                        if plot_gene_name_instead_transcriptID == 'yes':
                                            gene_label = self.get_gene_name()
                                        else:
                                            gene_label = gene_or_transcript_name
                                        self.write_to_file('<text text-anchor="start" font-family="Arial" x="' + str(self.x_start + self.total_width + 15) + '" y="' + str(y_genestart + 3) + '" font-size="8" >' + escape(gene_label) + '</text>')
                                        drawhight = 1.0
                                    elif self.line_split[2] == 'CDS':
                                        drawhight = 6.0
                                    elif self.line_split[2] == 'exon':
                                        drawhight = 3.0
                                    if draw_gene_or_transcript and display_genestart == 'yes':
                                        gene_start = self.line_split[3] if self.line_split[6] == '+' else self.line_split[4] - 1
                                        if self.region[1] <= gene_start < self.region[2]:
                                            self.draw_start_arrow(gene_start, self.line_split[6], y_genestart)
                                    if draw_gene_or_transcript == True or self.line_split[2] in ['CDS', 'exon']:
                                        self.write_to_file(self.draw_rect(self.x_start + (to_draw[0] - self.region[1]) * self.total_width / float(self.region[2] - self.region[1]), y_genestart - 0.3 + drawhight / 2, '#0B34FF', (to_draw[1] - to_draw[0]) * self.total_width / float(self.region[2] - self.region[1]), drawhight, 1))
                                    if self.line_split[2] == 'start_codon' and self.region[1] <= self.line_split[3] < self.region[2]:
                                        if draw_TSS == 'yes':
                                            plot = True
                                            if plot_transcripts == False:
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
                                                codon_start = self.line_split[3] if self.line_split[6] == '+' else self.line_split[4] - 1
                                                self.draw_start_arrow(codon_start, self.line_split[6], y_genestart)
                                                plotted_TSS = True
        self.write_to_file('</svg>')
        last_bottom = max((y + (track_hight if negative else 0) for _, y, negative in group_layouts))
        annotation_bottom = max(y_position_bed, locals().get('gene_row_bottom', y_position_bed))
        canvas_height = max(initial_canvas_height, last_bottom + 30, annotation_bottom + 25)
        self.finish_svg(max(320, self.x_start + self.total_width + 120), canvas_height)
        self.log('Saved ' + self.output_filename)
        return self.output_filename

def build_parser():
    parser = argparse.ArgumentParser(description='SparK: genomic signal plots as SVG')
    parser.add_argument('-pt', '--plot_type', help='choices: standard, STD, sine', choices=['standard', 'STD', 'sine'], required=False, type=str, default='standard')
    parser.add_argument('-ps', '--show_plots', help='choices: all, averages', choices=['all', 'averages'], required=False, type=str, default='all')
    parser.add_argument('-pr', '--region', help='zero-based half-open interval, e.g. chr1:10000-20000', required=False, type=str)
    parser.add_argument('-cf', '--control_files', help='separate by space', required=False, nargs='+', type=str)
    parser.add_argument('-tf', '--treat_files', help='separate by space', required=False, nargs='+', type=str, default=[])
    parser.add_argument('-cg', '--control_groups', help='group numbers separate by space', required=False, nargs='+', type=int, default=[])
    parser.add_argument('-tg', '--treat_groups', help='group numbers separate by space', required=False, nargs='+', type=int, default=[])
    parser.add_argument('-gl', '--group_labels', help='set group labels', required=False, nargs='+', type=str)
    parser.add_argument('-l', '--labels', help='set labels for controls and treatment', required=False, nargs='+', type=str)
    parser.add_argument('-gs', '--group_autoscale', help='set to "yes" to autoscale all groups(plots), except the groups excluded in -eg', choices=['yes', 'no'], required=False, type=str)
    parser.add_argument('-es', '--exclude_from_group_autoscale', help='group numbers of groups to be excluded from autoscale', required=False, nargs='+', type=int)
    parser.add_argument('-eg', '--exclude_groups', help='Exclude groups from the analysis', required=False, nargs='+', type=int)
    parser.add_argument('-f', '--fills', help='one color, control/treatment colors, a preset, or more than two per-input colors', required=False, nargs='+', type=str, default=None)
    parser.add_argument('-cs', '--custom_y_axis_scales', help='Enter y-axis values for all groups(plots). All groups need a value. Enter "D" for each group no value should be assigned, e.g. to keep autoscaling functionality for some groups', required=False, nargs='+', type=str)
    parser.add_argument('-dc', '--display_chrom_location', help='set to "no" if you do not want to plot the chromosomal coordinates', choices=['no', 'top_left', 'top_right', 'bottom_left', 'bottom_right'], required=False, type=str, default='top_left')
    parser.add_argument('-gtf', '--gtffile', help='link gtf file for drawing genes here', required=False, type=str)
    parser.add_argument('-tss', '--drawtss', help='set to no to hide translation start-codon arrows (legacy -tss behavior)', choices=['yes', 'no'], required=False, type=str, default='yes')
    parser.add_argument('-genestart', '--draw_genestart', help='set to yes to indicate gene/transcript starts (TSS)', choices=['yes', 'no'], required=False, type=str, default='no')
    parser.add_argument('-sp', '--spark', help='display significant change "yes"', choices=['yes', 'no'], required=False, type=str)
    parser.add_argument('-sc', '--spark_color', help='spark color', required=False, type=str, nargs='+')
    parser.add_argument('-sm', '--smoothen', help='smoothen tracks, int', required=False, type=int)
    parser.add_argument('-o', '--output_name', help='output graph name, str', required=False, type=str)
    parser.add_argument('-bed', '--bed_files', help='bed files to be plotted', required=False, type=str, nargs='+')
    parser.add_argument('-bedcol', '--bed_color', help='colors of bed files in hex', required=False, type=str, nargs='+')
    parser.add_argument('-bedlab', '--bed_labels', help='set labels for bed tracks', required=False, nargs='+', type=str)
    parser.add_argument('-w', '--track_width', help='width of the track, default = 150, int', required=False, type=int, default=150)
    parser.add_argument('-dg', '--display_genes', help='genes to display from the gtf file', nargs='+', required=False, type=str)
    parser.add_argument('--max_points', help='maximum datapoints per plot', required=False, type=int, default=2000)
    parser.add_argument('--point_stat', help='per-bin statistic for downsampling: average or max', choices=['average', 'max'], required=False, type=str, default='average')
    parser.add_argument('--y_scale', help='scale factor for y-axis', required=False, type=float, default=1.0)
    parser.add_argument('--x_scale', help='scale factor for x-axis', required=False, type=float, default=1.0)
    parser.add_argument('-dt', '--display_transcripts', help='display custom transcripts. By default, all transcripts annotated in the gtf file will be merged and displayed as one gene. Alternatively all can be plotted seperatelly by setting this to "all". Further, Transcript IDs can be listed to plot only certain transcripts', nargs='+', required=False, type=str, default=['mergeall'])
    parser.add_argument('-wg', '--write_genenames', help='write genename instead of transcript ID when transcripts are plotted. Set to "yes".', choices=['yes', 'no'], required=False, type=str, default='no')
    parser.add_argument('-scale', '--display_scalebar', help='set to "no" to remove scalebar', choices=['yes', 'no'], required=False, type=str, default='yes')
    parser.add_argument('--track_colors', nargs='+', help='one hex color per input, controls first then treatments; unambiguous for one or two inputs')
    parser.add_argument('--version', action='version', version=SparK_Version)
    parser.add_argument('--config', help='load saved plot settings from JSON')
    parser.add_argument('--save-config', help='save effective settings as reusable JSON')
    parser.add_argument('--track', action='append', help='named condition input: CONDITION=FILE; repeat for replicates')
    parser.add_argument('--gene', help='plot a gene name or gene ID using -gtf')
    parser.add_argument('--genes-file', help='plot every gene in a text file, one per line')
    parser.add_argument('--regions-file', help='batch BED file or text file of chr:start-end regions')
    parser.add_argument('--flank', type=int, default=10000, help='bases on each side of a gene; default 10000')
    parser.add_argument('--batch-shared-scale', action='store_true', help='use consistent scales for corresponding groups across the batch')
    parser.add_argument('--compare', nargs=2, metavar=('TARGET','REFERENCE'), help='condition pair for a transform or sparks')
    parser.add_argument('--transform', choices=['none','difference','ratio','log2ratio'], default='none')
    parser.add_argument('--pseudocount', type=float, default=1.0, help='added to both condition means for ratios; default 1')
    parser.add_argument('--scale-factors', nargs='+', type=float, help='explicit per-input multipliers, controls first then treatments')
    parser.add_argument('--formats', nargs='+', default=['svg','jpg'], help='additional exports: svg jpg png pdf; SVG is always retained')
    parser.add_argument('--no-jpg', action='store_true', help='disable the automatic JPG copy')
    parser.add_argument('--dpi', type=int, default=300, help='raster export resolution; default 300 DPI')
    parser.add_argument('--jpg-quality', type=int, default=95, help='JPEG quality from 1 to 100; default 95')
    parser.add_argument('--hash-inputs', action='store_true', help='include full input SHA-256 hashes in run records (reads each full file once)')
    parser.add_argument('--no-run-record', action='store_true', help='disable the automatic .run.json record')
    parser.add_argument('--no-update-check', action='store_true', help='disable the best-effort background version check')
    parser.add_argument('--gui', action='store_true', help='open the local visual interface')
    parser.add_argument('--data-dir', default='.', help='starting folder for the local interface file picker')
    parser.add_argument('--port', type=int, default=0, help='local interface port; default chooses a free port')
    parser.add_argument('--no-open', action='store_true', help='print the local interface URL without opening a browser')

    return parser

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
    controls, treatments = (args['control_files'], args['treat_files'])
    cg, tg = (args['control_groups'], args['treat_groups'])
    if cg or tg:
        if len(cg) != len(controls) or len(tg) != len(treatments):
            raise ValueError('Provide exactly one group assignment per file in -cg and -tg; no inputs may be omitted')
    elif args['show_plots'] == 'averages' or args['plot_type'] in ('STD', 'sine'):
        cg, tg = ([1] * len(controls), [1] * len(treatments))
    else:
        cg = list(range(1, len(controls) + 1))
        tg = list(range(len(controls) + 1, len(controls) + len(treatments) + 1))
    if any((group <= 0 for group in cg + tg)):
        raise ValueError('Group numbers must be positive integers')
    args['control_groups'], args['treat_groups'] = (cg, tg)
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
        nc, nt = (cg.count(group), tg.count(group))
        if not args.get('_named_tracks') and args['plot_type'] == 'sine' and (not nc or not nt):
            raise ValueError('Sine plots require control and treatment inputs in every plotted group')
        if not args.get('_named_tracks') and args['plot_type'] == 'STD' and any((n == 1 for n in (nc, nt))):
            raise ValueError('STD plots require at least two replicates for each plotted condition')
        if not args.get('_named_tracks') and args['spark'] == 'yes' and (nc < 2 or nt < 2):
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
    if not args.get('_named_tracks') and colors[-1] and (args['show_plots'] == 'averages' or args['plot_type'] == 'STD'):
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
    return (region, active, colors, output)

PALETTE = ['#FF1800', '#005CFF', '#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#666666', '#F0E442']
PATH_OPTIONS = ('control_files', 'treat_files', 'bed_files', 'gtffile', 'regions_file', 'genes_file')


def _absolute(path, base):
    return str((Path(base) / os.path.expanduser(str(path))).resolve())


def load_config(path):
    path = Path(path).expanduser().resolve()
    with path.open(encoding='utf-8') as handle:
        config = json.load(handle)
    if not isinstance(config, dict):
        raise ValueError('A configuration must be a JSON object')
    config = copy.deepcopy(config)
    if 'gtf' in config:
        config['gtffile'] = config.pop('gtf')
    for key in PATH_OPTIONS:
        if config.get(key):
            value = config[key]
            config[key] = [_absolute(p, path.parent) for p in value] if isinstance(value, list) else _absolute(value, path.parent)
    if isinstance(config.get('tracks'), list):
        for record in config['tracks']:
            if isinstance(record, dict):
                for key in ('file', 'path'):
                    if record.get(key):
                        record[key] = _absolute(record[key], path.parent)
    elif isinstance(config.get('tracks'), dict):
        config['tracks'] = {name: [_absolute(p, path.parent) for p in ([paths] if isinstance(paths, str) else paths)]
                            for name, paths in config['tracks'].items()}
    for key in ('output_name', 'output'):
        if config.get(key):
            config[key] = _absolute(config[key], path.parent)
    return config


def default_options():
    return vars(build_parser().parse_args([]))


def normalize_config(config):
    """Validate reusable settings without sharing mutable renderer state."""
    if not isinstance(config, dict):
        raise ValueError('Plot settings must be a dictionary or JSON configuration file')
    config = copy.deepcopy(config)
    for short, full in {'output': 'output_name', 'gtf': 'gtffile'}.items():
        if short in config:
            if full in config and config[full] != config[short]:
                raise ValueError('Use only one of {} and {}'.format(short, full))
            config[full] = config.pop(short)
    args = default_options()
    extra = {'tracks', 'condition_colors', 'condition_labels', 'regions', 'genes'}
    unknown = set(config) - set(args) - extra
    if unknown:
        raise ValueError('Unknown setting(s): ' + ', '.join(sorted(unknown)))
    args.update(config)
    for key in ('condition_colors','condition_labels'):
        if args.get(key) is not None and not isinstance(args[key],dict):
            raise ValueError('{} must be an object keyed by condition name'.format(key))
    for key in ('regions','genes'):
        if args.get(key) is not None and (not isinstance(args[key],list) or not args[key]):
            raise ValueError('{} must be a nonempty list'.format(key))
    for action in build_parser()._actions:
        value = args.get(action.dest)
        if action.choices and value is not None and value not in action.choices:
            raise ValueError('{} must be one of {}'.format(action.dest, ', '.join(action.choices)))
    for key in ('max_points', 'track_width', 'dpi', 'jpg_quality', 'flank', 'port'):
        if not isinstance(args[key], int) or isinstance(args[key], bool):
            raise ValueError('{} must be an integer'.format(key))
    if args['smoothen'] is not None and (not isinstance(args['smoothen'],int) or isinstance(args['smoothen'],bool)):
        raise ValueError('smoothen must be an integer')
    for key in ('no_jpg','hash_inputs','no_run_record','no_update_check','batch_shared_scale','gui','no_open'):
        if not isinstance(args[key],bool):
            raise ValueError('{} must be true or false'.format(key))
    if not 0 <= args['port'] <= 65535:
        raise ValueError('port must be between 0 and 65535')
    for key in ('x_scale', 'y_scale', 'pseudocount'):
        if not isinstance(args[key], (int, float)) or isinstance(args[key],bool) or not math.isfinite(args[key]):
            raise ValueError('{} must be a finite number'.format(key))
    if args['flank'] < 0 or args['pseudocount'] < 0:
        raise ValueError('Flank and pseudocount must not be negative')
    if not 1 <= args['jpg_quality'] <= 100 or not 1 <= args['dpi'] <= 2400:
        raise ValueError('JPG quality must be 1–100 and DPI must be 1–2400')
    for key in ('control_files', 'treat_files', 'bed_files', 'control_groups', 'treat_groups',
                'labels', 'group_labels', 'fills', 'track_colors', 'custom_y_axis_scales',
                'exclude_groups','exclude_from_group_autoscale','bed_color','bed_labels','spark_color',
                'display_transcripts','display_genes','scale_factors','track'):
        if args[key] is not None and not isinstance(args[key], list):
            raise ValueError('{} must be a JSON list'.format(key))
    for key in ('control_groups','treat_groups','exclude_groups','exclude_from_group_autoscale'):
        if args.get(key) is not None and any(not isinstance(value,int) or isinstance(value,bool) for value in args[key]):
            raise ValueError('{} must contain integers'.format(key))
    if args.get('scale_factors') is not None and any(not isinstance(v,(int,float)) or isinstance(v,bool) or not math.isfinite(v) or v<=0 for v in args['scale_factors']):
        raise ValueError('scale_factors must contain finite positive numbers')
    sources = [key for key in ('region','gene','genes_file','regions_file','regions','genes') if args.get(key)]
    if not args.get('gui') and len(sources) != 1:
        raise ValueError('Choose one region source: -pr, --gene, --genes-file, or --regions-file')
    if any(args.get(key) for key in ('gene','genes_file','genes')) and not args.get('gtffile'):
        raise ValueError('Gene plotting requires a GTF file (-gtf)')
    formats = args['formats']
    if isinstance(formats, str):
        formats = formats.split(',')
    if not isinstance(formats, list) or not all(isinstance(f, str) for f in formats):
        raise ValueError('formats must be a list such as ["svg", "jpg", "pdf"]')
    formats = [part.lower() for f in formats for part in f.split(',')]
    formats = ['jpg' if f == 'jpeg' else f for f in formats]
    if set(formats) - {'svg','jpg','png','pdf'}:
        raise ValueError('Export formats are svg, jpg, png and pdf')
    args['formats'] = list(dict.fromkeys(['svg'] + ([] if args['no_jpg'] else ['jpg']) + [f for f in formats if f != 'jpg' or not args['no_jpg']]))
    if args.get('tracks') and (args.get('control_files') or args.get('treat_files') or args.get('track')):
        raise ValueError('Use tracks metadata or legacy file arguments, not both')
    if args.get('track'):
        records = []
        for value in args['track']:
            if '=' not in value:
                raise ValueError('--track uses CONDITION=FILE; repeat it for each input')
            condition, filename = value.split('=', 1)
            records.append({'condition':condition, 'file':filename})
        args['tracks'] = records
    args['_named_tracks'] = bool(args.get('tracks'))
    if args['_named_tracks']:
        if args.get('scale_factors') is not None:
            raise ValueError('Put named-track scale factors in each track entry')
        tracks = args['tracks']
        if isinstance(tracks, dict):
            tracks = [dict(file=p, condition=name) for name, paths in tracks.items()
                      for p in ([paths] if isinstance(paths, str) else paths)]
        if not isinstance(tracks, list) or not tracks:
            raise ValueError('tracks must be a nonempty list or condition-to-files mapping')
        records, conditions, group_names = [], [], []
        condition_colors = args.get('condition_colors') or {}
        for raw in tracks:
            if not isinstance(raw, dict) or set(raw) - {'file','path','condition','group','label','color','scale_factor'}:
                raise ValueError('Each track uses file, condition, group, label, color and/or scale_factor')
            record = copy.deepcopy(raw)
            filename = record.pop('path', None) or record.get('file')
            if not isinstance(filename, (str,os.PathLike)) or not str(filename):
                raise ValueError('Every track requires a file path')
            filename = str(filename)
            condition = record.get('condition', 'Signal')
            if not isinstance(condition, str) or not condition.strip():
                raise ValueError('Every condition requires a nonempty name')
            if condition not in conditions:
                conditions.append(condition)
            group = str(record.get('group', 'Signal'))
            if group not in group_names:
                group_names.append(group)
            factor = record.get('scale_factor', 1.0)
            if not isinstance(factor, (int,float)) or isinstance(factor,bool) or not math.isfinite(factor) or factor <= 0:
                raise ValueError('Track scale_factor must be finite and positive')
            color = record.get('color') or condition_colors.get(condition) or PALETTE[conditions.index(condition) % len(PALETTE)]
            records.append(dict(file=_absolute(filename, Path.cwd()), condition=condition, group=group_names.index(group)+1,
                                label=str(record.get('label',Path(filename).name)), color=hex_color(color), scale_factor=float(factor)))
        args['_conditions_order'] = conditions
        first = conditions[0]
        controls = [r for r in records if r['condition']==first]
        treatments = [r for r in records if r['condition']!=first]
        args['control_files'], args['treat_files'] = [r['file'] for r in controls], [r['file'] for r in treatments]
        args['control_groups'], args['treat_groups'] = [r['group'] for r in controls], [r['group'] for r in treatments]
        if args['fills'] or args['track_colors']:
            raise ValueError('Put named-track colors in tracks or condition_colors')
        args['track_colors'] = [r['color'] for r in controls+treatments]
        args['group_labels'] = args['group_labels'] or group_names
        args['_track_records'] = records
    else:
        if not args.get('control_files'):
            raise ValueError('Provide -cf inputs, --track inputs, or tracks in a configuration')
        args['control_files'] = [_absolute(p, Path.cwd()) for p in args['control_files']]
        args['treat_files'] = [_absolute(p, Path.cwd()) for p in args['treat_files']]
        args['_conditions_order'] = ['Control'] + (['Treatment'] if args['treat_files'] else [])
    conditions = args['_conditions_order']
    if args.get('compare'):
        if not isinstance(args['compare'], (list,tuple)) or len(args['compare']) != 2:
            raise ValueError('compare requires [target, reference] condition names')
        target, reference = args['compare']
        if target not in conditions or reference not in conditions or target==reference:
            raise ValueError('Compare conditions must be distinct names present in the inputs')
    else:
        if len(conditions)>2 and (args['transform']!='none' or args['spark']=='yes'):
            raise ValueError('Choose --compare TARGET REFERENCE for comparisons involving more than two conditions')
        reference = conditions[0]
        target = conditions[1] if len(conditions)>1 else None
    args['_compare_reference'], args['_compare_target'] = reference, target
    if args['transform']!='none':
        if not target:
            raise ValueError('Difference and ratio tracks require two conditions')
        if args['plot_type']!='standard' or args['spark']=='yes':
            raise ValueError('Difference and ratio tracks use standard plots without sparks; they compare condition means')
    if args['plot_type']=='sine' and len(conditions)!=2:
        raise ValueError('Sine plots require exactly two conditions; use standard or STD for more conditions')
    for key in ('gtffile','regions_file','genes_file'):
        if args.get(key):
            args[key] = _absolute(args[key], Path.cwd())
    if args.get('bed_files'):
        if args.get('bed_labels') is None:
            args['bed_labels'] = list(args['bed_files'])
        args['bed_files'] = [_absolute(p, Path.cwd()) for p in args['bed_files']]
    output = str(args.get('output_name') or ('SparK_batch' if any(args.get(key) for key in ('genes_file','regions_file','genes','regions')) else 'SparK_graph'))
    if Path(output).suffix.lower() in {'.svg','.jpg','.jpeg','.png','.pdf'}:
        output = str(Path(output).with_suffix(''))
    args['output_name'] = _absolute(output, Path.cwd())
    return args


def finalize_arguments(args):
    """Finish validation once a gene/batch request has been resolved to coordinates."""
    region, active, colors, output = validate_arguments(args)
    if not args['_named_tracks']:
        factors = args.get('scale_factors') or [1.0]*len(args['control_files']+args['treat_files'])
        if len(factors)!=len(args['control_files']+args['treat_files']) or any(not math.isfinite(v) or v<=0 for v in factors):
            raise ValueError('Provide one finite positive scale factor per input')
        args['_track_records'] = []
        index = 0
        for name, files, groups, palette in [('Control',args['control_files'],args['control_groups'],colors[2]),
                                             ('Treatment',args['treat_files'],args['treat_groups'],colors[3])]:
            for filename, group, color in zip(files,groups,palette):
                args['_track_records'].append(dict(file=filename,condition=name,group=group,color=color,
                                                   label=Path(filename).name,scale_factor=float(factors[index])))
                index += 1
    for group in active:
        counts = {name:sum(r['group']==group and r['condition']==name for r in args['_track_records'])
                  for name in args['_conditions_order']}
        if args['plot_type']=='sine' and any(counts[name]==0 for name in args['_conditions_order']):
            raise ValueError('Sine plots require both conditions in every plotted group')
        if args['plot_type']=='STD' and any(n==1 for n in counts.values()):
            raise ValueError('STD plots require at least two replicates for every plotted condition')
        if args['spark']=='yes' and any(counts.get(name,0)<2 for name in (args['_compare_reference'],args['_compare_target'])):
            raise ValueError('Sparks require two replicates of both comparison conditions in every group')
        if args['transform']!='none' and any(counts.get(name,0)==0 for name in (args['_compare_reference'],args['_compare_target'])):
            raise ValueError('Each transformed group requires both comparison conditions')
        if args['show_plots']=='averages' or args['plot_type']=='STD':
            for name in counts:
                palette = {r['color'] for r in args['_track_records'] if r['group']==group and r['condition']==name}
                if len(palette)>1:
                    raise ValueError('Replicates averaged into one trace must share a color')
    return region, active, colors, output


def build_group_series(args, loaded):
    groups = {}
    for record in args['_track_records']:
        group = record['group']
        if group in (args.get('exclude_groups') or []):
            continue
        entries = groups.setdefault(group, [])
        entry = next((e for e in entries if e['name']==record['condition']), None)
        if entry is None:
            entry = dict(name=record['condition'], data=[], colors=[])
            entries.append(entry)
        data = loaded[record['file']]
        entry['data'].append(data if record['scale_factor']==1 else [v*record['scale_factor'] for v in data])
        entry['colors'].append(record['color'])
    for entries in groups.values():
        entries.sort(key=lambda entry:args['_conditions_order'].index(entry['name']))
    if args['transform']!='none':
        for group, entries in groups.items():
            ref, target = comparison_data(args, entries)
            left, right = np.mean(target,axis=0), np.mean(ref,axis=0)
            if args['transform']=='difference':
                values = left-right
                label = args['_compare_target']+' − '+args['_compare_reference']
            else:
                numerator, denominator = left+args['pseudocount'], right+args['pseudocount']
                if np.any(numerator<0) or np.any(denominator<=0) or (args['transform']=='log2ratio' and np.any(numerator<=0)):
                    raise ValueError('Ratio inputs plus pseudocount must be nonnegative (strictly positive for log2); increase the pseudocount or use difference')
                values = numerator/denominator
                if args['transform']=='log2ratio':
                    values = np.log2(values)
                label = ('log₂ ' if args['transform']=='log2ratio' else '')+args['_compare_target']+' / '+args['_compare_reference']
            if not np.all(np.isfinite(values)):
                raise ValueError('The requested transformation produced non-finite values')
            groups[group] = [dict(name=label,data=[values.tolist()],colors=['#005CFF'])]
    return groups


def comparison_data(args, series):
    reference = next((entry['data'] for entry in series if entry['name']==args['_compare_reference']), [])
    target = next((entry['data'] for entry in series if entry['name']==args['_compare_target']), [])
    return reference, target


def display_conditions(args, groups):
    names = []
    labels = args.get('condition_labels') or {}
    for entries in groups.values():
        for entry in entries:
            if entry['name'] not in [name for name,color in names]:
                names.append((entry['name'],entry['colors'][0]))
    return [(str(labels.get(name,name)),color) for name,color in names]


def reusable_config(args):
    ignored = {'config','save_config','gui','data_dir','port','no_open','track'}
    result = {key:copy.deepcopy(value) for key,value in args.items()
              if not key.startswith('_') and key not in ignored and value is not None}
    if args['_named_tracks']:
        result['tracks'] = copy.deepcopy(args['_track_records'])
        for record in result['tracks']:
            record['group'] = (args.get('group_labels') or [str(i) for i in range(1,max(r['group'] for r in result['tracks'])+1)])[record['group']-1]
        for key in ('control_files','treat_files','control_groups','treat_groups','track_colors','fills','scale_factors'):
            result.pop(key,None)
    return result


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w',encoding='utf-8',dir=path.parent,prefix='.spark-',suffix='.json',delete=False) as handle:
            temporary = Path(handle.name)
            json.dump(value,handle,indent=2,sort_keys=True,ensure_ascii=False,allow_nan=False)
            handle.write('\n')
        os.replace(temporary,path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def require_exports(formats):
    try:
        if set(formats) & {'jpg','png'}:
            import resvg_py
            from PIL import Image
        if 'pdf' in formats:
            from svglib.svglib import svg2rlg
            from reportlab.graphics import renderPDF
    except ImportError as exc:
        raise ImportError('Image export needs the updated dependencies. Run: python -m pip install -r requirements.txt (or install Pillow resvg-py reportlab svglib)') from exc


def export_figures(svg, formats, dpi=300, quality=95):
    """Render JPG/PNG on white and retain vector geometry in PDF exports."""
    paths = {'svg':Path(svg)}
    if set(formats) & {'jpg','png'}:
        import resvg_py
        from PIL import Image
        import xml.etree.ElementTree as ET
        geometry = ET.parse(svg).getroot().get('viewBox').split()
        width, height = float(geometry[2]), float(geometry[3])
        if width*height*(dpi/96)**2 > 100_000_000:
            raise ValueError('Raster output would exceed 100 million pixels; reduce --dpi or use --no-jpg --formats svg/pdf')
        raster = resvg_py.svg_to_bytes(svg_path=str(svg),zoom=dpi/96,background='white')
        with Image.open(io.BytesIO(raster)) as original:
            image = original.convert('RGB')
            for ext in ('jpg','png'):
                if ext in formats:
                    destination = Path(svg).with_suffix('.'+ext)
                    if ext=='jpg':
                        image.save(destination,format='JPEG',quality=quality,subsampling=0,dpi=(dpi,dpi))
                    else:
                        image.save(destination,format='PNG',dpi=(dpi,dpi))
                    paths[ext] = destination
    if 'pdf' in formats:
        from svglib.svglib import svg2rlg
        from reportlab.graphics import renderPDF
        from reportlab.pdfgen import canvas
        drawing = svg2rlg(str(svg))
        if drawing is None:
            raise ValueError('Could not convert the SVG to vector PDF')
        destination = Path(svg).with_suffix('.pdf')
        document = canvas.Canvas(str(destination),pagesize=(drawing.width,drawing.height),invariant=1)
        document.setTitle(Path(svg).stem)
        document.setCreator('SparK '+SparK_Version)
        renderPDF.draw(drawing,document,0,0)
        document.showPage()
        document.save()
        paths['pdf'] = destination
    return paths


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda:handle.read(8*1024*1024),b''):
            digest.update(chunk)
    return digest.hexdigest()


def input_fingerprint(path, full_hash=False):
    path = Path(path).resolve()
    info = path.stat()
    record = dict(path=str(path),size=info.st_size,mtime_ns=info.st_mtime_ns,ctime_ns=info.st_ctime_ns,inode=info.st_ino,
                  verification='sha256' if full_hash else 'metadata')
    if full_hash:
        signature = [str(path),info.st_size,info.st_mtime_ns,info.st_ctime_ns,info.st_ino]
        key = hashlib.sha256(json.dumps(signature).encode()).hexdigest()
        cache = Path(tempfile.gettempdir())/('spark-hashes-v1-'+str(getattr(os,'getuid',lambda:0)()))/(key+'.json')
        try:
            cached = json.loads(cache.read_text())
            if cached.get('signature') != signature:
                raise ValueError('Stale hash cache')
            record['sha256'] = cached['sha256']
        except (OSError,ValueError,KeyError):
            record['sha256'] = file_sha256(path)
            after = path.stat()
            if (after.st_size,after.st_mtime_ns,after.st_ctime_ns,after.st_ino)!=tuple(signature[1:]):
                raise OSError('Input changed while hashing: '+str(path))
            try:
                cache.parent.mkdir(mode=0o700,parents=True,exist_ok=True)
                atomic_json(cache,dict(signature=signature,sha256=record['sha256']))
            except OSError:
                pass
    return record


def run_record(args, staged, destinations):
    sources = list(dict.fromkeys(args['control_files']+args['treat_files']+(args.get('bed_files') or [])+([args['gtffile']] if args.get('gtffile') else [])))
    packages = {}
    for name in ('numpy','pyBigWig','resvg-py','Pillow','reportlab','svglib'):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            pass
    return dict(schema_version=1,spark_version=SparK_Version,created_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                python=sys.version,packages=packages,code_sha256=file_sha256(__file__),settings=reusable_config(args),
                resolved_region=args['region'],gene=args.get('_resolved_gene'),
                inputs=[input_fingerprint(path,args['hash_inputs']) for path in sources],
                outputs={name:dict(path=str(destinations[name]),size=Path(path).stat().st_size,sha256=file_sha256(path)) for name,path in staged.items()})


def publish_outputs(staged, destinations, work_dir):
    """Stage all formats first; restore prior files if publication fails."""
    backups, written = {}, []
    try:
        for name, destination in destinations.items():
            destination = Path(destination)
            if destination.exists():
                if not destination.is_file():
                    raise ValueError('Output destination is not a file: '+str(destination))
                backup = Path(work_dir)/('previous-'+destination.name)
                os.replace(destination,backup)
                backups[destination] = backup
            os.replace(staged[name],destination)
            written.append(destination)
    except Exception:
        for destination in written:
            if destination.exists():
                destination.unlink()
        for destination,backup in backups.items():
            os.replace(backup,destination)
        raise


VERSION_URL = 'https://api.github.com/repos/StefanKurtenbach/SparK/releases/latest'
RELEASE_URL = 'https://github.com/StefanKurtenbach/SparK/releases/latest'


def version_tuple(value):
    match = re.fullmatch(r'v?(\d+)\.(\d+)\.(\d+)', str(value).strip())
    return tuple(map(int,match.groups())) if match else None


def _update_cache_path():
    return Path(tempfile.gettempdir())/('spark-update-v1-'+str(getattr(os,'getuid',lambda:0)()))/'latest.json'


def _cached_update(path):
    try:
        value = json.loads(Path(path).read_text())
        age = time.time()-float(value['checked_at'])
        ttl = 86400 if version_tuple(value.get('latest_version')) else 600
        return value if 0 <= age < ttl else None
    except (OSError,ValueError,KeyError,TypeError):
        return None


def fetch_latest_version():
    """Check the latest published stable release; never send plot/input information."""
    request = urllib.request.Request(VERSION_URL,headers={'User-Agent':'SparK/'+SparK_Version,'Accept':'application/vnd.github+json'})
    with urllib.request.urlopen(request,timeout=1.5) as response:
        release = json.loads(response.read(65536).decode('utf-8'))
    if release.get('draft') or release.get('prerelease'):
        raise ValueError('Only published stable releases are checked')
    value = release.get('tag_name','')
    if not version_tuple(value):
        raise ValueError('The version response is not a stable semantic version')
    return value.lstrip('v')


def _print_update(latest, current, stream):
    newer, installed = version_tuple(latest), version_tuple(current)
    if newer and installed and newer > installed:
        print('\nSparK {} is available (installed: {}).\n{}'.format(latest,current,RELEASE_URL),file=stream,flush=True)


def _claim_update_notice(work):
    try:
        descriptor = os.open(Path(work)/'notified',os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o600)
        os.close(descriptor)
        return True
    except OSError:
        return False


def update_worker(cache_path, work_dir, current, terminal):
    """A short-lived worker that outlives a fast plot without holding its pipes open."""
    import shutil
    cache_path, work = Path(cache_path), Path(work_dir)
    latest = None
    try:
        latest = fetch_latest_version()
    except Exception:
        pass  # Update availability must never affect plotting or exit status.
    try:
        cache_path.parent.mkdir(mode=0o700,parents=True,exist_ok=True)
        atomic_json(cache_path,dict(checked_at=time.time(),latest_version=latest))
        atomic_json(work/'result.json',dict(latest_version=latest))
        marker = work/'done.json'
        if marker.exists():
            success = json.loads(marker.read_text()).get('success',False)
            if terminal and success and _claim_update_notice(work):
                _print_update(latest,current,sys.stderr)
            shutil.rmtree(work,ignore_errors=True)
    except Exception:
        pass
    return 0


class UpdateCheck:
    """Best-effort check: no prompts, automatic installs, joins or plotting delays."""
    def __init__(self, enabled=True, stream=None):
        self.stream = stream or sys.stderr
        self.cached = None
        self.work = None
        self.finished = False
        if not enabled or os.environ.get('SPARK_NO_UPDATE_CHECK')=='1':
            return
        self.cached = _cached_update(_update_cache_path())
        if self.cached:
            return
        try:
            self.work = Path(tempfile.mkdtemp(prefix='spark-version-check-'))
            terminal = bool(getattr(self.stream,'isatty',lambda:False)())
            subprocess.Popen([sys.executable,str(Path(__file__).resolve()),'--_update-worker',
                              str(_update_cache_path()),str(self.work),SparK_Version,'1' if terminal else '0'],
                             stdin=subprocess.DEVNULL,stdout=subprocess.DEVNULL,
                             stderr=self.stream if terminal else subprocess.DEVNULL,
                             close_fds=True,start_new_session=(os.name=='posix'))
        except Exception:
            if self.work:
                import shutil
                shutil.rmtree(self.work,ignore_errors=True)
            self.work = None

    def finish(self, success=True):
        if self.finished:
            return
        self.finished = True
        if self.cached:
            if success:
                _print_update(self.cached.get('latest_version'),SparK_Version,self.stream)
            return
        if self.work is None:
            return
        import shutil
        try:
            atomic_json(self.work/'done.json',dict(success=success))
            result = self.work/'result.json'
            if result.exists():
                latest = json.loads(result.read_text()).get('latest_version')
                if success and _claim_update_notice(self.work):
                    _print_update(latest,SparK_Version,self.stream)
                shutil.rmtree(self.work,ignore_errors=True)
        except Exception:
            pass


def gene_index(gtf):
    """Cache gene coordinates, indexing only gene records in plain-text GTFs."""
    gtf = Path(gtf).resolve()
    stat = gtf.stat()
    signature = [str(gtf),stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns,stat.st_ino,1]
    key = hashlib.sha256(json.dumps(signature).encode()).hexdigest()
    cache = Path(tempfile.gettempdir())/('spark-genes-v1-'+str(getattr(os,'getuid',lambda:0)()))/(key+'.json')
    try:
        saved = json.loads(cache.read_text())
        if saved['signature']==signature:
            return saved['genes']
    except (OSError,ValueError,KeyError,TypeError):
        pass
    genes = {}

    def add(line, fallback=False):
        fields = line.rstrip().split('\t')
        if len(fields)!=9 or fields[2] not in (('gene','transcript','exon') if fallback else ('gene',)):
            return
        attributes = dict(re.findall(r'(?:^|;)\s*(\S+)\s+"([^"]*)"',fields[8]))
        name, ident = attributes.get('gene_name'), attributes.get('gene_id')
        if not name and not ident:
            return
        start, stop = int(fields[3])-1,int(fields[4])
        if start<0 or stop<=start:
            raise ValueError('Invalid gene coordinates in '+str(gtf))
        record = dict(chromosome=fields[0],start=start,end=stop,strand=fields[6],gene_name=name or ident,gene_id=ident or name)
        for alias in set(filter(None,(name,ident))):
            values = genes.setdefault(alias,[])
            existing = next((v for v in values if v['gene_id']==record['gene_id'] and v['chromosome']==record['chromosome']),None)
            if existing:
                existing['start'],existing['end'] = min(existing['start'],start),max(existing['end'],stop)
            else:
                values.append(dict(record))
    if gtf.name.endswith('.gz'):
        with gzip.open(gtf,'rt') as handle:
            for line in handle:
                if not line.startswith('#'):
                    add(line)
    elif stat.st_size:
        with gtf.open('rb') as handle, mmap.mmap(handle.fileno(),0,access=mmap.ACCESS_READ) as mapped:
            if re.search(rb'\r(?!\n)',mapped):
                with gtf.open() as lines:
                    for line in lines:
                        if not line.startswith('#'):
                            add(line)
            else:
                position = mapped.find(b'\tgene\t')
                while position>=0:
                    start = mapped.rfind(b'\n',0,position)+1
                    stop = mapped.find(b'\n',position)
                    stop = len(mapped) if stop<0 else stop+1
                    add(mapped[start:stop].decode('utf-8'))
                    position = mapped.find(b'\tgene\t',stop)
    if not genes:
        with (gzip.open(gtf,'rt') if gtf.name.endswith('.gz') else gtf.open()) as handle:
            for line in handle:
                if not line.startswith('#'):
                    add(line,fallback=True)
    after = gtf.stat()
    if (after.st_size,after.st_mtime_ns,after.st_ctime_ns,after.st_ino)!=tuple(signature[1:5]):
        raise OSError('The GTF changed while the gene index was being built')
    try:
        cache.parent.mkdir(mode=0o700,parents=True,exist_ok=True)
        atomic_json(cache,dict(signature=signature,genes=genes))
    except OSError:
        pass
    return genes


def resolve_gene(name, gtf, flank=10000, index=None):
    matches = (index if index is not None else gene_index(gtf)).get(name,[])
    if not matches:
        raise ValueError('Gene {!r} was not found in the GTF; use its exact name or gene ID'.format(name))
    if len(matches)!=1:
        raise ValueError('Gene {!r} has multiple loci; use a unique gene ID or explicit coordinates'.format(name))
    record = dict(matches[0])
    record['query'] = name
    start, stop = max(0,record['start']-flank),record['end']+flank
    return '{}:{}-{}'.format(record['chromosome'],start,stop),record


def read_regions(path):
    jobs = []
    with Path(path).open() as handle:
        for number,line in enumerate(handle,1):
            line = line.strip()
            if not line or line.startswith(('#','track ','browser ')):
                continue
            fields = line.split()
            if ':' in fields[0]:
                region = fields[0]
                name = fields[1] if len(fields)>1 else region
            elif len(fields)>=3:
                region = '{}:{}-{}'.format(fields[0],fields[1],fields[2])
                name = fields[3] if len(fields)>3 else region
            else:
                raise ValueError('Invalid region on line {} of {}'.format(number,path))
            parse_region(region)
            jobs.append(dict(name=name,region=region))
    if not jobs:
        raise ValueError('The regions file contains no regions')
    return jobs


def _read_genes(path):
    with Path(path).open() as handle:
        names = [line.strip() for line in handle if line.strip() and not line.lstrip().startswith('#')]
    if not names:
        raise ValueError('The gene list is empty')
    return names


def _single_args(args, job):
    current = copy.deepcopy(args)
    current.update(region=job['region'],gene=None,genes_file=None,regions_file=None,regions=None,genes=None,save_config=None)
    if job.get('gene'):
        current['_resolved_gene'] = job['gene']
    return current


def _plot_one(args, quiet=True):
    finalize_arguments(args)
    require_exports(args['formats'])
    base = Path(args['output_name'])
    base.parent.mkdir(parents=True,exist_ok=True)
    input_paths = {Path(path).resolve() for path in args['control_files']+args['treat_files']+(args.get('bed_files') or [])+([args['gtffile']] if args.get('gtffile') else [])}
    names = {ext:Path(str(base)+'.'+ext) for ext in args['formats']}
    if any(path.resolve() in input_paths for path in names.values()):
        raise ValueError('An output must not overwrite an input')
    with tempfile.TemporaryDirectory(prefix='.spark-output-',dir=base.parent) as directory:
        staged_args = copy.deepcopy(args)
        staged_args['output_name'] = str(Path(directory)/(base.name+'.svg'))
        renderer = PlotRenderer(quiet=True)
        try:
            svg = renderer.render(staged_args)
        finally:
            renderer.discard_incomplete_svg()
        staged = export_figures(svg,args['formats'],args['dpi'],args['jpg_quality'])
        if not args['no_run_record']:
            staged['record'] = Path(directory)/(base.name+'.run.json')
            record = run_record(args,{key:value for key,value in staged.items() if key!='record'},names)
            atomic_json(staged['record'],record)
            names['record'] = Path(str(base)+'.run.json')
        publish_outputs(staged,names,directory)
    result = dict(version=SparK_Version,region=args['region'],outputs={key:str(path) for key,path in names.items()})
    if args.get('_resolved_gene'):
        result['gene'] = args['_resolved_gene']
    if not quiet:
        for name,path in result['outputs'].items():
            print('Saved '+path,flush=True)
    return result


def _make_config(config, options):
    value = load_config(config) if isinstance(config,(str,os.PathLike)) else copy.deepcopy(config or {})
    selectors = ('region','gene','regions_file','genes_file','regions','genes')
    for selected in selectors:
        if options.get(selected) is not None:
            for other in selectors:
                if other!=selected and other not in options:
                    value.pop(other,None)
    value.update(options)
    return value


def plot(region=None, tracks=None, output=None, config=None, quiet=True, **options):
    """Create a figure from explicit tracks or saved settings; return its output paths.

    Example: plot('chr1:0-1000', tracks={'Control':['a.bw'], 'Drug':['b.bw']}, output='figure')
    No update checks or browser windows are started by this Python API.
    """
    if region is not None:
        options['region'] = region
    if tracks is not None:
        options['tracks'] = tracks
    if output is not None:
        options['output_name'] = output
    args = normalize_config(_make_config(config,options))
    if any(args.get(key) for key in ('genes_file','regions_file','regions','genes')):
        return _plot_batch(args,quiet=quiet)
    original = copy.deepcopy(args)
    if args.get('gene'):
        resolved,record = resolve_gene(args['gene'],args['gtffile'],args['flank'])
        args = _single_args(args,dict(region=resolved,gene=record))
    result = _plot_one(args,quiet)
    if original.get('save_config'):
        atomic_json(original['save_config'],reusable_config(original))
    return result


def plot_many(regions=None, genes=None, tracks=None, output=None, config=None, quiet=True, **options):
    """Create a batch with an HTML index and per-figure run records."""
    if regions is not None and genes is not None:
        raise ValueError('Choose regions or genes, not both')
    if tracks is not None:
        options['tracks'] = tracks
    if output is not None:
        options['output_name'] = output
    config_value = _make_config(config,options)
    if regions is not None or genes is not None:
        for key in ('region','gene','regions_file','genes_file','regions','genes'):
            config_value.pop(key,None)
        if regions is not None:
            config_value['regions'] = [regions] if isinstance(regions,str) else list(regions)
        else:
            config_value['genes'] = [genes] if isinstance(genes,str) else list(genes)
    args = normalize_config(config_value)
    return _plot_batch(args,quiet=quiet)


def _batch_jobs(args, regions=None, genes=None):
    regions = args.get('regions') if regions is None else regions
    genes = args.get('genes') if genes is None else genes
    if regions is None and args.get('regions_file'):
        return read_regions(args['regions_file'])
    if regions is not None:
        jobs = []
        for item in regions:
            job = dict(item) if isinstance(item,dict) else dict(name=str(item),region=str(item))
            parse_region(job['region'])
            jobs.append(job)
        if not jobs:
            raise ValueError('The region list is empty')
        return jobs
    if genes is None and not args.get('genes_file'):
        raise ValueError('Supply a regions file, gene file, or a list of regions/genes for batch plotting')
    genes = list(genes) if genes is not None else _read_genes(args['genes_file'])
    if not args.get('gtffile'):
        raise ValueError('Gene batches require a GTF file')
    index = gene_index(args['gtffile'])
    jobs = []
    for name in genes:
        try:
            region,record = resolve_gene(str(name),args['gtffile'],args['flank'],index)
            jobs.append(dict(name=str(name),region=region,gene=record))
        except (ValueError,OSError) as exc:
            jobs.append(dict(name=str(name),error=str(exc)))
    return jobs


def _batch_scales(args, jobs):
    scales = {}
    for job in jobs:
        if job.get('error'):
            continue
        current = _single_args(args,job)
        try:
            region,active,_,_ = finalize_arguments(current)
            loaded = {filename:read_signals(region,[filename],0,current['max_points'],current['point_stat'],current['smoothen'])[0]
                      for filename in dict.fromkeys(current['control_files']+current['treat_files'])}
            groups = build_group_series(current,loaded)
            renderer = PlotRenderer()
            renderer.plot_type,renderer.show_plots = current['plot_type'],current['show_plots']
            for group in active:
                maximum = renderer.get_max_value(*[entry['data'] for entry in groups[group]])
                scales[group] = max(scales.get(group,0),maximum)
        except (OSError,ValueError,ImportError,RuntimeError) as exc:
            job['error'] = str(exc)
    return scales


def _plot_batch(args, regions=None, genes=None, quiet=True):
    require_exports(args['formats'])
    jobs = _batch_jobs(args,regions,genes)
    directory = Path(args['output_name'])
    directory.mkdir(parents=True,exist_ok=True)
    shared = _batch_scales(args,jobs) if args['batch_shared_scale'] else {}
    results,seen = [],set()
    for number,job in enumerate(jobs,1):
        name = str(job.get('name') or job.get('region') or number)
        if job.get('error'):
            results.append(dict(name=name,error=job['error']))
            continue
        stem = re.sub(r'[^A-Za-z0-9_.-]+','_',name).strip('._') or 'region'
        candidate = stem
        suffix = 2
        while candidate in seen:
            candidate = '{}_{}'.format(stem,suffix)
            suffix += 1
        seen.add(candidate)
        current = _single_args(args,job)
        current['output_name'] = str(directory/candidate)
        if shared:
            finalize_arguments(current)
            maximum_group = max(current['control_groups']+current['treat_groups'])
            scales = list(current.get('custom_y_axis_scales') or ['D']*maximum_group)
            for group,value in shared.items():
                if scales[group-1]=='D' and value>0:
                    scales[group-1] = str(value)
            current['custom_y_axis_scales'] = scales
        try:
            result = _plot_one(current,quiet)
            result['name'] = name
            results.append(result)
        except (OSError,ValueError,ImportError,RuntimeError) as exc:
            results.append(dict(name=name,error=str(exc)))
            if not quiet:
                print('Failed {}: {}'.format(name,exc),file=sys.stderr,flush=True)
    passed = sum('error' not in item for item in results)
    summary = dict(version=SparK_Version,total=len(results),passed=passed,failed=len(results)-passed,plots=results,
                   shared_scales=shared,index=str(directory/'index.html'),record=str(directory/'batch.run.json'))
    atomic_json(directory/'batch.run.json',summary)
    cards = []
    for item in results:
        if item.get('error'):
            cards.append('<article><h2>{}</h2><p class="error">{}</p></article>'.format(escape(item['name']),escape(item['error'])))
            continue
        outputs = item['outputs']
        preview = Path(outputs.get('jpg') or outputs.get('png') or outputs['svg']).name
        links = ' · '.join('<a href="{}">{}</a>'.format(escape(Path(path).name,quote=True),escape(ext.upper())) for ext,path in outputs.items())
        cards.append('<article><h2>{}</h2><p>{}</p><a href="{}"><img src="{}" alt="{}"></a><p>{}</p></article>'.format(
            escape(item['name']),escape(item['region']),escape(Path(outputs['svg']).name,quote=True),escape(preview,quote=True),escape(item['name'],quote=True),links))
    page = '<!doctype html><meta charset="utf-8"><title>SparK batch</title><style>body{font:16px system-ui;max-width:1200px;margin:32px auto;padding:20px;background:#f5f7fa}main{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:20px}article{background:white;border:1px solid #d8dfeb;padding:18px;border-radius:10px}img{max-width:100%;max-height:650px}h2{font-size:18px}.error{color:#a32738}</style>'
    page += '<h1>SparK {} batch: {} / {} plots</h1><main>{}</main>'.format(SparK_Version,passed,len(results),''.join(cards))
    destination = directory/'index.html'
    with tempfile.NamedTemporaryFile(mode='w',encoding='utf-8',dir=directory,delete=False) as handle:
        handle.write(page)
        temporary = handle.name
    os.replace(temporary,destination)
    if args.get('save_config'):
        atomic_json(args['save_config'],reusable_config(args))
    if not quiet:
        print('Batch: {}/{} plots. {}'.format(passed,len(results),destination),flush=True)
    return summary


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1]==['--_update-worker']:
        if len(argv)!=5:
            return 0
        return update_worker(argv[1],argv[2],argv[3],argv[4]=='1')
    parser = build_parser()
    parsed = vars(parser.parse_args(argv))
    options = {flag:action.dest for action in parser._actions for flag in action.option_strings}
    provided = {options[word.split('=',1)[0]] for word in argv if word.split('=',1)[0] in options}
    try:
        config = load_config(parsed['config']) if parsed.get('config') else {}
    except (OSError,ValueError) as exc:
        parser.exit(2,'SparK error: {}\n'.format(exc))
    config.update({key:value for key,value in parsed.items() if key in provided and key not in {'config'}})
    for selected in ('region','gene','genes_file','regions_file'):
        if selected in provided:
            for other in ('region','gene','genes_file','regions_file','regions','genes'):
                if other!=selected and other not in provided:
                    config.pop(other,None)
    checker = UpdateCheck(enabled=not config.get('no_update_check',False))
    success = False
    try:
        if config.get('gui'):
            serve_gui(config,checker)
            success = True
            return 0
        print('SparK Version '+SparK_Version,flush=True)
        result = plot(config=config,quiet=False)
        success = not result.get('failed',0)
        return 0 if success else 1
    except (ValueError,OSError,ImportError,RuntimeError,subprocess.SubprocessError) as exc:
        parser.exit(2,'SparK error: {}\n'.format(exc))
    finally:
        checker.finish(success)


GUI_HTML = r'''<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>SparK — local plot builder</title>
<style>
:root{font:15px system-ui;color:#172536;background:#f1f4f8}*{box-sizing:border-box}body{margin:0}header{padding:22px 28px;background:#162c43;color:white;display:flex;justify-content:space-between;align-items:center}h1{font-size:25px;margin:0}header p{margin:5px 0 0;color:#c5d3e3}main{display:grid;grid-template-columns:440px 1fr;gap:22px;max-width:1500px;margin:24px auto;padding:0 24px}section{background:white;border:1px solid #d8e0ea;border-radius:12px;padding:22px;margin-bottom:18px}h2{font-size:17px;margin:0 0 14px}label{display:block;font-weight:600;margin:12px 0 5px}input,select,textarea,button{font:inherit}input,select,textarea{width:100%;padding:9px;border:1px solid #b9c6d6;border-radius:6px;background:#fff}button{border:1px solid #a6b7ca;border-radius:6px;background:white;color:#172536;padding:9px 13px;cursor:pointer}button:hover{background:#e9f0f8}button.primary{background:#136a9d;color:white;border-color:#136a9d;font-weight:650}button:disabled{opacity:.5;cursor:wait}.row{display:flex;gap:9px;align-items:end}.row>div{flex:1}.track{background:#f6f8fb;border:1px solid #e0e6ee;padding:12px;border-radius:8px;margin-bottom:10px}.track .top{display:grid;grid-template-columns:1fr 42px 32px;gap:8px}.track input[type=color]{height:39px;padding:3px}.track .row{margin-top:8px}.track button.remove{padding:6px}.muted{color:#63748a;font-size:13px;line-height:1.5}.formats{display:flex;gap:14px;align-items:center;margin-top:12px}.formats label{font-weight:400;margin:0;display:flex;align-items:center;gap:5px}.formats input{width:auto}.actions{display:flex;gap:8px;flex-wrap:wrap;margin-top:16px}#status{white-space:pre-wrap;color:#344a64;line-height:1.5}#status.error{color:#a22638}#preview{max-width:100%;max-height:1000px;display:block;margin:auto;border:1px solid #e1e6ed}#empty{padding:90px 20px;text-align:center;color:#63748a;background:#f7f9fc;border:1px dashed #c1cddd;border-radius:8px}.downloads{display:flex;gap:10px;flex-wrap:wrap;margin-top:16px}.downloads a{display:block;padding:9px 15px;border:1px solid #b9c6d6;border-radius:6px;text-decoration:none;color:#135e90}details{margin-top:18px}summary{cursor:pointer;font-weight:600}textarea{font:12px ui-monospace,monospace;min-height:220px;margin-top:10px}dialog{width:min(700px,94vw);border:1px solid #b9c6d6;border-radius:12px;padding:24px}dialog::backdrop{background:#13263977}#file-list{max-height:55vh;overflow:auto;margin:15px 0}#file-list button{display:block;width:100%;text-align:left;margin:5px 0}#picker-location{overflow-wrap:anywhere}#initial-version{font-size:13px;border:1px solid #486178;padding:6px 10px;border-radius:20px}@media(max-width:900px){main{grid-template-columns:1fr;padding:0 14px}header{padding:18px}}
</style></head><body>
<header><div><h1>SparK</h1><p>Build genomic figures on your computer</p></div><span id="initial-version">Local interface</span></header>
<main><div><section><h2>1. Choose your tracks</h2><p class="muted">Files stay on this computer. Add the same condition name for its replicates.</p><div id="tracks"></div><button id="add-track">+ Add track</button>
<h2 style="margin-top:24px">2. Choose a region</h2><label for="region-mode">Find by</label><select id="region-mode"><option value="region">Coordinates</option><option value="gene">Gene name or gene ID</option></select>
<label for="region">Region</label><input id="region" value="chr1:0-1000" placeholder="chr12:6527512-6550512"><div id="gene-options" hidden><label for="gene">Gene</label><input id="gene" placeholder="GAPDH"><label for="flank">Flank on each side (bases)</label><input id="flank" type="number" value="10000" min="0"></div>
<label for="gtf">GTF annotation (optional for coordinates)</label><div class="row"><input id="gtf" placeholder="Choose a GTF file"><button id="pick-gtf">Browse</button></div>
<div class="row"><div><label for="plot-type">Plot style</label><select id="plot-type"><option value="standard">Standard</option><option value="STD">Mean ± SD</option><option value="sine">Sine (two conditions)</option></select></div><div><label for="show-plots">Replicates</label><select id="show-plots"><option value="all">Show all</option><option value="averages">Show means</option></select></div></div>
<label for="output-name">Figure name</label><input id="output-name" value="SparK_figure"><div class="formats"><span class="muted">SVG + JPG included</span><label><input type="checkbox" id="pdf">PDF</label><label><input type="checkbox" id="png">PNG</label></div>
<div class="actions"><button id="render" class="primary">Preview &amp; save</button><button id="demo">Try demo</button></div>
<details><summary>Saved settings and advanced options</summary><p class="muted">Import a small JSON settings file or edit settings below. When this box is nonempty it supplies the full plot settings.</p><input id="import-config" type="file" accept=".json,application/json"><textarea id="config" aria-label="Saved plot settings" placeholder="Optional JSON configuration"></textarea><div class="actions"><button id="save-config">Download settings</button><button id="clear-config">Use the form</button></div></details>
</section></div><div><section><h2>Figure preview</h2><p id="status" role="status">Choose files and a region, or try the small synthetic demo.</p><div id="empty">Your editable figure and automatic JPG will appear here.</div><img id="preview" hidden alt="Generated SparK figure"><div id="downloads" class="downloads"></div></section><p class="muted" id="output-location"></p></div></main>
<dialog id="picker"><h2>Choose a local file</h2><p class="muted" id="picker-location"></p><button id="parent-folder">↑ Parent folder</button><div id="file-list"></div><button id="close-picker">Cancel</button></dialog>
<script>
const token=new URLSearchParams(location.search).get('token');
const api=async(path,options={})=>{const response=await fetch(path,{...options,headers:{'X-SparK-Token':token,'Content-Type':'application/json',...(options.headers||{})}});const data=await response.json();if(!response.ok)throw Error(data.error||'Request failed');return data;};
const $=id=>document.getElementById(id);let pickerTarget=null,pickerPath='';
function addTrack(track={}){const div=document.createElement('div');div.className='track';div.innerHTML='<div class="top"><input class="condition" aria-label="Condition" placeholder="Condition"><input class="color" type="color" aria-label="Track color"><button class="remove" aria-label="Remove track">×</button></div><div class="row"><input class="file" aria-label="Signal file" placeholder="Path to a BigWig or bedGraph"><button class="browse">Browse</button></div><input class="group" aria-label="Track group" placeholder="Group (optional)" style="margin-top:8px">';div.querySelector('.condition').value=track.condition||'Control';div.querySelector('.file').value=track.file||track.path||'';div.querySelector('.color').value=track.color||'#005cff';div.querySelector('.group').value=track.group||'';div.dataset.factor=track.scale_factor||1;div.querySelector('.remove').onclick=()=>div.remove();div.querySelector('.browse').onclick=()=>openPicker(div.querySelector('.file'));$('tracks').appendChild(div);}
function formConfig(){const tracks=[...document.querySelectorAll('.track')].map(el=>({file:el.querySelector('.file').value.trim(),condition:el.querySelector('.condition').value.trim(),color:el.querySelector('.color').value,group:el.querySelector('.group').value||'Signal',scale_factor:Number(el.dataset.factor)}));const result={tracks,plot_type:$('plot-type').value,show_plots:$('show-plots').value,output_name:$('output-name').value,formats:['svg','jpg',...($('pdf').checked?['pdf']:[]),...($('png').checked?['png']:[])]};if($('region-mode').value==='gene'){result.gene=$('gene').value.trim();result.flank=Number($('flank').value);}else result.region=$('region').value.trim();if($('gtf').value.trim())result.gtffile=$('gtf').value.trim();return result;}
function getConfig(){return $('config').value.trim()?JSON.parse($('config').value):formConfig();}
function status(text,error=false){$('status').textContent=text;$('status').className=error?'error':'';}
$('add-track').onclick=()=>addTrack({condition:'Treatment',color:'#dd5500'});
$('pick-gtf').onclick=()=>openPicker($('gtf'));
$('region-mode').onchange=()=>{const gene=$('region-mode').value==='gene';$('gene-options').hidden=!gene;$('region').hidden=gene;document.querySelector('label[for=region]').hidden=gene;};
$('render').onclick=async()=>{try{$('render').disabled=true;status('Reading local data and creating your figure…');const data=await api('/api/plot',{method:'POST',body:JSON.stringify(getConfig())});$('empty').hidden=true;$('preview').hidden=false;$('preview').src=data.preview;$('downloads').replaceChildren();for(const [name,url] of Object.entries(data.links)){const a=document.createElement('a');a.href=url;a.textContent=name.toUpperCase();a.download='';$('downloads').appendChild(a);}status('Saved '+data.region+'. SVG, JPG and run settings are ready.');$('output-location').textContent='Saved in '+data.directory;}catch(error){status(error.message,true);}finally{$('render').disabled=false;}};
$('demo').onclick=async()=>{try{const data=await api('/api/demo');$('config').value='';$('tracks').replaceChildren();for(const track of data.tracks)addTrack(track);$('region-mode').value='region';$('region-mode').onchange();$('region').value=data.region;$('gtf').value='';$('plot-type').value='standard';$('show-plots').value=data.show_plots;$('output-name').value=data.output_name;$('pdf').checked=data.formats.includes('pdf');$('png').checked=data.formats.includes('png');status('Synthetic demo loaded. Edit the form or click Preview & save.');}catch(error){status(error.message,true);}};
$('save-config').onclick=()=>{try{const a=document.createElement('a');const url=URL.createObjectURL(new Blob([JSON.stringify(getConfig(),null,2)+'\n'],{type:'application/json'}));a.href=url;a.download='SparK_settings.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}catch(error){status(error.message,true);}};
$('clear-config').onclick=()=>{$('config').value='';status('Using the form settings.');};
$('import-config').onchange=async event=>{try{const file=event.target.files[0];if(!file)return;if(file.size>1048576)throw Error('Choose a small JSON settings file, not a genomic data file.');const text=await file.text();JSON.parse(text);$('config').value=text;document.querySelector('details').open=true;status('Saved settings loaded. Advanced JSON supplies the settings; click Use the form to switch back.');}catch(error){status(error.message,true);}};
async function openPicker(target){pickerTarget=target;pickerPath='';$('picker').showModal();await listFiles();}
async function listFiles(){try{const data=await api('/api/files?path='+encodeURIComponent(pickerPath));$('picker-location').textContent=data.directory;$('file-list').replaceChildren();for(const item of data.files){const button=document.createElement('button');button.textContent=(item.directory?'📁 ':'')+item.name;button.onclick=()=>{if(item.directory){pickerPath=item.relative;listFiles();}else{pickerTarget.value=item.path;$('picker').close();}};$('file-list').appendChild(button);}}catch(error){status(error.message,true);}}
$('parent-folder').onclick=()=>{pickerPath=pickerPath.split('/').slice(0,-1).join('/');listFiles();};$('close-picker').onclick=()=>$('picker').close();
api('/api/config').then(data=>{$('initial-version').textContent='Version '+data.version+' · local only';if(data.config&&Object.keys(data.config).some(key=>!['output_name','no_update_check'].includes(key))){$('config').value=JSON.stringify(data.config,null,2);document.querySelector('details').open=true;status('Advanced JSON supplies the saved settings. Click Use the form to switch back.');}$('output-location').textContent='Outputs: '+data.directory;}).catch(error=>status(error.message,true));
addTrack({condition:'Control',color:'#005cff'});addTrack({condition:'Treatment',color:'#dd5500'});
</script></body></html>'''


def create_gui_server(config=None, data_dir='.', output_dir=None, port=0):
    """Create a token-protected loopback server; no files are uploaded to a cloud service."""
    from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
    from urllib.parse import urlsplit, parse_qs, quote
    import secrets
    root = Path(data_dir).expanduser().resolve()
    if not root.is_dir():
        raise ValueError('The interface data directory does not exist: '+str(root))
    output = Path(output_dir or root/'SparK_GUI_outputs').expanduser().resolve()
    output.mkdir(parents=True,exist_ok=True)
    token = secrets.token_urlsafe(24)
    initial = copy.deepcopy(config or {})
    for key in ('gui','data_dir','port','no_open','config'):
        initial.pop(key,None)

    class Handler(BaseHTTPRequestHandler):
        def log_message(self,*args):
            pass

        def allowed(self):
            parsed = urlsplit(self.path)
            host = self.headers.get('Host','')
            allowed = {'127.0.0.1:'+str(self.server.server_port),'localhost:'+str(self.server.server_port)}
            origin = self.headers.get('Origin')
            if host not in allowed or (origin and origin not in {'http://'+value for value in allowed}):
                return False
            supplied = self.headers.get('X-SparK-Token') or parse_qs(parsed.query).get('token',[''])[0]
            return secrets.compare_digest(str(supplied),token)

        def send(self,status,data,content_type='application/json'):
            payload = json.dumps(data).encode() if content_type=='application/json' else data
            self.send_response(status)
            self.send_header('Content-Type',content_type)
            self.send_header('Content-Length',str(len(payload)))
            self.send_header('Cache-Control','no-store')
            self.send_header('Referrer-Policy','no-referrer')
            self.send_header('X-Content-Type-Options','nosniff')
            self.end_headers()
            self.wfile.write(payload)

        def do_GET(self):
            if not self.allowed():
                self.send(403,dict(error='Open the local URL printed by SparK to access this interface'))
                return
            from urllib.parse import unquote
            parsed = urlsplit(self.path)
            try:
                if parsed.path=='/':
                    self.send(200,GUI_HTML.encode(),'text/html; charset=utf-8')
                elif parsed.path=='/api/config':
                    self.send(200,dict(version=SparK_Version,config=initial,directory=str(output)))
                elif parsed.path=='/api/files':
                    relative = parse_qs(parsed.query).get('path',[''])[0]
                    directory = (root/relative).resolve()
                    if not directory.is_relative_to(root) or not directory.is_dir():
                        raise ValueError('Choose a folder inside the starting data directory')
                    files = []
                    for path in sorted(directory.iterdir(),key=lambda p:(not p.is_dir(),p.name.lower())):
                        if path.name.startswith('.') or not path.resolve().is_relative_to(root):
                            continue
                        files.append(dict(name=path.name,path=str(path),relative=str(path.relative_to(root)),directory=path.is_dir()))
                        if len(files)>=1000:
                            break
                    self.send(200,dict(directory=str(directory),files=files))
                elif parsed.path=='/api/demo':
                    self.send(200,self.server.demo_config())
                elif parsed.path.startswith('/files/'):
                    file = (output/unquote(parsed.path[len('/files/'):])).resolve()
                    if not file.is_relative_to(output) or not file.is_file():
                        raise ValueError('The output file was not found')
                    types = {'.svg':'image/svg+xml','.jpg':'image/jpeg','.png':'image/png','.pdf':'application/pdf','.json':'application/json'}
                    if file.suffix not in types:
                        raise ValueError('This file type is not available for download')
                    self.send(200,file.read_bytes(),types[file.suffix])
                else:
                    self.send(404,dict(error='Not found'))
            except (OSError,ValueError) as exc:
                self.send(400,dict(error=str(exc)))

        def do_POST(self):
            if not self.allowed():
                self.send(403,dict(error='This request is not authorized for the local interface'))
                return
            try:
                if urlsplit(self.path).path!='/api/plot':
                    self.send(404,dict(error='Not found'))
                    return
                length = int(self.headers.get('Content-Length','0'))
                if self.headers.get('Content-Type','').split(';')[0]!='application/json' or not 0<length<=1048576:
                    raise ValueError('Send a small JSON settings object; genomic files are read directly from disk')
                settings = json.loads(self.rfile.read(length))
                if not isinstance(settings,dict):
                    raise ValueError('Settings must be an object')
                if any(settings.get(key) for key in ('genes_file','regions_file','genes','regions')):
                    raise ValueError('Use the command line or Python API for a batch; this preview shows one region')
                name = Path(str(settings.get('output_name') or settings.pop('output','SparK_figure'))).name
                name = re.sub(r'[^A-Za-z0-9_.-]+','_',name).strip('._') or 'SparK_figure'
                settings['output_name'] = str(output/name)
                for key in ('gui','port','data_dir','no_open','save_config','output'):
                    settings.pop(key,None)
                settings['no_update_check'] = True
                result = plot(config=settings)
                links = {key:'/files/'+quote(str(Path(path).relative_to(output)))+'?token='+token for key,path in result['outputs'].items()}
                self.send(200,dict(region=result['region'],links=links,preview=links.get('jpg') or links.get('png') or links['svg'],directory=str(output)))
            except (OSError,ValueError,ImportError,RuntimeError) as exc:
                self.send(400,dict(error=str(exc)))

    server = ThreadingHTTPServer(('127.0.0.1',port),Handler)
    server.daemon_threads = True
    server.token,server.output_dir,server.data_dir = token,output,root
    server.demo_directory = None

    def demo_config():
        if server.demo_directory is None:
            server.demo_directory = Path(tempfile.mkdtemp(prefix='spark-demo-'))
        tracks = []
        for c,(condition,color) in enumerate([('Control','#005CFF'),('Drug A','#DD5500'),('Drug B','#009E73')]):
            for replicate in range(2):
                file = server.demo_directory/('condition{}_rep{}.bdg'.format(c,replicate+1))
                if not file.exists():
                    with file.open('w') as handle:
                        for start in range(0,1000,10):
                            value = 2+(18+4*c)*math.exp(-((start-300-80*c)/90)**2)+(8+2*replicate)*math.exp(-((start-700)/70)**2)
                            handle.write('chr1\t{}\t{}\t{:.6f}\n'.format(start,start+10,value))
                tracks.append(dict(file=str(file),condition=condition,color=color,group='Synthetic signal'))
        return dict(tracks=tracks,region='chr1:0-1000',show_plots='averages',formats=['svg','jpg','pdf'],output_name='SparK_demo')
    server.demo_config = demo_config
    return server


def serve_gui(config, checker=None):
    import webbrowser
    import shutil
    defaults = default_options()
    defaults.update(config)
    server = create_gui_server(config,data_dir=defaults['data_dir'],output_dir=defaults.get('output_name'),port=defaults['port'])
    url = 'http://127.0.0.1:{}/?token={}'.format(server.server_port,server.token)
    print('SparK {} local interface: {}'.format(SparK_Version,url),flush=True)
    print('Files stay on this computer. Stop with Ctrl+C.',flush=True)
    if checker:
        checker.finish(True)
    if not defaults['no_open']:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
        if server.demo_directory:
            shutil.rmtree(server.demo_directory,ignore_errors=True)


if __name__ == '__main__':
    raise SystemExit(main())
