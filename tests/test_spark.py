import contextlib
from fractions import Fraction
import gzip
import importlib.util
import math
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('spark', ROOT / 'SparK.py')
spark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(spark)


def polygons(tree):
    return tree.findall('{http://www.w3.org/2000/svg}polygon')


def points(polygon):
    return [tuple(map(float, p.split(','))) for p in polygon.get('points').split()]


class ReaderTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        spark.max_points, spark.point_stat, spark.smoothen_tracks = 2000, 'average', None

    def read(self, content, start=0, stop=10, offset=0):
        path = self.root / 'input.bdg'
        path.write_text(content)
        return spark.make_raw_data_filled(['1', start, stop, 'chr1'], [str(path)], offset)[0]

    def test_left_boundary_does_not_wrap(self):
        self.assertEqual(self.read('chr1\t0\t5\t10\n', 3, 13), [10,10]+[0]*8)

    def test_right_boundary_is_clipped(self):
        self.assertEqual(self.read('chr1\t8\t15\t4\n'), [0]*8+[4]*2)

    def test_offset_has_same_sign_for_query_and_assignment(self):
        self.assertEqual(self.read('chr1\t4\t7\t9\n', 3, 8, 1), [9,9,9,0,0])

    def test_noncontiguous_chromosomes_and_overlapping_records(self):
        text = 'chr1\t0\t4\t1\nchr2\t0\t10\t99\n1\t2\t6\t3\nChr1\t8\t10\t5\n'
        self.assertEqual(self.read(text), [1,1,3,3,3,3,0,0,5,5])

    def test_headers_and_empty_lines(self):
        self.assertEqual(self.read('track type=bedGraph\n# comment\nbrowser position chr1\n\nchr1\t0\t10\t1\n'), [1]*10)

    def test_invalid_records(self):
        for row in ['chr1\t0\t10\tnan\n','chr1\t0\t10\tinf\n','chr1\t-1\t10\t1\n', 'chr1\t4\t4\t1\n', 'chr1\tbad\t10\t1\n']:
            with self.subTest(row=row), self.assertRaises(ValueError):
                self.read(row)

    def test_cached_and_streamed_readers_match_numeric_oracle(self):
        rng = random.Random(17)
        records = [(rng.choice(['chr1','chr2']), rng.randrange(20), rng.randrange(1,7), rng.randrange(-4,5)) for _ in range(80)]
        text = ''.join(f'{chrom}\t{start}\t{start+length}\t{value}\n' for chrom,start,length,value in records)
        path = self.root/'random.bdg'; path.write_text(text)
        for start, stop in [(0,10),(3,19),(15,27)]:
            expected = [0]*(stop-start)
            for chrom, a, length, value in records:
                if chrom == 'chr1':
                    for x in range(max(start,a), min(stop,a+length)):
                        expected[x-start] = value
            for disabled in ['0','1']:
                with patch.dict(spark.os.environ, {'SPARK_DISABLE_TEXT_INDEX':disabled}):
                    self.assertEqual(spark.make_raw_data_filled(['1',start,stop],[str(path)],0)[0], expected)

    def test_gzip_and_universal_newlines(self):
        for separator in ['\n','\r','\r\n']:
            text = separator.join(['chr1\t0\t5\t1','chr1\t5\t10\t2'])
            path = self.root/'input.bdg'; path.write_bytes(text.encode())
            self.assertEqual(spark.make_raw_data_filled(['1',0,10],[str(path)],0)[0], [1]*5+[2]*5)
        path = self.root/'input.bdg.gz'
        with gzip.open(path,'wt') as handle:
            handle.write('chr1\t0\t10\t7\n')
        self.assertEqual(spark.make_raw_data_filled(['1',0,10],[str(path)],0)[0], [7]*10)

    def test_smoothing_has_no_padding_bins(self):
        spark.smoothen_tracks = 3
        result = self.read('chr1\t0\t1\t3\nchr1\t1\t2\t6\nchr1\t2\t5\t9\n',0,5)
        self.assertEqual(result, [3,4,6,8,9])

    def test_smoothing_window_larger_than_region(self):
        spark.smoothen_tracks = 10
        self.assertEqual(self.read('chr1\t0\t1\t1\nchr1\t1\t2\t3\n',0,2), [1,1.2])

    def test_tabix_query_conversion_and_failure(self):
        with patch.object(spark.subprocess, 'run', return_value=subprocess.CompletedProcess([],0,'chr1\t0\t1\t1\n','')) as run:
            self.assertEqual(list(spark.get_lines_tabix('a.bdg.gz',['1',0,1,'chr1'])), ['chr1\t0\t1\t1'])
            self.assertEqual(run.call_args.args[0], ['tabix','a.bdg.gz','chr1:1-1'])
            self.assertTrue(run.call_args.kwargs['check'])

    def test_large_index_blocks_keep_spanning_intervals(self):
        text = 'chr1\t0\t10000\t9\n' + ''.join(f'chr1\t{i}\t{i+1}\t{i%7}\n' for i in range(9000))
        self.assertEqual(self.read(text,8998,9003), [8998%7,8999%7,9,9,9])

    def test_cache_corruption_and_content_change(self):
        cache = self.root/'cache'
        with patch.object(spark.tempfile,'gettempdir',return_value=str(cache)):
            self.assertEqual(self.read('chr1\t0\t10\t1\n'),[1]*10)
            for p in cache.rglob('*.json'):
                p.write_text('invalid JSON')
            self.assertEqual(self.read('chr1\t0\t10\t7\n'),[7]*10)

    def test_unwritable_cache_falls_back(self):
        with patch.object(spark.os,'makedirs',side_effect=PermissionError('read only')):
            self.assertEqual(self.read('chr1\t0\t10\t2\n'),[2]*10)

    def test_empty_file_is_zero_signal(self):
        self.assertEqual(self.read(''),[0]*10)

    def test_bigwig_missing_values_boundaries_and_chromosome(self):
        if spark.pyBigWig is None:
            self.skipTest('pyBigWig not installed')
        path = self.root/'test.bw'
        with contextlib.closing(spark.pyBigWig.open(str(path),'w')) as bw:
            bw.addHeader([('chr1',10)])
            bw.addEntries(['chr1','chr1'],[0,5],ends=[3,10],values=[2.0,-4.0])
        result = spark.make_raw_data_filled(['1',0,12],[str(path)],0)[0]
        self.assertEqual(result, [2]*3+[0]*2+[-4]*5+[0]*2)
        with self.assertRaisesRegex(ValueError,'absent'):
            spark.make_raw_data_filled(['2',0,10],[str(path)],0)


class DownsamplingTests(unittest.TestCase):
    def test_fractional_bins_against_exact_overlap_oracle(self):
        # Both sides of each bin boundary are represented as rational numbers.
        for size in range(2,32):
            values = [((x*7)%13)-6 for x in range(size)]
            for bins in range(1,size):
                expected, peaks = [], []
                for b in range(bins):
                    left, right = Fraction(b*size,bins), Fraction((b+1)*size,bins)
                    weights = [max(Fraction(0), min(right,i+1)-max(left,i)) for i in range(size)]
                    expected.append(float(sum(v*w for v,w in zip(values,weights))/(right-left)))
                    peaks.append(max(v for v,w in zip(values,weights) if w))
                actual = spark.downsample(values,bins,'average')
                self.assertEqual(len(actual),bins)
                for x,y in zip(actual,expected):
                    self.assertAlmostEqual(x,y,places=12)
                self.assertEqual(spark.downsample(values,bins,'max'),peaks)

    def test_realistic_fractional_final_bin(self):
        self.assertEqual(spark.downsample([1]*4001,2000), [1]*2000)

    def test_max_preserves_spikes_and_average_preserves_integral(self):
        values = [0,0,9,0,0]
        self.assertEqual(max(spark.downsample(values,3,'max')),9)
        self.assertAlmostEqual(sum(spark.downsample(values,3))*5/3,9)


class CommandTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.files = []
        for n, value in enumerate([1,2,3,100]):
            path = self.root/f'track{n}.bdg'
            path.write_text(f'chr1\t0\t100\t{value}\n')
            self.files.append(str(path))

    def run_plot(self, options=(), files=None, success=True):
        output = self.root/'plot.svg'
        command = [sys.executable,str(ROOT/'SparK.py'),'-pr','chr1:0-5','-cf',*(files or self.files),'-o',str(output),*options]
        result = subprocess.run(command,capture_output=True,text=True)
        if success:
            self.assertEqual(result.returncode,0,result.stderr)
            return ET.parse(output).getroot()
        self.assertNotEqual(result.returncode,0)
        self.assertNotIn('Traceback',result.stderr)
        return result

    def test_four_individual_colors(self):
        tree = self.run_plot(['-f','000000','FF0000','00FF00','0000FF'])
        self.assertEqual([p.get('fill') for p in polygons(tree)], ['#000000','#FF0000','#00FF00','#0000FF'])

    def test_explicit_two_track_colors(self):
        tree = self.run_plot(['--track_colors','#abc','123456'],self.files[:2])
        self.assertEqual([p.get('fill') for p in polygons(tree)], ['#abc','#123456'])

    def test_legacy_two_colors_control_then_treatment(self):
        tree = self.run_plot(['-tf',self.files[1],'-f','FF0000','0000FF'],self.files[:1])
        self.assertEqual([p.get('fill') for p in polygons(tree)], ['#FF0000','#0000FF'])

    def test_default_treatment_fill_is_valid(self):
        tree = self.run_plot(['-tf',self.files[1]],self.files[:1])
        self.assertTrue(all(p.get('fill').startswith('#') for p in polygons(tree)))

    def test_invalid_colors_are_clean_errors(self):
        for colors in [['abcde'],['red,blue'],['ff0000','00ff00','0000ff']]:
            with self.subTest(colors=colors):
                self.run_plot(['-f',*colors],success=False)

    def test_missing_and_extra_group_assignments(self):
        for values in [['1','1','2'],['1','1','2','2','2']]:
            self.run_plot(['-cg',*values],success=False)

    def test_excluded_and_sparse_groups(self):
        for options in [['-cg','1','1','2','2','-eg','1'],['-cg','1','1','3','3','-eg','1']]:
            tree = self.run_plot(options)
            self.assertEqual(len(polygons(tree)),2)
            self.assertEqual(max(y for x,y in points(polygons(tree)[0])),100)
            self.assertEqual(min(y for x,y in points(polygons(tree)[1])),73)

    def test_excluding_every_group_is_error(self):
        self.run_plot(['-cg','1','1','2','2','-eg','1','2'],success=False)

    def test_treatment_only_average_group(self):
        tree = self.run_plot(['-tf',self.files[1],'-cg','1','-tg','2','-ps','averages'], self.files[:1])
        self.assertEqual(len(polygons(tree)),2)

    def test_sine_average_count_and_shared_scale(self):
        tree = self.run_plot(['-tf',*self.files,'-pt','sine','-ps','averages'])
        self.assertEqual(len(polygons(tree)),2)
        low, high = (points(p) for p in polygons(tree))
        self.assertAlmostEqual(max(y for x,y in low)-100,100-min(y for x,y in high))

    def test_sine_shared_scale_preserves_unequal_amplitudes(self):
        tree = self.run_plot(['-tf',self.files[1],'-pt','sine'],self.files[:1])
        low, high = (points(p) for p in polygons(tree))
        self.assertAlmostEqual((max(y for x,y in low)-100)/(100-min(y for x,y in high)),0.5)

    def test_zero_sine_signal_still_has_axis_labels(self):
        for path in self.files[:2]:
            Path(path).write_text('chr1\t0\t100\t0\n')
        tree = self.run_plot(['-tf',self.files[1],'-pt','sine'],self.files[:1])
        self.assertTrue(tree.findall('{http://www.w3.org/2000/svg}text'))

    def test_control_only_standard_deviation(self):
        tree = self.run_plot(['-pt','STD'],self.files[:2])
        self.assertEqual(len(polygons(tree)),1)

    def test_no_endpoint_signal_is_erased(self):
        tree = self.run_plot(files=self.files[:1])
        coordinates = points(polygons(tree)[0])
        self.assertIn((50.0,73.0),coordinates)
        self.assertIn((200.0,73.0),coordinates)

    def test_averages_with_conflicting_individual_colors_are_rejected(self):
        self.run_plot(['--track_colors','fff','000','f00','00f','-ps','averages'],success=False)

    def test_negative_shared_autoscale_preserves_direction(self):
        for path in self.files:
            Path(path).write_text('chr1\t0\t100\t-2\n')
        tree = self.run_plot(['-gs','yes'])
        self.assertGreater(max(y for x,y in points(polygons(tree)[0])),100)
        labels = [e.text for e in tree.findall('{http://www.w3.org/2000/svg}text')]
        self.assertIn('2.2',labels)
        self.assertIn('-2.2',labels)

    def test_all_groups_excluded_from_shared_scaling(self):
        tree = self.run_plot(['-cg','1','1','2','2','-gs','yes','-es','1','2'])
        self.assertEqual(len(polygons(tree)),4)

    def test_invalid_argument_values(self):
        for options in [['-pr','chr1-10'],['-pr','chr1:5-5'],['-pr','chr1:5-1'],['--max_points','0'],
                        ['-sm','0'],['--x_scale','nan'],['--y_scale','-1'],['-cs','0','D','D','D'],
                        ['-cg','0','1','1','1'],['-dc','wrong'],['--point_stat','wrong'],['-sp','yes']]:
            with self.subTest(options=options):
                self.run_plot(options,success=False)

    def test_error_does_not_replace_previous_figure(self):
        destination = self.root/'plot.svg'; destination.write_text('previous figure')
        Path(self.files[0]).write_text('chr1\t0\t100\tnan\n')
        self.run_plot(success=False)
        self.assertEqual(destination.read_text(),'previous figure')
        self.assertFalse(list(self.root.glob('.spark-*')))

    def test_annotation_error_does_not_publish_partial_svg(self):
        destination = self.root/'plot.svg'; destination.write_text('previous figure')
        gtf = self.root/'bad.gtf'; gtf.write_text('chr1\ttest\tgene\tbad\t5\t.\t+\t.\tgene_name "x";\n')
        self.run_plot(['-gtf',str(gtf)],success=False)
        self.assertEqual(destination.read_text(),'previous figure')
        self.assertFalse(list(self.root.glob('.spark-*')))

    def test_invalid_bed_is_clean_error(self):
        bed = self.root/'bad.bed'; bed.write_text('chr1\t3\n')
        self.run_plot(['-bed',str(bed)],success=False)

    def test_translation_and_transcription_arrows_are_optional(self):
        gtf = self.root/'arrows.gtf'
        gtf.write_text('chr1\ttest\tgene\t1\t5\t.\t-\t.\tgene_name "x";\n'
                       'chr1\ttest\tstart_codon\t2\t4\t.\t-\t.\tgene_name "x";\n')
        tree = self.run_plot(['-gtf',str(gtf),'-genestart','yes'],self.files[:1])
        self.assertEqual(len(polygons(tree)),3)  # signal + codon + gene start
        tree = self.run_plot(['-gtf',str(gtf),'-genestart','no','-tss','no'],self.files[:1])
        self.assertEqual(len(polygons(tree)),1)

    def test_sine_highlights_use_shared_signal_scale(self):
        for path,value in zip(self.files,[4,4,1,1]):
            Path(path).write_text(f'chr1\t0\t5\t{value}\n')
        tree = self.run_plot(['-tf',*self.files[2:],'-pt','sine','-sp','yes'],self.files[:2])
        highlight = next(p for p in polygons(tree) if p.get('fill')=='#00FF00')
        self.assertAlmostEqual(max(y for x,y in points(highlight)),120.25)

    def test_labels_are_escaped(self):
        tree = self.run_plot(['-gl','A&B','C<D','E>F','G"H','-l','C&T','<T>'])
        labels = [e.text for e in tree.findall('{http://www.w3.org/2000/svg}text')]
        self.assertIn('A&B',labels)
        self.assertIn('<T>',labels)

    def test_gtf_inclusive_single_base_and_gene_id_fallback(self):
        gtf = self.root/'genes.gtf'
        gtf.write_text('chr1\ttest\tgene\t2\t2\t.\t+\t.\tgene_id "A&B";\n')
        tree = self.run_plot(['-gtf',str(gtf)],self.files[:1])
        annotations = [r for r in tree.findall('{http://www.w3.org/2000/svg}rect') if r.get('fill')=='#0B34FF']
        self.assertEqual(len(annotations),1)
        self.assertEqual(float(annotations[0].get('x')),80)
        self.assertEqual(float(annotations[0].get('width')),30)

    def test_canvas_expands_for_wide_tracks_and_many_genes(self):
        gtf = self.root/'genes.gtf'
        gtf.write_text(''.join(f'chr1\ttest\tgene\t1\t5\t.\t+\t.\tgene_name "gene{n}";\n' for n in range(40)))
        tree = self.run_plot(['-gtf',str(gtf),'--x_scale','4'],self.files[:1])
        _,_,width,height = map(float,tree.get('viewBox').split())
        self.assertGreaterEqual(width,770)
        self.assertGreater(height,600)

    def test_sparks_split_at_nonsignificant_gap_and_keep_direction_color(self):
        for path in self.files[:2]:
            Path(path).write_text('chr1\t0\t2\t4\nchr1\t2\t3\t1\nchr1\t3\t5\t4\n')
        for path in self.files[2:]:
            Path(path).write_text('chr1\t0\t5\t1\n')
        tree = self.run_plot(['-tf',*self.files[2:],'-cg','1','1','-tg','1','1','-sp','yes'],self.files[:2])
        highlights = [p for p in polygons(tree) if p.get('fill')=='#00FF00']
        self.assertEqual(len(highlights),2)
        self.assertLess(max(x for x,y in points(highlights[0])),min(x for x,y in points(highlights[1])))


if __name__ == '__main__':
    unittest.main()
