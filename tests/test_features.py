import concurrent.futures
import ast
import contextlib
import gzip
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from unittest.mock import patch
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('spark_features',ROOT/'SparK.py')
spark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(spark)
NS = '{http://www.w3.org/2000/svg}'


class FeatureTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.files = []
        for index,value in enumerate([1,3,5,7,9,11]):
            path = self.root/('track{}.bdg'.format(index))
            path.write_text('chr1\t0\t10\t{}\nchr1\t10\t20\t{}\n'.format(value,value*2))
            self.files.append(str(path))
        self.gtf = self.root/'genes.gtf'
        self.gtf.write_text('chr1\ttest\tgene\t3\t6\t.\t+\t.\tgene_id "ID_A"; gene_name "A";\n'
                            'chr1\ttest\tgene\t12\t18\t.\t-\t.\tgene_id "ID_B"; gene_name "B";\n')

    def config(self, **options):
        return dict(control_files=self.files[:1],region='chr1:0-20',output_name=str(self.root/'figure'),**options)

    def named(self):
        return {'Control':self.files[:2],'Drug A':self.files[2:4],'Drug B':self.files[4:]}

    def test_default_jpg_and_run_record(self):
        from PIL import Image
        result = spark.plot(config=self.config())
        self.assertEqual(set(result['outputs']),{'svg','jpg','record'})
        with Image.open(result['outputs']['jpg']) as image:
            self.assertEqual(image.format,'JPEG')
            self.assertEqual(image.mode,'RGB')
            self.assertEqual(image.width,1000)
            self.assertAlmostEqual(image.info['dpi'][0],300,delta=1)
        record = json.loads(Path(result['outputs']['record']).read_text())
        self.assertEqual(record['spark_version'],spark.SparK_Version)
        self.assertEqual(record['inputs'][0]['verification'],'metadata')
        self.assertEqual(record['outputs']['svg']['sha256'],spark.file_sha256(result['outputs']['svg']))
        self.assertFalse(list(self.root.glob('.spark-output-*')))

    def test_all_export_formats_and_vector_pdf(self):
        result = spark.plot(config=self.config(formats=['svg','jpg','png','pdf']))
        self.assertEqual(set(result['outputs']),{'svg','jpg','png','pdf','record'})
        pdf = Path(result['outputs']['pdf']).read_bytes()
        self.assertTrue(pdf.startswith(b'%PDF'))
        self.assertNotIn(b'/Subtype /Image',pdf)
        self.assertIn(b'/MediaBox',pdf)

    def test_no_jpg_opt_out(self):
        result = spark.plot(config=self.config(no_jpg=True))
        self.assertNotIn('jpg',result['outputs'])
        self.assertFalse((self.root/'figure.jpg').exists())

    def test_no_record_opt_out(self):
        result = spark.plot(config=self.config(no_run_record=True))
        self.assertNotIn('record',result['outputs'])

    def test_export_failure_preserves_prior_outputs(self):
        svg,jpg = self.root/'figure.svg',self.root/'figure.jpg'
        svg.write_text('old SVG');jpg.write_bytes(b'old JPG')
        with patch.object(spark,'export_figures',side_effect=ValueError('export failed')):
            with self.assertRaisesRegex(ValueError,'export failed'):
                spark.plot(config=self.config())
        self.assertEqual(svg.read_text(),'old SVG')
        self.assertEqual(jpg.read_bytes(),b'old JPG')

    def test_input_hashes_are_optional_and_correct(self):
        result = spark.plot(config=self.config(hash_inputs=True))
        record = json.loads(Path(result['outputs']['record']).read_text())
        self.assertEqual(record['inputs'][0]['sha256'],spark.file_sha256(self.files[0]))
        self.assertEqual(record['inputs'][0]['verification'],'sha256')

    def test_config_paths_resolve_relative_to_config(self):
        config = self.root/'settings.json'
        config.write_text(json.dumps({'control_files':['track0.bdg'],'region':'chr1:0-20','output_name':'relative'}))
        result = spark.plot(config=config)
        self.assertTrue((self.root/'relative.jpg').exists())
        self.assertEqual(result['region'],'chr1:0-20')

    def test_save_and_reload_config_reproduce_svg(self):
        config = self.root/'saved.json'
        first = spark.plot(config=self.config(save_config=str(config)))
        expected = Path(first['outputs']['svg']).read_bytes()
        second = spark.plot(config=config,output=self.root/'second')
        self.assertEqual(Path(second['outputs']['svg']).read_bytes(),expected)

    def test_named_config_round_trip(self):
        config = self.root/'named.json'
        result = spark.plot('chr1:0-20',tracks=self.named(),output=self.root/'named',show_plots='averages',save_config=str(config))
        again = spark.plot(config=config,output=self.root/'again')
        self.assertEqual(Path(result['outputs']['svg']).read_bytes(),Path(again['outputs']['svg']).read_bytes())

    def test_unknown_and_invalid_config_are_errors(self):
        for config in [dict(region='chr1:0-20',control_files=self.files[:1],typo=True),
                       dict(region='chr1:0-20',control_files=self.files[:1],control_groups=['1']),
                       dict(region='chr1:0-20',control_files=self.files[:1],formats=['gif'])]:
            with self.subTest(config=config),self.assertRaises(ValueError):
                spark.plot(config=config)

    def test_three_conditions_with_means(self):
        result = spark.plot('chr1:0-20',tracks=self.named(),output=self.root/'three',show_plots='averages')
        tree = ET.parse(result['outputs']['svg']).getroot()
        self.assertEqual(len(tree.findall(NS+'polygon')),3)
        labels = [element.text for element in tree.findall(NS+'text')]
        for condition in self.named():
            self.assertIn(condition,labels)

    def test_three_conditions_with_sd(self):
        result = spark.plot('chr1:0-20',tracks=self.named(),output=self.root/'sd',plot_type='STD')
        tree = ET.parse(result['outputs']['svg']).getroot()
        self.assertEqual(len(tree.findall(NS+'polygon')),3)
        self.assertEqual(len(tree.findall(NS+'path')),3)

    def test_interleaved_replicates_have_consistent_colors(self):
        tracks = [dict(file=self.files[0],condition='Control'),dict(file=self.files[2],condition='Drug'),
                  dict(file=self.files[1],condition='Control'),dict(file=self.files[3],condition='Drug')]
        result = spark.plot('chr1:0-20',tracks=tracks,output=self.root/'interleaved',show_plots='averages')
        self.assertEqual(len(ET.parse(result['outputs']['svg']).getroot().findall(NS+'polygon')),2)

    def test_more_than_two_sine_conditions_is_error(self):
        with self.assertRaisesRegex(ValueError,'exactly two'):
            spark.plot('chr1:0-20',tracks=self.named(),output=self.root/'bad',plot_type='sine')

    def transformed(self, kind, pseudocount=1):
        args = spark.normalize_config(dict(tracks=self.named(),region='chr1:0-20',transform=kind,compare=['Drug A','Control'],pseudocount=pseudocount))
        spark.finalize_arguments(args)
        loaded = {record['file']:spark.read_signals(['1',0,20],[record['file']])[0] for record in args['_track_records']}
        return spark.build_group_series(args,loaded)[1][0]['data'][0]

    def test_difference_uses_condition_means(self):
        self.assertEqual(self.transformed('difference'),[4]*10+[8]*10)

    def test_ratio_and_log2ratio(self):
        import math
        self.assertEqual(self.transformed('ratio'),[7/3]*10+[13/5]*10)
        self.assertAlmostEqual(self.transformed('log2ratio')[0],math.log2(7/3))

    def test_multi_condition_transform_requires_pair(self):
        with self.assertRaisesRegex(ValueError,'--compare'):
            spark.plot('chr1:0-20',tracks=self.named(),output=self.root/'bad',transform='difference')

    def test_signed_ratio_with_invalid_denominator_is_error(self):
        for path in self.files[:2]:
            Path(path).write_text('chr1\t0\t20\t-5\n')
        with self.assertRaisesRegex(ValueError,'pseudocount'):
            self.transformed('log2ratio')

    def test_transformed_figure_has_comparison_label(self):
        result = spark.plot('chr1:0-20',tracks=self.named(),output=self.root/'delta',transform='difference',compare=['Drug A','Control'])
        labels = [e.text for e in ET.parse(result['outputs']['svg']).getroot().findall(NS+'text')]
        self.assertIn('Drug A − Control',labels)

    def test_explicit_scale_factors(self):
        args = spark.normalize_config(dict(control_files=self.files[:1],region='chr1:0-20',scale_factors=[2]))
        spark.finalize_arguments(args)
        data = spark.build_group_series(args,{args['_track_records'][0]['file']:[1,2]})[1][0]['data'][0]
        self.assertEqual(data,[2,4])

    def test_gene_lookup_by_name_or_id(self):
        self.assertEqual(spark.resolve_gene('A',self.gtf,2)[0],'chr1:0-8')
        self.assertEqual(spark.resolve_gene('ID_B',self.gtf,0)[0],'chr1:11-18')
        result = spark.plot(tracks={'Control':self.files[:1]},gene='A',gtffile=str(self.gtf),flank=2,output=self.root/'gene')
        self.assertEqual(result['region'],'chr1:0-8')
        self.assertEqual(result['gene']['strand'],'+')

    def test_gene_gzip_and_cache_invalidation(self):
        zipped = self.root/'genes.gtf.gz'
        with gzip.open(zipped,'wt') as handle:
            handle.write(self.gtf.read_text())
        self.assertEqual(spark.resolve_gene('A',zipped,0)[0],'chr1:2-6')
        self.gtf.write_text(self.gtf.read_text().replace('gene_name "A"','gene_name "C"'))
        with self.assertRaisesRegex(ValueError,'not found'):
            spark.resolve_gene('A',self.gtf)
        self.assertEqual(spark.resolve_gene('C',self.gtf,0)[0],'chr1:2-6')

    def test_ambiguous_and_missing_genes_are_errors(self):
        self.gtf.write_text(self.gtf.read_text()+'chr2\tt\tgene\t1\t5\t.\t+\t.\tgene_id "ID_C"; gene_name "A";\n')
        with self.assertRaisesRegex(ValueError,'multiple loci'):
            spark.resolve_gene('A',self.gtf)
        with self.assertRaisesRegex(ValueError,'not found'):
            spark.resolve_gene('missing',self.gtf)

    def test_batch_regions_and_report(self):
        result = spark.plot_many(['chr1:0-10','chr1:10-20'],tracks={'Control':self.files[:1]},output=self.root/'batch')
        self.assertEqual((result['passed'],result['failed']),(2,0))
        self.assertTrue(Path(result['index']).is_file())
        for item in result['plots']:
            self.assertTrue(Path(item['outputs']['jpg']).is_file())
            record = json.loads(Path(item['outputs']['record']).read_text())
            spark.normalize_config(record['settings'])

    def test_bed_batch_cli(self):
        bed = self.root/'regions.bed';bed.write_text('chr1\t0\t10\tfirst\nchr1\t10\t20\tsecond\n')
        result = subprocess.run([sys.executable,str(ROOT/'SparK.py'),'--regions-file',str(bed),'-cf',self.files[0],'-o',str(self.root/'bed-batch'),'--no-update-check'],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertTrue((self.root/'bed-batch/first.jpg').is_file())

    def test_batch_gene_failures_are_reported_and_escaped(self):
        result = spark.plot_many(genes=['A','missing<&'],tracks={'Control':self.files[:1]},gtffile=str(self.gtf),flank=0,output=self.root/'genes')
        self.assertEqual((result['passed'],result['failed']),(1,1))
        page = Path(result['index']).read_text()
        self.assertIn('missing&lt;&amp;',page)
        self.assertNotIn('<script',page)

    def test_batch_shared_scales(self):
        result = spark.plot_many(['chr1:0-10','chr1:10-20'],tracks={'Control':self.files[:1]},output=self.root/'shared',batch_shared_scale=True)
        self.assertEqual(result['shared_scales'],{1:2})
        for item in result['plots']:
            labels = [e.text for e in ET.parse(item['outputs']['svg']).getroot().findall(NS+'text')]
            self.assertIn('2.2',labels)

    def test_batch_duplicate_names_are_unique(self):
        result = spark.plot_many([dict(name='same',region='chr1:0-10'),dict(name='same',region='chr1:10-20')],tracks={'Control':self.files[:1]},output=self.root/'duplicate')
        self.assertNotEqual(result['plots'][0]['outputs']['svg'],result['plots'][1]['outputs']['svg'])

    def test_inline_batch_config_round_trip(self):
        config = self.root/'batch.json'
        first = spark.plot_many(['chr1:0-10','chr1:10-20'],tracks={'Control':self.files[:1]},output=self.root/'batch-one',save_config=str(config))
        again = spark.plot(config=config,output=self.root/'batch-two')
        self.assertEqual(again['passed'],2)
        for a,b in zip(first['plots'],again['plots']):
            self.assertEqual(Path(a['outputs']['svg']).read_bytes(),Path(b['outputs']['svg']).read_bytes())

    def test_batch_scale_with_highest_group_excluded(self):
        result = spark.plot_many(['chr1:0-10','chr1:10-20'],output=self.root/'excluded',control_files=self.files[:2],control_groups=[1,2],exclude_groups=[2],batch_shared_scale=True)
        self.assertEqual(result['passed'],2)
        self.assertEqual(result['shared_scales'],{1:2})

    def test_cli_region_overrides_saved_gene_selector(self):
        config = self.root/'gene.json'
        config.write_text(json.dumps(dict(control_files=self.files[:1],gene='A',gtffile=str(self.gtf))))
        result = subprocess.run([sys.executable,str(ROOT/'SparK.py'),'--config',str(config),'-pr','chr1:0-20','-o',str(self.root/'override'),'--no-update-check','--no-jpg'],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        record = json.loads((self.root/'override.run.json').read_text())
        self.assertEqual(record['settings']['region'],'chr1:0-20')

    def test_invalid_numeric_config_exits_cleanly(self):
        for extra in [dict(smoothen='ten'),dict(exclude_groups=2),dict(scale_factors=['two']),dict(no_jpg='false')]:
            with self.subTest(extra=extra),self.assertRaises(ValueError):
                spark.plot(config=self.config(**extra))

    def test_release_versions_agree(self):
        version = spark.SparK_Version
        self.assertEqual((ROOT/'VERSION').read_text().strip(),version)
        self.assertIn('version = "'+version+'"',(ROOT/'pyproject.toml').read_text())
        self.assertTrue((ROOT/'README.md').read_text().startswith('# SparK '+version+'\n'))
        self.assertIn('## '+version+'\n',(ROOT/'CHANGELOG.md').read_text())

    def test_concurrent_python_plots_have_independent_state(self):
        def make(index):
            return spark.plot('chr1:0-20',tracks=[dict(file=self.files[index],condition='C',color=['#FF0000','#0000FF'][index])],output=self.root/('parallel'+str(index)),y_scale=index+1)
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            results = list(pool.map(make,[0,1]))
        for index,result in enumerate(results):
            polygon = ET.parse(result['outputs']['svg']).getroot().find(NS+'polygon')
            self.assertEqual(polygon.get('fill'),['#FF0000','#0000FF'][index])
            record = json.loads(Path(result['outputs']['record']).read_text())
            self.assertEqual(record['settings']['y_scale'],index+1)


class UpdateTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.cache = self.root/'cache/latest.json'
        self.env = patch.dict(os.environ,{'SPARK_NO_UPDATE_CHECK':'0'});self.env.start();self.addCleanup(self.env.stop)
        self.cache_patch = patch.object(spark,'_update_cache_path',return_value=self.cache);self.cache_patch.start();self.addCleanup(self.cache_patch.stop)

    def test_numeric_version_comparison(self):
        self.assertGreater(spark.version_tuple('3.10.0'),spark.version_tuple('3.9.0'))
        self.assertEqual(spark.version_tuple('v3.2.0'),(3,2,0))
        self.assertIsNone(spark.version_tuple('3.3.0-rc1'))
        self.assertIsNone(spark.version_tuple('hello\nInstall something'))

    def test_cached_new_version_notifies_only_after_finish(self):
        spark.atomic_json(self.cache,dict(checked_at=time.time(),latest_version='99.0.0'))
        stream = io.StringIO()
        with patch.object(spark.subprocess,'Popen') as launch:
            check = spark.UpdateCheck(stream=stream)
            self.assertEqual(stream.getvalue(),'')
            check.finish(True);check.finish(True)
            self.assertEqual(stream.getvalue().count('is available'),1)
            launch.assert_not_called()

    def test_same_older_and_failed_checks_are_silent(self):
        for latest in [spark.SparK_Version,'1.0.0',None]:
            spark.atomic_json(self.cache,dict(checked_at=time.time(),latest_version=latest))
            stream = io.StringIO();check = spark.UpdateCheck(stream=stream);check.finish(True)
            self.assertEqual(stream.getvalue(),'')

    def test_opt_out_does_not_start_network_worker(self):
        with patch.dict(os.environ,{'SPARK_NO_UPDATE_CHECK':'1'}),patch.object(spark.subprocess,'Popen') as launch:
            spark.UpdateCheck().finish(True)
            launch.assert_not_called()

    def test_worker_failure_does_not_raise_or_notify(self):
        work = self.root/'job';work.mkdir();spark.atomic_json(work/'done.json',dict(success=True))
        stream=io.StringIO()
        with patch.object(spark,'fetch_latest_version',side_effect=TimeoutError('offline')),patch.object(spark.sys,'stderr',stream):
            self.assertEqual(spark.update_worker(self.cache,work,spark.SparK_Version,True),0)
        self.assertEqual(stream.getvalue(),'')
        self.assertIsNone(json.loads(self.cache.read_text())['latest_version'])

    def test_notice_after_fast_parent_finishes(self):
        work = self.root/'job';work.mkdir()
        stream = io.StringIO()
        with patch.object(spark.subprocess,'Popen') as launch:
            check = spark.UpdateCheck(stream=stream)
            worker_dir = check.work
            check.finish(True)
            self.assertTrue((worker_dir/'done.json').is_file())
            self.assertEqual(stream.getvalue(),'')
        with patch.object(spark,'fetch_latest_version',return_value='99.0.0'),patch.object(spark.sys,'stderr',stream):
            spark.update_worker(self.cache,worker_dir,spark.SparK_Version,True)
        self.assertIn('99.0.0 is available',stream.getvalue())
        self.assertFalse(worker_dir.exists())

    def test_piped_stderr_is_not_kept_open_by_background_worker(self):
        with patch.object(spark.subprocess,'Popen') as launch:
            check = spark.UpdateCheck(stream=io.StringIO())
            self.assertEqual(launch.call_args.kwargs['stderr'],subprocess.DEVNULL)
            self.assertEqual(launch.call_args.kwargs['stdout'],subprocess.DEVNULL)
            self.assertTrue(launch.call_args.kwargs['close_fds'])
            spark.atomic_json(check.work/'result.json',dict(latest_version=None))
            check.finish(True)

    def test_fetch_accepts_only_published_stable_releases(self):
        for release,valid in [(dict(tag_name='v99.0.0',draft=False,prerelease=False),True),
                              (dict(tag_name='v99.0.0',draft=True),False),
                              (dict(tag_name='v99.0.0-rc1'),False)]:
            with self.subTest(release=release),patch.object(spark.urllib.request,'urlopen',return_value=io.BytesIO(json.dumps(release).encode())):
                if valid:
                    self.assertEqual(spark.fetch_latest_version(),'99.0.0')
                else:
                    with self.assertRaises(ValueError):spark.fetch_latest_version()

    @unittest.skipUnless(os.name=='posix','Terminal integration uses a POSIX pseudo-terminal')
    def test_slow_server_does_not_hold_up_plot_and_notice_follows_saved_output(self):
        import pty
        import select
        started,released = threading.Event(),threading.Event()
        class Server(BaseHTTPRequestHandler):
            def log_message(self,*args):pass
            def do_GET(self):
                started.set();released.wait(timeout=4)
                payload=b'{"tag_name":"v99.0.0","draft":false,"prerelease":false}'
                try:
                    self.send_response(200);self.send_header('Content-Length',str(len(payload)));self.end_headers();self.wfile.write(payload)
                except OSError:pass
        server=ThreadingHTTPServer(('127.0.0.1',0),Server)
        thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
        script=self.root/'SparK.py'
        source=(ROOT/'SparK.py').read_text().replace(spark.VERSION_URL,'http://127.0.0.1:{}/latest'.format(server.server_port),1)
        tree=ast.parse(source);node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='_update_cache_path')
        lines=source.splitlines(keepends=True)
        lines[node.lineno-1:node.end_lineno]=['def _update_cache_path():\n    return Path('+repr(str(self.cache))+')\n']
        script.write_text(''.join(lines))
        signal=self.root/'input.bdg';signal.write_text('chr1\t0\t10\t2\n')
        master,slave=pty.openpty()
        process=subprocess.Popen([sys.executable,str(script),'-cf',str(signal),'-pr','chr1:0-10','-o',str(self.root/'plot'),'--no-jpg'],
                                 stdin=subprocess.DEVNULL,stdout=slave,stderr=slave,env=dict(os.environ,SPARK_NO_UPDATE_CHECK='0'))
        os.close(slave)
        try:
            # The HTTP response remains withheld until the plotting process exits.
            self.assertEqual(process.wait(timeout=8),0)
            self.assertTrue(started.wait(timeout=5))
            self.assertTrue((self.root/'plot.svg').exists())
            released.set()
            output=b'';deadline=time.monotonic()+5
            while time.monotonic()<deadline and b'99.0.0 is available' not in output:
                if select.select([master],[],[],0.2)[0]:
                    try:output+=os.read(master,65536)
                    except OSError:break
            text=output.decode(errors='replace')
            self.assertIn('99.0.0 is available',text)
            self.assertLess(text.index('Saved '),text.index('99.0.0 is available'))
        finally:
            released.set()
            if process.poll() is None:process.kill();process.wait()
            os.close(master);server.shutdown();server.server_close();thread.join(timeout=3)


class GuiTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.server = spark.create_gui_server(data_dir=self.root)
        self.thread = threading.Thread(target=self.server.serve_forever,daemon=True);self.thread.start()
        self.base = 'http://127.0.0.1:'+str(self.server.server_port)
        self.addCleanup(self.stop)

    def stop(self):
        self.server.shutdown();self.server.server_close();self.thread.join(timeout=3)
        if self.server.demo_directory:
            shutil.rmtree(self.server.demo_directory,ignore_errors=True)

    def request(self,path,data=None,headers=None,token=True):
        settings={'X-SparK-Token':self.server.token} if token else {}
        if data is not None:settings['Content-Type']='application/json'
        settings.update(headers or {})
        request=urllib.request.Request(self.base+path,data=json.dumps(data).encode() if data is not None else None,headers=settings)
        return urllib.request.urlopen(request,timeout=10)

    def test_interface_is_local_and_token_protected(self):
        self.assertEqual(self.server.server_address[0],'127.0.0.1')
        with self.assertRaises(urllib.error.HTTPError) as error:self.request('/api/config',token=False)
        self.assertEqual(error.exception.code,403)
        with self.request('/') as response:
            self.assertIn(b'Preview &amp; save',response.read())

    def test_foreign_origin_and_path_escape_are_rejected(self):
        with self.assertRaises(urllib.error.HTTPError) as error:self.request('/api/config',headers={'Origin':'https://example.com'})
        self.assertEqual(error.exception.code,403)
        with self.assertRaises(urllib.error.HTTPError) as error:self.request('/api/files?path=..')
        self.assertEqual(error.exception.code,400)

    def test_demo_preview_and_downloads(self):
        with self.request('/api/demo') as response:config=json.load(response)
        with self.request('/api/plot',config) as response:result=json.load(response)
        self.assertIn('jpg',result['links']);self.assertIn('pdf',result['links'])
        with self.request(result['preview']) as response:
            self.assertEqual(response.headers['Content-Type'],'image/jpeg')
            self.assertTrue(response.read().startswith(b'\xff\xd8'))

    def test_gui_rejects_inline_batches(self):
        with self.request('/api/demo') as response:config=json.load(response)
        config.pop('region');config['regions']=['chr1:0-10']
        with self.assertRaises(urllib.error.HTTPError) as error:self.request('/api/plot',config)
        self.assertEqual(error.exception.code,400)
        self.assertIn('batch',error.exception.read().decode())

    def test_gui_outputs_stay_in_selected_output_directory(self):
        with self.request('/api/demo') as response:config=json.load(response)
        config['output_name']='../../outside'
        with self.request('/api/plot',config) as response:json.load(response)
        self.assertTrue((self.server.output_dir/'outside.svg').is_file())
        self.assertFalse((self.root.parent/'outside.svg').exists())


if __name__=='__main__':
    unittest.main()
