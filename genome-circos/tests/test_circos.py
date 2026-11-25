"""Run with python -I -m unittest discover -s tests."""
from collections import Counter
import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest
from xml.etree import ElementTree

SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts'
SPEC = importlib.util.spec_from_file_location('draw_circos', SCRIPTS / 'draw_circos.py')
PLOT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PLOT)


class CircosTests(unittest.TestCase):
    def test_window_boundary_and_repeat_union(self):
        with TemporaryDirectory() as path:
            genes, repeats = (Path(path) / key for key in ('genes.gff3', 'te.gff3'))
            genes.write_text('chr1\t.\tgene\t10\t10\t.\t+\t.\tID=g1\n'
                             'chr1\t.\tgene\t11\t11\t.\t+\t.\tID=g2\n'
                             '##FASTA\n>chr1\nACGT\n', encoding='utf-8')
            repeats.write_text('chr1\t.\trepeat_region\t8\t13\t.\t.\t.\tClassification=LTR/Gypsy\n'
                               'chr1\t.\trepeat_region\t10\t18\t.\t.\t.\tClassification=LTR/Gypsy\n', encoding='utf-8')
            value = PLOT.densities({'genes': genes, 'repeats': repeats}, {'chr1': 24}, 10)['chr1']
            self.assertEqual(value[0], [1, 1, 0])
            self.assertEqual(value[1], [3, 8, 0])

    def test_paf_chromosome_pair_span_and_strand(self):
        with TemporaryDirectory() as path:
            paf = Path(path) / 'demo.paf'
            paf.write_text('b\t100\t10\t21\t-\ta\t100\t0\t12\t10\t12\t60\ttp:A:P\n'
                           'b\t100\t10\t19\t+\ta\t100\t0\t11\t9\t11\t60\n'
                           'b\t100\t10\t21\t+\ta\t100\t0\t12\t10\t12\t0\ttp:A:S\n'
                           'b\t100\t10\t21\t+\tz\t100\t0\t12\t10\t12\t60\n', encoding='utf-8')
            rows = list(PLOT.paf_links(paf, [{'a': 100}, {'b': 100}], [('chr1', 'a', 'b')], 10))
            self.assertEqual(rows, [(0, 'a', 0, 12, 'b', 10, 21, '-')])
            with self.assertRaisesRegex(ValueError, 'lengths do not match'):
                list(PLOT.paf_links(paf, [{'a': 110}, {'b': 100}], [('chr1', 'a', 'b')], 10))

    def test_synthetic_cli_and_missing_optional_data(self):
        with TemporaryDirectory() as folder:
            demo, output, second = (Path(folder) / name for name in ('demo', 'output', 'second'))
            subprocess.run([sys.executable, '-I', str(SCRIPTS / 'make_demo.py'), str(demo)], check=True)
            config_path = demo / 'manifest.json'
            cmd = [sys.executable, '-I', str(SCRIPTS / 'draw_circos.py'), str(config_path)]
            first = subprocess.run(cmd + [str(output)], check=True, capture_output=True, text=True)
            self.assertIn('11 chromosome pairs, 44 alignment ribbons', first.stdout)
            labels = {item.text for item in ElementTree.parse(output / 'circos.svg').iter() if item.tag.endswith('}text')}
            self.assertIn('chr11', labels)
            self.assertIn('Second: Assembly B', labels)
            caption = (output / 'legend.en.md').read_text(encoding='utf-8')
            self.assertIn('they do not represent structural-variant calls', caption)
            with (output / 'window_values.tsv').open(encoding='utf-8') as table:
                sectors = Counter((r['assembly'], r['chromosome']) for r in csv.DictReader(table, delimiter='\t'))
            self.assertEqual(len(sectors), 22)
            self.assertNotEqual(subprocess.run(cmd + [str(output)], capture_output=True).returncode, 0)
            manifest = json.loads(config_path.read_text(encoding='utf-8'))
            del manifest['paf']
            for assembly in manifest['assemblies']:
                del assembly['telomeres'], assembly['rdna']
            config_path.write_text(json.dumps(manifest), encoding='utf-8')
            subprocess.run(cmd + [str(second)], check=True, capture_output=True)
            self.assertNotIn('Central ribbons', (second / 'legend.en.md').read_text(encoding='utf-8'))
            svg = (second / 'circos.svg').read_text(encoding='utf-8')
            self.assertEqual(svg.count(f'fill="{PLOT.TEL["not_assessed"]}"'), 44)
