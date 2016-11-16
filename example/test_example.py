import os
import bin.stringrnautils as utils
import gzip
import unittest
import logging
import subprocess

logger = logging.getLogger(os.path.basename(__file__))
logging.basicConfig(level=logging.INFO)


class ExampleTestCase(unittest.TestCase):

    def setUp(self):
        logger.debug("Texting example.")
        example_dir = 'example'
        self.example_master_files = [os.path.join(example_dir, fname)
                                     for fname in ['experiments_example.tsv',
                                                   'predictions_example.tsv',
                                                   'textmining_example.tsv', ]]
        self.expected_integrated_score_file = os.path.join(example_dir,
                                                'integrated_score_expected.tsv')
        self.output_integrated_score_file = os.path.join(example_dir,
                                            'integrated_score_output.tsv.gz')
        subprocess.call('cat ' + ' '.join(self.example_master_files) +
                        '|  python bin/integrate_evidence_channels.py 0.15 '
                        '--STRING_prior 0.01 '
                        '--integrated_score_file_path {} >/dev/null'.format(
                         self.output_integrated_score_file), shell=True)

    def tearDown(self):
        if os.path.isfile(self.output_integrated_score_file):
            os.remove(self.output_integrated_score_file)
        logger.debug('Done.' + os.linesep + os.linesep)

    def test_expected_output_existing(self):
        self.assertTrue(os.path.isfile(self.expected_integrated_score_file),
                        'Expected integrated score file does not exist: ' +
                        self.expected_integrated_score_file)

    def test_expected_output_has_one_line(self):
        num_lines = sum(1 for line in open(self.expected_integrated_score_file))
        self.assertEqual(num_lines, 1,
                         'Expected integrated score file does not have '
                         'a single line but {:d}: {}'.format(num_lines,
                                        self.expected_integrated_score_file))

    def test_expected_output_content(self):
        with open(self.expected_integrated_score_file, 'r') as fin:
            expected_content = os.linesep.join(fin.readlines())
        with gzip.open(self.output_integrated_score_file, 'r') as fin:
            found_content = os.linesep.join(fin.readlines())
        self.assertMultiLineEqual(found_content, expected_content,
                'Found ({}) and expected ({}) integrated score files '
                'differ.'.format(self.output_integrated_score_file,
                                 self.expected_integrated_score_file))


if __name__ == '__main__':
    unittest.main()
