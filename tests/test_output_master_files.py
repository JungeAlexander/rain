import os
from bin.stringrnautils import Interaction
import unittest
import logging
logger = logging.getLogger(os.path.basename(__file__))
logging.basicConfig(level=logging.INFO)


class MyTestCase(unittest.TestCase):
    def setUp(self):
        logger.debug("Testing output.")
        master_dir = 'master_files'
        self.expected_master_files = [os.path.join(master_dir, fname)
                                      for fname in ['database.tsv', 'textmining.tsv', 'experiments.tsv',
                                                    'predictions.tsv']]
        self.expected_min_master_files_line_count = [267, 30000, 600000, 30000000]  # ordered as self.expected_master_files
        self.expected_master_files_to_line_count = {}
        for mf, count in zip(self.expected_master_files, self.expected_min_master_files_line_count):
            self.expected_master_files_to_line_count[mf] = count
        self.existing_master_files = [os.path.join(master_dir, fname)
                                      for fname in os.listdir(master_dir)
                                      if os.path.isfile(os.path.join(master_dir, fname))]
        self.aggregated_master_file = 'data/all.tsv'

    def tearDown(self):
        logger.debug('Done.' + os.linesep + os.linesep)

    def test_expected_master_files_existing(self):
        for exp_master_file in self.expected_master_files:
            self.assertTrue(os.path.isfile(exp_master_file), 'Master file {} does not exist.'.format(exp_master_file))

    def test_expected_master_files_have_five_lines(self):
        for exp_master_file in self.expected_master_files:
            num_lines = sum(1 for line in open(exp_master_file))
            self.assertGreaterEqual(num_lines, 5, 'Master file {} is too small.'.format(exp_master_file))

    def test_existance_aggregated_master_file(self):
        self.assertTrue(os.path.isfile(self.aggregated_master_file),
                        'Aggregated master file {} does not exist.'.format(self.aggregated_master_file))

    def test_scores_aggregated_master_file(self):
        with open(self.aggregated_master_file, 'r') as mf:
            for line in mf:
                org, ent1, ent2, directed, channel, score, sources, url, comment = line.rstrip('\n\r').split('\t')
                interact = Interaction(org, ent1, ent2, directed, channel, score, sources, url, comment)
                self.assertGreater(float(interact._score), 0,
                                   'Score of following interaction was not greater '
                                   'than {:f}: {}'.format(0, interact))
                self.assertLess(float(interact._score), 1,
                                'Score of following interaction was not less '
                                'than {:f}: {}'.format(1, interact))

    def test_size_aggregated_master_file(self):
        for exp_master_file in self.expected_master_files:
            found_num_lines = sum(1 for line in open(exp_master_file))
            exp_num_lines = self.expected_master_files_to_line_count[exp_master_file]
            self.assertGreaterEqual(found_num_lines, exp_num_lines,
                             'Master file {} is too small.'.format(exp_master_file))

if __name__ == '__main__':
    unittest.main()
