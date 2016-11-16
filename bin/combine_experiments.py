import stringrnautils
import argparse
import os
import logging
logger = logging.getLogger(os.path.basename(__file__))
logging.basicConfig(level=logging.INFO)

def is_valid_path(arg):
    if not os.path.exists(arg):
        parser.error("Given master file path %s does not exist." % arg)
    return arg


def main():
    stringrnautils.combine_masterfiles(("miRTarBase_NPInter_SBexcluded.tsv", "starbase.tsv"),
                                        'experiments.tsv', args.gold_standard_file, 'experiments', 40, 
                                         negative_evidence=False,rebenchmark_everything=True,
                                         ignore_fraction=0.3)


if __name__ == '__main__':
    logger.info("Integrating Experiments channel.")
    parser = argparse.ArgumentParser()
    parser.add_argument('gold_standard_file', type=is_valid_path,
                        help='The gold standard file to benchmark against.')
    args = parser.parse_args()
    main()
    logger.info('Done.' + os.linesep + os.linesep)

