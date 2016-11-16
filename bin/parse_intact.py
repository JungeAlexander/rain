#!/usr/bin/env python

__author__ = 'Alexander Junge (alexander.junge@gmail.com)'

import argparse
from ftplib import FTP
import os
import sys
import zipfile
import gzip

parser = argparse.ArgumentParser(description="""Extracts gold standard interactions from IntAct (http://www.ebi.ac.uk/intact/).

                                                Starts with downloading and extracting the complete database content.""")
args = parser.parse_args()

data_dir = 'data'
intact_version = '2015-08-21'
intact_interaction_file = os.path.join(data_dir, '_'.join(['intact', intact_version, 'psimitab', 'intact.zip']))
intact_interaction_dir = os.path.join(data_dir, '_'.join(['intact', intact_version, 'psimitab', 'intact']))

if not os.path.exists(intact_interaction_dir):
    if not os.path.exists(intact_interaction_file):
        # we want to download: ftp://ftp.ebi.ac.uk/pub/databases/intact/2015-08-21/psimitab/intact.zip
        ebi_ftp = 'ftp.ebi.ac.uk'
        intact_ftp_dir = 'pub/databases/intact/' + intact_version + '/psimitab/'
        intact_ftp_file = 'intact.zip'
        sys.stdout.write('Downloading %s%s from %s.\n' % (intact_ftp_dir, intact_ftp_file, ebi_ftp))

        ftp = FTP(ebi_ftp)
        ftp.login()
        ftp.cwd(intact_ftp_dir)

        with open(intact_interaction_file, 'wb') as out_file:
            ftp.retrbinary('RETR %s' % intact_ftp_file, out_file.write)
            sys.stdout.write('Download finished.\n')

        ftp.quit()
    else:
        sys.stdout.write('File %s exists, hence no Intact files are downloaded.\n' % intact_interaction_file)
    
    # extract previously downloaded intact zip archive and gzip
    sys.stdout.write('Extracting %s.\n' % intact_interaction_file)
    with zipfile.ZipFile(intact_interaction_file, 'r') as zfile:
            zfile.extractall(intact_interaction_dir)
    sys.stdout.write('Done extracting files.\n')

    files_to_gzip = [os.path.join(intact_interaction_dir, f) 
        for f in os.listdir(intact_interaction_dir) if os.path.isfile(os.path.join(intact_interaction_dir, f))]

    for src_path in files_to_gzip:
        sys.stdout.write('Gzipping %s.\n' % src_path)
        with open(src_path) as src, gzip.open(src_path + '.gzip', 'wb') as dst:        
            dst.writelines(src)
        os.remove(src_path)
        sys.stdout.write('Done gzipping.\n')

    os.remove(intact_interaction_file)
else:
    sys.stdout.write('Directory %s exists, hence no Intact files are changed.\n' % intact_interaction_dir)
