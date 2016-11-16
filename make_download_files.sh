#!/usr/bin/env bash
# Generates files that are available for download on the RAIN website.
# Files are written to download_files/.
# Author: Alexander Junge, RTH, University of Copenhagen, <alexander.junge@gmail.com>
set -o nounset
set -o errexit

DOWNLOAD_DIR=download_files/
RAIN_VERSION=1

cd master_files/
# files per evidence channel
gzip *.tsv
cp croft.tsv.gz ../"$DOWNLOAD_DIR"/v"$RAIN_VERSION".croft.tsv.gz
cp predictions.tsv.gz ../"$DOWNLOAD_DIR"/v"$RAIN_VERSION".predictions.tsv.gz
cp experiments.tsv.gz ../"$DOWNLOAD_DIR"/v"$RAIN_VERSION".experiments.tsv.gz
cp database.tsv.gz ../"$DOWNLOAD_DIR"/v"$RAIN_VERSION".database.tsv.gz
cp textmining.tsv.gz ../"$DOWNLOAD_DIR"/v"$RAIN_VERSION".textmining.tsv.gz
gunzip *.tsv.gz
cd ..

# create species-specific download files (interactions and aliases)
cd "$DOWNLOAD_DIR"
taxids=( 9606 10090 10116 4932 )
channels=( predictions experiments database textmining )
for taxid in "${taxids[@]}"
do
  for chan in "${channels[@]}"
  do
    zcat v"$RAIN_VERSION"."$chan".tsv.gz | awk -v taxid="$taxid" '$1 == taxid {print $0}' | gzip > "$taxid".v"$RAIN_VERSION"."$chan".tsv.gz
  done
done
cd ..

cd data/
# files with combined scores
cp 10090.combined.tsv.gz ../"$DOWNLOAD_DIR"/10090.v"$RAIN_VERSION".combined.tsv.gz
cp 10116.combined.tsv.gz ../"$DOWNLOAD_DIR"/10116.v"$RAIN_VERSION".combined.tsv.gz
cp 4932.combined.tsv.gz ../"$DOWNLOAD_DIR"/4932.v"$RAIN_VERSION".combined.tsv.gz
cp 9606.combined.tsv.gz ../"$DOWNLOAD_DIR"/9606.v"$RAIN_VERSION".combined.tsv.gz
cp combined.tsv.gz ../"$DOWNLOAD_DIR"/v"$RAIN_VERSION".combined.tsv.gz
# extended gold std
gzip extended_gold_standard.tsv
cp extended_gold_standard.tsv.gz ../"$DOWNLOAD_DIR"/v"$RAIN_VERSION".extended.gold.standard.tsv.gz
gunzip extended_gold_standard.tsv.gz
cd ..

cd id_dictionaries/
# move alias file to webserver
cp rna.aliases.search.interface.v"$RAIN_VERSION".txt.gz ../"$DOWNLOAD_DIR"/v"$RAIN_VERSION".rna.aliases.search.interface.txt.gz
cp rna.aliases.complete.v"$RAIN_VERSION".txt.gz ../"$DOWNLOAD_DIR"/v"$RAIN_VERSION".rna.aliases.complete.txt.gz
cp rna.aliases.universe.v"$RAIN_VERSION".txt.gz ../"$DOWNLOAD_DIR"/v"$RAIN_VERSION".rna.aliases.universe.txt.gz
cd ..
