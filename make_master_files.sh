#!/usr/bin/env bash
# Integrates and benchmarks RAIN interaction files.
# Author: Alexander Junge, RTH, University of Copenhagen, <alexander.junge@gmail.com>
set -o errexit
set -o nounset
# remove temporary files in case something goes wrong
trap 'rm -f master_files/* download_files/* data/all.tsv' ERR INT TERM

# set to true to skip benchmarking and integrating miRNA-target prediction tools and instead use a pre-computed set of
# predicted interactions and confidence scores
USE_PRECOMPUTED_PREDICTION_CHANNEL=false

# set to false to keep previously produced interaction master files, if they exist. Note that the textmining interaction
# master file is always produced, irrespective of the parameter below.
OVERWRITE_MASTER_FILES=true

# The RAIN version number to use for naming output files
VERSION=1

# set to true to use an extended gold standard consisting of Croft et al. data + high confidence interactions from miRTarBase/NPInter
USE_EXTENDED_GOLD_STD=true

master_file_dir=master_files/
data_file_dir=data/
id_dict_dir=id_dictionaries/
log_file=make_master_files.log

python bin/check_download_data.py 2>"$log_file"

# generate manually cleaned ncRNA alias file
Rscript bin/GenerateNcrnaAlias.R id_dictionaries/ManualAliasDict.tsv id_dictionaries/CircBase.tsv.gz \
id_dictionaries/non-coding_RNA.txt.gz id_dictionaries/ncRNAaliasfile.tsv.gz &>>"$log_file"

extended_gold_std_master_file="$data_file_dir"extended_gold_standard.tsv  # write to data dir to avoid messing with master file aggregation at the end

# check if all output master files exist in case they are to be overwritten
if [ -f "$master_file_dir"croft.tsv -a -f "$master_file_dir"predictions.tsv -a -f "$master_file_dir"experiments.tsv \
-a -f "$master_file_dir"database.tsv -a -f "$master_file_dir"textmining.tsv -a -f "$extended_gold_std_master_file" ]; then
    master_files_exist=true
else
    master_files_exist=false
fi

croft_master_file="$master_file_dir"croft.tsv
if [[ "$master_files_exist" = false || "$OVERWRITE_MASTER_FILES" = true ]]; then
    python bin/generate_ncRNAorth.py 2>>"$log_file"

    python bin/retrieve_knowledge_channel.py 2>>"$log_file"

    cat data/croft_mir_ensp.tsv | python bin/map_croft.py 1>"$croft_master_file" 2>>"$log_file"
    gold_standard_file="$croft_master_file"

    miRTarBase_master_file="$master_file_dir"miRTarBase_excluding_StarBase.tsv
    python bin/parse-miRTarBase.py "$miRTarBase_master_file" "$gold_standard_file" --exclude_starbase 2>>"$log_file"

    NPInter_SBexcluded_master_file="$master_file_dir"NPInter_excluding_StarBase.tsv
    python bin/parse-NPInter-database.py "$NPInter_SBexcluded_master_file" "$gold_standard_file" --exclude_starbase 2>>"$log_file"

    miRTarBase_NPInter_SBexcluded_master_file="$master_file_dir"miRTarBase_NPInter_SBexcluded.tsv
    python bin/combine-miRTarBase-NPInter.py "$miRTarBase_master_file" "$NPInter_SBexcluded_master_file" \
    "$miRTarBase_NPInter_SBexcluded_master_file" "$gold_standard_file" 2>>"$log_file"

    if [ "$USE_EXTENDED_GOLD_STD" = true ]; then

        extended_gold_std_master_file="$data_file_dir"extended_gold_standard.tsv  # write to data dir to avoid messing with master file aggregation at the end
        low_throughput_pids_file="$data_file_dir"low_throughoput_pids.tsv


        python bin/generate_extended_gold_std.py --human_only "$croft_master_file" "$miRTarBase_NPInter_SBexcluded_master_file" \
            --low-throughput-threshold 5 --n-low-throughput 2 --n-total 2 \
            --low-throughput-pids "$low_throughput_pids_file" \
            1>"$extended_gold_std_master_file" 2>> "$log_file"


        tmp_miRTarBase_NPInter_SBexcluded_master_file="$miRTarBase_NPInter_SBexcluded_master_file".tmp
        mv "$miRTarBase_NPInter_SBexcluded_master_file" "$tmp_miRTarBase_NPInter_SBexcluded_master_file"

        python bin/filter_extended_gold_std_rebenchmark_miRTarBase_NPInter.py  \
            "$extended_gold_std_master_file" \
            "$tmp_miRTarBase_NPInter_SBexcluded_master_file" \
            --low-throughput-pid-file "$low_throughput_pids_file" \
            1> "$miRTarBase_NPInter_SBexcluded_master_file" 2>>"$log_file"

        rm "$tmp_miRTarBase_NPInter_SBexcluded_master_file"

        # Replacing Croft et al. gold standard with extended gold standard..
        #rm "$croft_master_file"
        #mv "$extended_gold_std_master_file" "$master_file_dir"

        #gold_standard_file="$master_file_dir"extended_gold_standard.tsv
        gold_standard_file="$extended_gold_std_master_file"
    fi

    python bin/parse-starbase.py "$gold_standard_file" 2>>"$log_file"

    python bin/combine_experiments.py "$gold_standard_file" 2>>"$log_file"

    if [ "$USE_PRECOMPUTED_PREDICTION_CHANNEL" = true ]; then
        printf "Loading pre-computed predicted interaction file.\n\n" >>"$log_file"
        cp "$data_file_dir"/predictions_precomputed.tsv.gz "$master_file_dir"/predictions.tsv.gz
        gunzip "$master_file_dir"/predictions.tsv.gz
    else
        python bin/integratemiRNAPredictionTools.py "$gold_standard_file" -master_path "$master_file_dir" 2>>"$log_file"
    fi
else
    printf "Retaining previously computed master files.\n\n" >>"$log_file"
    if [ "$USE_EXTENDED_GOLD_STD" = true ]; then
        gold_standard_file="$extended_gold_std_master_file"
    else
        gold_standard_file="$croft_master_file"
    fi
fi

python bin/parse-textmining.py -disallow_mirpairs -disallow_sno_mir_pairs "$gold_standard_file" 2>>"$log_file"
# Keep all scores that have an integration score >= 0.15
grep -hv '#' "$master_file_dir"* | python bin/integrate_evidence_channels.py 0.15 > "$data_file_dir"all.tsv 2>>"$log_file"

min_interaction_count=500
# Filtering species such that only those with a minimum of $min_interaction_count are retained.
taxids="$(python bin/organism-statistics.py "$min_interaction_count" 2>>"$log_file")"
printf "Taxonomy IDs are kept: $taxids\n\n" >>"$log_file"

tax_arr=($taxids) # split at whitespace into array of tax ids
comb_file="$data_file_dir"combined.tsv
rm -f "$comb_file".gz
for tax in "${tax_arr[@]}"
do
    org_file="$data_file_dir$tax".combined.tsv.gz
    zcat data/integrated_scores.tsv.gz | grep "^$tax" | gzip > "$org_file"
    zcat "$org_file" >> "$comb_file"
done
gzip "$comb_file"

## getting statistics for  all available taxonomy identifiers.
#all_taxids="$(python bin/organism-statistics.py 1 2>/dev/null)"
#printf "Taxonomy IDs are kept: $all_taxids\n\n" >>"$log_file"

rna_aliases_file="$id_dict_dir"rna.aliases.search.interface.v"$VERSION".txt
rm -f "$rna_aliases_file".gz
python bin/create_rna_aliases_file.py "$data_file_dir"all.tsv $taxids 1>"$rna_aliases_file" 2>>"$log_file"
gzip "$rna_aliases_file"

# Aggregating all master files with no prior score filtering.
grep -hv '#' "$master_file_dir"* > "$data_file_dir"all_not_filtered.tsv

complete_rna_aliases_file="$id_dict_dir"rna.aliases.complete.v"$VERSION".txt
rm -f "$complete_rna_aliases_file".gz
python bin/create_rna_aliases_file.py "$data_file_dir"all_not_filtered.tsv $taxids 1>"$complete_rna_aliases_file" 2>>"$log_file"
gzip "$complete_rna_aliases_file"

## alias file containing all known aliases, i.e., the 'universe' of aliases
## this file is not of interest to RAIN itself but rather to other resources developed within RTH, e.g., WP7
universe_rna_aliases_file="$id_dict_dir"rna.aliases.universe.v"$VERSION".txt
rm -f "$universe_rna_aliases_file".gz
python bin/create_rna_aliases_file.py "$data_file_dir"all_not_filtered.tsv $taxids --print_universe 1>"$universe_rna_aliases_file" 2>>"$log_file"
gzip "$universe_rna_aliases_file"

python bin/find_cross_species_interactions.py "$complete_rna_aliases_file".gz "$master_file_dir"* 2>>"$log_file"

node_info_file=id_dictionaries/node_info.tsv
if [ -f "$node_info_file" ]
then
    printf "Loading pre-processed node information file.\n\n" >>"$log_file"
else
    python bin/parse-nodeinfo.py -alias "$complete_rna_aliases_file".gz -info "$node_info_file" 2>>"$log_file"
fi

python bin/create_payload.py -data_path "$data_file_dir" -payload_server "http://download.jensenlab.org/STRING-RNA/" \
-project_name "string-rna payload" -ncrna_info_file "$node_info_file" -tax_identifiers $taxids 2>>"$log_file"

python -m unittest discover tests/

./make_download_files.sh

echo "RAIN pipeline ran successfully."
exit 0
