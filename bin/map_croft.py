#!/usr/bin/env python
import sys
import os
import stringrnautils
import logging

logger = logging.getLogger(os.path.basename(__file__))
logging.basicConfig(level=logging.INFO)

# the file uses string8 identifiers
# need to map to string10
# ensp8 -> ensg -> ensp10


def main():
    ############################################################################
    # create string mapping dictionaries
    ############################################################################
    ensp8_to_ensg = stringrnautils.get_string_to_alias_mapper('9606', 'ENSP', 'ENSG', 8)['9606']
    ensg_to_ensp10 = stringrnautils.get_alias_to_string_mapper('9606', 'ENSP', 'ENSG', 10)['9606']
    ensp_to_ensp10 = stringrnautils.get_alias_to_string_mapper('9606', 'ENSP', 'ENSP', 10)['9606']
    mir_hash = stringrnautils.get_unique_mir_mapper()

    # cheating :) because one gene died between versions
    # ensp8_to_ensg['ENSP00000308970'] = 'ENSG00000231924'
    # ensg_to_ensp9['ENSG00000231924'] = 'ENSP00000244296'
    black_list = set()  # croft contains duplicates!!
    unmappable_mir_triples = set()
    unmappable_ensp_triples = set()

    for line in sys.stdin:
        mir, ensp, _ = line.rstrip().split('\t')

        new_mir, new_ensp = None, None

        # Map miRNA
        if mir in mir_hash:
            new_mir = mir_hash[mir]
        else:
            unmappable_mir_triples.add('({}, {})'.format(ensp, mir))
            # if mir not in mir_hash and mir in mir_alias_hash:
            #     sys.stderr.write(' - mir is one of these %s' % str(mir_alias_hash[mir]))

        # Map ENSP
        if ensp in ensp8_to_ensg and ensp8_to_ensg[ensp] in ensg_to_ensp10:
            ensg = ensp8_to_ensg[ensp]
            new_ensp = ensg_to_ensp10[ensg]
        elif ensp in ensp_to_ensp10:
            new_ensp = ensp_to_ensp10[ensp]
        else:
            unmappable_ensp_triples.add('({}, {})'.format(ensp, mir))

        if new_ensp and new_mir and (ensp, mir) not in black_list:
            out_line = '\t'.join(("9606", new_mir, new_ensp, "0", "database", "0.900", "Croft", "", ""))
            sys.stdout.write('%s\n' % out_line)
            black_list.add((ensp, mir))

    if len(unmappable_mir_triples) > 0:
        logger.warning("Could not map " + str(len(unmappable_mir_triples)) +
                       " miRNAs to IDs used in miRBase. Respective interactions were: " +
                       ', '.join(unmappable_mir_triples))

    if len(unmappable_ensp_triples) > 0:
        logger.warning("Could not map " + str(len(unmappable_ensp_triples)) +
                       " proteins to ENSPs used in STRING 10. Respective interactions were: " +
                       ', '.join(unmappable_ensp_triples))

if __name__ == '__main__':
    if len(sys.argv) != 1:
        raise ValueError('Run this script without arguments.')
    else:
        logger.info('Generating Croft gold standard set.')
        main()
        logger.info('Done.' + os.linesep + os.linesep)
