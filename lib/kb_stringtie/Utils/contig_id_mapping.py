#!/usr/bin/env python

import sys
import os
import re
import shutil

def create_sanitized_contig_ids(fasta_filename=None):
    """Create a tab delimited file containing a column for original contig ids
    and a column for sanitized contig ids that only contain alphanumeric characters"""

    if fasta_filename is None or not os.path.isfile(fasta_filename):
        raise IOError("Invalid FASTA file given: {}".format(fasta_filename))

    contig_ids = []
    # key = modified_id, value = original contig_id
    contig_id_mapping = {}

    # save all contig_ids
    with open(fasta_filename, 'r') as data:
        line = data.readline()
        while line:
            if line.startswith(">"):
                contig_ids.append(line.split()[0][1:])

            line = data.readline()

    for x in contig_ids:
        modified_id = re.sub('[^0-9a-zA-Z]+', '', x)

        while modified_id in contig_id_mapping:
            modified_id = modified_id + 'a'

        contig_id_mapping[modified_id] = x

    id_filename = fasta_filename.split('.')[0] + "_mapping.tab"
    with open(id_filename, 'w') as id_file:
        # write the header
        id_file.write("original\tmodified\n")

        # write the ids
        for x in contig_id_mapping:
            id_file.write("{}\t{}\n".format(contig_id_mapping[x], x)) 

    return id_filename

def _parse_mapping(mapping_filename, to_modified=True):
    contig_id_mapping = {}
    with open(mapping_filename, 'r') as id_mapping:
        header = id_mapping.readline()
        line = id_mapping.readline()
        while line:
            original_id, modified_id = line.strip().split('\t')
            if to_modified:
                from_id = original_id
                to_id = modified_id
            else:
                from_id = modified_id
                to_id = original_id

            contig_id_mapping[from_id] = to_id
            line = id_mapping.readline()
    return contig_id_mapping

def replace_fasta_contig_ids(fasta_filename=None, mapping_filename=None, to_modified=True):
    """Replace FASTA contig id strings with modified id strings (to_modified=True) or vice versa (to_modified=False) 
    from a tab delimited file created by create_sanitized_contig_ids()"""

    if fasta_filename is None or not os.path.isfile(fasta_filename):
        raise IOError("Invalid FASTA file given: {}".format(fasta_filename))

    # key = from, value = to
    contig_id_mapping = _parse_mapping(mapping_filename, to_modified)

    # create temp file for new contents
    temp_name = fasta_filename.split('.')[0] + "_temp.fa"
    temp_fasta = open(temp_name, 'w')

    with open(fasta_filename, 'r') as fasta_data:
        line = fasta_data.readline()
        while line:
            if line.startswith(">"):
                contig_id = line.split()[0][1:]
                modified_line = line.replace(contig_id, contig_id_mapping[contig_id])
                temp_fasta.write(modified_line)
            else:
                temp_fasta.write(line)
                
            line = fasta_data.readline()
    
    temp_fasta.close()
    shutil.copyfile(temp_name, fasta_filename)

def replace_gff_contig_ids(gff_filename=None, mapping_filename=None, to_modified=True):
    """Replace FASTA contig id strings with modified id strings (to_modified=True) or vice versa (to_modified=False) 
    from a tab delimited file created by create_sanitized_contig_ids()"""
    
    # key = from, value = to
    contig_id_mapping = _parse_mapping(mapping_filename, to_modified)

    # create temp file for new contents
    temp_name = gff_filename.split('.')[0] + "_temp.gff"
    temp_gff = open(temp_name, 'w')

    with open(gff_filename, 'r') as gff_data:
        line = gff_data.readline()
        while line:
            if line.startswith("##sequence-region"):
                contig_id = line.split()[1]
                modified_line = line.replace(contig_id, contig_id_mapping[contig_id])
                temp_gff.write(modified_line)
            elif line.startswith("#"):
                temp_gff.write(line)
            else:
                contig_id = line.split()[0]
                #modified_line = line.replace(contig_id, contig_id_mapping[contig_id])
                #temp_gff.write(modified_line)
                temp_gff.write(line)
                
            line = gff_data.readline()
    
    temp_gff.close()
    shutil.copyfile(temp_name, gff_filename)

if __name__ == "__main__":
    #fasta_filename = "test_fasta_original.fa"
    #gff_filename = "test_gff_original.gff"

    fasta_filename = "t/data/contig_id_mapping/arab_fasta.fasta"
    gff_filename = "t/data/contig_id_mapping/arab_gff.gff"

    #mapping_filename = create_sanitized_contig_ids(fasta_filename)
    #replace_fasta_contig_ids(fasta_filename, mapping_filename, to_modified=True)

    #print "Modified FASTA file"
    #with open(fasta_filename, 'r') as f:
    #    print f.read()

    #print "*" * 80

    #replace_fasta_contig_ids(fasta_filename, mapping_filename, to_modified=False)

    #print "Original FASTA file"
    #with open(fasta_filename, 'r') as f:
    #    print f.read()  
    #print "Modified GFF file"
    #replace_gff_contig_ids(gff_filename, mapping_filename, to_modified=True)

    #with open(gff_filename, 'r') as f:
    #    print f.read()

    #print "*" * 80

    #print "Original GFF file"
    #replace_gff_contig_ids(gff_filename, mapping_filename, to_modified=False)

    #with open(gff_filename, 'r') as f:
    #    print f.read()
