import csv
import datetime
import os
import time


def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


def _make_gff(file_path, append_file, novel_prefix='MSTRG.'):
    """Works on the very narrow case of gtf files from stringtie."""
    type_idx = 2
    written_genes = set()
    timestamp = str(datetime.datetime.now()).split('.')[0]

    if not os.path.isfile(file_path):
        raise ValueError('{} is not a file'.format(file_path))

    with open(append_file, 'a') as output_file:
        with open(file_path, 'r') as input_file:
            for line in input_file:
                if line[0] == "#":
                    continue
                sl = line.split('\t')
                gene_id = sl[-1].split('gene_id \"')[1].split('"')[0]
                transcript_id = line.split('transcript_id \"')[1].split('"')[0]
                if 'ref_gene_id \"' in sl[-1]:
                    gene_id = sl[-1].split('ref_gene_id \"')[1].split('"')[0]

                # write dummy genes
                elif sl[type_idx] == 'transcript' and gene_id not in written_genes:
                    sl[type_idx] = 'gene'
                    sl[-1] = 'ID={}; note="Spoofed gene for a RNASeq transcript\n'.format(gene_id)
                    output_file.write("\t".join(sl))
                    written_genes.add(gene_id)

                # write transcripts and exons
                if novel_prefix in transcript_id:
                    if sl[type_idx] == 'exon':
                        sl[-1] = "Parent={}".format(transcript_id)
                    elif gene_id:
                        sl[type_idx] = 'transcript'
                        sl[-1] = "ID={}; Parent={}".format(transcript_id, gene_id)
                    sl[-1] += "; note=Predicted transcript from RNASeq run on {}\n".format(
                        timestamp)
                    output_file.write("\t".join(sl))
    return append_file


def exchange_gene_ids(result_directory):
    """
    exchange_gene_ids: exchange gene_ids with gene_name
    """

    log('starting exchanging gene_ids with gene_name')

    result_files = os.listdir(result_directory)

    if 'transcripts.gtf' in result_files:
        pass
        _update_transcripts(result_directory)

    if 't_data.ctab' in result_files:
        _update_t_data(result_directory)


def _update_transcripts(result_directory):
    log('updating transcripts.gtf ')
    os.rename(os.path.join(result_directory, 'transcripts.gtf'),
              os.path.join(result_directory, 'original_transcripts.gtf'))
    original_transcript_path = os.path.join(result_directory,
                                            'original_transcripts.gtf')
    exchange_transcript_path = os.path.join(result_directory,
                                            'transcripts.gtf')
    with open(exchange_transcript_path, 'w') as output_file:
        with open(original_transcript_path, 'r') as input_file:
            for line in input_file:
                if 'transcript_id \"' in line and 'reference_id \"' in line:
                    trans_id = line.split('transcript_id \"')[1].split('"')[0]
                    ref_id = line.split('reference_id \"')[1].split('"')[0]
                    line = line.replace(trans_id, ref_id)
                if 'gene_id \"' in line and 'ref_gene_id \"' in line:
                    gene_id = line.split('gene_id \"')[1].split('"')[0]
                    ref_id = line.split('ref_gene_id \"')[1].split('"')[0]
                    line = line.replace(gene_id, ref_id)

                output_file.write(line)


def _update_t_data(result_directory):
    log('updating t_data.ctab gene_ids')
    os.rename(os.path.join(result_directory, 't_data.ctab'),
              os.path.join(result_directory, 'original_t_data.ctab'))
    original_tdata_path = os.path.join(result_directory,
                                       'original_t_data.ctab')
    exchange_tdata_path = os.path.join(result_directory,
                                       't_data.ctab')
    reader = csv.DictReader(open(original_tdata_path), dialect='excel-tab')
    writer = csv.DictWriter(open(exchange_tdata_path, 'w'), reader.fieldnames,
                            dialect='excel-tab')
    writer.writeheader()
    for line in reader:
        if "gene_name" in line:
            line['gene_id'] = line['gene_name']
        writer.writerow(line)


def _filter_merge_file(gtf_file):
    """
    _filter_merge_file: remove lines with no gene_name
    """

    log('start filtering merged gft file')

    dir_name = os.path.dirname(gtf_file)

    filtered_file_name = 'filtered_stringtie_merge.gtf'
    filtered_file_path = os.path.join(dir_name, filtered_file_name)

    with open(filtered_file_path, 'w') as output_file:
        with open(gtf_file, 'r') as input_file:
            for line in input_file:
                if line.startswith('#') or 'gene_name' in line:
                    output_file.write(line)
                else:
                    log('skipping line:\n{}'.format(line))

    return filtered_file_path
