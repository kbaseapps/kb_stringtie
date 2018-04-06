import csv
import datetime
import os
import time


def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


def exchange_gene_ids(result_directory, trans_id_dict=None):
    """
    exchange_gene_ids: exchange gene_ids with gene_name
    """

    log('starting exchanging gene_ids with gene_name')

    result_files = os.listdir(result_directory)

    if 'transcripts.gtf' in result_files:
        _update_transcripts(result_directory)

    if 't_data.ctab' in result_files:
        _update_t_data(result_directory, trans_id_dict)


def _make_gff(file_path):
    """Works on the very narrow case of gtf files from stringtie.
       Makes stringtie transcript IDs into timestamp IDs"""
    type_idx = 2
    timestamp = datetime.datetime.now().isoformat().split('.')[0]
    gene_id_dict = {}
    trans_id_dict = {}
    written_genes = set()

    if not os.path.isfile(file_path):
        raise ValueError('{} is not a file'.format(file_path))
    new_file_path = os.path.splitext(file_path)[0] + ".gff"
    if new_file_path == file_path:
        raise ValueError('{} appears to be a GFF file'.format(file_path))

    with open(new_file_path, 'w') as output_file:
        with open(file_path, 'r') as input_file:
            # On first pass, just collect gene IDs and write dummy genes
            for line in input_file:
                if line[0] == "#":
                    continue
                if 'gene_id \"' in line and 'ref_gene_id \"' in line:
                    gene_id = line.split('gene_id \"')[1].split('"')[0]
                    ref_gene = line.split('ref_gene_id \"')[1].split('"')[0]
                    gene_id_dict[gene_id] = ref_gene
                    if ref_gene not in written_genes:
                        sl = line.split('\t')
                        sl[type_idx] = 'gene'
                        sl[-1] = 'ID={}\n'.format(ref_gene)
                        output_file.write("\t".join(sl))
                        written_genes.add(ref_gene)

            # now we write the transcripts.
            input_file.seek(0)
            for line in input_file:
                if line[0] == "#":
                    continue
                sl = line.split('\t')
                gene_id = line.split('gene_id \"')[1].split('"')[0]
                gene_id = gene_id_dict.get(gene_id)
                transcript_id = line.split('transcript_id \"')[1].split('"')[0]
                if "MSTRG." in transcript_id:
                    if transcript_id not in trans_id_dict:
                        trans_id_dict[transcript_id] = "_".join(
                            [timestamp, str(len(trans_id_dict)+1)])
                    transcript_id = trans_id_dict[transcript_id]

                    if sl[type_idx] == 'exon':
                        sl[-1] = "Parent={}".format(transcript_id)
                    elif gene_id:
                        sl[-1] = "ID={}; Parent={}".format(transcript_id, gene_id)
                    else:
                        sl[-1] = "ID={}".format(transcript_id)
                    output_file.write("\t".join(sl+['\n']))
    return trans_id_dict, new_file_path


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
                if 'gene_id \"' in line and 'ref_gene_name \"' in line:
                    gene_id = line.split('gene_id \"')[1].split('"')[0]
                    gene_name = line.split('ref_gene_name \"')[1].split('"')[0]
                    line = line.replace(gene_id, gene_name)
                    output_file.write(line)
                elif 'gene_id \"' in line and 'gene_name \"' in line:
                    gene_id = line.split('gene_id \"')[1].split('"')[0]
                    gene_name = line.split('gene_name \"')[1].split('"')[0]
                    line = line.replace(gene_id, gene_name)
                    output_file.write(line)
                else:
                    output_file.write(line)


def _update_t_data(result_directory, trans_id_dict=None):
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
        if trans_id_dict and line.get('t_name') in trans_id_dict:
            line['t_name'] = trans_id_dict[line['t_name']]
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
