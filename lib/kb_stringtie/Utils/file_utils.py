import os
import time


def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


def _exchange_gene_ids(result_directory):
    """
    _exchange_gene_ids: exchange gene_ids with gene_name
    """

    log('starting exchanging gene_ids with gene_name')

    result_files = os.listdir(result_directory)

    for result_file_name in result_files:
        if result_file_name == 'transcripts.gtf':
            log('updating transcripts.gtf gene_ids')
            os.rename(os.path.join(result_directory, 'transcripts.gtf'),
                      os.path.join(result_directory, 'original_transcripts.gtf'))
            original_transcript_path = os.path.join(result_directory,
                                                    'original_transcripts.gtf')
            exchange_transcript_path = os.path.join(result_directory,
                                                    result_file_name)

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
        elif result_file_name == 't_data.ctab':
            log('updating t_data.ctab gene_ids')
            os.rename(os.path.join(result_directory, 't_data.ctab'),
                      os.path.join(result_directory, 'original_t_data.ctab'))
            original_tdata_path = os.path.join(result_directory,
                                               'original_t_data.ctab')
            exchange_tdata_path = os.path.join(result_directory,
                                               result_file_name)

            first_line = True
            with open(exchange_tdata_path, 'w') as output_file:
                with open(original_tdata_path, 'r') as input_file:
                    for line in input_file:
                        if first_line:
                            gene_id_index = line.split('\t').index('gene_id')
                            gene_name_index = line.split('\t').index('gene_name')
                            first_line = False
                            output_file.write(line)
                        else:
                            line_list = line.split('\t')
                            if len(line_list) >= max(gene_id_index, gene_name_index):
                                line_list[gene_id_index] = line_list[gene_name_index]
                                line = '\t'.join(line_list)
                                output_file.write(line)
                            else:
                                output_file.write(line)


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