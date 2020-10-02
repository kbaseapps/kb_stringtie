import errno
import json
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import time
import traceback
import uuid
import zipfile

from pathos.multiprocessing import ProcessingPool as Pool

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.ExpressionUtilsClient import ExpressionUtils
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.ReadsAlignmentUtilsClient import ReadsAlignmentUtils
from installed_clients.SetAPIServiceClient import SetAPI
from installed_clients.WorkspaceClient import Workspace as Workspace
from .file_utils import exchange_gene_ids, _update_merge_file, _make_gff


def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class StringTieUtil:
    STRINGTIE_TOOLKIT_PATH = '/kb/deployment/bin/StringTie'
    GFFREAD_TOOLKIT_PATH = '/kb/deployment/bin/gffread'
    GFFCOMPARE_TOOLKIT_PATH = '/kb/deployment/bin/gffcompare'

    OPTIONS_MAP = {'output_transcripts': '-o',
                   'gene_abundances_file': '-A',
                   'num_threads': '-p',
                   'fr_firststrand': '--rf',
                   'fr_secondstrand': '--fr',
                   'cov_refs_file': '-C',
                   'junction_base': '-a',
                   'junction_coverage': '-j',
                   'disable_trimming': '-t',
                   'min_locus_gap_sep_value': '-g',
                   'ballgown_mode': '-B',
                   'skip_reads_with_no_ref': '-e',
                   'maximum_fraction': '-M',
                   'label': '-l',
                   'gtf_file': '-G',
                   'min_length': '-m',
                   'min_read_coverage': '-c',
                   'min_isoform_abundance': '-f'}

    BOOLEAN_OPTIONS = ['disable_trimming', 'ballgown_mode', 'skip_reads_with_no_ref']

    def _validate_run_stringtie_params(self, params):
        """
        _validate_run_stringtie_params:
                validates params passed to run_stringtie method
        """

        log('start validating run_stringtie params')

        # check for required parameters
        for p in ['alignment_object_ref', 'workspace_name', 'expression_suffix', 
                  'expression_set_suffix']:
            if p not in params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

    def _mkdir_p(self, path):
        """
        _mkdir_p: make directory for given path
        """
        if not path:
            return
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    def _generate_command(self, params):
        """
        _generate_command: generate stringtie command
        """

        command = self.STRINGTIE_TOOLKIT_PATH + '/stringtie '

        for key, option in self.OPTIONS_MAP.items():
            option_value = params.get(key)
            if key in self.BOOLEAN_OPTIONS and option_value:
                option_value = ' '
            if option_value:
                command += '{} {} '.format(option, option_value)

        command += '{} '.format(params.get('input_file'))

        log('generated stringtie command: {}'.format(command))

        return command

    def _run_command(self, command):
        """
        _run_command: run command and print result
        """

        log('start executing command:\n{}'.format(command))

        pipe = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        output = pipe.communicate()[0]
        exitCode = pipe.returncode

        if (exitCode == 0):
            log('Executed command:\n{}\n'.format(command) +
                'Exit Code: {}\nOutput:\n{}'.format(exitCode, output))
        else:
            error_msg = 'Error running command:\n{}\n'.format(command)
            error_msg += 'Exit Code: {}\nOutput:\n{}'.format(exitCode, output)
            raise ValueError(error_msg)

    def _run_gffread(self, gff_path, gtf_path):
        """
        _run_gffread: run gffread script

        ref: http://ccb.jhu.edu/software/stringtie/gff.shtml
        """

        log('converting gff to gtf')

        command = self.GFFREAD_TOOLKIT_PATH + '/gffread '
        command += "-E {0} -T -o {1}".format(gff_path, gtf_path)

        self._run_command(command)

    def _run_gffcompare(self, gff_path, gtf_path):
        """
        _run_gffcompare: run gffcompare script

        ref: http://ccb.jhu.edu/software/stringtie/gff.shtml
        """

        log('converting gff to gtf')
        output = os.path.dirname(gtf_path) + "/gffcmp"

        command = self.GFFCOMPARE_TOOLKIT_PATH + '/gffcompare '
        command += "-r {} -G -o {} {}".format(gff_path, output, gtf_path)

        self._run_command(command)

    def _get_input_file(self, alignment_ref):
        """
        _get_input_file: get input  SAM/BAM file from Alignment object
        """

        log('getting bam file from alignment')

        bam_file_dir = self.rau.download_alignment({'source_ref': 
                                                    alignment_ref})['destination_dir']

        files = os.listdir(bam_file_dir)
        bam_file_list = [file for file in files if re.match(r'.*\_sorted\.bam', file)]
        if not bam_file_list:
            bam_file_list = [file for file in files if re.match(r'.*(?<!sorted)\.bam', file)]

        if not bam_file_list:
            raise ValueError('Cannot find .bam file from alignment {}'.format(alignment_ref))

        bam_file_name = bam_file_list[0]

        bam_file = os.path.join(bam_file_dir, bam_file_name)

        return bam_file

    def _get_gtf_file(self, alignment_ref, result_directory):
        """
        _get_gtf_file: get the reference annotation file (in GTF format)
        """

        alignment_data = self.ws.get_objects2({'objects':
                                               [{'ref': alignment_ref}]})['data'][0]['data']

        genome_ref = alignment_data.get('genome_id')

        # annotation_file = self._create_gtf_file(genome_ref, result_directory)
        annotation_file = self._create_gtf_annotation_from_genome(genome_ref, result_directory)

        gene_name_annotation_file = annotation_file.split('.gtf')[0] + '_append_name.gtf'

        with open(gene_name_annotation_file, 'w') as output_file:
            with open(annotation_file, 'r') as input_file:
                for line in input_file:
                    if ('gene_id \"' in line) and ('gene_name \"' not in line):
                        line = line.replace("\n", "")
                        gene_id = line.split('gene_id \"')[1].split('"')[0]
                        line += ' gene_name \"{}\";\n'.format(gene_id)
                        output_file.write(line)
                    else:
                        output_file.write(line)

        return gene_name_annotation_file

    def _create_gtf_annotation_from_genome(self, genome_ref, result_directory):
        """
         Create reference annotation file from genome
        """
        ref = self.ws.get_object_subset(
            [{'ref': genome_ref, 'included': ['contigset_ref', 'assembly_ref']}])
        contig_id = None
        if 'contigset_ref' in ref[0]['data']:
            contig_id = ref[0]['data']['contigset_ref']
        elif 'assembly_ref' in ref[0]['data']:
            contig_id = ref[0]['data']['assembly_ref']
        if contig_id is None:
            raise ValueError(
                "Genome at {0} does not have reference to the assembly object".format(
                    genome_ref))
        print(contig_id)
        log("Generating GFF file from Genome")
        try:
            ret = self.au.get_assembly_as_fasta({'ref': genome_ref + ";" + contig_id})
            fa_output_file = ret['path']

            if os.path.dirname(fa_output_file) != result_directory:
                shutil.copy(fa_output_file, result_directory)

            # get the GFF
            ret = self.gfu.genome_to_gff({'genome_ref': genome_ref,
                                          'target_dir': result_directory})
            genome_gff_file = ret['file_path']
            gtf_ext = ".gtf"

            if not genome_gff_file.endswith(gtf_ext):
                gtf_path = os.path.splitext(genome_gff_file)[0] + '.gtf'
                self._run_gffread(genome_gff_file, gtf_path)
            else:
                gtf_path = genome_gff_file

            log("gtf file : " + gtf_path)
        except Exception:
            raise ValueError(
                "Generating GTF file from Genome Annotation object Failed :  {}".format(
                    "".join(traceback.format_exc())))

        return gtf_path

    def _create_gtf_file(self, genome_ref, result_directory):
        """
        _create_gtf_file: create reference annotation file from genome
        """

        log('start generating reference annotation file')

        genome_gtf_file = self.gfu.genome_to_gff({'genome_ref': genome_ref,
                                                  'target_dir': result_directory,
                                                  'is_gtf': True})['file_path']

        return genome_gtf_file

    def _save_expression(self, result_directory, alignment_ref, workspace_name,
                         expression_suffix, genome_ref='', transcripts=0):
        """
        _save_expression: save Expression object to workspace
        """

        log('start saving Expression object')

        alignment_data_object = self.ws.get_objects2({'objects':
                                                     [{'ref': alignment_ref}]})['data'][0]

        alignment_name = alignment_data_object['info'][1]
        if re.match('.*_*[Aa]lignment', alignment_name):
            expression_obj_name = re.sub('_*[Aa]lignment', expression_suffix, alignment_name)
        else:
            expression_obj_name = alignment_name + expression_suffix

        destination_ref = workspace_name + '/' + expression_obj_name
        upload_expression_params = {'destination_ref': destination_ref,
                                    'source_dir': result_directory,
                                    'alignment_ref': alignment_ref,
                                    'tool_used': 'StringTie',
                                    'tool_version': '1.3.3',
                                    'genome_ref': genome_ref,
                                    'transcripts': transcripts
                                    }

        expression_ref = self.eu.upload_expression(upload_expression_params)['obj_ref']

        return expression_ref

    def _save_expression_set(self, alignment_expression_map, alignment_set_ref,
                             workspace_name, expression_set_suffix, genome_ref=None):
        """
        _save_expression_set: save ExpressionSet object to workspace
        """

        log('start saving ExpressionSet object')

        items = []
        for alignment_expression in alignment_expression_map:
            items.append({'ref': alignment_expression.get('expression_obj_ref'),
                          'label': alignment_expression.get('alignment_label')})

        expression_set_data = {'description': 'ExpressionSet using StringTie',
                               'items': items}

        alignment_set_data_object = self.ws.get_objects2({'objects':
                                                         [{'ref': alignment_set_ref}]})['data'][0]

        alignment_set_name = alignment_set_data_object['info'][1]
        if re.match('.*_*[Aa]lignment_*[Ss]et', alignment_set_name):
            expression_set_name = re.sub('_*[Aa]lignment_*[Ss]et',
                                         expression_set_suffix,
                                         alignment_set_name)
        else:
            expression_set_name = alignment_set_name + expression_set_suffix

        expression_set_save_params = {'data': expression_set_data,
                                      'workspace': workspace_name,
                                      'genome_ref': genome_ref,
                                      'output_object_name': expression_set_name}

        save_result = self.set_client.save_expression_set_v1(expression_set_save_params)
        expression_set_ref = save_result['set_ref']

        return expression_set_ref

    def _save_expression_matrix(self, expressionset_ref, workspace_name):
        """
        _save_expression_matrix: save FPKM and TPM ExpressionMatrix
        """

        log('start saving ExpressionMatrix object')

        expression_set_name = self.ws.get_object_info([{"ref": expressionset_ref}],
                                                      includeMetadata=None)[0][1]

        output_obj_name_prefix = re.sub('_*[Ee]xpression_*[Ss]et',
                                        '',
                                        expression_set_name)

        upload_expression_matrix_params = {'expressionset_ref': expressionset_ref,
                                           'output_obj_name': output_obj_name_prefix,
                                           'workspace_name': workspace_name}

        expression_matrix_refs = self.eu.get_expressionMatrix(upload_expression_matrix_params)

        return expression_matrix_refs

    def _generate_output_file_list(self, result_directory):
        """
        _generate_output_file_list: zip result files and generate file_links for report
        """

        log('start packing result files')

        output_files = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(output_directory)
        result_file = os.path.join(output_directory, 'stringtie_result.zip')

        with zipfile.ZipFile(result_file, 'w',
                             zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zip_file:
            for root, dirs, files in os.walk(result_directory):
                for file in files:
                    if not file.endswith('.DS_Store'):
                        zip_file.write(os.path.join(root, file),
                                       os.path.join(os.path.basename(root), file))

        output_files.append({'path': result_file,
                             'name': os.path.basename(result_file),
                             'label': os.path.basename(result_file),
                             'description': 'File(s) generated by StringTie App'})

        result_dirs = os.listdir(result_directory)
        if 'merge_result' in result_dirs:
            merge_file = os.path.join(result_directory, 'merge_result', 'stringtie_merge.gtf')
            output_files.append({'path': merge_file,
                                 'name': os.path.basename(merge_file),
                                 'label': os.path.basename(merge_file),
                                 'description': 'merge file generated by StringTie App'})

        return output_files

    def _generate_merge_html_report(self, result_directory):
        """
        _generate_html_report: generate html summary report
        """

        log('start generating merge html report')
        html_report = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(output_directory)
        result_file_path = os.path.join(output_directory, 'report.html')

        result_dirs = os.listdir(result_directory)

        Overview_Content = ''
        Overview_Content += '<br/><table><tr><th>Generated Files</th>'
        Overview_Content += '<th></th></tr>'
        Overview_Content += '<tr><th>Directory</th><th>File Name</th></tr>'
        for result_dir in result_dirs:
            result_files = os.listdir(os.path.join(result_directory, result_dir))
            result_files.sort()
            first_file = True
            for file_name in result_files:
                if first_file:
                    Overview_Content += '<tr><td>{}</td>'.format(result_dir)
                    Overview_Content += '<td>{}</td></tr>'.format(file_name)
                    first_file = False
                else:
                    Overview_Content += '<tr><td>{}</td>'.format('')
                    Overview_Content += '<td>{}</td></tr>'.format(file_name)
        Overview_Content += '</table>'

        with open(result_file_path, 'w') as result_file:
            with open(os.path.join(os.path.dirname(__file__), 'report_template.html'),
                      'r') as report_template_file:
                report_template = report_template_file.read()
                report_template = report_template.replace('<p>Overview_Content</p>',
                                                          Overview_Content)
                result_file.write(report_template)

        html_report.append({'path': result_file_path,
                            'name': os.path.basename(result_file_path),
                            'label': os.path.basename(result_file_path),
                            'description': 'HTML summary report for StringTie App'})
        return html_report

    def _generate_html_report(self, result_directory, obj_ref):
        """
        _generate_html_report: generate html summary report
        """

        log('start generating html report')
        html_report = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(output_directory)
        result_file_path = os.path.join(output_directory, 'report.html')

        expression_object = self.ws.get_objects2({'objects':
                                                 [{'ref': obj_ref}]})['data'][0]
        expression_info = expression_object['info']
        expression_data = expression_object['data']

        expression_object_type = expression_info[2]
        Overview_Content = ''
        if re.match('KBaseRNASeq.RNASeqExpression-\d.\d', expression_object_type):
            Overview_Content += '<br/><table><tr><th>Generated Expression Object</th>'
            Overview_Content += '<th></th></tr>'
            Overview_Content += '<tr><th>Expression Name</th><th>Condition</th></tr>'
            Overview_Content += '<tr><td>{} ({})</td>'.format(expression_info[1],
                                                              obj_ref)
            Overview_Content += '<td>{}</td></tr>'.format(expression_data['condition'])
            Overview_Content += '</table>'
        elif re.match('KBaseSets.ExpressionSet-\d.\d', expression_object_type):
            Overview_Content += '<br/><table><tr><th>Generated ExpressionSet Object</th></tr>'
            Overview_Content += '<tr><td>{} ({})'.format(expression_info[1],
                                                         obj_ref)
            Overview_Content += '</td></tr></table>'
            Overview_Content += '<p><br/></p>'
            Overview_Content += '<table><tr><th>Generated Expression Objects</th>'
            Overview_Content += '<th></th></tr>'
            Overview_Content += '<tr><th>Expression Name</th><th>Condition</th></tr>'
            for item in expression_data['items']:
                item_expression_object = self.ws.get_objects2({'objects':
                                                              [{'ref': item['ref']}]})['data'][0]
                item_expression_info = item_expression_object['info']
                item_expression_data = item_expression_object['data']
                expression_name = item_expression_info[1]
                Overview_Content += '<tr><td>{} ({})</td>'.format(expression_name,
                                                                  item['ref'])
                Overview_Content += '<td>{}</td>'.format(item_expression_data['condition'])
                Overview_Content += '</tr>'
            Overview_Content += '</table>'
        with open(result_file_path, 'w') as result_file:
            with open(os.path.join(os.path.dirname(__file__), 'report_template.html'),
                      'r') as report_template_file:
                report_template = report_template_file.read()
                report_template = report_template.replace('<p>Overview_Content</p>',
                                                          Overview_Content)
                result_file.write(report_template)

        html_report.append({'path': result_file_path,
                            'name': os.path.basename(result_file_path),
                            'label': os.path.basename(result_file_path),
                            'description': 'HTML summary report for StringTie App'})
        return html_report

    def _generate_merge_report(self, workspace_name, result_directory):
        """
        _generate_merge_report: generate summary report
        """

        log('creating merge report')

        output_files = self._generate_output_file_list(result_directory)
        output_html_files = self._generate_merge_html_report(result_directory)

        report_params = {'message': '',
                         'workspace_name': workspace_name,
                         'file_links': output_files,
                         'html_links': output_html_files,
                         'direct_html_link_index': 0,
                         'html_window_height': 366,
                         'report_object_name': 'kb_stringtie_report_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def _generate_report(self, obj_ref, workspace_name, result_directory,
                         exprMatrix_FPKM_ref=None, exprMatrix_TPM_ref=None, genome_ref=None):
        """
        _generate_report: generate summary report
        """

        log('creating report')

        output_files = self._generate_output_file_list(result_directory)
        output_html_files = self._generate_html_report(result_directory,
                                                       obj_ref)

        expression_object = self.ws.get_objects2({'objects':
                                                 [{'ref': obj_ref}]})['data'][0]
        expression_info = expression_object['info']
        expression_data = expression_object['data']
        objects_created = []

        expression_object_type = expression_info[2]
        if re.match('KBaseRNASeq.RNASeqExpression-\d+.\d+', expression_object_type):
            objects_created.append({'ref': obj_ref,
                                    'description': 'Expression generated by StringTie'})
        elif re.match('KBaseSets.ExpressionSet-\d+.\d+', expression_object_type):
            objects_created.append({'ref': obj_ref,
                                    'description': 'ExpressionSet generated by StringTie'})
            items = expression_data['items']
            for item in items:
                objects_created.append({'ref': item['ref'],
                                        'description': 'Expression generated by StringTie'})
            objects_created.append({'ref': exprMatrix_FPKM_ref,
                                    'description': 'FPKM ExpressionMatrix generated by StringTie'})
            objects_created.append({'ref': exprMatrix_TPM_ref,
                                    'description': 'TPM ExpressionMatrix generated by StringTie'})
        if genome_ref:
            objects_created.append({'ref': genome_ref,
                                    'description': 'Genome containing novel transcripts generated '
                                                   'by StringTie'})

        report_params = {'message': '',
                         'workspace_name': workspace_name,
                         'file_links': output_files,
                         'objects_created': objects_created,
                         'html_links': output_html_files,
                         'direct_html_link_index': 0,
                         'html_window_height': 366,
                         'report_object_name': 'kb_stringtie_report_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def _process_alignment_object(self, params):
        """
        _process_alignment_object: process KBaseRNASeq.RNASeqAlignment type input object
        """

        try:

            log('start processing RNASeqAlignment object\n')
            log('params:\n{}'.format(json.dumps(params, indent=1)))
            alignment_ref = params.get('alignment_ref')

            alignment_set_object = self.ws.get_objects2({'objects':
                                                        [{'ref': alignment_ref}]})['data'][0]

            alignment_info = alignment_set_object['info']
            alignment_data = alignment_set_object['data']

            alignment_name = alignment_info[1]
            alignment_label = alignment_data['condition']

            result_directory = os.path.join(self.scratch,
                                            alignment_name + '_' + str(uuid.uuid4()))
            self._mkdir_p(result_directory)

            # input files
            if not params.get('gtf_file'):
                params['gtf_file'] = self._get_gtf_file(alignment_ref, result_directory)
                if params.get('label'):
                    if params['label'] in open(params['gtf_file']).read():
                        raise ValueError("Provided prefix for transcripts matches an existing "
                                         "feature ID. Please select a different label for "
                                         "transcripts.")
            else:
                shutil.copy(params.get('gtf_file'), result_directory)
            params['input_file'] = self._get_input_file(alignment_ref)
            log('using {} as reference annotation file.'.format(params.get('gtf_file')))

            # output files
            self.output_transcripts = 'transcripts.gtf'
            params['output_transcripts'] = os.path.join(result_directory, self.output_transcripts)

            self.gene_abundances_file = 'genes.fpkm_tracking'
            params['gene_abundances_file'] = os.path.join(result_directory,
                                                          self.gene_abundances_file)

            command = self._generate_command(params)
            self._run_command(command)

            if params.get('exchange_gene_ids'):
                exchange_gene_ids(result_directory)

            if ('generate_ws_object' in params and not params.get('generate_ws_object')):
                log('skip generating expression object')
                expression_obj_ref = ''
            else:
                expression_obj_ref = self._save_expression(
                    result_directory,
                    alignment_ref,
                    params.get('workspace_name'),
                    params['expression_suffix'],
                    params.get('genome_ref'),
                    params.get('novel_isoforms', 0))

            returnVal = {'result_directory': result_directory,
                         'expression_obj_ref': expression_obj_ref,
                         'alignment_ref': alignment_ref,
                         'annotation_file': params['gtf_file'],
                         'alignment_label': alignment_label}
        except:
            log('caught exception in worker')
            exctype, value = sys.exc_info()[:2]

            returnVal = {'exception': '{}: {}'.format(exctype, value)}
        finally:
            return returnVal

    def _process_alignment_set_object(self, params):
        """
        _process_alignment_set_object: process KBaseRNASeq.RNASeqAlignmentSet type input object
        """

        log('start processing AlignmentSet object\nparams:\n{}'.format(json.dumps(params,
                                                                                  indent=1)))
        alignment_set_ref = params.get('alignment_set_ref')

        alignment_set = self.set_client.get_reads_alignment_set_v1({
                                                                    'ref': alignment_set_ref,
                                                                    'include_item_info': 0,
                                                                    'include_set_item_ref_paths': 1
                                                                    })
        # pull down the genome once so as to avoid duplicate effort
        if not params.get('gtf_file'):
            alignment_ref = alignment_set["data"]["items"][0]["ref_path"]
            params['gtf_file'] = self._get_gtf_file(alignment_ref, self.scratch)
            if params.get('label'):
                if params['label'] in open(params['gtf_file']).read():
                    raise ValueError("Provided prefix for transcripts matches an existing "
                                     "feature ID. Please select a different label for "
                                     "transcripts.")

        mul_processor_params = []
        for alignment in alignment_set["data"]["items"]:
            alignment_ref = alignment['ref_path']
            alignment_upload_params = params.copy()
            alignment_upload_params['alignment_ref'] = alignment_ref
            mul_processor_params.append(alignment_upload_params)

        cpus = min(params.get('num_threads'), multiprocessing.cpu_count())
        pool = Pool(ncpus=cpus)
        log('running _process_alignment_object with {} cpus'.format(cpus))
        alignment_expression_map = pool.map(self._process_alignment_object,
                                            mul_processor_params)

        for proc_alignment_return in alignment_expression_map:
            if 'exception' in proc_alignment_return:
                error_msg = 'Caught exception in worker\n'
                error_msg += 'Exception: {}'.format(proc_alignment_return['exception'])
                raise ValueError(error_msg)

        result_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(result_directory)

        for proc_alignment_return in alignment_expression_map:
            alignment_ref = proc_alignment_return.get('alignment_ref')
            alignment_info = self.ws.get_object_info3({'objects': [{"ref": alignment_ref}]})
            alignment_name = alignment_info['infos'][0][1]
            self._run_command('cp -R {} {}'.format(proc_alignment_return.get('result_directory'),
                                                   os.path.join(result_directory,
                                                                alignment_name)))
        if ('generate_ws_object' in params and not params.get('generate_ws_object')):
            log('skip generating expression set object')
            expression_obj_ref = ''
            expression_matrix_refs = {}
        else:
            expression_obj_ref = self._save_expression_set(alignment_expression_map,
                                                           alignment_set_ref,
                                                           params.get('workspace_name'),
                                                           params['expression_set_suffix'],
                                                           params.get('genome_ref'))
            expression_matrix_refs = self._save_expression_matrix(expression_obj_ref,
                                                                  params.get('workspace_name'))

        annotation_file_name = os.path.basename(alignment_expression_map[0]['annotation_file'])
        annotation_file_path = os.path.join(result_directory,
                                            os.listdir(result_directory)[0],
                                            annotation_file_name)

        returnVal = {'result_directory': result_directory,
                     'expression_obj_ref': expression_obj_ref,
                     'annotation_file': annotation_file_path,
                     'exprMatrix_FPKM_ref': expression_matrix_refs.get('exprMatrix_FPKM_ref'),
                     'exprMatrix_TPM_ref': expression_matrix_refs.get('exprMatrix_TPM_ref')}

        return returnVal

    def _run_merge_option(self, result_directory, params, annotation_file):

        log('start running stringtie merge')

        result_dirs = os.listdir(result_directory)

        merge_directory = os.path.join(result_directory, 'merge_result')
        self._mkdir_p(merge_directory)

        option_params = params.copy()

        option_params.pop('num_threads', None)
        option_params.pop('ballgown_mode', None)
        option_params.pop('skip_reads_with_no_ref', None)
        option_params.pop('junction_coverage', None)
        option_params.pop('junction_base', None)
        option_params.pop('min_read_coverage', None)
        option_params.pop('min_locus_gap_sep_value', None)

        output_merge = 'stringtie_merge.gtf'
        option_params['output_transcripts'] = os.path.join(merge_directory, output_merge)

        command = self.STRINGTIE_TOOLKIT_PATH + '/stringtie '
        command += '--merge '
        command += '-G {} '.format(annotation_file)

        for key, option in self.OPTIONS_MAP.items():
            option_value = option_params.get(key)
            if key in self.BOOLEAN_OPTIONS and option_value:
                option_value = ' '
            if option_value:
                command += '{} {} '.format(option, option_value)

        for result_dir in result_dirs:
            gtf_file = os.path.join(result_directory, result_dir, 'transcripts.gtf')
            command += '{} '.format(gtf_file)

        self._run_command(command)

    def _get_genome_ref(self, alignment_set_ref):
        """Get a genome ref from an alignment set"""
        alignment_set_data = self.dfu.get_objects(
            {"object_refs": [alignment_set_ref]})['data'][0]['data']

        for alignment in alignment_set_data['items']:
            alignment_data = self.dfu.get_objects(
                {"object_refs":[alignment['ref']]})['data'][0]['data']
            return alignment_data['genome_id']

    def _save_genome_with_novel_isoforms(self, workspace, genome_ref,
                                         gff_file, new_genome_name=None):
        """"""
        log('Saving genome with novel isoforms')
        genome_data = self.dfu.get_objects(
            {"object_refs": [genome_ref]})['data'][0]['data']
        if 'assembly_ref' in genome_data:
            assembly_ref = genome_data['assembly_ref']
        elif 'contigset_ref' in genome_data:
            assembly_ref = genome_data['contigset_ref']
        else:
            raise ValueError("Genome missing assembly")
        fasta_file = self.au.get_assembly_as_fasta(
            {'ref': assembly_ref})['path']
        if not new_genome_name:
            new_genome_name = genome_data['id'] + "_stringtie"
        ret = self.gfu.fasta_gff_to_genome({
                'workspace_name': workspace,
                'genome_name': new_genome_name,
                'fasta_file': {'path': fasta_file},
                'gff_file': {'path': gff_file},
                'source': 'StringTie'
            })
        return ret['genome_ref']

    def _novel_isoform_mode(self, alignment_object_ref, params):
        """This is a three step process: First, run StringTie on all the alignments individually
          which will produce novel transcripts. Next, merge the resulting transcripts together.
          Finally, rerun StringTie with the merged GTF file as the reference genome.
          """
        log('running Stringtie the 1st time')
        params.update({'ballgown_mode': 0,
                       'skip_reads_with_no_ref': 0,
                       'generate_ws_object': False,
                       'exchange_gene_ids': 1})
        returnVal = self._process_alignment_set_object(params)
        first_run_result_dir = returnVal.get('result_directory')
        annotation_file = returnVal['annotation_file']

        log('running StringTie merge')
        self._run_merge_option(first_run_result_dir, params, annotation_file)
        merge_file = os.path.join(first_run_result_dir,
                                  'merge_result',
                                  'stringtie_merge.gtf')

        old_genome_ref = self._get_genome_ref(alignment_object_ref)
        ret = self.gfu.genome_to_gff({'genome_ref': old_genome_ref,
                                      'target_dir': first_run_result_dir})
        self._run_gffcompare(ret['file_path'], merge_file)
        comp_file = os.path.join(first_run_result_dir, 'merge_result', 'gffcmp.annotated.gtf')
        upload_file = _make_gff(comp_file, ret['file_path'], params.get('label', 'MSTRG.'))
        params['genome_ref'] = self._save_genome_with_novel_isoforms(
            params['workspace_name'], old_genome_ref, upload_file,
            params.get('novel_isoforms', {}).get('stringtie_genome_name'))
        _update_merge_file(merge_file)

        log('running StringTie the 3rd time with merged gtf')
        params.update({'gtf_file': merge_file,
                       'generate_ws_object': True,
                       'exchange_gene_ids': 0,
                       'ballgown_mode': 1,
                       'skip_reads_with_no_ref': 1})
        returnVal = self._process_alignment_set_object(params)

        shutil.move(os.path.join(first_run_result_dir, 'merge_result'),
                    returnVal.get('result_directory'))

        report_output = self._generate_report(returnVal.get('expression_obj_ref'),
                                              params.get('workspace_name'),
                                              returnVal.get('result_directory'),
                                              returnVal.get('exprMatrix_FPKM_ref'),
                                              returnVal.get('exprMatrix_TPM_ref'),
                                              params['genome_ref'])
        return report_output, returnVal

    def __init__(self, config, parent_ref=None):
        self.ws_url = config["workspace-url"]
        self.callback_url = config['SDK_CALLBACK_URL']
        self.token = config['KB_AUTH_TOKEN']
        self.shock_url = config['shock-url']
        self.srv_wiz_url = config['srv-wiz-url']
        self.scratch = config['scratch']
        self.dfu = DataFileUtil(self.callback_url)
        self.gfu = GenomeFileUtil(self.callback_url)
        self.rau = ReadsAlignmentUtils(self.callback_url)
        self.au = AssemblyUtil(self.callback_url)
        self.eu = ExpressionUtils(self.callback_url)
        self.ws = Workspace(self.ws_url, token=self.token)
        self.set_client = SetAPI(self.srv_wiz_url, service_ver='dev')
        self.parent_ref = parent_ref

    def get_ref_with_parent(self, ref):
        """ Use a reference chain if the parent is defined """
        if self.parent_ref and self.parent_ref != ref:
            return self.parent_ref + ";" + ref
        return ref

    def run_stringtie_app(self, params):
        """
        run_stringtie_app: run StringTie app
        (http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)

        required params:
        alignment_object_ref: Alignment or AlignmentSet object reference
        workspace_name: the name of the workspace it gets saved to
        expression_set_suffix: suffix append to expression set object name
        expression_suffix: suffix append to expression object name
        mode: one of ['normal', 'merge', 'novel_isoform']

        optional params:
        num_threads: number of processing threads
        junction_base: junctions that don't have spliced reads
        junction_coverage: junction coverage
        disable_trimming: disables trimming at the ends of the assembled transcripts
        min_locus_gap_sep_value: minimum locus gap separation value
        ballgown_mode: enables the output of Ballgown input table files
        skip_reads_with_no_ref: reads with no reference will be skipped
        maximum_fraction: maximum fraction of muliple-location-mapped reads
        label: prefix for the name of the output transcripts
        min_length: minimum length allowed for the predicted transcripts
        min_read_coverage: minimum input transcript coverage
        min_isoform_abundance: minimum isoform abundance

        return:
        result_directory: folder path that holds all files generated by run_stringtie_app
        expression_obj_ref: generated Expression/ExpressionSet object reference
        report_name: report name generated by KBaseReport
        report_ref: report reference generated by KBaseReport
        """
        log('--->\nrunning StringTieUtil.run_stringtie\n' +
            'params:\n{}'.format(json.dumps(params, indent=1)))

        self._validate_run_stringtie_params(params)
        if isinstance(params.get('novel_isoforms'), dict) and \
                "transcript_label" in params['novel_isoforms']:
            params['label'] = params['novel_isoforms']['transcript_label']

        alignment_object_ref = self.get_ref_with_parent(ref=params.get('alignment_object_ref'))
        print("About to get info for ", alignment_object_ref)
        alignment_object_info = self.ws.get_object_info3({"objects":
                                                         [{"ref": alignment_object_ref}]}
                                                         )['infos'][0]
        alignment_object_type = alignment_object_info[2]

        if re.match('KBaseRNASeq.RNASeqAlignment-\d.\d', alignment_object_type):
            params.update({'alignment_ref': alignment_object_ref})
            returnVal = self._process_alignment_object(params)
            report_output = self._generate_report(returnVal.get('expression_obj_ref'),
                                                  params.get('workspace_name'),
                                                  returnVal.get('result_directory'))
            returnVal.update(report_output)
        elif (re.match('KBaseRNASeq.RNASeqAlignmentSet-\d.\d', alignment_object_type) or
              re.match('KBaseSets.ReadsAlignmentSet-\d.\d', alignment_object_type)):
            params.update({'alignment_set_ref': alignment_object_ref})
            if params.get('novel_isoforms'):
                report_output, returnVal = self._novel_isoform_mode(alignment_object_ref, params)
            else:
                params.update({'ballgown_mode': 1,
                               'skip_reads_with_no_ref': 1,
                               'exchange_gene_ids': 0})

                returnVal = self._process_alignment_set_object(params)

                report_output = self._generate_report(returnVal.get('expression_obj_ref'),
                                                      params.get('workspace_name'),
                                                      returnVal.get('result_directory'),
                                                      returnVal.get('exprMatrix_FPKM_ref'),
                                                      returnVal.get('exprMatrix_TPM_ref'))
            returnVal.update(report_output)
        else:
            error_msg = 'Invalid input object type\nObject info:\n{}'.format(alignment_object_info)
            raise ValueError(error_msg)

        return returnVal
