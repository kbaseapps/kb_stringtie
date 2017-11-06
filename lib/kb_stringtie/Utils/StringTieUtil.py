import time
import json
import os
import uuid
import errno
import subprocess
import re
from pathos.multiprocessing import ProcessingPool as Pool
import multiprocessing
import zipfile
import shutil
import sys
import traceback
import contig_id_mapping as c_mapping
from pprint import pprint

from DataFileUtil.DataFileUtilClient import DataFileUtil
from Workspace.WorkspaceClient import Workspace as Workspace
from KBaseReport.KBaseReportClient import KBaseReport
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
from ReadsAlignmentUtils.ReadsAlignmentUtilsClient import ReadsAlignmentUtils
from ExpressionUtils. ExpressionUtilsClient import ExpressionUtils
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from SetAPI.SetAPIServiceClient import SetAPI


def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class StringTieUtil:
    STRINGTIE_TOOLKIT_PATH = '/kb/deployment/bin/StringTie'
    GFFREAD_TOOLKIT_PATH = '/kb/deployment/bin/gffread'

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
            log('Executed commend:\n{}\n'.format(command) +
                'Exit Code: {}\nOutput:\n{}'.format(exitCode, output))
        else:
            error_msg = 'Error running commend:\n{}\n'.format(command)
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
        _get_gtf_file: get the reference annotation file (in GTF or GFF3 format)
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
        if 'contigset_ref' in ref[0]['data']:
            contig_id = ref[0]['data']['contigset_ref']
        elif 'assembly_ref' in ref[0]['data']:
            contig_id = ref[0]['data']['assembly_ref']
        if contig_id is None:
            raise ValueError(
                "Genome at {0} does not have reference to the assembly object".format(
                    genome_ref))
        print contig_id
        log("Generating GFF file from Genome")
        try:
            ret = self.au.get_assembly_as_fasta({'ref': contig_id})
            fa_output_file = ret['path']

            shutil.copy(fa_output_file, result_directory)
            fa_output_name = os.path.basename(fa_output_file)
            fa_output_file = os.path.join(result_directory, fa_output_name)

            mapping_filename = c_mapping.create_sanitized_contig_ids(fa_output_file)

            # get the GFF
            ret = self.gfu.genome_to_gff({'genome_ref': genome_ref,
                                          'target_dir': result_directory})
            genome_gff_file = ret['file_path']
            c_mapping.replace_gff_contig_ids(genome_gff_file, mapping_filename, to_modified=True)
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

    def _save_expression(self, result_directory, alignment_ref, workspace_name, gtf_file, 
                         expression_suffix):
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
                                    'tool_version': '1.3.3'}

        expression_ref = self.eu.upload_expression(upload_expression_params)['obj_ref']

        return expression_ref

    def _save_expression_set(self, alignment_expression_map, alignment_set_ref, workspace_name,
                             expression_set_suffix):
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
                    if not (file.endswith('.DS_Store') or 
                            os.path.basename(root) == 'merge_result'):
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

            filtered_merge_file = os.path.join(result_directory, 
                                               'merge_result', 
                                               'filtered_stringtie_merge.gtf')
            output_files.append({'path': filtered_merge_file,
                                 'name': os.path.basename(filtered_merge_file),
                                 'label': os.path.basename(filtered_merge_file),
                                 'description': 'filtered merge file generated by StringTie App'})

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
                         exprMatrix_FPKM_ref=None, exprMatrix_TPM_ref=None):
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

        expression_object_type = expression_info[2]
        if re.match('KBaseRNASeq.RNASeqExpression-\d+.\d+', expression_object_type):
            objects_created = [{'ref': obj_ref,
                                'description': 'Expression generated by StringTie'}]
        elif re.match('KBaseSets.ExpressionSet-\d+.\d+', expression_object_type):
            objects_created = [{'ref': obj_ref,
                                'description': 'ExpressionSet generated by StringTie'}]
            items = expression_data['items']
            for item in items:
                objects_created.append({'ref': item['ref'],
                                        'description': 'Expression generated by StringTie'})
            objects_created.append({'ref': exprMatrix_FPKM_ref,
                                    'description': 'FPKM ExpressionMatrix generated by StringTie'})
            objects_created.append({'ref': exprMatrix_TPM_ref,
                                    'description': 'TPM ExpressionMatrix generated by StringTie'})

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
            params['input_file'] = self._get_input_file(alignment_ref)
            if not params.get('gtf_file'):
                params['gtf_file'] = self._get_gtf_file(alignment_ref, result_directory)
            else:
                shutil.copy(params.get('gtf_file'), result_directory)
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
                self._exchange_gene_ids(result_directory)

            if ('generate_ws_object' in params and not params.get('generate_ws_object')):
                log('skip generating expression object')
                expression_obj_ref = ''
            else:
                expression_obj_ref = self._save_expression(result_directory,
                                                           alignment_ref,
                                                           params.get('workspace_name'),
                                                           params['gtf_file'],
                                                           params['expression_suffix'])

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
                                                           params['expression_set_suffix'])
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

    def _filter_merge_file(self, gtf_file):
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

    def _exchange_gene_ids(self, result_directory):
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

    def __init__(self, config):
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

        alignment_object_ref = params.get('alignment_object_ref')
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
            if params.get('mode') in ['merge', 'novel_isoform']:

                log('running Stringtie the 1st time')
                if params.get('mode') == 'novel_isoform':
                    params.update({'ballgown_mode': 0})
                    params.update({'skip_reads_with_no_ref': 0})
                elif params.get('mode') == 'merge':
                    params.update({'ballgown_mode': 1})
                    params.update({'skip_reads_with_no_ref': 1})

                params['generate_ws_object'] = False
                params['exchange_gene_ids'] = False
                returnVal = self._process_alignment_set_object(params)
                first_run_result_dir = returnVal.get('result_directory')
                annotation_file = returnVal['annotation_file']

                log('running StringTie merge')
                self._run_merge_option(first_run_result_dir, params, annotation_file)

                merge_file = os.path.join(first_run_result_dir, 
                                          'merge_result', 
                                          'stringtie_merge.gtf')

                log('running StringTie the 3rd time with merged gtf')
                if params.get('mode') == 'novel_isoform':
                    params.update({'gtf_file': merge_file})
                elif params.get('mode') == 'merge':
                    filtered_merge_file = self._filter_merge_file(merge_file)
                    params.update({'gtf_file': filtered_merge_file})
                    params.update({'generate_ws_object': True})
                    params.update({'exchange_gene_ids': 1})

                params.update({'ballgown_mode': 1})
                params.update({'skip_reads_with_no_ref': 1})
                returnVal = self._process_alignment_set_object(params)

                self._run_command('cp -R {} {}'.format(os.path.join(first_run_result_dir,
                                                       'merge_result'),
                                                       returnVal.get('result_directory')))

                if params.get('mode') == 'novel_isoform':
                    report_output = self._generate_merge_report(params.get('workspace_name'),
                                                                returnVal.get('result_directory'))
                elif params.get('mode') == 'merge':
                    report_output = self._generate_report(returnVal.get('expression_obj_ref'),
                                                          params.get('workspace_name'),
                                                          returnVal.get('result_directory'),
                                                          returnVal.get('exprMatrix_FPKM_ref'),
                                                          returnVal.get('exprMatrix_TPM_ref'))
            else:
                params.update({'ballgown_mode': 1})
                params.update({'skip_reads_with_no_ref': 1})
                params.update({'exchange_gene_ids': 0})

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
