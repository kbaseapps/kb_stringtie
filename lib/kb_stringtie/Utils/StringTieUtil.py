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

from DataFileUtil.DataFileUtilClient import DataFileUtil
from Workspace.WorkspaceClient import Workspace as Workspace
from KBaseReport.KBaseReportClient import KBaseReport
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
from ReadsAlignmentUtils.ReadsAlignmentUtilsClient import ReadsAlignmentUtils
from ExpressionUtils. ExpressionUtilsClient import ExpressionUtils
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

        if params.get('merge'):
            command += '--merge '

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
        # genome_name = self.ws.get_object_info([{"ref": genome_ref}], includeMetadata=None)[0][1]
        # ws_gtf = genome_name+"_GTF_Annotation"

        genome_data = self.ws.get_objects2({'objects':
                                            [{'ref': genome_ref}]})['data'][0]['data']

        gff_handle_ref = genome_data.get('gff_handle_ref')

        if gff_handle_ref:
            log('getting reference annotation file from genome')
            annotation_file = self.dfu.shock_to_file({'handle_id': gff_handle_ref,
                                                      'file_path': result_directory,
                                                      'unpack': 'unpack'})['file_path']
        else:
            annotation_file = self._create_gtf_file(genome_ref, result_directory)

        return annotation_file

    def _create_gtf_file(self, genome_ref, result_directory):
        """
        _create_gtf_file: create reference annotation file from genome
        """

        log('start generating reference annotation file')

        genome_gff_file = self.gfu.genome_to_gff({'genome_ref': genome_ref,
                                                  'target_dir': result_directory})['file_path']

        gtf_ext = '.gtf'
        if not genome_gff_file.endswith(gtf_ext):
            gtf_path = os.path.splitext(genome_gff_file)[0] + '.gtf'
            self._run_gffread(genome_gff_file, gtf_path)
        else:
            gtf_path = genome_gff_file

        return gtf_path

    def _save_expression(self, result_directory, alignment_ref, workspace_name, gtf_file, 
                         expression_suffix):
        """
        _save_expression: save Expression object to workspace
        """

        log('start saving Expression object')

        alignment_data_object = self.ws.get_objects2({'objects':
                                                     [{'ref': alignment_ref}]})['data'][0]

        alignment_name = alignment_data_object['info'][1]
        expression_obj_name = re.sub('_[Aa]lignment', expression_suffix, alignment_name)
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
            items.append({'ref': alignment_expression.get('expression_obj_ref')})

        expression_set_data = {'description': 'ExpressionSet using DESeq2', 
                               'items': items}

        alignment_set_data_object = self.ws.get_objects2({'objects':
                                                         [{'ref': alignment_set_ref}]})['data'][0]

        alignment_set_name = alignment_set_data_object['info'][1]
        expression_set_name = re.sub('_[Aa]lignment_*[Ss]et',
                                     expression_set_suffix,
                                     alignment_set_name)

        expression_set_save_params = {'data': expression_set_data,
                                      'workspace': workspace_name,
                                      'output_object_name': expression_set_name}

        save_result = self.set_client.save_expression_set_v1(expression_set_save_params)
        expression_set_ref = save_result['set_ref']

        return expression_set_ref

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
                    if not (file.endswith('.DS_Store')):
                        zip_file.write(os.path.join(root, file),
                                       os.path.join(os.path.basename(root), file))

        output_files.append({'path': result_file,
                             'name': os.path.basename(result_file),
                             'label': os.path.basename(result_file),
                             'description': 'File(s) generated by StringTie App'})

        return output_files

    def _generate_html_report(self, result_directory, obj_ref):
        """
        _generate_html_report: generate html summary report
        """

        log('Start generating html report')
        html_report = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(output_directory)
        result_file_path = os.path.join(output_directory, 'report.html')

        expression_object = self.ws.get_objects2({'objects':
                                                 [{'ref': obj_ref}]})['data'][0]

        expression_object_type = expression_object.get('info')[2]

        Overview_Content = ''
        if re.match('KBaseRNASeq.RNASeqExpression-\d.\d', expression_object_type):
            Overview_Content += '<p>Generated Expression Object:</p>'
            Overview_Content += '<p>{}</p>'.format(expression_object.get('info')[1])
        elif re.match('KBaseSets.ExpressionSet-\d.\d', expression_object_type):
            Overview_Content += '<p>Generated Expression Set Object:</p>'
            Overview_Content += '<p>{}</p>'.format(expression_object.get('info')[1])
            Overview_Content += '<br><p>Generated Expression Object:</p>'
            for item in expression_object['data']['items']:
                expression_name = self.ws.get_object_info([{"ref": item['ref']}],
                                                          includeMetadata=None)[0][1]
                Overview_Content += '<p>{}</p>'.format(expression_name)

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

    def _generate_report(self, obj_ref, workspace_name, result_directory):
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
        if re.match('KBaseRNASeq.RNASeqExpression-\d.\d', expression_object_type):
            objects_created = [{'ref': obj_ref,
                                'description': 'Expression generated by StringTie'}]
        elif re.match('KBaseSets.ExpressionSet-\d.\d', expression_object_type):
            objects_created = [{'ref': obj_ref,
                                'description': 'ExpressionSet generated by StringTie'}]
            items = expression_data['items']
            for item in items:
                objects_created.append({'ref': item['ref'],
                                        'description': 'Expression generated by StringTie'})

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
        log('start processing RNASeqAlignment object\nparams:\n{}'.format(json.dumps(params, 
                                                                                     indent=1)))
        alignment_ref = params.get('alignment_ref')

        alignment_object_info = self.ws.get_object_info3({"objects": 
                                                         [{"ref": alignment_ref}]}
                                                         )['infos'][0]
        alignment_name = alignment_object_info[1]

        result_directory = os.path.join(self.scratch, 
                                        alignment_name + '_' + str(int(time.time() * 100)))
        self._mkdir_p(result_directory)

        # input files
        params['input_file'] = self._get_input_file(alignment_ref)
        if not params.get('gtf_file'):
            params['gtf_file'] = self._get_gtf_file(alignment_ref, result_directory)
        else:
            shutil.copy(params.get('gtf_file'), result_directory)

        # output files
        self.output_transcripts = 'transcripts.gtf'
        params['output_transcripts'] = os.path.join(result_directory, self.output_transcripts)

        self.gene_abundances_file = 'genes.fpkm_tracking'
        params['gene_abundances_file'] = os.path.join(result_directory, self.gene_abundances_file)

        command = self._generate_command(params)
        self._run_command(command)

        expression_obj_ref = self._save_expression(result_directory,
                                                   alignment_ref,
                                                   params.get('workspace_name'),
                                                   params['gtf_file'],
                                                   params['expression_suffix'])

        returnVal = {'result_directory': result_directory,
                     'expression_obj_ref': expression_obj_ref,
                     'alignment_ref': alignment_ref}

        return returnVal

    def _process_alignment_set_object(self, params):
        """
        _process_alignment_set_object: process KBaseRNASeq.RNASeqAlignmentSet type input object
        """

        log('start processing AlignmentSet object\nparams:\n{}'.format(json.dumps(params, 
                                                                                  indent=1)))

        alignment_set_ref = params.get('alignment_set_ref')
        alignment_set_object = self.ws.get_objects2({'objects':
                                                    [{'ref': alignment_set_ref}]}
                                                    )['data'][0]

        alignment_set_info = alignment_set_object['info']
        alignment_set_data = alignment_set_object['data']

        alignment_set_type = alignment_set_info[2]

        mul_processor_params = []
        if re.match('KBaseRNASeq.RNASeqAlignmentSet-\d.\d', alignment_set_type):
            mapped_alignment_ids = alignment_set_data['mapped_alignments_ids']
            for i in mapped_alignment_ids:
                for sample_name, alignment_id in i.items():
                    aliment_upload_params = params.copy()
                    aliment_upload_params['alignment_ref'] = alignment_id
                    mul_processor_params.append(aliment_upload_params)
        elif re.match('KBaseSets.ReadsAlignmentSet-\d.\d', alignment_set_type):
            items = alignment_set_data['items']
            for item in items:
                alignment_ref = item['ref']
                aliment_upload_params = params.copy()
                aliment_upload_params['alignment_ref'] = alignment_ref
                mul_processor_params.append(aliment_upload_params)

        cpus = min(params.get('num_threads'), multiprocessing.cpu_count())
        pool = Pool(ncpus=cpus)
        log('running _process_alignment_object with {} cpus'.format(cpus))
        alignment_expression_map = pool.map(self._process_alignment_object, mul_processor_params)

        result_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(result_directory)

        for proc_alignment_return in alignment_expression_map:
            alignment_ref = proc_alignment_return.get('alignment_ref')
            alignment_name = self.ws.get_object_info([{"ref": alignment_ref}],
                                                     includeMetadata=None)[0][1]
            self._run_command('cp -R {} {}'.format(proc_alignment_return.get('result_directory'),
                                                   os.path.join(result_directory, 
                                                                alignment_name)))

        expression_obj_ref = self._save_expression_set(alignment_expression_map,
                                                       alignment_set_ref,
                                                       params.get('workspace_name'),
                                                       params['expression_set_suffix'])

        returnVal = {'result_directory': result_directory,
                     'expression_obj_ref': expression_obj_ref}

        report_output = self._generate_report(expression_obj_ref,
                                              params.get('workspace_name'),
                                              result_directory)
        returnVal.update(report_output)

        return returnVal

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
        self.eu = ExpressionUtils(self.callback_url, service_ver='dev')
        self.ws = Workspace(self.ws_url, token=self.token)
        self.set_client = SetAPI(self.srv_wiz_url)
       
    def run_stringtie_app(self, params):
        """
        run_stringtie_app: run StringTie app
        (http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)

        required params:
        alignment_object_ref: Alignment or AlignmentSet object reference
        workspace_name: the name of the workspace it gets saved to
        expression_set_suffix: suffix append to expression set object name
        expression_suffix: suffix append to expression object name

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
        merge: set transcript merge mode

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
            returnVal = self._process_alignment_set_object(params)
        else:
            error_msg = 'Invalid input object type\nObject info:\n{}'.format(alignment_object_info)
            raise ValueError(error_msg)

        return returnVal
