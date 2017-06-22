import time
import json
import os
import uuid
import errno
import subprocess
import math
import re
from pathos.multiprocessing import ProcessingPool as Pool
import multiprocessing
import zipfile

from DataFileUtil.DataFileUtilClient import DataFileUtil
from Workspace.WorkspaceClient import Workspace as Workspace
from KBaseReport.KBaseReportClient import KBaseReport
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
from ReadsAlignmentUtils.ReadsAlignmentUtilsClient import ReadsAlignmentUtils


def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class StringTieUtil:
    STRINGTIE_TOOLKIT_PATH = '/kb/deployment/bin/StringTie'
    PREPDE_TOOLKIT_PATH = '/kb/deployment/bin/prepDE'
    GFFREAD_TOOLKIT_PATH = '/kb/deployment/bin/gffread'

    OPTIONS_MAP = {
                    'output_transcripts': '-o',
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
                    'min_isoform_abundance': '-f'
                   }

    BOOLEAN_OPTIONS = ['disable_trimming', 'ballgown_mode', 'skip_reads_with_no_ref']

    def _validate_run_stringtie_params(self, params):
        """
        _validate_run_stringtie_params:
                validates params passed to run_stringtie method
        """

        log('Start validating run_stringtie params')

        # check for required parameters
        for p in ['alignment_object_ref', 'workspace_name']:
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

        log('Generated stringtie command: {}'.format(command))

        return command

    def _run_command(self, command):
        """
        _run_command: run command and print result
        """

        log('Start executing command:\n{}'.format(command))
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

    def _run_prepDE(self, result_directory):
        """
        _run_prepDE: run prepDE.py script

        ref: http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual#deseq
        """

        log('generating matrix of read counts')
        command = self.PREPDE_TOOLKIT_PATH + '/prepDE.py '
        command += '-i {} '.format(os.path.dirname(result_directory))
        command += '-g {} '.format(os.path.join(result_directory, 'gene_count_matrix.csv'))
        command += '-t {} '.format(os.path.join(result_directory, 'transcript_count_matrix.csv'))

        self._run_command(command)

    def _get_input_file(self, alignment_ref):
        """
        _get_input_file: get input  SAM/BAM file from Alignment object
        """

        bam_file_archive = self.rau.download_alignment({'source_ref': alignment_ref})['bam_file']
        self.dfu.unpack_file({'file_path': bam_file_archive})

        bam_file = os.path.join(os.path.dirname(bam_file_archive), 'accepted_hits.bam')

        return bam_file

    def _get_gtf_file(self, alignment_ref):
        """
        _get_gtf_file: get the reference annotation file (in GTF or GFF3 format)
        """
        result_directory = self.scratch
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
            annotation_file = self._create_gtf_file(genome_ref)

        return annotation_file

    def _save_gff_annotation(self, genome_id, gtf_file, workspace_name):
        """
        _save_gff_annotation: save GFFAnnotation object to workspace
        """
        log('start saving GffAnnotation object')

        if isinstance(workspace_name, int) or workspace_name.isdigit():
            workspace_id = workspace_name
        else:
            workspace_id = self.dfu.ws_name_to_id(workspace_name)

        genome_data = self.ws.get_objects2({'objects':
                                            [{'ref': genome_id}]})['data'][0]['data']
        genome_name = genome_data.get('id')
        genome_scientific_name = genome_data.get('scientific_name')
        gff_annotation_name = genome_name+"_GTF_Annotation"
        file_to_shock_result = self.dfu.file_to_shock({'file_path': gtf_file,
                                                       'make_handle': True})
        gff_annotation_data = {'handle': file_to_shock_result['handle'],
                               'size': file_to_shock_result['size'],
                               'genome_id': genome_id,
                               'genome_scientific_name': genome_scientific_name}

        object_type = 'KBaseRNASeq.GFFAnnotation'

        save_object_params = {
            'id': workspace_id,
            'objects': [{
                            'type': object_type,
                            'data': gff_annotation_data,
                            'name': gff_annotation_name
                        }]
        }

        dfu_oi = self.dfu.save_objects(save_object_params)[0]
        gff_annotation_obj_ref = str(dfu_oi[6]) + '/' + str(dfu_oi[0]) + '/' + str(dfu_oi[4])

        return gff_annotation_obj_ref

    def _create_gtf_file(self, genome_ref):
        """
        _create_gtf_file: create reference annotation file from genome
        """
        log('start generating reference annotation file')
        result_directory = self.scratch

        genome_gff_file = self.gfu.genome_to_gff({'genome_ref': genome_ref,
                                                  'target_dir': result_directory})['file_path']

        gtf_ext = '.gtf'
        if not genome_gff_file.endswith(gtf_ext):
            gtf_path = os.path.splitext(genome_gff_file)[0] + '.gtf'
            self._run_gffread(genome_gff_file, gtf_path)
        else:
            gtf_path = genome_gff_file

        return gtf_path

    def _generate_expression_set_data(self, alignment_expression_map, alignment_set_ref,
                                      workspace_name):
        """
        _generate_expression_set_data: generate ExpressionSet object with stringtie output files
        """
        alignment_set_data_object = self.ws.get_objects2({'objects':
                                                         [{'ref': alignment_set_ref}]})['data'][0]

        alignment_set_name = alignment_set_data_object['info'][1]
        expression_set_name = re.sub('_[Aa]lignment_*[Ss]et',
                                     '_stringtie_expression_set',
                                     alignment_set_name)

        alignment_set_data = alignment_set_data_object['data']

        expression_set_data = {
            'tool_used': 'StringTie',
            'tool_version': '1.3.3b',
            'id': expression_set_name,
            'alignmentSet_id': alignment_set_ref,
            'genome_id': alignment_set_data.get('genome_id'),
            'sampleset_id': alignment_set_data.get('sampleset_id')
        }

        sample_expression_ids = []
        mapped_expression_objects = []
        mapped_expression_ids = []

        for alignment_expression in alignment_expression_map:
            alignment_ref = alignment_expression.get('alignment_ref')
            expression_ref = alignment_expression.get('expression_obj_ref')
            sample_expression_ids.append(expression_ref)
            mapped_expression_ids.append({alignment_ref: expression_ref})
            alignment_name = self.ws.get_object_info([{"ref": alignment_ref}],
                                                     includeMetadata=None)[0][1]
            expression_name = self.ws.get_object_info([{"ref": expression_ref}],
                                                      includeMetadata=None)[0][1]
            mapped_expression_objects.append({alignment_name: expression_name})

        expression_set_data['sample_expression_ids'] = sample_expression_ids
        expression_set_data['mapped_expression_objects'] = mapped_expression_objects
        expression_set_data['mapped_expression_ids'] = mapped_expression_ids

        return expression_set_data

    def _generate_expression_data(self, result_directory, alignment_ref,
                                  gtf_file, workspace_name):
        """
        _generate_expression_data: generate Expression object with stringtie output files
        """
        alignment_data_object = self.ws.get_objects2({'objects':
                                                     [{'ref': alignment_ref}]})['data'][0]

        alignment_name = alignment_data_object['info'][1]
        expression_name = re.sub('_[Aa]lignment',
                                 '_stringtie_expression',
                                 alignment_name)

        expression_data = {
            'id': expression_name,
            'type': 'RNA-Seq',
            'numerical_interpretation': 'FPKM',
            'processing_comments': 'log2 Normalized',
            'tool_used': 'StringTie',
            'tool_version': '1.3.3b'
        }
        alignment_data = alignment_data_object['data']

        condition = alignment_data.get('condition')
        expression_data.update({'condition': condition})

        genome_id = alignment_data.get('genome_id')
        expression_data.update({'genome_id': genome_id})

        gff_annotation_obj_ref = self._save_gff_annotation(genome_id, gtf_file, workspace_name)
        expression_data.update({'annotation_id': gff_annotation_obj_ref})

        read_sample_id = alignment_data.get('read_sample_id')
        expression_data.update({'mapped_rnaseq_alignment': {read_sample_id: alignment_ref}})

        exp_dict = self._parse_FPKMtracking(os.path.join(result_directory,
                                                         self.gene_abundances_file), 'FPKM')
        expression_data.update({'expression_levels': exp_dict})

        tpm_exp_dict = self._parse_FPKMtracking(os.path.join(result_directory,
                                                             self.gene_abundances_file), 'TPM')
        expression_data.update({'tpm_expression_levels': tpm_exp_dict})

        handle = self.dfu.file_to_shock({'file_path': result_directory,
                                         'pack': 'zip',
                                         'make_handle': True})['handle']
        expression_data.update({'file': handle})

        return expression_data

    def _parse_FPKMtracking(self, filename, metric):
        result = {}
        pos1 = 0
        if metric == 'FPKM':
            pos2 = 7
        if metric == 'TPM':
            pos2 = 8

        with open(filename) as f:
            next(f)
            for line in f:
                larr = line.split("\t")
                if larr[pos1] != "":
                    result[larr[pos1]] = math.log(float(larr[pos2]) + 1, 2)
        return result

    def _save_expression(self, result_directory, alignment_ref, workspace_name, gtf_file):
        """
        _save_expression: save Expression object to workspace
        """
        log('start saving Expression object')
        if isinstance(workspace_name, int) or workspace_name.isdigit():
            workspace_id = workspace_name
        else:
            workspace_id = self.dfu.ws_name_to_id(workspace_name)

        expression_data = self._generate_expression_data(result_directory,
                                                         alignment_ref,
                                                         gtf_file,
                                                         workspace_name)

        object_type = 'KBaseRNASeq.RNASeqExpression'
        save_object_params = {
            'id': workspace_id,
            'objects': [{
                            'type': object_type,
                            'data': expression_data,
                            'name': expression_data.get('id')
                        }]
        }

        dfu_oi = self.dfu.save_objects(save_object_params)[0]
        expression_ref = str(dfu_oi[6]) + '/' + str(dfu_oi[0]) + '/' + str(dfu_oi[4])

        return expression_ref

    def _save_expression_set(self, alignment_expression_map, alignment_set_ref, workspace_name):
        """
        _save_expression_set: save ExpressionSet object to workspace
        """
        log('start saving ExpressionSet object')
        if isinstance(workspace_name, int) or workspace_name.isdigit():
            workspace_id = workspace_name
        else:
            workspace_id = self.dfu.ws_name_to_id(workspace_name)

        expression_set_data = self._generate_expression_set_data(alignment_expression_map,
                                                                 alignment_set_ref,
                                                                 workspace_name)

        object_type = 'KBaseRNASeq.RNASeqExpressionSet'
        save_object_params = {
            'id': workspace_id,
            'objects': [{
                            'type': object_type,
                            'data': expression_set_data,
                            'name': expression_set_data.get('id')
                        }]
        }

        dfu_oi = self.dfu.save_objects(save_object_params)[0]
        expression_set_ref = str(dfu_oi[6]) + '/' + str(dfu_oi[0]) + '/' + str(dfu_oi[4])

        return expression_set_ref

    def _generate_output_file_list(self, result_directory):
        """
        _generate_output_file_list: zip result files and generate file_links for report
        """
        log('Start packing result files')
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
            Overview_Content += '<p>Generated Expression Object:</p><p>{}</p>'.format(
                                                                expression_object.get('info')[1])
        elif re.match('KBaseRNASeq.RNASeqExpressionSet-\d.\d', expression_object_type):
            Overview_Content += '<p>Generated Expression Set Object:</p><p>{}</p>'.format(
                                                                expression_object.get('info')[1])
            Overview_Content += '<br><p>Generated Expression Object:</p>'
            for expression_ref in expression_object['data']['sample_expression_ids']:
                expression_name = self.ws.get_object_info([{"ref": expression_ref}],
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
                            'description': 'HTML summary report for MaxBin2 App'})
        return html_report

    def _generate_report(self, obj_ref, workspace_name, result_directory):
        """
        _generate_report: generate summary report
        """
        log('creating report')

        output_files = self._generate_output_file_list(result_directory)
        output_html_files = self._generate_html_report(result_directory,
                                                       obj_ref)

        report_params = {
              'message': '',
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

    def _process_alignment_object(self, params):
        """
        _process_alignment_object: process KBaseRNASeq.RNASeqAlignment type input object
        """
        log('start processing RNASeqAlignment object\nparams:\n{}'.format(
                                                                json.dumps(params, indent=1)))
        alignment_ref = params.get('alignment_ref')

        result_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(result_directory)

        # input files
        params['input_file'] = self._get_input_file(alignment_ref)
        if not params.get('gtf_file'):
            params['gtf_file'] = self._get_gtf_file(alignment_ref)

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
                                                   params['gtf_file'])

        returnVal = {'result_directory': result_directory,
                     'expression_obj_ref': expression_obj_ref,
                     'alignment_ref': alignment_ref}

        return returnVal

    def _process_alignment_set_object(self, params):
        """
        _process_alignment_set_object: process KBaseRNASeq.RNASeqAlignmentSet type input object
        """
        log('start processing RNASeqAlignmentSet object\nparams:\n{}'.format(
                                                                json.dumps(params, indent=1)))

        alignment_set_ref = params.get('alignment_set_ref')

        params['gtf_file'] = self._get_gtf_file(alignment_set_ref)

        alignment_set_data = self.ws.get_objects2({'objects':
                                                   [{'ref': alignment_set_ref}]})['data'][0]['data']

        mapped_alignment_ids = alignment_set_data['mapped_alignments_ids']
        mul_processor_params = []
        for i in mapped_alignment_ids:
            for sample_name, alignment_id in i.items():
                aliment_upload_params = params.copy()
                aliment_upload_params['alignment_ref'] = alignment_id
                mul_processor_params.append(aliment_upload_params)

        cpus = min(params.get('num_threads'), multiprocessing.cpu_count())
        pool = Pool(ncpus=cpus)
        log('running _process_alignment_object with {} cpus'.format(cpus))
        alignment_expression_map = pool.map(self._process_alignment_object, mul_processor_params)

        result_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(result_directory)

        for proc_alignment_return in alignment_expression_map:
            expression_obj_ref = proc_alignment_return.get('expression_obj_ref')
            expression_name = self.ws.get_object_info([{"ref": expression_obj_ref}],
                                                      includeMetadata=None)[0][1]
            self._run_command('cp -R {} {}'.format(proc_alignment_return.get('result_directory'),
                                                   os.path.join(result_directory, expression_name)))

        expression_obj_ref = self._save_expression_set(alignment_expression_map,
                                                       alignment_set_ref,
                                                       params.get('workspace_name'))

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
        self.dfu = DataFileUtil(self.callback_url)
        self.gfu = GenomeFileUtil(self.callback_url)
        self.rau = ReadsAlignmentUtils(self.callback_url)
        self.ws = Workspace(self.ws_url, token=self.token)

        self.scratch = os.path.join(config['scratch'], str(uuid.uuid4()))
        self._mkdir_p(self.scratch)

    def run_stringtie_app(self, params):
        """
        run_stringtie_app: run StringTie app
        (http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)

        required params:
        alignment_object_ref: Alignment or AlignmentSet object reference
        workspace_name: the name of the workspace it gets saved to.

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
        alignment_object_info = self.ws.get_object_info3({
                                        "objects": [{"ref": alignment_object_ref}]})['infos'][0]
        alignment_object_type = alignment_object_info[2]

        if re.match('KBaseRNASeq.RNASeqAlignment-\d.\d', alignment_object_type):
            params.update({'alignment_ref': alignment_object_ref})
            returnVal = self._process_alignment_object(params)
            report_output = self._generate_report(returnVal.get('expression_obj_ref'),
                                                  params.get('workspace_name'),
                                                  returnVal.get('result_directory'))
            returnVal.update(report_output)
        elif re.match('KBaseRNASeq.RNASeqAlignmentSet-\d.\d', alignment_object_type):
            params.update({'alignment_set_ref': alignment_object_ref})
            returnVal = self._process_alignment_set_object(params)
        else:
            raise ValueError('None RNASeqAlignment type\nObject info:\n{}'.format(
                                                                    alignment_object_info))

        return returnVal
