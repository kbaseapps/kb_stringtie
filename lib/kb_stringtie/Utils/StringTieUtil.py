import time
import json
import os
import uuid
import errno
import subprocess

from DataFileUtil.DataFileUtilClient import DataFileUtil
from Workspace.WorkspaceClient import Workspace as Workspace
# from KBaseReport.KBaseReportClient import KBaseReport
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil


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

    def _validate_run_stringtie_params(self, params):
        """
        _validate_run_stringtie_params:
                validates params passed to run_stringtie method
        """

        log('Start validating run_stringtie params')

        # check for required parameters
        for p in ['alignment_ref', 'expression_object_name', 'workspace_name']:
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
            if isinstance(option_value, bool) and option_value:
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

    def _get_input_file(self, alignment_ref):
        """
        _get_input_file: get input  SAM/BAM file from Alignment object
        """
        result_directory = self.scratch

        return 'input_file_path'

    def _get_gtf_file(self, alignment_ref):
        """
        _get_gtf_file: get the reference annotation file (in GTF or GFF3 format)
        """
        result_directory = self.scratch
        alignment_info = self.ws.get_objects2({'objects':
                                               [{'ref': alignment_ref}]})['data'][0]['info']

        genome_ref = alignment_info[-1].get('genome_id')

        # genome_ref = '15206/137/1'

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

    def _create_gtf_file(self, genome_ref):
        """
        _create_gtf_file: create reference annotation file from genome
        """

        result_directory = self.scratch

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

    def __init__(self, config):
        self.ws_url = config["workspace-url"]
        self.callback_url = config['SDK_CALLBACK_URL']
        self.token = config['KB_AUTH_TOKEN']
        self.shock_url = config['shock-url']
        self.dfu = DataFileUtil(self.callback_url)
        self.gfu = GenomeFileUtil(self.callback_url)
        self.au = AssemblyUtil(self.callback_url)
        self.ws = Workspace(self.ws_url, token=self.token)

        self.scratch = os.path.join(config['scratch'], str(uuid.uuid4()))
        self._mkdir_p(self.scratch)

    def run_stringtie_app(self, params):
        """
        run_stringtie_app: run StringTie app
        (http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)

        required params:
        alignment_ref: Alignment object reference
        expression_set_name: ExpressionSet object name and output file header
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

        result_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(result_directory)

        # input files
        alignment_ref = params.get('alignment_ref')
        params['input_file'] = self._get_input_file(alignment_ref)
        params['gtf_file'] = self._get_gtf_file(alignment_ref)

        # output files
        output_transcripts = 'transcripts.gtf'
        params['output_transcripts'] = os.path.join(result_directory, output_transcripts)

        gene_abundances_file = 'genes.fpkm_tracking'
        params['gene_abundances_file'] = os.path.join(result_directory, gene_abundances_file)

        command = self._generate_command(params)

        self._run_command(command)

        if params.get('run_matrix_count'):
            self._run_prepDE(result_directory)

        returnVal = {'result_directory': result_directory}
        return returnVal
