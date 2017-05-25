import time
import json
import os
import uuid
import errno
import subprocess
import shutil
import sys
import zipfile

# from DataFileUtil.DataFileUtilClient import DataFileUtil
# from KBaseReport.KBaseReportClient import KBaseReport
# from AssemblyUtil.AssemblyUtilClient import AssemblyUtil


def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class StringTieUtil:
    STRINGTIE_TOOLKIT_PATH = '/kb/deployment/bin/StringTie'

    MERGE_SPECIFIC_OPTION_MAP = {'min_fpkm': '-F',
                                 'min_tpm': '-T',
                                 'keep_introns': '-i'}

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

        log('Generated run_MaxBin command: {}'.format(command))

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

    def _get_input_file(self):
        return 'input_file_path'

    def _get_gtf_file(self):
        return ''

    def __init__(self, config):
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.shock_url = config['shock-url']
        # self.dfu = DataFileUtil(self.callback_url)
        # self.au = AssemblyUtil(self.callback_url)

    def run_stringtie_app(self, params):
        """
        run_stringtie_app: run StringTie app
        (http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)

        required params:
        assembly_ref: Alignment object reference
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
        input_file = self._get_input_file()
        params['input_file'] = input_file

        gtf_file = self._get_gtf_file()
        params['gtf_file'] = gtf_file

        # output files
        output_transcripts = 'transcripts.gtf'
        params['output_transcripts'] = os.path.join(result_directory, output_transcripts)

        gene_abundances_file = 'genes.fpkm_tracking'
        params['gene_abundances_file'] = os.path.join(result_directory, gene_abundances_file)

        command = self._generate_command(params)

        self._run_command(command)

        returnVal = {'a': 'a'}
        return returnVal
