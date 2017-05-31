# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import json  # noqa: F401
import time
import requests  # noqa: F401
import shutil
from mock import patch

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint  # noqa: F401

from biokbase.workspace.client import Workspace as workspaceService
from kb_stringtie.kb_stringtieImpl import kb_stringtie
from kb_stringtie.kb_stringtieServer import MethodContext
from kb_stringtie.authclient import KBaseAuth as _KBaseAuth
from kb_stringtie.Utils.StringTieUtil import StringTieUtil


class kb_stringtieTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_stringtie'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_stringtie',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL)
        cls.serviceImpl = kb_stringtie(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']

        cls.stringtie_runner = StringTieUtil(cls.cfg)
        cls.prepare_data()

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    @classmethod
    def prepare_data(cls):
        input_file = 'samExample.sam'
        input_file_path = os.path.join(cls.scratch, input_file)
        shutil.copy(os.path.join("data", input_file), input_file_path)

        cls.alignment_ref = '21746/5/12'

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_kb_stringtie_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def mock_get_input_file(alignment_ref):
        print 'Mocking StringTieUtil._get_input_file'

        return '/kb/module/work/tmp/samExample.sam'

    def test_bad_run_stringtie_app_params(self):
        invalidate_input_params = {
          'missing_alignment_ref': 'alignment_ref',
          'expression_object_name': 'expression_object_name',
          'workspace_name': 'workspace_name'
        }
        with self.assertRaisesRegexp(
                    ValueError, '"alignment_ref" parameter is required, but missing'):
            self.getImpl().run_stringtie_app(self.getContext(), invalidate_input_params)

        invalidate_input_params = {
          'alignment_ref': 'alignment_ref',
          'missing_expression_object_name': 'expression_object_name',
          'workspace_name': 'workspace_name'
        }
        with self.assertRaisesRegexp(
                    ValueError, '"expression_object_name" parameter is required, but missing'):
            self.getImpl().run_stringtie_app(self.getContext(), invalidate_input_params)

        invalidate_input_params = {
          'alignment_ref': 'alignment_ref',
          'expression_object_name': 'expression_object_name',
          'missing_workspace_name': 'workspace_name'
        }
        with self.assertRaisesRegexp(
                    ValueError, '"workspace_name" parameter is required, but missing'):
            self.getImpl().run_stringtie_app(self.getContext(), invalidate_input_params)

    def test_StringTieUtil_generate_command(self):
        command_params = {
            'num_threads': 4,
            'junction_base': 8,
            'junction_coverage': 0.8,
            'disable_trimming': True,
            'min_locus_gap_sep_value': 60,
            'maximum_fraction': 0.8,
            'label': 'Lable',
            'min_length': 100,
            'min_read_coverage': 1.6,
            'min_isoform_abundance': 0.6,
            'output_transcripts': 'output_transcripts_file',
            'gene_abundances_file': 'gene_abundances_file',
            'input_file': 'input_file',
            'ballgown_mode': True,
            'skip_reads_with_no_ref': True,
            'gtf_file': 'gtf_file',
            'cov_refs_file': 'cov_refs_file',
            }

        expect_command = '/kb/deployment/bin/StringTie/stringtie '
        expect_command += '-p 4 -B   -C cov_refs_file -e   -G gtf_file -m 100 '
        expect_command += '-o output_transcripts_file -A gene_abundances_file '
        expect_command += '-f 0.6 -j 0.8 -a 8 -t   -c 1.6 -l Lable -g 60 -M 0.8 input_file '

        command = self.stringtie_runner._generate_command(command_params)
        self.assertEquals(command, expect_command)

    @patch.object(StringTieUtil, "_get_input_file", side_effect=mock_get_input_file)
    def test_run_stringtie_app_single_file(self, _get_input_file):

        input_params = {
            'alignment_ref': self.alignment_ref,
            'expression_object_name': 'MyExpression',
            'workspace_name': self.getWsName(),
            'run_matrix_count': True,
            "min_read_coverage": 2.5,
            "junction_base": 10,
            "num_threads": 2,
            "min_isoform_abundance": 0.1,
            "min_length": 200,
            "skip_reads_with_no_ref": 1,
            "merge": 0,
            "junction_coverage": 1,
            "ballgown_mode": 1,
            "min_locus_gap_sep_value": 50
        }

        result = self.getImpl().run_stringtie_app(self.getContext(), input_params)[0]

        self.assertTrue('result_directory' in result)
        result_files = os.listdir(result['result_directory'])
        expect_result_files = ['gene_count_matrix.csv', 'genes.fpkm_tracking',
                               'transcript_count_matrix.csv', 'transcripts.gtf']
        self.assertTrue(all(x in result_files for x in expect_result_files))
        self.assertTrue('expression_obj_ref' in result)
        self.assertTrue('report_name' in result)
        self.assertTrue('report_ref' in result)
