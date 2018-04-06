# -*- coding: utf-8 -*-
import shutil
import unittest
import os  # noqa: F401

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from biokbase.workspace.client import Workspace as workspaceService
from Workspace.WorkspaceClient import Workspace as Workspace
from kb_stringtie.kb_stringtieImpl import kb_stringtie
from kb_stringtie.kb_stringtieServer import MethodContext
from kb_stringtie.authclient import KBaseAuth as _KBaseAuth
from kb_stringtie.Utils.StringTieUtil import StringTieUtil
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
from ReadsUtils.ReadsUtilsClient import ReadsUtils
from ReadsAlignmentUtils.ReadsAlignmentUtilsClient import ReadsAlignmentUtils
from DataFileUtil.DataFileUtilClient import DataFileUtil


class kb_stringtie_ciTest(unittest.TestCase):

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
        cls.ws = Workspace(cls.wsURL, token=token)
        cls.serviceImpl = kb_stringtie(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']

        cls.gfu = GenomeFileUtil(cls.callback_url)
        cls.dfu = DataFileUtil(cls.callback_url)
        cls.ru = ReadsUtils(cls.callback_url)
        cls.rau = ReadsAlignmentUtils(cls.callback_url)

        cls.stringtie_runner = StringTieUtil(cls.cfg)

        cls.prepare_data()

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    @classmethod
    def prepare_data(cls):
        """
        cls.genome_ref_1 = "30957/3/1"
        cls.condition_1 = 'Ecoli_WT'
        cls.alignment_ref_1 = "30957/13/2"
        cls.sample_set_ref_1 = "30957/12/1"
        cls.reads_alignment_set_ref_1 = "30957/19/1"
        """
        cls.genome_ref_1 = "30996/15/1"
        cls.condition_1 = 'WT_R1'
        cls.alignment_ref_1 = "30996/17/1"
        cls.sample_set_ref_1 = "30996/9/1"
        cls.reads_alignment_set_ref_1 = "30996/21/1"


    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        #return "jjeffryes:narrative_1522765762916"
        return "sunita:narrative_1522869877663"

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def test_update_genome_with_novel_isoforms(self):
        return
        gff_path = self.scratch+'/example_stringtie_merge.gff'
        shutil.copy('data/example_stringtie_merge.gff', gff_path)
        print self.stringtie_runner._update_genome_with_novel_isoforms(
            self.getWsName(), self.genome_ref, gff_path)

    def test_run_stringtie_app_alignment(self):
        return
        input_params = {
            'alignment_object_ref': self.alignment_ref_1,
            'workspace_name': self.getWsName(),
            'expression_suffix': '_stringtie_expression',
            'expression_set_suffix': '_stringtie_expression_set',

            "min_read_coverage": 2.5,
            "junction_base": 10,
            "num_threads": 2,
            "min_isoform_abundance": 0.1,
            "min_length": 200,
            "skip_reads_with_no_ref": 1,
            "merge": 0,
            "junction_coverage": 1,
            "ballgown_mode": 1,
            "min_locus_gap_sep_value": 50,
            "disable_trimming": 1
        }

        result = self.getImpl().run_stringtie_app(self.getContext(), input_params)[0]

        self.assertTrue('result_directory' in result)
        result_files = os.listdir(result['result_directory'])
        print(result_files)
        expect_result_files = ['genes.fpkm_tracking', 'transcripts.gtf',
                               'e2t.ctab', 'e_data.ctab', 'i2t.ctab', 'i_data.ctab', 't_data.ctab']
        self.assertTrue(all(x in result_files for x in expect_result_files))
        self.assertTrue('expression_obj_ref' in result)
        self.assertTrue('report_name' in result)
        self.assertTrue('report_ref' in result)
        expression_data = self.ws.get_objects2({'objects': 
                                               [{'ref': result.get('expression_obj_ref')}]}
                                               )['data'][0]['data']
        self.assertEqual(expression_data.get('genome_id'), self.genome_ref)
        self.assertEqual(expression_data.get('condition'), self.condition_1)

    def test_run_stringtie_app_novel_isoform(self):

        input_params = {
            'alignment_object_ref': self.reads_alignment_set_ref_1,
            'workspace_name': self.getWsName(),
            'expression_suffix': '_stringtie_expression',
            'expression_set_suffix': '_stringtie_expression_set',
            'mode': 'novel_isoform',

            "min_read_coverage": 2.5,
            "junction_base": 10,
            "num_threads": 2,
            "min_isoform_abundance": 0.1,
            "min_length": 200,
            "junction_coverage": 1,
            "min_locus_gap_sep_value": 50,
            "disable_trimming": 1,
            "ballgown_mode": 1,
            "merge": 1
        }

        result = self.getImpl().run_stringtie_app(self.getContext(), input_params)[0]

        self.assertTrue('result_directory' in result)
        result_dirs = os.listdir(result['result_directory'])
        print(result_dirs)
        self.assertTrue('merge_result' in result_dirs)
        self.assertTrue('expression_obj_ref' in result)
        self.assertTrue('' == result['expression_obj_ref'])
        self.assertTrue('report_name' in result)
        self.assertTrue('report_ref' in result)
