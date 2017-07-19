# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import json  # noqa: F401
import time
import requests  # noqa: F401
import shutil

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint  # noqa: F401

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
        cls.ws = Workspace(cls.wsURL, token=token)
        cls.serviceImpl = kb_stringtie(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']

        cls.gfu = GenomeFileUtil(cls.callback_url)
        cls.dfu = DataFileUtil(cls.callback_url)
        cls.ru = ReadsUtils(cls.callback_url)
        cls.rau = ReadsAlignmentUtils(cls.callback_url, service_ver='dev')

        cls.stringtie_runner = StringTieUtil(cls.cfg)

        suffix = int(time.time() * 1000)
        cls.wsName = "test_kb_stringtie_" + str(suffix)
        cls.wsClient.create_workspace({'workspace': cls.wsName})

        cls.prepare_data()

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    @classmethod
    def prepare_data(cls):
        # upload genome object
        genbank_file_name = 'minimal.gbff'
        genbank_file_path = os.path.join(cls.scratch, genbank_file_name)
        shutil.copy(os.path.join('data', genbank_file_name), genbank_file_path)

        genome_object_name = 'test_Genome'
        cls.genome_ref = cls.gfu.genbank_to_genome({'file': {'path': genbank_file_path},
                                                    'workspace_name': cls.wsName,
                                                    'genome_name': genome_object_name
                                                    })['genome_ref']

        # upload reads object
        reads_file_name = 'Sample1.fastq'
        reads_file_path = os.path.join(cls.scratch, reads_file_name)
        shutil.copy(os.path.join('data', reads_file_name), reads_file_path)

        reads_object_name_1 = 'test_Reads_1'
        cls.reads_ref_1 = cls.ru.upload_reads({'fwd_file': reads_file_path,
                                               'wsname': cls.wsName,
                                               'sequencing_tech': 'Unknown',
                                               'interleaved': 0,
                                               'name': reads_object_name_1
                                               })['obj_ref']

        reads_object_name_2 = 'test_Reads_2'
        cls.reads_ref_2 = cls.ru.upload_reads({'fwd_file': reads_file_path,
                                               'wsname': cls.wsName,
                                               'sequencing_tech': 'Unknown',
                                               'interleaved': 0,
                                               'name': reads_object_name_2
                                               })['obj_ref']

        # upload alignment object
        alignment_file_name = 'accepted_hits.bam'
        alignment_file_path = os.path.join(cls.scratch, alignment_file_name)
        shutil.copy(os.path.join('data', alignment_file_name), alignment_file_path)

        alignment_object_name_1 = 'test_Alignment_1'
        cls.condition_1 = 'test_condition_1'
        destination_ref = cls.wsName + '/' + alignment_object_name_1
        cls.alignment_ref_1 = cls.rau.upload_alignment({'file_path': alignment_file_path,
                                                        'destination_ref': destination_ref,
                                                        'read_library_ref': cls.reads_ref_1,
                                                        'condition': cls.condition_1,
                                                        'library_type': 'single_end',
                                                        'assembly_or_genome_ref': cls.genome_ref
                                                        })['obj_ref']

        alignment_object_name_2 = 'test_Alignment_2'
        cls.condition_2 = 'test_condition_2'
        destination_ref = cls.wsName + '/' + alignment_object_name_2
        cls.alignment_ref_2 = cls.rau.upload_alignment({'file_path': alignment_file_path,
                                                        'destination_ref': destination_ref,
                                                        'read_library_ref': cls.reads_ref_2,
                                                        'condition': cls.condition_2,
                                                        'library_type': 'single_end',
                                                        'assembly_or_genome_ref': cls.genome_ref
                                                        })['obj_ref']

        # upload sample_set object
        workspace_id = cls.dfu.ws_name_to_id(cls.wsName)
        sample_set_object_name = 'test_Sample_Set'
        sample_set_data = {'sampleset_id': sample_set_object_name,
                           'sampleset_desc': 'test sampleset object',
                           'Library_type': 'SingleEnd',
                           'condition': [cls.condition_1, cls.condition_2],
                           'domain': 'Unknown',
                           'num_samples': 2,
                           'platform': 'Unknown'}
        save_object_params = {
            'id': workspace_id,
            'objects': [{'type': 'KBaseRNASeq.RNASeqSampleSet',
                         'data': sample_set_data,
                         'name': sample_set_object_name}]
        }

        dfu_oi = cls.dfu.save_objects(save_object_params)[0]
        cls.sample_set_ref = str(dfu_oi[6]) + '/' + str(dfu_oi[0]) + '/' + str(dfu_oi[4])

        # upload RNASeqAlignmentSet object
        object_type = 'KBaseRNASeq.RNASeqAlignmentSet'
        alignment_set_object_name = 'test_RNASeq_Alignment_Set'
        alignment_set_data = {'genome_id': cls.genome_ref,
                              'read_sample_ids': [reads_object_name_1, reads_object_name_2],
                              'mapped_rnaseq_alignments': [{reads_object_name_1: 
                                                            alignment_object_name_1},
                                                           {reads_object_name_2: 
                                                            alignment_object_name_2}],
                              'mapped_alignments_ids': [{reads_object_name_1: cls.alignment_ref_1},
                                                        {reads_object_name_2: cls.alignment_ref_2}
                                                        ],
                              'sample_alignments': [cls.alignment_ref_1, cls.alignment_ref_2],
                              'sampleset_id': cls.sample_set_ref}
        save_object_params = {
            'id': workspace_id,
            'objects': [{'type': object_type,
                         'data': alignment_set_data,
                         'name': alignment_set_object_name}]
        }

        dfu_oi = cls.dfu.save_objects(save_object_params)[0]
        cls.rnaseq_alignment_set_ref = str(dfu_oi[6]) + '/' + str(dfu_oi[0]) + '/' + str(dfu_oi[4])

        # upload ReadsAlignmentSet object
        object_type = 'KBaseSets.ReadsAlignmentSet'
        alignment_set_object_name = 'test_reads_Alignment_Set'
        alignment_set_data = {'description': 'test ReadsAlignmentSet object',
                              'items': [{'ref': cls.alignment_ref_1},
                                        {'ref': cls.alignment_ref_2}]}
        save_object_params = {
            'id': workspace_id,
            'objects': [{'type': object_type,
                         'data': alignment_set_data,
                         'name': alignment_set_object_name}]
        }

        dfu_oi = cls.dfu.save_objects(save_object_params)[0]
        cls.reads_alignment_set_ref = str(dfu_oi[6]) + '/' + str(dfu_oi[0]) + '/' + str(dfu_oi[4])

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        return self.__class__.wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    # def test_bad_run_stringtie_app_params(self):
    #     invalidate_input_params = {'missing_alignment_object_ref': 'alignment_object_ref',
    #                                'workspace_name': 'workspace_name'}
    #     with self.assertRaisesRegexp(ValueError, 
    #                                  '"alignment_object_ref" parameter is required, but missing'):
    #         self.getImpl().run_stringtie_app(self.getContext(), invalidate_input_params)

    #     invalidate_input_params = {'alignment_object_ref': 'alignment_object_ref',
    #                                'missing_workspace_name': 'workspace_name'}
    #     with self.assertRaisesRegexp(ValueError, 
    #                                  '"workspace_name" parameter is required, but missing'):
    #         self.getImpl().run_stringtie_app(self.getContext(), invalidate_input_params)

    # def test_StringTieUtil_generate_command(self):
    #     command_params = {
    #         'num_threads': 4,
    #         'junction_base': 8,
    #         'junction_coverage': 0.8,
    #         'disable_trimming': True,
    #         'min_locus_gap_sep_value': 60,
    #         'maximum_fraction': 0.8,
    #         'label': 'Lable',
    #         'min_length': 100,
    #         'min_read_coverage': 1.6,
    #         'min_isoform_abundance': 0.6,
    #         'output_transcripts': 'output_transcripts_file',
    #         'gene_abundances_file': 'gene_abundances_file',
    #         'input_file': 'input_file',
    #         'ballgown_mode': True,
    #         'skip_reads_with_no_ref': True,
    #         'gtf_file': 'gtf_file',
    #         'cov_refs_file': 'cov_refs_file'
    #     }

    #     expect_command = '/kb/deployment/bin/StringTie/stringtie '
    #     expect_command += '-p 4 -B   -C cov_refs_file -e   -G gtf_file -m 100 '
    #     expect_command += '-o output_transcripts_file -A gene_abundances_file '
    #     expect_command += '-f 0.6 -j 0.8 -a 8 -t   -c 1.6 -l Lable -g 60 -M 0.8 input_file '

    #     command = self.stringtie_runner._generate_command(command_params)
    #     self.assertEquals(command, expect_command)

    # def test_run_stringtie_app_alignment(self):

    #     input_params = {
    #         'alignment_object_ref': self.alignment_ref_1,
    #         'workspace_name': self.getWsName(),
    #         'expression_suffix': '_stringtie_expression',
    #         'expression_set_suffix': '_stringtie_expression_set',

    #         "min_read_coverage": 2.5,
    #         "junction_base": 10,
    #         "num_threads": 4,
    #         "min_isoform_abundance": 0.1,
    #         "min_length": 200,
    #         "skip_reads_with_no_ref": 1,
    #         "merge": 0,
    #         "junction_coverage": 1,
    #         "ballgown_mode": 1,
    #         "min_locus_gap_sep_value": 50,
    #         "disable_trimming": 1
    #     }

    #     result = self.getImpl().run_stringtie_app(self.getContext(), input_params)[0]

    #     self.assertTrue('result_directory' in result)
    #     result_files = os.listdir(result['result_directory'])
    #     print result_files
    #     expect_result_files = ['genes.fpkm_tracking', 'transcripts.gtf',
    #                            'e2t.ctab', 'e_data.ctab', 'i2t.ctab', 'i_data.ctab', 't_data.ctab']
    #     self.assertTrue(all(x in result_files for x in expect_result_files))
    #     self.assertTrue('expression_obj_ref' in result)
    #     self.assertTrue('report_name' in result)
    #     self.assertTrue('report_ref' in result)
    #     expression_data = self.ws.get_objects2({'objects': 
    #                                            [{'ref': result.get('expression_obj_ref')}]}
    #                                            )['data'][0]['data']
    #     self.assertEqual(expression_data.get('genome_id'), self.genome_ref)
    #     self.assertEqual(expression_data.get('condition'), self.condition_1)
    #     self.assertEqual(expression_data.get('id'), 'test_stringtie_expression_1')

    # def test_run_stringtie_app_rnaseq_alignment_set(self):

    #     input_params = {
    #         'alignment_object_ref': self.rnaseq_alignment_set_ref,
    #         'workspace_name': self.getWsName(),
    #         'expression_suffix': '_stringtie_expression',
    #         'expression_set_suffix': '_stringtie_expression_set',

    #         "min_read_coverage": 2.5,
    #         "junction_base": 10,
    #         "num_threads": 4,
    #         "min_isoform_abundance": 0.1,
    #         "min_length": 200,
    #         "skip_reads_with_no_ref": 1,
    #         "merge": 0,
    #         "junction_coverage": 1,
    #         "ballgown_mode": 1,
    #         "min_locus_gap_sep_value": 50,
    #         "disable_trimming": 1
    #     }

    #     result = self.getImpl().run_stringtie_app(self.getContext(), input_params)[0]

    #     self.assertTrue('result_directory' in result)
    #     result_files = os.listdir(result['result_directory'])
    #     print result_files
    #     self.assertTrue('expression_obj_ref' in result)
    #     self.assertTrue('report_name' in result)
    #     self.assertTrue('report_ref' in result)
    #     expression_data = self.ws.get_objects2({'objects': 
    #                                            [{'ref': result.get('expression_obj_ref')}]}
    #                                            )['data'][0]['data']
    #     self.assertTrue('items' in expression_data)
    #     self.assertTrue('description' in expression_data)

    # def test_run_stringtie_app_reads_alignment_set(self):

    #     input_params = {
    #         'alignment_object_ref': self.reads_alignment_set_ref,
    #         'workspace_name': self.getWsName(),
    #         'expression_suffix': '_stringtie_expression',
    #         'expression_set_suffix': '_stringtie_expression_set',

    #         "min_read_coverage": 2.5,
    #         "junction_base": 10,
    #         "num_threads": 4,
    #         "min_isoform_abundance": 0.1,
    #         "min_length": 200,
    #         "skip_reads_with_no_ref": 1,
    #         "merge": 0,
    #         "junction_coverage": 1,
    #         "ballgown_mode": 1,
    #         "min_locus_gap_sep_value": 50,
    #         "disable_trimming": 1
    #     }

    #     result = self.getImpl().run_stringtie_app(self.getContext(), input_params)[0]

    #     self.assertTrue('result_directory' in result)
    #     result_files = os.listdir(result['result_directory'])
    #     print result_files
    #     self.assertTrue('expression_obj_ref' in result)
    #     self.assertTrue('report_name' in result)
    #     self.assertTrue('report_ref' in result)
    #     expression_data = self.ws.get_objects2({'objects': 
    #                                            [{'ref': result.get('expression_obj_ref')}]}
    #                                            )['data'][0]['data']
    #     self.assertTrue('items' in expression_data)
    #     self.assertTrue('description' in expression_data)

    def test_run_stringtie_app_merge(self):

        input_params = {
            'alignment_object_ref': self.reads_alignment_set_ref,
            'workspace_name': self.getWsName(),
            'expression_suffix': '_stringtie_expression',
            'expression_set_suffix': '_stringtie_expression_set',
            'merge': 1,

            "min_read_coverage": 2.5,
            "junction_base": 10,
            "num_threads": 4,
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
        print result['result_directory']
        print result_files
        self.assertTrue('expression_obj_ref' in result)
        self.assertTrue('report_name' in result)
        self.assertTrue('report_ref' in result)
        expression_data = self.ws.get_objects2({'objects': 
                                               [{'ref': result.get('expression_obj_ref')}]}
                                               )['data'][0]['data']
        self.assertTrue('items' in expression_data)
        self.assertTrue('description' in expression_data)
