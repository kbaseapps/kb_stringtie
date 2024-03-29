# -*- coding: utf-8 -*-
############################################################
#
# Autogenerated by the KBase type compiler -
# any changes made here will be overwritten
#
############################################################

from __future__ import print_function
# the following is a hack to get the baseclient to import whether we're in a
# package or not. This makes pep8 unhappy hence the annotations.
try:
    # baseclient and this client are in a package
    from .baseclient import BaseClient as _BaseClient  # @UnusedImport
except ImportError:
    # no they aren't
    from baseclient import BaseClient as _BaseClient  # @Reimport


class kb_stringtie(object):

    def __init__(
            self, url=None, timeout=30 * 60, user_id=None,
            password=None, token=None, ignore_authrc=False,
            trust_all_ssl_certificates=False,
            auth_svc='https://ci.kbase.us/services/auth/api/legacy/KBase/Sessions/Login',
            service_ver='dev',
            async_job_check_time_ms=100, async_job_check_time_scale_percent=150, 
            async_job_check_max_time_ms=300000):
        if url is None:
            raise ValueError('A url is required')
        self._service_ver = service_ver
        self._client = _BaseClient(
            url, timeout=timeout, user_id=user_id, password=password,
            token=token, ignore_authrc=ignore_authrc,
            trust_all_ssl_certificates=trust_all_ssl_certificates,
            auth_svc=auth_svc,
            async_job_check_time_ms=async_job_check_time_ms,
            async_job_check_time_scale_percent=async_job_check_time_scale_percent,
            async_job_check_max_time_ms=async_job_check_max_time_ms)

    def run_stringtie_app(self, params, context=None):
        """
        run_stringtie_app: run StringTie app
        ref: http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
        :param params: instance of type "StringTieInput" (required params:
           alignment_object_ref: Alignment or AlignmentSet object reference
           workspace_name: the name of the workspace it gets saved to
           expression_set_suffix: suffix append to expression set object name
           expression_suffix: suffix append to expression object name mode:
           one of ['normal', 'merge', 'novel_isoform'] optional params:
           num_threads: number of processing threads junction_base: junctions
           that don't have spliced reads junction_coverage: junction coverage
           disable_trimming: disables trimming at the ends of the assembled
           transcripts min_locus_gap_sep_value: minimum locus gap separation
           value ballgown_mode: enables the output of Ballgown input table
           files skip_reads_with_no_ref: reads with no reference will be
           skipped novel_isoforms: output expression matrices with novel
           isoforms maximum_fraction: maximum fraction of
           muliple-location-mapped reads min_length: minimum length allowed
           for the predicted transcripts min_read_coverage: minimum input
           transcript coverage min_isoform_abundance: minimum isoform
           abundance ref:
           http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual) ->
           structure: parameter "alignment_object_ref" of type "obj_ref" (An
           X/Y/Z style reference), parameter "workspace_name" of String,
           parameter "expression_set_suffix" of String, parameter
           "expression_suffix" of String, parameter "num_threads" of Long,
           parameter "junction_base" of Long, parameter "junction_coverage"
           of Double, parameter "disable_trimming" of type "boolean" (A
           boolean - 0 for false, 1 for true. @range (0, 1)), parameter
           "min_locus_gap_sep_value" of Long, parameter "ballgown_mode" of
           type "boolean" (A boolean - 0 for false, 1 for true. @range (0,
           1)), parameter "skip_reads_with_no_ref" of type "boolean" (A
           boolean - 0 for false, 1 for true. @range (0, 1)), parameter
           "novel_isoforms" of type "NovelIsoformParams"
           (stringtie_genome_name: name for the new genome including novel
           transcripts transcript_label: prefix for the name of the output
           transcripts) -> structure: parameter "label" of String, parameter
           "stringtie_genome_name" of String, parameter "maximum_fraction" of
           Double, parameter "min_length" of Long, parameter
           "min_read_coverage" of Double, parameter "min_isoform_abundance"
           of Double
        :returns: instance of type "StringTieResult" (result_directory:
           folder path that holds all files generated by run_stringtie
           expression_obj_ref: generated Expression/ExpressionSet object
           reference exprMatrix_FPKM/TPM_ref: generated FPKM/TPM
           ExpressionMatrix object reference report_name: report name
           generated by KBaseReport report_ref: report reference generated by
           KBaseReport) -> structure: parameter "result_directory" of String,
           parameter "expression_obj_ref" of type "obj_ref" (An X/Y/Z style
           reference), parameter "exprMatrix_FPKM_ref" of type "obj_ref" (An
           X/Y/Z style reference), parameter "exprMatrix_TPM_ref" of type
           "obj_ref" (An X/Y/Z style reference), parameter "report_name" of
           String, parameter "report_ref" of String
        """
        return self._client.run_job('kb_stringtie.run_stringtie_app',
                                    [params], self._service_ver, context)

    def status(self, context=None):
        return self._client.run_job('kb_stringtie.status',
                                    [], self._service_ver, context)
