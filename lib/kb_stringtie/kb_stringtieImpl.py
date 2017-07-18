# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import json

from kb_stringtie.Utils.StringTieUtil import StringTieUtil
#END_HEADER


class kb_stringtie:
    '''
    Module Name:
    kb_stringtie

    Module Description:
    A KBase module: kb_stringtie
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "1.0.2"
    GIT_URL = "https://github.com/Tianhao-Gu/kb_stringtie.git"
    GIT_COMMIT_HASH = "01a9ae483a6b8973a2528dc73bf8eb97f5c3eb99"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.config['SDK_CALLBACK_URL'] = os.environ['SDK_CALLBACK_URL']
        self.config['KB_AUTH_TOKEN'] = os.environ['KB_AUTH_TOKEN']
        #END_CONSTRUCTOR
        pass


    def run_stringtie_app(self, ctx, params):
        """
        run_stringtie_app: run StringTie app
        ref: http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
        :param params: instance of type "StringTieInput" (required params:
           alignment_object_ref: Alignment or AlignmentSet object reference
           workspace_name: the name of the workspace it gets saved to
           expression_set_suffix: suffix append to expression set object name
           expression_suffix: suffix append to expression object name
           optional params: num_threads: number of processing threads
           junction_base: junctions that don't have spliced reads
           junction_coverage: junction coverage disable_trimming: disables
           trimming at the ends of the assembled transcripts
           min_locus_gap_sep_value: minimum locus gap separation value
           ballgown_mode: enables the output of Ballgown input table files
           skip_reads_with_no_ref: reads with no reference will be skipped
           maximum_fraction: maximum fraction of muliple-location-mapped
           reads label: prefix for the name of the output transcripts
           min_length: minimum length allowed for the predicted transcripts
           min_read_coverage: minimum input transcript coverage
           min_isoform_abundance: minimum isoform abundance merge: set
           transcript merge mode ref:
           http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual) ->
           structure: parameter "alignment_object_ref" of type "obj_ref" (An
           X/Y/Z style reference), parameter "workspace_name" of String,
           parameter "expression_set_suffix" of String, parameter
           "expression_suffix" of String, parameter "merge" of type "boolean"
           (A boolean - 0 for false, 1 for true. @range (0, 1)), parameter
           "num_threads" of Long, parameter "junction_base" of Long,
           parameter "junction_coverage" of Double, parameter
           "disable_trimming" of type "boolean" (A boolean - 0 for false, 1
           for true. @range (0, 1)), parameter "min_locus_gap_sep_value" of
           Long, parameter "ballgown_mode" of type "boolean" (A boolean - 0
           for false, 1 for true. @range (0, 1)), parameter
           "skip_reads_with_no_ref" of type "boolean" (A boolean - 0 for
           false, 1 for true. @range (0, 1)), parameter "maximum_fraction" of
           Double, parameter "label" of String, parameter "min_length" of
           Long, parameter "min_read_coverage" of Double, parameter
           "min_isoform_abundance" of Double
        :returns: instance of type "StringTieResult" (result_directory:
           folder path that holds all files generated by run_stringtie
           expression_obj_ref: generated Expression/ExpressionSet object
           reference report_name: report name generated by KBaseReport
           report_ref: report reference generated by KBaseReport) ->
           structure: parameter "result_directory" of String, parameter
           "expression_obj_ref" of type "obj_ref" (An X/Y/Z style reference),
           parameter "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN run_stringtie_app
        print '--->\nRunning kb_stringtie.run_stringtie\nparams:'
        print json.dumps(params, indent=1)

        for key, value in params.iteritems():
            if isinstance(value, basestring):
                params[key] = value.strip()

        stringtie_runner = StringTieUtil(self.config)
        returnVal = stringtie_runner.run_stringtie_app(params)
        #END run_stringtie_app

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method run_stringtie_app return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
