# -*- coding: utf-8 -*-
############################################################
#
# Autogenerated by the KBase type compiler -
# any changes made here will be overwritten
#
############################################################


# the following is a hack to get the baseclient to import whether we're in a
# package or not. This makes pep8 unhappy hence the annotations.
try:
    # baseclient and this client are in a package
    from .baseclient import BaseClient as _BaseClient  # @UnusedImport
except ImportError:
    # no they aren't
    from .baseclient import BaseClient as _BaseClient  # @Reimport


class KBaseReport(object):

    def __init__(
            self, url=None, timeout=30 * 60, user_id=None,
            password=None, token=None, ignore_authrc=False,
            trust_all_ssl_certificates=False,
            auth_svc='https://ci.kbase.us/services/auth/api/legacy/KBase/Sessions/Login',
            service_ver='release',
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

    def create(self, params, context=None):
        """
        Create a KBaseReport with a brief summary of an App run.
        :param params: instance of type "CreateParams" (Provide the report
           information.  The structure is: params = { report: { text_message:
           '', warnings: ['w1'], objects_created: [ { ref: 'ws/objid',
           description: '' }] }, workspace_name: 'ws' }) -> structure:
           parameter "report" of type "Report" (A simple Report of a method
           run in KBase. It only provides for now a way to display a fixed
           width text output summary message, a list of warnings, and a list
           of objects created (each with descriptions). @optional warnings
           file_links html_links direct_html direct_html_link_index @metadata
           ws length(warnings) as Warnings @metadata ws length(text_message)
           as Size(characters) @metadata ws length(objects_created) as
           Objects Created) -> structure: parameter "text_message" of String,
           parameter "warnings" of list of String, parameter
           "objects_created" of list of type "WorkspaceObject" (Represents a
           Workspace object with some brief description text that can be
           associated with the object. @optional description) -> structure:
           parameter "ref" of type "ws_id" (@id ws), parameter "description"
           of String, parameter "file_links" of list of type "LinkedFile"
           (Represents a file or html archive that the report should like to
           @optional description label) -> structure: parameter "handle" of
           type "handle_ref" (Reference to a handle @id handle), parameter
           "description" of String, parameter "name" of String, parameter
           "label" of String, parameter "URL" of String, parameter
           "html_links" of list of type "LinkedFile" (Represents a file or
           html archive that the report should like to @optional description
           label) -> structure: parameter "handle" of type "handle_ref"
           (Reference to a handle @id handle), parameter "description" of
           String, parameter "name" of String, parameter "label" of String,
           parameter "URL" of String, parameter "direct_html" of String,
           parameter "direct_html_link_index" of Long, parameter
           "workspace_name" of String
        :returns: instance of type "ReportInfo" (The reference to the saved
           KBaseReport.  The structure is: reportInfo = { ref:
           'ws/objid/ver', name: 'myreport.2262323452' }) -> structure:
           parameter "ref" of type "ws_id" (@id ws), parameter "name" of
           String
        """
        return self._client.run_job('KBaseReport.create',
                                    [params], self._service_ver, context)

    def create_extended_report(self, params, context=None):
        """
        A more complex function to create a report that enables the user to specify files and html view that the report should link to
        :param params: instance of type "CreateExtendedReportParams"
           (Parameters used to create a more complex report with file and
           html links The following arguments allow the user to specify the
           classical data fields in the report object: string message -
           simple text message to store in report object list
           <WorkspaceObject> objects_created; list <string> warnings - a list
           of warning messages in simple text The following argument allows
           the user to specify the location of html files/directories that
           the report widget will render <or> link to: list <fileRef>
           html_links - a list of paths or shock node IDs pointing to a
           single flat html file or to the top level directory of a website
           The report widget can render one html view directly. Set one of
           the following fields to decide which view to render: string
           direct_html - simple html text that will be rendered within the
           report widget int  direct_html_link_index - use this to specify
           the index of the page in html_links to view directly in the report
           widget (ignored if html_string is set) The following argument
           allows the user to specify the location of files that the report
           widget should link for download: list <fileRef> file_links - a
           list of paths or shock node IDs pointing to a single flat file The
           following parameters indicate where the report object should be
           saved in the workspace: string report_object_name - name to use
           for the report object (job ID is used if left unspecified)
           html_window_height - height of the html window in the narrative
           output widget summary_window_height - height of summary window in
           the narrative output widget string workspace_name - name of
           workspace where object should be saved) -> structure: parameter
           "message" of String, parameter "objects_created" of list of type
           "WorkspaceObject" (Represents a Workspace object with some brief
           description text that can be associated with the object. @optional
           description) -> structure: parameter "ref" of type "ws_id" (@id
           ws), parameter "description" of String, parameter "warnings" of
           list of String, parameter "html_links" of list of type "File" ->
           structure: parameter "path" of String, parameter "shock_id" of
           String, parameter "name" of String, parameter "description" of
           String, parameter "direct_html" of String, parameter
           "direct_html_link_index" of Long, parameter "file_links" of list
           of type "File" -> structure: parameter "path" of String, parameter
           "shock_id" of String, parameter "name" of String, parameter
           "description" of String, parameter "report_object_name" of String,
           parameter "html_window_height" of Double, parameter
           "summary_window_height" of Double, parameter "workspace_name" of
           String
        :returns: instance of type "ReportInfo" (The reference to the saved
           KBaseReport.  The structure is: reportInfo = { ref:
           'ws/objid/ver', name: 'myreport.2262323452' }) -> structure:
           parameter "ref" of type "ws_id" (@id ws), parameter "name" of
           String
        """
        return self._client.run_job('KBaseReport.create_extended_report',
                                    [params], self._service_ver, context)

    def status(self, context=None):
        return self._client.run_job('KBaseReport.status',
                                    [], self._service_ver, context)
