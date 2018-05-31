import filecmp
import shutil
import unittest

from kb_stringtie.Utils import file_utils

temp_data = '/kb/module/work/tmp/temp_data'


class KBFileUtilsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        shutil.copytree('data', temp_data)

    def test_make_gff(self):
        file_utils._make_gff(temp_data+"/gffcmp.annotated.gtf",
                             temp_data + '/GCF_000005845.2_ASM584v2_genomic.gff',
                             'novel_isoform')

    def test_update_t_data(self):
        file_utils._update_t_data(temp_data)
        self.assertMultiLineEqual(open(temp_data+"/t_data.ctab").read(),
                                  open(temp_data+"/expected_t_data.ctab").read())

    def test_update_transcripts(self):
        file_utils._update_transcripts(temp_data)
        self.assertMultiLineEqual(open(temp_data+"/transcripts.gtf").read(),
                                  open(temp_data+"/expected_transcripts.gtf").read())
