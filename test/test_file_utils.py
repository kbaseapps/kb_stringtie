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
        file_utils._make_gff(temp_data+"/stringtie_merge.gtf")

    def test_update_t_data(self):
        file_utils._update_t_data(temp_data)
        assert filecmp.cmp(temp_data+"/t_data.ctab",
                           temp_data+"/expected_t_data.ctab")

    def test_update_transcripts(self):
        file_utils._update_transcripts(temp_data)
        assert filecmp.cmp(temp_data+"/transcripts.gtf",
                           temp_data+"/expected_transcripts.gtf")
