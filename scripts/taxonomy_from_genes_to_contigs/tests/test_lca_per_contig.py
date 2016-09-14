#!/usr/bin/env python
from nose.tools import assert_equal, assert_true, assert_almost_equal, nottest, assert_false
from os.path import isdir,isfile
from os import listdir
import os
import sys
import subprocess
import pandas as pd

file_path = os.path.realpath(__file__)
file_dir_path = os.path.dirname(file_path)
data_path = os.path.abspath(os.path.join(file_dir_path,"test_data/"))
tmp_dir_path = file_dir_path + '/nose_tmp_output'

CWD = os.getcwd()

class TestCMD(object):
    def setUp(self):
        """Create temporary dir if necessary,
        otherwise clear contents of it"""
        if not isdir(tmp_dir_path):
            os.mkdir(tmp_dir_path)
        self.tearDown()
        os.chdir(file_dir_path)
        self.data_path = data_path
        self.tmp_dir_path = tmp_dir_path

    def tearDown(self):
        """remove temporary output files"""
        for d in os.listdir(tmp_dir_path):
            d_path = os.path.join(tmp_dir_path,d)
            try:
                os.remove(d_path)
            except:
                for f in os.listdir(d_path):
                    f_path = os.path.join(d_path,f)
                    os.remove(f_path)
                os.rmdir(d_path)
        assert os.listdir(tmp_dir_path) == []


    def run_command(self):
        ncbi_tree = os.path.join(self.data_path, "ncbi.tre")
        ncbi_map = os.path.join(self.data_path, "test_ncbi.map")
        ncbi_nodes = os.path.join(self.data_path, "test_nodes.dmp")
        megan_result = os.path.join(self.data_path, "megan_result.csv")
        bed_file = os.path.join(self.data_path, "PROKKA.bed")
        no_lca_out = os.path.join(self.tmp_dir_path, "no_lca_out.csv")
        lca_out = os.path.join(self.tmp_dir_path, "lca_out.csv")
        call_string = ("python ../lca_per_contig.py --ncbi_tree_file {0} --ncbi_map_file {1} "
                    "--nodes_file {2} --megan_result {3} --bed_file {4} --no_lca_out {5} "
                    "--lca_out {6}").format(ncbi_tree, ncbi_map, ncbi_nodes, megan_result, \
                            bed_file, no_lca_out, lca_out)


        self.c = 0 # Exit code
        try:
            self.op = subprocess.check_output(
                call_string,
                shell=True)
        except subprocess.CalledProcessError as exc:
            self.c = exc.returncode

    def test_everything(self):
        ## Actinobacteria, phylum, id: 201174
        #ncbi_map.ix[201174]
        ## Proteobacteria, phylum, id: 1224
        #ncbi_map.ix[1224]['name']
        ## Alphaproteobacteria, class, id: 28211
        #ncbi_map.ix[28211]['name']
        ## Betaproteobacteria, class, id: 28216
        #ncbi_map.ix[28216]['name']
        ## Deltavirus, genus, id: 39759
        #ncbi_map.ix[39759]

        self.run_command()
        assert_equal(self.c,0,
                     msg = "Command exited with nonzero status")


        # File creation
        no_lca_file = os.path.join(self.tmp_dir_path, "no_lca_out.csv")
        assert_true(isfile(no_lca_file),
                    msg = "no lca out file is not created")

        lca_file = os.path.join(self.tmp_dir_path, "lca_out.csv")
        assert_true(isfile(lca_file),
                    msg = "lca out file is not created")

        lca_df = pd.read_table(lca_file, sep=',')
        # A case where no gene is assigned to any taxa
        # k99_100001
        assert_true(lca_df[lca_df['contig_id'] == 'k99_10001']['superkingdom'].isnull().all())

        # A case where only one gene is assigned to a taxa and no other genes on that contig
        # k99_100002
        assert_false(lca_df[lca_df['contig_id'] == 'k99_10002']['superkingdom'].isnull().any())
        assert_true((lca_df[lca_df['contig_id'] == 'k99_10002']['superkingdom'] == 'Bacteria').all())
        assert_true((lca_df[lca_df['contig_id'] == 'k99_10002']['phylum'] == 'Proteobacteria').all())

        # A case where only one gene is assigned to a taxa and but one other gene on that contig
        # k99_100003
        assert_false(lca_df[lca_df['contig_id'] == 'k99_10003']['superkingdom'].isnull().any())
        assert_true((lca_df[lca_df['contig_id'] == 'k99_10003']['superkingdom'] == 'Bacteria').all())
        assert_true((lca_df[lca_df['contig_id'] == 'k99_10003']['phylum'] == 'Proteobacteria').all())

        # A case where two genes are assigned to a taxa and they have the same taxa
        # k99_100004
        assert_false(lca_df[lca_df['contig_id'] == 'k99_10004']['superkingdom'].isnull().any())
        assert_true((lca_df[lca_df['contig_id'] == 'k99_10004']['superkingdom'] == 'Bacteria').all())
        assert_true((lca_df[lca_df['contig_id'] == 'k99_10004']['phylum'] == 'Proteobacteria').all())
        assert_true((lca_df[lca_df['contig_id'] == 'k99_10004']['class'] == 'Alphaproteobacteria').all())

        # A case where two genes are assigned to a taxa and they have a common ancestor
        # k99_10005
        assert_false(lca_df[lca_df['contig_id'] == 'k99_10005']['superkingdom'].isnull().any())
        assert_true((lca_df[lca_df['contig_id'] == 'k99_10005']['superkingdom'] == 'Bacteria').all())
        assert_true((lca_df[lca_df['contig_id'] == 'k99_10005']['phylum'] == 'Proteobacteria').all())
        assert_false(lca_df[lca_df['contig_id'] == 'k99_10005']['class'].isnull().any())

        # A case where two genes are assigned to a taxa and they have no common ancestor
        # k99_10006
        assert_true(lca_df[lca_df['contig_id'] == 'k99_100006']['superkingdom'].isnull().all())
        assert_true(lca_df[lca_df['contig_id'] == 'k99_100006']['phylum'].isnull().all())
        assert_true(lca_df[lca_df['contig_id'] == 'k99_100006']['class'].isnull().all())
        assert_true(lca_df[lca_df['contig_id'] == 'k99_100006']['order'].isnull().all())
        assert_true(lca_df[lca_df['contig_id'] == 'k99_100006']['family'].isnull().all())
        assert_true(lca_df[lca_df['contig_id'] == 'k99_100006']['genus'].isnull().all())
        assert_true(lca_df[lca_df['contig_id'] == 'k99_100006']['species'].isnull().all())

        # A case where ten out of eleven genes are assigned to a taxa and 9 have the same taxa and the tenth is a relative
        # k99_10007

        assert_false(lca_df[lca_df['contig_id'] == 'k99_100007']['superkingdom'].isnull().any())
        assert_true((lca_df[lca_df['contig_id'] == 'k99_100007']['superkingdom'] == 'Bacteria').all())
        assert_true((lca_df[lca_df['contig_id'] == 'k99_100007']['phylum'] == 'Proteobacteria').all())
        assert_true((lca_df[lca_df['contig_id'] == 'k99_100007']['class'] == 'Alphaproteobacteria').all())

        # A case where ten out of ten genes are assigned to a taxa and 9 share a common ancestor closer than the tenth.
        # k99_10008
        # Bacteria
        # |
        # --201174
        # |
        # --1224
        #   |
        #   |---28211
        #   |---28216

        assert_false(lca_df[lca_df['contig_id'] == 'k99_100008']['superkingdom'].isnull().any())
        assert_true((lca_df[lca_df['contig_id'] == 'k99_100008']['superkingdom'] == 'Bacteria').all())
        assert_true((lca_df[lca_df['contig_id'] == 'k99_100008']['phylum'] == 'Proteobacteria').all())
        assert_true(lca_df[lca_df['contig_id'] == 'k99_100008']['class'].isnull().all())

        # A case where ten out of ten genes are assigned to a taxa but only 8 share a common ancestor.
        # k99_10009
        # Virus
        # Bacteria
        # |
        # --201174
        # |
        # --1224
        #   |
        #   |---28211
        #   |---28216

        assert_true(lca_df[lca_df['contig_id'] == 'k99_100009']['superkingdom'].isnull().all())
        assert_true(lca_df[lca_df['contig_id'] == 'k99_100009']['phylum'].isnull().all())
        assert_true(lca_df[lca_df['contig_id'] == 'k99_100009']['class'].isnull().all())
        assert_true(lca_df[lca_df['contig_id'] == 'k99_100009']['order'].isnull().all())
        assert_true(lca_df[lca_df['contig_id'] == 'k99_100009']['family'].isnull().all())
        assert_true(lca_df[lca_df['contig_id'] == 'k99_100009']['genus'].isnull().all())
        assert_true(lca_df[lca_df['contig_id'] == 'k99_100009']['species'].isnull().all())

        # A case with a strange taxon where certain levels are skipped
        # k99_100010
        assert_false(lca_df[lca_df['contig_id'] == 'k99_100010']['superkingdom'].isnull().any())
        assert_true((lca_df[lca_df['contig_id'] == 'k99_100010']['superkingdom'] == 'Viruses').all())
        assert_true((lca_df[lca_df['contig_id'] == 'k99_100010']['phylum'].isnull()).all())
        assert_true((lca_df[lca_df['contig_id'] == 'k99_100010']['class'].isnull()).all())
        assert_true((lca_df[lca_df['contig_id'] == 'k99_100010']['order'].isnull()).all())
        assert_true((lca_df[lca_df['contig_id'] == 'k99_100010']['family'].isnull()).all())
        assert_true((lca_df[lca_df['contig_id'] == 'k99_100010']['genus'] == 'Deltavirus').all(),
                msg="{}".format(lca_df[lca_df['contig_id'] == 'k99_100010']))

        # A case where root is present, which should be considered as unclassified
        # k99_100011
        assert_false(lca_df[lca_df['contig_id'] == 'k99_100010']['superkingdom'].isnull().any())
        assert_true((lca_df[lca_df['contig_id'] == 'k99_100011']['superkingdom'] == 'Bacteria').all())

