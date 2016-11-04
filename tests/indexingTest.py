# =============================================================================
# AUSTRALIAN NATIONAL UNIVERSITY OPEN SOURCE LICENSE (ANUOS LICENSE)
# VERSION 1.3
# 
# The contents of this file are subject to the ANUOS License Version 1.3
# (the "License"); you may not use this file except in compliance with
# the License. You may obtain a copy of the License at:
# 
#   https://sourceforge.net/projects/febrl/
# 
# Software distributed under the License is distributed on an "AS IS"
# basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See
# the License for the specific language governing rights and limitations
# under the License.
# 
# The Original Software is: "indexingTest.py"
# 
# The Initial Developer of the Original Software is:
#   Dr Peter Christen (Research School of Computer Science, The Australian
#                      National University)
# 
# Copyright (C) 2002 - 2011 the Australian National University and
# others. All Rights Reserved.
# 
# Contributors:
# 
# Alternatively, the contents of this file may be used under the terms
# of the GNU General Public License Version 2 or later (the "GPL"), in
# which case the provisions of the GPL are applicable instead of those
# above. The GPL is available at the following URL: http://www.gnu.org/
# If you wish to allow use of your version of this file only under the
# terms of the GPL, and not to allow others to use your version of this
# file under the terms of the ANUOS License, indicate your decision by
# deleting the provisions above and replace them with the notice and
# other provisions required by the GPL. If you do not delete the
# provisions above, a recipient may use your version of this file under
# the terms of any one of the ANUOS License or the GPL.
# =============================================================================
#
# Freely extensible biomedical record linkage (Febrl) - Version 0.4.2
#
# See: http://datamining.anu.edu.au/linkage.html
#
# =============================================================================

"""Test module for indexing.py.
"""

# =============================================================================
# Import necessary modules (Python standard modules first, then Febrl modules)

import sets
import sys
import unittest
sys.path.append('..')

import comparison  # Assumed to have been tested successfully
import dataset     # Assumed to have been tested successfully
import stringcmp

import indexing

def print_log(x):  # Function to be used as 'log_funct'
  print x

# =============================================================================

class TestCase(unittest.TestCase):

  # Initialise test case  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def setUp(self):

    self.dataset1 = dataset.DataSetCSV(description='First test CSV data set',
                                       access_mode='read',
                                       rec_ident='rec_id',
                                       header_line=True,
                                       file_name='./test-data.csv')

    self.dataset2 = dataset.DataSetCSV(description='Second test CSV data set',
                                       access_mode='read',
                                       rec_ident='rec_id',
                                       header_line=True,
                                       file_name='./test-data.csv')

    self.rec_ident1 = []  # Get lists of all record identifiers
    self.rec_ident2 = []
    for (rec_ident,rec) in self.dataset1.readall():
      self.rec_ident1.append(rec_ident)
    for (rec_ident,rec) in self.dataset2.readall():
      self.rec_ident2.append(rec_ident)

    gn_jfc = comparison.FieldComparatorJaro(threshold = 0.75,
                                            desc = 'Givenname Jaro')
    sn_wfcc = comparison.FieldComparatorWinkler(threshold = 0.5,
                                                desc = 'Surname Winkler',
                                         do_cache=True)
    sub_pqfcc = comparison.FieldComparatorPosQGram(threshold = 0.6,
                                                   q = 2,
                                                   max_dist = 3,
                                                   common_div = 'average',
                                                   padded = True,
                                                   desc = 'Suburb PosQGram',
                                                   do_cache=True)
    pc_kfc = comparison.FieldComparatorKeyDiff(max_key_di = 2,
                                               desc = 'Postcode KeyDiff')

    field_comp_list = [(gn_jfc,    'given_name', 'given_name'),
                       (sn_wfcc,   'surname', 'surname'),
                       (sub_pqfcc, 'suburb', 'suburb'),
                       (pc_kfc,    'postcode', 'postcode')]

    self.rec_comp_link = comparison.RecordComparator(self.dataset1,
                                                     self.dataset2,
                                                     field_comp_list,
                                          'Test record comparator for linkage')

    self.rec_comp_dedupl = comparison.RecordComparator(self.dataset1,
                                                     self.dataset1,
                                                     field_comp_list,
                                    'Test record comparator for deduplication')

  # Clean up test case  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def tearDown(self):

    self.dataset1.finalise()
    self.dataset2.finalise()

    del self.rec_comp_link
    del self.rec_comp_dedupl

  # ---------------------------------------------------------------------------
  # Start test cases

  def testFullIndexLinkage(self):  # - - - - - - - - - - - - - - - - - - - - -
    """Test FullIndex linkage"""

    return

    full_index = indexing.FullIndex(description = 'Test full index',
                                    dataset1 = self.dataset1,
                                    dataset2 = self.dataset2,
                                    rec_comparator = self.rec_comp_link,
                                    progress=2,
                                    index_sep_str = ' ',
                                    index_def = [])

    assert isinstance(full_index.index1, dict)
    assert isinstance(full_index.index2, dict)
    assert isinstance(full_index.description, str)
    assert isinstance(full_index.rec_cache1, dict)
    assert isinstance(full_index.rec_cache2, dict)
    assert full_index.progress_report == 2
    assert full_index.skip_missing == True
    assert full_index.do_deduplication == False
    assert full_index.index_def == []
    assert full_index.num_rec_pairs ==len(self.rec_ident1)*len(self.rec_ident2)

    assert full_index.status == 'initialised'

    full_index.build()  # - - - - - - - - - - - - - - - - - - - - - - - - - - -

    assert full_index.num_rec_pairs ==len(self.rec_ident1)*len(self.rec_ident2)
    assert isinstance(full_index.index1, dict)
    assert isinstance(full_index.index2, dict)
    assert isinstance(full_index.description, str)
    assert isinstance(full_index.rec_cache1, dict)
    assert isinstance(full_index.rec_cache2, dict)
    assert full_index.progress_report == 2
    assert full_index.skip_missing == True
    assert full_index.do_deduplication == False
    assert full_index.index_def == []
    assert full_index.status == 'built'
    assert len(full_index.small_data_set_dict)==full_index.dataset1.num_records

    full_index.compact()  # - - - - - - - - - - - - - - - - - - - - - - - - - -

    assert isinstance(full_index.index1, dict)
    assert isinstance(full_index.index2, dict)
    assert isinstance(full_index.description, str)
    assert isinstance(full_index.rec_cache1, dict)
    assert isinstance(full_index.rec_cache2, dict)
    assert full_index.progress_report == 2
    assert full_index.skip_missing == True
    assert full_index.do_deduplication == False
    assert full_index.index_def == []
    assert full_index.num_rec_pairs ==len(self.rec_ident1)*len(self.rec_ident2)
    assert full_index.status == 'compacted'

    w_vec_dict_list = []  # A list of the dictionaries returned from run() - -

    for i in range(10):

      [field_names_list, weight_vec_dict] = full_index.run()
      w_vec_dict_list.append(weight_vec_dict)

      assert isinstance(weight_vec_dict, dict)
      assert len(weight_vec_dict) == full_index.num_rec_pairs
      assert full_index.num_rec_pairs == len(self.rec_ident1) * \
                                         len(self.rec_ident2)

      for rec_ident1 in self.rec_ident1:
        for rec_ident2 in self.rec_ident2:
          assert (rec_ident1,rec_ident2) in weight_vec_dict
          assert rec_ident1 in self.rec_ident1
          assert rec_ident2 in self.rec_ident2
          assert isinstance(weight_vec_dict[(rec_ident1,rec_ident2)], list)
          assert len(weight_vec_dict[(rec_ident1,rec_ident2)]) == \
                 len(full_index.rec_comparator.field_comparison_list)

    for i in range(9):  # Check if all the results are the same

      assert len(w_vec_dict_list[i]) == len(w_vec_dict_list[i+1])

      rec_keys1 = w_vec_dict_list[i].keys()
      rec_keys1.sort()
      rec_keys2 = w_vec_dict_list[i+1].keys()
      rec_keys2.sort()

      rec_values1 = w_vec_dict_list[i].values()
      rec_values1.sort()
      rec_values2 = w_vec_dict_list[i+1].values()
      rec_values2.sort()

      assert rec_values1 == rec_values2

    # Test length filter (which is not used)
    #
    [field_names_list, prev_w_vec_dict] = full_index.run()

    for lf in [99,50,30,20,15,10,5,1,0.5]:

      [field_names_list, this_w_vec_dict] = \
                                        full_index.run(length_filter_perc = lf)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

    # Test cut-off threshold
    #
    [field_names_list, prev_w_vec_dict] = full_index.run()

    for cot in [0.05, 0.1, 0.2, 0.5, 0.7, 0.8, 0.9]:

      [field_names_list, this_w_vec_dict] = \
                                        full_index.run(cut_off_threshold = cot)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

  def testFullIndexDedupl(self):  # - - - - - - - - - - - - - - - - - - - - - -
    """Test FullIndex deduplication"""

    return

    full_index = indexing.FullIndex(description = 'Test full index',
                                    dataset1 = self.dataset1,
                                    dataset2 = self.dataset1,
                                    log_funct = print_log,
                                    rec_comparator = self.rec_comp_dedupl,
                                    skip_m = False,
                                    index_sep_str = '*',
                                    index_def = [])

    assert isinstance(full_index.index1, dict)
    assert isinstance(full_index.index2, dict)
    assert isinstance(full_index.description, str)
    assert isinstance(full_index.rec_cache1, dict)
    assert isinstance(full_index.rec_cache2, dict)
    assert full_index.progress_report == 10
    assert full_index.skip_missing == False
    assert full_index.do_deduplication == True
    assert full_index.index_def == []
    assert full_index.num_rec_pairs == \
           len(self.rec_ident1)*(len(self.rec_ident1)-1)/2

    assert full_index.status == 'initialised'

    full_index.build()  # - - - - - - - - - - - - - - - - - - - - - - - - - - -

    assert full_index.num_rec_pairs == \
           len(self.rec_ident1)*(len(self.rec_ident1)-1)/2
    assert isinstance(full_index.index1, dict)
    assert isinstance(full_index.index2, dict)
    assert isinstance(full_index.description, str)
    assert isinstance(full_index.rec_cache1, dict)
    assert isinstance(full_index.rec_cache2, dict)
    assert full_index.progress_report == 10
    assert full_index.skip_missing == False
    assert full_index.do_deduplication == True
    assert full_index.index_def == []
    assert full_index.status == 'built'
    assert len(full_index.small_data_set_dict)==full_index.dataset1.num_records

    full_index.compact()  # - - - - - - - - - - - - - - - - - - - - - - - - - -

    assert isinstance(full_index.index1, dict)
    assert isinstance(full_index.index2, dict)
    assert isinstance(full_index.description, str)
    assert isinstance(full_index.rec_cache1, dict)
    assert isinstance(full_index.rec_cache2, dict)
    assert full_index.progress_report == 10
    assert full_index.skip_missing == False
    assert full_index.do_deduplication == True
    assert full_index.index_def == []
    assert full_index.num_rec_pairs == \
           len(self.rec_ident1)*(len(self.rec_ident1)-1)/2

    assert full_index.status == 'compacted'

    w_vec_dict_list = []  # A list of the dictionaries returned from run() - -

    for i in range(10):

      [field_names_list, weight_vec_dict] = full_index.run()
      w_vec_dict_list.append(weight_vec_dict)

      assert isinstance(weight_vec_dict, dict)
      assert full_index.num_rec_pairs == \
             len(self.rec_ident1)*(len(self.rec_ident1)-1)/2
      assert len(weight_vec_dict) == full_index.num_rec_pairs

      tmp_rec_ident1 = self.rec_ident1[:]
      tmp_rec_ident1.sort()
      tmp_rec_ident2 = tmp_rec_ident1[:]  # Copy of list

      for rec_ident1 in tmp_rec_ident1:
        tmp_rec_ident2.remove(rec_ident1)

        for rec_ident2 in tmp_rec_ident2:
          assert (rec_ident1,rec_ident2) in weight_vec_dict
          assert rec_ident1 in self.rec_ident1
          assert rec_ident2 in self.rec_ident2
          assert isinstance(weight_vec_dict[(rec_ident1,rec_ident2)], list)
          assert len(weight_vec_dict[(rec_ident1,rec_ident2)]) == \
                 len(full_index.rec_comparator.field_comparison_list)

    for i in range(9):  # Check if all the results are the same

      assert len(w_vec_dict_list[i]) == len(w_vec_dict_list[i+1])

      rec_keys1 = w_vec_dict_list[i].keys()
      rec_keys1.sort()
      rec_keys2 = w_vec_dict_list[i+1].keys()
      rec_keys2.sort()

      rec_values1 = w_vec_dict_list[i].values()
      rec_values1.sort()
      rec_values2 = w_vec_dict_list[i+1].values()
      rec_values2.sort()

      assert rec_values1 == rec_values2

    # Test length filter
    #
    [field_names_list, prev_w_vec_dict] = full_index.run()

    for lf in [99,50,30,20,15,10,5,1,0.5]:

      [field_names_list, this_w_vec_dict] = \
                                        full_index.run(length_filter_perc = lf)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

    # Test cut-off threshold
    #
    [field_names_list, prev_w_vec_dict] = full_index.run()

    for cot in [0.05, 0.1, 0.2, 0.5, 0.7, 0.8, 0.9]:

      [field_names_list, this_w_vec_dict] = \
                                        full_index.run(cut_off_threshold = cot)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

  # ---------------------------------------------------------------------------

  def testBlockingIndexLinkage(self):  # - - - - - - - - - - - - - - - - - - -
    """Test BlockingIndex linkage"""

    index_def1 = [['surname','surname',False,False,None,[]]]
    index_def2 = [['given_name','given_name',True,True,4,[]],
                  ['postcode','postcode',True,False,2,[]]]

    # Compare with full index, assume it is correct
    #
    full_index = indexing.FullIndex(description = 'Test full index',
                                    dataset1 = self.dataset1,
                                    dataset2 = self.dataset2,
                                    rec_comparator = self.rec_comp_link,
                                    progress=2,
                                    index_def = [])
    full_index.build()
    full_index.compact()
    [field_names_list, full_weight_vec_dict] = full_index.run()

    block_index = indexing.BlockingIndex(description = 'Test blocking index',
                                         dataset1 = self.dataset1,
                                         dataset2 = self.dataset2,
                                         rec_comparator = self.rec_comp_link,
                                         progress=2,
                                         index_def = [index_def1,index_def2])

    assert isinstance(block_index.index1, dict)
    assert isinstance(block_index.index2, dict)
    assert isinstance(block_index.description, str)
    assert isinstance(block_index.rec_cache1, dict)
    assert isinstance(block_index.rec_cache2, dict)
    assert block_index.progress_report == 2
    assert block_index.skip_missing == True
    assert block_index.do_deduplication == False
    assert len(block_index.index_def) == 2
    assert block_index.num_rec_pairs == None

    assert block_index.status == 'initialised'

    block_index.build()  # - - - - - - - - - - - - - - - - - - - - - - - - - -

    assert block_index.num_rec_pairs > 0
    assert block_index.num_rec_pairs <= \
           len(self.rec_ident1)*len(self.rec_ident2)
    assert isinstance(block_index.index1, dict)
    assert isinstance(block_index.index2, dict)
    assert isinstance(block_index.description, str)
    assert isinstance(block_index.rec_cache1, dict)
    assert isinstance(block_index.rec_cache2, dict)
    assert block_index.progress_report == 2
    assert block_index.skip_missing == True
    assert block_index.do_deduplication == False
    assert len(block_index.index_def) == 2

    assert block_index.status == 'built'

    block_index.compact()  # - - - - - - - - - - - - - - - - - - - - - - - - -

    assert block_index.num_rec_pairs > 0
    assert block_index.num_rec_pairs <= \
           len(self.rec_ident1)*len(self.rec_ident2)
    assert isinstance(block_index.index1, dict)
    assert isinstance(block_index.index2, dict)
    assert isinstance(block_index.description, str)
    assert isinstance(block_index.rec_cache1, dict)
    assert isinstance(block_index.rec_cache2, dict)
    assert block_index.progress_report == 2
    assert block_index.skip_missing == True
    assert block_index.do_deduplication == False
    assert len(block_index.index_def) == 2
    assert block_index.status == 'compacted'

    [field_names_list, weight_vec_dict] = block_index.run()  # - - - - - - - -

    assert isinstance(weight_vec_dict, dict)
    assert len(weight_vec_dict) == block_index.num_rec_pairs
    assert block_index.num_rec_pairs > 0
    assert block_index.num_rec_pairs <= \
           len(self.rec_ident1)*len(self.rec_ident2)

    for rec_ident1 in self.rec_ident1:
      for rec_ident2 in self.rec_ident2:
        if (rec_ident1 == rec_ident2):
          assert (rec_ident1,rec_ident2) in weight_vec_dict

    for (rec_ident1,rec_ident2) in weight_vec_dict:
      assert rec_ident1 in self.rec_ident1
      assert rec_ident2 in self.rec_ident2
      assert isinstance(weight_vec_dict[(rec_ident1,rec_ident2)], list)

      assert len(weight_vec_dict[(rec_ident1,rec_ident2)]) == \
             len(block_index.rec_comparator.field_comparison_list)

      assert (rec_ident1,rec_ident2) in full_weight_vec_dict

    # Test length filter
    #
    [field_names_list, prev_w_vec_dict] = block_index.run()

    for lf in [99,50,30,20,15,10,5,1,0.5]:

      [field_names_list, this_w_vec_dict] = \
                                       block_index.run(length_filter_perc = lf)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

    # Test cut-off threshold
    #
    [field_names_list, prev_w_vec_dict] = block_index.run()

    for cot in [0.05, 0.1, 0.2, 0.5, 0.7, 0.8, 0.9]:

      [field_names_list, this_w_vec_dict] = \
                                       block_index.run(cut_off_threshold = cot)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

  def testBlockIndexDedupl(self):  # - - - - - - - - - - - - - - - - - - - - -
    """Test BlockingIndex deduplication"""

    index_def1 = [['surname','surname',False,False,None,[]]]
    index_def2 = [['given_name','given_name',True,True,4,[]],
                  ['postcode','postcode',True,False,2,[]]]

    # Compare with full index, assume it is correct
    #
    full_index = indexing.FullIndex(description = 'Test full index',
                                    dataset1 = self.dataset1,
                                    dataset2 = self.dataset1,
                                    rec_comparator = self.rec_comp_dedupl,
                                    progress=2,
                                    index_sep_str = '$#',
                                    index_def = [])
    full_index.build()
    full_index.compact()
    [field_names_list, full_weight_vec_dict] = full_index.run()

    block_index = indexing.BlockingIndex(description = 'Test blocking index',
                                         dataset1 = self.dataset1,
                                         dataset2 = self.dataset1,
                                         rec_comparator = self.rec_comp_dedupl,
                                         log_funct = print_log,
                                         skip_mi = False,
                                         index_sep_str = '$#',
                                         index_def = [index_def1,index_def2])

    assert isinstance(block_index.index1, dict)
    assert isinstance(block_index.index2, dict)
    assert isinstance(block_index.description, str)
    assert isinstance(block_index.rec_cache1, dict)
    assert isinstance(block_index.rec_cache2, dict)
    assert block_index.progress_report == 10
    assert block_index.skip_missing == False
    assert block_index.do_deduplication == True
    assert len(block_index.index_def) == 2
    assert block_index.num_rec_pairs == None

    assert block_index.status == 'initialised'

    block_index.build()  # - - - - - - - - - - - - - - - - - - - - - - - - - -

    assert block_index.num_rec_pairs > 0
    assert block_index.num_rec_pairs <= \
           len(self.rec_ident1)*len(self.rec_ident2)
    assert isinstance(block_index.index1, dict)
    assert isinstance(block_index.index2, dict)
    assert isinstance(block_index.description, str)
    assert isinstance(block_index.rec_cache1, dict)
    assert isinstance(block_index.rec_cache2, dict)
    assert block_index.progress_report == 10
    assert block_index.skip_missing == False
    assert block_index.do_deduplication == True
    assert len(block_index.index_def) == 2
    assert block_index.status == 'built'

    block_index.compact()  # - - - - - - - - - - - - - - - - - - - - - - - - -

    assert block_index.num_rec_pairs > 0
    assert block_index.num_rec_pairs <= \
           len(self.rec_ident1)*len(self.rec_ident2)
    assert isinstance(block_index.index1, dict)
    assert isinstance(block_index.index2, dict)
    assert isinstance(block_index.description, str)
    assert isinstance(block_index.rec_cache1, dict)
    assert isinstance(block_index.rec_cache2, dict)
    assert block_index.progress_report == 10
    assert block_index.skip_missing == False
    assert block_index.do_deduplication == True
    assert len(block_index.index_def) == 2
    assert block_index.status == 'compacted'

    [field_names_list, weight_vec_dict] = block_index.run()  # - - - - - - - -

    assert isinstance(weight_vec_dict, dict)
    assert len(weight_vec_dict) == block_index.num_rec_pairs
    assert block_index.num_rec_pairs > 0
    assert block_index.num_rec_pairs <= \
           len(self.rec_ident1)*len(self.rec_ident2)

    for (rec_ident1,rec_ident2) in weight_vec_dict:
      assert rec_ident1 in self.rec_ident1
      assert rec_ident2 in self.rec_ident1
      assert isinstance(weight_vec_dict[(rec_ident1,rec_ident2)], list)

      assert len(weight_vec_dict[(rec_ident1,rec_ident2)]) == \
             len(block_index.rec_comparator.field_comparison_list)

      assert (rec_ident1,rec_ident2) in full_weight_vec_dict

    # Test length filter
    #
    [field_names_list, prev_w_vec_dict] = block_index.run()

    for lf in [99,50,30,20,15,10,5,1,0.5]:

      [field_names_list, this_w_vec_dict] = \
                                       block_index.run(length_filter_perc = lf)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

    # Test cut-off threshold
    #
    [field_names_list, prev_w_vec_dict] = block_index.run()

    for cot in [0.05, 0.1, 0.2, 0.5, 0.7, 0.8, 0.9]:

      [field_names_list, this_w_vec_dict] = \
                                       block_index.run(cut_off_threshold = cot)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

  # ---------------------------------------------------------------------------

  def testSortingIndexLinkage(self):  # - - - - - - - - - - - - - - - - - - - -
    """Test SortingIndex linkage"""

    return

    index_def1 = [['surname','surname',False,False,None,[]]]
    index_def2 = [['given_name','given_name',True,True,4,[]],
                  ['postcode','postcode',True,False,2,[]]]

    # Assume blocking index is OK, so get the record pairs from this index for
    # comparison
    #
    block_index = indexing.BlockingIndex(description = 'Test blocking index',
                                         dataset1 = self.dataset1,
                                         dataset2 = self.dataset2,
                                         rec_comparator = self.rec_comp_link,
                                         progress=2,
                                         index_def = [index_def1,index_def2])
    block_index.build()
    block_index.compact()
    [field_names_list, block_index_weight_vec_dict] = block_index.run()

    # Make sure record identifiers in each pair are sorted (smaller id first)
    #
    for ((rec_ident1,rec_ident2),w) in block_index_weight_vec_dict.items():
      if (rec_ident1 > rec_ident2):
        del block_index_weight_vec_dict[(rec_ident1,rec_ident2)]
        block_index_weight_vec_dict[(rec_ident2,rec_ident1)] = w

    prev_rec_dict = {}

    for w in [1,2,3,4,5,6,7,9,11,12,15,20]:

      sort_index = indexing.SortingIndex(description = 'Test Sorting index',
                                         dataset1 = self.dataset1,
                                         dataset2 = self.dataset2,
                                         rec_comparator = self.rec_comp_link,
                                         progress=2,
                                         window_s = w,
                                         index_def = [index_def1,index_def2])

      assert isinstance(sort_index.index1, dict)
      assert isinstance(sort_index.index2, dict)
      assert isinstance(sort_index.description, str)
      assert isinstance(sort_index.rec_cache1, dict)
      assert isinstance(sort_index.rec_cache2, dict)
      assert sort_index.progress_report == 2
      assert sort_index.skip_missing == True
      assert sort_index.do_deduplication == False
      assert len(sort_index.index_def) == 2
      assert sort_index.window_size == w
      assert sort_index.num_rec_pairs == None

      assert sort_index.status == 'initialised'

      sort_index.build()  # - - - - - - - - - - - - - - - - - - - - - - - - - -

      assert sort_index.num_rec_pairs == None
      assert isinstance(sort_index.index1, dict)
      assert isinstance(sort_index.index2, dict)
      assert isinstance(sort_index.description, str)
      assert isinstance(sort_index.rec_cache1, dict)
      assert isinstance(sort_index.rec_cache2, dict)
      assert sort_index.progress_report == 2
      assert sort_index.skip_missing == True
      assert sort_index.do_deduplication == False
      assert len(sort_index.index_def) == 2
      assert sort_index.window_size == w

      assert sort_index.status == 'built'

      sort_index.compact()  # - - - - - - - - - - - - - - - - - - - - - - - - -

      assert sort_index.num_rec_pairs > 0
      assert sort_index.num_rec_pairs <= \
             len(self.rec_ident1)*len(self.rec_ident2)
      assert isinstance(sort_index.index1, dict)
      assert isinstance(sort_index.index2, dict)
      assert isinstance(sort_index.description, str)
      assert isinstance(sort_index.rec_cache1, dict)
      assert isinstance(sort_index.rec_cache2, dict)
      assert sort_index.progress_report == 2
      assert sort_index.skip_missing == True
      assert sort_index.do_deduplication == False
      assert len(sort_index.index_def) == 2
      assert sort_index.window_size == w

      assert sort_index.status == 'compacted'

      [field_names_list, weight_vec_dict] = sort_index.run()  # - - - - - - - -

      assert isinstance(weight_vec_dict, dict)
      assert len(weight_vec_dict) == sort_index.num_rec_pairs
      assert sort_index.num_rec_pairs > 0
      assert sort_index.num_rec_pairs <= \
             len(self.rec_ident1)*len(self.rec_ident2)

      for rec_ident1 in self.rec_ident1:
        for rec_ident2 in self.rec_ident2:
          if (rec_ident1 == rec_ident2):
            assert (rec_ident1,rec_ident2) in weight_vec_dict, \
                   (rec_ident1,rec_ident2)

      for (rec_ident1,rec_ident2) in weight_vec_dict:
        assert rec_ident1 in self.rec_ident1
        assert rec_ident2 in self.rec_ident2
        assert isinstance(weight_vec_dict[(rec_ident1,rec_ident2)], list)

        assert len(weight_vec_dict[(rec_ident1,rec_ident2)]) == \
               len(sort_index.rec_comparator.field_comparison_list)

      # All record pairs from the blocking index should be in the sorting index
      # as well
      #
      assert len(block_index_weight_vec_dict) <= len(weight_vec_dict), \
             ('Number of record pairs in blocking index: %d, and in ' % \
             (len(block_index_weight_vec_dict))+' sorting index: %d' % \
             (len(weight_vec_dict)))

      # Build a sorting weight vector dict with pair identifiers sorted
      #
      for ((rec_ident1,rec_ident2), weight_list) in weight_vec_dict.items():
        if (rec_ident1 > rec_ident2):
          del weight_vec_dict[(rec_ident1,rec_ident2)]
          weight_vec_dict[(rec_ident2,rec_ident1)] = weight_list

      for (rec_ident1,rec_ident2) in block_index_weight_vec_dict:
        assert (rec_ident1,rec_ident2) in weight_vec_dict, \
               'Record pair %s from blocking index not in sorting index' % \
               (str((rec_ident1,rec_ident2)))

      for rec_ident_pair in prev_rec_dict.iterkeys():
        assert rec_ident_pair in weight_vec_dict, \
               'Record pair %s not in weight vector dict with window = %d' % \
               (str(rec_ident_pair), w)

      prev_rec_dict = weight_vec_dict

    # Sorting index with window size 1 should be same as blocking index - - - -
    #
    sort_index = indexing.SortingIndex(description = 'Test Sorting index',
                                       dataset1 = self.dataset1,
                                       dataset2 = self.dataset2,
                                       rec_comparator = self.rec_comp_link,
                                       log_funct = print_log,
                                       progress=2,
                                       window_s = 1,
                                       index_def = [index_def1,index_def2])
    sort_index.build()
    sort_index.compact()
    [field_names_list, weight_vec_dict] = sort_index.run()

    # Make sure record identifiers in each pair are sorted (smaller id first)
    #
    for ((rec_ident1,rec_ident2),w) in weight_vec_dict.items():
      if (rec_ident1 > rec_ident2):
        del weight_vec_dict[(rec_ident1,rec_ident2)]
        weight_vec_dict[(rec_ident2,rec_ident1)] = w

    for (rec_ident1,rec_ident2) in weight_vec_dict.keys():
      assert (rec_ident1,rec_ident2) in block_index_weight_vec_dict, \
             'Not in blocking weight vectors: %s' % \
             (str((rec_ident1,rec_ident2)))

    for (rec_ident1,rec_ident2) in block_index_weight_vec_dict.keys():
      assert (rec_ident1,rec_ident2) in weight_vec_dict, \
             'Not in sorting weight vectors: %s' % \
             (str((rec_ident1,rec_ident2)))

    # Test length filter
    #
    [field_names_list, prev_w_vec_dict] = sort_index.run()

    for lf in [99,50,30,20,15,10,5,1,0.5]:

      [field_names_list, this_w_vec_dict] = \
                                        sort_index.run(length_filter_perc = lf)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

    # Test cut-off threshold
    #
    [field_names_list, prev_w_vec_dict] = sort_index.run()

    for cot in [0.05, 0.1, 0.2, 0.5, 0.7, 0.8, 0.9]:

      [field_names_list, this_w_vec_dict] = \
                                        sort_index.run(cut_off_threshold = cot)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

  def testSortingIndexDedupl(self):  # - - - - - - - - - - - - - - - - - - - -
    """Test SortingIndex deduplication"""

    return

    index_def1 = [['surname','surname',False,False,None,[]]]
    index_def2 = [['given_name','given_name',True,True,4,[]],
                  ['postcode','postcode',True,False,2,[]]]

    # Assume blocking index is OK, so get the record pairs from this index for
    # comparison
    #
    block_index = indexing.BlockingIndex(description = 'Test blocking index',
                                         dataset1 = self.dataset1,
                                         dataset2 = self.dataset1,
                                         rec_comparator = self.rec_comp_dedupl,
                                         progress=4,
                                         index_def = [index_def1,index_def2])
    block_index.build()
    block_index.compact()
    [field_names_list, block_index_weight_vec_dict] = block_index.run()

    # Make sure record identifiers in each pair are sorted (smaller id first)
    #
    for ((rec_ident1,rec_ident2),w) in block_index_weight_vec_dict.items():
      if (rec_ident1 > rec_ident2):
        del block_index_weight_vec_dict[(rec_ident1,rec_ident2)]
        block_index_weight_vec_dict[(rec_ident2,rec_ident1)] = w

    prev_rec_dict = {}

    for w in [1, 2, 3, 4, 5, 6, 7, 9, 11, 12, 15, 20]:

      sort_index = indexing.SortingIndex(description = 'Test Sorting index',
                                         dataset1 = self.dataset1,
                                         dataset2 = self.dataset1,
                                         rec_comparator = self.rec_comp_dedupl,
                                         progress=4,
                                         window_s = w,
                                         skip_m = False,
                                         index_def = [index_def1,index_def2])

      assert isinstance(sort_index.index1, dict)
      assert isinstance(sort_index.index2, dict)
      assert isinstance(sort_index.description, str)
      assert isinstance(sort_index.rec_cache1, dict)
      assert isinstance(sort_index.rec_cache2, dict)
      assert sort_index.progress_report == 4
      assert sort_index.skip_missing == False
      assert sort_index.do_deduplication == True
      assert len(sort_index.index_def) == 2
      assert sort_index.window_size == w
      assert sort_index.num_rec_pairs == None

      assert sort_index.status == 'initialised'

      sort_index.build()  # - - - - - - - - - - - - - - - - - - - - - - - - - -

      assert sort_index.num_rec_pairs == None
      assert isinstance(sort_index.index1, dict)
      assert isinstance(sort_index.index2, dict)
      assert isinstance(sort_index.description, str)
      assert isinstance(sort_index.rec_cache1, dict)
      assert isinstance(sort_index.rec_cache2, dict)
      assert sort_index.progress_report == 4
      assert sort_index.skip_missing == False
      assert sort_index.do_deduplication == True
      assert len(sort_index.index_def) == 2
      assert sort_index.window_size == w

      assert sort_index.status == 'built'

      sort_index.compact()  # - - - - - - - - - - - - - - - - - - - - - - - - -

      assert sort_index.num_rec_pairs > 0
      assert sort_index.num_rec_pairs <= \
             len(self.rec_ident1)*len(self.rec_ident2)
      assert isinstance(sort_index.index1, dict)
      assert isinstance(sort_index.index2, dict)
      assert isinstance(sort_index.description, str)
      assert isinstance(sort_index.rec_cache1, dict)
      assert isinstance(sort_index.rec_cache2, dict)
      assert sort_index.progress_report == 4
      assert sort_index.skip_missing == False
      assert sort_index.do_deduplication == True
      assert len(sort_index.index_def) == 2
      assert sort_index.window_size == w

      assert sort_index.status == 'compacted'

      [field_names_list, weight_vec_dict] = sort_index.run()  # - - - - - - - -

      assert isinstance(weight_vec_dict, dict)
      assert len(weight_vec_dict) == sort_index.num_rec_pairs
      assert sort_index.num_rec_pairs > 0
      assert sort_index.num_rec_pairs <= \
             len(self.rec_ident1)*len(self.rec_ident2)

      for rec_ident1 in self.rec_ident1:
        assert (rec_ident1, rec_ident1) not in weight_vec_dict, \
               (rec_ident1, rec_ident1)

      for (rec_ident1,rec_ident2) in weight_vec_dict:
        assert rec_ident1 in self.rec_ident1
        assert rec_ident2 in self.rec_ident1
        assert isinstance(weight_vec_dict[(rec_ident1,rec_ident2)], list)

        assert len(weight_vec_dict[(rec_ident1,rec_ident2)]) == \
               len(sort_index.rec_comparator.field_comparison_list)

      # All record pairs from the blocking index should be in the sorting index
      # as well
      #
      assert len(block_index_weight_vec_dict) <= len(weight_vec_dict), \
             ('Number of record pairs in blocking index: %d, and in ' % \
             (len(block_index_weight_vec_dict))+' sorting index: %d' % \
             (len(weight_vec_dict)))

      # Build a sorting weight vector dict with pair identifiers sorted
      #
      for ((rec_ident1,rec_ident2), weight_list) in weight_vec_dict.items():
        if (rec_ident1 > rec_ident2):
          del weight_vec_dict[(rec_ident1,rec_ident2)]
          weight_vec_dict[(rec_ident2,rec_ident1)] = weight_list

      for (rec_ident1,rec_ident2) in block_index_weight_vec_dict:
        assert (rec_ident1,rec_ident2) in weight_vec_dict, \
               'Record pair %s from blocking index not in sorting index' % \
               (str((rec_ident1,rec_ident2)))

      for rec_ident_pair in prev_rec_dict.iterkeys():
        assert rec_ident_pair in weight_vec_dict, \
               'Record pair %s not in weight vector dict with window = %d' % \
               (str(rec_ident_pair), w)

      prev_rec_dict = weight_vec_dict

    # Sorting index with window size 1 should be same as blocking index - - - -
    #
    sort_index = indexing.SortingIndex(description = 'Test Sorting index',
                                       dataset1 = self.dataset1,
                                       dataset2 = self.dataset1,
                                       rec_comparator = self.rec_comp_dedupl,
                                       progress=2,
                                       window_s = 1,
                                       index_def = [index_def1,index_def2])
    sort_index.build()
    sort_index.compact()
    [field_names_list, weight_vec_dict] = sort_index.run()

    # Make sure record identifiers in each pair are sorted (smaller id first)
    #
    for ((rec_ident1,rec_ident2),w) in weight_vec_dict.items():
      if (rec_ident1 > rec_ident2):
        del weight_vec_dict[(rec_ident1,rec_ident2)]
        weight_vec_dict[(rec_ident2,rec_ident1)] = w

    for (rec_ident1,rec_ident2) in weight_vec_dict.keys():
      assert (rec_ident1,rec_ident2) in block_index_weight_vec_dict, \
             'Not in blocking weight vectors: %s' % \
             (str((rec_ident1,rec_ident2)))

    for (rec_ident1,rec_ident2) in block_index_weight_vec_dict.keys():
      assert (rec_ident1,rec_ident2) in weight_vec_dict, \
             'Not in sorting weight vectors: %s' % \
             (str((rec_ident1,rec_ident2)))

    # Test length filter
    #
    [field_names_list, prev_w_vec_dict] = sort_index.run()

    for lf in [99,50,30,20,15,10,5,1,0.5]:

      [field_names_list, this_w_vec_dict] = \
                                        sort_index.run(length_filter_perc = lf)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

    # Test cut-off threshold
    #
    [field_names_list, prev_w_vec_dict] = sort_index.run()

    for cot in [0.05, 0.1, 0.2, 0.5, 0.7, 0.8, 0.9]:

      [field_names_list, this_w_vec_dict] = \
                                        sort_index.run(cut_off_threshold = cot)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict


  # ---------------------------------------------------------------------------

  def testAdaptSortIndexLinkage(self):  # - - - - - - - - - - - - - - - - - - -
    """Test AdaptSortingIndex linkage"""

    #return

    index_def1 = [['surname','surname',False,False,None,[]]]
    index_def2 = [['given_name','given_name',True,True,4,[]],
                  ['postcode','postcode',True,False,2,[]]]

    # Assume blocking index is OK, so get the record pairs from this index for
    # comparison
    #
    block_index = indexing.BlockingIndex(description = 'Test blocking index',
                                         dataset1 = self.dataset1,
                                         dataset2 = self.dataset2,
                                         rec_comparator = self.rec_comp_link,
                                         progress=2,
                                         index_def = [index_def1,index_def2])
    block_index.build()
    block_index.compact()
    [field_names_list, block_index_weight_vec_dict] = block_index.run()

    # Make sure record identifiers in each pair are sorted (smaller id first)
    #
    for ((rec_ident1,rec_ident2),w) in block_index_weight_vec_dict.items():
      if (rec_ident1 > rec_ident2):
        del block_index_weight_vec_dict[(rec_ident1,rec_ident2)]
        block_index_weight_vec_dict[(rec_ident2,rec_ident1)] = w

    prev_rec_dict = {}

    for (str_cmp_funct, str_cmp_thres) in \
         [(stringcmp.jaro,0.9), (stringcmp.jaro,0.8), (stringcmp.jaro,0.7),
          (stringcmp.bigram,0.7),(stringcmp.bigram,0.8),
          (stringcmp.bigram,0.9),
          (stringcmp.lcs, 0.7), (stringcmp.lcs, 0.8), (stringcmp.lcs, 0.9)]:

      adsort_index = indexing.AdaptSortingIndex(description = \
                                              'Test Adapt Sorting index',
                                         dataset1 = self.dataset1,
                                         dataset2 = self.dataset2,
                                         rec_comparator = self.rec_comp_link,
                                         progress=2,
                                         str_cmp_funct = str_cmp_funct,
                                         str_cmp_thres = str_cmp_thres,
                                         index_def = [index_def1,index_def2])


      assert isinstance(adsort_index.index1, dict)
      assert isinstance(adsort_index.index2, dict)
      assert isinstance(adsort_index.description, str)
      assert isinstance(adsort_index.rec_cache1, dict)
      assert isinstance(adsort_index.rec_cache2, dict)
      assert adsort_index.progress_report == 2
      assert adsort_index.skip_missing == True
      assert adsort_index.do_deduplication == False
      assert len(adsort_index.index_def) == 2
      assert adsort_index.str_cmp_thres >= 0.0
      assert adsort_index.str_cmp_thres <= 1.0
      assert adsort_index.num_rec_pairs == None

      assert adsort_index.status == 'initialised'

      adsort_index.build()  # - - - - - - - - - - - - - - - - - - - - - - - - -

      assert adsort_index.num_rec_pairs == None
      assert isinstance(adsort_index.index1, dict)
      assert isinstance(adsort_index.index2, dict)
      assert isinstance(adsort_index.description, str)
      assert isinstance(adsort_index.rec_cache1, dict)
      assert isinstance(adsort_index.rec_cache2, dict)
      assert adsort_index.progress_report == 2
      assert adsort_index.skip_missing == True
      assert adsort_index.do_deduplication == False
      assert len(adsort_index.index_def) == 2
      assert adsort_index.str_cmp_thres >= 0.0
      assert adsort_index.str_cmp_thres <= 1.0

      assert adsort_index.status == 'built'

      adsort_index.compact()  # - - - - - - - - - - - - - - - - - - - - - - - -

      assert adsort_index.num_rec_pairs > 0
      assert adsort_index.num_rec_pairs <= \
             len(self.rec_ident1)*len(self.rec_ident2)
      assert isinstance(adsort_index.index1, dict)
      assert isinstance(adsort_index.index2, dict)
      assert isinstance(adsort_index.description, str)
      assert isinstance(adsort_index.rec_cache1, dict)
      assert isinstance(adsort_index.rec_cache2, dict)
      assert adsort_index.progress_report == 2
      assert adsort_index.skip_missing == True
      assert adsort_index.do_deduplication == False
      assert len(adsort_index.index_def) == 2
      assert adsort_index.str_cmp_thres >= 0.0
      assert adsort_index.str_cmp_thres <= 1.0

      assert adsort_index.status == 'compacted'

      [field_names_list, weight_vec_dict] = adsort_index.run()  # - - - - - - -

      assert isinstance(weight_vec_dict, dict)
      assert len(weight_vec_dict) == adsort_index.num_rec_pairs
      assert adsort_index.num_rec_pairs > 0
      assert adsort_index.num_rec_pairs <= \
             len(self.rec_ident1)*len(self.rec_ident2)

      for rec_ident1 in self.rec_ident1:
        for rec_ident2 in self.rec_ident2:
          if (rec_ident1 == rec_ident2):
            assert (rec_ident1,rec_ident2) in weight_vec_dict, \
                   (rec_ident1,rec_ident2)

      for (rec_ident1,rec_ident2) in weight_vec_dict:
        assert rec_ident1 in self.rec_ident1
        assert rec_ident2 in self.rec_ident2
        assert isinstance(weight_vec_dict[(rec_ident1,rec_ident2)], list)

        assert len(weight_vec_dict[(rec_ident1,rec_ident2)]) == \
               len(adsort_index.rec_comparator.field_comparison_list)

      # All record pairs from the blocking index should be in the sorting index
      # as well
      #
      assert len(block_index_weight_vec_dict) <= len(weight_vec_dict), \
             ('Number of record pairs in blocking index: %d, and in ' % \
             (len(block_index_weight_vec_dict))+' adapt sorting index: %d' % \
             (len(weight_vec_dict)))

      # Build a sorting weight vector dict with pair identifiers sorted
      #
      for ((rec_ident1,rec_ident2), weight_list) in weight_vec_dict.items():
        if (rec_ident1 > rec_ident2):
          del weight_vec_dict[(rec_ident1,rec_ident2)]
          weight_vec_dict[(rec_ident2,rec_ident1)] = weight_list

      for (rec_ident1,rec_ident2) in block_index_weight_vec_dict:
        assert (rec_ident1,rec_ident2) in weight_vec_dict, \
               'Record pair %s from blocking index not in sorting index' % \
               (str((rec_ident1,rec_ident2)))

      prev_rec_dict = weight_vec_dict

      # Make sure record identifiers in each pair are sorted (smaller id first)
      #
      for ((rec_ident1,rec_ident2),w) in weight_vec_dict.items():
        if (rec_ident1 > rec_ident2):
          del weight_vec_dict[(rec_ident1,rec_ident2)]
          weight_vec_dict[(rec_ident2,rec_ident1)] = w

      #for (rec_ident1,rec_ident2) in weight_vec_dict.keys():
      #  assert (rec_ident1,rec_ident2) in block_index_weight_vec_dict, \
      #         'Not in blocking weight vectors: %s' % \
      #         (str((rec_ident1,rec_ident2)))

      for (rec_ident1,rec_ident2) in block_index_weight_vec_dict.keys():
        assert (rec_ident1,rec_ident2) in weight_vec_dict, \
               'Not in sorting weight vectors: %s' % \
               (str((rec_ident1,rec_ident2)))

    # Test length filter
    #
    [field_names_list, prev_w_vec_dict] = adsort_index.run()

    for lf in [99,50,30,20,15,10,5,1,0.5]:

      [field_names_list, this_w_vec_dict] = \
                                      adsort_index.run(length_filter_perc = lf)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

    # Test cut-off threshold
    #
    [field_names_list, prev_w_vec_dict] = adsort_index.run()

    for cot in [0.05, 0.1, 0.2, 0.5, 0.7, 0.8, 0.9]:

      [field_names_list, this_w_vec_dict] = \
                                      adsort_index.run(cut_off_threshold = cot)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

  def testAdaptSortIndexDedupl(self):  # - - - - - - - - - - - - - - - - - - -
    """Test AdaptSortingIndex deduplication"""

    #return

    index_def1 = [['surname','surname',False,False,None,[]]]
    index_def2 = [['given_name','given_name',True,True,4,[]],
                  ['postcode','postcode',True,False,2,[]]]

    # Assume blocking index is OK, so get the record pairs from this index for
    # comparison
    #
    block_index = indexing.BlockingIndex(description = 'Test blocking index',
                                         dataset1 = self.dataset1,
                                         dataset2 = self.dataset1,
                                         rec_comparator = self.rec_comp_dedupl,
                                         progress=4,
                                         index_def = [index_def1,index_def2])
    block_index.build()
    block_index.compact()
    [field_names_list, block_index_weight_vec_dict] = block_index.run()

    # Make sure record identifiers in each pair are sorted (smaller id first)
    #
    for ((rec_ident1,rec_ident2),w) in block_index_weight_vec_dict.items():
      if (rec_ident1 > rec_ident2):
        del block_index_weight_vec_dict[(rec_ident1,rec_ident2)]
        block_index_weight_vec_dict[(rec_ident2,rec_ident1)] = w

    prev_rec_dict = {}

    for (str_cmp_funct, str_cmp_thres) in \
         [(stringcmp.jaro,0.9), (stringcmp.jaro,0.8), (stringcmp.jaro,0.7),
          (stringcmp.bigram,0.7),(stringcmp.bigram,0.8),
          (stringcmp.bigram,0.9),
          (stringcmp.lcs, 0.7), (stringcmp.lcs, 0.8), (stringcmp.lcs, 0.9)]:

      adsort_index = indexing.AdaptSortingIndex(description = \
                                              'Test Adapt Sorting index',
                                         dataset1 = self.dataset1,
                                         dataset2 = self.dataset1,
                                         rec_comparator = self.rec_comp_dedupl,
                                         progress=4,
                                         str_cmp_funct = str_cmp_funct,
                                         str_cmp_thres = str_cmp_thres,
                                         skip_m = False,
                                         index_def = [index_def1,index_def2])

      assert isinstance(adsort_index.index1, dict)
      assert isinstance(adsort_index.index2, dict)
      assert isinstance(adsort_index.description, str)
      assert isinstance(adsort_index.rec_cache1, dict)
      assert isinstance(adsort_index.rec_cache2, dict)
      assert adsort_index.progress_report == 4
      assert adsort_index.skip_missing == False
      assert adsort_index.do_deduplication == True
      assert adsort_index.str_cmp_thres >= 0.0
      assert adsort_index.str_cmp_thres <= 1.0
      assert len(adsort_index.index_def) == 2
      assert adsort_index.num_rec_pairs == None

      assert adsort_index.status == 'initialised'

      adsort_index.build()  # - - - - - - - - - - - - - - - - - - - - - - - - -

      assert adsort_index.num_rec_pairs == None
      assert isinstance(adsort_index.index1, dict)
      assert isinstance(adsort_index.index2, dict)
      assert isinstance(adsort_index.description, str)
      assert isinstance(adsort_index.rec_cache1, dict)
      assert isinstance(adsort_index.rec_cache2, dict)
      assert adsort_index.progress_report == 4
      assert adsort_index.skip_missing == False
      assert adsort_index.do_deduplication == True
      assert len(adsort_index.index_def) == 2
      assert adsort_index.str_cmp_thres >= 0.0
      assert adsort_index.str_cmp_thres <= 1.0

      assert adsort_index.status == 'built'

      adsort_index.compact()  # - - - - - - - - - - - - - - - - - - - - - - - -

      assert adsort_index.num_rec_pairs > 0
      assert adsort_index.num_rec_pairs <= \
             len(self.rec_ident1)*len(self.rec_ident2)
      assert isinstance(adsort_index.index1, dict)
      assert isinstance(adsort_index.index2, dict)
      assert isinstance(adsort_index.description, str)
      assert isinstance(adsort_index.rec_cache1, dict)
      assert isinstance(adsort_index.rec_cache2, dict)
      assert adsort_index.progress_report == 4
      assert adsort_index.skip_missing == False
      assert adsort_index.do_deduplication == True
      assert len(adsort_index.index_def) == 2
      assert adsort_index.str_cmp_thres >= 0.0
      assert adsort_index.str_cmp_thres <= 1.0

      assert adsort_index.status == 'compacted'

      [field_names_list, weight_vec_dict] = adsort_index.run()  # - - - - - - -

      assert isinstance(weight_vec_dict, dict)
      assert len(weight_vec_dict) == adsort_index.num_rec_pairs
      assert adsort_index.num_rec_pairs > 0
      assert adsort_index.num_rec_pairs <= \
             len(self.rec_ident1)*len(self.rec_ident2)

      for rec_ident1 in self.rec_ident1:
        assert (rec_ident1, rec_ident1) not in weight_vec_dict, \
               (rec_ident1, rec_ident1)

      for (rec_ident1,rec_ident2) in weight_vec_dict:
        assert rec_ident1 in self.rec_ident1
        assert rec_ident2 in self.rec_ident1
        assert isinstance(weight_vec_dict[(rec_ident1,rec_ident2)], list)

        assert len(weight_vec_dict[(rec_ident1,rec_ident2)]) == \
               len(adsort_index.rec_comparator.field_comparison_list)

      # All record pairs from the blocking index should be in the sorting index
      # as well
      #
      assert len(block_index_weight_vec_dict) <= len(weight_vec_dict), \
             ('Number of record pairs in blocking index: %d, and in ' % \
             (len(block_index_weight_vec_dict))+' adapt sorting index: %d' % \
             (len(weight_vec_dict)))

      # Build a sorting weight vector dict with pair identifiers sorted
      #
      for ((rec_ident1,rec_ident2), weight_list) in weight_vec_dict.items():
        if (rec_ident1 > rec_ident2):
          del weight_vec_dict[(rec_ident1,rec_ident2)]
          weight_vec_dict[(rec_ident2,rec_ident1)] = weight_list

      for (rec_ident1,rec_ident2) in block_index_weight_vec_dict:
        assert (rec_ident1,rec_ident2) in weight_vec_dict, \
               'Record pair %s from blocking index not in adapt sort index' % \
               (str((rec_ident1,rec_ident2)))

      prev_rec_dict = weight_vec_dict

      # Make sure record identifiers in each pair are sorted (smaller id first)
      #
      for ((rec_ident1,rec_ident2),w) in weight_vec_dict.items():
        if (rec_ident1 > rec_ident2):
          del weight_vec_dict[(rec_ident1,rec_ident2)]
          weight_vec_dict[(rec_ident2,rec_ident1)] = w

      #for (rec_ident1,rec_ident2) in weight_vec_dict.keys():
      #  assert (rec_ident1,rec_ident2) in block_index_weight_vec_dict, \
      #         'Not in blocking weight vectors: %s' % \
      #         (str((rec_ident1,rec_ident2)))

      for (rec_ident1,rec_ident2) in block_index_weight_vec_dict.keys():
        assert (rec_ident1,rec_ident2) in weight_vec_dict, \
               'Not in sorting weight vectors: %s' % \
               (str((rec_ident1,rec_ident2)))

    # Test length filter
    #
    [field_names_list, prev_w_vec_dict] = adsort_index.run()

    for lf in [99,50,30,20,15,10,5,1,0.5]:

      [field_names_list, this_w_vec_dict] = \
                                      adsort_index.run(length_filter_perc = lf)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

    # Test cut-off threshold
    #
    [field_names_list, prev_w_vec_dict] = adsort_index.run()

    for cot in [0.05, 0.1, 0.2, 0.5, 0.7, 0.8, 0.9]:

      [field_names_list, this_w_vec_dict] = \
                                      adsort_index.run(cut_off_threshold = cot)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

  # ---------------------------------------------------------------------------

  def testSortArrayIndexLinkage(self):  # - - - - - - - - - - - - - - - - - - -
    """Test SortingArrayIndex linkage"""

    return

    index_def1 = [['surname','surname',False,False,None,[]]]
    index_def2 = [['given_name','given_name',True,True,4,[]],
                  ['postcode','postcode',True,False,2,[]]]

    prev_rec_dict = {}

    for w in [2,3,4,5,6,7,9,11,12,15,20]:

      sort_index = indexing.SortingArrayIndex(description = \
                                              'Test SortingArray index',
                                         dataset1 = self.dataset1,
                                         dataset2 = self.dataset2,
                                         rec_comparator = self.rec_comp_link,
                                         progress=2,
                                         window_s = w,
                                         index_def = [index_def1,index_def2])

      assert isinstance(sort_index.index1, dict)
      assert isinstance(sort_index.index2, dict)
      assert isinstance(sort_index.description, str)
      assert isinstance(sort_index.rec_cache1, dict)
      assert isinstance(sort_index.rec_cache2, dict)
      assert sort_index.progress_report == 2
      assert sort_index.skip_missing == True
      assert sort_index.do_deduplication == False
      assert len(sort_index.index_def) == 2
      assert sort_index.window_size == w
      assert sort_index.num_rec_pairs == None

      assert sort_index.status == 'initialised'

      sort_index.build()  # - - - - - - - - - - - - - - - - - - - - - - - - - -

      assert sort_index.num_rec_pairs == None
      assert isinstance(sort_index.index1, dict)
      assert isinstance(sort_index.index2, dict)
      assert isinstance(sort_index.description, str)
      assert isinstance(sort_index.rec_cache1, dict)
      assert isinstance(sort_index.rec_cache2, dict)
      assert sort_index.progress_report == 2
      assert sort_index.skip_missing == True
      assert sort_index.do_deduplication == False
      assert len(sort_index.index_def) == 2
      assert sort_index.window_size == w

      assert sort_index.status == 'built'

      sort_index.compact()  # - - - - - - - - - - - - - - - - - - - - - - - - -

      assert sort_index.num_rec_pairs > 0
      assert sort_index.num_rec_pairs <= \
             len(self.rec_ident1)*len(self.rec_ident2)
      assert isinstance(sort_index.index1, dict)
      assert isinstance(sort_index.index2, dict)
      assert isinstance(sort_index.description, str)
      assert isinstance(sort_index.rec_cache1, dict)
      assert isinstance(sort_index.rec_cache2, dict)
      assert sort_index.progress_report == 2
      assert sort_index.skip_missing == True
      assert sort_index.do_deduplication == False
      assert len(sort_index.index_def) == 2
      assert sort_index.window_size == w

      assert sort_index.status == 'compacted'

      [field_names_list, weight_vec_dict] = sort_index.run()  # - - - - - - - -

      assert isinstance(weight_vec_dict, dict)
      assert len(weight_vec_dict) == sort_index.num_rec_pairs
      assert sort_index.num_rec_pairs > 0
      assert sort_index.num_rec_pairs <= \
             len(self.rec_ident1)*len(self.rec_ident2)

#      for rec_ident1 in self.rec_ident1:
#        for rec_ident2 in self.rec_ident2:
#          if (rec_ident1 == rec_ident2):
#            assert (rec_ident1,rec_ident2) in weight_vec_dict, \
#                   (rec_ident1,rec_ident2)

      for (rec_ident1,rec_ident2) in weight_vec_dict:
        assert rec_ident1 in self.rec_ident1
        assert rec_ident2 in self.rec_ident2
        assert isinstance(weight_vec_dict[(rec_ident1,rec_ident2)], list)

        assert len(weight_vec_dict[(rec_ident1,rec_ident2)]) == \
               len(sort_index.rec_comparator.field_comparison_list)

      # Build a sorting weight vector dict with pair identifiers sorted
      #
      for ((rec_ident1,rec_ident2), weight_list) in weight_vec_dict.items():
        if (rec_ident1 > rec_ident2):
          del weight_vec_dict[(rec_ident1,rec_ident2)]
          weight_vec_dict[(rec_ident2,rec_ident1)] = weight_list

      if (prev_rec_dict != {}):  # Only if there is a previous one
        for rec_ident_pair in prev_rec_dict.iterkeys():
          assert rec_ident_pair in weight_vec_dict, \
                 'Record pair %s not in weight vector dict with window = %d' % \
                 (str(rec_ident_pair), w)

      prev_rec_dict = weight_vec_dict

    # Test length filter
    #
    [field_names_list, prev_w_vec_dict] = sort_index.run()

    for lf in [99,50,30,20,15,10,5,1,0.5]:

      [field_names_list, this_w_vec_dict] = \
                                        sort_index.run(length_filter_perc = lf)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

    # Test cut-off threshold
    #
    [field_names_list, prev_w_vec_dict] = sort_index.run()

    for cot in [0.05, 0.1, 0.2, 0.5, 0.7, 0.8, 0.9]:

      [field_names_list, this_w_vec_dict] = \
                                        sort_index.run(cut_off_threshold = cot)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict


  def testSortArrayIndexDedupl(self):  # - - - - - - - - - - - - - - - - - - - -
    """Test SortingArrayIndex deduplication"""

    return

    index_def1 = [['surname','surname',False,False,None,[]]]
    index_def2 = [['given_name','given_name',True,True,4,[]],
                  ['postcode','postcode',True,False,2,[]]]

    prev_rec_dict = {}

    for w in [2, 3, 4, 5, 6, 7, 9, 11, 12, 15, 20]:

      sort_index = indexing.SortingArrayIndex(description = \
                                              'Test SortingArray index',
                                         dataset1 = self.dataset1,
                                         dataset2 = self.dataset1,
                                         rec_comparator = self.rec_comp_dedupl,
                                         progress=4,
                                         window_s = w,
                                         skip_m = False,
                                         index_def = [index_def1,index_def2])

      assert isinstance(sort_index.index1, dict)
      assert isinstance(sort_index.index2, dict)
      assert isinstance(sort_index.description, str)
      assert isinstance(sort_index.rec_cache1, dict)
      assert isinstance(sort_index.rec_cache2, dict)
      assert sort_index.progress_report == 4
      assert sort_index.skip_missing == False
      assert sort_index.do_deduplication == True
      assert len(sort_index.index_def) == 2
      assert sort_index.window_size == w
      assert sort_index.num_rec_pairs == None

      assert sort_index.status == 'initialised'

      sort_index.build()  # - - - - - - - - - - - - - - - - - - - - - - - - - -

      assert sort_index.num_rec_pairs == None
      assert isinstance(sort_index.index1, dict)
      assert isinstance(sort_index.index2, dict)
      assert isinstance(sort_index.description, str)
      assert isinstance(sort_index.rec_cache1, dict)
      assert isinstance(sort_index.rec_cache2, dict)
      assert sort_index.progress_report == 4
      assert sort_index.skip_missing == False
      assert sort_index.do_deduplication == True
      assert len(sort_index.index_def) == 2
      assert sort_index.window_size == w

      assert sort_index.status == 'built'

      sort_index.compact()  # - - - - - - - - - - - - - - - - - - - - - - - - -

      assert sort_index.num_rec_pairs > 0
      assert sort_index.num_rec_pairs <= \
             len(self.rec_ident1)*len(self.rec_ident2)
      assert isinstance(sort_index.index1, dict)
      assert isinstance(sort_index.index2, dict)
      assert isinstance(sort_index.description, str)
      assert isinstance(sort_index.rec_cache1, dict)
      assert isinstance(sort_index.rec_cache2, dict)
      assert sort_index.progress_report == 4
      assert sort_index.skip_missing == False
      assert sort_index.do_deduplication == True
      assert len(sort_index.index_def) == 2
      assert sort_index.window_size == w

      assert sort_index.status == 'compacted'

      [field_names_list, weight_vec_dict] = sort_index.run()  # - - - - - - - -

      assert isinstance(weight_vec_dict, dict)
      assert len(weight_vec_dict) == sort_index.num_rec_pairs
      assert sort_index.num_rec_pairs > 0
      assert sort_index.num_rec_pairs <= \
             len(self.rec_ident1)*len(self.rec_ident2)

      for rec_ident1 in self.rec_ident1:
        assert (rec_ident1, rec_ident1) not in weight_vec_dict, \
               (rec_ident1, rec_ident1)

      for (rec_ident1,rec_ident2) in weight_vec_dict:
        assert rec_ident1 in self.rec_ident1
        assert rec_ident2 in self.rec_ident1
        assert isinstance(weight_vec_dict[(rec_ident1,rec_ident2)], list)

        assert len(weight_vec_dict[(rec_ident1,rec_ident2)]) == \
               len(sort_index.rec_comparator.field_comparison_list)

      # Build a sorting weight vector dict with pair identifiers sorted
      #
      for ((rec_ident1,rec_ident2), weight_list) in weight_vec_dict.items():
        if (rec_ident1 > rec_ident2):
          del weight_vec_dict[(rec_ident1,rec_ident2)]
          weight_vec_dict[(rec_ident2,rec_ident1)] = weight_list

      if (prev_rec_dict != {}):  # Only if there is a previous one
        for rec_ident_pair in prev_rec_dict.iterkeys():
          assert rec_ident_pair in weight_vec_dict, \
                 'Record pair %s not in weight vector dict with window = %d' % \
                 (str(rec_ident_pair), w)

      prev_rec_dict = weight_vec_dict

    # Test length filter
    #
    [field_names_list, prev_w_vec_dict] = sort_index.run()

    for lf in [99,50,30,20,15,10,5,1,0.5]:

      [field_names_list, this_w_vec_dict] = \
                                        sort_index.run(length_filter_perc = lf)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

    # Test cut-off threshold
    #
    [field_names_list, prev_w_vec_dict] = sort_index.run()

    for cot in [0.05, 0.1, 0.2, 0.5, 0.7, 0.8, 0.9]:

      [field_names_list, this_w_vec_dict] = \
                                        sort_index.run(cut_off_threshold = cot)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

  # ---------------------------------------------------------------------------

  def testQGramIndexLinkage(self):  # - - - - - - - - - - - - - - - - - - - -
    """Test QGramIndex linkage"""

    return

    index_def1 = [['surname','surname',False,False,None,[]]]
    index_def2 = [['given_name','given_name',True,True,4,[]],
                  ['postcode','postcode',True,False,2,[]]]

    # Assume blocking index is OK, so get the record pairs from this index for
    # comparison
    #
    block_index = indexing.BlockingIndex(description = 'Test blocking index',
                                         dataset1 = self.dataset1,
                                         dataset2 = self.dataset2,
                                         rec_comparator = self.rec_comp_link,
                                         progress=2,
                                         index_sep_str = 'xXx',
                                         index_def = [index_def1,index_def2])
    block_index.build()
    block_index.compact()
    [field_names_list, block_index_weight_vec_dict] = block_index.run()

    # Make sure record identifiers in each pair are sorted (smaller id first)
    #
    for ((rec_ident1,rec_ident2),w) in block_index_weight_vec_dict.items():
      if (rec_ident1 > rec_ident2):
        del block_index_weight_vec_dict[(rec_ident1,rec_ident2)]
        block_index_weight_vec_dict[(rec_ident2,rec_ident1)] = w

    prev_rec_dict = {}

    for t in [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.3, 0.1]:

      qgram_index = indexing.QGramIndex(description = 'Test Q-Gram index',
                                        dataset1 = self.dataset1,
                                        dataset2 = self.dataset2,
                                        rec_comparator = self.rec_comp_link,
                                        progress=2,
                                        padded=True,
                                        q= 2,
                                        thresh=t,
                                        index_sep_str = 'xXx',
                                        index_def = [index_def1,index_def2])

      assert isinstance(qgram_index.index1, dict)
      assert isinstance(qgram_index.index2, dict)
      assert isinstance(qgram_index.description, str)
      assert isinstance(qgram_index.rec_cache1, dict)
      assert isinstance(qgram_index.rec_cache2, dict)
      assert qgram_index.progress_report == 2
      assert qgram_index.skip_missing == True
      assert qgram_index.do_deduplication == False
      assert len(qgram_index.index_def) == 2
      assert qgram_index.index_sep_str == 'xXx'
      assert qgram_index.padded == True
      assert qgram_index.q == 2
      assert qgram_index.threshold == t
      assert qgram_index.num_rec_pairs == None

      assert qgram_index.status == 'initialised'

      qgram_index.build()  # - - - - - - - - - - - - - - - - - - - - - - - - -

      assert qgram_index.num_rec_pairs == None
      assert isinstance(qgram_index.index1, dict)
      assert isinstance(qgram_index.index2, dict)
      assert isinstance(qgram_index.description, str)
      assert isinstance(qgram_index.rec_cache1, dict)
      assert isinstance(qgram_index.rec_cache2, dict)
      assert qgram_index.progress_report == 2
      assert qgram_index.skip_missing == True
      assert qgram_index.do_deduplication == False
      assert len(qgram_index.index_def) == 2
      assert qgram_index.index_sep_str == 'xXx'
      assert qgram_index.padded == True
      assert qgram_index.q == 2
      assert qgram_index.threshold == t

      assert qgram_index.status == 'built'

      qgram_index.compact()  # - - - - - - - - - - - - - - - - - - - - - - - -

      assert qgram_index.num_rec_pairs > 0
      assert qgram_index.num_rec_pairs <= \
             len(self.rec_ident1)*len(self.rec_ident2)
      assert isinstance(qgram_index.index1, dict)
      assert isinstance(qgram_index.index2, dict)
      assert isinstance(qgram_index.description, str)
      assert isinstance(qgram_index.rec_cache1, dict)
      assert isinstance(qgram_index.rec_cache2, dict)
      assert qgram_index.progress_report == 2
      assert qgram_index.skip_missing == True
      assert qgram_index.do_deduplication == False
      assert len(qgram_index.index_def) == 2
      assert qgram_index.index_sep_str == 'xXx'
      assert qgram_index.padded == True
      assert qgram_index.q == 2
      assert qgram_index.threshold == t

      assert qgram_index.status == 'compacted'

      [field_names_list, weight_vec_dict] = qgram_index.run()  # - - - - - - -

      assert isinstance(weight_vec_dict, dict)
      assert len(weight_vec_dict) == qgram_index.num_rec_pairs
      assert qgram_index.num_rec_pairs > 0
      assert qgram_index.num_rec_pairs <= \
             len(self.rec_ident1)*len(self.rec_ident2)

      for rec_ident1 in self.rec_ident1:
        for rec_ident2 in self.rec_ident2:
          if (rec_ident1 == rec_ident2):
            assert (rec_ident1,rec_ident2) in weight_vec_dict

      for (rec_ident1,rec_ident2) in weight_vec_dict:
        assert rec_ident1 in self.rec_ident1
        assert rec_ident2 in self.rec_ident2
        assert isinstance(weight_vec_dict[(rec_ident1,rec_ident2)], list)

        assert len(weight_vec_dict[(rec_ident1,rec_ident2)]) == \
               len(qgram_index.rec_comparator.field_comparison_list)

      # All record pairs from the blocking index should be in the q-gram index
      # as well
      #
      assert len(block_index_weight_vec_dict) <= len(weight_vec_dict), \
             ('Number of record pairs in blocking index: %d, and in ' % \
             (len(block_index_weight_vec_dict))+' q-gram index: %d' % \
             (len(weight_vec_dict)))

      # Build a sorting weight vector dict with pair identifiers sorted
      #
      for ((rec_ident1,rec_ident2), weight_list) in weight_vec_dict.items():
        if (rec_ident1 > rec_ident2):
          del weight_vec_dict[(rec_ident1,rec_ident2)]
          weight_vec_dict[(rec_ident2,rec_ident1)] = weight_list

      for (rec_ident1,rec_ident2) in block_index_weight_vec_dict:
        assert (rec_ident1,rec_ident2) in weight_vec_dict, \
               'Record pair %s from blocking index not in q-gram index' % \
               (str((rec_ident1,rec_ident2)))

      for rec_ident_pair in prev_rec_dict.iterkeys():

        assert rec_ident_pair in weight_vec_dict, \
               'Record pair %s not in weight vector dict with thresh = %f' % \
               (str(rec_ident_pair), t)

      prev_rec_dict = weight_vec_dict

    # Q-gram index with threshold 1.0 should be same as blocking index - - - -
    #
    qgram_index = indexing.QGramIndex(description = 'Test Q-Gram index',
                                      dataset1 = self.dataset1,
                                      dataset2 = self.dataset2,
                                      rec_comparator = self.rec_comp_link,
                                      progress=4,
                                      padded=False,
                                      q= 3,
                                      log_funct = print_log,
                                      thresh=1.0,
                                      skip_m = False,
                                      index_sep_str = 'xXx',
                                      index_def = [index_def1,index_def2])
    qgram_index.build()
    qgram_index.compact()
    [field_names_list, weight_vec_dict] = qgram_index.run()

    # Make sure record identifiers in each pair are sorted (smaller id first)
    #
    for ((rec_ident1,rec_ident2),w) in weight_vec_dict.items():
      if (rec_ident1 > rec_ident2):
        del weight_vec_dict[(rec_ident1,rec_ident2)]
        weight_vec_dict[(rec_ident2,rec_ident1)] = w

    for (rec_ident1,rec_ident2) in weight_vec_dict.keys():
      assert (rec_ident1,rec_ident2) in block_index_weight_vec_dict, \
             'Not in blocking weight vectors: %s' % \
             (str((rec_ident1,rec_ident2)))

    for (rec_ident1,rec_ident2) in block_index_weight_vec_dict.keys():
      assert (rec_ident1,rec_ident2) in weight_vec_dict, \
             'Not in q-gram weight vectors: %s' % \
             (str((rec_ident1,rec_ident2)))

    # Test length filter
    #
    [field_names_list, prev_w_vec_dict] = qgram_index.run()

    for lf in [99,50,30,20,15,10,5,1,0.5]:

      [field_names_list, this_w_vec_dict] = \
                                       qgram_index.run(length_filter_perc = lf)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

    # Test cut-off threshold
    #
    [field_names_list, prev_w_vec_dict] = qgram_index.run()

    for cot in [0.05, 0.1, 0.2, 0.5, 0.7, 0.8, 0.9]:

      [field_names_list, this_w_vec_dict] = \
                                       qgram_index.run(cut_off_threshold = cot)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

  def testQGramIndexDedupl(self):  # - - - - - - - - - - - - - - - - - - - - -
    """Test QGramIndex deduplication"""

    return

    index_def1 = [['surname','surname',False,False,None,[]]]
    index_def2 = [['given_name','given_name',True,True,4,[]],
                  ['postcode','postcode',True,False,2,[]]]

    # Assume blocking index is OK, so get the record pairs from this index for
    # comparison
    #
    block_index = indexing.BlockingIndex(description = 'Test blocking index',
                                         dataset1 = self.dataset1,
                                         dataset2 = self.dataset1,
                                         rec_comparator = self.rec_comp_dedupl,
                                         progress=4,
                                         index_sep_str = '_',
                                         index_def = [index_def1,index_def2])
    block_index.build()
    block_index.compact()
    [field_names_list, block_index_weight_vec_dict] = block_index.run()

    # Make sure record identifiers in each pair are sorted (smaller id first)
    #
    for ((rec_ident1,rec_ident2),w) in block_index_weight_vec_dict.items():

      if (rec_ident1 > rec_ident2):
        del block_index_weight_vec_dict[(rec_ident1,rec_ident2)]
        block_index_weight_vec_dict[(rec_ident2,rec_ident1)] = w

    prev_rec_dict = {}

    for t in [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.3, 0.1]:

      qgram_index = indexing.QGramIndex(description = 'Test Q-Gram index',
                                        dataset1 = self.dataset1,
                                        dataset2 = self.dataset1,
                                        rec_comparator = self.rec_comp_dedupl,
                                        progress=4,
                                        padded=False,
                                        q= 3,
                                        thresh=t,
                                        skip_m = False,
                                        index_sep_str = '_',
                                        index_def = [index_def1,index_def2])

      assert isinstance(qgram_index.index1, dict)
      assert isinstance(qgram_index.index2, dict)
      assert isinstance(qgram_index.description, str)
      assert isinstance(qgram_index.rec_cache1, dict)
      assert isinstance(qgram_index.rec_cache2, dict)
      assert qgram_index.progress_report == 4
      assert qgram_index.skip_missing == False
      assert qgram_index.do_deduplication == True
      assert len(qgram_index.index_def) == 2
      assert qgram_index.index_sep_str == '_'
      assert qgram_index.padded == False
      assert qgram_index.q == 3
      assert qgram_index.threshold == t
      assert qgram_index.num_rec_pairs == None

      assert qgram_index.status == 'initialised'

      qgram_index.build()  # - - - - - - - - - - - - - - - - - - - - - - - - -

      assert qgram_index.num_rec_pairs == None
      assert isinstance(qgram_index.index1, dict)
      assert isinstance(qgram_index.index2, dict)
      assert isinstance(qgram_index.description, str)
      assert isinstance(qgram_index.rec_cache1, dict)
      assert isinstance(qgram_index.rec_cache2, dict)
      assert qgram_index.progress_report == 4
      assert qgram_index.skip_missing == False
      assert qgram_index.do_deduplication == True
      assert len(qgram_index.index_def) == 2
      assert qgram_index.index_sep_str == '_'
      assert qgram_index.padded == False
      assert qgram_index.q == 3
      assert qgram_index.threshold == t

      assert qgram_index.status == 'built'

      qgram_index.compact()  # - - - - - - - - - - - - - - - - - - - - - - - -

      assert qgram_index.num_rec_pairs > 0
      assert qgram_index.num_rec_pairs <= \
             len(self.rec_ident1)*len(self.rec_ident2)
      assert isinstance(qgram_index.index1, dict)
      assert isinstance(qgram_index.index2, dict)
      assert isinstance(qgram_index.description, str)
      assert isinstance(qgram_index.rec_cache1, dict)
      assert isinstance(qgram_index.rec_cache2, dict)
      assert qgram_index.progress_report == 4
      assert qgram_index.skip_missing == False
      assert qgram_index.do_deduplication == True
      assert len(qgram_index.index_def) == 2
      assert qgram_index.index_sep_str == '_'
      assert qgram_index.padded == False
      assert qgram_index.q == 3
      assert qgram_index.threshold == t

      assert qgram_index.status == 'compacted'

      [field_names_list, weight_vec_dict] = qgram_index.run()  # - - - - - - -

      assert isinstance(weight_vec_dict, dict)
      assert len(weight_vec_dict) == qgram_index.num_rec_pairs
      assert qgram_index.num_rec_pairs > 0
      assert qgram_index.num_rec_pairs <= \
             len(self.rec_ident1)*len(self.rec_ident2)

      for rec_ident1 in self.rec_ident1:
        assert (rec_ident1, rec_ident1) not in weight_vec_dict, \
               (rec_ident1, rec_ident1)

      for (rec_ident1,rec_ident2) in weight_vec_dict:
        assert rec_ident1 in self.rec_ident1
        assert rec_ident2 in self.rec_ident1
        assert isinstance(weight_vec_dict[(rec_ident1,rec_ident2)], list)

        assert len(weight_vec_dict[(rec_ident1,rec_ident2)]) == \
               len(qgram_index.rec_comparator.field_comparison_list)

      # All record pairs from the blocking index should be in the q-gram index
      # as well
      #
      assert len(block_index_weight_vec_dict) <= len(weight_vec_dict), \
             ('Number of record pairs in blocking index: %d, and in ' % \
             (len(block_index_weight_vec_dict))+' q-gram index: %d' % \
             (len(weight_vec_dict)))

      # Build a q-gram weight vector dict with pair identifiers sorted
      #
      for ((rec_ident1,rec_ident2), weight_list) in weight_vec_dict.items():
        if (rec_ident1 > rec_ident2):
          del weight_vec_dict[(rec_ident1,rec_ident2)]
          weight_vec_dict[(rec_ident2,rec_ident1)] = weight_list

      for (rec_ident1,rec_ident2) in block_index_weight_vec_dict:
        assert (rec_ident1,rec_ident2) in weight_vec_dict, \
               'Record pair %s from blocking index not in q-gram index' % \
               (str((rec_ident1,rec_ident2)))

      for rec_ident_pair in prev_rec_dict.iterkeys():
        assert rec_ident_pair in weight_vec_dict, \
               'Record pair %s not in weight vector dict with thresh = %f' % \
               (str(rec_ident_pair), t)

      prev_rec_dict = weight_vec_dict

    # Q-gram index with threshold 1.0 should be same as blocking index - - - -
    #
    qgram_index = indexing.QGramIndex(description = 'Test Q-Gram index',
                                      dataset1 = self.dataset1,
                                      dataset2 = self.dataset1,
                                      rec_comparator = self.rec_comp_dedupl,
                                      progress=4,
                                      padded=False,
                                      q= 3,
                                      log_funct = print_log,
                                      thresh=1.0,
                                      skip_m = False,
                                      index_sep_str = '_',
                                      index_def = [index_def1,index_def2])
    qgram_index.build()
    qgram_index.compact()
    [field_names_list, weight_vec_dict] = qgram_index.run()

    # Make sure record identifiers in each pair are sorted (smaller id first)
    #
    for ((rec_ident1,rec_ident2),w) in weight_vec_dict.items():
      if (rec_ident1 > rec_ident2):
        del weight_vec_dict[(rec_ident1,rec_ident2)]
        weight_vec_dict[(rec_ident2,rec_ident1)] = w

    for (rec_ident1,rec_ident2) in weight_vec_dict.keys():
      assert (rec_ident1,rec_ident2) in block_index_weight_vec_dict, \
             'Not in blocking weight vectors: %s' % \
             (str((rec_ident1,rec_ident2)))

    for (rec_ident1,rec_ident2) in block_index_weight_vec_dict.keys():
      assert (rec_ident1,rec_ident2) in weight_vec_dict, \
             'Not in q-gram weight vectors: %s' % \
             (str((rec_ident1,rec_ident2)))

    # Test length filter
    #
    [field_names_list, prev_w_vec_dict] = qgram_index.run()

    for lf in [99,50,30,20,15,10,5,1,0.5]:

      [field_names_list, this_w_vec_dict] = \
                                      qgram_index.run(length_filter_perc = lf)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

    # Test cut-off threshold
    #
    [field_names_list, prev_w_vec_dict] = qgram_index.run()

    for cot in [0.05, 0.1, 0.2, 0.5, 0.7, 0.8, 0.9]:

      [field_names_list, this_w_vec_dict] = \
                                       qgram_index.run(cut_off_threshold = cot)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

  # ---------------------------------------------------------------------------

  def testCanopyIndexLinkage(self):  # - - - - - - - - - - - - - - - - - - - -
    """Test CanopyIndex linkage"""

    return

    index_def1 = [['surname','surname',False,False,None,[]]]
    index_def2 = [['given_name','given_name',True,False,4,[]],
                  ['postcode','postcode',True,True,2,[]]]

    # Assume blocking index is OK, so get the record pairs from this index for
    # comparison
    #
    block_index = indexing.BlockingIndex(description = 'Test blocking index',
                                         dataset1 = self.dataset1,
                                         dataset2 = self.dataset2,
                                         rec_comparator = self.rec_comp_link,
                                         progress=2,
                                         index_def = [index_def1,index_def2])
    block_index.build()
    block_index.compact()
    [field_names_list, block_index_weight_vec_dict] = block_index.run()

    # Make sure record identifiers in each pair are sorted (smaller id first)
    #
    for ((rec_ident1,rec_ident2),w) in block_index_weight_vec_dict.items():
      if (rec_ident1 > rec_ident2):
        del block_index_weight_vec_dict[(rec_ident1,rec_ident2)]
        block_index_weight_vec_dict[(rec_ident2,rec_ident1)] = w

    for canopy_method_list in [[('tfidf', 'threshold', 0.999999, 0.999999),
                                ('tfidf', 'threshold', 0.999999, 0.9),
                                ('tfidf', 'threshold', 0.9, 0.8),
                                ('tfidf', 'threshold', 0.9, 0.7),
                                ('tfidf', 'threshold', 0.9, 0.5),
                                ('tfidf', 'threshold', 0.8, 0.5),
                                ('tfidf', 'threshold', 0.7, 0.5)],
                               [('tfidf', 'nearest', 1, 1),
                                ('tfidf', 'nearest', 1, 2),
                                ('tfidf', 'nearest', 1, 3),
                                ('tfidf', 'nearest', 1, 5)],
                               [('tfidf', 'nearest', 2, 2),
                                ('tfidf', 'nearest', 2, 3),
                                ('tfidf', 'nearest', 2, 4),
                                ('tfidf', 'nearest', 2, 5)],
                               [('tfidf', 'nearest', 4, 4),
                                ('tfidf', 'nearest', 4, 6)],
                               [('jaccard', 'threshold', 1.0, 1.0),
                                ('jaccard', 'threshold', 1.0, 0.9),
                                ('jaccard', 'threshold', 0.9, 0.8),
                                ('jaccard', 'threshold', 0.9, 0.7),
                                ('jaccard', 'threshold', 0.9, 0.5),
                                ('jaccard', 'threshold', 0.8, 0.5),
                                ('jaccard', 'threshold', 0.7, 0.5),
                                ('jaccard', 'threshold', 0.7, 0.4)],
                               [('jaccard', 'nearest', 1, 1),
                                ('jaccard', 'nearest', 1, 2),
                                ('jaccard', 'nearest', 1, 3),
                                ('jaccard', 'nearest', 1, 5)],
                               [('jaccard', 'nearest', 2, 2),
                                ('jaccard', 'nearest', 2, 3),
                                ('jaccard', 'nearest', 2, 4),
                                ('jaccard', 'nearest', 2, 5)],
                               [('jaccard', 'nearest', 4, 4),
                                ('jaccard', 'nearest', 4, 6)]]:

      prev_rec_dict = {}

      for canopy_method in canopy_method_list:

        canopy_index = indexing.CanopyIndex(description = 'Test canopy index',
                                            dataset1 = self.dataset1,
                                            dataset2 = self.dataset2,
                                            rec_compar = self.rec_comp_link,
                                            progress=2,
                                            padded=True,
                                            q= 2,
                                            canopy_me = canopy_method,
                                            delete_perc = 80,
                                            index_def = [index_def1,
                                                         index_def2])

        assert isinstance(canopy_index.index1, dict)
        assert isinstance(canopy_index.index2, dict)
        assert isinstance(canopy_index.description, str)
        assert isinstance(canopy_index.rec_cache1, dict)
        assert isinstance(canopy_index.rec_cache2, dict)
        assert canopy_index.progress_report == 2
        assert canopy_index.skip_missing == True
        assert canopy_index.do_deduplication == False
        assert len(canopy_index.index_def) == 2
        assert canopy_index.padded == True
        assert canopy_index.q == 2
        assert isinstance(canopy_index.canopy_method, tuple)
        assert canopy_index.delete_perc == 80
        assert canopy_index.num_rec_pairs == None

        assert canopy_index.status == 'initialised'

        canopy_index.build()  # - - - - - - - - - - - - - - - - - - - - - - - -

        assert canopy_index.num_rec_pairs == None
        assert isinstance(canopy_index.index1, dict)
        assert isinstance(canopy_index.index2, dict)
        assert isinstance(canopy_index.description, str)
        assert isinstance(canopy_index.rec_cache1, dict)
        assert isinstance(canopy_index.rec_cache2, dict)
        assert canopy_index.progress_report == 2
        assert canopy_index.skip_missing == True
        assert canopy_index.do_deduplication == False
        assert len(canopy_index.index_def) == 2
        assert canopy_index.padded == True
        assert canopy_index.q == 2
        assert isinstance(canopy_index.canopy_method, tuple)
        assert canopy_index.delete_perc == 80
        assert canopy_index.num_rec_pairs == None

        assert canopy_index.status == 'built'

        canopy_index.compact()  # - - - - - - - - - - - - - - - - - - - - - - -

        assert canopy_index.num_rec_pairs > 0
        assert canopy_index.num_rec_pairs <= \
               len(self.rec_ident1)*len(self.rec_ident2)
        assert isinstance(canopy_index.index1, dict)
        assert isinstance(canopy_index.index2, dict)
        assert isinstance(canopy_index.description, str)
        assert isinstance(canopy_index.rec_cache1, dict)
        assert isinstance(canopy_index.rec_cache2, dict)
        assert canopy_index.progress_report == 2
        assert canopy_index.skip_missing == True
        assert canopy_index.do_deduplication == False
        assert len(canopy_index.index_def) == 2
        assert canopy_index.padded == True
        assert canopy_index.q == 2
        assert isinstance(canopy_index.canopy_method, tuple)
        assert canopy_index.delete_perc == 80
        assert canopy_index.status == 'compacted'

        [field_names_list, weight_vec_dict] = canopy_index.run()  # - - - - - -

        assert isinstance(weight_vec_dict, dict)
        assert len(weight_vec_dict) == canopy_index.num_rec_pairs
        assert canopy_index.num_rec_pairs > 0
        assert canopy_index.num_rec_pairs <= \
               len(self.rec_ident1)*len(self.rec_ident2)

        for rec_ident1 in self.rec_ident1:
          for rec_ident2 in self.rec_ident2:
            if (rec_ident1 == rec_ident2):
              assert (rec_ident1,rec_ident2) in weight_vec_dict, \
                     (rec_ident1,rec_ident2)

        for (rec_ident1,rec_ident2) in weight_vec_dict:
          assert rec_ident1 in self.rec_ident1
          assert rec_ident2 in self.rec_ident2
          assert isinstance(weight_vec_dict[(rec_ident1,rec_ident2)], list)

          assert len(weight_vec_dict[(rec_ident1,rec_ident2)]) == \
                 len(canopy_index.rec_comparator.field_comparison_list)

        # All record pairs from the blocking index should be in the canopy
        # index as well
        #
        assert len(block_index_weight_vec_dict) <= len(weight_vec_dict), \
               ('Number of record pairs in blocking index: %d, and in ' % \
               (len(block_index_weight_vec_dict))+' canopy index: %d' % \
               (len(weight_vec_dict)))

        # Build a canopy weight vector dict with pair identifiers sorted
        #
        for ((rec_ident1,rec_ident2), weight_list) in weight_vec_dict.items():
          if (rec_ident1 > rec_ident2):
            del weight_vec_dict[(rec_ident1,rec_ident2)]
            weight_vec_dict[(rec_ident2,rec_ident1)] = weight_list

        for (rec_ident1,rec_ident2) in block_index_weight_vec_dict:
          assert (rec_ident1,rec_ident2) in weight_vec_dict, \
                 'Record pair %s from blocking index not in canopy index' % \
                 (str((rec_ident1,rec_ident2)))

        for rec_ident_pair in prev_rec_dict.iterkeys():

          assert rec_ident_pair in weight_vec_dict, \
                 'Record pair %s not in weight vector dict with canopy ' % \
                 (str(rec_ident_pair)) + 'method: %s' % (str(canopy_method))

        prev_rec_dict = weight_vec_dict

    # Canopy index with threshold 1.0 / nearest 1 should be same as blocking -
    #
    for canopy_method in [('tfidf', 'threshold', 0.999999, 0.999999),
                          ('tfidf', 'nearest', 1, 1),
                          ('jaccard', 'threshold', 1.0, 1.0),
                          ('jaccard', 'nearest', 1, 1)]:

      canopy_index = indexing.CanopyIndex(description = 'Test canopy index',
                                          dataset1 = self.dataset1,
                                          dataset2 = self.dataset2,
                                          rec_compar = self.rec_comp_link,
                                          progress=2,
                                          padded=True,
                                          q= 2,
                                          canopy_me = canopy_method,
                                          delete_perc = 90,
                                          index_def = [index_def1,index_def2])
      canopy_index.build()
      canopy_index.compact()
      [field_names_list, weight_vec_dict] = canopy_index.run()

      # Make sure record identifiers in each pair are sorted (smaller id first)
      #
      for ((rec_ident1,rec_ident2),w) in weight_vec_dict.items():
        if (rec_ident1 > rec_ident2):
          del weight_vec_dict[(rec_ident1,rec_ident2)]
          weight_vec_dict[(rec_ident2,rec_ident1)] = w

      for (rec_ident1,rec_ident2) in weight_vec_dict.keys():
        assert (rec_ident1,rec_ident2) in block_index_weight_vec_dict, \
               'Not in blocking weight vectors: %s' % \
               (str((rec_ident1,rec_ident2)))

      for (rec_ident1,rec_ident2) in block_index_weight_vec_dict.keys():
        assert (rec_ident1,rec_ident2) in weight_vec_dict, \
               'Not in canopy weight vectors: %s' % \
              (str((rec_ident1,rec_ident2)))

    # Test length filter
    #
    [field_names_list, prev_w_vec_dict] = canopy_index.run()

    for lf in [99,50,30,20,15,10,5,1,0.5]:

      [field_names_list, this_w_vec_dict] = \
                                      canopy_index.run(length_filter_perc = lf)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

    # Test cut-off threshold
    #
    [field_names_list, prev_w_vec_dict] = canopy_index.run()

    for cot in [0.05, 0.1, 0.2, 0.5, 0.7, 0.8, 0.9]:

      [field_names_list, this_w_vec_dict] = \
                                      canopy_index.run(cut_off_threshold = cot)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

  def testCanopyIndexDedupl(self):  # - - - - - - - - - - - - - - - - - - - - -
    """Test CanopyIndex deduplication"""

    return

    index_def1 = [['surname','surname',False,False,None,[]]]
    index_def2 = [['given_name','given_name',True,False,4,[]],
                  ['postcode','postcode',True,True,2,[]]]

    # Assume blocking index is OK, so get the record pairs from this index for
    # comparison
    #
    block_index = indexing.BlockingIndex(description = 'Test blocking index',
                                         dataset1 = self.dataset1,
                                         dataset2 = self.dataset1,
                                         rec_comparator = self.rec_comp_dedupl,
                                         progress=4,
                                         index_sep_str = '123',
                                         index_def = [index_def1,index_def2])
    block_index.build()
    block_index.compact()
    [field_names_list, block_index_weight_vec_dict] = block_index.run()

    # Make sure record identifiers in each pair are sorted (smaller id first)
    #
    for ((rec_ident1,rec_ident2),w) in block_index_weight_vec_dict.items():

      if (rec_ident1 > rec_ident2):
        del block_index_weight_vec_dict[(rec_ident1,rec_ident2)]
        block_index_weight_vec_dict[(rec_ident2,rec_ident1)] = w

    for canopy_method_list in [[('tfidf', 'threshold', 0.999999, 0.999999),
                                ('tfidf', 'threshold', 0.999999, 0.9),
                                ('tfidf', 'threshold', 0.9, 0.8),
                                ('tfidf', 'threshold', 0.9, 0.7),
                                ('tfidf', 'threshold', 0.9, 0.5),
                                ('tfidf', 'threshold', 0.8, 0.5),
                                ('tfidf', 'threshold', 0.7, 0.5)],
                               [('tfidf', 'nearest', 1, 1),
                                ('tfidf', 'nearest', 1, 2),
                                ('tfidf', 'nearest', 1, 3),
                                ('tfidf', 'nearest', 1, 5)],
                               [('tfidf', 'nearest', 2, 2),
                                ('tfidf', 'nearest', 2, 3),
                                ('tfidf', 'nearest', 2, 4),
                                ('tfidf', 'nearest', 2, 5)],
                               [('tfidf', 'nearest', 4, 4),
                                ('tfidf', 'nearest', 4, 6)],
                               [('jaccard', 'threshold', 1.0, 1.0),
                                ('jaccard', 'threshold', 1.0, 0.9),
                                ('jaccard', 'threshold', 0.9, 0.8),
                                ('jaccard', 'threshold', 0.9, 0.7),
                                ('jaccard', 'threshold', 0.9, 0.5),
                                ('jaccard', 'threshold', 0.8, 0.5),
                                ('jaccard', 'threshold', 0.7, 0.5),
                                ('jaccard', 'threshold', 0.7, 0.4)],
                               [('jaccard', 'nearest', 1, 1),
                                ('jaccard', 'nearest', 1, 2),
                                ('jaccard', 'nearest', 1, 3),
                                ('jaccard', 'nearest', 1, 5)],
                               [('jaccard', 'nearest', 2, 2),
                                ('jaccard', 'nearest', 2, 3),
                                ('jaccard', 'nearest', 2, 4),
                                ('jaccard', 'nearest', 2, 5)],
                               [('jaccard', 'nearest', 4, 4),
                                ('jaccard', 'nearest', 4, 6)]]:

      prev_rec_dict = {}

      for canopy_method in canopy_method_list:

        canopy_index = indexing.CanopyIndex(description = 'Test canopy index',
                                            dataset1 = self.dataset1,
                                            dataset2 = self.dataset1,
                                            rec_compar =self.rec_comp_dedupl,
                                            progress=4,
                                            padded=False,
                                            q= 3,
                                            log_funct = print_log,
                                            canopy_me = canopy_method,
                                            delete_perc = 60,
                                            index_sep_str = '123',
                                            index_def = [index_def1,
                                                         index_def2])

        assert isinstance(canopy_index.index1, dict)
        assert isinstance(canopy_index.index2, dict)
        assert isinstance(canopy_index.description, str)
        assert isinstance(canopy_index.rec_cache1, dict)
        assert isinstance(canopy_index.rec_cache2, dict)
        assert canopy_index.progress_report == 4
        assert canopy_index.skip_missing == True
        assert canopy_index.do_deduplication == True
        assert len(canopy_index.index_def) == 2
        assert canopy_index.index_sep_str == '123'
        assert canopy_index.padded == False
        assert canopy_index.q == 3
        assert isinstance(canopy_index.canopy_method, tuple)
        assert canopy_index.delete_perc == 60
        assert canopy_index.num_rec_pairs == None

        assert canopy_index.status == 'initialised'

        canopy_index.build()  # - - - - - - - - - - - - - - - - - - - - - - - -

        assert isinstance(canopy_index.index1, dict)
        assert isinstance(canopy_index.index2, dict)
        assert isinstance(canopy_index.description, str)
        assert isinstance(canopy_index.rec_cache1, dict)
        assert isinstance(canopy_index.rec_cache2, dict)
        assert canopy_index.progress_report == 4
        assert canopy_index.skip_missing == True
        assert canopy_index.do_deduplication == True
        assert len(canopy_index.index_def) == 2
        assert canopy_index.index_sep_str == '123'
        assert canopy_index.padded == False
        assert canopy_index.q == 3
        assert isinstance(canopy_index.canopy_method, tuple)
        assert canopy_index.delete_perc == 60
        assert canopy_index.num_rec_pairs == None

        assert canopy_index.status == 'built'

        canopy_index.compact()  # - - - - - - - - - - - - - - - - - - - - - - -

        assert canopy_index.num_rec_pairs > 0
        assert canopy_index.num_rec_pairs <= \
               len(self.rec_ident1)*len(self.rec_ident2)
        assert isinstance(canopy_index.index1, dict)
        assert isinstance(canopy_index.index2, dict)
        assert isinstance(canopy_index.description, str)
        assert isinstance(canopy_index.rec_cache1, dict)
        assert isinstance(canopy_index.rec_cache2, dict)
        assert canopy_index.progress_report == 4
        assert canopy_index.skip_missing == True
        assert canopy_index.do_deduplication == True
        assert len(canopy_index.index_def) == 2
        assert canopy_index.index_sep_str == '123'
        assert canopy_index.padded == False
        assert canopy_index.q == 3
        assert isinstance(canopy_index.canopy_method, tuple)
        assert canopy_index.delete_perc == 60
        assert canopy_index.status == 'compacted'

        [field_names_list, weight_vec_dict] = canopy_index.run()  # - - - - - -

        assert isinstance(weight_vec_dict, dict)
        assert len(weight_vec_dict) == canopy_index.num_rec_pairs
        assert canopy_index.num_rec_pairs > 0
        assert canopy_index.num_rec_pairs <= \
               len(self.rec_ident1)*len(self.rec_ident2)

        for rec_ident1 in self.rec_ident1:
          assert (rec_ident1, rec_ident1) not in weight_vec_dict, \
                 (rec_ident1, rec_ident1)

        for (rec_ident1,rec_ident2) in weight_vec_dict:
          assert rec_ident1 in self.rec_ident1
          assert rec_ident2 in self.rec_ident1
          assert isinstance(weight_vec_dict[(rec_ident1,rec_ident2)], list)

          assert len(weight_vec_dict[(rec_ident1,rec_ident2)]) == \
                 len(canopy_index.rec_comparator.field_comparison_list)

        # All record pairs from the blocking index should be in canopy index
        # as well
        #
        assert len(block_index_weight_vec_dict) <= len(weight_vec_dict), \
               ('Number of record pairs in blocking index: %d, and in ' % \
               (len(block_index_weight_vec_dict))+' canopy index: %d' % \
               (len(weight_vec_dict)))

        # Build a canopy weight vector dict with pair identifiers sorted
        #
        for ((rec_ident1,rec_ident2), weight_list) in weight_vec_dict.items():
          if (rec_ident1 > rec_ident2):
            del weight_vec_dict[(rec_ident1,rec_ident2)]
            weight_vec_dict[(rec_ident2,rec_ident1)] = weight_list

        for (rec_ident1,rec_ident2) in block_index_weight_vec_dict:
          assert (rec_ident1,rec_ident2) in weight_vec_dict, \
                 'Record pair %s from blocking index not in canopy index' % \
                 (str((rec_ident1,rec_ident2)))

        for rec_ident_pair in prev_rec_dict.iterkeys():
          assert rec_ident_pair in weight_vec_dict, \
                 'Record pair %s not in weight vector dict with canopy ' % \
                 (str(rec_ident_pair)) + 'method: %s' % (str(canopy_method))

        prev_rec_dict = weight_vec_dict

    # Canopy index with threshold 1.0 / nearest 1 should be same as blocking -
    #
    for canopy_method in [('tfidf', 'threshold', 0.999999, 0.999999),
                          ('tfidf', 'nearest', 1, 1),
                          ('jaccard', 'threshold', 1.0, 1.0),
                          ('jaccard', 'nearest', 1, 1)]:

      canopy_index = indexing.CanopyIndex(description = 'Test canopy index',
                                          dataset1 = self.dataset1,
                                          dataset2 = self.dataset1,
                                          rec_compar = self.rec_comp_dedupl,
                                          progress=4,
                                          padded=False,
                                          q= 3,
                                          canopy_me = canopy_method,
                                          delete_perc = 90,
                                          index_sep_str = '123',
                                          index_def = [index_def1,index_def2])
      canopy_index.build()
      canopy_index.compact()
      [field_names_list, weight_vec_dict] = canopy_index.run()

      # Make sure record identifiers in each pair are sorted (smaller id first)
      #
      for ((rec_ident1,rec_ident2),w) in weight_vec_dict.items():
        if (rec_ident1 > rec_ident2):
          del weight_vec_dict[(rec_ident1,rec_ident2)]
          weight_vec_dict[(rec_ident2,rec_ident1)] = w

      for (rec_ident1,rec_ident2) in weight_vec_dict.keys():
        assert (rec_ident1,rec_ident2) in block_index_weight_vec_dict, \
               'Not in blocking weight vectors: %s' % \
               (str((rec_ident1,rec_ident2)))

      for (rec_ident1,rec_ident2) in block_index_weight_vec_dict.keys():
        assert (rec_ident1,rec_ident2) in weight_vec_dict, \
               'Not in canopy weight vectors: %s' % \
              (str((rec_ident1,rec_ident2)))

    # Test length filter
    #
    [field_names_list, prev_w_vec_dict] = canopy_index.run()

    for lf in [99,50,30,20,15,10,5,1,0.5]:

      [field_names_list, this_w_vec_dict] = \
                                      canopy_index.run(length_filter_perc = lf)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

    # Test cut-off threshold
    #
    [field_names_list, prev_w_vec_dict] = canopy_index.run()

    for cot in [0.05, 0.1, 0.2, 0.5, 0.7, 0.8, 0.9]:

      [field_names_list, this_w_vec_dict] = \
                                      canopy_index.run(cut_off_threshold = cot)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

  # ---------------------------------------------------------------------------

  def testStringMapIndexLinkage(self):  # - - - - - - - - - - - - - - - - - - -
    """Test StringMapIndex linkage"""

    return

    index_def1 = [['surname','surname',False,False,None,[]]]
    index_def2 = [['given_name','given_name',True,True,4,[]],
                  ['postcode','postcode',True,False,2,[]]]

    # Assume blocking index is OK, so get the record pairs from this index for
    # comparison
    #
    block_index = indexing.BlockingIndex(description = 'Test blocking index',
                                         dataset1 = self.dataset1,
                                         dataset2 = self.dataset2,
                                         rec_comparator = self.rec_comp_link,
                                         progress=2,
                                         index_sep_str = chr(3),
                                         index_def = [index_def1,index_def2])
    block_index.build()
    block_index.compact()
    [field_names_list, block_index_weight_vec_dict] = block_index.run()

    # Make sure record identifiers in each pair are sorted (smaller id first)
    #
    for ((rec_ident1,rec_ident2),w) in block_index_weight_vec_dict.items():
      if (rec_ident1 > rec_ident2):
        del block_index_weight_vec_dict[(rec_ident1,rec_ident2)]
        block_index_weight_vec_dict[(rec_ident2,rec_ident1)] = w

    for canopy_method_list in [[('nearest', 1, 1),
                                ('nearest', 1, 2),
                                ('nearest', 1, 3),
                                ('nearest', 1, 5)],
                               [('nearest', 2, 2),
                                ('nearest', 2, 3),
                                ('nearest', 2, 4)],
                               [('nearest', 3, 3),
                                ('nearest', 4, 6)],
                               [('threshold', 1.0, 1.0),
                                ('threshold', 1.0, 0.9),
                                ('threshold', 1.0, 0.7)],
                               [('threshold', 0.9, 0.9),
                                ('threshold', 0.9, 0.7),
                                ('threshold', 0.9, 0.5)],
                               [('threshold', 0.7, 0.7),
                                ('threshold', 0.7, 0.2)]]:

      for (dim, sub_dim) in [(10,1),(10,2),(15,1),(15,2),(15,3),(20,2),(20,3)]:

        prev_rec_dict = {}

        for canopy_method in canopy_method_list:

          strmap_index = indexing.StringMapIndex(descri = 'Test str-map index',
                                                 dataset1 = self.dataset1,
                                                 dataset2 = self.dataset2,
                                                 rec_comp = self.rec_comp_link,
                                                 progress=2,
                                                 canopy_me = canopy_method,
                                                 dim = dim,
                                                 sub_d = sub_dim,
                                                 grid_resol = 10,
                                                 sim_func = stringcmp.editdist,
                                                 index_sep_str = chr(3),
                                                 index_def = [index_def1,
                                                              index_def2])

          assert isinstance(strmap_index.index1, dict)
          assert isinstance(strmap_index.index2, dict)
          assert isinstance(strmap_index.description, str)
          assert isinstance(strmap_index.rec_cache1, dict)
          assert isinstance(strmap_index.rec_cache2, dict)
          assert strmap_index.progress_report == 2
          assert strmap_index.skip_missing == True
          assert strmap_index.do_deduplication == False
          assert len(strmap_index.index_def) == 2
          assert strmap_index.index_sep_str == chr(3)
          assert strmap_index.dim == dim
          assert strmap_index.sub_dim == sub_dim
          assert strmap_index.grid_resolution == 10
          assert isinstance(strmap_index.canopy_method, tuple)
          assert strmap_index.sim_funct == stringcmp.editdist
          assert strmap_index.num_rec_pairs == None

          assert strmap_index.status == 'initialised'

          strmap_index.build()  # - - - - - - - - - - - - - - - - - - - - - - -

          assert strmap_index.num_rec_pairs == None
          assert isinstance(strmap_index.index1, dict)
          assert isinstance(strmap_index.index2, dict)
          assert isinstance(strmap_index.description, str)
          assert isinstance(strmap_index.rec_cache1, dict)
          assert isinstance(strmap_index.rec_cache2, dict)
          assert strmap_index.progress_report == 2
          assert strmap_index.skip_missing == True
          assert strmap_index.do_deduplication == False
          assert len(strmap_index.index_def) == 2
          assert strmap_index.index_sep_str == chr(3)
          assert strmap_index.dim == dim
          assert strmap_index.sub_dim == sub_dim
          assert strmap_index.grid_resolution == 10
          assert isinstance(strmap_index.canopy_method, tuple)
          assert strmap_index.sim_funct == stringcmp.editdist
          assert strmap_index.num_rec_pairs == None

          assert strmap_index.status == 'built'

          strmap_index.compact()  # - - - - - - - - - - - - - - - - - - - - - -

          assert strmap_index.num_rec_pairs > 0
          assert strmap_index.num_rec_pairs <= \
                len(self.rec_ident1)*len(self.rec_ident2)
          assert isinstance(strmap_index.index1, dict)
          assert isinstance(strmap_index.index2, dict)
          assert isinstance(strmap_index.description, str)
          assert isinstance(strmap_index.rec_cache1, dict)
          assert isinstance(strmap_index.rec_cache2, dict)
          assert strmap_index.progress_report == 2
          assert strmap_index.skip_missing == True
          assert strmap_index.do_deduplication == False
          assert len(strmap_index.index_def) == 2
          assert strmap_index.index_sep_str == chr(3)
          assert strmap_index.dim == dim
          assert strmap_index.sub_dim == sub_dim
          assert strmap_index.grid_resolution == 10
          assert isinstance(strmap_index.canopy_method, tuple)
          assert strmap_index.sim_funct == stringcmp.editdist
          assert strmap_index.status == 'compacted'

          [field_names_list, weight_vec_dict] = strmap_index.run()  # - - - - -

          assert isinstance(weight_vec_dict, dict)
          assert len(weight_vec_dict) == strmap_index.num_rec_pairs
          assert strmap_index.num_rec_pairs > 0
          assert strmap_index.num_rec_pairs <= \
                 len(self.rec_ident1)*len(self.rec_ident2)

          for rec_ident1 in self.rec_ident1:
            for rec_ident2 in self.rec_ident2:
              if (rec_ident1 == rec_ident2):
                assert (rec_ident1,rec_ident2) in weight_vec_dict, \
                       (rec_ident1,rec_ident2)

          for (rec_ident1,rec_ident2) in weight_vec_dict:
            assert rec_ident1 in self.rec_ident1
            assert rec_ident2 in self.rec_ident2
            assert isinstance(weight_vec_dict[(rec_ident1,rec_ident2)], list)

            assert len(weight_vec_dict[(rec_ident1,rec_ident2)]) == \
                   len(strmap_index.rec_comparator.field_comparison_list)

          # All record pairs from the blocking index should be in the string-
          # map index as well
          #
          assert len(block_index_weight_vec_dict) <= len(weight_vec_dict), \
                 ('Number of record pairs in blocking index: %d, and in ' % \
                 (len(block_index_weight_vec_dict))+' string map index: %d' % \
                 (len(weight_vec_dict)))

          # Build a string-map weight vector dict with pair identifiers sorted
          #
          for ((rec_ident1,rec_ident2),weight_list) in weight_vec_dict.items():
            if (rec_ident1 > rec_ident2):
              del weight_vec_dict[(rec_ident1,rec_ident2)]
              weight_vec_dict[(rec_ident2,rec_ident1)] = weight_list

          for (rec_ident1,rec_ident2) in block_index_weight_vec_dict:
            assert (rec_ident1,rec_ident2) in weight_vec_dict, \
                   'Record pair %s from blocking index not in StrMap index' % \
                   (str((rec_ident1,rec_ident2)))

          for rec_ident_pair in prev_rec_dict.iterkeys():
            assert rec_ident_pair in weight_vec_dict, \
                   'Record pair %s not in weight vector dict with canopy ' % \
                   (str(rec_ident_pair)) + 'method: %s' % (str(canopy_method))

          prev_rec_dict = weight_vec_dict

    # Str-map index with threshold 1.0 / nearest 1 should be same as blocking
    #
    for canopy_method in [('threshold', 1.0, 1.0),
                          ('nearest', 1, 1)]:

      for (dim, sub_dim) in [(10,3), (10,2), (15,3), (15,2), (20,2), (20,3)]:

        strmap_index = indexing.StringMapIndex(descri = 'Test str-map index',
                                               dataset1 = self.dataset1,
                                               dataset2 = self.dataset2,
                                               rec_compa = self.rec_comp_link,
                                               progress=2,
                                               canopy_me = canopy_method,
                                               dim = dim,
                                               log_funct = print_log,
                                               sub_d = sub_dim,
                                               grid_resol = 10,
                                               sim_funct = stringcmp.editdist,
                                               index_sep_str = chr(3),
                                               index_def = [index_def1,
                                               index_def2])
        strmap_index.build()
        strmap_index.compact()
        [field_names_list, weight_vec_dict] = strmap_index.run()

        # Make sure record ident. in each pair are sorted (smaller id first)
        #
        for ((rec_ident1,rec_ident2),w) in weight_vec_dict.items():
          if (rec_ident1 > rec_ident2):
            del weight_vec_dict[(rec_ident1,rec_ident2)]
            weight_vec_dict[(rec_ident2,rec_ident1)] = w

        for (rec_ident1,rec_ident2) in weight_vec_dict.keys():
          assert (rec_ident1,rec_ident2) in block_index_weight_vec_dict, \
                 'Not in blocking weight vectors: %s' % \
                 (str((rec_ident1,rec_ident2)))

        for (rec_ident1,rec_ident2) in block_index_weight_vec_dict.keys():
          assert (rec_ident1,rec_ident2) in weight_vec_dict, \
                 'Not in string-map weight vectors: %s' % \
                (str((rec_ident1,rec_ident2)))

    # Test length filter
    #
    [field_names_list, prev_w_vec_dict] = strmap_index.run()

    for lf in [99,50,30,20,15,10,5,1,0.5]:

      [field_names_list, this_w_vec_dict] = \
                                      strmap_index.run(length_filter_perc = lf)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

    # Test cut-off threshold
    #
    [field_names_list, prev_w_vec_dict] = strmap_index.run()

    for cot in [0.05, 0.1, 0.2, 0.5, 0.7, 0.8, 0.9]:

      [field_names_list, this_w_vec_dict] = \
                                      strmap_index.run(cut_off_threshold = cot)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

  def testStringMapIndexDedupl(self):  # - - - - - - - - - - - - - - - - - - -
    """Test StringMapIndex deduplication"""

    return

    index_def1 = [['surname','surname',False,False,None,[]]]
    index_def2 = [['given_name','given_name',True,True,4,[]],
                  ['postcode','postcode',True,False,2,[]]]

    # Assume blocking index is OK, so get the record pairs from this index for
    # comparison
    #
    block_index = indexing.BlockingIndex(description = 'Test blocking index',
                                         dataset1 = self.dataset1,
                                         dataset2 = self.dataset1,
                                         rec_comparator = self.rec_comp_dedupl,
                                         progress=4,
                                         index_sep_str = '-',
                                         index_def = [index_def1,index_def2])
    block_index.build()
    block_index.compact()
    [field_names_list, block_index_weight_vec_dict] = block_index.run()

    # Make sure record identifiers in each pair are sorted (smaller id first)
    #
    for ((rec_ident1,rec_ident2),w) in block_index_weight_vec_dict.items():

      if (rec_ident1 > rec_ident2):
        del block_index_weight_vec_dict[(rec_ident1,rec_ident2)]
        block_index_weight_vec_dict[(rec_ident2,rec_ident1)] = w

    for canopy_method_list in [[('nearest', 1, 1),
                                ('nearest', 1, 2),
                                ('nearest', 1, 3),
                                ('nearest', 1, 5)],
                               [('nearest', 2, 2),
                                ('nearest', 2, 3),
                                ('nearest', 2, 4)],
                               [('nearest', 3, 3),
                                ('nearest', 4, 6)],
                               [('threshold', 1.0, 1.0),
                                ('threshold', 1.0, 0.9),
                                ('threshold', 1.0, 0.7)],
                               [('threshold', 0.9, 0.9),
                                ('threshold', 0.9, 0.7),
                                ('threshold', 0.9, 0.5)],
                               [('threshold', 0.7, 0.7),
                                ('threshold', 0.7, 0.3)]]:

      for (dim, sub_dim) in [(10,2),(10,3),(15,2),(15,3),(20,2),(20,3)]:

        prev_rec_dict = {}

        for canopy_method in canopy_method_list:

          strmap_index = indexing.StringMapIndex(descri = 'Test str-map index',
                                                 dataset1 = self.dataset1,
                                                 dataset2 = self.dataset1,
                                                 rec_comp=self.rec_comp_dedupl,
                                                 progress= 4,
                                                 canopy_me = canopy_method,
                                                 dim = dim,
                                                 sub_d = sub_dim,
                                                 grid_resol = 10,
                                                 sim_func = stringcmp.editdist,
                                                 index_sep_str = '-',
                                                 index_def = [index_def1,
                                                              index_def2])

          assert isinstance(strmap_index.index1, dict)
          assert isinstance(strmap_index.index2, dict)
          assert isinstance(strmap_index.description, str)
          assert isinstance(strmap_index.rec_cache1, dict)
          assert isinstance(strmap_index.rec_cache2, dict)
          assert strmap_index.progress_report == 4
          assert strmap_index.skip_missing == True
          assert strmap_index.do_deduplication == True
          assert len(strmap_index.index_def) == 2
          assert strmap_index.index_sep_str == '-'
          assert strmap_index.dim == dim
          assert strmap_index.sub_dim == sub_dim
          assert strmap_index.grid_resolution == 10
          assert isinstance(strmap_index.canopy_method, tuple)
          assert strmap_index.sim_funct == stringcmp.editdist
          assert strmap_index.num_rec_pairs == None

          assert strmap_index.status == 'initialised'

          strmap_index.build()  # - - - - - - - - - - - - - - - - - - - - - - -

          assert strmap_index.num_rec_pairs == None
          assert isinstance(strmap_index.index1, dict)
          assert isinstance(strmap_index.index2, dict)
          assert isinstance(strmap_index.description, str)
          assert isinstance(strmap_index.rec_cache1, dict)
          assert isinstance(strmap_index.rec_cache2, dict)
          assert strmap_index.progress_report == 4
          assert strmap_index.skip_missing == True
          assert strmap_index.do_deduplication == True
          assert len(strmap_index.index_def) == 2
          assert strmap_index.index_sep_str == '-'
          assert strmap_index.dim == dim
          assert strmap_index.sub_dim == sub_dim
          assert strmap_index.grid_resolution == 10
          assert isinstance(strmap_index.canopy_method, tuple)
          assert strmap_index.sim_funct == stringcmp.editdist
          assert strmap_index.num_rec_pairs == None

          assert strmap_index.status == 'built'

          strmap_index.compact()  # - - - - - - - - - - - - - - - - - - - - - -

          assert strmap_index.num_rec_pairs > 0
          assert strmap_index.num_rec_pairs <= \
                len(self.rec_ident1)*len(self.rec_ident2)
          assert isinstance(strmap_index.index1, dict)
          assert isinstance(strmap_index.index2, dict)
          assert isinstance(strmap_index.description, str)
          assert isinstance(strmap_index.rec_cache1, dict)
          assert isinstance(strmap_index.rec_cache2, dict)
          assert strmap_index.progress_report == 4
          assert strmap_index.skip_missing == True
          assert strmap_index.do_deduplication == True
          assert len(strmap_index.index_def) == 2
          assert strmap_index.index_sep_str == '-'
          assert strmap_index.dim == dim
          assert strmap_index.sub_dim == sub_dim
          assert strmap_index.grid_resolution == 10
          assert isinstance(strmap_index.canopy_method, tuple)
          assert strmap_index.sim_funct == stringcmp.editdist
          assert strmap_index.status == 'compacted'

          [field_names_list, weight_vec_dict] = strmap_index.run()  # - - - - -

          assert isinstance(weight_vec_dict, dict)
          assert len(weight_vec_dict) == strmap_index.num_rec_pairs
          assert strmap_index.num_rec_pairs > 0
          assert strmap_index.num_rec_pairs <= \
                 len(self.rec_ident1)*(len(self.rec_ident1)-1)/2   ##########

          for rec_ident1 in self.rec_ident1:
            assert (rec_ident1, rec_ident1) not in weight_vec_dict, \
                   (rec_ident1, rec_ident1)

          for (rec_ident1,rec_ident2) in weight_vec_dict:
            assert rec_ident1 in self.rec_ident1
            assert rec_ident2 in self.rec_ident1
            assert isinstance(weight_vec_dict[(rec_ident1,rec_ident2)], list)

            assert len(weight_vec_dict[(rec_ident1,rec_ident2)]) == \
                   len(strmap_index.rec_comparator.field_comparison_list)

          # All record pairs from the blocking index should be in str-map index
          # as well
          #
          assert len(block_index_weight_vec_dict) <= len(weight_vec_dict), \
                 ('Number of record pairs in blocking index: %d, and in ' % \
                 (len(block_index_weight_vec_dict))+' canopy index: %d' % \
                 (len(weight_vec_dict)))

          # Build a string-map weight vector dict with pair identifiers sorted
          #
          for ((rec_ident1,rec_ident2),weight_list) in weight_vec_dict.items():
            if (rec_ident1 > rec_ident2):
              del weight_vec_dict[(rec_ident1,rec_ident2)]
              weight_vec_dict[(rec_ident2,rec_ident1)] = weight_list

          for (rec_ident1,rec_ident2) in block_index_weight_vec_dict:
            assert (rec_ident1,rec_ident2) in weight_vec_dict, \
                   'Record pair %s from blocking index not in str-map index' \
                   % (str((rec_ident1,rec_ident2)))

          # Due to heuristic nature we can't test each record pair if it was
          # in previous weight vector dictionary

#          for rec_ident_pair in prev_rec_dict.iterkeys():
#            assert rec_ident_pair in weight_vec_dict, \
#                   'Record pair %s not in weight vector dict with canopy ' % \
#                   (str(rec_ident_pair)) + 'method: %s' % (str(canopy_method))

          assert (len(prev_rec_dict) <= len(weight_vec_dict)), \
                 (len(prev_rec_dict), len(weight_vec_dict), canopy_method,
                 dim, sub_dim)

          prev_rec_dict = weight_vec_dict

    # Str-map index with threshold 1.0 / nearest 1 should be same as blocking
    #
    for canopy_method in [('threshold', 1.0, 1.0),
                          ('nearest', 1, 1)]:

      for (dim, sub_dim) in [(10,2),(10,3),(15,2),(15,3),(20,2),(20,3)]:

        strmap_index = indexing.StringMapIndex(descri = 'Test str-map index',
                                               dataset1 = self.dataset1,
                                               dataset2 = self.dataset1,
                                               rec_comp = self.rec_comp_dedupl,
                                               progress=4,
                                               canopy_me = canopy_method,
                                               dim = dim,
                                               log_funct = print_log,
                                               sub_d = sub_dim,
                                               grid_resol = 10,
                                               sim_funct = stringcmp.editdist,
                                               index_sep_str = '-',
                                               index_def = [index_def1,
                                               index_def2])
        strmap_index.build()
        strmap_index.compact()
        [field_names_list, weight_vec_dict] = strmap_index.run()

        # Make sure record ident. in each pair are sorted (smaller id first)
        #
        for ((rec_ident1,rec_ident2),w) in weight_vec_dict.items():
          if (rec_ident1 > rec_ident2):
            del weight_vec_dict[(rec_ident1,rec_ident2)]
            weight_vec_dict[(rec_ident2,rec_ident1)] = w

        for (rec_ident1,rec_ident2) in weight_vec_dict.keys():
          assert (rec_ident1,rec_ident2) in block_index_weight_vec_dict, \
                 'Not in blocking weight vectors: %s' % \
                 (str((rec_ident1,rec_ident2)))

        for (rec_ident1,rec_ident2) in block_index_weight_vec_dict.keys():
          assert (rec_ident1,rec_ident2) in weight_vec_dict, \
                 'Not in canopy weight vectors: %s' % \
                (str((rec_ident1,rec_ident2)))

    # Test length filter
    #
    [field_names_list, prev_w_vec_dict] = strmap_index.run()

    for lf in [99,50,30,20,15,10,5,1,0.5]:

      [field_names_list, this_w_vec_dict] = \
                                      strmap_index.run(length_filter_perc = lf)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

    # Test cut-off threshold
    #
    [field_names_list, prev_w_vec_dict] = strmap_index.run()

    for cot in [0.05, 0.1, 0.2, 0.5, 0.7, 0.8, 0.9]:

      [field_names_list, this_w_vec_dict] = \
                                      strmap_index.run(cut_off_threshold = cot)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

  # ---------------------------------------------------------------------------

  def testBigMatchIndexLinkage(self):  # - - - - - - - - - - - - - - - - - - -
    """Test BigMatchIndex linkage"""

    return

    index_def1 = [['surname','surname',False,False,None,[]]]
    index_def2 = [['given_name','given_name',True,True,4,[]],
                  ['postcode','postcode',True,False,2,[]]]

    # Assume blocking index is OK, so get the record pairs from this index for
    # comparison
    #
    block_index = indexing.BlockingIndex(description = 'Test blocking index',
                                         dataset1 = self.dataset1,
                                         dataset2 = self.dataset2,
                                         rec_comparator = self.rec_comp_link,
                                         progress=2,
                                         index_sep_str = ' ',
                                         index_def = [index_def1,index_def2])
    block_index.build()
    block_index.compact()
    [field_names_list, block_index_weight_vec_dict] = block_index.run()

    # Make sure record identifiers in each pair are sorted (smaller id first)
    #
    for ((rec_ident1,rec_ident2),w) in block_index_weight_vec_dict.items():
      if (rec_ident1 > rec_ident2):
        del block_index_weight_vec_dict[(rec_ident1,rec_ident2)]
        block_index_weight_vec_dict[(rec_ident2,rec_ident1)] = w

    for block_method_list in [[('block',)],
                              [('sort',1),('sort',2),('sort',3),('sort',4),
                               ('sort',5),('sort',6),('sort',8)],
                              [('qgram',1,True,1.0),('qgram',1,True,0.9),
                               ('qgram',1,True,0.8),('qgram',1,True,0.7),
                               ('qgram',1,True,0.6),('qgram',1,True,0.4)],
                              [('qgram',1,False,1.0),('qgram',1,False,0.9),
                               ('qgram',1,False,0.8),('qgram',1,False,0.7),
                               ('qgram',1,False,0.6),('qgram',1,False,0.4)],
                              [('qgram',2,True,1.0),('qgram',2,True,0.9),
                               ('qgram',2,True,0.8),('qgram',2,True,0.7),
                               ('qgram',2,True,0.6),('qgram',2,True,0.4)],
                              [('qgram',2,False,1.0),('qgram',2,False,0.9),
                               ('qgram',2,False,0.8),('qgram',2,False,0.7),
                               ('qgram',2,False,0.6),('qgram',2,False,0.4)],
                              [('qgram',3,True,1.0),('qgram',3,True,0.9),
                               ('qgram',3,True,0.8),('qgram',3,True,0.7),
                               ('qgram',3,True,0.6),('qgram',3,True,0.4)],
                              [('qgram',3,False,1.0),('qgram',3,False,0.9),
                               ('qgram',3,False,0.8),('qgram',3,False,0.7),
                               ('qgram',3,False,0.6),('qgram',3,False,0.4)]]:

      prev_rec_dict = {}

      for block_method in block_method_list:

        bigmatch_index= indexing.BigMatchIndex(descrip = 'Test BigMatch index',
                                               dataset1 = self.dataset1,
                                               dataset2 = self.dataset2,
                                               block_method = block_method,
                                               rec_compar = self.rec_comp_link,
                                               progress=2,
                                               index_sep_str = ' ',
                                               index_d = [index_def1,
                                                          index_def2])

        assert isinstance(bigmatch_index.index1, dict)
        assert isinstance(bigmatch_index.index2, dict)
        assert isinstance(bigmatch_index.description, str)
        assert isinstance(bigmatch_index.rec_cache1, dict)
        assert isinstance(bigmatch_index.rec_cache2, dict)
        assert bigmatch_index.progress_report == 2
        assert bigmatch_index.skip_missing == True
        assert bigmatch_index.do_deduplication == False
        assert len(bigmatch_index.index_def) == 2
        assert bigmatch_index.index_sep_str == ' '
        assert bigmatch_index.num_rec_pairs == None

        assert bigmatch_index.status == 'initialised'

        bigmatch_index.build()  # - - - - - - - - - - - - - - - - - - - - - - -

        assert bigmatch_index.num_rec_pairs == None
        assert isinstance(bigmatch_index.index1, dict)
        assert isinstance(bigmatch_index.index2, dict)
        assert isinstance(bigmatch_index.description, str)
        assert isinstance(bigmatch_index.rec_cache1, dict)
        assert isinstance(bigmatch_index.rec_cache2, dict)
        assert bigmatch_index.progress_report == 2
        assert bigmatch_index.skip_missing == True
        assert bigmatch_index.do_deduplication == False
        assert len(bigmatch_index.index_def) == 2
        assert bigmatch_index.index_sep_str == ' '

        assert bigmatch_index.status == 'built'

        bigmatch_index.compact()  # - - - - - - - - - - - - - - - - - - - - - -

        assert bigmatch_index.num_rec_pairs == None
        assert isinstance(bigmatch_index.index1, dict)
        assert isinstance(bigmatch_index.index2, dict)
        assert isinstance(bigmatch_index.description, str)
        assert isinstance(bigmatch_index.rec_cache1, dict)
        assert isinstance(bigmatch_index.rec_cache2, dict)
        assert bigmatch_index.progress_report == 2
        assert bigmatch_index.skip_missing == True
        assert bigmatch_index.do_deduplication == False
        assert len(bigmatch_index.index_def) == 2
        assert bigmatch_index.index_sep_str == ' '
        assert bigmatch_index.status == 'compacted'

        [field_names_list, weight_vec_dict] = bigmatch_index.run()  # - - - - -

        assert isinstance(weight_vec_dict, dict)
        assert len(weight_vec_dict) == bigmatch_index.num_rec_pairs, \
               (len(weight_vec_dict), bigmatch_index.num_rec_pairs)
        assert bigmatch_index.num_rec_pairs > 0
        assert bigmatch_index.num_rec_pairs <= \
               len(self.rec_ident1)*len(self.rec_ident2)

        for rec_ident1 in self.rec_ident1:
          for rec_ident2 in self.rec_ident2:
            if (rec_ident1 == rec_ident2):
              assert (rec_ident1,rec_ident2) in weight_vec_dict

        for (rec_ident1,rec_ident2) in weight_vec_dict:
          assert rec_ident1 in self.rec_ident1
          assert rec_ident2 in self.rec_ident2
          assert isinstance(weight_vec_dict[(rec_ident1,rec_ident2)], list)

          assert len(weight_vec_dict[(rec_ident1,rec_ident2)]) == \
                 len(bigmatch_index.rec_comparator.field_comparison_list)

        # All record pairs from the blocking index should be in the BigMatch
        # map index as well
        #
        assert len(block_index_weight_vec_dict) <= len(weight_vec_dict), \
               ('Number of record pairs in blocking index: %d, and in ' % \
               (len(block_index_weight_vec_dict))+' BigMatch index: %d' % \
               (len(weight_vec_dict)))

        # Build a BigMatch weight vector dict with pair identifiers sorted
        #
        for ((rec_ident1,rec_ident2),weight_list) in weight_vec_dict.items():
          if (rec_ident1 > rec_ident2):
            del weight_vec_dict[(rec_ident1,rec_ident2)]
            weight_vec_dict[(rec_ident2,rec_ident1)] = weight_list

        for (rec_ident1,rec_ident2) in block_index_weight_vec_dict:
          assert (rec_ident1,rec_ident2) in weight_vec_dict, \
                 'Record pair %s from blocking index not in BigMatch index' % \
                 (str((rec_ident1,rec_ident2)))+' with block method: %s' % \
                 (block_method)

        for rec_ident_pair in prev_rec_dict.iterkeys():
          assert rec_ident_pair in weight_vec_dict, \
                 'Record pair %s not in previous weight vector dict: %s'  % \
                 (str(rec_ident_pair), str(block_method))

        prev_rec_dict = weight_vec_dict

    # BigMatch index with method block or sort window size 1 or q-gram
    # threshold 1 should be same as blocking
    #
    for block_method in [('block',), ('sort',1),
                         ('qgram',1,True,1.0),('qgram',1,False,1.0),
                         ('qgram',2,True,1.0),('qgram',2,False,1.0),
                         ('qgram',3,True,1.0),('qgram',3,False,1.0)]:

      bigmatch_index = indexing.BigMatchIndex(descrip = 'Test BigMatch index',
                                              dataset1 = self.dataset1,
                                              dataset2 = self.dataset2,
                                              block_method = block_method,
                                              rec_compar = self.rec_comp_link,
                                              log_funct = print_log,
                                              progress=2,
                                              index_sep_str = ' ',
                                              index_d = [index_def1,
                                                         index_def2])
      bigmatch_index.build()
      bigmatch_index.compact()
      [field_names_list, weight_vec_dict] = bigmatch_index.run()

      # Make sure record ident. in each pair are sorted (smaller id first)
      #
      for ((rec_ident1,rec_ident2),w) in weight_vec_dict.items():
        if (rec_ident1 > rec_ident2):
          del weight_vec_dict[(rec_ident1,rec_ident2)]
          weight_vec_dict[(rec_ident2,rec_ident1)] = w

      for (rec_ident1,rec_ident2) in weight_vec_dict.keys():
        assert (rec_ident1,rec_ident2) in block_index_weight_vec_dict, \
               'Not in blocking weight vectors: %s' % \
               (str((rec_ident1,rec_ident2)))

      for (rec_ident1,rec_ident2) in block_index_weight_vec_dict.keys():
        assert (rec_ident1,rec_ident2) in weight_vec_dict, \
               'Not in BigMatch weight vectors: %s' % \
              (str((rec_ident1,rec_ident2)))

    # Test length filter
    #
    [field_names_list, prev_w_vec_dict] = bigmatch_index.run()

    for lf in [99,50,30,20,15,10,5,1,0.5]:

      [field_names_list, this_w_vec_dict] = \
                                    bigmatch_index.run(length_filter_perc = lf)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

    # Test cut-off threshold
    #
    [field_names_list, prev_w_vec_dict] = bigmatch_index.run()

    for cot in [0.05, 0.1, 0.2, 0.5, 0.7, 0.8, 0.9]:

      [field_names_list, this_w_vec_dict] = \
                                    bigmatch_index.run(cut_off_threshold = cot)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

  # ---------------------------------------------------------------------------

  def testDedupIndexLinkage(self):  # - - - - - - - - - - - - - - - - - - - - -
    """Test DedupIndex deduplication"""

    return

    index_def1 = [['surname','surname',False,False,None,[]]]
    index_def2 = [['given_name','given_name',True,True,4,[]],
                  ['postcode','postcode',True,False,2,[]]]

    # Assume blocking index is OK, so get the record pairs from this index for
    # comparison
    #
    block_index = indexing.BlockingIndex(description = 'Test blocking index',
                                         dataset1 = self.dataset1,
                                         dataset2 = self.dataset1,
                                         rec_comparator = self.rec_comp_dedupl,
                                         progress=2,
                                         index_sep_str = chr(1),
                                         index_def = [index_def1,index_def2])
    block_index.build()
    block_index.compact()
    [field_names_list, block_index_weight_vec_dict] = block_index.run()

    # Make sure record identifiers in each pair are sorted (smaller id first)
    #
    for ((rec_ident1,rec_ident2),w) in block_index_weight_vec_dict.items():
      if (rec_ident1 > rec_ident2):
        del block_index_weight_vec_dict[(rec_ident1,rec_ident2)]
        block_index_weight_vec_dict[(rec_ident2,rec_ident1)] = w

    for block_method_list in [[('block',)],
                              [('sort',1),('sort',2),('sort',3),('sort',4),
                               ('sort',5),('sort',6),('sort',8)],
                              [('qgram',1,True,1.0),('qgram',1,True,0.9),
                               ('qgram',1,True,0.8),('qgram',1,True,0.7),
                               ('qgram',1,True,0.6),('qgram',1,True,0.4)],
                              [('qgram',1,False,1.0),('qgram',1,False,0.9),
                               ('qgram',1,False,0.8),('qgram',1,False,0.7),
                               ('qgram',1,False,0.6),('qgram',1,False,0.4)],
                              [('qgram',2,True,1.0),('qgram',2,True,0.9),
                               ('qgram',2,True,0.8),('qgram',2,True,0.7),
                               ('qgram',2,True,0.6),('qgram',2,True,0.4)],
                              [('qgram',2,False,1.0),('qgram',2,False,0.9),
                               ('qgram',2,False,0.8),('qgram',2,False,0.7),
                               ('qgram',2,False,0.6),('qgram',2,False,0.4)],
                              [('qgram',3,True,1.0),('qgram',3,True,0.9),
                               ('qgram',3,True,0.8),('qgram',3,True,0.7),
                               ('qgram',3,True,0.6),('qgram',3,True,0.4)],
                              [('qgram',3,False,1.0),('qgram',3,False,0.9),
                               ('qgram',3,False,0.8),('qgram',3,False,0.7),
                               ('qgram',3,False,0.6),('qgram',3,False,0.4)]]:
      prev_rec_dict = {}

      for block_method in block_method_list:

        dedup_index= indexing.DedupIndex(descrip = 'Test Dedup index',
                                               dataset1 = self.dataset1,
                                               dataset2 = self.dataset1,
                                               block_method = block_method,
                                               rec_comp = self.rec_comp_dedupl,
                                               progress=2,
                                               index_sep_str = chr(1),
                                               index_d = [index_def1,
                                                          index_def2])

        assert isinstance(dedup_index.index1, dict)
        assert isinstance(dedup_index.index2, dict)
        assert isinstance(dedup_index.description, str)
        assert isinstance(dedup_index.rec_cache1, dict)
        assert isinstance(dedup_index.rec_cache2, dict)
        assert dedup_index.progress_report == 2
        assert dedup_index.skip_missing == True
        assert dedup_index.do_deduplication == True
        assert len(dedup_index.index_def) == 2
        assert dedup_index.index_sep_str == chr(1)
        assert dedup_index.num_rec_pairs == None

        assert dedup_index.status == 'initialised'

        dedup_index.build()  # - - - - - - - - - - - - - - - - - - - - - - -

        assert dedup_index.num_rec_pairs == None
        assert isinstance(dedup_index.index1, dict)
        assert isinstance(dedup_index.index2, dict)
        assert isinstance(dedup_index.description, str)
        assert isinstance(dedup_index.rec_cache1, dict)
        assert isinstance(dedup_index.rec_cache2, dict)
        assert dedup_index.progress_report == 2
        assert dedup_index.skip_missing == True
        assert dedup_index.do_deduplication == True
        assert len(dedup_index.index_def) == 2
        assert dedup_index.index_sep_str == chr(1)

        assert dedup_index.status == 'built'

        dedup_index.compact()  # - - - - - - - - - - - - - - - - - - - - - -

        assert dedup_index.num_rec_pairs == None
        assert isinstance(dedup_index.index1, dict)
        assert isinstance(dedup_index.index2, dict)
        assert isinstance(dedup_index.description, str)
        assert isinstance(dedup_index.rec_cache1, dict)
        assert isinstance(dedup_index.rec_cache2, dict)
        assert dedup_index.progress_report == 2
        assert dedup_index.skip_missing == True
        assert dedup_index.do_deduplication == True
        assert len(dedup_index.index_def) == 2
        assert dedup_index.index_sep_str == chr(1)
        assert dedup_index.status == 'compacted'

        [field_names_list, weight_vec_dict] = dedup_index.run()  # - - - - - -

        assert isinstance(weight_vec_dict, dict)
        assert len(weight_vec_dict) == dedup_index.num_rec_pairs, \
               (len(weight_vec_dict), dedup_index.num_rec_pairs)
        assert dedup_index.num_rec_pairs > 0
        assert dedup_index.num_rec_pairs <= \
               len(self.rec_ident1)*(len(self.rec_ident1)-1)/2

        for rec_ident1 in self.rec_ident1:
          assert (rec_ident1, rec_ident1) not in weight_vec_dict, \
                 (rec_ident1, rec_ident1)

        for (rec_ident1,rec_ident2) in weight_vec_dict:
          assert rec_ident1 in self.rec_ident1
          assert rec_ident2 in self.rec_ident2
          assert isinstance(weight_vec_dict[(rec_ident1,rec_ident2)], list)

          assert len(weight_vec_dict[(rec_ident1,rec_ident2)]) == \
                 len(dedup_index.rec_comparator.field_comparison_list)

        # All record pairs from the blocking index should be in the Dedup index
        # as well
        #
        assert len(block_index_weight_vec_dict) <= len(weight_vec_dict), \
               ('Number of record pairs in blocking index: %d, and in ' % \
               (len(block_index_weight_vec_dict))+' Dedup index: %d' % \
               (len(weight_vec_dict)))

        # Build a Dedup weight vector dict with pair identifiers sorted
        #
        for ((rec_ident1,rec_ident2),weight_list) in weight_vec_dict.items():
          if (rec_ident1 > rec_ident2):
            del weight_vec_dict[(rec_ident1,rec_ident2)]
            weight_vec_dict[(rec_ident2,rec_ident1)] = weight_list

        for (rec_ident1,rec_ident2) in block_index_weight_vec_dict:
          assert (rec_ident1,rec_ident2) in weight_vec_dict, \
                 'Record pair %s from blocking index not in Dedup index' % \
                 (str((rec_ident1,rec_ident2)))+' with block method: %s' % \
                 (block_method)

        for rec_ident_pair in prev_rec_dict.iterkeys():
          assert rec_ident_pair in weight_vec_dict, \
                 'Record pair %s not in previous weight vector dict: %s'  % \
                 (str(rec_ident_pair), str(block_method))

        prev_rec_dict = weight_vec_dict

    # Dedup index with method block or q-gram threshold 1 should be same as
    # blocking
    #
    for block_method in [('block',), ##('sort',1),
                         ('qgram',1,True,1.0),('qgram',1,False,1.0),
                         ('qgram',2,True,1.0),('qgram',2,False,1.0),
                         ('qgram',3,True,1.0),('qgram',3,False,1.0)]:

      dedup_index = indexing.DedupIndex(descrip = 'Test Dedup index',
                                              dataset1 = self.dataset1,
                                              dataset2 = self.dataset1,
                                              block_method = block_method,
                                              rec_compa = self.rec_comp_dedupl,
                                              progress=2,
                                              log_funct = print_log,
                                              index_sep_str = chr(1),
                                              index_d = [index_def1,
                                                         index_def2])
      dedup_index.build()
      dedup_index.compact()
      [field_names_list, weight_vec_dict] = dedup_index.run()

      # Make sure record ident. in each pair are sorted (smaller id first)
      #
      for ((rec_ident1,rec_ident2),w) in weight_vec_dict.items():
        if (rec_ident1 > rec_ident2):
          del weight_vec_dict[(rec_ident1,rec_ident2)]
          weight_vec_dict[(rec_ident2,rec_ident1)] = w

      for (rec_ident1,rec_ident2) in weight_vec_dict.keys():
        assert (rec_ident1,rec_ident2) in block_index_weight_vec_dict, \
               'Not in blocking weight vectors: %s' % \
               (str((rec_ident1,rec_ident2)))

      for (rec_ident1,rec_ident2) in block_index_weight_vec_dict.keys():
        assert (rec_ident1,rec_ident2) in weight_vec_dict, \
               'Not in Dedup weight vectors: %s' % \
              (str((rec_ident1,rec_ident2)))

    # Test length filter
    #
    [field_names_list, prev_w_vec_dict] = dedup_index.run()

    for lf in [99,50,30,20,15,10,5,1,0.5]:

      [field_names_list, this_w_vec_dict] = \
                                       dedup_index.run(length_filter_perc = lf)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

    # Test cut-off threshold
    #
    [field_names_list, prev_w_vec_dict] = dedup_index.run()

    for cot in [0.05, 0.1, 0.2, 0.5, 0.7, 0.8, 0.9]:

      [field_names_list, this_w_vec_dict] = \
                                       dedup_index.run(cut_off_threshold = cot)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

  # ---------------------------------------------------------------------------

  def testSuffixArrayIndexLinkage(self):  # - - - - - - - - - - - - - - - - - -
    """Test SuffixArrayIndex linkage"""

    return

    index_def1 = [['surname','surname',False,False,None,[]]]
    index_def2 = [['given_name','given_name',True,True,4,[]],
                  ['postcode','postcode',True,False,2,[]]]

    # Assume blocking index is OK, so get the record pairs from this index for
    # comparison
    #
    block_index = indexing.BlockingIndex(description = 'Test blocking index',
                                         dataset1 = self.dataset1,
                                         dataset2 = self.dataset2,
                                         rec_comparator = self.rec_comp_link,
                                         progress=2,
                                         index_def = [index_def1,index_def2])
    block_index.build()
    block_index.compact()
    [field_names_list, block_index_weight_vec_dict] = block_index.run()

    # Make sure record identifiers in each pair are sorted (smaller id first)
    #
    for ((rec_ident1,rec_ident2),w) in block_index_weight_vec_dict.items():
      if (rec_ident1 > rec_ident2):
        del block_index_weight_vec_dict[(rec_ident1,rec_ident2)]
        block_index_weight_vec_dict[(rec_ident2,rec_ident1)] = w

    block_method_lists = [[(5,5,False),(4,5,False),(3,5,False),
                           (2,5,False),(1,5,False)],
                          [(5,5,True),(4,5,True),(3,5,True),
                           (2,5,True),(1,5,True)],
                          [(5,10,False),(4,10,False),(3,10,False),
                           (2,10,False),(1,10,False)],
                          [(5,10,True),(4,10,True),(3,10,True),
                           (2,10,True),(1,10,True)]]

    for block_method_list in block_method_lists:
      for (min_q_gram_len,max_block_size, padded) in block_method_list:

        for suff_method in ['suffixonly', 'allsubstr']:

          prev_rec_dict = {}

          sarray_index = indexing.SuffixArrayIndex(desc='Test suff-arr index',
                                                 dataset1 = self.dataset1,
                                                 dataset2 = self.dataset2,
                                                 rec_compar=self.rec_comp_link,
                                                 progress=2,
                                                 suffix_m = suff_method,
                                                 padd = padded,
                                                 block_method =(min_q_gram_len,
                                                               max_block_size),
                                                 index_def = [index_def1,
                                                              index_def2])
          assert isinstance(sarray_index.index1, dict)
          assert isinstance(sarray_index.index2, dict)
          assert isinstance(sarray_index.description, str)
          assert isinstance(sarray_index.rec_cache1, dict)
          assert isinstance(sarray_index.rec_cache2, dict)
          assert sarray_index.progress_report == 2
          assert sarray_index.skip_missing == True
          assert sarray_index.do_deduplication == False
          assert len(sarray_index.index_def) == 2
          assert sarray_index.padded == padded
          assert sarray_index.suffix_method == suff_method
          assert sarray_index.block_method == (min_q_gram_len, max_block_size)
          assert sarray_index.num_rec_pairs == None

          assert sarray_index.status == 'initialised'

          sarray_index.build()  # - - - - - - - - - - - - - - - - - - - - - - -

          assert sarray_index.num_rec_pairs == None
          assert isinstance(sarray_index.index1, dict)
          assert isinstance(sarray_index.index2, dict)
          assert isinstance(sarray_index.description, str)
          assert isinstance(sarray_index.rec_cache1, dict)
          assert isinstance(sarray_index.rec_cache2, dict)
          assert sarray_index.progress_report == 2
          assert sarray_index.skip_missing == True
          assert sarray_index.do_deduplication == False
          assert len(sarray_index.index_def) == 2
          assert sarray_index.padded == padded
          assert sarray_index.suffix_method == suff_method
          assert sarray_index.block_method == (min_q_gram_len, max_block_size)
          assert sarray_index.status == 'built'

          sarray_index.compact()  # - - - - - - - - - - - - - - - - - - - - - -

          assert sarray_index.num_rec_pairs > 0
          assert sarray_index.num_rec_pairs <= \
                 len(self.rec_ident1)*len(self.rec_ident2)
          assert isinstance(sarray_index.index1, dict)
          assert isinstance(sarray_index.index2, dict)
          assert isinstance(sarray_index.description, str)
          assert isinstance(sarray_index.rec_cache1, dict)
          assert isinstance(sarray_index.rec_cache2, dict)
          assert sarray_index.progress_report == 2
          assert sarray_index.skip_missing == True
          assert sarray_index.do_deduplication == False
          assert len(sarray_index.index_def) == 2
          assert sarray_index.padded == padded
          assert sarray_index.suffix_method == suff_method
          assert sarray_index.block_method == (min_q_gram_len, max_block_size)
          assert sarray_index.status == 'compacted'

          [field_names_list, weight_vec_dict] = sarray_index.run()  # - - - - -

          assert isinstance(weight_vec_dict, dict)
          assert len(weight_vec_dict) == sarray_index.num_rec_pairs
          assert sarray_index.num_rec_pairs > 0
          assert sarray_index.num_rec_pairs <= \
                 len(self.rec_ident1)*len(self.rec_ident2)

          for rec_ident1 in self.rec_ident1:
            for rec_ident2 in self.rec_ident2:
              if (rec_ident1 == rec_ident2):
                assert (rec_ident1,rec_ident2) in weight_vec_dict, \
                       (rec_ident1,rec_ident2)

          for (rec_ident1,rec_ident2) in weight_vec_dict:
            assert rec_ident1 in self.rec_ident1
            assert rec_ident2 in self.rec_ident2
            assert isinstance(weight_vec_dict[(rec_ident1,rec_ident2)], list)

            assert len(weight_vec_dict[(rec_ident1,rec_ident2)]) == \
                   len(sarray_index.rec_comparator.field_comparison_list)

          # All record pairs from the blocking index should be in the suffix
          # array index as well (only when min_g_gram_len = 1)
          #
          if (min_q_gram_len == 1):

            assert len(block_index_weight_vec_dict) <= len(weight_vec_dict), \
                 ('Number of record pairs in blocking index: %d, and in ' % \
                 (len(block_index_weight_vec_dict))+' suffix array index: %d' \
                 % (len(weight_vec_dict)))

            # Build suffix array weight vector dict with pair identifier sorted
            #
            for ((rec_ident1,rec_ident2),w_list) in weight_vec_dict.items():
              if (rec_ident1 > rec_ident2):
                del weight_vec_dict[(rec_ident1,rec_ident2)]
                weight_vec_dict[(rec_ident2,rec_ident1)] = w_list

            for (rec_ident1,rec_ident2) in block_index_weight_vec_dict:
              assert (rec_ident1,rec_ident2) in weight_vec_dict, \
                   'Record pair %s from blocking index not in suffix array' % \
                   (str((rec_ident1,rec_ident2)))+' index'

            for rec_ident_pair in prev_rec_dict.iterkeys():
              assert rec_ident_pair in weight_vec_dict, \
                 'Record pair %s not in weight vector dict with method: %s' % \
                 (str(rec_ident_pair), (str(sarray_index.block_method)))

          prev_rec_dict = weight_vec_dict

    # Suffix array index with min_q_gram_len 1 and max_block_size very large
    # should contain all pairs that are in blocking index
    #
    sarray_index = indexing.SuffixArrayIndex(desc='Test suffix array index',
                                             dataset1 = self.dataset1,
                                             dataset2 = self.dataset2,
                                             rec_compar = self.rec_comp_link,
                                             progress=2,
                                             padd = False,
                                             log_funct = print_log,
                                             suffix_method = 'suffixonly',
                                             block_method = (1,999),
                                             index_def = [index_def1,
                                                          index_def2])
    sarray_index.build()
    sarray_index.compact()
    [field_names_list, weight_vec_dict] = sarray_index.run()

    # Make sure record identifiers in each pair are sorted (smaller id first)
    #
    for ((rec_ident1,rec_ident2),w) in weight_vec_dict.items():
      if (rec_ident1 > rec_ident2):
        del weight_vec_dict[(rec_ident1,rec_ident2)]
        weight_vec_dict[(rec_ident2,rec_ident1)] = w

#    for (rec_ident1,rec_ident2) in weight_vec_dict.keys():
#      assert (rec_ident1,rec_ident2) in block_index_weight_vec_dict, \
#             'Not in blocking weight vectors: %s' % \
#             (str((rec_ident1,rec_ident2)))

    for (rec_ident1,rec_ident2) in block_index_weight_vec_dict.keys():
      assert (rec_ident1,rec_ident2) in weight_vec_dict, \
             'Not in suffix array weight vectors: %s' % \
             (str((rec_ident1,rec_ident2)))

    # Test length filter
    #
    [field_names_list, prev_w_vec_dict] = sarray_index.run()

    for lf in [99,50,30,20,15,10,5,1,0.5]:

      [field_names_list, this_w_vec_dict] = \
                                      sarray_index.run(length_filter_perc = lf)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

    # Test cut-off threshold
    #
    [field_names_list, prev_w_vec_dict] = sarray_index.run()

    for cot in [0.05, 0.1, 0.2, 0.5, 0.7, 0.8, 0.9]:

      [field_names_list, this_w_vec_dict] = \
                                      sarray_index.run(cut_off_threshold = cot)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

  def testSuffixArrayIndexDedupl(self):  # - - - - - - - - - - - - - - - - - - -
    """Test SuffixArrayIndex deduplication"""

    return

    index_def1 = [['surname','surname',False,False,None,[]]]
    index_def2 = [['given_name','given_name',True,True,4,[]],
                  ['postcode','postcode',True,False,2,[]]]

    # Assume blocking index is OK, so get the record pairs from this index for
    # comparison
    #
    block_index = indexing.BlockingIndex(description = 'Test blocking index',
                                         dataset1 = self.dataset1,
                                         dataset2 = self.dataset1,
                                         rec_comparator = self.rec_comp_dedupl,
                                         progress=4,
                                         index_def = [index_def1,index_def2])
    block_index.build()
    block_index.compact()
    [field_names_list, block_index_weight_vec_dict] = block_index.run()

    # Make sure record identifiers in each pair are sorted (smaller id first)
    #
    for ((rec_ident1,rec_ident2),w) in block_index_weight_vec_dict.items():
      if (rec_ident1 > rec_ident2):
        del block_index_weight_vec_dict[(rec_ident1,rec_ident2)]
        block_index_weight_vec_dict[(rec_ident2,rec_ident1)] = w

    block_method_lists = [[(5,5,False),(4,5,False),(3,5,False),
                           (2,5,False),(1,5,False)],
                          [(5,5,True),(4,5,True),(3,5,True),
                           (2,5,True),(1,5,True)],
                          [(5,10,False),(4,10,False),(3,10,False),
                           (2,10,False),(1,10,False)],
                          [(5,10,True),(4,10,True),(3,10,True),
                           (2,10,True),(1,10,True)]]

    for block_method_list in block_method_lists:
      for (min_q_gram_len, max_block_size, padded) in block_method_list:

        for suff_method in ['suffixonly', 'allsubstr']:

          prev_rec_dict = {}

          sarray_index = indexing.SuffixArrayIndex(desc='Test suff-arr index',
                                               dataset1 = self.dataset1,
                                               dataset2 = self.dataset1,
                                               rec_comp = self.rec_comp_dedupl,
                                               progress=4,
                                               suffix_m = suff_method,
                                               padd = padded,
                                               skip_m = False,
                                               block_method = (min_q_gram_len,
                                                               max_block_size),
                                               index_def = [index_def1,
                                                            index_def2])

          assert isinstance(sarray_index.index1, dict)
          assert isinstance(sarray_index.index2, dict)
          assert isinstance(sarray_index.description, str)
          assert isinstance(sarray_index.rec_cache1, dict)
          assert isinstance(sarray_index.rec_cache2, dict)
          assert sarray_index.progress_report == 4
          assert sarray_index.skip_missing == False
          assert sarray_index.do_deduplication == True
          assert len(sarray_index.index_def) == 2
          assert sarray_index.padded == padded
          assert sarray_index.suffix_method == suff_method
          assert sarray_index.block_method == (min_q_gram_len, max_block_size)
          assert sarray_index.num_rec_pairs == None

          assert sarray_index.status == 'initialised'

          sarray_index.build()  # - - - - - - - - - - - - - - - - - - - - - - -

          assert sarray_index.num_rec_pairs == None
          assert isinstance(sarray_index.index1, dict)
          assert isinstance(sarray_index.index2, dict)
          assert isinstance(sarray_index.description, str)
          assert isinstance(sarray_index.rec_cache1, dict)
          assert isinstance(sarray_index.rec_cache2, dict)
          assert sarray_index.progress_report == 4
          assert sarray_index.skip_missing == False
          assert sarray_index.do_deduplication == True
          assert len(sarray_index.index_def) == 2
          assert sarray_index.padded == padded
          assert sarray_index.suffix_method == suff_method
          assert sarray_index.block_method == (min_q_gram_len, max_block_size)

          assert sarray_index.status == 'built'

          sarray_index.compact()  # - - - - - - - - - - - - - - - - - - - - - -

          assert sarray_index.num_rec_pairs > 0
          assert sarray_index.num_rec_pairs <= \
                 len(self.rec_ident1)*len(self.rec_ident2)
          assert isinstance(sarray_index.index1, dict)
          assert isinstance(sarray_index.index2, dict)
          assert isinstance(sarray_index.description, str)
          assert isinstance(sarray_index.rec_cache1, dict)
          assert isinstance(sarray_index.rec_cache2, dict)
          assert sarray_index.progress_report == 4
          assert sarray_index.skip_missing == False
          assert sarray_index.do_deduplication == True
          assert len(sarray_index.index_def) == 2
          assert sarray_index.suffix_method == suff_method
          assert sarray_index.padded == padded
          assert sarray_index.block_method == (min_q_gram_len, max_block_size)

          assert sarray_index.status == 'compacted'

          [field_names_list, weight_vec_dict] = sarray_index.run()  # - - - - -

          assert isinstance(weight_vec_dict, dict)
          assert len(weight_vec_dict) == sarray_index.num_rec_pairs
          assert sarray_index.num_rec_pairs > 0
          assert sarray_index.num_rec_pairs <= \
                 len(self.rec_ident1)*len(self.rec_ident2)

          for rec_ident1 in self.rec_ident1:
            assert (rec_ident1, rec_ident1) not in weight_vec_dict, \
                   (rec_ident1, rec_ident1)

          for (rec_ident1,rec_ident2) in weight_vec_dict:
            assert rec_ident1 in self.rec_ident1
            assert rec_ident2 in self.rec_ident1
            assert isinstance(weight_vec_dict[(rec_ident1,rec_ident2)], list)

            assert len(weight_vec_dict[(rec_ident1,rec_ident2)]) == \
                   len(sarray_index.rec_comparator.field_comparison_list)

          # All record pairs from the blocking index should be in the suffix
          # array index as well (only when min_g_gram_len = 1)
          #
          if (min_q_gram_len == 1):

            assert len(block_index_weight_vec_dict) <= len(weight_vec_dict), \
                 ('Number of record pairs in blocking index: %d, and in ' % \
                 (len(block_index_weight_vec_dict))+' suffix array index: %d' \
                 % (len(weight_vec_dict)))

            # Build suffix array weight vector dict with pair identifier sorted
            #
            for ((rec_ident1,rec_ident2),w_list) in weight_vec_dict.items():
              if (rec_ident1 > rec_ident2):
                del weight_vec_dict[(rec_ident1,rec_ident2)]
                weight_vec_dict[(rec_ident2,rec_ident1)] = w_list

            for (rec_ident1,rec_ident2) in block_index_weight_vec_dict:
              assert (rec_ident1,rec_ident2) in weight_vec_dict, \
                   'Record pair %s from blocking index not in suffix array ' \
                    % (str((rec_ident1,rec_ident2)))+'index'

            for rec_ident_pair in prev_rec_dict.iterkeys():
              assert rec_ident_pair in weight_vec_dict, \
                   'Record pair %s not in weight vector dict with method: %s' \
                   % (str(rec_ident_pair), (str(sarray_index.block_method)))

          prev_rec_dict = weight_vec_dict

    # Suffix array index with min_q_gram_len 1 and max_block_size very large
    # should contain all pairs that are in blocking index
    #
    sarray_index = indexing.SuffixArrayIndex(desc='Test suffix array index',
                                             dataset1 = self.dataset1,
                                             dataset2 = self.dataset1,
                                             rec_compar = self.rec_comp_dedupl,
                                             progress=4,
                                             log_funct = print_log,
                                             skip_m = False,
                                             padd = False,
                                             suffix_method = 'suffixonly',
                                             block_method = (1,999),
                                             index_def = [index_def1,
                                                          index_def2])

    sarray_index.build()
    sarray_index.compact()
    [field_names_list, weight_vec_dict] = sarray_index.run()

    # Make sure record identifiers in each pair are sorted (smaller id first)
    #
    for ((rec_ident1,rec_ident2),w) in weight_vec_dict.items():
      if (rec_ident1 > rec_ident2):
        del weight_vec_dict[(rec_ident1,rec_ident2)]
        weight_vec_dict[(rec_ident2,rec_ident1)] = w

#    for (rec_ident1,rec_ident2) in weight_vec_dict.keys():
#      assert (rec_ident1,rec_ident2) in block_index_weight_vec_dict, \
#             'Not in blocking weight vectors: %s' % \
#             (str((rec_ident1,rec_ident2)))

    for (rec_ident1,rec_ident2) in block_index_weight_vec_dict.keys():
      assert (rec_ident1,rec_ident2) in weight_vec_dict, \
             'Not in suffix array weight vectors: %s' % \
             (str((rec_ident1,rec_ident2)))

    # Test length filter
    #
    [field_names_list, prev_w_vec_dict] = sarray_index.run()

    for lf in [99,50,30,20,15,10,5,1,0.5]:

      [field_names_list, this_w_vec_dict] = \
                                      sarray_index.run(length_filter_perc = lf)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

    # Test cut-off threshold
    #
    [field_names_list, prev_w_vec_dict] = sarray_index.run()

    for cot in [0.05, 0.1, 0.2, 0.5, 0.7, 0.8, 0.9]:

      [field_names_list, this_w_vec_dict] = \
                                      sarray_index.run(cut_off_threshold = cot)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict


  # ---------------------------------------------------------------------------

  def testRobustSuffixArrayIndexLinkage(self):  # - - - - - - - - - - - - - - -
    """Test RobustSuffixArrayIndex linkage"""

    return

    index_def1 = [['surname','surname',False,False,None,[]]]
    index_def2 = [['given_name','given_name',True,True,4,[]],
                  ['postcode','postcode',True,False,2,[]]]

    # Assume blocking index is OK, so get the record pairs from this index for
    # comparison
    #
    block_index = indexing.BlockingIndex(description = 'Test blocking index',
                                         dataset1 = self.dataset1,
                                         dataset2 = self.dataset2,
                                         rec_comparator = self.rec_comp_link,
                                         progress=2,
                                         index_def = [index_def1,index_def2])
    block_index.build()
    block_index.compact()
    [field_names_list, block_index_weight_vec_dict] = block_index.run()

    # Make sure record identifiers in each pair are sorted (smaller id first)
    #
    for ((rec_ident1,rec_ident2),w) in block_index_weight_vec_dict.items():
      if (rec_ident1 > rec_ident2):
        del block_index_weight_vec_dict[(rec_ident1,rec_ident2)]
        block_index_weight_vec_dict[(rec_ident2,rec_ident1)] = w

    block_method_lists = [[(5,5,False),(4,5,False),(3,5,False),
                           (2,5,False),(1,5,False)],
                          [(5,5,True),(4,5,True),(3,5,True),
                           (2,5,True),(1,5,True)],
                          [(5,10,False),(4,10,False),(3,10,False),
                           (2,10,False),(1,10,False)],
                          [(5,10,True),(4,10,True),(3,10,True),
                           (2,10,True),(1,10,True)]]

    for block_method_list in block_method_lists:
      for (min_q_gram_len,max_block_size, padded) in block_method_list:

        for (str_cmp_funct, str_cmp_thres) in \
            [(stringcmp.jaro,0.7), (stringcmp.jaro,0.8), (stringcmp.jaro,0.9),
             (stringcmp.bigram,0.7),(stringcmp.bigram,0.8),
             (stringcmp.bigram,0.9),
             (stringcmp.lcs, 0.7), (stringcmp.lcs, 0.8), (stringcmp.lcs, 0.9)]:

          prev_rec_dict = {}

          rsarray_index = indexing.RobustSuffixArrayIndex(desc = \
                                            'Test robust suff-arr index',
                                                 dataset1 = self.dataset1,
                                                 dataset2 = self.dataset2,
                                                 rec_compar=self.rec_comp_link,
                                                 progress=2,
                                                 str_cmp_funct = str_cmp_funct,
                                                 str_cmp_thres = str_cmp_thres,
                                                 skip_m = True,
                                                 padd = padded,
                                                 block_method =(min_q_gram_len,
                                                               max_block_size),
                                                 index_def = [index_def1,
                                                              index_def2])
          assert isinstance(rsarray_index.index1, dict)
          assert isinstance(rsarray_index.index2, dict)
          assert isinstance(rsarray_index.description, str)
          assert isinstance(rsarray_index.rec_cache1, dict)
          assert isinstance(rsarray_index.rec_cache2, dict)
          assert rsarray_index.progress_report == 2
          assert rsarray_index.skip_missing == True
          assert rsarray_index.do_deduplication == False
          assert len(rsarray_index.index_def) == 2
          assert rsarray_index.padded == padded
          assert rsarray_index.str_cmp_thres >= 0.0
          assert rsarray_index.str_cmp_thres <= 1.0
          assert rsarray_index.block_method == (min_q_gram_len, max_block_size)
          assert rsarray_index.num_rec_pairs == None

          assert rsarray_index.status == 'initialised'

          rsarray_index.build()  # - - - - - - - - - - - - - - - - - - - - - -

          assert rsarray_index.num_rec_pairs == None
          assert isinstance(rsarray_index.index1, dict)
          assert isinstance(rsarray_index.index2, dict)
          assert isinstance(rsarray_index.description, str)
          assert isinstance(rsarray_index.rec_cache1, dict)
          assert isinstance(rsarray_index.rec_cache2, dict)
          assert rsarray_index.progress_report == 2
          assert rsarray_index.skip_missing == True
          assert rsarray_index.do_deduplication == False
          assert len(rsarray_index.index_def) == 2
          assert rsarray_index.padded == padded
          assert rsarray_index.str_cmp_thres >= 0.0
          assert rsarray_index.str_cmp_thres <= 1.0
          assert rsarray_index.block_method == (min_q_gram_len, max_block_size)
          assert rsarray_index.status == 'built'

          rsarray_index.compact()  # - - - - - - - - - - - - - - - - - - - - -

          assert rsarray_index.num_rec_pairs > 0
          assert rsarray_index.num_rec_pairs <= \
                 len(self.rec_ident1)*len(self.rec_ident2)
          assert isinstance(rsarray_index.index1, dict)
          assert isinstance(rsarray_index.index2, dict)
          assert isinstance(rsarray_index.description, str)
          assert isinstance(rsarray_index.rec_cache1, dict)
          assert isinstance(rsarray_index.rec_cache2, dict)
          assert rsarray_index.progress_report == 2
          assert rsarray_index.skip_missing == True
          assert rsarray_index.do_deduplication == False
          assert len(rsarray_index.index_def) == 2
          assert rsarray_index.padded == padded
          assert rsarray_index.str_cmp_thres >= 0.0
          assert rsarray_index.str_cmp_thres <= 1.0
          assert rsarray_index.block_method == (min_q_gram_len, max_block_size)

          assert rsarray_index.status == 'compacted'

          [field_names_list, weight_vec_dict] = rsarray_index.run()  # - - - -

          assert isinstance(weight_vec_dict, dict)
          assert len(weight_vec_dict) == rsarray_index.num_rec_pairs
          assert rsarray_index.num_rec_pairs > 0
          assert rsarray_index.num_rec_pairs <= \
                 len(self.rec_ident1)*len(self.rec_ident2)

          for rec_ident1 in self.rec_ident1:
            for rec_ident2 in self.rec_ident2:
              if (rec_ident1 == rec_ident2):
                assert (rec_ident1,rec_ident2) in weight_vec_dict, \
                       (rec_ident1,rec_ident2)

          for (rec_ident1,rec_ident2) in weight_vec_dict:
            assert rec_ident1 in self.rec_ident1
            assert rec_ident2 in self.rec_ident2
            assert isinstance(weight_vec_dict[(rec_ident1,rec_ident2)], list)

            assert len(weight_vec_dict[(rec_ident1,rec_ident2)]) == \
                   len(rsarray_index.rec_comparator.field_comparison_list)

          # All record pairs from the blocking index should be in the suffix
          # array index as well (only when min_g_gram_len = 1)
          #
          if (min_q_gram_len == 1):

            assert len(block_index_weight_vec_dict) <= len(weight_vec_dict), \
                 ('Number of record pairs in blocking index: %d, and in ' % \
                 (len(block_index_weight_vec_dict))+ \
                 'robust  suffix array index: %d' % (len(weight_vec_dict)))

            # Build suffix array weight vector dict with pair identifier sorted
            #
            for ((rec_ident1,rec_ident2),w_list) in weight_vec_dict.items():
              if (rec_ident1 > rec_ident2):
                del weight_vec_dict[(rec_ident1,rec_ident2)]
                weight_vec_dict[(rec_ident2,rec_ident1)] = w_list

            for (rec_ident1,rec_ident2) in block_index_weight_vec_dict:
              assert (rec_ident1,rec_ident2) in weight_vec_dict, \
                   'Record pair %s from blocking index not in suffix array' % \
                   (str((rec_ident1,rec_ident2)))+' index'

            for rec_ident_pair in prev_rec_dict.iterkeys():
              assert rec_ident_pair in weight_vec_dict, \
                 'Record pair %s not in weight vector dict with method: %s' % \
                 (str(rec_ident_pair), (str(rsarray_index.block_method)))

          prev_rec_dict = weight_vec_dict

    # Suffix array index with min_q_gram_len 1 and max_block_size very large
    # should contain all pairs that are in blocking index
    #
    rsarray_index = indexing.RobustSuffixArrayIndex(desc= \
                                             'Test robust suffix array index',
                                             dataset1 = self.dataset1,
                                             dataset2 = self.dataset2,
                                             rec_compar = self.rec_comp_link,
                                             progress=2,
                                             skip_m = False,
                                             padd = False,
                                             log_funct = print_log,
                                             str_cmp_funct = str_cmp_funct,
                                             str_cmp_thres = str_cmp_thres,
                                             block_method = (1,999),
                                             index_def = [index_def1,
                                                          index_def2])
    rsarray_index.build()
    rsarray_index.compact()
    [field_names_list, weight_vec_dict] = rsarray_index.run()

    # Make sure record identifiers in each pair are sorted (smaller id first)
    #
    for ((rec_ident1,rec_ident2),w) in weight_vec_dict.items():
      if (rec_ident1 > rec_ident2):
        del weight_vec_dict[(rec_ident1,rec_ident2)]
        weight_vec_dict[(rec_ident2,rec_ident1)] = w

#    for (rec_ident1,rec_ident2) in weight_vec_dict.keys():
#      assert (rec_ident1,rec_ident2) in block_index_weight_vec_dict, \
#             'Not in blocking weight vectors: %s' % \
#             (str((rec_ident1,rec_ident2)))

    for (rec_ident1,rec_ident2) in block_index_weight_vec_dict.keys():
      assert (rec_ident1,rec_ident2) in weight_vec_dict, \
             'Not in robust suffix array weight vectors: %s' % \
             (str((rec_ident1,rec_ident2)))

    # Test length filter
    #
    [field_names_list, prev_w_vec_dict] = rsarray_index.run()

    for lf in [99,50,30,20,15,10,5,1,0.5]:

      [field_names_list, this_w_vec_dict] = \
                                     rsarray_index.run(length_filter_perc = lf)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

    # Test cut-off threshold
    #
    [field_names_list, prev_w_vec_dict] = rsarray_index.run()

    for cot in [0.05, 0.1, 0.2, 0.5, 0.7, 0.8, 0.9]:

      [field_names_list, this_w_vec_dict] = \
                                     rsarray_index.run(cut_off_threshold = cot)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict


  def testRobustSuffixArrayIndexDedupl(self):  # - - - - - - - - - - - - - - - -
    """Test RobustSuffixArrayIndex deduplication"""

    return

    index_def1 = [['surname','surname',False,False,None,[]]]
    index_def2 = [['given_name','given_name',True,True,4,[]],
                  ['postcode','postcode',True,False,2,[]]]

    # Assume blocking index is OK, so get the record pairs from this index for
    # comparison
    #
    block_index = indexing.BlockingIndex(description = 'Test blocking index',
                                         dataset1 = self.dataset1,
                                         dataset2 = self.dataset1,
                                         rec_comparator = self.rec_comp_dedupl,
                                         progress=4,
                                         index_def = [index_def1,index_def2])
    block_index.build()
    block_index.compact()
    [field_names_list, block_index_weight_vec_dict] = block_index.run()

    # Make sure record identifiers in each pair are sorted (smaller id first)
    #
    for ((rec_ident1,rec_ident2),w) in block_index_weight_vec_dict.items():
      if (rec_ident1 > rec_ident2):
        del block_index_weight_vec_dict[(rec_ident1,rec_ident2)]
        block_index_weight_vec_dict[(rec_ident2,rec_ident1)] = w

    block_method_lists = [[(5,5,False),(4,5,False),(3,5,False),
                           (2,5,False),(1,5,False)],
                          [(5,5,True),(4,5,True),(3,5,True),
                           (2,5,True),(1,5,True)],
                          [(5,10,False),(4,10,False),(3,10,False),
                           (2,10,False),(1,10,False)],
                          [(5,10,True),(4,10,True),(3,10,True),
                           (2,10,True),(1,10,True)]]

    for block_method_list in block_method_lists:
      for (min_q_gram_len, max_block_size, padded) in block_method_list:

        for (str_cmp_funct, str_cmp_thres) in \
            [(stringcmp.jaro,0.7), (stringcmp.jaro,0.8), (stringcmp.jaro,0.9),
             (stringcmp.bigram,0.7),(stringcmp.bigram,0.8),
             (stringcmp.bigram,0.9),
             (stringcmp.lcs, 0.7), (stringcmp.lcs, 0.8), (stringcmp.lcs, 0.9)]:

          prev_rec_dict = {}

          rsarray_index = indexing.RobustSuffixArrayIndex(desc= \
                                      'Test robust suffix-array index',
                                               dataset1 = self.dataset1,
                                               dataset2 = self.dataset1,
                                               rec_comp = self.rec_comp_dedupl,
                                               progress=4,
                                               str_cmp_funct = str_cmp_funct,
                                               str_cmp_thres = str_cmp_thres,
                                               padd = padded,
                                               skip_m = False,
                                               block_method = (min_q_gram_len,
                                                               max_block_size),
                                               index_def = [index_def1,
                                                            index_def2])

          assert isinstance(rsarray_index.index1, dict)
          assert isinstance(rsarray_index.index2, dict)
          assert isinstance(rsarray_index.description, str)
          assert isinstance(rsarray_index.rec_cache1, dict)
          assert isinstance(rsarray_index.rec_cache2, dict)
          assert rsarray_index.progress_report == 4
          assert rsarray_index.skip_missing == False
          assert rsarray_index.do_deduplication == True
          assert len(rsarray_index.index_def) == 2
          assert rsarray_index.padded == padded
          assert isinstance(rsarray_index.str_cmp_thres, float)
          assert rsarray_index.str_cmp_thres >= 0.0
          assert rsarray_index.str_cmp_thres <= 1.0
          assert rsarray_index.block_method == (min_q_gram_len, max_block_size)
          assert rsarray_index.num_rec_pairs == None

          assert rsarray_index.status == 'initialised'

          rsarray_index.build()  # - - - - - - - - - - - - - - - - - - - - - -

          assert rsarray_index.num_rec_pairs == None
          assert isinstance(rsarray_index.index1, dict)
          assert isinstance(rsarray_index.index2, dict)
          assert isinstance(rsarray_index.description, str)
          assert isinstance(rsarray_index.rec_cache1, dict)
          assert isinstance(rsarray_index.rec_cache2, dict)
          assert rsarray_index.progress_report == 4
          assert rsarray_index.skip_missing == False
          assert rsarray_index.do_deduplication == True
          assert len(rsarray_index.index_def) == 2
          assert rsarray_index.padded == padded
          assert isinstance(rsarray_index.str_cmp_thres, float)
          assert rsarray_index.str_cmp_thres >= 0.0
          assert rsarray_index.str_cmp_thres <= 1.0
          assert rsarray_index.block_method == (min_q_gram_len, max_block_size)

          assert rsarray_index.status == 'built'

          rsarray_index.compact()  # - - - - - - - - - - - - - - - - - - - - -

          assert rsarray_index.num_rec_pairs > 0
          assert rsarray_index.num_rec_pairs <= \
                 len(self.rec_ident1)*len(self.rec_ident2)
          assert isinstance(rsarray_index.index1, dict)
          assert isinstance(rsarray_index.index2, dict)
          assert isinstance(rsarray_index.description, str)
          assert isinstance(rsarray_index.rec_cache1, dict)
          assert isinstance(rsarray_index.rec_cache2, dict)
          assert rsarray_index.progress_report == 4
          assert rsarray_index.skip_missing == False
          assert rsarray_index.do_deduplication == True
          assert len(rsarray_index.index_def) == 2
          assert isinstance(rsarray_index.str_cmp_thres, float)
          assert rsarray_index.str_cmp_thres >= 0.0
          assert rsarray_index.str_cmp_thres <= 1.0
          assert rsarray_index.padded == padded
          assert rsarray_index.block_method == (min_q_gram_len, max_block_size)

          assert rsarray_index.status == 'compacted'

          [field_names_list, weight_vec_dict] = rsarray_index.run()  # - - - - -

          assert isinstance(weight_vec_dict, dict)
          assert len(weight_vec_dict) == rsarray_index.num_rec_pairs
          assert rsarray_index.num_rec_pairs > 0
          assert rsarray_index.num_rec_pairs <= \
                 len(self.rec_ident1)*len(self.rec_ident2)

          for rec_ident1 in self.rec_ident1:
            assert (rec_ident1, rec_ident1) not in weight_vec_dict, \
                   (rec_ident1, rec_ident1)

          for (rec_ident1,rec_ident2) in weight_vec_dict:
            assert rec_ident1 in self.rec_ident1
            assert rec_ident2 in self.rec_ident1
            assert isinstance(weight_vec_dict[(rec_ident1,rec_ident2)], list)

            assert len(weight_vec_dict[(rec_ident1,rec_ident2)]) == \
                   len(rsarray_index.rec_comparator.field_comparison_list)

          # All record pairs from the blocking index should be in the suffix
          # array index as well (only when min_g_gram_len = 1)
          #
          if (min_q_gram_len == 1):

            assert len(block_index_weight_vec_dict) <= len(weight_vec_dict), \
                 ('Number of record pairs in blocking index: %d, and in ' % \
                 (len(block_index_weight_vec_dict))+ \
                 ' robust suffix array index: %d' % (len(weight_vec_dict)))

            # Build suffix array weight vector dict with pair identifier sorted
            #
            for ((rec_ident1,rec_ident2),w_list) in weight_vec_dict.items():
              if (rec_ident1 > rec_ident2):
                del weight_vec_dict[(rec_ident1,rec_ident2)]
                weight_vec_dict[(rec_ident2,rec_ident1)] = w_list

            for (rec_ident1,rec_ident2) in block_index_weight_vec_dict:
              assert (rec_ident1,rec_ident2) in weight_vec_dict, \
                   'Record pair %s from blocking index not in suffix array ' \
                    % (str((rec_ident1,rec_ident2)))+'index'

            for rec_ident_pair in prev_rec_dict.iterkeys():
              assert rec_ident_pair in weight_vec_dict, \
                   'Record pair %s not in weight vector dict with method: %s' \
                   % (str(rec_ident_pair), (str(rsarray_index.block_method)))

          prev_rec_dict = weight_vec_dict

    # Suffix array index with min_q_gram_len 1 and max_block_size very large
    # should contain all pairs that are in blocking index
    #
    rsarray_index = indexing.RobustSuffixArrayIndex(desc= \
                                             'Test robust suffix-array index',
                                             dataset1 = self.dataset1,
                                             dataset2 = self.dataset1,
                                             rec_compar = self.rec_comp_dedupl,
                                             progress=4,
                                             str_cmp_funct = str_cmp_funct,
                                             str_cmp_thres = str_cmp_thres,
                                             log_funct = print_log,
                                             skip_m = False,
                                             padd = False,
                                             block_method = (1,999),
                                             index_def = [index_def1,
                                                          index_def2])

    rsarray_index.build()
    rsarray_index.compact()
    [field_names_list, weight_vec_dict] = rsarray_index.run()

    # Make sure record identifiers in each pair are sorted (smaller id first)
    #
    for ((rec_ident1,rec_ident2),w) in weight_vec_dict.items():
      if (rec_ident1 > rec_ident2):
        del weight_vec_dict[(rec_ident1,rec_ident2)]
        weight_vec_dict[(rec_ident2,rec_ident1)] = w

#    for (rec_ident1,rec_ident2) in weight_vec_dict.keys():
#      assert (rec_ident1,rec_ident2) in block_index_weight_vec_dict, \
#             'Not in blocking weight vectors: %s' % \
#             (str((rec_ident1,rec_ident2)))

    for (rec_ident1,rec_ident2) in block_index_weight_vec_dict.keys():
      assert (rec_ident1,rec_ident2) in weight_vec_dict, \
             'Not in robust suffix array weight vectors: %s' % \
             (str((rec_ident1,rec_ident2)))

    # Test length filter
    #
    [field_names_list, prev_w_vec_dict] = rsarray_index.run()

    for lf in [99,50,30,20,15,10,5,1,0.5]:

      [field_names_list, this_w_vec_dict] = \
                                      rsarray_index.run(length_filter_perc = lf)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict

    # Test cut-off threshold
    #
    [field_names_list, prev_w_vec_dict] = rsarray_index.run()

    for cot in [0.05, 0.1, 0.2, 0.5, 0.7, 0.8, 0.9]:

      [field_names_list, this_w_vec_dict] = \
                                     rsarray_index.run(cut_off_threshold = cot)

      assert len(this_w_vec_dict) <= len(prev_w_vec_dict)

      # Check that all weight vectors are also in previous weight vector dict
      #
      for rec_id_pair in this_w_vec_dict:
        assert rec_id_pair in prev_w_vec_dict

        assert this_w_vec_dict[rec_id_pair] == prev_w_vec_dict[rec_id_pair]

      prev_w_vec_dict = this_w_vec_dict


# =============================================================================
# Start tests when called from command line

if (__name__ == "__main__"):
  unittest.main()  # Run all test

# =============================================================================
