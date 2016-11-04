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
# The Original Software is: "datasetTest.py"
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

"""Test module for dataset.py.
"""

# =============================================================================
# Import necessary modules (Python standard modules first, then Febrl modules)

import logging
import os
import string
import sys
import unittest
sys.path.append('..')

import dataset

log_level = logging.WARNING  # logging.INFO

# =============================================================================

class TestCase(unittest.TestCase):

  # Initialise test case  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def setUp(self):
    pass  # Nothing to initialise

  # Clean up test case  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def tearDown(self):
    pass  # Nothing to clean up

  # ---------------------------------------------------------------------------
  # Start test cases

  def testCSV(self):   # - - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test CSV data set"""

    for test_file in ['./test-data.csv','./test-data.csv.gz']:

      # Initialise data set for reading
      #
      test_ds = dataset.DataSetCSV(description='A test CSV data set',
                                   access_mode='read',
                                   field_list=[('rec-id',0),('gname',1),
                                               ('surname',2),('streetnumb',3),
                                               ('streetname_type',4),
                                               ('suburb',5),('postcode',6)],
                                   rec_ident='rec-id',
                                   header_line=False,
                                   write_header=True,
                                   file_name=test_file)

      assert test_ds.dataset_type == 'CSV', \
             'Test data set has wrong type (should be "CSV"): "%s"' % \
             (str(test_ds.dataset_type))
      assert isinstance(test_ds.field_list,list), \
             'CSV data set field list is not of a list: %s' % \
             (str(test_ds.field_list))
      assert test_ds.access_mode == 'read', \
             'CSV data set has wrong access mode (should be "read"): %s' % \
             (str(test_ds.access_mode))
      assert isinstance(test_ds.file_name, str), \
             'CSV data set file name is not a string: %s, %s' % \
             (type(test_ds.file_name), str(test_ds.file_name))
      assert test_ds.num_records == 21, \
             'CSV data set has wrong number of records (should be 21): %d' % \
             (test_ds.num_records)
      assert test_ds.next_rec_num == 0, \
             'CSV data set has wrong next record number (should be 0): %d' % \
             (test_ds.next_rec_num)

      # Read a single record
      #
      test_rec_dict = test_ds.read()

      assert test_ds.next_rec_num == 1, \
             'CSV data set has wrong next record number (should be 1): %d' % \
             (test_ds.next_rec_num)
      assert isinstance(test_rec_dict, dict), \
             'Record returned is not of type dictionary: %s, %s' % \
             (type(test_rec_dict), str(test_rec_dict))
      assert len(test_rec_dict) == 1, \
             'More or less than one record returned: %d, %s' % \
             (len(test_rec_dict), str(test_rec_dict))

      # Read another single record
      #
      test_rec_dict = test_ds.read()

      assert test_ds.next_rec_num == 2, \
             'CSV data set has wrong next record number (should be 2): %d' % \
             (test_ds.next_rec_num)
      assert isinstance(test_rec_dict, dict), \
             'Record returned is not of type dictionary: %s, %s' % \
             (type(test_rec_dict), str(test_rec_dict))
      assert len(test_rec_dict) == 1, \
             'More or less than one record returned: %d, %s' % \
             (len(test_rec_dict), str(test_rec_dict))

      # Read three records
      #
      test_rec_dict = test_ds.read(3)
      assert test_ds.next_rec_num == 5, \
             'CSV data set has wrong next record number (should be 5): %d' % \
             (test_ds.next_rec_num)
      assert isinstance(test_rec_dict, dict), \
             'Record returned is not of type dictionary: %s, %s' % \
             (type(test_rec_dict), str(test_rec_dict))
      assert len(test_rec_dict) == 3, \
             'More or less than three record returned: %d, %s' % \
             (len(test_rec_dict), str(test_rec_dict))

      # Read twenty records (as there are only 21 in the data set only 16
      # should be returned, as 5 have been read before)
      #
      test_rec_dict = test_ds.read(20)

      assert test_ds.next_rec_num == 21, \
             'CSV data set has wrong next record number (should be 21): %d' % \
             (test_ds.next_rec_num)
      assert isinstance(test_rec_dict, dict), \
             'Record returned is not of type dictionary: %s, %s' % \
             (type(test_rec_dict), str(test_rec_dict))
      assert len(test_rec_dict) == 16, \
             'More or less than sixteen record returned: %d, %s' % \
             (len(test_rec_dict), str(test_rec_dict))

      # Read records 3 records from record number 15 onwards
      #
      test_rec_dict = test_ds.read(15,3)

      assert test_ds.next_rec_num == 18, \
             'CSV data set has wrong next record number (should be 18): %d' % \
             (test_ds.next_rec_num)
      assert isinstance(test_rec_dict, dict), \
             'Record returned is not of type dictionary: %s, %s' % \
             (type(test_rec_dict), str(test_rec_dict))
      assert len(test_rec_dict) == 3, \
             'More or less than three record returned: %d, %s' % \
             (len(test_rec_dict), str(test_rec_dict))

      # Read records 5 records from record number 2 onwards
      #
      test_rec_dict = test_ds.read(2,5)

      assert test_ds.next_rec_num == 7, \
             'CSV data set has wrong next record number (should be 7): %d' % \
             (test_ds.next_rec_num)
      assert isinstance(test_rec_dict, dict), \
             'Record returned is not of type dictionary: %s, %s' % \
             (type(test_rec_dict), str(test_rec_dict))
      assert len(test_rec_dict) == 5, \
             'More or less than five record returned: %d, %s' % \
             (len(test_rec_dict), str(test_rec_dict))

      # Read all records
      #
      all_test_rec_dict = test_ds.read(0,30)

      assert test_ds.next_rec_num == 21, \
             'CSV data set has wrong next record number (should be 21): %d' % \
             (test_ds.next_rec_num)
      assert isinstance(all_test_rec_dict, dict), \
             'Record returned is not of type dictionary: %s, %s' % \
             (type(all_test_rec_dict), str(all_test_rec_dict))
      assert len(all_test_rec_dict) == 21, \
             'More or less than tenty-one record returned: %d, %s' % \
             (len(all_test_rec_dict), str(all_test_rec_dict))

      test_ds.finalise()

    # Initialise data set for writing (with quotes set) - - - - - - - - - - - -
    #
    try:
      os.remove('./test-data2.csv')  # Remove previously written test data set
    except:
      pass

    test_ds = dataset.DataSetCSV(description='A test CSV data set',
                                 access_mode='write',
                                 field_list=[('rec-id',0),('gname',1),
                                             ('surname',2),('streetnumb',3),
                                             ('streetname_type',4),
                                             ('suburb',5),('postcode',6)],
                                 rec_ident='rec-id',
                                 header_line=True,
                                 write_header=True,
                                 write_quote="'",
                                 file_name='./test-data2.csv')

    assert test_ds.dataset_type == 'CSV', \
           'Test data set has wrong type (should be "CSV"): "%s"' % \
           (str(test_ds.dataset_type))
    assert isinstance(test_ds.field_list,list), \
           'CSV data set field list is not of a list: %s' % \
           (str(test_ds.field_list))
    assert test_ds.access_mode == 'write', \
           'CSV data set has wrong access mode (should be "write"): %s' % \
           (str(test_ds.access_mode))
    assert isinstance(test_ds.file_name, str), \
           'CSV data set file name is not a string: %s, %s' % \
           (type(test_ds.file_name), str(test_ds.file_name))
    assert test_ds.num_records == 0, \
           'CSV data set has wrong number of records (should be 0): %d' % \
           (test_ds.num_records)
    assert test_ds.next_rec_num == 0, \
           'CSV data set has wrong next record number (should be 0): %d' % \
           (test_ds.next_rec_num)

    # Write a single record
    #
    test_rec_0 = {'00':all_test_rec_dict['00']}

    test_ds.write(test_rec_0)

    assert test_ds.num_records == 1, \
           'CSV data set has wrong number of records: %d (should be 1)' % \
           (test_ds.num_records)
    assert test_ds.next_rec_num == 1, \
           'CSV data set has wrong next record number: %d (should be 1)' % \
           (test_ds.next_rec_num)

    # Write another single record
    #
    test_rec_1 = {'10':all_test_rec_dict['10']}

    test_ds.write(test_rec_1)

    assert test_ds.num_records == 2, \
           'CSV data set has wrong number of records: %d (should be 2)' % \
           (test_ds.num_records)
    assert test_ds.next_rec_num == 2, \
           'CSV data set has wrong next record number: %d (should be 2)' % \
           (test_ds.next_rec_num)

    # Write all records
    #
    test_ds.write(all_test_rec_dict)

    assert test_ds.num_records == 23, \
           'CSV data set has wrong number of records: %d (should be 23)' % \
           (test_ds.num_records)
    assert test_ds.next_rec_num == 23, \
           'CSV data set has wrong next record number: %d (should be 23)' % \
           (test_ds.next_rec_num)

    test_ds.finalise()
    test_ds = None

    # Initialise data set for appending (with quotes set) - - - - - - - - - - -
    #
    test_ds = dataset.DataSetCSV(description='A test CSV data set',
                                 access_mode='append',
                                 field_list=[('rec-id',0),('gname',1),
                                             ('surname',2),('streetnumb',3),
                                             ('streetname_type',4),
                                             ('suburb',5),('postcode',6)],
                                 rec_ident='rec-id',
                                 header_line=True,
                                 write_header=True,
                                 write_quote='"',
                                 file_name='./test-data2.csv')

    assert test_ds.dataset_type == 'CSV', \
           'Test data set has wrong type (should be "CSV"): "%s"' % \
           (str(test_ds.dataset_type))
    assert isinstance(test_ds.field_list,list), \
           'CSV data set field list is not of a list: %s' % \
           (str(test_ds.field_list))
    assert test_ds.access_mode == 'append', \
           'CSV data set has wrong access mode (should be "append"): %s' % \
           (str(test_ds.access_mode))
    assert isinstance(test_ds.file_name, str), \
           'CSV data set file name is not a string: %s, %s' % \
           (type(test_ds.file_name), str(test_ds.file_name))
    assert test_ds.num_records == 23, \
           'CSV data set has wrong number of records (should be 23): %d' % \
           (test_ds.num_records)
    assert test_ds.next_rec_num == 23, \
           'CSV data set has wrong next record number (should be 23): %d' % \
           (test_ds.next_rec_num)

    # Write a single record
    #
    test_rec_0 = {'00':all_test_rec_dict['00']}

    test_ds.write(test_rec_0)

    assert test_ds.num_records == 24, \
           'CSV data set has wrong number of records: %d (should be 24)' % \
           (test_ds.num_records)
    assert test_ds.next_rec_num == 24, \
           'CSV data set has wrong next record number: %d (should be 24)' % \
           (test_ds.next_rec_num)

    # Write another single record
    #
    test_rec_1 = {'10':all_test_rec_dict['10']}

    test_ds.write(test_rec_1)

    assert test_ds.num_records == 25, \
           'CSV data set has wrong number of records: %d (should be 25)' % \
           (test_ds.num_records)
    assert test_ds.next_rec_num == 25, \
           'CSV data set has wrong next record number: %d (should be 25)' % \
           (test_ds.next_rec_num)

    # Write all records
    #
    test_ds.write(all_test_rec_dict)

    assert test_ds.num_records == 46, \
           'CSV data set has wrong number of records: %d (should be 46)' % \
           (test_ds.num_records)
    assert test_ds.next_rec_num == 46, \
           'CSV data set has wrong next record number: %d (should be 46)' % \
           (test_ds.next_rec_num)

    test_ds.finalise()
    test_ds = None

    # Initialise data set for reading - - - - - - - - - - - - - - - - - - - - -
    #
    test_ds = dataset.DataSetCSV(description='A test CSV data set',
                                 access_mode='read',
                                 rec_ident='rec-id',
                                 header_line=True,
                                 write_header=True,
                                 file_name='./test-data2.csv')

    assert test_ds.dataset_type == 'CSV', \
           'Test data set has wrong type (should be "CSV"): "%s"' % \
           (str(test_ds.dataset_type))
    assert isinstance(test_ds.field_list,list), \
           'CSV data set field list is not of a list: %s' % \
           (str(test_ds.field_list))
    assert test_ds.access_mode == 'read', \
           'CSV data set has wrong access mode (should be "read"): %s' % \
           (str(test_ds.access_mode))
    assert isinstance(test_ds.file_name, str), \
           'CSV data set file name is not a string: %s, %s' % \
           (type(test_ds.file_name), str(test_ds.file_name))
    assert test_ds.num_records == 46, \
           'CSV data set has wrong number of records (should be 46): %d' % \
           (test_ds.num_records)
    assert test_ds.next_rec_num == 0, \
           'CSV data set has wrong next record number (should be 0): %d' % \
           (test_ds.next_rec_num)

    # Read a single record
    #
    test_rec_dict = test_ds.read()

    assert test_ds.next_rec_num == 1, \
           'CSV data set has wrong next record number (should be 1): %d' % \
           (test_ds.next_rec_num)
    assert isinstance(test_rec_dict, dict), \
           'Record returned is not of type dictionary: %s, %s' % \
           (type(test_rec_dict), str(test_rec_dict))
    assert len(test_rec_dict) == 1, \
           'More or less than one record returned: %d, %s' % \
           (len(test_rec_dict), str(test_rec_dict))

    # Read another single record
    #
    test_rec_dict = test_ds.read()

    assert test_ds.next_rec_num == 2, \
           'CSV data set has wrong next record number (should be 2): %d' % \
           (test_ds.next_rec_num)
    assert isinstance(test_rec_dict, dict), \
           'Record returned is not of type dictionary: %s, %s' % \
           (type(test_rec_dict), str(test_rec_dict))
    assert len(test_rec_dict) == 1, \
           'More or less than one record returned: %d, %s' % \
           (len(test_rec_dict), str(test_rec_dict))

    # Read three records
    #
    test_rec_dict = test_ds.read(3)
    assert test_ds.next_rec_num == 5, \
           'CSV data set has wrong next record number (should be 5): %d' % \
           (test_ds.next_rec_num)
    assert isinstance(test_rec_dict, dict), \
           'Record returned is not of type dictionary: %s, %s' % \
           (type(test_rec_dict), str(test_rec_dict))
    assert len(test_rec_dict) == 3, \
           'More or less than three record returned: %d, %s' % \
           (len(test_rec_dict), str(test_rec_dict))

    # Readall() iterator
    #
    rec_cnt = 0
    for test_rec_tuple in test_ds.readall():

      assert isinstance(test_rec_tuple, tuple), \
             'Record returned is not of type "tuple": %s, %s' % \
             (type(test_rec_tuple), str(test_rec_tuple))
      assert len(test_rec_tuple) == 2, \
             'More or less than one record returned: %d, %s' % \
             (len(test_rec_tuple), str(test_rec_tuple))

      rec_cnt += 1

      assert test_ds.next_rec_num == rec_cnt, \
             'CSV data set has wrong next record number (should be %d): %d' % \
             (rec_cnt, test_ds.next_rec_num)

    assert test_ds.num_records == rec_cnt, \
           'CSV data set has wrong number of records (should be %d): %d' % \
           (rec_cnt, test_ds.num_records)

    test_ds.finalise()
    test_ds = None

  def testCOL(self):   # - - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test COL data set"""

    test_data_col_width = [6,10,10,13,19,21,11,8]

    for test_file in ['./test-data.col','./test-data.col.gz']:

      # Initialise data set for reading, not using header line
      #
      test_ds = dataset.DataSetCOL(description='A test COL data set',
                                   access_mode='read',
                                   field_list=[('rec-id',6),('gname',10),
                                               ('surname',10),
                                               ('streetnumber',13),
                                               ('address_1',19),
                                               ('address_2',21),
                                               ('suburb',11),('postcode',8)],
                                   rec_ident='rec-id',
                                   header_line=False,
                                   write_header=True,
                                   file_name=test_file)

      assert test_ds.dataset_type == 'COL', \
             'Test data set has wrong type (should be "COL"): "%s"' % \
             (str(test_ds.dataset_type))
      assert isinstance(test_ds.field_list,list), \
             'COL data set field list is not of a list: %s' % \
             (str(test_ds.field_list))
      assert test_ds.access_mode == 'read', \
             'COL data set has wrong access mode (should be "read"): %s' % \
             (str(test_ds.access_mode))
      assert isinstance(test_ds.file_name, str), \
             'COL data set file name is not a string: %s, %s' % \
             (type(test_ds.file_name), str(test_ds.file_name))
      assert test_ds.num_records == 21, \
             'COL data set has wrong number of records (should be 21): %d' % \
             (test_ds.num_records)
      assert test_ds.next_rec_num == 0, \
             'COL data set has wrong next record number (should be 0): %d' % \
             (test_ds.next_rec_num)

      # Read a single record
      #
      test_rec_dict = test_ds.read()

      assert test_ds.next_rec_num == 1, \
             'COL data set has wrong next record number (should be 1): %d' % \
             (test_ds.next_rec_num)
      assert isinstance(test_rec_dict, dict), \
             'Record returned is not of type dictionary: %s, %s' % \
             (type(test_rec_dict), str(test_rec_dict))
      assert len(test_rec_dict) == 1, \
             'More or less than one record returned: %d, %s' % \
             (len(test_rec_dict), str(test_rec_dict))

      # Read another single record
      #
      test_rec_dict = test_ds.read()

      assert test_ds.next_rec_num == 2, \
             'COL data set has wrong next record number (should be 2): %d' % \
             (test_ds.next_rec_num)
      assert isinstance(test_rec_dict, dict), \
             'Record returned is not of type dictionary: %s, %s' % \
             (type(test_rec_dict), str(test_rec_dict))
      assert len(test_rec_dict) == 1, \
             'More or less than one record returned: %d, %s' % \
             (len(test_rec_dict), str(test_rec_dict))

      # Read three records
      #
      test_rec_dict = test_ds.read(3)
      assert test_ds.next_rec_num == 5, \
             'COL data set has wrong next record number (should be 5): %d' % \
             (test_ds.next_rec_num)
      assert isinstance(test_rec_dict, dict), \
             'Record returned is not of type dictionary: %s, %s' % \
             (type(test_rec_dict), str(test_rec_dict))
      assert len(test_rec_dict) == 3, \
             'More or less than three record returned: %d, %s' % \
             (len(test_rec_dict), str(test_rec_dict))

      # Read twenty records (as there are only 21 in the data set only 16
      # should be returned, as 5 have been read before)
      #
      test_rec_dict = test_ds.read(20)

      assert test_ds.next_rec_num == 21, \
             'COL data set has wrong next record number (should be 21): %d' % \
             (test_ds.next_rec_num)
      assert isinstance(test_rec_dict, dict), \
             'Record returned is not of type dictionary: %s, %s' % \
             (type(test_rec_dict), str(test_rec_dict))
      assert len(test_rec_dict) == 16, \
             'More or less than sixteen record returned: %d, %s' % \
             (len(test_rec_dict), str(test_rec_dict))

      # Read records 3 records from record number 15 onwards
      #
      test_rec_dict = test_ds.read(15,3)

      assert test_ds.next_rec_num == 18, \
             'COL data set has wrong next record number (should be 18): %d' % \
             (test_ds.next_rec_num)
      assert isinstance(test_rec_dict, dict), \
             'Record returned is not of type dictionary: %s, %s' % \
             (type(test_rec_dict), str(test_rec_dict))
      assert len(test_rec_dict) == 3, \
             'More or less than three record returned: %d, %s' % \
             (len(test_rec_dict), str(test_rec_dict))

      # Read records 5 records from record number 2 onwards
      #
      test_rec_dict = test_ds.read(2,5)

      assert test_ds.next_rec_num == 7, \
             'COL data set has wrong next record number (should be 7): %d' % \
             (test_ds.next_rec_num)
      assert isinstance(test_rec_dict, dict), \
             'Record returned is not of type dictionary: %s, %s' % \
             (type(test_rec_dict), str(test_rec_dict))
      assert len(test_rec_dict) == 5, \
             'More or less than five record returned: %d, %s' % \
             (len(test_rec_dict), str(test_rec_dict))

      # Read all records
      #
      all_test_rec_dict = test_ds.read(0,30)

      assert test_ds.next_rec_num == 21, \
             'COL data set has wrong next record number (should be 21): %d' % \
             (test_ds.next_rec_num)
      assert isinstance(all_test_rec_dict, dict), \
             'Record returned is not of type dictionary: %s, %s' % \
             (type(all_test_rec_dict), str(all_test_rec_dict))
      assert len(all_test_rec_dict) == 21, \
             'More or less than tenty-one record returned: %d, %s' % \
             (len(all_test_rec_dict), str(all_test_rec_dict))

      test_ds.finalise()

    for test_file in ['./test-data.col','./test-data.col.gz']:

      # Initialise data set for reading, now using header line
      #
      test_ds = dataset.DataSetCOL(description='A test COL data set',
                                   access_mode='read',
                                   field_list=test_data_col_width,
                                   rec_ident='rec_id',
                                   header_line=True,
                                   write_header=True,
                                   file_name=test_file)

      assert test_ds.dataset_type == 'COL', \
             'Test data set has wrong type (should be "COL"): "%s"' % \
             (str(test_ds.dataset_type))
      assert isinstance(test_ds.field_list,list), \
             'COL data set field list is not of a list: %s' % \
             (str(test_ds.field_list))
      assert test_ds.access_mode == 'read', \
             'COL data set has wrong access mode (should be "read"): %s' % \
             (str(test_ds.access_mode))
      assert isinstance(test_ds.file_name, str), \
             'COL data set file name is not a string: %s, %s' % \
             (type(test_ds.file_name), str(test_ds.file_name))
      assert test_ds.num_records == 20, \
             'COL data set has wrong number of records (should be 20): %d' % \
             (test_ds.num_records)
      assert test_ds.next_rec_num == 0, \
             'COL data set has wrong next record number (should be 0): %d' % \
             (test_ds.next_rec_num)

      # Read three records
      #
      test_rec_dict = test_ds.read(3)
      assert test_ds.next_rec_num == 3, \
             'COL data set has wrong next record number (should be 3): %d' % \
             (test_ds.next_rec_num)
      assert isinstance(test_rec_dict, dict), \
             'Record returned is not of type dictionary: %s, %s' % \
             (type(test_rec_dict), str(test_rec_dict))
      assert len(test_rec_dict) == 3, \
             'More or less than three record returned: %d, %s' % \
             (len(test_rec_dict), str(test_rec_dict))

      # Read records 7 records from record number 10 onwards
      #
      test_rec_dict = test_ds.read(10,7)

      assert test_ds.next_rec_num == 17, \
             'COL data set has wrong next record number (should be 17): %d' % \
             (test_ds.next_rec_num)
      assert isinstance(test_rec_dict, dict), \
             'Record returned is not of type dictionary: %s, %s' % \
             (type(test_rec_dict), str(test_rec_dict))
      assert len(test_rec_dict) == 7, \
             'More or less than seven record returned: %d, %s' % \
             (len(test_rec_dict), str(test_rec_dict))

      # Read all records
      #
      all_test_rec_dict = test_ds.read(0,30)

      assert test_ds.next_rec_num == 20, \
             'COL data set has wrong next record number (should be 20): %d' % \
             (test_ds.next_rec_num)
      assert isinstance(all_test_rec_dict, dict), \
             'Record returned is not of type dictionary: %s, %s' % \
             (type(all_test_rec_dict), str(all_test_rec_dict))
      assert len(all_test_rec_dict) == 20, \
             'More or less than tenty-one record returned: %d, %s' % \
             (len(all_test_rec_dict), str(all_test_rec_dict))

      test_ds.finalise()

    # Initialise data set for writing (with header lines) - - - - - - - - - - -
    #
    try:
      os.remove('./test-data2.col')  # Remove previously written test data set
    except:
      pass

    test_ds = dataset.DataSetCOL(description='A test COL data set',
                                 access_mode='write',
                                 field_list=[('rec-id',6),('gname',10),
                                               ('surname',10),
                                               ('streetnumber',13),
                                               ('address_1',19),
                                               ('address_2',21),
                                               ('suburb',11),('postcode',8)],
                                 rec_ident='rec-id',
                                 header_line=True,
                                 write_header=True,
                                 file_name='./test-data2.col')

    assert test_ds.dataset_type == 'COL', \
           'Test data set has wrong type (should be "COL"): "%s"' % \
           (str(test_ds.dataset_type))
    assert isinstance(test_ds.field_list,list), \
           'COL data set field list is not of a list: %s' % \
           (str(test_ds.field_list))
    assert test_ds.access_mode == 'write', \
           'COL data set has wrong access mode (should be "write"): %s' % \
           (str(test_ds.access_mode))
    assert isinstance(test_ds.file_name, str), \
           'COL data set file name is not a string: %s, %s' % \
           (type(test_ds.file_name), str(test_ds.file_name))
    assert test_ds.num_records == 0, \
           'COL data set has wrong number of records (should be 0): %d' % \
           (test_ds.num_records)
    assert test_ds.next_rec_num == 0, \
           'COL data set has wrong next record number (should be 0): %d' % \
           (test_ds.next_rec_num)

    # Write a single record
    #
    test_rec_0 = {'00':all_test_rec_dict['00']}

    test_ds.write(test_rec_0)

    assert test_ds.num_records == 1, \
           'COL data set has wrong number of records: %d (should be 1)' % \
           (test_ds.num_records)
    assert test_ds.next_rec_num == 1, \
           'COL data set has wrong next record number: %d (should be 1)' % \
           (test_ds.next_rec_num)

    # Write another single record
    #
    test_rec_1 = {'10':all_test_rec_dict['10']}

    test_ds.write(test_rec_1)

    assert test_ds.num_records == 2, \
           'COL data set has wrong number of records: %d (should be 2)' % \
           (test_ds.num_records)
    assert test_ds.next_rec_num == 2, \
           'COL data set has wrong next record number: %d (should be 2)' % \
           (test_ds.next_rec_num)

    # Write all records
    #
    test_ds.write(all_test_rec_dict)

    assert test_ds.num_records == 22, \
           'COL data set has wrong number of records: %d (should be 22)' % \
           (test_ds.num_records)
    assert test_ds.next_rec_num == 22, \
           'COL data set has wrong next record number: %d (should be 22)' % \
           (test_ds.next_rec_num)

    test_ds.finalise()
    test_ds = None

    # Initialise data set for appending (no header line) - - - - - - - - - - -
    #
    test_ds = dataset.DataSetCOL(description='A test COL data set',
                                 access_mode='append',
                                 field_list=[('rec-id',6),('gname',10),
                                               ('surname',10),
                                               ('streetnumber',13),
                                               ('address_1',19),
                                               ('address_2',21),
                                               ('suburb',11),('postcode',8)],
                                 rec_ident='rec-id',
                                 header_line=True,
                                 write_header=False,
                                 file_name='./test-data2.col')

    assert test_ds.dataset_type == 'COL', \
           'Test data set has wrong type (should be "COL"): "%s"' % \
           (str(test_ds.dataset_type))
    assert isinstance(test_ds.field_list,list), \
           'COL data set field list is not of a list: %s' % \
           (str(test_ds.field_list))
    assert test_ds.access_mode == 'append', \
           'COL data set has wrong access mode (should be "append"): %s' % \
           (str(test_ds.access_mode))
    assert isinstance(test_ds.file_name, str), \
           'COL data set file name is not a string: %s, %s' % \
           (type(test_ds.file_name), str(test_ds.file_name))
    assert test_ds.num_records == 22, \
           'COL data set has wrong number of records (should be 22): %d' % \
           (test_ds.num_records)
    assert test_ds.next_rec_num == 22, \
           'COL data set has wrong next record number (should be 22): %d' % \
           (test_ds.next_rec_num)

    # Write a single record
    #
    test_rec_0 = {'00':all_test_rec_dict['00']}

    test_ds.write(test_rec_0)

    assert test_ds.num_records == 23, \
           'COL data set has wrong number of records: %d (should be 23)' % \
           (test_ds.num_records)
    assert test_ds.next_rec_num == 23, \
           'COL data set has wrong next record number: %d (should be 23)' % \
           (test_ds.next_rec_num)

    # Write another single record
    #
    test_rec_1 = {'10':all_test_rec_dict['10']}

    test_ds.write(test_rec_1)

    assert test_ds.num_records == 24, \
           'COL data set has wrong number of records: %d (should be 24)' % \
           (test_ds.num_records)
    assert test_ds.next_rec_num == 24, \
           'COL data set has wrong next record number: %d (should be 24)' % \
           (test_ds.next_rec_num)

    # Write all records
    #
    test_ds.write(all_test_rec_dict)

    assert test_ds.num_records == 44, \
           'COL data set has wrong number of records: %d (should be 44)' % \
           (test_ds.num_records)
    assert test_ds.next_rec_num == 44, \
           'COL data set has wrong next record number: %d (should be 44)' % \
           (test_ds.next_rec_num)

    test_ds.finalise()
    test_ds = None

    # Initialise data set for reading - - - - - - - - - - - - - - - - - - - - -
    #
    test_ds = dataset.DataSetCOL(description='A test COL data set',
                                 access_mode='read',
                                 rec_ident='rec_id',
                                 field_list = test_data_col_width,
                                 header_line=True,
                                 write_header=True,
                                 file_name='./test-data2.col')

    assert test_ds.dataset_type == 'COL', \
           'Test data set has wrong type (should be "COL"): "%s"' % \
           (str(test_ds.dataset_type))
    assert isinstance(test_ds.field_list,list), \
           'COL data set field list is not of a list: %s' % \
           (str(test_ds.field_list))
    assert test_ds.access_mode == 'read', \
           'COL data set has wrong access mode (should be "read"): %s' % \
           (str(test_ds.access_mode))
    assert isinstance(test_ds.file_name, str), \
           'COL data set file name is not a string: %s, %s' % \
           (type(test_ds.file_name), str(test_ds.file_name))
    assert test_ds.num_records == 44, \
           'COL data set has wrong number of records (should be 44): %d' % \
           (test_ds.num_records)
    assert test_ds.next_rec_num == 0, \
           'COL data set has wrong next record number (should be 0): %d' % \
           (test_ds.next_rec_num)

    # Read a single record
    #
    test_rec_dict = test_ds.read()

    assert test_ds.next_rec_num == 1, \
           'COL data set has wrong next record number (should be 1): %d' % \
           (test_ds.next_rec_num)
    assert isinstance(test_rec_dict, dict), \
           'Record returned is not of type dictionary: %s, %s' % \
           (type(test_rec_dict), str(test_rec_dict))
    assert len(test_rec_dict) == 1, \
           'More or less than one record returned: %d, %s' % \
           (len(test_rec_dict), str(test_rec_dict))

    # Read another single record
    #
    test_rec_dict = test_ds.read()

    assert test_ds.next_rec_num == 2, \
           'COL data set has wrong next record number (should be 2): %d' % \
           (test_ds.next_rec_num)
    assert isinstance(test_rec_dict, dict), \
           'Record returned is not of type dictionary: %s, %s' % \
           (type(test_rec_dict), str(test_rec_dict))
    assert len(test_rec_dict) == 1, \
           'More or less than one record returned: %d, %s' % \
           (len(test_rec_dict), str(test_rec_dict))

    # Read three records
    #
    test_rec_dict = test_ds.read(3)
    assert test_ds.next_rec_num == 5, \
           'COL data set has wrong next record number (should be 5): %d' % \
           (test_ds.next_rec_num)
    assert isinstance(test_rec_dict, dict), \
           'Record returned is not of type dictionary: %s, %s' % \
           (type(test_rec_dict), str(test_rec_dict))
    assert len(test_rec_dict) == 3, \
           'More or less than three record returned: %d, %s' % \
           (len(test_rec_dict), str(test_rec_dict))

    # Readall() iterator
    #
    rec_cnt = 0
    for test_rec_tuple in test_ds.readall():

      assert isinstance(test_rec_tuple, tuple), \
             'Record returned is not of type "tuple": %s, %s' % \
             (type(test_rec_tuple), str(test_rec_tuple))
      assert len(test_rec_tuple) == 2, \
             'More or less than one record returned: %d, %s' % \
             (len(test_rec_tuple), str(test_rec_tuple))

      rec_cnt += 1

      assert test_ds.next_rec_num == rec_cnt, \
             'COL data set has wrong next record number (should be %d): %d' % \
             (rec_cnt, test_ds.next_rec_num)

    assert test_ds.num_records == rec_cnt, \
           'COL data set has wrong number of records (should be %d): %d' % \
           (rec_cnt, test_ds.num_records)

    test_ds.finalise()
    test_ds = None

  def testMemory(self):   # ---------------------------------------------------
    """Test Memory data set"""

    # First load records from CSV data set
    #
    csv_ds = dataset.DataSetCSV(description='A test CSV data set',
                                access_mode='read',
                                rec_ident='rec-id',
                                header_line=False,
                                field_list=[('rec-id',0),('gname',1),
                                            ('surname',2),('streetnumb',3),
                                            ('streetname_type',4),
                                            ('suburb',5),('postcode',6)],
                                write_header=True,
                                file_name='./test-data.csv')
    all_csv_rec_dict = csv_ds.read(21)  # Read all records
    assert isinstance(all_csv_rec_dict, dict), \
           'Record returned is not of type dictionary: %s, %s' % \
           (type(all_csv_rec_dict), str(all_csv_rec_dict))
    assert len(all_csv_rec_dict) == 21, \
           'CSV data set has wrong next record number (should be 21): %d' % \
           (len(all_csv_rec_dict))
    csv_ds.finalise()

    # Initialise memory data set for read-writing
    #
    test_ds = dataset.DataSetMemory(description='A test Memory data set',
                                    access_mode='readwrite',
                                    field_list=[('rec-id',''),('gname',''),
                                                ('surname',''),
                                                ('streetnumb',''),
                                                ('streetname_type',''),
                                                ('suburb',''),('postcode','')],
                                    rec_ident='rec-id')

    assert test_ds.dataset_type == 'MEMORY', \
           'Test data set has wrong type (should be "MEMORY"): "%s"' % \
           (str(test_ds.dataset_type))
    assert isinstance(test_ds.field_list,list), \
           'Memory data set field list is not of a list: %s' % \
           (str(test_ds.field_list))
    assert test_ds.access_mode == 'readwrite', \
           'Memory data set has wrong access mode (should be "readwrite"):' + \
           '%s' % (str(test_ds.access_mode))
    assert test_ds.num_records == 0, \
           'Memory data set has wrong number of records (should be 0): %d' % \
           (test_ds.num_records)

    # Write a single record
    #
    test_rec_0 = {'00':all_csv_rec_dict['00']}

    test_ds.write(test_rec_0)

    assert test_ds.num_records == 1, \
           'Memory data set has wrong number of records: %d (should be 1)' % \
           (test_ds.num_records)

    # Write another single record
    #
    test_rec_1 = {'10':all_csv_rec_dict['10']}

    test_ds.write(test_rec_1)

    assert test_ds.num_records == 2, \
           'Memory data set has wrong number of records: %d (should be 2)' % \
           (test_ds.num_records)

    # Write all records - should result in 2 warnings of duplicate identifiers
    #
    test_ds.write(all_csv_rec_dict)

    assert test_ds.num_records == 21, \
           'Memory data set has wrong number of records: %d (should be 21)' % \
           (test_ds.num_records)

    # Read a single record
    #
    test_rec_dict = test_ds.read('00')

    assert isinstance(test_rec_dict, dict), \
           'Record returned is not of type dictionary: %s, %s' % \
           (type(test_rec_dict), str(test_rec_dict))
    assert len(test_rec_dict) == 1, \
           'More or less than one record returned: %d, %s' % \
           (len(test_rec_dict), str(test_rec_dict))
    assert '00' in test_rec_dict, \
           'Record identifier "00" is not key in record dictionary:' + \
           '%s' % (str(test_rec_dict))

    # Read a single record
    #
    test_rec_dict = test_ds.read('73')

    assert isinstance(test_rec_dict, dict), \
           'Record returned is not of type dictionary: %s, %s' % \
           (type(test_rec_dict), str(test_rec_dict))
    assert len(test_rec_dict) == 1, \
           'More or less than one record returned: %d, %s' % \
           (len(test_rec_dict), str(test_rec_dict))
    assert '73' in test_rec_dict, \
           'Record identifier "73" is not key in record dictionary:' \
           + '%s' % (str(test_rec_dict))

    # Read three records
    #
    test_rec_dict = test_ds.read(['60','44','41'])

    assert isinstance(test_rec_dict, dict), \
           'Record returned is not of type dictionary: %s, %s' % \
           (type(test_rec_dict), str(test_rec_dict))
    assert len(test_rec_dict) == 3, \
           'More or less than three record returned: %d, %s' % \
           (len(test_rec_dict), str(test_rec_dict))

    # Read another three records (with one duplicate identifier given)
    #
    test_rec_dict = test_ds.read(['20','72','20'])

    assert isinstance(test_rec_dict, dict), \
           'Record returned is not of type dictionary: %s, %s' % \
           (type(test_rec_dict), str(test_rec_dict))
    assert len(test_rec_dict) == 2, \
           'More or less than three record returned: %d, %s' % \
           (len(test_rec_dict), str(test_rec_dict))

    # Readall() iterator
    #
    rec_cnt = 0
    for test_rec_tuple in test_ds.readall():

      assert isinstance(test_rec_tuple, tuple), \
             'Record returned is not of type "tuple": %s, %s' % \
             (type(test_rec_tuple), str(test_rec_tuple))
      assert len(test_rec_tuple) == 2, \
             'More or less than one record returned: %d, %s' % \
             (len(test_rec_tuple), str(test_rec_euple))

      rec_cnt += 1

    assert test_ds.num_records == rec_cnt, \
           'CSV data set has wrong number of records (should be %d): %d' % \
           (rec_cnt, test_ds.num_records)

    test_ds.finalise()
    test_ds = None

  def testShelve(self):   # ---------------------------------------------------
    """Test Shelve data set"""

    # First load records from CSV data set
    #
    csv_ds = dataset.DataSetCSV(description='A test CSV data set',
                                access_mode='read',
                                rec_ident='rec-id',
                                header_line=False,
                                field_list=[('rec-id',0),('gname',1),
                                            ('surname',2),('streetnumb',3),
                                            ('streetname_type',4),
                                            ('suburb',5),('postcode',6)],
                                write_header=True,
                                file_name='./test-data.csv')
    all_csv_rec_dict = csv_ds.read(21)  # Read all records
    assert isinstance(all_csv_rec_dict, dict), \
           'Record returned is not of type dictionary: %s, %s' % \
           (type(all_csv_rec_dict), str(all_csv_rec_dict))
    assert len(all_csv_rec_dict) == 21, \
           'CSV data set has wrong next record number (should be 21): %d' % \
           (len(all_csv_rec_dict))
    csv_ds.finalise()

    # Initialise shelve data set for read-writing
    #
    test_ds = dataset.DataSetShelve(description='A test Shelve data set',
                                    file_name = 'test-data.slv',
                                    clear = True,  # False,
                                    access_mode='readwrite',
                                    field_list=[('rec-id',''),('gname',''),
                                                ('surname',''),
                                                ('streetnumb',''),
                                                ('streetname_type',''),
                                                ('suburb',''),('postcode','')],
                                    rec_ident='rec-id')

    assert test_ds.dataset_type == 'SHELVE', \
           'Test data set has wrong type (should be "SHELVE"): "%s"' % \
           (str(test_ds.dataset_type))
    assert isinstance(test_ds.field_list,list), \
           'Shelve data set field list is not of a list: %s' % \
           (str(test_ds.field_list))
    assert test_ds.access_mode == 'readwrite', \
           'Shelve data set has wrong access mode (should be "readwrite"):' + \
           '%s' % (str(test_ds.access_mode))
    assert test_ds.num_records == 0, \
           'Shelve data set has wrong number of records (should be 0): %d' % \
           (test_ds.num_records)

    # Write a single record
    #
    test_rec_0 = {'00':all_csv_rec_dict['00']}

    test_ds.write(test_rec_0)

    assert test_ds.num_records == 1, \
           'Shelve data set has wrong number of records: %d (should be 1)' % \
           (test_ds.num_records)

    # Write another single record
    #
    test_rec_1 = {'10':all_csv_rec_dict['10']}

    test_ds.write(test_rec_1)

    assert test_ds.num_records == 2, \
           'Shelve data set has wrong number of records: %d (should be 2)' % \
           (test_ds.num_records)

    # Write all records - should result in 2 warnings of duplicate identifiers
    #
    test_ds.write(all_csv_rec_dict)

    assert test_ds.num_records == 21, \
           'Shelve data set has wrong number of records: %d (should be 21)' % \
           (test_ds.num_records)

    # Read a single record
    #
    test_rec_dict = test_ds.read('00')

    assert isinstance(test_rec_dict, dict), \
           'Record returned is not of type dictionary: %s, %s' % \
           (type(test_rec_dict), str(test_rec_dict))
    assert len(test_rec_dict) == 1, \
           'More or less than one record returned: %d, %s' % \
           (len(test_rec_dict), str(test_rec_dict))
    assert '00' in test_rec_dict, \
           'Record identifier "00" is not key in record dictionary:' + \
           '%s' % (str(test_rec_dict))

    # Read a single record
    #
    test_rec_dict = test_ds.read('73')

    assert isinstance(test_rec_dict, dict), \
           'Record returned is not of type dictionary: %s, %s' % \
           (type(test_rec_dict), str(test_rec_dict))
    assert len(test_rec_dict) == 1, \
           'More or less than one record returned: %d, %s' % \
           (len(test_rec_dict), str(test_rec_dict))
    assert '73' in test_rec_dict, \
           'Record identifier "73" is not key in record dictionary:' \
           + '%s' % (str(test_rec_dict))

    # Read three records
    #
    test_rec_dict = test_ds.read(['60','44','41'])

    assert isinstance(test_rec_dict, dict), \
           'Record returned is not of type dictionary: %s, %s' % \
           (type(test_rec_dict), str(test_rec_dict))
    assert len(test_rec_dict) == 3, \
           'More or less than three record returned: %d, %s' % \
           (len(test_rec_dict), str(test_rec_dict))

    # Read another three records (with one duplicate identifier given)
    #
    test_rec_dict = test_ds.read(['20','72','20'])

    assert isinstance(test_rec_dict, dict), \
           'Record returned is not of type dictionary: %s, %s' % \
           (type(test_rec_dict), str(test_rec_dict))
    assert len(test_rec_dict) == 2, \
           'More or less than three record returned: %d, %s' % \
           (len(test_rec_dict), str(test_rec_dict))

    # Readall() iterator
    #
    rec_cnt = 0
    for test_rec_tuple in test_ds.readall():

      assert isinstance(test_rec_tuple, tuple), \
             'Record returned is not of type "tuple": %s, %s' % \
             (type(test_rec_tuple), str(test_rec_tuple))
      assert len(test_rec_tuple) == 2, \
             'More or less than one record returned: %d, %s' % \
             (len(test_rec_tuple), str(test_rec_tuple))

      rec_cnt += 1

    assert test_ds.num_records == rec_cnt, \
           'CSV data set has wrong number of records (should be %d): %d' % \
           (rec_cnt, test_ds.num_records)

    test_ds.finalise()
    test_ds = None

  def testCSVdelimiter(self):   # - - - - - - - - - - - - - - - - - - - - - - -
    """Test CSV data set with different delimiters"""

    for (test_file, delim, num_col) in [('./test-comma.txt', ',', 6),
                                        ('./test-comma.txt', 'x', 1),
                                        ('./test-period.txt', ':', 6),
                                        ('./test-period.txt', ',', 1),
                                        ('./test-tabulator.txt', '\t', 6),
                                        ('./test-tabulator.txt', ',', 1),
                                        ('./test-tabulator-comma.txt', ',', 6),
                                        ('./test-tabulator-comma.txt', ',', 6),
                                        ('./test-tabulator-comma.txt',';', 1),
                                        ('./test-tabulator.txt', chr(9), 6)]:
      # Create a field list
      #
      field_name_list = []
      for i in range(num_col):
        field_name_list.append(('field-%d' % (i), i))

      # Initialise data set for reading
      #
      test_ds = dataset.DataSetCSV(description='A test CSV data set',
                                   access_mode='read',
                                   field_list = field_name_list,
                                   delimiter = delim,
                                   rec_ident='rec-id',
                                   header_line=False,
                                   strip_fields=False,
                                   write_header=True,
                                   file_name=test_file)

      assert test_ds.dataset_type == 'CSV', \
             'Test data set has wrong type (should be "CSV"): "%s"' % \
             (str(test_ds.dataset_type))
      assert isinstance(test_ds.field_list,list), \
             'CSV data set field list is not of a list: %s' % \
             (str(test_ds.field_list))
      assert test_ds.access_mode == 'read', \
             'CSV data set has wrong access mode (should be "read"): %s' % \
             (str(test_ds.access_mode))
      assert test_ds.delimiter == delim, \
             'CSV data set has wrong delimiter (should be "%s"): %s' % \
             (delim, str(test_ds.delimiter))
      assert len(test_ds.field_list) == num_col, \
             'CSV data set has wrong number of fields (should be %d): %d' % \
             (num_col, len(test_ds.field_list))

      # Now check for delimiters in the data set
      #
      rec_cnt = 0
      for (rec_id, rec_data) in test_ds.readall():
        for col_val in rec_data:
          assert delim not in col_val, \
                 'Delimiter character "%s" in data set values: "%s"' % \
                 (delim, col_val)

        if (('tabulator-comma' in test_file) and (delim not in [',', '\t'])):
          assert len(rec_data) == 1
        elif (('tabulator-comma' in test_file) and (delim in [',', '\t'])):
          assert len(rec_data) == 6

        elif (('comma' in test_file) and (delim == 'x')):
          assert len(rec_data) == 1
        elif (('comma' in test_file) and (delim == ',')):
          assert len(rec_data) == 6

        elif (('tabulator' in test_file) and (delim != '\t')):
          assert len(rec_data) == 1
        elif (('tabulator' in test_file) and (delim == '\t')):
          assert len(rec_data) == 6

        elif (('period' in test_file) and (delim != ':')):
          assert len(rec_data) == 1
        elif (('period' in test_file) and (delim == ':')):
          assert len(rec_data) == 6

        # Check the values in the current record
        #
        if (num_col == 6):
          col_cnt = 0

          for col_val in rec_data:
            this_num = (rec_cnt+col_cnt) % 10
            check_col_val = str(this_num)*(col_cnt+1)

            if (('tabulator-comma' in test_file) and (delim == '\t') and \
                (col_cnt > 0)):
              check_col_val = ','+check_col_val

            elif (('tabulator-comma' in test_file) and (delim == ',') and \
                  (col_cnt <5)):
              check_col_val = check_col_val+'\t'

            assert col_val == check_col_val, 'Different column value than ' + \
                   'expected in row %d: "%s" (should be "%s"' % \
                   (rec_cnt, col_val, check_col_val)

            col_cnt += 1

        rec_cnt += 1

      test_ds.finalise()
      test_ds = None

  def testMissValStripFields(self):   # - - - - - - - - - - - - - - - - - - - -
    """Test data sets regarding missing values and strip fields arguments"""

    for (miss_val, strip_fields_flag) in [('15', False),('15', True),
                                          (['3'], False),(['3'], True),
                                          (['3', '2906'], False),
                                          (['3', '2906'], True),
                                          (['xyz', '1'], False),
                                          (['xyz', '1'], True),
                                          ('rivett', False),
                                          ('rivett', True)]:

      test_ds = dataset.DataSetCSV(description='A test CSV data set',
                                   access_mode='read',
                                   field_list=[('rec-id',0),('gname',1),
                                               ('surname',2),('streetnumb',3),
                                               ('streetname_type',4),
                                               ('suburb',5),('postcode',6)],
                                   rec_ident='rec-id',
                                   header_line=False,
                                   write_header=True,
                                   miss_v=miss_val,
                                   strip_f=strip_fields_flag,
                                   file_name='./test-data.csv')

      assert test_ds.strip_fields in [True, False], \
             'CSV data set strip fields is not a boolean: %s' % \
             (str(test_ds.strip_fields))
      assert test_ds.strip_fields == strip_fields_flag, \
             'CSV data set has wrong strip fields values: %s (should be: %s)' \
             % (str(test_ds.strip_fields), str(strip_fields_flag))

      assert isinstance(test_ds.miss_val, list), \
             'CSV data set missing values variable is not a list: %s' % \
             (str(test_ds.miss_val))
      if (isinstance(miss_val, str)):
        tmp_miss_val = [miss_val]
      else:
        tmp_miss_val = miss_val
      assert test_ds.miss_val == tmp_miss_val, \
             'CSV data set has wrong missing values: %s (should be: %s)' % \
             (str(test_ds.miss_val), str(tmp_miss_val))

      # Now check all values in the records from the data set
      #
      rec_cnt = 0
      for (rec_id, rec_data) in test_ds.readall():
        for col_val in rec_data:
          if ((strip_fields_flag == True) and (len(col_val) > 0)):
            assert col_val[0] not in string.whitespace, \
                   'Field value not stripped of whitespaces: "%s"' % (col_val)
            assert col_val[-1] not in string.whitespace, \
                   'Field value not stripped of whitespaces: "%s"' % (col_val)

          assert col_val not in tmp_miss_val, \
                 'Missing value not removed: "%s" (from list: %s' % \
                 (col_val, str(tmp_miss_val))

      test_ds.finalise()
      test_ds = None

      test_ds = dataset.DataSetCOL(description='A test COL data set',
                                   access_mode='read',
                                   field_list=[('rec-id',6),('gname',10),
                                               ('surname',10),
                                               ('streetnumber',13),
                                               ('address_1',19),
                                               ('address_2',21),
                                               ('suburb',11),('postcode',8)],
                                   rec_ident='rec-id',
                                   header_line=False,
                                   write_header=True,
                                   miss_v=miss_val,
                                   strip_f=strip_fields_flag,
                                   file_name='./test-data.col')

      assert test_ds.strip_fields in [True, False], \
             'COL data set strip fields is not a boolean: %s' % \
             (str(test_ds.strip_fields))
      assert test_ds.strip_fields == strip_fields_flag, \
             'COL data set has wrong strip fields values: %s (should be: %s)' \
             % (str(test_ds.strip_fields), str(strip_fields_flag))

      assert isinstance(test_ds.miss_val, list), \
             'COL data set missing values variable is not a list: %s' % \
             (str(test_ds.miss_val))
      if (isinstance(miss_val, str)):
        tmp_miss_val = [miss_val]
      else:
        tmp_miss_val = miss_val
      assert test_ds.miss_val == tmp_miss_val, \
             'COL data set has wrong missing values: %s (should be: %s)' % \
             (str(test_ds.miss_val), str(tmp_miss_val))

      # Now check all values in the records from the data set
      #
      rec_cnt = 0
      for (rec_id, rec_data) in test_ds.readall():
        for col_val in rec_data:
          if ((strip_fields_flag == True) and (len(col_val) > 0)):
            assert col_val[0] not in string.whitespace, \
                   'Field value not stripped of whitespaces: "%s"' % (col_val)
            assert col_val[-1] not in string.whitespace, \
                   'Field value not stripped of whitespaces: "%s"' % (col_val)

          assert col_val not in tmp_miss_val, \
                 'Missing value not removed: "%s" (from list: %s' % \
                 (col_val, str(tmp_miss_val))

      test_ds.finalise()
      test_ds = None

      test_ds = dataset.DataSetMemory(description='A test Memory data set',
                                      access_mode='readwrite',
                                      field_list=[('rec-id',''),('gname',''),
                                                  ('surname',''),
                                                  ('streetnumb',''),
                                                  ('streetname_type',''),
                                                  ('suburb',''),
                                                  ('postcode','')],
                                      miss_v=miss_val,
                                      strip_f=strip_fields_flag,
                                      rec_ident='rec-id')

      assert test_ds.strip_fields in [True, False], \
             'MEM data set strip fields is not a boolean: %s' % \
             (str(test_ds.strip_fields))
      assert test_ds.strip_fields == strip_fields_flag, \
             'MEM data set has wrong strip fields values: %s (should be: %s)' \
             % (str(test_ds.strip_fields), str(strip_fields_flag))

      assert isinstance(test_ds.miss_val, list), \
             'MEM data set missing values variable is not a list: %s' % \
             (str(test_ds.miss_val))
      if (isinstance(miss_val, str)):
        tmp_miss_val = [miss_val]
      else:
        tmp_miss_val = miss_val
      assert test_ds.miss_val == tmp_miss_val, \
             'MEM data set has wrong missing values: %s (should be: %s)' % \
             (str(test_ds.miss_val), str(tmp_miss_val))

      # Write some records
      #
      test_data = {'00':['50','mitchell','polmear','341','fitchett street',
                         '',"    o'connor  ",'  2906   '],
                   '01':['61','isaad','whte  ','  15','  rivett circuit','',
                         'rivett','   2906'],
                   '02':['  62  ','isaac','wiglht','  15  ','tyrrell circuit',
                         'rivett  ','2906   ']}
      test_ds.write(test_data)

      # Now check all values in the records from the data set
      #
      rec_cnt = 0
      for (rec_id, rec_data) in test_ds.readall():
        for col_val in rec_data:
          if ((strip_fields_flag == True) and (len(col_val) > 0)):
            assert col_val[0] not in string.whitespace, \
                   'Field value not stripped of whitespaces: "%s"' % (col_val)
            assert col_val[-1] not in string.whitespace, \
                   'Field value not stripped of whitespaces: "%s"' % (col_val)

          assert col_val not in tmp_miss_val, \
                 'Missing value not removed: "%s" (from list: %s' % \
                 (col_val, str(tmp_miss_val))

      test_ds.finalise()
      test_ds = None



      test_ds = dataset.DataSetShelve(description='A test Shelve data set',
                                      file_name = 'test-data.slv',
                                      clear = True,  # False,
                                      access_mode='readwrite',
                                      field_list=[('rec-id',''),('gname',''),
                                                  ('surname',''),
                                                  ('streetnumb',''),
                                                  ('streetname_type',''),
                                                  ('suburb',''),
                                                  ('postcode','')],
                                      miss_v=miss_val,
                                      strip_f=strip_fields_flag,
                                      rec_ident='rec-id')

      assert test_ds.strip_fields in [True, False], \
             'SHL data set strip fields is not a boolean: %s' % \
             (str(test_ds.strip_fields))
      assert test_ds.strip_fields == strip_fields_flag, \
             'SHL data set has wrong strip fields values: %s (should be: %s)' \
             % (str(test_ds.strip_fields), str(strip_fields_flag))

      assert isinstance(test_ds.miss_val, list), \
             'SHL data set missing values variable is not a list: %s' % \
             (str(test_ds.miss_val))
      if (isinstance(miss_val, str)):
        tmp_miss_val = [miss_val]
      else:
        tmp_miss_val = miss_val
      assert test_ds.miss_val == tmp_miss_val, \
             'SHL data set has wrong missing values: %s (should be: %s)' % \
             (str(test_ds.miss_val), str(tmp_miss_val))

      # Write some records
      #
      test_data = {'00':['50','mitchell','polmear','341','fitchett street',
                         '',"    o'connor  ",'  2906   '],
                   '01':['61','isaad','whte  ','  15','  rivett circuit','',
                         'rivett','   2906'],
                   '02':['  62  ','isaac','wiglht','  15  ','tyrrell circuit',
                         'rivett  ','2906   ']}
      test_ds.write(test_data)

      # Now check all values in the records from the data set
      #
      rec_cnt = 0
      for (rec_id, rec_data) in test_ds.readall():
        for col_val in rec_data:
          if ((strip_fields_flag == True) and (len(col_val) > 0)):
            assert col_val[0] not in string.whitespace, \
                   'Field value not stripped of whitespaces: "%s"' % (col_val)
            assert col_val[-1] not in string.whitespace, \
                   'Field value not stripped of whitespaces: "%s"' % (col_val)

          assert col_val not in tmp_miss_val, \
                 'Missing value not removed: "%s" (from list: %s' % \
                 (col_val, str(tmp_miss_val))

      test_ds.finalise()
      test_ds = None

# =============================================================================
# Start tests when called from command line

if (__name__ == "__main__"):

  # Intialise a logger, set level to info
  #
  my_logger = logging.getLogger()  # New logger at root level
  my_logger.setLevel(log_level)

  unittest.main()  # Run all test

# =============================================================================
