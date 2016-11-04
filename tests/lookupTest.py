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
# The Original Software is: "lookupTest.py"
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

"""Test module for lookup.py.
"""

# =============================================================================
# Import necessary modules (Python standard modules first, then Febrl modules)

import os
import sys
import unittest
sys.path.append('..')

import lookup

# =============================================================================

class TestCase(unittest.TestCase):

  # Initialise test case  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def setUp(self):

    src_data_dir = '..'+os.sep+'data'+os.sep+'lookup'+os.sep

    self.tag_lookup_files = [src_data_dir+'post_address.tbl',
                             src_data_dir+'postcode.tbl',
                             src_data_dir+'address_misc.tbl',
                             src_data_dir+'saints.tbl',
                             src_data_dir+'givenname_f.tbl',
                             src_data_dir+'givenname_m.tbl',
                             src_data_dir+'surname.tbl']

    self.frequency_lookup_files = [src_data_dir+'givenname_f_freq.csv',
                                   src_data_dir+'suburb_nsw_freq.csv',
                                   src_data_dir+'postcode_act_freq.csv']

    self.geocode_lookup_files = [src_data_dir+'postcode_centroids.csv']

    self.correction_list_files = [src_data_dir+'address_corr.lst',
                                  src_data_dir+'name_corr.lst']

  # Clean up test case  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def tearDown(self):
    pass  # Nothing to clean up

  # ---------------------------------------------------------------------------
  #
  # Start test cases

  def testTagLookupTables(self):  # - - - - - - - - - - - - - - - - - - - - - -
    """Test tag look-up tables"""

    # First load all files separately
    #
    for f in self.tag_lookup_files:
      lookup_table = lookup.TagLookupTable(descr=f, default = '')

      assert (lookup_table.description == f), \
             'Look-up table "'+f+'" does not have correct description: "'+ \
             lookup_table.description+'"'

      assert (lookup_table.file_names == []), \
             'Look-up table "'+f+'" has non-empty file list: '+ \
             str(lookup_table.file_names)

      assert (isinstance(lookup_table,dict)), \
             'Look-up table "'+f+'" is not a dictionary'

      lookup_table.load(f)  # Load the table

      assert (len(lookup_table) > 0), \
             'Look-up table "'+f+'" is empty after load'

      assert (lookup_table['xyz1234zyx'] == lookup_table.default), \
             'Look-up table "'+f+'" returns wrond default: '+ \
             lookup_table.default

      assert (lookup_table.get('xyz1234zyx') == lookup_table.default), \
             'Look-up table "'+f+'" returns wrond default: '+ \
             lookup_table.default

      assert (lookup_table.max_key_length > 0), \
             'Look-up table "'+f+'" has illegal key length: '+ \
             str(lookup_table.max_key_length)

      for (key,value) in lookup_table.items():
        assert (isinstance(key,tuple)), \
               'Key in look-up table "'+f+'" is not a tuple: '+str(key)
        assert (len(value) == 2), \
               'Value in look-up table "'+f+'" does not contain two '+ \
               'elements: '+str(values)

    # Now load all files into one look-up table
    #
    lookup_table = lookup.TagLookupTable(descr=self.tag_lookup_files[0], \
                   default = '')

    assert (lookup_table.description == self.tag_lookup_files[0]), \
           'Combined look-up table does not have correct description: "'+ \
           str(self.tag_lookup_files[0])+'"'

    assert (lookup_table.file_names == []), \
           'Combined look-up table has non-empty file list: '+ \
           str(lookup_table.file_names)

    assert (isinstance(lookup_table,dict)), \
           'Combined lookup-table is not a '+'dictionary'

    lookup_table.load(self.tag_lookup_files)  # Load the table

    assert (len(lookup_table) > 0), \
           'Combined look-up table is empty after load'

    assert (lookup_table['xyz1234zyx'] == lookup_table.default), \
           'Combined look-up table returns wrond default: '+ \
           lookup_table.default

    assert (lookup_table.get('xyz1234zyx') == lookup_table.default), \
           'Combined look-up table returns wrond default: '+ \
           lookup_table.default

    assert (lookup_table.max_key_length > 0), \
           'Combined look-up table has illegal key length: '+ \
           str(lookup_table.max_key_length)

    for (key,value) in lookup_table.items():
      assert (isinstance(key,tuple)), \
             'Key in combined look-up table is not a tuple: '+str(key)
      assert (len(value) == 2), \
             'Value in combined look-up table does not contain two '+ \
             'elements: '+str(values)

  def testFrequencyLookupTables(self):  # - - - - - - - - - - - - - - - - - - -
    """Test frequency look-up tables"""

    # First load all files separately
    #
    for f in self.frequency_lookup_files:
      lookup_table = lookup.FrequencyLookupTable(descr = f, default = '')

      assert (lookup_table.description == f), \
             'Look-up table "'+f+'" does not have correct description: "'+ \
             lookup_table.description+'"'

      assert (lookup_table.file_names == []), \
             'Look-up table "'+f+'" has non-empty file list: '+ \
             str(lookup_table.file_names)

      assert (isinstance(lookup_table,dict)), \
             'Look-up table "'+f+'" is not a dictionary'

      assert (lookup_table.sum == None), \
             'Look-up table "'+f+'" has non-None sum'

      lookup_table.load(f)  # Load the table

      assert (len(lookup_table) > 0), \
             'Look-up table "'+f+'" is empty after load'

      assert (lookup_table.sum > 0), \
             'Look-up table "'+f+'" has zero sum'

      assert (lookup_table.sum > lookup_table.length), \
             'Look-up table "'+f+'" has sum smaller than number of entries: '+\
             str(self.length)

      assert (lookup_table['xyz1234zyx'] == lookup_table.default), \
             'Look-up table "'+f+'" returns wrond default: '+ \
             lookup_table.default

      assert (lookup_table.get('xyz1234zyx') == lookup_table.default), \
             'Look-up table "'+f+'" returns wrond default: '+ \
             lookup_table.default

      for (key,value) in lookup_table.items():
        assert (isinstance(key,str)), \
               'Key in look-up table "'+f+'" is not a string: '+str(key)
        assert (isinstance(value,int) and (value > 0)), \
               'Value in look-up table "'+f+'" is not an integer: '+str(value)

    # Now load all files into one look-up table
    #
    lookup_table = lookup.FrequencyLookupTable( \
                   descr=self.frequency_lookup_files[0], \
                   default = '')

    assert (lookup_table.description == self.frequency_lookup_files[0]), \
           'Combined look-up table does not have correct description: "'+ \
           str(self.frequency_lookup_files[0])+'"'

    assert (lookup_table.file_names == []), \
           'Combined look-up table has non-empty file list: '+ \
           str(lookup_table.file_names)

    assert (isinstance(lookup_table,dict)), \
           'Combined lookup-table is not a '+'dictionary'

    assert (lookup_table.sum == None), \
           'Combined look-up table has non-None sum'

    lookup_table.load(self.frequency_lookup_files)  # Load the table

    assert (len(lookup_table) > 0), \
           'Combined look-up table is empty after load'

    assert (lookup_table.sum > 0), \
           'Combined look-up table has zero sum'

    assert (lookup_table['xyz1234zyx'] == lookup_table.default), \
           'Combined look-up table returns wrond default: '+ \
           lookup_table.default

    assert (lookup_table.get('xyz1234zyx') == lookup_table.default), \
           'Combined look-up table returns wrond default: '+ \
           lookup_table.default

    for (key,value) in lookup_table.items():
      assert (isinstance(key,str)), \
             'Key in combined look-up table is not a string: '+str(key)
      assert (isinstance(value,int) and (value > 0)), \
             'Value in combined look-up table is not an integer: '+str(value)


  def testGeocodeLookupTables(self):  # - - - - - - - - - - - - - - - - - - - -
    """Test geocode look-up tables"""

    # First load all files separately
    #
    for f in self.geocode_lookup_files:
      lookup_table = lookup.GeocodeLookupTable(descrp = f, default = '')

      assert (lookup_table.description == f), \
             'Look-up table "'+f+'" does not have correct description: "'+ \
             lookup_table.description+'"'

      assert (lookup_table.file_names == []), \
             'Look-up table "'+f+'" has non-empty file list: '+ \
             str(lookup_table.file_names)

      assert (isinstance(lookup_table,dict)), \
             'Look-up table "'+f+'" is not a dictionary'

      lookup_table.load(f)  # Load the table

      assert (len(lookup_table) > 0), \
             'Look-up table "'+f+'" is empty after load'

      assert (lookup_table['xyz1234zyx'] == lookup_table.default), \
             'Look-up table "'+f+'" returns wrond default: '+ \
             lookup_table.default

      assert (lookup_table.get('xyz1234zyx') == lookup_table.default), \
             'Look-up table "'+f+'" returns wrond default: '+ \
             lookup_table.default

      for (key,value) in lookup_table.items():
        assert (isinstance(key,str)), \
               'Key in look-up table "'+f+'" is not a string: '+str(key)
        assert (len(value) == 2) and \
                isinstance(value[0],float) and isinstance(value[1],float), \
               'Value in look-up table "'+f+'" is not an location: '+str(value)

        assert (value[0] >= -180.0) and (value[0] <= 180.0), \
               'Location in look-up table "'+f+'" has illegal longitude: '+ \
               str(value[0])
        assert (value[1] >= -90.0) and (value[1] <= 90.0), \
               'Location in look-up table "'+f+'" has illegal latitude: '+ \
               str(value[1])

    # Now load all files into one look-up table
    #
    lookup_table = lookup.GeocodeLookupTable( \
                   descr = self.geocode_lookup_files[0], \
                   default = '')

    assert (lookup_table.description == self.geocode_lookup_files[0]), \
           'Combined look-up table does not have correct description: "'+ \
           str(self.geocode_lookup_files[0])+'"'

    assert (lookup_table.file_names == []), \
           'Combined look-up table has non-empty file list: '+ \
           str(lookup_table.file_names)

    assert (isinstance(lookup_table,dict)), \
           'Combined lookup-table is not a '+'dictionary'

    lookup_table.load(self.geocode_lookup_files)  # Load the table

    assert (len(lookup_table) > 0), \
           'Combined look-up table is empty after load'

    assert (lookup_table['xyz1234zyx'] == lookup_table.default), \
           'Combined look-up table returns wrond default: '+ \
           lookup_table.default

    assert (lookup_table.get('xyz1234zyx') == lookup_table.default), \
           'Combined look-up table returns wrond default: '+ \
           lookup_table.default

    for (key,value) in lookup_table.items():
      assert (isinstance(key,str)), \
             'Key in combined look-up table is not a string: '+str(key)
      assert (len(value) == 2) and \
              isinstance(value[0],float) and isinstance(value[1],float), \
             'Value in combined look-up table is not an location: '+str(value)

      assert (value[0] >= -180.0) and (value[0] <= 180.0), \
             'Location in combined look-up table has illegal longitude: '+ \
             str(value[0])
      assert (value[1] >= -90.0) and (value[1] <= 90.0), \
             'Location in combined look-up table has illegal latitude: '+ \
             str(value[1])

  def testCorrectionLists(self):  # - - - - - - - - - - - - - - - - - - - - - -
    """Test correction lists"""

    for f in self.correction_list_files:
      corr_list = lookup.CorrectionList(descr = f)

      assert (corr_list.description == f), \
             'Correction list "'+f+'" does not have correct description: "'+ \
             corr_list.description+'"'

      assert (corr_list.file_name == ''), \
             'Correction list "'+f+'" has non-empty file name: '+ \
             str(lookup_table.file_name)

      assert (isinstance(corr_list,list)), \
             'Correction list "'+f+'" is not a list'

      corr_list.load(f)  # Load the table

      assert (len(corr_list) > 0), \
             'Correction list "'+f+'" is empty after load'

      elem_len = 9999  # Long elements first
      old_key = ''

      for (key,value) in corr_list:
        assert (isinstance(key,str)), \
               'Key in correction list "'+f+'" is not a string: '+str(key)
        assert (isinstance(value,str)), \
               'Value in correction list "'+f+'" is not a string: '+str(value)

        assert (len(key) <= elem_len), \
               'Correction list "'+f+'" element key "'+value+ \
               '" is longer than previous key: "'+old_key+'"'

        old_key = key
        elem_len = len(key)

# =============================================================================
# Start tests when called from command line

if (__name__ == "__main__"):
  unittest.main()  # Run all test

# =============================================================================
