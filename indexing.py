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
# The Original Software is: "indexing.py"
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

"""Module with classes for indexing/blocking.

   This module provides classes for building blocks and indices that will be
   used in the linkage process for comparing record pairs in an efficient way.

   Various derived classes are provided that implement different indexing
   techniques:

     FullIndex               Performs no indexing but does the full comparison
                             of all record pairs (quadratic complexity in the
                             size of the number of records in the two data
                             sets).
     BlockingIndex           The 'standard' blocking index used for record
                             linkage.
     SortingIndex            Based on a sliding window over the sorted values
                             of the index variable definitions - uses an
                             inverted index approach where keys are unique
                             index variable values.
     SortingArrayIndex       Based on a sliding window over the sorted values
                             of the index variable definitions - uses an array
                             based approach where all index variable values
                             (including duplicates) are stored.
     AdaptSortingIndex       An adaptive version of the array based sorted
                             index.
     QGramIndex              Allows for fuzzy indexing with 'overlapping'
                             blocks, like clustering.
     CanopyIndex             Based on TF-IDF/Jaccard and canopy clustering.
     StringMapIndex          Based on the string-map multi-dimensional mapping
                             algorithm combined with canopy clustering.
     SuffixArrayIndex        Based on a suffix array, resulting in a similar
                             approach as the SortingIndex (as blocks are
                             created by going through the sorted suffix array).
     RobustSuffixArrayIndex  Combines the suffix array approach with the sorted
                             window approach to overcome variations in the
                             index key values.
     BigMatchIndex           Based on the BigMatch program developed by the US
                             Census Bureau. This index can only handle linkages
                             (not deduplications). Based on the idea to only
                             index the smaller data set into an inverted index,
                             and then process each record from the larger data
                             set as it is read from file.
     DedupIndex              A index specialised for deduplications. It
                             performs the build(), compact() and run() in one
                             routine by reading of the file, building of the
                             index and comparisons of record pairs are all done
                             in the run() routine.

   When initialising an index its index variables have to be defined using the
   attribute 'index_def' (see more details below).

   Creating and using an index normally consists of the following three steps:
   - build    Read the records from the data sets and build the indexing data
              structures.
   - compact  Make the indexing data structures more compact (more efficient
              for accessing the record pairs).
   - run      Run the comparison step (i.e. compare record pairs) on the index.

   Main bottlenecks in the implemented indices are:
   - QGramIndex:    Creating of the sub-lists (recursively), especially for
                    long index values and low thresholds. Two different
                    functions are implemented and used for different threshold
                    values.
   - SortingIndex:  The way the inverted index data is combined when a sliding
                    window is created.
   - CanopyIndex:   Creating canopies, again especially for index values having
                    many q-grams and large data set (resulting in long inverted
                    index lists).
"""

# =============================================================================
# Import necessary modules (Python standard modules first, then Febrl modules)

import csv
import heapq
import gc
import logging
import math
import random
import shelve
import time

import auxiliary
import dataset
import encode

# =============================================================================

class Indexing:
  """Base class for indexing. Handles index initialisation, as well as saving
     and loading of indices to/from files.

     Indexing (or blocking as it is called in traditional record linkage) is
     used to reduce the large number of possible record pair comparisons by
     grouping together records into 'blocks' that are similar according to some
     criteria.

     Indices have to be defined when an index is initialised. Basically, an
     index is specific to two fields (one per data set). An index definition is
     made of a list containing one or more sub-lists each made of:

     1) The name of a field from data set 1 which is used for the index.

     2) The name of a field from data set 2 which is used for the index.

     3) A flag, if set to True and the index value contains more than one word
        (i.e. contains whitespace characters), then these words will be sorted
        first (before being possibly reversed and then encoded). For example,
        a value 'peter a. christen' would be sorted into 'a. christen peter'.
        If set to False, no such sorting will be done.

     4) A flag, if set to True the index variable (assumed to be a string) will
        be reversed before it is used within an index definition, if set to
        False it is not reversed.

     5) A positive integer specifying the maximum length (in characters) of the
        index variable; or None (no maximum length of the index variable).

     6) A list made of a function and its input parameters that will be applied
        to the data set field values to create an index variable. If this list
        is empty (or no list is given but None instead), then the field values
        will be taken directly as index
        variables.
        The function can be any function that has a string as its first input
        argument and returns a string.

     The final index variable values are the concatenated values (possibly with
     a separator string between as detailed below) of each of the index
     definitions.

     Here is an example index definition:

       index_def_1 = [['surname','sname',True,True,None,[encode.soundex]],
                      ['postcode','zipcode',False,False,3,[]]]

     For this index, index variables will be the concatenations of the Soundex
     encoding of the sorted (if containing more than one word) reversed
     surnames (fields 'surname' in data set 1 and field 'sname' in data set 2)
     and the first three characters (digits) from the fields 'postcode' from
     data set 1 and 'zipcode' from data set 2.

     So for a record ['peter','meier','2000','sydney'] in data set 1 and a
     record ['paul','meyer','2010','sydney south'] the following two index
     variable values will be constructed (and inserted into the index data
     structure):

       ['peter','meier','2000',sydney'] ->       'r500'+'200' -> 'r500200'
       ['paul','meyer','2010','sydney south'] -> 'r500'+'201' -> 'r500201'

     So these two records will (depending upon the index implementation) likely
     be put into different 'blocks' and therefore not being compared.

     Here is another example:

       index_def_2 = \
          [['year_of_birth','yob',False,False,None,[encode.get_substring,3,5]],
           ['locality_name','loc_name',True,True,None,[encode.dmetaphone, 4]]]

     In this example, the third and fourth year digits extracted from the
     fields 'year_of_birth' (from data set 1) and 'yob' (from data set 2) will
     be concatenated with the Double-Metaphone encoding (first 4 characters of
     the encoding only) of the sorted and reversed values from the fields
     'locality_name' (from data set 1) and 'loc_name' (from data set 2).

     The final index definition given to an index when it is initialised is a
     list of index definitions, for example:

       full_index = indexing.BlockingIndex(description = 'Example full index',
                                           dataset1 = example_csv_data_set,
                                           dataset2 = example_csv_data_set,
                                           index_def = [index_def_1,
                                                        index_def_2],
                                           index_sep_str = ' ',
                                           rec_comparator = example_rec_comp)

     In this example, two indices will be built according to the two index
     definitions provided.

     All index classed have the following instance variables, which can be set
     when an index is initialised:

       description      A string describing the index
       dataset1         The first data set the index is built for
       dataset2         The second data set the index is built for
       rec_comparator   The record comparator
       index_def        Definitions of the indexes as described above
       index_sep_str    A separator string which will be inserted if an index
                        is made of more than one value. Default is the empty
                        string ''.
       skip_missing     A flag, if set to True records which have empty index
                        variable values will be skipped over, if set to False a
                        record with an empty indexing variable value will be
                        included in the index as well. Default value is True
       log_funct        This can be a Python function or method which will log
                        (print or save to a file) a progress report message. It
                        is assumed that this function or method has one input
                        argument of type string (the message to be printed).
                        Default is None, in which case the normal logging
                        module will be used.
       progress_report  Can be set to a percentage number between 1 and 50 in
                        which case a progress report is logged during the
                        record pair comparison stage (in the run() method)
                        every selected percentage number. Default values is 10.
                        If set to None no progress report will be logged.
       weight_vec_file  If this is set to a string (assumed to be a file name)
                        then the weight vectors will be written into this file
                        and NOT stored in the weight_vector dictionary. Each
                        weight vector is a comma separated line containing:
                          rec_id1, rec_id2, weight1, weight2, ... weightN
                        Default value is None, in which case the weight vectors
                        will not be written into a file but returned as a
                        dictionary.

     Note that skip_missing cannot be set to False for certain index methods,
     see their documentation for more details.

     Both the data sets and the index definitions must be provided when a index
     is initialised.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, base_kwargs):
    """Constructor
    """

    # General attributes for all index implementations
    #
    self.description =     ''
    self.dataset1 =        None
    self.dataset2 =        None
    self.rec_comparator =  None  # A record comparator
    self.index_def =       None  # The index definitions
    self.index_sep_str =   ''    # Separator string for in-between index values
    self.skip_missing =    True
    self.progress_report = 10
    self.log_funct =       None
    self.weight_vec_file = None

    self.index_def_proc = None        # Processed version of the index
                                      # definition for faster access to field
                                      # values
    self.index1 = {}                  # The index data structure for data set 1
    self.index2 = {}                  # The index data structure for data set 2
    self.index1_shelve_name = None    # If the index data structure for data
                                      # set 1 is to be file (shelve) based this
                                      # will be it's file name
    self.index2_shelve_name = None    # Same for data set 2
    self.rec_cache1 = {}              # A dictionary containing all records
                                      # from data set 1 with only the fields
                                      # needed for field comparisons
    self.rec_cache2 = {}              # Same for data set 2
    self.rec_cache1_file_name = None  # If the record cache for data sets 1
                                      # should be file (shelve) based this will
                                      # be it's file name
    self.rec_cache2_file_name = None  # Same for data sets 2
    self.num_rec_pairs = None         # The number of record pairs that will be
                                      # compared when the run() method is
                                      # called
    self.do_deduplication = None      # A Flag, set to True if both data sets
                                      # are the same
    self.comp_field_used1 = []        # A list of the fields in data set 1 as
                                      # used by the record comparator (column
                                      # indices)
    self.comp_field_used2 = []        # Same for data set 2
    self.rec_length_cache = {}        # Used in lenth filtering in run() method

    # Process base keyword arguments (all data set specific keywords were
    # processed in the derived class constructor)
    #
    for (keyword, value) in base_kwargs.items():

      if (keyword.startswith('desc')):
        auxiliary.check_is_string('description', value)
        self.description = value

      elif (keyword == 'dataset1'):
        auxiliary.check_is_list('Dataset 1 field list', value.field_list)
        self.dataset1 = value
      elif (keyword == 'dataset2'):
        auxiliary.check_is_list('Dataset 2 field list', value.field_list)
        self.dataset2 = value

      elif (keyword.startswith('index_d')):
        auxiliary.check_is_list('index_def', value)
        self.index_def = value

      elif (keyword.startswith('index_sep')):
        auxiliary.check_is_string('index_sep_str', value)
        self.index_sep_str = value

      elif (keyword.startswith('index1_she')):
        auxiliary.check_is_string('index1_shelve_name', value)
        self.index1_shelve_name = value
      elif (keyword.startswith('index2_she')):
        auxiliary.check_is_string('index2_shelve_name', value)
        self.index2_shelve_name = value

      elif (keyword.startswith('rec_cache1_f')):
        auxiliary.check_is_string('rec_cache1_file_name', value)
        self.rec_cache1_file_name = value
      elif (keyword.startswith('rec_cache2_f')):
        auxiliary.check_is_string('rec_cache2_file_name', value)
        self.rec_cache2_file_name = value

      elif (keyword.startswith('rec_com')):
        self.rec_comparator = value

      elif (keyword.startswith('skip')):
        auxiliary.check_is_flag('skip_missing', value)
        self.skip_missing = value

      elif (keyword.startswith('progr')):
        if (value == None):
          self.progress_report = value
        else:
          auxiliary.check_is_integer('progress_report', value)
          if ((value < 1) or (value > 50)):
            logging.exception('Illegal value for progress report, must be ' + \
                              'False or between 1 and 50: "%s"' % \
                              (str(value)))
            raise Exception
          self.progress_report = value

      elif (keyword.startswith('log_f')):
        auxiliary.check_is_function_or_method('log_funct', value)
        self.log_funct =  value

      elif (keyword.startswith('weight_v')):
        if (value != None):
          auxiliary.check_is_string('weight_vec_file', value)
        self.weight_vec_file = value

      else:
        logging.exception('Illegal constructor argument keyword: '+keyword)
        raise Exception

    # Check if the needed attributes are set  - - - - - - - - - - - - - - - - -
    #
    auxiliary.check_is_list('index_def', self.index_def)
    auxiliary.check_is_not_none('rec_comparator', self.rec_comparator)
    auxiliary.check_is_not_none('Dataset 1', self.dataset1)
    auxiliary.check_is_not_none('Dataset 2', self.dataset2)
    auxiliary.check_is_list('Dataset 1 field list', self.dataset1.field_list)
    auxiliary.check_is_list('Dataset 2 field list', self.dataset2.field_list)

    # Check if the data sets in the record comparator are the same as the ones
    # give in the index
    #
    if (self.rec_comparator.dataset1 != self.dataset1):
      logging.exception('Data set 1 in record comparator is different from' + \
                        ' data set 1 in index: "%s" / "%s"' % \
                        (self.rec_comparator.dataset1, self.dataset1))
      raise Exception
    if (self.rec_comparator.dataset2 != self.dataset2):
      logging.exception('Data set 2 in record comparator is different from' + \
                        ' data set 2 in index: "%s" / "%s"' % \
                        (self.rec_comparator.dataset2, self.dataset2))
      raise Exception

    # Check if a deduplication will be done or a linkage - - - - - - - - - - -
    #
    if (self.dataset1 == self.dataset2):
      self.do_deduplication = True
    else:
      self.do_deduplication = False

    # If indices or record caches are file (shelve) based open these shelves -
    #
    if (self.index1_shelve_name != None):
      self.index1 = self.__open_shelve_file__(self.index1_shelve_name)
    if (self.index2_shelve_name != None):
      self.index2 = self.__open_shelve_file__(self.index2_shelve_name)
    if (self.rec_cache1_file_name != None):
      self.rec_cache1 = self.__open_shelve_file__(self.rec_cache1_file_name)
    if (self.rec_cache2_file_name != None):
      self.rec_cache2 = self.__open_shelve_file__(self.rec_cache2_file_name)

    # Extract the field names from the two data set field name lists - - - - -
    #
    dataset1_field_names = []
    dataset2_field_names = []

    for (field_name, field_data) in self.dataset1.field_list:
      dataset1_field_names.append(field_name)

    for (field_name, field_data) in self.dataset2.field_list:
      dataset2_field_names.append(field_name)

    # Get the column indices of the fields used in the record comparator - - -
    #
    for (comp,f_ind1,f_ind2) in self.rec_comparator.field_comparison_list:
      if (f_ind1 not in self.comp_field_used1):
        self.comp_field_used1.append(f_ind1)
      if (f_ind2 not in self.comp_field_used2):
        self.comp_field_used2.append(f_ind2)
    self.comp_field_used1.sort()
    self.comp_field_used2.sort()

    # Check if definition of indices is correct and fields are in the data sets
    #
    self.index_def_proc = []  # Checked and processed index definitions will be
                              # added here

    for index_def_list in self.index_def:

      auxiliary.check_is_list('Index definition list "%s"' % \
                              (str(index_def_list)), index_def_list)

      index_def_list_proc = []

      for index_def in index_def_list:

        auxiliary.check_is_list('Index definition "%s"' % \
                                (str(index_def)), index_def)
        if (index_def == []):
          logging.info('Empty index definition given: %s' % \
                       (str(self.index_def)))

        # Check the two field names
        #
        field_name1 = index_def[0]
        field_name2 = index_def[1]

        if (field_name1 not in dataset1_field_names):
          logging.exception('Field "%s" is not in data set 1 field name ' \
                            % (field_name1) + 'list: %s' % \
                            (str(self.dataset1.field_list)))
          raise Exception
        field_index1 = dataset1_field_names.index(field_name1)

        if (field_name2 not in dataset2_field_names):
          logging.exception('Field "%s" is not in data set 2 field name ' \
                            % (field_name2) + 'list: %s' % \
                            (str(self.dataset2.field_list)))
          raise Exception
        field_index2 = dataset2_field_names.index(field_name2)

        index_def_proc = [field_index1,field_index2] # Processed index def.

        # Check if sort words flag is True or False
        #
        auxiliary.check_is_flag('Sort words flag', index_def[2])
        index_def_proc.append(index_def[2])

        # Check if reverse flag is True or False
        #
        auxiliary.check_is_flag('Reverse flag', index_def[3])
        index_def_proc.append(index_def[3])

        # Check maximum length is a positive integer or None
        #
        if (index_def[4] == None):
          index_def_proc.append(None)
        else:
          auxiliary.check_is_integer('Maximum length', index_def[4])
          auxiliary.check_is_positive('Maximum length', index_def[4])
          index_def_proc.append(index_def[4])

        # Check function definition
        #
        if ((index_def[5] != None) and (len(index_def[5]) > 0)):
          index_funct_def = index_def[5]
          auxiliary.check_is_function_or_method('Function "%s"' % \
                           (index_funct_def[0]), index_funct_def[0])
          index_def_proc.append(index_funct_def)
        else:
          index_def_proc.append(None)

        index_def_list_proc.append(index_def_proc)

      self.index_def_proc.append(index_def_list_proc)

    assert len(self.index_def) == len(self.index_def_proc)

    self.status = 'initialised'  # Status of the index (used by save and load
                                 # methods)

  # ---------------------------------------------------------------------------

  def build(self):
    """Method to build an index data structure.
       See implementations in derived classes for details.
    """

    logging.exception('Override abstract method in derived class')
    raise Exception

  # ---------------------------------------------------------------------------

  def __records_into_inv_index__(self):
    """Load the records from the data sets and put them into an inverted index
       data structure.

       If both data sets are the same (a deduplication) only insert records
       from one data set.

       This method builds an inverted index (one per index definition) as a
       Python dictionary with the keys being the indexing values (as returned
       by the _get_index_values__() method.
    """

    logging.info('Started to build inverted index:')

    num_indices = len(self.index_def)

    # Index data structure for blocks is one dictionary per index - - - - - - -
    #
    for i in range(num_indices):
      self.index1[i] = {}  # Index for data set 1
      self.index2[i] = {}  # Index for data set 2

    get_index_values_funct = self.__get_index_values__  # Shorthands
    skip_missing =           self.skip_missing

    # A list of data structures needed for the build process:
    # - the index data structure (dictionary)
    # - the record cache
    # - the data set to be read
    # - the comparison fields which are used
    # - a list index (0 for data set 1, 1 for data set 2)
    #
    build_list = [(self.index1, self.rec_cache1, self.dataset1,
                   self.comp_field_used1, 0)] # For data set 1

    if (self.do_deduplication == False):  # If linkage append data set 2
      build_list.append((self.index2, self.rec_cache2, self.dataset2,
                   self.comp_field_used2, 1))

    # Reading loop over all records in one or both data set(s) - - - - - - - -
    #
    for (index,rec_cache,dataset,comp_field_used_list,ds_index) in build_list:

      # Calculate a counter for the progress report
      #
      if (self.progress_report != None):
        progress_report_cnt = max(1, int(dataset.num_records / \
                                     (100.0 / self.progress_report)))
      else:  # So no progress report is being logged
        progress_report_cnt = dataset.num_records + 1

      start_time = time.time()

      rec_read = 0  # Number of records read from data set

      for (rec_ident, rec) in dataset.readall(): # Read all records in data set

        # Extract record fields needed for comparisons (set all others to '')
        #
        comp_rec = []

        field_ind = 0
        for field in rec:
          if (field_ind in comp_field_used_list):
            comp_rec.append(field.lower())  # Make them lower case
          else:
            comp_rec.append('')
          field_ind += 1

        rec_cache[rec_ident] = comp_rec  # Put into record cache

        # Now get the index variable values for this record - - - - - - - - - -
        #
        rec_index_val_list = get_index_values_funct(rec, ds_index)

        for i in range(num_indices):  # Put record identifier into all indices

          this_index = index[i]  # Shorthand

          block_val = rec_index_val_list[i]

          if ((block_val != '') or (skip_missing == False)):
            block_val_rec_list = this_index.get(block_val, [])
            block_val_rec_list.append(rec_ident)
            this_index[block_val] = block_val_rec_list

        rec_read += 1

        if ((rec_read % progress_report_cnt) == 0):
          self.__log_build_progress__(rec_read,dataset.num_records,start_time)

      used_sec_str = auxiliary.time_string(time.time()-start_time)
      rec_time_str = auxiliary.time_string((time.time()-start_time) / \
                                           dataset.num_records)
      logging.info('Read and indexed %d records in %s (%s per record)' % \
                   (dataset.num_records, used_sec_str, rec_time_str))
      logging.info('')

  # ---------------------------------------------------------------------------
  # Get sub-list functions are used for the q-gram and BigMatch index

  def __get_sublists1__(self, in_list, min_len):
    """An iterative method to compute all combinations of sub-lists of the list
       'in_list' of lengths down to 'min_len'.

       This routine seems to be faster for larger 'min_len' values (compared to
       length of the input list), i.e. which would correspond to less deep
       recursive levels (that correspond to higher threshold values) compared
       to the recursive method using sets given below.
    """

    all_list = [in_list]
    in_len =   len(in_list)

    if (in_len > min_len):

      for i in xrange(in_len):

        sub_list = in_list[:i]+in_list[i+1:]

        if (sub_list not in all_list):
          all_list.append(sub_list)

        if ((in_len-1) > min_len):
          for sub_sub_list in self.__get_sublists1__(sub_list, min_len):

            if (sub_sub_list not in all_list):
              all_list.append(sub_sub_list)

    return all_list

  # ---------------------------------------------------------------------------

  def __get_sublists2__(self, in_list, min_len):
    """Method to recursively compute all combinations of sub-lists of the list
       'in_list' of lengths down to 'min_len' using a Python generator.

       This routine seems to be faster for deeper recursive levels (that
       correspond to lower threshold values) compared to the iterative method
       using sets given above.
    """

    unique_combinations_funct = self.__unique_combinations2__  # Shorthand


    all_list = [in_list]
    in_len =   len(in_list)

    for i in xrange(in_len-1, min_len-1, -1):

      for sub_list in unique_combinations_funct(in_list, i):
        if (sub_list not in all_list):
          all_list.append(sub_list)

    return all_list

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  def __unique_combinations__(self, in_list, n):
    """Get unique combinations of length 'n' of the given input list. Based on
       a recipe from Activestate Python cookbook, see:

       http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/190465
    """

    unique_combinations_funct = self.__unique_combinations__  # Shorthand

    if (n == 0):
      yield []
    else:
      for i in xrange(len(in_list)):
        for sub_list in unique_combinations_funct(in_list[i+1:],n-1):
          yield [in_list[i]]+sub_list

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  def __unique_combinations2__(self, in_list, n):
    """Another version from Activestate Python cookbook, see:

       http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/66465

       Seems to be of similar speed than the previous function.
    """

    unique_combinations_funct = self.__unique_combinations2__  # Shorthand

    for i in xrange(len(in_list)):
      if (n == 1):
        yield in_list[i]
      else:
        for sub_list in unique_combinations_funct(in_list[i+1:], n-1):
          yield in_list[i] + sub_list

  # ---------------------------------------------------------------------------

  def compact(self):
    """Method to compact an index data structure.
       See implementations in derived classes for details.
    """

    logging.exception('Override abstract method in derived class')
    raise Exception

  # ---------------------------------------------------------------------------

  def __dedup_rec_pairs__(self, rec_id_list, rec_pair_dict):
    """Create record pairs for a deduplication using the given record
       identifier list and insert them into the given record pair dictionary.

       This version does not modify the input record identifier list. It does
       create a local copy of the record identifer list which is then sorted.
    """

    rec_cnt = 1  # Counter for second record identifier

    this_rec_id_list = rec_id_list[:]
    this_rec_id_list.sort()
    for rec_ident1 in this_rec_id_list:

      rec_ident2_set = rec_pair_dict.get(rec_ident1, set())

      for rec_ident2 in this_rec_id_list[rec_cnt:]:

        assert rec_ident1 != rec_ident2

        rec_ident2_set.add(rec_ident2)

      rec_pair_dict[rec_ident1] = rec_ident2_set

      rec_cnt += 1

    del this_rec_id_list

  # ---------------------------------------------------------------------------

  def __link_rec_pairs__(self, rec_id_list1, rec_id_list2, rec_pair_dict):
    """Create record pairs for a linkage using the given two record identifier
       lists and insert them into the given record pair dictionary.
    """

    for rec_ident1 in rec_id_list1:
      for rec_ident2 in rec_id_list2:

        rec_ident2_set = rec_pair_dict.get(rec_ident1, set())
        rec_ident2_set.add(rec_ident2)
        rec_pair_dict[rec_ident1] = rec_ident2_set

  # ---------------------------------------------------------------------------

  def run(self):
    """Run the record pair comparison accoding to the index.
       See implementations in derived classes for details.
    """

    logging.exception('Override abstract method in derived class')
    raise Exception

  # ---------------------------------------------------------------------------

  def __get_field_names_list__(self):
    """Returns the list of all field comparison descriptions, to be used when
       a weight vector file is written.
    """

    field_names_list = []
    for field_comp_tuple in self.rec_comparator.field_comparator_list:
        field_names_list.append(field_comp_tuple[0].description)
    return field_names_list

  # ---------------------------------------------------------------------------

  def __compare_rec_pairs_from_dict__(self, length_filter_perc = None,
                                      cut_off_threshold = None):
    """This method compares all the records pairs in the record pair dictionary
       and puts the resulting weight vectors into a dictionary which is then
       returned.

       If the argument 'length_filter_perc' is set to a percentage value
       between 0.0 and 100.0 the lengths of the concatenated record values in
       the fields used for comparison will be calculated, and if the difference
       is larger than the given percentage value the record pair will not be
       compared. For example, if a record A has values with a total length of
       12 characters and a record B values with a length of 16 characters, then
       their percentage difference will be:

          |12-16| / (max(12,16) = 4 / 16 = 0.25 (25%)

       So if the value of 'length_filter_perc' is set to 20 this record pair
       will be filtered out and not compared.

       Default value for 'length_filter_perc' is None, which means no length
       filtering will be performed. For more information about filtering in the
       record linkage process please refer to:

       - Adaptive Filtering for Efficient Record Linkage
         Lifang Gu and Rohan Baxter,
         SIAM Data Mining conference, 2004.

       The second argument 'cut_off_threshold' can be set to a numerical value,
       in which case all compared record pairs with a summed weight vector
       value less than this threshold will not be stored in the weight vector
       dictionary. Default value for 'cut_off_threshold' is None, which means
       all compared record pairs will be stored in the weight vector
       dictionary.
    """

    # Check if weight vector file should be written - - - - - - - - - - - - - -
    #
    if (self.weight_vec_file != None):
      try:
        weight_vec_fp = open(self.weight_vec_file, 'w')
      except:
        logging.exception('Cannot write weight vector file: %s' % \
                          (self.weight_vec_file))
        raise Exception
      weight_vec_writer = csv.writer(weight_vec_fp)

      # Write header line with descriptions of field comparisons
      #
      weight_vec_header_line = ['rec_id1', 'rec_id2'] + \
                                self.__get_field_names_list__()
      weight_vec_writer.writerow(weight_vec_header_line)

    # Calculate a counter for the progress report - - - - - - - - - - - - - - -
    #
    if (self.progress_report != None):
      progress_report_cnt = max(1, int(self.num_rec_pairs / \
                                   (100.0 / self.progress_report)))
    else:  # So no progress report is being logged
      progress_report_cnt = self.num_rec_pairs + 1

    weight_vec_dict = {}  # Dictionary with calculated weight vectors
    comp_done =       0   # Number of comparisons done

    rec_cache1 =       self.rec_cache1  # Shorthands to make program faster
    rec_pair_dict =    self.rec_pair_dict
    rec_comp =         self.rec_comparator.compare
    rec_length_cache = self.rec_length_cache

    # Check length filter and cut-off threshold arguments - - - - - - - - - - -
    #
    if (length_filter_perc != None):
      auxiliary.check_is_percentage('Length filter percentage',
                                    length_filter_perc)
      logging.info('  Length filtering set to %.1f%%' % (length_filter_perc))
      length_filter_perc /= 100.0  # Normalise

    if (cut_off_threshold != None):
      auxiliary.check_is_number('Cut-off threshold', cut_off_threshold)
      logging.info('  Cut-off threshold set to: %.2f' % (cut_off_threshold))

    num_rec_pairs_filtered =    0  # Count number of removed record pairs
    num_rec_pairs_below_thres = 0

    # Set shorthand depending upon deduplication or linkage - - - - - - - - - -
    #
    if (self.do_deduplication == True):  # A deduplication run
      rec_cache2 = self.rec_cache1
    else:
      rec_cache2 = self.rec_cache2

    start_time = time.time()

    for rec_ident1 in rec_pair_dict:

      rec1 = rec_cache1[rec_ident1]  # Get the actual first record

      if (length_filter_perc != None):
        rec1_len = len(''.join(rec1))  # Get length in characters for record

      for rec_ident2 in rec_pair_dict[rec_ident1]:

        rec2 = rec_cache2[rec_ident2]  # Get actual second record

        do_comp = True  # Flag, specify if comparison should be done

        if (length_filter_perc != None):
          if (rec_ident2 in rec_length_cache):  # Length is cached
            rec2_len = rec_length_cache[rec_ident2]
          else:
            rec2_len = len(''.join(rec2))
            rec_length_cache[rec_ident2] = rec2_len

          perc_diff = float(abs(rec1_len - rec2_len)) / max(rec1_len, rec2_len)

          if (perc_diff > length_filter_perc):
            do_comp = False  # Difference too large, don't do comparison
            num_rec_pairs_filtered += 1

        if (do_comp == True):
          w_vec = rec_comp(rec1, rec2)  # Compare them

          if (cut_off_threshold == None) or (sum(w_vec) >= cut_off_threshold):

            # Put result into weight vector dictionary
            #
            if (self.weight_vec_file == None):
              weight_vec_dict[(rec_ident1, rec_ident2)] = w_vec
            else:
              weight_vec_writer.writerow([rec_ident1, rec_ident2]+w_vec)

          else:
            num_rec_pairs_below_thres += 1

        comp_done += 1  # Count all record pair comparisons (even if not done)

        if ((comp_done % progress_report_cnt) == 0):
          self.__log_comparison_progress__(comp_done, start_time)

    used_sec_str = auxiliary.time_string(time.time()-start_time)
    rec_time_str = auxiliary.time_string((time.time()-start_time) / \
                                         self.num_rec_pairs)
    logging.info('Compared %d record pairs in %s (%s per pair)' % \
                 (self.num_rec_pairs, used_sec_str,rec_time_str))
    if (length_filter_perc != None):
      logging.info('  Length filtering (set to %.1f%%) filtered %d record ' % \
                   (length_filter_perc*100, num_rec_pairs_filtered) + 'pairs')
    if (cut_off_threshold != None):
      logging.info('  %d record pairs had summed weights below threshold ' % \
                   (num_rec_pairs_below_thres) + '%.2f' % (cut_off_threshold))

    memory_usage_str = auxiliary.get_memory_usage()
    if (memory_usage_str != None):
      logging.info('  '+memory_usage_str)

    if (self.weight_vec_file == None):
      return [self.__get_field_names_list__(), weight_vec_dict]
    else:
      weight_vec_fp.close()
      return None

  # ---------------------------------------------------------------------------

  def __find_closest__(self, sorted_list, elem):
    """Binary search of the given element 'elem' in the given sorted list, and
       return index of exact match or closest match (before where the element
       would be).

       Used in BigMatch and Dedup indexing methods.
    """

    start = 0
    end =   len(sorted_list)-1

    while 1:
      if (end < start):
        return end  # Not found

      middle = (start + end) / 2
      middle_val = sorted_list[middle]  # Get the element from the middle

      if (middle_val == elem):  # Found the key value
        return middle
      elif (elem < middle_val):
        end = middle - 1
      else:
        start = middle + 1

  # ---------------------------------------------------------------------------

  def load(self, index_file_name):
    """Load a previously saved index from a binary file.
    """

    pass  # To be done

  # ---------------------------------------------------------------------------

  def save(self, index_file_name):
    """Save an index into a binary file.
    """

    pass  # To be done

    # have to check that data sets are the same (have same field lists)

  # ---------------------------------------------------------------------------

  def __get_index_values__(self, rec, data_set_num):
    """For the given record (list of fields) extract and produce the indexing
       values. Returns a list with the indexing variable values (one per index
       definition).

       The data set number can be 0 (if the record is from the first data set)
       or 1 (if it is from the second data set).
    """

    index_var_values = []

    sep_str = self.index_sep_str

    assert (data_set_num == 0) or (data_set_num == 1)

    # Go through the index definitions and extract and process field values - -
    #
    for index_def_list in self.index_def_proc:

      index_val_list = []

      for index_def in index_def_list:  # Loop over all the index' definitions

        field_col = index_def[data_set_num]  # Column of the field to extract

        if (field_col >= len(rec)):
          field_val = ''
        else:
          field_val = rec[field_col].lower()

        if (field_val != ''):  # Field value is not empty

          # Check for sorting of words
          #
          if ((' ' in field_val) and (index_def[2] == True)):
            word_list = field_val.split()
            word_list.sort()
            field_val = ' '.join(word_list)

          if (index_def[3] == True):  # Reverse the index value
            field_val = field_val[::-1]

          funct_def = index_def[5]

          if (funct_def != None):  # There is a function defined for this index
            funct_call =    funct_def[0]  # The function itself
            num_funct_arg = len(funct_def)

            if (num_funct_arg == 1):  # No arguments
              funct_val = funct_call(field_val)
            elif (num_funct_arg == 2):  # One argument
              funct_val = funct_call(field_val, funct_def[1])
            elif (num_funct_arg == 3):  # Two arguments
              funct_val = funct_call(field_val, funct_def[1], funct_def[2])
            elif (num_funct_arg == 4):  # Three arguments
              funct_val = funct_call(field_val, funct_def[1], funct_def[2],
                                     funct_def[3])
            else:
              logging.exception('Too many arguments for function call: %s' % \
                                (str(funct_def)))
              raise Exception
          else:
           funct_val = field_val  # No function applied to the field value

          if (index_def[4] != None):  # There is  maximum length
            funct_val = funct_val[:index_def[4]]

          index_val_list.append(funct_val)

      # Make it a string and add to list of index values
      #
      index_val = sep_str.join(index_val_list)

      index_var_values.append(index_val)

    assert len(index_var_values) == len(self.index_def_proc)

    return index_var_values

  # ---------------------------------------------------------------------------

  def __open_shelve_file__(self, shelve_file_name):
    """Open a shelve with the given file name, and clear all it's content.

       Return the shelve.
    """

    try:
      this_shelve = shelve.open(shelve_file_name)
    except:
      logging.exception('Cannot open shelve: "%s"' % str(shelve_file_name))
      raise Exception

    for rec_key in this_shelve:  # Now clear the shelve and remove all entries

        del this_shelve[rec_key]

    this_shelve.sync()  # And make sure the database is updated

    return this_shelve

  # ---------------------------------------------------------------------------

  def __log_build_progress__(self, records_read, num_records, start_time):
    """Create a log message for the number of records read and indexed so far,
       the time used, and an estimation of much longer it will take.
    """

    used_time = time.time() - start_time
    perc_done = 100.0 * records_read / num_records
    rec_time  = used_time / records_read  # Time per record read and indexed
    togo_time = (num_records - records_read) * rec_time

    used_sec_str = auxiliary.time_string(used_time)
    rec_sec_str =  auxiliary.time_string(rec_time)
    togo_sec_str = auxiliary.time_string(togo_time)

    log_str = 'Read and indexed %d of %d records (%d%%) in %s (%s per' % \
              (records_read, num_records, round(perc_done), used_sec_str,
              rec_sec_str)+' record), estimated %s until finished.' % \
              (togo_sec_str)
    logging.info(log_str)

    if (self.log_funct != None):
      self.log_funct(log_str)

    memory_usage_str = auxiliary.get_memory_usage()
    if (memory_usage_str != None):
      logging.info('    '+memory_usage_str)

  # ---------------------------------------------------------------------------

  def __log_comparison_progress__(self, comparison_count, start_time):
    """Create a log message for the number of comparison done so far, the time
       used, and an estimation of much longer it will take.
    """

    used_time = time.time() - start_time
    perc_done = int(round(100.0 * comparison_count / self.num_rec_pairs))
    rec_time  = used_time / comparison_count  # Time per record pair comparison
    togo_time = (self.num_rec_pairs - comparison_count) * rec_time

    used_sec_str = auxiliary.time_string(used_time)
    rec_sec_str =  auxiliary.time_string(rec_time)
    togo_sec_str = auxiliary.time_string(togo_time)

    log_str = 'Compared %d of %d record pairs (%d%%) in %s (%s per ' % \
              (comparison_count, self.num_rec_pairs, perc_done, used_sec_str,
              rec_sec_str)+'pair), estimated %s until finished.' % \
              (togo_sec_str)
    logging.info(log_str)

    if (self.log_funct != None):
      self.log_funct(log_str)

    memory_usage_str = auxiliary.get_memory_usage()
    if (memory_usage_str != None):
      logging.info('    '+memory_usage_str)

  # ---------------------------------------------------------------------------

  def log(self, instance_var_list = None):
    """Write a log message with the basic index instance variables plus the
       instance variable provided in the given input list (assumed to contain
       pairs of names (strings) and values).
    """

    logging.info('')
    logging.info('Index:                    "%s"' % (self.description))
    logging.info('  Index status:           %s' % (self.status))
    logging.info('  Data set 1:             %s' % (self.dataset1.description))
    logging.info('  Data set 2:             %s' % (self.dataset2.description))
    if (self.do_deduplication == True):
      logging.info('  Data sets are the same: Deduplication')
    else:
      logging.info('  Data sets differ:       Linkage')
    logging.info('  Record comparator:      %s' % \
                 (self.rec_comparator.description))
    logging.info('    Used fields indices from data set 1: %s' % \
                 (str(self.comp_field_used1)))
    logging.info('    Used fields indices from data set 2: %s' % \
                 (str(self.comp_field_used2)))
    logging.info('  Skip missing:           %s' % (str(self.skip_missing)))
    logging.info('  Index separator string: "%s"' % (self.index_sep_str))

    if (self.num_rec_pairs == None):
      logging.info('  Number of record pairs: Not known yet')
    else:
      logging.info('  Number of record pairs: %d' % (self.num_rec_pairs))
    if (self.progress_report == None):
      logging.info('  No progress reported.')
    else:
      logging.info('  Progress report every:  %d%%' % (self.progress_report))
    if (self.weight_vec_file != None):
      logging.info('  Weight vectors will be written into: %s' % \
                   (self.weight_vec_file))

    if (instance_var_list != None):
      logging.info('  Index specific variables:')

      max_name_len = 0
      for (name, value) in instance_var_list:
        max_name_len = max(max_name_len, len(name))

      for (name, value) in instance_var_list:
        pad_spaces = (max_name_len-len(name))*' '
        logging.info('    %s %s' % (name+':'+pad_spaces, str(value)))

    if (self.index_def != []):  # Index definitions are not empty - - - - - - -
      logging.info('  Index definitions:')
      for i in range(len(self.index_def)):
        logging.info('    Index %d:' % (i))

        for j in range(len(self.index_def[i])):
          index_def = self.index_def[i][j]
          index_def_proc = self.index_def_proc[i][j]

          logging.info('      Data set 1 and 2 field names: "%s" (at column' \
                       % (index_def[0]) + ' %d) and "%s" (at column %d)' % \
                       (index_def_proc[0], index_def[1], index_def_proc[1]))
          logging.info('        Sort words flag set to %s' % \
                       (str(index_def[2])))
          logging.info('        Reverse flag set to %s' % (str(index_def[3])))

          if (index_def_proc[4] == None):
            logging.info('        No maximum length for index variable ' + \
                         'specified')
          else:
            logging.info('        Maximum length for index variable: %d' % \
                         (index_def_proc[4]))
          if (index_def_proc[5] == None):
            logging.info('        No function specified')
          else:
            funct_str = '        Function: %s' % (str(index_def_proc[5][0]))
            param_str = ''
            if (len(index_def_proc[5]) > 1):
              param_str += ' with parameter(s): '
              for p in index_def_proc[5][1:]:
                param_str += str(p)+', '
            logging.info(funct_str+param_str[:-2])

  # ---------------------------------------------------------------------------

  def get_index_stats(self):
    """Extract and log information about the index.
    """

    pass

# =============================================================================

class FullIndex(Indexing):
  """Class that implements a 'full' index which does not generate an index but
     will compare all record pairs from the two data sets.

     Note this will be very very slow for large data sets. The number of record
     pairs to be compared corresponds to the number of records in the first
     data set times the number of records in the second data set.

     Note that the index definition will be ignored, as all record pairs will
     be compared. Neither will length filtering nor the cut-off threshold be
     applied.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """No additional attributes needed, just call base class constructor.
    """

    Indexing.__init__(self, kwargs)  # Initialise base class

    num_rec1 = self.dataset1.num_records
    num_rec2 = self.dataset2.num_records

    if (num_rec1 < num_rec2):  # Make references to small and large data sets
      self.ds_swapped = False
      self.small_dataset = self.dataset1
      self.large_dataset = self.dataset2
    else:  # This will also hold for a deduplication
      self.ds_swapped = True
      self.small_dataset = self.dataset2
      self.large_dataset = self.dataset1

    if (self.do_deduplication == True):
      assert num_rec1 == num_rec2
      self.num_rec_pairs = num_rec1*(num_rec2-1)/2

    else:
      self.num_rec_pairs = num_rec1*num_rec2

    self.log()

  # ---------------------------------------------------------------------------

  def build(self):
    """Method to build an index data structure.

       If the smaller data set is not memory based it will be loaded into
       memory (as a dictionary) to make the comparisons faster.

       Otherwise nothing needs to be done.
    """

    logging.info('')
    logging.info('Built full index: "%s"' % (self.description))

    start_time = time.time()

    small_data_set_dict = {}

    # Copy all records into the memory based dictionary - - - - - - - - - - - -
    #
    for (rec_ident, rec_list) in self.small_dataset.readall():
      if (rec_ident in small_data_set_dict):
        logging.exception('Duplicate record identifier "%s" in small data ' % \
                          (rec_ident) + 'set')

      rec_list_lower = []  # Make all values lowercase
      for rec_val in rec_list:
        rec_list_lower.append(rec_val.lower())

      small_data_set_dict[rec_ident] = rec_list_lower

    # assert self.small_dataset.num_records == len(small_data_set_dict)

    self.small_data_set_dict = small_data_set_dict

    logging.info('Built full index in %s' % \
                 (auxiliary.time_string(time.time()-start_time)))

    memory_usage_str = auxiliary.get_memory_usage()
    if (memory_usage_str != None):
      logging.info('  '+memory_usage_str)

    self.status = 'built'  # Update index status

  # ---------------------------------------------------------------------------

  def compact(self):
    """Method to compact an index data structure.

       Nothing needs to be done.
    """

    logging.info('')
    logging.info('Compact full index: "%s"' % (self.description))
    logging.info('  Nothing needs to be done.')

    self.status = 'compacted'  # Update index status

  # ---------------------------------------------------------------------------

  def run(self, length_filter_perc = None, cut_off_threshold = None):
    """Iterate over all blocks in the index.

       Compare all record pairs (quadratic complexity), and return a weight
       vector dictionary with keys made of a tuple (record identifier 1,
       record identifier 2), and corresponding values the comparison weights.

       Length filtering and cut-off threshold will both be ignored.
    """

    logging.info('Started comparison of %d record pairs' % \
                 (self.num_rec_pairs))
    if (self.log_funct != None):
      self.log_funct('Started comparison of %d record pairs' % \
                     (self.num_rec_pairs))

    # Check if index has been compacted - - - - - - - - - - - - - - - - - - - -
    #
    if (self.status != 'compacted'):  # Not necessary, but keep code consistent
      logging.exception('Index "%s" has not been compacted, running ' % \
                        (self.description)+'comparisons not possible')
      raise Exception

    # Check if weight vector file should be written - - - - - - - - - - - - - -
    #
    if (self.weight_vec_file != None):
      try:
        weight_vec_fp = open(self.weight_vec_file, 'w')
      except:
        logging.exception('Cannot write weight vector file: %s' % \
                          (self.weight_vec_file))
        raise Exception
      weight_vec_writer = csv.writer(weight_vec_fp)

      # Write header line with descriptions of field comparisons
      #
      weight_vec_header_line = ['rec_id1', 'rec_id2'] + \
                                self.__get_field_names_list__()
      weight_vec_writer.writerow(weight_vec_header_line)

    start_time = time.time()

    # Calculate a counter for the progress report - - - - - - - - - - - - - - -
    #
    if (self.progress_report != None):
      progress_report_cnt = max(1, int(self.num_rec_pairs / \
                                   (100.0 / self.progress_report)))
    else:  # So no progress report is being logged
      progress_report_cnt = self.num_rec_pairs + 1

    weight_vec_dict = {}  # Dictionary with calculated weight vectors

    comp_done = 0  # Counter for the number of comparisons done so far

    compare_funct =       self.rec_comparator.compare  # Shorthands
    small_data_set_dict = self.small_data_set_dict
    rec_length_cache =    self.rec_length_cache

    if (self.do_deduplication == True):  # A deduplication run - - - - - - - -

      # Need a sorted list of all record identifiers
      #
      small_data_set_rec_id_list = small_data_set_dict.keys()
      small_data_set_rec_id_list.sort()

      rec_cnt = 1  # Counter for second record identifier

      for rec_ident1 in small_data_set_rec_id_list:
        rec1 = small_data_set_dict[rec_ident1]  # Get values of first record

        for rec_ident2 in small_data_set_rec_id_list[rec_cnt:]:

          assert rec_ident1 != rec_ident2  # Make sure they are different

          rec2 = small_data_set_dict[rec_ident2]  # Get values of second record

          w_vec = compare_funct(rec1, rec2)  # Compare them

          # Put result into weight vector dictionary
          #
          if (self.weight_vec_file == None):
            weight_vec_dict[(rec_ident1, rec_ident2)] = w_vec
          else:
            weight_vec_writer.writerow([rec_ident1, rec_ident2]+w_vec)

          comp_done += 1

          if ((comp_done % progress_report_cnt) == 0):
            self.__log_comparison_progress__(comp_done, start_time)

        rec_cnt += 1

    else: # A linkage run - - - - - - - - - - - - - - - - - - - - - - - - - - -

      for (rec_ident1, rec1) in self.large_dataset.readall():

        rec1_lower = []  # Make all values lowercase
        for rec_val in rec1:
          rec1_lower.append(rec_val.lower())
        rec1 = rec1_lower

        for (rec_ident2, rec2) in small_data_set_dict.iteritems():

          w_vec = compare_funct(rec1, rec2)  # Compare them

          # Put result into weight vector dictionary
          #
          if (self.weight_vec_file == None):
            if (self.ds_swapped == True):
              weight_vec_dict[(rec_ident1, rec_ident2)] = w_vec
            else:
              weight_vec_dict[(rec_ident2, rec_ident1)] = w_vec
          else:
            if (self.ds_swapped == True):
              weight_vec_writer.writerow([rec_ident1, rec_ident2]+w_vec)
            else:
              weight_vec_writer.writerow([rec_ident2, rec_ident1]+w_vec)

          comp_done += 1

          if ((comp_done % progress_report_cnt) == 0):
            self.__log_comparison_progress__(comp_done, start_time)

    used_sec_str = auxiliary.time_string(time.time()-start_time)
    rec_time_str = auxiliary.time_string((time.time()-start_time) / \
                                         self.num_rec_pairs)
    logging.info('Compared %d record pairs in %s (%s per pair)' % \
                 (self.num_rec_pairs, used_sec_str,rec_time_str))

    memory_usage_str = auxiliary.get_memory_usage()
    if (memory_usage_str != None):
      logging.info('  '+memory_usage_str)

    if (self.weight_vec_file == None):
      return [self.__get_field_names_list__(), weight_vec_dict]
    else:
      weight_vec_fp.close()
      return None

# =============================================================================

class BlockingIndex(Indexing):
  """Class that implements the 'classical' blocking used for record linkage.

     Records that have the same index variable values for an index are put into
     the same blocks, and only records within a block are then compared.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """No additional attributes needed, just call base class constructor.

       Note that number of record pairs will not be known after initialisation
       (so it is left at value None).
    """

    Indexing.__init__(self, kwargs)  # Initialise base class

    self.log()  # Log a message

  # ---------------------------------------------------------------------------

  def build(self):
    """Method to build an index data structure.

       Read all records from both files, extract blocking variables and then
       insert records into blocks.
    """

    logging.info('')
    logging.info('Build blocking index: "%s"' % (self.description))

    start_time = time.time()

    self.__records_into_inv_index__()

    num_indices = len(self.index_def)

    # Now calculate number of record pairs - - - - - - - - - - - - - - - - - -
    #
    self.num_rec_pairs = 0

    largest_block_num_rec = -1  # Also keep the largest block size and value
    largest_block_val = ''
    largest_block_index_num = -1  # Index number of the largest block

    for i in range(num_indices):

      if (self.do_deduplication == True):  # A deduplication - - - - - - - - -

        logging.info('  Index %d for data set 1 contains %d blocks' % \
                     (i, len(self.index1[i])))

        for block_val in self.index1[i]: # Loop over all block values in index
          block_num_recs = len(self.index1[i][block_val])

          if (block_num_recs > largest_block_num_rec):  # New largest block
            largest_block_num_rec =   block_num_recs
            largest_block_val =       block_val
            largest_block_index_num = i

          self.num_rec_pairs += block_num_recs*(block_num_recs-1)/2

      else:  # A linkage - - - - - - - - - - - - - - - - - - - - - - - - - - -

        logging.info('  Index %d for data set 1 contains %d blocks, and ' % \
                     (i, len(self.index1[i]))+'for data set 2 contains %d' % \
                     (len(self.index2[i]))+' blocks')

        for block_val in self.index1[i]: # Loop over all block values in index
          block_num_recs1 = len(self.index1[i][block_val])

          if block_val in self.index2[i]:  # Blocking values is both data sets'
                                           # index
            block_num_recs2 = len(self.index2[i][block_val])

            if (block_num_recs1 > largest_block_num_rec):  # New largest block
              largest_block_num_rec =   block_num_recs1
              largest_block_val =       block_val
              largest_block_index_num = i

            if (block_num_recs2 > largest_block_num_rec):  # New largest block
              largest_block_num_rec =   block_num_recs2
              largest_block_val =       block_val
              largest_block_index_num = i

            self.num_rec_pairs += block_num_recs1*block_num_recs2

    logging.info('Built blocking index in %s' % \
                 (auxiliary.time_string(time.time()-start_time)))

    logging.info('  Number of record pairs: %d' % (self.num_rec_pairs))
    logging.info('  Largest block with value "%s" contains %d records' % \
                 (largest_block_val, largest_block_num_rec)+' (in index %d)' \
                 % (largest_block_index_num))

    memory_usage_str = auxiliary.get_memory_usage()
    if (memory_usage_str != None):
      logging.info('  '+memory_usage_str)

    self.status = 'built'  # Update index status

  # ---------------------------------------------------------------------------

  def compact(self):
    """Method to compact an index data structure.

       Make a dictionary of all record pairs over all indices, which removes
       duplicate record pairs.
    """

    NUM_BLOCK_PROGRESS_REPORT = 1000

    logging.info('')
    logging.info('Compact blocking index: "%s"' % (self.description))

    start_time = time.time()

    num_indices = len(self.index_def)

    # Check if index has been built - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.status != 'built'):
      logging.exception('Index "%s" has not been built, compacting is not ' % \
                        (self.description)+'possible')
      raise Exception

    old_num_rec_pairs = self.num_rec_pairs  # Keep old number of record pairs

    rec_pair_dict = {}  # A dictionary with record identifiers from data set 1
                        # as keys and sets of identifiers from data set 2 as
                        # values

    for i in range(num_indices):

      istart_time = time.time()

      num_blocks_done = 0

      if (self.do_deduplication == True):  # A deduplication - - - - - - - - -

        this_index = self.index1[i]  # Shorthand

        for block_val in this_index: # Loop over all block values in index

          block_recs = this_index[block_val]  # All records in this block

          if (len(block_recs) > 1):

            self.__dedup_rec_pairs__(block_recs, rec_pair_dict)

          num_blocks_done += 1

          # Log progress report every XXX blocks processed - - - - - - - - - -
          #
          if ((num_blocks_done % NUM_BLOCK_PROGRESS_REPORT) == 0):
            logging.info('    Processed %d of %d blocks' % \
                         (num_blocks_done, len(this_index)))
            memory_usage_str = auxiliary.get_memory_usage()
            if (memory_usage_str != None):
              logging.info('      '+memory_usage_str)

      else:  # A linkage - - - - - - - - - - - - - - - - - - - - - - - - - - -

        this_index1 = self.index1[i]  # Shorthands
        this_index2 = self.index2[i]

        for block_val in this_index1: # Loop over all block values in index

          block_recs1 = this_index1[block_val]  # All records in this block
                                                # from data set 1
          if block_val in this_index2:  # Blocking values is both data sets'
                                        # index

            block_recs2 = this_index2[block_val] # All records in this block
                                                 # from data set 2

            self.__link_rec_pairs__(block_recs1, block_recs2, rec_pair_dict)

          num_blocks_done += 1

          # Log progress report every XXX blocks processed - - - - - - - - - -
          #
          if ((num_blocks_done % NUM_BLOCK_PROGRESS_REPORT) == 0):
            logging.info('    Processed %d of %d blocks' % \
                         (num_blocks_done, len(this_index1)))
            memory_usage_str = auxiliary.get_memory_usage()
            if (memory_usage_str != None):
              logging.info('      '+memory_usage_str)

      logging.info('  Compacted blocking index %d in %s' % \
                   (i, auxiliary.time_string(time.time()-istart_time)))

      self.index1[i].clear()  # Not needed anymore
      self.index2[i].clear()

      logging.info('    Explicitly run garbage collection')
      gc.collect()

      memory_usage_str = auxiliary.get_memory_usage()
      if (memory_usage_str != None):
        logging.info('      '+memory_usage_str)

    self.rec_pair_dict = rec_pair_dict

    self.num_rec_pairs = 0  # Count lengths of all record identifier sets - - -

    for rec_ident2_set in self.rec_pair_dict.itervalues():
      self.num_rec_pairs += len(rec_ident2_set)

    logging.info('Compacted blocking index in %s' % \
                 (auxiliary.time_string(time.time()-start_time)))
    logging.info('  Old number of record pairs: %d' % (old_num_rec_pairs))
    logging.info('  New number of record pairs: %d' % (self.num_rec_pairs))

    self.status = 'compacted'  # Update index status

  # ---------------------------------------------------------------------------

  def run(self, length_filter_perc = None, cut_off_threshold = None):
    """Iterate over all blocks in the index.

       Compare the record pairs as produced by the blocking process, and return
       a weight vector dictionary with keys made of a tuple (record identifier
       1, record identifier 2), and corresponding values the comparison
       weights.
    """

    logging.info('')
    logging.info('Started comparison of %d record pairs' % \
                 (self.num_rec_pairs))
    if (self.log_funct != None):
      self.log_funct('Started comparison of %d record pairs' % \
                     (self.num_rec_pairs))

    # Check if index has been compacted - - - - - - - - - - - - - - - - - - - -
    #
    if (self.status != 'compacted'):
      logging.exception('Index "%s" has not been compacted, running ' % \
                        (self.description)+'comparisons not possible')
      raise Exception

    # Compare the records
    #
    return self.__compare_rec_pairs_from_dict__(length_filter_perc,
                                                cut_off_threshold)

# =============================================================================

class SortingIndex(Indexing):
  """Class that implements the 'sorted neighbourhood' indexing approach based
     on an inverted index.

     The index variable values are sorted alphabetically and then a window
     (with user specified size) is moved over these sorted values. All records
     in the blocks covered by the current window position are then compared to
     each other.

     The additional argument (besides the base class arguments) which has to be
     set when this index is initialised is:

       window_size  A positive integer that gives the size of the moving window
                    in number of index variable values.

     Note that a window_size of 1 will result in the same records being
     compared as with the standard blocking approach (as a window of size 1
     does not cover any neighbouring index variable values).
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'window_size' argument first, then call the
       base class constructor.
    """

    self.window_size = None  # Set the window size to not defined

    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('window')):
        auxiliary.check_is_integer('window_size', value)
        auxiliary.check_is_positive('window_size', value)
        self.window_size = value

      else:
        base_kwargs[keyword] = value

    Indexing.__init__(self, base_kwargs)  # Initialise base class

    # Make sure 'window_size' attribute is set - - - - - - - - - - - - - - - -
    #
    auxiliary.check_is_integer('window_size', self.window_size)
    auxiliary.check_is_positive('window_size', self.window_size)

    self.log([('Window size', self.window_size)])  # Log a message

  # ---------------------------------------------------------------------------

  def build(self):
    """Method to build an index data structure.

       Read all records from both files, extract blocking variables and then
       insert records into blocks.
    """

    logging.info('')
    logging.info('Build sorting index: "%s"' % (self.description))

    start_time = time.time()

    self.__records_into_inv_index__()  # Read records and put into index

    logging.info('Built sorting index in %s' % \
                 (auxiliary.time_string(time.time()-start_time)))

    memory_usage_str = auxiliary.get_memory_usage()
    if (memory_usage_str != None):
      logging.info('  '+memory_usage_str)

    self.status = 'built'  # Update index status

  # ---------------------------------------------------------------------------

  def compact(self):
    """Method to compact an index data structure.

       Make a dictionary of all record pairs over all indices, which removes
       duplicate record pairs.
    """

    NUM_BLOCK_PROGRESS_REPORT = 1000

    logging.info('')
    logging.info('Compact sorting index: "%s"' % (self.description))

    start_time = time.time()

    # Check if index has been built - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.status != 'built'):
      logging.exception('Index "%s" has not been built, compacting is not ' % \
                        (self.description)+'possible')
      raise Exception

    num_indices = len(self.index_def)

    dedup_rec_pair_funct = self.__dedup_rec_pairs__  # Shorthands
    link_rec_pair_funct =  self.__link_rec_pairs__
    w =                    self.window_size

    rec_pair_dict = {}  # A dictionary with record identifiers from data set 1
                        # as keys and sets of identifiers from data set 2 as
                        # values

    for i in range(num_indices):

      istart_time = time.time()

      num_blocks_done = 0

      if (self.do_deduplication == True):  # A deduplication - - - - - - - - -

        this_index = self.index1[i]  # Shorthand

        w_block_len = [0]*w  # Counts of number of records in blocks in window

        curr_window_recs = []  # List of record identifiers in current window

        block_val_list = this_index.keys()  # Get all blocking values
        block_val_list.sort()

        num_block_vals = len(block_val_list)

        for j in xrange(num_block_vals):  # Loop over all blocks

          w_j = j % w  # Modulo window size

          # Remove record identifiers from previous block
          #
          curr_window_recs = curr_window_recs[w_block_len[w_j]:]

          new_window_recs = set()

          # Insert record identifiers from this block that are not in yet
          #
          for rec_ident in this_index[block_val_list[j]]:

            if (rec_ident not in curr_window_recs):
              new_window_recs.add(rec_ident)

          # Number of new records (previously not in window) from this block
          #
          w_block_len[w_j] = len(new_window_recs)

          curr_window_recs += list(new_window_recs)

          if (len(curr_window_recs) > 1):

            dedup_rec_pair_funct(curr_window_recs, rec_pair_dict)

          num_blocks_done += 1

          # Log progress report every XXX blocks processed - - - - - - - - - -
          #
          if ((num_blocks_done % NUM_BLOCK_PROGRESS_REPORT) == 0):
            logging.info('    Processed %d of %d blocks' % \
                         (num_blocks_done, num_block_vals))
            memory_usage_str = auxiliary.get_memory_usage()
            if (memory_usage_str != None):
              logging.info('      '+memory_usage_str)

        del curr_window_recs

      else:  # A linkage - - - - - - - - - - - - - - - - - - - - - - - - - - -

        this_index1 = self.index1[i]  # Shorthands
        this_index2 = self.index2[i]

        w_block_len1 = [0]*w  # Counts of number of records in blocks in window
        w_block_len2 = [0]*w

        curr_window_recs1 = []  # List of record identifiers in current window
        curr_window_recs2 = []

        block_val_list1 = this_index1.keys()  # Get all blocking values
        block_val_list1.sort()

        block_val_list2 = this_index2.keys()
        block_val_list2.sort()

        num_block_vals1 = len(block_val_list1)
        num_block_vals2 = len(block_val_list2)

        # Create a sorted combined list of all blocking values
        #
        comb_block_values = list(set(block_val_list1+block_val_list2))
        comb_block_values.sort()

        j = 0  # Iteration counters in block values 1 and 2
        k = 0

        do_delete1 = False  # Flags, signalling if current window block
        do_delete2 = False

        # Loop over combined list of block values - - - - - - - - - - - - - - -
        #
        for block_val in comb_block_values:

          if (block_val in block_val_list1):  # Advance window for data set 1

            w_j = j % w  # Modulo window size

            # Remove record identifiers from previous block
            #
            curr_window_recs1 = curr_window_recs1[w_block_len1[w_j]:]

            new_window_recs = set()

            # Insert record identifiers from this block that are not in yet
            #
            for rec_ident in this_index1[block_val]:

              if (rec_ident not in curr_window_recs1):
                new_window_recs.add(rec_ident)

            # Number of new records (previously not in window) from this block
            #
            w_block_len1[w_j] = len(new_window_recs)

            curr_window_recs1 += list(new_window_recs)

            j += 1  # Advance pointer into data set 1 blocks

          if (block_val in block_val_list2):  # Advance window for data set 2 -

            w_k = k % w  # Modulo window size

            # Remove record identifiers from previous block
            #
            curr_window_recs2 = curr_window_recs2[w_block_len2[w_k]:]

            new_window_recs = set()

            # Insert record identifiers from this block that are not in yet
            #
            for rec_ident in this_index2[block_val]:

              if (rec_ident not in curr_window_recs2):
                new_window_recs.add(rec_ident)

            # Number of new records (previously not in window) from this block
            #
            w_block_len2[w_k] = len(new_window_recs)

            curr_window_recs2 += list(new_window_recs)

            k += 1  # Advance pointer into data set 2 blocks

          link_rec_pair_funct(curr_window_recs1, curr_window_recs2,
                              rec_pair_dict)

          num_blocks_done += 1

          # Log progress report every XXX blocks processed - - - - - - - - - -
          #
          if ((num_blocks_done % NUM_BLOCK_PROGRESS_REPORT) == 0):
            logging.info('    Processed %d of %d blocks' % \
                         (num_blocks_done, len(comb_block_values)))
            memory_usage_str = auxiliary.get_memory_usage()
            if (memory_usage_str != None):
              logging.info('      '+memory_usage_str)

        del curr_window_recs1
        del curr_window_recs2

      logging.info('  Compacted sorting index %d in %s' % \
                   (i, auxiliary.time_string(time.time()-istart_time)))

      self.index1[i].clear()  # Not needed anymore
      self.index2[i].clear()

      logging.info('    Explicitly run garbage collection')
      gc.collect()

      memory_usage_str = auxiliary.get_memory_usage()
      if (memory_usage_str != None):
        logging.info('      '+memory_usage_str)

    self.rec_pair_dict = rec_pair_dict

    self.num_rec_pairs = 0  # Count lengths of all record identifier sets - - -

    for rec_ident2_set in self.rec_pair_dict.itervalues():
      self.num_rec_pairs += len(rec_ident2_set)

    logging.info('Compacted sorting index in %s' % \
                 (auxiliary.time_string(time.time()-start_time)))
    logging.info('  Number of record pairs: %d' % (self.num_rec_pairs))

    self.status = 'compacted'  # Update index status

  # ---------------------------------------------------------------------------

  def run(self, length_filter_perc = None, cut_off_threshold = None):
    """Iterate over all blocks in the index.

       Compare the record pairs as produced by the sorting indexing process,
       and return a weight vector dictionary with keys made of a tuple (record
       identifier 1, record identifier 2), and corresponding values the
       comparison weights.
    """

    logging.info('')
    logging.info('Started comparison of %d record pairs' % \
                 (self.num_rec_pairs))
    if (self.log_funct != None):
      self.log_funct('Started comparison of %d record pairs' % \
                     (self.num_rec_pairs))

    # Check if index has been compacted - - - - - - - - - - - - - - - - - - - -
    #
    if (self.status != 'compacted'):
      logging.exception('Index "%s" has not been compacted, running ' % \
                        (self.description)+'comparisons not possible')
      raise Exception

    # Compare the records
    #
    return self.__compare_rec_pairs_from_dict__(length_filter_perc,
                                                cut_off_threshold)


# =============================================================================

class SortingArrayIndex(Indexing):
  """Class that implements the 'sorted neighbourhood' indexing approach based
     on a sorted array.

     The index variable values are sorted alphabetically and then a window
     (with user specified size) is moved over these sorted values. All records
     in the blocks covered by the current window position are then compared to
     each other.

     The additional argument (besides the base class arguments) which has to be
     set when this index is initialised is:

       window_size  A positive integer that gives the size of the moving window
                    in number of index variable values.

     Note that a window_size of 1 will result in no records being compared with
     each other (this is different from the previous SortingIndex above).
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'window_size' argument first, then call the
       base class constructor.
    """

    self.window_size = None  # Set the window size to not defined

    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('window')):
        auxiliary.check_is_integer('window_size', value)
        auxiliary.check_is_positive('window_size', value)
        if (value == 1):
          logging.exception('Window size must be larger than 1.')
          raise Exception

        self.window_size = value

      else:
        base_kwargs[keyword] = value

    Indexing.__init__(self, base_kwargs)  # Initialise base class

    # Make sure 'window_size' attribute is set - - - - - - - - - - - - - - - -
    #
    auxiliary.check_is_integer('window_size', self.window_size)
    auxiliary.check_is_positive('window_size', self.window_size)

    self.log([('Window size', self.window_size)])  # Log a message

  # ---------------------------------------------------------------------------

  def build(self):
    """Method to build an index data structure.

       Read all records from both files, extract blocking variables and then
       insert records into blocks.
    """

    logging.info('')
    logging.info('Build sorted array index: "%s"' % (self.description))

    start_time = time.time()

    self.__records_into_inv_index__()  # Read records and put into index

    logging.info('Built sorted array index in %s' % \
                 (auxiliary.time_string(time.time()-start_time)))

    memory_usage_str = auxiliary.get_memory_usage()
    if (memory_usage_str != None):
      logging.info('  '+memory_usage_str)

    self.status = 'built'  # Update index status

  # ---------------------------------------------------------------------------

  def compact(self):
    """Method to compact an index data structure.

       Make a dictionary of all record pairs over all indices, which removes
       duplicate record pairs.
    """

    NUM_BLOCK_PROGRESS_REPORT = 1000

    logging.info('')
    logging.info('Compact sorted array index: "%s"' % (self.description))

    start_time = time.time()

    # Check if index has been built - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.status != 'built'):
      logging.exception('Index "%s" has not been built, compacting is not ' % \
                        (self.description)+'possible')
      raise Exception

    num_indices = len(self.index_def)

    dedup_rec_pair_funct = self.__dedup_rec_pairs__  # Shorthands
    link_rec_pair_funct =  self.__link_rec_pairs__
    w =                    self.window_size

    rec_pair_dict = {}  # A dictionary with record identifiers from data set 1
                        # as keys and sets of identifiers from data set 2 as
                        # values

    for i in range(num_indices):

      istart_time = time.time()

      num_blocks_done = 0

      if (self.do_deduplication == True):  # A deduplication - - - - - - - - -

        this_index = self.index1[i]  # Shorthand

        rec_sorted_array = []  # The sorted array with the record identifiers
                               # of all records

        # Sorting the unique blocking key values is faster than sorting all
        # values once they are in the sorted array
        #
        block_val_list = this_index.keys()  # Get all blocking values
        block_val_list.sort()

        num_block_vals = len(block_val_list)

        # Loop over all blocks in the inverted index
        #
        for block_key_val in block_val_list:

          rec_id_list = this_index[block_key_val]
          rec_sorted_array += rec_id_list

        # Can be shorter if empty blocking key values occur that
        #
        assert len(rec_sorted_array) <= self.dataset1.num_records

        # Now generate record pairs from the sliding window
        #
        for j in range(0, len(rec_sorted_array)-w+1):

          # Get record identifiers in the current window
          #
          win_rec_id_list = rec_sorted_array[j:j+w][:]
          assert len(win_rec_id_list) == w, \
                 (j,w,self.dataset1.num_records,win_rec_id_list)
          win_rec_id_list.sort()

          ## Possibly improvement: Win pos 1 and all others separate code
          rec_cnt = 1
          for rec_ident1 in win_rec_id_list:
            rec_ident2_set = rec_pair_dict.get(rec_ident1, set())
            for rec_ident2 in win_rec_id_list[rec_cnt:]:
              assert rec_ident1 != rec_ident2
              rec_ident2_set.add(rec_ident2)
            rec_pair_dict[rec_ident1] = rec_ident2_set
            rec_cnt += 1

          num_blocks_done += 1

          # Log progress report every XXX blocks processed - - - - - - - - - -
          #
          if ((num_blocks_done % NUM_BLOCK_PROGRESS_REPORT) == 0):
            logging.info('    Processed %d of %d blocks' % \
                         (num_blocks_done, num_block_vals))
            memory_usage_str = auxiliary.get_memory_usage()
            if (memory_usage_str != None):
              logging.info('      '+memory_usage_str)

        del rec_sorted_array

      else:  # A linkage - - - - - - - - - - - - - - - - - - - - - - - - - - -

        this_index1 = self.index1[i]  # Shorthands
        this_index2 = self.index2[i]

        rec_sorted_array = []  # The sorted array with the record identifiers
                               # of all records

        # Get all unique blocking key values from both indices
        #
        block_val_set = set(this_index1.keys()+this_index2.keys())
        block_val_list = list(block_val_set)
        block_val_list.sort()

        num_block_vals = len(block_val_list)

        # Loop over all blocks in the inverted indices
        #
        for block_key_val in block_val_list:

          rec_id_list1 = this_index1.get(block_key_val, [])
          rec_id_list2 = this_index2.get(block_key_val, [])

          # Merge lists of record identifiers
          #
          if (rec_id_list1 != []) and (rec_id_list2 != []):

            # Split 0-1 into equal intervals and give each record a
            # corresponding floating point number, then sort.
            # For example:
            # rec_id_list1=[a,b,c,d,e]
            #   => [(0.167,a), (0.333,b), (0.5,c), (0.663,d), (0.833,e)]
            # rec_id_list2=[x,y,z]
            #   => [(0.25,x), (0.5,y), (0.75,z)]
            # Merged list: [(0.167,a), (0.25,x), (0.333,b), (0.5,y), (0.5,c),
            #               (0.663,d), (0.75,z), (0.833,e)]

            merge_list = []

            interval1 = 1.0 / (len(rec_id_list1)+1.0)
            j = 1
            for rec_ident in rec_id_list1:
              merge_list.append((j*interval1, rec_ident, '1'))
              assert j*interval1 > 0 and j*interval1 < 1
              j += 1
            interval2 = 1.0 / (len(rec_id_list2)+1.0)
            j = 1
            for rec_ident in rec_id_list2:
              merge_list.append((j*interval2, rec_ident, '2'))
              assert j*interval2 > 0 and j*interval2 < 1
              j += 1

            merge_list.sort()

            assert len(merge_list) == len(rec_id_list1)+len(rec_id_list2)

            for (val, rec_ident, src_index) in merge_list:
              rec_sorted_array.append((rec_ident, src_index))

          elif (rec_id_list1 == []):
            for rec_ident in rec_id_list2:
              rec_sorted_array.append((rec_ident, '2'))  # From index 2

          elif (rec_id_list2 == []):
            for rec_ident in rec_id_list1:
              rec_sorted_array.append((rec_ident, '1'))  # From index 1

        # Can be shorter if empty blocking key values occur
        #
        assert len(rec_sorted_array) <= (self.dataset1.num_records + \
                                         self.dataset2.num_records), \
          (len(rec_sorted_array), (self.dataset1.num_records + \
                                         self.dataset2.num_records))

        # Now generate record pairs from the sliding window
        #
        for j in range(0, len(rec_sorted_array)-w+1):

          # Get record identifiers in the current window
          #
          win_rec_id_list = rec_sorted_array[j:j+w][:]
          assert len(win_rec_id_list) == w, \
                 (j,w,self.dataset1.num_records,win_rec_id_list)

          rec_id_list1 = []  # Record identifiers from index 1
          rec_id_list2 = []  # Record identifiers from index 2

          for (rec_ident, source_index) in win_rec_id_list:
            if (source_index == '1'):
              rec_id_list1.append(rec_ident)
            else:
              rec_id_list2.append(rec_ident)

          for rec_ident1 in rec_id_list1:
            for rec_ident2 in rec_id_list2:

              rec_ident2_set = rec_pair_dict.get(rec_ident1, set())
              rec_ident2_set.add(rec_ident2)
              rec_pair_dict[rec_ident1] = rec_ident2_set

          num_blocks_done += 1

          # Log progress report every XXX blocks processed - - - - - - - - - -
          #
          if ((num_blocks_done % NUM_BLOCK_PROGRESS_REPORT) == 0):
            logging.info('    Processed %d of %d blocks' % \
                         (num_blocks_done, num_block_vals))
            memory_usage_str = auxiliary.get_memory_usage()
            if (memory_usage_str != None):
              logging.info('      '+memory_usage_str)

        del rec_sorted_array

      logging.info('  Compacted sorting index %d in %s' % \
                   (i, auxiliary.time_string(time.time()-istart_time)))

      self.index1[i].clear()  # Not needed anymore
      self.index2[i].clear()

      logging.info('    Explicitly run garbage collection')
      gc.collect()

      memory_usage_str = auxiliary.get_memory_usage()
      if (memory_usage_str != None):
        logging.info('      '+memory_usage_str)

    self.rec_pair_dict = rec_pair_dict

    self.num_rec_pairs = 0  # Count lengths of all record identifier sets - - -

    for rec_ident2_set in self.rec_pair_dict.itervalues():
      self.num_rec_pairs += len(rec_ident2_set)

    logging.info('Compacted sorting index in %s' % \
                 (auxiliary.time_string(time.time()-start_time)))
    logging.info('  Number of record pairs: %d' % (self.num_rec_pairs))

    self.status = 'compacted'  # Update index status

  # ---------------------------------------------------------------------------

  def run(self, length_filter_perc = None, cut_off_threshold = None):
    """Iterate over all blocks in the index.

       Compare the record pairs as produced by the sorting indexing process,
       and return a weight vector dictionary with keys made of a tuple (record
       identifier 1, record identifier 2), and corresponding values the
       comparison weights.
    """

    logging.info('')
    logging.info('Started comparison of %d record pairs' % \
                 (self.num_rec_pairs))
    if (self.log_funct != None):
      self.log_funct('Started comparison of %d record pairs' % \
                     (self.num_rec_pairs))

    # Check if index has been compacted - - - - - - - - - - - - - - - - - - - -
    #
    if (self.status != 'compacted'):
      logging.exception('Index "%s" has not been compacted, running ' % \
                        (self.description)+'comparisons not possible')
      raise Exception

    # Compare the records
    #
    return self.__compare_rec_pairs_from_dict__(length_filter_perc,
                                                cut_off_threshold)

# =============================================================================

class AdaptSortingIndex(Indexing):
  """Class that implements the adaptive sorted neighbourhood indexing approach
     based on an inverted index and an adaptive construction of blocks based on
     the similarities of neighbouring index key values.

     For details see the following paper:

     - Adaptive sorted neighborhood methods for efficient record linkage
       S Yan, D Lee, M.Y Kan and L.C Giles.
       Proceedings of the 7th ACM/IEEE-CS joint conference on Digital libraries
       2007.

     The index variable values are sorted alphabetically first, then a search
     algorithm with adaptive window sizes is moved over the sorted index key
     values, and blocks are generated where two adjacent values have an
     approximate string similarity below a given threshold.

     The initial implementation is a simplified  version of the above
     description that compares all adjacent index key values in a liner fashion.

     The additional argument (besides the base class arguments) which has to be
     set when this index is initialised is:

       str_cmp_funct  A function to compare two strings (as implemented in the
                      stringcmp module).
       str_cmp_thres  The threshold for the string comparison function, must
                      be in (0..1).
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the index specific arguments first, then call the
       base class constructor.
    """

    self.str_cmp_funct = None
    self.str_cmp_thres = None

    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('str_cmp_f')):
        auxiliary.check_is_function_or_method('str_cmp_funct', value)
        self.str_cmp_funct = value

      elif (keyword.startswith('str_cmp_t')):
        auxiliary.check_is_normalised('str_cmp_thres', value)
        self.str_cmp_thres = value

      else:
        base_kwargs[keyword] = value

    Indexing.__init__(self, base_kwargs)  # Initialise base class

    # Make sure parameters have been set - - - - - - - - - - - - - - - - - - -
    #
    auxiliary.check_is_function_or_method('str_cmp_funct', self.str_cmp_funct)
    auxiliary.check_is_normalised('str_cmp_thres', self.str_cmp_thres)

    # A log a message
    #
    self.log([('String comparison function',  self.str_cmp_funct),
              ('String comparison threshold', self.str_cmp_thres)])

  # ---------------------------------------------------------------------------

  def build(self):
    """Method to build an index data structure.

       Read all records from both files, extract blocking variables and then
       insert records into blocks.
    """

    logging.info('')
    logging.info('Build adaptive sorted index: "%s"' % \
                 (self.description))

    start_time = time.time()

    self.__records_into_inv_index__()  # Read records and put into index

    logging.info('Built adaptive sorted index in %s' % \
                 (auxiliary.time_string(time.time()-start_time)))

    memory_usage_str = auxiliary.get_memory_usage()
    if (memory_usage_str != None):
      logging.info('  '+memory_usage_str)

    self.status = 'built'  # Update index status

  # ---------------------------------------------------------------------------

  def compact(self):
    """Method to compact an index data structure.

       Adaptively generate blocks according to the sorted index key values.

       Make a dictionary of all record pairs over all indices, which removes
       duplicate record pairs.
    """

    NUM_BLOCK_PROGRESS_REPORT = 1000

    logging.info('')
    logging.info('Compact adaptive sorted index: "%s"' % \
                 (self.description))

    start_time = time.time()

    # Check if index has been built - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.status != 'built'):
      logging.exception('Index "%s" has not been built, compacting is not ' % \
                        (self.description)+'possible')
      raise Exception

    num_indices = len(self.index_def)

    dedup_rec_pair_funct = self.__dedup_rec_pairs__  # Shorthands
    link_rec_pair_funct =  self.__link_rec_pairs__

    str_cmp_funct =  self.str_cmp_funct
    str_cmp_thres =  self.str_cmp_thres

    rec_pair_dict = {}  # A dictionary with record identifiers from data set 1
                        # as keys and sets of identifiers from data set 2 as
                        # values

    for i in range(num_indices):

      istart_time = time.time()

      num_blocks_done = 0

      if (self.do_deduplication == True):  # A deduplication - - - - - - - - -

        this_index = self.index1[i]  # Shorthand

        index_val_list = this_index.keys()  # Sort all unique index key values
        index_val_list.sort()

        num_index_vals = len(index_val_list)

#        win_start_index = 0  # Position in the index key value list of the
#                             # start of the current window
#        win_end_index =   0  # Make sure it works even if there is one index
#                             # key value only
#
#        while (win_end_index < num_index_vals):
#
#          first_val = index_val_list[win_start_index]
#          last_val =  index_val_list[win_end_index]
#
#          # Increase window size as long as index values are similar
#          #
#          while ((str_cmp_funct(first_val, last_val) > str_cmp_thres) and \
#                 (win_end_index < num_index_vals-1)):
#            win_end_index += 1
#            last_val =  index_val_list[win_end_index]
#
#          # Generate the list of record identifiers from this window
#          #
#          curr_win_record_set = set()
#          for this_val in index_val_list[win_start_index:win_end_index+1]:
#
#            curr_rec_id_list = this_index[this_val]
#
#            curr_win_record_set = \
#                               curr_win_record_set.union(set(curr_rec_id_list))
#
#          if (len(curr_win_record_set) > 1):
#
#            dedup_rec_pair_funct(list(curr_win_record_set), rec_pair_dict)
#
#          win_end_index += 1
#          win_start_index = win_end_index

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # The code below follows algorithm 2 from the above mentioned paper. It
        # is more efficient than the above simple linear scan.
        #
        # The notation given in this algorithm will be used as much as possible.
        #
        block_start_index = 0

        while (block_start_index < num_index_vals):

          w_first = block_start_index
          w = 1  # To make sure this works even with one single index key value
          w_last = w_first + w

          # Get the first and last index key values in the current window
          #
          first_val = index_val_list[w_first]
          last_val =  index_val_list[w_last-1]

          # Enlargement phase: Move the window forward as long as index key
          # values are similar
          #
          while ((str_cmp_funct(first_val, last_val) > str_cmp_thres) and
                 (w_last < num_index_vals)):
            w_first = w_last-1  # Make sure the windows overlap
            w *= 2              # Geometric increase in the window size
            w_last += w-1       # Adjust for overlap
            if (w_last > num_index_vals):
              w_last = num_index_vals  # Reached end of array

            first_val = last_val
            last_val =  index_val_list[w_last-1]

          # Retrenchment phase: Find the boundary pair (use simple linear scan)
          # (the retrenchment phase as described in the above mentioned paper is
          # yet to be implementd)
          #
          tmp_pos = block_start_index

          # Take care of special case where last block is 1 index key value only
          #
          if (tmp_pos+1 == num_index_vals):
            block_end_index = tmp_pos+1

          else:
            first_val =  index_val_list[tmp_pos]
            second_val = index_val_list[tmp_pos+1]

            while ((str_cmp_funct(first_val, second_val) > str_cmp_thres) and
                   ((tmp_pos+1) < num_index_vals)):
              tmp_pos += 1
              first_val =  second_val
              second_val = index_val_list[tmp_pos]

            block_end_index = tmp_pos+1

          # Generate the list of record identifiers from this block
          #
          curr_win_record_set = set()
          for this_val in index_val_list[block_start_index:block_end_index]:
            curr_rec_id_list = this_index[this_val]

            curr_win_record_set = \
                                curr_win_record_set.union(set(curr_rec_id_list))

          if (len(curr_win_record_set) > 1):
            dedup_rec_pair_funct(list(curr_win_record_set), rec_pair_dict)

          block_start_index = block_end_index

          num_blocks_done += 1

          del curr_win_record_set

          # Log progress report every XXX blocks processed - - - - - - - - - -
          #
          if ((num_blocks_done % NUM_BLOCK_PROGRESS_REPORT) == 0):
            logging.info('    Processed %d of %d blocks' % \
                         (num_blocks_done, num_index_vals))
            memory_usage_str = auxiliary.get_memory_usage()
            if (memory_usage_str != None):
              logging.info('      '+memory_usage_str)

      else:  # A linkage - - - - - - - - - - - - - - - - - - - - - - - - - - -

        this_index1 = self.index1[i]  # Shorthands
        this_index2 = self.index2[i]

        # Get and sort all unique index key values from both indices
        #
        index_val_set = set(this_index1.keys()+this_index2.keys())
        index_val_list = list(index_val_set)
        index_val_list.sort()

        num_index_vals = len(index_val_list)

#        win_start_index = 0  # Position in the index key value list of the
#                             # start of the current window
#        win_end_index =   0  # Make sure it works even if there is one index
#                             # key value only
#
#        while (win_end_index < num_index_vals):
#
#          first_val = index_val_list[win_start_index]
#          last_val =  index_val_list[win_end_index]
#
#          # Increase window size as long as index values are similar
#          #
#          while ((str_cmp_funct(first_val, last_val) > str_cmp_thres) and \
#                 (win_end_index < num_index_vals-1)):
#            win_end_index += 1
#            last_val =  index_val_list[win_end_index]
#
#          # Generate the lists of record identifiers from this window
#          #
#          curr_win_record_set1 = set()
#          curr_win_record_set2 = set()
#
#          for this_val in index_val_list[win_start_index:win_end_index+1]:
#
#            if (this_val in this_index1):
#              curr_rec_id_list1 = this_index1[this_val]
#
#              curr_win_record_set1 = \
#                             curr_win_record_set1.union(set(curr_rec_id_list1))
#
#            if (this_val in this_index2):
#              curr_rec_id_list2 = this_index2[this_val]
#
#              curr_win_record_set2 = \
#                             curr_win_record_set2.union(set(curr_rec_id_list2))
#
#          if ((len(curr_win_record_set1) > 0) and \
#              (len(curr_win_record_set2) > 0)):
#
#            link_rec_pair_funct(list(curr_win_record_set1),
#                                list(curr_win_record_set2), rec_pair_dict)
#
#          win_end_index += 1
#          win_start_index = win_end_index

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # The code below follows algorithm 2 from the above mentioned paper. It
        # is more efficient than the above simple linear scan.
        #
        # The notation given in this algorithm will be used as much as possible.
        #
        block_start_index = 0

        while (block_start_index < num_index_vals):

          w_first = block_start_index
          w = 1  # To make sure this works even with one single index key value
          w_last = w_first + w

          # Get the first and last index key values in the current window
          #
          first_val = index_val_list[w_first]
          last_val =  index_val_list[w_last-1]

          # Enlargement phase: Move the window forward as long as index key
          # values are similar
          #
          while ((str_cmp_funct(first_val, last_val) > str_cmp_thres) and
                 (w_last < num_index_vals)):
            w_first = w_last-1  # Make sure the windows overlap
            w *= 2              # Geometric increase in the window size
            w_last += w-1       # Adjust for overlap
            if (w_last > num_index_vals):
              w_last = num_index_vals  # Reached end of array

            first_val = last_val
            last_val =  index_val_list[w_last-1]

          # Retrenchment phase: Find the boundary pair (use simple linear scan)
          # (the retrenchment phase as described in the above mentioned paper is
          # yet to be implementd)
          #
          #
          tmp_pos = block_start_index

          # Take care of special case where last block is 1 index key value only
          #
          if (tmp_pos+1 == num_index_vals):
            block_end_index = tmp_pos+1

          else:
            first_val =  index_val_list[tmp_pos]
            second_val = index_val_list[tmp_pos+1]

            while ((str_cmp_funct(first_val, second_val) > str_cmp_thres) and
                   ((tmp_pos+1) < num_index_vals)):
              tmp_pos += 1
              first_val =  second_val
              second_val = index_val_list[tmp_pos]

            block_end_index = tmp_pos+1

          # Generate the list of record identifiers from this block
          #
          curr_win_record_set1 = set()
          curr_win_record_set2 = set()

          for this_val in index_val_list[block_start_index:block_end_index]:

            if (this_val in this_index1):
              curr_rec_id_list1 = this_index1[this_val]

              curr_win_record_set1 = \
                             curr_win_record_set1.union(set(curr_rec_id_list1))

            if (this_val in this_index2):
              curr_rec_id_list2 = this_index2[this_val]

              curr_win_record_set2 = \
                             curr_win_record_set2.union(set(curr_rec_id_list2))

          if ((len(curr_win_record_set1) > 0) and \
              (len(curr_win_record_set2) > 0)):
            link_rec_pair_funct(list(curr_win_record_set1),
                                list(curr_win_record_set2), rec_pair_dict)

          block_start_index = block_end_index

          num_blocks_done += 1

          del curr_win_record_set1, curr_win_record_set2

          # Log progress report every XXX blocks processed - - - - - - - - - -
          #
          if ((num_blocks_done % NUM_BLOCK_PROGRESS_REPORT) == 0):
            logging.info('    Processed %d of %d blocks' % \
                         (num_blocks_done, num_index_vals))
            memory_usage_str = auxiliary.get_memory_usage()
            if (memory_usage_str != None):
              logging.info('      '+memory_usage_str)

      logging.info('  Compacted sorting index %d in %s' % \
                   (i, auxiliary.time_string(time.time()-istart_time)))

      self.index1[i].clear()  # Not needed anymore
      self.index2[i].clear()

      logging.info('    Explicitly run garbage collection')
      gc.collect()

      memory_usage_str = auxiliary.get_memory_usage()
      if (memory_usage_str != None):
        logging.info('      '+memory_usage_str)

    self.rec_pair_dict = rec_pair_dict

    self.num_rec_pairs = 0  # Count lengths of all record identifier sets - - -

    for rec_ident2_set in self.rec_pair_dict.itervalues():
      self.num_rec_pairs += len(rec_ident2_set)

    logging.info('Compacted sorting index in %s' % \
                 (auxiliary.time_string(time.time()-start_time)))
    logging.info('  Number of record pairs: %d' % (self.num_rec_pairs))

    self.status = 'compacted'  # Update index status

  # ---------------------------------------------------------------------------

  def run(self, length_filter_perc = None, cut_off_threshold = None):
    """Iterate over all blocks in the index.

       Compare the record pairs as produced by the sorting indexing process,
       and return a weight vector dictionary with keys made of a tuple (record
       identifier 1, record identifier 2), and corresponding values the
       comparison weights.
    """

    logging.info('')
    logging.info('Started comparison of %d record pairs' % \
                 (self.num_rec_pairs))
    if (self.log_funct != None):
      self.log_funct('Started comparison of %d record pairs' % \
                     (self.num_rec_pairs))

    # Check if index has been compacted - - - - - - - - - - - - - - - - - - - -
    #
    if (self.status != 'compacted'):
      logging.exception('Index "%s" has not been compacted, running ' % \
                        (self.description)+'comparisons not possible')
      raise Exception

    # Compare the records
    #
    return self.__compare_rec_pairs_from_dict__(length_filter_perc,
                                                cut_off_threshold)






# =============================================================================

class QGramIndex(Indexing):
  """Class that implements an indexing structure based on q-grams and allows
     for fuzzy 'blocking'.

     This index implements a data structure based on q-grams (sub-strings of
     length q, e.g. bigrams (q=2): 'peter' -> 'pe','et','te','er') and allows
     for fuzzy blocking.

     The basic idea is that the index variable values will be converted into a
     list of q-grams, and sub-lists will be built using a user provided
     threshold (a number between 0.0 and 1.0) of all possible combinations.
     The resulting q-gram sub-lists will be inserted into an inverted index,
     i.e. record identifiers in the blocks will be inserted into Python
     dictionaries for each q-gram sub-list. This inverted index will then be
     used to retrieve the blocks.

     The additional argument (besides the base class arguments) which has to be
     set when this index is initialised is:

       q               The length of the q-grams to be used (must be at least
                       1). The default value is 2 (i.e. bigrams)
       padded          If set to True (default), the beginning and end of the
                       strings will be padded with (q-1) special characters, if
                       False no padding will be done.
       threshold       A number between 0.0 (not included) and 1.0, according
                       to which the q-gram sub-lists will be calculated.

     For example, assume an indexing definition contains:

       index_def_1 = [['surname','sname',None,[]], ... ]

     and with q=2 and the threshold set to 0.8. If a record value 'peter' is
     given in the 'surname' field, the corresponding q-gram list will be
     ['pe','et','te','er'] with four elements. With a threshold of 0.8, 4*0.8 =
     3.2, rounded to 3, all sub-list combinations of length 3 are calculated.
     For the given example they are:

       ['pe','et','te']
       ['pe','et','er']
       ['pe','te','er']
       ['et','te','er']

     Converted back into strings, the all records with a given name value
     'peter' will be inserted into the inverted index blocks with keys
     'peette', 'peeter', 'peteer', and 'etteer'.

     The lower the threshold, the shorter the sub-lists, but also the more
     sub-lists there will be per field value, resulting in more (smaller
     blocks) in the inverted index.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'q', padded' and 'threshold' arguments first,
       then call the base class constructor.
    """

    self.padded =    True
    self.q =         2
    self.threshold = None

    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('padd')):
        auxiliary.check_is_flag('padded', value)
        self.padded = value

      elif (keyword == 'q'):
        auxiliary.check_is_integer('q', value)
        auxiliary.check_is_positive('q', value)
        self.q = value

      elif (keyword.startswith('thres')):
        auxiliary.check_is_normalised('threshold', value)
        self.threshold = value

      else:
        base_kwargs[keyword] = value

    Indexing.__init__(self, base_kwargs)  # Initialise base class

    # Make sure 'threshold' attribute is set - - - - - - - - - - - - - - - - -
    #
    auxiliary.check_is_normalised('threshold', self.threshold)

    self.log([('Threshold', self.threshold),
              ('q', self.q),
              ('Padded flag', self.padded)])  # Log a message

    self.QGRAM_START_CHAR = chr(1)
    self.QGRAM_END_CHAR =   chr(2)

  # ---------------------------------------------------------------------------

  def build(self):
    """Method to build an index data structure.

       Read all records from both files, extract blocking variables and then
       insert records into blocks.
    """

    logging.info('')
    logging.info('Build %d-gram index: "%s"' % (self.q, self.description))

    start_time = time.time()

    # First build the basic (blocking) inverted index
    #
    self.__records_into_inv_index__()  # Read records and put into index

    num_indices = len(self.index_def)

    # Next create an index with keys being q-gram sub-list values - - - - - - -
    #
    self.qgram_index1 = {}
    self.qgram_index2 = {}

    for i in range(num_indices):  # Similar to basic index
      self.qgram_index1[i] = {}  # Index for data set 1
      self.qgram_index2[i] = {}  # Index for data set 2

    q = self.q  # Shorthands
    padded = self.padded
    threshold = self.threshold

    if (padded == True):  # More shorthands
      start_str = (q-1)*self.QGRAM_START_CHAR
      end_str =   (q-1)*self.QGRAM_END_CHAR

    # Depending upon threshold values use one of two sub-list methods
    # (the switch value of 0.75 was found experimentally)
    #
    if (self.threshold > 0.75):
      qgram_sublist_funct = self.__get_sublists1__
    else:
      qgram_sublist_funct = self.__get_sublists2__

    logging.info('Convert basic inverted index into %d-gram index' % (q))

    num_blocks =       0
    num_qgram_blocks = 0

    for i in range(num_indices):

      # Build a list of data structures needed for the build process
      #
      index_list = [(self.index1[i], self.qgram_index1[i],0)]  # For data set 1

      if (self.do_deduplication == False):  # If linkage append data set 2
        index_list.append((self.index2[i], self.qgram_index2[i],1))

      for (basic_index, qgram_index, ds_index) in index_list:

        qstart_time = time.time()

        num_blocks += len(basic_index)

        for index_val in basic_index:  # Loop over all the values in this index

          if (padded == True):
            qgram_str = '%s%s%s' % (start_str, index_val, end_str)
          else:
            qgram_str = index_val

          # Create the q-gram list
          #
          qgram_list= [qgram_str[j:j+q] for j in xrange(len(qgram_str)-(q-1))]

          # Calculate length of sub-lists needed and then create them
          #
          num_qgrams = len(qgram_list)

          min_num_qgrams = max(1, int(num_qgrams*threshold))

          qgram_sublists = qgram_sublist_funct(qgram_list, min_num_qgrams)

          # Convert q-gram sub-lists into strings and insert into q-gram index
          #
          for qgram_sublist in qgram_sublists:

            qgram_substr = ''.join(qgram_sublist)

            qgram_index_set = qgram_index.get(qgram_substr, set())
            qgram_index_set.add(index_val)
            qgram_index[qgram_substr] = qgram_index_set

        num_qgram_blocks += len(qgram_index)

        logging.info('  Built %d-gram index %d for data set %d in %s' % \
               (q,i,ds_index+1,auxiliary.time_string(time.time()-qstart_time)))

    logging.info('Built %d-gram index in %s' % \
                 (q, auxiliary.time_string(time.time()-start_time)))

    logging.info('  Number of basic index blocks (number of different ' + \
                 'index variable values): %d' % (num_blocks))
    logging.info('  Number of %d-gram index blocks (number of different ' % \
                 (q) + '%d-gram index values):  %d' % (q, num_qgram_blocks))

    memory_usage_str = auxiliary.get_memory_usage()
    if (memory_usage_str != None):
      logging.info('  '+memory_usage_str)

    self.status = 'built'  # Update index status

  # ---------------------------------------------------------------------------

  def compact(self):
    """Method to compact an index data structure.

       Make a dictionary of all record pairs over all indices, which removes
       duplicate record pairs.

       Finally calculate the total number of record pairs.
    """

    NUM_QGRAM_BLOCK_PROGRESS_REPORT = 1000

    logging.info('')
    logging.info('Compact %d-gram index: "%s"' % (self.q, self.description))

    start_time = time.time()

    num_indices = len(self.index_def)

    # Check if index has been built - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.status != 'built'):
      logging.exception('Index "%s" has not been built, compacting is not ' % \
                        (self.description)+'possible')
      raise Exception

    rec_pair_dict = {}  # A dictionary with record identifiers from data set 1
                        # as keys and sets of identifiers from data set 2 as
                        # values

    for i in range(num_indices):

      istart_time = time.time()

      num_qgram_blocks_done = 0

      if (self.do_deduplication == True):  # A deduplication - - - - - - - - -

        this_qgram_index = self.qgram_index1[i]
        this_index =       self.index1[i]

        for qgram_val in this_qgram_index:  # Loop over all q-gram values

          block_recs = []  # Combined list of all record identifiers

          # All record identifiers from the basic index for this q-gram value
          #
          index_val_set = this_qgram_index[qgram_val]

          # Get all the indexing variable values for this gqram value
          #
          for index_val in index_val_set:

            for block_rec in this_index[index_val]:
              if (block_rec not in block_recs):
                block_recs.append(block_rec)

          if (len(block_recs) > 1):

            self.__dedup_rec_pairs__(block_recs, rec_pair_dict)

          num_qgram_blocks_done += 1

          # Log progress report every XXX q-gram blocks processed
          #
          if ((num_qgram_blocks_done % NUM_QGRAM_BLOCK_PROGRESS_REPORT) == 0):
            logging.info('    Processed %d of %d %d-gram blocks' % \
                         (num_qgram_blocks_done, len(this_qgram_index),
                         self.q))
            memory_usage_str = auxiliary.get_memory_usage()
            if (memory_usage_str != None):
              logging.info('      '+memory_usage_str)

      else:  # A linkage - - - - - - - - - - - - - - - - - - - - - - - - - - -

        this_qgram_index1 = self.qgram_index1[i]
        this_qgram_index2 = self.qgram_index2[i]
        this_index1 =       self.index1[i]
        this_index2 =       self.index2[i]

        for qgram_val in this_qgram_index1:  # Loop over all q-gram values

          if (qgram_val in this_qgram_index2):  # A matching q-gram value

            block_recs1 = []  # Combined list of all record identifiers
            block_recs2 = []

            # All record identifiers from the basic indices for this q-gram
            #
            index_val_set1 = this_qgram_index1[qgram_val]

            # Get all the indexing variable values for this gqram value
            #
            for index_val1 in index_val_set1:

              for block_rec in this_index1[index_val1]:
                if (block_rec not in block_recs1):
                  block_recs1.append(block_rec)

            index_val_set2 = this_qgram_index2[qgram_val]

            for index_val2 in index_val_set2:

              for block_rec in this_index2[index_val2]:
                if (block_rec not in block_recs2):
                  block_recs2.append(block_rec)

            self.__link_rec_pairs__(block_recs1, block_recs2, rec_pair_dict)

          num_qgram_blocks_done += 1

          # Log progress report every XXX q-gram blocks processed
          #
          if ((num_qgram_blocks_done % NUM_QGRAM_BLOCK_PROGRESS_REPORT) == 0):
            logging.info('    Processed %d of %d %d-gram blocks' % \
                         (num_qgram_blocks_done, len(self.qgram_index1[i]),
                         self.q))
            memory_usage_str = auxiliary.get_memory_usage()
            if (memory_usage_str != None):
              logging.info('      '+memory_usage_str)

      logging.info('  Compacted %d-gram index %d in %s' % \
                   (self.q, i, auxiliary.time_string(time.time()-istart_time)))

      self.qgram_index1[i].clear()  # Not needed anymore
      self.qgram_index2[i].clear()
      self.index1[i].clear()
      self.index2[i].clear()

      logging.info('    Explicitly run garbage collection')
      gc.collect()

      memory_usage_str = auxiliary.get_memory_usage()
      if (memory_usage_str != None):
        logging.info('    '+memory_usage_str)

    self.rec_pair_dict = rec_pair_dict

    self.num_rec_pairs = 0  # Count lengths of all record identifier sets - - -

    for rec_ident2_set in self.rec_pair_dict.itervalues():
      self.num_rec_pairs += len(rec_ident2_set)

    logging.info('Compacted %d-gram index in %s' % \
                 (self.q, auxiliary.time_string(time.time()-start_time)))
    logging.info('  Number of record pairs: %d' % (self.num_rec_pairs))

    self.status = 'compacted'  # Update index status

  # ---------------------------------------------------------------------------

  def run(self, length_filter_perc = None, cut_off_threshold = None):
    """Iterate over all blocks in the index.

       Compare the record pairs as produced by the q-gram indexing process, and
       return a weight vector dictionary with keys made of a tuple (record
       identifier 1, record identifier 2), and corresponding values the
       comparison weights.
    """

    logging.info('')
    logging.info('Started comparison of %d record pairs' % \
                 (self.num_rec_pairs))
    if (self.log_funct != None):
      self.log_funct('Started comparison of %d record pairs' % \
                     (self.num_rec_pairs))

    # Check if index has been compacted - - - - - - - - - - - - - - - - - - - -
    #
    if (self.status != 'compacted'):
      logging.exception('Index "%s" has not been compacted, running ' % \
                        (self.description)+'comparisons not possible')
      raise Exception

    # Compare the records
    #
    return self.__compare_rec_pairs_from_dict__(length_filter_perc,
                                                cut_off_threshold)

# =============================================================================

class CanopyIndex(Indexing):
  """Class that implements the canopy clustering based indexing.

     For details see the following papers:

     - A comparison of fast blocking methods for record linkage
       Rohan Baxter, Peter Christen and Tim Churches,
       Workshop on Data Cleaning, Record Linkage and Object Consolidation, at
       the 9th ACM SIGKDD, 2003.

     - Learning to Match and cluster high-dimensional data sets
       William W. and Jacob Richman,
       8th ACM SIGKDD, 2002.

     The basic idea is to form overlapping clusters and to assign records into
     one or more of these clusters. These clusters will be formed using q-grams
     and either TF-IDF or simply calculating a proportion of common q-grams
     (an approximation to the Jaccard similarity).

     The additional argument (besides the base class arguments) which has to be
     set when this index is initialised is:

       canopy_method  Determines the method used to form and extract the canopy
                      clusters, as well its parameters. Possible are:
                        ('tfidf','threshold',tight_threshold,loose_threshold)
                        ('tfidf','nearest',remove_nearest,cluster_nearest)
                        ('jaccard','threshold',tight_threshold,loose_threshold)
                        ('jaccard','nearest',remove_nearest,cluster_nearest)
                      with parameters:
                      - tight_threshold   All records with a similarity equal
                                          to or larger than this will be
                                          removed from the pool of records
                                          after inserted into a canopy cluster.
                                          This is a normalised number 0..1 and
                                          must be larger than or equal to the
                                          loose threshold.
                      - loose_threshold   All records with a similarity equal
                                          to or larger than this will be added
                                          into a canopy cluster. This is a
                                          normalised number 0..1 and must be
                                          smaller or equal to the tight
                                          threshold.
                      - remove_nearest    Number of records nearest to the
                                          cluster center that will be removed
                                          from the pool of records after
                                          inserted into a canopy cluster. Must
                                          be a positive integer and smaller
                                          than or equal to the value of cluster
                                          nearest.
                      - cluster_nearest   Number of records nearest to the
                                          cluster center that will be added
                                          into a canopy cluster. Must be a
                                          positive integer and larger than or
                                          equal to the value of remove nearest.
       q                 The length of the q-grams used when creating the
                         TF-IDF data structure. Default is q=2.
       padded            If set to True (default), the beginning and end of the
                         strings will be padded with (q-1) special characters,
                         if False no padding will be done.
       delete_perc       Threshold for deleting common q-grams (if they appear
                         in more than this percentage of all records). Default
                         is None, in which case no q-grams will be deleted.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the arguments 'canopy_method', 'q', 'padded' and
       'delete_perc' first, then call the base class constructor.
    """

    self.canopy_method =  None
    self.q =              2
    self.padded =         True
    self.delete_perc =    None

    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('canopy_m')):
        auxiliary.check_is_tuple('canopy_method', value)
        self.canopy_method = value

      elif (keyword.startswith('padd')):
        auxiliary.check_is_flag('padded', value)
        self.padded = value

      elif (keyword == 'q'):
        auxiliary.check_is_integer('q', value)
        auxiliary.check_is_positive('q', value)
        self.q = value

      elif (keyword.startswith('skip')):
        auxiliary.check_is_flag('skip_missing', value)
        if (value == False):
          logging.exception('Value of argument"skip_missing" cannot be set' + \
                            ' to False for the CanopyIndex')
          raise Exception
        self.skip_missing = value

      elif (keyword.startswith('delete_p')):
        auxiliary.check_is_percentage('delete_perc', value)
        self.delete_perc = value

      else:
        base_kwargs[keyword] = value

    Indexing.__init__(self, base_kwargs)  # Initialise base class

    # Check if canopy method and parameters given are OK - - - - - - - - - - -
    #
    auxiliary.check_is_not_none('canopy_method', self.canopy_method)
    auxiliary.check_is_tuple('canopy_method', self.canopy_method)

    if (len(self.canopy_method) != 4):
      logging.exception('Canopy method tuple needs four elements: %s' % \
                        (str(self.canopy_method)))
      raise Exception

    if (self.canopy_method[0] not in ['tfidf', 'jaccard']):
      logging.exception('Illegal canopy method given (must be "tfidf" or ' + \
                        '"jaccard"): %s' % (str(self.canopy_method)))
      raise Exception

    if (self.canopy_method[1] not in ['threshold', 'nearest']):
      logging.exception('Illegal canopy retrieval method given (must be ' + \
                        '"threshold" or "nearest"): %s' % \
                        (str(self.canopy_method)))
      raise Exception

    if (self.canopy_method[1] == 'threshold'):
      auxiliary.check_is_normalised('tight_threshold', self.canopy_method[2])
      auxiliary.check_is_normalised('loose_threshold', self.canopy_method[3])

      if (self.canopy_method[2] < self.canopy_method[3]):
        logging.exception('Tight threshold is smaller than loose threshold:' \
                          + '%.3f / %.3f' % (self.canopy_method[2],
                                             self.canopy_method[3]))
        raise Exception

    else:  # Is 'nearest'
      auxiliary.check_is_integer('remove_nearest',   self.canopy_method[2])
      auxiliary.check_is_positive('remove_nearest',  self.canopy_method[2])
      auxiliary.check_is_integer('cluster_nearest',  self.canopy_method[3])
      auxiliary.check_is_positive('cluster_nearest', self.canopy_method[3])

      if (self.canopy_method[2] > self.canopy_method[3]):
        logging.exception('Remove nearest is larger than cluster nearest:' \
                          + '%d / %d' % (self.canopy_method[2],
                                         self.canopy_method[3]))
        raise Exception

    self.log([('Canopy method', self.canopy_method),
              ('q', self.q),
              ('Padded flag', self.padded),
              ('Delete percentage', self.delete_perc)])  # Log a message

    self.QGRAM_START_CHAR = chr(1)
    self.QGRAM_END_CHAR =   chr(2)

    if (self.padded == True):
      self.start_str = (self.q-1)*self.QGRAM_START_CHAR
      self.end_str =   (self.q-1)*self.QGRAM_END_CHAR
    else:
      self.start_str = ''
      self.end_str =   ''

  # ---------------------------------------------------------------------------

  def build(self):
    """Method to build an index data structure.

       Load the records from the data sets and put them into an inverted index
       data structure with q-grams as keys.

       If the canopy method uses TF-IDF, also builds a dictionary with the
       index (block) values for each record as well as a dictionary with the
       normalised Euclidean length (the weight W_d) for each index value (see
       book: Managing Gigabytes, page 187)

       For the Jaccard canopy method, the number of q-grams (incl. paddings)
       will be stored for each record.

       Note that for a linkage the records from both data sets will be inserted
       into the same index with a modified record identifier: '0' or '1' is
       added at the beginning of each identifier to distinguish between records
       coming from data set 1 or 2.
    """

    logging.info('')
    logging.info('Build canopy index: "%s"' % (self.description))

    start_time = time.time()

    num_indices = len(self.index_def)

    if (self.canopy_method[0] == 'tfidf'):
      canopy_method_tfidf = True
    else:
      canopy_method_tfidf = False  # For Jaccard methods

    # Initialise the index data structures needed - - - - - - - - - - - - - - -

    # Cache the index (block) value for each record so they can be used later
    # to extract canopies
    #
    self.index_val_cache = {}

    # Inverse document frequency for each q-gram (TF-IDF method only)
    #
    self.qgram_inv_doc_freq_cache = {}

    # Number of q-grams for each record identifier (i.e. its index value)
    # (Jaccard method only)
    #
    self.index_val_num_qgram = {}

    # Initialise dictionaries for each index - - - - - - - - - - - - - - - - -
    #
    for i in range(num_indices):
      self.index1[i] =                   {}  # Main inverted index
      self.index_val_cache[i] =          {}
      self.qgram_inv_doc_freq_cache[i] = {}
      self.index_val_num_qgram[i] =      {}

    # Keep maximum number of q-grams for every index (needed to get TF-IDF)
    #
    max_qgram_count = [-1]*num_indices

    max_qgram_count_val = [('','')]*num_indices  # Keep most frequent q-grams

    qgram_inv_doc_freq_cache = self.qgram_inv_doc_freq_cache  # Shorthands
    index_val_num_qgram =      self.index_val_num_qgram
    index_val_cache =          self.index_val_cache
    index =                    self.index1
    do_dedup =                 self.do_deduplication
    get_index_values_funct =   self.__get_index_values__
    get_qgram_list_funct =     self.__get_qgram_list__
    qgram_list_to_dict_funct = self.__qgram_list_to_dict__

    # Reference to data set 1 and record cache 1
    #
    build_list = [(self.dataset1, self.rec_cache1, self.comp_field_used1, 0)]
    if (do_dedup == False):  # If linkage append data set 2
      build_list.append((self.dataset2, self.rec_cache2,
                          self.comp_field_used2, 1))

    # Step 1: Read data set(s) and build basic inverted index - - - - - - - - -
    #
    for (dataset, rec_cache, comp_field_used, ds_index) in build_list:

      # Calculate a counter for the progress report
      #
      if (self.progress_report != None):
        progress_report_cnt = max(1, int(dataset.num_records / \
                                     (100.0 / self.progress_report)))
      else:  # So no progress report is being logged
        progress_report_cnt = dataset.num_records + 1

      rec_read = 0   # Number of records read from data set

      rstart_time = time.time()  # Start time reading data set

      for (rec_ident, rec) in dataset.readall(): # Read all records in data set

        # Extract record fields needed for comparisons - - - - - - - - - - - -
        #
        comp_rec = []

        field_ind = 0
        for field in rec:
          if (field_ind in comp_field_used):
            comp_rec.append(field.lower())
          else:
            comp_rec.append('')  # Set not needed fields to ''
          field_ind += 1

        rec_cache[rec_ident] = comp_rec  # Put into record cache

        # Now get the index variable values for this record - - - - - - - - - -
        #
        rec_index_val_list = get_index_values_funct(rec, ds_index)

        if (do_dedup == False):
          ds_rec_ident = str(ds_index)+rec_ident  # Add data set identifier
        else:
          ds_rec_ident = rec_ident  # Not needed for deduplication

        for i in range(num_indices):  # Put record identifier into all indices

          index_val = rec_index_val_list[i]

          if (index_val != ''):  # Missing values will be skipped

            # Put index (block) value for this record into index value cache
            #
            index_val_cache[i][ds_rec_ident] = index_val

            qgram_list = get_qgram_list_funct(index_val)

            if (canopy_method_tfidf == False):  # For Jaccard
              index_val_num_qgram[i][ds_rec_ident] = len(set(qgram_list))

            for qgram in qgram_list:

              # Get and update dictionary of record identifiers for this q-gram
              #
              qgram_rec_dict = index[i].get(qgram, {})

              if (ds_rec_ident not in qgram_rec_dict):
                qgram_rec_dict[ds_rec_ident] = 1
              else:
                qgram_count = qgram_rec_dict[ds_rec_ident] + 1
                qgram_rec_dict[ds_rec_ident] = qgram_count

                # Keep maximum count value for each index
                #
                if (qgram_count > max_qgram_count[i]):
                  max_qgram_count[i] =     qgram_count
                  max_qgram_count_val[i] = (qgram, index_val)

              index[i][qgram] = qgram_rec_dict  # Back into inverted index

        rec_read += 1

        if ((rec_read % progress_report_cnt) == 0):
          self.__log_build_progress__(rec_read,dataset.num_records,rstart_time)

      logging.info('  Explicitly run garbage collection')
      gc.collect()

      used_sec_str = auxiliary.time_string(time.time()-rstart_time)
      rec_time_str = auxiliary.time_string((time.time()-rstart_time) / \
                                           dataset.num_records)
      logging.info('Read and indexed %d records in %s (%s per record)' % \
                   (dataset.num_records, used_sec_str, rec_time_str))
      memory_usage_str = auxiliary.get_memory_usage()
      if (memory_usage_str != None):
        logging.info('  '+memory_usage_str)

      logging.info('  Maximum %d-gram counts and example values:', (self.q))
      for i in range(num_indices):
        if (max_qgram_count[i]>-1):  # Some q-grams occured more than once
          logging.info('    Index %d: %d-gram "%s" with count %d in value ' % \
                       (i, self.q, max_qgram_count_val[i][0],
                       max_qgram_count[i]) + '"%s"' % \
                       (max_qgram_count_val[i][1]))
        else:
          logging.info('    Index %d: No %d-gram appeared more than once' % \
                       (i, self.q))

    self.max_qgram_count = max_qgram_count

    # Step 2: Calculate q-gram inverse document frequencies and Euclidean - - -
    # length (weight W_d) for all index values.
    #
    if (canopy_method_tfidf == True):  # Only needed for TF-IDF

      for i in range(num_indices):  # Put record identifier into all indices

        this_index = index[i]

        this_max_qgram_count = float(max_qgram_count[i])

        index_val_Wd_cache = {}

        logging.info('  Calculate %d-gram IDF and Euclidean length (weight' % \
                     (self.q)+' W_d) for index %d' % (i))

        total_num_rec = float(len(index_val_cache[i]))

        for rec_ident in index_val_cache[i]:
          index_val = index_val_cache[i][rec_ident]

          # Check if the W_d weight has not been calculated yet for this value
          #
          if (index_val not in index_val_Wd_cache):

            qgram_list = get_qgram_list_funct(index_val)
            qgram_dict = qgram_list_to_dict_funct(qgram_list)

            W_d = 0.0  # See Managing Gigabytes, page 187

            # Loop over the q-grams and their counts in the given index value
            #
            for (qgram, qgram_count) in qgram_dict.iteritems():

              # Check if the q-gram IDF has not yet been calculated and cached
              #
              if (qgram not in qgram_inv_doc_freq_cache[i]):

                # Number of records with this q-gram
                #
                num_qgram_rec = len(this_index[qgram])

                # Inverse document frequency for this q-gram
                #
                inv_doc_freq = math.log(total_num_rec / num_qgram_rec, 2)

                qgram_inv_doc_freq_cache[i][qgram] = inv_doc_freq

                # Also normalise q-gram counts in the inverted index
                #
                for (tmp_rec_id, tmp_rec_co) in this_index[qgram].iteritems():
                  norm_rec_count = tmp_rec_co / this_max_qgram_count
                  this_index[qgram][tmp_rec_id] = norm_rec_count

              else:  # Get IDF for this g-gram from the cache

                inv_doc_freq = qgram_inv_doc_freq_cache[i][qgram]

              # Calculate TF-IDF weight W_dt for this q-gram in the index value
              # (have floating point value first to avoid integer division)
              #
              #W_dt = qgram_count/this_max_qgram_count * inv_doc_freq
              W_dt = inv_doc_freq * qgram_count / this_max_qgram_count

              # Add to Euclidean length (the weight W_d) of the index value
              #
              W_d += W_dt*W_dt

            # Put square-root of weight back
            #
            index_val_Wd_cache[index_val] = math.sqrt(W_d)

        logging.info('    Number of different index values: %d' % \
                     (len(index_val_Wd_cache)))

        # Divide all q-gram counts by their W_d so it does not have to be done
        # in the compact() routine (where it would have to be done many times)
        #
        for qgram in this_index:

          for tmp_rec_ident in this_index[qgram]:
            tmp_qgram_count = this_index[qgram][tmp_rec_ident]

            tmp_index_val = index_val_cache[i][tmp_rec_ident]
            tmp_W_d =       index_val_Wd_cache[tmp_index_val]

            this_index[qgram][tmp_rec_ident] = tmp_qgram_count / tmp_W_d

        index_val_Wd_cache.clear()  # Entries not needed anymore

        memory_usage_str = auxiliary.get_memory_usage()
        if (memory_usage_str != None):
          logging.info('    '+memory_usage_str)

        del index_val_Wd_cache

    # Log index details and delete q-grams occuring too frequent - - - - - - -

    # Get total number of records in data set(s)
    #
    num_records = self.dataset1.num_records
    if (do_dedup == False):
      num_records += self.dataset2.num_records

    for i in range(num_indices):
      logging.info('  Index %d contains %d different %d-grams' % \
                   (i, len(index[i]), self.q))
      max_len = -1
      max_val = None
      delete_qgram_list = []  # List of common q-grams to be deleted

      for (qgram, qgram_recs) in index[i].iteritems():
        qgram_rec_count = len(qgram_recs)

        if (qgram_rec_count > max_len):
          max_len = qgram_rec_count
          max_val = qgram

        if (self.delete_perc != None):  # Check if q-gram is too common
          if ((100.0*qgram_rec_count/num_records) > self.delete_perc):
            delete_qgram_list.append(qgram)

      logging.info('    Most common %d-gram: "%s" in %d of %d records' % \
                   (self.q, max_val, max_len, num_records))

      for qgram in delete_qgram_list:  # Delete all too frequent q-grams
        del index[i][qgram]

      if (canopy_method_tfidf == False):  # For Jaccard adjust q-gram numbers

        for ds_rec_ident in index_val_cache[i]:

          index_val = index_val_cache[i][ds_rec_ident]
          qgram_set = set(get_qgram_list_funct(index_val))

          changed_len = False
          for del_qgram in delete_qgram_list:

            if (del_qgram in qgram_set):
              qgram_set.remove(del_qgram)
              changed_len = True

          # Update to new length
          #
          if (changed_len == True):
            index_val_num_qgram[i][ds_rec_ident] = len(qgram_set)

      if (delete_qgram_list != []):
        logging.info('    Deleted the following %d-grams (as ' % (self.q) + \
                     'they occured in more than %.2f%% of all records):' % \
                     (self.delete_perc))
        logging.info('        %s' % (str(delete_qgram_list)))

    logging.info('Built canopy index in %s' % \
                 (auxiliary.time_string(time.time()-start_time)))

    logging.info('  Explicitly run garbage collection')
    gc.collect()

    memory_usage_str = auxiliary.get_memory_usage()
    if (memory_usage_str != None):
      logging.info('  '+memory_usage_str)

    self.status = 'built'  # Update index status

  # ---------------------------------------------------------------------------

  def __get_qgram_list__(self, in_str):
    """Convert the given input string into q-grams and return a list of these
       q-grams. Also pad the input string if pedding is activated.
    """

    q =  self.q  # Shorthands
    q1 = q-1

    if (self.padded == True):
      qgram_str = '%s%s%s' % (self.start_str, in_str, self.end_str)
    else:
      qgram_str = in_str

    qgram_list = [qgram_str[j:j+q] for j in xrange(len(qgram_str)-q1)]

    if (qgram_list == []):  # Make sure a string value is returned
      qgram_list = [qgram_str]

    return qgram_list

  # ---------------------------------------------------------------------------

  def __qgram_list_to_dict__(self, qgram_list):
    """Converts the input list of q-grams into a dictionary with the q-grams as
       keys and their frequencies as values.
    """

    qgram_dict = {}

    for qgram in qgram_list:
      qgram_count = qgram_dict.get(qgram, 0) + 1
      qgram_dict[qgram] = qgram_count

    return qgram_dict

  # ---------------------------------------------------------------------------

  def compact(self):
    """Method to compact an index data structure.

       Extract canopies from inverted index and put all resulting record pairs
       into a dictionary.

       Finally calculate the total number of record pairs.
    """

    NUM_CANOPY_PROGRESS_REPORT = 100  # Log a message every XXX canopies

    logging.info('')
    logging.info('Compact canopy index: "%s"' % (self.description))

    start_time = time.time()

    num_indices = len(self.index_def)

    dedup_rec_pairs_funct = self.__dedup_rec_pairs__  # Shorthands
    link_rec_pair_funct =   self.__link_rec_pairs__
    tfidf_canopy_funct =    self.__tfidf_canopy__
    jaccard_canopy_funct =  self.__jaccard_canopy__

    # Check if index has been built - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.status != 'built'):
      logging.exception('Index "%s" has not been built, compacting is not ' % \
                        (self.description)+'possible')
      raise Exception

    rec_pair_dict = {}  # A dictionary with record identifiers from data set 1
                        # as keys and sets of identifiers from data set 2 as
                        # values

    # Select get canopy function according to canopy method
    #
    if (self.canopy_method[0] == 'tfidf'):  # Only needed for TF-IDF
      do_tfidf = True
    else:
      do_tfidf = False

    for i in range(num_indices):

      istart_time = time.time()

      this_index_val_cache = self.index_val_cache[i]  # Shorthand

      num_canopies = 0  # Count the number of canopies created

      smallest_canopy_size =   999999
      smallest_canopy_center = ''    # Indexing value of the smallest canopy
      largest_canopy_size =    -99999
      largest_canopy_center =  ''    # Indexing value of the largest canopy

      # Total number of records in this index
      #
      total_num_rec = float(len(this_index_val_cache))

      logging.info('  Compacting index %d containing %d records and %d ' % \
                   (i, total_num_rec, len(self.index1[i]))+'%d-grams' % \
                   (self.q))

      # Loop over all values, extract canopies and delete records from values
      # cache that are within the tight threshold of a canopy
      #
      while(len(this_index_val_cache) > 0):

        # Get arbitrary record identifier and value from the values cache
        #
        (rec_ident, index_val) = this_index_val_cache.popitem()
        this_index_val_cache[rec_ident] = index_val  # Put back in

        # Get all records in this canopy - - - - - - - - - - - - - - - - - - -
        #
        if (do_tfidf == True):
          canopy_recs = tfidf_canopy_funct(self.index1[i], index_val,
                                          this_index_val_cache,
                                          self.qgram_inv_doc_freq_cache[i],
                                          self.max_qgram_count[i])
        else:
          canopy_recs = jaccard_canopy_funct(self.index1[i], index_val,
                                             this_index_val_cache,
                                             self.index_val_num_qgram[i])

        # Make sure center record is in its canopy
        #
        assert rec_ident in canopy_recs, (rec_ident,index_val,canopy_recs)

        num_canopy_rec = len(canopy_recs)
        num_canopies += 1

        if (num_canopy_rec < smallest_canopy_size):
          smallest_canopy_size =   num_canopy_rec
          smallest_canopy_center = index_val
        elif (num_canopy_rec > largest_canopy_size):
          largest_canopy_size =   num_canopy_rec
          largest_canopy_center = index_val

        # Process record list depending upon deduplication or linkage - - - - -
        #
        if (self.do_deduplication == True):

          if (num_canopy_rec > 1):  # For deduplication at least two records

            # Build record pairs from record identifiers in this canopy
            #
            dedup_rec_pairs_funct(canopy_recs, rec_pair_dict)

        else:  # A linkage - - - - - - - - - - - - - - - - - - - - - - - - - -

          canopy_recs1 = []  # Need to separate records
          canopy_recs2 = []

          for ds_rec_ident in canopy_recs:
            if (ds_rec_ident[0] == '0'):
              canopy_recs1.append(ds_rec_ident[1:])
            else:
              canopy_recs2.append(ds_rec_ident[1:])

            link_rec_pair_funct(canopy_recs1, canopy_recs2, rec_pair_dict)
          del canopy_recs1
          del canopy_recs2

        del canopy_recs

        # Log progress report every XXX canopies - - - - - - - - - - - - - - -
        #
        if ((num_canopies % NUM_CANOPY_PROGRESS_REPORT) == 0):
          logging.info('    Created %d canopies; %d records and ' % \
                       (num_canopies, len(self.index_val_cache[i])) + \
                       '%d %d-grams' % (len(self.index1[i]), self.q)+' left')
          memory_usage_str = auxiliary.get_memory_usage()
          if (memory_usage_str != None):
            logging.info('      '+memory_usage_str)

      # Delete not needed index data to free-up memory - - - - - - - - - - - -
      #
      self.index1[i].clear()  # Not needed anymore
      this_index_val_cache.clear()
      self.qgram_inv_doc_freq_cache[i].clear()

      logging.info('  Compacted canopy index %d in %s' % \
                   (i, auxiliary.time_string(time.time()-istart_time)))
      logging.info('    Produced %d canopies' % (num_canopies))
      logging.info('      Smallest canopy with %d records and center ' % \
                   (smallest_canopy_size)+'index value: "%s"' % \
                   (smallest_canopy_center))
      logging.info('      Largest canopy with %d records and center ' % \
                   (largest_canopy_size)+'index value: "%s"' % \
                   (largest_canopy_center))

      logging.info('  Explicitly run garbage collection')
      gc.collect()

      memory_usage_str = auxiliary.get_memory_usage()
      if (memory_usage_str != None):
        logging.info('    '+memory_usage_str)

    num_rec_pairs = 0  # Count lengths of all record identifier sets - - -

    for rec_ident2_set in rec_pair_dict.itervalues():
      num_rec_pairs += len(rec_ident2_set)

    self.rec_pair_dict = rec_pair_dict  # Save for later used in run()
    self.num_rec_pairs = num_rec_pairs

    logging.info('Compacted canopy index in %s' % \
                 (auxiliary.time_string(time.time()-start_time)))
    logging.info('  Number of record pairs: %d' % (num_rec_pairs))

    memory_usage_str = auxiliary.get_memory_usage()
    if (memory_usage_str != None):
      logging.info('  '+memory_usage_str)

    self.status = 'compacted'  # Update index status

  # ---------------------------------------------------------------------------

  def __tfidf_canopy__(self, index, index_val, index_val_cache,
                       qgram_inv_doc_freq_cache, max_qgram_count):
    """Returns a list of record identifiers in the given inverted index
       according to the canopy method:
       - tfidf,threshold  Returns records that have a TF-IDF similarity larger
                          than the loose threshold to the given index value,
                          and remove records that have a similarity value
                          larger than tight threshold.
       - tfidf,nearest    Returns the 'cluster_nearest' closest records, and
                          removes the 'remove_nearest' nearest records
                          according to TF-IDF similarity.

       Arguments provided:
         index                     The inverted index with q-grams as keys and
                                   dictionaries with record identifiers as
                                   values.
         index_val                 Value of the cluster center.
         index_val_cache           A dictionary with record identifiers as keys
                                   and their index value as values.
         qgram_inv_doc_freq_cache  A dictionary with q-grams as keys and their
                                   inverse document frequency (IDF) as values.
         max_qgram_count           The maximum count of a q-gram in one record
                                   in this index, used to calculate the
                                   normalised document frequency (DF).

       The cosine similarity calculation as provided in the book Managing
       Gigabytes, page 187 is used, as well as it's variable naming. The
       center of the canopy cluster is assumed to be the query document q.

       Note that due to numerical issues (rounding etc.) the cosine similarity
       values in the case of nearest neighbour canopy will be limited (rounded)
       to 10 digits.
    """

    get_qgram_list_funct = self.__get_qgram_list__

    qgram_list = get_qgram_list_funct(index_val)
    qgram_dict = self.__qgram_list_to_dict__(qgram_list)

    # Build a dictionary of the record identifiers that have at least one
    # q-gram in common with the given index value. Values in the dictionary are
    # the cosine similarity values between the chosen canopy center records and
    # the other records.
    #
    cos_sim_dict = {}

    # Euclidean length (the weight) of the index value (Manag GBytes, page 187)
    # (W_q in the book on page 187)
    #
    W_q = 0.0

    # Loop over the q-grams and their counts in the given index value - - - - -
    #
    for (qgram, qgram_count) in qgram_dict.iteritems():

      # Make sure the q-gram is in the index as it might have been deleted in
      # the build process if it occured too frequent and the self.delete_perc
      # was set
      #
      if (qgram in index):

        qgram_rec_dict = index[qgram]  # All records containing this q-gram

        num_qgram_rec = len(qgram_rec_dict)

        # Get the inverse document frequency of this q-gram from cache
        #
        inv_doc_freq = qgram_inv_doc_freq_cache[qgram]

        # Calculate TF-IDF weight for this q-gram in the index value
        # (W_qt in MG, page 187)
        #
        W_qt = inv_doc_freq * qgram_count / max_qgram_count

        # Add to Euclidean length (the weight) of the index value (W_q)
        #
        W_q += W_qt*W_qt

        # Loop over the record identifiers and their count for this q-gram
        # (note that these counts have already been normalised during the build
        # process)
        #
        for (rec_ident, rec_qgram_count) in qgram_rec_dict.iteritems():

          # Calculate TF-IDF weight for this q-gram in given record (W_dt)
          #
          W_dt = inv_doc_freq * rec_qgram_count

          # Update sum of W_qt*W_dt for this record
          #
          cos_sim = cos_sim_dict.get(rec_ident, 0.0) + W_qt*W_dt
          cos_sim_dict[rec_ident] = cos_sim

    W_q = math.sqrt(W_q)  # Calculate final W_q

    return_list = []  # Record identifiers to be returned
    delete_list = []  # Record identifers to be deleted

    # Now process list of record identifers and their similarities - - - - - -
    # according to canopy method
    #
    if (self.canopy_method[1] == 'threshold'):

      # Get thresholds and precompute scaling
      #
      t_tight = self.canopy_method[2]*W_q  # So not needed in loop below
      t_loose = self.canopy_method[3]*W_q

      # Loop over all record identifiers and their cosine similarities
      #
      for (rec_ident, cos_sim) in cos_sim_dict.iteritems():

        # Make sure the similarity is in the possible range
        #
        assert cos_sim / W_q >= -0.000000001, \
               (cos_sim / W_q, index_val, index_val_cache[rec_ident])
        assert cos_sim / W_q <=  1.000000001, \
               (cos_sim / W_q, index_val, index_val_cache[rec_ident])

        if (cos_sim >= t_loose):  # Within loose threshold, so keep it
          return_list.append(rec_ident)

          if (cos_sim >= t_tight):  # Within tight threshold
            delete_list.append(rec_ident)

      # If return list is empty make sure at least the record with the
      # highest similarity are returned
      #
      if (return_list == []):

        max_cos_val =     -1
        max_cos_val_rec = ''

        for (rec_ident, cos_sim) in cos_sim_dict.iteritems():
          if (cos_sim > max_cos_val):
            max_cos_val =     cos_sim
            max_cos_val_rec = rec_ident

        return_list.append(max_cos_val_rec)

    else:  # TF-IDF nearest

      remove_nearest =  self.canopy_method[2]
      cluster_nearest = self.canopy_method[3]

      # Build a dictionary with cosine similarities as keys
      #
      sim_dict = {}

      # Loop over all record identifiers and their cosine similarities
      #
      for (rec_ident, cos_sim) in cos_sim_dict.iteritems():

        # Due to numerical issues we round the final cosine similarities to 10
        # digits
        #
        round_cos_sim = round(cos_sim,10)  # Round to 10 digits accuracy #####

        sim_rec_ident_list = sim_dict.get(round_cos_sim, [])
        sim_rec_ident_list.append(rec_ident)
        sim_dict[round_cos_sim] = sim_rec_ident_list

      sim_values = sim_dict.keys()
      sim_values.sort(reverse=True)  # Largest values first

      for sim_val in sim_values:
        sim_val_rec_list = sim_dict[sim_val]

        if ((len(return_list) + len(sim_val_rec_list)) <= cluster_nearest):
          return_list += sim_val_rec_list

          if ((len(delete_list) + len(sim_val_rec_list)) <= remove_nearest) \
             or (delete_list == []):  # Make sure delete is not empty
            delete_list += sim_val_rec_list

        else:
          if (return_list == []):  # Make sure at least nearest neighbours are
                                   # returned and deleted
            return_list += sim_val_rec_list
            delete_list += sim_val_rec_list

          break  # Exit loop, enough nearest neighbours found

      del sim_values
      del sim_dict

    cos_sim_dict.clear()

    assert len(return_list) >= len(delete_list)

    # Delete records in the delete list - - - - - - - - - - - - - - - - - - - -
    #
    for rec_ident in delete_list:

      del_index_val = index_val_cache[rec_ident]  # Get value of this record

      # Delete the record in the inverted index for all q-grams of its value
      #
      del_qgram_list = get_qgram_list_funct(del_index_val)

      for qgram in del_qgram_list:
        if (qgram in index):
          if (rec_ident in index[qgram]):
            del index[qgram][rec_ident]

            if (len(index[qgram]) == 0):
              del index[qgram]  # Not needed anymore
              del qgram_inv_doc_freq_cache[qgram]

      # Delete the record identifier in indexing values cache
      #
      index_val_cache.pop(rec_ident)

    del delete_list

    return return_list

  # ---------------------------------------------------------------------------

  def __jaccard_canopy__(self, index, index_val, index_val_cache,
                         index_val_num_qgram):
    """Returns a list of record identifiers in the given inverted index
       according to the Jaccard canopy method. Uses thresholds to calculate how
       many of the index values q-grams have to be in common with other values.

       Arguments provided:
         index                     The inverted index with q-grams as keys and
                                   dictionaries with record identifiers as
                                   values.
         index_val                 Value of the cluster center.
         index_val_cache           A dictionary with record identifiers as keys
                                   and their index value as values.
         index_val_num_qgram       The number of q-grams in the index value for
                                   each record.
    """

    get_qgram_list_funct = self.__get_qgram_list__

    # First get a set with all q-grams in the index value
    #
    qgram_set =  set()

    for qgram in get_qgram_list_funct(index_val):
      if (qgram in index):    # Q-gram might have been deleted in build()
        qgram_set.add(qgram)  # because it was too common

    num_qgrams = len(qgram_set)

    # Get the number of record identifiers in the index for each q-gram - - - -
    #
    qgram_len_list = []

    for qgram in qgram_set:
      qgram_len_list.append((len(index[qgram]), qgram))

    qgram_len_list.sort()  # Sort so fewest q-gram numbers come first in list

    # Depending if thresholds or nearest neighbours are used, get values
    #
    if (self.canopy_method[1] == 'threshold'):
      do_threshold = True

      remove_threshold =  self.canopy_method[2]
      cluster_threshold = self.canopy_method[3]

      phase_switch_threshold = max(1,
                                 1+int((1.0-cluster_threshold)*num_qgrams))

    else:  # Jaccard nearest
      do_threshold = False

      remove_nearest =  self.canopy_method[2]
      cluster_nearest = self.canopy_method[3]

      qgram_count_dict = {}  # Number of records having q-grams in common

      for i in range(1,len(qgram_len_list)+1):
        qgram_count_dict[i] = 0  # Initialisation needed
###      for i in range(len(qgram_len_list)):
###        qgram_count_dict[i+1] = 0  # Initialisation needed

      # The phase switch will be calculated when ebough close records are found
      #
      phase_switch_threshold = len(qgram_len_list)

      jacc_sim_dict = {}  # Will contain jaccard similarity values as keys and
                          # lists of record identifiers as values

    # Create main dictionary with record identifiers and their counts (how many
    # q-grams in common with index value)
    #
    # This is done in 2 phases:
    # - In the first phase (counter below cluster threshold) new record
    #   identifers are added to this dictionary
    # - In the second phase (counter above cluster threshold) no new record
    #   identifers are added, as they would be able to get a count above the
    #   cluster threshold
    #
    qgram_rec_count_dict = {}

    # Main loop over all q-grams and their record identifer lists - - - - - - -
    #
    i = 0  # Count number of q-grams processed

    for (qgram_count, qgram) in qgram_len_list:
      qgram_rec_dict = index[qgram]  # All records containing this q-gram

      if (i < phase_switch_threshold):  # Consider new record identifers

        for rec_ident in qgram_rec_dict:
          rec_count = qgram_rec_count_dict.get(rec_ident, 0) + 1
          qgram_rec_count_dict[rec_ident] = rec_count

          if (do_threshold == False):
            qgram_count_dict[rec_count] = qgram_count_dict[rec_count] + 1

        # For nearest neighbor, need to check if we can switch the phase
        #
        if (do_threshold == False):

          if ((i+1) > len(qgram_len_list)/2):

            # Switch if there are enough q-grams with large count
            #
            if (qgram_count_dict[i+1] > cluster_nearest):  ### PC 23/01: i->i+1
              phase_switch_threshold = i

      # Don't consider new record identifiers as they won't reach the cluster
      # threshold or become nearest neighbours
      #
      else:
        for rec_ident in qgram_rec_count_dict:

          if (rec_ident in qgram_rec_dict):
            rec_ident_count = qgram_rec_count_dict[rec_ident] + 1
            qgram_rec_count_dict[rec_ident] = rec_ident_count

      i += 1

    return_list = []  # Record identifiers to be returned
    delete_list = []  # Record identifiers to be deleted

    # Calculate final thresholds and put into corresponding lists - - - - - - -
    #
    for (rec_ident, rec_ident_count) in qgram_rec_count_dict.iteritems():
      rec_index_val_num_qgrams = index_val_num_qgram[rec_ident]

      # Calculate the Jaccard similarity:
      # - intersection is: rec_ident_count (number of q-grams in common)
      # - union is: |q-grams in index val| + |q-grams in record index val| -
      #             |q-grams in common|
      qgram_union = rec_index_val_num_qgrams + num_qgrams - rec_ident_count

      jacc_sim = float(rec_ident_count) / qgram_union

##        print ' -> ', jacc_sim, index_val_cache[rec_ident], rec_ident
##        print '    ', rec_ident_count, rec_index_val_num_qgrams, num_qgrams

      if (do_threshold == True):

        if (jacc_sim >= cluster_threshold):
          return_list.append(rec_ident)

          if (jacc_sim >= remove_threshold):
           delete_list.append(rec_ident)

      else:  # Nearest, put into a similarity dictionary (sim. values as keys)

        rec_ident_list = jacc_sim_dict.get(jacc_sim, [])
        rec_ident_list.append(rec_ident)
        jacc_sim_dict[jacc_sim] = rec_ident_list

    if (do_threshold == False):  # Extract nearest from similariy dictionary -

      sim_values = jacc_sim_dict.keys()
      sim_values.sort(reverse=True)  # Largest values first

      for sim_val in sim_values:
        sim_val_rec_list = jacc_sim_dict[sim_val]

        if ((len(return_list) + len(sim_val_rec_list)) <= cluster_nearest):
          return_list += sim_val_rec_list

          if ((len(delete_list) + len(sim_val_rec_list)) <= remove_nearest) \
             or (delete_list == []):  # Make sure delete is not empty
            delete_list += sim_val_rec_list

        else:
          if (return_list == []):  # Make sure at least nearest neighbours are
                                   # returned and deleted
            return_list += sim_val_rec_list
            delete_list += sim_val_rec_list

          break  # Exit loop, enough nearest neighbours found

      del jacc_sim_dict
      del sim_values

    # Delete records in the delete list - - - - - - - - - - - - - - - - - - - -
    #
    for rec_ident in delete_list:

      del_index_val = index_val_cache[rec_ident]  # Get value of this record

      # Delete the record in the inverted index for all q-grams of its value
      #
      del_qgram_list = get_qgram_list_funct(del_index_val)

      for qgram in del_qgram_list:
        if (qgram in index):
          if (rec_ident in index[qgram]):
            del index[qgram][rec_ident]

            if (len(index[qgram]) == 0):
              del index[qgram]  # Not needed anymore

      # Delete the record identifier in indexing values cache
      #
      index_val_cache.pop(rec_ident)

    del delete_list

    return return_list

  # ---------------------------------------------------------------------------

  def run(self, length_filter_perc = None, cut_off_threshold = None):
    """Iterate over all blocks in the index.

       Compare the record pairs as produced by the canopy clustering indexing
       process, and return a weight vector dictionary with keys made of a tuple
       (record identifier 1, record identifier 2), and corresponding values the
       comparison weights.
    """

    logging.info('')
    logging.info('Started comparison of %d record pairs' % \
                 (self.num_rec_pairs))
    if (self.log_funct != None):
      self.log_funct('Started comparison of %d record pairs' % \
                     (self.num_rec_pairs))

    # Check if index has been compacted - - - - - - - - - - - - - - - - - - - -
    #
    if (self.status != 'compacted'):
      logging.exception('Index "%s" has not been compacted, running ' % \
                        (self.description)+'comparisons not possible')
      raise Exception

    # Compare the records
    #
    return self.__compare_rec_pairs_from_dict__(length_filter_perc,
                                                cut_off_threshold)

# =============================================================================

class StringMapIndex(Indexing):
  """Class that implements the string-map multidimensional mapping algorithm.

     For details see the following papers:

     - Efficient similarity string joins in large data sets
       Liang Jin, Chen Li and Sharad Methrotra,
       28th VLDB, 2002.

     - Efficient record linkage in large data sets
       Liang Jin, Chen Li and Sharad Methrotra,
       8th International Conference on Database Systems for Advanced
       Applications (DASFAA '03), 2003.

     The basic idea is to map strings into a multi-dimensional space that
     preserves the distances over strings, and to use an efficient tree-based
     data structure to retrieve similar pairs of strings.

     The data structure used to perform the nearest neighbour searching in the
     above papers is a R-tree. We choose an inverted index based approach
     similar to the iGrid index described in:

     - The iGrid Index: Reversing the Dimensionality Curse For Similarity
       Indexing in High Dimensional Space
       Charu C. Aggarwal and Philip S. Yu,
       KDD 2000.

     The additional argument (besides the base class arguments) which has to be
     set when this index is initialised is:

       dim              Dimension of the space into which strings will be
                        mapped.
       sub_dim          Sub-dimension (must be equal to or smaller than dim)
                        which will be the dimension used in the compact()
                        method to extract records from the inverted grid index
                        (record identifier sets from sub_dim dimensions are
                        intersected, and the union of these interested sets
                        will then be used to find nearest strings).
       sim_funct        The string distance function to be used. It is assumed
                        that this function has two input arguments (str1, str2)
                        and returns a normalised numerical similarity between
                        0.0 (totally different strings) and 1.0 (strings are
                        equal). Distances are then calculated as
                        (1.0-similarity value).
       cache_dist       A flag, if set to True then all distance calculations
                        done will be cached. This will speed up the index
                        building process but use more memory. If set to False,
                        distance calculations will not be cached. Default value
                        is True.
       grid_resolution  The inverted grid resolution in each dimensions, has to
                        be a power of 10 number (e.g. 10,100,1000,etc.)
       canopy_method    Determines how the nearest (most similar) strings are
                        extracted into clusters. Possible are:
                          ('threshold', tight_threshold, loose_threshold)
                          ('nearest',   remove_nearest,  cluster_nearest)
                        with parameters:
                        - tight_threshold   All records with a similarity equal
                                            to or larger than this will be
                                            removed from the pool of records
                                            after inserted into a canopy
                                            cluster. This is a normalised
                                            number 0..1 and must be larger than
                                            or equal to the loose threshold.
                        - loose_threshold   All records with a similarity equal
                                            to or larger than this will be
                                            added into a canopy cluster. This
                                            is a normalised number 0..1 and
                                            must be smaller or equal to the
                                            tight threshold.
                        - remove_nearest    Number of strings nearest to the
                                            cluster center string that will be
                                            removed from the pool of strings
                                            after inserted into a canopy
                                            cluster. Must be a positive integer
                                            and smaller than or equal to the
                                            value of cluster nearest.
                        - cluster_nearest   Number of strings nearest to the
                                            cluster center string that will be
                                            added into a canopy cluster. Must
                                            be a positive integer and larger
                                            than or equal to the value of
                                            remove nearest.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'dim', 'sub_dim', 'sim_funct', 'cache_dist' and
       'grid_resolution' arguments first, then call the base class constructor.
    """

    self.dim =              None
    self.sub_dim =          None
    self.sim_funct =        None
    self.cache_dist =       True
    self.grid_resolution =  None
    self.canopy_method =    None

    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword == 'dim'):
        auxiliary.check_is_integer('dim', value)
        auxiliary.check_is_positive('dim', value)
        self.dim = value

      elif (keyword.startswith('sub_d')):
        auxiliary.check_is_integer('sub_dim', value)
        auxiliary.check_is_positive('sub_dim', value)
        self.sub_dim = value

      elif (keyword.startswith('sim_f')):
        auxiliary.check_is_function_or_method('sim_funct', value)
        self.sim_funct = value

      elif (keyword.startswith('cache_d')):
        auxiliary.check_is_flag('cache_dist', value)
        self.cache_dist = value

      elif (keyword.startswith('grid_r')):
        auxiliary.check_is_integer('grid_resolution', value)
        auxiliary.check_is_positive('grid_resolution', value)
        self.grid_resolution = value

      elif (keyword.startswith('canopy_m')):
        auxiliary.check_is_tuple('canopy_method', value)
        self.canopy_method = value

      else:
        base_kwargs[keyword] = value

    Indexing.__init__(self, base_kwargs)  # Initialise base class

    # Make sure necessary attributes are set - - - - - - - - - - - - - - - - -
    #
    auxiliary.check_is_positive('dim', self.dim)
    auxiliary.check_is_positive('sub_dim', self.sub_dim)
    if (self.sub_dim > self.dim):
      logging.exception('Argument "sub_dim" is larger than "dim": %d / %d' % \
                        (self.sub_dim, self.dim))
      raise Exception

    auxiliary.check_is_function_or_method('sim_funct', self.sim_funct)
    auxiliary.check_is_integer('grid_resolution', self.grid_resolution)
    auxiliary.check_is_positive('grid_resolution', self.grid_resolution)
    if (self.grid_resolution not in [10,100,1000,10000]):
      logging.exception('Argument "grid_resolution" is not a power of 10 ' + \
                        'number: %d' % (elf.grid_resolution))
      raise Exception

    # Check if canopy method and parameters given are OK - - - - - - - - - - -
    #
    auxiliary.check_is_not_none('canopy_method', self.canopy_method)
    auxiliary.check_is_tuple('canopy_method', self.canopy_method)

    if (len(self.canopy_method) != 3):
      logging.exception('Canopy method tuple needs three elements: %s' % \
                        (str(self.canopy_method)))
      raise Exception

    if (self.canopy_method[0] not in ['threshold', 'nearest']):
      logging.exception('Illegal canopy retrieval method given (must be ' + \
                        '"threshold" or "nearest"): %s' % \
                        (str(self.canopy_method)))
      raise Exception

    if (self.canopy_method[0] == 'threshold'):
      auxiliary.check_is_normalised('tight_threshold', self.canopy_method[1])
      auxiliary.check_is_normalised('loose_threshold', self.canopy_method[2])

      if (self.canopy_method[1] < self.canopy_method[2]):
        logging.exception('Tight threshold is smaller than loose threshold:' \
                          + '%.3f / %.3f' % (self.canopy_method[1],
                                             self.canopy_method[2]))
        raise Exception

    else:  # Is 'nearest'
      auxiliary.check_is_integer('remove_nearest',   self.canopy_method[1])
      auxiliary.check_is_positive('remove_nearest',  self.canopy_method[1])
      auxiliary.check_is_integer('cluster_nearest',  self.canopy_method[2])
      auxiliary.check_is_positive('cluster_nearest', self.canopy_method[2])

      if (self.canopy_method[1] > self.canopy_method[2]):
        logging.exception('Remove nearest is larger than cluster nearest:' \
                          + '%d / %d' % (self.canopy_method[1],
                                         self.canopy_method[2]))
        raise Exception

    # Create necessary data structures (one per data set) - - - - - - - - - - -
    #
    self.string_list = {}  # Dictionary with lists with strings from data sets
    self.coord =       {}  # String object coordinates, dim x (num. of strings)
    self.grid_index =  {}

    self.m = 5  # Number of iterations in __choose_pivot__() to get two strings

    self.log([('Dimension', self.dim),
              ('Sub-space dimension', self.sub_dim),
              ('Cache distance calculations', self.cache_dist),
              ('Inverted grid resolution', self.grid_resolution),
              ('Canopy method', self.canopy_method),
              ('Similarity function', self.sim_funct)])  # Log a message

  # ---------------------------------------------------------------------------

  def __choose_pivot__(self, h, num_string, string_list, coord):
    """Method to choose two pivot strings on the h-th dimension.
       Returns the indices of the two pivots in the string list.
    """

    get_distance_funct = self.__get_distance__  # Shorthand

    # ind_a = random.randrange(num_string)
    ind_a = 0  # Simply take first

    for i in xrange(self.m):

      max_dist = -1
      ind_b =    -1  # Index of string at maximum distance from str_a

      for j in xrange(num_string):  # Calculate distance to all other strings
        dist = get_distance_funct(h, num_string, string_list, ind_a, j, coord)

        if (dist > max_dist):
          max_dist = dist
          ind_b =    j

      max_dist = -1
      ind_a =    -1

      for j in xrange(num_string):  # Calculate distance to all other strings
        dist = get_distance_funct(h, num_string, string_list, ind_b, j, coord)
        if (dist > max_dist):
          max_dist = dist
          ind_a =    j

    return (ind_a, ind_b)

  # ---------------------------------------------------------------------------

  def __get_distance__(self, h, num_string, string_list, ind1, ind2, coord):
    """Get distance of the two given strings (after strings have been projected
       onto the first h-1 axis).
    """

    self.num_dist_calc += 1

    if ((self.cache_dist == True) and ((h,ind1,ind2) in self.comp_dist_cache)):
      dist = self.comp_dist_cache[(h,ind1,ind2)]

    else:

      str1 = string_list[ind1]
      str2 = string_list[ind2]

      dist = 1.0 - self.sim_funct(str1, str2)  # Calculate string distance

      num_str_count = 0
      for i in xrange(h):
##        w = coord[i*num_string+ind1] - coord[i*num_string+ind2]
        w = coord[num_str_count + ind1] - coord[num_str_count + ind2]
        num_str_count += num_string

        dist = math.sqrt(abs(dist*dist - w*w))

      if (self.cache_dist == True):
        self.comp_dist_cache[(h,ind1,ind2)] = dist  # Save into cache

    return dist

  # ---------------------------------------------------------------------------

  def build(self):
    """Method to build an index data structure.

       First load the records from the data sets and put them into a basic
       inverted index data structure with blocking variables as keys. Next
       map them into a multi-dimensional space.

       If both data sets are the same (a deduplication) only insert records
       from one data set.
    """

    logging.info('')
    logging.info('Build string-map index: "%s"' % (self.description))

    start_time = time.time()

    # First build the basic (blocking) inverted index
    #
    self.__records_into_inv_index__()  # Read records and put into index

    num_indices = len(self.index_def)

    # Initialise the index data structures needed - - - - - - - - - - - - - - -
    #
    for i in range(num_indices):
      self.string_list[i] = {}
      self.coord[i] =       {}
      self.grid_index[i] =  {}

    choose_pivot_funct = self.__choose_pivot__  # Shorthands
    get_distance_funct = self.__get_distance__
    dim =                self.dim

    grid_round_digit = {10:1, 100:2, 1000:3, 10000:4}[self.grid_resolution]

    for i in range(num_indices):

      istart_time = time.time()

      # Get all string values from basic inverted indxex
      #
      string_list = self.index1[i].keys()  # From data set 1

      # For a linkage, build combined list of all string values
      #
      if (self.do_deduplication == False):
        string1_set = set(string_list)
        string2_set = set(self.index2[i].keys())  # From data set 2
        string_list = list(string1_set.union(string2_set))

      self.string_list[i] = string_list  # Save list of all strings

      num_string = len(string_list)  # Number of strings in this index

      coord = self.coord[i]  # Shorthand

      self.comp_dist_cache = {}  # Cache of calculated distances
      self.num_dist_calc =   0   # Count the number of distance calculations

      logging.info('  Map index %d containing %d string values into a ' % \
                   (i, num_string) + '%d-dimensional space' % (dim))

      coord = [0.0]*num_string*dim  # A list interpreted as matrix

      # Start string-map routine - - - - - - - - - - - - - - - - - - - - - - -
      #
      for h in xrange(dim):  # Loop over dimensions

        (p1, p2) = choose_pivot_funct(h, num_string, string_list, coord)
        dist = get_distance_funct(h, num_string, string_list, p1, p2, coord)

        if (dist == 0.0):  # All coordinates in the h-th dimension are 0
          break

        h_num_string = h*num_string  # Shorthands
        dist_square =  dist*dist
        dist_two =     2*dist

        # Calculate coordinates of all strings on this axis
        #
        for j in xrange(num_string):
          x = get_distance_funct(h, num_string, string_list, j, p1, coord)
          y = get_distance_funct(h, num_string, string_list, j, p2, coord)
          coord[h_num_string+j] = (x*x + dist_square - y*y) / dist_two

        logging.info('    Processed dimension %d' % (h))

      # Put coordinates into a data structure for efficient nearest neighbor -
      # search
      #
      index_grid = self.grid_index[i]  # Shorthand

      for h in xrange(dim):  # One dictionary per dimension
        index_grid[h] = {}

      logging.info('  Convert index %d into an inverted grid index' % (i))

      str_coord_dict = {}  # Convert coordinates array into a dictionary with
                           # strings as keys and their coordinates as values

      for j in xrange(num_string):

        str_val =   string_list[j]  # Get the string value
        str_coord = []

        num_str_count = 0
        for h in xrange(dim):  # Loop over dimensions
 ##         this_coord_val = coord[h*num_string+j]
          this_coord_val = coord[num_str_count + j]
          num_str_count += num_string

          str_coord.append(this_coord_val)

          # Put into inverted grid index
          #
          round_coord_val = round(this_coord_val, grid_round_digit)

          # Each grid cell contains a set of string values in this cell
          #
          grid_str_set = index_grid[h].get(round_coord_val, set())
          grid_str_set.add(str_val)
          index_grid[h][round_coord_val] = grid_str_set

        str_coord_dict[str_val] = str_coord

      del coord  # Not needed anymore
      self.coord[i] = str_coord_dict

      logging.info('  Built string-map index %d with %d strings values ' % \
                   (i, num_string)+'in %s' % \
                   (auxiliary.time_string(time.time()-istart_time)))
      logging.info('    Number of distance calculations done: %d' % \
                   (self.num_dist_calc))

      if (self.cache_dist == True):
        logging.info('    Length of cache: %d (in average %.2f distance ' % \
                     (len(self.comp_dist_cache),
                     float(self.num_dist_calc)/len(self.comp_dist_cache)) + \
                     'calculations per cache entry)')
        del self.comp_dist_cache

      logging.info('    Explicitly run garbage collection')
      gc.collect()

      memory_usage_str = auxiliary.get_memory_usage()
      if (memory_usage_str != None):
        logging.info('      '+memory_usage_str)

    logging.info('Built string-map index in %s' % \
                 (auxiliary.time_string(time.time()-start_time)))

    self.status = 'built'  # Update index status

  # ---------------------------------------------------------------------------

  def compact(self):
    """Method to compact an index data structure.

       An approach similar to canopy clustering is used by randomly picking a
       string and then extracting its k nearest neighbours.

       Finally calculate the total number of record pairs.
    """

    NUM_CANOPY_PROGRESS_REPORT = 100  # Log a message every XXX canopies

    logging.info('')
    logging.info('Compact string-map index: "%s"' % (self.description))

    start_time = time.time()

    num_indices = len(self.index_def)

    grid_res = self.grid_resolution
    grid_round_digit = {10:1, 100:2, 1000:3, 10000:4}[grid_res]

    interval_size = 1.0/self.grid_resolution  # To get neighbouring grid cells

    dedup_rec_pairs_funct = self.__dedup_rec_pairs__  # Shorthands
    link_rec_pair_funct =   self.__link_rec_pairs__
    sub_dim =               self.sub_dim
    dim =                   self.dim
    do_dedup =              self.do_deduplication

    if (self.canopy_method[0] == 'nearest'):  # Shorthands
      do_nearest = True
      remove_nearest =  self.canopy_method[1]
      cluster_nearest = self.canopy_method[2]

    else:  # Re-scale similarity measure and make it a distance
      do_nearest = False
      tight_threshold = (1.0-self.canopy_method[1])*math.sqrt(dim)
      loose_threshold = (1.0-self.canopy_method[2])*math.sqrt(dim)

    # Check if index has been built - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.status != 'built'):
      logging.exception('Index "%s" has not been built, compacting is not ' % \
                        (self.description)+'possible')
      raise Exception

    rec_pair_dict = {}  # A dictionary with record identifiers from data set 1
                        # as keys and sets of identifiers from data set 2 as
                        # values

    for i in range(num_indices):

      istart_time = time.time()

      num_canopies = 0  # Count the number of canopies created

      smallest_canopy_size =   999999
      smallest_canopy_center = ''    # Indexing value of the smallest canopy
      largest_canopy_size =    -99999
      largest_canopy_center =  ''    # Indexing value of the largest canopy

      index_grid = self.grid_index[i]  # Shorthand to inverted grid index

      this_index1 = self.index1[i]  # Shorthands to basic inverted index
      if (do_dedup == False):  # A linkage
        this_index2 = self.index2[i]

      # Get the dictionary with strings in this index and their coordinates
      #
      str_coord_dict = self.coord[i]

      str_val_list = str_coord_dict.keys()  # List of the string values only

      # Total number of strings in this index
      #
      total_num_str = len(str_val_list)

      logging.info('  Compacting index %d containing %d strings' % \
                   (i, total_num_str))

      # Loop over all string values, extract canopies and delete strings from
      # string dicionary
      #
      while(len(str_val_list) > 0):

        # Get arbitrary (the first) string and its coordinates
        #
        # str_ind = random.randrange(len(str_val_list))
        center_str_val =   str_val_list[0]
        center_str_coord = str_coord_dict[center_str_val]

        # In case not enough records can be extracted from the selected - - - -
        # string's grid cell and its direct neighbour cells, the search
        # has to be extended to more neighbouring cells.
        # The following flag will be set to True if enough records have been
        # extracted.
        #
        got_enough_records_for_canopy = False

        # Start with only one neighbouring grid cell in each direction
        #
        neighbour_offset = 0

        # Get all close neighbouring strings to the chosen center
        #
        final_candidate_set = set()  # Union of all candiate sets

        dist_comp_cache = {}  # Cache distance comparison made for this canopy

        # Until enough records extracted - - - - - - - - - - - - - - - - - - -
        #
        while (got_enough_records_for_canopy == False) and \
              (neighbour_offset < grid_res/2):

          neighbour_offset += 1

          for h in xrange(dim):  # Loop over dimensions

            # Get the key value of the grid cell where this string is
            #
            this_cell_val = round(center_str_coord[h], grid_round_digit)

            # Get strings in the grid cell and check if center string is there
            # (this is also the initial candidate set)
            #
            this_cand_set = index_grid[h][this_cell_val]
            assert center_str_val in this_cand_set

            # Get neighbouring grid cells
            #
            previous_cell_val = this_cell_val - neighbour_offset*interval_size
            next_cell_val =     this_cell_val + neighbour_offset*interval_size

            if (previous_cell_val in index_grid[h]):
              previous_cell_str_set = index_grid[h][previous_cell_val]
              this_cand_set = this_cand_set.union(previous_cell_str_set)

            if (next_cell_val in index_grid[h]):
              next_cell_str_set = index_grid[h][next_cell_val]
              this_cand_set = this_cand_set.union(next_cell_str_set)

            if ((h % sub_dim) == 0):  # Start new candidate set
              candidate_set = this_cand_set
            else:
              candidate_set = candidate_set.intersection(this_cand_set)

              if (((h+1) % sub_dim) == 0):

                # Union current candidate set with final candidate set
                #
                final_candidate_set = final_candidate_set.union(candidate_set)
                del candidate_set
                candidate_set = set()

          if (len(candidate_set) > 0):
            final_candidate_set = final_candidate_set.union(candidate_set)

          final_candidate_set_len = len(final_candidate_set)

          # Calculate proper distances to each of these strings - - - - - - - -
          #
          dist_dict = {}  # Distances as keys and lists of strings as values

          # Number of string pairs that differ but have distance zero
          #
          num_same_str_loc = 0

          for str_candidate_val in final_candidate_set:

            # Check if distance has already been calculated
            #
            if ((neighbour_offset > 1) and \
                ((center_str_val,str_candidate_val) in dist_comp_cache)):
              edist = dist_comp_cache[(center_str_val,str_candidate_val)]

            else:  # Calculate Euclidean distance

              # Get string coordinates first
              #
              str_candidate_coord = str_coord_dict[str_candidate_val]

              edist = 0.0  # Calculate Euclidean distance

              for h in xrange(dim):  # Loop over dimensions
                dim_diff = abs(str_candidate_coord[h] - center_str_coord[h])
                edist +=   dim_diff*dim_diff
              edist = math.sqrt(edist)

              # Check if strings that differ have non-zero distance
              #
              if ((edist == 0.0) and (center_str_val != str_candidate_val)):
                num_same_str_loc += 1

                # Change distance to a very small value, as otherwise lots of
                # records will be put into the canopy
                #
                edist = 0.0001

              # Cache distance
              #
              dist_comp_cache[(center_str_val,str_candidate_val)] = edist

            # Put distance and string into distance dictionary
            #
            dist_str_list = dist_dict.get(edist, [])
            dist_str_list.append(str_candidate_val)
            dist_dict[edist] = dist_str_list

          if (num_same_str_loc > 0):
            logging.warning('%d string pairs that differ had Euclidean ' % \
                            (num_same_str_loc) + 'distance of 0.0 (set to' + \
                            ' a very small distance)')

          dist_heap = dist_dict.keys()
          heapq.heapify(dist_heap)  # Heap of distances, easier to part. sort

          remove_str_list = []  # String values to be deleted

          canopy_recs1 = [] # All record identifiers from data set 1 in canopy
          canopy_recs2 = [] # Record identifiers from data set 2, linkage only

          # Get and remove strings according to canopy method and - - - - - - -
          # add corresponding record identifiers into canopies
          #
          while (dist_heap != []):

            smallest_dist = heapq.heappop(dist_heap) # Smallest distance value
            smallest_str_list = dist_dict[smallest_dist]  # Its string values

            this_str_val_recs1 = []  # All record identifiers for these strings
            this_str_val_recs2 = []

            for str_val in smallest_str_list:  # Get record ident. of strings

              if (str_val in this_index1):
                this_str_val_recs1 += this_index1[str_val]
              if (do_dedup == False) and (str_val in this_index2):
                this_str_val_recs2 += this_index2[str_val]

            if (do_nearest == True):

              comb_list_len = len(canopy_recs1) + len(canopy_recs2) + \
                              len(this_str_val_recs1) + len(this_str_val_recs2)

              if (comb_list_len <= cluster_nearest):  # Add more to canopy
                canopy_recs1 += this_str_val_recs1
                if (do_dedup == False):
                  canopy_recs2 += this_str_val_recs2

                # Remove string (make sure at least closest will be removed)
                #
                if ((comb_list_len <= remove_nearest) or
                    (len(remove_str_list) == 0)):
                  remove_str_list += smallest_str_list

              else:
                if ((canopy_recs1 == []) and (canopy_recs2 == [])):
                  canopy_recs1 = this_str_val_recs1
                  if (do_dedup == False):
                    canopy_recs2 = this_str_val_recs2
                  remove_str_list = smallest_str_list

                break

            else:  # Thresholds - - - - - - - - - - - - - - - - - - - - - - - -

              if (smallest_dist <= loose_threshold):  # Add to canopies

                canopy_recs1 += this_str_val_recs1
                if (do_dedup == False):  # A linkage
                  canopy_recs2 += this_str_val_recs2

                # Remove string (make sure at least closest will be removed)
                #
                if ((smallest_dist < tight_threshold) or
                    (len(remove_str_list) == 0)):
                  remove_str_list += smallest_str_list

              else:
                if ((canopy_recs1 == []) and (canopy_recs2 == [])):
                  canopy_recs1 = this_str_val_recs1
                  if (do_dedup == False):
                    canopy_recs2 = this_str_val_recs2
                  remove_str_list = smallest_str_list

                break  # Leave loop as threshold is reached

          # Check if the candidate string set was big enough - - - - - - - - -
          #
          if ((len(dist_heap) == 0) and \
              (((do_nearest == True) and (len(canopy_recs1)+len(canopy_recs2) \
                                          < cluster_nearest)) or \
               ((do_nearest == False) and (smallest_dist < loose_threshold)))):
            #logging.warning('Extracted candidate set with %d string(s) ' % \
            #                (final_candidate_set_len)+'was not enough to ' + \
            #                'build canopy (using %d neighbouring cells)' % \
            #                (neighbour_offset))
            pass
          else:
            got_enough_records_for_canopy = True  # Enough records extracted

        del this_str_val_recs1
        del this_str_val_recs2
        del dist_dict
        del dist_heap
        del dist_comp_cache
        del final_candidate_set

        # Make sure center string is in canopy and will be removed - - - - - -
        #
        assert center_str_val in remove_str_list

        # Delete strings from remove list from index and string list - - - - -
        #
        for str_val in remove_str_list:
          str_val_list.remove(str_val)

          if (do_dedup == True):
            del this_index1[str_val]
          else:
            if (str_val in this_index1):
              del this_index1[str_val]
            if (str_val in this_index2):
              del this_index2[str_val]

          str_coord = str_coord_dict[str_val]  # Also remove from grid index

          for h in xrange(dim):  # Loop over dimensions

            this_grid_cell_val = round(str_coord[h], grid_round_digit)
            this_grid_cell_str_set = index_grid[h][this_grid_cell_val]
            this_grid_cell_str_set.remove(str_val)
            if (len(this_grid_cell_str_set) > 0):
              index_grid[h][this_grid_cell_val] = this_grid_cell_str_set
            else:
              del index_grid[h][this_grid_cell_val]  # All strings removed

          # Finally remove from string coordinates dictionary
          #
          del str_coord_dict[str_val]

        del remove_str_list

        num_canopy_rec = len(canopy_recs1+canopy_recs2)
        num_canopies += 1

        if (num_canopy_rec < smallest_canopy_size):
          smallest_canopy_size =   num_canopy_rec
          smallest_canopy_center = center_str_val
        elif (num_canopy_rec > largest_canopy_size):
          largest_canopy_size =   num_canopy_rec
          largest_canopy_center = center_str_val

        # Retrieve all record identifiers for the cluster strings - - - - - - -
        #
        if (do_dedup == True):

          if (len(canopy_recs1) > 1):  # For deduplication at least two records

            # Build record pairs from record identifiers in this canopy
            #
            dedup_rec_pairs_funct(canopy_recs1, rec_pair_dict)

        else:  # A linkage - - - - - - - - - - - - - - - - - - - - - - - - - -

          if ((len(canopy_recs1) > 0) and (len(canopy_recs2) > 0)):
            link_rec_pair_funct(canopy_recs1, canopy_recs2, rec_pair_dict)

        del canopy_recs1
        del canopy_recs2

        # Log progress report every XXX canopies - - - - - - - - - - - - - - -
        #
        if ((num_canopies % NUM_CANOPY_PROGRESS_REPORT) == 0):
          logging.info('    Created %d canopies; %d strings left' % \
                       (num_canopies, len(str_val_list)))
          memory_usage_str = auxiliary.get_memory_usage()
          if (memory_usage_str != None):
            logging.info('      '+memory_usage_str)

      # Delete not needed index data to free-up memory - - - - - - - - - - - -
      #
      this_index1.clear()  # Not needed anymore
      if (do_dedup == False):
        this_index2.clear()

      logging.info('  Compacted canopy index %d in %s' % \
                   (i, auxiliary.time_string(time.time()-istart_time)))
      logging.info('    Produced %d canopies' % (num_canopies))
      logging.info('      Smallest canopy with %d strings and center ' % \
                   (smallest_canopy_size)+'index value: "%s"' % \
                   (smallest_canopy_center))
      logging.info('      Largest canopy with %d strings and center ' % \
                   (largest_canopy_size)+'index value: "%s"' % \
                   (largest_canopy_center))

      logging.info('  Explicitly run garbage collection')
      gc.collect()

      memory_usage_str = auxiliary.get_memory_usage()
      if (memory_usage_str != None):
        logging.info('    '+memory_usage_str)

    num_rec_pairs = 0  # Count lengths of all record identifier sets - - - - -

    for rec_ident2_set in rec_pair_dict.itervalues():
      num_rec_pairs += len(rec_ident2_set)

    self.rec_pair_dict = rec_pair_dict  # Save for later used in run()
    self.num_rec_pairs = num_rec_pairs

    logging.info('Compacted canopy index in %s' % \
                 (auxiliary.time_string(time.time()-start_time)))
    logging.info('  Number of record pairs: %d' % (num_rec_pairs))

    memory_usage_str = auxiliary.get_memory_usage()
    if (memory_usage_str != None):
      logging.info('  '+memory_usage_str)

    self.status = 'compacted'  # Update index status

  # ---------------------------------------------------------------------------

  def run(self, length_filter_perc = None, cut_off_threshold = None):
    """Iterate over all blocks in the index.

       Compare the record pairs as produced by the string map canopy clustering
       indexing process, and return a weight vector dictionary with keys made
       of a tuple (record identifier 1, record identifier 2), and corresponding
       values the comparison weights.
    """

    logging.info('')
    logging.info('Started comparison of %d record pairs' % \
                 (self.num_rec_pairs))
    if (self.log_funct != None):
      self.log_funct('Started comparison of %d record pairs' % \
                     (self.num_rec_pairs))

    # Check if index has been compacted - - - - - - - - - - - - - - - - - - - -
    #
    if (self.status != 'compacted'):
      logging.exception('Index "%s" has not been compacted, running ' % \
                        (self.description)+'comparisons not possible')
      raise Exception

    # Compare the records
    #
    return self.__compare_rec_pairs_from_dict__(length_filter_perc,
                                                cut_off_threshold)

# =============================================================================

class SuffixArrayIndex(Indexing):
  """Class that builds a suffix array on the values in the blocking variables,
     which can then efficiently be processed with different q-gram criterias.

     For details see the following paper:

     - A Fast Linkage Detection Scheme for Multi-Source Information Integration
       Akiko Aizawa and Keizo Oyama
       Proceedings of the 2005 International Workshop on Challenges in Web
       Information Retrieval and Integration (WIRI'05), 2005.

     The additional argument (besides the base class arguments) which has to be
     set when this index is initialised is:

       block_method   A tuple containing the details of how to form the blocks
                      using the strings in the suffix array. This tuple is of
                      the form: (min_suffix_len, max_block_size), with
                      parameters:
                      - min_suffix_len  The minimum length of sub-strings to be
                                        stored in the suffix-array and to be
                                        used as indexing (blocking) values.
                      - max_block_size  The maximum records to be put into a
                                        block.
       padded         If set to True (default), the beginning and end of the
                      indexing values taken from records will be padded with a
                      special start and end characters. If set to False no such
                      padding will be done.
       suffix_method  The way suffix values are generated from the blocking
                      key values. Possible are:
                      - 'suffixonly'  Only the true suffixes are generated, for
                                      example for 'peter', the values 'eter',
                                      'ter', etc. will be generated.
                      - 'allsubstr'   All sub-strings are generated, for
                                      example for 'peter', the values 'pete',
                                      'eter', 'pet', 'ete', 'ter', etc. will
                                      be generated.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'block_method' argument first, then call the
       base class constructor.

       Note that number of record pairs will not be known after initialisation
       (so it is left at value None).
    """

    self.block_method =  None
    self.padded =        True
    self.suffix_method = None

    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('block_m')):
        auxiliary.check_is_tuple('block_method', value)
        self.block_method = value

      elif (keyword.startswith('padd')):
        auxiliary.check_is_flag('padded', value)
        self.padded = value

      elif (keyword.startswith('suffix_m')):
        auxiliary.check_is_string('suffix_method', value)
        if (value not in ['suffixonly', 'allsubstr']):
          logging.exception('Illegal value for "suffix_method": %s ' % \
                            (value) + ' (has to be either "suffixonly" or' + \
                            '"allsubstr")')
          raise Exception
        self.suffix_method = value

      else:
        base_kwargs[keyword] = value

    Indexing.__init__(self, base_kwargs)  # Initialise base class

    # Check if block method and parameters given are OK - - - - - - - - - - - -
    #
    auxiliary.check_is_not_none('block_method', self.block_method)
    auxiliary.check_is_string('suffix_method', self.suffix_method)

    if (len(self.block_method) != 2):
      logging.exception('Blocking method tuple needs two elements: %s' % \
                        (str(self.block_method)))
      raise Exception

    auxiliary.check_is_integer('min_suffix_len', self.block_method[0])
    auxiliary.check_is_positive('min_suffix_len', self.block_method[0])
    auxiliary.check_is_integer('max_block_size', self.block_method[1])
    auxiliary.check_is_positive('max_block_size', self.block_method[1])

    self.log([('Blocking method', self.block_method),
              ('Suffix method', self.suffix_method),
              ('Padded flag', self.padded)])

    self.START_CHAR = chr(1)
    self.END_CHAR =   chr(2)

  # ---------------------------------------------------------------------------

  def build(self):
    """Method to build an index data structure.

       Read the data set(s) from file(s) and insert into suffix arrays (one per
       index defintion), then remove unneeded suffix array entries.
    """

    logging.info('')
    logging.info('Build suffix array index: "%s"' % (self.description))

    start_time = time.time()

    num_indices = len(self.index_def)

    # Index data structure for blocks is one dictionary per index - - - - - - -
    # (a suffix array, with suffix strings as keys and record
    # identifiers as lists)
    #
    for i in range(num_indices):
      self.index1[i] = {}  # Index for data set 1
      self.index2[i] = {}  # Index for data set 2

    get_index_values_funct = self.__get_index_values__  # Shorthands
    skip_missing =           self.skip_missing
    start_char =             self.START_CHAR
    end_char =               self.END_CHAR
    padded =                 self.padded
    min_suffix_len =         self.block_method[0]
    max_block_size =         self.block_method[1]

    # Set a flag for the way suffix array values are generated
    #
    if (self.suffix_method == 'allsubstr'):
      do_all_suffix_str = True
    else:
      do_all_suffix_str = False

    max_suff_str_len = [0]*num_indices  # Record longest suffix strings

    # A list of data structures needed for the build process:
    # - the index data structure (dictionary)
    # - the record cache
    # - the data set to be read
    # - the comparison fields which are used
    # - a list index (0 for data set 1, 1 for data set 2)
    #
    build_list = [(self.index1, self.rec_cache1, self.dataset1,
                   self.comp_field_used1, 0)] # For data set 1

    if (self.do_deduplication == False):  # If linkage append data set 2
      build_list.append((self.index2, self.rec_cache2, self.dataset2,
                   self.comp_field_used2, 1))

    # Reading loop over all records in one or both data set(s) - - - - - - - -
    #
    for (index,rec_cache,dataset,comp_field_used_list,ds_index) in build_list:

      # Calculate a counter for the progress report
      #
      if (self.progress_report != None):
        progress_report_cnt = max(1, int(dataset.num_records / \
                                     (100.0 / self.progress_report)))
      else:  # So no progress report is being logged
        progress_report_cnt = dataset.num_records + 1

      istart_time = time.time()

      rec_read = 0  # Number of records read from data set

      for (rec_ident, rec) in dataset.readall(): # Read all records in data set

        # Extract record fields needed for comparisons (set all others to '')
        #
        comp_rec = []

        field_ind = 0
        for field in rec:
          if (field_ind in comp_field_used_list):
            comp_rec.append(field.lower())
          else:
            comp_rec.append('')
          field_ind += 1

        rec_cache[rec_ident] = comp_rec  # Put into record cache

        # Now get the index variable values for this record - - - - - - - - - -
        #
        rec_index_val_list = get_index_values_funct(rec, ds_index)

        for i in range(num_indices):  # Put record identifier into all indices

          this_index = index[i]  # Shorthand

          index_val = rec_index_val_list[i]

          if ((index_val != '') or (skip_missing == False)):

            if (padded == True):  # Add start and end characters
              index_val = '%s%s%s' % (start_char, index_val, end_char)

            index_val_len = len(index_val)

            max_suff_str_len[i] = max(max_suff_str_len[i], index_val_len)

            # Create all suffix strings (up to minimum length) and insert them
            # into the index

            # A set of all suffix string values for this index value
            #
            this_suffix_str_set = set()
            this_suffix_str_set.add(index_val)  # Even if shorter than min_len

            if (do_all_suffix_str == True):  # Create all sub-string

              # According to Akiko Aizawa (e-mail 8/03/2007) not only suffix
              # strings are generated in their approach, but all sub-strings
              # down to length min_suffix_len

              # Outer loop over all sub-string length
              #
              #for s1 in range(index_val_len-1, min_suffix_len-1, -1):
              for s1 in range(min_suffix_len, index_val_len+1):
                # Inner loop over all possible sub-strings with this length
                #
                for s2 in range(index_val_len-s1+1):
                  this_suffix_str = index_val[s2:s2+s1]
                  this_suffix_str_set.add(this_suffix_str)
            else:  # Only create the true suffixes of the index value string

              for s in range(index_val_len-min_suffix_len+1):
                this_suffix_str = index_val[s:]
                this_suffix_str_set.add(this_suffix_str)

            # Now insert all suffix strings into index
            #
            for this_suffix_str in this_suffix_str_set:

              suffix_str_rec_list = this_index.get(this_suffix_str, [])

              if (suffix_str_rec_list != -1):  # Not too many records yet

                # Check if max_block records had this string value before
                #
                if (len(suffix_str_rec_list) < max_block_size):
                  suffix_str_rec_list.append(rec_ident)
                  this_index[this_suffix_str] = suffix_str_rec_list
                else:  # Too many records have this value
                  this_index[this_suffix_str] = -1 # Mark as being too frequent

            del this_suffix_str_set

        rec_read += 1

        if ((rec_read % progress_report_cnt) == 0):
          self.__log_build_progress__(rec_read,dataset.num_records,start_time)

      used_sec_str = auxiliary.time_string(time.time()-istart_time)
      rec_time_str = auxiliary.time_string((time.time()-istart_time) / \
                                           dataset.num_records)
      logging.info('Read and indexed %d records in %s (%s per record)' % \
                   (dataset.num_records, used_sec_str, rec_time_str))
      logging.info('')

    # Now remove unneeded entries in suffix array strings - - - - - - - - - - -
    #
    self.suffix_array_strings1 = []
    self.suffix_array_strings2 = []

    logging.info('Removing unneeded suffix array strings:')

    for i in range(num_indices):

      # For a deduplication, all suffix array string values that have only one
      # record can be removed, as are all entries containing more than
      # 'max_block_size' records)
      #
      this_index = self.index1[i]

      num_original = len(this_index)
      num_removed_1 =      0
      num_removed_large =  0

      this_suff_array_strings = this_index.keys()

      for s in this_suff_array_strings:

        if (this_index[s] == -1):  # This entry was too frequent
          del this_index[s]
          num_removed_large += 1
        else:
          this_block_len = len(this_index[s])

          # Check if length is 1 for a deduplication only
          #
          if (self.do_deduplication == True) and (this_block_len == 1):
            del this_index[s]
            num_removed_1 += 1

      if (self.do_deduplication == True):
        logging.info('  Removed %d (of %d) strings from suffix array 1 as ' % \
                     (num_removed_1,num_original) + 'they only had one record')
      logging.info('  Removed %d (of %d) strings from suffix array 1 as ' % \
                   (num_removed_large, num_original) + 'they had more ' + \
                   'than %d (maximum block size) records' % (max_block_size))

      this_suff_array_strings = this_index.keys()
      self.suffix_array_strings1.append(this_suff_array_strings)
      logging.info('  Suffix array in index %d for data set 1 contains %d ' \
                   % (i, len(this_suff_array_strings)) + 'strings')

      # Linkage (only remove entries with too many records) - - - - - - - - - -
      #
      if (self.do_deduplication == False):
        this_index = self.index2[i]

        num_original = len(this_index)
        num_removed_large = 0

        this_suff_array_strings = this_index.keys()

        for s in this_suff_array_strings:

          if (this_index[s] == -1):  # This entry was too frequent
            del this_index[s]
            num_removed_large += 1

        logging.info('  Removed %d (of %d) strings from suffix array 2 as ' % \
                     (num_removed_large, num_original) + 'they had more ' + \
                     'than %d (maximum block size) records' % (max_block_size))

        this_suff_array_strings = this_index.keys()
        self.suffix_array_strings2.append(this_suff_array_strings)
        logging.info('  Suffix array in index %d for data set 2 contains %d ' \
                     % (i, len(this_suff_array_strings)) + 'strings')

      logging.info('    Longest string in this suffix array is %d characters' \
                   % (max_suff_str_len[i]) + ' long')

    logging.info('Built suffix array index in %s' % \
                 (auxiliary.time_string(time.time()-start_time)))

    memory_usage_str = auxiliary.get_memory_usage()
    if (memory_usage_str != None):
      logging.info('  '+memory_usage_str)

    self.status = 'built'  # Update index status

  # ---------------------------------------------------------------------------

  def compact(self):
    """Method to compact an index data structure.

       Go through all entries in the sorted suffix array, and build blocks of
       records according to the parameter 'max_block_size'. Put all resulting
       record pairs into a dictionary.

       Finally calculate the total number of record pairs.
    """

    logging.info('')
    logging.info('Compact suffix array index: "%s"' % (self.description))

    start_time = time.time()

    # Check if index has been built - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.status != 'built'):
      logging.exception('Index "%s" has not been built, compacting is not ' % \
                        (self.description)+'possible')
      raise Exception

    dedup_rec_pair_funct = self.__dedup_rec_pairs__  # Shorthands
    link_rec_pair_funct =  self.__link_rec_pairs__

    num_indices = len(self.index_def)

    rec_pair_dict = {}  # A dictionary with record identifiers from data set 1
                        # as keys and sets of identifiers from data set 2 as
                        # values

    for i in range(num_indices):

      istart_time = time.time()

      num_strings_done = 0
      largest_block =    0

      if (self.do_deduplication == True):  # A deduplication - - - - - - - - -

        this_str_list = self.suffix_array_strings1[i]  # Shorthands
        this_index =    self.index1[i]

        this_str_list_len = len(this_str_list)

        # Calculate a counter for the progress report
        #
        if (self.progress_report != None):
          progress_report_cnt = max(1, int(this_str_list_len / \
                                       (100.0 / self.progress_report)))
        else:  # So no progress report is being logged
          progress_report_cnt = this_str_list_len + 1

        for s in this_str_list:

          curr_str_val_records = this_index[s]  # Get records for this string
          curr_block_size =      len(curr_str_val_records)
          largest_block =        max(largest_block, curr_block_size)

          dedup_rec_pair_funct(curr_str_val_records, rec_pair_dict)

          num_strings_done += 1

          if ((num_strings_done % progress_report_cnt) == 0):
            perc_done = int(round(100.0*num_strings_done / this_str_list_len))
            logging.info('    Processed %d of %d (%d%%) strings' % \
                         (num_strings_done, this_str_list_len, perc_done))
            memory_usage_str = auxiliary.get_memory_usage()
            if (memory_usage_str != None):
              logging.info('      '+memory_usage_str)

          del curr_str_val_records

      else:  # A linkage - - - - - - - - - - - - - - - - - - - - - - - - - - -

        this_str_list1 = self.suffix_array_strings1[i]  # Shorthands
        this_str_list2 = self.suffix_array_strings2[i]
        this_index1 =    self.index1[i]
        this_index2 =    self.index2[i]

        this_str_list_len = len(this_str_list1)

        # Calculate a counter for the progress report
        #
        if (self.progress_report != None):
          progress_report_cnt = max(1, int(this_str_list_len / \
                                       (100.0 / self.progress_report)))
        else:  # So no progress report is being logged
          progress_report_cnt = this_str_list_len + 1

        for s in this_str_list1:

          curr_str_val_records1 = this_index1[s] # Suffix array from data set 1
          curr_block_size1 =      len(curr_str_val_records1)

          if s in this_str_list2:  # String is in both suffix arrays

            curr_str_val_records2 = this_index2[s]  # From data set 2
            curr_block_size2 =      len(curr_str_val_records2)
            largest_block =         max(largest_block, \
                                        curr_block_size1 + curr_block_size2)

            link_rec_pair_funct(curr_str_val_records1, curr_str_val_records2,
                                rec_pair_dict)

            del curr_str_val_records2

          num_strings_done += 1

          if ((num_strings_done % progress_report_cnt) == 0):
            perc_done = int(round(100.0*num_strings_done / this_str_list_len))
            logging.info('    Processed %d of %d (%d%%) strings' % \
                         (num_strings_done, this_str_list_len, perc_done))
            memory_usage_str = auxiliary.get_memory_usage()
            if (memory_usage_str != None):
              logging.info('      '+memory_usage_str)

          del curr_str_val_records1

      logging.info('  Compacted suffix array index %d in %s' % \
                   (i, auxiliary.time_string(time.time()-istart_time)))
      logging.info('    Largest block contained %d records' % (largest_block))

      self.index1[i].clear()  # Not needed anymore
      self.index2[i].clear()

      logging.info('    Explicitly run garbage collection')
      gc.collect()

      memory_usage_str = auxiliary.get_memory_usage()
      if (memory_usage_str != None):
        logging.info('      '+memory_usage_str)

    self.rec_pair_dict = rec_pair_dict

    self.num_rec_pairs = 0  # Count lengths of all record identifier sets - - -

    for rec_ident2_set in self.rec_pair_dict.itervalues():
      self.num_rec_pairs += len(rec_ident2_set)

    logging.info('Compacted suffix array index in %s' % \
                 (auxiliary.time_string(time.time()-start_time)))
    logging.info('  Number of record pairs: %d' % (self.num_rec_pairs))

    self.status = 'compacted'  # Update index status

  # ---------------------------------------------------------------------------

  def run(self, length_filter_perc = None, cut_off_threshold = None):
    """Iterate over all blocks in the index.

       Compare the record pairs as produced by the suffix array indexing
       process, and return a weight vector dictionary with keys made of a tuple
       (record identifier 1, record identifier 2), and corresponding values the
       comparison weights.
    """

    logging.info('')
    logging.info('Started comparison of %d record pairs' % \
                 (self.num_rec_pairs))
    if (self.log_funct != None):
      self.log_funct('Started comparison of %d record pairs' % \
                     (self.num_rec_pairs))

    # Check if index has been compacted - - - - - - - - - - - - - - - - - - - -
    #
    if (self.status != 'compacted'):
      logging.exception('Index "%s" has not been compacted, running ' % \
                        (self.description)+'comparisons not possible')
      raise Exception

    # Compare the records
    #
    return self.__compare_rec_pairs_from_dict__(length_filter_perc,
                                                cut_off_threshold)


# =============================================================================

class RobustSuffixArrayIndex(Indexing):
  """Class that builds a suffix array based index on the values in the blocking
     variables, which can then efficiently be processed with different q-gram
     criterias.

     For details see the following paper:

     - Robust Record Linkage Blocking using Suffix Arrays
       Timothy de Vries, Hui Ke, Sanjay Chawla and Peter Christen.
       Proceedings of the ACM Conference on Information and Knowledge
       Management (CIKM), Hong Kong, November 2009.

     The additional argument (besides the base class arguments) which has to be
     set when this index is initialised is:

       block_method   A tuple containing the details of how to form the blocks
                      using the strings in the suffix array. This tuple is of
                      the form: (min_suffix_len, max_block_size), with
                      parameters:
                      - min_suffix_len  The minimum length of sub-strings to
                                        be stored in the suffix-array and to
                                        be used as indexing (blocking) values.
                      - max_block_size  The maximum records to be put into a
                                        block.
       padded         If set to True (default), the beginning and end of the
                      indexing values taken from records will be padded with a
                      special start and end characters. If set to False no
                      such padding will be done.
       str_cmp_funct  A function to compare two strings (as implemented in the
                      stringcmp module).
       str_cmp_thres  The threshold for the string comparison function, must
                      be in (0..1).
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the index specific arguments first, then call the
       base class constructor.

       Note that number of record pairs will not be known after initialisation
       (so it is left at value None).
    """

    self.block_method =  None
    self.padded =        True
    self.str_cmp_funct = None
    self.str_cmp_thres = None

    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('block_m')):
        auxiliary.check_is_tuple('block_method', value)
        self.block_method = value

      elif (keyword.startswith('padd')):
        auxiliary.check_is_flag('padded', value)
        self.padded = value

      elif (keyword.startswith('str_cmp_f')):
        auxiliary.check_is_function_or_method('str_cmp_funct', value)
        self.str_cmp_funct = value

      elif (keyword.startswith('str_cmp_t')):
        auxiliary.check_is_normalised('str_cmp_thres', value)
        self.str_cmp_thres = value

      else:
        base_kwargs[keyword] = value

    Indexing.__init__(self, base_kwargs)  # Initialise base class

    # Check if block method and parameters given are OK - - - - - - - - - - - -
    #
    auxiliary.check_is_not_none('block_method', self.block_method)
    auxiliary.check_is_function_or_method('str_cmp_funct', self.str_cmp_funct)
    auxiliary.check_is_normalised('str_cmp_thres', self.str_cmp_thres)

    if (len(self.block_method) != 2):
      logging.exception('Blocking method tuple needs two elements: %s' % \
                        (str(self.block_method)))
      raise Exception

    auxiliary.check_is_integer('max_block_size', self.block_method[1])
    auxiliary.check_is_positive('max_block_size', self.block_method[1])

    self.log([('Blocking method',             self.block_method),
              ('String comparison function',  self.str_cmp_funct),
              ('String comparison threshold', self.str_cmp_thres),
              ('Padded flag',                 self.padded)])

    self.START_CHAR = chr(1)
    self.END_CHAR =   chr(2)

  # ---------------------------------------------------------------------------

  def build(self):
    """Method to build an index data structure.

       This is the same as the build method for the SuffixArrayIndex.
    """

    logging.info('')
    logging.info('Build robust suffix array index: "%s"' % (self.description))

    start_time = time.time()

    num_indices = len(self.index_def)

    # Index data structure for blocks is one dictionary per index - - - - - - -
    # (a suffix array, with suffix strings as keys and record
    # identifiers as lists)
    #
    for i in range(num_indices):
      self.index1[i] = {}  # Index for data set 1
      self.index2[i] = {}  # Index for data set 2

    get_index_values_funct = self.__get_index_values__  # Shorthands
    skip_missing =           self.skip_missing
    start_char =             self.START_CHAR
    end_char =               self.END_CHAR
    padded =                 self.padded
    min_suffix_len =         self.block_method[0]
    max_block_size =         self.block_method[1]

    max_suff_str_len = [0]*num_indices  # Record longest suffix strings

    # A list of data structures needed for the build process:
    # - the index data structure (dictionary)
    # - the record cache
    # - the data set to be read
    # - the comparison fields which are used
    # - a list index (0 for data set 1, 1 for data set 2)
    #
    build_list = [(self.index1, self.rec_cache1, self.dataset1,
                   self.comp_field_used1, 0)] # For data set 1

    if (self.do_deduplication == False):  # If linkage append data set 2
      build_list.append((self.index2, self.rec_cache2, self.dataset2,
                   self.comp_field_used2, 1))

    # Reading loop over all records in one or both data set(s) - - - - - - - -
    #
    for (index,rec_cache,dataset,comp_field_used_list,ds_index) in build_list:

      # Calculate a counter for the progress report
      #
      if (self.progress_report != None):
        progress_report_cnt = max(1, int(dataset.num_records / \
                                     (100.0 / self.progress_report)))
      else:  # So no progress report is being logged
        progress_report_cnt = dataset.num_records + 1

      istart_time = time.time()

      rec_read = 0  # Number of records read from data set

      for (rec_ident, rec) in dataset.readall(): # Read all records in data set

        # Extract record fields needed for comparisons (set all others to '')
        #
        comp_rec = []

        field_ind = 0
        for field in rec:
          if (field_ind in comp_field_used_list):
            comp_rec.append(field.lower())
          else:
            comp_rec.append('')
          field_ind += 1

        rec_cache[rec_ident] = comp_rec  # Put into record cache

        # Now get the index variable values for this record - - - - - - - - - -
        #
        rec_index_val_list = get_index_values_funct(rec, ds_index)

        for i in range(num_indices):  # Put record identifier into all indices

          this_index = index[i]  # Shorthand

          index_val = rec_index_val_list[i]

          if ((index_val != '') or (skip_missing == False)):

            if (padded == True):  # Add start and end characters
              index_val = '%s%s%s' % (start_char, index_val, end_char)

            index_val_len = len(index_val)

            max_suff_str_len[i] = max(max_suff_str_len[i], index_val_len)

            # Create all suffix strings (up to minimum length) and insert them
            # into the index

            # A set of all suffix string values for this index value
            #
            this_suffix_str_set = set()
            this_suffix_str_set.add(index_val)  # Even if shorter than min_len

            # Create the true suffixes of the index value string
            #
            for s in range(index_val_len-min_suffix_len+1):
              this_suffix_str = index_val[s:]
              this_suffix_str_set.add(this_suffix_str)

            # Now insert all suffix strings into index
            #
            for this_suffix_str in this_suffix_str_set:

              suffix_str_rec_list = this_index.get(this_suffix_str, [])

              if (suffix_str_rec_list != -1):  # Not too many records yet

                # Check if max_block records had this string value before
                #
                if (len(suffix_str_rec_list) < max_block_size):
                  suffix_str_rec_list.append(rec_ident)
                  this_index[this_suffix_str] = suffix_str_rec_list
                else:  # Too many records have this value
                  this_index[this_suffix_str] = -1 # Mark as being too frequent

            del this_suffix_str_set

        rec_read += 1

        if ((rec_read % progress_report_cnt) == 0):
          self.__log_build_progress__(rec_read,dataset.num_records,start_time)

      used_sec_str = auxiliary.time_string(time.time()-istart_time)
      rec_time_str = auxiliary.time_string((time.time()-istart_time) / \
                                           dataset.num_records)
      logging.info('Read and indexed %d records in %s (%s per record)' % \
                   (dataset.num_records, used_sec_str, rec_time_str))
      logging.info('')

    # Now remove unneeded entries in suffix array strings - - - - - - - - - - -
    #
    self.suffix_array_strings1 = []
    self.suffix_array_strings2 = []

    logging.info('Removing unneeded suffix array strings:')

    for i in range(num_indices):

      # For a deduplication, all suffix array string values that have only one
      # record can be removed, as are all entries containing more than
      # 'max_block_size' records)
      #
      this_index = self.index1[i]

      num_original = len(this_index)
      num_removed_1 =      0
      num_removed_large =  0

      this_suff_array_strings = this_index.keys()

      for s in this_suff_array_strings:

        if (this_index[s] == -1):  # This entry was too frequent
          del this_index[s]
          num_removed_large += 1
        else:
          this_block_len = len(this_index[s])

          # Check if length is 1 for a deduplication only
          #
          if (self.do_deduplication == True) and (this_block_len == 1):
            del this_index[s]
            num_removed_1 += 1

      if (self.do_deduplication == True):
        logging.info('  Removed %d (of %d) strings from suffix array 1 as ' % \
                     (num_removed_1,num_original) + 'they only had one record')
      logging.info('  Removed %d (of %d) strings from suffix array 1 as ' % \
                   (num_removed_large, num_original) + 'they had more ' + \
                   'than %d (maximum block size) records' % (max_block_size))

      this_suff_array_strings = this_index.keys()
      self.suffix_array_strings1.append(this_suff_array_strings)
      logging.info('  Suffix array in index %d for data set 1 contains %d ' \
                   % (i, len(this_suff_array_strings)) + 'strings')

      # Linkage (only remove entries with too many records) - - - - - - - - - -
      #
      if (self.do_deduplication == False):
        this_index = self.index2[i]

        num_original = len(this_index)
        num_removed_large = 0

        this_suff_array_strings = this_index.keys()

        for s in this_suff_array_strings:

          if (this_index[s] == -1):  # This entry was too frequent
            del this_index[s]
            num_removed_large += 1

        logging.info('  Removed %d (of %d) strings from suffix array 2 as ' % \
                     (num_removed_large, num_original) + 'they had more ' + \
                     'than %d (maximum block size) records' % (max_block_size))

        this_suff_array_strings = this_index.keys()
        self.suffix_array_strings2.append(this_suff_array_strings)
        logging.info('  Suffix array in index %d for data set 2 contains %d ' \
                     % (i, len(this_suff_array_strings)) + 'strings')

      logging.info('    Longest string in this suffix array is %d characters' \
                   % (max_suff_str_len[i]) + ' long')

    logging.info('Built suffix array index in %s' % \
                 (auxiliary.time_string(time.time()-start_time)))

    memory_usage_str = auxiliary.get_memory_usage()
    if (memory_usage_str != None):
      logging.info('  '+memory_usage_str)

    self.status = 'built'  # Update index status

  # ---------------------------------------------------------------------------

  def compact(self):
    """Method to compact an index data structure.

       Sort suffix array, then do approximate string comparisons and merge
       lists if suffixarray values are similar.

       Also remove resulting blocks that are larger than 'max_block_size'.

       Then put all resulting record pairs into a dictionary.

       Finally calculate the total number of record pairs.
    """

    logging.info('')
    logging.info('Compact robust suffix array index: "%s"' % \
                 (self.description))

    start_time = time.time()

    # Check if index has been built - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.status != 'built'):
      logging.exception('Index "%s" has not been built, compacting is not ' % \
                        (self.description)+'possible')
      raise Exception

    dedup_rec_pair_funct = self.__dedup_rec_pairs__  # Shorthands
    link_rec_pair_funct =  self.__link_rec_pairs__

    max_block_size = self.block_method[1]
    str_cmp_funct =  self.str_cmp_funct
    str_cmp_thres =  self.str_cmp_thres

    num_indices = len(self.index_def)

    rec_pair_dict = {}  # A dictionary with record identifiers from data set 1
                        # as keys and sets of identifiers from data set 2 as
                        # values

    for i in range(num_indices):

      istart_time = time.time()

      num_strings_done = 0
      largest_block =    0

      if (self.do_deduplication == True):  # A deduplication - - - - - - - - -

        this_str_list = self.suffix_array_strings1[i]  # Shorthands
        this_index =    self.index1[i]

        this_str_list_len = len(this_str_list)
        this_str_list.sort()  # Sort all suffixes alphabetically

        # Create a new index and a new string list where similar strings (and
        # their record lists) are merged
        #
        merged_index =    {}
        merged_str_list = []

        j = 0
        while (j < (this_str_list_len-1)):
          this_str = this_str_list[j]  # Get current string
          k = j+1  # Compare with following strings until string similarity
                   # is below given threshold

          while ((k < this_str_list_len) and \
                 (str_cmp_funct(this_str, this_str_list[k]) >= str_cmp_thres)):
            k += 1

          if ((j+1) == k):  # Case 1: No merger, simply copy into merged index
            merged_index[this_str] = this_index[this_str]
            merged_str_list.append(this_str)

          else:  # Merge strings
            merged_str = this_str+'-'+this_str_list[k-1]  # New merged string
                                              # (first and last string value)
            merged_rec_set = set(this_index[this_str])

            for l in range(j+1,k):
              other_str = this_str_list[l]
              merged_rec_set = merged_rec_set.union(set(this_index[other_str]))

            # Do not store if too large (check with Tim deVries!) *************
            #
##            if (len(merged_rec_set) <= max_block_size):
            merged_index[merged_str] = list(merged_rec_set)
            merged_str_list.append(merged_str)

          j = k  # k is the first not merged string

        merged_str_list_len = len(merged_str_list)

        # Calculate a counter for the progress report
        #
        if (self.progress_report != None):
          progress_report_cnt = max(1, int(merged_str_list_len / \
                                       (100.0 / self.progress_report)))
        else:  # So no progress report is being logged
          progress_report_cnt = merged_str_list_len + 1

        for s in merged_str_list:

          curr_str_val_records = merged_index[s]  # Get records for this string
          curr_block_size =      len(curr_str_val_records)
          largest_block =        max(largest_block, curr_block_size)

          dedup_rec_pair_funct(curr_str_val_records, rec_pair_dict)

          num_strings_done += 1

          if ((num_strings_done % progress_report_cnt) == 0):
            perc_done = int(round(100.0*num_strings_done / merged_str_list_len))
            logging.info('    Processed %d of %d (%d%%) strings' % \
                         (num_strings_done, merged_str_list_len, perc_done))
            memory_usage_str = auxiliary.get_memory_usage()
            if (memory_usage_str != None):
              logging.info('      '+memory_usage_str)

          del curr_str_val_records

      else:  # A linkage - - - - - - - - - - - - - - - - - - - - - - - - - - -

        this_str_list1 = self.suffix_array_strings1[i]  # Shorthands
        this_str_list2 = self.suffix_array_strings2[i]
        this_index1 =    self.index1[i]
        this_index2 =    self.index2[i]

        # Create a combined sorted string list
        #
        comb_str_list = list(set(this_str_list1).union(set(this_str_list2)))

        comb_str_list_len = len(comb_str_list)
        comb_str_list.sort()

        # Create a new index and a new string list where similar strings (and
        # their record lists) are merged
        #
        merged_index =    {}  # Will have two list, one each per data set
        merged_str_list = []

        j = 0
        while (j < (comb_str_list_len-1)):
          this_str = comb_str_list[j]  # Get current string

          k = j+1  # Compare with following strings until string similarity
                   # is below given threshold

          while ((k < comb_str_list_len) and \
                 (str_cmp_funct(this_str, comb_str_list[k]) >= str_cmp_thres)):
            k += 1

          if ((j+1) == k):  # Case 1: No merger, simply copy into merged index

            # Only insert into merged index if there are records in this entry
            # from both data sets
            #
            this_rec_list1 = this_index1.get(this_str, [])
            this_rec_list2 = this_index2.get(this_str, [])
            if ((this_rec_list1 != []) and (this_rec_list2 != [])):
              merged_index[this_str] = (this_rec_list1, this_rec_list2)
              merged_str_list.append(this_str)

          else:  # Merge strings
            merged_str = this_str+'-'+comb_str_list[k-1]  # New merged string
                                              # (first and last string value)
            merged_rec_set1 = set(this_index1.get(this_str, []))
            merged_rec_set2 = set(this_index2.get(this_str, []))

            for l in range(j+1,k):
              other_str = comb_str_list[l]
              new_rec_set1 = set(this_index1.get(other_str, []))
              new_rec_set2 = set(this_index2.get(other_str, []))
              merged_rec_set1 = merged_rec_set1.union(new_rec_set1)
              merged_rec_set2 = merged_rec_set2.union(new_rec_set2)

            # Do not store if too large (check with Tim deVries!) *************
            #
##            if (len(merged_rec_set) <= max_block_size):

            # Only insert into merged index if there are records in this entry
            # from both data sets
            #
            if ((len(merged_rec_set1) > 0) and (len(merged_rec_set2) > 0)):
              merged_index[merged_str] = (list(merged_rec_set1),
                                          list(merged_rec_set2))
              merged_str_list.append(merged_str)

          j = k  # k is the first not merged string

        merged_str_list_len = len(merged_str_list)

        # Calculate a counter for the progress report
        #
        if (self.progress_report != None):
          progress_report_cnt = max(1, int(merged_str_list_len / \
                                       (100.0 / self.progress_report)))
        else:  # So no progress report is being logged
          progress_report_cnt = merged_str_list_len + 1

        for s in merged_str_list:

          curr_str_val_records = merged_index[s]  # Get records for this string
          curr_str_val_records1 = curr_str_val_records[0]
          curr_str_val_records2 = curr_str_val_records[1]
          largest_block = max(largest_block, \
                       len(curr_str_val_records1) + len(curr_str_val_records2))

          link_rec_pair_funct(curr_str_val_records1, curr_str_val_records2,
                              rec_pair_dict)

          num_strings_done += 1

          if ((num_strings_done % progress_report_cnt) == 0):
            perc_done = int(round(100.0*num_strings_done / merged_str_list_len))
            logging.info('    Processed %d of %d (%d%%) strings' % \
                         (num_strings_done, merged_str_list_len, perc_done))
            memory_usage_str = auxiliary.get_memory_usage()
            if (memory_usage_str != None):
              logging.info('      '+memory_usage_str)

          del curr_str_val_records1, curr_str_val_records2

      logging.info('  Compacted suffix array index %d in %s' % \
                   (i, auxiliary.time_string(time.time()-istart_time)))
      logging.info('    Largest block contained %d records' % (largest_block))

      self.index1[i].clear()  # Not needed anymore
      self.index2[i].clear()

      logging.info('    Explicitly run garbage collection')
      gc.collect()

      memory_usage_str = auxiliary.get_memory_usage()
      if (memory_usage_str != None):
        logging.info('      '+memory_usage_str)

    self.rec_pair_dict = rec_pair_dict

    self.num_rec_pairs = 0  # Count lengths of all record identifier sets - - -

    for rec_ident2_set in self.rec_pair_dict.itervalues():
      self.num_rec_pairs += len(rec_ident2_set)

    logging.info('Compacted suffix array index in %s' % \
                 (auxiliary.time_string(time.time()-start_time)))
    logging.info('  Number of record pairs: %d' % (self.num_rec_pairs))

    self.status = 'compacted'  # Update index status

  # ---------------------------------------------------------------------------

  def run(self, length_filter_perc = None, cut_off_threshold = None):
    """Iterate over all blocks in the index.

       Compare the record pairs as produced by the suffix array indexing
       process, and return a weight vector dictionary with keys made of a tuple
       (record identifier 1, record identifier 2), and corresponding values the
       comparison weights.
    """

    logging.info('')
    logging.info('Started comparison of %d record pairs' % \
                 (self.num_rec_pairs))
    if (self.log_funct != None):
      self.log_funct('Started comparison of %d record pairs' % \
                     (self.num_rec_pairs))

    # Check if index has been compacted - - - - - - - - - - - - - - - - - - - -
    #
    if (self.status != 'compacted'):
      logging.exception('Index "%s" has not been compacted, running ' % \
                        (self.description)+'comparisons not possible')
      raise Exception

    # Compare the records
    #
    return self.__compare_rec_pairs_from_dict__(length_filter_perc,
                                                cut_off_threshold)


# =============================================================================

class BigMatchIndex(Indexing):
  """Class that uses an idea implemented in the US Census Bureau linkage
     program BigMatch: Only index the smaller file, and process all records in
     the larger file as they are read from file.

     This index can only be used for linkages but not for deduplications.

     Records that have the same index variable values for an index are put into
     the same blocks, and only records within a block are then compared.

     The additional argument (besides the base class arguments) which has to be
     set when this index is initialised is:

       block_method  Determines the blocking method used. It can be set to:
                       ('block'),
                       ('sort',window_size), or
                       ('qgram',q,padded,threshold)
                     with parameters similar to the corresponding indexing
                     methods.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'block_method' argument first, then call the
       base class constructor.

       Note that number of record pairs will not be known after initialisation
       (so it is left at value None).
    """

    self.block_method = None

    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('block_m')):
        auxiliary.check_is_tuple('block_method', value)
        self.block_method = value

      else:
        base_kwargs[keyword] = value

    Indexing.__init__(self, base_kwargs)  # Initialise base class

    if (self.do_deduplication == True):
      logging.exception('BigMatchIndex can only be used for linkages, but ' + \
                        'not deduplications')
      raise Exception

    # Check if block method and parameters given are OK - - - - - - - - - - - -
    #
    auxiliary.check_is_not_none('block_method', self.block_method)

    if (self.block_method[0] == 'block'):
      if (len(self.block_method) != 1):
        logging.exception('Blocking method "block" does not need parameters:' \
                          + ' %s' % (str(self.block_method)))
        raise Exception

    elif (self.block_method[0] == 'sort'):
      if (len(self.block_method) != 2):
        logging.exception('Blocking method "sort" needs one parameter: %s' % \
                          (str(self.block_method)))
        raise Exception
      auxiliary.check_is_integer('window size', self.block_method[1])
      auxiliary.check_is_positive('window size', self.block_method[1])

    elif (self.block_method[0] == 'qgram'):
      if (len(self.block_method) != 4):
        logging.exception('Blocking method "qgram" needs three parameters: ' \
                          + '%s' % (str(self.block_method)))
        raise Exception
      auxiliary.check_is_integer('q', self.block_method[1])
      auxiliary.check_is_positive('q', self.block_method[1])
      auxiliary.check_is_flag('padded', self.block_method[2])
      auxiliary.check_is_normalised('threshold', self.block_method[3])

    else:
      logging.exception('Illegal blocking method given: %s' %
                        (str(self.block_method)))
      raise Exception

    num_rec1 = self.dataset1.num_records
    num_rec2 = self.dataset2.num_records

    if (num_rec1 < num_rec2):  # Make references to small and large data sets
      self.small_data_set_no =  0
      self.large_data_set_no =  1
      self.small_dataset = self.dataset1
      self.large_dataset = self.dataset2
      self.small_rec_cache = self.rec_cache1
      self.large_rec_cache = self.rec_cache2
      self.small_comp_field = self.comp_field_used1
      self.large_comp_field = self.comp_field_used2

    else:  # This will also hold for a deduplication
      self.small_data_set_no =  1
      self.large_data_set_no =  0
      self.small_dataset = self.dataset2
      self.large_dataset = self.dataset1
      self.small_rec_cache = self.rec_cache2
      self.large_rec_cache = self.rec_cache1
      self.small_comp_field = self.comp_field_used2
      self.large_comp_field = self.comp_field_used1

    self.log([('Small data set', self.small_dataset.description),
              ('Large data set',  self.large_dataset.description),
              ('Blocking method', self.block_method)])

    self.QGRAM_START_CHAR = chr(1)
    self.QGRAM_END_CHAR =   chr(2)

  # ---------------------------------------------------------------------------

  def build(self):
    """Method to build an index data structure.

       Read all records from small file, extract blocking variables and then
       insert records into blocks.
    """

    logging.info('')
    logging.info('Build BigMatch index: "%s"' % (self.description))

    start_time = time.time()

    num_indices = len(self.index_def)

    # Index data structure for blocks is one dictionary per index
    #
    for i in range(num_indices):
      self.index1[i] = {}  # Index for data set 1

    # Calculate a counter for the progress report
    #
    if (self.progress_report != None):
      progress_report_cnt = max(1, int(self.small_dataset.num_records / \
                                   (100.0 / self.progress_report)))
    else:  # So no progress report is being logged
      progress_report_cnt = self.small_dataset.num_records + 1

    rec_read = 0  # Number of records read from the small data set

    max_block_size = 0  # Keep size of largest block

    get_index_values_funct = self.__get_index_values__  # Shorthands
    small_rec_cache =        self.small_rec_cache
    small_data_set_no =      self.small_data_set_no
    skip_missing =           self.skip_missing
    this_index =             self.index1
    small_comp_field =       self.small_comp_field

    # Reading loop over all records in the small data set - - - - - - - - - - -
    #
    for (rec_ident, rec) in self.small_dataset.readall():

      # Extract record fields needed for comparisons (set all others to '')
      #
      comp_rec = []

      field_ind = 0

      for field in rec:
        if (field_ind in small_comp_field):
          comp_rec.append(field.lower())  # Make lowercase
        else:
          comp_rec.append('')
        field_ind += 1

      small_rec_cache[rec_ident] = comp_rec  # Put into record cache

      # Now get the index variable values for this record - - - - - - - - - - -
      #
      rec_index_val_list = get_index_values_funct(rec, small_data_set_no)

      for i in range(num_indices):  # Put record identifier into all indices

        index_val = rec_index_val_list[i]

        if ((index_val != '') or (skip_missing == False)):

          # Get and update list of record identifiers for this index value
          #
          block_val_list = this_index[i].get(index_val, [])
          block_val_list.append(rec_ident)
          max_block_size = max(max_block_size, len(block_val_list))

          # Put it back into the inverted index (dictionary)
          #
          this_index[i][index_val] = block_val_list

      rec_read += 1

      if ((rec_read % progress_report_cnt) == 0):
        self.__log_build_progress__(rec_read, self.small_dataset.num_records,
                                    start_time)

    # For sort block method a sorted list of all index values is needed - - - -
    #
    if (self.block_method[0] == 'sort'):

      self.sorted_index_val_list = []

      for i in range(num_indices):
        this_index_vals = this_index[i].keys()
        this_index_vals.sort()

        self.sorted_index_val_list.append(this_index_vals)

    # For q-gram block method, convert index values into q-gram sub-lists - - -
    # and create indices based on these sub-lists converted back to strings
    #
    elif (self.block_method[0] == 'qgram'):

      q =         self.block_method[1]  # Get parameter values for q-gram
      q1 =        q-1
      padded =    self.block_method[2]
      threshold = self.block_method[3]

      if (padded == True):  # More shorthands
        start_str = q1*self.QGRAM_START_CHAR
        end_str =   q1*self.QGRAM_END_CHAR

      # Depending upon threshold values use one of two sub-list methods
      # (the switch value of 0.75 was found experimentally)
      #
      if (threshold > 0.75):
        qgram_sublist_funct = self.__get_sublists1__
      else:
        qgram_sublist_funct = self.__get_sublists2__

      logging.info('Convert basic inverted index into %d-gram index' % (q))

      num_blocks =       0
      num_qgram_blocks = 0

      for i in range(num_indices):

        qstart_time = time.time()

        basic_index = self.index1[i]
        qgram_index = {}

        num_blocks += len(basic_index)

        for index_val in basic_index:  # Loop over all the values in this index

          # Get all records with this index value
          #
          basic_rec_ident_list = basic_index[index_val]

          if (padded == True):
            qgram_str = '%s%s%s' % (start_str, index_val, end_str)
          else:
            qgram_str = index_val

          # Create the q-gram list
          #
          qgram_list= [qgram_str[j:j+q] for j in xrange(len(qgram_str)-q1)]

          # Calculate length of sub-lists needed and then create them
          #
          num_qgrams = len(qgram_list)

          min_num_qgrams = max(1, int(num_qgrams*threshold))

          qgram_sublists = qgram_sublist_funct(qgram_list, min_num_qgrams)

          # Convert q-gram sub-lists into strings and insert into q-gram index
          #
          for qgram_sublist in qgram_sublists:

            qgram_substr = ''.join(qgram_sublist)
            qgram_rec_ident_list = qgram_index.get(qgram_substr, [])

            for rec_ident in basic_rec_ident_list:
              if (rec_ident not in qgram_rec_ident_list):
                qgram_rec_ident_list.append(rec_ident)
            qgram_index[qgram_substr] = qgram_rec_ident_list

        num_qgram_blocks += len(qgram_index)

        logging.info('  Built %d-gram index %d with %d blocks in %s' % \
               (q, i, len(qgram_index),
               auxiliary.time_string(time.time()-qstart_time)))

        this_index[i] = qgram_index

    total_num_blocks = 0
    for i in range(num_indices):
      total_num_blocks += len(self.index1[i])

    used_sec_str = auxiliary.time_string(time.time()-start_time)
    rec_time_str = auxiliary.time_string((time.time()-start_time) / \
                                           self.small_dataset.num_records)
    logging.info('Read and indexed %d records in %s (%s per record)' % \
                 (self.small_dataset.num_records, used_sec_str, rec_time_str))
    logging.info('')

    logging.info('Built BigMatch index containing %d blocks in %s' % \
                 (total_num_blocks,
                  auxiliary.time_string(time.time()-start_time)))
    logging.info('  Largest block contains %d records' % (max_block_size))

    memory_usage_str = auxiliary.get_memory_usage()
    if (memory_usage_str != None):
      logging.info('  '+memory_usage_str)

    self.status = 'built'  # Update index status

  # ---------------------------------------------------------------------------

  def compact(self):
    """Method to compact an index data structure.

       Nothing has to be done here.
    """

    logging.info('')
    logging.info('Compact BigMatch index: "%s"' % (self.description))
    logging.info('  Nothing needs to be done.')

    self.status = 'compacted'  # Update index status

  # ---------------------------------------------------------------------------

  def run(self, length_filter_perc = None, cut_off_threshold = None):
    """Iterate over all blocks in the index.

       Read the large data set and compare each record with the corresponding
       records from the small data sets as stored in the index.

       Compare the record pairs and return a weight vector dictionary with keys
       made of a tuple (record identifier 1, record identifier 2), and
       corresponding values the comparison weights.
    """

    logging.info('')
    logging.info('Read large data set and start comparison of record pairs')
    if (self.log_funct != None):
      self.log_funct('Read large data set and start comparison of record ' + \
                     'pairs')

    # Check if index has been compacted - - - - - - - - - - - - - - - - - - - -
    #
    if (self.status != 'compacted'):
      logging.exception('Index "%s" has not been compacted, running ' % \
                        (self.description)+'comparisons not possible')
      raise Exception

    # Check if weight vector file should be written - - - - - - - - - - - - - -
    #
    if (self.weight_vec_file != None):
      try:
        weight_vec_fp = open(self.weight_vec_file, 'w')
      except:
        logging.exception('Cannot write weight vector file: %s' % \
                          (self.weight_vec_file))
        raise Exception
      weight_vec_writer = csv.writer(weight_vec_fp)

      # Write header line with descriptions of field comparisons
      #
      weight_vec_header_line = ['rec_id1', 'rec_id2'] + \
                                self.__get_field_names_list__()
      weight_vec_writer.writerow(weight_vec_header_line)

    start_time = time.time()

    num_indices = len(self.index_def)

    # Check length filter and cut-off threshold arguments - - - - - - - - - - -
    #
    if (length_filter_perc != None):
      auxiliary.check_is_percentage('Length filter percentage',
                                    length_filter_perc)
      logging.info('  Length filtering set to %.1f%%' % (length_filter_perc))
      length_filter_perc /= 100.0  # Normalise

    if (cut_off_threshold != None):
      auxiliary.check_is_number('Cut-off threshold', cut_off_threshold)
      logging.info('  Cut-off threshold set to: %.2f' % (cut_off_threshold))

    num_rec_pairs_filtered =    0  # Count number of removed record pairs
    num_rec_pairs_below_thres = 0

    # For sort and q-gram block methods get their parameter values - - - - - -
    #
    block_method = self.block_method[0]

    if (block_method == 'sort'):
      win_size =  self.block_method[1]
      win_size_after =  win_size/2  # Number of index values after and before
      win_size_before = win_size - win_size_after - 1
      sorted_index_val_list = self.sorted_index_val_list

      num_sorted_index_val = []  # Get number of values for each index
      for i in range(num_indices):
        num_sorted_index_val.append(len(sorted_index_val_list[i]))

    elif (block_method == 'qgram'):
      q =         self.block_method[1]
      q1 =        q-1
      padded =    self.block_method[2]
      threshold = self.block_method[3]

      if (padded == True):  # More shorthands
        start_str = q1*self.QGRAM_START_CHAR
        end_str =   q1*self.QGRAM_END_CHAR

      # Depending upon threshold values use one of two sub-list methods
      # (the switch value of 0.75 was found experimentally)
      #
      if (threshold > 0.75):
        qgram_sublist_funct = self.__get_sublists1__
      else:
        qgram_sublist_funct = self.__get_sublists2__

    compare_funct =          self.rec_comparator.compare  # Shorthands
    get_index_values_funct = self.__get_index_values__
    skip_missing =           self.skip_missing
    large_data_set_no =      self.large_data_set_no
    this_index =             self.index1
    rec_length_cache =       self.rec_length_cache
    comp_field_used_list =   self.large_comp_field
    find_closest_funct =     self.__find_closest__
    small_rec_cache =        self.small_rec_cache
    small_data_set_no =      self.small_data_set_no

    # Calculate a counter for the progress report
    #
    if (self.progress_report != None):
      progress_report_cnt = max(1, int(self.large_dataset.num_records / \
                                   (100.0 / self.progress_report)))
    else:  # So no progress report is being logged
      progress_report_cnt = self.large_dataset.num_records + 1

    weight_vec_dict = {}  # Dictionary with calculated weight vectors

    rec_read =  0  # Number of records read from the large data set
    comp_done = 0  # Number of comparisons done

    # Set of all records from the small data set in a block
    #
    small_block_rec_set = set()

    # Reading loop over all records in the large data set - - - - - - - - - - -
    #
    for (large_rec_ident, large_rec) in self.large_dataset.readall():

      if (length_filter_perc != None):  # Get length of record in characters
                                        # (only for fields used in matching)
        large_rec_len = 0

        field_ind = 0
        for field in large_rec:
          if (field_ind in comp_field_used_list):
            large_rec_len += len(field)
          field_ind += 1

      large_rec_lower = []  # Make all values lowercase
      for field_val in large_rec:
        large_rec_lower.append(field_val.lower())
      large_rec = large_rec_lower

      # Get the index variable values for this record
      #
      rec_index_val_list = get_index_values_funct(large_rec, large_data_set_no)

      for i in range(num_indices):  # Put record identifier into all indices

        index_val = rec_index_val_list[i]

        if ((index_val != '') or (skip_missing == False)):

          # Process index value depending upon the block method
          #
          if (block_method == 'block'):  # Simply get record identifiers

            small_rec_ident_list = this_index[i].get(index_val, [])

          elif (block_method == 'sort'):  # Get record identifiers in window -

            # Find the list index of this index value or value before
            #
            list_index = find_closest_funct(sorted_index_val_list[i],index_val)

            # Get indices in list for the window and then the index values
            #
            win_start = max(0, list_index-win_size_before)
            win_end = min(num_sorted_index_val[i], list_index+win_size_after+1)

            win_index_val_list = sorted_index_val_list[i][win_start:win_end]

            small_rec_ident_list = []

            for sort_index_val in win_index_val_list:

              for small_rec_ident in this_index[i].get(sort_index_val, []):

                if (small_rec_ident not in small_rec_ident_list):
                  small_rec_ident_list.append(small_rec_ident)

          elif (block_method == 'qgram'):  # Make q-gram sub-lists - - - - - -

            if (padded == True):
              qgram_str = '%s%s%s' % (start_str, index_val, end_str)
            else:
              qgram_str = index_val

            # Create the q-gram list
            #
            qgram_list= [qgram_str[j:j+q] for j in xrange(len(qgram_str)-q1)]

            # Calculate length of sub-lists needed and then create them
            #
            num_qgrams = len(qgram_list)

            min_num_qgrams = max(1, int(num_qgrams*threshold))

            qgram_sublists = qgram_sublist_funct(qgram_list, min_num_qgrams)

            small_rec_ident_list = []

            # Convert q-gram sub-lists into strings and get record identifiers
            # from small data set
            #
            for qgram_sublist in qgram_sublists:

              qgram_substr = ''.join(qgram_sublist)
              for small_rec_ident in this_index[i].get(qgram_substr, []):

                if (small_rec_ident not in small_rec_ident_list):
                  small_rec_ident_list.append(small_rec_ident)

          else:
            logging.exception('Illegal blocking method given: %s' %
                              (str(self.block_method)))
            raise Exception

          # Now loop over all record identifiers from small data set - - - - -
          # with this index value
          #
          for small_rec_ident in small_rec_ident_list:

            # Only compare record if it hasn't been compared before
            #
            if (small_rec_ident not in small_block_rec_set):
              small_block_rec_set.add(small_rec_ident)

              # Get record values from cache
              #
              small_rec = small_rec_cache[small_rec_ident]

              do_comp = True  # Flag, specify if comparison should be done

              if (length_filter_perc != None):
                if (small_rec_ident in rec_length_cache):  # Length is cached
                  small_rec_len = rec_length_cache[small_rec_ident]
                else:
                  small_rec_len = len(''.join(small_rec))  # Length in chars
                  rec_length_cache[small_rec_ident] = small_rec_len

                perc_diff = float(abs(large_rec_len - small_rec_len)) / \
                            max(large_rec_len, small_rec_len)

                if (perc_diff > length_filter_perc):
                  do_comp = False  # Difference too large, don't do comparison
                  num_rec_pairs_filtered += 1

              if (do_comp == True):

                if (small_data_set_no == 0):
                  w_vec = compare_funct(small_rec, large_rec)

                  if (cut_off_threshold == None) or \
                     (sum(w_vec) >= cut_off_threshold):
                    if (self.weight_vec_file == None):
                      weight_vec_dict[(small_rec_ident,large_rec_ident)]= w_vec
                    else:
                      weight_vec_writer.writerow([small_rec_ident,
                                                  large_rec_ident]+w_vec)
                  else:
                    num_rec_pairs_below_thres += 1

                else:
                  w_vec = compare_funct(large_rec, small_rec)

                  if (cut_off_threshold == None) or \
                     (sum(w_vec) >= cut_off_threshold):
                    if (self.weight_vec_file == None):
                      weight_vec_dict[(large_rec_ident,small_rec_ident)]= w_vec
                    else:
                      weight_vec_writer.writerow([large_rec_ident,
                                                  small_rec_ident]+w_vec)
                  else:
                    num_rec_pairs_below_thres += 1

              comp_done += 1

      small_block_rec_set.clear()

      rec_read += 1

      if ((rec_read % progress_report_cnt) == 0):
        self.__log_build_progress__(rec_read, self.large_dataset.num_records,
                                    start_time)
        logging.info('    Number of comparisons done so far: %d (%.1f in ' % \
                     (comp_done, float(comp_done)/rec_read)+'average per ' + \
                     'record from the large data set)')

    self.num_rec_pairs = comp_done

    used_sec_str = auxiliary.time_string(time.time()-start_time)
    rec_read_time_str = auxiliary.time_string((time.time()-start_time) / \
                                              self.large_dataset.num_records)
    if (comp_done > 0):
      rec_comp_time_str = auxiliary.time_string((time.time()-start_time) / \
                                                comp_done)
    else:
      rec_comp_time_str = 0
    logging.info('Read %d records in %s (%s per record)' % \
                 (self.large_dataset.num_records, used_sec_str,
                  rec_read_time_str))
    logging.info('  Compared %d record pairs in %s (%s per pair)' % \
                 (comp_done, used_sec_str, rec_comp_time_str))
    if (length_filter_perc != None):
      logging.info('  Length filtering (set to %.1f%%) filtered %d record ' % \
                   (length_filter_perc*100, num_rec_pairs_filtered) + 'pairs')
    if (cut_off_threshold != None):
      logging.info('  %d record pairs had summed weights below threshold ' % \
                   (num_rec_pairs_below_thres) + '%.2f' % (cut_off_threshold))

    memory_usage_str = auxiliary.get_memory_usage()
    if (memory_usage_str != None):
      logging.info('  '+memory_usage_str)

    if (self.weight_vec_file == None):
      return [self.__get_field_names_list__(), weight_vec_dict]
    else:
      weight_vec_fp.close()
      return None

# =============================================================================

class DedupIndex(Indexing):
  """Class that implements an index method optimised for deduplication: The
     index is built while the data set is loaded, and the current record is
     compared to previous records according to a blocking method.

     This index can only be used for deduplications but not for linkages.

     Records that have the same index variable values for an index are put into
     the same blocks, and only records within a block are then compared.

     The additional argument (besides the base class arguments) which has to be
     set when this index is initialised is:

       block_method  Determines the blocking method used. It can be set to:
                       ('block'),
                       ('sort', window_size), or
                       ('qgram', q, padded, threshold)
                     with parameters similar to the corresponding indexing
                     methods.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'block_method' argument first, then call the
       base class constructor.

       Note that number of record pairs will not be known after initialisation
       (so it is left at value None).
    """

    self.block_method = None

    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('block_m')):
        auxiliary.check_is_tuple('block_method', value)
        self.block_method = value

      else:
        base_kwargs[keyword] = value

    Indexing.__init__(self, base_kwargs)  # Initialise base class

    if (self.do_deduplication == False):
      logging.exception('DedupIndex can only be used for deduplications, ' + \
                        'but not linkages')
      raise Exception

    # Check if block method and parameters given are OK - - - - - - - - - - - -
    #
    auxiliary.check_is_not_none('block_method', self.block_method)

    if (self.block_method[0] == 'block'):
      if (len(self.block_method) != 1):
        logging.exception('Blocking method "block" does not need parameters:' \
                          + ' %s' % (str(self.block_method)))
        raise Exception

    elif (self.block_method[0] == 'sort'):
      if (len(self.block_method) != 2):
        logging.exception('Blocking method "sort" needs one parameter: %s' % \
                          (str(self.block_method)))
        raise Exception
      auxiliary.check_is_integer('window size', self.block_method[1])
      auxiliary.check_is_positive('window size', self.block_method[1])

    elif (self.block_method[0] == 'qgram'):
      if (len(self.block_method) != 4):
        logging.exception('Blocking method "qgram" needs three parameters: ' \
                          + '%s' % (str(self.block_method)))
        raise Exception
      auxiliary.check_is_integer('q', self.block_method[1])
      auxiliary.check_is_positive('q', self.block_method[1])
      auxiliary.check_is_flag('padded', self.block_method[2])
      auxiliary.check_is_normalised('threshold', self.block_method[3])

    else:
      logging.exception('Illegal blocking method given: %s' %
                        (str(self.block_method)))
      raise Exception

    self.log([('Blocking method', self.block_method)])

    self.QGRAM_START_CHAR = chr(1)
    self.QGRAM_END_CHAR =   chr(2)

  # ---------------------------------------------------------------------------

  def build(self):
    """Method to build an index data structure.

       Nothing has to be done here.
    """

    logging.info('')
    logging.info('Build Dedup index: "%s"' % (self.description))
    logging.info('  Nothing needs to be done.')

    self.status = 'built'  # Update index status

  # ---------------------------------------------------------------------------

  def compact(self):
    """Method to compact an index data structure.

       Nothing has to be done here.
    """

    logging.info('')
    logging.info('Compact Dedup index: "%s"' % (self.description))
    logging.info('  Nothing needs to be done.')

    self.status = 'compacted'  # Update index status

  # ---------------------------------------------------------------------------

  def run(self, length_filter_perc = None, cut_off_threshold = None):
    """Iterate over all blocks in the index.

       Read the data set and compare each record with the earlier read records
       that have been inserted into the index. Finally insert record itself.

       Compare the record pairs and return a weight vector dictionary with keys
       made of a tuple (record identifier 1, record identifier 2), and
       corresponding values the comparison weights.
    """

    logging.info('')
    logging.info('Read data set and compare record pairs')
    if (self.log_funct != None):
      self.log_funct('Read data set and compare record pairs')

    # Check if index has been compacted - - - - - - - - - - - - - - - - - - - -
    #
    if (self.status != 'compacted'):
      logging.exception('Index "%s" has not been compacted, running ' % \
                        (self.description)+'comparisons not possible')
      raise Exception

    # Check if weight vector file should be written - - - - - - - - - - - - - -
    #
    if (self.weight_vec_file != None):
      try:
        weight_vec_fp = open(self.weight_vec_file, 'w')
      except:
        logging.exception('Cannot write weight vector file: %s' % \
                          (self.weight_vec_file))
        raise Exception
      weight_vec_writer = csv.writer(weight_vec_fp)

      # Write header line with descriptions of field comparisons
      #
      weight_vec_header_line = ['rec_id1', 'rec_id2'] + \
                                self.__get_field_names_list__()
      weight_vec_writer.writerow(weight_vec_header_line)

    start_time = time.time()

    num_indices = len(self.index_def)

    index = {}  # Index data structure for blocks is one dictionary per index
    for i in range(num_indices):
      index[i] = {}  # Index for data set 1

    # Check length filter and cut-off threshold arguments - - - - - - - - - - -
    #
    if (length_filter_perc != None):
      auxiliary.check_is_percentage('Length filter percentage',
                                    length_filter_perc)
      logging.info('  Length filtering set to %.1f%%' % (length_filter_perc))
      length_filter_perc /= 100.0  # Normalise

    if (cut_off_threshold != None):
      auxiliary.check_is_number('Cut-off threshold', cut_off_threshold)
      logging.info('  Cut-off threshold set to: %.2f' % (cut_off_threshold))

    num_rec_pairs_filtered =    0  # Count number of removed record pairs
    num_rec_pairs_below_thres = 0

    # For sort and q-gram block methods get their parameter values - - - - - -
    #
    block_method = self.block_method[0]

    if (block_method == 'sort'):
      win_size =  self.block_method[1]
      win_size_after =  win_size/2  # Number of index values after and before
      win_size_before = win_size - win_size_after - 1

      sorted_index_val_list = {}  # have a list per index of the sorted values
      for i in range(num_indices):
        sorted_index_val_list[i] = []

    elif (block_method == 'qgram'):
      q =         self.block_method[1]
      q1 =        q-1
      padded =    self.block_method[2]
      threshold = self.block_method[3]

      if (padded == True):  # More shorthands
        start_str = q1*self.QGRAM_START_CHAR
        end_str =   q1*self.QGRAM_END_CHAR

      # Depending upon threshold values use one of two sub-list methods
      # (the switch value of 0.75 was found experimentally)
      #
      if (threshold > 0.75):
        qgram_sublist_funct = self.__get_sublists1__
      else:
        qgram_sublist_funct = self.__get_sublists2__

    compare_funct =          self.rec_comparator.compare  # Shorthands
    find_closest_funct =     self.__find_closest__
    get_index_values_funct = self.__get_index_values__
    skip_missing =           self.skip_missing
    rec_cache =              self.rec_cache1
    rec_length_cache =       self.rec_length_cache
    comp_field_used_list =   self.comp_field_used1

    # Calculate a counter for the progress report
    #
    if (self.progress_report != None):
      progress_report_cnt = max(1, int(self.dataset1.num_records / \
                                   (100.0 / self.progress_report)))
    else:  # So no progress report is being logged
      progress_report_cnt = self.dataset1.num_records + 1

    weight_vec_dict = {}  # Dictionary with calculated weight vectors

    rec_read =  0  # Number of records read from the data set
    comp_done = 0  # Number of comparisons done

    # Reading loop over all records in the data set - - - - - - - - - - - - - -
    #
    for (rec_ident1, rec1) in self.dataset1.readall():

      if (rec_ident1 in rec_cache):
        logging.warn('Record with identifier "%s" appears more than once ' % \
                     (rec_ident1) + 'in data set')

      # Extract record fields needed for comparisons (set all others to '')
      # (also count length of field values for length filtering)
      #
      comp_rec =  []
      rec_len1 =  0

      for field_ind in range(len(rec1)):
        rec1[field_ind] = rec1[field_ind].lower()
        field = rec1[field_ind]

        if (field_ind in comp_field_used_list):
          comp_rec.append(field)  # Make lowecase
          rec_len1 += len(field)
        else:
          comp_rec.append('')

      rec_cache[rec_ident1] = comp_rec  # Put into record cache

      if (length_filter_perc != None):  # Cache length for length filtering
        rec_length_cache[rec_ident1] = rec_len1

      # Get the index variable values for this record
      #
      rec_index_val_list = get_index_values_funct(rec1, 0)

      # List of record identifiers for this record over all indices
      #
      this_rec_block_rec_list = []

      for i in range(num_indices):  # Put record identifier into all indices

        index_val = rec_index_val_list[i]

        if ((index_val != '') or (skip_missing == False)):

          # Get records in the same block as current records out of indices,
          # then insert record into index (depending upon blocking method)
          #
          if (block_method == 'block'):  # Simply get record identifiers

            block_rec_ident_list = index[i].get(index_val, [])

            for rec_ident2 in block_rec_ident_list:
              if (rec_ident2 not in this_rec_block_rec_list):
                this_rec_block_rec_list.append(rec_ident2)

            block_rec_ident_list.append(rec_ident1)  # Append this record

            # Put it back into the inverted index (dictionary)
            #
            index[i][index_val] = block_rec_ident_list

          elif (block_method == 'sort'):  # Get record identifiers in window -

            # First find the list index of this index value or the value before
            #
            list_index = find_closest_funct(sorted_index_val_list[i],index_val)

            # Get indices in list for the window and then the index values
            #
            win_start = max(0, list_index-win_size_before)
            win_end =  min(len(sorted_index_val_list[i]),
                               list_index+win_size_after+1)

            win_index_val_list = sorted_index_val_list[i][win_start:win_end]

            # Get all record identifiers from index for the values in window
            #
            for sort_index_val in win_index_val_list:

              for rec_ident2 in index[i][sort_index_val]:
                if (rec_ident2 not in this_rec_block_rec_list):
                  this_rec_block_rec_list.append(rec_ident2)

            # Get block for this record from index and put it into index
            #
            if (index_val in index[i]):
              block_rec_ident_list = index[i][index_val]
              block_rec_ident_list.append(rec_ident1)  # Append this record

              # Put it back into the inverted index (dictionary)
              #
              index[i][index_val] = block_rec_ident_list

            else:  # A new value, insert into index and sorted list of values
              index[i][index_val] = [rec_ident1]
              sorted_index_val_list[i].append(index_val)
              sorted_index_val_list[i].sort()

          elif (block_method == 'qgram'):  # Make q-gram sub-lists - - - - - -

            if (padded == True):
              qgram_str = '%s%s%s' % (start_str, index_val, end_str)
            else:
              qgram_str = index_val

            # Create the q-gram list
            #
            qgram_list= [qgram_str[j:j+q] for j in xrange(len(qgram_str)-q1)]

            # Calculate length of sub-lists needed and then create them
            #
            num_qgrams = len(qgram_list)

            min_num_qgrams = max(1, int(num_qgrams*threshold))

            qgram_sublists = qgram_sublist_funct(qgram_list, min_num_qgrams)

            # Convert q-gram sub-lists into strings and get record identifiers
            # from index for q-gram values
            #
            for qgram_sublist in qgram_sublists:

              qgram_substr = ''.join(qgram_sublist)  # Make it a string

              for rec_ident2 in index[i].get(qgram_substr, []):
                if (rec_ident2 not in this_rec_block_rec_list):
                  this_rec_block_rec_list.append(rec_ident2)

              # Add this record to block for this q-gram and put it into index
              #
              qgram_rec_ident_list = index[i].get(qgram_substr, [])
              qgram_rec_ident_list.append(rec_ident1)

              # Put it back into the inverted index (dictionary)
              #
              index[i][qgram_substr] = qgram_rec_ident_list

          else:
            logging.exception('Illegal blocking method given: %s' %
                              (str(self.block_method)))
            raise Exception

      # Now compare current record with all records in blocking set - - - - - -
      #
      for rec_ident2 in this_rec_block_rec_list:

        rec2 = rec_cache[rec_ident2]  # Get record values from cache

        do_comp = True  # Flag, specify if comparison should be done

        if (length_filter_perc != None):
          rec_len2 = rec_length_cache[rec_ident2]

          perc_diff = float(abs(rec_len1 - rec_len2)) / max(rec_len1, rec_len2)

          if (perc_diff > length_filter_perc):
            do_comp = False  # Difference too large, don't do comparison
            num_rec_pairs_filtered += 1

        if (do_comp == True):

          w_vec = compare_funct(rec1, rec2)

          if (cut_off_threshold == None) or (sum(w_vec) >= cut_off_threshold):

            # Make sure record identifiers are sorted
            #
            if (rec_ident1 < rec_ident2):
              if (self.weight_vec_file == None):
                weight_vec_dict[(rec_ident1, rec_ident2)] = w_vec
              else:
                weight_vec_writer.writerow([rec_ident1, rec_ident2]+w_vec)

            else:
              if (self.weight_vec_file == None):
                weight_vec_dict[(rec_ident2, rec_ident1)] = w_vec
              else:
                weight_vec_writer.writerow([rec_ident1, rec_ident2]+w_vec)
          else:
            num_rec_pairs_below_thres += 1

          comp_done += 1

      del this_rec_block_rec_list

      rec_read += 1

      if ((rec_read % progress_report_cnt) == 0):
        self.__log_build_progress__(rec_read, self.dataset1.num_records,
                                    start_time)
        logging.info('    Number of comparisons done so far: %d (%.1f in ' % \
                     (comp_done, float(comp_done)/rec_read)+'average per ' + \
                     'record)')

    self.num_rec_pairs = comp_done

    rec_cache.clear()
    rec_length_cache.clear()

    used_sec_str = auxiliary.time_string(time.time()-start_time)
    rec_read_time_str = auxiliary.time_string((time.time()-start_time) / \
                                              self.dataset1.num_records)
    if (comp_done > 0):
      rec_comp_time_str = auxiliary.time_string((time.time()-start_time) / \
                                                comp_done)
    else:
      rec_comp_time_str = 0
    logging.info('Read %d records in %s (%s per record)' % \
                 (self.dataset1.num_records, used_sec_str, rec_read_time_str))
    logging.info('  Compared %d record pairs in %s (%s per pair)' % \
                 (comp_done, used_sec_str, rec_comp_time_str))
    if (length_filter_perc != None):
      logging.info('  Length filtering (set to %.1f%%) filtered %d record ' % \
                   (length_filter_perc*100, num_rec_pairs_filtered) + 'pairs')
    if (cut_off_threshold != None):
      logging.info('  %d record pairs had summed weights below threshold ' % \
                   (num_rec_pairs_below_thres) + '%.2f' % (cut_off_threshold))

    memory_usage_str = auxiliary.get_memory_usage()
    if (memory_usage_str != None):
      logging.info('  '+memory_usage_str)

    if (self.weight_vec_file == None):
      return [self.__get_field_names_list__(), weight_vec_dict]
    else:
      weight_vec_fp.close()
      return None

# =============================================================================
