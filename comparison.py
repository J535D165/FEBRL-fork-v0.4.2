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
# The Original Software is: "comparison.py"
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

"""Module with classes for field (attribute) and record comparisons.

   This module provides classes for record and field comparisons that can be
   used for the linkage process.

   TODO:
   - improve cache: how to remove oldest count-1 pair efficiently?
   - do caching timing test -> comparisonTiming.py module
   - improve value frequency based weight calculations

   TODO (old version):
   - What do we do if we want to link two data sets, and we know that we have
     different distributions for e.g. surnames, but only a lookup table for the
     surnames from one data set?
   - Add frequency table capabilities for numeric, data and age field
     comparators.
   - Do a proper date/time comparison function (maybe using mxDateTime module)
     that takes both date and time into consideration.
"""

# =============================================================================
# Import necessary modules (Python standard modules first, then Febrl modules)

import bz2
import datetime
import difflib
import logging
import math
import time
import zlib

import auxiliary
import encode
import mymath

# =============================================================================

class RecordComparator:
  """Class that implements a record comparator to compare two records and
     compute (and return) a weight vector.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, dataset1, dataset2, field_comparator_list, descr = ''):
    """Constructor.

       Has as arguments two data set objects and a list of field comparators,
       that each are made of a tuple:

         (field comparator, field name in dataset1, field name in dataset2)

       The constructor checks if the field names are available in the two data
       sets and then generates an efficient list of comparison methods and
       field columns (into the data sets).
    """

    auxiliary.check_is_string('description', descr)
    self.description = descr

    # Check if the input objects needed are lists
    #
    auxiliary.check_is_list('dataset1.field_list', dataset1.field_list)
    auxiliary.check_is_list('dataset2.field_list', dataset2.field_list)
    auxiliary.check_is_list('field_comparator_list', field_comparator_list)

    self.dataset1 = dataset1
    self.dataset2 = dataset2

    self.field_comparator_list = field_comparator_list
    self.field_comparison_list = []  # Only compare methods and field columns

    # Extract field names from the two data set field name lists
    #
    dataset1_field_names = []
    dataset2_field_names = []

    for (field_name, field_data) in dataset1.field_list:
      dataset1_field_names.append(field_name)

    for (field_name, field_data) in dataset2.field_list:
      dataset2_field_names.append(field_name)

    # Go through field list and check if available
    #
    for (field_comp, field_name1, field_name2) in self.field_comparator_list:

      if (field_name1 not in  dataset1_field_names):
        logging.exception('Field "%s" is not in data set 1 field name list: ' \
                          % (field_name1) + '%s' % (str(dataset1.field_list)))
        raise Exception
      field_index1 = dataset1_field_names.index(field_name1)

      if (field_name2 not in dataset2_field_names):
        logging.exception('Field "%s" is not in data set 2 field name list: ' \
                          % (field_name2) + '%s' % (str(dataset2.field_list)))
        raise Exception
      field_index2 = dataset2_field_names.index(field_name2)

      field_tuple = (field_comp.compare, field_index1, field_index2)

      self.field_comparison_list.append(field_tuple)

    assert len(self.field_comparison_list) == len(self.field_comparator_list)

    # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.info('')
    logging.info('Initialised record comparator')
    logging.info('  Data set 1:    %s' % (dataset1.description))
    logging.info('    Field names: %s' % (str(dataset1_field_names)))
    logging.info('  Data set 2:    %s' % (dataset2.description))
    logging.info('    Field names: %s' % (str(dataset2_field_names)))

    logging.info('  Field comparators:')

    for i in range(len(self.field_comparison_list)):
      (field_comp, field_name1, field_name2) = self.field_comparator_list[i]
      (comp_method, field_index1, field_index2) = self.field_comparison_list[i]

      logging.info('    Description:  %s' % (field_comp.description))
      logging.info('    Field name 1: %s (at column %d)' % \
                   (field_name1, field_index1))
      logging.info('    Field name 2: %s (at column %d)' % \
                   (field_name2, field_index2))

  # ---------------------------------------------------------------------------

  def compare(self, rec1, rec2):
    """Compare two records (list of fields) and return a vector with weight
       values (floating-point numbers)
    """

    weight_vector = []

    # Compute a weight for each field comparator
    #
    for (comp_method,field_index1,field_index2) in self.field_comparison_list:

      if (field_index1 >= len(rec1)):
        val1 = ''
      else:
        val1 = rec1[field_index1]

      if (field_index2 >= len(rec2)):
        val2 = ''
      else:
        val2 = rec2[field_index2]

      val1 = val1.lower()
      val2 = val2.lower()

      w = comp_method(val1,val2)
      weight_vector.append(w)

    return weight_vector

  # ---------------------------------------------------------------------------

  def get_cache_stats(self):
    """Extract information about the cache size, maximum and average counts for
       all the field comparators that have an activated cache.
    """

    logging.info('Caching statistics for record comparator "%s"' % \
                 (self.description))

    # Go through field list and check if available
    #
    for (field_comp, field_name1, field_name2) in self.field_comparator_list:
      logging.info('  Field comparator: "%s"' % (field_comp.description))

      if (field_comp.do_caching == False):
        logging.info('    Caching is not activated')

      elif (field_comp.cache == {}):
        logging.info('    Caching is activated but cache is empty')

      else:
        cache_size = len(field_comp.cache)

        cache_max_count = -1
        cache_count_sum = 0.0

        for cache_key in field_comp.cache:  # Get all cache entries and counts
          (cache_weight, access_count) = field_comp.cache[cache_key]

          cache_max_count = max(cache_max_count, access_count)
          cache_count_sum += access_count

        cache_avrg_count = cache_count_sum / cache_size

        cache_size_str = '    Number of cache entries: %s' % (cache_size)
        if (field_comp.max_cache_size == None):
          cache_size_str += ' (Cache size is unlimited)'
        else:
          cache_size_str += ' (Cache size is limited to %d entries)' % \
                            (field_comp.max_cache_size)
        logging.info(cache_size_str)
        logging.info('    Maximum and average cache entry count: %d / %.2f' % \
                     (cache_max_count, cache_avrg_count))

# =============================================================================

class FieldComparator:
  """Base class for field comparators.

     All field comparators have the following instance variables, which can be
     set when a field comparator is initialised:

       description      A string describing the field comparator.
       do_caching       A flag, True or False, to enable or disable caching.
       max_cache_size   The maximum number of comparisons to be cached.
       cache            A dictionary with cached comparisons.
       missing_values   A list of one or more strings that correspond to
                        missing values.
       missing_weight   Numerical weight when values are missing.
       agree_weight     Numerical weight when two values agree.
       disagree_weight  Numerical weight when two values totally disagree.
       val_freq_table   A dictionary with values (as keys) and their counts
                        (as values). If provided, the comparison weight will be
                        frequency adjusted. Default is None (not provided).
  """

  # ---------------------------------------------------------------------------

  def __init__(self, base_kwargs):
    """Constructor.
    """

    # General attributes for all field comparators
    #
    self.description =     ''

    self.do_caching =           False        # Caching disabled by default
    self.max_cache_size =       None         # None - no maximum cache size
    self.cache =                {}           # Dictionary with cached values
    self.cache_num_not_cached = 0            # Number of comparisons not cached
    self.cache_warn_counts =    [2,5,10,50]  # List of when warnings should be
                                             # given (minimum counts)
    self.missing_values =  ['']
    self.missing_weight =  0.0
    self.agree_weight =    1.0
    self.disagree_weight = 0.0


    self.val_freq_table = None
    self.val_freq_sum =   None   # If a frequency table is provided, the sum of
                                 # all counts will be calculated and stored
    self.freq_max_weight = None  # A maximum weight value for frequency based
                                 # agreement values. If not provided it will be
                                 # set to the general agreement value

    # Process base keyword arguments (all data set specific keywords were
    # processed in the derived class constructor)
    #
    for (keyword, value) in base_kwargs.items():

      if (keyword.startswith('desc')):
        auxiliary.check_is_string('description', value)
        self.description = value

      elif (keyword.startswith('do_c')):
        auxiliary.check_is_flag('do_caching', value)
        self.do_caching = value

      elif (keyword.startswith('max_cach')):
        if (value == None):  # No maximum cache size (use with care!)
          self.max_cache_size = value
        else:
          auxiliary.check_is_integer('max_cache_size', value)
          auxiliary.check_is_not_negative('max_cache_size', value)
          self.max_cache_size = value

      elif (keyword.startswith('missing_v')):
        auxiliary.check_is_list('missing_values', value)
        self.missing_values = value

      elif (keyword.startswith('missing_w')):
        auxiliary.check_is_number('missing_weight', value)
        self.missing_weight = value

      elif (keyword.startswith('agr')):
        auxiliary.check_is_number('agree_weight', value)
        self.agree_weight = value

      elif (keyword.startswith('disa')):
        auxiliary.check_is_number('disagree_weight', value)
        self.disagree_weight = value

      elif (keyword.startswith('val_fr')):
        auxiliary.check_is_dictionary('val_freq_table', value)
        self.val_freq_table = value

      elif (keyword.startswith('freq_m')):
        auxiliary.check_is_number('freq_max_weight', value)
        self.freq_max_weight = value

      else:
        logging.exception('Illegal constructor argument keyword: %s' % \
                          (str(keyword)))
        raise Exception

    # Create a dictionary with the warning counts - - - - - - - - - - - - - - -
    #
    self.cache_warn_dict_counts = {}

    for i in self.cache_warn_counts:
      self.cache_warn_dict_counts[i] = 0  # No value pairs with count i so far

    # If a frequency table is given calculate the sum of all counts
    #
    if (self.val_freq_table != None):
      self.val_freq_sum = 0

      for (k,v) in self.val_freq_table.items():
        auxiliary.check_is_integer('frequency table entry "%s"' % (k), v)
        auxiliary.check_is_positive('frequency table entry "%s"' % (k), v)
        self.val_freq_sum += v

      # Make sure a maximum frequency agreement weight is set
      #
      if (self.freq_max_weight == None):  # Has not be given as argument
        self.freq_max_weight = self.agree_weight

        logging.warning('Setting maximum frequency agreement weight to ' + \
                        'general agreement weight')

    self.__check_weights__()

  # ---------------------------------------------------------------------------

  def __check_weights__(self):
    """Check the values of weights. Should not be used from outside the module.

       Make sure the agreement weight is larger than the disagreement weight,
       and that the missing weight is somewhere in between.
    """

    if (self.agree_weight < self.disagree_weight):
      logging.exception('Agreement weight set to a value smaller than the ' + \
                        'disagreement weight')
      raise Exception

    if (self.agree_weight < self.missing_weight):
      logging.exception('Agreement weight set to a value smaller than the ' + \
                        'missing weight')
      raise Exception

    if (self.disagree_weight > self.missing_weight):
      logging.exception('Disagreement weight set to a value larger than ' + \
                        'the missing weight')
      raise Exception

    if (self.val_freq_table != None):
      auxiliary.check_is_number('freq_max_weight', self.freq_max_weight)

      if (self.freq_max_weight < self.disagree_weight):
        logging.exception('Maximum frequency agreement weight set to a ' + \
                          'value smaller than the disagreement weight')
        raise Exception

      if (self.freq_max_weight < self.missing_weight):
        logging.exception('Maximum frequency agreement weight set to a ' + \
                          'value smaller than the missing weight')
        raise Exception

  # ---------------------------------------------------------------------------

  def __get_from_cache__(self, val1, val2):
    """Check if the given pair of values is in the cache, if so return cached
       similarity weight. Otherwise return None.

       If found in the cache, the pair's count is increased by one.

       warning messages are logged if the cache size is limited and all entries
       have certain counts (see numbers: self.cache_warn_counts).

       If caching is disabled, returned None.
    """

    if (self.do_caching == False):
      return None

    # Comparisons have to be symmetric: Only one of the pairs (val1,val2) and
    # (val2,val1) should be stored in the cache, so sort them
    #
    if (val1 < val2):
      cache_key = (val1, val2)
    else:
      cache_key = (val2, val1)

    if cache_key not in self.cache:  # The values pair is not in the cache
      return None

    # Get weight and access count from cache, increase count, and put it back
    #
    (cache_weight, access_count) = self.cache[cache_key]
    access_count += 1
    self.cache[cache_key] = (cache_weight, access_count)

    if (self.max_cache_size == None):
      return cache_weight  # Unlimited cache size, simply return

    if (access_count in self.cache_warn_counts):  # Check if warning needed

      num_count_pairs = self.cache_warn_dict_counts[access_count]
      self.cache_warn_dict_counts[access_count] = num_count_pairs+1

      if (num_count_pairs == self.max_cache_size):  # All pairs have this count
        logging.warning('All cache entries have a count of %d' % \
                        (access_count))

    return cache_weight

  # ---------------------------------------------------------------------------

  def __put_into_cache__(self, val1, val2, weight):
    """If caching is enabled and there is room in the cache put the given pair
       of values into the cache with the given similarity weight.

       If the cache is full don't insert values pair but increase the number of
       non-cached comparisons.

       If caching is disabled do nothing.
    """

    if (self.do_caching == False):
      return

    # Check if the cache is full
    #
    if ((self.max_cache_size != None) and \
        (len(self.cache) == self.max_cache_size)):

      self.cache_num_not_cached += 1  # One more pair not cached

      if (self.cache_num_not_cached in [100, 1000, 10000, 100000]):
        logging.warning('Cache is full, %d comparisons cannot be cached' % \
                        (self.cache_num_not_cached))

      return

    # Comparisons have to be symmetric: Only one of the pairs (val1,val2) and
    # (val2,val1) should be stored in the cache, so sort them
    #
    if (val1 < val2):
      cache_key = (val1, val2)
    else:
      cache_key = (val2, val1)

    # Insert new pair into cache and list of pairs with count 1
    #
    self.cache[cache_key] = (weight, 1)

  # ---------------------------------------------------------------------------

  def __calc_freq_agree_weight__(self, val):
    """Check if a frequency table is given and if so if the given value is in
       there - in which case a frequency based agreement weight is calculated
       and returned. Otherwise the general agreement weight is returned.
    """

    if ((self.val_freq_table != None) and (val in self.val_freq_table)):

      val_count = self.val_freq_table[val]
      val_freq = float(val_count) / float(self.val_freq_sum)

      # Agreement weight, computed according to L. Gill (2001), page 66
      #
      freq_weight = math.log(1.0 / val_freq, 2)  # log_2

      # Make sure frequency weight is not larger than maximum frequency
      # agreement weight
      #
      return min(freq_weight, self.freq_max_weight)

    else:
      return self.agree_weight

  # ---------------------------------------------------------------------------

  def __calc_freq_weights__(self, val1, val2):
    """Check if a frequency table is given and if so if the given values are in
       there - in which case frequency based agreement weights are calculated.
       The minimum frequency based weight is then returned. Otherwise the
       general agreement weight is returned.
    """

    if (self.val_freq_table == None):  # No frequency table given
      return self.agree_weight

    if (val1 in self.val_freq_table):

      val1_count = self.val_freq_table[val1]
      val1_freq = float(val1_count) / float(self.val_freq_sum)

      # Agreement weight, computed according to L. Gill (2001), page 66
      #
      freq1_weight = math.log(1.0 / val1_freq, 2)  # log_2

      # Make sure frequency weight is not larger than maximum frequency
      # agreement weight
      #
      freq1_weight = min(freq1_weight, self.freq_max_weight)

    else:
      freq1_weight =  self.agree_weight

    if (val2 in self.val_freq_table):
      val2_count = self.val_freq_table[val2]
      val2_freq = float(val2_count) / float(self.val_freq_sum)

      freq2_weight = math.log(1.0 / val2_freq, 2)  # log_2
      freq2_weight = min(freq2_weight, self.freq_max_weight)

    else:
      freq2_weight = self.agree_weight

    return min(freq1_weight, freq2_weight)

  # ---------------------------------------------------------------------------

  def set_weights(self, **kwargs):
    """Provide new values for the four possible weights.
    """

    # Process base keyword arguments
    #
    for (keyword, value) in kwargs.items():

      if (keyword.startswith('missing_w')):
        auxiliary.check_is_number('missing_weight', value)
        self.missing_weight = value

      elif (keyword.startswith('agr')):
        auxiliary.check_is_number('agree_weight', value)
        self.agree_weight = value

      elif (keyword.startswith('disa')):
        auxiliary.check_is_number('disagree_weight', value)
        self.disagree_weight = value

      elif (keyword.startswith('freq_m')):
        if (self.val_freq_table == None):
          logging.warning('No frequency table given, so maximum frequency ' + \
                          'agreement weight will not be used.')
        else:
          auxiliary.check_is_number('freq_max_weight', value)
          self.freq_max_weight = value

      else:
        logging.exception('Illegal constructor argument keyword: %s' % \
                          (str(keyword)))
        raise Exception

    self.__check_weights__()

  # ---------------------------------------------------------------------------

  def train(self):
    """Method which allows training of a field comparator before using it.

       Most field comparators do not need to be trained, see implementations in
       derived classes for details.
    """

    logging.info('No training needed for this comparator')

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two fields values, compute and return a numerical weight. See
       implementations in derived classes for details.
    """

    logging.exception('Override abstract method in derived class')
    raise Exception

  # ---------------------------------------------------------------------------

  def log(self, instance_var_list = None):
    """Write a log message with the basic field comparator instance variables
       plus the instance variable provided in the given input list (assumed to
       contain pairs of names (strings) and values).
    """

    logging.info('')
    logging.info('Field comparator: "%s"' % (self.description))
    logging.info('  Do caching:          %s' % (str(self.do_caching)))
    if (self.max_cache_size != None):
      logging.info('  Maximum cache size:  %s' % (str(self.max_cache_size)))
    else:
      logging.info('  Unlimited cache size')
    logging.info('    Warnings will be given once all cache entries have ' + \
                 'counts: %s' % (str(self.cache_warn_counts)))
    logging.info('  Missing values:      %s' % (str(self.missing_weight)))
    logging.info('  Missing weight:      %f' % (self.missing_weight))
    logging.info('  Agreement weight:    %f' % (self.agree_weight))
    logging.info('  Disagreement weight: %f' % (self.disagree_weight))

    if (self.val_freq_table != None):
      logging.info('  Number of entires in frequency table: %d' % \
                   (len(self.val_freq_table)))
      logging.info('  Frequency table sum:                  %f' % \
                   (self.val_freq_sum))
      logging.info('  Maximum frequency agreement weight:   %f' % \
                   (self.freq_max_weight))

    if (instance_var_list != None):
      logging.info('  Comparator specific variables:')

      max_name_len = 0
      for (name, value) in instance_var_list:
        max_name_len = max(max_name_len, len(name))

      for (name, value) in instance_var_list:
        pad_spaces = (max_name_len-len(name))*' '
        logging.info('    %s %s' % (name+':'+pad_spaces, str(value)))

  # ---------------------------------------------------------------------------

  def get_cache_stats(self):
    """Extract information about the cache size, maximum and average counts.
    """

    logging.info('Field comparator: "%s"' % (self.description))

    if (self.do_caching == False):
      logging.info('  Caching is not activated')

    elif (self.cache == {}):
      logging.info('  Caching is activated but cache is empty')

    else:
      cache_size = len(self.cache)

      cache_max_count = -1
      cache_count_sum = 0.0

      for cache_key in self.cache:  # Get all cache entries and extract counts

        (cache_weight, access_count) = self.cache[cache_key]

        cache_max_count = max(cache_max_count, access_count)
        cache_count_sum += access_count

      cache_avrg_count = cache_count_sum / cache_size

      cache_size_str = '  Number of cache entries: %s' % (cache_size)
      if (self.max_cache_size == None):
        cache_size_str += ' (Cache size is unlimited)'
      else:
        cache_size_str += ' (Cache size is limited to %d entries)' % \
                          (self.max_cache_size)
      logging.info(cache_size_str)
      logging.info('  Maximum and average cache entry count: %d / %.2f' % \
                   (cache_max_count, cache_avrg_count))

# =============================================================================

class FieldComparatorExactString(FieldComparator):
  """A field comparator based on exact string comparison.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """No additional attributes needed, just call base class constructor.
    """

    FieldComparator.__init__(self, kwargs)

    self.log()  # Log a message

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values using exact string comparator.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    elif (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    else:
      return self.disagree_weight

# =============================================================================

class FieldComparatorContainsString(FieldComparator):
  """A field comparator that checks if the shorter of two strings is contained
     in the longer string.

     Returns the agreement weight if the shorter string is contained in the
     longer string and the dis-agreement weight otherwise.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """No additional attributes needed, just call base class constructor.
    """

    FieldComparator.__init__(self, kwargs)

    self.log()  # Log a message

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values checking if the shorter value is contained in
       the longer value (both assumed to be strings).
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    len1 = len(val1)
    len2 = len(val2)

    if (len1 < len2):
      is_contained = (val1 in val2)
    else:
      is_contained = (val2 in val1)

    if (is_contained == True):
      return self.__calc_freq_agree_weight__(val1)

    else:
      return self.disagree_weight

# =============================================================================

class FieldComparatorTruncateString(FieldComparator):
  """A field comparator based on exact string comparison with the limitation in
     number of characters that are compared (strings are truncated).

     The additional argument (besides the base class arguments) which has to be
     set when this field comparator is initialised is:

       num_char_compared  Positive integer that gives the number of characters
                          that are compared.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'num_char_compared' argument first, then call
       the base class constructor.
    """

    self.num_char_compared = None

    # Process all keyword arguments
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('num_char')):
        auxiliary.check_is_integer('num_char_compared', value)
        auxiliary.check_is_not_negative('num_char_compared', value)
        self.num_char_compared = value

      else:
        base_kwargs[keyword] = value

    FieldComparator.__init__(self, base_kwargs)  # Process base arguments

    # Make sure 'num_char_compared' attribute is set  - - - - - - - - - - - - -
    #
    auxiliary.check_is_integer('num_char_compared', self.num_char_compared)
    auxiliary.check_is_not_negative('num_char_compared', \
                                     self.num_char_compared)

    self.log([('Maximum number of characters compared',
               self.num_char_compared)])  # Log a message

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values using exact string comparator with truncated
       strings.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    elif (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    str1 = str(val1)   # Make sure values are strings
    str2 = str(val2)

    if (str1[:self.num_char_compared] == str2[:self.num_char_compared]):
      return self.__calc_freq_weights__(val1, val2)

    else:
      return self.disagree_weight

# =============================================================================

class FieldComparatorKeyDiff(FieldComparator):
  """A field comparator that compares the fields character-wise with a maximum
     number of errors (different characters) being tolerated. This comparator
     can be used to compare numerical fields, such as date of birth, telephone
     numbers, etc.

     The additional argument (besides the base class arguments) which has to be
     set when this field comparator is initialised is:

       max_key_diff  Positive integer that gives the maximum number of
                     different characters tolerated.

     For example, if 'max_key_diff' is set to X and the number of different
     characters counted between two strings is Y (with 0 < Y <= X), then the
     resulting partial agreement weight will be computed as:

       weight = agree_weight - (Y/(X+1))*(agree_weight+abs(disagree_weight))

     If the number of different characters counted is larger than X, the
     disagreement weight will be returned.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'max_key_diff' argument first, then call the
       base class constructor.
    """

    self.max_key_diff =    None

    # Process all keyword arguments
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('max_k')):
        auxiliary.check_is_integer('max_key_diff', value)
        auxiliary.check_is_not_negative('max_key_diff', value)
        self.max_key_diff = value

      else:
        base_kwargs[keyword] = value

    FieldComparator.__init__(self, base_kwargs)  # Process base arguments

    # Make sure 'max_key_diff' attribute is set - - - - - - - - - - - - - - - -
    #
    auxiliary.check_is_integer('max_key_diff', self.max_key_diff)
    auxiliary.check_is_not_negative('max_key_diff', self.max_key_diff)

    self.log([('Maximum key difference', self.max_key_diff)])  # Log a message

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values using the key difference field comparator.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    elif (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # Calculate the key difference - - - - - - - - - - - - - - - - - - - - - -
    #
    str1 = str(val1)   # Make sure values are strings
    str2 = str(val2)
    len1 = len(str1)
    len2 = len(str2)

    # The initial number of errors is the difference in the string lengths
    #
    num_err = abs(len1 - len2)

    if (num_err > self.max_key_diff):
      return self.disagree_weight  # Too many different characters

    check_len = min(len1, len2)

    for i in range(check_len):  # Loop over positions in strings
      if (str1[i] != str2[i]):
        num_err += 1

    if (num_err > self.max_key_diff):
      return self.disagree_weight  # Too many different characters

    # Get general or frequency based agreement weight
    #
    agree_weight = self.__calc_freq_weights__(val1, val2)

    # Calculate partial agreement weight
    #
    return agree_weight - (float(num_err)/(self.max_key_diff+1.0)) * \
           (agree_weight + abs(self.disagree_weight))

# =============================================================================

class FieldComparatorNumericPerc(FieldComparator):
  """A field comparator for numeric fields, where a given percentage difference
     can be tolerated. The agreement weight is returned if the numbers are the
     same, and the disagreement weight if the percentage difference is larger
     than a maximum tolerated percentage value.

     The additional argument (besides the base class arguments) which has to be
     set when this field comparator is initialised is:

       max_perc_diff  A floating-point number between 0.0 and 100.0 giving the
                      maximum percentage difference tolerated.

     If the percentage difference is smaller than 'max_perc_diff' the resulting
     partial agreement weight is calculated according to the following formula:

       weight = agree_weight - (perc_diff / (max_perc_diff + 1.0)) *
                                (agree_weight+abs(disagree_weight))

     where the percentage difference is calculated as:

       perc_diff = 100.0 *
                   abs(value_1 - value_2) / max(abs(value_1), abs(value_2))
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'max_perc_diff' argument first, then call the
       base class constructor.
    """

    self.max_perc_diff =   None

    # Process all keyword arguments - - - - - - - - - - - - - - - - - - - - - -
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('max_p')):
        auxiliary.check_is_percentage('max_perc_diff', value)
        self.max_perc_diff = float(value)

      else:
        base_kwargs[keyword] = value

    FieldComparator.__init__(self, base_kwargs)  # Process base arguments

    # Make sure 'max_perc_diff' attribute is set - - - - - - - - - - - - - - -
    #
    auxiliary.check_is_percentage('max_perc_diff', self.max_perc_diff)

    self.log([('Maximum percentage difference', self.max_perc_diff)])

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two numerical field values and tolerate a percentage difference.

       If at least one of the two values to be compared is not a number, then
       the disagreement value is returned.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    elif (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # Calculate the percentage difference - - - - - - - - - - - - - - - - - - -
    #
    try:  # Check if the field values are numbers, if not return disagreement
      float_val1 = float(val1)
    except:
      return self.disagree_weight
    try:
      float_val2 = float(val2)
    except:
      return self.disagree_weight

    if (float_val1 == float_val2):
      return self.__calc_freq_agree_weight__(val1)

    elif (self.max_perc_diff == 0.0):  # No percentage difference tolerated
      return self.disagree_weight    # Because values are different

    # Calculate percentage difference and weight  - - - - - - - - - - - - - - -
    #
    perc_diff = 100.0 * abs(float_val1 - float_val2) / \
                        max(abs(float_val1), abs(float_val2))

    if (perc_diff > self.max_perc_diff):  # Percentage difference too large
      return self.disagree_weight

    # Get general or frequency based agreement weight
    #
    agree_weight = self.__calc_freq_weights__(val1, val2)

    # Calculate partial agreement weight
    #
    return agree_weight - (perc_diff / (self.max_perc_diff+1.0)) * \
           (agree_weight + abs(self.disagree_weight))

# =============================================================================

class FieldComparatorNumericAbs(FieldComparator):
  """A field comparator for numeric fields, where a given absolute difference
     can be tolerated. The agreement weight is returned if the numbers are the
     same, and the disagreement weight if the absolute difference is larger
     than a maximum tolerated absolute value.

     The disagreement value will be returned if either of the two number
     strings given for a comparison are not valid numbers.

     The additional argument (besides the base class arguments) which has to be
     set when this field comparator is initialised is:

       max_abs_diff  A non-negative floating-point giving the maximum absolute
                     difference tolerated.

     If the absolute difference is equal or smaller than 'max_abs_diff' the
     resulting partial agreement weight is calculated according to the
     following formula:

       weight = agree_weight - (abs_diff / (max_abs_diff + 1.0)) *
                                (agree_weight+abs(disagree_weight))

     where the absolute difference is calculated as:

       abs_diff = abs(value_a - value_b)
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'max_abs_diff' argument first, then call the
       base class constructor.
    """

    self.max_abs_diff = None

    # Process all keyword arguments - - - - - - - - - - - - - - - - - - - - - -
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('max_a')):
        auxiliary.check_is_number('max_abs_diff', value)
        auxiliary.check_is_not_negative('max_abs_diff', value)
        self.max_abs_diff = float(value)

      else:
        base_kwargs[keyword] = value

    FieldComparator.__init__(self, base_kwargs)  # Process base arguments

    # Make sure 'max_abs_diff' attribute is set - - - - - - - - - - - - - - - -
    #
    auxiliary.check_is_number('max_abs_diff', self.max_abs_diff)
    auxiliary.check_is_not_negative('max_abs_diff', self.max_abs_diff)

    self.log([('Maximum absolute difference', self.max_abs_diff)])

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two numerical field values and tolerate an absolute difference.

       If at least one of the two values to be compared is not a number, then
       the disagreement value is returned.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    elif (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # Calculate the absolute difference - - - - - - - - - - - - - - - - - - - -
    #
    try:  # Check if the field values are numbers
      float_val1 = float(val1)
    except:
      return self.disagree_weight
    try:
      float_val2 = float(val2)
    except:
      return self.disagree_weight

    if (float_val1 == float_val2):
      return self.__calc_freq_agree_weight__(val1)

    elif (self.max_abs_diff == 0.0):  # No absolute difference tolerated
      return self.disagree_weight   # Because values are different

    # Calculate absolute difference and weight  - - - - - - - - - - - - - - - -
    #
    abs_diff = abs(float_val1 - float_val2)

    if (abs_diff > self.max_abs_diff):  # Absolute difference too large
      return self.disagree_weight

    # Get general or frequency based agreement weight
    #
    agree_weight = self.__calc_freq_weights__(val1, val2)

    # Calculate partial agreement weight
    #
    return agree_weight - (abs_diff / (self.max_abs_diff+1.0)) * \
           (agree_weight + abs(self.disagree_weight))

# =============================================================================

class FieldComparatorEncodeString(FieldComparator):
  """A field comparator for string fields which are encoded using a phonetic
     encoding method before being compared (exact comparison of the encodings).

     Phonetic encoding routines are implemented in the Febrl module 'encode.py'
     (for more details on the encoding methods please see this module).

     The additional arguments (besides the base class arguments) which have to
     be set when this field comparator is initialised are:

       encode_method    The encoding method, possible are:
                        None            No encoding, use strings directly
                                        (same as exact or truncate string
                                        comparison)
                        'soundex'       Soundex
                        'mod_soundex'   Modified Soundex
                        'phonex'        Phonex
                        'phonix'        Phonix
                        'nysiis'        NYSIIS
                        'dmetaphone'    Double-Metaphone
                        'fuzzysoundex'  Fuzzy Soundex
       reverse          A flag, if set to True the two input strings are first
                        reversed before they are encoded. Default is False.
       max_code_length  Can be used to set the maximal length (in characters)
                        of the codes. Default is 4.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'encode_method', 'reverse', and
       'max_code_length' arguments first, then call the base class constructor.
    """

    self.encode_method =   None
    self.reverse =         False
    self.max_code_length = 4

    # Process all keyword arguments - - - - - - - - - - - - - - - - - - - - - -
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('encode')):
        if (value not in ['', None, 'None', 'none']):
          auxiliary.check_is_string('encode_method', value)
          if (value not in ['soundex','mod_soundex','phonex','phonix','nysiis',
                            'dmetaphone','fuzzysoundex']):
            logging.exception('Value of argument "encode_method" is not ' + \
                              'one of: "soundex", "mod_soundex", "phonex",' + \
                              ' "phonix", "nysiis", "dmetaphone", or ' + \
                              '"fuzzysoundex", or None: %s' % (value))
            raise Exception
          self.encode_method = value

      elif (keyword.startswith('rev')):
        auxiliary.check_is_flag('reverse', value)
        self.reverse = value

      elif (keyword.startswith('max_code')):
        auxiliary.check_is_positive('max_code_length', value)
        self.max_code_length = value

      else:
        base_kwargs[keyword] = value

    FieldComparator.__init__(self, base_kwargs)  # Process base arguments

    # Make sure 'encode_method' attribute is set - - - - - - - - - - - - - - -
    #
    if (self.encode_method != None):
      auxiliary.check_is_string('encode_method', self.encode_method)

    self.log([('Encoding method', str(self.encode_method)),
              ('Reverse flag',self.reverse),
              ('Maximum code length',self.max_code_length)])  # Log a message

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two string values using a phonetic encoding method. If the two
       strings or the encodings of the two strings are the same then return the
       agreement weight, otherwise the disagreement weight.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    if (self.do_caching == True):  # Check if values pair is in the cache
      cache_weight = self.__get_from_cache__(val1, val2)
      if (cache_weight != None):
        return cache_weight

    if (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # Compare the encodings - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.reverse == True):
      rev = list(val1)
      rev.reverse()
      str1 = ''.join(rev)
      rev = list(val2)
      rev.reverse()
      str2 = ''.join(rev)
    else:
      str1 = val1
      str2 = val2

    str1 = str1.lower()  # Encodings assume all lowercase
    str2 = str2.lower()

    max_len=self.max_code_length

    if (self.encode_method == None):
      code1 = str1[:max_len]
      code2 = str2[:max_len]
    elif (self.encode_method == 'soundex'):
      code1 = encode.soundex(str1, maxlen=max_len)
      code2 = encode.soundex(str2, maxlen=max_len)
    elif (self.encode_method == 'mod_soundex'):
      code1 = encode.mod_soundex(str1, maxlen=max_len)
      code2 = encode.mod_soundex(str2, maxlen=max_len)
    elif (self.encode_method == 'phonex'):
      code1 = encode.phonex(str1, maxlen=max_len)
      code2 = encode.phonex(str2, maxlen=max_len)
    elif (self.encode_method == 'phonix'):
      code1 = encode.phonix(str1, maxlen=max_len)
      code2 = encode.phonix(str2, maxlen=max_len)
    elif (self.encode_method == 'nysiis'):
      code1 = encode.nysiis(str1, maxlen=max_len)
      code2 = encode.nysiis(str2, maxlen=max_len)
    elif (self.encode_method == 'dmetaphone'):
      code1 = encode.dmetaphone(str1, maxlen=max_len)
      code2 = encode.dmetaphone(str2, maxlen=max_len)
    elif (self.encode_method == 'fuzzysoundex'):
      code1 = encode.fuzzy_soundex(str1, maxlen=max_len)
      code2 = encode.fuzzy_soundex(str2, maxlen=max_len)
    else:
      logging.exception('Illegal string encoding: %s' % (self.encode_method))
      raise Exception

    # Check if encodings are the same or different  - - - - - - - - - - - - - -
    #
    if (code1 == code2):

      # Get general or frequency based agreement weight
      #
      w = self.__calc_freq_weights__(val1, val2)

    else:
      w = self.disagree_weight

    if (self.do_caching == True):  # Put values pair into the cache
      self.__put_into_cache__(val1, val2, w)

    return w

# =============================================================================

class FieldComparatorDistance(FieldComparator):
  """A field comparator that computes the geographical distance between the two
     given fields. A geocode look-up table needs to be provided, which is a
     dictionary with values (like postcodes or suburb names) as keys and their
     geographic locations (longitude and latitude) as values.

     The additional arguments (besides the base class arguments) which have to
     be set when this field comparator is initialised are:

       geocode_table  A reference to the geocode look-up table (a dictionary)
       max_distance   A positive number that gives the maximum distance (in
                      kilometers) tolerated.

     If the computed distance between the two field values is smaller or equal
     to 'max_distance', the partial agreement weight is calculated using the
     following formula.

       weight = agree_weight - (distance/(max_distance + 1.0)) *
                               (agree_weight+abs(disagree_weight))

     For field values that are not found in the geocode look-up table the
     missing weight is returned.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'geocode_table' and 'max_distance' arguments
       first, then call the base class constructor.
    """

    self.geocode_table = None
    self.max_distance =  None

    # Process all keyword arguments - - - - - - - - - - - - - - - - - - - - - -
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('geoc')):
        auxiliary.check_is_dictionary('geocode_table', value)
        self.geocode_table = value

      elif (keyword.startswith('max_d')):
        auxiliary.check_is_not_negative('max_distance', value)
        self.max_distance = value

      else:
        base_kwargs[keyword] = value

    FieldComparator.__init__(self, base_kwargs)  # Process base arguments

    # Make sure 'geocode_table' and 'max_distance' arguments are set - - - - -
    #
    auxiliary.check_is_dictionary('geocode_table', self.geocode_table)
    auxiliary.check_is_not_negative('max_distance', self.max_distance)

    self.log([('Geocode table length',len(self.geocode_table)),
              ('Maximum distance tolerated', self.max_distance)])

    # Class constants
    #
    self.earth_radius =  6372.0  # Approximate radius of earth in kilometers
    self.deg2rad =       math.pi / 180.0  # Factor for degrees to radians

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values using the distance comparator.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    if (self.do_caching == True):  # Check if values pair is in the cache
      cache_weight = self.__get_from_cache__(val1, val2)
      if (cache_weight != None):
        return cache_weight

    if (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # Calculate distance comparison value - - - - - - - - - - - - - - - - - - -

    # Check if field values are in geocode look-up table
    #
    loc1 = self.geocode_table.get(val1, None)
    loc2 = self.geocode_table.get(val2, None)

    if (loc1 == None) or (loc2 == None):
      w = self.missing_weight  # One or both values are not in look-up table

    else:  # Check if both locations are the same

      if (loc1 == loc2):
        w = self.__calc_freq_agree_weight__(val1)

      else:  # Calculate distance on Earth surface
        long1, lati1 = loc1[0]*self.deg2rad, loc1[1]*self.deg2rad
        long2, lati2 = loc2[0]*self.deg2rad, loc2[1]*self.deg2rad

        alpha = math.cos(long1 - long2)
        x     = alpha * math.cos(lati1)*math.cos(lati2) + \
                        math.sin(lati1)*math.sin(lati2)
        dist  = self.earth_radius * math.acos(x)

        assert dist >= 0.0, 'Distance calculated is: %f' % (dist)

        # Get general or frequency based agreement weight
        #
        agree_weight = self.__calc_freq_weights__(val1, val2)

        # Check if distance is zero or too far
        #
        if (dist <= 0.00001):  # Allow for numerical issues
          w = agree_weight

        elif (dist > self.max_distance):
          w = self.disagree_weight

        else:  # Compute partial agrement weight
          w = agree_weight - (dist / (self.max_distance+1.0)) * \
                             (agree_weight + abs(self.disagree_weight))

    if (self.do_caching == True):  # Put values pair into the cache
      self.__put_into_cache__(val1, val2, w)

    return w

# =============================================================================

class FieldComparatorDate(FieldComparator):
  """A field comparator that compares two dates given as strings of one of
     several forms as detailed below (date_format')

     The disagreement value will be returned if either of the two dates given
     for a comparison are not a valid tuple of numbers.

     Separators (/:;,.) are automatically removed from the dates before they
     are parsed.

     The additional arguments (besides the base class arguments) which have to
     be set when this field comparator is initialised are:

       date_format           A format string that will be used to extract the
                             date values. Possible are:
                               'ddmmyyyy'
                               'mmddyyyy'
                               'yyyymmdd'
                               'ddmmyy'
                               'mmddyy'
       max_day1_before_day2  The maximum tolerated number of days that the
                             first day value can be before the second day value
       max_day2_before_day1  The maximum tolerated number of days that the
                             second day value can be before the first day value

     Default values for both arguments are 0, which results in exact day
     comparison only (i.e. if the dates are not the same the disagreement
     weight will be returned).

     If the first day value is equal to or less than 'max_day1_before_day2'
     days before the second day value then a partial agreement weight will be
     calculated as follows:

       weight = agree_weight - (date_diff/(max_day1_before_day2 + 1)) *
                               (agree_weight+abs(disagree_weight))

     and similar, if the second day value is equal to or less than
     'max_day2_before_day1' days before the first day value then a partial
     agreement weight will be calculated as follows:

       weight = agree_weight - (date_diff/(max_day2_before_day1 + 1)) *
                               (agree_weight+abs(disagree_weight))

     If the day difference is larger than 'max_day1_before_day2' or
     'max_day2_before_day1' then the disagreement weight will be returned.

     For two special cases, partial weights are calculated as follows
     (suggested by Mike Berry <berrym@hln.com>):

     1) A 'swap' partial agreement will be calculated when the day and month
        values are swapped (and the year values are the same) according to the
        following formula:

        weight = agree_weight - 0.5 * (agree_weight+abs(disagree_weight))

        For example: day1 = [12,10,1999] / day2 = [10,12,1999]

     2) If the difference between the days is larger than
        'max_day1_before_day2' or 'max_day2_before_day1', respectively, but
        both the day and year values are the same (i.e. only the month values
        are different), then a partial agreement weight will be calculated
        according to the following formula:

        weight = 0.75 * disagree_weight
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'max_day1_before_day2' and
       'max_day2_before_day1' arguments first, then call the base class
       constructor.
    """

    self.max_day1_before_day2 = 0
    self.max_day2_before_day1 = 0
    self.date_format =          None

    # Process all keyword arguments - - - - - - - - - - - - - - - - - - - - - -
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('max_day1')):
        auxiliary.check_is_not_negative('max_day1_before_day2', value)
        self.max_day1_before_day2 = value

      elif(keyword.startswith('max_day2')):
        auxiliary.check_is_not_negative('max_day2_before_day1', value)
        self.max_day2_before_day1 = value

      elif (keyword.startswith('date_f')):
        auxiliary.check_is_string('date_format', value)
        if (value not in ['ddmmyyyy','mmddyyyy','yyyymmdd','ddmmyy','mmddyy']):
          logging.exception('Illegal date format string: %s' % (value))
          raise Exception
        self.date_format = value

      else:
        base_kwargs[keyword] = value

    FieldComparator.__init__(self, base_kwargs)  # Process base arguments

    # Log a message
    #
    self.log([('Maximum day 1 before day 2 tolerated',
               self.max_day1_before_day2),
              ('Maximum day 2 before day 1 tolerated',
               self.max_day2_before_day1),
              ('Date format', self.date_format)])

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values - assumed to be dates made of triplets (tuples)
       (day,month,year) - using the date comparator.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    if (self.do_caching == True):  # Check if values pair is in the cache
      cache_weight = self.__get_from_cache__(val1, val2)
      if (cache_weight != None):
        return cache_weight

    if (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # Remove separator characters - - - - - - - - - - - - - - - - - - - - - - -
    #
    for c in '/:;,.\\':
      if c in val1:
        val1 = val1.replace(c,'')
      if c in val2:
        val2 = val2.replace(c,'')

    # Parse values to be compared - - - - - - - - - - - - - - - - - - - - - - -
    #
    if ((len(val1) not in [6,8]) or (len(val2) not in [6,8])):
      logging.warn('At least one of the field values is not of correct ' + \
                   'length: %s / %s' % (val1, val2))
      return self.disagree_weight

    if (self.date_format in ['ddmmyyyy', 'ddmmyy']):
      day1, month1, year1 = val1[:2],val1[2:4],val1[4:]
      day2, month2, year2 = val2[:2],val2[2:4],val2[4:]
    elif (self.date_format in ['mmddyyyy','mmddyy']):
      day1, month1, year1 = val1[2:4],val1[:2],val1[4:]
      day2, month2, year2 = val2[2:4],val2[:2],val2[4:]
    elif (self.date_format == 'yyyymmdd'):
      day1, month1, year1 = val1[6:],val1[4:6],val1[:4]
      day2, month2, year2 = val2[6:],val2[4:6],val2[:4]

    # Check values for validity - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if ((month1 < '01') or (month1 > '12') or \
        (month2 < '01') or (month2 > '12')):
      logging.warn('At least one of the month values is out of range: ' + \
                   '%s / %s' % (val1, val2))
      return self.disagree_weight

    days_ok = True

    if ((day1 < '01') or (day2 < '01')):
      days_ok = False
    else:

      if ((month1 in ['01','03','05','07','08','10','12']) and (day1 > '31')):
        days_ok = False
      elif ((month1 in ['04','06','09','11']) and (day1 > '30')):
        days_ok = False
      elif ((month1 == '02') and (day1 > '29')):
        days_ok = False

      if ((month2 in ['01','03','05','07','08','10','12']) and (day2 > '31')):
        days_ok = False
      elif ((month2 in ['04','06','09','11']) and (day2 > '30')):
        days_ok = False
      elif ((month2 == '02') and (day2 > '29')):
        days_ok = False

    if (days_ok == False):
      logging.warn('At least one of the day values is out of range: ' + \
                   '%s / %s' % (val1, val2))
      return self.disagree_weight

    # Convert into integer numbers - - - - - - - - - - - - - - - - - - - - - -
    #
    try:
      int_day1 = int(day1)
    except:
      return self.disagree_weight
    try:
      int_month1 = int(month1)
    except:
      return self.disagree_weight
    try:
      int_year1 = int(year1)
    except:
      return self.disagree_weight

    try:
      int_day2 = int(day2)
    except:
      return self.disagree_weight
    try:
      int_month2 = int(month2)
    except:
      return self.disagree_weight
    try:
      int_year2 = int(year2)
    except:
      return self.disagree_weight

    # Create date objects
    #
    date1 = datetime.date(int_year1, int_month1, int_day1)
    date2 = datetime.date(int_year2, int_month2, int_day2)

    if (date1 == date2):  # Same dates
      return self.__calc_freq_agree_weight__(val1)

    date_diff = date2 - date1
    day_diff =  date_diff.days

    # Get general or frequency based agreement weight
    #
    agree_weight = self.__calc_freq_weights__(val1, val2)

    # Check for swapped day and month values first - - - - - - - - - - - - - -
    # (berrym@hln.com: added "swap" agreement where month/day are swapped)
    #
    if ((date1.day == date2.month) and (date1.month == date2.day) and \
        (date1.year == date2.year)):
      w = agree_weight-0.5*(agree_weight+abs(self.disagree_weight))

    # Check if day difference is in the permitted tolerance range - - - - - - -
    #
    elif ((day_diff <= self.max_day1_before_day2) and \
          (-day_diff <= self.max_day2_before_day1)):

      if (day_diff > 0):
        w = agree_weight - \
            (float(day_diff) / (self.max_day1_before_day2+1.0)) * \
            (agree_weight + abs(self.disagree_weight))
      else:  # (day_diff < 0)
        w = agree_weight - \
            (float(-day_diff) / (self.max_day2_before_day1+1.0)) * \
            (agree_weight + abs(self.disagree_weight))

    # Day difference too large - - - - - - - - - - - - - - - - - - - - - - -
    #
    else:

      # berrym@hln.com: added "Day not month agreement", if day and year are
      # the same, don't penalize quite so badly
      #
      if ((date1.day == date2.day) and (date1.year == date2.year)):
        w = self.disagree_weight * 0.75
      else:
        w = self.disagree_weight

    if (self.do_caching == True):  # Put values pair into the cache
      self.__put_into_cache__(val1, val2, w)

    return w

# =============================================================================

class FieldComparatorTime(FieldComparator):
  """A field comparator for time fields, which must be given in 24 hours format
     (00:00 is midnight and 23:59 is 11:59 pm).

     The field values must be strings of the form 'HHMM' or 'HH:MM', otherwise
     an error is triggered.

     The disagreement value will be returned if either of the two time strings
     given for a comparison are not valid.

     The additional arguments (besides the base class arguments) which have to
     be set when this field comparator is initialised are:

       max_time1_before_time2  The maximum tolerated time in minutes that the
                               first time value can be before the second time
                               value
       max_time2_before_time1  The maximum tolerated time in minutes that the
                               second time value can be before the first time
                               value
       day_start               The time (given as a string 'HHMM') where a
                               24-hours period starts. It is then assumed that
                               both times values (when comparing them) are
                               within the same 24-hours period. The default
                               value is midnight (00:00).

     Default values for both arguments 'max_time1_before_time2' and
     'max_time2_before_time1' are 0, which results in exact time comparison
     only (i.e. if the times are not the same the disagreement weight will be
     returned).

     If the first time value is equal to or less than 'max_time1_before_time2'
     minutes before the second time value then a partial agreement weight will
     be calculated as follows:

       weight = agree_weight - (time_diff/(max_time1_before_time2 + 1)) *
                               (agree_weight+abs(disagree_weight))

     and similar, if the second time value is equal to or less than
     'max_time2_before_time1' minutes before the first time value then a
     partial agreement weight will be calculated as follows:

       weight = agree_weight - (time_diff/(max_time2_before_time1 + 1)) *
                               (agree_weight+abs(disagree_weight))

     If the time difference is larger than 'max_time1_before_time2' or
     'max_time2_before_time1' then the disagreement weight will be returned.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'max_time1_before_time2' and
       'max_time2_before_time1' arguments first, then call the base class
       constructor.
    """

    self.max_time1_before_time2 = 0
    self.max_time2_before_time1 = 0
    self.day_start =              0000  # Default value midnight

    # Process all keyword arguments - - - - - - - - - - - - - - - - - - - - - -
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('max_time1')):
        auxiliary.check_is_not_negative('max_time1_before_time2', value)
        self.max_time1_before_time2 = value

      elif (keyword.startswith('max_time2')):
        auxiliary.check_is_not_negative('max_time2_before_time1', value)
        self.max_time2_before_time1 = value

      elif (keyword.startswith('day_s')):
        self.day_start = value

      else:
        base_kwargs[keyword] = value

    FieldComparator.__init__(self, base_kwargs)  # Process base arguments

    # Process the day start - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (isinstance(self.day_start,str)):
      if (len(self.day_start) == 5) and (':' in self.day_start):
        daystart_str = self.day_start.replace(':','')
      elif (len(self.day_start) == 4):
        daystart_str = self.day_start
      else:
        logging.exception('Day start is not a string of length 4 or 5: "%s"' \
                          (str(self.day_start)))
        raise Exception

      if (not daystart_str.isdigit()):
        logging.exception('Day start is not a string made of digits: %s' % \
                          (str(daystart_str)))
        raise Exception

      hrs,min = int(daystart_str[:2]),int(daystart_str[2:])

      if (hrs < 0) or (hrs > 23) or (min < 0) or (min > 59):
        logging.exception('Day start time value out of range: "%s"' % \
                          (daystart_str))
        raise Exception

      self.day_start = (hrs * 60) + min

    elif (not isinstance(self.day_start,int)):
      logging.exception('Day start is not a string of length 4 or 5 nor an' + \
                        'integer number: "%s"' % (str(self.day_start)))
      raise Exception

    # Log a message
    #
    self.log([('Maximum time 1 before time 2 tolerated (in minutes)',
               self.max_time1_before_time2),
              ('Maximum time 2 before time 1 tolerated (in minutes)',
               self.max_time2_before_time1),
              ('Day start value (in minutes)', self.day_start)])

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values using the time comparator.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    if (self.do_caching == True):  # Check if values pair is in the cache
      cache_weight = self.__get_from_cache__(val1, val2)
      if (cache_weight != None):
        return cache_weight

    if (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # Calculate time comparison value - - - - - - - - - - - - - - - - - - - - -
    #
    if (isinstance(val1,str)):
      if (len(val1) == 5) and (':' in val1):
        time1str = val1.replace(':','')
      elif (len(val1) == 4):
        time1str = val1
      else:
        logging.warn('Value 1 is not a string of length 4 or 5: "%s"' % \
                     (str(val1)))
        return self.disagree_weight
    else:
      logging.warn('Value 1 is not a string: "%s"' % (str(val1)))
      return self.disagree_weight

    if (isinstance(val2,str)):
      if (len(val2) == 5) and (':' in val2):
        time2str = val2.replace(':','')
      elif (len(val2) == 4):
        time2str = val2
      else:
        logging.warn('Value 2 is not a string of length 4 or 5: "%s"' % \
                     (str(val2)))
        return self.disagree_weight
    else:
      logging.exception('Value 2 is not a string: "%s"' % (str(val2)))
      return self.disagree_weight

    if ((not time1str.isdigit()) or (not time2str.isdigit())):
      logging.warn('At least one of the field values is not a string made ' + \
                   'of digits only: %s / %s' % (str(time1str),str(time2str)))
      return self.disagree_weight

    hrs1,min1 = int(time1str[:2]),int(time1str[2:])
    hrs2,min2 = int(time2str[:2]),int(time2str[2:])

    if (hrs1 < 0) or (hrs1 > 23) or (min1 < 0) or (min1 > 59):
      logging.warn('Time 1 value out of range: "%s"' % (time1str))
      return self.disagree_weight
    if (hrs2 < 0) or (hrs2 > 23) or (min2 < 0) or (min2 > 59):
      logging.warn('Time 2 value out of range: "%s"' % (time2str))
      return self.disagree_weight

    time1 = (hrs1 * 60) + min1  # Convert into minute values
    time2 = (hrs2 * 60) + min2

    # Get general or frequency based agreement weight
    #
    agree_weight = self.__calc_freq_weights__(val1, val2)

    # Check if times are the same - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (time1 == time2):
        w = agree_weight

    elif ((self.max_time1_before_time2 == 0) and \
          (self.max_time2_before_time1 == 0)):
      w = self.disagree_weight  # Because values are different

    else:  # Check for time difference

      # Convert into times according to 'day_start' value first
      #
      time1 -= self.day_start
      if (time1 < 0):
        time1 += 1440  # Adjust into 24-hours period
      time2 -= self.day_start
      if (time2 < 0):
        time2 += 1440  # Adjust into 24-hours period

      # Times are now in a 24-hours period (values between 0 and 1439)  - - -
      #
      time_diff = time2 - time1

      # Check if time difference is in the permitted tolerance range - - - - -
      #
      if ((time_diff <= self.max_time1_before_time2) and \
          (-time_diff <= self.max_time2_before_time1)):

        if (time_diff > 0):
          w = agree_weight - \
                   (float(time_diff) / (self.max_time1_before_time2+1.0)) * \
                   (agree_weight + abs(self.disagree_weight))
        else:  # (time_diff < 0)
          w = agree_weight - \
                  (float(-time_diff) / (self.max_time2_before_time1+1.0)) * \
                   (agree_weight + abs(self.disagree_weight))

      # Day difference too large - - - - - - - - - - - - - - - - - - - - - - -
      #
      else:
        w = self.disagree_weight

    if (self.do_caching == True):  # Put values pair into the cache
      self.__put_into_cache__(val1, val2, w)

    return w

# =============================================================================

class FieldComparatorAge(FieldComparator):
  """A field comparator that compares ages given as strings of one of several
     forms as detailed below (date_format'). The 'age' of these dates is
     calculated according to a fix date.

     The disagreement value will be returned if either of the two dates given
     for a comparison are not a valid tuple of numbers.

     Separators (/:;,.) are automatically removed from the dates before they
     are parsed.

     The additional arguments (besides the base class arguments) which have to
     be set when this field comparator is initialised are:

       date_format    A format string that will be used to extract the date
                      values. Possible are:
                        'ddmmyyyy'
                        'mmddyyyy'
                        'yyyymmdd'
                        'ddmmyy'
                        'mmddyy'
       fix_date       A date relative to which ages are computed. This argument
                      can either be given as a date triplet (a tuple) of the
                      form (day,month,year), or as the string 'today' (in which
                      case the system date is taken as fix date). Default is
                      'today'.
       max_perc_diff  A floating-point number between 0.0 and 100.0 giving the
                      maximum age percentage difference tolerated.

     If the percentage difference is smaller than 'max_perc_diff' the resulting
     partial agreement weight is calculated according to the following formula:

       weight = agree_weight - (perc_diff / (max_perc_diff + 1.0)) *
                                (agree_weight+abs(disagree_weight))

     where the percentage difference is calculated as:

       perc_diff = 100.0 * abs(age1 - age2) / max(abs(age1), abs(age2))
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'fix_date' and 'max_perc_diff' arguments first,
       then call the base class constructor.
    """

    self.fix_date =      'today'
    self.date_format =   None
    self.max_perc_diff = None

    # Process all keyword arguments - - - - - - - - - - - - - - - - - - - - - -
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('fix')):
        self.fix_date = value

      elif (keyword.startswith('max_p')):
        auxiliary.check_is_percentage('max_perc_diff', value)
        self.max_perc_diff = float(value)

      elif (keyword.startswith('date_f')):
        auxiliary.check_is_string('date_format', value)
        if (value not in ['ddmmyyyy','mmddyyyy','yyyymmdd','ddmmyy','mmddyy']):
          logging.exception('Illegal date format string: %s' % (value))
          raise Exception
        self.date_format = value

      else:
        base_kwargs[keyword] = value

    FieldComparator.__init__(self, base_kwargs)  # Process base arguments

    # Make sure 'max_perc_diff' attribute is set - - - - - - - - - - - - - - -
    #
    auxiliary.check_is_percentage('max_perc_diff', self.max_perc_diff)

    # Process the fix date - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.fix_date == 'today'):
      self.fix_date_val = datetime.date.today()  # Today as (year,month,day)

    else:
      if (not isinstance(self.fix_date,tuple)):
        logging.exception('Fix date is not "today" nor a tuple: %s' % \
                          (str(self.fix_date)))
        raise Exception
      if (len(self.fix_date) != 3):
        logging.exception('Fix date is not a triplet: %s' % \
                          (str(self.fix_date)))
        raise Exception

      self.fix_date_val = datetime.date(self.fix_date[2], self.fix_date[1], \
                                        self.fix_date[0])

    # Log a message
    #
    self.log([('Fix date', self.fix_date),('Fix date value',self.fix_date_val),
              ('Date format', self.date_format),
              ('Maximum percentage difference', self.max_perc_diff)])

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values - assumed to be dates made of triplets (tuples)
       (day,month,year) -  using the age comparator.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    if (self.do_caching == True):  # Check if values pair is in the cache
      cache_weight = self.__get_from_cache__(val1, val2)
      if (cache_weight != None):
        return cache_weight

    if (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # Remove separator characters - - - - - - - - - - - - - - - - - - - - - - -
    #
    for c in '/:;,.\\':
      if c in val1:
        val1 = val1.replace(c,'')
      if c in val2:
        val2 = val2.replace(c,'')

    # Parse values to be compared - - - - - - - - - - - - - - - - - - - - - - -
    #
    if ((len(val1) not in [6,8]) or (len(val2) not in [6,8])):
      logging.warn('At least one of the field values is not of correct ' + \
                   'length: %s / %s' % (val1, val2))
      return self.disagree_weight

    if (self.date_format in ['ddmmyyyy', 'ddmmyy']):
      day1, month1, year1 = val1[:2],val1[2:4],val1[4:]
      day2, month2, year2 = val2[:2],val2[2:4],val2[4:]
    elif (self.date_format in ['mmddyyyy','mmddyy']):
      day1, month1, year1 = val1[2:4],val1[:2],val1[4:]
      day2, month2, year2 = val2[2:4],val2[:2],val2[4:]
    elif (self.date_format == 'yyyymmdd'):
      day1, month1, year1 = val1[6:],val1[4:6],val1[:4]
      day2, month2, year2 = val2[6:],val2[4:6],val2[:4]

    # Check values for validity - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if ((month1 < '01') or (month1 > '12') or \
        (month2 < '01') or (month2 > '12')):
      logging.warn('At least one of the month values is out of range: ' + \
                   '%s / %s' % (val1, val2))
      return self.disagree_weight

    days_ok = True

    if ((day1 < '01') or (day2 < '01')):
      days_ok = False
    else:

      if ((month1 in ['01','03','05','07','08','10','12']) and (day1 > '31')):
        days_ok = False
      elif ((month1 in ['04','06','09','11']) and (day1 > '30')):
        days_ok = False
      elif ((month1 == '02') and (day1 > '29')):
        days_ok = False

      if ((month2 in ['01','03','05','07','08','10','12']) and (day2 > '31')):
        days_ok = False
      elif ((month2 in ['04','06','09','11']) and (day2 > '30')):
        days_ok = False
      elif ((month2 == '02') and (day2 > '29')):
        days_ok = False

    if (days_ok == False):
      logging.warn('At least one of the day values is out of range: ' + \
                   '%s / %s' % (val1, val2))
      return self.disagree_weight

    # Convert into integer numbers - - - - - - - - - - - - - - - - - - - - - -
    #
    try:
      int_day1 = int(day1)
    except:
      return self.disagree_weight
    try:
      int_month1 = int(month1)
    except:
      return self.disagree_weight
    try:
      int_year1 = int(year1)
    except:
      return self.disagree_weight

    try:
      int_day2 = int(day2)
    except:
      return self.disagree_weight
    try:
      int_month2 = int(month2)
    except:
      return self.disagree_weight
    try:
      int_year2 = int(year2)
    except:
      return self.disagree_weight

    # Create date objects
    #
    date1 = datetime.date(int_year1, int_month1, int_day1)
    date2 = datetime.date(int_year2, int_month2, int_day2)

    if (date1 == date2):  # Same date
      return self.__calc_freq_agree_weight__(val1)

    age1 = (self.fix_date_val - date1).days  # Get ages in days
    age2 = (self.fix_date_val - date2).days

    # Check age values  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (age1 < 0) or (age1 > 54750):  # Number of days in 150 years
      logging.warning('Illegal value for age 1: %d days' % (age1) + \
                      ', set weight to missing value')
      w = self.disagree_weight

    elif (age2 < 0) or (age2 > 54750):
      logging.warning('Illegal value for age 2: %d days' % (age2) + \
                      ', set weight to missing value')
      w = self.disagree_weight

    else:

      # Calculate age percentage difference and weight  - - - - - - - - - - - -
      #
      if (self.max_perc_diff == 0.0):  # No percentage tolerance allowed
        w =  self.disagree_weight  # Because age values are different

      else:
        perc_diff = 100.0*float(abs(age1 - age2)) / max(abs(age1), abs(age2))

        if (perc_diff > self.max_perc_diff):  # Percentage diff too large
          w = self.disagree_weight

        else:

          # Get general or frequency based agreement weight
          #
          agree_weight = self.__calc_freq_weights__(val1, val2)

          w = agree_weight - (perc_diff / (self.max_perc_diff+1.0)) * \
              (agree_weight + abs(self.disagree_weight))

    if (self.do_caching == True):  # Put values pair into the cache
      self.__put_into_cache__(str(val1), str(val2), w)

    return w

# =============================================================================

# All the comparators below are approximate string comparators, and as such
# sub-classes of the class 'FieldComparatorApproxString' which is a sub-class
# of the base class 'FieldComparator'.

# =============================================================================

class FieldComparatorApproxString(FieldComparator):
  """A generic field comparator based on a approximate string comparator. A
     number of specific approximate comparison methods are based on this class
     (and implemented as sub-classes).

     The additional argument (besides the base class arguments) which has to be
     set when this field comparator is initialised is:

       threshold  A floating-point number between 0.0 and 1.0 giving the
                  minimum approximate value that is tolerated.

     If the approximate string comparator calculates a similarity value between
     1.0 and 'threshold', then the partial agreement weight is calculated using
     the following formula:

       weight = agree_weight - (1.0 - approx_str_val) / (1.0 - threshold) *
                               (agree_weight+abs(disagree_weight))

     If the approximate approximate string comparator calculates a similarity
     value less than the 'threshold' then the disagreement weight is returned.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, kwargs):
    """Constructor. Process the 'threshold' argument first, then call the
       base class constructor.
    """

    self.threshold = None

    # Process all keyword arguments - - - - - - - - - - - - - - - - - - - - - -
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('thres')):
        auxiliary.check_is_number('threshold', value)
        auxiliary.check_is_normalised('threshold', value)
        if (value == 1.0):
          logging.exception('Value of argument "threshold" must be smaller' + \
                            'then 1.0: %s' % (value))
          raise Exception
        self.threshold = float(value)

      else:
        base_kwargs[keyword] = value

    FieldComparator.__init__(self, base_kwargs)  # Process base arguments

    # Make sure the 'threshold' attribute is set - - - - - - - - - - - - - - -
    #
    auxiliary.check_is_number('threshold', self.threshold)
    auxiliary.check_is_normalised('threshold', self.threshold)

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two fields values, compute and return a numerical weight. See
       implementations in derived sub-classes for details.
    """

    logging.exception('Override abstract method in derived class')
    raise Exception

  # ---------------------------------------------------------------------------

  def __calc_partagree_weight__(self, val1, val2, approx_sim_val):
    """Calculate the partial agreement weight. Should not be used from outside
       the module.
    """

    # Check if similarity values is below or above threshold
    #
    if (approx_sim_val < self.threshold):
      return self.disagree_weight

    # Get general or frequency based agreement weight
    #
    agree_weight = self.__calc_freq_weights__(val1, val2)

    # Compute final adjusted weight - - - - - - - - - - - - - - - - - - - - -
    # Modified after formula in Winkler and Thibaudeau, 1991, page 12
    #
    return agree_weight-(1.0-approx_sim_val) / (1.0-self.threshold) * \
           (agree_weight + abs(self.disagree_weight))

# =============================================================================

class FieldComparatorJaro(FieldComparatorApproxString):
  """A field comparator based on the Jaro approximate string comparator.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """No additional attributes needed, just call base class constructor.
    """

    FieldComparatorApproxString.__init__(self, kwargs)

    self.log([('Threshold', self.threshold)])  # Log a message

    self.JARO_MARKER_CHAR = chr(1)  # Special character used to mark assigned
                                    # characters

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values using the Jaro approximate string comparator.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    if (self.do_caching == True):  # Check if values pair is in the cache
      cache_weight = self.__get_from_cache__(val1, val2)
      if (cache_weight != None):
        return cache_weight

    if (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # Calculate Jaro similarity value - - - - - - - - - - - - - - - - - - - - -
    #
    len1, len2 = len(val1), len(val2)
    halflen = max(len1,len2) / 2 - 1  # Or + 1?? PC 12/03/2009

    ass1, ass2 = '', ''  # Characters assigned in string 1 and string 2
    workstr1, workstr2 = val1, val2  # Copies of the original strings

    common1, common2 = 0.0, 0.0  # Number of common characters

    for i in range(len1):  # Analyse the first string
      start = max(0,i-halflen)
      end   = min(i+halflen+1,len2)
      index = workstr2.find(val1[i],start,end)
      if (index > -1):  # Found common character
        common1 += 1
        ass1 = ass1 + val1[i]
        workstr2 = workstr2[:index]+self.JARO_MARKER_CHAR+workstr2[index+1:]

    for i in range(len2):  # Analyse the second string
      start = max(0,i-halflen)
      end   = min(i+halflen+1,len1)
      index = workstr1.find(val2[i],start,end)
      if (index > -1):  # Found common character
        common2 += 1
        ass2 = ass2 + val2[i]
        workstr1 = workstr1[:index]+self.JARO_MARKER_CHAR+workstr1[index+1:]

    assert (common1 == common2), 'Jaro: Different "common" values'

    if (common1 == 0.0):  # No characters in common
      w = self.disagree_weight

    else:  # Compute number of transpositions  - - - - - - - - - - - - - - - -

      transp = 0.0
      for i in range(len(ass1)):
        if (ass1[i] != ass2[i]):
          transp += 0.5

      w = 1./3.*(common1 / float(len1) + common1 / float(len2) + \
          (common1-transp) / common1)

      assert (w > 0.0), 'Jaro: Weight is smaller than 0.0: %f' % (w)
      assert (w < 1.0), 'Jaro: Weight is larger than 1.0: %f' % (w)

      w = self.__calc_partagree_weight__(val1, val2, w)

    if (self.do_caching == True):  # Put values pair into the cache
      self.__put_into_cache__(val1, val2, w)

    return w

# =============================================================================

class FieldComparatorWinkler(FieldComparatorApproxString):
  """A field comparator based on the Jaro-Winkler approximate string
     comparator.

     The comparison method is based on the C code 'strcmp95' as available in
     the US Census Bureau BigMatch program.

     The additional arguments (besides the base class arguments) which have to
     be set when this field comparator is initialised are:

       check_sim   Check for similar character pairs and modify similarity
                   measure accordingly.
       check_init  Check for same initial characters (up to 4) and modify
                   similarity measure accordingly.
       check_long  Check for more agreeing characters in long strings and
                   modify similarity measure accordingly.

     These three arguments are flags that can be set to True or False. Default
     values are True.

       multi_word  An extension allowing to handle strings containing multiple
                   words. Possible actions are:
                   'sort'  Words are first sorted before the Winkler similarity
                           measure is calculated
                   'perm'  All permutations of words are given to the  Winkler
                           comparator and the highest similarity value is
                           returned (this can be slow for values containing
                           many words)
                   None    Do nothing, simply give values to Winkler comparator
                           (this is the default).
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'check_sim', 'check_init' and 'check_long'
       arguments first, then call the base class constructor.
    """

    self.check_sim =  True
    self.check_init = True
    self.check_long = True
    self.multi_word = None

    # Process all keyword arguments - - - - - - - - - - - - - - - - - - - - - -
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('check_s')):
        auxiliary.check_is_flag('check_sim', value)
        self.check_sim = value

      elif (keyword.startswith('check_i')):
        auxiliary.check_is_flag('check_init', value)
        self.check_init = value

      elif (keyword.startswith('check_l')):
        auxiliary.check_is_flag('check_long', value)
        self.check_long = value

      elif (keyword.startswith('multi')):
        if (value is None):
          self.multi_word = None
        else:
          auxiliary.check_is_string('multi_word', value)
          if (value not in ['sort','perm']):
            logging.exception('Value of argument "multi_word" is not one' + \
                            'of: "sort" or "perm": %s' % (value))
            raise Exception
          self.multi_word = value

      else:
        base_kwargs[keyword] = value

    FieldComparatorApproxString.__init__(self, base_kwargs)

    self.log([('Threshold', self.threshold),
              ('Multi word handling', self.multi_word),
              ('Check similar characters flag', self.check_sim),
              ('Check initial same characters flag', self.check_init),
              ('Check long strings flag', self.check_long)])  # Log a message

    self.JARO_MARKER_CHAR = chr(1)  # Special character used to mark assigned
                                    # characters

    # Taken from US Census Bureau BigMatch C code 'stringcmp'
    #
    self.sim_char_pairs = frozenset([('a','e'),('e','a'),('a','i'),('i','a'),
                                     ('a','o'),('o','a'),('a','u'),('u','a'),
                                     ('b','v'),('v','b'),('e','i'),('i','e'),
                                     ('e','o'),('o','e'),('e','u'),('u','e'),
                                     ('i','o'),('o','i'),('i','u'),('u','i'),
                                     ('o','u'),('u','o'),('i','y'),('y','i'),
                                     ('e','y'),('y','e'),('c','g'),('g','c'),
                                     ('e','f'),('f','e'),('w','u'),('u','w'),
                                     ('w','v'),('v','w'),('x','k'),('k','x'),
                                     ('s','z'),('z','s'),('x','s'),('s','x'),
                                     ('q','c'),('c','q'),('u','v'),('v','u'),
                                     ('m','n'),('n','m'),('l','i'),('i','l'),
                                     ('q','o'),('o','q'),('p','r'),('r','p'),
                                     ('i','j'),('j','i'),('2','z'),('z','2'),
                                     ('5','s'),('s','5'),('8','b'),('b','8'),
                                     ('1','i'),('i','1'),('1','l'),('l','1'),
                                     ('0','o'),('o','0'),('0','q'),('q','o'),
                                     ('c','k'),('k','c'),('g','j'),('j','g'),
                                     ('e',' '),(' ','e'),('y',' '),(' ','y'),
                                     ('s',' '),(' ','s')])

  # ---------------------------------------------------------------------------

  def __do_winkler__(self, val1, val2):
    """Calculate basic Winkler similarity measure for two input strings.

       Should not be used from outside the module.
    """

    len1, len2 = len(val1), len(val2)

    if (len1 < 4) or (len2 < 4):  # Both strings must be at least 4 chars long
      return self.disagree_weight

    halflen = max(len1,len2) / 2 - 1  # Or + 1?? PC 12/03/2009

    ass1, ass2 = '', ''  # Characters assigned in string 1 and string 2
    workstr1, workstr2 = val1, val2  # Copies of the original strings

    common1, common2 = 0.0, 0.0  # Number of common characters

    # Analyse the first string  - - - - - - - - - - - - - - - - - - - - - - - -
    #
    for i in range(len1):
      start = max(0,i-halflen)
      end   = min(i+halflen+1,len2)
      index = workstr2.find(val1[i],start,end)
      if (index > -1):  # Found common character
        common1 += 1
        ass1 = ass1 + val1[i]
        workstr2 = workstr2[:index]+self.JARO_MARKER_CHAR+workstr2[index+1:]

    # Analyse the second string - - - - - - - - - - - - - - - - - - - - - - - -
    #
    for i in range(len2):
      start = max(0,i-halflen)
      end   = min(i+halflen+1,len1)
      index = workstr1.find(val2[i],start,end)
      if (index > -1):  # Found common character
        common2 += 1
        ass2 = ass2 + val2[i]
        workstr1 = workstr1[:index]+self.JARO_MARKER_CHAR+workstr1[index+1:]

    assert (common1 == common2), 'Winkler: Different "common" values'

    if (common1 == 0.0):  # No characters in common
      return self.disagree_weight

    # Compute number of transpositions  - - - - - - - - - - - - - - - - - - - -
    #
    transp = 0.0
    for i in range(len(ass1)):
      if (ass1[i] != ass2[i]):
        transp += 0.5

    # Check for similarities in non-matched characters - - - - - - - - - - - -
    #
    if (self.check_sim == True):

      sim_weight = 0.0

      workstr1 = workstr1.replace(self.JARO_MARKER_CHAR ,'')  # Remove assigned
      workstr2 = workstr2.replace(self.JARO_MARKER_CHAR ,'')  # characters

      for c1 in workstr1:
        for j in range(len(workstr2)):
          if (c1,workstr2[j]) in self.sim_char_pairs:
            sim_weight += 3
            workstr2 = workstr2[:j]+self.JARO_MARKER_CHAR+workstr2[j+1:]
            break       # Mark character as used

      common1 += sim_weight / 10.0

    # Calculate basic (Jaro) weight
    #
    w = 1./3.*(common1 / float(len1) + common1 / float(len2) + \
               (common1-transp) / common1)

    assert (w > 0.0), 'Winkler: Basic weight is smaller than 0.0: %f' % (w)
    assert (w < 1.0), 'Winkler: Basic weight is larger than 1.0: %f' % (w)

    # Check for same characters at the beginning - - - - - - - - - - - - - - -
    #
    same_init = 0  # Variable needed later on

    if (self.check_init == True):

      minlen = min(len1, len2, 4)

      for same_init in range(1, minlen+1):
        if (val1[:same_init] != val2[:same_init]):
          break
      same_init -= 1

      assert (same_init >= 0), 'Winkler: "same_init" value smaller than 0'
      assert (same_init <= 4), 'Winkler: "same_init" value larger than 4'

      w += same_init*0.1 * (1.0 - w)

    # Check for long strings and possibly adjust weight - - - - - - - - - - - -
    #
    if (self.check_long == True):

      if ((min(len1,len2) > 4) and (common1 > same_init+1) and \
          (common1 >= min(len1,len2)+same_init)):
        w_mod = w + (1.0-w) * (common1-same_init-1) / \
                    (float(len1)+float(len2)-same_init*2+2)
                                                # Fixed -2 => +2 PC 12/03/2009

        assert (w_mod >= w), 'Winkler: Long string adjustment decreases weight'
        assert (w_mod < 1.0), 'Winkler: Long strings adjustment weight > 1.0'

        w = w_mod

    return w

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values using the Winkler approximate string
       comparator.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    if (self.do_caching == True):  # Check if values pair is in the cache
      cache_weight = self.__get_from_cache__(val1, val2)
      if (cache_weight != None):
        return cache_weight

    if (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # Calculate Jaro-Winkler similarity value - - - - - - - - - - - - - - - - -

    # Check if multi-word handling is activated and then for whitespaces
    #
    if (self.multi_word != None) and ((' ' in val1) or (' 'in val2)):

      if (self.multi_word == 'sort'):  # Sort words first

        if (' ' in val1):  # Sort value 1
          val1_list = val1.split(' ')
          val1_list.sort()
          val1sorted = ' '.join(val1_list)
        else:
          val1sorted = val1

        if (' ' in val2):  # Sort value 2
          val2_list = val2.split(' ')
          val2_list.sort()
          val2sorted = ' '.join(val2_list)
        else:
          val2sorted = val1

        w = self.__do_winkler__(val1sorted, val2sorted)

      else:  # Create all permutations
        val1_list = val1.split(' ')
        val2_list = val2.split(' ')

        perm_list1 = mymath.permute(val1_list)
        perm_list2 = mymath.permute(val2_list)

        w = -1.0  # Maximal similarity measure
        max_perm = None

        for perm1 in perm_list1:
          for perm2 in perm_list2:

            # Calculate Winkler similarity measure for this permutation
            #
            this_w = self.__do_winkler__(perm1, perm2)

            if (this_w > w):
              w        = this_w
              max_perm = [perm1, perm2]

        logging.debug('Permutation Winkler best permutation: %s' % \
                      (str(max_perm)))

    else:  # No multi word handling or no whitespaces in values

      w = self.__do_winkler__(val1, val2)

    w = self.__calc_partagree_weight__(val1, val2, w)

    if (self.do_caching == True):  # Put values pair into the cache
      self.__put_into_cache__(val1, val2, w)

    return w

# =============================================================================

class FieldComparatorQGram(FieldComparatorApproxString):
  """A field comparator based on the q-gram approximate string comparator.

     q-grams are the sub-strings containing q characters in a string. For
     example, 'peter' contains the bigrams (q=2): ['pe','et','te','er'].

     This method counts the number of common q-grams and divides by either the
     minimum, average, or maximum number of q-grams. The resulting similarity
     weight is returned.

     The additional arguments (besides the base class arguments) which have to
     be set when this field comparator is initialised are:

       q               The length of the q-grams to be used (must be at least
                       1). The default value is 2 (i.e. bigrams)
       common_divisor  Method of how to calculate the divisor. Can be set to
                       'average','shortest', or 'longest', and is calculated
                       according to the lengths of the two input strings.
       padded          If set to True (default), the beginnng and end of the
                       strings will be padded with (q-1) special characters, if
                       False no padding will be done.

     Padding will result in specific q-grams at the beginning and end of a
     string, for example 'peter' converted into padded bigrams (q=2) will
     result in the following 2-gram list: ['*p','pe','et','te','er','r@'], with
     '*' illustrating the start and '@' the end character.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'q', 'common_divisor' and 'padded' arguments
       first, then call the base class constructor.
    """

    self.q =              2
    self.common_divisor = None
    self.padded =         True

    # Process all keyword arguments - - - - - - - - - - - - - - - - - - - - - -
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword == 'q'):
        auxiliary.check_is_positive('q', value)
        self.q = value

      elif (keyword.startswith('common_d')):
        auxiliary.check_is_string('common_divisor', value)
        if (value not in ['average','shortest','longest']):
          logging.exception('Value of argument "common_divisor" is not one' + \
                            'of: "average", "shortest", or "longest": %s' % \
                            (value))
          raise Exception
        self.common_divisor = value

      elif (keyword.startswith('pad')):
        auxiliary.check_is_flag('padded', value)
        self.padded = value

      else:
        base_kwargs[keyword] = value

    FieldComparatorApproxString.__init__(self, base_kwargs)

    # Make sure 'common_divisor' attribute is set  - - - - - - - - - - - - - -
    #
    auxiliary.check_is_string('common_divisor', self.common_divisor)

    self.log([('Threshold', self.threshold),
              ('q', self.q),
              ('Common divisor', self.common_divisor),
              ('Padded flag', self.padded)])  # Log a message

    self.QGRAM_START_CHAR = chr(1)
    self.QGRAM_END_CHAR =   chr(2)

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values using the q-gram approximate string comparator.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    if (self.do_caching == True):  # Check if values pair is in the cache
      cache_weight = self.__get_from_cache__(val1, val2)
      if (cache_weight != None):
        return cache_weight

    if (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # Calculate q-gram similarity value - - - - - - - - - - - - - - - - - - - -
    #
    q = self.q  # Faster access

    # Calculate number of q-grams in strings (plus start and end characters)
    #
    if (self.padded == True):
      num_qgram1 = len(val1)+q-1
      num_qgram2 = len(val2)+q-1
    else:
      num_qgram1 = max(len(val1)-(q-1), 0)  # Make sure its not negative
      num_qgram2 = max(len(val2)-(q-1), 0)

    # Check if there are q-grams at all from both strings - - - - - - - - - - -
    # (no q-grams if length of a string is less than q)
    #
    if ((self.padded == False) and (min(num_qgram1, num_qgram2) == 0)):
      w = self.disagree_weight

    else:
      if (self.common_divisor == 'average'):  # Calculate the divisor
        divisor = 0.5*(num_qgram1 + num_qgram2)
      elif (self.common_divisor == 'shortest'):
        divisor = min(num_qgram1, num_qgram2)
      else:  # Longest
        divisor = max(num_qgram1, num_qgram2)

      # Use number of q-grams to quickly check if below threshold - - - - - - -
      #
      max_common_qgram = min(num_qgram1, num_qgram2)  # Max possible q-grams
      w = float(max_common_qgram) / float(divisor)    #   in common

      if (w  < self.threshold):  # Similariy is smaller than threshold
        w = self.disagree_weight

      else:

        # Add start and end characters (padding) - - - - - - - - - - - - - - -
        #
        if (self.padded == True):
          qgram_str1 = (q-1)*self.QGRAM_START_CHAR+val1+(q-1)* \
                       self.QGRAM_END_CHAR
          qgram_str2 = (q-1)*self.QGRAM_START_CHAR+val2+(q-1)* \
                       self.QGRAM_END_CHAR
        else:
          qgram_str1 = val1
          qgram_str2 = val2

        # Make a list of q-grams for both strings - - - - - - - - - - - - - - -
        #
        qgram_list1 = [qgram_str1[i:i+q] for i in range(len(qgram_str1)-(q-1))]
        qgram_list2 = [qgram_str2[i:i+q] for i in range(len(qgram_str2)-(q-1))]

        # Get common q-grams  - - - - - - - - - - - - - - - - - - - - - - - - -
        #
        common = 0

        if (num_qgram1 < num_qgram2):  # Count using the shorter q-gram list
          short_qgram_list = qgram_list1
          long_qgram_list =  qgram_list2
        else:
          short_qgram_list = qgram_list2
          long_qgram_list =  qgram_list1

        for q_gram in short_qgram_list:
          if (q_gram in long_qgram_list):
            common += 1
            long_qgram_list.remove(q_gram)  # Remove the counted q-gram

        w = float(common) / float(divisor)

        assert (w >= 0.0), 'Q-gram: Similarity weight < 0.0'
        assert (w <= 1.0), 'Q-gram: Similarity weight > 1.0'

        w = self.__calc_partagree_weight__(val1, val2, w)

    if (self.do_caching == True):  # Put values pair into the cache
      self.__put_into_cache__(val1, val2, w)

    return w

# =============================================================================

class FieldComparatorPosQGram(FieldComparatorApproxString):
  """A field comparator based on the positional q-gram approximate comparator.

     q-grams are the sub-strings containing q characters in a string. For
     example, 'peter' contains the bigrams (q=2): ['pe','et','te','er'].

     Positional q-grams also contain the position within the string:
     [('pe',0),('et',1),('te',2),('er',3)].

     This method counts the number of common q-grams within a maximum distance
     and divides by either the minimum, average, or maximum number of q-grams.
     The resulting similarity weight is returned.

     The additional arguments (besides the base class arguments) which have to
     be set when this field comparator is initialised are:

       q               The length of the q-grams to be used (must be at least
                       1). The default value is 2 (i.e. bigrams)
       max_dist        Maximum distance allowed between two positional q-grams
                       (for example, with max_dist = 2 ('pe',6) and ('pe',8)
                       will be compared, however, ('pe',1) and ('pe',7) will
                       not be compared).
       common_divisor  Method of how to calculate the divisor. Can be set to
                       'average','shortest', or 'longest', and is calculated
                       according to the lengths of the two input strings.
       padded          If set to True (default), the beginnng and end of the
                       strings will be padded with (q-1) special characters, if
                       False no padding will be done.

     Padding will result in specific q-grams at the beginning and end of a
     string, for example 'peter' converted into padded bigrams (q=2) will
     result in the following 2-gram list:
     [('*p',0),('pe',1),('et',2),('te',3),('er',4),('r@',5)], with '*'
     illustrating the start and '@' the end character.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'q', 'max_dist', 'common_divisor' and 'padded'
       arguments first, then call the base class constructor.
    """

    self.q =              2
    self.max_dist =       None
    self.common_divisor = None
    self.padded =         True

    # Process all keyword arguments - - - - - - - - - - - - - - - - - - - - - -
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword == 'q'):
        auxiliary.check_is_positive('q', value)
        self.q = value

      elif (keyword.startswith('max_di')):
        auxiliary.check_is_positive('max_dist', value)
        self.max_dist = value

      elif (keyword.startswith('common_d')):
        auxiliary.check_is_string('common_divisor', value)
        if (value not in ['average','shortest','longest']):
          logging.exception('Value of argument "common_divisor" is not one' + \
                            'of: "average", "shortest", or "longest": %s' % \
                            (value))
          raise Exception
        self.common_divisor = value

      elif (keyword.startswith('pad')):
        auxiliary.check_is_flag('padded', value)
        self.padded = value

      else:
        base_kwargs[keyword] = value

    FieldComparatorApproxString.__init__(self, base_kwargs)

    # Make sure 'common_divisor' and 'max_dist' attributes are set  - - - - - -
    #
    auxiliary.check_is_string('common_divisor', self.common_divisor)
    auxiliary.check_is_positive('max_dist', self.max_dist)

    self.log([('Threshold', self.threshold),
              ('Maximum distance tolerated', self.max_dist),
              ('q', self.q),
              ('Common divisor', self.common_divisor),
              ('Padded flag', self.padded)])  # Log a message

    self.QGRAM_START_CHAR = chr(1)
    self.QGRAM_END_CHAR =   chr(2)

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values using the q-gram approximate string comparator.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    if (self.do_caching == True):  # Check if values pair is in the cache
      cache_weight = self.__get_from_cache__(val1, val2)
      if (cache_weight != None):
        return cache_weight

    if (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # Calculate positional q-gram similarity value - - - - - - - - - - - - - -
    #
    q = self.q  # Faster access

    # Calculate number of q-grams in strings (plus start and end characters)
    #
    if (self.padded == True):
      num_qgram1 = len(val1)+q-1
      num_qgram2 = len(val2)+q-1
    else:
      num_qgram1 = max(len(val1)-(q-1), 0)  # Make sure its not negative
      num_qgram2 = max(len(val2)-(q-1), 0)

    # Check if there are q-grams at all from both strings - - - - - - - - - - -
    # (no q-grams if length of a string is less than q)
    #
    if ((self.padded == False) and (min(num_qgram1, num_qgram2) == 0)):
      w = self.disagree_weight

    else:
      if (self.common_divisor == 'average'):  # Calculate the divisor
        divisor = 0.5*(num_qgram1+num_qgram2)
      elif (self.common_divisor == 'shortest'):
        divisor = min(num_qgram1, num_qgram2)
      else:  # Longest
        divisor = max(num_qgram1, num_qgram2)

      # Use number of q-grams to quickly check if below threshold - - - - - - -
      #
      max_common_qgram = min(num_qgram1, num_qgram2)  # Max possible q-grams
      w = float(max_common_qgram) / float(divisor)    # ... in common

      if (w  < self.threshold):  # Similariy is smaller than threshold
        w = self.disagree_weight

      else:

        # Add start and end characters (padding) - - - - - - - - - - - - - - -
        #
        if (self.padded == True):
          qgram_str1 = (q-1)*self.QGRAM_START_CHAR+val1+(q-1)* \
                       self.QGRAM_END_CHAR
          qgram_str2 = (q-1)*self.QGRAM_START_CHAR+val2+(q-1)* \
                       self.QGRAM_END_CHAR
        else:
          qgram_str1 = val1
          qgram_str2 = val2

        # Make a list of q-grams for both strings - - - - - - - - - - - - - - -
        #
        qgram_list1 = \
                  [(qgram_str1[i:i+q],i) for i in range(len(qgram_str1)-(q-1))]
        qgram_list2 = \
                  [(qgram_str2[i:i+q],i) for i in range(len(qgram_str2)-(q-1))]

        # Get common q-grams  - - - - - - - - - - - - - - - - - - - - - - - - -
        #
        common = 0

        if (num_qgram1 < num_qgram2):  # Count using the shorter q-gram list
          short_qgram_list = qgram_list1
          long_qgram_list =  qgram_list2
        else:
          short_qgram_list = qgram_list2
          long_qgram_list =  qgram_list1

        for pos_q_gram in short_qgram_list:
          (q_gram,pos) = pos_q_gram

          pos_range = range(max(pos-self.max_dist,0), pos+self.max_dist+1)

          for test_pos in pos_range:
            test_pos_q_gram = (q_gram,test_pos)
            if (test_pos_q_gram in long_qgram_list):
              common += 1
              long_qgram_list.remove(test_pos_q_gram) # Remove counted q-gram
              break

        w = float(common) / float(divisor)

        assert (w >= 0.0), 'Positional Q-gram: Similarity weight < 0.0'
        assert (w <= 1.0), 'Positional Q-gram: Similarity weight > 1.0'

        w = self.__calc_partagree_weight__(val1, val2, w)

    if (self.do_caching == True):  # Put values pair into the cache
      self.__put_into_cache__(val1, val2, w)

    return w

# =============================================================================

class FieldComparatorSGram(FieldComparatorApproxString):
  """A field comparator based on the skip-grams approximate string comparator.

     S-grams or skip-grams are based on bigrams (q=2) which can contain gaps.
     For more details please refer to:

     "Non-adjacent Digrams Improve Matching of Cross-Lingual Spelling Variants"
     by H. Keskustalo, A. Pirkola, K. Visala, E. Leppanen and J. Jarvelin,
     SPIRE 2003.

     This method counts the number of common s-grams and divides by either the
     minimum, average, or maximum number of q-grams. The resulting similarity
     weight is returned.

     The additional arguments (besides the base class arguments) which have to
     be set when this field comparator is initialised are:

       gram_class_list  Gram class list (see the mentioned publication above
                        for more details).
       common_divisor   Method of how to calculate the divisor. Can be set to
                        'average','shortest', or 'longest', and is calculated
                        according to the lengths of the two input strings.
       padded           If set to True (default), the beginning and end of the
                        strings will be padded with (q-1) special characters,
                        if False no padding will be done.

     Padding will result in special start and end characters being added at the
     beginning and the end of the character, similar as done with the q-gram
     comparator.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'gram_class_list', 'common_divisor' and
       'padded' arguments first, then call the base class constructor.
    """

    self.gram_class_list = None
    self.common_divisor =  None
    self.padded =          True

    # Process all keyword arguments - - - - - - - - - - - - - - - - - - - - - -
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('gram')):
        auxiliary.check_is_list('gram_class_list', value)
        self.gram_class_list = value

      elif (keyword.startswith('common_d')):
        auxiliary.check_is_string('common_divisor', value)
        if (value not in ['average','shortest','longest']):
          logging.exception('Value of argument "common_divisor" is not one' + \
                            'of: "average", "shortest", or "longest": %s' % \
                            (value))
          raise Exception
        self.common_divisor = value

      elif (keyword.startswith('pad')):
        auxiliary.check_is_flag('padded', value)
        self.padded = value

      else:
        base_kwargs[keyword] = value

    FieldComparatorApproxString.__init__(self, base_kwargs)

    # Make sure 'common_divisor' and 'gram_class_list' attributes are set  - -
    #
    auxiliary.check_is_string('common_divisor', self.common_divisor)
    auxiliary.check_is_list('gram_class_list', self.gram_class_list)

    self.log([('Threshold', self.threshold),
              ('Gram class list', self.gram_class_list),
              ('Common divisor', self.common_divisor),
              ('Padded flag', self.padded)])  # Log a message

    self.SGRAM_START_CHAR = chr(1)
    self.SGRAM_END_CHAR =   chr(2)

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values using the s-gram approximate string comparator.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    if (self.do_caching == True):  # Check if values pair is in the cache
      cache_weight = self.__get_from_cache__(val1, val2)
      if (cache_weight != None):
        return cache_weight

    if (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # Calculate s-gram similarity value - - - - - - - - - - - - - - - - - - - -

    # Extend strings with start and end characters
    #
    if (self.padded == True):
      tmp_str1 = self.SGRAM_START_CHAR+val1+self.SGRAM_END_CHAR
      tmp_str2 = self.SGRAM_START_CHAR+val2+self.SGRAM_END_CHAR
    else:
      tmp_str1 = val1
      tmp_str2 = val2

    len1 = len(tmp_str1)
    len2 = len(tmp_str2)

    common = 0.0   # Sum number of common s-grams over gram classes
    divisor = 0.0  # Sum of divisors over gram classes

    for c in self.gram_class_list:  # Loop over all gram classes given - - - -

      sgram_list1 = []
      sgram_list2 = []

      for s in c:  # Skip distances
        for i in range(0,len1-s-1):
          sgram_list1.append(tmp_str1[i]+tmp_str1[i+s+1])
        for i in range(0,len2-s-1):
          sgram_list2.append(tmp_str2[i]+tmp_str2[i+s+1])

      num_sgram1 = len(sgram_list1)
      num_sgram2 = len(sgram_list2)

      if (self.common_divisor == 'average'):
        this_divisor = 0.5*(num_sgram1+num_sgram2)  # Average num of s-grams
      elif (self.common_divisor == 'shortest'):
        this_divisor = min(num_sgram1,num_sgram2)
      else:  # Longest
        this_divisor = max(num_sgram1,num_sgram2)

      if (num_sgram1 < num_sgram2):  # Count using the shorter s-gram list
        short_sgram_list = sgram_list1
        long_sgram_list =  sgram_list2
      else:
        short_sgram_list = sgram_list2
        long_sgram_list =  sgram_list1

      this_common = 0  # Number of common s-grams for this gram class

      for s_gram in short_sgram_list:
        if (s_gram in long_sgram_list):
          this_common += 1
          long_sgram_list.remove(s_gram)  # Remove the counted s-gram

      common +=  this_common
      divisor += this_divisor

    if (divisor == 0):  # One string did not have any s-gram
      w = 0.0
    else:
      w = common / divisor

    assert (w >= 0.0), 'S-gram: Similarity weight < 0.0'
    assert (w <= 1.0), 'S-gram: Similarity weight > 1.0'

    w = self.__calc_partagree_weight__(val1, val2, w)

    if (self.do_caching == True):  # Put values pair into the cache
      self.__put_into_cache__(val1, val2, w)

    return w

# =============================================================================

class FieldComparatorEditDist(FieldComparatorApproxString):
  """A field comparator based on the edit distance (or Levenshtein) approximate
     string comparator.

     The edit distance is the minimal number of insertions, deletions and
     substitutions needed to transform one string into the other.

     For more information on edit distance see for example:

       http://www.nist.gov/dads/HTML/editdistance.html
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """No additional attributes needed, just call base class constructor.
    """

    FieldComparatorApproxString.__init__(self, kwargs)

    self.log([('Threshold', self.threshold)])  # Log a message

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values using the edit-distance (or Levenshtein)
       approximate string comparator.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    if (self.do_caching == True):  # Check if values pair is in the cache
      cache_weight = self.__get_from_cache__(val1, val2)
      if (cache_weight != None):
        return cache_weight

    if (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # Calculate edit distance similarity value - - - - - - - - - - - - - - - -
    #
    n = len(val1)
    m = len(val2)
    max_len = max(n,m)

    # Quick check if edit distance is below threshold - - - - - - - - - - - - -
    #
    len_diff = abs(n-m)
    w = 1.0 - float(len_diff) / float(max_len)

    if (w  < self.threshold):  # Similariy is smaller than threshold
      w = self.disagree_weight

    else: # Calculate the maximum distance possible with this threshold
      max_dist = (1.0-self.threshold)*max_len

      if (n > m):  # Make sure n <= m, to use O(min(n,m)) space
        str1 = val2
        str2 = val1
        n, m = m, n
      else:
        str1 = val1
        str2 = val2

      current = range(n+1)

      w = -1  # Set weight to an illegal value (so it can be chacked later)

      for i in range(1, m+1):
        previous = current
        current =  [i]+n*[0]
        str2char = str2[i-1]

        for j in range(1,n+1):
          substitute = previous[j-1]
          if (str1[j-1] != str2char):
            substitute += 1

          # Get minimum of insert, delete and substitute
          #
          current[j] = min(previous[j]+1, current[j-1]+1, substitute)

        if (min(current) > max_dist):  # Distance is already too large
          w = max(1.0 - float(max_dist+1) / float(max_len), 0.0)
          break  # Exit loop

      if (w == -1):  # Weight has not been calculated
        w = 1.0 - float(current[n]) / float(max_len)

      assert (w >= 0.0), 'Edit distance: Similarity weight < 0.0'
      assert (w <= 1.0), 'Edit distance: Similarity weight > 1.0'

      w = self.__calc_partagree_weight__(val1, val2, w)

    if (self.do_caching == True):  # Put values pair into the cache
      self.__put_into_cache__(val1, val2, w)

    return w

# =============================================================================

class FieldComparatorDaLeDist(FieldComparatorApproxString):
  """A field comparator based on a modified edit (or Levenshtein) distance that
     counts transpositions as elementary operations as well. This is also
     called the Damerau-Levenshtein distance.

     The DaLe distance is the minimal number of insertions, deletions,
     substitutions and transpositions needed to transform one string into the
     other.

     Compared to the original edit distance function, which handles a
     transposition (like: 'sydney' <-> 'sydeny' as 2 operations (two
     substitutions or one insert and one delete), this modified version handles
     this as 1 operation.

     Based on code from Justin Zobel's 'vrank'.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """No additional attributes needed, just call base class constructor.
    """

    FieldComparatorApproxString.__init__(self, kwargs)

    self.log([('Threshold', self.threshold)])  # Log a message

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values using the Damerau-Levenshtein distance
       approximate string comparator.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    if (self.do_caching == True):  # Check if values pair is in the cache
      cache_weight = self.__get_from_cache__(val1, val2)
      if (cache_weight != None):
        return cache_weight

    if (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # Calculate Damerau-Levenshtein distance similarity value - - - - - - - - -
    #
    n = len(val1)
    m = len(val2)
    max_len = max(n,m)

    # Quick check if Damerau-Levenshtein is below threshold - - - - - - - - - -
    #
    len_diff = abs(n-m)
    w = 1.0 - float(len_diff) / float(max_len)

    if (w  < self.threshold):  # Similariy is smaller than threshold
      w = self.disagree_weight

    else: # Calculate the maximum distance possible with this threshold
      max_dist = (1.0-self.threshold)*max_len

      if (n > m):  # Make sure n <= m, to use O(min(n,m)) space
        str1 = val2
        str2 = val1
        n, m = m, n
      else:
        str1 = val1
        str2 = val2

      d = []  # Table with the full distance matrix

      current = range(n+1)
      d.append(current)

      w = -1  # Set weight to an illegal value (so it can be chacked later)

      for i in range(1,m+1):

        previous = current
        current =  [i]+n*[0]
        str2char = str2[i-1]

        for j in range(1,n+1):
          substitute = previous[j-1]
          if (str1[j-1] != str2char):
            substitute += 1

          if (i == 1) or (j == 1):  # First characters, no transp possible

            # Get minimum of insert, delete and substitute
            #
            current[j] = min(previous[j]+1, current[j-1]+1, substitute)

          else:
            if (str1[j-2] == str2[i-1]) and (str1[j-1] == str2[i-2]):
              transpose = d[i-2][j-2] + 1
            else:
              transpose = d[i-2][j-2] + 3

            current[j] = min(previous[j]+1, current[j-1]+1, substitute, \
                             transpose)

        d.append(current)

        if (min(current) > max_dist):  # Distance is already too large
          w = max(1.0 - float(max_dist+1) / float(max_len), 0.0)
          break  # Exit loop

      if (w == -1):  # Weight has not been calculated

        w = 1.0 - float(current[n]) / float(max_len)

      assert (w >= 0.0), 'DaLe distance: Similarity weight < 0.0'
      assert (w <= 1.0), 'DaLe distance: Similarity weight > 1.0'

      w = self.__calc_partagree_weight__(val1, val2, w)

    if (self.do_caching == True):  # Put values pair into the cache
      self.__put_into_cache__(val1, val2, w)

    return w

# =============================================================================

class FieldComparatorBagDist(FieldComparatorApproxString):
  """A field comparator based on the bag distance approximate string
     comparator.

     Bag distance is a cheap method to calculate the distance between two
     strings. It is always smaller or equal to the edit distance, and therefore
     the similarity measure returned by the method is always larger than the
     edit distance similarity measure.

     For more details see for example:

     "String Matching with Metric Trees Using an Approximate Distance"
     Ilaria Bartolini, Paolo Ciaccia and Marco Patella,
     in Proceedings of the 9th International Symposium on String Processing
     and Information Retrieval, Lisbone, Purtugal, September 2002.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """No additional attributes needed, just call base class constructor.
    """

    FieldComparatorApproxString.__init__(self, kwargs)

    self.log([('Threshold', self.threshold)])  # Log a message

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values using the bag distance approximate string
       comparator.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    if (self.do_caching == True):  # Check if values pair is in the cache
      cache_weight = self.__get_from_cache__(val1, val2)
      if (cache_weight != None):
        return cache_weight

    if (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # Calculate bag distance similarity value - - - - - - - - - - - - - - - - -
    #
    n = len(val1)
    m = len(val2)

    list1 = list(val1)
    list2 = list(val2)

    for ch in val1:
      if (ch in list2):
        list2.remove(ch)

    for ch in val2:
      if (ch in list1):
        list1.remove(ch)

    b = max(len(list1),len(list2))

    w = 1.0 - float(b) / float(max(n,m))

    assert (w >= 0.0), 'Bag distance: Similarity weight < 0.0'
    assert (w <= 1.0), 'Bag distance: Similarity weight > 1.0'

    w = self.__calc_partagree_weight__(val1, val2, w)

    if (self.do_caching == True):  # Put values pair into the cache
      self.__put_into_cache__(val1, val2, w)

    return w

# =============================================================================

class FieldComparatorSWDist(FieldComparatorApproxString):
  """A field comparator based on the Smith-Waterman distance approximate string
     comparator.

     Smith-Waterman distance is commonly used in biological sequence alignment.

     Scores for matches, misses, gap and extension penalties are set to values
     described in:

     "The field matching problem: Algorithms and applications"
     by A.E. Monge and C.P. Elkan, 1996.

     The additional argument (besides the base class arguments) which has to be
     set when this field comparator is initialised is:

       common_divisor  Method of how to calculate the divisor. Can be set to
                       'average','shortest', or 'longest', and is calculated
                       according to the lengths of the two input strings.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'common_divisor' argument first, then call the
       base class constructor.
    """

    self.common_divisor =  None

    # Process all keyword arguments - - - - - - - - - - - - - - - - - - - - - -
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('common_d')):
        auxiliary.check_is_string('common_divisor', value)
        if (value not in ['average','shortest','longest']):
          logging.exception('Value of argument "common_divisor" is not one' + \
                            'of: "average", "shortest", or "longest": %s' % \
                            (value))
          raise Exception
        self.common_divisor = value

      else:
        base_kwargs[keyword] = value

    FieldComparatorApproxString.__init__(self, base_kwargs)

    # Make sure 'common_divisor' attribute is set  - - - - - - - - - - - - - -
    #
    auxiliary.check_is_string('common_divisor', self.common_divisor)

    self.log([('Threshold', self.threshold),
              ('Common divisor', self.common_divisor)])  # Log a message

    # Scores used for Smith-Waterman algorithm - - - - - - - - - - - - - - - -
    #
    self.match_score =       5
    self.approx_score =      2
    self.mismatch_score =   -5
    self.gap_penalty =       5
    self.extension_penalty = 1

    # Dictionary with approximate match characters mapped into numbers
    # {a,e,i,o,u} -> 0, {d,t} -> 1, {g,j} -> 2, {l,r} -> 3, {m,n} -> 4,
    # {b,p,v} -> 5
    #
    self.approx_matches = {'a':0,'b':5,'d':1,'e':0,'g':2,'i':0,'j':2,'l':3,
                           'm':4,'n':4,'o':0,'p':5,'r':3,'t':1,'u':0,'v':5}

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values using the Smith-Waterman distance approximate
       string comparator.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    if (self.do_caching == True):  # Check if values pair is in the cache
      cache_weight = self.__get_from_cache__(val1, val2)
      if (cache_weight != None):
        return cache_weight

    if (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # Calculate Smith-Waterman distance similarity value - - - - - - - - - - -
    #
    n = len(val1)
    m = len(val2)

    if (self.common_divisor == 'average'):
      divisor = 0.5*(n+m)*self.match_score  # Average maximum score
    elif (self.common_divisor == 'shortest'):
      divisor = min(n,m)*self.match_score
    else:  # Longest
      divisor = max(n,m)*self.match_score

    best_score = 0  # Keep the best score while calculating table

    d = []  # Table with the full distance matrix

    for i in range(n+1):  # Initalise table
      d.append([0.0]*(m+1))

    for i in range(1,n+1):
      vali1 = val1[i-1]
      approx_match1 = self.approx_matches.get(vali1,-1)

      for j in range(1,m+1):
        valj2 = val2[j-1]

        match = d[i-1][j-1]

        if (vali1 == valj2):
          match += self.match_score
        else:
          approx_match2 = self.approx_matches.get(valj2,-1)

          if (approx_match1 >= 0) and (approx_match2 >= 0) and \
             (approx_match1 == approx_match2):
            match += self.approx_score
          else:
            match += self.mismatch_score

        insert = 0
        for k in range(1,i):
          score = d[i-k][j] - self.gap_penalty - k*self.extension_penalty
          insert = max(insert, score)

        delete = 0
        for l in range(1,j):
          score = d[i][j-l] - self.gap_penalty - l*self.extension_penalty
          delete = max(delete, score)

        d[i][j] = max(match, insert, delete, 0)
        best_score = max(d[i][j], best_score)

    # best_score can be min(len(str1),len)str2))*match_score (if one string is
    # a sub-string of the other string)
    #
    # The lower best_score the less similar the sequences are.
    #
    w = float(best_score) / float(divisor)

    assert (w >= 0.0), 'Smith-Waterman distance: Similarity weight < 0.0'
    assert (w <= 1.0), 'Smith-Waterman distance: Similarity weight > 1.0'

    w = self.__calc_partagree_weight__(val1, val2, w)

    if (self.do_caching == True):  # Put values pair into the cache
      self.__put_into_cache__(val1, val2, w)

    return w

# =============================================================================

class FieldComparatorSyllAlDist(FieldComparatorApproxString):
  """A field comparator based on the syllable alignment distance approximate
     string comparator.

     The syllable alignment distance is based on syllables instead of
     characters and calculates a distance similar to edit distance.

     For more information see:
     "Syllable Alignment: A Novel Approach for Phonetic String Search"
     by Ruibin Gong and Tony k.Y. Chan, IEICE, 2006.

     The additional arguments (besides the base class arguments) which have to
     be set when this field comparator is initialised are:

       common_divisor  Method of how to calculate the divisor. Can be set to
                       'average','shortest', or 'longest', and is calculated
                       according to the number of syllables of the two input
                       strings.
       do_phonix       A flag, if set to True (the default) the Phonix sound
                       encoding transformation will be applied first to both
                       strings, otherwise the original strings will be used.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'common_divisor' and 'do_phonix' argument
       first, then call the base class constructor.
    """

    self.common_divisor = None
    self.do_phonix =      True

    # Process all keyword arguments - - - - - - - - - - - - - - - - - - - - - -
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('common_d')):
        auxiliary.check_is_string('common_divisor', value)
        if (value not in ['average','shortest','longest']):
          logging.exception('Value of argument "common_divisor" is not one' + \
                            'of: "average", "shortest", or "longest": %s' % \
                            (value))
          raise Exception
        self.common_divisor = value

      elif (keyword.startswith('do_p')):
        auxiliary.check_is_flag('do_phonix', value)
        self.do_phonix = value

      else:
        base_kwargs[keyword] = value

    FieldComparatorApproxString.__init__(self, base_kwargs)

    # Make sure 'common_divisor' attribute is set  - - - - - - - - - - - - - -
    #
    auxiliary.check_is_string('common_divisor', self.common_divisor)

    self.log([('Threshold', self.threshold),
              ('Do phonix transformation flag', self.do_phonix),
              ('Common divisor', self.common_divisor)])  # Log a message

    # Substitution and gap penalty weights
    #
    self.s1 = 1   # Aligning two chars (not syllable start) that are the same
    self.s2 = -1  # Aligning two chars (not syllable start) that are different
    self.s3 = -4  # Aligning a character with a syllable start
    self.s4 = 6   # Aligning two syllable starts that are the same
    self.s5 = -2  # Aligning two syllable starts that are different
    self.g1 = -1  # Aligning a gap with a character (not syllable start)
    self.g2 = -3  # Aligning a gap with a syllable start

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values using the syllable alignment distance
       approximate string comparator.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    if (self.do_caching == True):  # Check if values pair is in the cache
      cache_weight = self.__get_from_cache__(val1, val2)
      if (cache_weight != None):
        return cache_weight

    if (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # Calculate syllable alignment distance similarity value - - - - - - - - -
    #
    if (self.do_phonix == True):
      workstr1 = encode.phonix_transform(val1)
      workstr2 = encode.phonix_transform(val2)
    else:
      workstr1 = val1
      workstr2 = val2

    # Syllable scan, make beginning of each syllable an uppercase character
    #
    syll_str_list = []  # List for the two syllable strings

    for s in (workstr1, workstr2):
      str_list = list(s)
      str_list[0] = str_list[0].upper() # First char is start of 1st syllable
      str_len = len(s)

      for i in range(1, str_len):

        if (str_list[i] not in 'aeiouyAEIOUY'):

          if (i < (str_len-1)):  # Not last character
            if (str_list[i+1] in 'aeiouyAEIOUYhrw'):
              str_list[i] = str_list[i].upper()

          elif (str_list[i] not in 'aeiouyAEIOUY'):
            str_list[i] = str_list[i].upper()

          if (str_list[i] in 'HRW') and (str_list[i-1] <= 'Z'):
            str_list[i] = str_list[i].lower()

      syll_str_list.append(''.join(str_list))  # Convert back to string

    wstr1 = syll_str_list[0]
    wstr2 = syll_str_list[1]

    n, m = len(wstr1), len(wstr2)

    # Calculate maximum number of syllable starts and other characters to get
    # maximum possible alignment weight
    #
    max_w1 = 0
    for c in wstr1:
      if c.isupper():
        max_w1 += self.s4  # Syllable start
      else:
        max_w1 += self.s1  # Other characters

    max_w2 = 0
    for c in wstr2:
      if c.isupper():
        max_w2 += self.s4  # Syllable start
      else:
        max_w2 += self.s1  # Other characters

    # Calculate the divisor - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.common_divisor == 'average'):
      divisor = 0.5*(max_w1+max_w2)  # Average weight
    elif (self.common_divisor == 'shortest'):
      divisor = min(max_w1,max_w2)
    else:  # Longest
      divisor = max(max_w1,max_w2)

    d = []  # Table with the full distance matrix

    for i in range(n+1):  # Initalise table
      d.append([0.0]*(m+1))

    for i in range(1,n+1):  # First column
      if (wstr1[i-1].isupper()):
        d[i][0] = d[i-1][0]+self.g2
      else:
        d[i][0] = d[i-1][0]+self.g1

    for j in range(1,m+1):  # First row
      if (wstr2[j-1].isupper()):
        d[0][j] = d[0][j-1]+self.g2
      else:
        d[0][j] = d[0][j-1]+self.g1

    for j in range(1,m+1):  # Fill in rest of table
      c2 = wstr2[j-1]

      for i in range(1,n+1):
        c1 = wstr1[i-1]

        if (c1.isupper()):
          x = d[i-1][j]+self.g2
        else:
          x = d[i-1][j]+self.g1

        if (c2.isupper()):
          y = d[i][j-1]+self.g2
        else:
          y = d[i][j-1]+self.g1

        if (c1.isupper() and c2.isupper()):
          if (c1 == c2):
            z = d[i-1][j-1]+self.s4
          else:
            z = d[i-1][j-1]+self.s5
        elif (c1.islower() and c2.islower()):
          if (c1 == c2):
            z = d[i-1][j-1]+self.s1
          else:
            z = d[i-1][j-1]+self.s2
        else:
          z = d[i-1][j-1]+self.s3

        d[i][j] = max(x,y,z)

    w = max(float(d[i][j]) / float(divisor), 0.0)

    assert (w >= 0.0), 'Syllable-alignment distance: Similarity weight < 0.0'
    assert (w <= 1.0), 'Syllable-alignment distance: Similarity weight > 1.0'

    w = self.__calc_partagree_weight__(val1, val2, w)

    if (self.do_caching == True):  # Put values pair into the cache
      self.__put_into_cache__(val1, val2, w)

    return w

# =============================================================================

class FieldComparatorSeqMatch(FieldComparatorApproxString):
  """A field comparator based on the Python standard library 'difflib' sequence
     matcher.

     For more details see the Python module 'difflib' at:

       http://www.python.org/doc/current/lib/module-difflib.html

     Because the matches are not commutative, the pair and the swapped pair are
     compared and the average is taken.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """No additional attributes needed, just call base class constructor.
    """

    FieldComparatorApproxString.__init__(self, kwargs)

    self.log([('Threshold', self.threshold)])  # Log a message

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values using the Python standard library 'difflib'
       sequence matcher approximate string comparator.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    if (self.do_caching == True):  # Check if values pair is in the cache
      cache_weight = self.__get_from_cache__(val1, val2)
      if (cache_weight != None):
        return cache_weight

    if (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # Calculate sequence matcher similarity value - - - - - - - - - - - - - - -
    #
    seq_matcher_1 = difflib.SequenceMatcher(None, val1, val2)
    seq_matcher_2 = difflib.SequenceMatcher(None, val2, val1)

    w = (seq_matcher_1.ratio()+seq_matcher_2.ratio()) / 2.0 # Calc average

    assert (w >= 0.0), 'Python sequence matcher: Similarity weight < 0.0'
    assert (w <= 1.0), 'Python sequence matcher: Similarity weight > 1.0'

    w = self.__calc_partagree_weight__(val1, val2, w)

    if (self.do_caching == True):  # Put values pair into the cache
      self.__put_into_cache__(val1, val2, w)

    return w

# =============================================================================

class FieldComparatorEditex(FieldComparatorApproxString):
  """A field comparator based on the editex approximate string comparator.

     Based on ideas described in:

     "Phonetic String Matching: Lessons Learned from Information Retrieval"
     by Justin Zobel and Philip Dart, SIGIR 1995.

     Important: This function assumes that the input strings only contain
     letters and whitespace, but no other characters.

     All non-letter characters are handled like whitespaces (like a silent
     sound).
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """No additional attributes needed, just call base class constructor.
    """

    FieldComparatorApproxString.__init__(self, kwargs)

    self.log([('Threshold', self.threshold)])  # Log a message

    # Values for edit costs - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    self.BIG_COSTS = 3  # If characters are not in same group
    self.SML_COSTS = 2  # If characters are in same group

    # Mappings of letters into groups - - - - - - - - - - - - - - - - - - - - -
    #
    self.groupsof_dict = {'a':0,'b':1,'c':2,'d':3,'e':0,'f':1,'g':2,'h':7,
                          'i':0,'j':2,'k':2,'l':4,'m':5,'n':5,'o':0,'p':1,
                          'q':2,'r':6,'s':2,'t':3,'u':0,'v':1,'w':7,'x':2,
                          'y':0,'z':2,'{':7}

  # ---------------------------------------------------------------------------

  def __delete_cost__(self, char1, char2):

    if (char1 == char2):
      return 0

    code1 = self.groupsof_dict.get(char1,-1)  # -1 is not a char
    code2 = self.groupsof_dict.get(char2,-2)  # -2 if not a char

    if (code1 == code2) or (code2 == 7):  # Same or silent
      return self.SML_COSTS  # Small difference costs
    else:
      return self.BIG_COSTS

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values using the editex approximate string comparator.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    if (self.do_caching == True):  # Check if values pair is in the cache
      cache_weight = self.__get_from_cache__(val1, val2)
      if (cache_weight != None):
        return cache_weight

    if (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # Calculate editex similarity value - - - - - - - - - - - - - - - - - - - -
    #
    n, m = len(val1), len(val2)

    if (n > m):  # Make sure n <= m, to use O(min(n,m)) space
      str1 = val2.lower()
      str2 = val1.lower()
      n, m = m, n
    else:
      str1 = val1.lower()
      str2 = val2.lower()

    if (' ' in str1):
      str1 = str1.replace(' ','{')
    if (' ' in str2):
      str2 = str2.replace(' ','{')

    row = [0]*(m+1)  # Generate empty cost matrix
    F = []
    for i in range(n+1):
      F.append(row[:])

    F[1][0] = self.BIG_COSTS  # Initialise first row and first column of
    F[0][1] = self.BIG_COSTS  #   cost matrix

    sum = self.BIG_COSTS
    for i in range(2,n+1):
      sum += self.__delete_cost__(str1[i-2], str1[i-1])
      F[i][0] = sum

    sum = self.BIG_COSTS
    for j in range(2,m+1):
      sum += self.__delete_cost__(str2[j-2], str2[j-1])
      F[0][j] = sum

    for i in range(1,n+1):

      if (i == 1):
        inc1 = self.BIG_COSTS
      else:
        inc1 = self.__delete_cost__(str1[i-2], str1[i-1])

      for j in range(1,m+1):

        if (j == 1):
          inc2 = self.BIG_COSTS
        else:
          inc2 = self.__delete_cost__(str2[j-2], str2[j-1])

        if (str1[i-1] == str2[j-1]):
          diag = 0
        else:
          code1 = self.groupsof_dict.get(str1[i-1],-1)  # -1 is not a char
          code2 = self.groupsof_dict.get(str2[j-1],-2)  # -2 if not a char

          if (code1 == code2):  # Same phonetic group
            diag = self.SML_COSTS
          else:
            diag = self.BIG_COSTS

        F[i][j] = min(F[i-1][j]+inc1, F[i][j-1]+inc2, F[i-1][j-1]+diag)

    w = 1.0 - float(F[n][m]) / float(max(F[0][m],F[n][0]))

    if (w < 0.0):
      w = 0.0

    assert (w >= 0.0), 'Editex: Similarity weight < 0.0'
    assert (w <= 1.0), 'Editex: Similarity weight > 1.0'

    w = self.__calc_partagree_weight__(val1, val2, w)

    if (self.do_caching == True):  # Put values pair into the cache
      self.__put_into_cache__(val1, val2, w)

    return w

# =============================================================================

class FieldComparatorLCS(FieldComparatorApproxString):
  """A field comparator based on the longest common substring approximate
     string comparator.

     Based on a dynamic programming algorithm, see for example:

       http://www.ics.uci.edu/~dan/class/161/notes/6/Dynamic.html
       http://www.unixuser.org/~euske/python/index.html
       http://en.wikipedia.org/wiki/Longest_common_substring_problem

     The algorithm extracts common substrings until no more are found with a
     minimum common length and then calculates a similairy measure.

     Note that the repeated lcs method is not symmetric, i.e. string pairs:
       'prap' / 'papr' -> 1.0  ('ap' is extracted first, leaving 'pr' / 'pr')
       'papr' / 'prap' -> 0.5  ('pr' is extracted first, leaving 'pa' / 'ap')
     (assuming minimum common length is set to 2). Therefore, lcs is run twice
     with input strings swapped and the similarity value averaged.

     The additional arguments (besides the base class arguments) which have to
     be set when this field comparator is initialised are:

       common_divisor  Method of how to calculate the divisor. Can be set to
                       'average','shortest', or 'longest', and is calculated
                       according to the number of syllables of the two input
                       strings.
       min_common_len  The minimum length of a common substring. Default value
                       is 2.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'common_divisor' and 'min_common_len' arguments
       first, then call the base class constructor.
    """

    self.common_divisor = None
    self.min_common_len = 2

    # Process all keyword arguments - - - - - - - - - - - - - - - - - - - - - -
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('common_d')):
        auxiliary.check_is_string('common_divisor', value)
        if (value not in ['average','shortest','longest']):
          logging.exception('Value of argument "common_divisor" is not one' + \
                            'of: "average", "shortest", or "longest": %s' % \
                            (value))
          raise Exception
        self.common_divisor = value

      elif (keyword.startswith('min_c')):
        auxiliary.check_is_positive('min_common_len', value)
        self.min_common_len = value

      else:
        base_kwargs[keyword] = value

    FieldComparatorApproxString.__init__(self, base_kwargs)

    # Make sure 'common_divisor' attribute is set  - - - - - - - - - - - - - -
    #
    auxiliary.check_is_string('common_divisor', self.common_divisor)

    self.log([('Threshold', self.threshold),
              ('Minimum common length', self.min_common_len),
              ('Common divisor', self.common_divisor)])  # Log a message

  # ---------------------------------------------------------------------------

  def __do_lcs__(self, str1, str2):
    """Method to extract longest common substring from the two input strings.
       Returns the common substring, its length, and the two input strings with
       the common substring removed.

       Should not be used from outside the module.
    """

    n = len(str1)
    m = len(str2)

    if (n > m):  # Make sure n <= m, to use O(min(n,m)) space
      str1, str2 = str2, str1
      n, m =       m, n
      swapped = True
    else:
      swapped = False

    current = (n+1)*[0]

    com_len = 0
    com_ans1 = -1
    com_ans2 = -1

    for i in range(m):
      previous = current
      current =  (n+1)*[0]

      for j in range(n):
        if (str1[j] != str2[i]):
          current[j] = 0
        else:
          current[j] = previous[j-1]+1
          if (current[j] > com_len):
            com_len = current[j]
            com_ans1 = j
            com_ans2 = i

    com1 = str1[com_ans1-com_len+1:com_ans1+1]
    com2 = str2[com_ans2-com_len+1:com_ans2+1]

    if (com1 != com2):
      logging.exception('LCS: Different common substrings: %s / %s in ' % \
                        (com1, com2) + 'original strings: %s / %s' % \
                        (str1, str2))
      raise Exception

    # Remove common substring from input strings
    #
    str1 = str1[:com_ans1-com_len+1] + str1[1+com_ans1:]
    str2 = str2[:com_ans2-com_len+1] + str2[1+com_ans2:]

    if (swapped == True):
      return com1, com_len, str2, str1
    else:
      return com1, com_len, str1, str2

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values using the longest common substring approximate
       string comparator.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    if (self.do_caching == True):  # Check if values pair is in the cache
      cache_weight = self.__get_from_cache__(val1, val2)
      if (cache_weight != None):
        return cache_weight

    if (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # Calculate longest common substring similarity value - - - - - - - - - - -
    #
    len1, len2 = len(val1), len(val2)

    if (self.common_divisor == 'average'):
      divisor = 0.5*(len1+len2)  # Compute average string length
    elif (self.common_divisor == 'shortest'):
      divisor = min(len1,len2)
    else:  # Longest
      divisor = max(len1,len2)

    # Quick check if below threshold - - - - - - - - - - - - - - - - - - - -
    #
    max_common_len = min(len1,len2)

    w = float(max_common_len) / float(divisor)

    if (w  < self.threshold):  # Similariy is smaller than threshold
      w = self.disagree_weight

    else:
      # Iterative calculation of longest common substring until strings to s
      #
      w = 0.0

      for (s1,s2) in [(val1,val2), (val2,val1)]:

        com_str, com_len, s1, s2 = self.__do_lcs__(s1, s2) # Find initial LCS

        total_com_str = com_str
        total_com_len = com_len

        while (com_len >= self.min_common_len):
          com_str, com_len, s1n, s2n = self.__do_lcs__(s1, s2)

          if (com_len >= self.min_common_len):
            total_com_str += com_str
            total_com_len += com_len
            s1, s2 = s1n, s2n

        w += float(total_com_len) / float(divisor)

      w /= 2.0

      assert (w >= 0.0), 'Longest common substring: Similarity weight < 0.0'
      assert (w <= 1.0), 'Longest common substring: Similarity weight > 1.0'

      w = self.__calc_partagree_weight__(val1, val2, w)

    if (self.do_caching == True):  # Put values pair into the cache
      self.__put_into_cache__(val1, val2, w)

    return w

# =============================================================================

class FieldComparatorOntoLCS(FieldComparatorApproxString):
  """A field comparator based on the ontology longest common substring
     approximate string comparator as described in:

       Giorgos Stoilos, Giorgos Stamou and Stefanos Kollinas:
       A String Metric for Ontology Alignment
       ISWC 2005, Springer LNCS 3729, pp 624-637, 2005.

     It starts calculating the basic longest common substrings, then calculates
     a difference measure using the unmatched characers, and finally applies
     the Winkler modification (increase wirght for same initial characters).

     The additional arguments (besides the base class arguments) which have to
     be set when this field comparator is initialised are:

       common_divisor  Method of how to calculate the divisor. Can be set to
                       'average','shortest', or 'longest', and is calculated
                       according to the number of syllables of the two input
                       strings.
       min_common_len  The minimum length of a common substring. Default value
                       is 2.
       p               Constant for Hamacher product difference, see above
                       mentioned paper, can be in [0,1]. Default value is 0.6
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'common_divisor' and 'min_common_len' arguments
       first, then call the base class constructor.
    """

    self.common_divisor = None
    self.min_common_len = 2
    self.p =              0.6

    # Process all keyword arguments - - - - - - - - - - - - - - - - - - - - - -
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('common_d')):
        auxiliary.check_is_string('common_divisor', value)
        if (value not in ['average','shortest','longest']):
          logging.exception('Value of argument "common_divisor" is not one' + \
                            'of: "average", "shortest", or "longest": %s' % \
                            (value))
          raise Exception
        self.common_divisor = value

      elif (keyword.startswith('min_c')):
        auxiliary.check_is_positive('min_common_len', value)
        self.min_common_len = value

      elif (keyword == 'p'):
        auxiliary.check_is_normalised('p', value)
        self.p = value

      else:
        base_kwargs[keyword] = value

    FieldComparatorApproxString.__init__(self, base_kwargs)

    # Make sure 'common_divisor' attribute is set  - - - - - - - - - - - - - -
    #
    auxiliary.check_is_string('common_divisor', self.common_divisor)

    self.log([('Threshold', self.threshold),
              ('Minimum common length', self.min_common_len),
              ('p (Hamacher difference constant)', self.p),
              ('Common divisor', self.common_divisor)])  # Log a message

  # ---------------------------------------------------------------------------

  def __do_lcs__(self, str1, str2):
    """Method to extract longest common substring from the two input strings.
       Returns the common substring, its length, and the two input strings with
       the common substring removed.

       Should not be used from outside the module.
    """

    n = len(str1)
    m = len(str2)

    if (n > m):  # Make sure n <= m, to use O(min(n,m)) space
      str1, str2 = str2, str1
      n, m =       m, n
      swapped = True
    else:
      swapped = False

    current = (n+1)*[0]

    com_len = 0
    com_ans1 = -1
    com_ans2 = -1

    for i in range(m):
      previous = current
      current =  (n+1)*[0]

      for j in range(n):
        if (str1[j] != str2[i]):
          current[j] = 0
        else:
          current[j] = previous[j-1]+1
          if (current[j] > com_len):
            com_len = current[j]
            com_ans1 = j
            com_ans2 = i

    com1 = str1[com_ans1-com_len+1:com_ans1+1]
    com2 = str2[com_ans2-com_len+1:com_ans2+1]

    if (com1 != com2):
      logging.exception('LCS: Different common substrings: %s / %s in ' % \
                        (com1, com2) + 'original strings: %s / %s' % \
                        (str1, str2))
      raise Exception

    # Remove common substring from input strings
    #
    str1 = str1[:com_ans1-com_len+1] + str1[1+com_ans1:]
    str2 = str2[:com_ans2-com_len+1] + str2[1+com_ans2:]

    if (swapped == True):
      return com1, com_len, str2, str1
    else:
      return com1, com_len, str1, str2

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values using the longest common substring approximate
       string comparator.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    if (self.do_caching == True):  # Check if values pair is in the cache
      cache_weight = self.__get_from_cache__(val1, val2)
      if (cache_weight != None):
        return cache_weight

    if (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # Calculate longest common substring similarity value - - - - - - - - - - -
    #
    len1, len2 = len(val1), len(val2)

    if (self.common_divisor == 'average'):
      divisor = 0.5*(len1+len2)  # Compute average string length
    elif (self.common_divisor == 'shortest'):
      divisor = min(len1,len2)
    else:  # Longest
      divisor = max(len1,len2)

    # Iterative calculation of longest common sub-strings - - - - - - - - - - -
    #
    w_lcs =  0.0  # Basic longest common sub-string weight
    h_diff = 0.0  # Hamacher product difference

    for (s1,s2) in [(val1,val2), (val2,val1)]:

      com_str, com_len, s1, s2 = self.__do_lcs__(s1, s2) # Find initial LCS

      total_com_str = com_str
      total_com_len = com_len

      while (com_len >= self.min_common_len):
        com_str, com_len, s1n, s2n = self.__do_lcs__(s1, s2)

        if (com_len >= self.min_common_len):
          total_com_str += com_str
          total_com_len += com_len
          s1, s2 = s1n, s2n

      w_lcs += float(total_com_len) / float(divisor)

      # Calculate Hamacher product difference for sub-strings left
      #
      s1_len = float(len(s1)) / len1
      s2_len = float(len(s2)) / len2

      h_diff += s1_len*s2_len / (self.p + (1-self.p) * \
                                          (s1_len + s2_len - s1_len*s2_len))

    w_lcs /=  2.0
    h_diff /= 2.0

    assert (w_lcs >= 0.0) and (w_lcs <= 1.0), \
           'Basic LCS similarity weight outside [0,1]: %f' % (w_lcs)
    assert (h_diff >= 0.0) and (h_diff <= 1.0), \
           'Hamacher product difference outside [0,1]: %f' % (h_diff)

    # Apply Winkler modification (check for same characters at the beginning) -
    #
    same_init = 0  # Number of same characters at beginning

    minlen = min(len1, len2, 4)

    for same_init in range(1, minlen+1):
      if (val1[:same_init] != val2[:same_init]):
        break
    same_init -= 1

    assert (same_init >= 0), 'Winkler: "same_init" value smaller than 0'
    assert (same_init <= 4), 'Winkler: "same_init" value larger than 4'

    w_lcs_wink = w_lcs + 0.1 * same_init * (1.0 - w_lcs)

    assert (w_lcs_wink >= 0.0) and (w_lcs_wink <= 1.0), \
           'Winkler modified LCS similarity weight outside [0,1]: %f' % \
           (w_lcs_wink)

    w = w_lcs_wink - h_diff  # A weight in interval [-1,1]

    w = w/2.0 + 0.5  # Scale into [0,1]

    assert (w >= 0.0) and (w <= 1.0), \
           'Ontology LCS similarity weight outside [0,1]: %f' % (w)

    w = self.__calc_partagree_weight__(val1, val2, w)

    if (self.do_caching == True):  # Put values pair into the cache
      self.__put_into_cache__(val1, val2, w)

    return w

# =============================================================================

class FieldComparatorCompress(FieldComparatorApproxString):
  """A field comparator based on the compression approximate string comparator.

     For more information about using compression for similarity measures see:

     - Cilibrasi, R. and Vitanyi, P.: Clustering by compression. IEEE Trans.
       Infomat. Th. Submitted, 2004. See: http://arxiv.org/abs/cs.CV/0312044

     - Keogh, E., Lonardi, S. and Ratanamahatana, C.A.: Towards parameter-free
       data mining. Proceedings of the 2004 ACM SIGKDD international conference
       on Knowledge discovery and data mining, pp. 206-215, Seattle, 2004.

     - http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/306626
       for details about the arithmetic coder.

     The additional arguments (besides the base class arguments) which have to
     be set when this field comparator is initialised are:

       compressor  The compressor to be used, currently supported are:
                   'zlib' (default) using the Python standard libray zlib.py
                                    compressor
                   'bz2' using the Python standard library bz2.py compressor
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'compressor' argument first, then call the base
       class constructor.
    """

    self.compressor = None

    # Process all keyword arguments - - - - - - - - - - - - - - - - - - - - - -
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('comp')):
        auxiliary.check_is_string('compressor', value)
        if (value not in ['zlib','bz2']):
          logging.exception('Value of argument "compressor" is not one' + \
                            'of: "zlib" or "bz2": %s' % (value))
          raise Exception
        self.compressor = value

      else:
        base_kwargs[keyword] = value

    FieldComparatorApproxString.__init__(self, base_kwargs)

    # Make sure 'compressor' attribute is set  - - - - - - - - - - - - - - - -
    #
    auxiliary.check_is_string('compressor', self.compressor)

    self.log([('Threshold', self.threshold),
              ('Compression method', self.compressor)])  # Log a message

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values using the compressor approximate string
       comparator.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    if (self.do_caching == True):  # Check if values pair is in the cache
      cache_weight = self.__get_from_cache__(val1, val2)
      if (cache_weight != None):
        return cache_weight

    if (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # Calculate the compressor similarity value - - - - - - - - - - - - - - - -
    #
    if (self.compressor == 'zlib'):
      c1 =  float(len(zlib.compress(val1)))
      c2 =  float(len(zlib.compress(val2)))
      c12 = 0.5*(len(zlib.compress(val1+val2))+len(zlib.compress(val2+val1)))

    elif (self.compressor == 'bz2'):
      c1 =  float(len(bz2.compress(val1)))
      c2 =  float(len(bz2.compress(val2)))
      c12 = 0.5*(len(bz2.compress(val1+val2)) + len(bz2.compress(val2+val1)))

    # else:  # More to be added later

    if (c12 == 0.0):
      w = self.disagree_weight  # Maximal distance

    else:  # Calculate normalised compression distance
      w = 1.0 - (c12 - min(c1,c2)) / max(c1,c2)

      if (w < 0.0):
        logging.warning('%s Compression based comparison smaller than 0.0 ' % \
                        (self.compressor)+'with strings "%s" and "%s": %.3f' \
                        % (val1, val2, w))
        w = 0.0

      assert (w >= 0.0), 'Compression: Similarity weight < 0.0'
      assert (w <= 1.0), 'Compression: Similarity weight > 1.0'

      w = self.__calc_partagree_weight__(val1, val2, w)

    if (self.do_caching == True):  # Put values pair into the cache
      self.__put_into_cache__(val1, val2, w)

    return w

# =============================================================================

class FieldComparatorTokenSet(FieldComparatorApproxString):
  """A field comparator based on token sets (i.e. words and numbers separated
     by whitespaces) and allowing for stop-words to be removed before the
     comparison.

     The additional arguments (besides the base class arguments) which have to
     be set when this field comparator is initialised are:

       stop_word_list  A list containing strings that will be removed if found
                       in the input strings before the comparison.
       common_divisor  Method of how to calculate the divisor. Can be set to
                       'average','shortest', or 'longest', and is calculated
                       according to the lengths of the two input strings.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'stop_word_list' and 'common_divisor' arguments
       first, then call the base class constructor.
    """

    self.stop_word_list = []
    self.common_divisor = None

    # Process all keyword arguments - - - - - - - - - - - - - - - - - - - - - -
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('stop_w')):
        auxiliary.check_is_list('stop_word_list', value)
        for i in range(len(value)):  # Check all elements in list are strings
          auxiliary.check_is_string('stop_word_list[%d]' % (i), value[i])
        self.stop_word_list = value

      elif (keyword.startswith('common_d')):
        auxiliary.check_is_string('common_divisor', value)
        if (value not in ['average','shortest','longest']):
          logging.exception('Value of argument "common_divisor" is not one' + \
                            'of: "average", "shortest", or "longest": %s' % \
                            (value))
          raise Exception
        self.common_divisor = value

      else:
        base_kwargs[keyword] = value

    FieldComparatorApproxString.__init__(self, base_kwargs)

    # Make sure 'common_divisor' attribute is set  - - - - - - - - - - - - - -
    #
    auxiliary.check_is_string('common_divisor', self.common_divisor)

    self.log([('Threshold', self.threshold),
              ('Stop word list', self.stop_word_list),
              ('Common divisor', self.common_divisor)])  # Log a message

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values using the token approximate string comparator.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    if (self.do_caching == True):  # Check if values pair is in the cache
      cache_weight = self.__get_from_cache__(val1, val2)
      if (cache_weight != None):
        return cache_weight

    if (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # Remove all stop words - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    clean_val1 = val1
    clean_val2 = val2

    for stop_word in self.stop_word_list:
      if stop_word in clean_val1:
        clean_val1 = clean_val1.replace(stop_word, '')
      if stop_word in clean_val2:
        clean_val2 = clean_val2.replace(stop_word, '')

    set1 = set(clean_val1.split())  # Make the cleaned values a set of tokens
    set2 = set(clean_val2.split())

    num_token1 = len(set1)
    num_token2 = len(set2)

    if ((num_token1 == 0) or (num_token2 == 0)): # No tokens in one of the sets
      w = self.disagree_weight

    else:
      num_common_token = len(set1.intersection(set2))

      if (num_common_token == 0):
        w = self.disagree_weight  # No common tokens

      else:
        if (self.common_divisor == 'average'):  # Calculate the divisor
          divisor = 0.5*(num_token1 + num_token2)
        elif (self.common_divisor == 'shortest'):
          divisor = min(num_token1, num_token2)
        else:  # Longest
          divisor = max(num_token1, num_token2)

        w = float(num_common_token) / float(divisor)

        assert (w >= 0.0), 'Token comparison: Similarity weight < 0.0'
        assert (w <= 1.0), 'Token comparison: Similarity weight > 1.0'

        w = self.__calc_partagree_weight__(val1, val2, w)

    if (self.do_caching == True):  # Put values pair into the cache
      self.__put_into_cache__(val1, val2, w)

    return w

# =============================================================================

class FieldComparatorCharHistogram(FieldComparatorApproxString):
  """A field comparator based on counting all letters, digits and whitespaces
     in two strings and building two corresponding histrograms. The similarity
     is then calculated using the cosine similarity measure between these two
     histogram vectors.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """No additional attributes needed, just call base class constructor.
    """

    FieldComparatorApproxString.__init__(self, kwargs)

    self.log([('Threshold', self.threshold)])  # Log a message

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values using the character histogram approximate
       string comparator.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    if (self.do_caching == True):  # Check if values pair is in the cache
      cache_weight = self.__get_from_cache__(val1, val2)
      if (cache_weight != None):
        return cache_weight

    if (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    histo1 = [0]*37
    histo2 = [0]*37

    workval1 = val1.lower()
    workval2 = val2.lower()

    for c in workval1:
      if (c == ' '):
        histo1[0] += 1
      elif ((c >= 'a') and (c <= 'z')):  # Count characters
        histo1[ord(c)-96] += 1
      elif ((c >= '0') and (c <= '9')):  # Count digits
        histo1[ord(c)-21] += 1

    for c in workval2:
      if (c == ' '):
        histo2[0] += 1
      elif ((c >= 'a') and (c <= 'z')):
        histo2[ord(c)-96] += 1
      elif ((c >= '0') and (c <= '9')):  # Count digits
        histo2[ord(c)-21] += 1

    vec1sum =  0.0
    vec2sum =  0.0
    vec12sum = 0.0

    for i in range(27):
      vec1sum +=  histo1[i]*histo1[i]
      vec2sum +=  histo2[i]*histo2[i]
      vec12sum += histo1[i]*histo2[i]

    if (vec1sum*vec2sum == 0.0):
      cos_sim = 0.0  # At least one vector is all zeros

    else:
      vec1sum = math.sqrt(vec1sum)
      vec2sum = math.sqrt(vec2sum)

      cos_sim = vec12sum / (vec1sum * vec2sum)

      # Due to rounding errors the similarity can be slightly larger than 1.0
      #
      cos_sim = min(cos_sim, 1.0)

    assert (cos_sim >= 0.0) and (cos_sim <= 1.0), (cos_sim, vec1, vec2)

    if (cos_sim == 0.0):
      w = self.disagree_weight
    else:
      w = self.__calc_partagree_weight__(val1, val2, cos_sim)

    if (self.do_caching == True):  # Put values pair into the cache
      self.__put_into_cache__(val1, val2, w)

    return w

# =============================================================================

class FieldComparatorTwoLevelJaro(FieldComparatorApproxString):
  """A field comparator that applies the Jaro comparator at word level, and
     additionally allows the comparison of individual words to be done using an
     approximate comparison function.

     The additional arguments (besides the base class arguments) which have to
     be set when this field comparator is initialised are:

       comp_funct     The function used to compare individual words. Either the
                      string 'equal' (default) or one of the string comparison
                      functions available in the 'stringcmp' module (i.e. a
                      function which takes two strings as input and returns a
                      similarity value between 0 and 1)
       min_threshold  Minimum threshold between 0 and 1. If an approximate
                      string comparison function is used for 'comp_funct' then
                      this 'min_threshold' needs to be set as well in order to
                      select the words that can match in the current window -
                      otherwise the 'best' match will be selected, even if it
                      has a very low similarity value.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'comp_funct' and 'min_threshold' arguments
       first, then call the base class constructor.
    """

    self.comp_funct =   'equal'
    self.min_threshold = None

    self.JARO_MARKER_CHAR = chr(1)  # Special character used to mark assigned
                                    # characters

    # Process all keyword arguments - - - - - - - - - - - - - - - - - - - - - -
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('comp_f')):
        if (isinstance(value, str)):
          if (value == 'equal'):
            self.comp_funct = 'equal'
          else:
            logging.exception('Value of "comp_funct" must be a function or' + \
                              ' the string "equal"')
            raise Exception
        else:
          auxiliary.check_is_function_or_method('comp_funct', value)
          self.comp_funct = value

      elif (keyword.startswith('min_t')):
        auxiliary.check_is_normalised('min_threshold', value)
        self.min_threshold = value

      else:
        base_kwargs[keyword] = value

    FieldComparatorApproxString.__init__(self, base_kwargs)

    # Make sure 'min_threshold' attribute is set  - - - - - - - - - - - - - -
    #
    if ((self.comp_funct != 'equal') and (self.min_threshold == None)):
      logging.exception('Comparison function is given but no minimal ' + \
                        'threshold')
      raise Exception

    self.log([('Threshold', self.threshold),
              ('Comparison function', self.comp_funct),
              ('Minimum threshold', self.min_threshold)])  # Log a message

  # ---------------------------------------------------------------------------

  def compare(self, val1, val2):
    """Compare two field values using the two-level Jaro approximate string
       comparator.
    """

    # Check if one of the values is a missing value
    #
    if (val1 in self.missing_values) or (val2 in self.missing_values):
      return self.missing_weight

    if (self.do_caching == True):  # Check if values pair is in the cache
      cache_weight = self.__get_from_cache__(val1, val2)
      if (cache_weight != None):
        return cache_weight

    if (val1 == val2):
      return self.__calc_freq_agree_weight__(val1)

    # If neither value contains a space (i.e. both are only one word) then use
    # the given word level comparison function
    #
    if (' ' not in val1) and (' ' not in val2):
      if (self.comp_funct == 'equal'):
        # Already tested if strings are the same, so here they are not
        #
        return self.disagree_weight

      # Calculate simple similarity value
      #
      w = self.comp_funct(val1, val2)

      w = self.__calc_partagree_weight__(val1, val2, w)

      if (self.do_caching == True):  # Put values pair into the cache
        self.__put_into_cache__(val1, val2, w)

      return w

    # Convert values into lists of words (whitespace separated)
    #
    list1 = val1.split()
    list2 = val2.split()

    len1 = len(list1)
    len2 = len(list2)

    halflen = max(len1,len2) / 2

    ass_list1 = []  # Words assigned in list1
    ass_list2 = []  # Words assigned in list2

    work_list1 = list1[:]  # Copy of original lists
    work_list2 = list2[:]

    common1 = 0  # Number of common characters
    common2 = 0

    # If 'equal' comparison function is given, then Jaro can be - - - - - - - -
    # directly applied at word level
    #
    if (self.comp_funct == 'equal'):

      for i in range(len1):  # Analyse the first word list
        start = max(0,i-halflen)
        end   = min(i+halflen+1,len2)
        if (list1[i] in work_list2[start:end]):  # Found common word
          ind = work_list2[start:end].index(list1[i])
          common1 += 1
          ass_list1.append(list1[i])
          work_list2[ind+start] = self.JARO_MARKER_CHAR

      for i in range(len2):  # Analyse the second string
        start = max(0,i-halflen)
        end   = min(i+halflen+1,len1)
        if (list2[i] in work_list1[start:end]):  # Found common word
          ind = work_list1[start:end].index(list2[i])
          common2 += 1
          ass_list2.append(list2[i])
          work_list1[ind+start] = self.JARO_MARKER_CHAR

      if (common1 != common2):
        logging.error('Two-level-Jaro: Wrong common values for strings ' + \
                      '"%s" and "%s"' % (val1, val2) + \
                      ', common1: %i, common2: %i' % (common1, common2) + \
                      ', common should be the same.')
        common1 = float(common1+common2) / 2.0  ##### This is just a fix #####

    # For approximate comparison function, compare all words within current
    # 'window' and keep all matches above threshold, then select the best match
    #
    else:

      for i in range(len1):  # Analyse the first word list
        start = max(0,i-halflen)
        end   = min(i+halflen+1,len2)
        search_word = list1[i]
        ind = -1  # The index of the best match found
        best_match_sim = -1
        word_ind = 0
        for word in work_list2[start:end]:
          tmp_sim = self.comp_funct(search_word, word)
          if (tmp_sim >= self.min_threshold):
            if (tmp_sim > best_match_sim):
              ind = word_ind
              best_match_sim = tmp_sim
          word_ind += 1
        if (ind >= 0):  # Found common word
          common1 += 1
          ass_list1.append(list1[i])
          work_list2[ind+start] = self.JARO_MARKER_CHAR

      for i in range(len2):  # Analyse the second string
        start = max(0,i-halflen)
        end   = min(i+halflen+1,len1)
        search_word = list2[i]
        ind = -1  # The index of the best match found
        best_match_sim = -1
        word_ind = 0
        for word in work_list1[start:end]:
          tmp_sim = self.comp_funct(search_word, word)
          if (tmp_sim >= self.min_threshold):
            if (tmp_sim > best_match_sim):
              ind = word_ind
              best_match_sim = tmp_sim
          word_ind += 1
        if (ind >= 0):  # Found common word
          common2 += 1
          ass_list2.append(list2[i])
          work_list1[ind+start] = self.JARO_MARKER_CHAR

      # For approximate comparisons, the assignment can be asymmetric, and thus
      # the values of common can differ. For example consider the following two
      # article titles:
      # - synaptic activation of transient recepter potential channels by
      #   metabotropic glutamate receptors in the lateral amygdala
      # - synaptic activation of transient receptor potential channels by
      #   metabotropic glutamate receptors in the lateral amygdala
      # In the first assignment loop, 'recepter' will match with 'receptor',
      # while in the second loop 'receptor' will match with 'receptors' (if for
      # example the q-gram comparison function is used).
      #
      if (common1 != common2):
        logging.warning('Two-level-Jaro: Different common values for ' + \
                        'strings "%s" and "%s"' % (val1, val2) + \
                        ', common1: %i, common2: %i' % (common1, common2))
        common1 = float(common1+common2) / 2.0

    if (common1 == 0):
      w = self.disagree_weight

    else:

      # Compute number of transpositions  - - - - - - - - - - - - - - - - - - -
      #
      min_num_ass_words = min(len(ass_list1), len(ass_list2))
      transposition = 0
      for i in range(min_num_ass_words):
        if (self.comp_funct == 'equal'):  # Standard way like for in basic Jaro
          if (ass_list1[i] != ass_list2[i]):
            transposition += 1
        else:  # Again use approximate string comp. to calculate similarities
          tmp_sim = self.comp_funct(ass_list1[i], ass_list2[i])
          if (tmp_sim >= self.min_threshold):
            transposition += 1

      common1 = float(common1)
      w = 1./3.*(common1 / float(len1) + common1 / float(len2) + \
               (common1-transposition) / common1)

      assert (w > 0.0), 'Two-Step Jaro: Weight is smaller than 0.0: %f' % (w)
      assert (w < 1.0), 'Two-Step Jaro: Weight is larger than 1.0: %f' % (w)

      w = self.__calc_partagree_weight__(val1, val2, w)

    if (self.do_caching == True):  # Put values pair into the cache
      self.__put_into_cache__(val1, val2, w)

    return w

# =============================================================================
