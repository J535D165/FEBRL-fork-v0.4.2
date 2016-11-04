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
# The Original Software is: "dataset.py"
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

"""Module with several classes of data set implementations.

   This module provides classes for access to different types of data sets,
   including text files (coma separated values and column wise), databases
   (using the Python database API), binary files (using Python shelves), and
   memory based data (using Python dictionaries).

   A data set is generally defined as a set of fields (or attributes), with the
   possibility that one of these fields contains unique record identifiers. If
   such a field is not available unique record identifiers will be generated on
   the fly.

   Records read from or written into a data set are basically a dictionary
   entry made of a record identifier (dictionary key) and a list of the
   attribute values of the record, like:

     {'rec-id-42':['peter','miller','42','main','st','sydney','2000','nsw']}

   See the doc strings of individual classes and methods for detailed
   documentation.

   TODO:
   - Implement SQL data set
   - Implement handling of .ZIP and .BZ2 files for CSV and COL data sets
   - Allow multiple files or tables in one data set
"""

# =============================================================================
# Import necessary modules (Python standard modules first, then Febrl modules)

import csv
import gzip
import logging
import math
import os
import random
import shelve
import string
import sys
import time

import auxiliary
import mymath

# =============================================================================
# Some constants used by the analyse() method

QUANT_LIST = [0.0,0.05,0.25,0.5,0.75,0.95,1.0]  # List of quantiles for analyse
NUM_VALUES = 15  # Maximum number of values to be shown for a column in analyse
MISS_PERC_THRES = 5.0  # Threshold in percentage above which a column will not
                       # be classified suitable for blocking in analyse

# =============================================================================

class DataSet:
  """Base class for data set access.

     This class and all its derived classes provide methods for reading from
     and writing to various data set implementations.

     Currently a data set can consists of one underlying file or table only.

     Each record in a data set is made of fields (or attributes). A list of
     tuples with the field names needs to be given when a data set is
     initialised. The first element in each tuple is a field name while the
     second element depends upon the data set type (e.g. column numbers for CSV
     files, or the names of the attributes in an SQL table).

     Records are returned as dictionary entries containing a unique record
     identifier (dictionary key) and a list of the record field (attribute)
     values from the data set.

     All data sets have the following instance variables, which can be set when
     a data set is initialised:

       description    A string describing the data set.
       access_mode    The data set's access mode, which can be either 'read',
                      'write', 'append' or 'readwrite'. Not all access modes
                      are possible for all data set implementations.
       field_list     Field names and data set specific information about the
                      fields (or attributes) in the data set.
       rec_ident      Either the name of a field in the data set containing
                      unique record identifiers, or a string (different from
                      any field name in the data set) that will be used to
                      generate record identifiers (by adding unique numbers
                      to that string) on the fly.
       strip_fields   A flag (True/False) stating if whitespaces should be
                      stripped off fields before they are returned. Default is
                      True.
       miss_val       If not None (the default), this can be either a string or
                      a list of strings. This (or these) string(s), if found in
                      an attribute, will be removed and the attribute set to an
                      empty string. Note that leading and trailing whitespaces
                      will be removed from all strings in this missing values
                      list.
"""

  # ---------------------------------------------------------------------------

  def __init__(self, base_kwargs):
    """Constructor, set general attributes.
    """

    # General attributes for all data sets
    #
    self.type =         ''
    self.description =  ''
    self.access_mode =  None
    self.field_list =   None
    self.rec_ident =    None
    self.strip_fields = True
    self.miss_val =     None

    self.num_records =  None  # To be set when a data set is initialised

    # Process base keyword arguments (all data set specific keywords were
    # processed in the derived class constructor)
    #
    for (keyword, value) in base_kwargs.items():

      if (keyword.startswith('desc')):
        auxiliary.check_is_string('description', value)
        self.description = value

      elif (keyword.startswith('access')):
        auxiliary.check_is_string('access_mode', value)
        if (value not in ['read','write','append','readwrite']):
          logging.exception('Illegal "access_mode" value: %s' % (str(value)))
          raise Exception
        self.access_mode = value

      elif (keyword.startswith('field_l')):
        auxiliary.check_is_list('field_list',value)
        self.field_list = value

      elif (keyword.startswith('rec_id')):
        auxiliary.check_is_string('rec_ident',value)
        self.rec_ident = value

      elif (keyword.startswith('strip_f')):
        auxiliary.check_is_flag('strip_fields', value)
        self.strip_fields = value

      elif (keyword.startswith('miss_v')):
        if (value != None):
          if (isinstance(value, str)):
            value = [value]  # Make it a list
          auxiliary.check_is_list('miss_val',value)
        self.miss_val = value

      else:
        logging.exception('Illegal constructor argument keyword: "%s"' % \
              (str(keyword)))
        raise Exception

    # Check if access mode and record identifier are defined - - - - - - - - -
    #
    auxiliary.check_is_string('access_mode', self.access_mode)
    auxiliary.check_is_string('rec_ident', self.rec_ident)

    # Remove all leading and trailing whitespaces from missing values - - - - -
    #
    if (self.miss_val != None):

      clean_miss_val = []

      for miss_val in self.miss_val:

        auxiliary.check_is_string('miss_val entry', miss_val)

        stripped_miss_val = miss_val.strip()
        if (stripped_miss_val != ''):  # Not empty string
          clean_miss_val.append(stripped_miss_val)

      if (clean_miss_val != []):  # Check if there are non-empty missing values
        self.miss_val = clean_miss_val
      else:
        self.miss_val = None

  # ---------------------------------------------------------------------------

  def finalise(self):
    """Finalise a data set.
       See implementations in derived classes for details.
    """

    logging.exception('Override abstract method in derived class')
    raise Exception

  # ---------------------------------------------------------------------------

  def read(self, *recs):
    """Read and return one or more records.

       Calling arguments are not the same for all derived data set classes, see
       implementations in derived classes for details.
    """

    logging.exception('Override abstract method in derived class')
    raise Exception

  # ---------------------------------------------------------------------------

  def readall(self):
    """An iterator which will return one record per call as a tuple (record
       identifier, record field list).

       Use like:  for rec in dataset.readall():

       See implementations in derived classes for details.
    """

    logging.exception('Override abstract method in derived class')
    raise Exception

  # ---------------------------------------------------------------------------

  def write(self, rec_dict):
    """Write one or more records into the data set.
       See implementations in derived classes for details.
    """

    logging.exception('Override abstract method in derived class')
    raise Exception

  # ---------------------------------------------------------------------------

  def analyse(self, sample, word_analysis, log_funct=None, log_num_recs=None):
    """Read the data and analyse a sample (or all) of the records in it.

       The following arguments have to be set:

         sample         A number between 0 and 100 that gives the percentage of
                        records that will be randomly selected and analysed. If
                        set to 100 then all records in the data set will be
                        analysed.
         word_analysis  This flag determines if the values in fields should
                        be analysed as separate words or as complete values.
                        If set to True, then a value that contains more than
                        one word (i.e. contains at least one whitespace) will
                        be split into its words and they will be used in the
                        analysis for frequencies, rather than the original
                        value.
         log_funct      This can be a Python function or method which will log
                        (print or save to a file) a progress report message. It
                        is assumed that this function or method has one input
                        argument of type string (the message to be printed).
         log_num_recs   If 'log_funct' is defined, this must be a positive
                        integer number which interval in number of records read
                        from the data set between calls to the 'log_funct'
                        function.

       For each field (column / attribute) in the data set this method collects
       and generates the following basic statistics for the sampled records:
       - The number of unique values
       - The smallest and largest values (as strings)
       - The average frequency and standard deviation of the values
       - A list of quantiles of the values
       - The most and least frequent values and their frequencies
       - The minimum and maximum value length
       - If the field only contains digits, only contains letters, or only
         contains digits and letters
       - The maximum number of spaces in values
       - The number of records with missing value

       It then generates a table summarising the field statistics and a table
       summarising the quantiles statistics.

       Finally, it generates a table containing details on the suitability of
       fields for blocking (according to their number of values and proportion
       of missing values).

       It returns a list containing strings (each assumed to be a line of text)
       that can be printed or save into a file.
    """

    # Check input parameters
    #
    auxiliary.check_is_percentage('Sample percentage', sample)
    auxiliary.check_is_flag('Word analysis', word_analysis)

    if (log_funct != None):
      auxiliary.check_is_function_or_method('Log function', log_funct)
      auxiliary.check_is_integer('Log number of records', log_num_recs)
      auxiliary.check_is_positive('Log number of records', log_num_recs)

    test_rec_dict = self.read()  # Read one record to get the number of fields
    num_fields = len(test_rec_dict.values()[0])

    # Define list of data structures for information to be collected
    #
    num_missing_list =    []
    values_dict_list =    []
    min_length_list =     []
    max_length_list =     []
    isdigit_list =        []
    isalpha_list =        []
    isalnum_list =        []
    max_num_spaces_list = []
    warn_message_dict = {} # Keys will be warning messages, values their counts

    for i in range(num_fields):
      num_missing_list.append(0)
      values_dict_list.append({})
      min_length_list.append(999)
      max_length_list.append(0)
      isdigit_list.append(True)
      isalpha_list.append(True)
      isalnum_list.append(True)
      max_num_spaces_list.append(0)

    num_warnings = 0  # Count all warnings

    num_records = 0        # Nunmber of records in data set
    num_recs_analysed = 0  # Number of records analysed (sampled)

    # Read and process data lines - - - - - - - - - - - - - - - - - - - - - - -
    #
    start_time = time.time()

    for (rec_id, rec_list) in self.readall():

      if (len(rec_list) != num_fields):
        warn_msg = 'Line does have %d fields (not %d as expected)' % \
                   (len(rec_list), num_fields)
        warn_msg_count = warn_message_dict.get(warn_msg, 0) + 1
        warn_message_dict[warn_msg] = warn_msg_count

        num_warnings += 1

        if (len(rec_list) < num_fields):  # Correct by adding empty fields
          rec_list += ['']*(num_fields-len(rec_list))

      if (random.random()*100 < sample):  # Randomly select a record

        c = 0
        # Loop over the field values in this record
        #
        for col_val in rec_list[:num_fields]:

          # Value has been stripped and/or missing values removed in .readall()

          if (col_val == ''):
            num_missing_list[c] += 1

          else:  # Value is not empty

            if (word_analysis == True): # Split into words and analyse separate
              col_val_list = col_val.split()  # Split at whitespaces
            else:
              col_val_list = [col_val]

            for col_val in col_val_list:
              col_dict = values_dict_list[c]
              col_val_freq = col_dict.get(col_val, 0) + 1  # Increase count
              col_dict[col_val] = col_val_freq
              values_dict_list[c] =   col_dict

            col_val_len = len(col_val)
            min_length_list[c] = min(min_length_list[c], col_val_len)
            max_length_list[c] = max(max_length_list[c], col_val_len)

            isdigit_list[c] = isdigit_list[c] and col_val.isdigit()
            isalpha_list[c] = isalpha_list[c] and col_val.isalpha()
            isalnum_list[c] = isalnum_list[c] and col_val.isalnum()

            if (' ' in col_val):  # Count number of whitespaces in value
              max_num_spaces_list[c] = max(max_num_spaces_list[c],
                                           (len(col_val.split(' '))-1))
          c += 1

        num_recs_analysed += 1

      num_records += 1

      if ((log_funct != None) and (log_num_recs != None) and \
          ((num_records % log_num_recs) == 0)):
        used_time = (time.time() - start_time)
        processed_perc = 100.0 * float(num_records) / self.num_records
        log_funct('Read %.2f%% of %d records in %.2f sec ' % \
                  (processed_perc, self.num_records, used_time))
        #log_funct('Processed a %.1f%% sample of %d records in %.2f sec' % \
        #          (sample, num_records, used_time))

    if ((log_funct != None) and (log_num_recs != None)):
      used_time = (time.time() - start_time)
      log_funct('Finished reading a %.1f%% sample of %d records in %.2f' % \
                (sample, num_records, used_time)+' sec')

    # Log warnings if there were some - - - - - - - - - - - - - - - - - - - - -
    #
    if (warn_message_dict != {}):
      warn_msg_tuples = warn_message_dict.items()
      warn_msg_tuples.sort()

      logging.warn('Warnings from "analyse" read records:')
      for (warn_msg, warn_count) in warn_msg_tuples:
        logging.warn('  %s occured %d times' % (warn_msg, warn_count))

    # Calculate and report final statistics - - - - - - - - - - - - - - - - - -
    #
    final_stats =   []  # Build a list of lines to be returned
    stddev_list =   []
    avrg_list =     []
    num_val_list =  []
    num_miss_list = []
    type_list =     []  # If field is digits only, letters only, etc.

    final_stats.append('')
    final_stats.append('Detailed statistics for data set: %s' % \
                       (self.file_name))
    final_stats.append('='*(34+len(self.file_name)))
    final_stats.append('')
    final_stats.append('Analysis conducted on '+time.asctime())
    final_stats.append('')
    final_stats.append('Analysed %d of %d records (%.2f%% sample)' % \
                       (num_recs_analysed, num_records, sample))
    final_stats.append('')

    if (word_analysis == True):
      final_stats.append('Frequency analysis based on separate words ' + \
                         'in fields')
    else:
      final_stats.append('Frequency analysis based on field values')

    final_stats.append('')

    for c in range(num_fields):

      header_line_str = 'Field %d: "%s"' % (c, self.field_list[c][0])
      final_stats.append(header_line_str)

      value_list = values_dict_list[c].keys()
      freq_list = values_dict_list[c].values()
      freq_list_len = len(freq_list)
      final_stats.append('  Number of unique values: %d' % (freq_list_len))
      num_val_list.append(freq_list_len)

      if (freq_list_len > 0):
        list_tuples = map(None, freq_list, value_list)
        list_tuples.sort()
        value_list.sort()
        final_stats.append('    Smallest and largest values (as strings): ' + \
                       '"%s" / "%s"' % (str(value_list[0]),
                       str(value_list[-1])))

        avrg = float(sum(freq_list)) / float(freq_list_len)
        final_stats.append('  Average frequency: %.2f' % (avrg))
        avrg_list.append(avrg)

        stddev = 0.0
        for v in freq_list:
          stddev += (v - avrg)*(v - avrg)
        stddev = math.sqrt(stddev / float(freq_list_len))
        final_stats.append('  Frequency stddev:  %.2f' % (stddev))
        stddev_list.append(stddev)

        quantiles_str = auxiliary.str_vector(mymath.quantiles(freq_list,
                                              QUANT_LIST), num_digits=2)
        final_stats.append('  Quantiles: %s' % (quantiles_str))

        if (freq_list_len < NUM_VALUES):
          final_stats.append('  All field values:')
          for lt in list_tuples:
            final_stats.append('    '+str(lt))
        else:
          final_stats.append('  Most frequent field values:')
          for j in [-1, -2, -3, -4, -5, -6]:
            final_stats.append('    '+str(list_tuples[j]))
          final_stats.append('  Least frequent field values:')
          for j in [5, 4, 3, 2, 1, 0]:
            final_stats.append('    '+str(list_tuples[j]))

        final_stats.append('  Minimum and maximum value lengths: %d / %d' % \
                           (min_length_list[c], max_length_list[c]))
        final_stats.append('  Is-digit: %s, is-alpha: %s, is-alnum: %s' % \
                           (str(isdigit_list[c]),str(isalpha_list[c]),
                           str(isalnum_list[c])))
        final_stats.append('  Maximum number of spaces in values: %d' % \
                           (max_num_spaces_list[c]))
        final_stats.append('  Number of records with missing value: %d' % \
                           (num_missing_list[c]))
        num_miss_list.append(num_missing_list[c])

        # Determine if field contents are digits only, letters only, etc.
        #
        if (isdigit_list[c] == True):
          type_list.append('Only digits')
        elif (isalpha_list[c] == True):
          type_list.append('Only letters')
        elif (isalnum_list[c] == True):
          type_list.append('Digits and letters')
        else:
          type_list.append('Various')

      else:  # freq_list_len = 0 (no values in a field at all)
        avrg_list.append(-1)
        stddev_list.append(-1)
        type_list.append('Empty')
        num_miss_list.append(num_records)
        num_val_list.append(0)
        final_stats.append('    Empty field!')

      final_stats.append('')

    # Give a summary of fields - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    final_stats.append('Summary of field statistics')
    final_stats.append('===========================')
    final_stats.append('')

    final_stats.append('                            Unique    Missing     ' + \
                       '   Frequencies')

    final_stats.append('Field names                 values     values     ' + \
                       '  Avrg    StdDev  Field type')
    final_stats.append('-'*86)

    for c in range(num_fields):
      hv = self.field_list[c][0]

      av = '%.2f' % (avrg_list[c])
      if (av != '-1.00'):
        av = av.rjust(10)
      else:
        av = '        --'

      sv = '%.2f' % (stddev_list[c])
      if (sv != '-1.00'):
        sv = sv.rjust(9)
      else:
        sv = '       --'

      tv = type_list[c]

      final_stats.append('%22s  %10d %10d %s %s  %s' % \
           (hv.ljust(22), num_val_list[c], num_miss_list[c], av, sv, tv))
    final_stats.append('')

    # Give a summary of quantiles - - - - - - - - - - - - - - - - - - - - - - -
    #
    final_stats.append('Field quantiles')
    final_stats.append('===============')
    final_stats.append('')

    quant_list_str = '['
    for v in QUANT_LIST:
      quant_list_str += '%d%%, ' % (int(v*100))
    quant_list_str = quant_list_str[:-2]+']'

    final_stats.append('Field names             Quantiles  %s' % \
                       (quant_list_str))
    final_stats.append('-'*86)

    for c in range(num_fields):
      hv = self.field_list[c][0]

      freq_list = values_dict_list[c].values()

      if (len(freq_list) > 0):

        qv = auxiliary.str_vector(mymath.quantiles(freq_list, QUANT_LIST),
                                  num_digits=2)
        final_stats.append('%22s  %s' % (hv.ljust(22), qv))
      else:  # An empty field
        final_stats.append('%22s  [Empty]' % (hv.ljust(22)))

    final_stats.append('')

    # Calculate number of blocks as well as percentage of missing values - - -
    #
    final_stats.append('Suitability of fields for blocking')
    final_stats.append('==================================')
    final_stats.append('(number of record pair comparisons is calculated ' + \
                       'for a deduplication)')
    final_stats.append('')

    final_stats.append('Field names             Suitability')
    final_stats.append('-'*86)

    for c in range(num_fields):
      hv = self.field_list[c][0]

      num_miss_rec = num_missing_list[c]

      if (num_miss_rec > 0):
        miss_perc = 100.0 * float(num_miss_rec) / float(num_records)

        if (miss_perc > MISS_PERC_THRES):  # Field not suitable for blocking
          final_stats.append('%22s  %.2f%% (%d) records with missing ' % \
                             (hv.ljust(22), miss_perc, num_miss_rec) + \
                             'values (not suitable)')
      else:
        miss_perc = 0.0

      if (miss_perc <= MISS_PERC_THRES):  # Field suitable for blocking

        freq_list = values_dict_list[c].values()
        num_val =   len(freq_list)

        # Start with missing value records
        #
        num_comp = num_miss_rec*(num_miss_rec-1)

        for j in freq_list:
          num_comp += j*(j-1)  # Add number of record pair comparisons

        final_stats.append('%22s  %d unique values, resulting in %d ' % \
                           (hv.ljust(22), num_val, num_comp) + \
                           'record pair comparisons')

        if (miss_perc > 0.0):
          final_stats.append(' '*24+'Note field contains %.2f%% (%d) ' % \
                         (miss_perc, num_miss_rec) + 'records with missing' + \
                         ' values')
      final_stats.append('')

    final_stats.append('')

    return final_stats

  # ---------------------------------------------------------------------------

  def log(self, instance_var_list = None):
    """Write a log message with the basic data set instance variables plus the
       instance variable provided in the given input list (assumed to contain
       pairs of names (strings) and values).
    """

    logging.info('')
    logging.info('Data set: "%s"' % (self.description))
    logging.info('  Data set type:     %s' % (self.dataset_type))
    logging.info('  Access mode:       %s' % (self.access_mode))
    logging.info('  Strip fields:      %s' % (str(self.strip_fields)))
    logging.info('  Number of records: %d' % (self.num_records))
    logging.info('  Record identifier: %s' % (self.rec_ident))
    if (self.miss_val != None):
      logging.info('  Missing values:    %s' % (str(self.miss_val)))
    else:
      logging.info('  Missing values:    None')
    logging.info('  Fields list:')
    for (field_name,val) in self.field_list:
      if (field_name == self.rec_ident):
        logging.info('    %s: %s (record identifier)' % \
                     (str(field_name), str(val)))
      else:
        logging.info('    %s: %s' % (str(field_name), str(val)))

    if (instance_var_list != None):
      logging.info('  Data set specific variables:')

      max_name_len = 0
      for (name, value) in instance_var_list:
        max_name_len = max(max_name_len, len(name))

      for (name, value) in instance_var_list:
        pad_spaces = (max_name_len-len(name))*' '
        logging.info('    %s %s' % (name+':'+pad_spaces, str(value)))

# =============================================================================

class DataSetCSV(DataSet):
  """Implementation of a CSV (comma separated values) data set class.

     The 'field_list' attribute can either be given when a CSV data set is
     initialised (in which case it must contain pairs of field names and column
     numbers, starting from zero), or it can be generated from the header line
     in the CSV data file (if the 'header_line' attribute is set to True and
     the 'access_mode' is 'read').

     Possible values for the 'access_mode' argument are: 'read', 'write', or
     'append' (but not 'readwrite').

     If the file name ends with '.gz' it is assumed the file is GZIP compressed
     and it will be opened using the gzip library. This will only be checked
     when opening a CSV data set for reading, writing and appending will always
     be into uncompressed files.

     The additional arguments (besides the base class arguments) which have to
     be set when this data set is initialised are:

       file_name         A string containing the name of the underlying CSV
                         file.
       delimiter         A one character string which designates the delimter
                         used to split a line into columns/attributes. Default
                         value is a comma (','), aother commen alternative is
                         a tabulator ('\t').
       header_line       A flag, if set to True (and access mode is "read") the
                         first line in the CSV data file is assumed to contain
                         the field (attribute) names and these will be used to
                         generate the 'field_list' attribute. If set to False
                         (default) then the 'field_list' attribute needs to be
                         provided.
       write_header      A flag, if set to "True" a header line with the field
                         names is written into the file when it is opened in
                         'write' or 'append' (only if file empty) mode.
       write_quote_char  A quote character, used when writing to file. Default
                         is no quote character (empty string '').

     Note that all values returned from a CSV data set (from it's read methods)
     are strings, while non-string values written to the data set will be
     stored as strings in the CSV data file.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the derived attributes first, then call the base
       class constructor.
    """

    self.dataset_type = 'CSV'

    self.file_name =        None   # The name of the CSV file
    self.header_line =      False  # Flag, default set to no header line
    self.write_header =     False  # Flag, set to not write header line
    self.file =             None   # File pointer to current file
    self.write_quote_char = ''     # The quote character for writing fields
    self.delimiter  =       ','    # The delimiter character

    self.next_rec_num =  None
    self.rec_ident_col = -1    # Column of the record identifier field

    # Process all keyword arguments
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('file')):
        auxiliary.check_is_string('file_name', value)
        self.file_name = value

      elif (keyword.startswith('header')):
        auxiliary.check_is_flag('header_line', value)
        self.header_line = value

      elif (keyword.startswith('write_he')):
        auxiliary.check_is_flag('write_header', value)
        self.write_header = value

      elif (keyword.startswith('write_qu')):
        auxiliary.check_is_string('write_quote_char', value)
        self.write_quote_char = value

      elif (keyword.startswith('delimi')):
        auxiliary.check_is_string('delimiter', value)
        if (len(value) != 1):
          logging.exception('Value of "delimiter" argument must be a one-' + \
                            'character string, but it is: "%s"' % (delimiter))
          raise Exception
        self.delimiter = value

      else:
        base_kwargs[keyword] = value

    DataSet.__init__(self, base_kwargs)  # Process base arguments

    if ((self.header_line == False) and (self.field_list == None)):
      logging.exception('Argument "field_list" must be given if field ' + \
                        'names are not taken from header line')
      raise Exception

    if ((self.header_line == True) and (self.access_mode == 'read')):
      self.field_list = []  # Will be generated from file header line


    # Make sure the 'file_name' attribute is set - - - - - - - - - - - - - - -
    #
    auxiliary.check_is_string('file_name', self.file_name)

    # Now perform various checks for each access mode and open file - - - - - -
    #
    if (self.access_mode == 'read'):

      if (self.file_name.endswith('.gz')) or (self.file_name.endswith('.GZ')):
        file_gzipped = True

        try:
          self.file = gzip.open(self.file_name) # Open gzipped file
        except:
          logging.exception('Cannot open gzipped CSV file "%s" for reading' % \
                            (self.file_name))
          raise IOError

      else:  # Open normal file for reading
        file_gzipped = False

        try:  # Try to open the file in read mode
          self.file = open(self.file_name,'r')
        except:
          logging.exception('Cannot open CSV file "%s" for reading' % \
                            (self.file_name))
          raise IOError

      # Initialise the CSV parser - - - - - - - - - - - - - - - - - - - - - - -
      #
      self.csv_parser = csv.reader(self.file, delimiter = self.delimiter)

      # If header line is set to True get field names
      #
      if (self.header_line == True):
        header_line = self.csv_parser.next()
        self.field_list = []  # Make sure field list is empty

        col_num = 0  # Create a list with field names and column numbers
        for field_name in header_line:
          if (self.strip_fields == True):
            field_name = field_name.strip()
          self.field_list.append((field_name,col_num))
          col_num += 1

      # Count number of records in the file
      #
      if ((sys.platform[0:5] in ['linux','sunos']) and \
          (file_gzipped == False)):  # Fast line counting
        if (' ' not in self.file_name):
          wc = os.popen('wc -l ' + self.file_name)
        else:
          wc = os.popen('wc -l "%s"' % (self.file_name))

        num_rows = int(string.split(wc.readline())[0])
        wc.close()
      else:  # Slow line counting method

        num_rows = 0

        if (file_gzipped == True):
          fp = gzip.open(self.file_name)
        else:
          fp = open(self.file_name,'r')
        for l in fp:
          num_rows += 1
        fp.close()

      self.num_records = num_rows
      if (self.header_line == True):
        self.num_records -= 1

      # Check that there are records in the data set
      #
      if (self.num_records == 0):
        logging.exception('No records in CSV data set opened for reading')
        raise Exception

      self.next_rec_num = 0

    elif (self.access_mode == 'write'):   # - - - - - - - - - - - - - - - - - -

      # Try to open the file in write mode
      #
      try:
        self.file = open(self.file_name,'w')
      except:
        logging.exception('Cannot open CSV file "%s" for writing' % \
                          (self.file_name))
        raise IOError

      # Initialise the CSV parser as writer - - - - - - - - - - - - - - - - - -
      #
      self.csv_parser = csv.writer(self.file, delimiter = self.delimiter)

      # Write the header line with field names if desired
      #
      if (self.write_header == True):
        header_list = []
        for (field_name,col) in self.field_list:
          if (self.strip_fields == True):
            field_name = field_name.strip()
          if (self.write_quote_char != ''):
            field_name = self.write_quote_char + field_name + \
                         self.write_quote_char
          header_list.append(field_name)

        self.csv_parser.writerow(header_list)
        self.file.flush()

      self.num_records =  0
      self.next_rec_num = 0

    elif (self.access_mode == 'append'):  # - - - - - - - - - - - - - - - - - -

      # Try to open the file in append mode
      #
      try:
        self.file = open(self.file_name,'a')
      except:
        logging.exception('Cannot open CSV file "%s" for appending' % \
                          (self.file_name))
        raise IOError

      # Count number of records in the file
      #
      if (sys.platform[0:5] in ['linux','sunos']):  # Fast line counting
        wc = os.popen('wc -l ' + self.file_name)
        num_rows = int(string.split(wc.readline())[0])
        wc.close()
      else:  # Slow line counting method
        num_rows = 0
        fp = open(self.file_name,'r')
        for l in fp:
          num_rows += 1
        fp.close()

      # Initialise the CSV parser as writer - - - - - - - - - - - - - - - - - -
      #
      self.csv_parser = csv.writer(self.file, delimiter = self.delimiter)

      # If no records are stored write header if desired
      #
      if (num_rows == 0) and (self.write_header == True):
        header_list = []
        for (field_name,col) in self.field_list:
          if (self.strip_fields == True):
            field_name = field_name.strip()
          if (self.write_quote_char != ''):
            field_name = self.write_quote_char + field_name + \
                         self.write_quote_char
          header_list.append(field_name)

        self.csv_parser.writerow(header_list)
        self.file.flush()

        self.num_records =  0
        self.next_rec_num = 0

      else:  # There are records in the data set already

        self.num_records = num_rows
        if (self.header_line == True):
          self.num_records -= 1

        self.next_rec_num = self.num_records

    else:  # Illegal data set access mode - - - - - - - - - - - - - - - - - - -

      logging.exception('Illegal data set access mode: "%s" (not allowed ' % \
            (str(self.access_mode)) + 'with CSV data set implementation).')
      raise Exception

    # Check if the values in the 'field_list' argument are all positive integer
    # numbers (column numbers), and also check if one is the record identifier
    #
    this_col_num = 0  # Check if column numbers are consecutive

    for (field_name, field_col) in self.field_list:

      auxiliary.check_is_string('Field name for column %d is not a string' % \
                                 (field_col), field_name)
      auxiliary.check_is_integer('Value of field column "%s" is not an ' % \
                                  (field_name) + ' integer number', field_col)
      auxiliary.check_is_not_negative('Value of field column "%s" is a ' % \
                                  (field_name) + 'negative number', field_col)
      if (this_col_num != field_col):
        logging.exception('Column numbers are not consecutive: %s' % \
                          (str(self.field_list)))
        raise Exception
      else:
        this_col_num += 1

      # Check if this is the record identifier field
      #
      if (self.rec_ident == field_name):
        self.rec_ident_col = field_col

    self.log([('CSV file name', self.file_name),
              ('Header line', self.header_line),
              ('Write header', self.write_header),
              ('Quote character', self.write_quote_char),
              ('Record identifier column', self.rec_ident_col),
              ('Delimiter', self.delimiter)])

  # ---------------------------------------------------------------------------

  def finalise(self):
    """Finalise a data set, i.e. close the file and set various attributes to
       None.
    """

    # Close file if it is open
    #
    if (self.file != None):

      self.file.close()
      self.file = None

    self.access_mode =  None
    self.file_name =    None
    self.num_records =  None
    self.next_rec_num = None

    # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.info('Finalised CSV data set "%s"' % (self.description))

  # ---------------------------------------------------------------------------

  def __read_one_record__(self):
    """Read and return the next record. Should not be used from outside the
       module.

       Updates 'next_rec_num' attribute.

       Returns a dictionary {record_identifier : record_field_list}, or {} once
       the end of the file is reached.
    """

    try:
      rec = self.csv_parser.next()
    except:
      return {}  # Reached end of file

    if (self.strip_fields == True):  # Strip leading and trailing whitespace -
      rec = map(string.strip, rec)

    if (self.miss_val != None):  # Check for missing values in record - - - - -
      clean_rec = []
      miss_val_list = self.miss_val  # Faster reference access

      for val in rec:
        if (val in miss_val_list):  # Found a missing value
          clean_rec.append('')  # Replace with empty string
        else:
          clean_rec.append(val)
      rec = clean_rec

    if (self.rec_ident_col == -1):  # Generate or get record identifier - - - -
      rec_ident = self.rec_ident+'-%d' % (self.next_rec_num)

    else:  # Get record identifier from the record itself
      rec_ident = rec[self.rec_ident_col]

    self.next_rec_num += 1

    return {rec_ident:rec}

  # ---------------------------------------------------------------------------

  def read(self, *recs):
    """Read and return one or more records.
       - If no argument is given return the next record in the data set.
       - If one argument (n) is given (must be a positive number) return the
         next n records (or less if there are not enough records in the data
         set).
       - If two arguments (s,n) are given (both must be numbers) return the
         next n records starting from record s (or less if there are not enough
         records in the data set).

       Returns an empty dictionary if no more records are available.
    """

    if (self.file == None):
      logging.exception('Data set not initialised')
      raise Exception

    if (self.access_mode != 'read'):
      logging.exception('Data set not initialised for "read" access')
      raise Exception

    len_recs_arg = len(recs)

    # Process the three different forms of calling
    #
    if (len_recs_arg == 0):  # No arguments, just read one record - - - - - - -

      return self.__read_one_record__()

    elif (len_recs_arg == 1):  # One argument, read n records - - - - - - - - -

      num_recs = recs[0]
      if (not isinstance(num_recs, int) or (num_recs < 1)):
        logging.exception('Number of records given is not a positive integer' \
                          + ' number: %s' % (str(num_recs)))
        raise Exception

      rec_dict = {}

      for i in range(num_recs):

        this_rec_dict = self.__read_one_record__()

        if (this_rec_dict == {}):  # No more records in data set
          return rec_dict

        if this_rec_dict.keys()[0] in rec_dict:  # Check for unique rec ident.
          logging.exception('Already existing record identifier returned: %s' \
                            % (str(this_rec_dict)))
          raise Exception

        rec_dict.update(this_rec_dict)

      return rec_dict

    elif (len_recs_arg == 2):  # Two arguments, read n records - - - - - - - -

      start_num = recs[0]
      num_recs =  recs[1]
      if (not isinstance(start_num, int) or (start_num < 0)):
        logging.exception('Start record number given is not a positive ' + \
                          'integer number: %s' % (str(start_num)))
        raise Exception
      if (not isinstance(num_recs, int) or (num_recs < 1)):
        logging.exception('Number of records given is not a positive integer' \
                          + ' number: %s' % (str(num_recs)))
        raise Exception
      if (start_num >= self.num_records):
        logging.exception('Start record number is larger than the number ' + \
                          'of records in the data set: %d' % (start_num))
        raise Exception

      # Check if the start record number is at the current position or not
      #
      if (start_num != self.next_rec_num):

        self.file.close()  # Close currently open file

        if (self.file_name.endswith('.gz')) or \
           (self.file_name.endswith('.GZ')):
          self.file = gzip.open(self.file_name) # Open gzipped file
        else:
          self.file = open(self.file_name,'r')  # Re-open file

        # Initialise the CSV parser as reader
        #
        self.csv_parser = csv.reader(self.file, delimiter = self.delimiter)

        # Skip over header (if there is one) and skip to start record
        #
        if (self.header_line == True):
          self.csv_parser.next()

        for i in range(start_num):
          self.csv_parser.next()

        self.next_rec_num = start_num  # Update record counter

      # Now read the desired number of records - - - - - - - - - - - - - - - -
      #
      rec_dict = {}

      for i in range(num_recs):

        this_rec_dict = self.__read_one_record__()

        if (this_rec_dict == {}):  # No more records in data set
          return rec_dict

        if this_rec_dict.keys()[0] in rec_dict:  # Check for unique rec ident.
          logging.exception('Already existing record identifier returned: %s' \
                            % (str(this_rec_dict)))
          raise Exception

        rec_dict.update(this_rec_dict)

      return rec_dict

    else:  # Illegal calling (3 or more arguments)
      logging.exception('Illegal call to read(): 3 or more arguments: %s' % \
                        (str(recs)))
      raise Exception

  # ---------------------------------------------------------------------------

  def readall(self):
    """An iterator which will return one record per call as a tuple (record
       identifier, record field list).

       The file is first closed and then re-opened returning the first record.
    """

    if (self.file == None):
      logging.exception('Data set not initialised')
      raise Exception

    if (self.access_mode != 'read'):
      logging.exception('Data set not initialised for "read" access')
      raise Exception

    self.file.close()  # Close currently open file

    if (self.file_name.endswith('.gz')) or (self.file_name.endswith('.GZ')):
      self.file = gzip.open(self.file_name) # Open gzipped file
    else:
      self.file = open(self.file_name,'r')  # Re-open file

    self.next_rec_num = 0   # Initialise next record counter

    # Initialise the CSV parser as reader
    #
    self.csv_parser = csv.reader(self.file, delimiter = self.delimiter)

    # Skip over header (if there is one) and skip to start record
    #
    if (self.header_line == True):
      self.csv_parser.next()

    for rec in self.csv_parser:

      if (self.strip_fields == True):  # Strip leading and trailing whitespace
        rec = map(string.strip, rec)

      if (self.miss_val != None):  # Check for missing values in record
        clean_rec = []
        miss_val_list = self.miss_val  # Faster reference access

        for val in rec:
          if (val in miss_val_list):  # Found a missing value
            clean_rec.append('')  # Replace with empty string
          else:
            clean_rec.append(val)
        rec = clean_rec

      if (self.rec_ident_col == -1):  # Generate record identifier
        rec_ident = self.rec_ident+'-%d' % (self.next_rec_num)

      else:  # Get record identifier from the record itself
        rec_ident = rec[self.rec_ident_col]

      self.next_rec_num += 1

      yield (rec_ident,rec)

  # ---------------------------------------------------------------------------

  def write(self, rec_dict):
    """Write one or more records into the data set.
       The input dictionary with records is first sorted (according to the
       record identifiers) and then written into the CSV file.
    """

    if (self.file == None):
      logging.exception('Data set not initialised')
      raise Exception

    if (self.access_mode not in ['write','append']):
      logging.exception('Data set not initialised for "write" or "append" ' + \
                        'access')
      raise Exception

    # Sort record identifiers
    #
    rec_ident_keys = rec_dict.keys()
    rec_ident_keys.sort()

    for rec_ident in rec_ident_keys:
      rec = rec_dict[rec_ident]

      if (self.strip_fields == True):  # Strip leading and trailing whitespace
        rec = map(string.strip,rec)

      if (self.miss_val != None):  # Check for missing values in record
        clean_rec = []
        miss_val_list = self.miss_val  # Faster reference access

        for val in rec:
          if (val in miss_val_list):  # Found a missing value
            clean_rec.append('')  # Replace with empty string
          else:
            clean_rec.append(val)
        rec = clean_rec

      if (self.write_quote_char != ''):
        rec = map(lambda s:self.write_quote_char+s+self.write_quote_char, rec)

      self.csv_parser.writerow(rec)  # Write record to file

      self.num_records +=  1
      self.next_rec_num += 1

    self.file.flush()

# =============================================================================

class DataSetCOL(DataSet):
  """Implementation of a column wise fields (with fixed width) data set class.

     The 'field_list' attribute can either be a list of positive integer
     numbers only (assumed to be the column width) if the 'header_line'
     attribute is set to True (so that the field names are taken from the first
     line in the data set) and the 'access_mode' is 'read', or (if
     'header_line' is set to False or access is not 'read') a list with tuples
     made of a field name (strings) and column width (positive integers).

     Possible values for the 'access_mode' argument are: 'read', 'write', or
     'append' (but not 'readwrite').

     If the file name ends with '.gz' it is assumed the file is GZIP compressed
     and it will be opened using the gzip library. This will only be checked
     when opening a COL data set for reading, writing and appending will always
     be into a uncompressed file.

     The additional arguments (besides the base class arguments) which have to
     be set when this data set is initialised are:

       file_name         A string containing the name of the underlying COL
                         file.
       header_line       A flag, if set to True (and access mode is "read") the
                         first line in the COL data file is assumed to contain
                         the field (attribute) names and these will be used to
                         generate the 'field_list' attribute. If set to False
                         (default) then the 'field_list' attribute needs to be
                         provided.
       write_header      A flag, if set to "True" a header line with the field
                         names is written into the file when it is opened in
                         'write' or 'append' (only if file empty) mode.

     Note that all values returned from a COL data set (from it's read methods)
     are strings, while non-string values written to the data set will be
     stored as strings in the COL data file.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the derived attributes first, then call the base
       class constructor.
    """

    self.dataset_type = 'COL'

    self.file_name =        None   # The name of the COL file
    self.header_line =      False  # Flag, default set to no header line
    self.write_header =     False  # Flag, set to not write header line
    self.file =             None   # File pointer to current file

    self.col_start_end = None  # List of tuples with column starts and ends
    self.next_rec_num =  None

    self.rec_ident_field = -1    # Number of the record identifier field

    # Process all keyword arguments
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('file')):
        auxiliary.check_is_string('file_name', value)
        self.file_name = value

      elif (keyword.startswith('header')):
        auxiliary.check_is_flag('header_line', value)
        self.header_line = value

      elif (keyword.startswith('write_he')):
        auxiliary.check_is_flag('write_header', value)
        self.write_header = value

      else:
        base_kwargs[keyword] = value

    DataSet.__init__(self, base_kwargs)  # Process base arguments

    # Check if field list argument is given and contains correct details
    #
    auxiliary.check_is_not_none('field_list', self.field_list)
    auxiliary.check_is_list('Field list', self.field_list)

    field_cnt = 0
    for field_data in self.field_list:  # Check all entries in field list

      if ((self.header_line == True) and (self.access_mode == 'read')):
        auxiliary.check_is_integer('field_list[%d]' % (field_cnt), field_data)
        auxiliary.check_is_positive('field_list[%d]' % (field_cnt), field_data)

      else:  # Check if field list details are given (names, column widths)

        auxiliary.check_is_tuple('field_list[%d]' % (field_cnt), field_data)
        if (len(field_data) != 2):
          logging.exception('field_list[%d] does not have 2 elements' % \
                            (field_cnt))
          raise Exception
        auxiliary.check_is_string('field_list[%d][0]' % (field_cnt),
                                  field_data[0])
        auxiliary.check_is_integer('field_list[%d][1]' % (field_cnt),
                                   field_data[1])
        auxiliary.check_is_positive('field_list[%d][1]' % (field_cnt),
                                    field_data[1])

      field_cnt += 1

    # Calculate the column starts and ends - - - - - - - - - - - - - - - - - -
    #
    self.col_start_end = []
    start_col = 0

    for field_data in self.field_list:
      if (isinstance(field_data, tuple)):
        col_width = field_data[1]
      else:
        col_width = field_data

      self.col_start_end.append((start_col, start_col+col_width))
      start_col += col_width

    auxiliary.check_is_string('file_name', self.file_name)

    # Now perform various checks for each access mode and open file - - - - - -
    #
    if (self.access_mode == 'read'):

      if (self.file_name.endswith('.gz')) or (self.file_name.endswith('.GZ')):
        file_gzipped = True

        try:
          self.file = gzip.open(self.file_name) # Open gzipped file
        except:
          logging.exception('Cannot open gzipped CSV file "%s" for reading' % \
                            (self.file_name))
          raise IOError

      else:  # Open normal file for reading
        file_gzipped = False

        try:  # Try to open the file in read mode
          self.file = open(self.file_name,'r')
        except:
          logging.exception('Cannot open CSV file "%s" for reading' % \
                            (self.file_name))
          raise IOError

      # If header line is set to True get field names from file
      #
      if (self.header_line == True):
        header_line = self.file.readline()

        new_field_list = []

        field_cnt = 0
        for (s,e) in self.col_start_end:  # Extract fields
          field_name = header_line[s:e]  # Extract a field name

          if (self.strip_fields == True):
            field_name = field_name.strip()

          new_field_list.append((field_name, self.field_list[field_cnt]))
          field_cnt += 1

        self.field_list = new_field_list  # Store back

      # Count number of records in the file
      #
      if ((sys.platform[0:5] in ['linux','sunos']) and \
          (file_gzipped == False)):  # Fast line counting
        wc = os.popen('wc -l ' + self.file_name)
        num_rows = int(string.split(wc.readline())[0])
        wc.close()
      else:  # Slow line counting method

        num_rows = 0

        if (file_gzipped == True):
          fp = gzip.open(self.file_name)
        else:
          fp = open(self.file_name,'r')
        for l in fp:
          num_rows += 1
        fp.close()

      self.num_records = num_rows
      if (self.header_line == True):
        self.num_records -= 1

      # Check that there are records in the data set
      #
      if (self.num_records == 0):
        logging.exception('No records in COL data set opened for reading')
        raise Exception

      self.next_rec_num = 0

    elif (self.access_mode == 'write'):   # - - - - - - - - - - - - - - - - - -

      # Try to open the file in write mode
      #
      try:
        self.file = open(self.file_name,'w')
      except:
        logging.exception('Cannot open COL file "%s" for writing' % \
                          (self.file_name))
        raise IOError

      # Write the header line with field names if desired
      #
      if (self.write_header == True):
        header_line = ''
        for (field_name, col_width) in self.field_list:
          if (self.strip_fields == True):
            field_name = field_name.strip()

          header_line += field_name[:col_width].ljust(col_width)

        self.file.write(header_line+os.linesep)
        self.file.flush()

      self.num_records =  0
      self.next_rec_num = 0

    elif (self.access_mode == 'append'):  # - - - - - - - - - - - - - - - - - -

      # Try to open the file in append mode
      #
      try:
        self.file = open(self.file_name,'a')
      except:
        logging.exception('Cannot open CSV file "%s" for appending' % \
                          (self.file_name))
        raise IOError

      # Count number of records in the file
      #
      if (sys.platform[0:5] in ['linux','sunos']):  # Fast line counting
        wc = os.popen('wc -l ' + self.file_name)
        num_rows = int(string.split(wc.readline())[0])
        wc.close()
      else:  # Slow line counting method
        num_rows = 0
        fp = open(self.file_name,'r')
        for l in fp:
          num_rows += 1
        fp.close()

      # If no records are stored write header if desired
      #
      if (num_rows == 0) and (self.write_header == True):
        header_line = ''
        for (field_name,col_width) in self.field_list:
          if (self.strip_fields == True):
            field_name = field_name.strip()

          header_line += field_name[:col_width].ljust(col_width)

        self.file.write(header_line+os.linesep)
        self.file.flush()

        self.num_records =  0
        self.next_rec_num = 0

      else:  # There are records in the data set already

        self.num_records = num_rows
        if (self.header_line == True):
          self.num_records -= 1

        self.next_rec_num = self.num_records

    else:  # Illegal data set access mode - - - - - - - - - - - - - - - - - - -

      logging.exception('Illegal data set access mode: "%s" (not allowed ' % \
            (str(self.access_mode)) + 'with COL data set implementation).')
      raise Exception

    # Check if one of the fields is the record identifier - - - - - - - - - -
    #
    field_cnt = 0

    for (field_name, field_col_width) in self.field_list:

      if (field_name == self.rec_ident):
        self.rec_ident_field = field_cnt
      field_cnt += 1

    self.log([('COL file name', self.file_name),
              ('Header line', self.header_line),
              ('Write header', self.write_header),
              ('Column start and ends', self.col_start_end),
              ('Record identifier field number', self.rec_ident_field)])

  # ---------------------------------------------------------------------------

  def finalise(self):
    """Finalise a data set, i.e. close the file and set various attributes to
       None.
    """

    # Close file if it is open
    #
    if (self.file != None):

      self.file.close()
      self.file = None

    self.access_mode =  None
    self.file_name =    None
    self.num_records =  None
    self.next_rec_num = None

    # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.info('Finalised COL data set "%s"' % (self.description))

  # ---------------------------------------------------------------------------

  def __read_one_record__(self):
    """Read and return the next record. Should not be used from outside the
       module.

       Updates 'next_rec_num' attribute.

       Returns a dictionary {record_identifier : record_field_list}, or {} once
       the end of the file is reached.
    """

    file_line = self.file.readline()

    if (file_line == ''):  # Reached end of file
      return {}

    rec = []  # Extract fields from file string

    for (s,e) in self.col_start_end:
      field_val = file_line[s:e]  # Extract a field name

      if (self.strip_fields == True):  # Strip leading and trailing whitespace
        field_val =  field_val.strip()

      if (self.miss_val != None):  # Check for missing values in record
        if (field_val in self.miss_val):
          field_val = ''  # Found a missing value

      rec.append(field_val)

    if (self.rec_ident_field == -1):  # Generate or get record identifier
      rec_ident = self.rec_ident+'-%d' % (self.next_rec_num)

    else:  # Get record identifier from the record itself
      rec_ident = rec[self.rec_ident_field]

    self.next_rec_num += 1

    return {rec_ident:rec}

  # ---------------------------------------------------------------------------

  def read(self, *recs):
    """Read and return one or more records.
       - If no argument is given return the next record in the data set.
       - If one argument (n) is given (must be a positive number) return the
         next n records (or less if there are not enough records in the data
         set).
       - If two arguments (s,n) are given (both must be numbers) return the
         next n records starting from record s (or less if there are not enough
         records in the data set).

       Returns an empty dictionary if no more records are available.
    """

    if (self.file == None):
      logging.exception('Data set not initialised')
      raise Exception

    if (self.access_mode != 'read'):
      logging.exception('Data set not initialised for "read" access')
      raise Exception

    len_recs_arg = len(recs)

    # Process the three different forms of calling
    #
    if (len_recs_arg == 0):  # No arguments, just read one record - - - - - - -

      return self.__read_one_record__()

    elif (len_recs_arg == 1):  # One argument, read n records - - - - - - - - -

      num_recs = recs[0]
      if (not isinstance(num_recs, int) or (num_recs < 1)):
        logging.exception('Number of records given is not a positive integer' \
                          + ' number: %s' % (str(num_recs)))
        raise Exception

      rec_dict = {}

      for i in range(num_recs):

        this_rec_dict = self.__read_one_record__()

        if (this_rec_dict == {}):  # No more records in data set
          return rec_dict

        if this_rec_dict.keys()[0] in rec_dict:  # Check for unique rec ident.
          logging.exception('Already existing record identifier returned: %s' \
                            % (str(this_rec_dict)))
          raise Exception

        rec_dict.update(this_rec_dict)

      return rec_dict

    elif (len_recs_arg == 2):  # Two arguments, read n records - - - - - - - -

      start_num = recs[0]
      num_recs =  recs[1]
      if (not isinstance(start_num, int) or (start_num < 0)):
        logging.exception('Start record number given is not a positive ' + \
                          'integer number: %s' % (str(start_num)))
        raise Exception
      if (not isinstance(num_recs, int) or (num_recs < 1)):
        logging.exception('Number of records given is not a positive integer' \
                          + ' number: %s' % (str(num_recs)))
        raise Exception
      if (start_num >= self.num_records):
        logging.exception('Start record number is larger than the number ' + \
                          'of records in the data set: %d' % (start_num))
        raise Exception

      # Check if the start record number is at the current position or not
      #
      if (start_num != self.next_rec_num):

        self.file.close()  # Close currently open file

        if (self.file_name.endswith('.gz')) or \
           (self.file_name.endswith('.GZ')):
          self.file = gzip.open(self.file_name) # Open gzipped file
        else:
          self.file = open(self.file_name,'r')  # Re-open file

        # Skip over header (if there is one) and skip to start record
        #
        if (self.header_line == True):
          self.file.readline()

        for i in range(start_num):
          self.file.readline()

        self.next_rec_num = start_num  # Update record counter

      # Now read the desired number of records - - - - - - - - - - - - - - - -
      #
      rec_dict = {}

      for i in range(num_recs):

        this_rec_dict = self.__read_one_record__()

        if (this_rec_dict == {}):  # No more records in data set
          return rec_dict

        if this_rec_dict.keys()[0] in rec_dict:  # Check for unique rec ident.
          logging.exception('Already existing record identifier returned: %s' \
                            % (str(this_rec_dict)))
          raise Exception

        rec_dict.update(this_rec_dict)

      return rec_dict

    else:  # Illegal calling (3 or more arguments)
      logging.exception('Illegal call to read(): 3 or more arguments: %s' % \
                        (str(recs)))
      raise Exception

  # ---------------------------------------------------------------------------

  def readall(self):
    """An iterator which will return one record per call as a tuple (record
       identifier, record field list).

       The file is first closed and then re-opened returning the first record.
    """

    if (self.file == None):
      logging.exception('Data set not initialised')
      raise Exception

    if (self.access_mode != 'read'):
      logging.exception('Data set not initialised for "read" access')
      raise Exception

    self.file.close()  # Close currently open file

    if (self.file_name.endswith('.gz')) or (self.file_name.endswith('.GZ')):
      self.file = gzip.open(self.file_name) # Open gzipped file
    else:
      self.file = open(self.file_name,'r')  # Re-open file

    self.next_rec_num = 0   # Initialise next record counter


    # Skip over header (if there is one) and skip to start record
    #
    if (self.header_line == True):
      self.file.readline()

    for file_line in self.file.xreadlines():

      rec = []  # Extract fields from file string

      for (s,e) in self.col_start_end:
        field_val = file_line[s:e]  # Extract a field name

        if (self.strip_fields == True):  # whitespace
          field_val =  field_val.strip()

        if (self.miss_val != None):  # Check for missing values in record
          if (field_val in self.miss_val):
            field_val = ''  # Found a missing value

        rec.append(field_val)

      if (self.rec_ident_field == -1):  # Generate or get record identifier
        rec_ident = self.rec_ident+'-%d' % (self.next_rec_num)

      else:  # Get record identifier from the record itself
        rec_ident = rec[self.rec_ident_field]

      self.next_rec_num += 1

      yield (rec_ident,rec)

  # ---------------------------------------------------------------------------

  def write(self, rec_dict):
    """Write one or more records into the data set.
       The input dictionary with records is first sorted (according to the
       record identifiers) and then written into the COL file.
    """

    if (self.file == None):
      logging.exception('Data set not initialised')
      raise Exception

    if (self.access_mode not in ['write','append']):
      logging.exception('Data set not initialised for "write" or "append" ' + \
                        'access')
      raise Exception

    # Sort record identifiers
    #
    rec_ident_keys = rec_dict.keys()
    rec_ident_keys.sort()

    for rec_ident in rec_ident_keys:
      rec = rec_dict[rec_ident]

      if (self.strip_fields == True):  # Strip leading and trailing whitespace
        rec = map(string.strip,rec)

      file_line = ''  # String to be written to file

      field_cnt = 0
      for (field_name, col_width) in self.field_list:

        field_val = rec[field_cnt]

        if (self.miss_val != None):  # Check for missing values in record
          if (field_val in self.miss_val):
            field_val = ''  # Found a missing value

        file_line += field_val[:col_width].ljust(col_width)
        field_cnt += 1

      self.file.write(file_line+os.linesep)

      self.num_records +=  1
      self.next_rec_num += 1

    self.file.flush()

# =============================================================================

class DataSetMemory(DataSet):
  """Implementation of a memory based data set class using Python dictionaries.

     The 'field_list' attribute must be given when a memory based data set is
     initialised. The field list must contain tuples where the first element is
     a field name and the second element can be an empty string (that will not
     be used), for example:

       field_list=[('rec-id',''),('title',''),('gname',''),('surname','')]

     The only possible value for the 'access_mode' argument is: 'readwrite'.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """No additional attributes needed, just call base class constructor.
    """

    self.dataset_type = 'MEMORY'

    self.dict = {}           # The dictionary containing the records
    self.rec_ident_col = -1  # Column of the record identifier field

    DataSet.__init__(self, kwargs)  # Process base arguments

    if (self.access_mode != 'readwrite'):
      logging.exception('Memory data set must be initialised in "readwrite" ' \
                        + ' access mode, not: "%s"' % (self.access_mode))
      raise Exception

    self.num_records =   0
    self.rec_ident_col = -1     # Column of the record identifier field

    # Check if the record identifier is one of the field names
    #
    field_col = 0

    for (field_name,not_used) in self.field_list:

      auxiliary.check_is_string('field name in column %d' % (field_col), \
                                field_name)

      # Check if this is the record identifier field
      #
      if (self.rec_ident == field_name):
        self.rec_ident_col = field_col

      field_col += 1

    self.log([('Record identifier column', self.rec_ident_col)])

  # ---------------------------------------------------------------------------

  def finalise(self):
    """Finalise a data set. All data will be lost.
    """

    self.dict.clear()  # Properly clean up memory
    self.dict = None

    self.access_mode =  None
    self.num_records =  None

    # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.info('Finalised Memory data set "%s"' % (self.description))

  # ---------------------------------------------------------------------------

  def read(self, recs):
    """Read and return one or more records.

       - If the argument is a string it is assumed to be a record identifier
         and the corresponding record (if it is in the data set) will be
         returned a a dictionary (otherwise an empty dictionary will be
         returned).
       - If the argument is a list or a set of strings (assumed to be record
         identifiers) then all corresponding records in the data set will be
         returned. if no record is found an empty dictionary will be returned.
    """

    if (self.dict == None):
      logging.exception('Data set not initialised')
      raise Exception

    if (isinstance(recs, str)):  # One record identifier only - - - - - - - - -

      if (recs in self.dict):
        rec = self.dict[recs]

        if (self.strip_fields == True):  # Strip whitespace
          rec = map(string.strip,rec)

        if (self.miss_val != None):  # Check for missing values in record
          clean_rec = []
          miss_val_list = self.miss_val  # Faster reference access

          for val in rec:
            if (val in miss_val_list):  # Found a missing value
              clean_rec.append('')  # Replace with empty string
            else:
              clean_rec.append(val)
          rec = clean_rec

        return {recs:rec}

      else:
        return {}

    # List or set of record identifiers - - - - - - - - - - - - - - - - - - - -
    #
    elif (isinstance(recs, list) or isinstance(recs, set) or \
          isinstance(recs, set)):

      rec_dict = {}

      for rec_ident in recs:

        if (not isinstance(rec_ident, str)):
          logging.exception('Record identifier is not a string: "%s"' % \
                            (str(rec_ident)))
          raise Exception

        if (rec_ident in self.dict) and (rec_ident not in rec_dict):

          rec = self.dict[rec_ident]  # Get the record as a list of values

          if (self.strip_fields == True):  # Strip whitespace
            rec = map(string.strip,rec)

          if (self.miss_val != None):  # Check for missing values in record
            clean_rec = []
            miss_val_list = self.miss_val  # Faster reference access

            for val in rec:
              if (val in miss_val_list):  # Found a missing value
                clean_rec.append('')  # Replace with empty string
              else:
                clean_rec.append(val)
            rec = clean_rec

          rec_dict[rec_ident] = rec

      return rec_dict

    else:
      logging.exception('Illegal argument given to read(): "%s" of type %s' % \
                        (str(recs), type(recs)))
      raise Exception

  # ---------------------------------------------------------------------------

  def readall(self):
    """An iterator which will return one record per call as a tuple (record
       identifier, record field list).

       It returns the records in an arbitrary order (not sorted according to
       record identifiers).
    """

    if (self.dict == None):
      logging.exception('Data set not initialised')
      raise Exception

    for rec_key in self.dict.keys():
      rec = self.dict[rec_key]

      if (self.strip_fields == True):  # Strip leading and trailing whitespace
        rec = map(string.strip, rec)

      if (self.miss_val != None):  # Check for missing values in record
        clean_rec = []
        miss_val_list = self.miss_val  # Faster reference access

        for val in rec:
          if (val in miss_val_list):  # Found a missing value
            clean_rec.append('')  # Replace with empty string
          else:
            clean_rec.append(val)
        rec = clean_rec

      if (self.rec_ident_col == -1):  # Use the dictionary key
        rec_ident = rec_key

      else:  # Get record identifier from the record itself
        rec_ident = rec[self.rec_ident_col]

      yield (rec_ident,rec)

  # ---------------------------------------------------------------------------

  def write(self, rec_dict):
    """Write one or more records into the data set.

       The given records are simply stored into the dictionary, but record keys
       (identifiers) are checked for duplicates - if found warnings are logged.
    """

    if (self.dict == None):
      logging.exception('Data set not initialised')
      raise Exception

    for rec_ident in rec_dict:
      if (rec_ident in self.dict):
        logging.warn('Record with identifer "%s" is already in the memory ' % \
                     (rec_ident)+'data set - overwrite old version.')
      else:
        self.num_records += 1  # This is a new record

      rec = rec_dict[rec_ident]  # Get record as a list of values

      if (self.strip_fields == True):  # Strip leading and trailing whitespace
        rec = map(string.strip, rec)

      if (self.miss_val != None):  # Check for missing values in record
        clean_rec = []
        miss_val_list = self.miss_val  # Faster reference access

        for val in rec:
          if (val in miss_val_list):  # Found a missing value
            clean_rec.append('')  # Replace with empty string
          else:
            clean_rec.append(val)
        rec = clean_rec

      self.dict[rec_ident] = rec

# =============================================================================

class DataSetShelve(DataSet):
  """Implementation of a disk based data set class using Python shelves.

     The 'field_list' attribute must be given when a memory based data set is
     initialised. The field list must contain tuples where the first element is
     a field name and the second element can be an empty string (that will not
     be used), for example:

       field_list=[('rec-id',''),('title',''),('gname',''),('surname','')]

     The only possible value for the 'access_mode' argument is: 'readwrite'.

     The additional arguments (besides the base class arguments) which have to
     be set when this data set is initialised are:

       file_name  A string containing the name of the underlying shelve file.
       clear      A flag (True or False), when True the content of the shelve
                  database file will be cleared when opened. Default value is
                  False.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the derived attributes first, then call the base
       class constructor.
    """

    self.dataset_type = 'SHELVE'

    self.file_name = None   # The name of the shelve file (without extensions)
    self.clear =     False  # Flag (True or False) for clearing the database
                            # when opening or not

    self.shelve =        None  # The 'shelve'fFile pointer to current file
    self.db =            None  # A reference to the underlying databas
    self.rec_ident_col = -1    # Column of the record identifier field

    # Process all keyword arguments
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('file')):
        auxiliary.check_is_string('file_name', value)
        self.file_name = value

      elif (keyword.startswith('cle')):
        auxiliary.check_is_flag('clear', value)
        self.clear = value

      else:
        base_kwargs[keyword] = value

    DataSet.__init__(self, base_kwargs)  # Process base arguments

    # Make sure the 'file_name' attribute is set and 'access_mode' is correct -
    #
    auxiliary.check_is_string('file_name', self.file_name)

    if (self.access_mode != 'readwrite'):
      logging.exception('Shelve data set must be initialised in "readwrite" ' \
                        + ' access mode, not: "%s"' % (self.access_mode))
      raise Exception

    # Check if the record identifier is one of the field names
    #
    field_col = 0

    for (field_name,not_used) in self.field_list:

      auxiliary.check_is_string('field name in column %d' % (field_col), \
                                field_name)

      # Check if this is the record identifier field
      #
      if (self.rec_ident == field_name):
        self.rec_ident_col = field_col

      field_col += 1

    # Now open the shelve - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    try:
      self.shelve = shelve.open(self.file_name)
    except:
      logging.exception('Cannot open shelve: "%s"' % str(self.file_name))
      raise Exception

    # Now clear the shelve if the 'clear' flag is set to True - - - - - - - - -
    #
    if (self.clear == True):

      for rec_key in self.shelve:  # Remove all entries

        del self.shelve[rec_key]

      self.shelve.sync()  # And make sure the database is updated

    # Get the number of records in the shelve - - - - - - - - - - - - - - - - -
    #
    self.num_records = len(self.shelve)

    self.log([('Shelve file name', self.file_name),
              ('Clear flag', self.clear),
              ('Record identifier column', self.rec_ident_col)])

  # ---------------------------------------------------------------------------

  def finalise(self):
    """Finalise a data set. Close the shelve file.
    """

    # Close shelve if it is open
    #
    if (self.shelve != None):

      self.shelve.close()

      self.shelve = None

    self.access_mode = None
    self.file_name =   None
    self.num_records = None

    # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.info('Finalised shelve data set "%s"' % (self.description))

  # ---------------------------------------------------------------------------

  def read(self, recs):
    """Read and return one or more records.

       - If the argument is a string it is assumed to be a record identifier
         and the corresponding record (if it is in the data set) will be
         returned a a dictionary (otherwise an empty dictionary will be
         returned).
       - If the argument is a list or a set of strings (assumed to be record
         identifiers) then all corresponding records in the data set will be
         returned. if no record is found an empty dictionary will be returned.
    """

    if (self.shelve == None):
      logging.exception('Data set not initialised')
      raise Exception

    if (isinstance(recs, str)):  # One record identifier only - - - - - - - - -

      if (recs in self.shelve):
        rec = self.shelve[recs]

        if (self.strip_fields == True):  # Strip whitespace
          rec = map(string.strip,rec)

        if (self.miss_val != None):  # Check for missing values in record
          clean_rec = []
          miss_val_list = self.miss_val  # Faster reference access

          for val in rec:
            if (val in miss_val_list):  # Found a missing value
              clean_rec.append('')  # Replace with empty string
            else:
              clean_rec.append(val)
          rec = clean_rec

        return {recs:rec}

      else:
        return {}

    # List or set of record identifiers - - - - - - - - - - - - - - - - - - - -
    #
    elif (isinstance(recs, list) or isinstance(recs, set) or \
          isinstance(recs, set)):

      rec_dict = {}

      for rec_ident in recs:

        if (not isinstance(rec_ident, str)):
          logging.exception('Record identifier is not a string: "%s"' % \
                            (str(rec_ident)))
          raise Exception

        if (rec_ident in self.shelve) and (rec_ident not in rec_dict):

          rec = self.shelve[rec_ident]  # Get the record as a list of values

          if (self.strip_fields == True):  # Strip whitespace
            rec = map(string.strip,rec)

          if (self.miss_val != None):  # Check for missing values in record
            clean_rec = []
            miss_val_list = self.miss_val  # Faster reference access

            for val in rec:
              if (val in miss_val_list):  # Found a missing value
                clean_rec.append('')  # Replace with empty string
              else:
                clean_rec.append(val)
            rec = clean_rec

          rec_dict[rec_ident] = rec

      return rec_dict

    else:
      logging.exception('Illegal argument given to read(): "%s" of type %s' % \
                        (str(recs), type(recs)))
      raise Exception

  # ---------------------------------------------------------------------------

  def readall(self):
    """An iterator which will return one record per call as a tuple (record
       identifier, record field list).

       It returns the records in an arbitrary order (not sorted according to
       record identifiers).
    """

    if (self.shelve == None):
      logging.exception('Data set not initialised')
      raise Exception

    for rec_key in self.shelve:
      rec = self.shelve[rec_key]

      if (self.strip_fields == True):
        rec = map(string.strip, rec)

      if (self.miss_val != None):  # Check for missing values in record
        clean_rec = []
        miss_val_list = self.miss_val  # Faster reference access

        for val in rec:
          if (val in miss_val_list):  # Found a missing value
            clean_rec.append('')  # Replace with empty string
          else:
            clean_rec.append(val)
        rec = clean_rec

      if (self.rec_ident_col == -1):  # Use the dictionary key
        rec_ident = rec_key

      else:  # Get record identifier from the record itself
        rec_ident = rec[self.rec_ident_col]

      yield (rec_ident,rec)

  # ---------------------------------------------------------------------------

  def write(self, rec_dict):
    """Write one or more records into the data set.

       The given records are simply stored into the shelve, but record keys
       (identifiers) are checked for duplicates - if found warnings are logged.
    """

    if (self.shelve == None):
      logging.exception('Data set not initialised')
      raise Exception

    for rec_ident in rec_dict:
      if (rec_ident in self.shelve):
        logging.warn('Record with identifer "%s" is already in the memory ' % \
                     (rec_ident)+'data set - overwrite old version.')
      else:
        self.num_records += 1  # This is a new record

      rec = rec_dict[rec_ident]  # Get record as a list of values

      if (self.strip_fields == True):  # Strip leading and trailing whitespace
        rec = map(string.strip, rec)

      if (self.miss_val != None):  # Check for missing values in record
        clean_rec = []
        miss_val_list = self.miss_val  # Faster reference access

        for val in rec:
          if (val in miss_val_list):  # Found a missing value
            clean_rec.append('')  # Replace with empty string
          else:
            clean_rec.append(val)
        rec = clean_rec

      self.shelve[rec_ident] = rec

    self.shelve.sync()  # And make sure the database is updated

# =============================================================================
