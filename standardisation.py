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
# The Original Software is: "standardisation.py"
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

"""Module standardisation.py - Classes for cleaning and standardisations.

   This module provides classes for record cleaning and standardisations,
   either based on rules or a machine learning approach (Hidden Markov models)

   TODO
   - PC 29/11/2002 Check for missing values in date fields (data standardiser)
                   And how to deal with them?
"""

# =============================================================================
# Import necessary modules (Python standard modules first, then Febrl modules)

import logging
import os
import string
import time

import auxiliary
import mymath
import phonenum

# =============================================================================

class RecordStandardiser:
  """Main class for the standardisation process. Implements a method to clean
     and standardise all records in a data set.

     A record standardiser will include one or more component standardisers.
     Currently four different types of components can be standardised: names,
     addresses, date and phone numbers. A list of component standardisers needs
     to be given when a record standardiser is initialised.

     A record standardiser needs the following instance variables to be set
     when it is initialised:

       input_dataset    Reference to the input data set. Has to be initalised
                        in read access mode.
       output_dataset   Reference to the output data set. Has to be initalised
                        in write access mode.
       comp_stand_list  List with one or more component standardisers.
       log_funct        This can be a Python function or method which will log
                        (print or save to a file) a progress report message.
                        It is assumed that this function or method has one
                        input argument of type string (the message to be
                        printed). Default is None, in which case the normal
                        logging module will be used.
       progress_report  Can be set to a percentage number between 1 and 50 in
                        which case a progress report is logged during the
                        record pair comparison stage (in the run() method)
                        every selected percentage number. Default values is
                        10. If set to None no progress report will be logged.
       pass_field_list  A list of tuples (input_field, output_field) that will
                        simply be passed from the input data set to the output
                        data set without being cleaned or standardised. Each
                        input_field must be in the field list of the input
                        data set and correspondingly each output_filed in the
                        list of field of the output data set. Default value for
                        this argument is an empty list.

     The standardise() method can be used to standardise all records from the
     input data set and write them into the output data set.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor.
    """

    self.description =     ''
    self.inp_dataset =     None
    self.out_dataset =     None
    self.comp_stand_list = []    # A list of one or more component standardiser
                                 # tuples (each made of the component
                                 # standardiser itself, a list with the field
                                 # indices from the input data set, and a list
                                 # of the field indices from the output data
                                 # set. Only the standardisers themselves will
                                 # be provided as input, the field index lists
                                 # will be generated from them during
                                 # initialisation.
    self.progress_report = 10
    self.log_funct =       None
    self.pass_field_list = []

    # The indices of the pass fields into the input and output data sets
    #
    self.pass_field_index_list = []

    # Process all keyword arguments
    #
    for (keyword, value) in kwargs.items():
      if (keyword.startswith('desc')):
        auxiliary.check_is_string('description', value)
        self.description = value

      elif (keyword.startswith('input_d')):
        auxiliary.check_is_list('Input dataset field list', value.field_list)
        self.in_dataset = value
      elif (keyword.startswith('output_d')):
        auxiliary.check_is_list('Output dataset field list', value.field_list)
        self.out_dataset = value

      elif (keyword.startswith('comp_stand')):
        auxiliary.check_is_list('comp_stand_list', value)
        self.comp_stand_list = value

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

      elif (keyword.startswith('pass_fi')):
        auxiliary.check_is_list('pass_field_list', value)
        self.pass_field_list = value

      else:
        logging.exception('Illegal constructor argument keyword: '+keyword)
        raise Exception

    # Check if the needed attributes are set  - - - - - - - - - - - - - - - - -
    #
    auxiliary.check_is_not_none('Input dataset', self.in_dataset)
    auxiliary.check_is_not_none('Output dataset', self.out_dataset)
    auxiliary.check_is_not_none('Component standardiser list',
                                self.comp_stand_list)

    if (self.in_dataset.access_mode != 'read'):
      logging.exception('Input data set "%s" must be in "read" access mode' % \
                        (in_dataset.description))
      raise Exception

    if (self.out_dataset.access_mode != 'write'):
      logging.exception('Output data set "%s" must be in "write" access ' % \
                        (out_dataset.description) + 'mode')
      raise Exception

    # Get a list of the fields from the input and output data sets
    #
    in_dataset_fields = []
    for f_data in self.in_dataset.field_list:
      in_dataset_fields.append(f_data[0])  # Append field name only

    out_dataset_fields = []
    for f_data in self.out_dataset.field_list:
      out_dataset_fields.append(f_data[0])

    # Check if there is no conflict in the output fields definition - - - - - -
    #
    output_fields_set = set()

    for cs in self.comp_stand_list:

      for field in cs.out_fields:
        if (field != None):  # Only check fields that are not set to None
          if (field in output_fields_set):
            logging.exception('Output fields definition conflict with ' + \
                              'field "%s"' % (str(field)))
            raise Exception
          output_fields_set.add(field)

    # Check that fields in pass field list are in input and output data sets -
    #
    pass_out_field_set = set()  # Check for duplicate output field names

    i = 0
    for field_tuple in self.pass_field_list:
      auxiliary.check_is_tuple('pass_field_list[%d]' % (i), field_tuple)
      if (len(field_tuple) != 2):
        logging.exception('Field tuple in "pass_field_list" does not ' + \
                          'contain two fields: %s' % (str(field_tuple)))
        raise Exception

      in_field, out_field = field_tuple

      if (in_field not in in_dataset_fields):
        logging.exception('Input field in "pass_field_list" is not in the ' + \
                          'input data set: %s' % (str(field_tuple)))
        raise Exception

      if (out_field not in out_dataset_fields):
        logging.exception('Output field in "pass_field_list" is not in the' + \
                          ' output data set: %s' % (str(field_tuple)))
        raise Exception

      if (out_field in pass_out_field_set):
        logging.exception('Output field in "pass_field_list" appears twice' + \
                          ' field tuples: %s' % (str(self.pass_field_list)))
        raise Exception

      pass_out_field_set.add(out_field)

      # Check for conflicts with component standardiser output fields
      #
      if (out_field in output_fields_set):
        logging.exception('Output field in "pass_field_list" is also an ' + \
                          'output field of a component standardiser: %s' % \
                          (str(self.pass_field_list)))
        raise Exception

      in_pass_field_ind =  in_dataset_fields.index(in_field)
      out_pass_field_ind = out_dataset_fields.index(out_field)
      self.pass_field_index_list.append((in_pass_field_ind,out_pass_field_ind))

      i += 1

    # Check all component standardisers and generate their lists of field - - -
    # indices to be used.
    #
    for i in range(len(self.comp_stand_list)):
      cs = self.comp_stand_list[i]
      in_field_indices = []

      in_not_none = 0
      for in_f in cs.in_fields:
        if (in_f != None):
          in_not_none += 1
          if in_f not in in_dataset_fields:
            logging.exception('Component standardiser "%s" contains a ' % \
                              (cs.description) + 'field that is not in ' + \
                              'the input data set: %s' % (in_f))
            raise Exception

          in_field_indices.append(in_dataset_fields.index(in_f))
        else:
          in_field_indices.append(None)

      if (in_not_none == 0):
        logging.exception('No input field, or only None input fields ' + \
                          'defined for component standardiser "%s"' % \
                          (cs.description))
        raise Exception

      # Now the same for output fields
      #
      out_field_indices = []

      out_not_none = 0
      for out_f in cs.out_fields:
        if (out_f != None):
          out_not_none += 1
          if out_f not in out_dataset_fields:
            logging.exception('Component standardiser "%s" contains a ' % \
                              (cs.description) + 'field that is not in the' + \
                              ' output data set: %s' % (out_f))
            raise Exception

          out_field_indices.append(out_dataset_fields.index(out_f))
        else:
          out_field_indices.append(None)

      if (out_not_none == 0):
        logging.exception('No output field, or only None output fields' + \
                          'defined for component standardiser "%s"' % \
                          (cs.description))
        raise Exception

      self.comp_stand_list[i] = (cs, in_field_indices, out_field_indices)

    # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.info('Initialised record standardiser: %s' % (self.description))
    logging.info('  Input data set:  %s' % (self.in_dataset.description))
    logging.info('  Output data set: %s' % (self.out_dataset.description))
    if (self.progress_report == None):
      logging.info('  No progress reported.')
    else:
      logging.info('  Progress report every:  %d%%' % (self.progress_report))
    logging.info('  Pass field list: %s' % (str(self.pass_field_list)))

    logging.info('  Component standardisers:')
    for cs in self.comp_stand_list:
      logging.info('    Name: %s' % (str(cs[0].description)))
      logging.info('      Input fields:         %s' % (str(cs[0].in_fields)))
      logging.info('      Input field indices:  %s' % (str(cs[1])))
      logging.info('      Output fields:        %s' % (str(cs[0].out_fields)))
      logging.info('      Output field indices: %s' % (str(cs[2])))

  # ---------------------------------------------------------------------------

  def standardise(self):
    """Standardise all records in the input data set and write them into the
       output data set.
    """

    num_output_fields = len(self.out_dataset.field_list)

    # Calculate a counter for the progress report
    #
    if (self.progress_report != None):
      progress_report_cnt = int(self.in_dataset.num_records / \
                                (100.0 / self.progress_report))
      progress_report_cnt = max(1, progress_report_cnt)  # Make it positive

    else:  # So no progress report is being logged
      progress_report_cnt = in_dataset.num_records + 1

    start_time = time.time()

    rec_read = 0  # Number of records read from data set

    # Loop over all records from input data set - - - - - - - - - - - - - - - -
    #
    for (rec_ident, in_rec) in self.in_dataset.readall():

      out_rec = ['']*num_output_fields # Create output record with empty fields

      # First copy pass fields from input to output record - - - - - - - - - -
      #
      for (in_field_ind, out_field_ind) in self.pass_field_index_list:
        out_rec[out_field_ind] = in_rec[in_field_ind]

      # Process one component standardiser after the other - - - - - - - - - -
      #
      for cs_details in self.comp_stand_list:
        cs =                   cs_details[0]
        in_field_index_list =  cs_details[1]
        out_field_index_list = cs_details[2]

        num_input_fields = len(in_field_index_list)

        # First get the input record field values
        #
        in_rec_val_list = []
        for in_index in in_field_index_list:
          if (in_index != None):
            in_rec_val_list.append(in_rec[in_index])

        # Check for word spilling if more than one field - - - - - - - - - - -
        #
        if ((cs.check_word_spill == True) and (num_input_fields > 1)):

          in_str = ''  # The input string for the component standardiser

          for in_val in in_rec_val_list:
            spill_flag = cs.check_field_spill(in_str, in_val)
            if (spill_flag == True):
              in_str = in_str + in_val
            else:  # Use given field separator
              in_str = in_str + cs.field_sep + in_val
        else:
          in_str = cs.field_sep.join(in_rec_val_list)

        # Clean the input string - - - - - - - - - - - - - - - - - - - - - - -
        #
        clean_in_str = cs.clean_component(in_str)

        out_field_list = cs.standardise(in_str, clean_in_str)
        assert len(out_field_list) == len(out_field_index_list), \
               (clean_in_str, out_field_list)

        i = 0
        for out_index in out_field_index_list:
          if (out_index != None):
            out_rec[out_index] = out_field_list[i]
          i += 1

      # Write the standardised record into the output data set
      #
      self.out_dataset.write({rec_ident:out_rec})

      rec_read += 1

      if ((rec_read % progress_report_cnt) == 0):  # Log progress - - - - - - -
        used_time = time.time() - start_time
        perc_done = 100.0 * rec_read / self.in_dataset.num_records
        rec_time  = used_time / rec_read  # Time per record read and indexed
        togo_time = (self.in_dataset.num_records - rec_read) * rec_time

        used_sec_str = auxiliary.time_string(used_time)
        rec_sec_str =  auxiliary.time_string(rec_time)
        togo_sec_str = auxiliary.time_string(togo_time)

        log_str = 'Read and standardised %d of %d records (%d%%) in %s (%s' % \
                  (rec_read, self.in_dataset.num_records, round(perc_done),
                  used_sec_str, rec_sec_str) + \
                  ' per record), estimated %s until finished.' % (togo_sec_str)
        logging.info(log_str)

        if (self.log_funct != None):
          self.log_funct(log_str)

        memory_usage_str = auxiliary.get_memory_usage()
        if (memory_usage_str != None):
          logging.info('    '+memory_usage_str)

    # Finalise all component standardisers() - - - - - - - - - - - - - - - - -
    #
    for cs_details in self.comp_stand_list:
      cs_details[0].finalise()

    used_sec_str = auxiliary.time_string(time.time()-start_time)
    rec_time_str = auxiliary.time_string((time.time()-start_time) / \
                                         self.in_dataset.num_records)
    logging.info('Read and standardised %d records in %s (%s per record)' % \
                 (self.in_dataset.num_records, used_sec_str, rec_time_str))
    logging.info('')

# =============================================================================

class ComponentStandardiser:
  """Base class for component standardisers.

     The following arguments must be given to the constructor of the base class
     and all derived classes:

       input_fields   A list containing field names (that must be in the input
                      data set)
       output_fields  A list containing field names (that must be in the output
                      data set)

     The following optional arguments can be used with all component
     standardisers:

       description       A string describing the component standardiser.
       field_separator   A character that is used to join together field values
                         in a component before cleaning and standardisation is
                         done. Default values is the empty string ''.
       check_word_spill  A Flag, if set to true and a component is made of
                         several fields and the field_separator is set to a
                         whitespace character (' ') then the field spill method
                         will be called. Default is False.
       tag_table         Reference to a tag lookup table.
       corr_list         Reference to a correction list.

      Note that for date and phone number standardisers the field separator
      will automatically be set to the empty string '' and word spilling will
      be set to False.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, base_kwargs):
    """Constructor
    """

    # General attributes for all component standardisers.
    #
    self.description = ''
    self.in_fields =   None  # A string with a list of field names form the
                             # input data set, that will be standardised. The
                             # values in the fields will be concatenated into
                             # one string with the field separator string
                             # in between the values.
    self.out_fields = None   # A list of field names form the output data set
                             # into which the standardised component values
                             # will be stored.
    self.corr_list =  []     # A correction list for the cleaning step.
    self.tag_table =  {}     # A lookup table containing tags and values with
                             # their corrections

    self.field_sep =  ''           # Field separator character
    self.check_word_spill = False  # Flag for checking word spilling

    self.alphanum = string.letters+string.digits  # For word spill method

    # Process base keyword arguments (all data set specific keywords were
    # processed in the derived class constructor)
    #
    for (keyword, value) in base_kwargs.items():
      if (keyword.startswith('desc')):
        auxiliary.check_is_string('description', value)
        self.description = value

      elif (keyword.startswith('input_f')):
        auxiliary.check_is_list('input_fields', value)
        self.in_fields = value

      elif (keyword.startswith('output_f')):
        auxiliary.check_is_list('output_fields', value)
        self.out_fields = value

      elif (keyword.startswith('field_s')):
        auxiliary.check_is_string('field_separator', value)
        self.field_sep = value

      elif (keyword.startswith('check_wo')):
        auxiliary.check_is_flag('check_word_spill', value)
        self.check_word_spill = value

      elif (keyword.startswith('corr_l')):
        auxiliary.check_is_list('corr_list', value)
        self.corr_list = value

      elif (keyword.startswith('tag_t')):
        auxiliary.check_is_dictionary('tag_table', value)
        self.tag_table = value

      else:
        logging.exception('Illegal constructor argument keyword: %s' % \
                          (str(keyword)))
        raise Exception

    # Check if fields are defined and not empty - - - - - - - - - - - - - - - -
    #
    auxiliary.check_is_list('input_fields', self.in_fields)
    auxiliary.check_is_list('output_fields', self.out_fields)

    if (len(self.in_fields) == 0) or (len(self.out_fields) == 0):
      logging.exception('At least one of the field lists is empty.')
      raise Exception

    if (self.field_sep == ''):  # No word spill checking with empty separator
      self.check_word_spill = False

  # ---------------------------------------------------------------------------

  def clean_component(self, in_str):
    """Clean the given input string using a correction list.

       This method cleans the input string by checking if any of the 'original'
       strings given in the correction list is in the string. If so it will
       replace such an original string with the corresponding replacement
       string. It also strips off all leading and trailing spaces and makes all
       letters lowercase. A cleaned string is returned.
    """

    if (in_str.strip() == ''):  # Check if the string only contains whitespaces
      return ''

    tmp_str = in_str.replace('\t',' ')  # Replace tabs with whitespaces

    # Make all lowercase and add leading and trailing whitespaces  - - - - - -
    # (this is to make sure replacement strings do match at beginning and end)
    #
    tmp_str = ' '+tmp_str.lower()+' '


    for (org, repl) in self.corr_list:  # Check for strings in correction list
      if (org in tmp_str):
        tmp_str = tmp_str.replace(org, repl)

    # Make sure commas are separated from words so they become list elements  -
    #
    tmp_str = tmp_str.replace(',', ' , ')

    while ('  ' in tmp_str):  # Make sure there are no repeated spaces
      tmp_str = tmp_str.replace('  ', ' ')

    out_str = tmp_str.strip()

    return out_str

  # --------------------------------------------------------------------------

  def check_field_spill(self, str1, str2):
    """Method to check if a known word is spilling from one string into
       another.

       Uses the keys the tag look-up table to check for known words.

       Return 'True' if a word spills over, otherwise 'False'
    """

    org_str1 = str1  # Keep copies of the original input string
    org_str2 = str2

    if ((str1 == '') or (str2 == '')):  # One string is empty
      return False

    # Check if the characters at the 'boundary' are either letters or digits
    #
    if (str1[-1] in self.alphanum) and (str2[0] in self.alphanum):

      min_ind = -len(str1)
      ind =-1  # Start with last character
      while (str1[ind] in self.alphanum) and (ind > min_ind):
        ind -= 1
      test_str1 = str1[ind:]  # Only keep last word in string 1

      max_ind = len(str2)
      ind =0  # Start with first character
      while (str1[ind] in self.alphanum) and (ind < max_ind):
        ind += 1
      test_str2 = str2[:ind]  # Only keep first word in string 2

      # Concatenate into one word and make a tuple so it conforms to tag table
      # keys
      #
      check_tuple = (test_str1.lower() + test_str2.lower(),)

      if (check_tuple in self.tag_table):  # Found a word spilling
        logging.info('  Found word spilling: "%s","%s" -> "%s"' % \
                        (org_str1, org_str2, check_word))
        return True

    else:
      return False

  # ---------------------------------------------------------------------------

  def standardise(self, in_str, clean_in_str):
    """Standardise the given input string.

       Will return a list with the standardised fields corresponding to the
       list of output fields defined in the component standardiser.

       See implementations in derived classes for details.
    """

    logging.exception('Override abstract method in derived class')
    raise Exception

  # ---------------------------------------------------------------------------

  def __tag_component__(self, in_str):
    """Tag an input component string using the tag look-up table and make it a
       list.

       This method cleans the input string and extracts words, numbers and
       separators into a list. Each element of this list is assigned one or
       more tags. A 'greedy tagger' is applied, which cheques sequences of list
       elements in the given lookup table (longer sequences first) and replaces
       them with the string and tag from the lookup table if found.

       Returns two lists:
         token_list  Contains the tokens (words, numbers, etc.) from the input
                     string separated at whitespace characters.
         tag_list    One or more tags for each of the tokens in the token list.
    """

    # First, split input string into elements at spaces - - - - - - - - - - - -
    #
    org_list = in_str.split()  # The original list from the input string

    token_list = []  # The initially empty list of tokens
    tag_list  =  []  # The initially empty list of tags

    max_key_len = self.tag_table.max_key_length

    while (org_list != []):  # As long as not all elements have been processed
      tmp_list = org_list[:max_key_len]
                                                  # Extract longest sub-list
      tmp_val = []  # Start with empty value
      tmp_key = tuple(tmp_list)

      while (tmp_key != ()): # As long as key not empty and not found in lookup
        if (self.tag_table.has_key(tmp_key)):
          tmp_val = self.tag_table[tmp_key]
          break
        tmp_key = tmp_key[:-1]  # Remove last element in key

      if (tmp_val != []):  # A value has been found in the dictionary - - - - -
        tmp_len = len(tmp_key)  # Length of found sequence

        if (tmp_val[0] != ''):  # It's not an empty value
          token_list.append(tmp_val[0])  # Append corrected word (or sequence)
          tag_list.append(tmp_val[1])    # Append tag or tags

      else:  # No value has been found in the lookup dictionary, try other tags

        tmp_val = org_list[0]  # Value is first element in the original list
        tmp_len = 1

        if (len(tmp_val) == 1) and (tmp_val.isalpha()):  # A 1-letter word - -
          token_list.append(tmp_val)
          tag_list.append('II')

        elif (tmp_val.isdigit()):  # Element is a number - - - - - - - - - - -
          token_list.append(tmp_val)
          if (len(tmp_val) == 4):
            tag_list.append('N4')
          else:
            tag_list.append('NU')

        elif (not tmp_val.isalpha()) and tmp_val.isalnum():  # Alpha-numeric -
          token_list.append(tmp_val)
          tag_list.append('AN')

        elif (tmp_val == '-'):  # Element is a hyphen - - - - - - - - - - - - -

          # Only append if previous element is not a hyphen, vertical bar, or
          # brackets
          #
          if ((token_list == []) or \
              (tag_list[-1] not in ['HY','VB','OB','CB'])):
            token_list.append(tmp_val)
            tag_list.append('HY')

        elif (tmp_val == ','):  # Element is a comma - - - - - - - - - - - - -
          pass  # Don't append

        elif (tmp_val == '|'):  # Element is a vertical bar - - - - - - - - - -

          # Remove a pair of vertical bars
          #
          if ((token_list != []) and (tag_list[-1] == 'VB')):
            tag_list =   tag_list[:-1]
            token_list = token_list[:-1]

          # Replace previous hyphen a vertical bar
          #
          elif ((token_list != []) and (tag_list[-1] == 'HY')):
            token_list[-1] = tmp_val
            tag_list[-1] = 'VB'

          else:  # Simple append vertical bar
            token_list.append(tmp_val)
            tag_list.append('VB')

        elif (tmp_val == '('):  # Element is an opening bracket - - - - - - - -

          # Replace previous hyphen with an opening bracket
          #
          if ((token_list != []) and (tag_list[-1] == 'HY')):
            token_list[-1] = tmp_val
            tag_list[-1] = 'OB'

          # Only append if previous element is not an opening bracket
          #
          elif ((token_list == []) or (tag_list[-1] != 'OB')):
            token_list.append(tmp_val)
            tag_list.append('OB')

        elif (tmp_val == ')'):  # Element is an closing bracket - - - - - - - -

          # Remove a pair of opening / closing brackets
          #
          if ((token_list != []) and (tag_list[-1] == 'OB')):
            tag_list =   tag_list[:-1]
            token_list = token_list[:-1]

          # Replace previous hyphen with a closing bracket
          #
          elif ((token_list != []) and (tag_list[-1] == 'HY')):
            token_list[-1] = tmp_val
            tag_list[-1] = 'CB'

          # Only append if previous element is not a closing bracket
          #
          elif ((token_list == []) or (tag_list[-1] != 'CB')):
            token_list.append(tmp_val)
            tag_list.append('CB')

        else:  # An unknown element - - - - - - - - - - - - - - - - - - - - - -
          token_list.append(tmp_val)
          tag_list.append('UN')

      # Finally remove the processed elements from the original element list
      #
      org_list = org_list[tmp_len:]  # Remove processed elements

    # Remove certain elements from start and end - - - - - - - - - - - - - - -
    #
    while ((len(tag_list) > 1) and (tag_list[0] in ['CO','HY','SP'])):
      tag_list =   tag_list[1:]
      token_list = token_list[1:]
    while ((len(tag_list) > 1) and (tag_list[-1] in ['CO','HY','SP'])):
      tag_list =   tag_list[:-1]
      token_list = token_list[:-1]

    # Remove vertical bars or brackets from start and end - - - - - - - - - - -
    #
    while ((len(tag_list) > 2) and \
           ((tag_list[0] == 'VB') and (tag_list[-1] == 'VB')) or \
           ((tag_list[0] == 'OB') and (tag_list[-1] == 'CB'))):
      tag_list =   tag_list[1:-1]
      token_list = token_list[1:-1]

    return [token_list, tag_list]

  # ---------------------------------------------------------------------------

  def __write_hmm_train_record__(self, in_str, tok_list, tag_perm_list,
                                 state_list, hmm_prob):
    """Write the given record, the list of its tag sequence permutations, and
       possibly its HMM state sequence into the training file.

       If the sequence frequency file is given also insert the sequences with
       their probabilities and input strings into the sequence frequency
       dictionary.
    """

    self.hmm_train_fp.write('# Input string: %s' % (in_str)+os.linesep)
    self.hmm_train_fp.write('# Token list:   %s' % (tok_list)+os.linesep)

    list_len = len(tok_list)

    if (state_list == None):  # Create an empty state sequenct list
      state_list = [' ']*list_len

    for tag_list in tag_perm_list:

      assert len(tag_list) == len(tok_list)

      tag_state_list = []

      for i in range(list_len):
        tag_state_list.append(tag_list[i]+':'+state_list[i])

      tag_state_str = ', '.join(tag_state_list)
      self.hmm_train_fp.write('%s' % (tag_state_str)+os.linesep)

      #  Add sequence into sequence probability dictionary if defined - - - - -
      #
      if (self.hmm_seq_prob_file != None):
        seq_list = self.hmm_seq_prob_dict.get(tag_state_str, [])
        seq_list.append((in_str, hmm_prob))
        self.hmm_seq_prob_dict[tag_state_str] = seq_list

    self.hmm_train_fp.write(os.linesep)

  # ---------------------------------------------------------------------------

  def __write_hmm_seq_prob_file__(self, hmm_name_str):
    """Write the dictionary with sequence frequencies into file.

       Sort according to frequency (count of occurrence) of a sequence, with
       high frequencies first (as these are the ones one wants to get correct).
    """

    try:
      self.hmm_seq_prob_fp = open(self.hmm_seq_prob_file, 'w')
    except:
      logging.exception('Cannot write to HMM sequence probability file: ' + \
                        '"%s"' % (self.hmm_seq_prob_file))

    self.hmm_seq_prob_fp.write('#'+'#'*70 + os.linesep)  # Write a header
    self.hmm_seq_prob_fp.write('# HMM sequence probabilities for: "%s"' % \
                            (hmm_name_str)+ os.linesep)
    self.hmm_seq_prob_fp.write('#' + os.linesep)
    self.hmm_seq_prob_fp.write('# Created '+time.ctime(time.time()) + \
                            os.linesep)
    self.hmm_seq_prob_fp.write('#' + os.linesep)

    sorted_seq_key_list = []  # Make a list for easy sorting

    num_rec = 0  # Count total number of records

    for (tag_state_str, seq_list) in self.hmm_seq_prob_dict.iteritems():
      hmm_prob = self.hmm_seq_prob_dict[tag_state_str][0][1]

      sorted_seq_key_list.append((len(seq_list), hmm_prob, tag_state_str))

      num_rec += len(seq_list)

    sorted_seq_key_list.sort(reverse=True)

    cum_freq = 0  # Cumulative frequency seen so far

    for (freq, hmm_prob, tag_state_str) in sorted_seq_key_list:
      self.hmm_seq_prob_fp.write('# Pattern: %s' % (tag_state_str) + \
                                 os.linesep)
      cum_freq += freq
      freq_perc = 100.0*freq / float(num_rec)
      cum_freq_perc = 100.0*cum_freq / float(num_rec)

      self.hmm_seq_prob_fp.write('#   Frequency: %d (%.2f%%, cumulative:' % \
                                 (freq, freq_perc) + ' %.2f%%)' % \
                                 (cum_freq_perc) + os.linesep)
      examples = self.hmm_seq_prob_dict[tag_state_str]
      self.hmm_seq_prob_fp.write('#   Normalised HMM Viterbi probability: ' + \
                                 '%.20f' % (hmm_prob) + os.linesep)
      self.hmm_seq_prob_fp.write('# Example input strings:' + os.linesep)
      for (example_str, prob) in examples[:10]:  # Only 10 first examples
        self.hmm_seq_prob_fp.write('#   %s' % (example_str) + os.linesep)

      self.hmm_seq_prob_fp.write('%s' % (tag_state_str) + os.linesep)
      self.hmm_seq_prob_fp.write(os.linesep)

    self.hmm_seq_prob_fp.write('# End.' + os.linesep)
    self.hmm_seq_prob_fp.close()

  # ---------------------------------------------------------------------------

  def finalise(self):
    """Method to do any final work, such as writing information to files.
    """

    return  # Nothing to do for general case.

  # ---------------------------------------------------------------------------

  def log(self, instance_var_list = None):
    """Write a log message with the basic component standardiser instance
       variables plus the instance variable provided in the given input list
       (assumed to contain pairs of names (strings) and values).
    """

    logging.info('')
    logging.info('Component standardiser: "%s"' % (self.description))
    logging.info('  Input fields:        %s' % (str(self.in_fields)))
    logging.info('  Output fields:       %s' % (str(self.out_fields)))
    logging.info('  Field separator:     "%s"' % (self.field_sep))
    logging.info('  Check word spilling: %s' % (str(self.check_word_spill)))
    if (self.corr_list != []):
      logging.info('  Length of correction list: %d' % (len(self.corr_list)))
    if (self.tag_table != {}):
      logging.info('  Length of tag lookup table: %d' % (len(self.tag_table)))

    if (instance_var_list != None):
      logging.info('  Standardiser specific variables:')
      max_name_len = 0
      for (name, value) in instance_var_list:
        max_name_len = max(max_name_len, len(name))

      for (name, value) in instance_var_list:
        pad_spaces = (max_name_len-len(name))*' '
        logging.info('    %s %s' % (name+':'+pad_spaces, str(value)))

# =============================================================================

class DateStandardiser(ComponentStandardiser):
  """Class for standardising dates into three field: day, month, year.

     The 'output_fields' must be a list of three fields, the first field is for
     day, the second for month and the third for year. Fields can be set to
     'None' if no output is to be written, as long as at least one field is
     set. For example, if one is only interested in the year, then the output
     fields can be set to: output_fields = [None, None, "year_field"].

     Note that field spill checking is not possible for this standardiser.

     The additional arguments (besides the base class arguments) that need to
     be set when this standardiser is initialised are:

       parse_formats  A list containing one or several strings with a date
                      parsing format. The parsing formats must each contain
                      three of the following format strings with a space in
                      between (e.g. '%d %m %Y'):
                        %b  Abbreviated month name (Jan, Feb, Mar, etc.)
                        %B  Full month name (January, February, etc.)
                        %d  Day of the month as a decimal number [01,31]
                        %m  Month as a decimal number [01,12]
                        %y  Year without century as a decimal number [00,99]
                        %Y  Year with century as a decimal number
                      The parsing routine tries one format after the other and
                      as soon as a format is successful returns the
                      corresponding parsed date. Date formats more commonly
                      used should therefore be at the top of the parse format
                      list.
                      Besides the spaces between the three parsing directives,
                      format strings must not contain any other characters,
                      such as '\/:-' etc. as they will be removed from the
                      input values before date parsing is attempted.
       pivot_year     Value of pivot year (between 00 and 99) that controls
                      expansion of two-digit year values into four-digit year
                      values. Two-digits years smaller than the pivot year will
                      be expanded into 20XX, years larger and equal than the
                      pivot year will be expanded into 19xx
                      For example: pivot_year = 03:  68 -> 1968
                                                     03 -> 1903
                                                     02 -> 2002
                      The default value is the current year plus 1.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor.
    """

    self.parse_formats = None  # Initialise attributes
    self.pivot_year =    None

    # Process all keyword arguments
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('parse_for')):
        auxiliary.check_is_list('parse_formats', value)
        self.parse_formats = value

      elif (keyword.startswith('pivot_y')):
        auxiliary.check_is_integer('pivot_year', value)
        if (value < 0) or (value > 99):
          logging.exception('Argument "pivot_year" is not an integer ' + \
                            'number between 0 and 99: %s' % (str(value)))
          raise Exception
        self.pivot_year = value

      else:
        base_kwargs[keyword] = value

    ComponentStandardiser.__init__(self, base_kwargs)  # Process base arguments

    # Check if the needed attributes are set  - - - - - - - - - - - - - - - - -
    #
    if (len(self.out_fields) != 3):
      logging.exception('Attribute "output_fields" is not a list with ' + \
                        'three elements: %s' % (str(self.out_fields)))
      raise Exception

    auxiliary.check_is_not_none('parse_formats', self.parse_formats)

    # Check if parse format strings are in a valid form - - - - - - - - - - - -
    #
    for format_str in self.parse_formats:
      auxiliary.check_is_string('Format string "%s"' % (str(format_str)),
                                format_str)
      if (len(format_str) != 8):
        logging.exception('Format string has wrong length (should be 8 ' + \
                          'characters long): "%s"' % (format_str))
        raise Exception

      tmp_str = format_str

      while (tmp_str != ''):  # Check if all directives are valid
        if (tmp_str[:2] not in ['%b','%B','%d','%m','%y','%Y']):
          logging.exception('Illegal directive in format string "%s"' % \
                          (format_str))
          raise Exception

        tmp_str = tmp_str[3:]  # Remove directive and following space

    # Check if pivot year is set, if not set to current year plus one - - - - -
    #
    if (self.pivot_year == None):
      this_year = time.localtime(time.time())[0]
      next_year = (int(this_year) % 100) + 1  # Get next year as integer number

      if (next_year < 0) or (next_year > 99):
        logging.exception('Illegal value for "next_year" (not between 0 ' + \
                          'and 99): '+ str(next_year))
        raise Exception

      self.pivot_year = next_year

    # Set the correction list specifically - - - - - - - - - - - - - - - - - -
    #
    self.corr_list = [('/',' '), ('\\',' '), ('.',' '), (',',' '), ('+',' '),
                      ('~',' '), (':',' '),  (';',' '), ('-',' '), ('=',' '),
                      ('~',' '), ("'",' '),  ('"',' ')]

    self.field_sep =        ''    # Set to a specific separator string
    self.check_word_spill = False  # No word spill checking for dates

    self.log([('Pivot year', self.pivot_year),
              ('Parse formats', self.parse_formats)])  # Log a message

  # ---------------------------------------------------------------------------

  def standardise(self, in_str, clean_in_str):
    """Standardise the date defined in the input string.

       Returns a list containing year, month, day values (as strings).
    """

    if (clean_in_str.strip() == ''):  # No date given
      return ['','','']

    date_try = None

    # Try one date format after the other until success or none worked  - - - -
    #
    for format_str in self.parse_formats:
      try:
        date_try = time.strptime(clean_in_str, format_str)
      except:
        pass

      if (date_try != None):  # Parsing was successful
        break

    if (date_try != None):  # Date successfully parsed - - - - - - - - - - - -
      day =   str(date_try[2])
      month = str(date_try[1])
      year =  str(date_try[0])

      if ('%y' in format_str):  # A two-digit year was parsed, check for pivot
                                # string
        if ((year.startswith('20')) and (int(year[2:]) > self.pivot_year)):
          year = '19'+year[2:]
    else:
      logging.warn('Could not parse date string: "%s"' % (in_str))
      return ['','','']

    return [day, month, year]


# =============================================================================

class PhoneNumStandardiser(ComponentStandardiser):
  """Class for standardising telephone numbers into lists with elements
     [country_code, country_name, area_code, number, extension].

     The 'output_fields' must be a list of five fields with the above given
     elements. Fields can be set to 'None' if no output is to be written into a
     corresponding field, as long as at least one field is set.

     The additional argument (besides the base class arguments) that can be set
     when this standardiser is initialised is:

       default_country  A string which sets the default country for parsing
                        telephone numbers (if no country code is available in
                        the input phone number). Currently 'australia' or
                        'canada/usa' are possible values.

     If no 'default_country' is given the value will be set to 'australia'.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor.
    """

    self.default_country = 'australia'  # Initialise attributes

    # Process all keyword arguments
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('default_c')):
        auxiliary.check_is_string('default_country', value)

        if (value.lower() not in ['australia','canada/usa']):
          logging.exception('Illegal value for argument "default_country"' + \
                            ' (not "australia" or "canada/usa"): %s' % (value))
        self.default_country = value.lower()

      else:
        base_kwargs[keyword] = value

    ComponentStandardiser.__init__(self, base_kwargs)  # Process base arguments

    # Check if the needed attributes are set  - - - - - - - - - - - - - - - - -
    #
    if (len(self.out_fields) != 5):
      logging.exception('Attribute "output_fields" is not a list with five '+ \
                        'elements: %s' % (str(self.out_fields)))
      raise Exception

    self.log([('Default country', self.default_country)])  # Log a message

  # ---------------------------------------------------------------------------

  def standardise(self, in_str, clean_in_str):
    """Standardise the phone number in the input string.

       Returns a list containing country_code, country_name, area_code, number,
       and extension.
    """

    if (clean_in_str.strip() == ''):  # No phone number given
      return ['','','','','']

    parsed_phone_number = phonenum.str_to_phonenum(clean_in_str,
                                         default_country=self.default_country)

    # Check for an empty parsed phone number
    #
    if (parsed_phone_number == []):
      logging.warn('Could not parse phone number string "%s"' % (in_str))
      return ['','','','','']

    return parsed_phone_number


# =============================================================================

class NameStandardiser(ComponentStandardiser):
  """Class for standardising names into lists with elements: [title, gender
     guess, given names, alternative given names, surnames, alternative
     surnames].

     The 'output_fields' must be a list of six field names with the above given
     elements. Fields can be set to 'None' if no output is to be written into a
     corresponding field, as long as at least one field is set.

     Simple names (made of one or two words) will be handled using a rule based
     approach, while more complex names will be standardised using a hidden
     Markov model (HMM) approach.

     The additional argument (besides the base class arguments) that can be set
     when this standardiser is initialised is:

       female_titles      A list of female titles (like 'ms'), used to guess a
                          gender.
       male_titles        A list of male titles (like 'mr'), used to guess a
                          gender.
       first_name_comp    The expected first component in the names, either
                          'gname' (for given names) or 'sname' (for surnames).
                          Default value is 'gname'.
       name_hmm           This can either be set to None (default) or to a HMM,
                          which will be used to parse the cleaned and tagged
                          input names. If set to None then the input names will
                          not be written into the output data set, but only
                          into the 'hmm_train' file (which has to be set to a
                          string).
       hmm_train_file     This can either be set to None (default) or to a
                          string which is assumed to be a file name. If a file
                          name is given, then the cleaned and tagged input
                          names will be written into a HMM training file for
                          later training of a corresponding HMM. If set to None
                          no such training file will be written.
       hmm_seq_prob_file  This can either be set to None (default) or to a
                          string which is assumed to be a file name. It can
                          only be set to a file name if the name HMM is given.
                          If a file name is provided, then the training
                          sequences will be written into file grouped and
                          sorted according to their HMM probabilities, with the
                          sequences with low probabilities coming first. This
                          allows for easy identification and correction of
                          names that were hard to standardise by the HMM.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor.
    """

    self.female_titles =     []
    self.male_titles =       []
    self.first_name_comp =   'gname'
    self.name_hmm =          None
    self.hmm_train_file =    None
    self.hmm_seq_prob_file = None
    self.hmm_seq_prob_dict = {}

    # Process all keyword arguments
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('male_t')):
        auxiliary.check_is_list('male_titles', value)
        i = 0
        for title in value:
          auxiliary.check_is_string('male_title[%d]' % (i), value[i])
          i += 1
        self.male_titles = value

      elif (keyword.startswith('female_t')):
        auxiliary.check_is_list('female_titles', value)
        i = 0
        for title in value:
          auxiliary.check_is_string('female_title[%d]' % (i), value[i])
          i += 1
        self.female_titles = value

      elif (keyword.startswith('first_na')):
        if (value not in ['gname','sname']):
          logging.exception('Argument "first_name_comp" must be either ' + \
                            '"gname" or "sname"')
          raise Exception
        self.first_name_comp = value

      elif (keyword.startswith('name_hm')):
        self.name_hmm = value

      elif (keyword.startswith('hmm_tr')):
        if (value != None):
          auxiliary.check_is_string('hmm_train_file', value)
        self.hmm_train_file = value

      elif (keyword.startswith('hmm_seq')):
        if (value != None):
          auxiliary.check_is_string('hmm_seq_prob_file', value)
        self.hmm_seq_prob_file = value

      else:
        base_kwargs[keyword] = value

    ComponentStandardiser.__init__(self, base_kwargs)  # Process base arguments

    # Initialise a dictionary with counters for how many names were
    # standardised into which name structure (G=given name, S=surname)
    #
    self.count_dict = {'empty':0, 'G':0, 'S':0, 'GS':0, 'SG':0, 'HMM':0}

    # Check if the needed attributes are set  - - - - - - - - - - - - - - - - -
    #
    if (len(self.out_fields) != 6):
      logging.exception('Attribute "output_fields" is not a list with six '+ \
                        'elements: %s' % (str(self.out_fields)))
      raise Exception

    # Check the male and female titles have no common values
    #
    if (set(self.female_titles).intersection(set(self.male_titles)) != set()):
      logging.exception('At least one title word is in both the male and ' + \
                        'female title list')
      raise Exception

    # Check HMM related attributes - - - - - - - - - - - - - - - - - - - - - -
    #
    if ((self.name_hmm == None) and (self.hmm_train_file == None)):
      logging.exception('Both name HMM and HMM training file set to None.' + \
                        ' At least one must be not None.')
      raise Exception

    if (self.hmm_train_file != None):  # Open HMM training file for writing
      try:
        self.hmm_train_fp = open(self.hmm_train_file, 'w')
      except:
        logging.exception('Cannot write to HMM training file: "%s"' % \
                          (self.hmm_train_file))

      self.hmm_train_fp.write('#'+'#'*70 + os.linesep)  # Write a header
      self.hmm_train_fp.write('# Tagged name data for HMM training.' + \
                         os.linesep)
      self.hmm_train_fp.write('#' + os.linesep)
      self.hmm_train_fp.write('# Created '+time.ctime(time.time()) + \
                              os.linesep)
      self.hmm_train_fp.write('#' + os.linesep)
      self.hmm_train_fp.write('# Possible HMM state for names are:' + \
                              os.linesep)
      self.hmm_train_fp.write('# baby, knwn, andor, gname1, gname2, ghyph,' + \
                              ' gopbr, gclbr, agname1,' + os.linesep)
      self.hmm_train_fp.write('# agname2, coma, sname1, sname2, shyph, ' + \
                              'sopbr, sclbr, asname1, asname2,' + os.linesep)
      self.hmm_train_fp.write('# pref1, pref2, rubb' + os.linesep)
      self.hmm_train_fp.write(os.linesep)

    if ((self.hmm_seq_prob_file != None) and (self.name_hmm == None)):
      logging.exception('HMM sequence probability file given but no HMM' + \
                        ' provided.')
      raise Exception

    if (self.name_hmm != None):
      name_hmm_name = self.name_hmm.description
    else:
      name_hmm_name = 'None'

    self.log([('Female title list',             self.female_titles),
              ('Male title list',               self.male_titles),
              ('First name component',          self.first_name_comp),
              ('Name HMM',                      name_hmm_name),
              ('HMM training file',             str(self.hmm_train_file)),
              ('HMM sequence probability file', str(self.hmm_seq_prob_file))])

  # ---------------------------------------------------------------------------

  def standardise(self, in_str, clean_in_str):
    """Standardise the name defined in the input string.

       Returns a list containing (as strings): title, gender guess, given
       names, alternative given names, surnames, alternative surnames.
    """

    if (clean_in_str == ''):  # No name given
      self.count_dict['empty'] = self.count_dict['empty'] + 1
      return ['','','','','','']

    # Tag the name string and split into elements - - - - - - - - - - - - - - -
    #
    [tok_list, tag_list] = self.__tag_component__(clean_in_str)

    # Get title words - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    [title_list, tok_list, tag_list] = self.__get_title__(tok_list, tag_list)
    title_str = ' '.join(title_list)

    # Guess gender - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    gender_guess = self.__get_gender_guess__(title_str, tag_list)

    # Check if the word 'nee' (with tag 'NE') is in the list, and if so decide
    # if it is a separator or a name word.
    #
    # Also check for pairs of vertical bars or opening/closing brackets with
    # nothing in between.
    #
    i = 0
    while (i < len(tag_list)):

      if (tag_list[i] == 'NE'):

        # 'nee' is not a separator if it:
        # - is at the beginning or at the end of the sequence
        # - if the word before 'nee' is either a name prefix or a 'baby of'
        #   element
        # - if the word before is a known saint name (but not the word 'saint'
        #   by itself
        # For example:
        # 'nee miller' -> 'nee' is given name
        # 'peter de nee' -> 'de nee' is surname
        # 'paul saint nee' -> 'saint nee' is surname
        # 'saint paul nee' -> 'nee' is surname
        # 'paula miller (nee jones)' -> 'nee' is separator ('jones' is the
        #                               alternative surname
        # 'tim jones nee miller' -> 'nee' is separator ('miller' is the
        #                           alternative surname
        # 'peter, son of nee miller' -> 'nee' is 'surname' of 'peter's' dad
        #
        if (((i == 0) or (i == len(tag_list)-1)) or \
            (tag_list[i-1] in ['PR','BO']) or \
            ((tag_list[i-1] == 'ST') and ('_' in tok_list[i-1]))):
          tag_list[i] = 'UN'  # Make it an unknown word tag
        else:
          tag_list[i] = 'SP'  # Make it a separator tag

      elif (tag_list[i] == 'VB'):  # Pair of vertical bars
        if ((i < len(tag_list)-1) and (tag_list[i+1] == 'VB')):
          tag_list = tag_list[:i]+tag_list[i+2:] # Remove pair of vertical bars
          tok_list = tok_list[:i]+tok_list[i+2:]

      elif (tag_list[i] == 'OB'):  # Pair of opening and closing brackets
        if ((i < len(tag_list)-1) and (tag_list[i+1] == 'CB')):
          tag_list = tag_list[:i]+tag_list[i+2:]  # Remove pair of brackets
          tok_list = tok_list[:i]+tok_list[i+2:]

      i += 1

    # First handle the various common cases specifically to speed up execution
    #
    if (len(tag_list) == 1):  # Only one element - - - - - - - - - - - - - - -

      if (tag_list[0] in ['GF','GM']):  # Known given name
        self.count_dict['G'] = self.count_dict['G'] + 1
        return [title_str, gender_guess, tok_list[0], '', '', '']
      else:  # Otherwise output as surname
        self.count_dict['S'] = self.count_dict['S'] + 1
        return [title_str, gender_guess, '', '', tok_list[0], '']

    elif (len(tag_list) == 2):  # Two elements, assume given and surname - - -

      if (self.first_name_comp == 'gname'):

        # Special case: First element is known surname and second element is
        # known given name (so assume they are swapped)
        #
        if ((tag_list[0] == 'SN') and (tag_list[1] in ['GM','GF'])):
          self.count_dict['SG'] = self.count_dict['SG'] + 1
          return [title_str, gender_guess, tok_list[1], '', tok_list[0], '']

        else:  # Assume first element is given name and second is surname
          self.count_dict['GS'] = self.count_dict['GS'] + 1
          return [title_str, gender_guess, tok_list[0], '', tok_list[1], '']

      else:  # Same for surname (assumed to be first name element)

        if ((tag_list[0] in ['GM','GF']) and (tag_list[1] == 'SN')):
          self.count_dict['GS'] = self.count_dict['GS'] + 1
          return [title_str, gender_guess, tok_list[0], '', tok_list[1], '']

        else: # Assume first element is surname and second is given name
          self.count_dict['SG'] = self.count_dict['SG'] + 1
          return [title_str, gender_guess, tok_list[1], '', tok_list[0], '']

    # Any other case that contains at least three elements - - - - - - - - - -

    # Create all permutations of the tag list
    #
    tag_perm_list = mymath.perm_tag_sequence(tag_list)

    if (self.name_hmm != None):  # Parse using the HMM
      (name_list, state_seq, hmm_prob) = self.__get_name_hmm__(tok_list,
                                                               tag_perm_list)
      self.count_dict['HMM'] = self.count_dict['HMM'] + 1

    else:
      name_list = ['','','','']
      state_seq =  None  # No HMM state sequence available
      hmm_prob = None

    if (self.hmm_train_file != None):
      self.__write_hmm_train_record__(in_str, tok_list, tag_perm_list,
                                      state_seq, hmm_prob)

    given_name_str =     ' '.join(name_list[0])
    alt_given_name_str = ' '.join(name_list[1])
    surname_str =        ' '.join(name_list[2])
    alt_surname_str =    ' '.join(name_list[3])

    return [title_str, gender_guess, given_name_str, alt_given_name_str,
            surname_str, alt_surname_str]

  # ---------------------------------------------------------------------------

  def __get_gender_guess__(self, title_str, tag_list):
    """Try to guess a gender from the given title string and given name lists.

       First, the title string is checked, and if no gender can be guessed from
       this the list of tags is searched for given name tags.

       The returned string is either 'female', 'male' or '' (if no gender has
       been found).
    """

    # Check title string - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (title_str in self.female_titles):
      return 'female'

    elif (title_str in self.male_titles):
      return 'male'

    male_count   = 0
    female_count = 0

    # Check given name tags only if no gender in title words found - - - - - -
    #
    for tag in tag_list:
      if ('GM' in tag):
        male_count += 1
      if ('GF' in tag):
        female_count += 1

    # Check if only male or only female gender found - - - - - - - - - - - - -
    #
    if (male_count > 0) and (female_count == 0):
      gender_guess = 'male'
    elif (male_count == 0) and (female_count > 0):
      gender_guess = 'female'
    else:
      gender_guess = ''

    return gender_guess

  # ---------------------------------------------------------------------------

  def __get_title__(self, tok_list, tag_list):
    """Extract the title component of a name.

       This method extracts all title tokens (and removes them from both the
       given token and tag lists).

       It returns three lists: the extracted title words, and the modified
       token and tag lists with all title words removed.
    """

    assert len(tok_list) == len(tag_list)

    title_list =   []
    mod_tok_list = []
    mod_tag_list = []

    list_len = len(tag_list)

    # Extract all title words (list elements with a TI tag) - - - - - - - - - -
    #
    for i in range(list_len):
      if (tag_list[i] == 'TI'):
        if (tok_list[i] not in title_list):  # Don't insert duplicates
          title_list.append(tok_list[i])
      else:
        mod_tag_list.append(tag_list[i])
        mod_tok_list.append(tok_list[i])

    return [title_list, mod_tok_list, mod_tag_list]

  # ---------------------------------------------------------------------------

  def __get_name_hmm__(self, tok_list, tag_perm_list):
    """This method parses the input token sequence using the tag sequence
       permutations based on a hidden Markov model and extracts given- and
       surnames into four lists (given names, alternative given names,
       surnames, alternative surnames) using the Viterbi alogrithm which for
       each tag sequence returns the best state sequence through the HMM. The
       overall best state sequence will be selected.

       Various post-processing steps are done after the HMM sequence is
       assigned.

       Returns a list with the four name sub-lists (containing the parsed given
       names, alternative given names, surnames, and alternative surnames), as
       well as the best HMM observation sequence and corresponding HMM
       probability.
    """

    # Give all tag sequences to the HMM and keep - - - - - - - - - - - - - - -
    # the one with highest probability
    #
    max_prob =       -1.0
    best_state_seq = []
    best_tag_list =  []

    for tag_list in tag_perm_list:

      assert len(tok_list) == len(tag_list)

      [state_seq, prob] = self.name_hmm.viterbi(tag_list)
      if (prob > max_prob):
        best_state_seq = state_seq
        best_tag_list =  tag_list
        max_prob =       prob

      logging.info('  Sequence: %s has Viterbi probability: %f' % \
                   (str(tag_list), prob))

    logging.info('Best state sequence: %s with tag sequence %s' % \
                 (str(best_state_seq), str(best_tag_list)))

    if (max_prob == 0.0):
      logging.warn('Probability is 0.0 for best state sequence: %s' % \
                   (str(best_state_seq)))

    tag_list = best_tag_list  # Shorter

    if (len(tag_list) != len(tok_list)):
      logging.exception('Length of token list and best tag list differ: ' + \
                        '%s, %s' % (str(tok_list), str(tag_list)))
      raise Exception

    if (tag_list == []):
      logging.warn('Empty tag list returned from HMM')
      return [[[],[],[],[]],[],0.0]  # Return empty name list

    tag_list_len = len(tag_list)

    # Normalise the HMM probability by the length of the tag sequence
    #
    norm_max_prob = max_prob / float(tag_list_len)

    # Post-process the state and token sequences - - - - - - - - - - - - - - -
    #
    name_list = [[],[],[],[],[]]  # Resulting output list

    i = 0

    while (i < tag_list_len):  # Loop over tags, tokens and HMM states - - - -
      token = tok_list[i]
      tag =   tag_list[i]
      state = best_state_seq[i]

      if (token in ['|','(',')']):  # Do not output vertical bars and brackets
        pass

      elif (state in ['knwn','andor']):  # Skip over' known as' or 'and/or' - -
        pass  # These should only be separators between names and altern. names

      elif (state == 'baby'):  # 'Baby of' sequence - - - - - - - - - - - - - -

        if (i < tag_list_len-1):
          if (best_state_seq[i+1] in ['gname1','gname2']):  # Add to givenname
            name_list[1].append(token)

          elif (best_state_seq[i+1] in ['agname1','agname2']): # Add to agname
            name_list[2].append(token)

          elif (best_state_seq[i+1] in ['sname1','sname2']):  # Add to surname
            name_list[3].append(token)

          elif (best_state_seq[i+1] in ['asname1','asname2']): # Add to asname
            name_list[4].append(token)

          else:
            logging.warn('Strange situation with "baby of": %s' % \
                         (str(tok_list)))
        else:
          logging.warn('Strange situation with "baby of": %s' % \
                       (str(tok_list)))

      elif (state in ['gname1','gname2']):  # Givennames  - - - - - - - - - - -
        name_list[0].append(token)  # Append to given name list

      elif (state in ['agname1','agname2']):  # Alternative givennames  - - - -
        name_list[1].append(token)  # Append to alternative given name list

      elif (state == 'ghyph'):  # Givenname hyphen  - - - - - - - - - - - - - -

        if (i > 0) and (best_state_seq[i-1] in ['gname1','gname2']):  # G-name
          if (name_list[1] != []):  # Only append hyphen after given name
            name_list[0].append(token)  # Append to given name list
          else:
            logging.warn('Strange hyphen situation in given name: %s' % \
                         (str(tok_list)))

        elif (i > 0) and (best_state_seq[i-1] in ['agname1','agname2']): # agn
          if (name_list[1] != []):  # Only append hyphen after alt. given name
            name_list[1].append(token)  # Append to alternative given name list
          else:
            logging.warn('Strange hyphen situation in alternative given ' + \
                         'name: %s' % (str(tok_list)))
        else:
          logging.warn('Strange hyphen situation in alternative given ' + \
                       'name: %s' % (str(tok_list)))

      elif (state in ['sname1','sname2']):  # Surnames  - - - - - - - - - - - -
        name_list[2].append(token)  # Append to surname list

      elif (state in ['asname1','asname2']):  # Alternative surnames  - - - - -
        name_list[3].append(token)  # Append to alternative surname list

      elif (state == 'shyph'):  # Surname hyphen  - - - - - - - - - - - - - - -
        if (i > 0) and (best_state_seq[i-1] in ['sname1','sname2']):  # Surname
          if (name_list[2] != []):  # Only append hyphen after a surname
            name_list[2].append(token)  # Append to surname list
          else:
            logging.warn('Strange hyphen situation in surname: %s' % \
                         (str(tok_list)))

        elif (i > 0) and (best_state_seq[i-1] in ['asname1','asname2']): # asn
          if (name_list[3] != []):  # Only append hyphen after a alt. surname
            name_list[3].append(token)  # Append to alternative surname list
          else:
            logging.warn('Strange hyphen situation in alternative ' + \
                         'surname: %s' % (str(tok_list)))
        else:
          logging.warn('Strange hyphen situation in surname: %s' % \
                       (str(tok_list)))

      elif (state in ['pref1','pref2']):  # Name prefix - - - - - - - - - - - -
        if (i < tag_list_len-1) and (best_state_seq[i+1] in ['pref1','pref2']):
                                          # Followed by another prefix
          token = token+' '+tok_list[i+1]  # Concatenate
          i = i+1
          state = best_state_seq[i]

        if ((i < tag_list_len-1) and \
            (best_state_seq[i+1] in ['gname1','gname2'])):
          if (name_list[0] != []):  # There is already a given name
            name_list[0].append(token)  # Append to given name list
            logging.warn('Strange name prefix situation in given name: %s' % \
                         (str(tok_list)))
          else:
            name_list[0].append(token)  # Append to given name list

        elif ((i < tag_list_len-1) and \
              (best_state_seq[i+1] in ['agname1','agname2'])):
          if (name_list[1] != []): # There is already an alternative given name
            name_list[1].append(token)  # Append to alternative given name list
            logging.warn('Strange name prefix situation in alternative ' + \
                         'given name: %s' % (str(tok_list)))
          else:
            name_list[1].append(token)  # Append to alternative given name list

        elif ((i < tag_list_len-1) and \
              (best_state_seq[i+1] in ['sname1','sname2'])):
          if (name_list[2] != []):  # There is already a surname
            name_list[2].append(token)  # Append to surname list
            logging.warn('Strange name prefix situation in surname: %s' % \
                         (str(tok_list)))
          else:
            name_list[2].append(token)  # Append to surname list

        elif ((i < tag_list_len-1) and \
              (best_state_seq[i+1] in ['asname1','asname2'])):
          if (name_list[3] != []):  # There is already an alternative  surname
            name_list[3].append(token)  # Append to alternative surname list
            logging.warn('Strange name prefix situation in alternative ' + \
                         'surname: %s' % (str(tok_list)))
          else:
            name_list[3].append(token)  # Append to alternative surname list
        else:
          logging.warn('Strange name prefix situation: %s' % (str(tok_list)))

      else:  # Should never happen
        logging.warn('This should never happen! State: %s, tag: %s, token' % \
                     (state, tag) + ': %s, token list: %s' % (token, tok_list))
      i +=1

    # Finally do some tests on the output fields  - - - - - - - - - - - - - - -

    # Check if a name sub-list has more than three elements, if so print out
    #
    for n_list in name_list:
      if (len(n_list) > 3):
        logging.warn('A name output field contains more than three ' + \
                     'elements: %s' % (str(n_list)))

    # Check if a name component is not allocated but its alternative name is
    #
    if (name_list[0] == []) and (name_list[1] != []):
      logging.warn('No given name but an alternative given name: %s' % \
                    (str(name_list[1])) + ' -> Corrected')
      name_list[0] = name_list[1][:]  # Move alt. given name to given name
      name_list[1] = []

    if (name_list[2] == []) and (name_list[3] != []):
      logging.warn('No surname but an alternative surname: %s' % \
                   (str(name_list[3])) + ' -> Corrected')
      name_list[2] = name_list[3][:]  # Move alternative surname to given name
      name_list[3] = []

    # Check if a givenname is given but no surname  - - - - - - - - - - - - - -
    #
    if (name_list[0] != []) and (name_list[2] == []):

      if (name_list[1] != []):  # An alternative given name is availabe
        logging.warn('No surname but an alternative given name: %s' % \
                     (str(name_list)) + ' -> Corrected')
        name_list[2] = name_list[1][:]  # Move alt. given name to surname
        name_list[1] = []

      elif (len(name_list[1]) > 1):  # More than one given name available
        logging.debug('No surname but more than one given name: %s' % \
                      (str(name_list)) + ' -> Corrected')
        name_list[2] = name_list[0][-1][:]  # Move last given name to surname
        name_list[0] = name_list[0][:-1][:]

    # Log current status of sub-lists extracted - - - - - - - - - - - - - - - -
    #
    logging.info('  Best HMM state sequence:   %s' % (str(best_state_seq)))
    logging.info('    Normalised HMM probability:   %f' % (norm_max_prob))
    logging.info('  Parsed name elements:')
    logging.info('    Given names:             %s' % (str(name_list[0])))
    logging.info('    Alternative given names: %s' % (str(name_list[1])))
    logging.info('    Surnames:                %s' % (str(name_list[2])))
    logging.info('    Alternative surnames:    %s' % (str(name_list[3])))

    return [name_list, best_state_seq, norm_max_prob]

  # ---------------------------------------------------------------------------

  def finalise(self):
    """If the HMM sequence probability file is given write it out.
    """

    std_type_dict = {'empty':'Empty address', 'G':'Given name only',
                     'S':'Surname only', 'GS':'Given name / surname',
                     'SG':'Surname / given name', 'HMM':'HMM standardisation'}

    logging.info('Standardisation statistics:')
    for t in self.count_dict:

      logging.info('  %20s: %d' % (std_type_dict[t], self.count_dict[t]))
      print '  %20s: %d' % (std_type_dict[t], self.count_dict[t])  ###########

    if (self.hmm_seq_prob_file != None):
      self.__write_hmm_seq_prob_file__('Name HMM')

# =============================================================================

class AddressStandardiser(ComponentStandardiser):
  """Class for standardising addresses into lists using a hidden Markov model
     (HMM) based approach.

     The 'output_fields' must be a list of field names from the output data set
     into which the standardised names will be written.
       building_name
       post_address_type
       post_address_number
       lot_number_prefix
       lot_number
       lot_number_suffix
       flat_number_prefix
       flat_number
       flat_number_suffix
       flat_type
       level_number_prefix
       level_number
       level_number_suffix
       level_type
       number_first_prefix
       number_first
       number_first_suffix
       number_last_prefix
       number_last
       number_last_suffix
       street_name
       street_suffix
       street_type
       locality_name
       postcode
       state_abbrev
       country

     The fields correspond to the state names given in the addres HMM supplied.

     Fields can be set to 'None' if no output is to be written into a
     corresponding field, as long as at least one field is set.

     The additional argument (besides the base class arguments) that can be set
     when this standardiser is initialised is:

     The additional arguments (besides the base class arguments) are

       address_hmm        This can either be set to None (default) or to a HMM,
                          which will be used to parse the cleaned and tagged
                          input addresses. If set to None then the input
                          addresses will not be written into the output data
                          set, but only into the 'hmm_train' file (which has to
                          be set to a string).
       hmm_train_file     This can either be set to None (default) or to a
                          string which is assumed to be a file name. If a file
                          name is given, then the cleaned and tagged input
                          addresses will be written into a HMM training file
                          for later training of a corresponding HMM. If set to
                          None no such training file will be written.
       hmm_seq_prob_file  This can either be set to None (default) or to a
                          string which is assumed to be a file name. It can
                          only be set to a file name if the name HMM is given.
                          If a file name is provided, then the training
                          sequences will be written into file grouped and
                          sorted according to their HMM probabilities, with the
                          sequences with low probabilities coming first. This
                          allows for easy identification and correction of
                          names that were hard to standardise by the HMM.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor.
    """

    self.address_hmm =       None
    self.hmm_train_file =    None
    self.hmm_seq_prob_file = None
    self.hmm_seq_prob_dict = {}

    # Process all keyword arguments
    #
    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('address_hm')):
        self.address_hmm = value

      elif (keyword.startswith('hmm_tr')):
        if (value != None):
          auxiliary.check_is_string('hmm_train_file', value)
        self.hmm_train_file = value

      elif (keyword.startswith('hmm_seq')):
        if (value != None):
          auxiliary.check_is_string('hmm_seq_prob_file', value)
        self.hmm_seq_prob_file = value

      else:
        base_kwargs[keyword] = value

    ComponentStandardiser.__init__(self, base_kwargs)  # Process base arguments

    # Check if the needed attributes are set  - - - - - - - - - - - - - - - - -
    #
    if (len(self.out_fields) != 27):
      logging.exception('Attribute "output_fields" is not a list with 27 '+ \
                        ' elements: %s' % (str(self.out_fields)))
      raise Exception

    # Check HMM related attributes - - - - - - - - - - - - - - - - - - - - - -
    #
    if ((self.address_hmm == None) and (self.hmm_train_file == None)):
      logging.exception('Both address HMM and HMM training file set to ' + \
                        'None. At least one must be not None.')
      raise Exception

    if (self.hmm_train_file != None):  # Open HMM training file for writing
      try:
        self.hmm_train_fp = open(self.hmm_train_file, 'w')
      except:
        logging.exception('Cannot write to HMM training file: "%s"' % \
                          (self.hmm_train_file))

      self.hmm_train_fp.write('#'+'#'*70 + os.linesep)  # Write a header
      self.hmm_train_fp.write('# Tagged address data for HMM training.' + \
                         os.linesep)
      self.hmm_train_fp.write('#' + os.linesep)
      self.hmm_train_fp.write('# Created '+time.ctime(time.time()) + \
                              os.linesep)
      self.hmm_train_fp.write('#' + os.linesep)
      self.hmm_train_fp.write('# Possible HMM state for addresses are:' + \
                              os.linesep)
      for hmm_state_line in [
        'building_name1, building_name2, building_name3, country,',
        'post_address_type, post_address_number, lot_number_prefix,',
        'lot_number, lot_number_suffix, flat_number_prefix, flat_number,',
        'flat_number_suffix, flat_type, level_number_prefix, level_number,',
        'level_number_suffix, level_type, number_first_prefix, number_first,',
        'number_first_suffix, number_last_prefix, number_last,',
        'number_last_suffix, street_name1, street_name2, street_name3,',
        'street_suffix, street_type, locality_name1, locality_name2,',
        'locality_name3, postcode, state_abbrev', 'slash', 'hyph', 'rubb']:
        self.hmm_train_fp.write('# %s' % (hmm_state_line) + os.linesep)
      self.hmm_train_fp.write(os.linesep)

    if ((self.hmm_seq_prob_file != None) and (self.address_hmm == None)):
      logging.exception('HMM sequence probability file given but no HMM' + \
                        ' provided.')
      raise Exception

    if (self.address_hmm != None):
      address_hmm_name = self.address_hmm.description
    else:
      address_hmm_name = 'None'

    self.log([('Address HMM',                   address_hmm_name),
              ('HMM training file',             str(self.hmm_train_file)),
              ('HMM sequence probability file', str(self.hmm_seq_prob_file))])

  # ---------------------------------------------------------------------------

  def standardise(self, in_str, clean_in_str):
    """Standardise the address defined in the input string.

       Returns a list containing (as strings) the 27 address elements listed at
       the beginning of this class.
    """

    if (clean_in_str == ''):  # No address given
      return 27*['']

    # Tag the address string and split into elements - - - - - - - - - - - - -
    #
    [tok_list, tag_list] = self.__tag_component__(clean_in_str)

    # Perform some post-processing on tag and token lists - - - - - - - - - - -
    #
    mod_tag_list = []
    mod_tok_list = []

    for i in range(len(tag_list)):
      if (tag_list[i] == 'AN'):  # An alpha-numeric token, try to take apart

        if (tok_list[i][0].isdigit() == True):  # Starts with a number
          last_num = 0
          while (tok_list[i][last_num].isdigit() == True):
            last_num += 1  # Find last number starting from beginning
          first_letter = len(tok_list[i])-1
          while (tok_list[i][first_letter].isalpha() == True):
            first_letter -= 1  # Find first letter starting from end

          if (last_num-1 == first_letter):  # Same position, no mix
            mod_tok_list.append(tok_list[i][:last_num])
            if (len(mod_tok_list[-1]) == 4):
              if (((mod_tok_list[-1],) in self.tag_table) and \
                  (self.tag_table[(mod_tok_list[-1],)][1] == 'PC')):
                mod_tag_list.append('PC')
              else:
                mod_tag_list.append('N4')
            else:
               mod_tag_list.append('NU')
            mod_tok_list.append(tok_list[i][first_letter+1:])
            if (len(mod_tok_list[-1]) == 1):
              mod_tag_list.append('II')
            else:
               mod_tag_list.append('UN')
          else:  # Mixed letters and digits, append original AN tag and token
            mod_tag_list.append(tag_list[i])
            mod_tok_list.append(tok_list[i])

        else:  # Starts with a letter
          last_letter = 0
          while (tok_list[i][last_letter].isalpha() == True):
            last_letter += 1  # Find last letter starting from beginning
          first_num = len(tok_list[i])-1
          while (tok_list[i][first_num].isdigit() == True):
            first_num -= 1  # Find first number starting from end

          if (last_letter-1 == first_num):  # Same position, no mix
            mod_tok_list.append(tok_list[i][:last_letter])
            if (len(mod_tok_list[-1]) == 1):
              mod_tag_list.append('II')
            else:
               mod_tag_list.append('UN')
            mod_tok_list.append(tok_list[i][first_num+1:])
            if (len(mod_tok_list[-1]) == 4):
              if (((mod_tok_list[-1],) in self.tag_table) and \
                  (self.tag_table[(mod_tok_list[-1],)][1] == 'PC')):
                mod_tag_list.append('PC')
              else:
                mod_tag_list.append('N4')

            else:
               mod_tag_list.append('NU')
          else:  # Mixed letters and digits, append original AN tag and token
            mod_tag_list.append(tag_list[i])
            mod_tok_list.append(tok_list[i])


      else:
        mod_tag_list.append(tag_list[i])
        mod_tok_list.append(tok_list[i])

    tag_list = mod_tag_list
    tok_list = mod_tok_list

    # Create all permutations of the tag list - - - - - - - - - - - - - - - - -
    #
    tag_perm_list = mymath.perm_tag_sequence(tag_list)

    if (self.address_hmm != None):  # Parse using the HMM
      (address_list, state_seq, hmm_prob) = self.__get_address_hmm__(tok_list,
                                                                 tag_perm_list)
    else:
      address_list = 27*['']
      state_seq =   None  # No HMM state sequence available
      hmm_prob = None

    if (self.hmm_train_file != None):
      self.__write_hmm_train_record__(in_str, tok_list, tag_perm_list,
                                      state_seq, hmm_prob)

    address_result_list = []  # Put address result list together

    for i in range(27):
      address_result_list.append(' '.join(address_list[i]))

    return address_result_list

  # ---------------------------------------------------------------------------

  def __get_address_hmm__(self, tok_list, tag_perm_list):
    """This method parses the input token sequence using the tag sequence
       permutations based on a hidden Markov model and extracts upto 27 address
       elements using the Viterbi alogrithm which for each tag sequence returns
       the best state sequence through the HMM. The overall best state sequence
       will be selected.

       Various post-processing steps are done after the HMM sequence is
       assigned.

       Returns a list with three sub-lists, the first containing the parsed
       address elements in 27 sub-lists, and the second and third contain the
       the best HMM observation sequence and corresponding HMM probability,
       respectively.
    """

    # Give all tag sequences to the HMM and keep - - - - - - - - - - - - - - -
    # the one with highest probability
    #
    max_prob =       -1.0
    best_state_seq = []
    best_tag_list  = []

    for tag_list in tag_perm_list:

      assert len(tok_list) == len(tag_list)

      [state_seq, prob] = self.address_hmm.viterbi(tag_list)
      if (prob > max_prob):
        best_state_seq = state_seq
        best_tag_list =  tag_list
        max_prob =       prob

      logging.info('  Sequence: %s has Viterbi probability: %f' % \
                   (str(tag_list), prob))

    logging.info('Best state sequence: %s with tag sequence %s' % \
                 (str(best_state_seq), str(best_tag_list)))

    if (max_prob == 0.0):
      logging.warn('Probability is 0.0 for best state sequence: %s' % \
                   (str(best_state_seq)))

    tag_list = best_tag_list  # Shorter

    if (len(tag_list) != len(tok_list)):
      logging.exception('Length of token list and best tag list differ: ' + \
                        '%s, %s' % (str(tok_list), str(tag_list)))
      raise Exception

    if (tag_list == []):
      logging.warn('Empty tag list returned from HMM')
      return [27*[[]], [], 0.0]  # Return empty address list

    tag_list_len = len(tag_list)

    # Normalise the HMM probability by the length of the tag sequence
    #
    norm_max_prob = max_prob / float(tag_list_len)

    # Post-process the state and token sequences - - - - - - - - - - - - - - -

    # Build a dictionary with HMM states as keys and tokens as values
    #
    address_dict = {}

    for i in range(tag_list_len):  # Loop over tokens and states
      token = tok_list[i]
      state = best_state_seq[i]

      # Do not output commas, brackets, vertical bars and hyphens  - - - - - -
      #
      if (token in [',', '(', ')', '|', '-', '/']):
        logging.info('Discard character "%s" from input (with tag: %s)' % \
                     (token, state))
      else:
        val_list = address_dict.get(state, [])
        val_list.append(token)
        address_dict[state] = val_list

    # Merge name states (building, street, locality) - - - - - - - - - - - - -
    #
    for name_state in ['building_name','street_name','locality_name']:
      if (name_state+'1' in address_dict) or (name_state+'2' in address_dict) \
         or (name_state+'3' in address_dict):
        val_list = address_dict.get(name_state+'1', [])
        val_list += address_dict.get(name_state+'2',[])
        val_list += address_dict.get(name_state+'3',[])

        address_dict[name_state] = val_list
        for i in ['1','2','3']:
          if (name_state+i in address_dict):
            del address_dict[name_state+i]

    # Check if concatenated locality, state or country words are in - - - - - -
    # lookup-table
    #
    for name_state in ['locality_name','state_abbrev','country']:
      if (name_state in address_dict):
        name_state_list = address_dict[name_state]
        if (len(name_state_list) > 1):  # Contains more than one word
          name_state_tuple = tuple(name_state_list)  # Make it a tuple
          if (self.tag_table.has_key(name_state_tuple)):
            new_name_state = self.tag_table[name_state_tuple][0]
            address_dict[name_state] = [new_name_state]

    # Correct special cases 'st' and 'mt' at beginning of locality names - - -
    # (if there is no street type - as 'mt' (mount) and 'st' (street/saint)
    # might be wrong due to greedy tagging approach)
    #
    if (('street_type' not in address_dict) and \
        ('locality_name' in address_dict)):
      org_loc_name_list = address_dict['locality_name']
      first_loc_name = org_loc_name_list[0]

      if ((first_loc_name.startswith('saint_')) or \
          (first_loc_name.startswith('mount_'))):
        loc_name_tuple = tuple(first_loc_name[6:].replace('_',' ').split())

        if (self.tag_table.has_key(loc_name_tuple)):
          lookup_val = self.tag_table[loc_name_tuple]  # Found in look-up table
          loc_name_list = [lookup_val[0]]+org_loc_name_list[1:]

          address_dict['locality_name'] = loc_name_list
          if (first_loc_name[:5] == 'saint'):
            address_dict['street_type'] = ['street']
          else:
            address_dict['street_type'] = ['mount']
          logging.info('Corrected locality name "%s" to "%s" and moved ' % \
                       (' '.join(org_loc_name_list),' '.join(loc_name_list)) \
                       + '"%s" to street type as "%s"' % \
                       (first_loc_name[:5],
                        ' '.join(address_dict['street_type'])))

    # Finally do some tests on the output fields  - - - - - - - - - - - - - - -
    #
    for (out_field, val_list) in address_dict.items():

      # Check if an output field value list has more than three elements
      #
      if (len(val_list) > 3):
        logging.warn('Output field "%s" contains more than ' % (out_field) + \
                   'three elements: %s' % (str(val_list)))

      # Check if a 'number' output field value list has more than one element,
      # and also check if it only contains numbers
      #
      if (out_field.endswith('number')):  # Numbers only, not pre- or suffix
        if (len(val_list) > 1):
          logging.warn('More than one element in output field "%s": %s' % \
                       (out_field, str(val_list)))
        for num in val_list:
          if (not num.isdigit()):  # Element contains not only digits
            logging.warn('Output field "%s" contains no number: %s' % \
                         (out_field, str(val_list)))

      # Check if type element contain one word only and if it's a known type
      # word
      #
      if (out_field.endswith('type')):  # street_type, flat_type, level_type
        if (len(val_list) > 1):
          logging.warn('More than element in output field "%s": %s' % \
                       (out_field, str(val_list)))
        for word in val_list:
          word_tuple = (word,)  # Make it a tuple
          if (self.tag_table.has_key(word_tuple)):
            word_tag = self.tag_table[word_tuple][1]  # Get it's tag

            if (out_field.startswith('street') and ('WT' not in word_tag)) or \
               (out_field.startswith('flat') and ('FT' not in word_tag)) or \
               (out_field.startswith('level') and ('LT' not in word_tag)):
              logging.warn('Type word "%s" in output field "%s" does not ' % \
                           (word, out_field) + 'have a correct "type" ' + \
                           'tag: %s' % (word_tag))

    # Check if 'qualifier' elements only contain known qualifier words  - - - -
    #
    if (address_dict.has_key('wayfare_qualifier')):  # Check wayfare qualifier
      wf_quali = address_dict['wayfare_qualifier']
      for wq in wf_quali:
        if ((not self.tag_table.has_key((wq,))) or \
            ((self.tag_table.has_key((wq,)) and \
             (self.tag_table[(wq,)][1].find('LQ') < 0)))):
          logging.warn('Wayfare qualifier word is not known: %s' % (str(wq)))
          break  # Exit for loop

    if (address_dict.has_key('locality_qualifier')): # Check locality qualifier
      loc_quali = address_dict['locality_qualifier']
      for lq in loc_quali:
        if ((not self.tag_table.has_key((lq,))) or \
            ((self.tag_table.has_key((lq,)) and \
             (self.tag_table[(lq,)][1].find('LQ') < 0)))):
          logging.warn('Locality qualifier word is not known: %s' % (str(lq)))
          break  # Exit for loop

    # Check if type, postcode state abbreviation states do not contain
    # duplicate words
    #
    for check_state in ['state_abbrev', 'postcode', 'street_type', 'flat_type',
                        'level_type']:
      if (address_dict.has_key(check_state)):
        if (len(address_dict[check_state]) > 1):
          address_dict[check_state] = list(set(address_dict[check_state]))

    # Finally build the address list - - - - - - - - - - - - - - - - - - - - -
    #
    address_list = 27*[[]]

    logging.info('Final segmented address:')
    i = 0
    for state_name in ['building_name','post_address_type',
                       'post_address_number','lot_number_prefix',
                       'lot_number','lot_number_suffix','flat_number_prefix',
                       'flat_number','flat_number_suffix','flat_type',
                       'level_number_prefix','level_number',
                       'level_number_suffix','level_type',
                       'number_first_prefix','number_first',
                       'number_first_suffix','number_last_prefix',
                       'number_last','number_last_suffix','street_name',
                       'street_suffix','street_type','locality_name',
                       'postcode','state_abbrev','country']:

      if (state_name in address_dict):
        address_list[i] = address_dict[state_name]

      logging.info('  %25s: %s' % (state_name, str(address_list[i])))
      i += 1

    logging.info('  Best state sequence: %s' % (str(best_state_seq)))
    logging.info('  Normalised HMM probability: %f' % (norm_max_prob))

    return [address_list, best_state_seq, norm_max_prob]

  # ---------------------------------------------------------------------------

  def finalise(self):
    """If the HMM sequence probability file is given write it out.
    """

    if (self.hmm_seq_prob_file != None):
      self.__write_hmm_seq_prob_file__('Address HMM')

# =============================================================================
