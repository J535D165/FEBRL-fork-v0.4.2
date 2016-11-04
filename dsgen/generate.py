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
# The Original Software is: "generate.py"
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

"""Module generate.py - Auxiliary program to create records using various
                        frequency tables and introduce duplicates with errors.

   USAGE:
     python generate.py [output_file] [num_originals] [num_duplicates]
                        [max_duplicate_per_record] [max_modification_per_field]
                        [max_modification_per_record] [distribution]

   ARGUMENTS:
     output_file                  Name of the output file (currently this is a
                                  CSV file).
     num_originals                Number of original records to be created.
     num_duplicates               Number of duplicate records to be created
                                  (maximum number is 9).
     max_duplicate_per_record     The maximal number of duplicates that can be
                                  created for one original record.
     max_modification_per_field   The maximum number of modifications per field
     max_modification_per_record  The maximum number of modifications per
                                  record.
     distribution                 The probability distribution used to create
                                  the duplicates (i.e the number of duplicates
                                  for one original).
                                  Possible are: - uniform
                                                - poisson
                                                - zipf

   DESCRIPTION:
     This program can be used to create a data set with records that contain
     randomly created names and addresses (using frequency files), dates,
     phone numbers, and identifier numbers. Duplicate records will then be
     created following a given probability distribution, with different single
     errors being introduced.

     Various parameters on how theses duplicates are created can be given
     within the program, see below.

     New: It is possible to load dictionaries (look-up table) with misspellings
          that will be used to replace a correct word with a randomly chosen
          misspelling. A user can easily customise this misspelling files.

   TODO:
     - add substitution matrix with character substitution probabilities
       (instead of keyboard based substitutions).

     - Improve performance (loading and creating frequency tables)

     - for each field have a counter num_modifcations in the field dictionary

     - do swap between field first (count as 2 rec. modifications)

     - Allow various probability distributions for fields of type 'date' and
       'iden' (using a new keyword in field dictionaries).

     - Try to find real world error distributions for typographical errors and
       integrate them into the random error creation

     - Add random word spilling between fields (similar to field swapping)
"""

# =============================================================================
# Imports go here

import math
import os
import random
import sets
import string
import sys
import time

# Set this flag to True for verbose output, otherwise to False - - - - - - - -
#
VERBOSE_OUTPUT = True

# =============================================================================
#
# For each field (attribute), a dictionary has to be defined with the following
# keys (probabilities can have values between 0.0 and 1.0, or they can be
# missing - in which case it is assumed they have a value of 0.0):
# - name           The field name to be used when a header is written into the
#                  output file.
# - type           The type of the field. Possible are:
#                  'freq'  (for fields that use a frequency table with field
#                           values)
#                  'date'  (for date fields in a certain range)
#                  'phone' (for phone numbers)
#                  'ident' (for numerical identifier fields in a certain range)
# - char_range     The range of random characters that can be introduced. Can
#                  be one of 'alpha', 'digit', or 'alphanum'.
#
# For fields of type 'freq' the following keys must be given:
# - freq_file      The name of a frequency file.
# - misspell_file  The name of a misspellings file.
#
# For fields of type 'date' the following keys must be given:
# - start_date     A start date, must be a tuple (day,month,year).
# - end_date       A end date, must be a tuple (day,month,year).
#
# For fields of type 'phone' the following keys must be given:
# - area_codes     A list with possible area codes (as strings).
# - num_digits     The number of digits in the phone numbers (without the area
#                  code).
#
# For fields of type 'ident' the following keys must be given:
# - start_id       A start identification number.
# - end_id         An end identification number.
#
# For all fields the following keys must be given:
# - select_prob    Probability of selecting a field for introducing one or
#                  more modifications (set this to 0.0 if no modifications
#                  should be introduced into this field ever). Note: The sum
#                  of these select probabilities over all defined fields must
#                  be 100.
# - misspell_prob  Probability to swap an original value with a randomly
#                  chosen misspelling from the corresponding misspelling
#                  dictionary (can only be set to larger than 0.0 if such a
#                  misspellings dictionary is defined for the given field).
# - ins_prob       Probability to insert a character into a field value.
# - del_prob       Probability to delete a character from a field value.
# - sub_prob       Probability to substitute a character in a field value with
#                  another character.
# - trans_prob     Probability to transpose two characters in a field value.
# - val_swap_prob  Probability to swap the value in a field with another
#                  (randomly selected) value for this field (taken from this
#                  field's look-up table).
# - wrd_swap_prob  Probability to swap two words in a field (given there are
#                  at least two words in a field).
# - spc_ins_prob   Probability to insert a space into a field value (thus
#                  splitting a word).
# - spc_del_prob   Probability to delete a space (if available) in a field (and
#                  thus merging two words).
# - miss_prob      Probability to set a field value to missing (empty).
# - new_val_prob   Probability to insert a new value given the original value
#                  was empty.
#
# Note: The sum over the probabilities ins_prob, del_prob, sub_prob,
#       trans_prob, val_swap_prob, wrd_swap_prob, spc_ins_prob, spc_del_prob,
#       and miss_prob for each defined field must be 1.0; or 0.0 if no
#       modification should be done at all on a given field.
#
# =============================================================================
# Comments about typographical errors and misspellings found in the literature:
#
# Damerau 1964: - 80% are single errors: insert, delete, substitute or
#                                        transpose
#               - Statistic given: 567/964 (59%) substitutions
#                                  153/964 (16%) deletions
#                                   23/964 ( 2%) transpositions
#                                   99/964 (10%) insertions
#                                  122/964 (13%) multiple errors
#
# Hall 1980: - OCR and other automatic devices introduce similar errors of
#              substitutions, deletions and insertions, but not transpositions;
#              frequency and type of errors are characteristics of the device.
#
# Pollock/Zamora 1984: - OCR output contains almost exclusively substitution
#                        errors which ordinarily account for less than 20% of
#                        key boarded misspellings.
#                      - 90-95% of misspellings in raw keyboarding typically
#                        only contain one error.
#                      - Only 7.8% of the first letter of misspellings were
#                        incorrect, compared to 11.7% of the second and 19.2%
#                        of the third.
#                      - Generally assumed that vowels are less important than
#                        consonants.
#                      - The frequency of a misspelling seems to be determined
#                        more by the frequency of it's parent word than by the
#                        difficulty of spelling it.
#                      - Most errors are mechanical (typos), not the result of
#                        poor spelling.
#                      - The more frequent a letter, the more likely it is to
#                        be miskeyed.
#                      - Deletions are similar frequent than transpositions,
#                        but more frequent than insertions and again more
#                        frequent than substitutions.
#
# Pollock/Zamora 1983: - Study of 50,000 nonword errors, 3-4 character
#                        misspellings constitute only 9.2% of total
#                        misspellings, but they generate 40% of miscorrections.
#
# Peterson 1986: In two studies:
#                - Transpose two letters:  2.6% / 13.1%
#                - Insert letter:         18.7% / 20.3%
#                - Delete letter:         31.6% / 34.4%
#                - Substitute letter:     40.0% / 26.9%
#
# Kukich 1990: - Over 63% of errors in TDD conversations occur in words of
#                length 2, 3 or 4.
#
# Kukich 1992: - 13% of non-word spelling errors in a 40,000 corpus of typed
#                conversations involved merging of two words, 2% splitting a
#                word (often at valid forms, "forgot" -> "for got").
#              - Most misspellings seem to be within two characters in length
#                of the correct spelling.
#
# =============================================================================
# Other comments:
#
# - Intuitively, one can assume that long and unfrequent words are more likely
#   to be misspelt than short and common words.
#
# =============================================================================

givenname_dict = {'name':'given_name',
                  'type':'freq',
            'char_range':'alpha',
#             'freq_file':'data'+os.sep+'givenname-freq.csv',
             'freq_file':'data-org'+os.sep+'givenname.csv',
           'select_prob':0.10,
         'misspell_file':'data'+os.sep+'givenname-misspell.tbl',
         'misspell_prob':0.30,
              'ins_prob':0.05,
              'del_prob':0.15,
              'sub_prob':0.35,
            'trans_prob':0.05,
         'val_swap_prob':0.02,
         'wrd_swap_prob':0.02,
          'spc_ins_prob':0.01,
          'spc_del_prob':0.01,
             'miss_prob':0.02,
          'new_val_prob':0.02}

surname_dict = {'name':'surname',
                'type':'freq',
          'char_range':'alpha',
#           'freq_file':'data'+os.sep+'surname-freq.csv',
           'freq_file':'data-org'+os.sep+'surname.csv',
         'select_prob':0.15,
       'misspell_file':'data'+os.sep+'surname-misspell.tbl',
       'misspell_prob':0.30,
            'ins_prob':0.10,
            'del_prob':0.10,
            'sub_prob':0.35,
          'trans_prob':0.04,
       'val_swap_prob':0.02,
       'wrd_swap_prob':0.02,
        'spc_ins_prob':0.01,
        'spc_del_prob':0.02,
           'miss_prob':0.02,
        'new_val_prob':0.02}

streetnumber_dict = {'name':'street_number',
                     'type':'freq',
               'char_range':'digit',
#                'freq_file':'data'+os.sep+'streetnumber-freq.csv',
                'freq_file':'data-org'+os.sep+'streetnumber.csv',
              'select_prob':0.10,
                 'ins_prob':0.10,
                 'del_prob':0.15,
                 'sub_prob':0.60,
               'trans_prob':0.05,
            'val_swap_prob':0.05,
            'wrd_swap_prob':0.01,
             'spc_ins_prob':0.00,
             'spc_del_prob':0.00,
                'miss_prob':0.02,
             'new_val_prob':0.02}

address1_dict = {'name':'address_1',
                 'type':'freq',
           'char_range':'alpha',
#            'freq_file':'data'+os.sep+'address1-freq.csv',
            'freq_file':'data-org'+os.sep+'address1.csv',
          'select_prob':0.10,
             'ins_prob':0.10,
             'del_prob':0.15,
             'sub_prob':0.55,
           'trans_prob':0.05,
        'val_swap_prob':0.02,
        'wrd_swap_prob':0.03,
         'spc_ins_prob':0.02,
         'spc_del_prob':0.03,
            'miss_prob':0.04,
         'new_val_prob':0.01}

# Address 2 contains property and institution names - only use rarely
# (set missing probability to a high value)
#
address2_dict = {'name':'address_2',
                 'type':'freq',
           'char_range':'alpha',
#            'freq_file':'data'+os.sep+'address2-freq.csv',
            'freq_file':'data-org'+os.sep+'address2.csv',
          'select_prob':0.10,
             'ins_prob':0.04,
             'del_prob':0.04,
             'sub_prob':0.10,
           'trans_prob':0.02,
        'val_swap_prob':0.03,
        'wrd_swap_prob':0.04,
         'spc_ins_prob':0.02,
         'spc_del_prob':0.01,
            'miss_prob':0.60,
         'new_val_prob':0.10}

suburb_dict = {'name':'suburb',
               'type':'freq',
         'char_range':'alpha',
#          'freq_file':'data'+os.sep+'suburb-freq.csv',
          'freq_file':'data-org'+os.sep+'suburb.csv',
        'select_prob':0.10,
      'misspell_file':'data'+os.sep+'suburb-misspell.tbl',
      'misspell_prob':0.40,
           'ins_prob':0.10,
           'del_prob':0.15,
           'sub_prob':0.22,
         'trans_prob':0.04,
      'val_swap_prob':0.01,
      'wrd_swap_prob':0.02,
       'spc_ins_prob':0.02,
       'spc_del_prob':0.02,
          'miss_prob':0.01,
       'new_val_prob':0.01}

postcode_dict = {'name':'postcode',
                 'type':'freq',
           'char_range':'digit',
#            'freq_file':'data'+os.sep+'postcode-freq.csv',
            'freq_file':'data-org'+os.sep+'postcode.csv',
          'select_prob':0.05,
             'ins_prob':0.00,
             'del_prob':0.00,
             'sub_prob':0.35,
           'trans_prob':0.60,
        'val_swap_prob':0.03,
        'wrd_swap_prob':0.00,
         'spc_ins_prob':0.00,
         'spc_del_prob':0.00,
            'miss_prob':0.01,
         'new_val_prob':0.01}

state_dict = {'name':'state',
              'type':'freq',
        'char_range':'alpha',
#         'freq_file':'data'+os.sep+'state-freq.csv',
         'freq_file':'data-org'+os.sep+'state.csv',
       'select_prob':0.05,
          'ins_prob':0.10,
          'del_prob':0.10,
          'sub_prob':0.55,
        'trans_prob':0.02,
     'val_swap_prob':0.03,
     'wrd_swap_prob':0.00,
      'spc_ins_prob':0.00,
      'spc_del_prob':0.00,
         'miss_prob':0.10,
      'new_val_prob':0.10}

dob_dict = {'name':'date_of_birth',
            'type':'date',
      'char_range':'digit',
      'start_date':(01,01,1900),
        'end_date':(31,12,1999),
     'select_prob':0.10,
        'ins_prob':0.00,
        'del_prob':0.00,
        'sub_prob':0.50,
      'trans_prob':0.30,
   'val_swap_prob':0.05,
   'wrd_swap_prob':0.00,
    'spc_ins_prob':0.00,
    'spc_del_prob':0.00,
       'miss_prob':0.10,
    'new_val_prob':0.05}

age_dict = {'name':'age',
            'type':'freq',
      'char_range':'digit',
       'freq_file':'data'+os.sep+'age-freq.csv',
     'select_prob':0.05,
        'ins_prob':0.00,
        'del_prob':0.00,
        'sub_prob':0.30,
      'trans_prob':0.20,
   'val_swap_prob':0.20,
   'wrd_swap_prob':0.00,
    'spc_ins_prob':0.00,
    'spc_del_prob':0.00,
       'miss_prob':0.20,
    'new_val_prob':0.10}

phonenum_dict = {'name':'phone_number',
                 'type':'phone',
           'char_range':'digit',
           'area_codes':['02','03','04','07','08'],  # Australian area codes
           'num_digits':8,                    # For Australian phone numbers
          'select_prob':0.05,
             'ins_prob':0.00,
             'del_prob':0.00,
             'sub_prob':0.40,
           'trans_prob':0.30,
        'val_swap_prob':0.15,
        'wrd_swap_prob':0.00,
         'spc_ins_prob':0.00,
         'spc_del_prob':0.00,
            'miss_prob':0.05,
         'new_val_prob':0.10}

ssid_dict = {'name':'soc_sec_id',
             'type':'ident',
       'char_range':'digit',
         'start_id':1000000,
           'end_id':9999999,
      'select_prob':0.05,
         'ins_prob':0.00,
         'del_prob':0.00,
         'sub_prob':0.50,
       'trans_prob':0.40,
    'val_swap_prob':0.10,
    'wrd_swap_prob':0.00,
     'spc_ins_prob':0.00,
     'spc_del_prob':0.00,
        'miss_prob':0.00,
     'new_val_prob':0.00}

# Create a field which can be used for blocking (generate values 0 to 9 which
# are not modified in duplicates).
#
blocking_dict = {'name':'blocking_number',
                 'type':'ident',
           'char_range':'digit',
             'start_id':0,
               'end_id':10,
          'select_prob':0.00,
             'ins_prob':0.00,
             'del_prob':0.00,
             'sub_prob':0.00,
           'trans_prob':0.00,
        'val_swap_prob':0.00,
        'wrd_swap_prob':0.00,
         'spc_ins_prob':0.00,
         'spc_del_prob':0.00,
            'miss_prob':0.00,
         'new_val_prob':0.00}

# -----------------------------------------------------------------------------
# Probabilities (between 0.0 and 1.0) for swapping values between two fields.
# Use field names as defined in the field directories (keys 'name').

field_swap_prob = {('address_1',  'address_2'):0.02,
                   ('given_name', 'surname'):  0.05,
                   ('postcode',   'suburb'):   0.01}

# -----------------------------------------------------------------------------
# Probabilities (between 0.0 and 1.0) for creating a typographical error (a new
# character) in the same row or the same column. This is used in the random
# selection of a new character in the 'sub_prob' (substitution of a character
# in a field).

single_typo_prob = {'same_row':0.40,
                    'same_col':0.30}

# -----------------------------------------------------------------------------
# Now add all field dictionaries into a list according to how they should be
# saved in the output file.

field_list = [givenname_dict, surname_dict, streetnumber_dict, address1_dict,
              address2_dict, suburb_dict, postcode_dict, state_dict,
              dob_dict, age_dict, phonenum_dict, ssid_dict, blocking_dict]

# -----------------------------------------------------------------------------
# Flag for writing a header line (keys 'name' of field dictionaries).

save_header = True  # Set to 'False' if no header should be written

# -----------------------------------------------------------------------------
# String to be inserted for missing values.

missing_value = ''

# =============================================================================
# Nothing to be changed below here
# =============================================================================

# Initialise random number generator  - - - - - - - - - - - - - - - - - - - - -
#
random.seed()

# =============================================================================
# Functions used by the main program come here

def error_position(input_string, len_offset):
  """A function that randomly calculates an error position within the given
     input string and returns the position as integer number 0 or larger.

     The argument 'len_offset' can be set to an integer (e.g. -1, 0, or 1) and
     will give an offset relative to the string length of the maximal error
     position that can be returned.

     Errors do not likely appear at the beginning of a word, so a gauss random
     distribution is used with the mean being one position behind half the
     string length (and standard deviation 1.0)
  """

  str_len = len(input_string)
  max_return_pos = str_len - 1 + len_offset  # Maximal position to be returned

  if (str_len == 0):
    return None  # Empty input string

  mid_pos = (str_len + len_offset) / 2 + 1

  random_pos = random.gauss(float(mid_pos), 1.0)
  random_pos = max(0,int(round(random_pos)))  # Make it integer and 0 or larger

  return min(random_pos, max_return_pos)

# -----------------------------------------------------------------------------

def error_character(input_char, char_range):
  """A function which returns a character created randomly. It uses row and
     column keyboard dictionaires.
  """

  # Keyboard substitutions gives two dictionaries with the neigbouring keys for
  # all letters both for rows and columns (based on ideas implemented by
  # Mauricio A. Hernandez in his dbgen).
  #
  rows = {'a':'s',  'b':'vn', 'c':'xv', 'd':'sf', 'e':'wr', 'f':'dg', 'g':'fh',
          'h':'gj', 'i':'uo', 'j':'hk', 'k':'jl', 'l':'k',  'm':'n',  'n':'bm',
          'o':'ip', 'p':'o',  'q':'w',  'r':'et', 's':'ad', 't':'ry', 'u':'yi',
          'v':'cb', 'w':'qe', 'x':'zc', 'y':'tu', 'z':'x',
          '1':'2',  '2':'13', '3':'24', '4':'35', '5':'46', '6':'57', '7':'68',
          '8':'79', '9':'80', '0':'9'}

  cols = {'a':'qzw','b':'gh', 'c':'df', 'd':'erc','e':'d', 'f':'rvc','g':'tbv',
          'h':'ybn','i':'k',  'j':'umn','k':'im', 'l':'o',  'm':'jk', 'n':'hj',
          'o':'l',  'p':'p',  'q':'a',  'r':'f',  's':'wxz','t':'gf',  'u':'j',
          'v':'fg', 'w':'s',  'x':'sd', 'y':'h',  'z':'as'}

  rand_num = random.random()  # Create a random number between 0 and 1

  if (char_range == 'digit'):

    # A randomly chosen neigbouring key in the same keyboard row
    #
    if (input_char.isdigit()) and (rand_num <= single_typo_prob['same_row']):
      output_char = random.choice(rows[input_char])
    else:
      choice_str =  string.replace(string.digits, input_char, '')
      output_char = random.choice(choice_str)  # A randomly choosen digit

  elif (char_range == 'alpha'):

    # A randomly chosen neigbouring key in the same keyboard row
    #
    if (input_char.isalpha()) and (rand_num <= single_typo_prob['same_row']):
      output_char = random.choice(rows[input_char])

    # A randomly chosen neigbouring key in the same keyboard column
    #
    elif (input_char.isalpha()) and \
       (rand_num <= (single_typo_prob['same_row'] + \
                     single_typo_prob['same_col'])):
      output_char = random.choice(cols[input_char])
    else:
      choice_str =  string.replace(string.lowercase, input_char, '')
      output_char = random.choice(choice_str)  # A randomly choosen letter

  else:  # Both letters and digits possible

    # A randomly chosen neigbouring key in the same keyboard row
    #
    if (rand_num <= single_typo_prob['same_row']):
      if (input_char in rows):
        output_char = random.choice(rows[input_char])
      else:
        choice_str =  string.replace(string.lowercase+string.digits, \
                                     input_char, '')
        output_char = random.choice(choice_str)  # A randomly choosen character

    # A randomly chosen neigbouring key in the same keyboard column
    #
    elif (rand_num <= (single_typo_prob['same_row'] + \
                       single_typo_prob['same_col'])):
      if (input_char in cols):
        output_char = random.choice(cols[input_char])
      else:
        choice_str =  string.replace(string.lowercase+string.digits, \
                                     input_char, '')
        output_char = random.choice(choice_str)  # A randomly choosen character

    else:
      choice_str =  string.replace(string.lowercase+string.digits, \
                                   input_char, '')
      output_char = random.choice(choice_str)  # A randomly choosen character

  return output_char

# -----------------------------------------------------------------------------

# Some simple funcions used for date conversions follow
# (based on functions from the 'normalDate.py' module by Jeff Bauer, see:
# http://starship.python.net/crew/jbauer/normalDate/)

days_in_month = [[31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], \
                 [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]]

def first_day_of_year(year):
  """Calculate the day number (relative to 1 january 1900) of the first day in
     the given year.
  """

  if (year == 0):
    print 'Error: A year value of 0 is not possible'
    raise Exception

  elif (year < 0):
    first_day = (year * 365) + int((year - 1) / 4) - 693596
  else:  # Positive year
    leap_adj = int ((year + 3) / 4)
    if (year > 1600):
      leap_adj = leap_adj - int((year + 99 - 1600) / 100) + \
                 int((year + 399 - 1600) / 400)

    first_day = year * 365 + leap_adj - 693963

    if (year > 1582):
      first_day -= 10

  return first_day

# -----------------------------------------------------------------------------

def is_leap_year(year):
  """Determine if the given year is a leap year. Returns 0 (no) or 1 (yes).
  """

  if (year < 1600):
    if ((year % 4) != 0):
      return 0
    else:
      return 1

  elif ((year % 4) != 0):
    return 0

  elif ((year % 100) != 0):
    return 1

  elif ((year % 400) != 0):
    return 0

  else:
    return 1

# -----------------------------------------------------------------------------

def epoch_to_date(daynum):
  """Convert an epoch day number into a date [day, month, year], with
     day, month and year being strings of length 2, 2, and 4, respectively.
     (based on a function from the 'normalDate.py' module by Jeff Bauer, see:
     http://starship.python.net/crew/jbauer/normalDate/)

  USAGE:
    [year, month, day] = epoch_to_date(daynum)

  ARGUMENTS:
    daynum  A integer giving the epoch day (0 = 1 January 1900)

  DESCRIPTION:
    Function for converting a number of days (integer value) since epoch time
    1 January 1900 (integer value) into a date tuple [day, month, year].

  EXAMPLES:
    [day, month, year] = epoch_to_date(0)      # returns ['01','01','1900']
    [day, month, year] = epoch_to_date(37734)  # returns ['25','04','2003']
  """

  if (not (isinstance(daynum, int) or isinstance(daynum, long))):
    print 'Error: Input value for "daynum" is not of integer type: %s' % \
          (str(daynum))
    raise Exception

  if (daynum >= -115860):
    year = 1600 + int(math.floor((daynum + 109573) / 365.2425))
  elif (daynum >= -693597):
    year = 4 + int(math.floor((daynum + 692502) / 365.2425))
  else:
    year = -4 + int(math.floor((daynum+695058) / 365.2425))

  days = daynum - first_day_of_year(year) + 1

  if (days <= 0):
    year -= 1
    days = daynum - first_day_of_year(year) + 1

  days_in_year = 365 + is_leap_year(year)  # Adjust for a leap year

  if (days > days_in_year):
    year += 1
    days = daynum - first_day_of_year(year) + 1

  # Add 10 days for dates between 15 October 1582 and 31 December 1582
  #
  if (daynum >= -115860) and (daynum <= -115783):
    days += 10

  day_count = 0
  month = 12
  leap_year_flag = is_leap_year(year)

  for m in range(12):
    day_count += days_in_month[leap_year_flag][m]
    if (day_count >= days):
      month = m + 1
      break

  # Add up the days in the prior months
  #
  prior_month_days = 0
  for m in range(month-1):
    prior_month_days += days_in_month[leap_year_flag][m]

  day = days - prior_month_days

  day_str =   string.zfill(str(day),2)  # Add '0' if necessary
  month_str = string.zfill(str(month),2)  # Add '0' if necessary
  year_str =  str(year)  # Is always four digits long

  return [day_str, month_str, year_str]

# -----------------------------------------------------------------------------

def date_to_epoch(day, month, year):
  """ Convert a date [day, month, year] into an epoch day number.
     (based on a function from the 'normalDate.py' module by Jeff Bauer, see:
     http://starship.python.net/crew/jbauer/normalDate/)

  USAGE:
    daynum = date_to_epoch(year, month, day)

  ARGUMENTS:
    day    Day value (string or integer number)
    month  Month value (string or integer number)
    year   Year value (string or integer number)

  DESCRIPTION:
    Function for converting a date into a epoch day number (integer value)
    since 1 january 1900.

  EXAMPLES:
    day = date_to_epoch('01', '01', '1900')  # returns 0
    day = date_to_epoch('25', '04', '2003')  # returns 37734
  """

  # Convert into integer values
  #
  try:
    day_int = int(day)
  except:
    print 'Error: "day" value is not an integer'
    raise Exception
  try:
    month_int = int(month)
  except:
    print 'Error: "month" value is not an integer'
    raise Exception
  try:
    year_int = int(year)
  except:
    print 'Error: "year" value is not an integer'
    raise Exception

  # Test if values are within range
  #
  if (year_int <= 1000):
    print 'Error: Input value for "year" is not a positive integer ' + \
          'number: %i' % (year)
    raise Exception

  leap_year_flag = is_leap_year(year_int)

  if (month_int <= 0) or (month_int > 12):
    print 'Error: Input value for "month" is not a possible day number: %i' % \
          (month)
    raise Exception

  if (day_int <= 0) or (day_int > days_in_month[leap_year_flag][month_int-1]):
    print 'Error: Input value for "day" is not a possible day number: %i' % \
          (day)
    raise Exception

  days = first_day_of_year(year_int) + day_int - 1

  for m in range(month_int-1):
    days += days_in_month[leap_year_flag][m]

  if (year_int == 1582):
    if (month_int > 10) or ((month_int == 10) and (day_int > 4)):
      days -= 10

  return days

# -----------------------------------------------------------------------------

def load_misspellings_dict(misspellings_file_name):
  """Load a look-up table containing misspellings for common words, which can
     be used to introduce realistic errors.

     Returns a dictionary where the keys are the correct spellings and the
     values are a list of one or more misspellings.
  """

  # Open file and read all lines into a list
  #
  try:
    f = open(misspellings_file_name, 'r')
  except:
    print 'Error: Can not read from misspellings file "%s"' % \
          (misspellings_file_name)
    raise IOError

  file_data = f.readlines()  # Read complete file
  f.close()

  misspell_dict = {}

  key = None  # Start with a non-existing eky word (correct word)

  # Now process all lines - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  for line in file_data:
    l = line.strip()  # Remove line separators

    if (len(l) > 0) and (l[0] != '#'):  # Not empty line and not comment

      ll = l.split(':')  # Separate key from values

      if (ll[0] == '') and (len(ll) > 1):
        ll = ll[1:]

      if (len(ll) == 2):  # Line contains a key - - - - - - - - - - - - - - - -

        key = ll[0].strip().lower()  # Get key, make lower and strip spaces

        if (key == ''):
          print 'This should not happen: "%s"' % (l)
          raise Exception

        vals = ll[1].strip().lower() # Get values in a string

        if (vals == ''):
          print 'Error: No misspellings given for "%s" in line: "%s"' % \
                (key, l)
          raise Exception

        val_list = vals.split(',')
        val_set = sets.Set()
        for val in val_list:
          if (val != ''):
            val_set.add(val.strip())  # Remove all spaces

        # Check that all misspellings are different from the original
        #
        if (key in val_set):
          print 'Error: A misspelling is the same as the original value' + \
                ' "%s" in line: "%s"' % (key, l)
          raise Exception

        # Now insert into misspellings dictionary
        #
        key_val_set = misspell_dict.get(key, sets.Set())
        key_val_set = key_val_set.union(val_set)
        misspell_dict[key] = key_val_set

      elif (len(ll) == 1):  # Line contains only values - - - - - - - - - - - -

        if (key == None):
          print 'Error: No key (correct word) defined in line: "%s"' % (l)
          raise Exception

        vals = ll[0].lower() # Get values in a string
        val_list = vals.split(',')

        val_set = sets.Set()
        for val in val_list:
          if (val != ''):
            val_set.add(val.strip())  # Remove all spaces

        # Check that all misspellings are different from the original
        #
        if (key in val_set):
          print 'Error: A misspelling is the same as the original value' + \
                ' "%s" in line: "%s"' % (key, l)
          raise Exception

        # Now insert into misspellings dictionary
        #
        key_val_set = misspell_dict.get(key, sets.Set())
        key_val_set = key_val_set.union(val_set)
        misspell_dict[key] = key_val_set

      else:
        print 'error:Illegal line format in line: "%s"' % (l)
        raise Exception

  # Now convert all sets into lists - - - - - - - - - - - - - - - - - - - - -
  #
  for k in misspell_dict:
    misspell_dict[k] = list(misspell_dict[k])

  # print '  Length of misspellings dictionary: %d' % (len(misspell_dict))

  return misspell_dict

# -----------------------------------------------------------------------------

def random_select(prob_dist_list):
  """Randomly select one of the list entries (tuples of value and probability
     values).
  """

  rand_num = random.random()  # Random number between 0.0 and 1.0

  ind = -1
  while (prob_dist_list[ind][1] > rand_num):
    ind -= 1

  return prob_dist_list[ind][0]

# =============================================================================
# Start main program

if (len(sys.argv) != 8):
  print 'Seven arguments needed with %s:' % (sys.argv[0])
  print '  - Output file name'
  print '  - Number of original records'
  print '  - Number of duplicate records'
  print '  - Maximal number of duplicate records for one original record'
  print '  - Maximum number of modifications per field'
  print '  - Maximum number of modifications per record'
  print '  - Probability distribution for duplicates (uniform, poisson, zipf)'
  print 'All other parameters have to be set within the code'
  sys.exit()

output_file =           sys.argv[1]
num_org_records =       int(sys.argv[2])
num_dup_records =       int(sys.argv[3])
max_num_dups =          int(sys.argv[4])
max_num_field_modifi =  int(sys.argv[5])
max_num_record_modifi = int(sys.argv[6])
prob_distribution =     sys.argv[7][:3]

if (num_org_records <= 0):
  print 'Error: Number of original records must be positive'
  sys.exit()

if (num_dup_records < 0):
  print 'Error: Number of duplicate records must be zero or positive'
  sys.exit()

if (max_num_dups <= 0) or (max_num_dups > 9):
  print 'Error: Maximal number of duplicates per record must be positive ' + \
        'and less than 10'
  sys.exit()

if (max_num_field_modifi <= 0):
  print 'Error: Maximal number of modifications per field must be positive'
  sys.exit()

if (max_num_record_modifi <= 0):
  print 'Error: Maximal number of modifications per record must be positive'
  sys.exit()

if (max_num_record_modifi < max_num_field_modifi):
  print 'Error: Maximal number of modifications per record must be equal to'
  print '       or larger than maximal number of modifications per field'
  sys.exit()

if (prob_distribution not in ['uni', 'poi', 'zip']):
  print 'Error: Illegal probability distribution: %s' % (sys.argv[7])
  print '       Must be one of: "uniform", "poisson", or "zipf"'
  sys.exit()

# -----------------------------------------------------------------------------
# Check all user options within generate.py for validity
#
field_names = []  # Make a list of all field names

# A list of all probabilities to check ('select_prob' is checked separately)
#
prob_names = ['ins_prob','del_prob','sub_prob','trans_prob','val_swap_prob',
              'wrd_swap_prob','spc_ins_prob','spc_del_prob','miss_prob',
              'misspell_prob','new_val_prob']

select_prob_sum = 0.0  # Sum over all select probabilities

# Check if all defined field dictionaries have the necessary keys
#
i = 0  # Loop counter
for field_dict in field_list:

  if ('name' not in field_dict):
    print 'Error: No field name given for field dictionary'
    raise Exception
  elif (field_dict['name'] == 'rec_id'):
    print 'Error: Illegal field name "rec_id" (used for record identifier)'
    raise Exception
  else:
    field_names.append(field_dict['name'])

  if (field_dict.get('type','') not in ['freq','date','phone','ident']):
    print 'Error: Illegal or no field type given for field "%s": %s' % \
          (field_dict['name'], field_dict.get('type',''))
    raise Exception

  if (field_dict.get('char_range','') not in ['alpha', 'alphanum','digit']):
    print 'Error: Illegal or no random character range given for ' + \
          'field "%s": %s' % (field_dict['name'], \
                              field_dict.get('char_range',''))
    raise Exception

  if (field_dict['type'] == 'freq'):
    if (not field_dict.has_key('freq_file')):
      print 'Error: Field of type "freq" has no file name given'
      raise Exception

  elif (field_dict['type'] == 'date'):
    if (not (field_dict.has_key('start_date') and \
             field_dict.has_key('end_date'))):
      print 'Error: Field of type "date" has no start and/or end date given'
      raise Exception

    else:  # Process start and end date
      start_date = field_dict['start_date']
      end_date =   field_dict['end_date']

      start_epoch = date_to_epoch(start_date[0], start_date[1], start_date[2])
      end_epoch =   date_to_epoch(end_date[0], end_date[1], end_date[2])
      field_dict['start_epoch'] = start_epoch
      field_dict['end_epoch'] =   end_epoch
      field_list[i] = field_dict

  elif (field_dict['type'] == 'phone'):
    if (not (field_dict.has_key('area_codes') and \
             field_dict.has_key('num_digits'))):
      print 'Error: Field of type "phone" has no area codes and/or number ' + \
            'of digits given'
      raise Exception

    else:  # Process area codes and number of digits
      if (isinstance(field_dict['area_codes'],str)):  # Only one area code
        field_dict['area_codes'] = [field_dict['area_codes']]  # Make it a list
      if (not isinstance(field_dict['area_codes'],list)):
        print 'Error: Area codes given are not a string or a list: %s' % \
              (str(field_dict['area_codes']))
        raise Exception

      if (not isinstance(field_dict['num_digits'],int)):
        print 'Error: Number of digits given is not an integer: %s (%s)' % \
              (str(field_dict['num_digits']), type(field_dict['num_digits']))
        raise Exception

      field_list[i] = field_dict

  elif (field_dict['type'] == 'ident'):
    if (not (field_dict.has_key('start_id') and \
             field_dict.has_key('end_id'))):
      print 'Error: Field of type "iden" has no start and/or end ' + \
            'identification number given'
      raise Exception

  # Check all the probabilities for this field
  #
  if ('select_prob' not in field_dict):
    field_dict['select_dict'] = 0.0
  elif (field_dict['select_prob'] < 0.0) or (field_dict['select_prob'] > 1.0):
    print 'Error: Illegal value for select probability in dictionary for ' + \
          'field "%s": %f' % (field_dict['name'], field_dict['select_prob'])
  else:
    select_prob_sum += field_dict['select_prob']

  field_prob_sum = 0.0

  for prob in prob_names:
    if (prob not in field_dict):
      field_dict[prob] = 0.0
    elif (field_dict[prob] < 0.0) or (field_dict[prob] > 1.0):
      print 'Error: Illegal value for "%s" probability in dictionary for ' % \
            (prob) + 'field "%s": %f' % (field_dict['name'], field_dict[prob])
      raise Exception
    else:
      field_prob_sum += field_dict[prob]

  if (field_prob_sum > 0.0) and (abs(field_prob_sum - 1.0) > 0.001):
      print 'Error: Sum of probabilities for field "%s" is not 1.0: %f' % \
            (field_dict['name'], field_prob_sum)
      raise Exception

  # Create a list of field probabilities and insert into field dictionary
  #
  prob_list = []
  prob_sum =  0.0

  for prob in prob_names:
    prob_list.append((prob, prob_sum))
    prob_sum += field_dict[prob]

  field_dict['prob_list'] = prob_list
  field_list[i] = field_dict  # Store dictionary back into dictionary list

  i += 1

if (abs(select_prob_sum - 1.0) > 0.001):
  print 'Error: Field select probabilities do not sum to 1.0: %f' % \
        (select_prob_sum)
  raise Exception

# Create list of select probabilities - - - - - - - - - - - - - - - - - - - - -
#
select_prob_list = []
prob_sum =         0.0

for field_dict in field_list:
  select_prob_list.append((field_dict, prob_sum))
  prob_sum += field_dict['select_prob']

# -----------------------------------------------------------------------------
# Create a distribution for the number of duplicates for an original record
#
num_dup =  1
prob_sum = 0.0
prob_dist_list = [(num_dup, prob_sum)]

if (prob_distribution == 'uni'):  # Uniform distribution of duplicates - - - -

  uniform_val = 1.0 / float(max_num_dups)

  for i in range(max_num_dups-1):
    num_dup += 1
    prob_dist_list.append((num_dup, uniform_val+prob_dist_list[-1][1]))

elif (prob_distribution == 'poi'):  # Poisson distribution of duplicates - - -

  def fac(n):  # Factorial of an integer number (recursive calculation)
    if (n > 1.0):
      return n*fac(n - 1.0)
    else:
      return 1.0

  poisson_num = []  # A list of poisson numbers
  poisson_sum = 0.0  # The sum of all poisson number

  # The mean (lambda) for the poisson numbers
  #
  mean = 1.0 + (float(num_dup_records) / float(num_org_records))

  for i in range(max_num_dups):
    poisson_num.append((math.exp(-mean) * (mean ** i)) / fac(i))
    poisson_sum += poisson_num[-1]

  for i in range(max_num_dups):  # Scale so they sum up to 1.0
    poisson_num[i] = poisson_num[i] / poisson_sum

  for i in range(max_num_dups-1):
    num_dup += 1
    prob_dist_list.append((num_dup, poisson_num[i]+prob_dist_list[-1][1]))

elif (prob_distribution == 'zip'):  # Zipf distribution of duplicates - - - - -
  zipf_theta = 0.5

  denom = 0.0
  for i in range(num_org_records):
    denom += (1.0 / (i+1) ** (1.0 - zipf_theta))

  zipf_c = 1.0 / denom
  zipf_num = []  # A list of Zipf numbers
  zipf_sum = 0.0  # The sum of all Zipf number

  for i in range(max_num_dups):
    zipf_num.append(zipf_c / ((i+1) ** (1.0 - zipf_theta)))
    zipf_sum += zipf_num[-1]

  for i in range(max_num_dups):  # Scale so they sum up to 1.0
    zipf_num[i] = zipf_num[i] / zipf_sum

  for i in range(max_num_dups-1):
    num_dup += 1
    prob_dist_list.append((num_dup, zipf_num[i]+prob_dist_list[-1][1]))

print
print 'Create %i original and %i duplicate records' % \
      (num_org_records, num_dup_records)
print '  Distribution of number of duplicates (maximal %i duplicates):' % \
      (max_num_dups)
print '  %s' % (prob_dist_list)

# -----------------------------------------------------------------------------
# Load frequency files and misspellings dictionaries
#
print
print 'Step 1: Load and process frequency tables and misspellings dictionaries'

freq_files = {}
freq_files_length = {}

i = 0  # Loop counter
for field_dict in field_list:
  field_name = field_dict['name']

  if (field_dict['type'] == 'freq'):  # Check for 'freq' field type

    file_name = field_dict['freq_file']  # Get the corresponding file name

    if (file_name != None):
      try:
        fin = open(file_name)  # Open file for reading
      except:
        print '  Error: Can not open frequency file %s' % (file_name)
        raise Exception
      value_list = []  # List with all values of the frequency file

      for line in fin:
        line = line.strip()
        line_list = line.split(',')
        if (len(line_list) != 2):
          print '  Error: Illegal format in  frequency file %s: %s' % \
                (file_name, line)
          raise Exception

        line_val =  line_list[0].strip()
        line_freq = int(line_list[1])

        # Append value as many times as given in frequency file
        #
        new_list = [line_val]* line_freq
        value_list += new_list

      random.shuffle(value_list)  # Randomly shuffle the list of values

      freq_files[field_name] = value_list
      freq_files_length[field_name] = len(value_list)

      if (VERBOSE_OUTPUT == True):
        print '  Loaded frequency file for field "%s" from file: %s' % \
              (field_dict['name'], file_name)
        print

    else:
      print '  Error: No file name defined for frequency field "%s"' % \
            (field_dict['name'])
      raise Exception

  if ('misspell_file' in field_dict):  # Load misspellings dictionary file
    misspell_file_name = field_dict['misspell_file']
    field_dict['misspell_dict'] = load_misspellings_dict(misspell_file_name)

    if (VERBOSE_OUTPUT == True):
      print '  Loaded misspellings dictionary for field "%s" from file: "%s' \
            % (field_dict['name'], misspell_file_name)
      print

    field_list[i] = field_dict  # Store dictionary back into dictionary list

  i += 1

# -----------------------------------------------------------------------------
# Create original records
#
print
print 'Step 2: Create original records'
print

org_rec = {}  # Dictionary for original records
all_rec_set = sets.Set()  # Set of all records (without identifier) used for
                          # checking that all records are different
rec_cnt = 0

while (rec_cnt < num_org_records):
  rec_id = 'rec-%i-org' % (rec_cnt)  # The records identifier

  rec_dict = {'rec_id':rec_id}  # Save record identifier

  # Now randomly create all the fields in a record  - - - - - - - - - - - - - -
  #
  for field_dict in field_list:
    field_name = field_dict['name']

    # Randomly set field values to missing
    #
    if (random.random() <= field_dict['miss_prob']):
      rand_val = missing_value

    elif (field_dict['type'] == 'freq'):  # A frequency file based field
      rand_num = random.randint(0, freq_files_length[field_name]-1)
      rand_val = freq_files[field_name][rand_num]

    elif (field_dict['type'] == 'date'):  # A date field
      rand_num = random.randint(field_dict['start_epoch'], \
                                field_dict['end_epoch']-1)
      rand_date = epoch_to_date(rand_num)
      rand_val = rand_date[2]+rand_date[1]+rand_date[0]  # ISO format: yyyymmdd

    elif (field_dict['type'] == 'phone'):  # A phone number field
      area_code = random.choice(field_dict['area_codes'])
      max_digit = int('9'*field_dict['num_digits'])
      min_digit = int('1'*(int(1+round(field_dict['num_digits']/2.))))
      rand_num = random.randint(min_digit, max_digit)
      rand_val = area_code+' '+str(rand_num).zfill(field_dict['num_digits'])

    elif (field_dict['type'] == 'ident'):  # A identification number field
      rand_num = random.randint(field_dict['start_id'], \
                                field_dict['end_id']-1)
      rand_val = str(rand_num)

    if (rand_val != missing_value):  # Don't save missing values
      rec_dict[field_name] = rand_val

  # Create a string representation which can be used to check for uniqueness
  #
  rec_data = rec_dict.copy()  # Make a copy of the record dictionary
  del(rec_data['rec_id'])     # Remove the record identifier
  rec_list = rec_data.items()
  rec_list.sort()
  rec_str = str(rec_list)

  if (rec_str not in all_rec_set):  # Check if same record already created
    all_rec_set.add(rec_str)
    org_rec[rec_id] = rec_dict  # Insert into original records
    rec_cnt += 1

    # Print original record - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (VERBOSE_OUTPUT == True):
      print '  Original:'
      print '    Record ID         : %-30s' % (rec_dict['rec_id'])
      for field_name in field_names:
        print '    %-18s: %-30s' % (field_name, \
                                    rec_dict.get(field_name, missing_value))
      print

  else:
    if (VERBOSE_OUTPUT == True):
      print '***** Record "%s" already crated' % (rec_str)

# -----------------------------------------------------------------------------
# Create duplicate records
#
print
print 'Step 2: Create duplicate records'
print

dup_rec = {}  # Dictionary for duplicate records

org_rec_used = {}  # Dictionary with record IDs of original records used to
                   # create duplicates

rec_cnt = 0  # Record counter

while (rec_cnt < num_dup_records):

  # Find an original record that has so far not been used to create - - - - - -
  # duplicates
  #
  rand_rec_num = random.randint(0, num_org_records)
  org_rec_id = 'rec-%i-org' % (rand_rec_num)

  while (org_rec_id in org_rec_used) or (org_rec_id not in org_rec):
    rand_rec_num = random.randint(0, num_org_records)  # Get new record number
    org_rec_id = 'rec-%i-org' % (rand_rec_num)

  # Randomly choose how many duplicates to create from this record
  #
  num_dups = random_select(prob_dist_list)

  if (VERBOSE_OUTPUT == True):
    print '  Use record %s to create %i duplicates' % (org_rec_id, num_dups)

  org_rec_dict = org_rec[org_rec_id]  # Get the original record

  d = 0  # Loop counter for duplicates for this record

  # Loop to create duplicate records - - - - - - - - - - - - - - - - - - - - -
  #
  while (d < num_dups) and (rec_cnt < num_dup_records):

    # Create a duplicate of the original record
    #
    dup_rec_dict = org_rec_dict.copy()  # Make a copy of the original record
    dup_rec_id =             'rec-%i-dup-%i' % (rand_rec_num, d)
    dup_rec_dict['rec_id'] = dup_rec_id

    num_modif_in_record = 0  # Count the number of modifications in this record

    # Set the field modification counters to zero for all fields
    #
    field_mod_count_dict = {}

    for field_dict in field_list:
      field_mod_count_dict[field_dict['name']] = 0

    # Do random swapping between fields if two or more modifications in record
    #
    if (max_num_record_modifi > 1):

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Random swapping of values between a pair of field values
      #
      field_swap_pair_list = field_swap_prob.keys()
      random.shuffle(field_swap_pair_list)

      for field_pair in field_swap_pair_list:

        if (random.random() <= field_swap_prob[field_pair]) and \
           (num_modif_in_record <= (max_num_record_modifi-2)):

          fname_a, fname_b = field_pair

          # Make sure both fields are in the record dictionary
          #
          if (fname_a in dup_rec_dict) and (fname_b in dup_rec_dict):
            fvalue_a = dup_rec_dict[fname_a]
            fvalue_b = dup_rec_dict[fname_b]

            dup_rec_dict[fname_a] = fvalue_b  # Swap field values
            dup_rec_dict[fname_b] = fvalue_a

            num_modif_in_record += 2

            field_mod_count_dict[fname_a] = field_mod_count_dict[fname_a] + 1
            field_mod_count_dict[fname_b] = field_mod_count_dict[fname_b] + 1

            if (VERBOSE_OUTPUT == True):
              print '    Swapped fields "%s" and "%s": "%s" <-> "%s"' % \
                    (fname_a, fname_b, fvalue_a, fvalue_b)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Now introduce modifications up to the given maximal number

    while (num_modif_in_record < max_num_record_modifi):

      # Randomly choose a field
      #
      field_dict = random_select(select_prob_list)
      field_name = field_dict['name']

      # Make sure this field hasn't been modified already
      #
      while (field_mod_count_dict[field_name] == max_num_field_modifi):
        field_dict = random_select(select_prob_list)
        field_name = field_dict['name']

      if (field_dict['char_range'] == 'digit'):
        field_range = string.digits
      elif (field_dict['char_range'] == 'alpha'):
        field_range = string.lowercase
      elif (field_dict['char_range'] == 'alphanum'):
        field_range = string.digits+string.lowercase

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Randomly select the number of modifications to be done in this field
      # (and make sure we don't too many modifications in the record)
      #
      if (max_num_field_modifi == 1):
        num_field_mod_to_do = 1
      else:
        num_field_mod_to_do = random.randint(1, max_num_field_modifi)

      num_rec_mod_to_do = max_num_record_modifi - num_modif_in_record

      if (num_field_mod_to_do > num_rec_mod_to_do):
        num_field_mod_to_do = num_rec_mod_to_do

      if (VERBOSE_OUTPUT == True):
        print '    Choose field "%s" for %d modification' % \
              (field_name, num_field_mod_to_do)

      num_modif_in_field = 0  # Count the number of modifications in this field

      org_field_val = org_rec_dict.get(field_name, None) # Get original value

      # Loop over chosen number of modifications - - - - - - - - - - - - - - -
      #
      for m in range(num_field_mod_to_do):

        # Randomly choose a modification
        #
        mod_op = random_select(field_dict['prob_list'])

        old_field_val = dup_rec_dict.get(field_name, None)
        dup_field_val = old_field_val  # Modify this value

        # ---------------------------------------------------------------------
        # Do the selected modification

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Randomly choose a misspelling if the field value is found in the
        # misspellings dictionary
        #
        if (mod_op == 'misspell_prob') and ('misspell_dict' in field_dict) \
           and (old_field_val in field_dict['misspell_dict']):

          misspell_list = field_dict['misspell_dict'][old_field_val]

          if (len(misspell_list) == 1):
            dup_field_val = misspell_list[0]

          else:  # Randomly choose a value
            dup_field_val = random.choice(misspell_list)

          if (VERBOSE_OUTPUT == True):
            print '      Exchanged value "%s" in field "%s" with "%s"' % \
                  (old_field_val, field_name, dup_field_val) + \
                  ' from misspellings dictionary'

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Randomly exchange of a field value with another value
        #
        elif (mod_op == 'val_swap_prob') and (old_field_val != None):

          if (field_dict['type'] == 'freq'):  # A frequency file based field
            rand_num = random.randint(0, freq_files_length[field_name]-1)
            dup_field_val = freq_files[field_name][rand_num]

          elif (field_dict['type'] == 'date'):  # A date field
            rand_num = random.randint(field_dict['start_epoch'], \
                                      field_dict['end_epoch']-1)
            rand_date = epoch_to_date(rand_num)
            dup_field_val = rand_date[2]+rand_date[1]+rand_date[0]

          elif (field_dict['type'] == 'phone'):  # A phone number field
            area_code = random.choice(field_dict['area_codes'])
            max_digit = int('9'*field_dict['num_digits'])
            min_digit = int('1'*(int(1+round(field_dict['num_digits']/2.))))
            rand_num = random.randint(min_digit, max_digit)
            dup_field_val = area_code+' '+ \
                            str(rand_num).zfill(field_dict['num_digits'])

          elif (field_dict['type'] == 'ident'):  # A identification numb. field
            rand_num = random.randint(field_dict['start_id'], \
                                      field_dict['end_id']-1)
            dup_field_val = str(rand_num)

          if (dup_field_val != old_field_val):

            if (VERBOSE_OUTPUT == True):
              print '      Exchanged value in field "%s": "%s" -> "%s"' % \
                         (field_name, old_field_val, dup_field_val)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Randomly set to missing value
        #
        elif (mod_op == 'miss_prob') and (old_field_val != None):

          dup_field_val = missing_value  # Set to a missing value

          if (VERBOSE_OUTPUT == True):
            print '      Set field "%s" to missing value: "%s" -> "%s"' % \
                      (field_name, old_field_val, dup_field_val)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Randomly swap two words if the value contains at least two words
        #
        elif (mod_op == 'wrd_swap_prob') and (old_field_val != None) and \
             (' ' in old_field_val):

          # Count number of words
          #
          word_list = old_field_val.split(' ')
          num_words = len(word_list)

          if (num_words == 2):  # If only 2 words given
            swap_index = 0
          else:  # If more words given select position randomly
            swap_index = random.randint(0, num_words-2)

          tmp_word =                word_list[swap_index]
          word_list[swap_index] =   word_list[swap_index+1]
          word_list[swap_index+1] = tmp_word

          dup_field_val = ' '.join(word_list)

          if (dup_field_val != old_field_val):

            if (VERBOSE_OUTPUT == True):
              print '      Swapped words in field "%s": "%s" -> "%s"' % \
                    (field_name, old_field_val, dup_field_val)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Randomly create a new value if the field value is empty (missing)
        #
        elif (mod_op == 'new_val_prob') and (old_field_val == None):

          if (field_dict['type'] == 'freq'):  # A frequency file based field
            rand_num = random.randint(0, freq_files_length[field_name]-1)
            dup_field_val = freq_files[field_name][rand_num]

          elif (field_dict['type'] == 'date'):  # A date field
            rand_num = random.randint(field_dict['start_epoch'], \
                                      field_dict['end_epoch']-1)
            rand_date = epoch_to_date(rand_num)
            dup_field_val = rand_date[2]+rand_date[1]+rand_date[0]

          elif (field_dict['type'] == 'phone'):  # A phone number field
            area_code = random.choice(field_dict['area_codes'])
            max_digit = int('9'*field_dict['num_digits'])
            min_digit = int('1'*(int(1+round(field_dict['num_digits']/2.))))
            rand_num = random.randint(min_digit, max_digit)
            dup_field_val = area_code+' '+ \
                            str(rand_num).zfill(field_dict['num_digits'])

          elif (field_dict['type'] == 'ident'):  # A identification number
            rand_num = random.randint(field_dict['start_id'], \
                                      field_dict['end_id']-1)
            dup_field_val = str(rand_num)

          if (VERBOSE_OUTPUT == True):
            print '      Exchanged missing value "%s" in field "%s" with "%s"'\
                  % (missing_value, field_name, dup_field_val)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Random substitution of a character
        #
        elif (mod_op == 'sub_prob') and (old_field_val != None):

          # Get an substitution position randomly
          #
          rand_sub_pos = error_position(dup_field_val, 0)

          if (rand_sub_pos != None):  # If a valid position was returned

            old_char = dup_field_val[rand_sub_pos]
            new_char = error_character(old_char, field_dict['char_range'])

            new_field_val = dup_field_val[:rand_sub_pos] + new_char + \
                            dup_field_val[rand_sub_pos+1:]

            if (new_field_val != dup_field_val):
              dup_field_val = new_field_val

              if (VERBOSE_OUTPUT == True):
                print '      Substituted character "%s" with "%s" in field ' \
                      % (old_char, new_char) + '"%s": "%s" -> "%s"' % \
                      (field_name, old_field_val, dup_field_val)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Random insertion of a character
        #
        elif (mod_op == 'ins_prob') and (old_field_val != None):

          # Get an insert position randomly
          #
          rand_ins_pos = error_position(dup_field_val, +1)
          rand_char =    random.choice(field_range)

          if (rand_ins_pos != None):  # If a valid position was returned
            dup_field_val = dup_field_val[:rand_ins_pos] + rand_char + \
                            dup_field_val[rand_ins_pos:]

            if (VERBOSE_OUTPUT == True):
              print '      Inserted char "%s" into field "%s": "%s" -> "%s"' \
                    % (rand_char, field_name, old_field_val, dup_field_val)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Random deletion of a character
        #
        elif (mod_op == 'del_prob') and (old_field_val != None) and \
             (len(old_field_val) > 1):  # Field must have at least 2 characters

          # Get a delete position randomly
          #
          rand_del_pos = error_position(dup_field_val, 0)

          del_char = dup_field_val[rand_del_pos]

          dup_field_val = dup_field_val[:rand_del_pos] + \
                          dup_field_val[rand_del_pos+1:]

          if (VERBOSE_OUTPUT == True):
            print '      Deleted character "%s" in field "%s": "%s" -> "%s"' \
                  % (del_char, field_name, old_field_val, dup_field_val)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Random transposition of two characters
        #
        elif (mod_op == 'trans_prob') and (old_field_val != None) and \
             (len(dup_field_val) > 1):  # Field must have at least 2 characters

          # Get a transposition position randomly
          #
          rand_trans_pos = error_position(dup_field_val, -1)

          trans_chars = dup_field_val[rand_trans_pos:rand_trans_pos+2]
          trans_chars2 = trans_chars[1] + trans_chars[0]  # Do transposition

          new_field_val = dup_field_val[:rand_trans_pos] + trans_chars2 + \
                          dup_field_val[rand_trans_pos+2:]

          if (new_field_val != dup_field_val):
            dup_field_val = new_field_val

            if (VERBOSE_OUTPUT == True):
              print '      Transposed characters "%s" in field "%s": "%s" ' % \
                    (trans_chars, field_name, old_field_val) + \
                    '-> "%s"' % (dup_field_val)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Random insertion of a space (thus splitting a word)
        #
        elif (mod_op == 'spc_ins_prob') and (old_field_val != None) and \
             (len(dup_field_val) > 1):  # Field must have at least 2 characters

          # Randomly select the place where to insert a space (make sure no
          # spaces are next to this place)
          #
          rand_ins_pos = error_position(dup_field_val, 0)
          while (dup_field_val[rand_ins_pos-1] == ' ') or \
                (dup_field_val[rand_ins_pos] == ' '):
            rand_ins_pos = error_position(dup_field_val, 0)

          new_field_val = dup_field_val[:rand_ins_pos] + ' ' + \
                          dup_field_val[rand_ins_pos:]

          if (new_field_val != dup_field_val):
            dup_field_val = new_field_val

            if (VERBOSE_OUTPUT == True):
              print '      Inserted space " " into field "%s": "%s" -> "%s"' \
                    % (field_name, old_field_val, dup_field_val)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Random deletion of a space (thus merging two words)
        #
        elif (mod_op == 'spc_del_prob') and (old_field_val != None) and \
             (' ' in dup_field_val):  # Field must contain a space character

          # Count number of spaces and randomly select one to be deleted
          #
          num_spaces = dup_field_val.count(' ')

          if (num_spaces == 1):
            space_ind = dup_field_val.index(' ')  # Get index of the space
          else:
            rand_space = random.randint(1, num_spaces-1)
            space_ind = dup_field_val.index(' ', 0)  # Get index of first space
            for i in range(rand_space):
              # Get index of following spaces
              space_ind = dup_field_val.index(' ', space_ind)

          new_field_val = dup_field_val[:space_ind] + \
                          dup_field_val[space_ind+1:]

          if (new_field_val != dup_field_val):
            dup_field_val = new_field_val

            if (VERBOSE_OUTPUT == True):
              print '      Deleted space " " from field "%s": "%s" -> "%s"' % \
                    (field_name, old_field_val, dup_field_val)

        # Now check if the modified field value is different - - - - - - - - -
        #
        if (old_field_val == org_field_val) and \
           (dup_field_val != old_field_val):  # The first field modification
          field_mod_count_dict[field_name] = 1
          num_modif_in_record += 1

        elif (old_field_val != org_field_val) and \
             (dup_field_val != old_field_val):  # Following field modifications
          field_mod_count_dict[field_name] += 1
          num_modif_in_record += 1

        if (dup_field_val != old_field_val):
          dup_rec_dict[field_name] = dup_field_val

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Now check if the duplicate record differs from the original
    #
    rec_data = dup_rec_dict.copy()  # Make a copy of the record dictionary
    del(rec_data['rec_id'])  # Remove the record identifier
    rec_list = rec_data.items()
    rec_list.sort()
    rec_str = str(rec_list)

    if (rec_str not in all_rec_set):  # Check if same record already created
      all_rec_set.add(rec_str)
      org_rec_used[org_rec_id] = 1

      dup_rec[dup_rec_id] = dup_rec_dict  # Insert into duplicate records
      rec_cnt += 1

      d += 1  # Duplicate counter (loop counter)

      # Print original and duplicate records field by field - - - - - - - - - -
      #
      if (VERBOSE_OUTPUT == True):
        print '  Original and duplicate records:'
        print '    Number of modifications in record: %d' % \
              (num_modif_in_record)
        print '    Record ID         : %-30s | %-30s' % \
              (org_rec_dict['rec_id'], dup_rec_dict['rec_id'])
        for field_name in field_names:
          print '    %-18s: %-30s | %-30s' % \
                (field_name, org_rec_dict.get(field_name, missing_value), \
                 dup_rec_dict.get(field_name, missing_value))
        print

    else:
      if (VERBOSE_OUTPUT == True):
        print '  No random modifications for record "%s" -> Choose another' % \
              (dup_rec_id)

    if (VERBOSE_OUTPUT == True):
      print

# -----------------------------------------------------------------------------
# Write output csv file
#
print
print 'Step 3: Write output file'

all_rec = org_rec  # Merge original and duplicate records
all_rec.update(dup_rec)

# Get all record IDs and shuffle them randomly
#
all_rec_ids = all_rec.keys()
random.shuffle(all_rec_ids)

# Make a list of field names and sort them according to column number
#

field_name_list = ['rec_id']+field_names

# Open output file
#
try:
  f_out = open(output_file,'w')
except:
  print 'Error: Can not write to output file "%s"' % (output_file)
  sys.exit()

# Write header line - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
if (save_header == True):
  header_line = ''
  for field_name in field_name_list:
    header_line = header_line + field_name+ ', '
  header_line = header_line[:-2]
  f_out.write(header_line+os.linesep)

# Loop over all record IDs
#
for rec_id in all_rec_ids:
  rec_dict = all_rec[rec_id]
  out_line = ''
  for field_name in field_name_list:
    out_line = out_line + rec_dict.get(field_name, missing_value) + ', '
    if (field_name == 'rec_id') and (out_line[-6:] == '-org, '):
      out_line += '  '

  # Remove last comma and space and add line separator
  #
  out_line = out_line[:-2]
  # print out_line
  f_out.write(out_line+os.linesep)

f_out.close()

print 'End.'

# =============================================================================
