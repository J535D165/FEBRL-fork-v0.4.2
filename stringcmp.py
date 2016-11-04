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
# The Original Software is: "stringcmp.py"
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

"""Module with various approximate string comparison methods.

Provides routines for various approximate string comparisons. All return a
similarity value between 0.0 (strings are totally different) to 1.0 (strings
are the same).

Comparison methods provided:

  exact          Exact comparison
  jaro           Jaro
  winkler        Winkler (based on Jaro)  (for backwards compatibility)
  qgram          q-gram based
  bigram         2-gram based  (for backwards compatibility)
  posqgram       Positional q-gram based
  sgram          Skip-gram based
  editdist       Edit-distance (or Levenshtein distance)
  mod_editdist   Modified edit-distance (with transposition cost 1, not 2)
  bagdist        Bag distance (cheap distance based method)
  swdist         Smith-Waternam distance
  syllaligndist  Syllable alignment distance
  seqmatch       Uses Python's standard library 'difflib'
  compression    Based on Zlib compression algorithm
  lcs            (Repeated) longest common substring, improves results for
                 swapped words
  ontolcs        Ontology alignment string comparison based on longest common
                 substring, Hamacher product and Winkler heuristics.
  permwinkler    Winkler combined with permutations of words, improves results
                 for swapped words
  sortwinkler    Winkler with sorted words (if more than one), improves results
                 for swapped words
  editex         Phonetic aware edit-distance (Zobel et al. 1996)
  twoleveljaro   Apply Jaro comparator at word level, with words being compared
                 using a selectable approximate string comparison function
  charhistogram  Get histogram of characters for both strings and calculate the
                 cosine similarity between the two histogram vectors

See doc strings of individual functions for detailed documentation.

If called from command line, a test routine is run which prints example
approximate string comparisons for various string pairs.
"""

# =============================================================================
# Imports go here

import bz2
import difflib
import logging
import math
import time
import zlib

import encode  # For Phonix transformation routine (used in syllable alignment
               # distance)
import mymath  # Contains arithmetic coder

# =============================================================================
# Special character used in the Jaro, Winkler and q-gram comparions functions.
# Thanks to Luca Montecchiani (luca.mon@aliceposta.it).
#
JARO_MARKER_CHAR = chr(1)
QGRAM_START_CHAR = chr(1)
QGRAM_END_CHAR =   chr(2)

# =============================================================================

def do_stringcmp(cmp_method, str1, str2, min_threshold = None):
  """A 'chooser' functions which performs the selected comparison method.

  For each approximate string comparison method, various callable versions are
  provided that set their parameters to commonly used values.

  For each method, appending the string '-winkler' will result in the Winkler
  modification being applied (increasing the similarity weight if the beginning
  of the two strings are the same, up to first four characters).

  Possible values for 'cmp_method' are:

    exact            Exact comparison
    jaro             Jaro's method
    winkler          Jaro's method with Winkler modification (same as calling
                     'jaro-winkler')
    qgram1short      q-grams of length 1, divisor is shortest string length
    qgram1avrg       q-grams of length 1, divisor is average string length
    qgram1long       q-grams of length 1, divisor is longest string length
    qgram2short      q-grams of length 2, divisor is shortest string length
    qgram2avrg       q-grams of length 2, divisor is average string length
    qgram2long       q-grams of length 2, divisor is longest string length
    qgram3short      q-grams of length 3, divisor is shortest string length
    qgram3avrg       q-grams of length 3, divisor is average string length
    qgram3long       q-grams of length 3, divisor is longest string length
    qgram1Pshort     Padded q-grams of length 1, divisor is shortest string
                     length
    qgram1Pavrg      Padded q-grams of length 1, divisor is average string
                     length
    qgram1Plong      Padded q-grams of length 1, divisor is longest string
                     length
    qgram2Pshort     Padded q-grams of length 2, divisor is shortest string
                     length
    qgram2Pavrg      Padded q-grams of length 2, divisor is average string
                     length
    qgram2Plong      Padded q-grams of length 2, divisor is longest string
                     length
    qgram3Pshort     Padded q-grams of length 3, divisor is shortest string
                     length
    qgram3Pavrg      Padded q-grams of length 3, divisor is average string
                     length
    qgram3Plong      Padded q-grams of length 3, divisor is longest string
                     length
    posqgram1short   Positional q-grams of length 1, divisor is shortest length
    posqgram1avrg    Positional q-grams of length 1, divisor is average length
    posqgram1long    Positional q-grams of length 1, divisor is longest string
    posqgram2short   Positional q-grams of length 2, divisor is shortest length
    posqgram2avrg    Positional q-grams of length 2, divisor is average length
    posqgram2long    Positional q-grams of length 2, divisor is longest string
    posqgram3short   Positional q-grams of length 3, divisor is shortest length
    posqgram3avrg    Positional q-grams of length 3, divisor is average length
    posqgram3long    Positional q-grams of length 3, divisor is longest string
    posqgram1Pshort  Padded positional q-grams of length 1, divisor is shortest
                     string length
    posqgram1Pavrg   Padded positional q-grams of length 1, divisor is average
                     string length
    posqgram1Plong   Padded positional q-grams of length 1, divisor is longest
                     string length
    posqgram2Pshort  Padded positional q-grams of length 2, divisor is shortest
                     string length
    posqgram2Pavrg   Padded positional q-grams of length 2, divisor is average
                     string length
    posqgram2Plong   Padded positional q-grams of length 2, divisor is longest
                     string length
    posqgram3Pshort  Padded positional q-grams of length 3, divisor is shortest
                     string length
    posqgram3Pavrg   Padded positional q-grams of length 3, divisor is average
                     string length
    posqgram3lPong   Padded positional q-grams of length 3, divisor is longest
                     string length
    sgramshort       Skip-grams, divisor is shortest string length
    sgramavrg        Skip-grams, divisor is average string length
    sgramlong        Skip-grams, divisor is longest string length
    sgramPshort      Padded skip-grams, divisor is shortest string length
    sgramPavrg       Padded skip-grams, divisor is average string length
    sgramPlong       Padded skip-grams, divisor is longest string length
    editdist         Edit-distance (or Levenshtein distance)
    mod_editdist     Modified edit-distance (with transposition cost 1, not 2)
    editex           Phonetic aware edit-distance (Zobel et al. 1996)
    bagdist          Bag distance (cheap distance based method)
    swdistshort      Smith-Waterman distance, divisor is shortest length
    swdistavrg       Smith-Waterman distance, divisor is average length
    swdistlong       Smith-Waterman distance, divisor is longest length
    syllaldistshort  Syllable alignment distance, divisor is shortest length
    syllaldistavrg   Syllable alignment distance, divisor is average length
    syllaldistlong   Syllable alignment distance, divisor is longest length
    seqmatch         Uses Python's standard library 'difflib'
    compressZLib     Based on Zlib compression algorithm
    compressBZ2      Based on BZ2 compression algorithm
    compressArith    Based on arithmetic compression algorithm
    lcs2short        Longest common substring with minimum length of substrings
                     2, and divisor is shortest string length
    lcs2avrg         Longest common substring with minimum length of substrings
                     2, and divisor is average string length
    lcs2long         Longest common substring with minimum length of substrings
                     2, and divisor is longest string length
    lcs3short        Longest common substring with minimum length of substrings
                     3, and divisor is shortest string length
    lcs3avrg         Longest common substring with minimum length of substrings
                     3, and divisor is average string length
    lcs3long         Longest common substring with minimum length of substrings
                     3, and divisor is longest string length
    ontolcs2short    Ontology longest common substring with minimum length of
                     substrings 2, and divisor is shortest string length
    ontolcs2avrg     Ontology longest common substring with minimum length of
                     substrings 2, and divisor is average string length
    ontolcs2long     Ontology longest common substring with minimum length of
                     substrings 2, and divisor is longest string length
    ontolcs3short    Ontology longest common substring with minimum length of
                     substrings 3, and divisor is shortest string length
    ontolcs3avrg     Ontology longest common substring with minimum length of
                     substrings 3, and divisor is average string length
    ontolcs3long     Ontology longest common substring with minimum length of
                     substrings 3, and divisor is longest string length
    permwinkler      Winkler combined with permutations of words, improves
                     results for swapped words
    sortwinkler      Winkler with sorted words (if more than one), improves
                     results for swapped words

  This functions returns the similarity value (between 0.0 and 1.0) as well as
  the time needed to compare the strings (as floating-point value in seconds).
  """

  # Check if there is a 'divisor' value given (needed for several methods)
  #
  if ('short' in cmp_method):
    divisor = 'shortest'
  elif ('long' in cmp_method):
    divisor = 'longest'
  elif ('avrg' in cmp_method):
    divisor = 'average'
  else:
    divisor = None

  # For q- and s-gram based methods check for a padding 'P'
  #
  if (('gram' in cmp_method) and ('P' in cmp_method)):
    padded = True
  else:
    padded = False

  if (cmp_method.startswith('exa')):
    start_time = time.time()
    sim_weight = exact(str1, str2)
    time_used = time.time() - start_time

  elif (cmp_method.startswith('jaro')):
    start_time = time.time()
    sim_weight = jaro(str1, str2, min_threshold)
    time_used = time.time() - start_time
  elif (cmp_method.startswith('winkler')):
    start_time = time.time()
    sim_weight = winkler(str1, str2, min_threshold)
    time_used = time.time() - start_time

  elif (cmp_method.startswith('qgram')):
    q = int(cmp_method[5])  # Length of q-grams
    start_time = time.time()
    sim_weight = qgram(str1, str2, q, divisor, min_threshold, padded)
    time_used = time.time() - start_time
  elif (cmp_method.startswith('posqgram')):
    q = int(cmp_method[8])  # Length of q-grams
    max_dist = 2
    start_time = time.time()
    sim_weight = posqgram(str1, str2, q, max_dist, divisor, min_threshold,
                          padded)
    time_used = time.time() - start_time
  elif (cmp_method.startswith('sgram')):
    start_time = time.time()
    sim_weight = sgram(str1, str2, [[0],[0,1],[1,2]], divisor, min_threshold,
                       padded)
    time_used = time.time() - start_time

  elif (cmp_method.startswith('editdist')):
    start_time = time.time()
    sim_weight = editdist(str1, str2, min_threshold)
    time_used = time.time() - start_time
  elif (cmp_method.startswith('mod_editdist')):
    start_time = time.time()
    sim_weight = mod_editdist(str1, str2, min_threshold)
    time_used = time.time() - start_time
  elif (cmp_method.startswith('swdist')):
    start_time = time.time()
    sim_weight = swdist(str1, str2, divisor, min_threshold)
    time_used = time.time() - start_time
  elif (cmp_method.startswith('syllaldist')):
    start_time = time.time()
    sim_weight = syllaligndist(str1, str2, divisor, min_threshold, \
                               do_phonix=False)
    time_used = time.time() - start_time
  elif (cmp_method.startswith('bagdist')):
    start_time = time.time()
    sim_weight = bagdist(str1, str2, min_threshold)
    time_used = time.time() - start_time

  elif (cmp_method.startswith('seqmatch')):
    start_time = time.time()
    sim_weight = seqmatch(str1, str2, min_threshold)
    time_used = time.time() - start_time

  elif (cmp_method.startswith('compress')):
    if ('ZLib' in cmp_method):
      compr_method = 'zlib'
    elif ('BZ2' in cmp_method):
      compr_method = 'bz2'
    elif ('Arith' in cmp_method):
      compr_method = 'arith'
    else:
      logging.exception('Illegal compression method given: %s' % (cmp_method))
      raise Exception
    start_time = time.time()
    sim_weight = compression(str1, str2, compr_method, min_threshold)
    time_used = time.time() - start_time

  elif (cmp_method.startswith('lcs')):
    m = int(cmp_method[3])
    start_time = time.time()
    sim_weight = lcs(str1, str2, m, divisor, min_threshold)
    time_used = time.time() - start_time

  elif (cmp_method.startswith('ontolcs')):
    m = int(cmp_method[7])
    start_time = time.time()
    sim_weight = ontolcs(str1, str2, m, divisor, min_threshold)
    time_used = time.time() - start_time

  elif (cmp_method.startswith('sortwinkler')):
    start_time = time.time()
    sim_weight = sortwinkler(str1, str2, min_threshold)
    time_used = time.time() - start_time
  elif (cmp_method.startswith('permwinkler')):
    start_time = time.time()
    sim_weight = permwinkler(str1, str2, min_threshold)
    time_used = time.time() - start_time

  elif (cmp_method.startswith('editex')):
    start_time = time.time()
    sim_weight = editex(str1, str2, min_threshold)
    time_used = time.time() - start_time

  elif (not cmp_method.endswith('-winkler')):
    logging.exception('Illegal approximate string comparison method: %s' \
                      % (cmp_method))
    raise Exception

  # Check if Winkler modification should be applied - - - - - - - - - - - - - -
  #
  if (cmp_method.endswith('-winkler') == True) and \
     (sim_weight > 0.0) and (sim_weight < 1.0):
    sim_weight = winklermod(str1, str2, sim_weight)

  return sim_weight, time_used

# =============================================================================

def exact(str1, str2):
  """Do exact comparison of two strings.
  """

  if (str1 == '') or (str2 == ''):
    return 0.0
  elif (str1 == str2):
    return 1.0
  else:
    return 0.0

# =============================================================================

def jaro(str1, str2, min_threshold = None):
  """Return approximate string comparator measure (between 0.0 and 1.0)

  USAGE:
    score = jaro(str1, str2, min_threshold)

  ARGUMENTS:
    str1           The first string
    str2           The second string
    min_threshold  Minimum threshold between 0 and 1 (currently not used)

  DESCRIPTION:
    As desribed in 'An Application of the Fellegi-Sunter Model of
    Record Linkage to the 1990 U.S. Decennial Census' by William E. Winkler
    and Yves Thibaudeau.
  """

  # Quick check if the strings are empty or the same - - - - - - - - - - - - -
  #
  if (str1 == '') or (str2 == ''):
    return 0.0
  elif (str1 == str2):
    return 1.0

  len1 = len(str1)
  len2 = len(str2)

  halflen = max(len1,len2) / 2 - 1  # Or + 1?? PC 12/03/2009

  ass1 = ''  # Characters assigned in str1
  ass2 = ''  # Characters assigned in str2

  workstr1 = str1  # Copy of original string
  workstr2 = str2

  common1 = 0  # Number of common characters
  common2 = 0

  # Analyse the first string  - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  for i in range(len1):
    start = max(0,i-halflen)
    end   = min(i+halflen+1,len2)
    index = workstr2.find(str1[i],start,end)
    if (index > -1):  # Found common character
      common1 += 1
      ass1 = ass1 + str1[i]
      workstr2 = workstr2[:index]+JARO_MARKER_CHAR+workstr2[index+1:]

  # Analyse the second string - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  for i in range(len2):
    start = max(0,i-halflen)
    end   = min(i+halflen+1,len1)
    index = workstr1.find(str2[i],start,end)
    if (index > -1):  # Found common character
      common2 += 1
      ass2 = ass2 + str2[i]
      workstr1 = workstr1[:index]+JARO_MARKER_CHAR+workstr1[index+1:]

  if (common1 != common2):
    logging.error('Jaro: Wrong common values for strings "%s" and "%s"' % \
          (str1, str2) + ', common1: %i, common2: %i' % (common1, common2) + \
          ', common should be the same.')
    common1 = float(common1+common2) / 2.0  ##### This is just a fix #####

  if (common1 == 0):
    return 0.0

  # Compute number of transpositions  - - - - - - - - - - - - - - - - - - - - -
  #
  transposition = 0
  for i in range(len(ass1)):
    if (ass1[i] != ass2[i]):
      transposition += 1
  transposition = transposition / 2.0

  common1 = float(common1)
  w = 1./3.*(common1 / float(len1) + common1 / float(len2) + \
           (common1-transposition) / common1)

  assert (w >= 0.0) and (w <= 1.0), 'Similarity weight outside 0-1: %f' % (w)

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Jaro comparator string "%s" with "%s" value: %.3f' % \
                (str1, str2, w))
  return w

# =============================================================================

def winklermod(str1, str2, in_weight):
  """Applies the Winkler modification if beginning of strings is the same.

  USAGE:
    score = winklermod(str1, str2, in_weight)

  ARGUMENTS:
    str1       The first string
    str2       The second string
    in_weight  The basic similariy weight calculated by a string comparison
               method

  DESCRIPTION:
    As desribed in 'An Application of the Fellegi-Sunter Model of
    Record Linkage to the 1990 U.S. Decennial Census' by William E. Winkler
    and Yves Thibaudeau.

    If the begining of the two strings (up to fisrt four characters) are the
    same, the similarity weight will be increased.
 """

  # Quick check if the strings are empty or the same - - - - - - - - - - - - -
  #
  if (str1 == '') or (str2 == ''):
    return 0.0
  elif (str1 == str2):
    return 1.0

  # Compute how many characters are common at beginning - - - - - - - - - - - -
  #
  minlen = min(len(str1), len(str2))

  for same in range(1,minlen+1):
    if (str1[:same] != str2[:same]):
      break
  same -= 1
  if (same > 4):
    same = 4

  assert (same >= 0)

  winkler_weight = in_weight + same*0.1 * (1.0 - in_weight)

  assert (winkler_weight >= in_weight), 'Winkler modification is negative'

  assert (winkler_weight >= 0.0) and (winkler_weight <= 1.0), \
         'Similarity weight outside 0-1: %f' % (w)

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Winkler modification for string "%s" and "%s": Input ' % \
                (str1, str2)+'weight %.3f modified to %.3f' % \
                (in_weight, winkler_weight))

  return winkler_weight

# =============================================================================

def winkler(str1, str2, min_threshold = None):
  """For backwards compatibility, call Jaro followed by Winkler modification.
  """

  jaro_weight = jaro(str1, str2, min_threshold)

  return winklermod(str1, str2, jaro_weight)

# =============================================================================

def qgram(str1, str2, q=2, common_divisor = 'average', min_threshold = None,
          padded=True):
  """Return approximate string comparator measure (between 0.0 and 1.0)
     using q-grams (with default bigrams: q = 2).

  USAGE:
    score = qgram(str1, str2, q, common_divisor, min_threshold, padded)

  ARGUMENTS:
    str1            The first string
    str2            The second string
    q               The length of the q-grams to be used. Must be at least 1.
    common_divisor  Method of how to calculate the divisor, it can be set to
                    'average','shortest', or 'longest' , and is calculated
                    according to the lengths of the two input strings
    min_threshold   Minimum threshold between 0 and 1
    padded          If set to True (default), the beginnng and end of the
                    strings will be padded with (q-1) special characters, if
                    False no padding will be done.

  DESCRIPTION:
    q-grams are q-character sub-strings contained in a string. For example,
    'peter' contains the bigrams (q=2): ['pe','et','te','er'].

    Padding will result in specific q-grams at the beginning and end of a
    string, for example 'peter' converted into padded bigrams (q=2) will result
    in the following 2-gram list: ['*p','pe','et','te','er','r@'], with '*'
    illustrating the start and '@' the end character.

    This routine counts the number of common q-grams and divides by the
    average number of q-grams. The resulting number is returned.
  """

  if (q < 1):
    logging.exception('Illegal value for q: %d (must be at least 1)' % (q))
    raise Exception

  # Quick check if the strings are empty or the same - - - - - - - - - - - - -
  #
  if (str1 == '') or (str2 == ''):
    return 0.0
  elif (str1 == str2):
    return 1.0

  # Calculate number of q-grams in strings (plus start and end characters) - -
  #
  if (padded == True):
    num_qgram1 = len(str1)+q-1
    num_qgram2 = len(str2)+q-1
  else:
    num_qgram1 = max(len(str1)-(q-1),0)  # Make sure its not negative
    num_qgram2 = max(len(str2)-(q-1),0)

  # Check if there are q-grams at all from both strings - - - - - - - - - - - -
  # (no q-grams if length of a string is less than q)
  #
  if ((padded == False) and (min(num_qgram1, num_qgram2) == 0)):
    return 0.0

  # Calculate the divisor - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  if (common_divisor not in ['average','shortest','longest']):
    logging.exception('Illegal value for common divisor: %s' % \
                      (common_divisor))
    raise Exception

  if (common_divisor == 'average'):
    divisor = 0.5*(num_qgram1+num_qgram2)  # Compute average number of q-grams
  elif (common_divisor == 'shortest'):
    divisor = min(num_qgram1,num_qgram2)
  else:  # Longest
    divisor = max(num_qgram1,num_qgram2)

  # Use number of q-grams to quickly check for minimum threshold - - - - - - -
  #
  if (min_threshold != None):
    if (isinstance(min_threshold, float)) and (min_threshold > 0.0) and \
       (min_threshold > 0.0):

      max_common_qgram = min(num_qgram1,num_qgram2)

      w = float(max_common_qgram) / float(divisor)

      if (w  < min_threshold):
        return 0.0  # Similariy is smaller than minimum threshold

    else:
      logging.exception('Illegal value for minimum threshold (not between' + \
                        ' 0 and 1): %f' % (min_threshold))
      raise Exception

  # Add start and end characters (padding) - - - - - - - - - - - - - - - - - -
  #
  if (padded == True):
    qgram_str1 = (q-1)*QGRAM_START_CHAR+str1+(q-1)*QGRAM_END_CHAR
    qgram_str2 = (q-1)*QGRAM_START_CHAR+str2+(q-1)*QGRAM_END_CHAR
  else:
    qgram_str1 = str1
    qgram_str2 = str2

  # Make a list of q-grams for both strings - - - - - - - - - - - - - - - - - -
  #
  qgram_list1 = [qgram_str1[i:i+q] for i in range(len(qgram_str1) - (q-1))]
  qgram_list2 = [qgram_str2[i:i+q] for i in range(len(qgram_str2) - (q-1))]

  # Get common q-grams  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

  assert (w >= 0.0) and (w <= 1.0), 'Similarity weight outside 0-1: %f' % (w)

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('%d-gram comparator string "%s" with "%s" value: %.3f' % \
                (q, str1, str2, w))
  return w

# =============================================================================

def bigram(str1, str2, min_threshold = None):
  """For backwards compatibility.
  """

  return qgram(str1, str2, 2, 'average', min_threshold)

# =============================================================================

def posqgram(str1, str2, q=2, max_dist = 2, common_divisor = 'average',
             min_threshold = None, padded=True):
  """Return approximate string comparator measure (between 0.0 and 1.0)
     using positional q-grams (with default bigrams: q = 2).

  USAGE:
    score = posqgram(str1, str2, q, max_dist, common_divisor, min_threshold,
                     padded)

  ARGUMENTS:
    str1            The first string
    str2            The second string
    q               The length of the q-grams to be used. Must be at least 1.
    max_dist        Maximum distance allowed between two positional q-grams
                    (for example, with max_dist = 2 ('pe',6) and ('pe',8) are
                    considered to be similar, however, ('pe',1) and ('pe',7)
                    are not).
    common_divisor  Method of how to calculate the divisor, it can be set to
                    'average','shortest', or 'longest' , and is calculated
                    according to the lengths of the two input strings
    min_threshold   Minimum threshold between 0 and 1
    padded          If set to True (default), the beginnng and end of the
                    strings will be padded with (q-1) special characters, if
                    False no padding will be done.

  DESCRIPTION:
    q-grams are q-character sub-strings contained in a string. For example,
    'peter' contains the bigrams (q=2): ['pe','et','te','er'].

    Positional q-grams also contain the position within the string:
    [('pe',0),('et',1),('te',2),('er',3)].

    Padding will result in specific q-grams at the beginning and end of a
    string, for example 'peter' converted into padded bigrams (q=2) will result
    in the following 2-gram list:
    [('*p',0),('pe',1),('et',2),('te',3),('er',4),('r@',5)], with '*'
    illustrating the start and '@' the end character.

    This routine counts the number of common q-grams within the maximum
    distance and divides by the average number of q-grams. The resulting number
    is returned.
  """

  if (q < 1):
    logging.exception('Illegal value for q: %d (must be at least 1)' % (q))
    raise Exception

  if (max_dist < 0):
    logging.exception('Illegal value for maximum distance:: %d (must be ' % \
                      (max_dist) + 'zero or positive)')
    raise Exception

  # Quick check if the strings are empty or the same - - - - - - - - - - - - -
  #
  if (str1 == '') or (str2 == ''):
    return 0.0
  elif (str1 == str2):
    return 1.0

  # Calculate number of q-grams in strings (plus start and end characters) - -
  #
  if (padded == True):
    num_qgram1 = len(str1)+q-1
    num_qgram2 = len(str2)+q-1
  else:
    num_qgram1 = max(len(str1)-(q-1),0)  # Make sure its not negative
    num_qgram2 = max(len(str2)-(q-1),0)

  # Check if there are q-grams at all from both strings - - - - - - - - - - - -
  # (no q-grams if length of a string is less than q)
  #
  if ((padded == False) and (min(num_qgram1, num_qgram2) == 0)):
    return 0.0

  # Calculate the divisor - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  if (common_divisor not in ['average','shortest','longest']):
    logging.exception('Illegal value for common divisor: %s' % \
                      (common_divisor))
    raise Exception

  if (common_divisor == 'average'):
    divisor = 0.5*(num_qgram1+num_qgram2)  # Compute average number of q-grams
  elif (common_divisor == 'shortest'):
    divisor = min(num_qgram1,num_qgram2)
  else:  # Longest
    divisor = max(num_qgram1,num_qgram2)

  # Use number of q-grams to quickly check for minimum threshold - - - - - - -
  #
  if (min_threshold != None):
    if (isinstance(min_threshold, float)) and (min_threshold > 0.0) and \
       (min_threshold > 0.0):

      max_common_qgram = min(num_qgram1,num_qgram2)

      w = float(max_common_qgram) / float(divisor)

      if (w  < min_threshold):
        return 0.0  # Similariy is smaller than minimum threshold

    else:
      logging.exception('Illegal value for minimum threshold (not between' + \
                        ' 0 and 1): %f' % (min_threshold))
      raise Exception

  # Add start and end characters (padding) - - - - - - - - - - - - - - - - - -
  #
  if (padded == True):
    qgram_str1 = (q-1)*QGRAM_START_CHAR+str1+(q-1)*QGRAM_END_CHAR
    qgram_str2 = (q-1)*QGRAM_START_CHAR+str2+(q-1)*QGRAM_END_CHAR
  else:
    qgram_str1 = str1
    qgram_str2 = str2

  # Make a list of q-grams for both strings - - - - - - - - - - - - - - - - - -
  #
  qgram_list1 = [(qgram_str1[i:i+q],i) for i in range(len(qgram_str1) - (q-1))]
  qgram_list2 = [(qgram_str2[i:i+q],i) for i in range(len(qgram_str2) - (q-1))]

  # Get common q-grams  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

    pos_range = range(max(pos-max_dist,0), pos+max_dist+1)

    for test_pos in pos_range:
      test_pos_q_gram = (q_gram,test_pos)
      if (test_pos_q_gram in long_qgram_list):
        common += 1
        long_qgram_list.remove(test_pos_q_gram)  # Remove the counted q-gram
        break

  w = float(common) / float(divisor)

  assert (w >= 0.0) and (w <= 1.0), 'Similarity weight outside 0-1: %f' % (w)

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Positional %d-gram (max distance=%d) comparator string ' % \
                (q,max_dist) + '"%s" with "%s" value: %.3f' % (str1, str2, w))
  return w

# =============================================================================

def sgram(str1, str2, gc, common_divisor = 'average', min_threshold = None,
          padded=True):
  """Return approximate string comparator measure (between 0.0 and 1.0)
     using s-grams (skip-grams) with bigrams.

  USAGE:
    score = sgram(str1, str2, gc, common_divisor, min_threshold, padded)

  ARGUMENTS:
    str1            The first string
    str2            The second string
    gc              Gram class list (see below).
    common_divisor  Method of how to calculate the divisor, it can be set to
                    'average','shortest', or 'longest' , and is calculated
                    according to the lengths of the two input strings
    min_threshold   Minimum threshold between 0 and 1
    padded          If set to True (default), the beginnng and end of the
                    strings will be padded with (q-1) special characters, if
                    False no padding will be done.

  DESCRIPTION:
    Uses s-grams as described in:

    "Non-adjacent Digrams Improve Matching of Cross-Lingual Spelling Variants"
    by H. Keskustalo, A. Pirkola, K. Visala, E. Leppanen and J. Jarvelin,
    SPIRE 2003.

    Padding will result in special start and end characters being added at the
    beginning and the end of the character, similar to as done for the qgram
    and posqgram routines.
  """

  # Quick check if the strings are empty or the same - - - - - - - - - - - - -
  #
  if (str1 == '') or (str2 == ''):
    return 0.0
  elif (str1 == str2):
    return 1.0

  # Check if divisor is OK - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  if (common_divisor not in ['average','shortest','longest']):
    logging.exception('Illegal value for common divisor: %s' % \
                      (common_divisor))
    raise Exception

  # Extend strings with start and end characters
  #
  if (padded == True):
    tmp_str1 = QGRAM_START_CHAR+str1+QGRAM_END_CHAR
    tmp_str2 = QGRAM_START_CHAR+str2+QGRAM_END_CHAR
  else:
    tmp_str1 = str1
    tmp_str2 = str2

  len1 = len(tmp_str1)
  len2 = len(tmp_str2)

  common = 0.0   # Sum number of common s-grams over gram classes
  divisor = 0.0  # Sum of divisors over gram classes

  # Loop over all gram classes given - - - - - - - - - - - - - - - - - - - - -
  #
  for c in gc:

    sgram_list1 = []
    sgram_list2 = []

    for s in c:  # Skip distances

      for i in range(0,len1-s-1):
        sgram_list1.append(tmp_str1[i]+tmp_str1[i+s+1])
      for i in range(0,len2-s-1):
        sgram_list2.append(tmp_str2[i]+tmp_str2[i+s+1])

    num_sgram1 = len(sgram_list1)
    num_sgram2 = len(sgram_list2)

    if (common_divisor == 'average'):
      this_divisor = 0.5*(num_sgram1+num_sgram2)  # Average number of s-grams
    elif (common_divisor == 'shortest'):
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

  if (divisor == 0):  # One string did not have any s-grams
    w = 0.0
  else:
    w = common / divisor

  assert (w >= 0.0) and (w <= 1.0), 'Similarity weight outside 0-1: %f' % (w)

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('s-gram comparator string "%s" with "%s" value: %.3f' % \
                (str1, str2, w))
  return w

# =============================================================================

def editdist(str1, str2, min_threshold = None):
  """Return approximate string comparator measure (between 0.0 and 1.0)
     using the edit (or Levenshtein) distance.

  USAGE:
    score = editdist(str1, str2, min_threshold)

  ARGUMENTS:
    str1           The first string
    str2           The second string
    min_threshold  Minimum threshold between 0 and 1

  DESCRIPTION:
    The edit distance is the minimal number of insertions, deletions and
    substitutions needed to make two strings equal.

    For more information on the modified Soundex see:
    - http://www.nist.gov/dads/HTML/editdistance.html
  """

  # Quick check if the strings are empty or the same - - - - - - - - - - - - -
  #
  if (str1 == '') or (str2 == ''):
    return 0.0
  elif (str1 == str2):
    return 1.0

  n = len(str1)
  m = len(str2)
  max_len = max(n,m)

  if (min_threshold != None):
    if (isinstance(min_threshold, float)) and (min_threshold > 0.0) and \
       (min_threshold > 0.0):

      len_diff = abs(n-m)
      w = 1.0 - float(len_diff) / float(max_len)

      if (w  < min_threshold):
        return 0.0  # Similariy is smaller than minimum threshold

      else: # Calculate the maximum distance possible with this threshold
        max_dist = (1.0-min_threshold)*max_len

    else:
      logging.exception('Illegal value for minimum threshold (not between' + \
                        ' 0 and 1): %f' % (min_threshold))
      raise Exception

  if (n > m):  # Make sure n <= m, to use O(min(n,m)) space
    str1, str2 = str2, str1
    n, m =       m, n

  current = range(n+1)

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

    if (min_threshold != None) and (min(current) > max_dist):
      return 1.0 - float(max_dist+1) / float(max_len)

  w = 1.0 - float(current[n]) / float(max_len)

  assert (w >= 0.0) and (w <= 1.0), 'Similarity weight outside 0-1: %f' % (w)

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Edit-distance comparator string "%s" with "%s" value: %.3f' \
                % (str1, str2, w))
  return w

# =============================================================================

def mod_editdist(str1, str2, min_threshold = None):
  """Return approximate string comparator measure (between 0.0 and 1.0)
     using a modified edit (or Levenshtein) distance that counts transpositions
     as elementary operations as well. This is also called the Damerau-
     Levenshtein distance.

  USAGE:
    score = mod_editdist(str1, str2, min_threshold)

  ARGUMENTS:
    str1           The first string
    str2           The second string
    min_threshold  Minimum threshold between 0 and 1

  DESCRIPTION:
    The edit distance is the minimal number of insertions, deletions,
    substitutions and transpositions needed to make two strings equal.

    Compared to the original editdist function, which handles a transposition
    (like: 'sydney' <-> 'sydeny' as 2 operations (two substitutions or one
    insert and one delet), this modified version handles this as 1 operation.

    Based on code from Justin Zobel's 'vrank'.
  """

  # Quick check if the strings are empty or the same - - - - - - - - - - - - -
  #
  if (str1 == '') or (str2 == ''):
    return 0.0
  elif (str1 == str2):
    return 1.0

  n = len(str1)
  m = len(str2)
  max_len = max(n,m)

  if (min_threshold != None):
    if (isinstance(min_threshold, float)) and (min_threshold > 0.0) and \
       (min_threshold > 0.0):

      len_diff = abs(n-m)
      w = 1.0 - float(len_diff) / float(max_len)

      if (w  < min_threshold):
        return 0.0  # Similariy is smaller than minimum threshold

      else: # Calculate the maximum distance possible with this threshold
        max_dist = (1.0-min_threshold)*max_len

    else:
      logging.exception('Illegal value for minimum threshold (not between' + \
                        ' 0 and 1): %f' % (min_threshold))
      raise Exception

  if (n > m):  # Make sure n <= m, to use O(min(n,m)) space
    str1, str2 = str2, str1
    n, m =       m, n

  d = []  # Table with the full distance matrix

  current = range(n+1)
  d.append(current)

  for i in range(1,m+1):

    previous = current
    current =  [i]+n*[0]
    str2char = str2[i-1]

    for j in range(1,n+1):
      substitute = previous[j-1]
      if (str1[j-1] != str2char):
        substitute += 1

      if (i == 1) or (j == 1):  # First characters, no transposition possible

        # Get minimum of insert, delete and substitute
        #
        current[j] = min(previous[j]+1, current[j-1]+1, substitute)

      else:
        if (str1[j-2] == str2[i-1]) and (str1[j-1] == str2[i-2]):
          transpose = d[i-2][j-2] + 1
        else:
          transpose = d[i-2][j-2] + 3

        current[j] = min(previous[j]+1, current[j-1]+1, substitute, transpose)

    d.append(current)

    if (min_threshold != None) and (min(current) > max_dist):
      return 1.0 - float(max_dist+1) / float(max_len)

  w = 1.0 - float(current[n]) / float(max_len)

  assert (w >= 0.0) and (w <= 1.0), 'Similarity weight outside 0-1: %f' % (w)

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Modified edit-distance comparator string "%s" with "%s" ' % \
                (str1, str2) + 'value: %.3f' % (w))

  return w

# =============================================================================

def editdist_edits(str1, str2):
  """Return approximate string comparator measure (between 0.0 and 1.0)
     using the edit (or Levenshtein) distance as well as a triplet with the
     counts of the actual edits (inserts, deletes and substitutions).

  USAGE:
    score, edit_counts = editdist_edits(str1, str2)

  ARGUMENTS:
    str1           The first string
    str2           The second string

  DESCRIPTION:
    The edit distance is the minimal number of insertions, deletions and
    substitutions needed to make two strings equal.

    edit_counts  is a list with three elements that contain the number of
                 inserts, deletes and substitutions that were needed to convert
                 str1 into str2.

    For more information on the modified Soundex see:
    - http://www.nist.gov/dads/HTML/editdistance.html
  """

  # Check if the strings are empty or the same - - - - - - - - - - - - - - - -
  #
  if (str1 == '') and (str2 == ''):
    return 0.0, [0,0,0]

  elif (str1 == '') or (str2 == ''):
    if (str1 == ''):
      return 0.0, [len(str2),0,0]    # Inserts needed to get from empty to str1
    else:
      return 0.0, [0,len(str1),0,0]  # Deletes nedded to get from str2 to empty

  elif (str1 == str2):
    return 1.0, [0,0,0]

  n = len(str1)
  m = len(str2)

  d = []  # Table with the full distance matrix

  current = range(n+1)
  d.append(current)

  for i in range(1,m+1):

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

    d.append(current)

  # Count the number of edits that were needed - - - - - - - - - - - - - - - -
  #
  num_edits = [0,0,0]  # Number of Inserts, deletes and substitutions

  d_curr = d[m][n]  # Start with final position in table
  j = n
  i = m

  while (d_curr > 0):
    if (d[i-1][j-1]+1 == d_curr):  # Substitution
      i -= 1
      j -= 1
      num_edits[2] += 1
    elif (d[i-1][j]+1 == d_curr):  # Delete
      i -= 1
      num_edits[1] += 1
    elif (d[i][j-1]+1 == d_curr):  # Insert
      j -= 1
      num_edits[0] += 1

    else:  # Current position not larger than any of the previous positions
      if (d[i-1][j-1] == d_curr):
        i -= 1
        j -= 1
      elif (d[i-1][j] == d_curr):
        i -= 1
      elif (d[i][j-1] == d_curr):
        j -= 1
    d_curr = d[i][j]  # Update current position in table

  w = 1.0 - float(d[m][n]) / float(max(n,m))

  assert (w >= 0.0) and (w <= 1.0), 'Similarity weight outside 0-1: %f' % (w)

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Edit-distance comparator string "%s" with "%s" value: %.3f' \
                % (str1, str2, w) + ' and edit operation counts %s' % \
                (str(num_edits)))

  return w, num_edits

# =============================================================================

def bagdist(str1, str2, min_threshold = None):
  """Return approximate string comparator measure (between 0.0 and 1.0)
     using the bag distance.

  USAGE:
    score = bagdist(str1, str2, min_threshold)

  ARGUMENTS:
    str1           The first string
    str2           The second string
    min_threshold  Minimum threshold between 0 and 1 (currently not used)

  DESCRIPTION:
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

  # Quick check if the strings are empty or the same - - - - - - - - - - - - -
  #
  if (str1 == '') or (str2 == ''):
    return 0.0
  elif (str1 == str2):
    return 1.0

  n = len(str1)
  m = len(str2)

  list1 = list(str1)
  list2 = list(str2)

  for ch in str1:
    if (ch in list2):
      list2.remove(ch)

  for ch in str2:
    if (ch in list1):
      list1.remove(ch)

  b = max(len(list1),len(list2))

  w = 1.0 - float(b) / float(max(n,m))

  assert (w >= 0.0) and (w <= 1.0), 'Similarity weight outside 0-1: %f' % (w)

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Bag-distance comparator string "%s" with "%s" value: %.3f' % \
                (str1, str2, w))
  return w

# =============================================================================

def swdist(str1, str2, common_divisor = 'average', min_threshold = None):
  """Return approximate string comparator measure (between 0.0 and 1.0)
     using the Smith-Waterman distance.

  USAGE:
    score = swdist(str1, str2, common_divisor, min_threshold)

  ARGUMENTS:
    str1            The first string
    str2            The second string
    common_divisor  Method of how to calculate the divisor, it can be set to
                    'average','shortest', or 'longest' , and is calculated
                    according to the lengths of the two input strings
    min_threshold   Minimum threshold between 0 and 1

  DESCRIPTION:
    Smith-Waterman distance is commonly used in biological sequence alignment.

    Scores for matches, misses, gap and extension penalties are set to values
    described in:

    "The field matching problem: Algorithms and applications"
    by A.E. Monge and C.P. Elkan, 1996.
  """

  # Quick check if the strings are empty or the same - - - - - - - - - - - - -
  #
  if (str1 == '') or (str2 == ''):
    return 0.0
  elif (str1 == str2):
    return 1.0

  n = len(str1)
  m = len(str2)

  # Scores used for Smith-Waterman algorithm - - - - - - - - - - - - - - - - -
  #
  match_score =       5
  approx_score =      2
  mismatch_score =   -5
  gap_penalty =       5
  extension_penalty = 1

  # Calculate the divisor - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  if (common_divisor not in ['average','shortest','longest']):
    logging.exception('Illegal value for common divisor: %s' % \
                      (common_divisor))
    raise Exception

  if (common_divisor == 'average'):
    divisor = 0.5*(n+m)*match_score  # Average maximum score
  elif (common_divisor == 'shortest'):
    divisor = min(n,m)*match_score
  else:  # Longest
    divisor = max(n,m)*match_score

  # Dictionary with approximate match characters mapped into numbers
  # {a,e,i,o,u} -> 0, {d,t} -> 1, {g,j} -> 2, {l,r} -> 3, {m,n} -> 4,
  # {b,p,v} -> 5
  #
  approx_matches = {'a':0, 'b':5, 'd':1, 'e':0, 'g':2, 'i':0, 'j':2, 'l':3,
                    'm':4, 'n':4, 'o':0, 'p':5, 'r':3, 't':1, 'u':0, 'v':5}

  best_score = 0  # Keep the best score while calculating table

  d = []  # Table with the full distance matrix

  for i in range(n+1):  # Initalise table
    d.append([0.0]*(m+1))

  for i in range(1,n+1):
    for j in range(1,m+1):

      match = d[i-1][j-1]

      if (str1[i-1] == str2[j-1]):
        match += match_score
      else:
        approx_match1 = approx_matches.get(str1[i-1],-1)
        approx_match2 = approx_matches.get(str2[j-1],-1)

        if (approx_match1 >= 0) and (approx_match2 >= 0) and \
           (approx_match1 == approx_match2):
          match += approx_score
        else:
          match += mismatch_score

      insert = 0
      for k in range(1,i):
        score = d[i-k][j] - gap_penalty - k*extension_penalty
        insert = max(insert, score)

      delete = 0
      for l in range(1,j):
        score = d[i][j-l] - gap_penalty - l*extension_penalty
        delete = max(delete, score)

      d[i][j] = max(match, insert, delete, 0)
      best_score = max(d[i][j], best_score)

  # best_score can be min(len(str1),len)str2))*match_score (if one string is
  # a sub-string ofd the other string).
  #
  # The lower best_score the less similar the sequences are.
  #
  w = float(best_score) / float(divisor)

  assert (w >= 0.0) and (w <= 1.0), 'Similarity weight outside 0-1: %f' % (w)

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Smith-Waterman distance comparator string "%s" with "%s" ' % \
                (str1, str2) + 'value: %.3f' % (w))

  return w

# =============================================================================

def syllaligndist(str1, str2, common_divisor = 'average', min_threshold = None,
                  do_phonix=True):
  """Return approximate string comparator measure (between 0.0 and 1.0)
     using the syllable alignment distance.

  USAGE:
    score = syllaligndist(str1, str2, common_divisor, min_threshold, do_phonix)

  ARGUMENTS:
    str1            The first string
    str2            The second string
    common_divisor  Method of how to calculate the divisor, it can be set to
                    'average','shortest', or 'longest' , and is calculated
                    according to the lengths and number of syllables of the two
                    input strings
    min_threshold   Minimum threshold between 0 and 1
    do_phonix       A flag, if set to True the Phonix transformation will be
                    applied first to poth strings, otherwise the original
                    strings will be used.

  DESCRIPTION:
    The syllable alignment distance is based on syllables instead of characters
    and calculates a distance similar to edit distance.

    For more information see:
    "Syllable Alignment: A Novel Approach for Phonetic String Search"
    by Ruibin Gong and Tony k.Y. Chan, IEICE, 2006.
  """

  # Quick check if the strings are empty or the same - - - - - - - - - - - - -
  #
  if (str1 == '') or (str2 == ''):
    return 0.0
  elif (str1 == str2):
    return 1.0

  if (do_phonix == True):
    workstr1 = encode.phonix_transform(str1)
    workstr2 = encode.phonix_transform(str2)
  else:
    workstr1 = str1
    workstr2 = str2

  # Substitution and gap penalty weights
  #
  s1 = 1   # Aligning two characters (not syllable start) that are the same
  s2 = -1  # Aligning two characters (not syllable start) that are different
  s3 = -4  # Aligning a character with a syllable start
  s4 = 6   # Aligning two syllable starts that are the same
  s5 = -2  # Aligning two syllable starts that are different
  g1 = -1  # Aligning a gap with a character (not syllable start)
  g2 = -3  # Aligning a gap with a syllable start

  # Syllable scan, make beginning of each syllable an uppercase character - - -
  #
  syll_str_list = []  # List for the two syllable strings

  for s in (workstr1, workstr2):
    str_list = list(s)
    str_list[0] = str_list[0].upper()  # First char is start of first syllable
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

  n = len(wstr1)
  m = len(wstr2)

  # Calculate maximum number of syllable starts and other characters to get
  # maximum possible alignment weight
  #
  max_w1 = 0
  for c in wstr1:
    if c.isupper():
      max_w1 += s4  # Syllable start
    else:
      max_w1 += s1  # Other characters

  max_w2 = 0
  for c in wstr2:
    if c.isupper():
      max_w2 += s4  # Syllable start
    else:
      max_w2 += s1  # Other characters

  # Calculate the divisor - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  if (common_divisor not in ['average','shortest','longest']):
    logging.exception('Illegal value for common divisor: %s' % \
                      (common_divisor))
    raise Exception

  if (common_divisor == 'average'):
    divisor = 0.5*(max_w1+max_w2)  # Average weight
  elif (common_divisor == 'shortest'):
    divisor = min(max_w1,max_w2)
  else:  # Longest
    divisor = max(max_w1,max_w2)

  d = []  # Table with the full distance matrix

  for i in range(n+1):  # Initalise table
    d.append([0.0]*(m+1))

  for i in range(1,n+1):  # First column
    if (wstr1[i-1].isupper()):
      d[i][0] = d[i-1][0]+g2
    else:
      d[i][0] = d[i-1][0]+g1

  for j in range(1,m+1):  # First row
    if (wstr2[j-1].isupper()):
      d[0][j] = d[0][j-1]+g2
    else:
      d[0][j] = d[0][j-1]+g1

  for j in range(1,m+1):  # Fill in rest of table
    c2 = wstr2[j-1]

    for i in range(1,n+1):
      c1 = wstr1[i-1]

      if (c1.isupper()):
        x = d[i-1][j]+g2
      else:
        x = d[i-1][j]+g1

      if (c2.isupper()):
        y = d[i][j-1]+g2
      else:
        y = d[i][j-1]+g1

      if (c1.isupper() and c2.isupper()):
        if (c1 == c2):
          z = d[i-1][j-1]+s4
        else:
          z = d[i-1][j-1]+s5
      elif (c1.islower() and c2.islower()):
        if (c1 == c2):
          z = d[i-1][j-1]+s1
        else:
          z = d[i-1][j-1]+s2
      else:
        z = d[i-1][j-1]+s3

      d[i][j] = max(x,y,z)

  w = float(d[i][j]) / float(divisor)

  if (w < 0.0):
    w = 0.0

  assert (w >= 0.0) and (w <= 1.0), 'Similarity weight outside 0-1: %f' % (w)

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Syllable alignment distance comparator string ' + \
                '"%s" with "%s" value: %.3f' % (str1, str2, w))
  return w

# =============================================================================

def seqmatch(str1, str2, min_threshold = None):
  """Return approximate string comparator measure (between 0.0 and 1.0)
     using the Python standard library 'difflib' sequence matcher.

     Because the matches are not commutative, the pair and the swapped pair are
     compared and the average is taken.

  USAGE:
    score = seqmatch(str1, str2, min_threshold)

  ARGUMENTS:
    str1           The first string
    str2           The second string
    min_threshold  Minimum threshold between 0 and 1 (currently not used)

  DESCRIPTION:
    For more information on Python's 'difflib' library see:

      http://www.python.org/doc/current/lib/module-difflib.html
  """

  # Quick check if the strings are empty or the same - - - - - - - - - - - - -
  #
  if (str1 == '') or (str2 == ''):
    return 0.0
  elif (str1 == str2):
    return 1.0

  seq_matcher_1 = difflib.SequenceMatcher(None, str1, str2)
  seq_matcher_2 = difflib.SequenceMatcher(None, str2, str1)

  w = (seq_matcher_1.ratio()+seq_matcher_2.ratio()) / 2.0  # Return average

  assert (w >= 0.0) and (w <= 1.0), 'Similarity weight outside 0-1: %f' % (w)

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Seq-match comparator string "%s" with "%s" value: %.3f' % \
                (str1, str2, w))
  return w

# =============================================================================

def compression(str1, str2, compressor='zlib', min_threshold = None):
  """Return approximate string comparator measure (between 0.0 and 1.0)
     using the zlib compression library.

  USAGE:
    score = compression(str1, str2, compressor, min_threshold)

  ARGUMENTS:
    str1           The first string
    str2           The second string
    compressor     The compressor to be used, currently supported are:
                   'zlib' (default) using the Python standard libray zlib.py
                   compressor
                   'arith' using the arithmetic compressor implemented in the
                   mymath.py module.
                   'bz2' using the Python standard library bz2.py compressor
    min_threshold  Minimum threshold between 0 and 1 (currently not used)

  DESCRIPTION:
    For more information about using compression for similarity measures see:

    - Cilibrasi, R. and Vitanyi, P.: Clustering by compression. IEEE Trans.
      Infomat. Th. Submitted, 2004. See: http://arxiv.org/abs/cs.CV/0312044

    - Keogh, E., Lonardi, S. and Ratanamahatana, C.A.: Towards parameter-free
      data mining. Proceedings of the 2004 ACM SIGKDD international conference
      on Knowledge discovery and data mining, pp. 206-215, Seattle, 2004.

    - http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/306626
      for details about the arithmetic coder.
  """

  if (compressor not in ['zlib','arith','bz2']):
    logging.exception('Compression: Illegal compressor selected: "%s"' % \
                      (compressor)+' (must be either "zlib", "bz" or "arith")')
    raise Exception

  # Quick check if the strings are empty or the same - - - - - - - - - - - - -
  #
  if (str1 == '') or (str2 == ''):
    return 0.0
  elif (str1 == str2):
    return 1.0

  if (compressor == 'zlib'):
    c1 =  float(len(zlib.compress(str1)))
    c2 =  float(len(zlib.compress(str2)))
    c12 = 0.5 * (len(zlib.compress(str1+str2)) + len(zlib.compress(str2+str1)))

  elif (compressor == 'bz2'):
    c1 =  float(len(bz2.compress(str1)))
    c2 =  float(len(bz2.compress(str2)))
    c12 = 0.5 * (len(bz2.compress(str1+str2)) + len(bz2.compress(str2+str1)))

  else:  # Arithmetic compressor
    traintext = str1+str2  # All text as training to get probabilities
    trainprob = mymath.arith_coder_train(traintext)

    #print str1, str2
    #for (k,v) in trainprob.items():
    #  print k, v[0],'/',v[1]

    c1 =  float(mymath.arith_coder_encode(str1, trainprob))
    c2 =  float(mymath.arith_coder_encode(str2, trainprob))
    c12 = float(mymath.arith_coder_encode(str1+str2, trainprob))
    c21 = float(mymath.arith_coder_encode(str2+str1, trainprob))

    #print c1, c2, c12, 1.0 - (c12 - min(c1,c2)) / max(c1,c2)
    #print

    if (c21 != c12):
      print str1, str2, c12, c21

  if (c12 == 0.0):
    return 0.0  # Maximal distance

  w = 1.0 - (c12 - min(c1,c2)) / max(c1,c2)

  if (w < 0.0):
    print 'warning:Compression based comparison smaller than 0.0 with ' + \
          'strings "%s" and "%s": %.3f (cap to 1.0)' % (str1, str2, w)
    w = 0.0

  assert (w >= 0.0) and (w <= 1.0), 'Similarity weight outside 0-1: %f' % (w)

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Compression comparator string "%s" with "%s" value: %.3f' % \
                (str1, str2, w))
  return w

# =============================================================================

def lcs(str1, str2, min_common_len = 2, common_divisor = 'average',
        min_threshold = None):
  """Return approximate string comparator measure (between 0.0 and 1.0) using
     repeated longest common substring extractions.

  USAGE:
    score = lcs(str1, str2, min_common_len, common_divisor, min_threshold)

  ARGUMENTS:
    str1            The first string
    str2            The second string
    min_common_len  The minimum length of a common substring
    common_divisor  Method of how to calculate the divisor, it can be set to
                    'average','shortest', or 'longest' , and is calculated
                    according to the lengths of the two input strings
    min_threshold   Minimum threshold between 0 and 1

  DESCRIPTION:
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
  """

  if (min_common_len < 1):
    logging.exception('Minimum common length must be at least 1: %d' % \
                      (min_common_len))
    raise Exception

  # Quick check if the strings are empty or the same - - - - - - - - - - - - -
  #
  if (str1 == '') or (str2 == ''):
    return 0.0
  elif (str1 == str2):
    return 1.0

  len1 = len(str1)
  len2 = len(str2)

  # Calculate the divisor - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  if (common_divisor not in ['average','shortest','longest']):
    logging.exception('Illegal value for common divisor: %s' % \
                      (common_divisor))
    raise Exception

  if (common_divisor == 'average'):
    divisor = 0.5*(len1+len2)  # Compute average string length
  elif (common_divisor == 'shortest'):
    divisor = min(len1,len2)
  else:  # Longest
    divisor = max(len1,len2)

  # Use string length to quickly check for minimum threshold - - - - - - - - -
  #
  if (min_threshold != None):
    if (isinstance(min_threshold, float)) and (min_threshold > 0.0) and \
       (min_threshold < 1.0):

      max_common_len = min(len1,len2)

      w = float(max_common_len) / float(divisor)

      if (w  < min_threshold):
        return 0.0  # Similariy is smaller than minimum threshold

    else:
      logging.exception('Illegal value for minimum threshold (not between' + \
                        ' 0 and 1): %f' % (min_threshold))
      raise Exception

  w = 0.0

  for (s1,s2) in [(str1,str2),(str2,str1)]:
    #print '0:', s1, s2

    com_str, com_len, s1, s2 = do_lcs(s1, s2)  # Find initial LCS on input
    #print ' 1:',com_str, com_len, (s1, s2)

    total_com_str = com_str
    total_com_len = com_len

    while (com_len >= min_common_len): # As long as there are common substrings
      com_str, com_len, s1n, s2n = do_lcs(s1, s2)

      if (com_len >= min_common_len):
        #print ' 2:',com_str, com_len, (s1, s2)

        total_com_str += com_str
        total_com_len += com_len
        s1,s2 = s1n, s2n

    w += float(total_com_len) / float(divisor)

    #print '3:', s1, s2

  w /= 2.0

  assert (w >= 0.0) and (w <= 1.0), 'Similarity weight outside 0-1: %f' % (w)

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Longest common substring comparator string "%s" with "%s"' % \
                (str1, str2) + ' value: %.3f' % (w))

  return w

# -----------------------------------------------------------------------------

def do_lcs(str1, str2):
  """Subroutine to extract longest common substring from the two input strings.
     Returns the common substring, its length, and the two input strings with
     the common substring removed.
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

# =============================================================================

def ontolcs(str1, str2, min_common_len = 2, common_divisor = 'average',
            min_threshold = None):
  """Return approximate string comparator measure (between 0.0 and 1.0) using
     repeated longest common substring extractions, Hamacher difference and the
     Winkler heuristic.

  USAGE:
    score = ontolcs(str1, str2, min_common_len, common_divisor, min_threshold)

  ARGUMENTS:
    str1            The first string
    str2            The second string
    min_common_len  The minimum length of a common substring
    common_divisor  Method of how to calculate the divisor, it can be set to
                    'average','shortest', or 'longest' , and is calculated
                    according to the lengths of the two input strings
    min_threshold   Minimum threshold between 0 and 1

  DESCRIPTION:
    For more information about the ontology similarity measures see:

    - Giorgos Stoilos, Giorgos Stamou and Stefanos Kollinas:
      A String Metric for Ontology Alignment
      ISWC 2005, Springer LNCS 3729, pp 624-637, 2005.
  """

  P = 0.6 # Constant for Hamacher product difference, see above mentioned paper

  if (min_common_len < 1):
    logging.exception('Minimum common length must be at least 1: %d' % \
                      (min_common_len))
    raise Exception

  if (str1 == '') or (str2 == ''):
    return 0.0
  elif (str1 == str2):
    return 1.0

  len1 = len(str1)
  len2 = len(str2)

  # Calculate the divisor - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  if (common_divisor not in ['average','shortest','longest']):
    logging.exception('Illegal value for common divisor: %s' % \
                      (common_divisor))
    raise Exception

  if (common_divisor == 'average'):
    divisor = 0.5*(len1+len2)  # Compute average string length
  elif (common_divisor == 'shortest'):
    divisor = min(len1,len2)
  else:  # Longest
    divisor = max(len1,len2)

  w_lcs =  0.0  # Basic longest common sub-string weight
  h_diff = 0.0  # Hamacher product difference

  for (s1,s2) in [(str1,str2),(str2,str1)]:

    com_str, com_len, s1, s2 = do_lcs(s1, s2)  # Find initial LCS on input

    total_com_str = com_str
    total_com_len = com_len

    while (com_len >= min_common_len): # As long as there are common substrings
      com_str, com_len, s1n, s2n = do_lcs(s1, s2)

      if (com_len >= min_common_len):
        total_com_str += com_str
        total_com_len += com_len
        s1,s2 = s1n, s2n

    w_lcs += float(total_com_len) / float(divisor)

    # Calculate Hamacher product difference for sub-strings left
    #
    s1_len = float(len(s1)) / len1
    s2_len = float(len(s2)) / len2

    h_diff += s1_len*s2_len / (P + (1-P) * (s1_len + s2_len - s1_len*s2_len))

  w_lcs /=  2.0
  h_diff /= 2.0

  assert (w_lcs >= 0.0) and (w_lcs <= 1.0), \
         'Basic LCS similarity weight outside 0-1: %f' % (w_lcs)
  assert (h_diff >= 0.0) and (h_diff <= 1.0), \
         'Hamacher product difference outside 0-1: %f' % (h_diff)

  w_lcs_wink = winklermod(str1, str2, w_lcs)

  w = w_lcs_wink - h_diff  # A weight in interval [-1,1]

  w = w/2.0 + 0.5  # Scale into [0,1]

  assert (w >= 0.0) and (w <= 1.0), \
         'Ontology LCS similarity weight outside 0-1: %f' % (w)

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Ontology longest common substring comparator string ' + \
                '"%s" with "%s"' % (str1, str2) + ' value: %.3f' % (w))

  return w

# =============================================================================

def permwinkler(str1, str2, min_threshold = None):
  """Return approximate string comparator measure (between 0.0 and 1.0) using
     a combination of the Winkler string comparator on all permutations of
     words (ifd there are more than one in the input strings), which improves
     the results for swapped words.

  USAGE:
    score = permwinkler(str1, str2, min_threshold)

  ARGUMENTS:
    str1            The first string
    str2            The second string
    min_threshold   Minimum threshold between 0 and 1 (currently not used)

  DESCRIPTION:
    If one or both of the input strings contain more than one words all
    possible permutations of are compared using the Winkler approximate string
    comparator, and the maximum value is returned.

    If both input strings contain one word only then the standard Winkler
    string comparator is used.
  """

  # Quick check if the strings are empty or the same - - - - - - - - - - - - -
  #
  if (str1 == '') or (str2 == ''):
    return 0.0
  elif (str1 == str2):
    return 1.0

  if (' ' not in str1) and (' ' not in str2):
    w = winkler(str1, str2, min_threshold)  # Standard Winkler

  else:  # At least one of the strings contains two words

    str_list1 = str1.split(' ')
    str_list2 = str2.split(' ')

    perm_list1 = mymath.permute(str_list1)
    perm_list2 = mymath.permute(str_list2)

    w =        -1.0  # Maximal similarity measure
    max_perm = None

    for perm1 in perm_list1:
      for perm2 in perm_list2:

        # Calculate standard winkler for this permutation
        #
        this_w = winkler(perm1, perm2)

        if (this_w > w):
          w        = this_w
          max_perm = [perm1, perm2]

    logging.debug('Permutation Winkler best permutation: %s' % (str(max_perm)))

  assert (w >= 0.0) and (w <= 1.0), 'Similarity weight outside 0-1: %f' % (w)

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Permutation Winkler comparator string "%s" with "%s" value:' \
                % (str1, str2) + ' %.3f' % (w))
  return w

# =============================================================================

def sortwinkler(str1, str2, min_threshold = None):
  """Return approximate string comparator measure (between 0.0 and 1.0) using
     the Winkler string comparator on the word-sorted input strings (if there
     are more than one in the input strings), which improves the results for
     swapped words.

  USAGE:
    score = sortwinkler(str1, str2, min_threshold)

  ARGUMENTS:
    str1            The first string
    str2            The second string
    min_threshold   Minimum threshold between 0 and 1 (currently not used)

  DESCRIPTION:
    If one or both of the input strings contain more than one words then the
    input string is word-sorted before the standard Winkler approximate string
    comparator is applied.

    If both input strings contain one word only then the standard Winkler
    string comparator is used.
  """

  # Quick check if the strings are empty or the same - - - - - - - - - - - - -
  #
  if (str1 == '') or (str2 == ''):
    return 0.0
  elif (str1 == str2):
    return 1.0

  if (' ' in str1):  # Sort string 1
    word_list = str1.split(' ')
    word_list.sort()
    str1 = ' '.join(word_list)

  if (' ' in str2):  # Sort string 2
    word_list = str2.split(' ')
    word_list.sort()
    str2 = ' '.join(word_list)

  w = winkler(str1, str2)  # Standard Winkler

  assert (w >= 0.0) and (w <= 1.0), 'Similarity weight outside 0-1: %f' % (w)

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Sorted Winkler comparator string "%s" with "%s" value: %.3f' \
                % (str1, str2, w))
  return w

# =============================================================================

def editex(str1, str2, min_threshold = None):
  """Return approximate string comparator measure (between 0.0 and 1.0)
     using the editex distance.

  USAGE:
    score = editex(str1, str2, min_threshold)

  ARGUMENTS:
    str1           The first string
    str2           The second string
    min_threshold  Minimum threshold between 0 and 1

  DESCRIPTION:
    Based on ideas described in:

    "Phonetic String Matching: Lessons Learned from Information Retrieval"
    by Justin Zobel and Philip Dart, SIGIR 1995.

    Important: This function assumes that the input strings only contain
    letters and whitespace, but no other characters. A whitespace is handled
    like a slient sounds.
  """

  # Quick check if the strings are empty or the same - - - - - - - - - - - - -
  #
  if (str1 == '') or (str2 == ''):
    return 0.0
  elif (str1 == str2):
    return 1.0

  n = len(str1)
  m = len(str2)

  # Values for edit costs - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  BIG_COSTS = 3  # If characters are not in same group
  SML_COSTS = 2  # If characters are in same group

  # Mappings of letters into groups - - - - - - - - - - - - - - - - - - - - - -
  #
  groupsof_dict = {'a':0, 'b':1, 'c':2, 'd':3, 'e':0, 'f':1, 'g':2, 'h':7,
                   'i':0, 'j':2, 'k':2, 'l':4, 'm':5, 'n':5, 'o':0, 'p':1,
                   'q':2, 'r':6, 's':2, 't':3, 'u':0, 'v':1, 'w':7, 'x':2,
                   'y':0, 'z':2, '{':7}

  # Function to calculate cost of a deletion - - - - - - - - - - - - - - - - -
  #
  def delcost(char1, char2, groupsof_dict):

    if (char1 == char2):
      return 0

    code1 = groupsof_dict.get(char1,-1)  # -1 is not a char
    code2 = groupsof_dict.get(char2,-2)  # -2 if not a char

    if (code1 == code2) or (code2 == 7):  # Same or silent
      return SML_COSTS  # Small difference costs
    else:
      return BIG_COSTS

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if (' ' in str1):
    str1 = str1.replace(' ','{')
  if (' ' in str2):
    str2 = str2.replace(' ','{')

  if (n > m):  # Make sure n <= m, to use O(min(n,m)) space
    str1, str2 = str2, str1
    n, m =       m, n

  row = [0]*(m+1)  # Generate empty cost matrix
  F = []
  for i in range(n+1):
    F.append(row[:])

  F[1][0] = BIG_COSTS   # Initialise first row and first column of cost matrix
  F[0][1] = BIG_COSTS

  sum = BIG_COSTS
  for i in range(2,n+1):
    sum += delcost(str1[i-2], str1[i-1], groupsof_dict)
    F[i][0] = sum

  sum = BIG_COSTS
  for j in range(2,m+1):
    sum += delcost(str2[j-2], str2[j-1], groupsof_dict)
    F[0][j] = sum

  for i in range(1,n+1):

    if (i == 1):
      inc1 = BIG_COSTS
    else:
      inc1 = delcost(str1[i-2], str1[i-1], groupsof_dict)

    for j in range(1,m+1):
      if (j == 1):
        inc2 = BIG_COSTS
      else:
        inc2 = delcost(str2[j-2], str2[j-1], groupsof_dict)

      if (str1[i-1] == str2[j-1]):
        diag = 0
      else:
        code1 = groupsof_dict.get(str1[i-1],-1)  # -1 is not a char
        code2 = groupsof_dict.get(str2[j-1],-2)  # -2 if not a char

        if (code1 == code2):  # Same phonetic group
          diag = SML_COSTS
        else:
          diag = BIG_COSTS

      F[i][j] = min(F[i-1][j]+inc1, F[i][j-1]+inc2, F[i-1][j-1]+diag)

  w = 1.0 - float(F[n][m]) / float(max(F[0][m],F[n][0]))

  if (w < 0.0):
    w = 0.0

  assert (w >= 0.0) and (w <= 1.0), 'Similarity weight outside 0-1: %f' % (w)

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Editex comparator string "%s" with "%s" value: %.3f' \
                % (str1, str2, w))
  return w

# =============================================================================

def twoleveljaro(str1, str2, comp_funct = 'equal', min_threshold = None):
  """Return approximate string comparator measure (between 0.0 and 1.0)

  USAGE:
    score = jaro(str1, str2, comp_funct, min_threshold)

  ARGUMENTS:
    str1           The first string
    str2           The second string
    comp_funct     The function used to compare individual words. Either the
                   string 'equal' (default) or one of the string comparison
                   functions available in this module (i.e. a function which
                   takes two strings as input and returns a similarity value
                   between 0 and 1)
    min_threshold  Minimum threshold between 0 and 1 (currently not used)

  DESCRIPTION:
    This function applies Jaro comparator at word level, and additionally
    allows the comparison of individual words to be done using an approximate
    comparison function.

    If an approximate string comparison function is used for 'comp_funct' then
    the 'min_threshold' needs to be set as well in order to select the number
    of words that can match in the current window - otherwise the 'best' match
    will be selected, even if it has a very low similarity value.

    For a description of the Jaro string comparator see 'An Application of the
    Fellegi-Sunter Model of Record Linkage to the 1990 U.S. Decennial Census'
    by William E. Winkler and Yves Thibaudeau.
  """

  # Quick check if the strings are empty or the same - - - - - - - - - - - - -
  #
  if (str1 == '') or (str2 == ''):
    return 0.0
  elif (str1 == str2):
    return 1.0

  # If neither string contains a space (i.e. both are only one word) then use
  # the given word level comparison function
  #
  if (' ' not in str1) and (' ' not in str2):
    if (comp_funct == 'equal'):
      return 0.0 # Already tested if strings are the same, so here they are not

    # Calculate simple similarity value
    #
    w = comp_funct(str1, str2)

    assert (w >= 0.0) and (w <= 1.0), 'Similarity weight outside 0-1: %f' % (w)

    return w

  # If a comparison function is given, a minimum threshold is also required
  #
  if ((comp_funct != 'equal') and (min_threshold == None)):
    logging.exception('Comparison function is given but no minimal threshold')
    raise Exception

  # Convert strings into lists of words (whitespace separated)
  #
  list1 = str1.split()
  list2 = str2.split()

  len1 = len(list1)
  len2 = len(list2)

  halflen = max(len1,len2) / 2

  ass_list1 = []  # Words assigned in list1
  ass_list2 = []  # Words assigned in list2

  work_list1 = list1[:]  # Copy of original lists
  work_list2 = list2[:]

  common1 = 0  # Number of common characters
  common2 = 0

#  print halflen
#  print 'word lists:'
#  print ' ', list1
#  print ' ', list2
#  print

  # If 'equal' comparison function is given, then Jaro can be - - - - - - - - -
  # directly applied at word level
  #
  if (comp_funct == 'equal'):

#    print 'equal: analyse word list 1:', list1
    for i in range(len1):  # Analyse the first word list
#      print i,   worklist1, asslist1
#      print ' ', worklist2, asslist2

      start = max(0,i-halflen)
      end   = min(i+halflen+1,len2)
#      print start, end, list1[i], worklist2[start:end]
      if (list1[i] in work_list2[start:end]):  # Found common word
        ind = work_list2[start:end].index(list1[i])
        common1 += 1
        ass_list1.append(list1[i])
        work_list2[ind+start] = JARO_MARKER_CHAR
#      print

#    print
#    print 'equal: analyse word list 2:', list2
    for i in range(len2):  # Analyse the second string
#      print i,   worklist1, asslist1
#      print ' ', worklist2, asslist2

      start = max(0,i-halflen)
      end   = min(i+halflen+1,len1)
#      print start, end, list2[i], worklist1[start:end]
      if (list2[i] in work_list1[start:end]):  # Found common word
        ind = work_list1[start:end].index(list2[i])
        common2 += 1
        ass_list2.append(list2[i])
        work_list1[ind+start] = JARO_MARKER_CHAR
#      print

#    print 'common:', common1
#    print 'assigned:'
#    print ass_list1
#    print ass_list2

    if (common1 != common2):
      logging.error('Two-level-Jaro: Wrong common values for strings ' + \
                    '"%s" and "%s"' % (str1, str2) + \
                    ', common1: %i, common2: %i' % (common1, common2) + \
                    ', common should be the same.')
      common1 = float(common1+common2) / 2.0  ##### This is just a fix #####

  # For approximate comparison function, compare all words within current
  # 'window' and keep all matches above threshold, then select the best match
  #
  else:

#    print 'approx: analyse word list 1:', list1
    for i in range(len1):  # Analyse the first word list
#      print i,   work_list1, ass_list1
#      print ' ', work_list2, ass_list2
      start = max(0,i-halflen)
      end   = min(i+halflen+1,len2)
#      print start, end, list1[i], work_list2[start:end]
      search_word = list1[i]
      ind = -1  # The index of the best match found
      best_match_sim = -1
      word_ind = 0
      for word in work_list2[start:end]:
        tmp_sim = comp_funct(search_word, word)
        if (tmp_sim >= min_threshold):
          if (tmp_sim > best_match_sim):
            ind = word_ind
            best_match_sim = tmp_sim
        word_ind += 1
      if (ind >= 0):  # Found common word
#        print '  found match:', search_word, work_list2[ind+start], best_match_sim
        common1 += 1
        ass_list1.append(list1[i])
        work_list2[ind+start] = JARO_MARKER_CHAR
#        print '*', work_list2
#      print

#    print
#    print 'approx: analyse word list 2:', list2
    for i in range(len2):  # Analyse the second string
#      print i,   work_list1, ass_list1
#      print ' ', work_list2, ass_list2
      start = max(0,i-halflen)
      end   = min(i+halflen+1,len1)
#      print start, end, list2[i], work_list1[start:end]
      search_word = list2[i]
      ind = -1  # The index of the best match found
      best_match_sim = -1
      word_ind = 0
      for word in work_list1[start:end]:
        tmp_sim = comp_funct(search_word, word)
        if (tmp_sim >= min_threshold):
          if (tmp_sim > best_match_sim):
            ind = word_ind
            best_match_sim = tmp_sim
        word_ind += 1
      if (ind >= 0):  # Found common word
#        print '  found match:', search_word, work_list1[ind+start], best_match_sim
        common2 += 1
        ass_list2.append(list2[i])
        work_list1[ind+start] = JARO_MARKER_CHAR
#        print '*', work_list1
#      print

#    print 'common:', common1
#    print 'assigned:'
#    print ass_list1
#    print ass_list2

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
      logging.warning('Two-level-Jaro: Different common values for strings' + \
                      ' "%s" and "%s"' % (str1, str2) + \
                      ', common1: %i, common2: %i' % (common1, common2))
      common1 = float(common1+common2) / 2.0

  if (common1 == 0):
    return 0.0

  # Compute number of transpositions  - - - - - - - - - - - - - - - - - - - - -
  #
  min_num_ass_words = min(len(ass_list1), len(ass_list2))
  transposition = 0
  for i in range(min_num_ass_words):
    if (comp_funct == 'equal'):  # Standard way like done in Jaro comparator
      if (ass_list1[i] != ass_list2[i]):
        transposition += 1

    else:  # Again use approximate stringcomparison to calculate similarities
      tmp_sim = comp_funct(ass_list1[i], ass_list2[i])
      if (tmp_sim >= min_threshold):
#        print tmp_sim, ass_list1[i], ass_list2[i]
        transposition += 1

#  print 'transpositions:', transposition

  common1 = float(common1)
  w = 1./3.*(common1 / float(len1) + common1 / float(len2) + \
           (common1-transposition) / common1)

  assert (w >= 0.0) and (w <= 1.0), 'Similarity weight outside 0-1: %f' % (w)

  # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  logging.debug('Two-Level-Jaro comparator string "%s" with "%s" value: %.3f' \
                % (str1, str2, w))
  return w

# =============================================================================

def charhistogram(str1, str2, min_threshold = None):
  """Return approximate string comparator measure (between 0.0 and 1.0)

  USAGE:
    score = charhistogram(str1, str2, min_threshold)

  ARGUMENTS:
    str1           The first string
    str2           The second string
    min_threshold  Minimum threshold between 0 and 1 (currently not used)

  DESCRIPTION:
    This function counts all characters (and whitespaces) in the two strings
    and builds histrograms of characters. It then calculates the cosine
    similarity measure between these two histogram vectors.
  """

  # Quick check if the strings are empty or the same - - - - - - - - - - - - -
  #
  if (str1 == '') or (str2 == ''):
    return 0.0
  elif (str1 == str2):
    return 1.0

  histo1 = [0]*37
  histo2 = [0]*37

  workstr1 = str1.lower()
  workstr2 = str2.lower()

  for c in workstr1:
    if (c == ' '):
      histo1[0] += 1
    elif ((c >= 'a') and (c <= 'z')):  # Count characters
      histo1[ord(c)-96] += 1
    elif ((c >= '0') and (c <= '9')):  # Count digits
      histo1[ord(c)-21] += 1


  for c in workstr2:
    if (c == ' '):
      histo2[0] += 1
    elif ((c >= 'a') and (c <= 'z')):
      histo2[ord(c)-96] += 1
    elif ((c >= '0') and (c <= '9')):  # Count digits
      histo2[ord(c)-21] += 1

  #print histo1
  #print histo2

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

  return cos_sim

# =============================================================================
#
# Do some tests if called from command line
#
# Most test strings are taken from:
#   Approximate String Comparison and its Effect on an Advanced Record
#   Linkage System, Edward H. Porter and William W. Winkler, Bureau of
#   Census, 1997. Research report RR97/02.
#
if (__name__ == '__main__'):

  msg = []

  msg.append('Febrl module "stringcmp.py"')
  msg.append('---------------------------')
  msg.append('')

  strings = [['shackleford','dunningham','nichleson','jones','massey', \
              'abroms','hardin','itman','jeraldine','marhta','michelle', \
              'julies','tanya','dwayne','sean','jon','jon','brookhaven', \
              'brook hallow','decatur','fitzrureiter','higbee','higbee', \
              'lacura','iowa','lst','peter','abcde','yz','cunningham', \
              'campell','galloway','frederick','michele','jesse', \
              'jonathon','julies','yvette','dickson','dixon','peter', \
              'gondiwindi','delfinni','ein','do','doe', \
              'louise marie', 'maria louisa', 'mighty joe', 'kim zhu', \
              'lim zhau kim','prap','papr','gail','vest','sydney',
              'tsetung'], \
             ['shackelford','cunnigham','nichulson','johnson','massie', \
              'abrams','martinez','smith','geraldine','martha','michael', \
              'julius','tonya','duane','susan','john','jan','brrokhaven', \
              'brook hllw','decatir','fitzenreiter','highee','higvee', \
              'locura','iona','ist','peter','fghij','abcdef', \
              'cunnigham','campbell','calloway','fredrick','michelle', \
              'jessie','jonathan','juluis','yevett','dixon','dickson', \
              'ole', 'gondiwindiro','delfini','eni','od','deo', \
              'marie louise', 'louisa marie', 'joe mighty', 'zhou kim', \
              'kim lim zhao','papr','prap','gayle','west','sydeny',
              'zedong']]

  msg.append('     String 1      String 2  Jaro  J-W  1gram 2gram 3gram' + \
             ' 1pqgr 2pqgr 3pqgr Sgram eDist meDis bDist Editx SeqMa ComBZ' + \
             ' ComZL ComAC LCS2  LCS3  OLCS2 OLCS3 P-Win S-Win SWDis SyADi' + \
             ' Histo 2LJaro 2LJaroA')

  for i in range(len(strings[0])):
    str1 = strings[0][i]
    str2 = strings[1][i]

    s = '%13s %13s' % (str1, str2)
    s += ' %.3f' % (jaro(str1,str2))
    s += ' %.3f' % (winkler(str1,str2))
    s += ' %.3f' % (qgram(str1,str2,1))
    s += ' %.3f' % (qgram(str1,str2,2))
    s += ' %.3f' % (qgram(str1,str2,3))
    s += ' %.3f' % (posqgram(str1,str2,1))
    s += ' %.3f' % (posqgram(str1,str2,2))
    s += ' %.3f' % (posqgram(str1,str2,3))
    s += ' %.3f' % (sgram(str1,str2,[[0],[0,1],[1,2]]))
    s += ' %.3f' % (editdist(str1,str2))
    s += ' %.3f' % (mod_editdist(str1,str2))
    s += ' %.3f' % (bagdist(str1,str2))
    s += ' %.3f' % (editex(str1,str2))
    s += ' %.3f' % (seqmatch(str1,str2))
    s += ' %.3f' % (compression(str1,str2,'bz2'))
    s += ' %.3f' % (compression(str1,str2,'zlib'))
    s += ' %.3f' % (compression(str1,str2,'arith'))
    s += ' %.3f' % (lcs(str1,str2,2))
    s += ' %.3f' % (lcs(str1,str2,3))
    s += ' %.3f' % (ontolcs(str1,str2,2))
    s += ' %.3f' % (ontolcs(str1,str2,3))
    s += ' %.3f' % (permwinkler(str1,str2))
    s += ' %.3f' % (sortwinkler(str1,str2))
    s += ' %.3f' % (swdist(str1,str2))
    s += ' %.3f' % (syllaligndist(str1,str2))
    s += ' %.3f' % (charhistogram(str1,str2))
    s += ' %.3f' % (twoleveljaro(str1,str2))
    s += ' %.3f' % (twoleveljaro(str1,str2, qgram, 0.8))
    msg.append(s)

    if (qgram(str1, str2, 2) != sgram(str1,str2,gc=[[0]])):
      msg.append('  Error: 2-gram != s-gram (with gc=[[0]])')

    if (editdist(str1, str2) > bagdist(str1,str2)):
      msg.append('  Error: BadD > EditD')

    if (lcs(str1,str2,1) < lcs(str1,str2,2)):
      msg.append('  Error: LCS1 < LCS2')

    if (lcs(str1,str2,2) < lcs(str1,str2,3)):
      msg.append('  Error: LCS2 < LCS3')

    if (editdist(str1, str2) > mod_editdist(str1,str2)):
      msg.append('  Error: EditD > Modified EditD')

  for m in msg:
    print m

# =============================================================================
