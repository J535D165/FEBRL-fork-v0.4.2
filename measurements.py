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
# The Original Software is: "measurements.py"
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

"""Module measurements.py - Functions to assess quality and complexity.

   Provides functions to calculate the following measures:

   - Pairs completeness
   - Pairs quality
   - Reduction ratio

   - Accuracy
   - Precision
   - Recall
   - F-measure

   For more details about these measures and a general discussion on measuring
   data linkage quality and complexity please refer to:

     "Quality and Complexity Measures for Data Linkage and Deduplication"
     Peter Christen and Karl Goiser
     Book chapter in "Quality Measures in Data Mining"
     F. Guillet and H. Hamilton (eds)

     Studies in Computational Intelligence, vol. 43, Springer, February 2007.

  The measurement routines assume that the weight vector dictionary contains
  two record pair identifiers as keys, and that these are also available in
  the two data sets.

  It is also assumed the record identifiers are made of the form 'X..XY' with
  X..X being the record identifier (one or more digits) and Y being one digit,
  with 0 for the original record, and 1..9 for up to 9 duplicates. For example
  '12345' is duplicate number '5' of record '1234', while '100010' is the
  original record number '10001'.

  TODO:
  - PC correct for both dedup and linkage??
  - possibly a classifier..? or a dict with tuples: (rec_id1, rec_id2, M/N/P)
    (at least for quality measures, RR and PC only need rec-id pairs
"""

# =============================================================================
# Import necessary modules (Febrl modules first, then Python standard modules)

import auxiliary

import logging

# =============================================================================

def pairs_completeness(weight_vec_dict, dataset1, dataset2, get_id_funct,
                       match_check_funct):
  """Pairs completeness is measured as

       pc = Nm / M

     with Nm (<= M) being the number of correctly classified truly matched
     record pairs in the blocked comparison space, and M the total number of
     true matches.

     If both data sets are the same a deduplication is assumed, otherwise a
     linkage.

     The arguments that have to be set when this method is called are:
       weight_vec_dict    A dictionary containing weight vectors.
       dataset1           The initialised first data set object.
       dataset2           The initialised second data set object.
       get_id_funct       This has to be a function (or method), assumed to
                          have argument a record (assumed to be a list fo field
                          values), and returns the record identifier from that
                          record.
       match_check_funct  This has to be a function (or method), assumed to
                          have as arguments the two record identifiers of a
                          record pair and its weight vector, and returns True
                          if the record pair is from a true match, or False
                          otherwise. Thus, 'match_check_funct' is of the form:

                            match_flag = match_check_funct(rec_id1, rec_id2,
                                                           weight_vec)
  """

  auxiliary.check_is_dictionary('weight_vec_dict', weight_vec_dict)
  auxiliary.check_is_not_none('dataset1', dataset1)
  auxiliary.check_is_not_none('dataset2', dataset2)
  auxiliary.check_is_function_or_method('get_id_funct', get_id_funct)
  auxiliary.check_is_function_or_method('match_check_funct', match_check_funct)

  # Check if a deduplication will be done or a linkage - - - - - - - - - - - -
  #
  if (dataset1 == dataset2):
    do_dedup = True
  else:
    do_dedup = False

  logging.info('')
  logging.info('Calculate pairs completeness:')
  logging.info('  Data set 1: %s (containing %d records)' % \
               (dataset1.description, dataset1.num_records))
  if (do_dedup == True):
    logging.info('  Data sets are the same: Deduplication')
  else:
    logging.info('  Data set 2: %s (containing %d records)' % \
                 (dataset2.description, dataset2.num_records))
    logging.info('  Data sets differ:       Linkage')
  logging.info('  Number of record pairs in weight vector dictionary: %d' % \
               (len(weight_vec_dict)))

  num_all_true_matches = 0  # Count the total number of all true matches

  # For a deduplication only process data set 1 - - - - - - - - - - - - - - - -
  #
  if (do_dedup == True):

    # Build a dictionary with entity identifiers as keys and a list of their
    # record identifier (rec_ident) as values
    #
    entity_ident_dict = {}

    for (rec_ident, rec) in dataset1.readall():
      ent_id = get_id_funct(rec)

      this_rec_list = entity_ident_dict.get(ent_id, [])
      this_rec_list.append(rec_ident)
      entity_ident_dict[ent_id] = this_rec_list

    logging.info('  Number of unique entity identifiers in data set 1: %d' % \
                 (len(entity_ident_dict)))

    for (ent_id, rec_list) in entity_ident_dict.iteritems():
      num_this_rec = len(rec_list)

      if (num_this_rec > 1):
        num_all_true_matches += num_this_rec*(num_this_rec-1)/2

    # More efficent version: Only count number of matches ber record don't
    # store them
    #
    entity_ident_dict2 = {}

    for (rec_ident, rec) in dataset1.readall():
      ent_id = get_id_funct(rec)
      ent_id_count = entity_ident_dict2.get(ent_id, 0) + 1
      entity_ident_dict2[ent_id] = ent_id_count

    assert sum(entity_ident_dict2.values()) == dataset1.num_records

    tm = 0  # Total number of true matches (without indexing)

    for (ent_id, ent_count) in entity_ident_dict2.iteritems():
      tm += ent_count*(ent_count-1)/2

    assert num_all_true_matches == tm

  else:  # For a linkage - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # Build two dictionaries with  entity identifiers as keys and a list of
    # their record identifier (rec_ident) as values
    #
    entity_ident_dict1 = {}
    entity_ident_dict2 = {}

    for (rec_ident, rec) in dataset1.readall():
      ent_id = get_id_funct(rec)

      this_rec_list = entity_ident_dict1.get(ent_id, [])
      this_rec_list.append(rec_ident)
      entity_ident_dict1[ent_id] = this_rec_list

    logging.info('  Number of unique entity identifiers in data set 1: %d' % \
                 (len(entity_ident_dict1)))

    for (rec_ident, rec) in dataset2.readall():
      ent_id = get_id_funct(rec)

      this_rec_list = entity_ident_dict2.get(ent_id, [])
      this_rec_list.append(rec_ident)
      entity_ident_dict2[ent_id] = this_rec_list

    logging.info('  Number of unique entity identifiers in data set 2: %d' % \
                 (len(entity_ident_dict2)))

    # Now calculate total true match number (loop over smaller dict)
    #
    if (len(entity_ident_dict1) < len(entity_ident_dict2)):
      for (ent_id1, rec_list1) in entity_ident_dict1.iteritems():

        if (ent_id1 in entity_ident_dict2):
          rec_list2 = entity_ident_dict2[ent_id1]

          num_all_true_matches += len(rec_list1) * len(rec_list2)
    else:
      for (ent_id2, rec_list2) in entity_ident_dict2.iteritems():

        if (ent_id2 in entity_ident_dict1):
          rec_list1 = entity_ident_dict1[ent_id2]

          num_all_true_matches += len(rec_list1) * len(rec_list2)

    # More efficent version: Only count number of matches ber record don't
    # store them
    #
    entity_ident_dict3 = {}
    entity_ident_dict4 = {}

    for (rec_ident, rec) in dataset1.readall():
      ent_id = get_id_funct(rec)
      ent_id_count = entity_ident_dict3.get(ent_id, 0) + 1
      entity_ident_dict3[ent_id] = ent_id_count

    for (rec_ident, rec) in dataset2.readall():
      ent_id = get_id_funct(rec)
      ent_id_count = entity_ident_dict4.get(ent_id, 0) + 1
      entity_ident_dict4[ent_id] = ent_id_count

    assert sum(entity_ident_dict3.values()) == dataset1.num_records
    assert sum(entity_ident_dict4.values()) == dataset2.num_records

    tm = 0  # Total number of true matches (without indexing)

    if (len(entity_ident_dict3) < len(entity_ident_dict4)):
      for (ent_id, ent_count) in entity_ident_dict3.iteritems():
        if ent_id in entity_ident_dict4:
          tm += ent_count*entity_ident_dict4[ent_id]
    else:
      for (ent_id, ent_count) in entity_ident_dict4.iteritems():
        if ent_id in entity_ident_dict3:
          tm += ent_count*entity_ident_dict3[ent_id]

    assert num_all_true_matches == tm

  logging.info('  Number of all true matches: %d' % (num_all_true_matches))

  # Get number of true matches in weight vector dictionary - - - - - - - - - -
  #
  num_true_matches =  0
  num_false_matches = 0

  for (rec_id_tuple, this_vec) in weight_vec_dict.iteritems():

    if (match_check_funct(rec_id_tuple[0], rec_id_tuple[1], this_vec) == True):
      num_true_matches += 1
    else:
      num_false_matches += 1

  assert len(weight_vec_dict) == num_true_matches+num_false_matches

  logging.info('  Number of true and false matches in weight vector ' + \
               'dictionary: %d / %d' % (num_true_matches,num_false_matches))

  if (num_all_true_matches > 0):

    pc = float(num_true_matches) / float(num_all_true_matches)

    logging.info('  Pairs completeness: %.4f%%' % (100.0*pc)) # As percentage

  else:

    pc = 0.0

    logging.info('  No true matches - cannot calculate pairs completeness')

  assert pc <= 1.0, pc

  return pc

# -----------------------------------------------------------------------------

def pairs_quality(weight_vec_dict, match_check_funct):
  """Pairs quality is the ratio of true matches divided by the total number of
     matches of the compared record pairs returned after blocking. It is
     measured as:

       pq = |TP| / all_matches

     with TP being the true positives, and all matches being the number of
     weight vectors given.

     The arguments that have to be set when this method is called are:
       weight_vec_dict    A dictionary containing weight vectors.
       match_check_funct  This has to be a function (or method), assumed to
                          have as arguments the two record identifiers of a
                          record pair and its weight vector, and returns True
                          if the record pair is from a true match, or False
                          otherwise. Thus, 'match_check_funct' is of the form:

                            match_flag = match_check_funct(rec_id1, rec_id2,
                                                           weight_vec)
  """

  auxiliary.check_is_dictionary('weight_vec_dict', weight_vec_dict)
  auxiliary.check_is_function_or_method('match_check_funct', match_check_funct)

  total_num_rec_pairs = len(weight_vec_dict)

  logging.info('')
  logging.info('Calculate pairs quality:')
  logging.info('  Number of record pairs in weight vector dictionary: %d' % \
               (total_num_rec_pairs))

  # Get number of true matches in weight vector dictionary - - - - - - - - - -
  #
  num_true_matches =  0
  num_false_matches = 0

  for (rec_id_tuple, this_vec) in weight_vec_dict.iteritems():

    if (match_check_funct(rec_id_tuple[0], rec_id_tuple[1], this_vec) == True):
      num_true_matches += 1
    else:
      num_false_matches += 1

  assert total_num_rec_pairs == (num_true_matches + num_false_matches)

  logging.info('  Number of true and false matches in weight vector ' + \
               'dictionary: %d / %d' % (num_true_matches, num_false_matches))

  pq = float(num_true_matches) / total_num_rec_pairs

  logging.info('  Pairs quality: %.4f%%' % (100.0*pq)) # As percentage

  assert pq <= 1.0

  return pq

# -----------------------------------------------------------------------------

def reduction_ratio(weight_vector_dict, dataset1, dataset2):
  """For a linkage, reduction ratio is measured as

       rr = 1 - Nb / (|A| x |B|)

     with Nb being the number of record pairs produced by a blocking algorithm
     (i.e. the number of record pairs not removed by blocking) and |A| and |B|
     the number of records in the first and second data sets, respectively.

     For a deduplication the reduction ratio is measured as:

       rr = 1 - Nb / ((|A| x (|A|-1))/2)

     The reduction ratio measures the relative reduction of the comparison
     space, without taking into account the quality of the reduction (i.e. how
     many true matches and how many false matches are removed by blocking).

     If both data sets are the same a deduplication is assumed, otherwise a
     linkage.
  """

  # Check if a deduplication will be done or a linkage - - - - - - - - - - - -
  #
  if (dataset1 == dataset2):
    do_deduplication = True
  else:
    do_deduplication = False

  logging.info('')
  logging.info('Calculate reduction ratio:')
  logging.info('  Data set 1: %s (containing %d records)' % \
               (dataset1.description, dataset1.num_records))
  if (do_deduplication == True):
    logging.info('  Data sets are the same: Deduplication')
    logging.info('  Total number of possible comparisons: %d' % \
                 (dataset1.num_records*(dataset1.num_records-1)/2))
  else:
    logging.info('  Data set 2: %s (containing %d records)' % \
                 (dataset2.description, dataset2.num_records))
    logging.info('  Data sets differ:       Linkage')
    logging.info('    Total number of possible comparisons: %d' % \
                 (dataset1.num_records*dataset2.num_records))
  logging.info('  Number of record pairs in weight vector dictionary: %d' % \
               (len(weight_vector_dict)))

  Nb = len(weight_vector_dict)
  A = dataset1.num_records

  if (do_deduplication == True):

    rr = 1 - float(Nb) / float(0.5*A*(A-1))

  else:  # A linkage

    B = dataset2.num_records

    rr = 1 - float(Nb) / (float(A)*float(B))

  logging.info('  Reduction ratio: %.4f%%' % (100.0*rr))  # Log as percentage

  assert rr <= 1.0

  return rr

# =============================================================================

def quality_measures(weight_vec_dict, match_set, non_match_set,
                     match_check_funct):
  """Calculate several quality measures based on the number of true positives,
     true negatives, false positives and false negatives in the given match
     and non-match sets and weight vector dictionary using the given match
     check function.

     The function calculates and returns:

     - Accuracy:       (|TP|+|TN|)
                  ---------------------
                  (|TP|+|TN|+|FP|+|FN|)

     - Precision:    |TP|
                  -----------
                  (|TP|+|FP|)

     - Recall:       |TP|
                  -----------
                  (|TP|+|FN|)

     - F-Measure:   2 * (Precision * Recall)
                  --------------------------
                     (Precision + Recall)

     With TP the True Positives, TN the True negatives, FP the False Positives
     and FN the False Negatives.

     For a discussion about measuring data linkage and deduplication quality
     please refer to:

       Quality and Complexity Measures for Data Linkage and Deduplication
       Peter Christen and Karl Goiser

       Book chapter in "Quality Measures in Data Mining"
                       Studies in Computational Intelligence, Vol. 43
                       F. Guillet and H. Hamilton (eds), Springer
                       March 2007.
  """

  auxiliary.check_is_dictionary('weight_vec_dict', weight_vec_dict)
  auxiliary.check_is_set('match set', match_set)
  auxiliary.check_is_set('non match set', non_match_set)
  auxiliary.check_is_function_or_method('match_check_funct', match_check_funct)

  if ((len(match_set) + len(non_match_set)) != len(weight_vec_dict)):
    logging.exception('Match and non-match set are not of same length as ' + \
                      'weight vector dictionary: %d, %d / %d' % \
                      (len(match_set),len(non_match_set),len(weight_vec_dict)))
    raise Exception

  tp = 0.0
  fp = 0.0
  tn = 0.0
  fn = 0.0

  for rec_id_tuple in match_set:
    w_vec = weight_vec_dict[rec_id_tuple]

    if (match_check_funct(rec_id_tuple[0], rec_id_tuple[1], w_vec) == True):
      tp += 1
    else:
      fp += 1

  for rec_id_tuple in non_match_set:
    w_vec = weight_vec_dict[rec_id_tuple]

    if (match_check_funct(rec_id_tuple[0], rec_id_tuple[1], w_vec) == False):
      tn += 1
    else:
      fn += 1

  logging.info('')
  logging.info('Classification results: TP=%d, FP=%d / TN=%d, FN=%d' % \
               (tp, fp, tn, fn))

  if ((tp != 0) or (fp != 0) or (tn != 0) or (fn != 0)):
    acc = (tp + tn) / (tp + fp + tn + fn)
  else:
    acc = 0.0

  if ((tp != 0) or (fp != 0)):
    prec = tp / (tp + fp)
  else:
    prec = 0.0

  if ((tp != 0) or (fn != 0)):
    reca = tp / (tp + fn)
  else:
    reca = 0.0
  if ((prec != 0.0) or (reca != 0.0)):
    fmeas = 2*(prec*reca) / (prec+reca)
  else:
    fmeas = 0.0

  logging.info('Quality measures:')
  logging.info('  Accuracy: %.6f  Precision:%.4f  Recall: %.4f  ' % \
               (acc, prec, reca)+'F-measure: %.4f' % (fmeas))

  return acc, prec, reca, fmeas

# =============================================================================


def get_examples_quality(match_set, true_match_set, non_match_set,
                         true_non_match_set):
  """Count number of true and false positives and negatives in the match and
     non-match sets based on the true matches and non-matches.

     Returns two percentage values (match and non-match set  accuracy), or 0
     is a match or non-match set is empty.
  """

  logging.info('Number of matches:          %d' % (len(match_set)))
  logging.info('Number of true matches:     %d' % (len(true_match_set)))
  logging.info('Number of non-matches:      %d' % (len(non_match_set)))
  logging.info('Number of true non-matches: %d' % (len(true_non_match_set)))

  tp = 0.0
  fp = 0.0
  tn = 0.0
  fn = 0.0

  for m in match_set:
    if m in true_match_set:
      tp += 1
    elif m in true_non_match_set:
      fp += 1
    else:
      print m

  for nm in non_match_set:
    if nm in true_non_match_set:
      tn += 1
    elif nm in true_match_set:
      fn += 1
    else:
      print nm

  assert (tp+fp) == len(match_set)
  assert (tn+fn) == len(non_match_set)

  logging.info('  TP=%.0f, FP=%.0f, TN=%.0f, FN=%.0f' % (tp, fp, tn, fn))

  if (len(match_set) > 0):
    match_acc = tp / (tp+fp)
  else:
    match_acc = 0.0

  if (len(non_match_set) > 0):
    non_match_acc = tn / (tn+fn)
  else:
    non_match_acc = 0.0

  logging.info('  Match accuracy:     %.2f' % (match_acc))
  logging.info('  Non-match accuracy: %.2f' % (non_match_acc))

  return match_acc, non_match_acc

# =============================================================================
