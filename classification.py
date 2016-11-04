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
# The Original Software is: "classification.py"
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

"""Module with classes for record pair classification.

   This module provides various classifiers, both based on supervised and
   unsupervised methods, that classify weight vectors into 'matches' and
   'non-matches', and possibly 'possible matches' (some classifiers don't
   classify weight vectors into this third class).

     FellegiSunter     The classical Fellegi and Sunter classifier with two
                       thresholds.
     OptimalThreshold  A classifiers that uses the true match and non-match
                       status to optimally set threshold values.
     KMeans            Unsupervised K-means clustering algorithm with.
     FarthestFirst     Unsupervised farthest first clustering algorithm.
     SuppVecMachine    Supervised support vector machine (SVM) classifier.
     TwoStep           Unsupervised two-step classifier.
     TAILOR            Unsupervised hybrid classifier as described in the paper
                       TAILOR: A record linkage toolbox (Elfeky MG, Verykios
                       VS, Elmagarmid AK, ICDE, San Jose, 2002.

##
TODO: DecisionTree    Supervised decision tree induction based classifier.
##

   Creating and using a classifier normally consists of the following steps:
   - initialise  The classifier is initialised and trained if training data is
                 provided (i.e. if a weight vector dictionary and a match and
                 non-match set is given, the training method will be called).
   - train       Train the classifier using training data.
   - test        Testing the trained classifier using test data (with known
                 match and non-match status).
   - classify    Use the trained classifier to classify weight vectors with
                 unknown match status.

   Each classifier also has a cross_validate() method that allows evaluation of
   the classifier by conducting a cross validation.

   Additional auxiliary functions in this module that are related to record
   pair classification are:

     get_true_matches_nonmatches      Checks for all weight vectors in a weight
                                      vector dictionary if they correspond to
                                      true matches or true non-matches.
     extract_collapse_weight_vectors  A function that allows manipulation of
                                      the weight vectors in a weight vector
                                      dictionary, such as summing weight vector
                                      elements or filtering them out. Returns a
                                      modified weight vector dictionary.

   TODO:
   - Decision Tree based -> improve, make faster
   - EM clustering

   - have an argument collapse_vector [0,0,1,2,0,1] of same lengths as weight
     vector, which will take weights and summ them according to index numbers
     e.g. new_w-vec[0] = w_vec[0]+w_vec[1]+w_vec[4], new_w_vec[1] = ...
"""

# =============================================================================
# Import necessary modules (Python standard modules first, then Febrl modules)

import auxiliary
import mymath

import heapq
import logging
import math
import os
import random

#try:
#  import Numeric
#  imp_numeric = True
#except:
#  imp_numeric = False
#
#if (imp_numeric == True):
#  try:
#    import PyML
#    import PyML.datafunc
#    import PyML.svm
#    imp_pyml = True
#  except:
#    imp_pyml = False
#else:
#  imp_pyml = False

try:
  import svm
  imp_svm = True
except:
  imp_svm = False

#if (imp_pyml == False):
#  logging.warn('Cannot import Numeric and PyML modules')
if (imp_svm == False):
  logging.warn('Cannot import svm module')

# =============================================================================

class Classifier:
  """Base class for classifiers.

     All classifiers have the following instance variables, which can be set
     when a classifier is initialised:

       description          A string describing the classifier.
       train_w_vec_dict     A weight vector dictionary that will be used for
                            training.
       train_match_set      A set with record identifier pairs which are
                            assumed to be the match training examples.
       train_non_match_set  A set with record identifier pairs which are
                            assumed to be the non-match training examples.

     If the last three arguments (train_w_vec_dict, train_match_set, and
     train_non_match_set) are given when a classifier is initialised, it is
     trained straight away (so the 'train' method does not have to be called).
     Default is that these three values are set to None, so no training is done
     when a classifier is initialised.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, base_kwargs):
    """Constructor.
    """

    # General attributes for all classifiers.
    #
    self.description =     ''        # A description of the classifier.

    self.train_w_vec_dict = None     # The dictionary containing weight vectors
                                     # used for training.
    self.train_match_set = None      # A set with record identifier pairs that
                                     # are matches used for training.
    self.train_non_match_set = None  # A set with record identifier pairs that
                                     # are non-matches used for training.

    # Process base keyword arguments (all data set specific keywords were
    # processed in the derived class constructor)
    #
    for (keyword, value) in base_kwargs.items():

      if (keyword.startswith('desc')):
        auxiliary.check_is_string('description', value)
        self.description = value

      elif (keyword.startswith('train_w_vec')):
        auxiliary.check_is_dictionary('train_w_vec_dict', value)
        self.train_w_vec_dict = value

      elif (keyword.startswith('train_mat')):
        auxiliary.check_is_set('train_match_set', value)
        self.train_match_set = value
      elif (keyword.startswith('train_non')):
        auxiliary.check_is_set('train_non_match_set', value)
        self.train_non_match_set = value

      else:
        logging.exception('Illegal constructor argument keyword: '+keyword)
        raise Exception

  # ---------------------------------------------------------------------------

  def train(self, w_vec_dict, match_set, non_match_set):
    """Method to train a classifier using the given weight vector dictionary
       and match and non-match sets of record identifier pairs.

       See implementations in derived classes for details.
    """

    logging.exception('Override abstract method in derived class')
    raise Exception

  # ---------------------------------------------------------------------------

  def test(self, w_vec_dict, match_set, non_match_set):
    """Method to test a classifier using the given weight vector dictionary and
       match and non-match sets of record identifier pairs.

       Will return a confusion matrix as a list of the form: [TP, FN, FP, TN].

       See implementations in derived classes for details.
    """

    logging.exception('Override abstract method in derived class')
    raise Exception

  # ---------------------------------------------------------------------------

  def cross_validate(self, w_vec_dict, match_set, non_match_set):
    """Method to conduct a cross validation using the given weight vector
       dictionary and match and non-match sets of record identifier pairs.

       Will return a confusion matrix as a list of the form: [TP, FN, FP, TN].

       See implementations in derived classes for details.
    """

    logging.exception('Override abstract method in derived class')
    raise Exception

  # ---------------------------------------------------------------------------

  def classify(self, w_vec_dict):
    """Method to classify the given weight vector dictionary using the trained
       classifier.

       Will return three sets with record identifier pairs:
       1) match set
       2) non-match set
       3) possible match set  (this will always be empty for certain
                               classifiers that only classify into matches and
                               non-matches)

       See implementations in derived classes for details.
    """

    logging.exception('Override abstract method in derived class')
    raise Exception

  # ---------------------------------------------------------------------------

  def log(self, instance_var_list = None):
    """Write a log message with the basic classifier instance variables plus
       the instance variable provided in the given input list (assumed to
       contain pairs of names (strings) and values).
    """

    logging.info('')
    logging.info('Classifier:              "%s"' % (self.description))
    if (self.train_w_vec_dict != None):
      logging.info('  Number of weight vectors provided:              %d' % \
                   (len(self.train_w_vec_dict)))
    if (self.train_match_set != None):
      logging.info('  Number of match training examples provided:     %d' % \
                   (len(self.train_match_set)))
    if (self.train_non_match_set != None):
      logging.info('  Number of non-match training examples provided: %d' % \
                   (len(self.train_non_match_set)))

    if (instance_var_list != None):
      logging.info('  Classifier specific variables:')

      max_name_len = 0
      for (name, value) in instance_var_list:
        max_name_len = max(max_name_len, len(name))

      for (name, value) in instance_var_list:
        pad_spaces = (max_name_len-len(name))*' '
        logging.info('    %s %s' % (name+':'+pad_spaces, str(value)))


# =============================================================================

class FellegiSunter(Classifier):
  """Implements the classical Fellegi and Sunter classifier.

     This classifier sums all weights in a weight vector into one matching
     weight and then uses two threshold to classify this weight vector as
     either a match, non-match or a possible match.

     The arguments that have to be set when this classifier is initialised are:

       lower_threshold  All weight vectors with a summed matching weight below
                        this threshold are classified as non-matches.
       upper_threshold  All weight vectors with a summed matching weight above
                        this threshold are classified as matches.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'lower_threshold' and 'upper_threshold'
       arguments first, then call the base class constructor.
    """

    self.lower_threshold = None
    self.upper_threshold = None

    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('lower_t')):
        auxiliary.check_is_number('lower_threshold', value)
        self.lower_threshold = value
      elif (keyword.startswith('upper_t')):
        auxiliary.check_is_number('upper_threshold', value)
        self.upper_threshold = value

      else:
        base_kwargs[keyword] = value

    Classifier.__init__(self, base_kwargs)  # Initialise base class

    # Check threshold values are set and valid - - - - - - - - - - - - - - - -
    #
    auxiliary.check_is_number('lower_threshold', self.lower_threshold)
    auxiliary.check_is_number('upper_threshold', self.upper_threshold)

    if (self.lower_threshold > self.upper_threshold):
      logging.exception('Lower threshold is larger than upper threshold: ' + \
                        '%.3f / %.3f' % (lower_threshold, upper_threshold))
      raise Exception

    self.log([('Lower threshold', self.lower_threshold),
              ('Upper threshold', self.upper_threshold)])  # Log a message

    # If the weight vector dictionary and both match and non-match sets - - - -
    # are given start the training process
    #
    if ((self.train_w_vec_dict != None) and (self.train_match_set != None) \
        and (self.train_non_match_set != None)):
      logging.info('Train Fellegi and Sunter classifier: "%s"' % \
                   (self.description))
      logging.info('  Nothing needs to be done.')

  # ---------------------------------------------------------------------------

  def train(self, w_vec_dict, match_set, non_match_set):
    """Method to train a classifier using the given weight vector dictionary
       and match and non-match sets of record identifier pairs.

       Nothing needs to be done.
    """

    logging.info('')
    logging.info('Train Fellegi and Sunter classifier: "%s"' % \
                 (self.description))
    logging.info('  Nothing needs to be done.')

  # ---------------------------------------------------------------------------

  def test(self, w_vec_dict, match_set, non_match_set):
    """Method to test a classifier using the given weight vector dictionary and
       match and non-match sets of record identifier pairs.

       Will return a confusion matrix as a list of the form: [TP, FN, FP, TN].

       The numbers returned in the confusion matrix will only consider the
       classified matches and non-matches, but not the possible matches.

       TODO:
       - Is this correct, does this makes sense?
    """

    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)
    auxiliary.check_is_set('match_set', match_set)
    auxiliary.check_is_set('non_match_set', non_match_set)

    # Check that match and non-match sets are separate and do cover all weight
    # vectors given
    #
    if (len(match_set.intersection(non_match_set)) > 0):
      logging.exception('Intersection of match and non-match set not empty')
      raise Exception
    if ((len(match_set)+len(non_match_set)) != len(w_vec_dict)):
      logging.exception('Weight vector dictionary of different length than' + \
                        ' summed lengths of match and non-match sets: ' + \
                        '%d / %d+%d=%d' % (len(w_vec_dict), len(match_set),
                        len(non_match_set), len(match_set)+len(non_match_set)))
      raise Exception

    logging.info('')
    logging.info('Testing Fellegi and Sunter classifier using %d weight ' % \
                 (len(w_vec_dict))+'vectors')
    logging.info('  Match and non-match sets with %d and %d entries' % \
                 (len(match_set), len(non_match_set)))

    num_true_m =   0
    num_false_m =  0
    num_true_nm =  0
    num_false_nm = 0
    num_poss_m =   0

    for (rec_id_tuple, w_vec) in w_vec_dict.iteritems():
      w_sum = sum(w_vec)

      if (w_sum > self.upper_threshold):
        if (rec_id_tuple in match_set):
          num_true_m += 1
        else:
          num_false_m += 1

      elif (w_sum < self.lower_threshold):
        if (rec_id_tuple in non_match_set):
          num_true_nm += 1
        else:
          num_false_nm += 1

      else:
        num_poss_m += 1

    assert (num_true_m+num_false_nm+num_false_m+num_true_nm+num_poss_m) == \
           len(w_vec_dict)

    logging.info('  Results: TP = %d, FN = %d, FP = %d, TN = %d; ' % \
                 (num_true_m,num_false_nm,num_false_m,num_true_nm) + \
                 'possible matches = %d' % (num_poss_m))

    return [num_true_m, num_false_nm, num_false_m, num_true_nm]

  # ---------------------------------------------------------------------------

  def cross_validate(self, w_vec_dict, match_set, non_match_set, n=10):
    """Method to conduct a cross validation using the given weight vector
       dictionary and match and non-match sets of record identifier pairs.

       Will return a confusion matrix as a list of the form: [TP, FN, FP, TN].

       The Fellegi and Sunter classifier cannot perform cross validation, as
       the thresholds are set by the user. Therefore, in this method the given
       weight vectors are tested once only by calling the 'test' method.

       See documentation of 'test' for more information.
    """

    logging.info('')
    logging.info('Cross validation for Fellegi and Sunter is the same as' + \
                 'testing.')

    return self.test(w_vec_dict, match_set, non_match_set)

  # ---------------------------------------------------------------------------

  def classify(self, w_vec_dict):
    """Method to classify the given weight vector dictionary using the trained
       classifier.

       Will return three sets with record identifier pairs: 1) match set,
       2) non-match set, and 3) possible match set
    """

    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)

    logging.info('')
    logging.info('Classify %d weight vectors using Fellegi and Sunter ' % \
                 (len(w_vec_dict))+'classifier')

    match_set =      set()
    non_match_set =  set()
    poss_match_set = set()

    for (rec_id_tuple, w_vec) in w_vec_dict.iteritems():
      w_sum = sum(w_vec)

      if (w_sum > self.upper_threshold):
        match_set.add(rec_id_tuple)

      elif (w_sum < self.lower_threshold):
        non_match_set.add(rec_id_tuple)

      else:
        poss_match_set.add(rec_id_tuple)

    assert (len(match_set) + len(non_match_set) + len(poss_match_set)) == \
           len(w_vec_dict)

    logging.info('Classified %d weight vectors: %d as matches, %d as ' % \
                 (len(w_vec_dict), len(match_set), len(non_match_set)) + \
                 'non-matches, and %d as possible matches' % \
                 (len(poss_match_set)))

    return match_set, non_match_set, poss_match_set


# =============================================================================

class OptimalThreshold(Classifier):
  """Implements a classifier that has access to the true matches and true
     non-matches, and can thus set an optimal threshold (one for each weight
     vector element / dimension) so that the sum of either the number of (a)
     false matches and false non-matches, (b) false matches only, or (c) false
     non-matches only, is minimised.

     The weight vector values are binned first (according to the value of the
     'bin_width' argument) and then the optimal threshold is used on the binned
     values in each vector element (dimension).

     The arguments that have to be set when this classifier is initialised are:

       bin_width   The numerical width of the bins to be used.
       min_method  A string which designates what to minimise. This can either
                   be 'pos-neg' (default) in which case the sum of both false
                   matches and false non-matches will be minimised, or 'pos' or
                   'neg', in which case only the corresponding false numbers
                   will be minimised. Default is 'pos-neg'.
    """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'bin_width' and 'min_method' arguments first,
       then call the base class constructor.
    """

    self.bin_width =  None
    self.min_method = 'pos-neg'

    self.opt_threshold_list = None  # Will be calculated in training phase

    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('bin_w')):
        auxiliary.check_is_number('bin_width', value)
        auxiliary.check_is_positive('bin_width', value)
        self.bin_width = value

      elif (keyword.startswith('min_m')):
        auxiliary.check_is_string('min_method', value)
        if (value not in ['pos-neg', 'pos', 'neg']):
          logging.exception('Value of "min_method" is not one of "pos-neg"' + \
                            ', "pos", or "neg": %s' % (value))
          raise Exception
        self.min_method = value

      else:
        base_kwargs[keyword] = value

    Classifier.__init__(self, base_kwargs)  # Initialise base class

    # Check the bin width value is set - - - - - - - - - - - - - - - - - - - -
    #
    auxiliary.check_is_number('bin_width', self.bin_width)
    auxiliary.check_is_positive('bin_width', self.bin_width)

    self.log([('Bin width', self.bin_width),
              ('Minimise method', self.min_method)])  # Log a message

    # If the weight vector dictionary and both match and non-match sets - - - -
    # are given start the training process
    #
    if ((self.train_w_vec_dict != None) and (self.train_match_set != None) \
        and (self.train_non_match_set != None)):
      self.train(self.train_w_vec_dict, self.train_match_set,
                 (self.train_non_match_set))

  # ---------------------------------------------------------------------------

  def train(self, w_vec_dict, match_set, non_match_set):
    """Method to train a classifier using the given weight vector dictionary
       and match and non-match sets of record identifier pairs.

       Note that all weight vectors must either be in the match or the
       non-match training sets.

       This method will calculate the optimal threshold for each vector element
       (dimension).
    """

    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)
    auxiliary.check_is_set('match_set', match_set)
    auxiliary.check_is_set('non_match_set', non_match_set)

    # Check that match and non-match sets are separate and do cover all weight
    # vectors given
    #
    if (len(match_set.intersection(non_match_set)) > 0):
      logging.exception('Intersection of match and non-match set not empty')
      raise Exception
    if ((len(match_set)+len(non_match_set)) != len(w_vec_dict)):
      logging.exception('Weight vector dictionary of different length than' + \
                        ' summed lengths of match and non-match sets: ' + \
                        '%d / %d+%d=%d' % (len(w_vec_dict), len(match_set),
                        len(non_match_set), len(match_set)+len(non_match_set)))
      raise Exception

    self.train_w_vec_dict =    w_vec_dict  # Save
    self.train_match_set =     match_set
    self.train_non_match_set = non_match_set

    # Get a random vector dictionary element to get dimensionality of vectors
    #
    (rec_id_tuple, w_vec) = w_vec_dict.popitem()
    v_dim = len(w_vec)
    w_vec_dict[rec_id_tuple] = w_vec  # Put back in

    logging.info('Train optimal threshold classifier using %d weight ' % \
                 (len(w_vec_dict))+'vectors')
    logging.info('  Match and non-match sets with %d and %d entries' % \
                 (len(match_set), len(non_match_set)))
    logging.info('  Dimensionality:   %d' % (v_dim))

    # One dictionary with binned weights and their counts per dimension - - - -
    #
    match_weight_dict_list =     []
    non_match_weight_dict_list = []

    for i in range(v_dim):
      match_weight_dict_list.append({})
      non_match_weight_dict_list.append({})

    # Go through all weight vectors and put them into match or non-match bins -
    #
    for (rec_id_tuple, w_vec) in w_vec_dict.iteritems():

      for i in range(v_dim):
        match_dict =     match_weight_dict_list[i]
        non_match_dict = non_match_weight_dict_list[i]

        # Bin by rounding values down
        #
        binned_w = w_vec[i] - (w_vec[i] % self.bin_width)

        if (rec_id_tuple in match_set):
          w_count = match_dict.get(binned_w, 0) + 1
          match_dict[binned_w] = w_count

        elif (rec_id_tuple in non_match_set):
          w_count = non_match_dict.get(binned_w, 0) + 1
          non_match_dict[binned_w] = w_count
        else:
          logging.exception('Record identifier tuple %s not in match sets!' % \
                            (str(rec_id_tuple)))
          raise Exception

    # Get minimum and maximum binned weights - - - - - - - - - - - - - - - - -
    #
    opt_threshold_list = []  # One optimal threshold per dimension

    for i in range(v_dim):
      match_dict =     match_weight_dict_list[i]
      non_match_dict = non_match_weight_dict_list[i]

      min_match_weight =  99999.99999
      max_match_weight = -99999.99999

      all_weights_set = set()

      for w in match_dict:
        min_match_weight = min(w, min_match_weight)
        max_match_weight = max(w, max_match_weight)
        all_weights_set.add(w)

      min_non_match_weight =  99999.99999
      max_non_match_weight = -99999.99999

      for w in non_match_dict:
        min_non_match_weight = min(w, min_non_match_weight)
        max_non_match_weight = max(w, max_non_match_weight)
        all_weights_set.add(w)

      all_weights_list = list(all_weights_set)
      all_weights_list.sort()

      assert min(min_match_weight,min_non_match_weight) == all_weights_list[0]
      assert max(max_match_weight,max_non_match_weight) == all_weights_list[-1]

      logging.info('  Minimum and maximum binnded weights in dimension' + \
                   ' %d: %.3f / %.3f' % (i, all_weights_list[0],
                                            all_weights_list[-1]))
      logging.info('      True match weights range:     %.3f to %.3f' % \
                   (min_match_weight, max_match_weight))
      logging.info('      True non-match weights range: %.3f to %.3f' % \
                   (min_non_match_weight, max_non_match_weight))

      # Go through weight count dictionaries to find optimal thresholds - - - -
      #
      tp = len(match_set)  # Set initial classification counts
      tn = 0
      fp = len(non_match_set)
      fn = 0

      if (self.min_method == 'pos-neg'):  # Init classification information
        min_num_wrong = fp+fn
      elif (self.min_method == 'pos'):
        min_num_wrong = fp
      else:
        min_num_wrong = fn

      opt_threshold = all_weights_list[0]

      for w in all_weights_list:
        m_count =  match_dict.get(w, 0)
        nm_count = non_match_dict.get(w, 0)

        tp -= m_count
        fn += m_count
        tn += nm_count
        fp -= nm_count

        if (self.min_method == 'pos-neg'):
          if ((fp+fn) < min_num_wrong):
            min_num_wrong = (fp+fn)
            opt_threshold = w

        elif (self.min_method == 'pos'):
          if (fp < min_num_wrong):
            min_num_wrong = fp
            opt_threshold = w

        else:  # Minimise false negatives
          if ((min_num_wrong == 0 ) and (fn > 0)):  # First time there are FM
            min_num_wrong = fn
            opt_threshold = w-self.bin_width

      opt_threshold_list.append(opt_threshold)

      if (self.min_method == 'neg'):
        min_num_wrong = 0  # Adjust, this is always possible with very low thr.

      logging.info('    Optimal threshold in dimension %d is %.3f' % \
                   (i, opt_threshold))
      logging.info('      Number of "%s" misclassifications: %d' % \
                   (self.min_method, min_num_wrong))

    self.opt_threshold_list = opt_threshold_list

  # ---------------------------------------------------------------------------

  def test(self, w_vec_dict, match_set, non_match_set):
    """Method to test a classifier using the given weight vector dictionary and
       match and non-match sets of record identifier pairs.

       Weight vectors will be assigned to matches or non-matches according to
       the summed distances for their values from the thresholds in each vector
       element (dimension).

       Will return a confusion matrix as a list of the form: [TP, FN, FP, TN].
    """

    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)
    auxiliary.check_is_set('match_set', match_set)
    auxiliary.check_is_set('non_match_set', non_match_set)

    # Check that match and non-match sets are separate and do cover all weight
    # vectors given
    #
    if (len(match_set.intersection(non_match_set)) > 0):
      logging.exception('Intersection of match and non-match set not empty')
      raise Exception
    if ((len(match_set)+len(non_match_set)) != len(w_vec_dict)):
      logging.exception('Weight vector dictionary of different length than' + \
                        ' summed lengths of match and non-match sets: ' + \
                        '%d / %d+%d=%d' % (len(w_vec_dict), len(match_set),
                        len(non_match_set), len(match_set)+len(non_match_set)))
      raise Exception

    # Get a random vector dictionary element to get dimensionality of vectors
    #
    (rec_id_tuple, w_vec) = w_vec_dict.popitem()
    v_dim = len(w_vec)
    w_vec_dict[rec_id_tuple] = w_vec  # Put back in

    logging.info('')
    logging.info('Testing optimal threshold classifier using %d weight ' % \
                 (len(w_vec_dict))+'vectors')
    logging.info('  Match and non-match sets with %d and %d entries' % \
                 (len(match_set), len(non_match_set)))
    logging.info('  Dimensionality:   %d' % (v_dim))

    num_true_m =   0
    num_false_m =  0
    num_true_nm =  0
    num_false_nm = 0

    for (rec_id_tuple, w_vec) in w_vec_dict.iteritems():
      w_sum = sum(w_vec)

      diff_sum = 0.0  # Sum of differences over vector elements (dimensions)

      for i in range(v_dim):

        # Get difference between this weight value and threshold
        #
        diff_sum += (w_vec[i] - self.opt_threshold_list[i])

      if (diff_sum >= 0.0):
        if (rec_id_tuple in match_set):
          num_true_m += 1
        else:
          num_false_m += 1
      else:
        if (rec_id_tuple in non_match_set):
          num_true_nm += 1
        else:
          num_false_nm += 1

    assert (num_true_m+num_false_nm+num_false_m+num_true_nm) == len(w_vec_dict)

    logging.info('  Results: TP = %d, FN = %d, FP = %d, TN = %d' % \
                 (num_true_m,num_false_nm,num_false_m,num_true_nm))

    return [num_true_m, num_false_nm, num_false_m, num_true_nm]

  # --------------------------------------------------------------------------

  def cross_validate(self, w_vec_dict, match_set, non_match_set, n=10):
    """Method to conduct a cross validation using the given weight vector
       dictionary and match and non-match sets of record identifier pairs.

       Will return a confusion matrix as a list of the form: [TP, FN, FP, TN].

       The cross validation approach randomly splits the weight vector
       dictionary into 'n' parts (and 'n' corresponding sub-set for matches and
       non-matches), and then generates 'n' optimal threshold classifiers,
       tests them and finally returns the average performance of these 'n'
       classifiers.

       At the end of the cross validation procedure the optimal thresholds will
       be set to the average values of the 'n' optimal thresholds (in each
       dimension).
    """

    auxiliary.check_is_integer('n', n)
    auxiliary.check_is_positive('n', n)
    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)
    auxiliary.check_is_set('match_set', match_set)
    auxiliary.check_is_set('non_match_set', non_match_set)

    # Check that match and non-match sets are separate and do cover all weight
    # vectors given
    #
    if (len(match_set.intersection(non_match_set)) > 0):
      logging.exception('Intersection of match and non-match set not empty')
      raise Exception
    if ((len(match_set)+len(non_match_set)) != len(w_vec_dict)):
      logging.exception('Weight vector dictionary of different length than' + \
                        ' summed lengths of match and non-match sets: ' + \
                        '%d / %d+%d=%d' % (len(w_vec_dict), len(match_set),
                        len(non_match_set), len(match_set)+len(non_match_set)))
      raise Exception

    # Get a random vector dictionary element to get dimensionality of vectors
    #
    (rec_id_tuple, w_vec) = w_vec_dict.popitem()
    v_dim = len(w_vec)
    w_vec_dict[rec_id_tuple] = w_vec  # Put back in

    logging.info('')
    logging.info('Conduct %d-fold cross validation on optimal threshold ' % \
                 (n) + 'classifier using %d weight vectors' % \
                 (len(w_vec_dict)))
    logging.info('  Match and non-match sets with %d and %d entries' % \
                 (len(match_set), len(non_match_set)))
    logging.info('  Dimensionality:   %d' % (v_dim))

    opt_thres = []  # Keep the threshold from all folds

    # Create the sub-sets of record identifier pairs for folds - - - - - - - -
    #
    rec_id_tuple_list = w_vec_dict.keys()
    random.shuffle(rec_id_tuple_list)
    fold_num_rec_id_tuple = max(1,int(round(float(len(rec_id_tuple_list))/n)))

    # Split the weight vector dictionary and match and non-match sets into
    # (lists containing one entry per fold) and only store test elements
    #
    w_vec_dict_test_list = []
    m_set_test_list =      []
    nm_set_test_list =     []

    for fold in range(n):
      w_vec_dict_test_list.append({})
      m_set_test_list.append(set())
      nm_set_test_list.append(set())

    for fold in range(n):

      # Calculate start and end indices for test elements for this fold
      #
      if (fold == (n-1)):  # The last fold, get remainder of list
        start = fold*fold_num_rec_id_tuple
        this_fold_test_ids = rec_id_tuple_list[start:]
      else:  # All other folds
        start = fold*fold_num_rec_id_tuple
        end = start+fold_num_rec_id_tuple
        this_fold_test_ids = rec_id_tuple_list[start:end]

      for rec_id_tuple in this_fold_test_ids:

        w_vec_dict_test_list[fold][rec_id_tuple] = w_vec_dict[rec_id_tuple]

        if (rec_id_tuple in match_set):
          m_set_test_list[fold].add(rec_id_tuple)
        else:
          nm_set_test_list[fold].add(rec_id_tuple)

      assert len(w_vec_dict_test_list[fold]) == len(this_fold_test_ids)
      assert len(m_set_test_list[fold]) + len(nm_set_test_list[fold]) == \
             len(this_fold_test_ids)

    # Loop over folds - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Generate training and test dictionaries and sets
    #
    for fold in range(n):  # First extract test record identifier pairs

      this_fold_test_m_set =       m_set_test_list[fold]
      this_fold_test_nm_set =      nm_set_test_list[fold]
      this_fold_test_w_vec_dict =  w_vec_dict_test_list[fold]

      this_fold_train_m_set =  match_set.difference(m_set_test_list[fold])
      this_fold_train_nm_set = non_match_set.difference(nm_set_test_list[fold])
      this_fold_train_w_vec_dict = {}
      for f2 in range(n):
        if (f2 != fold):
          this_fold_train_w_vec_dict.update(w_vec_dict_test_list[f2])

      assert len(this_fold_test_m_set) + len(this_fold_train_m_set) == \
             len(match_set)
      assert len(this_fold_test_nm_set) + len(this_fold_train_nm_set) == \
             len(non_match_set)
      assert len(this_fold_test_w_vec_dict) + \
             len(this_fold_train_w_vec_dict) == len(w_vec_dict)

     #assert this_fold_test_m_set.intersection(this_fold_train_m_set) == set()
     #assert this_fold_test_m_set.intersection(this_fold_test_nm_set) == set()
     #assert this_fold_test_m_set.intersection(this_fold_train_nm_set) == set()
     #assert this_fold_test_nm_set.intersection(this_fold_train_m_set) ==set()
     #assert this_fold_test_nm_set.intersection(this_fold_train_nm_set) ==set()
     #assert this_fold_train_m_set.intersection(this_fold_train_nm_set) ==set()

      # Train optimal thrshold classifier and save calculated thresholds
      #
      self.train(this_fold_train_w_vec_dict, this_fold_train_m_set,
                 this_fold_train_nm_set)
      opt_thres.append(self.opt_threshold_list)

      del this_fold_train_w_vec_dict

    # Calculate final averaged optimal thresholds - - - - - - - - - - - - - - -
    #
    self.opt_threshold_list = [0.0]*v_dim
    for fold in range(n):
      for i in range(v_dim):
        self.opt_threshold_list[i] += (opt_thres[fold][i]/float(n))

    logging.info('Optimal thresholds: %s' % \
                 (auxiliary.str_vector(self.opt_threshold_list)))

    # Test on complete weight vector dictionary
    #
    [num_true_m,num_false_nm,num_false_m,num_true_nm]= self.test(w_vec_dict,
                                                                 match_set,
                                                                 non_match_set)
    return [num_true_m, num_false_nm, num_false_m, num_true_nm]

  # ---------------------------------------------------------------------------

  def classify(self, w_vec_dict):
    """Method to classify the given weight vector dictionary using the trained
       classifier.

       Will return three sets with record identifier pairs: 1) match set,
       2) non-match set, and 3) possible match set.

       The possible match set will be empty, as this classifier classifies all
       weight vectors as either matches or non-matches.
    """

    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)

    # Get a random vector dictionary element to get dimensionality of vectors
    #
    (rec_id_tuple, w_vec) = w_vec_dict.popitem()
    v_dim = len(w_vec)
    w_vec_dict[rec_id_tuple] = w_vec  # Put back in

    logging.info('')
    logging.info('Classify %d weight vectors using optimal threshold ' % \
                 (len(w_vec_dict))+'classifier')

    match_set =      set()
    non_match_set =  set()
    poss_match_set = set()

    for (rec_id_tuple, w_vec) in w_vec_dict.iteritems():

      w_sum = sum(w_vec)

      diff_sum = 0.0  # Sum of differences over vector elements (dimensions)

      for i in range(v_dim):

        # Get difference between this weight value and threshold
        #
        diff_sum += (w_vec[i] - self.opt_threshold_list[i])

      if (diff_sum >= 0.0):
        match_set.add(rec_id_tuple)
      else:
        non_match_set.add(rec_id_tuple)

    assert (len(match_set) + len(non_match_set)) == len(w_vec_dict)

    logging.info('Classified %d weight vectors: %d as matches and %d as ' % \
                 (len(w_vec_dict), len(match_set), len(non_match_set)) + \
                 'non-matches')

    return match_set, non_match_set, poss_match_set


# =============================================================================

class KMeans(Classifier):
  """Implements the unsupervised K-means clustering algorithm with two options
     of how to set the initial cluster centroids.

     Two clusters will be generated, one for matches and one for non-matches.

     Cluster centroids can be initialised either by randomly picking two weight
     vectors, or by taking the maximum and the minimum values (for each vector
     element / dimension) from the given weight vectors.

     When clustering (training) is conducted, it is also possible to only use a
     fraction of the given weight vectors for the clustering process through
     sampling.

     The arguments that have to be set when this classifier is initialised are:

       max_iter_count  The maximum number of iterations allowed.
       dist_measure    A function that calculates a distance measure between
                       two vectors (see the Febrl mymath.py module for such
                       functions).
       sample          A number between 0 and 100 that gives the percentage of
                       weight vectors that will be randomly selected and used
                       for clustering in the training process. If set to 100
                       (the default) then all given weight vectors will be
                       used.
       centroid_init   Method on how to initialise the centroids, can be either
                       'random' or 'min/max' (the default).
       fuzz_reg_thres  The fuzzy region threshold as described in:
                         L. Gu and R. Baxter: Decision models for record
                         linkage, Selected Papers from AusDM, Springer LNCS
                         3755, 2006.
                       See equation (5). Here, it is implemented slightly
                       different, in that the division by 2 is not done.
                       If the distance between a weight vector and the match
                       centroid is denoted by 'm_dist', and the distance
                       between a weight vector and the non-match centroid is
                       denoted by 'nm_dist', then the relative distance will be
                       calculated as:
                          rel_dist = abs(m_dist - nm_dist) / (m_dist + nm_dist)
                       with: 0 <= rel_dist <= 1.
                       Thus, the value of 'fuzz_reg_thres' can be set to a
                       value between 0 and 1. All weight vectors that have a
                       'rel_dist' < 'fuzz_reg_thres' will be inserted into the
                       possible match set.
                       Default value for 'fuzz_reg_thres' is None, in which
                       case no fuzzy region calculation will be done and no
                       weight vectors will be inserted into the possible match
                       set.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the K-means specific arguments first, then call the
       base class constructor.
    """

    self.max_iter_count = None
    self.dist_measure =   None
    self.sample =         100.0
    self.centroid_init =  'min/max'
    self.fuzz_reg_thres = None

    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('max_it')):
        auxiliary.check_is_integer('max_iter_count', value)
        auxiliary.check_is_positive('max_iter_count', value)
        self.max_iter_count = value

      elif (keyword.startswith('dist_m')):
        auxiliary.check_is_function_or_method('dist_measure', value)
        self.dist_measure = value

      elif (keyword.startswith('samp')):
        if (value != None):
          auxiliary.check_is_percentage('sample', value)
          self.sample = value

      elif (keyword.startswith('centr')):
        auxiliary.check_is_string('centroid_init', value)
        if (value not in ['random', 'min/max']):
          logging.exception('Value of "centroid_init" is not one of ' + \
                            '"random" or "min/max": %s' % (value))
          raise Exception
        self.centroid_init = value

      elif (keyword.startswith('fuzz')):
        if (value != None):
          auxiliary.check_is_normalised('fuzz_reg_thres', value)
          self.fuzz_reg_thres = value

      else:
        base_kwargs[keyword] = value

    Classifier.__init__(self, base_kwargs)  # Initialise base class

    # Check attribute values are set and valid - - - - - - - - - - - - - - - -
    #
    auxiliary.check_is_integer('max_iter_count', self.max_iter_count)
    auxiliary.check_is_positive('max_iter_count', self.max_iter_count)
    auxiliary.check_is_function_or_method('dist_measure', self.dist_measure)

    self.log([('Maximum iteration count', self.max_iter_count),
              ('Distance measure function', self.dist_measure),
              ('Sampling rate', self.sample),
              ('Centroid initialisation', self.centroid_init),
              ('Fuzzy match threshold', self.fuzz_reg_thres)]) # Log a message

    # If the weight vector dictionary and both match and non-match sets - - - -
    # are given start the training process
    #
    if ((self.train_w_vec_dict != None) and (self.train_match_set != None) \
        and (self.train_non_match_set != None)):
      self.train(self.train_w_vec_dict, self.train_match_set,
                 (self.train_non_match_set))

  # ---------------------------------------------------------------------------

  def train(self, w_vec_dict, match_set, non_match_set):
    """Method to train a classifier using the given weight vector dictionary.
       Note that the given match and non-match sets of record identifier pairs
       will not be used (unsupervised training).

       This method will calculate two cluster centroids (one for matches and
       one for non-matches), possibly using a sampled sub-set of the weight
       vectors given.
    """

    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)

    self.train_w_vec_dict =    w_vec_dict  # Save
    self.train_match_set =     match_set
    self.train_non_match_set = non_match_set

    # Get a random vector dictionary element to get dimensionality of vectors
    #
    (rec_id_tuple, w_vec) = w_vec_dict.popitem()
    v_dim = len(w_vec)
    w_vec_dict[rec_id_tuple] = w_vec  # Put back in

    logging.info('Train K-means classifier using %d weight vectors' % \
                 (len(w_vec_dict)))
    logging.info('  Dimensionality:   %d' % (v_dim))

    # Sample the weight vectors - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.sample == 100.0):
      use_w_vec_dict = w_vec_dict

    else:
      num_w_vec_sample = max(2, int(len(w_vec_dict)*self.sample/100.0))

      use_w_vec_dict = {}  # Create a new weight vector dictionary with samples

      rec_id_tuple_sample = random.sample(w_vec_dict.keys(),num_w_vec_sample)
      assert len(rec_id_tuple_sample) == num_w_vec_sample

      for rec_id_tuple in rec_id_tuple_sample:
        use_w_vec_dict[rec_id_tuple] = w_vec_dict[rec_id_tuple]

    logging.info('  Number of weight vectors to be used for clustering: %d' % \
                 (len(use_w_vec_dict)))

    zero_w_vec = [0.0]*v_dim  # Weight vector with all zeros

    # Initialise the cluster centroid - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.centroid_init == 'random'):
      rec_id_tuple_list = use_w_vec_dict.keys()

      m_centroid_rec_id = random.choice(rec_id_tuple_list)
      rec_id_tuple_list.remove(m_centroid_rec_id)
      nm_centroid_rec_id = random.choice(rec_id_tuple_list)
      del rec_id_tuple_list

      m_centroid =  use_w_vec_dict[m_centroid_rec_id]
      nm_centroid = use_w_vec_dict[nm_centroid_rec_id]

      # Make sure match centroid has larger values than non-match centroid
      #
      m_centr_zero_dist =  self.dist_measure(zero_w_vec, m_centroid)
      nm_centr_zero_dist = self.dist_measure(zero_w_vec, nm_centroid)

      # If match-centroid closer to zero than nno-match-centroid then swap
      #
      if (m_centr_zero_dist < nm_centr_zero_dist):
        m_centroid =  use_w_vec_dict[nm_centroid_rec_id]
        nm_centroid = use_w_vec_dict[m_centroid_rec_id]

    else:  # Get the minimum and maximum values in each weight vector element
      m_centroid =  [-999.99]*v_dim
      nm_centroid = [999.99]*v_dim

      for w_vec in use_w_vec_dict.itervalues():
        for i in range(v_dim):
          m_centroid[i] =  max(w_vec[i], m_centroid[i])
          nm_centroid[i] = min(w_vec[i], nm_centroid[i])

    logging.info('Initial cluster centroids using method "%s":' % \
                 (self.centroid_init))
    logging.info('  Initial match centroid:     %s' % \
                 (auxiliary.str_vector(m_centroid)))
    logging.info('  Initial non-match centroid: %s' % \
                 (auxiliary.str_vector(nm_centroid)))

    # Start iterations - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    cluster_assign_dict = {}  # Dictionary with cluster assignments

    iter_cnt = 1  # Iteration counter

    num_changed = 1

    while (num_changed > 0) and (iter_cnt < self.max_iter_count):

      num_changed =     0
      new_m_centroid =  [0.0]*v_dim  # Summed new distances
      new_nm_centroid = [0.0]*v_dim

      # Step 1: Calculate cluster membership for each weight vector - - - - - -
      #
      num_m =  0  # Number of weight vectors assigned to matches
      num_nm = 0  # Number of weight vectors assigned to non-matches

      for (rec_id_tuple, w_vec) in use_w_vec_dict.iteritems():

        m_dist =  self.dist_measure(w_vec, m_centroid)
        nm_dist = self.dist_measure(w_vec, nm_centroid)

        if (m_dist < nm_dist):  # Assign to cluster M (matches)
          old_assign = cluster_assign_dict.get(rec_id_tuple, 'X')
          if (old_assign != 'M'):
            num_changed += 1
          cluster_assign_dict[rec_id_tuple] = 'M'
          num_m += 1

          for i in range(v_dim):  # Add to summed cluster distances
            new_m_centroid[i] += w_vec[i]

        else:  # Assign to cluster NM (non-matches)
          old_assign = cluster_assign_dict.get(rec_id_tuple, 'X')
          if (old_assign != 'NM'):
            num_changed += 1
          cluster_assign_dict[rec_id_tuple] = 'NM'
          num_nm += 1

          for i in range(v_dim):  # Add to summed cluster distances
            new_nm_centroid[i] += w_vec[i]

      num_all = len(cluster_assign_dict)

      if ((num_m + num_nm) != num_all):
        logging.exception('Not all %d weight vectors assigned: M=%d, U=%d' % \
                          (num_all, num_m, num_nm))
        raise Exception

      if (num_m == 0) or (num_nm == 0):
        logging.warn('One cluster is empty: matches=%d, non-matches=%d' % \
                     (num_m, num_nm))
        break  # Stop K-means iterations

      # Step 2: Calculate new cluster centroids - - - - - - - - - - - - - - - -
      #
      for i in range(v_dim):  # Normalise new centroids
        new_m_centroid[i] /=  float(num_m)
        new_nm_centroid[i] /= float(num_nm)

      m_centroid =  new_m_centroid  # Assign new centroids
      nm_centroid = new_nm_centroid

      logging.info('Iteration %d: %d vectors changed cluster assignment' % \
            (iter_cnt, num_changed))

      iter_cnt += 1

    self.m_centroid =  m_centroid  # Save for later use
    self.nm_centroid = nm_centroid

    logging.info('Final cluster centroids using method "%s":' % \
                 (self.centroid_init))
    logging.info('  Match centroid:     %s' % \
                 (auxiliary.str_vector(m_centroid)))
    logging.info('  Non-match centroid: %s' % \
                 (auxiliary.str_vector(nm_centroid)))
    logging.info('  Cluster sizes: M=%d, NM=%d' % (num_m, num_nm))

  # ---------------------------------------------------------------------------

  def test(self, w_vec_dict, match_set, non_match_set):
    """Method to test a classifier using the given weight vector dictionary and
       match and non-match sets of record identifier pairs.

       Weight vectors will be assigned to matches or non-matches according to
       their distances to the match and non-match centroids.

       Note that no weight vector will be assigned to the possible match set
       (even if the fuzzy region threshold 'fuzz_reg_thres' has been set.

       Will return a confusion matrix as a list of the form: [TP, FN, FP, TN].
    """

    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)
    auxiliary.check_is_set('match_set', match_set)
    auxiliary.check_is_set('non_match_set', non_match_set)

    # Check that match and non-match sets are separate and do cover all weight
    # vectors given
    #
    if (len(match_set.intersection(non_match_set)) > 0):
      logging.exception('Intersection of match and non-match set not empty')
      raise Exception
    if ((len(match_set)+len(non_match_set)) != len(w_vec_dict)):
      logging.exception('Weight vector dictionary of different length than' + \
                        ' summed lengths of match and non-match sets: ' + \
                        '%d / %d+%d=%d' % (len(w_vec_dict), len(match_set),
                        len(non_match_set), len(match_set)+len(non_match_set)))
      raise Exception

    logging.info('')
    logging.info('Testing K-means classifier using %d weight vectors' % \
                 (len(w_vec_dict)))
    logging.info('  Match and non-match sets with %d and %d entries' % \
                 (len(match_set), len(non_match_set)))

    num_true_m =   0
    num_false_m =  0
    num_true_nm =  0
    num_false_nm = 0

    for (rec_id_tuple, w_vec) in w_vec_dict.iteritems():

      m_dist =  self.dist_measure(w_vec, self.m_centroid)
      nm_dist = self.dist_measure(w_vec, self.nm_centroid)

      if (m_dist < nm_dist):  # Assign to match cluster
        if (rec_id_tuple in match_set):
          num_true_m += 1
        else:
          num_false_m += 1

      else:
        if (rec_id_tuple in non_match_set):
          num_true_nm += 1
        else:
          num_false_nm += 1

    assert (num_true_m+num_false_nm+num_false_m+num_true_nm) == len(w_vec_dict)

    logging.info('  Results: TP = %d, FN = %d, FP = %d, TN = %d' % \
                 (num_true_m,num_false_nm,num_false_m,num_true_nm))

    return [num_true_m, num_false_nm, num_false_m, num_true_nm]

  # --------------------------------------------------------------------------

  def cross_validate(self, w_vec_dict, match_set, non_match_set, n=10):
    """Method to conduct a cross validation using the given weight vector
       dictionary and match and non-match sets of record identifier pairs.

       Will return a confusion matrix as a list of the form: [TP, FN, FP, TN].

       The cross validation approach splits the weight vector dictionary into
       'n' parts (and 'n' corresponding sub-set for matches and non-matches),
       and then generates 'n' K-means clusterings, tests them and finally
       returns the average performance of these 'n' classifiers, i.e. the
       final centroids will be set to the average of all 'n' centroids.
    """

    auxiliary.check_is_integer('n', n)
    auxiliary.check_is_positive('n', n)
    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)
    auxiliary.check_is_set('match_set', match_set)
    auxiliary.check_is_set('non_match_set', non_match_set)

    # Check that match and non-match sets are separate and do cover all weight
    # vectors given
    #
    if (len(match_set.intersection(non_match_set)) > 0):
      logging.exception('Intersection of match and non-match set not empty')
      raise Exception
    if ((len(match_set)+len(non_match_set)) != len(w_vec_dict)):
      logging.exception('Weight vector dictionary of different length than' + \
                        ' summed lengths of match and non-match sets: ' + \
                        '%d / %d+%d=%d' % (len(w_vec_dict), len(match_set),
                        len(non_match_set), len(match_set)+len(non_match_set)))
      raise Exception

    # Get a random vector dictionary element to get dimensionality of vectors
    #
    (rec_id_tuple, w_vec) = w_vec_dict.popitem()
    v_dim = len(w_vec)
    w_vec_dict[rec_id_tuple] = w_vec  # Put back in

    logging.info('')
    logging.info('Conduct %d-fold cross validation on K-means classifier ' % \
                 (n) + 'using %d weight vectors' % (len(w_vec_dict)))
    logging.info('  Match and non-match sets with %d and %d entries' % \
                 (len(match_set), len(non_match_set)))
    logging.info('  Dimensionality:   %d' % (v_dim))

    m_centroids =  []  # Keep the centroids from all folds
    nm_centroids = []

    # Create the sub-sets of record identifier pairs for folds - - - - - - - -
    #
    rec_id_tuple_list = w_vec_dict.keys()
    random.shuffle(rec_id_tuple_list)
    fold_num_rec_id_tuple = max(1,int(round(float(len(rec_id_tuple_list))/n)))

    # Split the weight vector dictionary and match and non-match sets into
    # (lists containing one entry per fold) and only store test elements
    #
    w_vec_dict_test_list = []
    m_set_test_list =      []
    nm_set_test_list =     []

    for fold in range(n):
      w_vec_dict_test_list.append({})
      m_set_test_list.append(set())
      nm_set_test_list.append(set())

    for fold in range(n):

      # Calculate start and end indices for test elements for this fold
      #
      if (fold == (n-1)):  # The last fold, get remainder of list
        start = fold*fold_num_rec_id_tuple
        this_fold_test_ids = rec_id_tuple_list[start:]
      else:  # All other folds
        start = fold*fold_num_rec_id_tuple
        end = start+fold_num_rec_id_tuple
        this_fold_test_ids = rec_id_tuple_list[start:end]

      for rec_id_tuple in this_fold_test_ids:

        w_vec_dict_test_list[fold][rec_id_tuple] = w_vec_dict[rec_id_tuple]

        if (rec_id_tuple in match_set):
          m_set_test_list[fold].add(rec_id_tuple)
        else:
          nm_set_test_list[fold].add(rec_id_tuple)

      assert len(w_vec_dict_test_list[fold]) == len(this_fold_test_ids)
      assert len(m_set_test_list[fold]) + len(nm_set_test_list[fold]) == \
             len(this_fold_test_ids)

    # Loop over folds - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Generate training and test dictionaries and sets
    #
    for fold in range(n):  # First extract test record identifier pairs

      this_fold_test_m_set =       m_set_test_list[fold]
      this_fold_test_nm_set =      nm_set_test_list[fold]
      this_fold_test_w_vec_dict =  w_vec_dict_test_list[fold]

      this_fold_train_m_set =  match_set.difference(m_set_test_list[fold])
      this_fold_train_nm_set = non_match_set.difference(nm_set_test_list[fold])
      this_fold_train_w_vec_dict = {}
      for f2 in range(n):
        if (f2 != fold):
          this_fold_train_w_vec_dict.update(w_vec_dict_test_list[f2])

      assert len(this_fold_test_m_set) + len(this_fold_train_m_set) == \
             len(match_set)
      assert len(this_fold_test_nm_set) + len(this_fold_train_nm_set) == \
             len(non_match_set)
      assert len(this_fold_test_w_vec_dict) + \
             len(this_fold_train_w_vec_dict) == len(w_vec_dict)

     #assert this_fold_test_m_set.intersection(this_fold_train_m_set) == set()
     #assert this_fold_test_m_set.intersection(this_fold_test_nm_set) == set()
     #assert this_fold_test_m_set.intersection(this_fold_train_nm_set) == set()
     #assert this_fold_test_nm_set.intersection(this_fold_train_m_set) ==set()
     #assert this_fold_test_nm_set.intersection(this_fold_train_nm_set) ==set()
     #assert this_fold_train_m_set.intersection(this_fold_train_nm_set) ==set()

      # Train K-means classifier and save calculated centroids
      #
      self.train(this_fold_train_w_vec_dict, this_fold_train_m_set,
                 this_fold_train_nm_set)

      m_centroids.append(self.m_centroid)
      nm_centroids.append(self.nm_centroid)

      del this_fold_train_w_vec_dict

    # Calculate final averaged centroids - - - - - - - - - - - - - - - - - - -
    #
    self.m_centroid =  [0.0]*v_dim
    self.nm_centroid = [0.0]*v_dim

    for fold in range(n):
      for i in range(v_dim):
        self.m_centroid[i] +=  (m_centroids[fold][i]/float(n))
        self.nm_centroid[i] += (nm_centroids[fold][i]/float(n))

    logging.info('Final cluster centroids using method "%s":' % \
                 (self.centroid_init))
    logging.info('  Match centroid:     %s' % \
                 (auxiliary.str_vector(self.m_centroid)))
    logging.info('  Non-match centroid: %s' % \
                 (auxiliary.str_vector(self.nm_centroid)))

    # Test on complete weight vector dictionary
    #
    [num_true_m,num_false_nm,num_false_m,num_true_nm]= self.test(w_vec_dict,
                                                                 match_set,
                                                                 non_match_set)
    return [num_true_m, num_false_nm, num_false_m, num_true_nm]

  # ---------------------------------------------------------------------------

  def classify(self, w_vec_dict):
    """Method to classify the given weight vector dictionary using the trained
       classifier.

       Will return three sets with record identifier pairs: 1) match set,
       2) non-match set, and 3) possible match set.

       If the fuzzy region threshold 'fuzz_reg_thres' has been set, weight
       vectors within this fuzzy region will be assigned to the possible match
       set, otherwise the possible match set will be empty.
    """

    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)

    logging.info('')
    logging.info('Classify %d weight vectors using K-means classifier' % \
                 (len(w_vec_dict)))

    match_set =      set()
    non_match_set =  set()
    poss_match_set = set()

    for (rec_id_tuple, w_vec) in w_vec_dict.iteritems():

      m_dist =  self.dist_measure(w_vec, self.m_centroid)
      nm_dist = self.dist_measure(w_vec, self.nm_centroid)

      if (self.fuzz_reg_thres == None):
        if (m_dist < nm_dist):  # Assign to match cluster
          match_set.add(rec_id_tuple)
        else:
          non_match_set.add(rec_id_tuple)

      else:  # Check if weight vector is in fuzzy region
        rel_dict = abs(m_dist - nm_dist) / (m_dist + nm_dist)

        if (rel_dict < self.fuzz_reg_thres):  # Assign to possible matches
          poss_match_set.add(rec_id_tuple)
        elif (m_dist < nm_dist):  # Assign to matches
          match_set.add(rec_id_tuple)
        else:
          non_match_set.add(rec_id_tuple)

    assert (len(match_set) + len(non_match_set) + len(poss_match_set)) == \
            len(w_vec_dict)

    logging.info('Classified %d weight vectors: %d as matches, %d as ' % \
                 (len(w_vec_dict), len(match_set), len(non_match_set)) + \
                 'non-matches, and %d as possible matches' % \
                 (len(poss_match_set)))

    return match_set, non_match_set, poss_match_set


# =============================================================================

class FarthestFirst(Classifier):
  """Implements several variations of the farthest first clustering algorithm
     with different possibilities of how to calculate the centroids.

     Calculates two centroids in the training method, and uses similar methods
     as the K-means classifier for cross validation, testing and
     classification.

     Two clusters will be generated, one for matches and one for non-matches.

     Cluster centroids can be initialised either by (a) the traditional
     farthest first approach (number of times a farthest vector is chosen is
     set to 10, but this can be changed within the code), or (b) using the two
     weight vectors with maximum and minimum values as match and non-match
     centroids, or (c) using the weight vector with maximum values as match
     centroid and the mode of all weight vectors as the non-match centroid. For
     more details please see:
       Goiser K. and Christen, P: Towards automated record linkage,
       Australasian Data Mining Conference (AusDM'06), Sydney, Conferences in
       Research and Practice in Information Technology (CRPIT), vol. 61.

     When clustering (training) is conducted, it is also possible to only use a
     fraction of the given weight vectors for the clustering process through
     sampling.

     The arguments that have to be set when this classifier is initialised are:

       dist_measure    A function that calculates a distance measure between
                       two vectors (see the Febrl mymath.py module for such
                       functions).
       sample          A number between 0 and 100 that gives the percentage of
                       weight vectors that will be randomly selected and used
                       for clustering in the training process. If set to 100
                       (the default) then all given weight vectors will be
                       used.
       centroid_init   Method on how to initialise the centroids, can be either
                       'traditional' (the default), 'min/max', or 'mode/max'.
       fuzz_reg_thres  The fuzzy region threshold, see K-means classifier for
                       more detailed information
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the farthest first specific arguments first, then
       call the base class constructor.
    """

    self.dist_measure =   None
    self.sample =         100.0
    self.centroid_init =  'traditional'
    self.fuzz_reg_thres = None
    self.num_choices =    10  # Number of choices for the traditional approach

    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('dist_m')):
        auxiliary.check_is_function_or_method('dist_measure', value)
        self.dist_measure = value

      elif (keyword.startswith('samp')):
        if (value != None):
          auxiliary.check_is_percentage('sample', value)
          self.sample = value

      elif (keyword.startswith('centr')):
        auxiliary.check_is_string('centroid_init', value)
        if (value not in ['traditional', 'min/max', 'mode/max']):
          logging.exception('Value of "centroid_init" is not one of ' + \
                            '"traditional, "min/max" or "mode/max": %s' % \
                            (value))
          raise Exception
        self.centroid_init = value

      elif (keyword.startswith('fuzz')):
        if (value != None):
          auxiliary.check_is_normalised('fuzz_reg_thres', value)
          self.fuzz_reg_thres = value

      else:
        base_kwargs[keyword] = value

    Classifier.__init__(self, base_kwargs)  # Initialise base class

    # Check attribute values are set and valid - - - - - - - - - - - - - - - -
    #
    auxiliary.check_is_function_or_method('dist_measure', self.dist_measure)

    self.log([('Distance measure function', self.dist_measure),
              ('Sampling rate', self.sample),
              ('Centroid initialisation', self.centroid_init),
              ('Fuzzy match threshold', self.fuzz_reg_thres)]) # Log a message

    # If the weight vector dictionary and both match and non-match sets - - - -
    # are given start the training process
    #
    if ((self.train_w_vec_dict != None) and (self.train_match_set != None) \
        and (self.train_non_match_set != None)):
      self.train(self.train_w_vec_dict, self.train_match_set,
                 (self.train_non_match_set))

  # ---------------------------------------------------------------------------

  def train(self, w_vec_dict, match_set, non_match_set):
    """Method to train a classifier using the given weight vector dictionary.
       Note that the given match and non-match sets of record identifier pairs
       will not be used (unsupervised training).

       This method will calculate two cluster centroids (one for matches and
       one for non-matches), possibly using a sampled sub-set of the weight
       vectors given.
    """

    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)

    self.train_w_vec_dict =    w_vec_dict  # Save
    self.train_match_set =     match_set
    self.train_non_match_set = non_match_set

    # Get a random vector dictionary element to get dimensionality of vectors
    #
    (rec_id_tuple, w_vec) = w_vec_dict.popitem()
    v_dim = len(w_vec)
    w_vec_dict[rec_id_tuple] = w_vec  # Put back in

    logging.info('Train farthest first classifier using %d weight vectors' % \
                 (len(w_vec_dict)))
    logging.info('  Dimensionality:   %d' % (v_dim))

    # Sample the weight vectors - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.sample == 100.0):
      use_w_vec_dict = w_vec_dict

    else:
      num_w_vec_sample = max(2, int(len(w_vec_dict)*self.sample/100.0))

      use_w_vec_dict = {}  # Create a new weight vector dictionary with samples

      rec_id_tuple_sample = random.sample(w_vec_dict.keys(),num_w_vec_sample)
      assert len(rec_id_tuple_sample) == num_w_vec_sample

      for rec_id_tuple in rec_id_tuple_sample:
        use_w_vec_dict[rec_id_tuple] = w_vec_dict[rec_id_tuple]

    logging.info('  Number of weight vectors to be used for clustering: %d' % \
                 (len(use_w_vec_dict)))

    # Iniialise the cluster centroid - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.centroid_init == 'traditional'):

      # Select a weight vector as first centroid
      #
      (rec_id_tuple, w_vec) = use_w_vec_dict.popitem()
      use_w_vec_dict[rec_id_tuple] = w_vec  # Put back in

      centroid1 = w_vec

      # Search the farthest weight vector from the initial centroid
      #
      for i in range(self.num_choices):
        max_dist =          -1.0
        max_dist_centroid = None

        for (rec_id_tuple, w_vec) in use_w_vec_dict.iteritems():
          dist = self.dist_measure(centroid1, w_vec)

          if (dist > max_dist):
            max_dist =          dist
            max_dist_centroid = w_vec

        centroid2 = centroid1  # Update farthest away centroid
        centroid1 = max_dist_centroid

      # Determine which is the match and which non-match centroid
      # (assume higher weights correspond to matches!)
      #
      if (sum(centroid1) > sum(centroid2)):
        m_centroid =  centroid1
        nm_centroid = centroid2
      else:
        m_centroid =  centroid2
        nm_centroid = centroid1

    # Get the weight vectors with minimum and maximum values in each - - - - -
    # vector element
    #
    elif (self.centroid_init == 'min/max'):
      m_centroid_sum = -999999
      nm_centroid_sum = 999999

      # Assume a larger summed weight is a match, a lower summed weight a
      # non-match
      #
      for w_vec in use_w_vec_dict.itervalues():
        w_vec_sum = sum(w_vec)

        if (w_vec_sum > m_centroid_sum):
          m_centroid = w_vec
          m_centroid_sum = w_vec_sum

        if (w_vec_sum < nm_centroid_sum):
          nm_centroid = w_vec
          nm_centroid_sum = w_vec_sum

    else:  # Mode/max method - - - - - - - - - - - - - - - - - - - - - - - - -

      m_centroid_sum = -999

      nm_histograms = []  # Estimate mode with a histogram based apprach
      bin_width = 0.01  # Maximum 100 bins from 0.0 to 1.0
      for i in range(v_dim):  # One dictionary per dimension
        nm_histograms.append({})

      for w_vec in use_w_vec_dict.itervalues():
        w_vec_sum = sum(w_vec)

        if (w_vec_sum > m_centroid_sum):  # Get match weight vector
          m_centroid = w_vec
          m_centroid_sum = w_vec_sum

        for i in range(v_dim):
          binned_w = w_vec[i] - (w_vec[i] % bin_width)
          bin_count = nm_histograms[i].get(binned_w, 0) + 1
          nm_histograms[i][binned_w] = bin_count

      # Get bin with highest counts in each dimension
      #
      nm_centroid = []

      for i in range(v_dim):
        max_count = -1
        for (binned_w, count) in nm_histograms[i].iteritems():
          if (count > max_count):
           centroid_w = binned_w
           max_count = count
        nm_centroid.append(centroid_w)

    self.m_centroid =  m_centroid  # Save for later use
    self.nm_centroid = nm_centroid

    logging.info('Final cluster centroids using method "%s":' % \
                 (self.centroid_init))
    logging.info('  Match centroid:     %s' % \
                 (auxiliary.str_vector(m_centroid)))
    logging.info('  Non-match centroid: %s' % \
                 (auxiliary.str_vector(nm_centroid)))

  # ---------------------------------------------------------------------------

  def test(self, w_vec_dict, match_set, non_match_set):
    """Method to test a classifier using the given weight vector dictionary and
       match and non-match sets of record identifier pairs.

       Weight vectors will be assigned to matches or non-matches according to
       their distances to the match and non-match centroids.

       Note that no weight vector will be assigned to the possible match set
       (even if the fuzzy region threshold 'fuzz_reg_thres' has been set.

       Will return a confusion matrix as a list of the form: [TP, FN, FP, TN].
    """

    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)
    auxiliary.check_is_set('match_set', match_set)
    auxiliary.check_is_set('non_match_set', non_match_set)

    # Check that match and non-match sets are separate and do cover all weight
    # vectors given
    #
    if (len(match_set.intersection(non_match_set)) > 0):
      logging.exception('Intersection of match and non-match set not empty')
      raise Exception
    if ((len(match_set)+len(non_match_set)) != len(w_vec_dict)):
      logging.exception('Weight vector dictionary of different length than' + \
                        ' summed lengths of match and non-match sets: ' + \
                        '%d / %d+%d=%d' % (len(w_vec_dict), len(match_set),
                        len(non_match_set), len(match_set)+len(non_match_set)))
      raise Exception

    logging.info('')
    logging.info('Testing farthest first classifier using %d weight ' % \
                 (len(w_vec_dict))+'vectors')
    logging.info('  Match and non-match sets with %d and %d entries' % \
                 (len(match_set), len(non_match_set)))

    num_true_m =   0
    num_false_m =  0
    num_true_nm =  0
    num_false_nm = 0

    for (rec_id_tuple, w_vec) in w_vec_dict.iteritems():

      m_dist =  self.dist_measure(w_vec, self.m_centroid)
      nm_dist = self.dist_measure(w_vec, self.nm_centroid)

      if (m_dist < nm_dist):  # Assign to match cluster
        if (rec_id_tuple in match_set):
          num_true_m += 1
        else:
          num_false_m += 1

      else:
        if (rec_id_tuple in non_match_set):
          num_true_nm += 1
        else:
          num_false_nm += 1

    assert (num_true_m+num_false_nm+num_false_m+num_true_nm) == len(w_vec_dict)

    logging.info('  Results: TP = %d, FN = %d, FP = %d, TN = %d' % \
                 (num_true_m,num_false_nm,num_false_m,num_true_nm))

    return [num_true_m, num_false_nm, num_false_m, num_true_nm]

  # --------------------------------------------------------------------------

  def cross_validate(self, w_vec_dict, match_set, non_match_set, n=10):
    """Method to conduct a cross validation using the given weight vector
       dictionary and match and non-match sets of record identifier pairs.

       Will return a confusion matrix as a list of the form: [TP, FN, FP, TN].

       The cross validation approach splits the weight vector dictionary into
       'n' parts (and 'n' corresponding sub-set for matches and non-matches),
       and then generates 'n' farthest first clusterings, tests them and
       finally returns the average performance of these 'n' classifiers, i.e.
       the final centroids will be set to the average of all 'n' centroids.
    """

    auxiliary.check_is_integer('n', n)
    auxiliary.check_is_positive('n', n)
    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)
    auxiliary.check_is_set('match_set', match_set)
    auxiliary.check_is_set('non_match_set', non_match_set)

    # Check that match and non-match sets are separate and do cover all weight
    # vectors given
    #
    if (len(match_set.intersection(non_match_set)) > 0):
      logging.exception('Intersection of match and non-match set not empty')
      raise Exception
    if ((len(match_set)+len(non_match_set)) != len(w_vec_dict)):
      logging.exception('Weight vector dictionary of different length than' + \
                        ' summed lengths of match and non-match sets: ' + \
                        '%d / %d+%d=%d' % (len(w_vec_dict), len(match_set),
                        len(non_match_set), len(match_set)+len(non_match_set)))
      raise Exception

    # Get a random vector dictionary element to get dimensionality of vectors
    #
    (rec_id_tuple, w_vec) = w_vec_dict.popitem()
    v_dim = len(w_vec)
    w_vec_dict[rec_id_tuple] = w_vec  # Put back in

    logging.info('')
    logging.info('Conduct %d-fold cross validation on farthest first ' % (n) \
                 + 'classifier using %d weight vectors' % (len(w_vec_dict)))
    logging.info('  Match and non-match sets with %d and %d entries' % \
                 (len(match_set), len(non_match_set)))
    logging.info('  Dimensionality:   %d' % (v_dim))

    m_centroids =  []  # Keep the centroids from all folds
    nm_centroids = []

    # Create the sub-sets of record identifier pairs for folds - - - - - - - -
    #
    rec_id_tuple_list = w_vec_dict.keys()
    random.shuffle(rec_id_tuple_list)
    fold_num_rec_id_tuple = max(1,int(round(float(len(rec_id_tuple_list))/n)))

    # Split the weight vector dictionary and match and non-match sets into
    # (lists containing one entry per fold) and only store test elements
    #
    w_vec_dict_test_list = []
    m_set_test_list =      []
    nm_set_test_list =     []

    for fold in range(n):
      w_vec_dict_test_list.append({})
      m_set_test_list.append(set())
      nm_set_test_list.append(set())

    for fold in range(n):

      # Calculate start and end indices for test elements for this fold
      #
      if (fold == (n-1)):  # The last fold, get remainder of list
        start = fold*fold_num_rec_id_tuple
        this_fold_test_ids = rec_id_tuple_list[start:]
      else:  # All other folds
        start = fold*fold_num_rec_id_tuple
        end = start+fold_num_rec_id_tuple
        this_fold_test_ids = rec_id_tuple_list[start:end]

      for rec_id_tuple in this_fold_test_ids:

        w_vec_dict_test_list[fold][rec_id_tuple] = w_vec_dict[rec_id_tuple]

        if (rec_id_tuple in match_set):
          m_set_test_list[fold].add(rec_id_tuple)
        else:
          nm_set_test_list[fold].add(rec_id_tuple)

      assert len(w_vec_dict_test_list[fold]) == len(this_fold_test_ids)
      assert len(m_set_test_list[fold]) + len(nm_set_test_list[fold]) == \
             len(this_fold_test_ids)

    # Loop over folds - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Generate training and test dictionaries and sets
    #
    for fold in range(n):  # First extract test record identifier pairs

      this_fold_test_m_set =       m_set_test_list[fold]
      this_fold_test_nm_set =      nm_set_test_list[fold]
      this_fold_test_w_vec_dict =  w_vec_dict_test_list[fold]

      this_fold_train_m_set =  match_set.difference(m_set_test_list[fold])
      this_fold_train_nm_set = non_match_set.difference(nm_set_test_list[fold])
      this_fold_train_w_vec_dict = {}
      for f2 in range(n):
        if (f2 != fold):
          this_fold_train_w_vec_dict.update(w_vec_dict_test_list[f2])

      assert len(this_fold_test_m_set) + len(this_fold_train_m_set) == \
             len(match_set)
      assert len(this_fold_test_nm_set) + len(this_fold_train_nm_set) == \
             len(non_match_set)
      assert len(this_fold_test_w_vec_dict) + \
             len(this_fold_train_w_vec_dict) == len(w_vec_dict)

     #assert this_fold_test_m_set.intersection(this_fold_train_m_set) == set()
     #assert this_fold_test_m_set.intersection(this_fold_test_nm_set) == set()
     #assert this_fold_test_m_set.intersection(this_fold_train_nm_set) == set()
     #assert this_fold_test_nm_set.intersection(this_fold_train_m_set) ==set()
     #assert this_fold_test_nm_set.intersection(this_fold_train_nm_set) ==set()
     #assert this_fold_train_m_set.intersection(this_fold_train_nm_set) ==set()

      # Train fathest first classifier and save calculated centroids
      #
      self.train(this_fold_train_w_vec_dict, this_fold_train_m_set,
                 this_fold_train_nm_set)

      m_centroids.append(self.m_centroid)
      nm_centroids.append(self.nm_centroid)

    # Calculate final averaged centroids - - - - - - - - - - - - - - - - - - -
    #
    self.m_centroid =  [0.0]*v_dim
    self.nm_centroid = [0.0]*v_dim

    for fold in range(n):
      for i in range(v_dim):
        self.m_centroid[i] +=  (m_centroids[fold][i]/float(n))
        self.nm_centroid[i] += (nm_centroids[fold][i]/float(n))

    logging.info('Final cluster centroids using method "%s":' % \
                 (self.centroid_init))
    logging.info('  Match centroid:     %s' % \
                 (auxiliary.str_vector(self.m_centroid)))
    logging.info('  Non-match centroid: %s' % \
                 (auxiliary.str_vector(self.nm_centroid)))

    # Test on complete weight vector dictionary
    #
    [num_true_m,num_false_nm,num_false_m,num_true_nm]= self.test(w_vec_dict,
                                                                 match_set,
                                                                 non_match_set)
    return [num_true_m, num_false_nm, num_false_m, num_true_nm]

  # ---------------------------------------------------------------------------

  def classify(self, w_vec_dict):
    """Method to classify the given weight vector dictionary using the trained
       classifier.

       Will return three sets with record identifier pairs: 1) match set,
       2) non-match set, and 3) possible match set.

       If the fuzzy region threshold 'fuzz_reg_thres' has been set, weight
       vectors within this fuzzy region will be assigned to the possible match
       set, otherwise the possible match set will be empty.
    """

    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)

    logging.info('')
    logging.info('Classify %d weight vectors using farthest first ' % \
                 (len(w_vec_dict))+'classifier')

    match_set =      set()
    non_match_set =  set()
    poss_match_set = set()

    for (rec_id_tuple, w_vec) in w_vec_dict.iteritems():

      m_dist =  self.dist_measure(w_vec, self.m_centroid)
      nm_dist = self.dist_measure(w_vec, self.nm_centroid)

      if (self.fuzz_reg_thres == None):
        if (m_dist < nm_dist):  # Assign to match cluster
          match_set.add(rec_id_tuple)
        else:
          non_match_set.add(rec_id_tuple)

      else:  # Check if weight vector is in fuzzy region
        rel_dict = abs(m_dist - nm_dist) / (m_dist + nm_dist)

        if (rel_dict < self.fuzz_reg_thres):  # Assign to possible matches
          poss_match_set.add(rec_id_tuple)
        elif (m_dist < nm_dist):  # Assign to matches
          match_set.add(rec_id_tuple)
        else:
          non_match_set.add(rec_id_tuple)

    assert (len(match_set) + len(non_match_set) + len(poss_match_set)) == \
            len(w_vec_dict)

    logging.info('Classified %d weight vectors: %d as matches, %d as ' % \
                 (len(w_vec_dict), len(match_set), len(non_match_set)) + \
                 'non-matches, and %d as possible matches' % \
                 (len(poss_match_set)))

    return match_set, non_match_set, poss_match_set


# =============================================================================

class SuppVecMachine(Classifier):
  """Implements a classifier based on a support vector machine (SVM). The
     'libsvm' library and its Python interface (module svm.py) will be used.
     For more details and downloads please see:

       http://www.csie.ntu.edu.tw/~cjlin/libsvm

     If this module is not implemented this classifier will be be usable.

     Note that the cross_validation() method only provides performance measures
     but no trained SVM model that can be used for classifying weight vectors
     later on. The train() method needs to be used to get a trained SVM.

     It is possible to do random sampling of training data from all weight
     vectors using the 'sample' argument.

     The arguments that have to be set when this classifier is initialised are:
     (the kernel type will be mapped to corresponding svm.py argument).

       kernel_type  The kernel type from from libsvm. Default value LINEAR,
                    other possibilities are: POLY, RBF, SIGMOID.
       C            The 'C' parameter from libsvm. Default value is 10.
       sample       A number between 0 and 100 that gives the percentage of
                    weight vectors that will be randomly selected and used for
                    clustering in the training process. If set to 100 (the
                    default) then all given weight vectors will be used.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the 'svmlib' and 'sample' arguments first, then
       call the base class constructor.
    """

    # Check if svm module is installed or not
    #
    if (imp_svm == False):
      logging.exception('Module "svm.py" not installed, cannot use ' + \
                        'SuppVectorMach classifier')
      raise Exception

    self.svm_type =    svm.C_SVC
    self.kernel_type = 'LINEAR'
    self.C =           10
    self.svm_model =   None  # Will be set in train() method
    self.sample =      100.0

    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('kernel')):
        auxiliary.check_is_string('kernel_type', value)
        if (value not in ['LINEAR', 'POLY', 'RBF', 'SIGMOID']):
          logging.exception('Illegal value for kernel type: %s ' % (value) + \
                            '(possible are: LINEAR, POLY, RBF, SIGMOID)')
          raise Exception
        self.kernel_type = value

      elif (keyword == 'C'):
        auxiliary.check_is_number('C', value)
        auxiliary.check_is_not_negative('C', value)
        self.C = value

      elif (keyword.startswith('samp')):
        if (value != None):
          auxiliary.check_is_percentage('sample', value)
          self.sample = value

      else:
        base_kwargs[keyword] = value

    Classifier.__init__(self, base_kwargs)  # Initialise base class

    self.log([('SVM kernel type', self.kernel_type),
              ('C',               self.C),
              ('Sampling rate',   self.sample)])  # Log a message

    # If the weight vector dictionary and both match and non-match sets - - - -
    # are given start the training process
    #
    if ((self.train_w_vec_dict != None) and (self.train_match_set != None) \
        and (self.train_non_match_set != None)):
      self.train(self.train_w_vec_dict, self.train_match_set,
                 (self.train_non_match_set))

  # ---------------------------------------------------------------------------

  def train(self, w_vec_dict, match_set, non_match_set):
    """Method to train a classifier using the given weight vector dictionary
       and match and non-match sets of record identifier pairs.

       Note that all weight vectors must either be in the match or the
       non-match training sets.
    """

    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)
    auxiliary.check_is_set('match_set', match_set)
    auxiliary.check_is_set('non_match_set', non_match_set)

    # Check that match and non-match sets are separate and do cover all weight
    # vectors given
    #
    if (len(match_set.intersection(non_match_set)) > 0):
      logging.exception('Intersection of match and non-match set not empty')
      raise Exception
    if ((len(match_set)+len(non_match_set)) != len(w_vec_dict)):
      logging.exception('Weight vector dictionary of different length than' + \
                        ' summed lengths of match and non-match sets: ' + \
                        '%d / %d+%d=%d' % (len(w_vec_dict), len(match_set),
                        len(non_match_set), len(match_set)+len(non_match_set)))
      raise Exception

    self.train_w_vec_dict =    w_vec_dict  # Save
    self.train_match_set =     match_set
    self.train_non_match_set = non_match_set

    logging.info('Train SVM classifier using %d weight vectors' % \
                 (len(w_vec_dict)))
    logging.info('  Match and non-match sets with %d and %d entries' % \
                 (len(match_set), len(non_match_set)))

    # Sample the weight vectors - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.sample == 100.0):
      use_w_vec_dict = w_vec_dict

    else:
      num_w_vec_sample = max(2, int(len(w_vec_dict)*self.sample/100.0))

      use_w_vec_dict = {}  # Create a new weight vector dictionary with samples

      rec_id_tuple_sample = random.sample(w_vec_dict.keys(),num_w_vec_sample)
      assert len(rec_id_tuple_sample) == num_w_vec_sample

      for rec_id_tuple in rec_id_tuple_sample:
        use_w_vec_dict[rec_id_tuple] = w_vec_dict[rec_id_tuple]

    logging.info('  Number of weight vectors to be used for SVM ' + \
                 'classification: %d' % (len(use_w_vec_dict)))

    train_data =   []
    train_labels = []

    for (rec_id_tuple, w_vec) in use_w_vec_dict.iteritems():
      train_data.append(w_vec)
      if (rec_id_tuple in match_set):
        train_labels.append(1.0)  # Match class
      else:
        train_labels.append(-1.0)  # Non-match class

    assert len(train_data) == len(train_labels)
    assert len(train_data) == len(use_w_vec_dict)

    # Initialise and train the SVM - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.kernel_type == 'LINEAR'):
      svm_kernel = svm.LINEAR
    elif (self.kernel_type == 'POLY'):
      svm_kernel = svm.POLY
    elif (self.kernel_type == 'RBF'):
      svm_kernel = svm.RBF
    elif (self.kernel_type == 'SIGMOID'):
      svm_kernel = svm.SIGMOID

    svm_prob =  svm.svm_problem(train_labels, train_data)

    # Due to change in SVM parameter setting in svm module, we need to catch
    # possible error
    #
    try:
      svm_param = svm.svm_parameter(svm_type = svm.C_SVC, C=self.C,
                                    kernel_type=svm_kernel)
      self.svm_model = svm.svm_model(svm_prob, svm_param)
      self.svm_version = 'old'

    except:
      svm_param = svm.svm_parameter('-s %d -c %f -t %d' % \
                  (svm.C_SVC, self.C, svm_kernel))
      self.svm_model = svm.libsvm.svm_train(svm_prob, svm_param)
      self.svm_version = 'new'

    logging.info('Trained SVM with %d training examples' % \
                 (len(use_w_vec_dict)))

  # ---------------------------------------------------------------------------

  def test(self, w_vec_dict, match_set, non_match_set):
    """Method to test a classifier using the given weight vector dictionary and
       match and non-match sets of record identifier pairs.

       Weight vectors will be assigned to matches or non-matches according to
       the SVM classification. No weight vector will be assigned to the
       possible match set.

       Will return a confusion matrix as a list of the form: [TP, FN, FP, TN].
    """

    if (self.svm_model == None):
      logging.warn('SVM has not been trained, testing not possible')
      return [0,0,0,0]

    svm_version = self.svm_version  # Shortcut

    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)
    auxiliary.check_is_set('match_set', match_set)
    auxiliary.check_is_set('non_match_set', non_match_set)

    # Check that match and non-match sets are separate and do cover all weight
    # vectors given
    #
    if (len(match_set.intersection(non_match_set)) > 0):
      logging.exception('Intersection of match and non-match set not empty')
      raise Exception
    if ((len(match_set)+len(non_match_set)) != len(w_vec_dict)):
      logging.exception('Weight vector dictionary of different length than' + \
                        ' summed lengths of match and non-match sets: ' + \
                        '%d / %d+%d=%d' % (len(w_vec_dict), len(match_set),
                        len(non_match_set), len(match_set)+len(non_match_set)))
      raise Exception

    num_true_m =   0
    num_false_m =  0
    num_true_nm =  0
    num_false_nm = 0

    for (rec_id_tuple, w_vec) in w_vec_dict.iteritems():

      if (svm_version == 'old'):
        if (self.svm_model.predict(w_vec) == 1.0):  # Match prediction
          pred_match = True
        else:
          pred_match = False

      else:  # New SVM module version
        x0, max_idx = svm.gen_svm_nodearray(w_vec)

        if (svm.libsvm.svm_predict(self.svm_model, x0) == 1.0):  # Match
          pred_match = True
        else:
          pred_match = False

      if (pred_match == True):
        if (rec_id_tuple in match_set):
          num_true_m += 1
        else:
          num_false_m += 1
      else:  # Non-match prediction
        if (rec_id_tuple in non_match_set):
          num_true_nm += 1
        else:
          num_false_nm += 1

    assert (num_true_m+num_false_nm+num_false_m+num_true_nm) == len(w_vec_dict)

    logging.info('  Results: TP = %d, FN = %d, FP = %d, TN = %d' % \
                 (num_true_m,num_false_nm,num_false_m,num_true_nm))

    return [num_true_m, num_false_nm, num_false_m, num_true_nm]

  # --------------------------------------------------------------------------

  def cross_validate(self, w_vec_dict, match_set, non_match_set, n=10):
    """Method to conduct a cross validation using the given weight vector
       dictionary and match and non-match sets of record identifier pairs.

       Will return a confusion matrix as a list of the form: [TP, FN, FP, TN].

       The cross validation approach splits the weight vector dictionary into
       'n' parts (and 'n' corresponding sub-set for matches and non-matches),
       and then generates 'n' SVM classifications, tests them and finally
       returns the average performance of these 'n' classifiers.
    """

    auxiliary.check_is_integer('n', n)
    auxiliary.check_is_positive('n', n)
    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)
    auxiliary.check_is_set('match_set', match_set)
    auxiliary.check_is_set('non_match_set', non_match_set)

    # Check that match and non-match sets are separate and do cover all weight
    # vectors given
    #
    if (len(match_set.intersection(non_match_set)) > 0):
      logging.exception('Intersection of match and non-match set not empty')
      raise Exception
    if ((len(match_set)+len(non_match_set)) != len(w_vec_dict)):
      logging.exception('Weight vector dictionary of different length than' + \
                        ' summed lengths of match and non-match sets: ' + \
                        '%d / %d+%d=%d' % (len(w_vec_dict), len(match_set),
                        len(non_match_set), len(match_set)+len(non_match_set)))
      raise Exception

    # Get a random vector dictionary element to get dimensionality of vectors
    #
    (rec_id_tuple, w_vec) = w_vec_dict.popitem()
    v_dim = len(w_vec)
    w_vec_dict[rec_id_tuple] = w_vec  # Put back in

    logging.info('')
    logging.info('Conduct %d-fold cross validation on SVM classifier ' % \
                 (n) + 'using %d weight vectors' % (len(w_vec_dict)))
    logging.info('  Match and non-match sets with %d and %d entries' % \
                 (len(match_set), len(non_match_set)))
    logging.info('  Dimensionality:   %d' % (v_dim))

    # Create the sub-sets of record identifier pairs for folds - - - - - - - -
    #
    rec_id_tuple_list = w_vec_dict.keys()
    random.shuffle(rec_id_tuple_list)
    fold_num_rec_id_tuple = max(1,int(round(float(len(rec_id_tuple_list))/n)))

    # Split the weight vector dictionary and match and non-match sets into
    # (lists containing one entry per fold) and only store test elements
    #
    w_vec_dict_test_list = []
    m_set_test_list =      []
    nm_set_test_list =     []

    for fold in range(n):
      w_vec_dict_test_list.append({})
      m_set_test_list.append(set())
      nm_set_test_list.append(set())

    for fold in range(n):

      # Calculate start and end indices for test elements for this fold
      #
      if (fold == (n-1)):  # The last fold, get remainder of list
        start = fold*fold_num_rec_id_tuple
        this_fold_test_ids = rec_id_tuple_list[start:]
      else:  # All other folds
        start = fold*fold_num_rec_id_tuple
        end = start+fold_num_rec_id_tuple
        this_fold_test_ids = rec_id_tuple_list[start:end]

      for rec_id_tuple in this_fold_test_ids:

        w_vec_dict_test_list[fold][rec_id_tuple] = w_vec_dict[rec_id_tuple]

        if (rec_id_tuple in match_set):
          m_set_test_list[fold].add(rec_id_tuple)
        else:
          nm_set_test_list[fold].add(rec_id_tuple)

      assert len(w_vec_dict_test_list[fold]) == len(this_fold_test_ids)
      assert len(m_set_test_list[fold]) + len(nm_set_test_list[fold]) == \
             len(this_fold_test_ids)

    # Initialise the total classification results - - - - - - - - - - - - - - -
    #
    num_true_m =   0
    num_false_nm = 0
    num_false_m =  0
    num_true_nm =  0

    # Loop over folds - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Generate training and test dictionaries and sets
    #
    for fold in range(n):  # First extract test record identifier pairs

      this_fold_test_m_set =       m_set_test_list[fold]
      this_fold_test_nm_set =      nm_set_test_list[fold]
      this_fold_test_w_vec_dict =  w_vec_dict_test_list[fold]

      this_fold_train_m_set =  match_set.difference(m_set_test_list[fold])
      this_fold_train_nm_set = non_match_set.difference(nm_set_test_list[fold])
      this_fold_train_w_vec_dict = {}
      for f2 in range(n):
        if (f2 != fold):
          this_fold_train_w_vec_dict.update(w_vec_dict_test_list[f2])

      assert len(this_fold_test_m_set) + len(this_fold_train_m_set) == \
             len(match_set)
      assert len(this_fold_test_nm_set) + len(this_fold_train_nm_set) == \
             len(non_match_set)
      assert len(this_fold_test_w_vec_dict) + \
             len(this_fold_train_w_vec_dict) == len(w_vec_dict)

      assert this_fold_test_m_set.intersection(this_fold_train_m_set) == set()
      assert this_fold_test_m_set.intersection(this_fold_test_nm_set) == set()
      assert this_fold_test_m_set.intersection(this_fold_train_nm_set) == set()
      assert this_fold_test_nm_set.intersection(this_fold_train_m_set) ==set()
      assert this_fold_test_nm_set.intersection(this_fold_train_nm_set) ==set()
      assert this_fold_train_m_set.intersection(this_fold_train_nm_set) ==set()

      # Train and test SVM classifier on this fold's data
      #
      self.train(this_fold_train_w_vec_dict, this_fold_train_m_set,
                 this_fold_train_nm_set)

      [this_num_true_m,this_num_false_nm,this_num_false_m,this_num_true_nm] = \
                                           self.test(this_fold_test_w_vec_dict,
                                                     this_fold_test_m_set,
                                                     this_fold_test_nm_set)
      num_true_m +=   this_num_true_m
      num_false_nm += this_num_false_nm
      num_false_m +=  this_num_false_m
      num_true_nm +=  this_num_true_nm

    # Calculate final cross validation results - - - - - - - - - - - - - - - -
    #
    num_true_m /=   float(n)
    num_false_nm /= float(n)
    num_false_m /=  float(n)
    num_true_nm /=  float(n)

    logging.info('  Results: TP = %d, FN = %d, FP = %d, TN = %d' % \
                 (num_true_m,num_false_nm,num_false_m,num_true_nm))

    return [num_true_m, num_false_nm, num_false_m, num_true_nm]

  # ---------------------------------------------------------------------------

  def classify(self, w_vec_dict):
    """Method to classify the given weight vector dictionary using the trained
       classifier.

       Will return three sets with record identifier pairs: 1) match set,
       2) non-match set, and 3) possible match set.

       The possible match set will be empty, as this classifier classifies all
       weight vectors as either matches or non-matches.
    """

    if (self.svm_model == None):
      logging.warn('SVM has not been trained, classification not possible')
      return set(), set(), set()

    svm_version = self.svm_version  # Shortcut

    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)

    match_set =      set()
    non_match_set =  set()
    poss_match_set = set()

    for (rec_id_tuple, w_vec) in w_vec_dict.iteritems():

      if (svm_version == 'old'):
        if (self.svm_model.predict(w_vec) == 1.0):  # Match prediction
          match_set.add(rec_id_tuple)
        else:  # Non-match prediction
          non_match_set.add(rec_id_tuple)

      else:  # New SVM module version
        x0, max_idx = svm.gen_svm_nodearray(w_vec)

        if (svm.libsvm.svm_predict(self.svm_model, x0) == 1.0):  # Match
          match_set.add(rec_id_tuple)
        else:  # Non-match prediction
          non_match_set.add(rec_id_tuple)

    assert (len(match_set) + len(non_match_set) + len(poss_match_set)) == \
            len(w_vec_dict)

    logging.info('Classified %d weight vectors: %d as matches, %d as ' % \
                 (len(w_vec_dict), len(match_set), len(non_match_set)) + \
                 'non-matches, and %d as possible matches' % \
                 (len(poss_match_set)))

    return match_set, non_match_set, poss_match_set

# =============================================================================

class TwoStep(Classifier):
  """Implements a two-step classifier that in a first step extracts likely
     true matches and true non-matches from the weight vector dictionary, and
     in a second step the trains a binary classifier using these training
     examples.

     For more information on this approach please see:
       Christen, P: A Two-Step Classification Approach to Unsupervised Record
       Linkage, Australasian Data Mining Conference (AusDM'07), Gold Coast,
       Conferences in Research and Practice in Information Technology (CRPIT),
       vol. 70. To be published.

     The arguments that have to be set when this classifier is initialised are:

       s1_match_method      The method to be used in step one to extract the
                            match training examples. This has to be a tuple
                            containing three elements:
                            1) A numerical value which is the match values (for
                               each element in the weight vectors), normally
                               this is 1.0. If this is set to None, then 1.0
                               will be used.
                            2) The way examples are selected, can be either
                               'threshold' or 'nearest'.
                            3) A corresponding parameter; for 'threshold' a
                               number that gives the maximum difference between
                               the match value and what is included in the
                               match example set, for 'nearest' a positive
                               number giving the number of nearest vectors that
                               will be included in the match set.
                            4) Only for 'nearest', a flag 'unique' which if set
                               to True means the nearest unique weight vectors
                               will be selected, if set to False simply the
                               nearest weight vectors will be selected (if if
                               many of them are the same).
       s1_non_match_method  A similar three-element tuple used to determined
                            how the non-match training examples will be
                            extracted from the weight vector dictionary.
       random_selection     The method used to add random weight vectors to the
                            match and non-match training examples sets. This
                            has to be either None (so no randomly selected
                            weight vectors will be added), or a tuple made of
                            three elements, with the first element being the
                            random selection method (see possible methods
                            below), the second element the percentage of so far
                            not-assigned weight vectors to be added to the
                            match training example set, and the third element
                            the percentage of so far not-assigned weight
                            vectors to be added to the non-match training
                            example set. The sum of these two percentage values
                            can be maximum 80% (due to the random selection
                            process). Possible random selection methods are:
                            - 'uniform'      Uniform random distribution.
                            - 'linear'       Linear random distribution.
                            - 'exponential'  Exponential random distribution.
                            Note that the default values for this is None, so
                            no random insertion of weight vectors is done.
       s2_classifier        The classifier to be used for step two. This
                            argument has to be a tuple with the first element
                            being the name of the classifier and the following
                            elements parameters for the classifier. Possible
                            are currently:
                            - ('kmeans', dist_measure, max_iter_count)
                              A K-means clustering will be done on the match
                              and non-match training example sets, i.e. two
                              centroids will be calculated in the training
                              process, to be used for testing and
                              classification later. The parameter dist_measure
                              has to be a function that calculates a distance
                              between two vectors (see the Febrl mymath.py
                              module for such functions). If the value of
                              'max_iter_count' is set to 0, then the centroids
                              will be calculated only using the training sets,
                              but no iterations using all weight vectors will
                              be done.
                            - ('nn', dist_measure, k)
                              This classifier uses a nearest neighbour approach
                              (with k the number of nearest weight vectors to
                              be considered, this has to be an odd positive
                              integer number, e.g. 1,3,5 etc.) using the given
                              distance measure (which has to be a function that
                              calculates a distance between two vectors, please
                              see the Febrl mymath.py module for such
                              functions).
                            - ('svm', kernel_type, C, increment, train_perc)
                              A SVM will be trained using the match and
                              non-match and non-match training example sets.
                              See the SuppVecMachine documentation for more
                              information on the parameters. 'increment' is a
                              percentage number that will determine the
                              incremental inclusion of additional training
                              examples from the weight vectors not included in
                              the training sets from step 1 for refinement of
                              the training step (calculated as percentage of
                              the number of weight vectors classified as
                              matches or non-matches in an iteration). If
                              'increment' is set to 0, no additional weight
                              vectors will be included. The 'train_perc'
                              argument then gives the total percentage of
                              weight vectors to be included into the training
                              sets. For example, if set to 100 then all weight
                              vectors will be included into a training set,
                              while if set to 50 only half will be included.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the classifier specific arguments first, then call
       call the base class constructor.
    """

    self.s1_m_method =   None
    self.s1_nm_method =  None
    self.rand_sel =      None
    self.s2_classifier = None

    self.svm_model =   None  # Will be set in training method
    self.m_centroid =  None
    self.nm_centroid = None

    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('s1_match')):
        auxiliary.check_is_tuple('s1_match_method', value)
        if (len(value) not in [3,4]):
          logging.exception('Match method tuple is not of length 3 or 4: %s' \
                            % (str(value)))
          raise Exception
        self.s1_m_method = value

      elif (keyword.startswith('s1_non_m')):
        auxiliary.check_is_tuple('s1_non_match_method', value)
        if (len(value) not in [3,4]):
          logging.exception('Non-match method tuple is not of length 3 or ' + \
                            '4: %s' % (str(value)))
          raise Exception
        self.s1_nm_method = value

      elif (keyword.startswith('random')):
        if (value != None):
          auxiliary.check_is_tuple('random_selection', value)
          if (len(value) != 3):
            logging.exception('Random selection tuple is not of length 3:' + \
                              ' %s' % (str(value)))
            raise Exception
          self.rand_sel = value

      elif (keyword.startswith('s2_class')):
        auxiliary.check_is_tuple('s2_classifier', value)
        if (len(value) not in [3,5]):
          logging.exception('Value of "s2_classifier" is not a tuple with ' + \
                              'three or five elements: %s' % (str(value)))
          raise Exception
        if (value[0] not in ['kmeans', 'nn', 'svm']):
          logging.exception('Value of "s2_classifier[0]" is not one of ' + \
                            '"kmeans", "nn", or "svm": %s' % (str(value[0])))
          raise Exception
        self.s2_classifier = value

      else:
        base_kwargs[keyword] = value

    Classifier.__init__(self, base_kwargs)  # Initialise base class

    # Check the match and non-match method tuples - - - - - - - - - - - - - - -
    #
    if (self.s1_m_method[0] == None):
      self.s1_m_method[0] = 1.0  # Set to a default value
    if (self.s1_nm_method[0] == None):
      self.s1_nm_method[0] = 1.0
    auxiliary.check_is_number('s1_match_method[0]',     self.s1_m_method[0])
    auxiliary.check_is_number('s1_non_match_method[0]', self.s1_nm_method[0])
    auxiliary.check_is_string('s1_match_method[1]',     self.s1_m_method[1])
    auxiliary.check_is_string('s1_non_match_method[1]', self.s1_nm_method[1])
    if (self.s1_m_method[1] not in ['threshold','nearest']):
      logging.exception('Value of "s1_match_method[1]" is neither ' + \
                        '"threshold" nor "nearest"')
      raise Exception
    if (self.s1_nm_method[1] not in ['threshold','nearest']):
      logging.exception('Value of "s1_non_match_method[1]" is neither ' + \
                        '"threshold" nor "nearest"')
      raise Exception
    auxiliary.check_is_number('s1_match_method[2]',       self.s1_m_method[2])
    auxiliary.check_is_number('s1_non_match_method[2]',   self.s1_nm_method[2])
    auxiliary.check_is_positive('s1_match_method[2]',     self.s1_m_method[2])
    auxiliary.check_is_positive('s1_non_match_method[2]', self.s1_nm_method[2])
    if (self.s1_m_method[1] == 'nearest'):
      auxiliary.check_is_integer('s1_match_method[2] for "nearest"',
                                 self.s1_m_method[2])
      if (len(self.s1_m_method) != 4):
        logging.exception('Value of "s1_match_method" tuple with "nearest"' + \
                          ' must contain four elements')
        raise Exception
      auxiliary.check_is_flag('s1_match_method[3]', self.s1_m_method[3])
    if (self.s1_nm_method[1] == 'nearest'):
      auxiliary.check_is_integer('s1_non_match_method[2] for "nearest"',
                                 self.s1_nm_method[2])
      if (len(self.s1_nm_method) != 4):
        logging.exception('Value of "s1_non_match_method" tuple with ' + \
                          '"nearest" must contain four elements')
        raise Exception
      auxiliary.check_is_flag('s1_non_match_method[3]', self.s1_nm_method[3])

    # If set, check random selection method - - - - - - - - - - - - - - - - - -
    #
    if (self.rand_sel != None):
      auxiliary.check_is_string('random_selection[0]', self.rand_sel[0])
      if (self.rand_sel[0] not in ['uniform', 'linear','exponential']):
        logging.exception('Value of "random_selection[0]" is not one of ' + \
                          '"uniform", "linear" or "exponential": %s' % (value))
        raise Exception
      auxiliary.check_is_percentage('random_selection[1]', self.rand_sel[1])
      auxiliary.check_is_percentage('random_selection[2]', self.rand_sel[2])
      auxiliary.check_is_not_negative('random_selection[1]', self.rand_sel[1])
      auxiliary.check_is_not_negative('random_selection[2]', self.rand_sel[2])

      if (float(self.rand_sel[1]) + float(self.rand_sel[2]) > 80):
        logging.exception('Sum of random selection percentages must be ' + \
                          'less than or equal to 80%')
        raise Exception

    # Make sure a classifier is defined correctly - - - - - - - - - - - - - - -
    #
    auxiliary.check_is_tuple('s2_classifier', self.s2_classifier)

    if (self.s2_classifier[0] == 'kmeans'):
      if (len(self.s2_classifier) != 3):
        logging.exception('Step 2 classifer "kmeans" reguires two ' + \
                          'parameters (distance measure and maximum ' + \
                          'iteration count): %s' % (str(self.s2_classifier)))
        raise Exception
      auxiliary.check_is_function_or_method('dist_measure',
                                            self.s2_classifier[1])
      auxiliary.check_is_integer('max_iter_count', self.s2_classifier[2])
      auxiliary.check_is_not_negative('max_iter_count', self.s2_classifier[2])

    elif (self.s2_classifier[0] == 'nn'):
      if (len(self.s2_classifier) != 3):
        logging.exception('Step 2 classifer "nn" reguires two parameters ' + \
                          '(distance measure and number of nearest ' + \
                          'neighbours, k): %s' % (str(self.s2_classifier)))
        raise Exception
      auxiliary.check_is_function_or_method('dist_measure',
                                            self.s2_classifier[1])
      auxiliary.check_is_integer('k', self.s2_classifier[2])
      auxiliary.check_is_positive('k', self.s2_classifier[2])
      if ((self.s2_classifier[2] % 2) == 0):  # An even number
        logging.exception('Number of nearest neighbours "k" for step 2 ' + \
                          'classifier "nn" has to be a positive odd integer')
        raise Exception

    elif (self.s2_classifier[0] == 'svm'):
      if (len(self.s2_classifier) != 5):
        logging.exception('Step 2 classifer "svm" reguires four parameters' + \
                          ' (kernel type, C, increment, and training ' + \
                          'percentage): %s' % (str(self.s2_classifier)))
        raise Exception
      if (self.s2_classifier[1] not in ['LINEAR', 'POLY', 'RBF', 'SIGMOID']):
        logging.exception('Illegal value for kernel type: %s ' % \
                          (self.s2_classifier[1]) + \
                          '(possible are: LINEAR, POLY, RBF, SIGMOID)')
        raise Exception
      auxiliary.check_is_number('C', self.s2_classifier[2])
      auxiliary.check_is_not_negative('C', self.s2_classifier[2])
      auxiliary.check_is_integer('increment', self.s2_classifier[3])
      auxiliary.check_is_not_negative('increment', self.s2_classifier[3])
      auxiliary.check_is_integer('train_perc', self.s2_classifier[4])
      auxiliary.check_is_not_negative('train_perc', self.s2_classifier[4])

      if ((self.s2_classifier[3] > 0) and (self.s2_classifier[4] == 0)):
        logging.exception('Two-step SVM increment larger than 0 but ' + \
                          'training percentage set to 0 - not possible.')
        raise Exception
      if ((self.s2_classifier[3] == 0) and (self.s2_classifier[4] > 0)):
        logging.warning('Two-step SVM increment set to 0 but training ' + \
                        'percentage set to larger than zero - Not used.')

    # Check if step 2 classifier is SVM and svm.py module is installed or not
    #
    if ((self.s2_classifier == 'svm') and (imp_svm == False)):
      logging.exception('Module "svm.py" not installed, cannot use "svm" ' + \
                        'classifier in step 2')
      raise Exception

    self.log([('Step 1 match method',     str(self.s1_m_method)),
              ('Step 1 non-match method', str(self.s1_nm_method)),
              ('Random selection',        str(self.rand_sel)),
              ('Step 2 classifier',  str(self.s2_classifier))]) # Log a message

    # If the weight vector dictionary and both match and non-match sets - - - -
    # are given start the training process
    #
    if ((self.train_w_vec_dict != None) and (self.train_match_set != None) \
        and (self.train_non_match_set != None)):
      self.train(self.train_w_vec_dict, self.train_match_set,
                 (self.train_non_match_set))

  # ---------------------------------------------------------------------------

  def train(self, w_vec_dict, match_set, non_match_set):
    """Method to train a classifier using the given weight vector dictionary.
       Note that the given match and non-match sets of record identifier pairs
       will not be used (unsupervised training).

       This method extracts training weight vectors that are with high
       likelihood true matches and true non-matches, and then trains a
       classifier on them
    """

    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)

    self.train_w_vec_dict =    w_vec_dict  # Save
    self.train_match_set =     match_set
    self.train_non_match_set = non_match_set

    # Get a random vector to get dimensionality of vectors
    #
    (rec_id_tuple, w_vec) = w_vec_dict.popitem()
    v_dim = len(w_vec)
    w_vec_dict[rec_id_tuple] = w_vec  # Put back in

    logging.info('Train two-step classifier using %d weight vectors' % \
                 (len(w_vec_dict)))
    logging.info('  Match and non-match sets with %d and %d entries' % \
                 (len(match_set), len(non_match_set)))

    m_train_set =  set()  # Match and non-match training set to be extracted
    nm_train_set = set()

    # Define some short hands for faster access
    #
    m_val = self.s1_m_method[0]
    if (self.s1_m_method[1] == 'threshold'):
      m_do_threshold = True
      m_threshold =    self.s1_m_method[2]
    else:
      m_do_threshold = False  # Do nearest method
      m_nearest =      self.s1_m_method[2]
      m_unique =       self.s1_m_method[3]
      m_heap =         []
    nm_val = self.s1_nm_method[0]
    if (self.s1_nm_method[1] == 'threshold'):
      nm_do_threshold = True
      nm_threshold =    self.s1_nm_method[2]
    else:
      nm_do_threshold = False  # Do nearest method
      nm_nearest =      self.s1_nm_method[2]
      nm_unique =       self.s1_nm_method[3]
      nm_heap =         []

    # Main loop over all vectors in the weight vector dictionary - - - - - - -
    #
    for (rec_id_tuple, w_vec) in w_vec_dict.iteritems():

      # Check if this weight vector is a good match example - - - - - - - - - -
      #
      if (m_do_threshold == True):
        w_vec_ok = True
        for w in w_vec:  # Check if all weights are within threshold
          if (abs(w-m_val) > m_threshold):
            w_vec_ok = False
            break
        if (w_vec_ok == True):
          m_train_set.add(rec_id_tuple)

      else:
        vec_sum = 0.0
        for w in w_vec:
#          vec_sum += abs(w-m_val)         # Manhatten distance
          vec_sum += (w-m_val)*(w-m_val)  # Euclidean distance
        heapq.heappush(m_heap, (vec_sum, rec_id_tuple))

      # Check if this weight vector is a good non-match example - - - - - - - -
      #
      if (nm_do_threshold == True):
        w_vec_ok = True
        for w in w_vec:  # Check if all weights are within threshold
          if (abs(w-nm_val) > nm_threshold):
            w_vec_ok = False
            break
        if (w_vec_ok == True):
          nm_train_set.add(rec_id_tuple)

      else:
        vec_sum = 0.0
        for w in w_vec:
#          vec_sum += abs(w-nm_val)          # Manhatten distance
          vec_sum += (w-nm_val)*(w-nm_val)  # Euclidean distance
        heapq.heappush(nm_heap, (vec_sum, rec_id_tuple))

    # Nothing further to process for threshold methods

    # Process heaps for 'nearest' method - - - - - - - - - - - - - - - - - - -
    #
    if (self.s1_m_method[1] == 'nearest'):
      prev_tuple = -1  # Keep previous record tuple for checking uniqueness

      while ((len(m_train_set) < m_nearest) and (len(m_heap) > 0)):
        this_tuple = heapq.heappop(m_heap)[1]  # Get rec-IDs from first tuple

        if ((m_unique == True) and (prev_tuple != -1)):
          if (w_vec_dict[prev_tuple] != w_vec_dict[this_tuple]):
            m_train_set.add(this_tuple)  # A different weight vector, add tuple
            prev_tuple = this_tuple
        else:  # Don't check for uniquness, simple add
          m_train_set.add(this_tuple) # Add record identifier tuple
          prev_tuple = this_tuple

    if (self.s1_nm_method[1] == 'nearest'):
      prev_tuple = -1  # Keep previous record tuple for checking uniqueness

      while ((len(nm_train_set) < nm_nearest) and (len(nm_heap) > 0)):
        this_tuple = heapq.heappop(nm_heap)[1]

        if ((nm_unique == True) and (prev_tuple != -1)):
          if (w_vec_dict[prev_tuple] != w_vec_dict[this_tuple]):
            nm_train_set.add(this_tuple) # A different weight vector, add tuple
            prev_tuple = this_tuple
        else:  # Don't check for uniquness, simple add
          nm_train_set.add(this_tuple) # Add record identifier tuple
          prev_tuple = this_tuple

    logging.info('Set %d weight vectors as matches and %d as non-matches' % \
                 (len(m_train_set), len(nm_train_set)))
    logging.info('  %d weight vectors not included into set of examples' % \
                 (len(w_vec_dict)-len(m_train_set)-len(nm_train_set)))

    # Have to check if a weight vector is in both match and non-match set - - -
    #
    remove_tuple_set = set()

    for rec_id_tuple in m_train_set:
      if (rec_id_tuple in nm_train_set):  # Remove it from both sets
        remove_tuple_set.add(rec_id_tuple)

    if (len(remove_tuple_set) > 0):
      for rec_id_tuple in remove_tuple_set:
        m_train_set.remove(rec_id_tuple)
        nm_train_set.remove(rec_id_tuple)

      logging.warn('Removed %d weight vectors from both match and ' % \
                   (len(remove_tuple_set)) + 'non-match training sets')

    # Random selection of additional weight vectors into training example sets
    #
    if (self.rand_sel != None):
      rand_method = self.rand_sel[0]
      num_m_perc =  self.rand_sel[1]
      num_nm_perc = self.rand_sel[2]

      # Calculate how many weight vectors have not yet been assigend to
      # training sets
      #
      num_w_vec_left = len(w_vec_dict)-(len(m_train_set)+len(nm_train_set))

      # At least one random example
      #
      num_m_examples =  max(1, int(num_w_vec_left * num_m_perc/100.0))
      num_nm_examples = max(1, int(num_w_vec_left * num_nm_perc/100.0))

      # Check if enough weight vectors still unassigned for random selection
      #
      if ((num_w_vec_left - (num_m_examples+num_nm_examples)) < 0):

        logging.warn('Cannot select %d match and %d non-match examples, ' % \
                     (num_m_examples, num_nm_examples) + \
                     'not enough not-assigned weight vectors left: %d' % \
                     (num_w_vec_left))

        # Reduce the number of examples to be selected
        #
        while ((num_w_vec_left - (num_m_examples + num_nm_examples)) <= 0):
          if (num_m_examples > 0):
            num_m_examples -=  1
          if (num_nm_examples > 0):
            num_nm_examples -= 1

        logging.warn('  Ajusted number of random match examples to %d and' % \
                     (num_m_examples) + ' number of non-match examples to' + \
                     ' %d' % (num_nm_examples))

      assert (num_w_vec_left - (num_m_examples + num_nm_examples)) >= 0

      # Only insert random examples if possible - - - - - - - - - - - - - - - -
      #
      if ((num_m_examples > 0) or (num_nm_examples > 0)):

        logging.info('Randomly select %d match and %d non-match examples ' % \
                     (num_m_examples, num_nm_examples) + 'using random ' + \
                     'method: %s' % (rand_method))

        # Calculate distance from 0.0 for all weight vectors - - - - - - - - -
        # (i.e. assumes 0 is the non-match value)
        #
        w_vec_dist_dict = {}  # Distance as keys with lists of record identifer
                              # pairs as values
        w_vec_m_min_dist =  9999.99  # Also get minimum and maximum distances
        w_vec_m_max_dist =  -999.99  # for elements so far in match and
        w_vec_nm_min_dist = 9999.99  # non-match training sets
        w_vec_nm_max_dist = -999.99

        for (rec_id_tuple, w_vec) in w_vec_dict.iteritems():
          w_vec_sum = sum(w_vec)

          if ((rec_id_tuple not in m_train_set) and \
              (rec_id_tuple not in nm_train_set)):
            w_vec_list = w_vec_dist_dict.get(w_vec_sum, [])
            w_vec_list.append(rec_id_tuple)
            w_vec_dist_dict[w_vec_sum] = w_vec_list

          elif (rec_id_tuple in match_set):
            w_vec_m_min_dist = min(w_vec_m_min_dist, w_vec_sum)
            w_vec_m_max_dist = max(w_vec_m_max_dist, w_vec_sum)
          elif (rec_id_tuple in non_match_set):
            w_vec_nm_min_dist = min(w_vec_nm_min_dist, w_vec_sum)
            w_vec_nm_max_dist = max(w_vec_nm_max_dist, w_vec_sum)

        w_vec_dist_list = w_vec_dist_dict.keys()
        w_vec_dist_list.sort()  # Sort according to distance from 0.0

        logging.info('  Minimum and maximum distance from 0.0:           ' + \
                     '%.3f / %.3f' % (w_vec_dist_list[0], w_vec_dist_list[-1]))
        logging.info('    Minimum and maximum distance in match set:     ' + \
                     '%.3f / %.3f' % (w_vec_m_min_dist, w_vec_m_max_dist))
        logging.info('    Minimum and maximum distance in non-match set: ' + \
                     '%.3f / %.3f' % (w_vec_nm_min_dist, w_vec_nm_max_dist))

        num_w_vec_dist = len(w_vec_dist_list) # Number of unique distances

        num_w_vec_not_sel = 0
        for dist in w_vec_dist_list:
          num_w_vec_not_sel += len(w_vec_dist_dict[dist])

        logging.info('  Number of different distances: %d' % (num_w_vec_dist))
        logging.info('  %d weight vectors have so far not been selected' % \
                     (num_w_vec_not_sel))

        # Randomly select match and non-match examples - - - - - - - - - - - -
        #
        random_m_set =  set()
        random_nm_set = set()

        while ((len(random_m_set) < num_m_examples) or \
               (len(random_nm_set) < num_nm_examples)):

          if (len(random_m_set) < num_m_examples):  # Select a match example

            if (rand_method[:3] == 'uni'):
              r_ind = random.randint(0, num_w_vec_dist-1)
            elif (rand_method[:3] == 'lin'):
              r_ind = int(mymath.random_linear(num_w_vec_dist))
            else:  # Exponential
              r_ind = num_w_vec_dist-int(mymath.random_expo(num_w_vec_dist))-1

            assert (0 <= r_ind) and (r_ind < num_w_vec_dist)

            # Get list of record identifiers at this distance
            #
            w_vec_list = w_vec_dist_dict[w_vec_dist_list[r_ind]]

            # Check if there are still record identifiers available
            #
            if (w_vec_list != []):  # If so select one at random
              rec_id = random.choice(w_vec_list)
              random_m_set.add(rec_id)  # Add record identifiers to match set

              # Remove the chosen record identifier tuple and put it back
              #
              w_vec_list.remove(rec_id)
              w_vec_dist_dict[w_vec_dist_list[r_ind]] = w_vec_list

          if (len(random_nm_set) < num_nm_examples):  # Select a non-match ex.

            if (rand_method[:3] == 'uni'):
              r_ind = random.randint(0, num_w_vec_dist-1)
            elif (rand_method[:3] == 'lin'):
              r_ind = num_w_vec_dist - \
                                    int(mymath.random_linear(num_w_vec_dist))-1
            else:  # Exponential
              r_ind = int(mymath.random_expo(num_w_vec_dist))

            assert (0 <= r_ind) and (r_ind < num_w_vec_dist)

            # Get list of record identifiers at this distance
            #
            w_vec_list = w_vec_dist_dict[w_vec_dist_list[r_ind]]

            # Check if there are still record identifiers available
            #
            if (w_vec_list != []):
              rec_id = random.choice(w_vec_list)
              random_nm_set.add(rec_id) # Add record identifiers to non-matches

              # Remove the chosen record identifier tuple and put it back
              #
              w_vec_list.remove(rec_id)
              w_vec_dist_dict[w_vec_dist_list[r_ind]] = w_vec_list

        assert len(random_m_set) ==  num_m_examples
        assert len(random_nm_set) == num_nm_examples

        m_train_set =  m_train_set.union(random_m_set)  # Merge training sets
        nm_train_set = nm_train_set.union(random_nm_set)

    if ((len(m_train_set) == 0) or (len(nm_train_set) == 0)):
      logging.warn('At least one of the two training sets is empty: ' + \
                   'Number of matches: %d and non-matches: %d' % \
                   (len(m_train_set), len(nm_train_set)))

      if (self.s2_classifier[0] == 'svm'):  # Cannot train SVM classifier
        self.svm_model = None
      elif (self.s2_classifier[0] == 'kmeans'):
        self.m_centroid = None
        self.nm_centroid = None
      elif (self.s2_classifier[0] == 'nn'):
        self.nn_m_train_w_vec_set =  None
        self.nn_nm_train_w_vec_set = None

      return  # Return without model training

    self.m_train_set =  m_train_set  # Save for later
    self.nm_train_set = nm_train_set

    # End of step one, now classify training example sets - - - - - - - - - - -
    #
    logging.info('Start step two classification using method: %s ' % \
                 (str(self.s2_classifier)) + 'with %d match and %d non-match' \
                 % (len(m_train_set), len(nm_train_set)) +' training examples')

    if (self.s2_classifier[0] == 'svm'):  # - - - - - - - - - - - - - - - - - -

      svm_type =   svm.C_SVC
      if (self.s2_classifier[1] == 'LINEAR'):
        svm_kernel = svm.LINEAR
      elif (self.s2_classifier[1] == 'POLY'):
        svm_kernel = svm.POLY
      elif (self.s2_classifier[1] == 'RBF'):
        svm_kernel = svm.RBF
      elif (self.s2_classifier[1] == 'SIGMOID'):
        svm_kernel = svm.SIGMOID
      C =          self.s2_classifier[2]
      increment =  self.s2_classifier[3]
      train_perc = self.s2_classifier[4]

      train_data =   []  # Generate training data
      train_labels = []

      for rec_id_tuple in m_train_set:
        train_data.append(w_vec_dict[rec_id_tuple])
        train_labels.append(1.0)  # Match class
      for rec_id_tuple in nm_train_set:
        train_data.append(w_vec_dict[rec_id_tuple])
        train_labels.append(-1.0)  # Match class

      # Initialise and train the SVM - - - - - - - - - - - - - - - - - - - - -
      #
      svm_prob =  svm.svm_problem(train_labels, train_data)

      # Due to change in SVM parameter setting in svm module, we need to catch
      # possible error
      #
      try:
        svm_param = svm.svm_parameter(svm_type = svm.C_SVC, C=C,
                                      kernel_type=svm_kernel)
        self.svm_model = svm.svm_model(svm_prob, svm_param)
        self.svm_version = 'old'

      except:
        svm_param = svm.svm_parameter('-s %d -c %f -t %d' % \
                    (svm_type, C, svm_kernel))
        self.svm_model = svm.libsvm.svm_train(svm_prob, svm_param)
        self.svm_version = 'new'

      # Iterative refinement by inclusion of additional weight vectors
      #
      if (increment > 0):

        #print '-------------------------------------------------'
        #print 'Initial training sets size:',len(m_train_set),len(nm_train_set)

        # Total number of weight vectors to add to both training sets
        #
        total_train_num_w_vec = int(len(w_vec_dict) * train_perc / 100.0)

        # Make a dictionary of weight vectors not yet in training data
        #
        un_used_w_vec_dict = {}

        for (rec_id_tuple, w_vec) in w_vec_dict.iteritems():
          if ((rec_id_tuple not in m_train_set) and \
              (rec_id_tuple not in nm_train_set)):
            un_used_w_vec_dict[rec_id_tuple] = w_vec

        assert len(train_data) == len(m_train_set)+len(nm_train_set), \
               (len(train_data), len(m_train_set)+len(nm_train_set))
        assert len(w_vec_dict) == len(un_used_w_vec_dict)+len(train_data), \
               (len(w_vec_dict), len(un_used_w_vec_dict)+len(train_data))

        # Loop until enough weight vectors in training data
        #
        while (len(train_data) < total_train_num_w_vec):

          # Two lists with the classified weight vectors not yet in the
          # training sets
          #
          new_m_class_set_list =  []
          new_nm_class_set_list = []

          # Classify so far un-used weight vectors
          #
          for (rec_id_tuple, w_vec) in un_used_w_vec_dict.iteritems():

            if (self.svm_version == 'old'):
              c1 = self.svm_model.predict(w_vec)
              c2 = self.svm_model.predict_values(w_vec)

              if (c1 == 1.0):  # Classified as match
                new_m_class_set_list.append((c2[(1, -1)], rec_id_tuple))
              else:  # A non-match
                new_nm_class_set_list.append((c2[(-1, 1)], rec_id_tuple))

            else:  # new svm module version
              x0, max_idx = svm.gen_svm_nodearray(w_vec)
              c1 = svm.libsvm.svm_predict(self.svm_model, x0)
              # From: http://code.google.com/p/search/source/browse/...
              #       trunk/libsvm-2.91/python/svmutil.py?r=5
              nr_class = 2 #svm.toPyModel(self.svm_model).get_nr_class()
              nr_classifier = nr_class*(nr_class-1)//2
              dec_values = (svm.c_double * nr_classifier)()
              c2 = svm.libsvm.svm_predict_values(self.svm_model, x0, dec_values)
              assert c1==c2
              values = dec_values[:nr_classifier]

              if (c1 == 1.0):  # Classified as match
                new_m_class_set_list.append((values[0], rec_id_tuple))
              else:  # A non-match
                new_nm_class_set_list.append((values[0], rec_id_tuple))

          # If any of the two lists are empty (i.e. all weight vectors were
          # either classified as matches or non-matches) then exit the loop
          #
          if ((new_m_class_set_list == []) or (new_nm_class_set_list == [])):
            logging.warn('Two-step SVM classifier: One of the training ' + \
                         'sets is empty! Abort iterative refinement.')
            break  # Abort and keep previous model  ### CHECk IF THIS IS OK??

          # Sort so largest probability values from SCM predict_values first
          #
          new_m_class_set_list.sort(reverse=True)
          new_nm_class_set_list.sort(reverse=True)

          #print 'New classified set sizes: matches=%d, non-matches=%d' % \
          #      (len(new_m_class_set_list), len(new_nm_class_set_list))

          #print
          #print 'Best match training examples:'
          #for i in range(min(10, len(new_m_class_set_list))):
          #  print ' ', new_m_class_set_list[i][0], \
          #        sum(w_vec_dict[new_m_class_set_list[i][1]]), \
          #        w_vec_dict[new_m_class_set_list[i][1]]
          #print
          #print 'Best non-match training examples:'
          #for i in range(min(10, len(new_nm_class_set_list))):
          #  print ' ', new_nm_class_set_list[i][0], \
          #        sum(w_vec_dict[new_nm_class_set_list[i][1]]), \
          #        w_vec_dict[new_nm_class_set_list[i][1]]
          #print

          # Calculate the number of new weight vectors to add to training sets
          #
          add_num_m =  max(1, int(len(new_m_class_set_list)*increment / 100.0))
          add_num_nm = max(1, int(len(new_nm_class_set_list)*increment /100.0))

          #print 'Add number of weight vectors: %d / %d' % \
          #      (add_num_m, add_num_nm)

          for i in range(add_num_m):  # Add new match training records
            rec_id_tuple = new_m_class_set_list[i][1]
            m_train_set.add(rec_id_tuple)
            train_data.append(w_vec_dict[rec_id_tuple])
            train_labels.append(1.0)
            del un_used_w_vec_dict[rec_id_tuple]  # It is used in training now

          for i in range(add_num_nm):  # Add new non-match training records
            rec_id_tuple = new_nm_class_set_list[i][1]
            nm_train_set.add(rec_id_tuple)
            train_data.append(w_vec_dict[rec_id_tuple])
            train_labels.append(-1.0)
            del un_used_w_vec_dict[rec_id_tuple]  # It is used in training now

          print 'Size of new training sets:',len(m_train_set),len(nm_train_set)

          assert len(train_data) == len(m_train_set)+len(nm_train_set), \
                 (len(train_data), len(m_train_set)+len(nm_train_set))
          assert len(w_vec_dict) == len(un_used_w_vec_dict)+len(train_data), \
                 (len(w_vec_dict), len(un_used_w_vec_dict)+len(train_data))

          self.m_train_set =  m_train_set  # Save for later use
          self.nm_train_set = nm_train_set

          # Re-train SVM classifier
          #
          svm_prob = svm.svm_problem(train_labels, train_data)
          if (self.svm_version == 'old'):
            self.svm_model = svm.svm_model(svm_prob, svm_param)
          else:
            self.svm_model = svm.libsvm.svm_train(svm_prob, svm_param)

        print 'Final training sets size:',len(m_train_set),len(nm_train_set)
        print (len(m_train_set)+len(nm_train_set)) / float(len(w_vec_dict))
        print '  increment:', increment, 'train_perc', train_perc
        print

      logging.info('Trained SVM with %d training examples' % (len(train_data)))

    # K-means clustering - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    elif (self.s2_classifier[0] == 'kmeans'):

      max_iter_count = self.s2_classifier[2]
      dist_meas =      self.s2_classifier[1]

      # Initialise centroids as center of training example sets from step 1
      #
      m_centroid =  [0.0]*v_dim
      nm_centroid = [0.0]*v_dim

      for rec_id_tuple in m_train_set:
        m_w_vec = w_vec_dict[rec_id_tuple]
        for i in range(v_dim):
          m_centroid[i] += m_w_vec[i]

      for rec_id_tuple in nm_train_set:
        nm_w_vec = w_vec_dict[rec_id_tuple]
        for i in range(v_dim):
          nm_centroid[i] += nm_w_vec[i]

      num_m =  len(m_train_set)
      num_nm = len(nm_train_set)

      for i in range(v_dim):  # Normalise centroid weights

        if (num_m > 0):
          m_centroid[i] /=  float(num_m)
        else:
          m_centroid[i] = 1.0  # Set to exact match value

        if (num_nm > 0):
          nm_centroid[i] /= float(num_nm)
        else:
          nm_centroid[i] = 0.0  # Set to total dissimilarity value

      # If max_iter_count > 0 start normal k-means iterations - - - - - - - - -
      #
      if (max_iter_count > 0):

        cluster_assign_dict = {}  # Dictionary with cluster assignments

        iter_cnt =    1  # Iteration counter
        num_changed = 1

        while (num_changed > 0) and (iter_cnt < max_iter_count):

          num_changed =     0
          new_m_centroid =  [0.0]*v_dim  # Summed new distances
          new_nm_centroid = [0.0]*v_dim

          # Step 1: Calculate cluster membership for all weight vectors (not
          #         just the ones in the training example sets
          #
          num_m =  0  # Number of weight vectors assigned to matches
          num_nm = 0  # Number of weight vectors assigned to non-matches

          for (rec_id_tuple, w_vec) in w_vec_dict.iteritems():

            m_dist =  dist_meas(w_vec, m_centroid)
            nm_dist = dist_meas(w_vec, nm_centroid)

            if (m_dist < nm_dist):  # Assign to cluster M (matches)
              old_assign = cluster_assign_dict.get(rec_id_tuple, 'X')
              if (old_assign != 'M'):
                num_changed += 1
              cluster_assign_dict[rec_id_tuple] = 'M'
              num_m += 1

              for i in range(v_dim):  # Add to summed cluster distances
                new_m_centroid[i] += w_vec[i]

            else:  # Assign to cluster NM (non-matches)
              old_assign = cluster_assign_dict.get(rec_id_tuple, 'X')
              if (old_assign != 'NM'):
                num_changed += 1
              cluster_assign_dict[rec_id_tuple] = 'NM'
              num_nm += 1

              for i in range(v_dim):  # Add to summed cluster distances
                new_nm_centroid[i] += w_vec[i]

          num_all = len(cluster_assign_dict)

          if ((num_m + num_nm) != num_all):
            logging.exception('Not all %d weight vectors were assigned: ' % \
                              (num_all) + 'M=%d, U=%d' % (num_m, num_nm))
            raise Exception

          if (num_m == 0) or (num_nm == 0):
            logging.warn('One cluster is empty: ' + \
                         'matches=%d, non-matches=%d' % (num_m, num_nm))
            break  # Stop K-means iterations

          # Step 2: Calculate new cluster centroids
          #
          for i in range(v_dim):  # Normalise new centroids
            new_m_centroid[i] /=  float(num_m)
            new_nm_centroid[i] /= float(num_nm)

          m_centroid =  new_m_centroid  # Assign new centroids
          nm_centroid = new_nm_centroid

          logging.info('Iteration %d: %d vectors changed cluster ' % \
                       (iter_cnt, num_changed) + 'assignment')

          iter_cnt += 1

      self.m_centroid =  m_centroid  # Save for later use
      self.nm_centroid = nm_centroid

      logging.info('Initial non-match centroid: %s' % \
                   (auxiliary.str_vector(nm_centroid)))
      logging.info('Initials match centroid:     %s' % \
                   (auxiliary.str_vector(m_centroid)))

    # Nearest neighbour based classification - - - - - - - - - - - - - - - - -
    #
    elif (self.s2_classifier[0] == 'nn'):

      dist_meas = self.s2_classifier[1]
      k =         self.s2_classifier[2]

      k1 = k+1  # Number of nearest weight vectors to store in nearest lists

      # Two dictionaries with the match and non-match training examples (weight
      # vectors as tuples), values will be lists with their nearest not yet
      # classified weight vectors
      #
      nn_m_train_w_vec_dict =  {}
      nn_nm_train_w_vec_dict = {}

      # First insert all step 1 training weight vectors into match and
      # non-match nearest neighbour training stes
      #
      for rec_id_tuple in m_train_set:
        m_w_vec = tuple(w_vec_dict[rec_id_tuple])
        nn_m_train_w_vec_dict[m_w_vec] = []

      for rec_id_tuple in nm_train_set:
        nm_w_vec = tuple(w_vec_dict[rec_id_tuple])
        nn_nm_train_w_vec_dict[nm_w_vec] = []

      # A dictionary which will contain the weight vectors not in the training
      # sets and information about their closest k training weight vectors
      # (distance and if a match ('M') or non-match ('NM'), and the training
      # vector itself)
      #
      nn_to_classify_w_vec_dict = {}

      # Generate a heap with all weight vectors so far not classified
      # according to their smallest distances to either a match or a non-match
      # training example (for each keep distance, weight vector and either 'M'
      # or 'NM')
      #
      nearest_w_vec_heap = []

      # Step 1: For each (unique) non-training weight vector find its closest k
      # weight vectors from the training sets
      #
      for (rec_id_tuple, w_vec) in w_vec_dict.iteritems():

        if ((rec_id_tuple not in m_train_set) and \
            (rec_id_tuple not in nm_train_set)):

          this_w_vec = tuple(w_vec)  # Tuple can be used as dictionary key

          # Check if this specific weight vector has not yet been processed
          # (i.e. duplicate vectors with same weights) or if it is in the
          # match or non-match training sets
          #
          if ((this_w_vec not in nn_to_classify_w_vec_dict) and \
              (this_w_vec not in nn_m_train_w_vec_dict) and \
              (this_w_vec not in nn_nm_train_w_vec_dict)):
            nearest_list = []  # List of k nearest weight vectors
            largest_dist = 9999999.999

            # Check its distance to weight vectors in the match training set
            #
            for m_w_vec in nn_m_train_w_vec_dict:
              m_dist = dist_meas(this_w_vec, m_w_vec)

              # Insert into nearest list for this weight vector
              #
              if (m_dist < largest_dist):
                nearest_list.append((m_dist, 'M', m_w_vec))
                nearest_list.sort()  # Smallest distances first
                nearest_list = nearest_list[:k1]  # Only keep k+1 nearest elem.
                largest_dist = nearest_list[-1][0]

              # Insert weight vector into nearest list for M weight vector
              #
              m_nearest_list = nn_m_train_w_vec_dict[m_w_vec]
              if (m_nearest_list == []):
                nn_m_train_w_vec_dict[m_w_vec] = [(m_dist, this_w_vec)]
              elif (m_dist < m_nearest_list[-1][0]):
                m_nearest_list.append((m_dist, this_w_vec))
                m_nearest_list.sort()  # Smallest distances first
                nn_m_train_w_vec_dict[m_w_vec] = m_nearest_list[:k1]

            for nm_w_vec in nn_nm_train_w_vec_dict:
              nm_dist = dist_meas(this_w_vec, nm_w_vec)

              # Insert into nearest list for this weight vector
              #
              if (nm_dist < largest_dist):
                nearest_list.append((nm_dist, 'NM', nm_w_vec))
                nearest_list.sort()  # Smallest distances first
                nearest_list = nearest_list[:k1]  # Only keep k+1 nearest elem.
                largest_dist = nearest_list[-1][0]

              # Insert weight vector into nearest list for NM weight vector
              #
              nm_nearest_list = nn_nm_train_w_vec_dict[nm_w_vec]
              if (nm_nearest_list == []):
                nn_nm_train_w_vec_dict[nm_w_vec] = [(nm_dist, this_w_vec)]
              elif (nm_dist < nm_nearest_list[-1][0]):
                nm_nearest_list.append((nm_dist, this_w_vec))
                nm_nearest_list.sort()  # Smallest distances first
                nn_nm_train_w_vec_dict[nm_w_vec] = nm_nearest_list[:k1]

            # Now calculate sum of k nearest distances and insert into heap
            #
            dist_sum = 0.0
            for (dist_val, match_type, train_w_vec) in nearest_list[:k]:
              dist_sum += dist_val

            heapq.heappush(nearest_w_vec_heap, (dist_sum, this_w_vec))

            # Insert into dictionary of weight vectors to be classified
            #
            nn_to_classify_w_vec_dict[this_w_vec] = (nearest_list, dist_sum)

      assert len(nearest_w_vec_heap) == len(nn_to_classify_w_vec_dict), \
             (len(nearest_w_vec_heap), len(nn_to_classify_w_vec_dict))

      # Step 2: Insert element on top of heap (the overall nearest to either
      # training set) into one of the training sets and remove from dictionary
      # of weight vectors to be classified. Then re-calculate the distances to
      # its nesrest neighbours yet to be classified
      #
      while (nn_to_classify_w_vec_dict != {}):

        nearest_w_vec_info = heapq.heappop(nearest_w_vec_heap)
        nearest_w_vec =      nearest_w_vec_info[1]

        while (nearest_w_vec not in nn_to_classify_w_vec_dict):
          nearest_w_vec_info = heapq.heappop(nearest_w_vec_heap)
          nearest_w_vec =      nearest_w_vec_info[1]

        nearest_list = nn_to_classify_w_vec_dict[nearest_w_vec][0]

        # Remove from dictionary of weight vectors to classify
        #
        del nn_to_classify_w_vec_dict[nearest_w_vec]

        # Determine if this nearest is going to be inserted into the match or
        # non-match training set
        #
        num_m, num_nm = 0,0
        for (dist_val, match_type, train_w_vec) in nearest_list[:k]:
          if (match_type == 'M'):
            num_m += 1
          else:
            num_nm += 1

        if (num_m > num_nm):
          nearest_w_vec_type = 'M'
        else:
          nearest_w_vec_type = 'NM'

        # If there are more weight vectors to classify, update their nearest
        # lists
        #
        if (nn_to_classify_w_vec_dict != {}):

          # Get weight vectors from so far unclassified weight vectors via the
          # nearest lists of the weight vectors from the training sets that are
          # nearest to this weight vector
          #
          nearest_w_vec_set = set()  # Only add unique weight vectors

          for (dist_val, match_type, train_w_vec) in nearest_list:

            if (match_type == 'M'):  # Get nearest to match training example
              train_nearest_list = nn_m_train_w_vec_dict[train_w_vec]

            else:  # Nearest to non-match training example
              train_nearest_list = nn_nm_train_w_vec_dict[train_w_vec]

            for (dist, w_vec) in train_nearest_list:
              if (w_vec in nn_to_classify_w_vec_dict):  # Still to classify
                nearest_w_vec_set.add(w_vec)

                assert w_vec not in nn_m_train_w_vec_dict
                assert w_vec not in nn_nm_train_w_vec_dict

          # Generate the nearest list of the new training example
          #
          new_train_w_vec_nearest_list = []

          # For all selected so far not classified weight vector update their
          # nearest lists as well as the nearest lists of their corresponding
          # nearest weight vectors from the training sets
          #
          for this_w_vec in nearest_w_vec_set:

            this_nearest_list = nn_to_classify_w_vec_dict[this_w_vec][0]
            this_dist_sum =     nn_to_classify_w_vec_dict[this_w_vec][1]
            largest_dist =      this_nearest_list[-1][0]

            # Claculate distance to the new training weight vector
            #
            new_dist = dist_meas(nearest_w_vec, this_w_vec)

            if (new_dist < largest_dist):
              this_nearest_list.append((new_dist, nearest_w_vec_type,
                                        nearest_w_vec))
              this_nearest_list.sort()  # Smallest distances first
              this_nearest_list = this_nearest_list[:k1]

            # Insert into either the nearest list of the new training example
            #
            if (new_train_w_vec_nearest_list == []):
              new_train_w_vec_nearest_list.append((new_dist, this_w_vec))
            elif (new_dist < new_train_w_vec_nearest_list[-1][0]):
              new_train_w_vec_nearest_list.append((new_dist, this_w_vec))
              new_train_w_vec_nearest_list.sort()  # Smallest distances first
              new_train_w_vec_nearest_list = \
                   new_train_w_vec_nearest_list[:k1]

            # Calculate new distance sum, and if smaller insert into heap
            #
            dist_sum = 0.0
            for (dist_val, match_type, train_w_vec) in this_nearest_list[:k]:
              dist_sum += dist_val

            if (dist_sum < this_dist_sum):
              heapq.heappush(nearest_w_vec_heap, (dist_sum, this_w_vec))

            # Update in dictionary of weight vectors to be classified
            #
            nn_to_classify_w_vec_dict[this_w_vec] = (this_nearest_list,
                                                     dist_sum)

        else:  # Last step, no more weight vectors to be classified
          new_train_w_vec_nearest_list = []  # Not needed anymore

        # Finally insert the new training example into appropriate training set
        #
        if (nearest_w_vec_type == 'M'):
          nn_m_train_w_vec_dict[nearest_w_vec] = new_train_w_vec_nearest_list
        else:
          nn_nm_train_w_vec_dict[nearest_w_vec] = new_train_w_vec_nearest_list

      # Save final training sets for later use
      #
      self.nn_m_train_w_vec_set =  set(nn_m_train_w_vec_dict.keys())
      self.nn_nm_train_w_vec_set = set(nn_nm_train_w_vec_dict.keys())

    else:
      logging.exception('Illegal step classifier method: %s' % \
                        (str(self.s2_classifier)))
      raise Exception

  # ---------------------------------------------------------------------------

  def test(self, w_vec_dict, match_set, non_match_set):
    """Method to test a classifier using the given weight vector dictionary and
       match and non-match sets of record identifier pairs.

       Weight vectors will be assigned to matches or non-matches according to
       either the SVM or the K-means classification. No weight vector will be
       assigned to the possible match set.

       Will return a confusion matrix as a list of the form: [TP, FN, FP, TN].
    """

    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)
    auxiliary.check_is_set('match_set', match_set)
    auxiliary.check_is_set('non_match_set', non_match_set)

    # Check that match and non-match sets are separate and do cover all weight
    # vectors given
    #
    if (len(match_set.intersection(non_match_set)) > 0):
      logging.exception('Intersection of match and non-match set not empty')
      raise Exception
    if ((len(match_set)+len(non_match_set)) != len(w_vec_dict)):
      logging.exception('Weight vector dictionary of different length than' + \
                        ' summed lengths of match and non-match sets: ' + \
                        '%d / %d+%d=%d' % (len(w_vec_dict), len(match_set),
                        len(non_match_set), len(match_set)+len(non_match_set)))
      raise Exception

    logging.info('')
    logging.info('Testing two-step classifier using %d weight vectors' % \
                 (len(w_vec_dict)))
    logging.info('  Match and non-match sets with %d and %d entries' % \
                 (len(match_set), len(non_match_set)))

    num_true_m =   0
    num_false_m =  0
    num_true_nm =  0
    num_false_nm = 0

    if (self.s2_classifier[0] == 'svm'):  # SVM classifier - - - - - - - - - -

      if (self.svm_model == None):
        logging.warn('SVM has not been trained, testing not possible')
        return [0,0,0,0]

      for (rec_id_tuple, w_vec) in w_vec_dict.iteritems():

        if (self.svm_version == 'old'):
          if (self.svm_model.predict(w_vec) == 1.0):  # Match prediction
            pred_match = True
          else:
            pred_match = False
# PC
        else:  # New SVM module version
          x0, max_idx = svm.gen_svm_nodearray(w_vec)

          if (svm.libsvm.svm_predict(self.svm_model, x0) == 1.0):  # Match
            pred_match = True
          else:
            pred_match = False

        if (pred_match == True):
          if (rec_id_tuple in match_set):
            num_true_m += 1
          else:
            num_false_m += 1
        else:  # Non-match prediction
          if (rec_id_tuple in non_match_set):
            num_true_nm += 1
          else:
            num_false_nm += 1

    elif (self.s2_classifier[0] == 'kmeans'):  # K-means clustering - - - - - -

      if ((self.m_centroid == None) or (self.nm_centroid == None)):
        logging.warn('K-means centroids have not been calculated, testing ' + \
                     ' not possible')
        return [0,0,0,0]

      dist_meas = self.s2_classifier[1]

      for (rec_id_tuple, w_vec) in w_vec_dict.iteritems():

        m_dist =  dist_meas(w_vec, self.m_centroid)
        nm_dist = dist_meas(w_vec, self.nm_centroid)

        if (m_dist < nm_dist):  # Assign to match cluster
          if (rec_id_tuple in match_set):
            num_true_m += 1
          else:
            num_false_m += 1
        else:
          if (rec_id_tuple in non_match_set):
            num_true_nm += 1
          else:
            num_false_nm += 1

    elif (self.s2_classifier[0] == 'nn'):  # Nearest neighbour classifier - - -

      nn_m_train_w_vec_set =  self.nn_m_train_w_vec_set
      nn_nm_train_w_vec_set = self.nn_nm_train_w_vec_set

      if ((nn_m_train_w_vec_set == None) or (nn_nm_train_w_vec_set == None)):
        logging.warn('Nearest neighbour classifier has not been trained, ' + \
                     'testing not possible')
        return [0,0,0,0]

      dist_meas = self.s2_classifier[1]
      k =         self.s2_classifier[2]

      for (rec_id_tuple, w_vec) in w_vec_dict.iteritems():

        this_w_vec = tuple(w_vec)  # Tuple can be used as dictionary key

        # Check if this weight vector is in one of the training sets
        #
        if (this_w_vec in nn_m_train_w_vec_set):
          if (rec_id_tuple in match_set):
            num_true_m += 1
          else:
            num_false_m += 1

        elif (this_w_vec in nn_nm_train_w_vec_set):
          if (rec_id_tuple in non_match_set):
            num_true_nm += 1
          else:
            num_false_nm += 1

        else:  # Have to find its k nearest neighbours from training sets

          nearest_list = [(9999999, '')]  # List of k nearest weight vectors

          for m_w_vec in nn_m_train_w_vec_set:
            m_dist = dist_meas(this_w_vec, m_w_vec)
            nearest_list.append((m_dist,'M'))
            nearest_list.sort()  # Smallest distances first
            nearest_list = nearest_list[:k]  # Only keep k nearest elements

          for nm_w_vec in nn_nm_train_w_vec_set:
            nm_dist = dist_meas(this_w_vec, nm_w_vec)
            nearest_list.append((nm_dist,'NM'))
            nearest_list.sort()  # Smallest distances first
            nearest_list = nearest_list[:k]  # Only keep k nearest elements

          # Now determine if this the the overall new nearest weight vector to
          # either match or non-match training examples
          #
          dist_sum = 0.0  # Sum of its k distances
          num_m =    0    # Number o matches in k nearest
          num_nm =   0    # Number of non-matches in k nearest

          for (dist_val, set_val) in nearest_list:
            dist_sum += dist_val
            if (set_val == 'M'):
              num_m += 1
            else:
              num_nm += 1

          if (num_m > num_nm):  # NN classifies this as a match
            if (rec_id_tuple in match_set):
              num_true_m += 1
            else:
              num_false_m += 1

          else:  # NN classifies this as a non-match
            if (rec_id_tuple in non_match_set):
              num_true_nm += 1
            else:
              num_false_nm += 1

    else:
      logging.exception('Illegal step classifier method: %s' % \
                        (self.s2_classifier))
      raise Exception

    assert (num_true_m+num_false_nm+num_false_m+num_true_nm) == len(w_vec_dict)

    logging.info('  Results: TP = %d, FN = %d, FP = %d, TN = %d' % \
                 (num_true_m,num_false_nm,num_false_m,num_true_nm))

    return [num_true_m, num_false_nm, num_false_m, num_true_nm]

  # ---------------------------------------------------------------------------

  def cross_validate(self, w_vec_dict, match_set, non_match_set, n=10):
    """Method to conduct a cross validation using the given weight vector
       dictionary and match and non-match sets of record identifier pairs.

       Will return a confusion matrix as a list of the form: [TP, FN, FP, TN].

       The two-step classifier cannot perform cross validation, as the
       selection of match and non-match training examples is determined by the
       methods given by the user. Therefore, in this method the given weight
       vectors are trained and then tested once only by calling the 'train' and
       then the 'test' methods.

       See documentation of 'test' for more information.
    """

    logging.info('')
    logging.info('Cross validation for two-step classifier is the same as' + \
                 'testing.')

    self.train(w_vec_dict, match_set, non_match_set)

    return self.test(w_vec_dict, match_set, non_match_set)

  # ---------------------------------------------------------------------------

  def classify(self, w_vec_dict):
    """Method to classify the given weight vector dictionary using the trained
       classifier.

       Will return three sets with record identifier pairs: 1) match set,
       2) non-match set, and 3) possible match set.

       The possible match set will be empty, as this classifier classifies all
       weight vectors as either matches or non-matches.
    """

    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)

    match_set =      set()
    non_match_set =  set()
    poss_match_set = set()

    if (self.s2_classifier[0] == 'svm'):  # SVM classifier - - - - - - - - - -

      if (self.svm_model == None):
        logging.warn('SVM has not been trained, classification not possible')
        return set(), set(), set()

      for (rec_id_tuple, w_vec) in w_vec_dict.iteritems():

        if (self.svm_version == 'old'):
          if (self.svm_model.predict(w_vec) == 1.0):  # Match prediction
            match_set.add(rec_id_tuple)
          else:  # Non-match prediction
            non_match_set.add(rec_id_tuple)

        else:  # New SVM module version
          x0, max_idx = svm.gen_svm_nodearray(w_vec)

          if (svm.libsvm.svm_predict(self.svm_model, x0) == 1.0):  # Match
            match_set.add(rec_id_tuple)
          else:  # Non-match prediction
            non_match_set.add(rec_id_tuple)

    elif (self.s2_classifier[0] == 'kmeans'):  # K-means clustering - - - - - -

      if ((self.m_centroid == None) or (self.nm_centroid == None)):
        logging.warn('K-means centroids have not been calculated, ' + \
                     'classification not possible')
        return set(), set(), set()

      dist_meas = self.s2_classifier[1]

      for (rec_id_tuple, w_vec) in w_vec_dict.iteritems():

        m_dist =  dist_meas(w_vec, self.m_centroid)
        nm_dist = dist_meas(w_vec, self.nm_centroid)

        if (m_dist < nm_dist):  # Assign to match set
          match_set.add(rec_id_tuple)
        else:
          non_match_set.add(rec_id_tuple)

    elif (self.s2_classifier[0] == 'nn'):  # Nearest neighbour classifier - - -

      nn_m_train_w_vec_set =  self.nn_m_train_w_vec_set
      nn_nm_train_w_vec_set = self.nn_nm_train_w_vec_set

      if ((nn_m_train_w_vec_set == None) or (nn_nm_train_w_vec_set == None)):
        logging.warn('Nearest neighbour classifier has not been trained, ' + \
                     'classification not possible')
        return set(), set(), set()

      dist_meas = self.s2_classifier[1]
      k =         self.s2_classifier[2]

      for (rec_id_tuple, w_vec) in w_vec_dict.iteritems():

        this_w_vec = tuple(w_vec)  # Tuple can be used as dictionary key

        # Check if this weight vector is in one of the training sets
        #
        if (this_w_vec in nn_m_train_w_vec_set):
          match_set.add(rec_id_tuple)

        elif (this_w_vec in nn_nm_train_w_vec_set):
          non_match_set.add(rec_id_tuple)

        else:  # Have to find its k nearest neighbours from training sets

          nearest_list = []  # List of k nearest weight vectors

          for m_w_vec in nn_m_train_w_vec_set:
            m_dist = dist_meas(this_w_vec, m_w_vec)
            nearest_list.append((m_dist,'M'))
            nearest_list.sort()  # Smallest distances first
            nearest_list = nearest_list[:k]  # Only keep k nearest elements

          for nm_w_vec in nn_nm_train_w_vec_set:
            nm_dist = dist_meas(this_w_vec, nm_w_vec)
            nearest_list.append((nm_dist,'NM'))
            nearest_list.sort()  # Smallest distances first
            nearest_list = nearest_list[:k]  # Only keep k nearest elements

          # Now determine if this the the overall new nearest weight vector to
          # either match or non-match training examples
          #
          dist_sum = 0.0  # Sum of its k distances
          num_m =    0    # Number o matches in k nearest
          num_nm =   0    # Number of non-matches in k nearest

          for (dist_val, set_val) in nearest_list:
            dist_sum += dist_val
            if (set_val == 'M'):
              num_m += 1
            else:
              num_nm += 1

          if (num_m > num_nm):  # NN classifies this as a match
            match_set.add(rec_id_tuple)
          else:
            non_match_set.add(rec_id_tuple)

    else:
      logging.exception('Illegal step classifier method: %s' % \
                        (str(self.s2_classifier)))
      raise Exception

    assert (len(match_set) + len(non_match_set) + len(poss_match_set)) == \
            len(w_vec_dict)

    logging.info('Classified %d weight vectors: %d as matches, %d as ' % \
                 (len(w_vec_dict), len(match_set), len(non_match_set)) + \
                 'non-matches, and %d as possible matches' % \
                 (len(poss_match_set)))

    return match_set, non_match_set, poss_match_set

# =============================================================================

class TAILOR(Classifier):
  """Implements the unsupervised hybrid classifier (based on k-means clustering
     followed by SVM classification) as described in the paper:

     TAILOR: A record linkage toolbox (Elfeky MG, Verykios VS, Elmagarmid AK,
     ICDE, San Jose, 2002.

     Note that in TAILOR originally a decision tree was used for classification
     rather than an SVM.

     Three clusters will be generated in a first step, one for matches,
     non-matches and possible matches each, and in the second step the match
     and non-match clusters will be used to train a binary SVM classifier.

     The arguments that have to be set when this classifier is initialised are:

       max_iter_count  The maximum number of iterations allowed.
       dist_measure    A function that calculates a distance measure between
                       two vectors (see the Febrl mymath.py module for such
                       functions).
       sample          A number between 0 and 100 that gives the percentage of
                       weight vectors that will be randomly selected and used
                       for clustering in the training process. If set to 100
                       (the default) then all given weight vectors will be
                       used.
       kernel_type     The kernel type from from libsvm. Default value LINEAR,
                       other possibilities are: POLY, RBF, SIGMOID.
       C               The 'C' parameter from libsvm. Default value is 10.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor. Process the classifer specific arguments first, then call
       the base class constructor.
    """

    # Check if svm module is installed or not
    #
    if (imp_svm == False):
      logging.exception('Module "svm.py" not installed, cannot use ' + \
                        'SuppVectorMach classifier')
      raise Exception

    self.max_iter_count = None
    self.dist_measure =   None
    self.sample =         100.0
    self.svm_type =       svm.C_SVC
    self.kernel_type =    'LINEAR'
    self.C =              10
    self.svm_model =      None  # Will be set in train() method

    base_kwargs = {}  # Dictionary, will contain unprocessed arguments for base
                      # class constructor

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('max_it')):
        auxiliary.check_is_integer('max_iter_count', value)
        auxiliary.check_is_positive('max_iter_count', value)
        self.max_iter_count = value

      elif (keyword.startswith('dist_m')):
        auxiliary.check_is_function_or_method('dist_measure', value)
        self.dist_measure = value

      elif (keyword.startswith('samp')):
        if (value != None):
          auxiliary.check_is_percentage('sample', value)
          self.sample = value

      elif (keyword.startswith('kernel')):
        auxiliary.check_is_string('kernel_type', value)
        if (value not in ['LINEAR', 'POLY', 'RBF', 'SIGMOID']):
          logging.exception('Illegal value for kernel type: %s ' % (value) + \
                            '(possible are: LINEAR, POLY, RBF, SIGMOID)')
          raise Exception
        self.kernel_type = value

      elif (keyword == 'C'):
        auxiliary.check_is_number('C', value)
        auxiliary.check_is_not_negative('C', value)
        self.C = value

      else:
        base_kwargs[keyword] = value

    Classifier.__init__(self, base_kwargs)  # Initialise base class

    # Check attribute values are set and valid - - - - - - - - - - - - - - - -
    #
    auxiliary.check_is_integer('max_iter_count', self.max_iter_count)
    auxiliary.check_is_positive('max_iter_count', self.max_iter_count)
    auxiliary.check_is_function_or_method('dist_measure', self.dist_measure)

    self.log([('Maximum iteration count', self.max_iter_count),
              ('Distance measure function', self.dist_measure),
              ('Sampling rate', self.sample),
              ('SVM kernel type', self.kernel_type),
              ('C', self.C)]) # Log a message

    # If the weight vector dictionary and both match and non-match sets - - - -
    # are given start the training process
    #
    if ((self.train_w_vec_dict != None) and (self.train_match_set != None) \
        and (self.train_non_match_set != None)):
      self.train(self.train_w_vec_dict, self.train_match_set,
                 (self.train_non_match_set))

  # ---------------------------------------------------------------------------

  def train(self, w_vec_dict, match_set, non_match_set):
    """Method to train a classifier using the given weight vector dictionary.
       Note that the given match and non-match sets of record identifier pairs
       will not be used (unsupervised training).

       This method will calculate three cluster centroids (one for matches,
       non-matches and possible matches each), possibly using a sampled sub-set
       of the weight vectors given, and then use the match and non-match
       clusters to train a SVM.
    """

    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)

    self.train_w_vec_dict =    w_vec_dict  # Save
    self.train_match_set =     match_set
    self.train_non_match_set = non_match_set

    # Get a random vector dictionary element to get dimensionality of vectors
    #
    (rec_id_tuple, w_vec) = w_vec_dict.popitem()
    v_dim = len(w_vec)
    w_vec_dict[rec_id_tuple] = w_vec  # Put back in

    logging.info('Train TAILOR classifier using %d weight vectors' % \
                 (len(w_vec_dict)))
    logging.info('  Dimensionality:   %d' % (v_dim))

    # Sample the weight vectors - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.sample == 100.0):
      use_w_vec_dict = w_vec_dict

    else:
      num_w_vec_sample = max(2, int(len(w_vec_dict)*self.sample/100.0))

      use_w_vec_dict = {}  # Create a new weight vector dictionary with samples

      rec_id_tuple_sample = random.sample(w_vec_dict.keys(),num_w_vec_sample)
      assert len(rec_id_tuple_sample) == num_w_vec_sample

      for rec_id_tuple in rec_id_tuple_sample:
        use_w_vec_dict[rec_id_tuple] = w_vec_dict[rec_id_tuple]

    logging.info('  Number of weight vectors to be used for clustering: %d' % \
                 (len(use_w_vec_dict)))

    # Initialise the cluster centroid - - - - - - - - - - - - - - - - - - - - -
    #
    m_centroid =  [-999.99]*v_dim  # Get the minimum and maximum values in
    nm_centroid = [999.99]*v_dim   # each weight vector element

    for w_vec in use_w_vec_dict.itervalues():
      for i in range(v_dim):
        m_centroid[i] =  max(w_vec[i], m_centroid[i])
        nm_centroid[i] = min(w_vec[i], nm_centroid[i])

    # Set possible match centroid half-way in between
    #
    pm_centroid = []

    for i in range(v_dim):
      pm_centroid.append((m_centroid[i]+nm_centroid[i])/2.0)


    logging.info('Initial cluster centroids:')
    logging.info('  Initial match centroid:     %s' % \
                 (auxiliary.str_vector(m_centroid)))
    logging.info('  Initial non-match centroid: %s' % \
                 (auxiliary.str_vector(nm_centroid)))
    logging.info('  Initial possible match centroid: %s' % \
                 (auxiliary.str_vector(pm_centroid)))

    # Step 1: k-means clustering - - - - - - - - - - - - - - - - - - - - - - -
    #
    cluster_assign_dict = {}  # Dictionary with cluster assignments

    iter_cnt = 1  # Iteration counter

    num_changed = 1

    while (num_changed > 0) and (iter_cnt < self.max_iter_count):

      num_changed =     0
      new_m_centroid =  [0.0]*v_dim  # Summed new distances
      new_nm_centroid = [0.0]*v_dim
      new_pm_centroid = [0.0]*v_dim

      # Step 1a: Calculate cluster membership for each weight vector - - - - -
      #
      num_m =  0  # Number of weight vectors assigned to matches
      num_nm = 0  # Number of weight vectors assigned to non-matches
      num_pm = 0  # Number of weight vectors assigned to possible matches

      for (rec_id_tuple, w_vec) in use_w_vec_dict.iteritems():

        m_dist =  self.dist_measure(w_vec, m_centroid)
        nm_dist = self.dist_measure(w_vec, nm_centroid)
        pm_dist = self.dist_measure(w_vec, pm_centroid)

        if ((m_dist < nm_dist) and (m_dist < pm_dist)):  # Assign to matches
          old_assign = cluster_assign_dict.get(rec_id_tuple, 'X')
          if (old_assign != 'M'):
            num_changed += 1
          cluster_assign_dict[rec_id_tuple] = 'M'
          num_m += 1

          for i in range(v_dim):  # Add to summed cluster distances
            new_m_centroid[i] += w_vec[i]

        elif (nm_dist < pm_dist):  # Assign to non-matches
          old_assign = cluster_assign_dict.get(rec_id_tuple, 'X')
          if (old_assign != 'NM'):
            num_changed += 1
          cluster_assign_dict[rec_id_tuple] = 'NM'
          num_nm += 1

          for i in range(v_dim):  # Add to summed cluster distances
            new_nm_centroid[i] += w_vec[i]

        else:  # Add to possible matches
          old_assign = cluster_assign_dict.get(rec_id_tuple, 'X')
          if (old_assign != 'PM'):
            num_changed += 1
          cluster_assign_dict[rec_id_tuple] = 'PM'
          num_pm += 1

          for i in range(v_dim):  # Add to summed cluster distances
            new_pm_centroid[i] += w_vec[i]

      num_all = len(cluster_assign_dict)

      if ((num_m + num_nm + num_pm) != num_all):
        logging.exception('Not all %d weight vectors assigned: ' + \
                          'M=%d, NM=%d, PM=%d' % \
                          (num_all, num_m, num_nm, num_pm))
        raise Exception

      if (num_m == 0) or (num_nm == 0):
        logging.warn('One cluster is empty: matches=%d, non-matches=%d, ' % \
                     (num_m, num_nm) + ' possible matches=%d' % (num_pm))
        break  # Stop K-means iterations

      # Step 1b: Calculate new cluster centroids - - - - - - - - - - - - - - -
      #
      for i in range(v_dim):  # Normalise new centroids
        new_m_centroid[i] /=  float(num_m)
        new_nm_centroid[i] /= float(num_nm)
        new_pm_centroid[i] /= float(num_pm)

      m_centroid =  new_m_centroid  # Assign new centroids
      nm_centroid = new_nm_centroid
      pm_centroid = new_pm_centroid

      logging.info('Iteration %d: %d vectors changed cluster assignment' % \
            (iter_cnt, num_changed))

      iter_cnt += 1

    self.m_centroid =  m_centroid  # Save for later use
    self.nm_centroid = nm_centroid
    self.pm_centroid = pm_centroid

    logging.info('Final cluster centroids:')
    logging.info('  Match centroid:          %s' % \
                 (auxiliary.str_vector(m_centroid)))
    logging.info('  Non-match centroid:      %s' % \
                 (auxiliary.str_vector(nm_centroid)))
    logging.info('  Possible match centroid: %s' % \
                 (auxiliary.str_vector(pm_centroid)))
    logging.info('  Cluster sizes: M=%d, NM=%d, PM=%d' % \
                 (num_m, num_nm, num_pm))

    if ((num_m == 0) or (num_nm == 0)):
      logging.warn('One of the training clusters is empty - cannot train SVM')
      self.svm_model = None
      return

    # Step 2: SVM classification using match and non-match clusters - - - - - -
    #
    train_data =   []
    train_labels = []

    for (rec_id_tuple, w_vec) in use_w_vec_dict.iteritems():
      if (cluster_assign_dict[rec_id_tuple] == 'M'):
        train_data.append(w_vec)
        train_labels.append(1.0)  # Match class

      elif (cluster_assign_dict[rec_id_tuple] == 'NM'):
        train_data.append(w_vec)
        train_labels.append(-1.0)  # Non-match class

    assert len(train_data) == num_m + num_nm

    # Initialise and train the SVM - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.kernel_type == 'LINEAR'):
      svm_kernel = svm.LINEAR
    elif (self.kernel_type == 'POLY'):
      svm_kernel = svm.POLY
    elif (self.kernel_type == 'RBF'):
      svm_kernel = svm.RBF
    elif (self.kernel_type == 'SIGMOID'):
      svm_kernel = svm.SIGMOID

    svm_prob =  svm.svm_problem(train_labels, train_data)

    # Due to change in SVM parameter setting in svm module, we need to catch
    # possible error
    #
    try:
      svm_param = svm.svm_parameter(svm_type = svm.C_SVC, C=self.C,
                                    kernel_type=svm_kernel)
      self.svm_model = svm.svm_model(svm_prob, svm_param)
      self.svm_version = 'old'

    except:
      svm_param = svm.svm_parameter('-s %d -c %f -t %d' % \
                  (svm.C_SVC, self.C, svm_kernel))
      self.svm_model = svm.libsvm.svm_train(svm_prob, svm_param)
      self.svm_version = 'new'

    logging.info('Trained SVM with %d training examples' % \
                 (len(use_w_vec_dict)))

  # ---------------------------------------------------------------------------

  def test(self, w_vec_dict, match_set, non_match_set):
    """Method to test a classifier using the given weight vector dictionary and
       match and non-match sets of record identifier pairs.

       Weight vectors will be assigned to matches or non-matches according to
       the SVM classification. No weight vector will be assigned to the
       possible match set.

       Will return a confusion matrix as a list of the form: [TP, FN, FP, TN].
    """

    if (self.svm_model == None):
      logging.warn('SVM has not been trained, testing not possible')
      return [0,0,0,0]

    svm_version = self.svm_version  # Shortcut

    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)
    auxiliary.check_is_set('match_set', match_set)
    auxiliary.check_is_set('non_match_set', non_match_set)

    # Check that match and non-match sets are separate and do cover all weight
    # vectors given
    #
    if (len(match_set.intersection(non_match_set)) > 0):
      logging.exception('Intersection of match and non-match set not empty')
      raise Exception
    if ((len(match_set)+len(non_match_set)) != len(w_vec_dict)):
      logging.exception('Weight vector dictionary of different length than' + \
                        ' summed lengths of match and non-match sets: ' + \
                        '%d / %d+%d=%d' % (len(w_vec_dict), len(match_set),
                        len(non_match_set), len(match_set)+len(non_match_set)))
      raise Exception

    num_true_m =   0
    num_false_m =  0
    num_true_nm =  0
    num_false_nm = 0

    for (rec_id_tuple, w_vec) in w_vec_dict.iteritems():

      if (svm_version == 'old'):
        if (self.svm_model.predict(w_vec) == 1.0):  # Match prediction
          pred_match = True
        else:
          pred_match = False

      else:  # New SVM module version
        x0, max_idx = svm.gen_svm_nodearray(w_vec)

        if (svm.libsvm.svm_predict(self.svm_model, x0) == 1.0):  # Match
          pred_match = True
        else:
          pred_match = False

      if (pred_match == True):
        if (rec_id_tuple in match_set):
          num_true_m += 1
        else:
          num_false_m += 1
      else:  # Non-match prediction
        if (rec_id_tuple in non_match_set):
          num_true_nm += 1
        else:
          num_false_nm += 1

    assert (num_true_m+num_false_nm+num_false_m+num_true_nm) == len(w_vec_dict)

    logging.info('  Results: TP = %d, FN = %d, FP = %d, TN = %d' % \
                 (num_true_m,num_false_nm,num_false_m,num_true_nm))

    return [num_true_m, num_false_nm, num_false_m, num_true_nm]

  # --------------------------------------------------------------------------

  def cross_validate(self, w_vec_dict, match_set, non_match_set, n=10):
    """Method to conduct a cross validation using the given weight vector
       dictionary and match and non-match sets of record identifier pairs.

       Will return a confusion matrix as a list of the form: [TP, FN, FP, TN].

       The cross validation approach splits the weight vector dictionary into
       'n' parts (and 'n' corresponding sub-set for matches and non-matches),
       and then generates 'n' TAILOR classifications, tests them and finally
       returns the average performance of these 'n' classifiers.
    """

    auxiliary.check_is_integer('n', n)
    auxiliary.check_is_positive('n', n)
    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)
    auxiliary.check_is_set('match_set', match_set)
    auxiliary.check_is_set('non_match_set', non_match_set)

    # Check that match and non-match sets are separate and do cover all weight
    # vectors given
    #
    if (len(match_set.intersection(non_match_set)) > 0):
      logging.exception('Intersection of match and non-match set not empty')
      raise Exception
    if ((len(match_set)+len(non_match_set)) != len(w_vec_dict)):
      logging.exception('Weight vector dictionary of different length than' + \
                        ' summed lengths of match and non-match sets: ' + \
                        '%d / %d+%d=%d' % (len(w_vec_dict), len(match_set),
                        len(non_match_set), len(match_set)+len(non_match_set)))
      raise Exception

    # Get a random vector dictionary element to get dimensionality of vectors
    #
    (rec_id_tuple, w_vec) = w_vec_dict.popitem()
    v_dim = len(w_vec)
    w_vec_dict[rec_id_tuple] = w_vec  # Put back in

    logging.info('')
    logging.info('Conduct %d-fold cross validation on TAILOR classifier ' % \
                 (n) + 'using %d weight vectors' % (len(w_vec_dict)))
    logging.info('  Match and non-match sets with %d and %d entries' % \
                 (len(match_set), len(non_match_set)))
    logging.info('  Dimensionality:   %d' % (v_dim))

    m_centroids =  []  # Keep the centroids from all folds
    nm_centroids = []

    # Create the sub-sets of record identifier pairs for folds - - - - - - - -
    #
    rec_id_tuple_list = w_vec_dict.keys()
    random.shuffle(rec_id_tuple_list)
    fold_num_rec_id_tuple = max(1,int(round(float(len(rec_id_tuple_list))/n)))

    # Split the weight vector dictionary and match and non-match sets into
    # (lists containing one entry per fold) and only store test elements
    #
    w_vec_dict_test_list = []
    m_set_test_list =      []
    nm_set_test_list =     []

    for fold in range(n):
      w_vec_dict_test_list.append({})
      m_set_test_list.append(set())
      nm_set_test_list.append(set())

    for fold in range(n):

      # Calculate start and end indices for test elements for this fold
      #
      if (fold == (n-1)):  # The last fold, get remainder of list
        start = fold*fold_num_rec_id_tuple
        this_fold_test_ids = rec_id_tuple_list[start:]
      else:  # All other folds
        start = fold*fold_num_rec_id_tuple
        end = start+fold_num_rec_id_tuple
        this_fold_test_ids = rec_id_tuple_list[start:end]

      for rec_id_tuple in this_fold_test_ids:

        w_vec_dict_test_list[fold][rec_id_tuple] = w_vec_dict[rec_id_tuple]

        if (rec_id_tuple in match_set):
          m_set_test_list[fold].add(rec_id_tuple)
        else:
          nm_set_test_list[fold].add(rec_id_tuple)

      assert len(w_vec_dict_test_list[fold]) == len(this_fold_test_ids)
      assert len(m_set_test_list[fold]) + len(nm_set_test_list[fold]) == \
             len(this_fold_test_ids)

    # Initialise the total classification results - - - - - - - - - - - - - - -
    #
    num_true_m =   0
    num_false_nm = 0
    num_false_m =  0
    num_true_nm =  0

    # Loop over folds - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Generate training and test dictionaries and sets
    #
    for fold in range(n):  # First extract test record identifier pairs

      this_fold_test_m_set =       m_set_test_list[fold]
      this_fold_test_nm_set =      nm_set_test_list[fold]
      this_fold_test_w_vec_dict =  w_vec_dict_test_list[fold]

      this_fold_train_m_set =  match_set.difference(m_set_test_list[fold])
      this_fold_train_nm_set = non_match_set.difference(nm_set_test_list[fold])
      this_fold_train_w_vec_dict = {}
      for f2 in range(n):
        if (f2 != fold):
          this_fold_train_w_vec_dict.update(w_vec_dict_test_list[f2])

      assert len(this_fold_test_m_set) + len(this_fold_train_m_set) == \
             len(match_set)
      assert len(this_fold_test_nm_set) + len(this_fold_train_nm_set) == \
             len(non_match_set)
      assert len(this_fold_test_w_vec_dict) + \
             len(this_fold_train_w_vec_dict) == len(w_vec_dict)

      assert this_fold_test_m_set.intersection(this_fold_train_m_set) == set()
      assert this_fold_test_m_set.intersection(this_fold_test_nm_set) == set()
      assert this_fold_test_m_set.intersection(this_fold_train_nm_set) == set()
      assert this_fold_test_nm_set.intersection(this_fold_train_m_set) ==set()
      assert this_fold_test_nm_set.intersection(this_fold_train_nm_set) ==set()
      assert this_fold_train_m_set.intersection(this_fold_train_nm_set) ==set()

      # Train and test TAILOR classifier on this fold's data
      #
      self.train(this_fold_train_w_vec_dict, this_fold_train_m_set,
                 this_fold_train_nm_set)

      [this_num_true_m,this_num_false_nm,this_num_false_m,this_num_true_nm] = \
                                           self.test(this_fold_test_w_vec_dict,
                                                     this_fold_test_m_set,
                                                     this_fold_test_nm_set)
      num_true_m +=   this_num_true_m
      num_false_nm += this_num_false_nm
      num_false_m +=  this_num_false_m
      num_true_nm +=  this_num_true_nm

    # Calculate final cross validation results - - - - - - - - - - - - - - - -
    #
    num_true_m /=   float(n)
    num_false_nm /= float(n)
    num_false_m /=  float(n)
    num_true_nm /=  float(n)

    logging.info('  Results: TP = %d, FN = %d, FP = %d, TN = %d' % \
                 (num_true_m,num_false_nm,num_false_m,num_true_nm))

    return [num_true_m, num_false_nm, num_false_m, num_true_nm]

  # ---------------------------------------------------------------------------

  def classify(self, w_vec_dict):
    """Method to classify the given weight vector dictionary using the trained
       classifier.

       Will return three sets with record identifier pairs: 1) match set,
       2) non-match set, and 3) possible match set.

       The possible match set will be empty, as this classifier classifies all
       weight vectors as either matches or non-matches.
    """

    if (self.svm_model == None):
      logging.warn('SVM has not been trained, classification not possible')
      return set(), set(), set()

    svm_version = self.svm_version  # Shortcut

    auxiliary.check_is_dictionary('w_vec_dict', w_vec_dict)

    match_set =      set()
    non_match_set =  set()
    poss_match_set = set()

    for (rec_id_tuple, w_vec) in w_vec_dict.iteritems():

      if (svm_version == 'old'):
        if (self.svm_model.predict(w_vec) == 1.0):  # Match prediction
          match_set.add(rec_id_tuple)
        else:  # Non-match prediction
          non_match_set.add(rec_id_tuple)

      else:  # New SVM module version
        x0, max_idx = svm.gen_svm_nodearray(w_vec)

        if (svm.libsvm.svm_predict(self.svm_model, x0) == 1.0):  # Match
          match_set.add(rec_id_tuple)
        else:  # Non-match prediction
          non_match_set.add(rec_id_tuple)

    assert (len(match_set) + len(non_match_set) + len(poss_match_set)) == \
            len(w_vec_dict)

    logging.info('Classified %d weight vectors: %d as matches, %d as ' % \
                 (len(w_vec_dict), len(match_set), len(non_match_set)) + \
                 'non-matches, and %d as possible matches' % \
                 (len(poss_match_set)))

    return match_set, non_match_set, poss_match_set


# =============================================================================
# Following are several auxiliary functions that are helpful for classification

def get_true_matches_nonmatches(weight_vec_dict, match_check_funct):
  """Checks for all weight vectors in the given dictionary if they are from
     a true match or a true non-match, assuming this information is (somehow)
     available in the record identifiers or the weight vectors.

     The function returns two sets with matches and non-matches, respectively,
     whose elements are the record identifier tuples.

     Arguments:
       weight_vec_dict    A dictionary containing weight vectors, with the keys
                          in the dictionary being record identifier tuples and
                          the values being the actual vectors.
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

  true_match_set =     set()
  true_non_match_set = set()

  for (rec_id_tuple, this_vec) in weight_vec_dict.iteritems():

    if (match_check_funct(rec_id_tuple[0], rec_id_tuple[1], this_vec) == True):
      true_match_set.add(rec_id_tuple)
    else:
      true_non_match_set.add(rec_id_tuple)

  return true_match_set, true_non_match_set

# -----------------------------------------------------------------------------

def extract_collapse_weight_vectors(manipulate_list, weight_vec_dict,
                                    vec_weights=None):
  """This function allows to manipulate the vectors in the given weight vector
     dictionary by (a) removing vector elements, (b) summing several vector
     elements into one, and (c) giving vector elements different weights (if
     the 'vec_weights' argument is given).

     Here are some examples (assuming the weight vectors are of length 6):

       manipulate_list = [(0,)]
         Returns vectors of dimensionality 1 containing only the first element
         of the input weight vectors.

       manipulate_list = [(0,1,2)]
         Returns vectors of dimensionality 1 containing the summed values of
         the first three elements of the input weight vectors.

       manipulate_list = [(0,1),(2,3),(4,5)]
         Returns vectors of dimensionality 3 containing the summed values of
         the first two, middle two, and last two elements of the input weight
         vectors.

     Arguments:
       manipulate_list  This has to be a list with tuples that details which
                        elements of the input vector shall be summed.
       weight_vec_dict  A dictionary containing weight vectors, with the keys
                        in the dictionary being record identifier tuples and
                        the values being the actual vectors.
       vec_weights      A vector (of same lengths as the vectors in the weight
                        vector dictionary) giving weights for each element in
                        the weight vectors, or None (default) in which case
                        all elements in the weight vectors get a weight 1.0.
                        All weights given in 'vec_weights' have to be positive.
  """

  auxiliary.check_is_dictionary('weight_vec_dict', weight_vec_dict)
  auxiliary.check_is_list('manipulate_list', manipulate_list)

  # Get a random vector dictionary element to get dimensionality of vectors
  #
  (rec_id_tuple, w_vec) = weight_vec_dict.popitem()
  v_dim = len(w_vec)
  weight_vec_dict[rec_id_tuple] = w_vec  # Put back in

  if (vec_weights != None):
    auxiliary.check_is_list('vec_weights', vec_weights)
    if (len(vec_weights) != v_dim):
      logging.exception('Argument "vec_weights" given is of different ' + \
                        'length compared to weight vectors in dictionary: ' + \
                        '%d / %d' % (len(vec_weights), v_dim))
      raise Exception
    for i in range(v_dim):
      auxiliary.check_is_positive('vec_weights[%d]' % (i), vec_weights[i])

  if (vec_weights == None):
    use_vec_weights = [1.0]*v_dim  # Uniform weights
  else:
    use_vec_weights = vec_weights

  # Check the manipulate list for correctness - - - - - - - - - - - - - - - - -
  #
  for i in range(len(manipulate_list)):
    auxiliary.check_is_tuple('manipulate_list[%d]' % (i), manipulate_list[i])
    for e in manipulate_list[i]:
      if (e >= v_dim):
        logging.exception('Element in tuple %d of the manipulate list ' % \
                          (i+1)+'is out of range: %s' % (str(manipulate_list)))
        raise Exception

  out_vec_dict = {}

  for (rec_id_tuple, this_vec) in weight_vec_dict.iteritems():
    new_vec = []
    for coll_tuple in manipulate_list:
      w = 0.0
      for e in coll_tuple:
        w += this_vec[e] * use_vec_weights[e]
      new_vec.append(w)
    out_vec_dict[rec_id_tuple] = new_vec

    assert len(new_vec) <= len(this_vec)

  assert len(out_vec_dict) == len(weight_vec_dict)

  return out_vec_dict

# =============================================================================

# Old stuff below, PC 9/08/07

# =============================================================================

def DecisionTree(weight_vec_dict, match_set, non_match_set):
  """Classifier that has access to the true and matches and non-matches and
     uses the supervised learning algorithm of ID3 decision tree induction
     using all weight vectors provided.

     The arguments that have to be set when this classifier is called are:

       weight_vec_dict  A dictionary containing weight vectors
       match_set        A set containing the true matches (as tuples of record
                        identifiers)
       non_match_set    A set containing the true non-matches (as tuples of
                        record identifiers)
  """

  assert len(match_set)+len(non_match_set) == len(weight_vec_dict)

  auxiliary.check_is_dictionary('weight_vec_dict', weight_vec_dict)
  auxiliary.check_is_set('match set', match_set)
  auxiliary.check_is_set('non match set', non_match_set)

  logging.info('')
  logging.info('Classify %d weight vectors using the ID3 decision tree ' % \
               (len(weight_vec_dict))+'classifier')
  logging.info('  Match and non-match sets with %d and %d entries' % \
               (len(match_set), len(non_match_set)))

  # Functions needed for decision tree classifier -----------------------------

  def entropy(weight_vec_dict, match_set, non_match_set):
    """Calculate the entropy of the given weight vector dictionary according to
       the match and non-match sets.
    """

    m_count =  0.0  # Count number of matched and non-matched weight vectors
    nm_count = 0.0

    for rec_id_tuple in weight_vec_dict.iterkeys():
      if rec_id_tuple in match_set:
        m_count += 1.0
      elif rec_id_tuple in non_match_set:
        nm_count += 1.0
      else:
        print 'error!', rec_id_tuple

    num_weight_vec = float(len(weight_vec_dict))

    data_entropy = 0.0

    if (m_count > 0.0):
      data_entropy += -m_count/num_weight_vec * \
                      math.log(m_count/num_weight_vec, 2)
    if (nm_count > 0.0):
      data_entropy += -nm_count/num_weight_vec * \
                      math.log(nm_count/num_weight_vec, 2)

    return data_entropy

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  def fitness(weight_vec_dict, split_column, match_set, non_match_set):
    """Calculate the information gain (reduction in entropy) that would result
       by splitting the given weight vector dictionary on the chosen column.
    """

    val_count_dict = {}

    # Calculate the frequency of each value in the selected split cloumn
    #
    for this_w_vector in weight_vec_dict.itervalues():
      val = this_w_vector[split_column]

      val_count = val_count_dict.get(val, 0) + 1.0
      val_count_dict[val] = val_count

    num_weight_vec = float(len(weight_vec_dict))

    subset_entropy = 0.0

    # Calculate the sum of the entropy for each sub-set of weight vectors
    # weighted by their probability of occuring
    #
    for (val, val_count) in val_count_dict.iteritems():

        val_prob = val_count / num_weight_vec

        weight_vec_subset_dict = {}

        for (rec_id_tuple, this_w_vector) in weight_vec_dict.iteritems():

          if (this_w_vector[split_column] == val):
            weight_vec_subset_dict[rec_id_tuple] = this_w_vector

        subset_entropy += val_prob * entropy(weight_vec_subset_dict,
                                             match_set, non_match_set)

    # Subtract the entropy of the chosen column from the entropy of the whole
    # weight vector dictionary with respect to the match and non-match sets
    #
    return entropy(weight_vec_dict, match_set, non_match_set) - subset_entropy

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  def majority_value(weight_vec_dict, match_set, non_match_set):
    """Return 'M' if the given weight vector dictionary contains more matches
       than non-matches, or 'NM' otherwise.
    """

    m_count =  0
    nm_count = 0

    # Count the number of matches and non-matches in the given weight vectors
    #
    for rec_id_tuple in weight_vec_dict.iterkeys():
      if rec_id_tuple in match_set:
        m_count += 1
      elif rec_id_tuple in non_match_set:
        nm_count += 1
      else:
        print 'error!', rec_id_tuple

    if (m_count > nm_count):
      return ('M', m_count)
    else:
      return ('NM', nm_count)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  def get_examples(weight_vec_dict, column, val):
    """Return a weight vector dictionary of all the weight vectors that have a
       value 'val' in column 'column'.
    """

    new_weight_vec_dict = {}

    for (rec_id_tuple, w_vector) in weight_vec_dict.iteritems():

      if (w_vector[column] == val):
        new_weight_vec_dict[rec_id_tuple] = w_vector

    return new_weight_vec_dict

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  def print_tree(tree, str):
    """Function that recursively crawls through the decision tree and prints it
       in readable format.
    """

    if isinstance(tree, dict):
      logging.info('    %scolumn %s:' % (str, tree.keys()[0]))
      for item in tree.values()[0].keys():
        logging.info('    %s    %s' % (str, item))
        print_tree(tree.values()[0][item], str + '    ')
    else:
        logging.info('    %s    ->    %s' % (str, tree))

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # Main recursive function to construct a decision tree
  #
  def create_decision_tree(weight_vec_dict, column_list, match_set,
                           non_match_set):
    """Return a new decision tree based on the examples given.
    """

    # Get the most common value and its count
    #
    (most_common_class, most_common_count) = \
                        majority_value(weight_vec_dict,match_set,non_match_set)

    if ((most_common_count == len(weight_vec_dict)) or (column_list == [])):
      return most_common_class  # All values are in the same class

    # Choose the next best column to best classify the data
    #
    best_gain =   0.0
    best_column = None

    for column in column_list:
      gain = fitness(weight_vec_dict, column, match_set, non_match_set)
      if (gain >= best_gain):
        best_gain =   gain
        best_column = column

    # Create a list of unique values in the chosen best column in the given
    # dictionary data
    #
    best_column_val_list = []

    for w_vector in weight_vec_dict.itervalues():
      if (w_vector[best_column] not in best_column_val_list):
        best_column_val_list.append(w_vector[best_column])

    # Create a new decision tree/node with the best column and an empty
    # dictionary object -- we'll fill that up next.
    #
    tree = {best_column:{}}

    # Create a new decision tree/sub-node for each of the values in the best
    # column
    #
    for val in best_column_val_list:

      # Create a sub-tree for the current value under the "best" field
      #
      this_weight_vec_dict = get_examples(weight_vec_dict, best_column, val)
      new_column_list = column_list[:]
      new_column_list.remove(best_column)

      subtree = create_decision_tree(this_weight_vec_dict, new_column_list,
                                     match_set, non_match_set)

      # Add the new sub-tree to the empty dictionary object in our new
      # tree/node we just created.
      #
      tree[best_column][val] = subtree

    return tree

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  def get_classification(w_vector, tree):
    """This function recursively traverses the decision tree and returns a
       classification for the given weight vector.
    """

    # If the current node is a string, then we've reached a leaf node and we
    # can return it as our answer
    #
    if isinstance(tree, str):
      return tree

    else:  # Traverse the tree further until a leaf node is found.

      column = tree.keys()[0]
      t = tree[column][w_vector[column]]

      return get_classification(w_vector, t)

  # Start decision tree classifier --------------------------------------------

  # Select a weight vector to get the vector length
  #
  centroid = weight_vec_dict.popitem()  # Get arbitrary weight vector
  weight_vec_dict[centroid[0]] = centroid[1]  # Put it back into the dictionary

  num_weights = len(centroid[1])
  column_list = range(num_weights)

  # Create the decision tree using all weight vectors - - - - - - - - - - - - -
  #
  id3_tree = create_decision_tree(weight_vec_dict, column_list, match_set,
                                  non_match_set)
  logging.info('  Built ID3 decision tree:')

##  print_tree(id3_tree,'')

  # Assign all weight vectors to either the match or non-match set - - - - - -
  #
  match_set =     set()
  non_match_set = set()

  for (rec_id_tuple, this_w_vector) in weight_vec_dict.iteritems():

    decs_class = get_classification(this_w_vector, id3_tree)

    if (decs_class == 'M'):
      match_set.add(rec_id_tuple)
    elif (decs_class == 'NM'):
      non_match_set.add(rec_id_tuple)
    else:
      print 'XX:', decs_class

  logging.info('Classified %d weight vectors as matches, %d as non-matches' % \
               (len(match_set), len(non_match_set)))

  return match_set, non_match_set

# =============================================================================
