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
# The Original Software is: "classificationTest.py"
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

"""Module classificationTest.py - Test module for classification.py.
"""

# -----------------------------------------------------------------------------

import random
import sys
import unittest
sys.path.append('..')

import classification
import mymath  # For K-means distance measures

import logging
my_logger = logging.getLogger()  # New logger at root level
my_logger.setLevel(logging.WARNING) # INFO)

# -----------------------------------------------------------------------------

class TestCase(unittest.TestCase):

  # Initialise test case  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def setUp(self):

    # Define a test weight vector dictionary with the following convention:
    # - if numbers in record identifiers are the same it is a match
    # - if numbers in record identifiers differe it is a non-match
    # - Weight vectors given here have all very high or low similarity values
    #
    self.w_vec_dict = {('1','1'):[1.0,  0.95, 1.0,  0.98, 1.0 ],
                       ('1','2'):[0.0,  0.05, 0.02, 0.1,  0.0 ],
                       ('1','3'):[0.06, 0.01, 0.15, 0.02, 0.88],
                       ('1','4'):[0.04, 0.1,  0.08, 0.0,  0.06],
                       ('2','2'):[0.98, 1.0,  1.0,  1.0,  0.92],
                       ('2','3'):[0.01, 0.1,  0.05, 0.16, 0.0 ],
                       ('2','4'):[0.12, 0.0,  0.01, 0.08, 0.04],
                       ('3','3'):[0.98, 1.0,  0.93, 0.9,  0.89],
                       ('3','4'):[0.1,  0.02, 0.05, 0.04, 0.01],
                       ('4','3'):[0.03, 0.1,  0.02, 0.1,  0.05],
                       ('4','4'):[0.92, 0.98, 0.95, 1.0,  0.95]}

    self.m_set =  set([('1','1'), ('2','2'), ('3','3'), ('4','4')])
    self.nm_set = set([('1','2'), ('1','3'), ('1','4'), ('2','3'),
                       ('2','4'), ('3','4'), ('4','3')])

    self.test_w_vec_dict = {('5','5'):[0.81, 1.0, 0.89, 1.0,  0.94],
                            ('1','5'):[0.02, 0.1, 0.01, 0.09, 0.0],
                            ('2','5'):[0.06, 0.1, 0.12, 0.0,  0.1]}

    # Randomly generate more weight vectors
    #
    for i in range(10,110):
      rec_tuple = ('%s' % (i) ,'%s' % (i))  # Generate matches
      w_vec = []
      for j in range(5):
        w_vec.append(1.0-random.random()/1.8)
      self.w_vec_dict[rec_tuple] = w_vec
      self.m_set.add(rec_tuple)

    for i in range(500):  # Generate non-matches
      x = random.randrange(10,110)
      y = random.randrange(10,110)
      rec_tuple = ('%s' % (x) ,'%s' % (y))
      while rec_tuple in self.w_vec_dict:
        x = random.randrange(10,110)
        y = random.randrange(10,110)
        rec_tuple = ('%s' % (x) ,'%s' % (y))
      w_vec = []
      for j in range(5):
        w_vec.append(random.random()/1.8)
      self.w_vec_dict[rec_tuple] = w_vec
      self.nm_set.add(rec_tuple)

    for i in range(500):  # Weight vectors for testing
      x = random.randrange(110,300)
      y = random.randrange(110,300)
      rec_tuple = ('%s' % (x) ,'%s' % (y))
      while rec_tuple in self.test_w_vec_dict:
        x = random.randrange(0,100)
        y = random.randrange(0,100)
        rec_tuple = ('%s' % (x) ,'%s' % (y))
      w_vec = []
      for j in range(5):
        if (random.random() < 0.5):
          w_vec.append(random.random())
        else:
          w_vec.append(1.0-random.random())
      self.test_w_vec_dict[rec_tuple] = w_vec

  # Clean up test case  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def tearDown(self):
    pass  # Nothing to clean up

  # ---------------------------------------------------------------------------
  #
  # Start test cases

  def testFellegiSunterClassifier(self):  # - - - - - - - - - - - - - - - - - -
    """Test Fellegi and Sunter classifier"""

    # Loop over various threshold pair values
    #
    for (lt,ut) in [(0.0,1.0), (0.2,0.7),(0.4,0.6),(0.5,0.5)]:

#      pass
#    for x in []:

      fs_class = classification.FellegiSunter(descr = 'fell-sunter',
                                            lower_t = lt,
                                            upper_t = ut)
      assert fs_class.lower_threshold == lt
      assert fs_class.upper_threshold == ut
      assert fs_class.train_w_vec_dict == None
      assert fs_class.train_match_set == None
      assert fs_class.train_non_match_set == None

      fs_class.train(self.w_vec_dict, self.m_set, self.nm_set)
      assert fs_class.lower_threshold == lt
      assert fs_class.upper_threshold == ut
      assert fs_class.train_w_vec_dict == None
      assert fs_class.train_match_set == None
      assert fs_class.train_non_match_set == None

      test_res = fs_class.test(self.w_vec_dict, self.m_set, self.nm_set)
      assert len(test_res) == 4
      assert fs_class.lower_threshold == lt
      assert fs_class.upper_threshold == ut
      assert fs_class.train_w_vec_dict == None
      assert fs_class.train_match_set == None
      assert fs_class.train_non_match_set == None

      test_res = fs_class.cross_validate(self.w_vec_dict, self.m_set,
                                         self.nm_set)
      assert len(test_res) == 4
      assert fs_class.lower_threshold == lt
      assert fs_class.upper_threshold == ut
      assert fs_class.train_w_vec_dict == None
      assert fs_class.train_match_set == None
      assert fs_class.train_non_match_set == None

      class_res = fs_class.classify(self.test_w_vec_dict)
      assert len(class_res) == 3
      assert fs_class.lower_threshold == lt
      assert fs_class.upper_threshold == ut
      assert fs_class.train_w_vec_dict == None
      assert fs_class.train_match_set == None
      assert fs_class.train_non_match_set == None
      assert isinstance(class_res[0], set) == True
      assert isinstance(class_res[1], set) == True
      assert isinstance(class_res[2], set) == True
      assert len(class_res[0]) + len(class_res[1]) + len(class_res[2]) == \
             len(self.test_w_vec_dict)

      # Train classifier when initialising it - - - - - - - - - - - - - - - - -
      #
      fs_class2 = classification.FellegiSunter(descr = 'fell-sunter',
                                             lower_t = lt,
                                             upper_t = ut,
                                    train_w_vec_dict = self.w_vec_dict,
                                     train_match_set = self.m_set,
                                 train_non_match_set = self.nm_set)

      assert fs_class2.lower_threshold == lt
      assert fs_class2.upper_threshold == ut
      assert fs_class2.train_w_vec_dict == self.w_vec_dict
      assert fs_class2.train_match_set == self.m_set
      assert fs_class2.train_non_match_set == self.nm_set

      test_res = fs_class2.test(self.w_vec_dict, self.m_set, self.nm_set)
      assert len(test_res) == 4
      assert fs_class2.lower_threshold == lt
      assert fs_class2.upper_threshold == ut
      assert fs_class2.train_w_vec_dict == self.w_vec_dict
      assert fs_class2.train_match_set == self.m_set
      assert fs_class2.train_non_match_set == self.nm_set

      class_res = fs_class2.classify(self.test_w_vec_dict)
      assert len(class_res) == 3
      assert fs_class2.lower_threshold == lt
      assert fs_class2.upper_threshold == ut
      assert fs_class2.train_w_vec_dict == self.w_vec_dict
      assert fs_class2.train_match_set == self.m_set
      assert fs_class2.train_non_match_set == self.nm_set
      assert isinstance(class_res[0], set) == True
      assert isinstance(class_res[1], set) == True
      assert isinstance(class_res[2], set) == True
      assert len(class_res[0]) + len(class_res[1]) + len(class_res[2]) == \
             len(self.test_w_vec_dict)


  def testOptimalThresholdClassifier(self):  # - - - - - - - - - - - - - - - -
    """Test optimal threshold classifier"""

    # Loop over various bin width and minimise methods
    #
    for (bw,mm) in [(0.1,'pos-neg'),(0.1,'pos'),(0.1,'neg'),
                    (0.02,'pos-neg'),(0.02,'pos'),(0.02,'neg'),
                    (0.005,'pos-neg'),(0.005,'pos'),(0.005,'neg')]:

#      pass
#    for x in []:

      ot_class = classification.OptimalThreshold(descr = 'optimal-threshold',
                                             bin_width = bw,
                                            min_method = mm)
      assert ot_class.bin_width == bw
      assert ot_class.min_method == mm
      assert ot_class.train_w_vec_dict == None
      assert ot_class.train_match_set == None
      assert ot_class.train_non_match_set == None

      ot_class.train(self.w_vec_dict, self.m_set, self.nm_set)
      assert ot_class.bin_width == bw
      assert ot_class.min_method == mm
      assert ot_class.train_w_vec_dict == self.w_vec_dict
      assert ot_class.train_match_set == self.m_set
      assert ot_class.train_non_match_set == self.nm_set
      for i in range(5):
        assert ot_class.opt_threshold_list[i] >= 0.0
        assert ot_class.opt_threshold_list[i] <= 1.0

      test_res = ot_class.test(self.w_vec_dict,self.m_set,self.nm_set)
      assert len(test_res) == 4
      assert ot_class.bin_width == bw
      assert ot_class.min_method == mm
      assert ot_class.train_w_vec_dict == self.w_vec_dict
      assert ot_class.train_match_set == self.m_set
      assert ot_class.train_non_match_set == self.nm_set
      for i in range(5):
        assert ot_class.opt_threshold_list[i] >= 0.0
        assert ot_class.opt_threshold_list[i] <= 1.0

      test_res = ot_class.cross_validate(self.w_vec_dict, self.m_set,
                                         self.nm_set, 2)
      assert len(test_res) == 4
      assert ot_class.bin_width == bw
      assert ot_class.min_method == mm
      for i in range(5):
        assert ot_class.opt_threshold_list[i] >= 0.0
        assert ot_class.opt_threshold_list[i] <= 1.0

      test_res = ot_class.cross_validate(self.w_vec_dict, self.m_set,
                                         self.nm_set, 3)
      assert len(test_res) == 4
      assert ot_class.bin_width == bw
      assert ot_class.min_method == mm
      for i in range(5):
        assert ot_class.opt_threshold_list[i] >= 0.0
        assert ot_class.opt_threshold_list[i] <= 1.0

      test_res = ot_class.cross_validate(self.w_vec_dict, self.m_set,
                                         self.nm_set, 5)
      assert len(test_res) == 4
      assert ot_class.bin_width == bw
      assert ot_class.min_method == mm
      for i in range(5):
        assert ot_class.opt_threshold_list[i] >= 0.0
        assert ot_class.opt_threshold_list[i] <= 1.0

      class_res = ot_class.classify(self.test_w_vec_dict)
      assert len(class_res) == 3
      assert isinstance(class_res[0], set) == True
      assert isinstance(class_res[1], set) == True
      assert isinstance(class_res[2], set) == True
      assert len(class_res[0]) + len(class_res[1]) + len(class_res[2]) == \
             len(self.test_w_vec_dict)
      for i in range(5):
        assert ot_class.opt_threshold_list[i] >= 0.0
        assert ot_class.opt_threshold_list[i] <= 1.0

      # Train classifier when initialising it - - - - - - - - - - - - - - - - -
      #
      ot_class2 = classification.OptimalThreshold(descr = 'optimal-threshold',
                                              bin_width = bw,
                                             min_method = mm,
                                       train_w_vec_dict = self.w_vec_dict,
                                        train_match_set = self.m_set,
                                    train_non_match_set = self.nm_set)
      assert ot_class2.bin_width == bw
      assert ot_class2.min_method == mm
      assert ot_class2.train_w_vec_dict == self.w_vec_dict
      assert ot_class2.train_match_set == self.m_set
      assert ot_class2.train_non_match_set == self.nm_set
      for i in range(5):
        assert ot_class2.opt_threshold_list[i] >= 0.0
        assert ot_class2.opt_threshold_list[i] <= 1.0

      test_res = ot_class2.test(self.w_vec_dict,self.m_set,self.nm_set)
      assert len(test_res) == 4
      assert ot_class2.bin_width == bw
      assert ot_class2.min_method == mm
      assert ot_class2.train_w_vec_dict == self.w_vec_dict
      assert ot_class2.train_match_set == self.m_set
      assert ot_class2.train_non_match_set == self.nm_set

      test_res = ot_class2.cross_validate(self.w_vec_dict, self.m_set,
                                          self.nm_set, 2)
      assert len(test_res) == 4
      assert ot_class2.bin_width == bw
      assert ot_class2.min_method == mm
      for i in range(5):
        assert ot_class2.opt_threshold_list[i] >= 0.0
        assert ot_class2.opt_threshold_list[i] <= 1.0

      test_res = ot_class2.cross_validate(self.w_vec_dict, self.m_set,
                                         self.nm_set, 3)
      assert len(test_res) == 4
      assert ot_class2.bin_width == bw
      assert ot_class2.min_method == mm
      for i in range(5):
        assert ot_class2.opt_threshold_list[i] >= 0.0
        assert ot_class2.opt_threshold_list[i] <= 1.0

      test_res = ot_class2.cross_validate(self.w_vec_dict, self.m_set,
                                         self.nm_set, 5)
      assert len(test_res) == 4
      assert ot_class2.bin_width == bw
      assert ot_class2.min_method == mm
      for i in range(5):
        assert ot_class2.opt_threshold_list[i] >= 0.0
        assert ot_class2.opt_threshold_list[i] <= 1.0

      class_res = ot_class2.classify(self.test_w_vec_dict)
      assert len(class_res) == 3
      assert isinstance(class_res[0], set) == True
      assert isinstance(class_res[1], set) == True
      assert isinstance(class_res[2], set) == True
      assert len(class_res[0]) + len(class_res[1]) + len(class_res[2]) == \
             len(self.test_w_vec_dict)
      for i in range(5):
        assert ot_class2.opt_threshold_list[i] >= 0.0
        assert ot_class2.opt_threshold_list[i] <= 1.0


  def testKMeansClassifier(self):  # - - - - - - - - - - - - - - - - - - - - -
    """Test K-means classifier"""

    # Loop over various distance measure, sampling, centroid initialisation
    # methods and fuzzy region thresholds
    #
    for (dm, sr, ci,frt) in [(mymath.distL1,   100,  'min/max', None),
                             (mymath.distL1,   100,  'random',  None),
                             (mymath.distL2,   100,  'min/max', None),
                             (mymath.distL2,   100,  'random',  None),
                             (mymath.distCanberra, 100,  'min/max', None),
                             (mymath.distCanberra, 100,  'random',  None),
                             (mymath.distCosine, 100,  'min/max', None),
                             (mymath.distCosine, 100,  'random',  None),
                             (mymath.distL1,   10.0, 'min/max', None),
                             (mymath.distL1,   10.0, 'random',  None),
                             (mymath.distL2,   10.0, 'min/max', None),
                             (mymath.distL2,   10.0, 'random',  None),
                             (mymath.distCanberra, 10.0, 'min/max', None),
                             (mymath.distCanberra, 10.0, 'random',  None),
                             (mymath.distCosine, 10.0, 'min/max', None),
                             (mymath.distCosine, 10.0, 'random',  None),
                             (mymath.distL1,   99.9, 'min/max', None),
                             (mymath.distL1,   99.9, 'random',  None),
                             (mymath.distL2,   99.9, 'min/max', None),
                             (mymath.distL2,   99.9, 'random',  None),
                             (mymath.distCanberra, 99.9, 'min/max', None),
                             (mymath.distCanberra, 99.9, 'random',  None),
                             (mymath.distCosine, 99.9, 'min/max', None),
                             (mymath.distCosine, 99.9, 'random',  None),
                             (mymath.distL1,   100,  'min/max', 0.1),
                             (mymath.distL1,   100,  'random',  0.1),
                             (mymath.distL2,   100,  'min/max', 0.1),
                             (mymath.distL2,   100,  'random',  0.1),
                             (mymath.distCanberra, 100,  'min/max', 0.1),
                             (mymath.distCanberra, 100,  'random',  0.1),
                             (mymath.distCosine, 100,  'min/max', 0.1),
                             (mymath.distCosine, 100,  'random',  0.1),
                             (mymath.distL1,   10.0, 'min/max', 0.1),
                             (mymath.distL1,   10.0, 'random',  0.1),
                             (mymath.distL2,   10.0, 'min/max', 0.1),
                             (mymath.distL2,   10.0, 'random',  0.1),
                             (mymath.distCanberra, 10.0, 'min/max', 0.1),
                             (mymath.distCanberra, 10.0, 'random',  0.1),
                             (mymath.distCosine, 10.0, 'min/max', 0.1),
                             (mymath.distCosine, 10.0, 'random',  0.1),
                             (mymath.distL1,   99.9, 'min/max', 0.1),
                             (mymath.distL1,   99.9, 'random',  0.1),
                             (mymath.distL2,   99.9, 'min/max', 0.1),
                             (mymath.distL2,   99.9, 'random',  0.1),
                             (mymath.distCanberra, 99.9, 'min/max', 0.1),
                             (mymath.distCanberra, 99.9, 'random',  0.1),
                             (mymath.distCosine, 99.9, 'min/max', 0.1),
                             (mymath.distCosine, 99.9, 'random',  0.1),
                             (mymath.distL1,   100,  'min/max', 0.4),
                             (mymath.distL1,   100,  'random',  0.4),
                             (mymath.distL2,   100,  'min/max', 0.4),
                             (mymath.distL2,   100,  'random',  0.4),
                             (mymath.distCanberra, 100,  'min/max', 0.4),
                             (mymath.distCanberra, 100,  'random',  0.4),
                             (mymath.distCosine, 100,  'min/max', 0.4),
                             (mymath.distCosine, 100,  'random',  0.4),
                             (mymath.distL1,   10.0, 'min/max', 0.4),
                             (mymath.distL1,   10.0, 'random',  0.4),
                             (mymath.distL2,   10.0, 'min/max', 0.4),
                             (mymath.distL2,   10.0, 'random',  0.4),
                             (mymath.distCanberra, 10.0, 'min/max', 0.4),
                             (mymath.distCanberra, 10.0, 'random',  0.4),
                             (mymath.distCosine, 10.0, 'min/max', 0.4),
                             (mymath.distCosine, 10.0, 'random',  0.4),
                             (mymath.distL1,   99.9, 'min/max', 0.4),
                             (mymath.distL1,   99.9, 'random',  0.4),
                             (mymath.distL2,   99.9, 'min/max', 0.4),
                             (mymath.distL2,   99.9, 'random',  0.4),
                             (mymath.distCanberra, 99.9, 'min/max', 0.4),
                             (mymath.distCanberra, 99.9, 'random',  0.4),
                             (mymath.distCosine, 99.9, 'min/max', 0.4),
                             (mymath.distCosine, 99.9, 'random',  0.4)]:

#      pass
#    for x in []:

      km_class = classification.KMeans(descr = 'K-means',
                                       sample = sr,
                                       max_iter_count = 100,
                                       dist_me = dm,
                                       centroid_i = ci,
                                       fuzz_reg_thres = frt)
      assert km_class.sample == sr
      assert km_class.max_iter_count == 100
      assert km_class.dist_measure == dm
      assert km_class.centroid_init == ci
      assert km_class.fuzz_reg_thres == frt
      assert km_class.train_w_vec_dict == None
      assert km_class.train_match_set == None
      assert km_class.train_non_match_set == None

      km_class.train(self.w_vec_dict, self.m_set, self.nm_set)
      assert km_class.sample == sr
      assert km_class.max_iter_count == 100
      assert km_class.dist_measure == dm
      assert km_class.centroid_init == ci
      assert km_class.fuzz_reg_thres == frt
      assert km_class.train_w_vec_dict == self.w_vec_dict
      assert km_class.train_match_set == self.m_set
      assert km_class.train_non_match_set == self.nm_set

      for i in range(5):
        assert km_class.m_centroid[i] >= 0.0
        assert km_class.m_centroid[i] <= 1.0
        assert km_class.nm_centroid[i] >= 0.0
        assert km_class.nm_centroid[i] <= 1.0

      test_res = km_class.test(self.w_vec_dict,self.m_set,self.nm_set)
      assert len(test_res) == 4
      assert km_class.sample == sr
      assert km_class.max_iter_count == 100
      assert km_class.dist_measure == dm
      assert km_class.centroid_init == ci
      assert km_class.fuzz_reg_thres == frt
      assert km_class.train_w_vec_dict == self.w_vec_dict
      assert km_class.train_match_set == self.m_set
      assert km_class.train_non_match_set == self.nm_set

      test_res = km_class.cross_validate(self.w_vec_dict, self.m_set,
                                         self.nm_set, 2)
      assert len(test_res) == 4
      assert km_class.sample == sr
      assert km_class.fuzz_reg_thres == frt
      assert km_class.max_iter_count == 100
      assert km_class.dist_measure == dm
      assert km_class.centroid_init == ci
      for i in range(5):
        assert km_class.m_centroid[i] >= 0.0
        assert km_class.m_centroid[i] <= 1.0
        assert km_class.nm_centroid[i] >= 0.0
        assert km_class.nm_centroid[i] <= 1.0

      test_res = km_class.cross_validate(self.w_vec_dict, self.m_set,
                                         self.nm_set, 3)
      assert len(test_res) == 4
      assert km_class.sample == sr
      assert km_class.fuzz_reg_thres == frt
      assert km_class.max_iter_count == 100
      assert km_class.dist_measure == dm
      assert km_class.centroid_init == ci
      for i in range(5):
        assert km_class.m_centroid[i] >= 0.0
        assert km_class.m_centroid[i] <= 1.0
        assert km_class.nm_centroid[i] >= 0.0
        assert km_class.nm_centroid[i] <= 1.0

      test_res = km_class.cross_validate(self.w_vec_dict, self.m_set,
                                         self.nm_set, 5)
      assert len(test_res) == 4
      assert km_class.sample == sr
      assert km_class.fuzz_reg_thres == frt
      assert km_class.max_iter_count == 100
      assert km_class.dist_measure == dm
      assert km_class.centroid_init == ci
      for i in range(5):
        assert km_class.m_centroid[i] >= 0.0
        assert km_class.m_centroid[i] <= 1.0
        assert km_class.nm_centroid[i] >= 0.0
        assert km_class.nm_centroid[i] <= 1.0

      class_res = km_class.classify(self.test_w_vec_dict)
      assert len(class_res) == 3
      assert isinstance(class_res[0], set) == True
      assert isinstance(class_res[1], set) == True
      assert isinstance(class_res[2], set) == True
      assert len(class_res[0]) + len(class_res[1]) + len(class_res[2]) == \
             len(self.test_w_vec_dict)
      for i in range(5):
        assert km_class.m_centroid[i] >= 0.0
        assert km_class.m_centroid[i] <= 1.0
        assert km_class.nm_centroid[i] >= 0.0
        assert km_class.nm_centroid[i] <= 1.0

      # Train classifier when initialising it - - - - - - - - - - - - - - - - -
      #
      km_class2 = classification.KMeans(descr = 'K-means',
                                       sample = sr,
                                       max_iter_count = 10,
                                       dist_me = dm,
                                       centroid_i = ci,
                                       fuzz_reg_thres = frt,
                                       train_w_vec_dict = self.w_vec_dict,
                                       train_match_set = self.m_set,
                                       train_non_match_set = self.nm_set)
      assert km_class2.sample == sr
      assert km_class2.max_iter_count == 10
      assert km_class2.dist_measure == dm
      assert km_class2.centroid_init == ci
      assert km_class2.fuzz_reg_thres == frt
      assert km_class2.train_w_vec_dict == self.w_vec_dict
      assert km_class2.train_match_set == self.m_set
      assert km_class2.train_non_match_set == self.nm_set
      for i in range(5):
        assert km_class2.m_centroid[i] >= 0.0
        assert km_class2.m_centroid[i] <= 1.0
        assert km_class2.nm_centroid[i] >= 0.0
        assert km_class2.nm_centroid[i] <= 1.0

      test_res = km_class2.test(self.w_vec_dict,self.m_set,self.nm_set)
      assert len(test_res) == 4
      assert km_class2.sample == sr
      assert km_class2.max_iter_count == 10
      assert km_class2.dist_measure == dm
      assert km_class2.centroid_init == ci
      assert km_class2.fuzz_reg_thres == frt
      assert km_class2.train_w_vec_dict == self.w_vec_dict
      assert km_class2.train_match_set == self.m_set
      assert km_class2.train_non_match_set == self.nm_set

      test_res = km_class2.cross_validate(self.w_vec_dict, self.m_set,
                                          self.nm_set, 2)
      assert len(test_res) == 4
      assert km_class2.sample == sr
      assert km_class2.fuzz_reg_thres == frt
      assert km_class2.max_iter_count == 10
      assert km_class2.dist_measure == dm
      assert km_class2.centroid_init == ci
      for i in range(5):
        assert km_class2.m_centroid[i] >= 0.0
        assert km_class2.m_centroid[i] <= 1.0
        assert km_class2.nm_centroid[i] >= 0.0
        assert km_class2.nm_centroid[i] <= 1.0

      test_res = km_class2.cross_validate(self.w_vec_dict, self.m_set,
                                         self.nm_set, 3)
      assert len(test_res) == 4
      assert km_class2.sample == sr
      assert km_class2.fuzz_reg_thres == frt
      assert km_class2.max_iter_count == 10
      assert km_class2.dist_measure == dm
      assert km_class2.centroid_init == ci
      for i in range(5):
        assert km_class2.m_centroid[i] >= 0.0
        assert km_class2.m_centroid[i] <= 1.0
        assert km_class2.nm_centroid[i] >= 0.0
        assert km_class2.nm_centroid[i] <= 1.0

      test_res = km_class2.cross_validate(self.w_vec_dict, self.m_set,
                                         self.nm_set, 5)
      assert len(test_res) == 4
      assert km_class2.sample == sr
      assert km_class2.fuzz_reg_thres == frt
      assert km_class2.max_iter_count == 10
      assert km_class2.dist_measure == dm
      assert km_class2.centroid_init == ci
      for i in range(5):
        assert km_class2.m_centroid[i] >= 0.0
        assert km_class2.m_centroid[i] <= 1.0
        assert km_class2.nm_centroid[i] >= 0.0
        assert km_class2.nm_centroid[i] <= 1.0

      class_res = km_class2.classify(self.test_w_vec_dict)
      assert len(class_res) == 3
      assert isinstance(class_res[0], set) == True
      assert isinstance(class_res[1], set) == True
      assert isinstance(class_res[2], set) == True
      assert len(class_res[0]) + len(class_res[1]) + len(class_res[2]) == \
             len(self.test_w_vec_dict)
      for i in range(5):
        assert km_class2.m_centroid[i] >= 0.0
        assert km_class2.m_centroid[i] <= 1.0
        assert km_class2.nm_centroid[i] >= 0.0
        assert km_class2.nm_centroid[i] <= 1.0


  def testFarthestFirstClassifier(self):  # - - - - - - - - - - - - - - - - - -
    """Test Farthest First classifier"""

    # Loop over various distance measure, sampling, centroid initialisation
    # methods and fuzzy region thresholds
    #
    for (dm, sr, ci,frt) in [(mymath.distL1,   100,  'min/max', None),
                             (mymath.distL1,   100,  'mode/max',  None),
                             (mymath.distL1,   100,  'traditional',  None),
                             (mymath.distL2,   100,  'min/max', None),
                             (mymath.distL2,   100,  'mode/max',  None),
                             (mymath.distL2,   100,  'traditional',  None),
                             (mymath.distCanberra, 100,  'min/max', None),
                             (mymath.distCanberra, 100,  'mode/max',  None),
                             (mymath.distCanberra, 100,  'traditional',  None),
                             (mymath.distCosine, 100,  'min/max', None),
                             (mymath.distCosine, 100,  'mode/max',  None),
                             (mymath.distCosine, 100,  'traditional',  None),
                             (mymath.distL1,   10.0, 'min/max', None),
                             (mymath.distL1,   10.0, 'mode/max',  None),
                             (mymath.distL1,   10.0, 'traditional',  None),
                             (mymath.distL2,   10.0, 'min/max', None),
                             (mymath.distL2,   10.0, 'mode/max',  None),
                             (mymath.distL2,   10.0, 'traditional',  None),
                             (mymath.distCanberra, 10.0, 'min/max', None),
                             (mymath.distCanberra, 10.0, 'mode/max',  None),
                             (mymath.distCanberra, 10.0, 'traditional',  None),
                             (mymath.distCosine, 10.0, 'min/max', None),
                             (mymath.distCosine, 10.0, 'mode/max',  None),
                             (mymath.distCosine, 10.0, 'traditional',  None),
                             (mymath.distL1,   99.9, 'min/max', None),
                             (mymath.distL1,   99.9, 'mode/max',  None),
                             (mymath.distL1,   99.9, 'traditional',  None),
                             (mymath.distL2,   99.9, 'min/max', None),
                             (mymath.distL2,   99.9, 'mode/max',  None),
                             (mymath.distL2,   99.9, 'traditional',  None),
                             (mymath.distCanberra, 99.9, 'min/max', None),
                             (mymath.distCanberra, 99.9, 'mode/max',  None),
                             (mymath.distCanberra, 99.9, 'traditional',  None),
                             (mymath.distCosine, 99.9, 'min/max', None),
                             (mymath.distCosine, 99.9, 'mode/max',  None),
                             (mymath.distCosine, 99.9, 'traditional',  None),
                             (mymath.distL1,   100,  'min/max', 0.1),
                             (mymath.distL1,   100,  'mode/max',  0.1),
                             (mymath.distL1,   100,  'traditional',  0.1),
                             (mymath.distL2,   100,  'min/max', 0.1),
                             (mymath.distL2,   100,  'mode/max',  0.1),
                             (mymath.distL2,   100,  'traditional',  0.1),
                             (mymath.distCanberra, 100,  'min/max', 0.1),
                             (mymath.distCanberra, 100,  'mode/max',  0.1),
                             (mymath.distCanberra, 100,  'traditional',  0.1),
                             (mymath.distCosine, 100,  'min/max', 0.1),
                             (mymath.distCosine, 100,  'mode/max',  0.1),
                             (mymath.distCosine, 100,  'traditional',  0.1),
                             (mymath.distL1,   10.0, 'min/max', 0.1),
                             (mymath.distL1,   10.0, 'mode/max',  0.1),
                             (mymath.distL1,   10.0, 'traditional',  0.1),
                             (mymath.distL2,   10.0, 'min/max', 0.1),
                             (mymath.distL2,   10.0, 'mode/max',  0.1),
                             (mymath.distL2,   10.0, 'traditional',  0.1),
                             (mymath.distCanberra, 10.0, 'min/max', 0.1),
                             (mymath.distCanberra, 10.0, 'mode/max',  0.1),
                             (mymath.distCanberra, 10.0, 'traditional',  0.1),
                             (mymath.distCosine, 10.0, 'min/max', 0.1),
                             (mymath.distCosine, 10.0, 'mode/max',  0.1),
                             (mymath.distCosine, 10.0, 'traditional',  0.1),
                             (mymath.distL1,   99.9, 'min/max', 0.1),
                             (mymath.distL1,   99.9, 'mode/max',  0.1),
                             (mymath.distL1,   99.9, 'traditional',  0.1),
                             (mymath.distL2,   99.9, 'min/max', 0.1),
                             (mymath.distL2,   99.9, 'mode/max',  0.1),
                             (mymath.distL2,   99.9, 'traditional',  0.1),
                             (mymath.distCanberra, 99.9, 'min/max', 0.1),
                             (mymath.distCanberra, 99.9, 'mode/max',  0.1),
                             (mymath.distCanberra, 99.9, 'traditional',  0.1),
                             (mymath.distCosine, 99.9, 'min/max', 0.1),
                             (mymath.distCosine, 99.9, 'mode/max',  0.1),
                             (mymath.distCosine, 99.9, 'traditional',  0.1),
                             (mymath.distL1,   100,  'min/max', 0.4),
                             (mymath.distL1,   100,  'mode/max',  0.4),
                             (mymath.distL1,   100,  'traditional',  0.4),
                             (mymath.distL2,   100,  'min/max', 0.4),
                             (mymath.distL2,   100,  'mode/max',  0.4),
                             (mymath.distL2,   100,  'traditional',  0.4),
                             (mymath.distCanberra, 100,  'min/max', 0.4),
                             (mymath.distCanberra, 100,  'mode/max',  0.4),
                             (mymath.distCanberra, 100,  'traditional',  0.4),
                             (mymath.distCosine, 100,  'min/max', 0.4),
                             (mymath.distCosine, 100,  'mode/max',  0.4),
                             (mymath.distCosine, 100,  'traditional',  0.4),
                             (mymath.distL1,   10.0, 'min/max', 0.4),
                             (mymath.distL1,   10.0, 'mode/max',  0.4),
                             (mymath.distL1,   10.0, 'traditional',  0.4),
                             (mymath.distL2,   10.0, 'min/max', 0.4),
                             (mymath.distL2,   10.0, 'mode/max',  0.4),
                             (mymath.distL2,   10.0, 'traditional',  0.4),
                             (mymath.distCanberra, 10.0, 'min/max', 0.4),
                             (mymath.distCanberra, 10.0, 'mode/max',  0.4),
                             (mymath.distCanberra, 10.0, 'traditional',  0.4),
                             (mymath.distCosine, 10.0, 'min/max', 0.4),
                             (mymath.distCosine, 10.0, 'mode/max',  0.4),
                             (mymath.distCosine, 10.0, 'traditional',  0.4),
                             (mymath.distL1,   99.9, 'min/max', 0.4),
                             (mymath.distL1,   99.9, 'mode/max',  0.4),
                             (mymath.distL1,   99.9, 'traditional',  0.4),
                             (mymath.distL2,   99.9, 'min/max', 0.4),
                             (mymath.distL2,   99.9, 'mode/max',  0.4),
                             (mymath.distL2,   99.9, 'traditional',  0.4),
                             (mymath.distCanberra, 99.9, 'min/max', 0.4),
                             (mymath.distCanberra, 99.9, 'mode/max',  0.4),
                             (mymath.distCanberra, 99.9, 'traditional',  0.4),
                             (mymath.distCosine, 99.9, 'min/max', 0.4),
                             (mymath.distCosine, 99.9, 'mode/max',  0.4),
                             (mymath.distCosine, 99.9, 'traditional',  0.4)]:
#      pass
#    for x in []:

      ff_class = classification.FarthestFirst(descr = 'Farthest first',
                                              sample = sr,
                                              dist_me = dm,
                                              centroid_i = ci,
                                              fuzz_reg_thres = frt)
      assert ff_class.sample == sr
      assert ff_class.dist_measure == dm
      assert ff_class.centroid_init == ci
      assert ff_class.fuzz_reg_thres == frt
      assert ff_class.train_w_vec_dict == None
      assert ff_class.train_match_set == None
      assert ff_class.train_non_match_set == None

      ff_class.train(self.w_vec_dict, self.m_set, self.nm_set)
      assert ff_class.sample == sr
      assert ff_class.dist_measure == dm
      assert ff_class.centroid_init == ci
      assert ff_class.fuzz_reg_thres == frt
      assert ff_class.train_w_vec_dict == self.w_vec_dict
      assert ff_class.train_match_set == self.m_set
      assert ff_class.train_non_match_set == self.nm_set
      for i in range(5):
        assert ff_class.m_centroid[i] >= 0.0
        assert ff_class.m_centroid[i] <= 1.0
        assert ff_class.nm_centroid[i] >= 0.0
        assert ff_class.nm_centroid[i] <= 1.0

      test_res = ff_class.test(self.w_vec_dict,self.m_set,self.nm_set)
      assert len(test_res) == 4
      assert ff_class.sample == sr
      assert ff_class.dist_measure == dm
      assert ff_class.centroid_init == ci
      assert ff_class.fuzz_reg_thres == frt
      assert ff_class.train_w_vec_dict == self.w_vec_dict
      assert ff_class.train_match_set == self.m_set
      assert ff_class.train_non_match_set == self.nm_set

      test_res = ff_class.cross_validate(self.w_vec_dict, self.m_set,
                                         self.nm_set, 2)
      assert len(test_res) == 4
      assert ff_class.sample == sr
      assert ff_class.fuzz_reg_thres == frt
      assert ff_class.dist_measure == dm
      assert ff_class.centroid_init == ci
      for i in range(5):
        assert ff_class.m_centroid[i] >= 0.0
        assert ff_class.m_centroid[i] <= 1.0
        assert ff_class.nm_centroid[i] >= 0.0
        assert ff_class.nm_centroid[i] <= 1.0

      test_res = ff_class.cross_validate(self.w_vec_dict, self.m_set,
                                         self.nm_set, 3)
      assert len(test_res) == 4
      assert ff_class.sample == sr
      assert ff_class.fuzz_reg_thres == frt
      assert ff_class.dist_measure == dm
      assert ff_class.centroid_init == ci
      for i in range(5):
        assert ff_class.m_centroid[i] >= 0.0
        assert ff_class.m_centroid[i] <= 1.0
        assert ff_class.nm_centroid[i] >= 0.0
        assert ff_class.nm_centroid[i] <= 1.0

      test_res = ff_class.cross_validate(self.w_vec_dict, self.m_set,
                                         self.nm_set, 5)
      assert len(test_res) == 4
      assert ff_class.sample == sr
      assert ff_class.fuzz_reg_thres == frt
      assert ff_class.dist_measure == dm
      assert ff_class.centroid_init == ci
      for i in range(5):
        assert ff_class.m_centroid[i] >= 0.0
        assert ff_class.m_centroid[i] <= 1.0
        assert ff_class.nm_centroid[i] >= 0.0
        assert ff_class.nm_centroid[i] <= 1.0

      class_res = ff_class.classify(self.test_w_vec_dict)
      assert len(class_res) == 3
      assert isinstance(class_res[0], set) == True
      assert isinstance(class_res[1], set) == True
      assert isinstance(class_res[2], set) == True
      assert len(class_res[0]) + len(class_res[1]) + len(class_res[2]) == \
             len(self.test_w_vec_dict)
      for i in range(5):
        assert ff_class.m_centroid[i] >= 0.0
        assert ff_class.m_centroid[i] <= 1.0
        assert ff_class.nm_centroid[i] >= 0.0
        assert ff_class.nm_centroid[i] <= 1.0

      # Train classifier when initialising it - - - - - - - - - - - - - - - - -
      #
      ff_class2 = classification.FarthestFirst(descr = 'Farthest first',
                                       sample = sr,
                                       dist_me = dm,
                                       centroid_i = ci,
                                       fuzz_reg_thres = frt,
                                       train_w_vec_dict = self.w_vec_dict,
                                       train_match_set = self.m_set,
                                       train_non_match_set = self.nm_set)
      assert ff_class2.sample == sr
      assert ff_class2.dist_measure == dm
      assert ff_class2.centroid_init == ci
      assert ff_class2.fuzz_reg_thres == frt
      assert ff_class2.train_w_vec_dict == self.w_vec_dict
      assert ff_class2.train_match_set == self.m_set
      assert ff_class2.train_non_match_set == self.nm_set
      for i in range(5):
        assert ff_class2.m_centroid[i] >= 0.0
        assert ff_class2.m_centroid[i] <= 1.0
        assert ff_class2.nm_centroid[i] >= 0.0
        assert ff_class2.nm_centroid[i] <= 1.0

      test_res = ff_class2.test(self.w_vec_dict,self.m_set,self.nm_set)
      assert len(test_res) == 4
      assert ff_class2.sample == sr
      assert ff_class2.dist_measure == dm
      assert ff_class2.centroid_init == ci
      assert ff_class2.fuzz_reg_thres == frt
      assert ff_class2.train_w_vec_dict == self.w_vec_dict
      assert ff_class2.train_match_set == self.m_set
      assert ff_class2.train_non_match_set == self.nm_set

      test_res = ff_class2.cross_validate(self.w_vec_dict, self.m_set,
                                          self.nm_set, 2)
      assert len(test_res) == 4
      assert ff_class2.sample == sr
      assert ff_class2.fuzz_reg_thres == frt
      assert ff_class2.dist_measure == dm
      assert ff_class2.centroid_init == ci
      for i in range(5):
        assert ff_class2.m_centroid[i] >= 0.0
        assert ff_class2.m_centroid[i] <= 1.0
        assert ff_class2.nm_centroid[i] >= 0.0
        assert ff_class2.nm_centroid[i] <= 1.0

      test_res = ff_class2.cross_validate(self.w_vec_dict, self.m_set,
                                         self.nm_set, 3)
      assert len(test_res) == 4
      assert ff_class2.sample == sr
      assert ff_class2.fuzz_reg_thres == frt
      assert ff_class2.dist_measure == dm
      assert ff_class2.centroid_init == ci
      for i in range(5):
        assert ff_class2.m_centroid[i] >= 0.0
        assert ff_class2.m_centroid[i] <= 1.0
        assert ff_class2.nm_centroid[i] >= 0.0
        assert ff_class2.nm_centroid[i] <= 1.0

      assert len(test_res) == 4
      assert ff_class2.sample == sr
      assert ff_class2.fuzz_reg_thres == frt
      assert ff_class2.dist_measure == dm
      assert ff_class2.centroid_init == ci
      for i in range(5):
        assert ff_class2.m_centroid[i] >= 0.0
        assert ff_class2.m_centroid[i] <= 1.0
        assert ff_class2.nm_centroid[i] >= 0.0
        assert ff_class2.nm_centroid[i] <= 1.0

      class_res = ff_class2.classify(self.test_w_vec_dict)
      assert len(class_res) == 3
      assert isinstance(class_res[0], set) == True
      assert isinstance(class_res[1], set) == True
      assert isinstance(class_res[2], set) == True
      assert len(class_res[0]) + len(class_res[1]) + len(class_res[2]) == \
             len(self.test_w_vec_dict)
      for i in range(5):
        assert ff_class2.m_centroid[i] >= 0.0
        assert ff_class2.m_centroid[i] <= 1.0
        assert ff_class2.nm_centroid[i] >= 0.0
        assert ff_class2.nm_centroid[i] <= 1.0


  def testSuppVecMachineClassifier(self):  # - - - - - - - - - - - - - - - -
    """Test SVM classifier"""

    # Loop over various SVM and sample values
    #
    for (svmk, svmC, sr) in [('LINEAR',  10,  100),
                             ('LINEAR',  10,  10.0),
                             ('LINEAR',  10,  99.9),
                             ('LINEAR',  1,   100),
                             ('LINEAR',  1,   10.0),
                             ('LINEAR',  1,   99.9),
                             ('LINEAR',  0.1, 100),
                             ('LINEAR',  0.1, 10.0),
                             ('LINEAR',  0.1, 99.9),
                             ('POLY',    10,  100),
                             ('POLY',    10,  10.0),
                             ('POLY',    10,  99.9),
                             ('POLY',    1,   100),
                             ('POLY',    1,   10.0),
                             ('POLY',    1,   99.9),
                             ('POLY',    0.1, 100),
                             ('POLY',    0.1, 10.0),
                             ('POLY',    0.1, 99.9),
                             ('RBF',     10,  100),
                             ('RBF',     10,  10.0),
                             ('RBF',     10,  99.9),
                             ('RBF',     1,   100),
                             ('RBF',     1,   10.0),
                             ('RBF',     1,   99.9),
                             ('RBF',     0.1, 100),
                             ('RBF',     0.1, 10.0),
                             ('RBF',     0.1, 99.9),
                             ('SIGMOID', 10,  100),
                             ('SIGMOID', 10,  10.0),
                             ('SIGMOID', 10,  99.9),
                             ('SIGMOID', 1,   100),
                             ('SIGMOID', 1,   10.0),
                             ('SIGMOID', 1,   99.9),
                             ('SIGMOID', 0.1, 100),
                             ('SIGMOID', 0.1, 10.0),
                             ('SIGMOID', 0.1, 99.9)]:
#      pass
#    for x in []:

      svm_class = classification.SuppVecMachine(descr = 'SVM',
                                                sample = sr,
                                                kernel_type = svmk,
                                                C = svmC)
      assert svm_class.sample == sr
      assert svm_class.C == svmC
      assert svm_class.kernel_type == svmk
      assert svm_class.train_w_vec_dict == None
      assert svm_class.train_match_set == None
      assert svm_class.train_non_match_set == None

      svm_class.train(self.w_vec_dict, self.m_set, self.nm_set)
      assert svm_class.sample == sr
      assert svm_class.C == svmC
      assert svm_class.kernel_type == svmk
      assert svm_class.train_w_vec_dict == self.w_vec_dict
      assert svm_class.train_match_set == self.m_set
      assert svm_class.train_non_match_set == self.nm_set

      test_res = svm_class.test(self.w_vec_dict,self.m_set,self.nm_set)
      assert len(test_res) == 4
      assert svm_class.sample == sr
      assert svm_class.C == svmC
      assert svm_class.kernel_type == svmk
      assert svm_class.train_w_vec_dict == self.w_vec_dict
      assert svm_class.train_match_set == self.m_set
      assert svm_class.train_non_match_set == self.nm_set

      test_res = svm_class.cross_validate(self.w_vec_dict, self.m_set,
                                          self.nm_set, 2)
      assert len(test_res) == 4
      assert svm_class.sample == sr
      assert svm_class.C == svmC
      assert svm_class.kernel_type == svmk
      #assert svm_class.train_w_vec_dict == self.w_vec_dict
      #assert svm_class.train_match_set == self.m_set
      #assert svm_class.train_non_match_set == self.nm_set

      test_res = svm_class.cross_validate(self.w_vec_dict, self.m_set,
                                          self.nm_set, 3)
      assert len(test_res) == 4
      assert svm_class.sample == sr
      assert svm_class.C == svmC
      assert svm_class.kernel_type == svmk
      #assert svm_class.train_w_vec_dict == self.w_vec_dict
      #assert svm_class.train_match_set == self.m_set
      #assert svm_class.train_non_match_set == self.nm_set

      test_res = svm_class.cross_validate(self.w_vec_dict, self.m_set,
                                          self.nm_set, 5)
      assert len(test_res) == 4
      assert svm_class.sample == sr
      assert svm_class.C == svmC
      assert svm_class.kernel_type == svmk
      #assert svm_class.train_w_vec_dict == self.w_vec_dict
      #assert svm_class.train_match_set == self.m_set
      #assert svm_class.train_non_match_set == self.nm_set

      class_res = svm_class.classify(self.test_w_vec_dict)
      assert len(class_res) == 3
      assert isinstance(class_res[0], set) == True
      assert isinstance(class_res[1], set) == True
      assert isinstance(class_res[2], set) == True
      assert len(class_res[0]) + len(class_res[1]) + len(class_res[2]) == \
             len(self.test_w_vec_dict)

      # Train classifier when initialising it - - - - - - - - - - - - - - - - -
      #
      svm_class2 = classification.SuppVecMachine(descr = 'SVM',
                                       sample = sr,
                                       kernel_type = svmk,
                                       C = svmC,
                                       train_w_vec_dict = self.w_vec_dict,
                                       train_match_set = self.m_set,
                                       train_non_match_set = self.nm_set)
      assert svm_class2.sample == sr
      assert svm_class2.C == svmC
      assert svm_class2.kernel_type == svmk
      assert svm_class2.train_w_vec_dict == self.w_vec_dict
      assert svm_class2.train_match_set == self.m_set
      assert svm_class2.train_non_match_set == self.nm_set

      test_res = svm_class2.test(self.w_vec_dict,self.m_set,self.nm_set)
      assert len(test_res) == 4
      assert svm_class2.sample == sr
      assert svm_class2.C == svmC
      assert svm_class2.kernel_type == svmk
      assert svm_class2.train_w_vec_dict == self.w_vec_dict
      assert svm_class2.train_match_set == self.m_set
      assert svm_class2.train_non_match_set == self.nm_set

      test_res = svm_class2.cross_validate(self.w_vec_dict, self.m_set,
                                          self.nm_set, 2)
      assert len(test_res) == 4
      assert svm_class2.sample == sr
      assert svm_class2.C == svmC
      assert svm_class2.kernel_type == svmk
      #assert svm_class2.train_w_vec_dict == self.w_vec_dict
      #assert svm_class2.train_match_set == self.m_set
      #assert svm_class.train_non_match_set == self.nm_set

      test_res = svm_class2.cross_validate(self.w_vec_dict, self.m_set,
                                          self.nm_set, 3)
      assert len(test_res) == 4
      assert svm_class2.sample == sr
      assert svm_class2.C == svmC
      assert svm_class2.kernel_type == svmk
      #assert svm_class2.train_w_vec_dict == self.w_vec_dict
      #assert svm_class2.train_match_set == self.m_set
      #assert svm_class2.train_non_match_set == self.nm_set

      test_res = svm_class.cross_validate(self.w_vec_dict, self.m_set,
                                          self.nm_set, 5)
      assert len(test_res) == 4
      assert svm_class2.sample == sr
      assert svm_class2.C == svmC
      assert svm_class2.kernel_type == svmk
      #assert svm_class2.train_w_vec_dict == self.w_vec_dict
      #assert svm_class2.train_match_set == self.m_set
      #assert svm_class2.train_non_match_set == self.nm_set

      class_res = svm_class2.classify(self.test_w_vec_dict)
      assert len(class_res) == 3
      assert isinstance(class_res[0], set) == True
      assert isinstance(class_res[1], set) == True
      assert isinstance(class_res[2], set) == True
      assert len(class_res[0]) + len(class_res[1]) + len(class_res[2]) == \
             len(self.test_w_vec_dict)


  def testTwoStepClassifier(self):  # - - - - - - - - - - - - - - - - - - - - -
    """Test Two-Step classifier"""

    # Loop over various parameter values
    #
    for (s1mm, s1nmm) in [((1.0,'threshold',0.1),    (0.0,'threshold',0.1)),
                          ((1.0,'threshold',0.3),    (0.0,'threshold',0.3)),
                          ((1.0,'threshold',0.6),    (0.0,'threshold',0.6)),
                          ((1.0,'nearest',10,False), (0.0,'nearest',10,False)),
                          ((1.0,'nearest',20,False), (0.0,'nearest',20,False)),
                          ((1.0,'nearest',50,False), (0.0,'nearest',20,False)),
                          ((1.0,'nearest',10,True),  (0.0,'nearest',10,True)),
                          ((1.0,'nearest',20,True),  (0.0,'nearest',20,True)),
                          ((1.0,'nearest',50,True),  (0.0,'nearest',20,True)),
                          ((1.0,'threshold',0.4),    (0.0,'nearest',10,True)),
                          ((1.0,'nearest',20,True),  (0.0,'threshold',0.6)),
                          ((1.0,'threshold',0.4),    (0.0,'nearest',20,False)),
                          ((1.0,'nearest',20,False), (0.0,'threshold',0.6))]:

      for (rs, s2c) in [(None, ('svm', 'LINEAR', 5, 10, 50)),
                        (None, ('kmeans', mymath.distCanberra, 1000)),
                        (None, ('kmeans', mymath.distCosine, 1000)),
                        (None, ('nn', mymath.distL1, 1)),
                        (None, ('nn', mymath.distL1, 3)),
                        (None, ('nn', mymath.distL1, 7)),

                        (('uniform',10,10), ('svm', 'LINEAR', 10, 20, 60)),
                        (('uniform',10,10), ('kmeans', mymath.distL1, 1000)),
                        (('uniform',10,10), ('nn', mymath.distL1, 1)),
                        (('uniform',10,10), ('nn', mymath.distL2, 1)),
                        (('uniform',5,20), ('svm', 'POLY', 0.1, 10, 50)),
                        (('uniform',5,20), ('kmeans', mymath.distL2, 0)),
                        (('uniform',20,5), ('svm', 'RBF', 0.5, 0, 0)),
                        (('uniform',20,5), ('kmeans', mymath.distLInf, 1000)),
                        (('linear',10,10), ('svm', 'SIGMOID',1, 20, 60)),
                        (('linear',10,10), ('kmeans', mymath.distCanberra,
                                            1000)),
                        (('linear',10,10), ('kmeans', mymath.distCosine,1000)),
                        (('linear',10,50), ('svm', 'LINEAR', 0.01, 0, 0)),
                        (('linear',10,50), ('kmeans', mymath.distL1, 0)),
                        (('linear',10,50), ('nn', mymath.distLInf, 3)),
                        (('linear',10,50), ('nn', mymath.distCanberra, 3)),
                        (('linear',10,50), ('nn', mymath.distCosine, 3)),
                        (('linear',50,30), ('svm', 'POLY', 12, 25, 70)),
                        (('linear',40,40), ('kmeans', mymath.distCanberra,
                                            1000)),
                        (('linear',40,40), ('kmeans', mymath.distCosine,1000)),
                        (('exponential',10,10), ('svm', 'RBF', 99, 30, 80)),
                        (('exponential',10,10), ('kmeans', mymath.distLInf,
                                                 1000)),
                        (('exponential',10,50), ('svm', 'SIGMOID', 0.97,0,0)),
                        (('exponential',10,50), ('kmeans',mymath.distL1,1000)),
                        (('exponential',30,50), ('svm', 'LINEAR',1.0, 5, 20)),
                        (('exponential',40,40), ('kmeans',mymath.distL2,0)),
                        (('exponential',40,40), ('nn',mymath.distL2,3)),
                        (('exponential',40,40), ('nn',mymath.distLInf,1))]:
#        pass
#      for x in []:

        ts_class = classification.TwoStep(descr = '2-step',
                                          s1_match_method = s1mm,
                                          s1_non_match_method = s1nmm,
                                          random_selection = rs,
                                          s2_classifier = s2c)
        assert ts_class.s1_m_method == s1mm
        assert ts_class.s1_nm_method == s1nmm
        assert ts_class.rand_sel == rs
        assert ts_class.s2_classifier == s2c
        assert ts_class.train_w_vec_dict == None
        assert ts_class.train_match_set == None
        assert ts_class.train_non_match_set == None

        ts_class.train(self.w_vec_dict, self.m_set, self.nm_set)
        assert ts_class.s1_m_method == s1mm
        assert ts_class.s1_nm_method == s1nmm
        assert ts_class.rand_sel == rs
        assert ts_class.s2_classifier == s2c
        assert ts_class.train_w_vec_dict == self.w_vec_dict
        assert ts_class.train_match_set == self.m_set
        assert ts_class.train_non_match_set == self.nm_set

        test_res = ts_class.test(self.w_vec_dict,self.m_set,self.nm_set)
        assert ts_class.s1_m_method == s1mm
        assert ts_class.s1_nm_method == s1nmm
        assert ts_class.rand_sel == rs
        assert ts_class.s2_classifier == s2c
        assert ts_class.train_w_vec_dict == self.w_vec_dict
        assert ts_class.train_match_set == self.m_set
        assert ts_class.train_non_match_set == self.nm_set

        test_res = ts_class.cross_validate(self.w_vec_dict, self.m_set,
                                           self.nm_set, 2)
        assert ts_class.s1_m_method == s1mm
        assert ts_class.s1_nm_method == s1nmm
        assert ts_class.rand_sel == rs
        assert ts_class.s2_classifier == s2c
        assert ts_class.train_w_vec_dict == self.w_vec_dict
        assert ts_class.train_match_set == self.m_set
        assert ts_class.train_non_match_set == self.nm_set

        test_res = ts_class.cross_validate(self.w_vec_dict, self.m_set,
                                           self.nm_set, 3)
        assert ts_class.s1_m_method == s1mm
        assert ts_class.s1_nm_method == s1nmm
        assert ts_class.rand_sel == rs
        assert ts_class.s2_classifier == s2c
        assert ts_class.train_w_vec_dict == self.w_vec_dict
        assert ts_class.train_match_set == self.m_set
        assert ts_class.train_non_match_set == self.nm_set

        test_res = ts_class.cross_validate(self.w_vec_dict, self.m_set,
                                           self.nm_set, 5)
        assert ts_class.s1_m_method == s1mm
        assert ts_class.s1_nm_method == s1nmm
        assert ts_class.rand_sel == rs
        assert ts_class.s2_classifier == s2c
        assert ts_class.train_w_vec_dict == self.w_vec_dict
        assert ts_class.train_match_set == self.m_set
        assert ts_class.train_non_match_set == self.nm_set

        class_res = ts_class.classify(self.test_w_vec_dict)
        assert len(class_res) == 3
        assert isinstance(class_res[0], set) == True
        assert isinstance(class_res[1], set) == True
        assert isinstance(class_res[2], set) == True
        assert len(class_res[0]) + len(class_res[1]) + len(class_res[2]) == \
               len(self.test_w_vec_dict)

        # Train classifier when initialising it - - - - - - - - - - - - - - - -
        #
        ts_class2 = classification.TwoStep(descr = '2-step',
                                           s1_match_method = s1mm,
                                           s1_non_match_method = s1nmm,
                                           random_selection = rs,
                                           s2_classifier = s2c,
                                           train_w_vec_dict = self.w_vec_dict,
                                           train_match_set = self.m_set,
                                           train_non_match_set = self.nm_set)
        assert ts_class2.s1_m_method == s1mm
        assert ts_class2.s1_nm_method == s1nmm
        assert ts_class2.rand_sel == rs
        assert ts_class2.s2_classifier == s2c
        assert ts_class2.train_w_vec_dict == self.w_vec_dict
        assert ts_class2.train_match_set == self.m_set
        assert ts_class2.train_non_match_set == self.nm_set

        test_res = ts_class2.test(self.w_vec_dict,self.m_set,self.nm_set)
        assert ts_class2.s1_m_method == s1mm
        assert ts_class2.s1_nm_method == s1nmm
        assert ts_class2.rand_sel == rs
        assert ts_class2.s2_classifier == s2c
        assert ts_class2.train_w_vec_dict == self.w_vec_dict
        assert ts_class2.train_match_set == self.m_set
        assert ts_class2.train_non_match_set == self.nm_set

        test_res = ts_class2.cross_validate(self.w_vec_dict, self.m_set,
                                           self.nm_set, 2)
        assert ts_class2.s1_m_method == s1mm
        assert ts_class2.s1_nm_method == s1nmm
        assert ts_class2.rand_sel == rs
        assert ts_class2.s2_classifier == s2c
        assert ts_class2.train_w_vec_dict == self.w_vec_dict
        assert ts_class2.train_match_set == self.m_set
        assert ts_class2.train_non_match_set == self.nm_set

        test_res = ts_class2.cross_validate(self.w_vec_dict, self.m_set,
                                           self.nm_set, 3)
        assert ts_class2.s1_m_method == s1mm
        assert ts_class2.s1_nm_method == s1nmm
        assert ts_class2.rand_sel == rs
        assert ts_class2.s2_classifier == s2c
        assert ts_class2.train_w_vec_dict == self.w_vec_dict
        assert ts_class2.train_match_set == self.m_set
        assert ts_class2.train_non_match_set == self.nm_set

        test_res = ts_class2.cross_validate(self.w_vec_dict, self.m_set,
                                           self.nm_set, 5)
        assert ts_class2.s1_m_method == s1mm
        assert ts_class2.s1_nm_method == s1nmm
        assert ts_class2.rand_sel == rs
        assert ts_class2.s2_classifier == s2c
        assert ts_class2.train_w_vec_dict == self.w_vec_dict
        assert ts_class2.train_match_set == self.m_set
        assert ts_class2.train_non_match_set == self.nm_set

        class_res = ts_class2.classify(self.test_w_vec_dict)
        assert len(class_res) == 3
        assert isinstance(class_res[0], set) == True
        assert isinstance(class_res[1], set) == True
        assert isinstance(class_res[2], set) == True
        assert len(class_res[0]) + len(class_res[1]) + len(class_res[2]) == \
               len(self.test_w_vec_dict)


  def testTAILORClassifier(self):  # - - - - - - - - - - - - - - - - - - - - -
    """Test TAILOR classifier"""

    # Loop over various k-means/SVM parameters and sample values
    #
    for (svmk, svmC, sr) in [('LINEAR',  10,  100),
                             ('LINEAR',  10,  10.0),
                             ('LINEAR',  10,  99.9),
                             ('LINEAR',  1,   100),
                             ('LINEAR',  1,   10.0),
                             ('LINEAR',  1,   99.9),
                             ('LINEAR',  0.1, 100),
                             ('LINEAR',  0.1, 10.0),
                             ('LINEAR',  0.1, 99.9),
                             ('POLY',    10,  100),
                             ('POLY',    10,  10.0),
                             ('POLY',    10,  99.9),
                             ('POLY',    1,   100),
                             ('POLY',    1,   10.0),
                             ('POLY',    1,   99.9),
                             ('POLY',    0.1, 100),
                             ('POLY',    0.1, 10.0),
                             ('POLY',    0.1, 99.9),
                             ('RBF',     10,  100),
                             ('RBF',     10,  10.0),
                             ('RBF',     10,  99.9),
                             ('RBF',     1,   100),
                             ('RBF',     1,   10.0),
                             ('RBF',     1,   99.9),
                             ('RBF',     0.1, 100),
                             ('RBF',     0.1, 10.0),
                             ('RBF',     0.1, 99.9),
                             ('SIGMOID', 10,  100),
                             ('SIGMOID', 10,  10.0),
                             ('SIGMOID', 10,  99.9),
                             ('SIGMOID', 1,   100),
                             ('SIGMOID', 1,   10.0),
                             ('SIGMOID', 1,   99.9),
                             ('SIGMOID', 0.1, 100),
                             ('SIGMOID', 0.1, 10.0),
                             ('SIGMOID', 0.1, 99.9)]:
#      pass
#    for x in []:

      for (dist_meas, max_i_cnt) in [(mymath.distL1,    2),
                                     (mymath.distL1, 10000),
                                     (mymath.distL2,    2),
                                     (mymath.distL2, 10000),
                                     (mymath.distCanberra, 2),
                                     (mymath.distCanberra, 10000),
                                     (mymath.distCosine, 2),
                                     (mymath.distCosine, 10000)]:

        tailor_class = classification.TAILOR(descr = 'TAILOR',
                                             max_iter_count = max_i_cnt,
                                             dist_measure = dist_meas,
                                             sample = sr,
                                             kernel_type = svmk,
                                             C = svmC)

        assert tailor_class.sample == sr
        assert tailor_class.C == svmC
        assert tailor_class.kernel_type == svmk
        assert tailor_class.max_iter_count == max_i_cnt
        assert tailor_class.dist_measure == dist_meas
        assert tailor_class.train_w_vec_dict == None
        assert tailor_class.train_match_set == None
        assert tailor_class.train_non_match_set == None

        tailor_class.train(self.w_vec_dict, self.m_set, self.nm_set)
        assert tailor_class.sample == sr
        assert tailor_class.C == svmC
        assert tailor_class.kernel_type == svmk
        assert tailor_class.max_iter_count == max_i_cnt
        assert tailor_class.dist_measure == dist_meas
        assert tailor_class.train_w_vec_dict == self.w_vec_dict
        assert tailor_class.train_match_set == self.m_set
        assert tailor_class.train_non_match_set == self.nm_set

        test_res = tailor_class.test(self.w_vec_dict,self.m_set,self.nm_set)
        assert len(test_res) == 4
        assert tailor_class.sample == sr
        assert tailor_class.C == svmC
        assert tailor_class.kernel_type == svmk
        assert tailor_class.max_iter_count == max_i_cnt
        assert tailor_class.dist_measure == dist_meas
        assert tailor_class.train_w_vec_dict == self.w_vec_dict
        assert tailor_class.train_match_set == self.m_set
        assert tailor_class.train_non_match_set == self.nm_set

        test_res = tailor_class.cross_validate(self.w_vec_dict, self.m_set,
                                              self.nm_set, 2)
        assert len(test_res) == 4
        assert tailor_class.sample == sr
        assert tailor_class.C == svmC
        assert tailor_class.kernel_type == svmk
        assert tailor_class.max_iter_count == max_i_cnt
        assert tailor_class.dist_measure == dist_meas
        #assert tailor_class.train_w_vec_dict == self.w_vec_dict
        #assert tailor_class.train_match_set == self.m_set
        #assert tailor_class.train_non_match_set == self.nm_set

        test_res = tailor_class.cross_validate(self.w_vec_dict, self.m_set,
                                               self.nm_set, 3)
        assert len(test_res) == 4
        assert tailor_class.sample == sr
        assert tailor_class.C == svmC
        assert tailor_class.kernel_type == svmk
        assert tailor_class.max_iter_count == max_i_cnt
        assert tailor_class.dist_measure == dist_meas
        #assert tailor_class.train_w_vec_dict == self.w_vec_dict
        #assert tailor_class.train_match_set == self.m_set
        #assert tailor_class.train_non_match_set == self.nm_set

        test_res = tailor_class.cross_validate(self.w_vec_dict, self.m_set,
                                               self.nm_set, 5)
        assert len(test_res) == 4
        assert tailor_class.sample == sr
        assert tailor_class.C == svmC
        assert tailor_class.kernel_type == svmk
        assert tailor_class.max_iter_count == max_i_cnt
        assert tailor_class.dist_measure == dist_meas
        #assert tailor_class.train_w_vec_dict == self.w_vec_dict
        #assert tailor_class.train_match_set == self.m_set
        #assert tailor_class.train_non_match_set == self.nm_set

        class_res = tailor_class.classify(self.test_w_vec_dict)
        assert len(class_res) == 3
        assert isinstance(class_res[0], set) == True
        assert isinstance(class_res[1], set) == True
        assert isinstance(class_res[2], set) == True
        if ((len(class_res[0]) > 0) or (len(class_res[1]) > 0) or \
            (len(class_res[2]) > 0)):  # Only if classification was possible
          assert len(class_res[0]) + len(class_res[1]) + len(class_res[2]) == \
                 len(self.test_w_vec_dict)

        # Train classifier when initialising it - - - - - - - - - - - - - - - -
        #
        tailor_class2 = classification.TAILOR(descr = 'TAILOR',
                                              max_iter_count = max_i_cnt,
                                              dist_measure = dist_meas,
                                              sample = sr,
                                              kernel_type = svmk,
                                              C = svmC,
                                              train_w_vec_dict=self.w_vec_dict,
                                              train_match_set=self.m_set,
                                              train_non_match_set=self.nm_set)
        assert tailor_class2.sample == sr
        assert tailor_class2.C == svmC
        assert tailor_class2.kernel_type == svmk
        assert tailor_class2.max_iter_count == max_i_cnt
        assert tailor_class2.dist_measure == dist_meas
        assert tailor_class2.train_w_vec_dict == self.w_vec_dict
        assert tailor_class2.train_match_set == self.m_set
        assert tailor_class2.train_non_match_set == self.nm_set

        test_res = tailor_class2.test(self.w_vec_dict,self.m_set,self.nm_set)
        assert len(test_res) == 4
        assert tailor_class2.sample == sr
        assert tailor_class2.C == svmC
        assert tailor_class2.kernel_type == svmk
        assert tailor_class2.max_iter_count == max_i_cnt
        assert tailor_class2.dist_measure == dist_meas
        assert tailor_class2.train_w_vec_dict == self.w_vec_dict
        assert tailor_class2.train_match_set == self.m_set
        assert tailor_class2.train_non_match_set == self.nm_set

        test_res = tailor_class2.cross_validate(self.w_vec_dict, self.m_set,
                                                self.nm_set, 2)
        assert len(test_res) == 4
        assert tailor_class2.sample == sr
        assert tailor_class2.C == svmC
        assert tailor_class2.kernel_type == svmk
        assert tailor_class2.max_iter_count == max_i_cnt
        assert tailor_class2.dist_measure == dist_meas
        #assert tailor_class2.train_w_vec_dict == self.w_vec_dict
        #assert tailor_class2.train_match_set == self.m_set
        #assert tailor_class.train_non_match_set == self.nm_set

        test_res = tailor_class2.cross_validate(self.w_vec_dict, self.m_set,
                                                self.nm_set, 3)
        assert len(test_res) == 4
        assert tailor_class2.sample == sr
        assert tailor_class2.C == svmC
        assert tailor_class2.kernel_type == svmk
        assert tailor_class2.max_iter_count == max_i_cnt
        assert tailor_class2.dist_measure == dist_meas
        #assert tailor_class2.train_w_vec_dict == self.w_vec_dict
        #assert tailor_class2.train_match_set == self.m_set
        #assert tailor_class2.train_non_match_set == self.nm_set

        test_res = tailor_class.cross_validate(self.w_vec_dict, self.m_set,
                                               self.nm_set, 5)
        assert len(test_res) == 4
        assert tailor_class2.sample == sr
        assert tailor_class2.C == svmC
        assert tailor_class2.kernel_type == svmk
        assert tailor_class2.max_iter_count == max_i_cnt
        assert tailor_class2.dist_measure == dist_meas
        #assert tailor_class2.train_w_vec_dict == self.w_vec_dict
        #assert tailor_class2.train_match_set == self.m_set
        #assert tailor_class2.train_non_match_set == self.nm_set

        class_res = tailor_class2.classify(self.test_w_vec_dict)
        assert len(class_res) == 3
        assert isinstance(class_res[0], set) == True
        assert isinstance(class_res[1], set) == True
        assert isinstance(class_res[2], set) == True
        if ((len(class_res[0]) > 0) or (len(class_res[1]) > 0) or \
            (len(class_res[2]) > 0)):  # Only if classification was possible
          assert len(class_res[0]) + len(class_res[1]) + len(class_res[2]) == \
                 len(self.test_w_vec_dict)


  def testGetTrueMatchesNonMatches(self):  # - - - - - - - - - - - - - - - - -
    """Test get_true_matches_nonmatches function"""

    # True matches have lastdigit in record identifiers the same and first
    # weight 1.0

    w_vec_dict = {('r11','r21'):[1.0, 0.5, 0.0],
                  ('r12','r21'):[0.6, 0.5, 0.0],
                  ('r22','r32'):[1.0, 0.0, 7.0],
                  ('r13','r22'):[0.0, 0.5, 0.0],
                  ('r22','r30'):[0.9, 0.1, 7.0],
                  ('r10','r60'):[1.0, 1.0, 0.1]}

    def match_funct1(r_id1, r_id2, w_vec):
      if (w_vec[0] == 1.0):
        return True
      else:
        return False

    def match_funct2(r_id1, r_id2, w_vec):
      if (r_id1[-1] == r_id2[-1]):
        return True
      else:
        return False

    m_set, nm_set = classification.get_true_matches_nonmatches(w_vec_dict,
                                                               match_funct1)
    assert isinstance(m_set, set) == True
    assert isinstance(nm_set, set) == True
    assert len(m_set) + len(nm_set) == len(w_vec_dict)
    assert len(m_set.intersection(nm_set)) == 0

    m_set, nm_set = classification.get_true_matches_nonmatches(w_vec_dict,
                                                               match_funct2)
    assert isinstance(m_set, set) == True
    assert isinstance(nm_set, set) == True
    assert len(m_set) + len(nm_set) == len(w_vec_dict)
    assert len(m_set.intersection(nm_set)) == 0

  def testExtractCollapseWeightVectors(self):  # - - - - - - - - - - - - - - -
    """Test extract_collapse_weight_vectors function"""

    w_vec_dict = {('r11','r21'):[1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
                  ('r12','r21'):[0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
                  ('r22','r32'):[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                  ('r13','r22'):[2.0, 1.5, 1.0, 0.5, 1.0, 1.5],
                  ('r22','r30'):[0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
                  ('r10','r60'):[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]}


    mani_list = [(0,1,2,3,4,5)]  # Sum all into one

    res_w_vec_dict = classification.extract_collapse_weight_vectors(mani_list,
                                                                    w_vec_dict)
    assert isinstance(res_w_vec_dict, dict) == True
    assert len(res_w_vec_dict) == len(w_vec_dict)

    org_keys = w_vec_dict.keys()
    new_keys = res_w_vec_dict.keys()
    org_keys.sort()
    new_keys.sort()
    assert org_keys == new_keys

    for (k, v) in res_w_vec_dict.items():
      assert len(v) == 1
      assert sum(w_vec_dict[k]) == v[0]

    vec_w_list = [2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
    res_w_vec_dict = classification.extract_collapse_weight_vectors(mani_list,
                                                                    w_vec_dict,
                                                                    vec_w_list)
    assert isinstance(res_w_vec_dict, dict) == True
    assert len(res_w_vec_dict) == len(w_vec_dict)
    org_keys = w_vec_dict.keys()
    new_keys = res_w_vec_dict.keys()
    org_keys.sort()
    new_keys.sort()
    assert org_keys == new_keys

    for (k, v) in res_w_vec_dict.items():
      assert len(v) == 1
      assert 2.0*sum(w_vec_dict[k]) == v[0]

    mani_list = [(0,1),(2,3),(4,5)]  # Sum all into three

    res_w_vec_dict = classification.extract_collapse_weight_vectors(mani_list,
                                                                    w_vec_dict)
    assert isinstance(res_w_vec_dict, dict) == True
    assert len(res_w_vec_dict) == len(w_vec_dict)

    org_keys = w_vec_dict.keys()
    new_keys = res_w_vec_dict.keys()
    org_keys.sort()
    new_keys.sort()
    assert org_keys == new_keys

    for (k, v) in res_w_vec_dict.items():
      assert len(v) == 3
      assert sum(w_vec_dict[k]) == sum(v)
      assert sum(w_vec_dict[k][0:2]) == v[0]
      assert sum(w_vec_dict[k][2:4]) == v[1]
      assert sum(w_vec_dict[k][4:6]) == v[2]

    vec_w_list = [3.0, 3.0, 4.0, 4.0, 5.0, 5.0]
    res_w_vec_dict = classification.extract_collapse_weight_vectors(mani_list,
                                                                    w_vec_dict,
                                                                    vec_w_list)
    assert isinstance(res_w_vec_dict, dict) == True
    assert len(res_w_vec_dict) == len(w_vec_dict)
    org_keys = w_vec_dict.keys()
    new_keys = res_w_vec_dict.keys()
    org_keys.sort()
    new_keys.sort()
    assert org_keys == new_keys

    for (k, v) in res_w_vec_dict.items():
      assert len(v) == 3
      assert 3.0*sum(w_vec_dict[k][0:2]) == v[0]
      assert 4.0*sum(w_vec_dict[k][2:4]) == v[1]
      assert 5.0*sum(w_vec_dict[k][4:6]) == v[2]

    mani_list = [(0,),(4,),(5,)]  # Only first and last two elements

    res_w_vec_dict = classification.extract_collapse_weight_vectors(mani_list,
                                                                    w_vec_dict)
    assert isinstance(res_w_vec_dict, dict) == True
    assert len(res_w_vec_dict) == len(w_vec_dict)

    org_keys = w_vec_dict.keys()
    new_keys = res_w_vec_dict.keys()
    org_keys.sort()
    new_keys.sort()
    assert org_keys == new_keys

    for (k, v) in res_w_vec_dict.items():
      assert len(v) == 3
      assert w_vec_dict[k][0] == v[0]
      assert w_vec_dict[k][4] == v[1]
      assert w_vec_dict[k][5] == v[2]

    vec_w_list = [3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    res_w_vec_dict = classification.extract_collapse_weight_vectors(mani_list,
                                                                    w_vec_dict,
                                                                    vec_w_list)
    assert isinstance(res_w_vec_dict, dict) == True
    assert len(res_w_vec_dict) == len(w_vec_dict)
    org_keys = w_vec_dict.keys()
    new_keys = res_w_vec_dict.keys()
    org_keys.sort()
    new_keys.sort()
    assert org_keys == new_keys

    for (k, v) in res_w_vec_dict.items():
      assert len(v) == 3
      assert 3.0*w_vec_dict[k][0] == v[0]
      assert 7.0*w_vec_dict[k][4] == v[1]
      assert 8.0*w_vec_dict[k][5] == v[2]

# =============================================================================
# Start tests when called from command line

if (__name__ == "__main__"):
  unittest.main()  # Run all test

# =============================================================================
