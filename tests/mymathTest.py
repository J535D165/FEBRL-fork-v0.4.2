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
# The Original Software is: "mymathTest.py"
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

"""Test module for mymath.py.
"""

# =============================================================================
# Import necessary modules (Python standard modules first, then Febrl modules)

import sets
import sys
import unittest
sys.path.append('..')

import mymath

# =============================================================================

class TestCase(unittest.TestCase):

  # Initialise test case  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def setUp(self):

    # Vectors plus their mean and standard deviations
    #
    self.vectors = [([1],             1.0,  0.0),
                    ([1,2.0],         1.5,  0.5),
                    ([10,100],       55.0, 45.0),
                    ([10.0,100.0],   55.0, 45.0),
                    ([1,2,3,4,5,6,7], 4.0, 2.0)]

    # Numbers and their log2 values
    #
    self.log2numbers = [(1,0),(2,1),(4,2),(8,3),(16,4),(32,5),(64,6),(1024,10)]

    # A list of tag sequences and their permuations (number and permutatations)
    #
    self.tag_lists = [(['a','b'],1,[['a','b']]),
                      ([['a','c'],'b'],2,[['a','b'],['c','b']]),
                      (['a',['b','c']],2,[['a','b'],['a','c']]),
                      (['a','b','c','d'],1,[['a','b','c','d']]),
                      ([['a','b'],['c','d']],4,[['a','c'],['b','c'], \
                                                ['a','d'],['b','d']]),
                      ([['a','b'],['c','d','e'],'f'],6,[['a','c','f'], \
                                                        ['b','c','f'], \
                                                        ['a','d','f'], \
                                                        ['b','d','f'], \
                                                        ['a','e','f'], \
                                                        ['b','e','f']]),
                     ]

    # Data for quantiles routine
    #
    self.quant_test_list = [([1,2,3,4,5,6,7,8,9,10],[0.0,0.25,0.5,0.75,1.0]),
                            ([10,9,8,7,6,5,4,3,2,1],[0.0,0.25,0.5,0.75,1.0]),
                            ([10,2,6,5,8,9,1,3,4,7],[0.0,0.25,0.5,0.75,1.0]),
                            ([14,13,12,11,10,9,8,7,6,5,4,3],
                             [0.05,0.25,0.5,0.75,0.95]),
                            ([3,4,5,6,7,8,9,10,11,12,13,14],
                             [0.05,0.25,0.5,0.75,0.95]),
                            ([3,14,5,7,6,12,9,10,13,8,11,4],
                             [0.05,0.25,0.5,0.75,0.95])]
    self.quant_res_list = [[1.0,3.25,5.5,7.75,10.0],
                           [1.0,3.25,5.5,7.75,10.0],
                           [1.0,3.25,5.5,7.75,10.0],
                           [3.55,5.75,8.5,11.25,13.45],
                           [3.55,5.75,8.5,11.25,13.45],
                           [3.55,5.75,8.5,11.25,13.45]]

  # Clean up test case  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def tearDown(self):
    pass  # Nothing to clean up

  # ---------------------------------------------------------------------------
  #
  # Start test cases

  def testMean(self):   # - - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'mean' routine"""

    for i in self.vectors:

      m = mymath.mean(i[0])

      assert (isinstance(m,float)), \
             'Value returned from "mean" is not a float: '+str(m)

      assert m == i[1], \
             'Wrong "mean" with data: '+str(i[0])+' (should be: '+str(i[1])+ \
             '): '+str(m)

  def testStdDev(self):   # - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'stddev' routine"""

    for i in self.vectors:

      s = mymath.stddev(i[0])

      assert (isinstance(s, float)), \
             'Value returned from "stddev" is not a float: '+str(s)

      assert s == i[2], \
             'Wrong "stddev" with data: '+str(i[0])+' (should be: '+ \
             str(i[2])+'): '+str(s)

  def testLog2(self):   # - - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'log2' routine"""

    for n in self.log2numbers:

      l = mymath.log2(n[0])

      assert (isinstance(l, float)), \
             'Value returned from "log2" is not a float: '+str(l)

      assert l == n[1], \
             'Wrong "log2" with value: '+str(n[0])+' (should be: '+ \
             str(n[1])+'): '+str(l)

  def testPermTagSeq(self):   # - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'perm_tag_sequence' routine"""

    for l in self.tag_lists:

      t = mymath.perm_tag_sequence(l[0])

      assert len(t) == l[1], \
             '"perm_tag_sequence" returns wrong number of permutations with '+\
             'list: '+str(l[0])+' (should be: '+str(l[1])+'): '+str(len(t))

      for i in range(len(t)):
        assert t[i] == l[2][i], \
               '"perm_tag_sequence" returns wrong permutation: '+str(t[i])+ \
               ', should be: '+str(l[2][i])

  def testQuantiles(self):  # - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'quantiles' routine"""

    for i in range(len(self.quant_test_list)):
      (in_data, quant_list) = self.quant_test_list[i]
      exp_list = self.quant_res_list[i]

      val_list = mymath.quantiles(in_data, quant_list)

      assert exp_list == val_list, \
             '"quantiles" returns wrong value list: %s (should be: %s)' % \
             (str(val_list), str(exp_list))

  def testDistances(self):  # - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test distances routines"""

    test_vec_pair_list = [([1,2,3],[2,3,4]),
                          ([3,5,2,1],[4,5,3,2]),
                          ([1],[4]),
                          ([0.1,0.2,0.3,0.4,999.99],[0,0,0,1,88.01]),
                         ]

    for (v1,v2) in test_vec_pair_list:

      l1_dist =  mymath.distL1(v1,v2)
      l2_dist =  mymath.distL2(v1,v2)
      li_dist =  mymath.distLInf(v1,v2)
      cbr_dist = mymath.distCanberra(v1,v2)
      cos_dist = mymath.distCosine(v1,v2)

      assert isinstance(l1_dist, float)
      assert isinstance(l2_dist, float)
      assert isinstance(li_dist, float)
      assert isinstance(cbr_dist, float)
      assert isinstance(cos_dist, float)

      assert l1_dist >=  0
      assert l2_dist >=  0
      assert li_dist >=  0
      assert cbr_dist >= 0
      assert cos_dist >= 0

      l1_dist2 =  mymath.distL1(v2,v1)
      l2_dist2 =  mymath.distL2(v2,v1)
      li_dist2 =  mymath.distLInf(v2,v1)
      cbr_dist2 = mymath.distCanberra(v2,v1)
      cos_dist2 = mymath.distCosine(v2,v1)

      assert l1_dist == l1_dist2
      assert l2_dist == l2_dist2
      assert li_dist == li_dist2
      assert cbr_dist == cbr_dist2
      assert cos_dist == cos_dist2

      l1_dist = mymath.distL1(v1,v1)
      l2_dist = mymath.distL2(v2,v2)
      li_dist = mymath.distLInf(v1,v1)
      cbr_dist = mymath.distCanberra(v1,v1)
      cos_dist = mymath.distCosine(v1,v1)

      assert isinstance(l1_dist, float)
      assert isinstance(l2_dist, float)
      assert isinstance(li_dist, float)
      assert isinstance(cbr_dist, float)
      assert isinstance(cos_dist, float)

      assert l1_dist ==  0
      assert l2_dist ==  0
      assert li_dist ==  0
      assert cbr_dist == 0
      assert cos_dist == 0

  def testRandom(self):  # - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test random distributions routine"""

    for i in range(1,100000):

      r1 = mymath.random_linear(i)
      assert (r1 >= 0) and (r1 < i), (i, r1)

      r2 = mymath.random_expo(i)
      assert (r2 >= 0) and (r2 < i), (i, r2)

# =============================================================================
# Start tests when called from command line

if (__name__ == "__main__"):
  unittest.main()  # Run all test

# =============================================================================
