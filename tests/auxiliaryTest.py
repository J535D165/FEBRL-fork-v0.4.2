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
# The Original Software is: "auxiliaryTest.py"
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

"""Test module for auxiliary.py.
"""

# =============================================================================
# Import necessary modules (Python standard modules first, then Febrl modules)

import sets
import sys
import unittest
sys.path.append('..')

import auxiliary

# =============================================================================

class TestCase(unittest.TestCase):

  # Initialise test case  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def setUp(self):
    pass  # Nothing to initialise

  # Clean up test case  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def tearDown(self):
    pass  # Nothing to clean up

  # ---------------------------------------------------------------------------
  # Start test cases

  def testIsNotNone(self):  # - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'check_is_not_none' function."""

    assert (auxiliary.check_is_not_none('TestArgument','hello') == None)
    assert (auxiliary.check_is_not_none('TestArgument',1) ==       None)
    assert (auxiliary.check_is_not_none('TestArgument',0) ==       None)
    assert (auxiliary.check_is_not_none('TestArgument',-111) ==    None)
    assert (auxiliary.check_is_not_none('TestArgument',0.555) ==   None)
    assert (auxiliary.check_is_not_none('TestArgument',{}) ==      None)
    assert (auxiliary.check_is_not_none('TestArgument',[]) ==      None)

  def testIsString(self):  # - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'check_is_string' function."""

    assert (auxiliary.check_is_string('TestArgument','hello') == None)
    assert (auxiliary.check_is_string('TestArgument','')      == None)
    assert (auxiliary.check_is_string('TestArgument',"123")   == None)
    assert (auxiliary.check_is_string('TestArgument','-1.23') == None)
    assert (auxiliary.check_is_string('TestArgument',"HELlo") == None)
    assert (auxiliary.check_is_string('TestArgument',"'!?!'") == None)
    assert (auxiliary.check_is_string('TestArgument',"[..]")  == None)

  def testIsNumber(self):  # - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'check_is_number' function."""

    assert (auxiliary.check_is_number('TestArgument',0)       == None)
    assert (auxiliary.check_is_number('TestArgument',-0)      == None)
    assert (auxiliary.check_is_number('TestArgument',1.23)    == None)
    assert (auxiliary.check_is_number('TestArgument',-24.41)  == None)
    assert (auxiliary.check_is_number('TestArgument',1289837) == None)
    assert (auxiliary.check_is_number('TestArgument',-973293) == None)

  def testIsPositive(self):  # - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'check_is_positive' function."""

    assert (auxiliary.check_is_positive('TestArgument',1.)        == None)
    assert (auxiliary.check_is_positive('TestArgument',0.001)     == None)
    assert (auxiliary.check_is_positive('TestArgument',0.0000001) == None)
    assert (auxiliary.check_is_positive('TestArgument',1)         == None)
    assert (auxiliary.check_is_positive('TestArgument',1.474)     == None)
    assert (auxiliary.check_is_positive('TestArgument',1236967)   == None)
    assert (auxiliary.check_is_positive('TestArgument',17676.474) == None)

  def testIsNotNegative(self):  # - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'check_is_not_negative' function."""

    assert (auxiliary.check_is_not_negative('TestArgument',0)         == None)
    assert (auxiliary.check_is_not_negative('TestArgument',-0)        == None)
    assert (auxiliary.check_is_not_negative('TestArgument',1.)        == None)
    assert (auxiliary.check_is_not_negative('TestArgument',0.00)      == None)
    assert (auxiliary.check_is_not_negative('TestArgument',-0.000)    == None)
    assert (auxiliary.check_is_not_negative('TestArgument',1)         == None)
    assert (auxiliary.check_is_not_negative('TestArgument',1.474)     == None)
    assert (auxiliary.check_is_not_negative('TestArgument',1236967)   == None)
    assert (auxiliary.check_is_not_negative('TestArgument',17676.474) == None)

  def testIsNormalised(self):  # - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'check_is_normalised' function."""

    assert (auxiliary.check_is_normalised('TestArgument',0)       == None)
    assert (auxiliary.check_is_normalised('TestArgument',-0)      == None)
    assert (auxiliary.check_is_normalised('TestArgument',0.00)    == None)
    assert (auxiliary.check_is_normalised('TestArgument',-0.0)    == None)
    assert (auxiliary.check_is_normalised('TestArgument',1)       == None)
    assert (auxiliary.check_is_normalised('TestArgument',1.0)     == None)
    assert (auxiliary.check_is_normalised('TestArgument',0.0001)  == None)
    assert (auxiliary.check_is_normalised('TestArgument',1.00000) == None)
    assert (auxiliary.check_is_normalised('TestArgument',1)       == None)
    assert (auxiliary.check_is_normalised('TestArgument',0.5)     == None)
    assert (auxiliary.check_is_normalised('TestArgument',1)       == None)
    assert (auxiliary.check_is_normalised('TestArgument',0.99999) == None)

  def testIsPercentage(self):  # - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'check_is_percentage' function."""

    assert (auxiliary.check_is_percentage('TestArgument',0)        == None)
    assert (auxiliary.check_is_percentage('TestArgument',-0)       == None)
    assert (auxiliary.check_is_percentage('TestArgument',0.00)     == None)
    assert (auxiliary.check_is_percentage('TestArgument',-0.0)     == None)
    assert (auxiliary.check_is_percentage('TestArgument',1)        == None)
    assert (auxiliary.check_is_percentage('TestArgument',1.0)      == None)
    assert (auxiliary.check_is_percentage('TestArgument',0.0001)   == None)
    assert (auxiliary.check_is_percentage('TestArgument',99.00000) == None)
    assert (auxiliary.check_is_percentage('TestArgument',100)      == None)
    assert (auxiliary.check_is_percentage('TestArgument',0.5)      == None)
    assert (auxiliary.check_is_percentage('TestArgument',50)       == None)
    assert (auxiliary.check_is_percentage('TestArgument',50.0001)  == None)
    assert (auxiliary.check_is_percentage('TestArgument',100.0)    == None)
    assert (auxiliary.check_is_percentage('TestArgument',0.99999)  == None)


  def testIsInteger(self):  # - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'check_is_integer' function."""

    assert (auxiliary.check_is_integer('TestArgument',0)      == None)
    assert (auxiliary.check_is_integer('TestArgument',1234)   == None)
    assert (auxiliary.check_is_integer('TestArgument',-96234) == None)

  def testIsFloat(self):  # - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'check_is_float' function."""

    assert (auxiliary.check_is_float('TestArgument',0.0)      == None)
    assert (auxiliary.check_is_float('TestArgument',-0.0)     == None)
    assert (auxiliary.check_is_float('TestArgument',0.123)    == None)
    assert (auxiliary.check_is_float('TestArgument',-65.9203) == None)
    assert (auxiliary.check_is_float('TestArgument',42.123)   == None)

  def testIsDictionary(self):  # - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'check_is_dictionary' function."""

    assert (auxiliary.check_is_dictionary('TestArgument',{})         == None)
    assert (auxiliary.check_is_dictionary('TestArgument', {1:2,6:0}) == None)
    assert (auxiliary.check_is_dictionary('TestArgument', \
            {'a':4,'t':1,(1,4,6):'tr'}) == None)

  def testIsList(self):  # - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'check_is_list' function."""

    assert (auxiliary.check_is_list('TestArgument',[])              == None)
    assert (auxiliary.check_is_list('TestArgument',[1,3,5])         == None)
    assert (auxiliary.check_is_list('TestArgument',['a','56',1,{}]) == None)

  def testIsSet(self):  # - - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'check_is_set' function."""

    assert (auxiliary.check_is_set('TestArgument',sets.Set())          == None)
    assert (auxiliary.check_is_set('TestArgument',sets.Set([1,2,3]))   == None)
    assert (auxiliary.check_is_set('TestArgument',sets.Set(['a','a'])) == None)

    assert (auxiliary.check_is_set('TestArgument',set())          == None)
    assert (auxiliary.check_is_set('TestArgument',set([1,2,3]))   == None)
    assert (auxiliary.check_is_set('TestArgument',set(['a','a'])) == None)

  def testIsTuple(self):  # - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'check_is_tuple' function."""

    assert (auxiliary.check_is_tuple('TestArgument',())        == None)
    assert (auxiliary.check_is_tuple('TestArgument',(1,))      == None)
    assert (auxiliary.check_is_tuple('TestArgument',('a',))    == None)
    assert (auxiliary.check_is_tuple('TestArgument',('a','b')) == None)
    assert (auxiliary.check_is_tuple('TestArgument',(42,'b'))  == None)

  def testIsFlag(self):  # - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'check_is_flag' function."""

    assert (auxiliary.check_is_flag('TestArgument', True)  == None)
    assert (auxiliary.check_is_flag('TestArgument', False) == None)
    assert (auxiliary.check_is_flag('TestArgument', 0)     == None)
    assert (auxiliary.check_is_flag('TestArgument', 1)     == None)

  def testIsFunction(self):  # - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'check_is_function_or_method' function."""

    def f1(x):
      print x

    def f2():
      print 'hello'

    assert (auxiliary.check_is_function_or_method('TestArgument', f1)  == None)
    assert (auxiliary.check_is_function_or_method('TestArgument', f2)  == None)
    assert (auxiliary.check_is_function_or_method('TestArgument',
            self.setUp)  == None)

  def testIsTimeString(self):  # - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'time_string' function."""

    assert (isinstance(auxiliary.time_string(0), str) == True)
    assert (isinstance(auxiliary.time_string(0.0), str) == True)
    assert (isinstance(auxiliary.time_string(1), str) == True)
    assert (isinstance(auxiliary.time_string(1.0), str) == True)
    assert (isinstance(auxiliary.time_string(123456), str) == True)
    assert (isinstance(auxiliary.time_string(123456.65), str) == True)

  def testIsStrVector(self):  # - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'str_vector' function."""

    assert (isinstance(auxiliary.str_vector([]), str) == True)
    assert (isinstance(auxiliary.str_vector([1], 2), str) == True)
    assert (isinstance(auxiliary.str_vector([1], 5), str) == True)
    assert (isinstance(auxiliary.str_vector([1], 2, False), str) == True)
    assert (isinstance(auxiliary.str_vector([1], 5, False), str) == True)
    assert (isinstance(auxiliary.str_vector([1,2], 2), str) == True)
    assert (isinstance(auxiliary.str_vector([1,2], 5), str) == True)
    assert (isinstance(auxiliary.str_vector([1,2], 2, False), str) == True)
    assert (isinstance(auxiliary.str_vector([1,2], 5, False), str) == True)
    assert (isinstance(auxiliary.str_vector([-1,4.2], 2), str) == True)
    assert (isinstance(auxiliary.str_vector([-1,4.2], 5), str) == True)
    assert (isinstance(auxiliary.str_vector([-1,4.2], 2, False), str) == True)
    assert (isinstance(auxiliary.str_vector([-1,4.2], 5, False), str) == True)
    assert (isinstance(auxiliary.str_vector([1,2,123455673], 2), str) == True)
    assert (isinstance(auxiliary.str_vector([1,2,657,964543], 5), str) == True)
    assert (isinstance(auxiliary.str_vector([1,2,677], 2, False), str) == True)
    assert (isinstance(auxiliary.str_vector([1,2,8,2], 5, False), str) == True)

  def testGetMemoryUsage(self):  # - - - - - - - - - - - - - - - - - - - - - -
    """Test 'get_memory_usage' function."""

    x = auxiliary.get_memory_usage()
    assert (x == None) or (isinstance(x,str) == True)
    x = auxiliary.get_memory_usage()
    assert (x == None) or (isinstance(x,str) == True)


# =============================================================================
# Start tests when called from command line

if (__name__ == "__main__"):
  unittest.main()  # Run all test

# =============================================================================
