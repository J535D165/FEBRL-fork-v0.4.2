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
# The Original Software is: "phonenumTest.py"
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

"""Test module for phonenum.py.
"""

# =============================================================================
# Import necessary modules (Python standard modules first, then Febrl modules)

import sets
import sys
import unittest
sys.path.append('..')

import phonenum

# =============================================================================

class TestCase(unittest.TestCase):

  # Initialise test case  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def setUp(self):
    self.tests = [('++61 2 6125 5690',
                   ['61', 'Australia', '02', '6125-5690', '']),
                  ('0061 02 6125 5690',
                   ['61', 'Australia', '02', '6125-5690', '']),
                  ('0061   02    6125-5690',
                   ['61', 'Australia', '02', '6125-5690', '']),
                  ('41 312 17 84',
                   ['41', 'Switzerland', '', '312 17 84', '']),
                  ('6125 0010',
                   ['61', 'Australia', '', '6125-0010', '']),
                  ('1-800-764-0432',
                   ['1', 'USA/Canada', '800', '764-0432', '']),
                  ('02 6125 0010',
                   ['61', 'Australia', '02', '6125-0010', '']),
                  ('00 1 317-923 4523',
                   ['1', 'USA/Canada', '317', '923-4523', '']),
                  ('1 317-923 4523',
                   ['1', 'USA/Canada', '317', '923-4523', '']),
                  ('00111 41 312 17 84',
                   ['41', 'Switzerland', '', '312 17 84', '']),
                  ('00001 41 312 17 84',
                   ['41', 'Switzerland', '', '312 17 84', '']),
                  ('01 41 312 17 84',
                   ['41', 'Switzerland', '', '312 17 84', '']),
                  ('1-541-754-3010',
                   ['1', 'USA/Canada', '541', '754-3010', '']),
                  ('754-3010',
                   ['1', 'USA/Canada',   '', '754-3010', '']),
                  ('754-3010ext 42',
                   ['1', 'USA/Canada',   '', '754-3010', '42']),
                  ('754-3010x 42',
                   ['1', 'USA/Canada',   '', '754-3010', '42']),
                  ('754-3010 ext 42',
                   ['1', 'USA/Canada',   '', '754-3010', '42']),
                  ('754-3010 ext. 42',
                   ['1', 'USA/Canada',   '', '754-3010', '42']),
                  ('754-3010 x. 42',
                   ['1', 'USA/Canada',   '', '754-3010', '42']),
                  ('754-3010 x42',
                   ['1', 'USA/Canada',   '', '754-3010', '42']),
                  ('(541) 754-3010',
                   ['1', 'USA/Canada', '541', '754-3010', '']),
                  ('+1-541-754-3010',
                   ['1', 'USA/Canada', '541', '754-3010', '']),
                  ('191 541 754 3010',
                   ['', '', '', '915417543010', '']),
                  ('001-541-754-3010',
                   ['1', 'USA/Canada', '541', '754-3010', '']),
                  ('636-48018',
                   ['61', 'Australia', '', '6364-8018', '']),
                  ('(089) / 636-48018',
                   ['1', 'USA/Canada', '896', '364-8018', '']),
                  ('+49-89-636-48018',
                   ['49', 'Germany', '', '89 636 48018', '']),
                  ('19-49-89-636-48018',
                   ['', '', '', '9498963648018', '']),
                  ('+61 (02) 6125 0101',
                   ['61', 'Australia', '02', '6125-0101', '']),
                  ('++61 (02) 6125 0101',
                   ['61', 'Australia', '02', '6125-0101', '']),
                  ('++61 (2) 6125 0101',
                   ['61', 'Australia', '02', '6125-0101', '']),
                  ('11 +61 (2) 6125 0101',
                   ['', '', '', '161261250101', '']),
                  ('0011 ++61 (2) 6125 0101',
                   ['61', 'Australia', '02', '6125-0101', '']),
                  ('0111 ++61 (2) 6125 0101',
                   ['61', 'Australia', '02', '6125-0101', '']),
                  ('0111 61 02 6125 0101',
                   ['61', 'Australia', '02', '6125-0101', '']),
                  ('61 (2) 6125 0101',
                   ['61', 'Australia', '02', '6125-0101', '']),
                 ]

  # Clean up test case  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def tearDown(self):
    pass  # Nothing to clean up

  # ---------------------------------------------------------------------------
  #
  # Start test cases

  def testIt(self):   # - - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'str_to_phonenum' standardisation routine"""

    for (input_number, expected_result) in self.tests:
      this_result = phonenum.str_to_phonenum(input_number)

      assert (isinstance(this_result, list)), \
             '"str_to_phonenum" does not return a list: %s' % \
             (str(this_result))

      assert (len(this_result) == len(expected_result)), \
             '"str_to_phonenum" does not return a list of correct length: %s' \
             % (str(this_result)) + ' (should be: %s)' % (str(expected_result))

      for i in range(len(this_result)):

        assert (isinstance(this_result[i], str)), \
               '"str_to_phonenum" does not return a string at position ' + \
               '%d: "%s"' % (i, str(this_result[i]))

        assert (this_result[i] == expected_result[i]), \
               '"str_to_phonenum" returns a different value in element at ' + \
               'position %d: "%s" (should be: "%s")' % \
               (i, this_result[i], expected_result[i])

# =============================================================================
# Start tests when called from command line

if (__name__ == "__main__"):
  unittest.main()  # Run all test

# =============================================================================
