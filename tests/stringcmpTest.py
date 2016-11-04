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
# The Original Software is: "stringcmpTest.py"
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

"""Module stringcmpTest.py - Test module for stringcmp.py.
"""

# =============================================================================
# Import necessary modules (Python standard modules first, then Febrl modules)

import logging
import sys
import unittest
sys.path.append('..')

import stringcmp

log_level = logging.WARNING  # logging.INFO

# =============================================================================

class TestCase(unittest.TestCase):

  # Initialise test case  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def setUp(self):
    self.string_pairs = [['peter',            'peter'            ],
                         ['peter',            ' peter'           ],
                         ['peter',            'xanthalope"'      ],
                         ['christensen',      'christensen'      ],
                         ['prapasi asawakun', 'prapasi asawakun',],
                         ['prapasi asawakun', 'paprasi asawakun',],
                         ['shackleford',      'shackelford'      ],
                         ['dunningham',       'cunnigham'        ],
                         ['nichleson',        'nichulson'        ],
                         ['jones',            'johnson'          ],
                         ['massey',           'massie'           ],
                         ['abroms',           'abrams'           ],
                         ['hardin',           'martinez'         ],
                         ['aa',               'aaa'              ],
                         ['aab',              'aaa'              ],
                         ['aab  ',            'a  a  a'          ],
                         ['aaaaaaaaaa',       'aaaaaaaaab'       ],
                         ['aaaaaaaaaa',       'aaaaaaaa'         ],
                         ['$',                '#'                ],
                         ['1234567890',       '2345678901'       ],
                         ['re$ mkM )"-',      're$ mkM )"-'      ],
                         ['re$ mkM )"-',      're$ mkM )"-'      ],
                         ['re$% mkM %")-',    '" or4%~~!][{ .'   ],
                         ['kim zho',          'zhou kim'         ],
                         ['lim zhao',         'zho  kim'         ],
                         ['kim lim zhao',     'lim zhau kim'     ]]

  # Clean up test case  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def tearDown(self):
    pass  # Nothing to clean up

  # ---------------------------------------------------------------------------
  #
  # Start test cases

  def testJaro(self):   # - - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'Jaro' approximate string comparator"""

    for pair in self.string_pairs:

      approx_str_value = stringcmp.jaro(pair[0],pair[1])

      assert (isinstance(approx_str_value,float)), \
             '"Jaro" does not return a floating point number for: '+str(pair)

      assert (approx_str_value >= 0.0), \
             '"Jaro" returns a negative number for: '+str(pair)

      assert (approx_str_value <= 1.0), \
             '"Jaro" returns a number larger than 1.0 for: '+str(pair)

      approx_str_value_1 = stringcmp.jaro(pair[0],pair[1])
      approx_str_value_2 = stringcmp.jaro(pair[1],pair[0])

      assert (approx_str_value_1 == approx_str_value_2), \
             '"Jaro" returns different values for pair and swapped pair: '+ \
             str(pair)+': '+str(approx_str_value_1)+', '+ \
             str(approx_str_value_2)

      # Check for value 1.0 if the strings are the same
      #
      if (pair[0] == pair[1]):

        assert (approx_str_value == 1.0), \
               '"Jaro" does not return 1.0 if strings are equal: '+str(pair)


  def testWinkler(self):  # - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'Winkler' approximate string comparator"""

    for pair in self.string_pairs:

      approx_str_value = stringcmp.winkler(pair[0],pair[1])

      assert (isinstance(approx_str_value,float)), \
             '"Winkler" does not return a floating point number for:'+ \
             str(pair)

      assert (approx_str_value >= 0.0), \
             '"Winkler" returns a negative number for:'+str(pair)

      assert (approx_str_value <= 1.0), \
             '"Winkler" returns a number larger than 1.0 for:'+str(pair)

      approx_str_value_1 = stringcmp.winkler(pair[0],pair[1])
      approx_str_value_2 = stringcmp.winkler(pair[1],pair[0])

      assert (approx_str_value_1 == approx_str_value_2), \
             '"Winkler" returns different values for pair and swapped ' + \
             'pair: '+str(pair)+': '+str(approx_str_value_1)+', '+ \
             str(approx_str_value_2)

      # Check for value 1.0 if the strings are the same
      #
      if (pair[0] == pair[1]):

        assert (approx_str_value == 1.0), \
               '"Winkler" does not return 1.0 if strings are equal: '+str(pair)

      # Winkler should always return a value equal to or larger than Jaro
      #
      approx_str_value_winkler = stringcmp.winkler(pair[0],pair[1])
      approx_str_value_jaro =    stringcmp.jaro(pair[0],pair[1])

      assert (approx_str_value_winkler >= approx_str_value_jaro), \
             '"Winkler" value smaller than "Jaro" value for:'+str(pair)


  def testBigram(self):   # - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'Bigram' approximate string comparator"""

    for pair in self.string_pairs:

      approx_str_value = stringcmp.bigram(pair[0],pair[1])

      assert (isinstance(approx_str_value,float)), \
             '"Bigram" does not return a floating point number for: '+str(pair)

      assert (approx_str_value >= 0.0), \
             '"Bigram" returns a negative number for: '+str(pair)

      assert (approx_str_value <= 1.0), \
             '"Bigram" returns a number larger than 1.0 for: '+str(pair)

      approx_str_value_1 = stringcmp.bigram(pair[0],pair[1])
      approx_str_value_2 = stringcmp.bigram(pair[1],pair[0])

      assert (approx_str_value_1 == approx_str_value_2), \
             '"Bigram" returns different values for pair and swapped ' + \
             'pair: '+str(pair)+': '+str(approx_str_value_1)+', '+ \
             str(approx_str_value_2)

      # Check for value 1.0 if the strings are the same
      #
      if (pair[0] == pair[1]):

        assert (approx_str_value == 1.0), \
               '"Bigram" does not return 1.0 if strings are equal: '+str(pair)


  def testBagDist(self):   # - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'BagDist' approximate string comparator"""

    for pair in self.string_pairs:

      approx_str_value =   stringcmp.bagdist(pair[0],pair[1])

      assert (isinstance(approx_str_value,float)), \
             '"BagDist" does not return a floating point number for: '+ \
             str(pair)

      assert (approx_str_value >= 0.0), \
             '"BagDist" returns a negative number for: '+str(pair)

      assert (approx_str_value <= 1.0), \
             '"BagDist" returns a number larger than 1.0 for: '+str(pair)

      approx_str_value_1 = stringcmp.bagdist(pair[0],pair[1])
      approx_str_value_2 = stringcmp.bagdist(pair[1],pair[0])

      assert (approx_str_value_1 == approx_str_value_2), \
             '"BagDist" returns different values for pair and swapped ' + \
             'pair: '+str(pair)+': '+str(approx_str_value_1)+', '+ \
             str(approx_str_value_2)

      # Check for value 1.0 if the strings are the same
      #
      if (pair[0] == pair[1]):

        assert (approx_str_value == 1.0), \
               '"BagDist" does not return 1.0 if strings are equal: '+ \
               str(pair)

      # Check if bad distance is always larger than edit distance

      editdist_str_value = stringcmp.editdist(pair[0],pair[1])
      assert (approx_str_value >= editdist_str_value), \
             '"BagDist" value is smaller than "EditDist" value for: '+ \
             str(pair)

  def testEditDist(self):   # - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'EditDist' approximate string comparator"""

    for pair in self.string_pairs:

      approx_str_value = stringcmp.editdist(pair[0],pair[1])

      assert (isinstance(approx_str_value,float)), \
             '"EditDist" does not return a floating point number for: '+ \
             str(pair)

      assert (approx_str_value >= 0.0), \
             '"EditDist" returns a negative number for: '+str(pair)

      assert (approx_str_value <= 1.0), \
             '"EditDist" returns a number larger than 1.0 for: '+str(pair)

      approx_str_value_1 = stringcmp.editdist(pair[0],pair[1])
      approx_str_value_2 = stringcmp.editdist(pair[1],pair[0])

      assert (approx_str_value_1 == approx_str_value_2), \
             '"EditDist" returns different values for pair and swapped ' + \
             'pair: '+str(pair)+': '+str(approx_str_value_1)+', '+ \
             str(approx_str_value_2)

      # Check for value 1.0 if the strings are the same
      #
      if (pair[0] == pair[1]):

        assert (approx_str_value == 1.0), \
               '"EditDist" does not return 1.0 if strings are equal: '+ \
               str(pair)


  def testSeqMatch(self):   # - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'SeqMatch' approximate string comparator"""

    for pair in self.string_pairs:

      approx_str_value = stringcmp.seqmatch(pair[0],pair[1])

      assert (isinstance(approx_str_value,float)), \
             '"SeqMatch" does not return a floating point number for: '+ \
             str(pair)

      assert (approx_str_value >= 0.0), \
             '"SeqMatch" returns a negative number for: '+str(pair)

      assert (approx_str_value <= 1.0), \
             '"SeqMatch" returns a number larger than 1.0 for: '+str(pair)

      approx_str_value_1 = stringcmp.seqmatch(pair[0],pair[1])
      approx_str_value_2 = stringcmp.seqmatch(pair[1],pair[0])

      assert (approx_str_value_1 == approx_str_value_2), \
             '"SeqMatch" returns different values for pair and swapped ' + \
             'pair: '+str(pair)+': '+str(approx_str_value_1)+', '+ \
             str(approx_str_value_2)

      # Check for value 1.0 if the strings are the same
      #
      if (pair[0] == pair[1]):

        assert (approx_str_value == 1.0), \
               '"SeqMatch" does not return 1.0 if strings are equal: '+ \
               str(pair)

  def testCompression(self):   # - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'Compression' approximate string comparator"""

    for pair in self.string_pairs:

      approx_str_value = stringcmp.compression(pair[0],pair[1])

      assert (isinstance(approx_str_value,float)), \
             '"Compression" does not return a floating point number for: '+ \
             str(pair)

      assert (approx_str_value >= 0.0), \
             '"Compression" returns a negative number for: '+str(pair)

      assert (approx_str_value <= 1.0), \
             '"Compression" returns a number larger than 1.0 for: '+str(pair)

      approx_str_value_1 = stringcmp.compression(pair[0],pair[1])
      approx_str_value_2 = stringcmp.compression(pair[1],pair[0])

      assert (approx_str_value_1 == approx_str_value_2), \
             '"Compression" returns different values for pair and swapped ' + \
             'pair: '+str(pair)+': '+str(approx_str_value_1)+', '+ \
             str(approx_str_value_2)

      # Check for value 1.0 if the strings are the same
      #
      if (pair[0] == pair[1]):

        assert (approx_str_value == 1.0), \
               '"Compression" does not return 1.0 if strings are equal: '+ \
               str(pair)


  def testLCS(self):   # - - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'LCS' approximate string comparator"""

    for pair in self.string_pairs:

      approx_str_value = stringcmp.lcs(pair[0],pair[1])

      assert (isinstance(approx_str_value,float)), \
             '"LCS" does not return a floating point number for: '+ \
             str(pair)

      assert (approx_str_value >= 0.0), \
             '"LCS" returns a negative number for: '+str(pair)

      assert (approx_str_value <= 1.0), \
             '"LCS" returns a number larger than 1.0 for: '+str(pair)

      approx_str_value_1 = stringcmp.lcs(pair[0],pair[1])
      approx_str_value_2 = stringcmp.lcs(pair[1],pair[0])

      assert (approx_str_value_1 == approx_str_value_2), \
             '"LCS" returns different values for pair and swapped ' + \
             'pair: '+str(pair)+': '+str(approx_str_value_1)+', '+ \
             str(approx_str_value_2)

      # Check for value 1.0 if the strings are the same
      #
      if (pair[0] == pair[1]):

        assert (approx_str_value == 1.0), \
               '"LCS" does not return 1.0 if strings are equal: '+str(pair)


  def testOntoLCS(self):   # - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'OntoLCS' approximate string comparator"""

    for pair in self.string_pairs:

      approx_str_value = stringcmp.ontolcs(pair[0],pair[1])

      assert (isinstance(approx_str_value,float)), \
             '"OntoLCS" does not return a floating point number for: '+ \
             str(pair)

      assert (approx_str_value >= 0.0), \
             '"OntoLCS" returns a negative number for: '+str(pair)

      assert (approx_str_value <= 1.0), \
             '"OntoLCS" returns a number larger than 1.0 for: '+str(pair)

      approx_str_value_1 = stringcmp.ontolcs(pair[0],pair[1])
      approx_str_value_2 = stringcmp.ontolcs(pair[1],pair[0])

      assert (approx_str_value_1 == approx_str_value_2), \
             '"OntoLCS" returns different values for pair and swapped ' + \
             'pair: '+str(pair)+': '+str(approx_str_value_1)+', '+ \
             str(approx_str_value_2)

      # Check for value 1.0 if the strings are the same
      #
      if (pair[0] == pair[1]):

        assert (approx_str_value == 1.0), \
               '"OntoLCS" does not return 1.0 if strings are equal: '+str(pair)


  def testPermWinkler(self):  # - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'PermWinkler' approximate string comparator"""

    for pair in self.string_pairs:

      approx_str_value = stringcmp.permwinkler(pair[0],pair[1])

      assert (isinstance(approx_str_value,float)), \
             '"PermWinkler" does not return a floating point number for:'+ \
             str(pair)

      assert (approx_str_value >= 0.0), \
             '"PermWinkler" returns a negative number for:'+str(pair)

      assert (approx_str_value <= 1.0), \
             '"PermWinkler" returns a number larger than 1.0 for:'+str(pair)

      approx_str_value_1 = stringcmp.permwinkler(pair[0],pair[1])
      approx_str_value_2 = stringcmp.permwinkler(pair[1],pair[0])

      assert (approx_str_value_1 == approx_str_value_2), \
             '"PermWinkler" returns different values for pair and swapped ' + \
             'pair: '+str(pair)+': '+str(approx_str_value_1)+', '+ \
             str(approx_str_value_2)

      # Check for value 1.0 if the strings are the same
      #
      if (pair[0] == pair[1]):

        assert (approx_str_value == 1.0), \
               '"PermWinkler" does not return 1.0 if strings are equal: '+ \
               str(pair)

      # PermWinkler should always return a value equal to or larger than
      # Winkler
      #
      approx_str_value_permwinkler = stringcmp.permwinkler(pair[0],pair[1])
      approx_str_value_winkler =     stringcmp.winkler(pair[0],pair[1])

      assert (approx_str_value_permwinkler >= approx_str_value_winkler), \
             '"PermWinkler" value smaller than "Winkler" value for:'+str(pair)

# =============================================================================
# Start tests when called from command line

if (__name__ == "__main__"):

  # Intialise a logger, set level to info
  #
  my_logger = logging.getLogger()  # New logger at root level
  my_logger.setLevel(log_level)

  unittest.main()  # Run all test

# =============================================================================
