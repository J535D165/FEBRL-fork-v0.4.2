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
# The Original Software is: "encodeTest.py"
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

"""Module encodeTest.py - Test module for encode.py.
"""

# =============================================================================
# Import necessary modules (Python standard modules first, then Febrl modules)

import logging
import sys
import unittest
sys.path.append('..')

import encode

log_level = logging.WARNING  # logging.INFO

# =============================================================================

class TestCase(unittest.TestCase):

  # Initialise test case  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def setUp(self):
    self.strings = ['peter','christen','ole','nielsen','markus','hegland',
                    'stephen','steve','roberts','tim','churches','xiong',
                    'ng','miller','millar','foccachio','van de hooch',
                    'xiao ching','asawakun','prapasri','von der felde','vest',
                    'west','oioi','ohio','oihcca', 'nielsen', 'kim', 'lim',
                    'computer','record','linkage','probabilistic',
                    'aa','aaaa aaa','  x   ','x']

  # Clean up test case  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def tearDown(self):
    pass  # Nothing to clean up

  # ---------------------------------------------------------------------------
  #
  # Start test cases

  def testSoundex(self):  # - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'Soundex' string encoding"""

    c = encode.soundex('')  # Test with empty string
    assert (c == '0000'), \
           '"Soundex" of empty string is not "0000"'

    for s in self.strings:

      c = encode.soundex(s)

      assert (isinstance(c,str)), \
             '"Soundex" of string "'+s+'"does not return a string: '+ \
             str(type(c))

      assert (s[0] == c[0]), \
             'First character in "Soundex" code for string "'+s+ \
             '" differs from original string: '+str(c)

      assert (len(c) == 4), \
             'Length of "Soundex" code for string "'+s+'" is not four '+ \
             'characters: '+str(c)+' with length: '+str(len(c))

      if (len(c) > 1):
        assert (c[1:].isdigit() == 1), \
               'Characters after first in "Soundex" '+ \
               'code for string "'+s+'" are not digits: '+str(c)

      c = encode.soundex(s,maxlen=1)

      assert (isinstance(c,str)), \
             '"Soundex" of string "'+s+'"does not return a string: '+ \
             str(type(c))

      assert (s[0] == c[0]), \
             'First character in "Soundex" code for string "'+ \
             s+'" differs from original string: '+str(c)

      assert (len(c) == 1), \
             'Length of "Soundex" code for string "'+s+'" is '+ \
             'not one character: '+str(c)+' with length: '+str(len(c))

      c = encode.soundex(s,maxlen=6)

      assert (isinstance(c,str)), \
             '"Soundex" of string "'+s+'"does not return'+ \
             ' a string: '+str(type(c))

      assert (s[0] == c[0]), \
             'First character in "Soundex" code for string "'+ \
             s+'" differs from original string: '+str(c)

      assert (len(c) <= 6), \
             '"Soundex" code for string "'+s+'" is longer than'+ \
             ' six characters: '+str(c)+' with length: '+str(len(c))

      if (len(c) > 1):
        assert (c[1:].isdigit() == 1), \
               'Characters after first in "Soundex" '+ \
               'code for string "'+s+'" are not digits: '+str(c)


  def testModSoundex(self):   # - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'ModSoundex' string encoding"""

    c = encode.mod_soundex('')  # Test with empty string
    assert (c == '0000'), '"ModSoundex" of empty string is not "0000"'

    for s in self.strings:

      c = encode.mod_soundex(s)

      assert (isinstance(c,str)), '"ModSoundex" of string "'+s+'"does not '+ \
             'return a string: '+str(type(c))

      assert (s[0] == c[0]), 'First character in "ModSoundex" code for '+ \
             'string "'+s+'" differs from original string: '+str(c)

      assert (len(c) == 4), 'Length of "ModSoundex" code for string "'+s+'" '+\
             'is not four characters: '+str(c)+' with length: '+str(len(c))

      if (len(c) > 1):
        assert (c[1:].isdigit() == 1), \
               'Characters after first in "ModSoundex'+ \
               '" code for string "'+s+'" are not digits: '+str(c)

      c = encode.mod_soundex(s,maxlen=1)

      assert (isinstance(c,str)), '"ModSoundex" of string "'+s+'"does not '+ \
             'return a string: '+str(type(c))

      assert (s[0] == c[0]), 'First character in "ModSoundex" code for '+ \
             'string "'+s+'" differs from original string: '+str(c)

      assert (len(c) == 1), 'Length of "ModSoundex" code for string "'+s+'" '+\
             'is not one character: '+str(c)+' with length: '+str(len(c))

      c = encode.mod_soundex(s,maxlen=6)

      assert (isinstance(c,str)), '"ModSoundex" of string "'+s+'"does not '+ \
             'return a string: '+str(type(c))

      assert (s[0] == c[0]), 'First character in "ModSoundex" code for '+ \
             'string "'+s+'" differs from original string: '+str(c)

      assert (len(c) <= 6), '"ModSoundex" code for string "'+s+'" is longer '+\
             'than six characters: '+str(c)+' with length: '+str(len(c))

      if (len(c) > 1):
        assert (c[1:].isdigit() == 1), \
               'Characters after first in "ModSoundex'+ \
               '" code for string "'+s+'" are not digits: '+str(c)


  def testPhonex(self):  # - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'Phonex' string encoding"""

    c = encode.phonex('')  # Test with empty string
    assert (c == '0000'), '"Phonex" of empty string is not "0000"'

    for s in self.strings:

      c = encode.phonex(s)

      assert (isinstance(c,str)), '"Phonex" of string "'+s+'"does not return'+\
             ' a string: '+str(type(c))

      assert (len(c) == 4), 'Length of "Phonex" code for string "'+s+'" is '+ \
             'not four characters: '+str(c)+' with length: '+str(len(c))

      if (len(c) > 1):
        assert (c[1:].isdigit() == 1), 'Characters after first in "Phonex" '+ \
               'code for string "'+s+'" are not digits: '+str(c)

      c = encode.phonex(s,maxlen=1)

      assert (isinstance(c,str)), '"Phonex" of string "'+s+'"does not return'+\
             ' a string: '+str(type(c))

      assert (len(c) == 1), 'Length of "Phonex" code for string "'+s+'" is '+ \
             'not one character: '+str(c)+' with length: '+str(len(c))

      c = encode.phonex(s,maxlen=6)

      assert (isinstance(c,str)), '"Phonex" of string "'+s+'"does not return'+\
             ' a string: '+str(type(c))

      assert (len(c) <= 6), '"Phonex" code for string "'+s+'" is longer than'+\
             ' six characters: '+str(c)+' with length: '+str(len(c))

      if (len(c) > 1):
        assert (c[1:].isdigit() == 1), 'Characters after first in "Phonex" '+ \
               'code for string "'+s+'" are not digits: '+str(c)


  def testPhonix(self):  # - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'Phonix' string encoding"""

    c = encode.phonix('')  # Test with empty string
    assert (c == 'a000'), '"Phonix" of empty string is not "a000"'

    c = encode.phonix('', maxlen=3)  # Test with empty string
    assert (c == 'a00'), '"Phonix" of empty string with maxlen=3 is not "a00"'

    for s in self.strings:

      c = encode.phonix(s)

      assert (isinstance(c,str)), '"Phonix" of string "'+s+'"does not return'+\
             ' a string: '+str(type(c))

      assert (len(c) == 4), 'Length of "Phonix" code for string "'+s+'" is '+ \
             'not four characters: '+str(c)+' with length: '+str(len(c))

      if (len(c) > 1):
        assert (c[1:].isdigit() == 1), 'Characters after first in "Phonix" '+ \
               'code for string "'+s+'" are not digits: '+str(c)

      c = encode.phonex(s,maxlen=1)

      assert (isinstance(c,str)), '"Phonix" of string "'+s+'"does not return'+\
             ' a string: '+str(type(c))

      assert (len(c) == 1), 'Length of "Phonix" code for string "'+s+'" is '+ \
             'not one character: '+str(c)+' with length: '+str(len(c))

      c = encode.phonix(s,maxlen=6)

      assert (isinstance(c,str)), '"Phonix" of string "'+s+'"does not return'+\
             ' a string: '+str(type(c))

      assert (len(c) <= 6), '"Phonix" code for string "'+s+'" is longer than'+\
             ' six characters: '+str(c)+' with length: '+str(len(c))

      if (len(c) > 1):
        assert (c[1:].isdigit() == 1), 'Characters after first in "Phonix" '+ \
               'code for string "'+s+'" are not digits: '+str(c)

  def testNYSIIS(self):  # - - - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'NYSIIS' string encoding"""

    c = encode.nysiis('')  # Test with empty string
    assert (c == ''), '"NYSIIS" of empty string is not ""'

    for s in self.strings:

      c = encode.nysiis(s)

      assert (isinstance(c,str)), '"NYSIIS" of string "'+s+'"does not return'+\
             ' a string: '+str(type(c))

      assert (len(c) <= 4), 'Length of "NYSIIS" code for string "'+s+'" is '+ \
             'more than four characters: '+str(c)+' with length: '+str(len(c))

      if (len(c) > 1):
        assert (c[1:].isalpha() == 1), 'Characters after first in "NYSIIS" '+ \
               'code for string "'+s+'" are not letters: '+str(c)

      c = encode.nysiis(s,maxlen=1)

      assert (isinstance(c,str)), '"NYSIIS" of string "'+s+'"does not return'+\
             ' a string: '+str(type(c))

      assert (len(c) == 1), 'Length of "NYSIIS" code for string "'+s+'" is '+ \
             'not one character: '+str(c)+' with length: '+str(len(c))

      c = encode.nysiis(s,maxlen=6)

      assert (isinstance(c,str)), '"NYSIIS" of string "'+s+'"does not return'+\
             ' a string: '+str(type(c))

      assert (len(c) <= 6), '"NYSIIS" code for string "'+s+'" is longer than'+\
             ' six characters: '+str(c)+' with length: '+str(len(c))

      if (len(c) > 1):
        assert (c[1:].isalpha() == 1), 'Characters after first in "NYSIIS" '+ \
               'code for string "'+s+'" are not letters: '+str(c)


  def testDoubleMetaphone(self):  # - - - - - - - - - - - - - - - - - - - - - -
    """Test 'DoubleMetaphone' string encoding"""

    c = encode.dmetaphone('')  # Test with empty string
    assert (c == ''), '"DoubleMetaphone" of empty string is not ""'

    for s in self.strings:

      c = encode.dmetaphone(s)

      assert (isinstance(c,str)), '"DoubleMetaphone" of string "'+s+'"does '+ \
             'not return a string: '+str(type(c))

      assert (len(c) <= 4), 'Length of "DoubleMetaphone" code for string "'+s+\
             '" is more than four characters: '+str(c)+' with length: '+ \
             str(len(c))

      if (len(c) > 1):
        assert (c[1:].isalpha() == 1), 'Characters after first in '+ \
               '"DoubleMetaphone" code for string "'+s+'" are not letters: '+\
               str(c)

      c = encode.dmetaphone(s,maxlen=1)

      assert (isinstance(c,str)), '"DoubleMetaphone" of string "'+s+'"does '+ \
             'not return a string: '+str(type(c))

      assert (len(c) == 1), 'Length of "DoubleMetaphone" code for string "'+s+\
             '" is '+'not one character: '+str(c)+' with length: '+str(len(c))

      c = encode.dmetaphone(s,maxlen=6)

      assert (isinstance(c,str)), '"DoubleMetaphone" of string "'+s+'"does '+ \
             'not return a string: '+str(type(c))

      assert (len(c) <= 6), '"DoubleMetaphone" code for string "'+s+'" is '+ \
             'longer than six characters: '+str(c)+' with length: '+str(len(c))

      if (len(c) > 1):
        assert (c[1:].isalpha() == 1), 'Characters after first in '+ \
               '"DoubleMetaphone" code for string "'+s+'" are not letters: '+ \
               str(c)

  def testGetSubstring(self):  # - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'Get-Substring' string encoding"""

    # A list of input string, start and end values, and expected output string
    #
    test_data = [('peter',0,6,'peter'), ('peter',1,6,'eter'),
                 ('peter',1,6,'eter'), ('peter',2,6,'ter'), ('peter',3,6,'er'),
                 ('peter',4,6,'r'), ('peter',5,6,''), ('peter',6,6,''),
                 ('hello world',0,11,'hello world'),
                 ('hello world',0,10,'hello worl'),
                 ('hello world',0,9,'hello wor'),
                 ('hello world',0,8,'hello wo'), ('hello world',0,7,'hello w'),
                 ('hello world',0,6,'hello '), ('hello world',0,5,'hello'),
                 ('hello world',0,4,'hell'), ('hello world',0,3,'hel'),
                 ('hello world',0,2,'he'), ('hello world',0,1,'h'),
                 ('hello world',0,0,''),
                 ('hello world',0,11,'hello world'),
                 ('hello world',1,10,'ello worl'),
                 ('hello world',3,9,'lo wor'), ('hello world',2,8,'llo wo'),
                 ('hello world',4,7,'o w'), ('hello world',1,6,'ello '),
                 ('hello world',2,5,'llo'), ('hello world',3,4,'l'),
                 ('hello world',2,3,'l'), ('hello world',1,2,'e'),
                 ('hello world',1,1,''), ('hello world',5,5,'')]

    for (in_str, start_index, end_index, out_str) in test_data:

      test_out_str = encode.get_substring(in_str, start_index, end_index)

      assert test_out_str == out_str, \
             'Get-substring returns wrong result (should be "%s"): "%s"' % \
             (out_str, test_out_str)


  def testFreqVector(self):  # - - - - - - - - - - - - - - - - - - - - - - - -
    """Test 'freq_vector' string encoding"""

    c = encode.freq_vector('')  # Test with empty string
    assert (c == ''), \
           '"freq_vector" of empty string is not an empty string'

    for s in self.strings:

      c = encode.freq_vector(s)

      assert (isinstance(c,list)), \
             '"freq_vector" of string "'+s+'" does not return a list: '+ \
             str(type(c))

      assert (len(c) == 26), \
             'Length of frequency vector with no encoding is not 26 for ' + \
             'string "'+s+'": '+str(c)

      for i in c:
        assert (isinstance(i,int) and (i >= 0)), \
               '"freq_vector" of string "'+s+'" does not contain integers:' \
               + str(type(c))

      c = encode.freq_vector(s, 'phonix')

      assert (isinstance(c,list)), \
             '"freq_vector" of string "'+s+'" does not return a list: '+ \
             str(type(c))

      assert (len(c) == 9), \
             'Length of frequency vector with "phonix" encoding is not 9 ' + \
             'for string "'+s+'": '+str(c)

      for i in c:
        assert (isinstance(i,int) and (i >= 0)), \
               '"freq_vector" of string "'+s+'" does not contain integers:' \
               + str(type(c))

      c = encode.freq_vector(s, 'soundex')

      assert (isinstance(c,list)), \
             '"freq_vector" of string "'+s+'" does not return a list: '+ \
             str(type(c))

      assert (len(c) == 7), \
             'Length of frequency vector with "soundex" encoding is not 7 ' + \
             'for string "'+s+'": '+str(c)

      for i in c:
        assert (isinstance(i,int) and (i >= 0)), \
               '"freq_vector" of string "'+s+'" does not contain integers:' \
               + str(type(c))

      c = encode.freq_vector(s, 'mod_soundex')

      assert (isinstance(c,list)), \
             '"freq_vector" of string "'+s+'" does not return a list: '+ \
             str(type(c))

      assert (len(c) == 10), \
             'Length of frequency vector with "mod_soundex" encoding is ' + \
             'not 10 for string "'+s+'": '+str(c)

      for i in c:
        assert (isinstance(i,int) and (i >= 0)), \
               '"freq_vector" of string "'+s+'" does not contain integers:' \
               + str(type(c))

# =============================================================================
# Start tests when called from command line

if (__name__ == "__main__"):

  # Intialise a logger, set level to info
  #
  my_logger = logging.getLogger()  # New logger at root level
  my_logger.setLevel(log_level)

  unittest.main()  # Run all test

# =============================================================================
