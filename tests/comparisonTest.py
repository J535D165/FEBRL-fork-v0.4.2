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
# The Original Software is: "comparisonTest.py"
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

"""Test module for comparison.py.
"""

# =============================================================================
# Import necessary modules (Python standard modules first, then Febrl modules)

import logging
import sys
import unittest
sys.path.append('..')

import comparison
import stringcmp

log_level = logging.WARNING  # logging.INFO

# =============================================================================

class testdataset:  # Used for record comparator testing

  def __init__(self, field_names_list, descr):
    self.field_list = field_names_list
    self.description = descr

# -----------------------------------------------------------------------------

class TestCase(unittest.TestCase):

  # Initialise test case  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def setUp(self):

    # Triplets of weight values (missing, disagreement, agreement)
    #
    self.weight_values = [(0.0, 0.0, 1.0),(0.0,-2.0,2.0),(-1.0,-12.0,20.0),
                          (-7.5,-22.0,4.0)]

    # Define a small geographic look-up table with locations
    #
    self.geo_lookup_table = {'2000':(151.20710, -33.87060),
                             '2010':(151.21524, -33.88299),
                             '2020':(151.19482, -33.92858),
                             '2030':(151.27912, -33.86007),
                             '2880':(141.45874, -31.95809),
                             '3500':(150.84669, -34.56312),
                             '3501':(150.84669, -34.56312),
                             '3585':(143.71791, -35.27436),
                             '3644':(145.75410, -35.87677),
                             '3691':(147.03416, -36.10235),
                             '4385':(150.91642, -28.74569),
                             '4388':(150.74207, -28.67069)}

    # Define pairs of exactly the same values for strings, numbers, dates, etc.
    #
    self.exact_string_pairs = [('peter','peter'),('PAUL','PAUL'),('T','T'),
                               ('a very long example!','a very long example!'),
                               ('x','x'),('123','123')]

    self.exact_number_pairs = [(0,0),('0',0),(0,'0'),('0','0'),(-1,-1),
                               ('-1',-1),(-1,'-1'),('-1','-1'),(123,123),
                               (123,'123'),('123',123),('123','123'),
                               (1.91,1.91),('1.91',1.91),(1.91,'1.91'),
                               ('1.91','1.91'),(0.0,-0.0),('0.0',-0.0),
                               (0.0,'-0.0'),('0.0','-0.0'),(-99.99,-99.99),
                               ('-99.99',-99.99),(-99.99,'-99.99'),
                               ('-99.99','-99.99')]

    self.exact_distance_pairs = [('2000','2000'),('2020','2020'),
                                 ('3500','3501')]

    self.exact_date_pairs = [('09112006','09:11:2006'),
                             ('11092000','11,09,2000'),
                             ('01011919','01/01/1919'),
                             ('10122022','10\\12\\2022')]

    self.exact_age_pairs = [('09112006','09/11/2006'),
                            ('09:01:1890','09\\01\\1890'),
                            ('11031963','11,03,1963'),
                            ('01:08.1900','0108:1900')]

    self.exact_time_pairs = [('1200','1200'),('0000','0000'),('2359','2359'),
                             ('12:21','1221'),('0101','01:01'),
                             ('22:22','22:22')]

    # Define pairs with one value missing
    #
    self.missing_values_list = ['','n/a','missing']  # Set when initialising

    self.missing_string_pairs = [('','peter'),('n/a','peter'),
                                 ('missing','peter'),('peter',''),
                                 ('peter','n/a'),('peter','missing'),
                                 ('',''),('','n/a'),('n/a',''),('n/a','n/a'),
                                 ('','missing'),('n/a','missing'),
                                 ('missing','missing'),('missing','n/a'),
                                 ('missing',''),('','PAUL'),('n/a','PAUL'),
                                 ('missing','PAUL'),('PAUL',''),('PAUL','n/a'),
                                 ('missing','a very long example!'),
                                 ('a very long example!','n/a')]

    self.missing_number_pairs = [('',0),('','0'),(0,''),('0',''),('n/a',0),
                                 ('missing','0'),(0,'n/a'),('0','missing'),
                                 ('',1.0),('','1.0'),(1.0,''),('1.0',''),
                                 ('n/a',1.0),('missing','1.0'),(1.0,'n/a'),
                                 ('1.0','missing'),('',-99.99),('','-99.99'),
                                 (-99.99,''),('-99.99',''),('n/a',-99.99),
                                 ('missing','-99.99'),(-99.99,'n/a'),
                                 ('-99.99','missing')]

    self.missing_distance_pairs = [('2000','2009'),('2113','2020'),
                                   ('2000',''),('n/a','2020'),('3500','3591'),
                                   ('perth','sydney'),('missing','n/a'),
                                   ('',''),('','n/a'),('n/a',''),('n/a','n/a'),
                                   ('','missing'),('n/a','missing'),
                                   ('missing','missing'),('missing','n/a')]

    self.missing_date_pairs = [('','09112006'),('09:11:2006',''),
                               ('n/a','11.09.2000'),
                               ('11/09\\2000','n/a'),
                               ('missing','01011919'),
                               ('01\\01\\1919','missing'),
                               ('',''),('','n/a'),('n/a',''),('n/a','n/a'),
                               ('','missing'),('n/a','missing'),
                               ('missing','missing'),('missing','n/a')]

    self.missing_age_pairs = [('','09012006'),('09;01;2006',''),
                               ('n/a','11:09:2000'),
                               ('11092000','n/a'),
                               ('missing','01:011919'),
                               ('01011919','missing'),
                               ('',''),('','n/a'),('n/a',''),('n/a','n/a'),
                               ('','missing'),('n/a','missing'),
                               ('missing','missing'),('missing','n/a')]

    self.missing_time_pairs = [('1200',''),('','0000'),('',''),('12:21',''),
                               ('','01:01'),('22:22','missing'),
                               ('n/a','22:22'),('',''),('','n/a'),('n/a',''),
                               ('n/a','n/a'),('','missing'),('n/a','missing'),
                               ('missing','missing'),('missing','n/a')]

    # Define pairs with similar values
    #
    self.similar_string_pairs = [('A very long example','A very long exampl'),
                                 ('A multi word example','example word multi'),
                                 ('A long example','a long example?'),
                                 ('PETER','PETAR'),('PAUL','POUL'),
                                 ('123','121')]

    self.similar_number_pairs = [(0,0.1),('0',0.1),(0,'0.1'),('0','0.1'),
                                 (0.2,0),('0.2',0),(0.2,'0'),('0.2','0'),
                                 (-10,-12),('-10',-11),(-11,'-12'),
                                 ('-12','-10'),(121,123),
                                 (123,'126'),('126',123),('125','123'),
                                 (1.90,1.91),('1.91',1.90),(1.92,'1.91'),
                                 ('1.89','1.94'),(0.0,-0.1),('0.0',-0.2),
                                 (0.0,'-0.1'),('0.0','-0.2'),(-99.9,-99.99),
                                 ('-99.99',-99.9),(-99.8,'-99.99'),
                                 ('-99.99','-99.89')]

    self.similar_distance_pairs = [('2000','2010'),('2000','2020'),
                                   ('2000','2030'),('2020','2010'),
                                   ('2010','2030'),('2030','2000'),
                                   ('4385','4388'),('4388','4385')]

    self.similar_date_pairs = [('09112006','10/11/2006'),
                               ('11092000','09,09,2000'),
                               ('01121999','12.01.1999'),
                               ('12091978','12\\07\\1978'),
                               ('01011919','03:01;1919'),
                               ('01012022','12;12:2022')]

    self.similar_age_pairs = [('09\\01\\2006','10012006'),
                              ('11:09;2000','09:09:2000'),
                              ('01121999','12,01,1999'),
                              ('12:09:1978','12,07,1978'),
                              ('01;01,1919','03/01/1919'),
                              ('01,02/1899','12:12:1901')]

    self.similar_time_pairs = [('1200','1155'),('0003','0000'),
                               ('12:21','1222'),('0059','01:01'),
                               ('22:22','22:19')]

    # Define a sequence of pairs with similar values (start with a same values
    # pair)
    #
    self.similar_string_seq = [('peter miller','peter miller'),
                               ('peter miller','peter miler'),
                               ('peter miller','pedro miller'),
                               ('peter miller','pedro milier'),
                               ('peter miller','joseph seier')]

    self.similar_number_seq = [(42.24,42.24),
                               (42.24,42.20),
                               (42.24,42.00),
                               (42.24,40.20),
                               (42.24,10),
                               (42.24,-42.20)]

    self.similar_distance_seq = [('2000','2000'),
                                 ('2000','2010'),
                                 ('2000','2020'),
                                 ('2000','2880'),
                                 ('2000','3500'),
                                 ('2000','3585')]

    self.similar_date_seq = [('09:01:2006','09,11,2006'),
                             ('09:01:2006','08,11,2006'),
                             ('09:01:2006','07/11/2006'),
                             ('09:01:2006','06.11.2006'),
                             ('09:01:2006','05;11;2006'),
                             ('09:01:2006','01:11:2006')]

    self.similar_age_seq = [('09:01:2006','09,01,2006'),
                            ('09:01:2006','01,01,2006'),
                            ('09:01:2006','01,10,2005'),
                            ('09:01:2006','01,01,2005'),
                            ('09:01:2006','01,01,2000'),
                            ('09:01:2006','11,12,1999'),
                            ('09:01:2006','01,11,1968'),
                            ('09:01:2006','03,03,1919')]

    self.similar_time_seq = [('1200','1200'),
                             ('1200','1155'),
                             ('1200','1147'),
                             ('1200','1130'),
                             ('1200','1031'),
                             ('1200','0909'),
                             ('1200','1731')]

    # Define pairs with totally different values
    #
    self.different_string_pairs = [('XRZ','BADL'),('peter','zyc')]

    self.different_number_pairs = [(-100,9999),(99999.98,-99999.98),
                                   (0,99999.999),(99999.999,0),(0,-99999.999),
                                   (-99999.999,0)]

    self.different_distance_pairs = [('2000','4388'),('2880','3691'),
                                     ('3501','3585')]

    self.different_date_pairs = [('09,11,2006','10,11,1971'),
                                ('11,09,1963','01,04,2005'),
                                ('09,09,1919','01,01,1979'),
                                ('01,12,1881','02,12,2002')]

    self.different_age_pairs = [('09,11,2006','10,11,1971'),
                                ('11,09,1963','01,04,2005'),
                                ('09,09,1919','01,01,1979'),
                                ('01,12,1881','02,12,2002')]

    self.different_time_pairs = [('0000','12:01'),('23:11','1109'),
                                 ('0315','1646'),('14:14','01:59')]

    # Define two test 'data sets' for record comparator testing - - - - - - - -
    #
    self.test_data_set1 = testdataset([('rec-id',0),('gname',1),('surname',2),
                                       ('streetnumb',3),('streetname_type',4),
                                       ('suburb',5),('postcode',6)],
                                      'Test data set 1')
    self.recs1 = [['rec-0-org','james','whiteway','2','?','red hill','2611'],
                  ['rec-1-org','yasmin','weidenbach','51','?','deakin','2602'],
                  ['rec-2-org','anari','astley','204','?','yarralumla','2905'],
                  ['rec-3-org','mitchell','devin','7','?','holder','2606']]

    self.test_data_set2 = testdataset([('unused',0),('rec-id',1),('sname',2),
                                      ('given_name',3),('street',4),
                                      ('locality',5),('zipcode',6)],
                                      'Test data set 2')
    self.recs2 =[['','rec-0-org','whiteway','james','2','red hill','2611'],
                 ['','rec-0-dup','whitway','jams','z','red hill','26l1'],
                 ['','rec-1-dup','weidenbach','yasi','15','deakin','260'],
                 ['','rec-2-dup','astly','anary','20a','yaraluma','2905'],
                 ['','rec-3-dup','deivn','michell','8','holden','2066']]

    self.contain_string_pairs = [('peter','peter christen'),
                                 ('peter christen','peter'),
                                 ('peter','christen peter'),
                                 ('christen peter','peter'),
                                 ('1', '1234567890'),
                                 ('1234567890', '1'),
                                 ('4', '1234567890'),
                                 ('1234567890', '4'),
                                 ('0', '1234567890'),
                                 ('1234567890', '0'),
                                 ('90', '1234567890'),
                                 ('1234567890', '90'),
                                 ('one', 'one two three'),
                                 ('one two three', 'one'),
                                 ('two', 'one two three'),
                                 ('one two three', 'two'),
                                 ('three', 'one two three'),
                                 ('one two three', 'three'),
                                 (' ', 'one two three'),
                                 ('one two three', ' '),
                                 ('e t', 'one two three'),
                                 ('one two three', 'e t')]

    self.not_contain_string_pairs = [('pedro','peter christen'),
                                     ('peter christen','pedro'),
                                     ('peter','christen petra'),
                                     ('christen petra','peter'),
                                     ('x', '1234567890'),
                                     ('1234567890', 'x'),
                                     ('4', '123x567890'),
                                     ('123x567890', '4'),
                                     ('0', '123456789x'),
                                     ('123456789x', '0'),
                                     ('90', '1234567891'),
                                     ('1234567891', '90'),
                                     ('four', 'one two three'),
                                     ('one two three', 'four'),
                                     ('twoo', 'one two three'),
                                     ('one two three', 'twoo'),
                                     ('tree', 'one two three'),
                                     ('one two three', 'tree'),
                                     ('_', 'one two three'),
                                     ('one two three', '_'),
                                     ('e-t', 'one two three'),
                                     ('one two three', 'e-t')]

  # Clean up test case  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def tearDown(self):
    pass  # Nothing to clean up

  # ---------------------------------------------------------------------------
  # Define tests for exact, contains, string, number, date, age, etc. field
  # comparators
  #
  def doExactFieldComparisonTest(self, fc):

    for (mw, daw, aw) in self.weight_values:

      fc.cache = {}  # Clear old cache entries

      fc.set_weights(missing_w = mw, agree_w = aw, disagree_w = daw)

      for (val1, val2) in self.exact_string_pairs:

        aw_test = fc.compare(val1, val2)
        assert aw_test == aw, \
               'Wrong agreement weight for "%s" (should be: %f but is %f)' % \
               (fc.description, aw, aw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        aw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert aw_test == aw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (aw_test, aw_test2, str(val1), str(val2))

      for (val1, val2) in self.missing_string_pairs:

        mw_test = fc.compare(val1, val2)
        assert mw_test == mw, \
               'Wrong missing weight for "%s" (should be: %f but is %f)' % \
               (fc.description, mw, mw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        mw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert mw_test == mw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (mw_test, mw_test2, str(val1), str(val2))

      for (val1, val2) in self.similar_string_pairs:

        paw_test = fc.compare(val1, val2)
        assert paw_test == daw, \
               'Wrong disagreement weight for "%s" (should be: %f but is %f)' \
               % (fc.description, daw, daw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        paw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert paw_test == paw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (paw_test, paw_test2, str(val1), str(val2))

      for (val1, val2) in self.different_string_pairs:

        daw_test = fc.compare(val1, val2)
        assert daw_test == daw, \
               'Wrong disagreement weight for "%s" (should be: %f but is %f)' \
               % (fc.description, daw, daw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        daw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert daw_test == daw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (daw_test, daw_test2, str(val1), str(val2))

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  def doContainsFieldComparisonTest(self, fc):

    for (mw, daw, aw) in self.weight_values:

      fc.cache = {}  # Clear old cache entries

      fc.set_weights(missing_w = mw, agree_w = aw, disagree_w = daw)

      for (val1, val2) in self.contain_string_pairs:

        aw_test = fc.compare(val1, val2)
        assert aw_test == aw, \
               'Wrong agreement weight for "%s" (should be: %f but is %f)' % \
               (fc.description, aw, aw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        aw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert aw_test == aw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (aw_test, aw_test2, str(val1), str(val2))

      for (val1, val2) in self.missing_string_pairs:

        mw_test = fc.compare(val1, val2)
        assert mw_test == mw, \
               'Wrong missing weight for "%s" (should be: %f but is %f)' % \
               (fc.description, mw, mw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        mw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert mw_test == mw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (mw_test, mw_test2, str(val1), str(val2))

      for (val1, val2) in self.not_contain_string_pairs:

        paw_test = fc.compare(val1, val2)
        assert paw_test == daw, \
               'Wrong disagreement weight for "%s" (should be: %f but is %f)' \
               % (fc.description, daw, daw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        paw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert paw_test == paw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (paw_test, paw_test2, str(val1), str(val2))

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  def doStringFieldComparisonTest(self, fc):

    for (mw, daw, aw) in self.weight_values:

      fc.cache = {}  # Clear old cache entries

      fc.set_weights(missing_w = mw, agree_w = aw, disagree_w = daw)

      for (val1, val2) in self.exact_string_pairs:
        aw_test = fc.compare(val1, val2)
        assert aw_test == aw, \
               'Wrong agreement weight for "%s" (should be: %f but is %f)' % \
               (fc.description, aw, aw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        aw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert aw_test == aw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (aw_test, aw_test2, str(val1), str(val2))

      for (val1, val2) in self.missing_string_pairs:

        mw_test = fc.compare(val1, val2)
        assert mw_test == mw, \
               'Wrong missing weight for "%s" (should be: %f but is %f)' % \
               (fc.description, mw, mw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        mw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert mw_test == mw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (mw_test, mw_test2, str(val1), str(val2))

      for (val1, val2) in self.similar_string_pairs:

        paw_test = fc.compare(val1, val2)
        assert paw_test <= aw, \
               'Partial agreement weight for "%s" is larger than agreement' % \
               (fc.description)+' weight (should be equal to or smaller ' + \
               'than %f but is %f)' % (aw, paw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))
        assert paw_test >= daw, \
               'Partial agreement weight for "%s" is smaller than ' % \
               (fc.description)+'disagreement weight (should be equal to ' + \
               'or larger than %f but is %f)' % (daw, paw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        paw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert paw_test == paw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (paw_test, paw_test2, str(val1), str(val2))

      prev_test = aw+1
      old_vals = ('','')
      for (val1, val2) in self.similar_string_seq:
        paw_test = fc.compare(val1, val2)
        assert paw_test <= aw, \
               'Partial agreement weight for "%s" is larger than agreement' % \
               (fc.description)+' weight (should be equal to or smaller ' + \
               'than %f but is %f)' % (aw, paw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))
        assert paw_test >= daw, \
               'Partial agreement weight for "%s" is smaller than ' % \
               (fc.description)+'disagreement weight (should be equal to ' + \
               'or larger than %f but is %f)' % (daw, paw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))
        assert paw_test <= prev_test, \
               'Partial agreement weight for "%s" is larger than the ' % \
               (fc.description)+'partial agreement weight of previous pair' + \
               ' (%s,%s) / (%s,%s)' % (str(old_vals[0]), str(old_vals[1]), \
               str(val1), str(val2)) + ': %.3f / %.3f' % (prev_test,paw_test)

        prev_test = paw_test
        old_vals = (val1,val2)

      for (val1, val2) in self.different_string_pairs:

        daw_test = fc.compare(val1, val2)

        # For the compression comparator this cannot be guaranteed!
        #
        if (('Compress' not in fc.description) and \
            ('compress' not in fc.description)):
          assert daw_test == daw, \
               'Wrong disagreement weight for "%s" (should be: %f but is %f)' \
               % (fc.description, daw, daw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        daw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert daw_test == daw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (daw_test, daw_test2, str(val1), str(val2))

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  def doNumericFieldComparisonTest(self, fc):

    for (mw, daw, aw) in self.weight_values:

      fc.cache = {}  # Clear old cache entries

      fc.set_weights(missing_w = mw, agree_w = aw, disagree_w = daw)

      for (val1, val2) in self.exact_number_pairs:
        aw_test = fc.compare(val1, val2)
        assert aw_test == aw, \
               'Wrong agreement weight for "%s" (should be: %f but is %f)' % \
               (fc.description, aw, aw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        aw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert aw_test == aw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (aw_test, aw_test2, str(val1), str(val2))

      for (val1, val2) in self.missing_number_pairs:

        mw_test = fc.compare(val1, val2)
        assert mw_test == mw, \
               'Wrong missing weight for "%s" (should be: %f but is %f)' % \
               (fc.description, mw, mw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        mw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert mw_test == mw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (mw_test, mw_test2, str(val1), str(val2))

      for (val1, val2) in self.similar_number_pairs:

        paw_test = fc.compare(val1, val2)
        assert paw_test <= aw, \
               'Partial agreement weight for "%s" is larger than agreement' % \
               (fc.description)+' weight (should be equal to or smaller ' + \
               'than %f but is %f)' % (aw, paw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))
        assert paw_test >= daw, \
               'Partial agreement weight for "%s" is smaller than ' % \
               (fc.description)+'disagreement weight (should be equal to ' + \
               'or larger than %f but is %f)' % (daw, paw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        paw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert paw_test == paw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (paw_test, paw_test2, str(val1), str(val2))

      prev_test = aw+1
      old_vals = ('','')
      for (val1, val2) in self.similar_number_seq:
        paw_test = fc.compare(val1, val2)
        assert paw_test <= aw, \
               'Partial agreement weight for "%s" is larger than agreement' % \
               (fc.description)+' weight (should be equal to or smaller ' + \
               'than %f but is %f)' % (aw, paw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))
        assert paw_test >= daw, \
               'Partial agreement weight for "%s" is smaller than ' % \
               (fc.description)+'disagreement weight (should be equal to ' + \
               'or larger than %f but is %f)' % (daw, paw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        assert paw_test <= prev_test, \
               'Partial agreement weight for "%s" is larger than the ' % \
               (fc.description)+'partial agreement weight of previous pair' + \
               ' (%s,%s) / (%s,%s)' % (str(old_vals[0]), str(old_vals[1]), \
               str(val1), str(val2)) + ': %.3f / %.3f' % (prev_test,paw_test)

        prev_test = paw_test
        old_vals = (val1,val2)

      for (val1, val2) in self.different_number_pairs:

        daw_test = fc.compare(val1, val2)
        assert daw_test == daw, \
               'Wrong disagreement weight for "%s" (should be: %f but is %f)' \
               % (fc.description, daw, daw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        daw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert daw_test == daw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (daw_test, daw_test2, str(val1), str(val2))

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  def doDistanceFieldComparisonTest(self, fc):

    for (mw, daw, aw) in self.weight_values:

      fc.cache = {}  # Clear old cache entries

      fc.set_weights(missing_w = mw, agree_w = aw, disagree_w = daw)

      for (val1, val2) in self.exact_distance_pairs:

        aw_test = fc.compare(val1, val2)
        assert aw_test == aw, \
               'Wrong agreement weight for "%s" (should be: %f but is %f)' % \
               (fc.description, aw, aw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        aw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert abs(aw_test-aw_test2) < 0.00000001, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (aw_test, aw_test2, str(val1), str(val2))

      for (val1, val2) in self.missing_distance_pairs:

        mw_test = fc.compare(val1, val2)
        assert mw_test == mw, \
               'Wrong missing weight for "%s" (should be: %f but is %f)' % \
               (fc.description, mw, mw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        mw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert abs(mw_test-mw_test2) < 0.00000001, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (mw_test, mw_test2, str(val1), str(val2))

      for (val1, val2) in self.similar_distance_pairs:

        paw_test = fc.compare(val1, val2)
        assert paw_test <= aw, \
               'Partial agreement weight for "%s" is larger than agreement' % \
               (fc.description)+' weight (should be equal to or smaller ' + \
               'than %f but is %f)' % (aw, paw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))
        assert paw_test >= daw, \
               'Partial agreement weight for "%s" is smaller than ' % \
               (fc.description)+'disagreement weight (should be equal to ' + \
               'or larger than %f but is %f)' % (daw, paw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        paw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert abs(paw_test-paw_test2) < 0.00000001, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (paw_test, paw_test2, str(val1), str(val2))

      prev_test = aw+1
      old_vals = ('','')
      for (val1, val2) in self.similar_distance_seq:
        paw_test = fc.compare(val1, val2)
        assert paw_test <= aw, \
               'Partial agreement weight for "%s" is larger than agreement' % \
               (fc.description)+' weight (should be equal to or smaller ' + \
               'than %f but is %f)' % (aw, paw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))
        assert paw_test >= daw, \
               'Partial agreement weight for "%s" is smaller than ' % \
               (fc.description)+'disagreement weight (should be equal to ' + \
               'or larger than %f but is %f)' % (daw, paw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        assert paw_test <= prev_test, \
               'Partial agreement weight for "%s" is larger than the ' % \
               (fc.description)+'partial agreement weight of previous pair' + \
               ' (%s,%s) / (%s,%s)' % (str(old_vals[0]), str(old_vals[1]), \
               str(val1), str(val2)) + ': %.3f / %.3f' % (prev_test,paw_test)

        prev_test = paw_test
        old_vals = (val1,val2)

      for (val1, val2) in self.different_distance_pairs:

        daw_test = fc.compare(val1, val2)
        assert daw_test == daw, \
               'Wrong disagreement weight for "%s" (should be: %f but is %f)' \
               % (fc.description, daw, daw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        daw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert abs(daw_test-daw_test2) < 0.00000001, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (daw_test, daw_test2, str(val1), str(val2))

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  def doDateFieldComparisonTest(self, fc):

    for (mw, daw, aw) in self.weight_values:

      fc.cache = {}  # Clear old cache entries

      fc.set_weights(missing_w = mw, agree_w = aw, disagree_w = daw)

      for (val1, val2) in self.exact_date_pairs:

        aw_test = fc.compare(val1, val2)
        assert aw_test == aw, \
               'Wrong agreement weight for "%s" (should be: %f but is %f)' % \
               (fc.description, aw, aw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        aw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert aw_test == aw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (aw_test, aw_test2, str(val1), str(val2))

      for (val1, val2) in self.missing_date_pairs:

        mw_test = fc.compare(val1, val2)
        assert mw_test == mw, \
               'Wrong missing weight for "%s" (should be: %f but is %f)' % \
               (fc.description, mw, mw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        mw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert mw_test == mw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (mw_test, mw_test2, str(val1), str(val2))

      for (val1, val2) in self.similar_date_pairs:

        paw_test = fc.compare(val1, val2)
        assert paw_test <= aw, \
               'Partial agreement weight for "%s" is larger than agreement' % \
               (fc.description)+' weight (should be equal to or smaller ' + \
               'than %f but is %f)' % (aw, paw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))
        assert paw_test >= daw, \
               'Partial agreement weight for "%s" is smaller than ' % \
               (fc.description)+'disagreement weight (should be equal to ' + \
               'or larger than %f but is %f)' % (daw, paw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        paw_test = fc.compare(val2, val1)
        assert paw_test <= aw, \
               'Partial agreement weight for "%s" is larger than agreement' % \
               (fc.description)+' weight (should be equal to or smaller ' + \
               'than %f but is %f)' % (aw, paw_test) + \
               ' with values (%s,%s)' % (str(val2), str(val1))
        assert paw_test >= daw, \
               'Partial agreement weight for "%s" is smaller than ' % \
               (fc.description)+'disagreement weight (should be equal to ' + \
               'or larger than %f but is %f)' % (daw, paw_test) + \
               ' with values (%s,%s)' % (str(val2), str(val1))

      prev_test = aw+1
      old_vals = ('','')
      for (val1, val2) in self.similar_date_seq:
        paw_test = fc.compare(val1, val2)
        assert paw_test <= aw, \
               'Partial agreement weight for "%s" is larger than agreement' % \
               (fc.description)+' weight (should be equal to or smaller ' + \
               'than %f but is %f)' % (aw, paw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))
        assert paw_test >= daw, \
               'Partial agreement weight for "%s" is smaller than ' % \
               (fc.description)+'disagreement weight (should be equal to ' + \
               'or larger than %f but is %f)' % (daw, paw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        assert paw_test <= prev_test, \
               'Partial agreement weight for "%s" is larger than the ' % \
               (fc.description)+'partial agreement weight of previous pair' + \
               ' (%s,%s) / (%s,%s)' % (str(old_vals[0]), str(old_vals[1]), \
               str(val1), str(val2)) + ': %.3f / %.3f' % (prev_test,paw_test)+\
               ' %d/%d' %(fc.max_day1_before_day2,fc.max_day2_before_day1)

        prev_test = paw_test
        old_vals = (val1,val2)

      for (val1, val2) in self.different_date_pairs:

        daw_test = fc.compare(val1, val2)
        assert daw_test == daw, \
               'Wrong disagreement weight for "%s" (should be: %f but is %f)' \
               % (fc.description, daw, daw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        daw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert daw_test == daw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (daw_test, daw_test2, str(val1), str(val2))


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  def doTimeFieldComparisonTest(self, fc):

    for (mw, daw, aw) in self.weight_values:

      fc.cache = {}  # Clear old cache entries

      fc.set_weights(missing_w = mw, agree_w = aw, disagree_w = daw)

      for (val1, val2) in self.exact_time_pairs:

        aw_test = fc.compare(val1, val2)
        assert aw_test == aw, \
               'Wrong agreement weight for "%s" (should be: %f but is %f)' % \
               (fc.description, aw, aw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        aw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert aw_test == aw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (aw_test, aw_test2, str(val1), str(val2))

      for (val1, val2) in self.missing_time_pairs:

        mw_test = fc.compare(val1, val2)
        assert mw_test == mw, \
               'Wrong missing weight for "%s" (should be: %f but is %f)' % \
               (fc.description, mw, mw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        mw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert mw_test == mw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (mw_test, mw_test2, str(val1), str(val2))

      for (val1, val2) in self.similar_time_pairs:

        paw_test = fc.compare(val1, val2)
        assert paw_test <= aw, \
               'Partial agreement weight for "%s" is larger than agreement' % \
               (fc.description)+' weight (should be equal to or smaller ' + \
               'than %f but is %f)' % (aw, paw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))
        assert paw_test >= daw, \
               'Partial agreement weight for "%s" is smaller than ' % \
               (fc.description)+'disagreement weight (should be equal to ' + \
               'or larger than %f but is %f)' % (daw, paw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        paw_test = fc.compare(val2, val1)
        assert paw_test <= aw, \
               'Partial agreement weight for "%s" is larger than agreement' % \
               (fc.description)+' weight (should be equal to or smaller ' + \
               'than %f but is %f)' % (aw, paw_test) + \
               ' with values (%s,%s)' % (str(val2), str(val1))
        assert paw_test >= daw, \
               'Partial agreement weight for "%s" is smaller than ' % \
               (fc.description)+'disagreement weight (should be equal to ' + \
               'or larger than %f but is %f)' % (daw, paw_test) + \
               ' with values (%s,%s)' % (str(val2), str(val1))

      prev_test = aw+1
      old_vals = ('','')
      for (val1, val2) in self.similar_time_seq:
        paw_test = fc.compare(val1, val2)
        assert paw_test <= aw, \
               'Partial agreement weight for "%s" is larger than agreement' % \
               (fc.description)+' weight (should be equal to or smaller ' + \
               'than %f but is %f)' % (aw, paw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))
        assert paw_test >= daw, \
               'Partial agreement weight for "%s" is smaller than ' % \
               (fc.description)+'disagreement weight (should be equal to ' + \
               'or larger than %f but is %f)' % (daw, paw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        assert paw_test <= prev_test, \
               'Partial agreement weight for "%s" is larger than the ' % \
               (fc.description)+'partial agreement weight of previous pair' + \
               ' (%s,%s) / (%s,%s)' % (str(old_vals[0]), str(old_vals[1]), \
               str(val1), str(val2)) + ': %.3f / %.3f' % (prev_test,paw_test)

        prev_test = paw_test
        old_vals = (val1,val2)

      for (val1, val2) in self.different_time_pairs:

        daw_test = fc.compare(val1, val2)
        assert daw_test == daw, \
               'Wrong disagreement weight for "%s" (should be: %f but is %f)' \
               % (fc.description, daw, daw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        daw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert daw_test == daw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (daw_test, daw_test2, str(val1), str(val2))

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  def doAgeFieldComparisonTest(self, fc):

    for (mw, daw, aw) in self.weight_values:

      fc.cache = {}  # Clear old cache entries

      fc.set_weights(missing_w = mw, agree_w = aw, disagree_w = daw)

      for (val1, val2) in self.exact_age_pairs:

        aw_test = fc.compare(val1, val2)
        assert aw_test == aw, \
               'Wrong agreement weight for "%s" (should be: %f but is %f)' % \
               (fc.description, aw, aw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        aw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert aw_test == aw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (aw_test, aw_test2, str(val1), str(val2))

      for (val1, val2) in self.missing_age_pairs:

        mw_test = fc.compare(val1, val2)
        assert mw_test == mw, \
               'Wrong missing weight for "%s" (should be: %f but is %f)' % \
               (fc.description, mw, mw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        mw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert mw_test == mw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (mw_test, mw_test2, str(val1), str(val2))

      for (val1, val2) in self.similar_age_pairs:

        paw_test = fc.compare(val1, val2)
        assert paw_test <= aw, \
               'Partial agreement weight for "%s" is larger than agreement' % \
               (fc.description)+' weight (should be equal to or smaller ' + \
               'than %f but is %f)' % (aw, paw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))
        assert paw_test >= daw, \
               'Partial agreement weight for "%s" is smaller than ' % \
               (fc.description)+'disagreement weight (should be equal to ' + \
               'or larger than %f but is %f)' % (daw, paw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        paw_test = fc.compare(val2, val1)
        assert paw_test <= aw, \
               'Partial agreement weight for "%s" is larger than agreement' % \
               (fc.description)+' weight (should be equal to or smaller ' + \
               'than %f but is %f)' % (aw, paw_test) + \
               ' with values (%s,%s)' % (str(val2), str(val1))
        assert paw_test >= daw, \
               'Partial agreement weight for "%s" is smaller than ' % \
               (fc.description)+'disagreement weight (should be equal to ' + \
               'or larger than %f but is %f)' % (daw, paw_test) + \
               ' with values (%s,%s)' % (str(val2), str(val1))

        paw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert paw_test == paw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (paw_test, paw_test2, str(val1), str(val2))

      prev_test = aw+1
      old_vals = ('','')
      for (val1, val2) in self.similar_age_seq:
        paw_test = fc.compare(val1, val2)
        assert paw_test <= aw, \
               'Partial agreement weight for "%s" is larger than agreement' % \
               (fc.description)+' weight (should be equal to or smaller ' + \
               'than %f but is %f)' % (aw, paw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))
        assert paw_test >= daw, \
               'Partial agreement weight for "%s" is smaller than ' % \
               (fc.description)+'disagreement weight (should be equal to ' + \
               'or larger than %f but is %f)' % (daw, paw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        assert paw_test <= prev_test, \
               'Partial agreement weight for "%s" is larger than the ' % \
               (fc.description)+'partial agreement weight of previous pair' + \
               ' (%s,%s) / (%s,%s)' % (str(old_vals[0]), str(old_vals[1]), \
               str(val1), str(val2)) + ': %.3f / %.3f' % (prev_test,paw_test)

        prev_test = paw_test
        old_vals = (val1,val2)

      for (val1, val2) in self.different_age_pairs:

        daw_test = fc.compare(val1, val2)
        assert daw_test == daw, \
               'Wrong disagreement weight for "%s" (should be: %f but is %f)' \
               % (fc.description, daw, daw_test) + \
               ' with values (%s,%s)' % (str(val1), str(val2))

        daw_test2 = fc.compare(val2, val1)  # Make sure comparison is symmetric
        assert daw_test == daw_test2, \
               'Different, non-symmetric comparison results for "%s": ' % \
               (fc.description)+' %f and %f with values (%s,%s)' % \
               (daw_test, daw_test2, str(val1), str(val2))


  # ---------------------------------------------------------------------------
  # Test calling functions for exact, string, number, etc. field comparators

  def testExactComparison(self):  # - - - - - - - - - - - - - - - - - - - - - -

    efc = comparison.FieldComparatorExactString(missing_v = \
                                                      self.missing_values_list,
                                           desc = 'FieldComparatorExactString')
    efcc = comparison.FieldComparatorExactString(missing_v = \
                                                      self.missing_values_list,
                                           desc = 'FieldComparatorExactString',
                                           do_cache = True)
    self.doExactFieldComparisonTest(efc)
    self.doExactFieldComparisonTest(efcc)  # Use cache

  def testContainsComparison(self):  # - - - - - - - - - - - - - - - - - - - -

    cfc = comparison.FieldComparatorContainsString(missing_v = \
                                                      self.missing_values_list,
                                           desc = 'FieldComparatorExactString')
    cfcc = comparison.FieldComparatorContainsString(missing_v = \
                                                      self.missing_values_list,
                                           desc = 'FieldComparatorExactString',
                                           do_cache = True)
    self.doContainsFieldComparisonTest(cfc)
    self.doContainsFieldComparisonTest(cfcc)  # Use cache

  def testApproxStringComparison(self):  # - - - - - - - - - - - - - - - - - -

    for c in range(1,6):

      tfc = comparison.FieldComparatorTruncateString(num_char_co = c,
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorTruncateString')
      tfcc = comparison.FieldComparatorTruncateString(num_char_co = c,
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorTruncateString',
                                        do_cache = True)
      self.doStringFieldComparisonTest(tfc)
      self.doStringFieldComparisonTest(tfcc)  # Use cache

    for k in range(1,4):

      kfc = comparison.FieldComparatorKeyDiff(max_key_di = k,
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorKeyDiff')
      kfcc = comparison.FieldComparatorKeyDiff(max_key_di = k,
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorKeyDiff',
                                        do_cache = True)
      self.doStringFieldComparisonTest(kfc)
      self.doStringFieldComparisonTest(kfcc)

    for em in ['soundex','mod_soundex','phonex','phonix','nysiis',
               'dmetaphone','fuzzysoundex']:
      for mcl in [1,2,3,4,5,6,7]:

        efc = comparison.FieldComparatorEncodeString(encode_method = em,
                                        max_code_l = mcl,
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorEncodeString')
        efcr = comparison.FieldComparatorEncodeString(encode_method = em,
                                        max_code_l = mcl,
                                        reverse = True,
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorEncodeString')

        efcc = comparison.FieldComparatorEncodeString(encode_method = em,
                                        max_code_l = mcl,
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorEncodeString',
                                        do_cache = True)
        efcrc = comparison.FieldComparatorEncodeString(encode_method = em,
                                        max_code_l = mcl,
                                        reverse = True,
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorEncodeString',
                                        do_cache = True)

        self.doStringFieldComparisonTest(efc)
        self.doStringFieldComparisonTest(efcr)
        self.doStringFieldComparisonTest(efcc)
        self.doStringFieldComparisonTest(efcrc)

    for t in [0.1, 0.25, 0.5, 0.6666, 0.8, 0.9]:

      chfc = comparison.FieldComparatorCharHistogram(threshold = t,
                                         missing_v = self.missing_values_list,
                                         desc = 'FieldComparatorCharHistogram')
      chfcc = comparison.FieldComparatorCharHistogram(threshold = t,
                                         missing_v = self.missing_values_list,
                                         desc = 'FieldComparatorCharHistogram',
                                         do_cache = True)
      self.doStringFieldComparisonTest(chfc)
      self.doStringFieldComparisonTest(chfcc)

      for comp_f in ['equal', stringcmp.jaro, stringcmp.qgram, stringcmp.lcs]:

        tljfc = comparison.FieldComparatorTwoLevelJaro(threshold = t,
                                          missing_v = self.missing_values_list,
                                          comp_funct = comp_f,
                                          min_thresh = 0.75,
                                          desc = 'FieldComparatorTwoLevelJaro')
        tljfcc = comparison.FieldComparatorTwoLevelJaro(threshold = t,
                                          missing_v = self.missing_values_list,
                                          comp_funct = comp_f,
                                          min_thresh = 0.75,
                                          desc = 'FieldComparatorTwoLevelJaro',
                                          do_cache = True)
        self.doStringFieldComparisonTest(tljfc)
        self.doStringFieldComparisonTest(tljfcc)

      jfc = comparison.FieldComparatorJaro(threshold = t,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorJaro')
      jfcc = comparison.FieldComparatorJaro(threshold = t,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorJaro',
                                          do_cache = True)
      self.doStringFieldComparisonTest(jfc)
      self.doStringFieldComparisonTest(jfcc)

      wfc = comparison.FieldComparatorWinkler(threshold = t,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorWinkler')
      wfcc = comparison.FieldComparatorWinkler(threshold = t,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorWinkler',
                                          do_cache = True)

      wfcl = comparison.FieldComparatorWinkler(threshold = t,
                                          check_long = False,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorWinkler')
      wfccl = comparison.FieldComparatorWinkler(threshold = t,
                                          check_long = False,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorWinkler',
                                          do_cache = True)

      wfci = comparison.FieldComparatorWinkler(threshold = t,
                                          check_init = False,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorWinkler')
      wfcci = comparison.FieldComparatorWinkler(threshold = t,
                                          check_init = False,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorWinkler',
                                          do_cache = True)

      wfcs = comparison.FieldComparatorWinkler(threshold = t,
                                          check_sim = False,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorWinkler')
      wfccs = comparison.FieldComparatorWinkler(threshold = t,
                                          check_sim = False,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorWinkler',
                                          do_cache = True)

      wfcli = comparison.FieldComparatorWinkler(threshold = t,
                                          check_long = False,
                                          check_init = False,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorWinkler')
      wfccli = comparison.FieldComparatorWinkler(threshold = t,
                                          check_long = False,
                                          check_init = False,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorWinkler',
                                          do_cache = True)

      wfclis = comparison.FieldComparatorWinkler(threshold = t,
                                          check_long = False,
                                          check_init = False,
                                          check_sim = False,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorWinkler')
      wfcclis = comparison.FieldComparatorWinkler(threshold = t,
                                          check_long = False,
                                          check_init = False,
                                          check_sim = False,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorWinkler',
                                          do_cache = True)

      wfcs = comparison.FieldComparatorWinkler(threshold = t,
                                          multi_w = 'sort',
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorWinkler')
      wfccs = comparison.FieldComparatorWinkler(threshold = t,
                                          multi_w = 'sort',
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorWinkler',
                                          do_cache = True)
      wfcp = comparison.FieldComparatorWinkler(threshold = t,
                                          multi_w = 'perm',
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorWinkler')
      wfccp = comparison.FieldComparatorWinkler(threshold = t,
                                          multi_w = 'perm',
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorWinkler',
                                          do_cache = True)

      self.doStringFieldComparisonTest(wfc)
      self.doStringFieldComparisonTest(wfcc)
      self.doStringFieldComparisonTest(wfcl)
      self.doStringFieldComparisonTest(wfccl)
      self.doStringFieldComparisonTest(wfci)
      self.doStringFieldComparisonTest(wfcci)
      self.doStringFieldComparisonTest(wfcs)
      self.doStringFieldComparisonTest(wfccs)
      self.doStringFieldComparisonTest(wfcli)
      self.doStringFieldComparisonTest(wfccli)
      self.doStringFieldComparisonTest(wfclis)
      self.doStringFieldComparisonTest(wfcclis)

      self.doStringFieldComparisonTest(wfcs)
      self.doStringFieldComparisonTest(wfccs)
      self.doStringFieldComparisonTest(wfcp)
      self.doStringFieldComparisonTest(wfccp)

      for q in [1,2,3,4,5]:

        for c in ['average','shortest','longest']:

          for stop_word_list in [[],['peter'],['and','or','peter'],
                                 ['PAUL','very','x','123','2000']]:

            tfc = comparison.FieldComparatorTokenSet(threshold = t,
                                          common_div = c,
                                          stop_word_list = stop_word_list,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorTokenSet')

          tfcc = comparison.FieldComparatorTokenSet(threshold = t,
                                          common_div = c,
                                          stop_word_list = stop_word_list,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorTokenSet',
                                          do_cache = True)

          self.doStringFieldComparisonTest(tfc)
          self.doStringFieldComparisonTest(tfcc)

          qfc = comparison.FieldComparatorQGram(threshold = t,
                                          q = q,
                                          common_div = c,
                                          padded = False,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorQGram')

          qfcp = comparison.FieldComparatorQGram(threshold = t,
                                          q = q,
                                          common_div = c,
                                          padded = True,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorQGram')

          qfcc = comparison.FieldComparatorQGram(threshold = t,
                                          q = q,
                                          common_div = c,
                                          padded = False,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorQGram',
                                          do_cache = True)

          qfccp = comparison.FieldComparatorQGram(threshold = t,
                                          q = q,
                                          common_div = c,
                                          padded = True,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorQGram',
                                          do_cache = True)

          self.doStringFieldComparisonTest(qfc)
          self.doStringFieldComparisonTest(qfcp)
          self.doStringFieldComparisonTest(qfcc)
          self.doStringFieldComparisonTest(qfccp)

          for max_d in [1,2,3,4]:

            pqfc = comparison.FieldComparatorPosQGram(threshold = t,
                                          q = q,
                                          max_dist = max_d,
                                          common_div = c,
                                          padded = False,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorPosQGram')

            pqfcp = comparison.FieldComparatorPosQGram(threshold = t,
                                          q = q,
                                          max_dist = max_d,
                                          common_div = c,
                                          padded = True,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorPosQGram')

            pqfcc = comparison.FieldComparatorPosQGram(threshold = t,
                                          q = q,
                                          max_dist = max_d,
                                          common_div = c,
                                          padded = False,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorPosQGram',
                                          do_cache = True)

            pqfccp = comparison.FieldComparatorPosQGram(threshold = t,
                                          q = q,
                                          max_dist = max_d,
                                          common_div = c,
                                          padded = True,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorPosQGram',
                                          do_cache = True)

            self.doStringFieldComparisonTest(pqfc)
            self.doStringFieldComparisonTest(pqfcp)
            self.doStringFieldComparisonTest(pqfcc)
            self.doStringFieldComparisonTest(pqfccp)

      for c in ['average','shortest','longest']:

        # Gram class [[0]] should be the same as q=2 grams
        #
        for gc in [[[0]],
                   [[1]],
                   [[2]],
                   [[0],[1]],
                   [[0,1,2]],
                   [[0],[0,1],[1,2]]]:

          sfc = comparison.FieldComparatorSGram(threshold = t,
                                          gram_class = gc,
                                          common_div = c,
                                          padded = False,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorSGram')

          sfcp = comparison.FieldComparatorSGram(threshold = t,
                                          gram_class = gc,
                                          common_div = c,
                                          padded = True,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorSGram')

          sfcc = comparison.FieldComparatorSGram(threshold = t,
                                          gram_class = gc,
                                          common_div = c,
                                          padded = False,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorSGram',
                                          do_cache = True)

          sfccp = comparison.FieldComparatorSGram(threshold = t,
                                          gram_class = gc,
                                          common_div = c,
                                          padded = True,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorSGram',
                                          do_cache = True)

          self.doStringFieldComparisonTest(sfc)
          self.doStringFieldComparisonTest(sfcp)
          self.doStringFieldComparisonTest(sfcc)
          self.doStringFieldComparisonTest(sfccp)

      edfc = comparison.FieldComparatorEditDist(threshold = t,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorEditDist')
      edfcc = comparison.FieldComparatorEditDist(threshold = t,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorEditDist',
                                          do_cache = True)
      self.doStringFieldComparisonTest(edfc)
      self.doStringFieldComparisonTest(edfcc)

      dldfc = comparison.FieldComparatorDaLeDist(threshold = t,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorDaLeDist')

      dldfcc = comparison.FieldComparatorDaLeDist(threshold = t,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorDaLeDist',
                                          do_cache = True)
      self.doStringFieldComparisonTest(dldfc)
      self.doStringFieldComparisonTest(dldfcc)

      bdfc = comparison.FieldComparatorBagDist(threshold = t,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorBagDist')

      bdfcc = comparison.FieldComparatorBagDist(threshold = t,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorBagDist',
                                          do_cache = True)
      self.doStringFieldComparisonTest(bdfc)
      self.doStringFieldComparisonTest(bdfcc)

      zcfc = comparison.FieldComparatorCompress(threshold = t,
                                          compr = 'zlib',
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorCompress')
      zcfcc = comparison.FieldComparatorCompress(threshold = t,
                                          compr = 'zlib',
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorCompress',
                                          do_cache = True)
      self.doStringFieldComparisonTest(zcfc)
      self.doStringFieldComparisonTest(zcfcc)

      bcfc = comparison.FieldComparatorCompress(threshold = t,
                                          compr = 'bz2',
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorCompress')
      bcfcc = comparison.FieldComparatorCompress(threshold = t,
                                          compr = 'bz2',
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorCompress',
                                          do_cache = True)
      self.doStringFieldComparisonTest(bcfc)
      self.doStringFieldComparisonTest(bcfcc)

      for c in ['average','shortest','longest']:

        swdfc = comparison.FieldComparatorSWDist(threshold = t,
                                          common_d = c,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorSWDist')

        swdfcc = comparison.FieldComparatorSWDist(threshold = t,
                                          common_d = c,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorSWDist',
                                          do_cache = True)
        self.doStringFieldComparisonTest(swdfc)
        self.doStringFieldComparisonTest(swdfcc)

        sadfcp = comparison.FieldComparatorSyllAlDist(threshold = t,
                                          common_d = c,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorSyllAlDist')

        sadfccp = comparison.FieldComparatorSyllAlDist(threshold = t,
                                          common_d = c,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorSyllAlDist',
                                          do_cache = True)

        sadfc = comparison.FieldComparatorSyllAlDist(threshold = t,
                                          common_d = c,
                                          do_ph = False,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorSyllAlDist')

        sadfcc = comparison.FieldComparatorSyllAlDist(threshold = t,
                                          common_d = c,
                                          missing_v = self.missing_values_list,
                                          do_ph = False,
                                          desc = 'FieldComparatorSyllAlDist',
                                          do_cache = True)

        self.doStringFieldComparisonTest(sadfcp)
        self.doStringFieldComparisonTest(sadfccp)
        self.doStringFieldComparisonTest(sadfc)
        self.doStringFieldComparisonTest(sadfcc)

        for l in [1,2,3,4]:

          lcsfc = comparison.FieldComparatorLCS(threshold = t,
                                          missing_v = self.missing_values_list,
                                          common_d = c,
                                          min_co = l,
                                          desc = 'FieldComparatorLCS')

          lcsfcc = comparison.FieldComparatorLCS(threshold = t,
                                          missing_v = self.missing_values_list,
                                          common_d = c,
                                          min_c = l,
                                          desc = 'FieldComparatorLCS',
                                          do_cache = True)
          self.doStringFieldComparisonTest(lcsfc)
          self.doStringFieldComparisonTest(lcsfcc)

          ontolcsfc = comparison.FieldComparatorOntoLCS(threshold = t,
                                          missing_v = self.missing_values_list,
                                          common_d = c,
                                          min_co = l,
                                          desc = 'FieldComparatorOntoLCS')

          ontolcsfcc = comparison.FieldComparatorOntoLCS(threshold = t,
                                          missing_v = self.missing_values_list,
                                          common_d = c,
                                          min_c = l,
                                          desc = 'FieldComparatorOntoLCS',
                                          do_cache = True)
          self.doStringFieldComparisonTest(lcsfc)
          self.doStringFieldComparisonTest(lcsfcc)

      smdfc = comparison.FieldComparatorSeqMatch(threshold = t,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorSeqMatch')
      smdfcc = comparison.FieldComparatorSeqMatch(threshold = t,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorSeqMatch',
                                          do_cache = True)
      self.doStringFieldComparisonTest(smdfc)
      self.doStringFieldComparisonTest(smdfcc)

      edfc = comparison.FieldComparatorEditex(threshold = t,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorEditex')
      edfcc = comparison.FieldComparatorEditex(threshold = t,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorEditex',
                                          do_cache = True)
      self.doStringFieldComparisonTest(edfc)
      self.doStringFieldComparisonTest(edfcc)


  def testNumericComparison(self):  # - - - - - - - - - - - - - - - - - - - - -

    for p in [0, 1, 2, 5, 10, 20, 50]:

      pfc = comparison.FieldComparatorNumericPerc(max_p = p,
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorNumericPerc')
      pfcc = comparison.FieldComparatorNumericPerc(max_p = p,
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorNumericPerc',
                                        do_cache = True)
      self.doNumericFieldComparisonTest(pfc)
      self.doNumericFieldComparisonTest(pfcc)  # Use cache

    for a in [0, 1, 2, 5, 10, 20, 50, 100, 1000]:

      afc = comparison.FieldComparatorNumericAbs(max_a = a,
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorNumericAbs')
      afcc = comparison.FieldComparatorNumericAbs(max_a = a,
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorNumericAbs',
                                        do_cache = True)
      self.doNumericFieldComparisonTest(afc)
      self.doNumericFieldComparisonTest(afcc)  # Use cache

  def testDistanceComparison(self):  # - - - - - - - - - - - - - - - - - - - -

    for d in [0.1, 1.0, 1, 2, 10, 11.11, 50]:

      dfc = comparison.FieldComparatorDistance(max_d = d,
                                        geocode = self.geo_lookup_table,
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorDistance')
      dfcc = comparison.FieldComparatorDistance(max_d = d,
                                        geocode = self.geo_lookup_table,
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorDistance',
                                        do_cache = True)
      self.doDistanceFieldComparisonTest(dfc)
      self.doDistanceFieldComparisonTest(dfcc)  # Use cache

  def testDateComparison(self):  # - - - - - - - - - - - - - - - - - - - - - -

    for dp in [(3,3),(5,5),(7,9),(7,5),(1,3),(3,1)]:

      for d_format in ['ddmmyyyy','mmddyyyy','ddmmyy','mmddyy']:

        dfc = comparison.FieldComparatorDate(max_day1=dp[0], max_day2 = dp[1],
                                        date_format = d_format,
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorDate')
        dfcc = comparison.FieldComparatorDate(max_day1=dp[0], max_day2 = dp[1],
                                        date_format = d_format,
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorDate',
                                        do_cache = True)
        self.doDateFieldComparisonTest(dfc)
        self.doDateFieldComparisonTest(dfcc)  # Use cache

  def testTimeComparison(self):  # - - - - - - - - - - - - - - - - - - - - - -

    for tp in [(1,1),(1,3),(3,1),(3,3),(15,15),(7,9),(17,5),(10,23),(23,10)]:

      tfc = comparison.FieldComparatorTime(max_time1=tp[0], max_time2 = tp[1],
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorTime')
      tfcc = comparison.FieldComparatorTime(max_time1=tp[0], max_time2 = tp[1],
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorTime',
                                        do_cache = True)
      tfcd = comparison.FieldComparatorTime(max_time1=tp[0], max_time2 = tp[1],
                                        day_start = '0030',
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorTime')
      tfccd = comparison.FieldComparatorTime(max_time1=tp[0], max_time2 =tp[1],
                                        day_start = '20:03',
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorTime',
                                        do_cache = True)
      self.doTimeFieldComparisonTest(tfc)
      self.doTimeFieldComparisonTest(tfcc)  # Use cache
      self.doTimeFieldComparisonTest(tfcd)
      self.doTimeFieldComparisonTest(tfccd)  # Use cache

  def testAgeComparison(self):  # - - - - - - - - - - - - - - - - - - - - - - -

    for p in [0, 1, 2, 5, 10, 20, 50]:

      afc = comparison.FieldComparatorAge(max_p = p,
                                        date_format = 'ddmmyyyy',
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorAge')
      afcc = comparison.FieldComparatorAge(max_p = p,
                                        date_format = 'ddmmyyyy',
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorAge',
                                        do_cache = True)
      afct = comparison.FieldComparatorAge(max_p = p,
                                        date_format = 'ddmmyyyy',
                                        fix_d = 'today',
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorAge')
      afcct = comparison.FieldComparatorAge(max_p = p,
                                        date_format = 'ddmmyyyy',
                                        fix_d = 'today',
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorAge',
                                        do_cache = True)
      afco = comparison.FieldComparatorAge(max_p = p,
                                        date_format = 'ddmmyyyy',
                                        fix_d = (10,11,2006),
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorAge')
      afcco = comparison.FieldComparatorAge(max_p = p,
                                        date_format = 'ddmmyyyy',
                                        fix_d = (11,11,2006),
                                        missing_v = self.missing_values_list,
                                        desc = 'FieldComparatorAge',
                                        do_cache = True)
      self.doAgeFieldComparisonTest(afc)
      self.doAgeFieldComparisonTest(afcc)  # Use cache
      self.doAgeFieldComparisonTest(afct)
      self.doAgeFieldComparisonTest(afcct)  # Use cache
      self.doAgeFieldComparisonTest(afco)
      self.doAgeFieldComparisonTest(afcco)  # Use cache

  # ---------------------------------------------------------------------------
  # Test frequency tables

  def testFreqTables(self):

    self.freq_table = {'miller':10,'smith':20,'johns':2,'west':5,
                       'dijkstra':2,'meyer':25,'smyth':16,'john':4,'meier':22}

    for t in [0.1, 0.25, 0.5, 0.6666, 0.8, 0.9]:

      jfc = comparison.FieldComparatorJaro(threshold = t,
                                          val_freq_table = self.freq_table,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorJaro')
      wfcc = comparison.FieldComparatorWinkler(threshold = t,
                                          val_freq_table = self.freq_table,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorWinkler',
                                          do_cache = True)

      self.doStringFieldComparisonTest(jfc)
      self.doStringFieldComparisonTest(wfcc)

      for (mw, daw, aw, faw) in ((0.0,-20.0,20.0,20.0),(0.0,-10.0,10.0,10.0),
                                 (0.0,-20.0,20.0,15.0),(0.0,-10.0,10.0,15.0)):

        jfc.set_weights(missing_w=mw, agr=aw,disag=daw,freq_m=faw)
        wfcc.set_weights(missing_w=mw, agr=aw,disag=daw,freq_m=faw)

        wfcc.cache = {}  # Clear old cache entries

        for val in self.freq_table:

          fajw = jfc.compare(val,val)
          faww = wfcc.compare(val,val)

          assert fajw > daw, \
                 'Jaro frequency agreement weight not larger than ' + \
                 'disagreement weight for value %s' % (val)
          assert faww > daw, \
                 'Winkler frequency agreement weight not larger than ' + \
                 'disagreement weight for value %s' % (val)
          assert fajw < aw, \
                 'Jaro frequency agreement weight not smaller than general' + \
                 ' agreement weight for value %s: %f / %f' % (val, fajw, aw)
          assert faww < aw, \
                 'Winkler frequency agreement weight not smaller than ' + \
                 'general agreement weight for value: %f / %f' % (val,faww,aw)
          assert fajw <= faww, \
                 'Jaro frequency agreement weight is larger than Winkler ' + \
                 'frequency agreement weight for value: %s' % (val) + \
                 ' (%f / %f)' %(fajw, faww)

        # First value of each pair is in frequency table, second values are not
        #
        for (val1,val2) in [('miller','miler'),('john','johnny'),
                            ('meyer','meier'),('dijkstra','dijkstras'),
                            ('smith','smyth')]:

          fajw1 = jfc.compare(val2,val2)
          faww1 = wfcc.compare(val2,val2)
          fajw2 = jfc.compare(val1,val2)
          faww2 = wfcc.compare(val1,val2)

          assert fajw1 >= daw, \
                 'Jaro frequency agreement weight not larger than ' + \
                 'disagreement weight for value %s' % (val1)
          assert faww1 >= daw, \
                 'Winkler frequency agreement weight not larger than ' + \
                 'disagreement weight for value %s' % (val1)
          assert fajw2 >= daw, \
                 'Jaro frequency agreement weight not larger than ' + \
                 'disagreement weight for values %s / %s' % (val1,val2)
          assert faww2 >= daw, \
                 'Winkler frequency agreement weight not larger than ' + \
                 'disagreement weight for value %s / %s' % (val1,val2) + \
                 'thres=%f, agr=%f, disag=%f, winkl w=%f' % \
                 (t,aw,daw,faww2)

          assert fajw1 > fajw2, \
                 'Jaro wrong frequency weight calculations for values: ' + \
                 '%s / %s' % (val1,val2)
          assert faww1 > faww2, \
                 'Winkler wrong frequency weight calculations for values: ' + \
                 '%s / %s' % (val1,val2)

  # ---------------------------------------------------------------------------
  # Test caching

  def testCaching(self):

    # Create very simple small test data
    #
    rec_values_list_easy1 = ['pete','99','98','88','98','88','peter']
    rec_values_list_easy2 = rec_values_list_easy1[:]
    rec_values_list_easy2.reverse()

    # Create a list of 60 integers, each value 6 times
    #
    int_list_short = range(1000,1010)  # List of 10 integers
    rec_values_list_short1 = int_list_short[:] * 6
    rec_values_list_short2 = int_list_short[:] * 6
    rec_values_list_short2.reverse()

    # Create a list of 500 integers, each value 10 times
    #
    int_list_long =  range(2050,2100)  # List of 50 integers
    rec_values_list_long1 = int_list_long[:] * 10
    rec_values_list_long1.reverse()
    rec_values_list_long2 = int_list_long[:] * 10

    rec_lists = [(rec_values_list_easy1, rec_values_list_easy2),
                 (rec_values_list_short1, rec_values_list_short2),
                 (rec_values_list_long1,  rec_values_list_long2)]

    for (l1, l2) in rec_lists:  # Loop over pairs of lists - - - - - - - - - -

      # Test with disabled caching
      #
      jfc = comparison.FieldComparatorJaro(threshold = 0.5,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorJaro')

      assert jfc.do_caching == False, \
             'Wrong value for caching for "%s" (should be %s but is %s)' % \
             (jfc.description, jfc.do_caching, False)

      jfc.get_cache_stats()

      for i in l1:
        for j in l2:

          w = jfc.compare(str(i),str(j))

          assert jfc.do_caching == False, \
             'Wrong value for caching for "%s" (should be %s but is %s)' % \
             (jfc.description, jfc.do_caching, False)

          assert jfc.max_cache_size == None, \
             'Wrong maximum cache size for "%s" (should be %s but is None)' % \
             (jfc.description, jfc.max_cache_size)

          assert jfc.cache == {}, \
             'Cache is not empty for "%s" (Length should be 0 but is %d)' % \
             (jfc.description, len(jfc.cache))

          for i in jfc.cache_warn_dict_counts:
            assert jfc.cache_warn_dict_counts[i] == 0, \
            'Cache count %d list is not empty for "%s" (Length should be ' % \
             (i, jfc.description) + '0 but is %d)' % \
             (len(jfc.fc.cache_warn_dict_counts[i]))

      jfc.get_cache_stats()

      # Test with enabled caching (no maximum size)
      #
      jfc = comparison.FieldComparatorJaro(threshold = 0.5,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorJaro',
                                          do_cache = True)

      assert jfc.do_caching == True, \
             'Wrong value for caching for "%s" (should be %s but is %s)' % \
             (jfc.description, jfc.do_caching, True)

      jfc.get_cache_stats()

      for i in l1:
        for j in l2:

          w = jfc.compare(str(i), str(j))

          assert jfc.do_caching == True, \
             'Wrong value for caching for "%s" (should be %s but is %s)' % \
             (jfc.description, jfc.do_caching, True)

          assert jfc.max_cache_size == None, \
             'Wrong maximum cache size for "%s" (should be %s but is None)' % \
             (jfc.description, jfc.max_cache_size)

          assert jfc.cache != {}, \
             'Cache is empty for "%s" (Length should be %d but is 0)' % \
             (jfc.description, len(jfc.cache))

      jfc.get_cache_stats()

      # Unlimited cache size does not update cache count lists
      #
      for i in jfc.cache_warn_dict_counts:
        assert jfc.cache_warn_dict_counts[i] == 0, \
        'Cache count %d list is not empty for "%s" (Length should be ' % \
         (i, jfc.description) + '0 but is %d)' % \
         (len(jfc.fc.cache_warn_dict_counts[i]))

      jfc.get_cache_stats()

      # Test with enabled caching (with maximum size 10)
      #
      jfc = comparison.FieldComparatorJaro(threshold = 0.5,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorJaro',
                                          do_cache = True,
                                          max_cache_size = 10)

      assert jfc.do_caching == True, \
             'Wrong value for caching for "%s" (should be %s but is %s)' % \
             (jfc.description, jfc.do_caching, True)

      jfc.get_cache_stats()

      for i in l1:
        for j in l2:

          w = jfc.compare(str(i), str(j))

          assert jfc.do_caching == True, \
             'Wrong value for caching for "%s" (should be %s but is %s)' % \
             (jfc.description, jfc.do_caching, True)

          assert jfc.max_cache_size == 10, \
             'Wrong maximum cache size for "%s" (should be %s but is None)' % \
             (jfc.description, jfc.max_cache_size)

          assert jfc.cache != {}, \
             'Cache is empty for "%s" (Length should be %d but is 0)' % \
             (jfc.description, len(jfc.cache))

          assert len(jfc.cache) <= 10, \
             'Cache is too large for "%s" (Length should be <= 10 but is' % \
             (jfc.description)+' %d)' % (len(jfc.cache))

      assert jfc.cache_warn_dict_counts[2] != 0, \
             'Cache count 2 list is empty for "%s"' % (jfc.description)

      jfc.get_cache_stats()

      # Test with enabled caching (with maximum size 501)
      #
      jfc = comparison.FieldComparatorJaro(threshold = 0.5,
                                          missing_v = self.missing_values_list,
                                          desc = 'FieldComparatorJaro',
                                          do_cache = True,
                                          max_cache_size = 1000)

      assert jfc.do_caching == True, \
             'Wrong value for caching for "%s" (should be %s but is %s)' % \
             (jfc.description, jfc.do_caching, True)

      jfc.get_cache_stats()

      for i in l1:
        for j in l2:

          w = jfc.compare(str(i), str(j))

          assert jfc.do_caching == True, \
             'Wrong value for caching for "%s" (should be %s but is %s)' % \
             (jfc.description, jfc.do_caching, True)

          assert jfc.max_cache_size == 1000, \
             'Wrong maximum cache size for "%s" (should be %s but is None)' % \
             (jfc.description, jfc.max_cache_size)

          assert jfc.cache != {}, \
             'Cache is empty for "%s" (Length should be %d but is 0)' % \
             (jfc.description, len(jfc.cache))

          assert len(jfc.cache) <= 1000, \
             'Cache is too large for "%s" (Length should be <= 10 but is' % \
             (jfc.description)+' %d)' % (len(jfc.cache))

      assert jfc.cache_warn_dict_counts[2] != 0, \
             'Cache count 2 list is empty for "%s"' % (jfc.description)

      jfc.get_cache_stats()

  # ---------------------------------------------------------------------------
  # Test record comparator
  #
  def testRecordComparator(self):

    gn_jfc = comparison.FieldComparatorJaro(threshold = 0.6,
                                          missing_v = self.missing_values_list,
                                          desc = 'Givenname Jaro')
    sn_wfcc = comparison.FieldComparatorWinkler(threshold = 0.5,
                                         missing_v = self.missing_values_list,
                                         desc = 'Surname Winkler',
                                         do_cache=True)
    ln_pqfcc = comparison.FieldComparatorPosQGram(threshold = 0.6,
                                          q = 2,
                                          max_dist = 3,
                                          common_div = 'average',
                                          padded = True,
                                          missing_v = self.missing_values_list,
                                          desc = 'Locality PosQGram',
                                          do_cache=True)
    pc_kfc = comparison.FieldComparatorKeyDiff(max_key_di = 2,
                                          missing_v = self.missing_values_list,
                                          desc = 'Postcode KeyDiff')

    # Do comparisons with different weights - - - - - - - - - - - - - - - - - -
    #
    for (mw, daw, aw) in self.weight_values:

      sn_wfcc.cache =  {}  # Clear old cache entries
      ln_pqfcc.cache = {}

      gn_jfc.set_weights(missing_w = mw, agree_w = aw, disagree_w = daw)
      sn_wfcc.set_weights(missing_w = mw, agree_w = aw, disagree_w = daw)
      ln_pqfcc.set_weights(missing_w = mw, agree_w = aw, disagree_w = daw)
      pc_kfc.set_weights(missing_w = mw, agree_w = aw, disagree_w = daw)


      field_comp_list = [(gn_jfc,   'gname', 'given_name'),
                         (sn_wfcc,  'surname', 'sname'),
                         (ln_pqfcc, 'suburb', 'locality'),
                         (pc_kfc,   'postcode', 'zipcode')]

      rc = comparison.RecordComparator(self.test_data_set1,self.test_data_set2,
                                     field_comp_list, 'Test record comparator')
      rc.get_cache_stats()

      # Do some comparisons with example records - - - - - - - - - - - - - - -
      #
      for r1 in self.recs1:
        for r2 in self.recs2:

          w_vec = rc.compare(r1,r2)

          assert isinstance(w_vec, list), \
                 'Weight vector returned from record comparator is not a ' + \
                 'list: %s' % (str(w_vec))

          assert len(w_vec) == len(field_comp_list), \
                 'Weight vector returned has wronglength (should be %d but' % \
                 (len(w_vec))+' is %d: %s' % (len(field_comp_list), str(w_vec))

          for w in w_vec:
            assert isinstance(w,float), \
                   'Weight returned is not a floating-point number: %s ' % \
                   (str(w))+' (in weight vector: %s)' % (str(w_vec))
            assert w <= aw, \
                   'Weight returned is larger than agreement weight (%f): %f' \
                   % (aw, w)+' (in weight vector: %s)' % (str(w_vec))
            assert w >= daw, \
                   'Weight returned is smaller than dis-agreement weight ' + \
                   '(%f): %f' % (aw, w)+' (in weight vector: %s)' % \
                   (str(w_vec))

      rc.get_cache_stats()

# =============================================================================
# Start tests when called from command line

if (__name__ == "__main__"):

  # Intialise a logger, set level to info
  #
  my_logger = logging.getLogger()  # New logger at root level
  my_logger.setLevel(log_level)

  unittest.main()  # Run all test

# =============================================================================
