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
# The Original Software is: "standardisationTest.py"
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

"""Test module for standardisation.py.
"""

# =============================================================================
# Import necessary modules (Python standard modules first, then Febrl modules)

import os
import sets
import sys
import unittest
sys.path.append('..')

import dataset
import lookup
import standardisation

# Set the logging level to warnings - - - - - - - - - - - - - - - - - - - - - -
#
import logging
my_logger = logging.getLogger()  # New logger at root level
my_logger.setLevel(logging.WARNING)

# =============================================================================

class TestCase(unittest.TestCase):

  # Initialise test case  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def setUp(self):

    self.in_ds = dataset.DataSetCSV(descri='A standardisation test data set',
                                    access_mode='read',
                                    rec_ident = 'rec_id',
                                    field_list = [],
                                    header_line=True,
                                    write_header=True,
                                    file_name='test-standard-dataset.csv')

    self.out_ds = dataset.DataSetCSV(descrip='The standardised test data set',
                                    access_mode='write',
                                    rec_ident = 'rec_id',
                                    field_list = [('day',0),
                                                  ('month',1),
                                                  ('year',2),
                                                  ('country_code',3),
                                                  ('country_name',4),
                                                  ('area_code',5),
                                                  ('number',6),
                                                  ('extension',7),
                                                  ('title',8),
                                                  ('gender_guess',9),
                                                  ('given_name',10),
                                                  ('alt_given_name',11),
                                                  ('surname',12),
                                                  ('alt_surname',13),
                                                  ('out_pass1',14),
                                                  ('out_pass2',15)],
                                    header_line=True,
                                    write_header=True,
                                    file_name='test-standardised-dataset.csv')


    self.dates = [['Sep 1, 68',        ['1', '9', '1968']],
                  ['18 Jan 2002',      ['18','1', '2002']],
                  ['17:2:2002',        ['17','2', '2002']],
                  ['2002-02-25',       ['25','2', '2002']],
                  ['18,03,2001',       ['18','3', '2001']],
                  ['21.12.1999',       ['21','12','1999']],
                  ['February 18,19',   ['18','2', '1919']],
                  ['23\\July\\1968',   ['23','7', '1968']],
                  ['18-02-2002',       ['18','2', '2002']],
                  ['5/03/01',          ['5', '3', '2001']],
                  ['19680429',         ['29','4', '1968']],
                  ['600810',           ['10','8', '1960']],
                  ['3:05:2000',        ['3', '5', '2000']],
                  ['30.11.1989',       ['30','11','1989']],
                  ["1. January '70",   ['1', '1', '1970']],
                  ['01011970',         ['1', '1', '1970']],
                  ['10011970',         ['10','1', '1970']],
                  ['31 dec  1969',     ['31','12','1969']],
                  ['30 december  69',  ['30','12','1969']],
                  ['01011970',         ['1', '1', '1970']],
                  ['13 Feb 1945',      ['13','2', '1945']],
                  ['Feb 13, \'45',     ['13','2', '1945']],
                  ['April 29 1968',    ['29','4', '1968']],
                  ['29-4=68',          ['29','4', '1968']],
                  ['11-01-1972',       ['11','1', '1972']],
                  ['January 10. 1972', ['10','1', '1972']],
                  ['29 Feb 1932',      ['29','2', '1932']],
                  ['29 Feb 32',        ['29','2', '1932']],
                  ['11 Jun 1902',      ['11','6', '1902']],
                  ['11 Jul 1989',      ['11','7', '1989']],
                  ['12111968',         ['12','11','1968']],
                  ['      21111969  ', ['21','11','1969']]]

    self.date_parse_formats = ['%d %m %Y',   # 24 04 2002  or  24 4 2002
                               '%d %B %Y',   # 24 April 2002
                               '%d %b %Y',   # 24 Apr 2002
                               '%m %d %Y',   # 04 24 2002  or  4 24 2002
                               '%B %d %Y',   # April 24 2002
                               '%b %d %Y',   # Apr 24 2002
                               '%Y %m %d',   # 2002 04 24  or  2002 4 24
                               '%Y %B %d',   # 2002 April 24
                               '%Y %b %d',   # 2002 Apr 24
                               '%d %m %y',   # 24 04 02    or  24 4 02
                               '%d %B %y',   # 24 April 02
                               '%d %b %y',   # 24 Apr 02
                               '%y %m %d',   # 02 04 24    or  02 4 24
                               '%y %B %d',   # 02 April 24
                               '%y %b %d',   # 02 Apr 24
                               '%m %d %y',   # 04 24 02    or  4 24 02
                               '%B %d %y',   # April 24 02
                               '%b %d %y']   # Apr 24 02

    self.phonenums = \
      [('++61 2 6125 5690',        ['61', 'Australia', '02', '6125-5690', '']),
       ('0061 02 6125 5690',       ['61', 'Australia', '02', '6125-5690', '']),
       ('0061   02    6125-5690',  ['61', 'Australia', '02', '6125-5690', '']),
       ('41 312 17 84',            ['41', 'Switzerland', '', '312 17 84', '']),
       ('6125 0010',               ['61', 'Australia', '', '6125-0010', '']),
       ('1-800-764-0432',          ['1', 'USA/Canada', '800', '764-0432', '']),
       ('02 6125 0010',            ['61', 'Australia', '02', '6125-0010', '']),
       ('00 1 317-923 4523',       ['1', 'USA/Canada', '317', '923-4523', '']),
       ('1 317-923 4523',          ['1', 'USA/Canada', '317', '923-4523', '']),
       ('00111 41 312 17 84',      ['41', 'Switzerland', '', '312 17 84', '']),
       ('00001 41 312 17 84',      ['41', 'Switzerland', '', '312 17 84', '']),
       ('01 41 312 17 84',         ['41', 'Switzerland', '', '312 17 84', '']),
       ('1-541-754-3010',          ['1', 'USA/Canada', '541', '754-3010', '']),
       ('754-3010',                ['1', 'USA/Canada',   '', '754-3010', '']),
       ('754-3010ext 42',          ['1', 'USA/Canada',   '', '754-3010','42']),
       ('754-3010x 42',            ['1', 'USA/Canada',   '', '754-3010','42']),
       ('754-3010 ext 42',         ['1', 'USA/Canada',   '', '754-3010','42']),
       ('754-3010 ext. 42',        ['1', 'USA/Canada',   '', '754-3010','42']),
       ('754-3010 x. 42',          ['1', 'USA/Canada',   '', '754-3010','42']),
       ('754-3010 x42',            ['1', 'USA/Canada',   '', '754-3010','42']),
       ('(541) 754-3010',          ['1', 'USA/Canada', '541', '754-3010', '']),
       ('+1-541-754-3010',         ['1', 'USA/Canada', '541', '754-3010', '']),
       ('191 541 754 3010',        ['', '', '', '915417543010', '']),
       ('001-541-754-3010',        ['1', 'USA/Canada', '541', '754-3010', '']),
       ('636-48018',               ['61', 'Australia', '', '6364-8018', '']),
       ('(089) / 636-48018',       ['1', 'USA/Canada', '896', '364-8018', '']),
       ('+49-89-636-48018',        ['49', 'Germany', '', '89 636 48018', '']),
       ('19-49-89-636-48018',      ['', '', '', '9498963648018', '']),
       ('+61 (02) 6125 0101',      ['61', 'Australia', '02', '6125-0101', '']),
       ('++61 (02) 6125 0101',     ['61', 'Australia', '02', '6125-0101', '']),
       ('++61 (2) 6125 0101',      ['61', 'Australia', '02', '6125-0101', '']),
       ('11 +61 (2) 6125 0101',    ['', '', '', '161261250101', '']),
       ('0011 ++61 (2) 6125 0101', ['61', 'Australia', '02', '6125-0101', '']),
       ('0111 ++61 (2) 6125 0101', ['61', 'Australia', '02', '6125-0101', '']),
       ('0111 61 02 6125 0101',    ['61', 'Australia', '02', '6125-0101', '']),
       ('61 (2) 6125 0101',        ['61', 'Australia', '02', '6125-0101', ''])]

    # Names with given names first
    #
    self.names_gnames = \
      [('',                        ['','','','','','']),
      ('Peter Christen',           ['male', '','peter', '','christen', '']),
      ('"DR" Peter Christen',      ['male', 'dr', 'peter', '','christen', '']),
      ('<mr> Peter Christen',      ['male', 'mr', 'peter', '','christen', '']),
      ('{ Dr > Peter Christen',    ['male', 'dr', 'peter', '','christen', '']),
      (' " Dr Peter Christen',     ['male', 'dr', 'peter', '','christen', '']),
      ('Peter () Christen',        ['male', 'dr', 'peter', '','christen', '']),
      ('Peter Christen(DR]]',      ['male', 'dr', 'peter', '','christen', '']),
      ('Peter Christen (mister',   ['male', 'mr', 'peter', '','christen', '']),
      ('Peter Christen " mr',      ['male', 'mr', 'peter', '','christen', '']),
      ('Peter Christen {mr } ',    ['male', 'mr', 'peter', '','christen', '']),
      ('Peter Christen "dr"',      ['male', 'dr', 'peter', '','christen', '']),
      (' ( ) Peter Christen',      ['male', '','peter', '','christen', '']),
      ('Peter " " Christen',       ['male', '','peter', '','christen', '']),
      ('Peter (> Christen',        ['male', '','peter', '','christen', '']),
      (',Peter Christen--',        ['male', '','peter', '','christen', '']),
      ('-,- Peter Christen-,-',    ['male', '','peter', '','christen', '']),
      (' //  Peter Christen//',    ['male', '','peter', '','christen', '']),
      ('(Peter,Christen  )  ',     ['male', '','peter', '','christen', '']),
      ('[Peter   Christen]',       ['male', '','peter', '','christen', '']),
      ('<<Peter ,  Christen>>',    ['male', '','peter', '','christen', '']),
      ('{  Peter Christen }',      ['male', '','peter', '','christen', '']),
      ('"Peter Christen"',         ['male', '','peter', '','christen', '']),
      ("''Peter ; Christen''",     ['male', '','peter', '','christen', '']),
      ("'|Peter ?: Christen'|",    ['male', '','peter', '','christen', '']),
      ('Mr peter Christen',        ['male','mr','peter', '','christen', '']),
      ('Mister Peter CHRISTEN',    ['male','mr','peter', '','christen', '']),
      ('Petra~ Christen',          ['female', '','petra', '','christen', '']),
      ('Ms petra Christen',        ['female','ms','petra', '','christen', '']),
      ('Misses Petra CHRISTEN',    ['female','ms','petra', '','christen', '']),
      ('Peter Marco Jones',        ['male','','peter','mark','jones','']),
      ('peter almond',             ['male','','peter','','almond','']),
      ('almond peter',             ['male','','peter','','almond','']),
      ('Peter',                    ['male','','peter','','','']),
      ('alison de francesco',      ['','','','','','']),
      ('alison de-francesco',      ['','','','','','']),
      ('peter de la placa',        ['','','','','','']),
      ('peter marco de la placa',  ['','','','','','']),
      ('maria petra de la placa-miller', ['','','','','','']),
      ('maria petra vonder felde', ['','','','','','']),
      ('Christen',                 ['','','','','christen','']),
      ('Jane',                     ['female','','jane','','','']),
      ('miss anita',               ['female','ms','anita','','','']),
      ('mr p. christen',           ['male','mr','p','','christen','']),
      ('Peter mary jones',         ['','','peter','mary','jones','']),
      ('mr Peter mary jones',      ['male','mr','peter','mary','jones','']),
      ('mister Paul PETER jones-miller',
                                ['male','mr','paul','peter','jones','miller']),
      ('peter known as pete',      ['male','','peter','peter','','']),
      ('nee miller',               ['','','nee','','miller','']),
      ('peter de nee',             ['','','peter','','de nee','']),
      ('paul saint nee',           ['','','paul','','saint nee','']),
      ('saint paul nee',           ['','','saint paul','','nee','']),
      ('paula miller (nee jones)', ['','','','','','']),
      ('peter, son of nee miller', ['','','','','','']),
      ('peter (known as  pete) christen',  ['','','','','','']),
      ('peter (known as "pete") christen', ['','','','','','']),
      ('peter christen miller',       ['','','','','','']),
      ('peter christen-miller',       ['','','','','','']),
      ('peter joe christen-miller',       ['','','','','','']),
      ("peter 'joe' christen-miller",       ['','','','','','']),
      ('"sharky" peter miller',       ['','','','','','']),
      ("'barbie' sue smith-jones",       ['','','','','','']),
      ('sue "barbie" smith meyer',       ['','','','','','']),
      ('sue known as "barbie" smith meyer',       ['','','','','','']),
      ("sue 'barbie' smith-jones",       ['','','','','','']),
      ("sue 'barbie' smith jones",       ['','','','','','']),
      ('sue baby of maria jones',       ['','','','','','']),

      ('jane co lo-schiavo',       ['','','','','','']),
      ('martina louis barber',       ['','','','','','']),
      ('lisa-anne hennessy',       ['','','','','','']),
      ('michelle southam-byrnes',       ['','','','','','']),
      ('nicole win jordan',       ['','','','','','']),
      ('caroline and clarke',       ['','','','','','']),
      ('jocelyn or buskens',       ['','','','','','']),
      ('yee fung nee cheng',       ['','','','','','']),
      ('jenny khaw nee yii',       ['','','','','','']),
      ('roslyn kay sta maria',       ['','','','','','']),
      ('shelley lee di stefano',       ['','','','','','']),
      ('li qing van huisstede',       ['','','','','','']),
      ('patricia ann van den hurk',       ['','','','','','']),
      ('kim maree nguyen su',       ['','','','','','']),
      ('adriana haile de lange',       ['','','','','','']),
      ("jodene akke op't land",       ['','','','','','']),
      ('cleo ann di blasio',       ['','','','','','']),
      ('debbie saphire st quintin',       ['','','','','','']),
      ('nehmat e el chaar',       ['','','','','','']),
      ('yan chen ping yang',       ['','','','','','']),
      ('sharon leoni van ant werpen',       ['','','','','','']),
      ('nicole maria de oliveira',       ['','','','','','']),
      ('sonia denni de arman',       ['','','','','','']),
      ('nicole dan de arman',       ['','','','','','']),
      ('johdy louise dal santo',       ['','','','','','']),
      ('tamara lou st. john-morton',       ['','','','','','']),
      ('mercy jacq john peter',       ['','','','','','']),
      ('carly evelyn de st germain',       ['','','','','','']),
      ('rachael jane van buuren',       ['','','','','','']),
      ('joanna lilli van ryswyk',       ['','','','','','']),
      ('melissa ma romijn-van stey',       ['','','','','','']),
      ('wong jing ling huang',       ['','','','','','']),
      ('julie  maree mackenzie - hun',       ['','','','','','']),
      ('joanne agnes righettli (dr)',       ['','dr','','','','']),
      ('siu har ng (hung)',       ['','','','','','']),
      ('anne-maree lawrence-franks',       ['','','','','','']),
      ('mao-yao rong-fong',       ['','','','','','']),
      ('wai-fun wheeler-smith',       ['','','','','','']),
      ('lee-anne westerbrook-sim',       ['','','','','','']),
      ('kasey-lee so-chan',       ['','','','','','']),
      ('sherri-anne hilder-penningt',       ['','','','','','']),
      ('yoon-sun ahn-wu',       ['','','','','','']),
      ('ying-xia yu-guo',       ['','','','','','']),
      ('hee-jing hyde-page',       ['','','','','','']),
      ('mary-anne chung-kwon',       ['','','','','','']),
      ('marie-reine attallah-boulos',       ['','','','','','']),
      ('tracy-lea zanco-hinds',       ['','','','','','']),
      ('tracy-maria beardow-brooks',       ['','','','','','']),
      ('el-masri sheehan-hill',       ['','','','','','']),
      ('vicki-maree cheryle-anne',       ['','','','','','']),
      ('vicki-mare sheehan-anna',       ['','','','','','']),
      ('cindy-lou mckie-bailey',       ['','','','','','']),
      ('jo-ann bakoss-parson',       ['','','','','','']),
      ('wan-ching tsui-chan',       ['','','','','','']),
      ('sue-ellen bruechert-reich',       ['','','','','','']),
      ('anna-marie vearing-brown',       ['','','','','','']),
      ("lisa-jane o'connor",       ['','','','','','']),
      ("julie-anne o'malley",       ['','','','','','']),
      ("mary-jane o'doherty",       ['','','','','','']),
      ("jose-carol o'leary",       ['','','','','','']),
      ("rose-merrie o'kane",       ['','','','','','']),
      ("ymeka-emily o'neill",       ['','','','','','']),
                        ]

## check field spill - have 2 input fields


    # Names with surnames first
    #
    self.names_snames = [('Christen Peter',
                          ['male', '','peter', '','christen', '']),
                         ('Christen, Peter',
                          ['male', '','peter', '','christen', '']),
                         ('Mr Christen Peter',
                          ['male','mr','peter', '','christen', '']),
                         ('Mister CHRISTEN, Peter',
                          ['male','mr','peter', '','christen', '']),
                         ('Christen Petra',
                          ['female', '','petra', '','christen', '']),
                         ('Ms Christen, petra',
                          ['female','ms','petra', '','christen', '']),
                         ('Misses CHRISTEN, PETRA',
                          ['female','ms','petra', '','christen', '']),
                         ('peter almond',
                          ['male','','peter','','almond','']),
                         ('almond peter',
                          ['male','','peter','','almond','']),
                         ('',
                          ['','','','','','']),
                         ('Peter',
                          ['male','','peter','','','']),
                         ('Christen',
                          ['','','','','christen','']),
                         ('Jane',
                          ['female','','jane','','','']),
                         ('miss anita',
                          ['female','ms','anita','','','']),
                         ('mr p. christen',
                          ['male','mr','p','','christen','']),
                         ('jones, Peter mary',
                          ['','','peter','mary','jones','']),
                         ('mr jones Peter mary',
                          ['male','mr','peter','mary','jones','']),
                         ('mister jones-miller, Paul PETER',
                          ['male','mr','paul','peter','jones','miller']),

                        ]

    self.name_male_titles = ['mr']
    self.name_female_titles = ['ms']

    src_data_dir = '..'+os.sep+'data'+os.sep+'lookup'+os.sep

    self.name_tag_table = lookup.TagLookupTable(descr='Name tag test table')
    self.name_tag_table.load([src_data_dir+'givenname_f.tbl',
                              src_data_dir+'givenname_m.tbl',
                              src_data_dir+'name_misc.tbl',
                              src_data_dir+'name_prefix.tbl',
                              src_data_dir+'name_title.tbl',
                              src_data_dir+'saints.tbl',
                              src_data_dir+'surname.tbl'])

    self.name_corr_list = lookup.CorrectionList(descr = 'Name corr test list')
    self.name_corr_list.load(src_data_dir+'name_corr.lst')

  # Clean up test case  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  def tearDown(self):
    pass  # Nothing to clean up

  # ---------------------------------------------------------------------------
  #
  # Start test cases

  def testDateStandardiser(self):  # - - - - - - - - - - - - - - - - - - - - -
    """Test date standardiser routines"""

    return

    ds = standardisation.DateStandardiser(descript = 'Test date standardiser',
                                          parse_form = self.date_parse_formats,
                                          input_fields = ['in_date'],
                                          output_fiel = ['day','month','year'])

    rs = standardisation.RecordStandardiser(descr = 'Test record standardiser',
                                            input_dataset = self.in_ds,
                                            output_dataset = self.out_ds,
                                            comp_stand_list = [ds],
                                            pass_fiel=[('pass1','out_pass1'),
                                                       ('pass2','out_pass2')])

    for (date_str, date_res) in self.dates:

      clean_date_str = ds.clean_component(date_str)
      test_date_res = ds.standardise(date_str, clean_date_str)

      assert date_res == test_date_res, \
             'Wrong date standardisation: %s, should be: %s' % \
             (str(test_date_res), str(date_res))

    rs.standardise()  # Use record standardiser and write output file

    # Test the content of the output data set
    #
    test_ds = dataset.DataSetCSV(description='Test standardised data set',
                                 access_mode='read',
                                 rec_ident = 'rec_id',
                                 field_list = [],
                                 header_line=True,
                                 write_header=True,
                                 file_name='test-standardised-dataset.csv')

    i = 0
    for (rec_id, rec_list) in test_ds.readall():
      test_day =   rec_list[0]
      test_month = rec_list[1]
      test_year =  rec_list[2]

      true_day =   self.dates[i][1][0]
      true_month = self.dates[i][1][1]
      true_year =  self.dates[i][1][2]

      assert test_day   == true_day, (i, rec_list[0:3], self.dates[i][1])
      assert test_month == true_month, (i, rec_list[0:3], self.dates[i][1])
      assert test_year  == true_year, (i, rec_list[0:3], self.dates[i][1])

      i += 1

  # Now another date standardiser with day and year set to None - - - - - - -
  #
  def testDateStandardiserNone(self):
    """Test date standardiser routines"""

    return

    ds = standardisation.DateStandardiser(descript = 'Test date standardiser',
                                          parse_form = self.date_parse_formats,
                                          input_fields = ['in_date'],
                                          pivot_year = 11,
                                          output_fiel = [None,'month',None])

    rs = standardisation.RecordStandardiser(descr = 'Test record standardiser',
                                            input_dataset = self.in_ds,
                                            output_dataset = self.out_ds,
                                            comp_stand_list = [ds],
                                            progress_report = 1)

    for (date_str, date_res) in self.dates:

      clean_date_str = ds.clean_component(date_str)
      test_date_res = ds.standardise(date_str, clean_date_str)

      assert date_res == test_date_res, \
             'Wrong date standardisation: %s, should be: %s' % \
             (str(test_date_res), str(date_res))

    rs.standardise()  # Use record standardiser and write output file

    # Test the content of the output data set
    #
    test_ds = dataset.DataSetCSV(description='Test standardised data set',
                                 access_mode='read',
                                 rec_ident = 'rec_id',
                                 field_list = [],
                                 header_line=True,
                                 write_header=True,
                                 file_name='test-standardised-dataset.csv')

    i = 0
    for (rec_id, rec_list) in test_ds.readall():
      test_day =   rec_list[0]
      test_month = rec_list[1]
      test_year =  rec_list[2]

      true_month = self.dates[i][1][1]

      assert test_day   == '', (i, rec_list[0:3], test_day)
      assert test_month == true_month, (i, rec_list[0:3], self.dates[i][1])
      assert test_year  == '', (i, rec_list[0:3], test_year)

      i += 1

  def testPhoneNumStandardiser(self):  # --------------------------------------
    """Test phone number standardiser routines"""

    return

    ps = standardisation.PhoneNumStandardiser(descript = \
                                              'Test phone number standardiser',
                                          input_fields = ['in_phonenum'],
                                          output_fiel = ['country_code',
                                                         'country_name',
                                                         'area_code', 'number',
                                                         'extension'])

    rs = standardisation.RecordStandardiser(descr = 'Test record standardiser',
                                            input_dataset = self.in_ds,
                                            output_dataset = self.out_ds,
                                            comp_stand_list = [ps])

    for (phonenum_str, phonenum_res) in self.phonenums:

      clean_phonenum_str = ps.clean_component(phonenum_str)
      test_phonenum_res =  ps.standardise(phonenum_str, clean_phonenum_str)

      assert phonenum_res == test_phonenum_res, \
             'Wrong phone number standardisation: %s, should be: %s' % \
             (str(test_phonenum_res), str(phonenum_res))

    rs.standardise()  # Use record standardiser and write output file

    # Test the content of the output data set
    #
    test_ds = dataset.DataSetCSV(description='Test standardised data set',
                                 access_mode='read',
                                 rec_ident = 'rec_id',
                                 field_list = [],
                                 header_line=True,
                                 write_header=True,
                                 file_name='test-standardised-dataset.csv')

    i = 0
    for (rec_id, rec_list) in test_ds.readall():
      test_country_code = rec_list[3]
      test_country_name = rec_list[4]
      test_area_code =    rec_list[5]
      test_number =       rec_list[6]
      test_extension =    rec_list[7]

      true_country_code = self.phonenums[i][1][0]
      true_country_name = self.phonenums[i][1][1]
      true_area_code =    self.phonenums[i][1][2]
      true_number =       self.phonenums[i][1][3]
      true_extension =    self.phonenums[i][1][4]

      assert test_country_code == true_country_code, \
             (i, rec_list[3:8], self.phonenums[i][1])
      assert test_country_name == true_country_name, \
             (i, rec_list[3:8], self.phonenums[i][1])
      assert test_area_code == true_area_code, \
             (i, rec_list[3:8], self.phonenums[i][1])
      assert test_number == true_number, \
             (i, rec_list[3:8], self.phonenums[i][1])
      assert test_extension == true_extension, \
             (i, rec_list[3:8], self.phonenums[i][1])

      i += 1

  # Now another phone number standardiser with two components set to None - - -
  #
  def testPhoneNumStandardiserNone(self):
    """Test phone number standardiser routines"""

    return

    ps = standardisation.PhoneNumStandardiser(descript = \
                                              'Test phone number standardiser',
                                          input_fields = ['in_phonenum'],
                                          output_fiel = ['country_code',
                                                         None,
                                                         'area_code', 'number',
                                                         None])

    rs = standardisation.RecordStandardiser(descr = 'Test record standardiser',
                                            input_dataset = self.in_ds,
                                            output_dataset = self.out_ds,
                                            comp_stand_list = [ps])

    for (phonenum_str, phonenum_res) in self.phonenums:

      clean_phonenum_str = ps.clean_component(phonenum_str)
      test_phonenum_res =  ps.standardise(phonenum_str, clean_phonenum_str)

      assert phonenum_res == test_phonenum_res, \
             'Wrong phone number standardisation: %s, should be: %s' % \
             (str(test_phonenum_res), str(phonenum_res))

    rs.standardise()  # Use record standardiser and write output file

    # Test the content of the output data set
    #
    test_ds = dataset.DataSetCSV(description='Test standardised data set',
                                 access_mode='read',
                                 rec_ident = 'rec_id',
                                 field_list = [],
                                 header_line=True,
                                 write_header=True,
                                 file_name='test-standardised-dataset.csv')

    i = 0
    for (rec_id, rec_list) in test_ds.readall():
      test_country_code = rec_list[3]
      test_country_name = rec_list[4]
      test_area_code =    rec_list[5]
      test_number =       rec_list[6]
      test_extension =    rec_list[7]

      true_country_code = self.phonenums[i][1][0]
      true_area_code =    self.phonenums[i][1][2]
      true_number =       self.phonenums[i][1][3]

      assert test_country_code == true_country_code, \
             (i, rec_list[3:8], self.phonenums[i][1])
      assert test_country_name == '', \
             (i, rec_list[3:8], self.phonenums[i][1])
      assert test_area_code == true_area_code, \
             (i, rec_list[3:8], self.phonenums[i][1])
      assert test_number == true_number, \
             (i, rec_list[3:8], self.phonenums[i][1])
      assert test_extension == '', \
             (i, rec_list[3:8], self.phonenums[i][1])

      i += 1

  def testGNameStandardiser(self):  # -----------------------------------------
    """Test name standardiser routines (given name first)"""

#    return

    ns = standardisation.NameStandardiser(descript = 'Test name standardiser',
                                          input_fields = ['in_gname'],
                                          output_fiel = ['title',
                                                         'gender_guess',
                                                         'given_name',
                                                         'alt_given_name',
                                                         'surname',
                                                         'alt_surname'],
                                          female_t = self.name_female_titles,
                                          male_t = self.name_male_titles,
                                          tag_t=self.name_tag_table,
                                          corr_l=self.name_corr_list,
                                          hmm_train_fil = 'test-hmm-train.txt')

    rs = standardisation.RecordStandardiser(descr = 'Test record standardiser',
                                            input_dataset = self.in_ds,
                                            output_dataset = self.out_ds,
                                            comp_stand_list = [ns])

    for (name_str, name_res) in self.names_gnames:

      clean_name_str = ns.clean_component(name_str)
      test_name_res =  ns.standardise(name_str, clean_name_str)

#      assert name_res == test_name_res, \
#             'Wrong given name first standardisation: %s, should be: %s' % \
#             (str(test_name_res), str(name_res))

#    rs.standardise()  # Use record standardiser and write output file

    print 'Count dict:', ns.count_dict

  # Now another name standardiser that assumes surnames before given names - -
  #
  def testSNameStandardiser(self):
    """Test name standardiser routines (surname first)"""

#    return

    ns = standardisation.NameStandardiser(descript = 'Test name standardiser',
                                          input_fields = ['in_sname'],
                                          output_fiel = ['title',
                                                         'gender_guess',
                                                         'given_name',
                                                         'alt_given_name',
                                                         'surname',
                                                         'alt_surname'],
                                          female_t = self.name_female_titles,
                                          male_t = self.name_male_titles,
                                          tag_t=self.name_tag_table,
                                          corr_l=self.name_corr_list,
                                          first_name_c = 'sname',
                                          hmm_train_fi = 'test-hmm-train.txt')

    rs = standardisation.RecordStandardiser(descr = 'Test record standardiser',
                                            input_dataset = self.in_ds,
                                            output_dataset = self.out_ds,
                                            comp_stand_list =[ns])

    for (name_str, name_res) in self.names_snames:

      clean_name_str = ns.clean_component(name_str)
      test_name_res =  ns.standardise(name_str, clean_name_str)

#      assert name_res == test_name_res, \
#             'Wrong surname first standardisation: %s, should be: %s' % \
#             (str(test_name_res), str(name_res))

    print 'Count dict:', ns.count_dict

#    rs.standardise()  # Use record standardiser and write output file


  # -------------------------

  ### Finally test all CS in one RS

# =============================================================================
# Start tests when called from command line

if (__name__ == "__main__"):
  unittest.main()  # Run all test

# =============================================================================
