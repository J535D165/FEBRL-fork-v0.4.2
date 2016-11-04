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
# The Original Software is: "evalClassification.py"
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

"""Module to evaluate the various classification techniques implemented in the
   Febrl module 'classification.py'.

   Uses Febrl-0.4.

   Peter Christen, 2008/01/30

   Defines a series of data sets, indices and comparisons, then runs them and
   evaluates their linkage quality using various classifiers.

   TODO:
   - maybe also add sampling for training to the various classifiers
"""

# =============================================================================
# Imports go here

import logging
import os
import time

#import Gnuplot

import classification
import comparison
import dataset
import encode
import indexing
import measurements
import mymath
import output

# =============================================================================
# Various settings

do_plotting =   False
do_complexity = True

n = 10  # Number of cross validation folds

do_opt_thres =      True  # Flags, set to True to do certain classifiers
do_svm =            True
do_kmeans =         True
do_ffirst =         True
do_tailor =         True
do_two_step_thres = True
do_two_step_near =  True
do_timing =         Truee

# File name for writing results to (incl. date and time)
#
result_file_name = './results/evalClassification-' + \
                   time.strftime('%Y%m%d-%H%M')+'.res'

progress_precentage = 10

num_random_select_iterations = 10

# Various possible deduplication indexing (blocking) techniques
#
index_dedup_block_method = ('block',)
#index_dedup_block_method = ('sort', 3)
#index_dedup_block_method = ('qgram', 2, False, 0.8)

weight_vect_dir = './weight-vectors/'  # Where weight vector files are stored

# =============================================================================

def get_measures(result_list):
  """Function which calculates quality measures from raw classification
     counts.
     Returns accuracy, precision, recall, and f-measure values.
  """

  tp, fn, fp, tn = result_list
  tp = float(tp)
  tn = float(tn)
  fp = float(fp)
  fn = float(fn)

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

  return acc, prec, reca, fmeas

# =============================================================================
# Define a project logger

my_logger = logging.getLogger()  # New logger at root level
my_logger.setLevel(logging.WARNING)
#my_logger.setLevel(logging.INFO)

# =============================================================================
# Open the results file
#
res_file = open(result_file_name, 'w')

# =============================================================================
# Generate a list with information for each data set:
# 1) A tuple with the two data set objects (will be the same for deduplication)
# 2) The data set index object
# 3) A list with tuples on a field comparison selection, each having:
#    a) A comment string
#    b) A list with the fields used in one experiment for the function
#       classification.extract_collapse()
# 4) The function to be used to check for true matches and non-matches
#
experiment_list = []

# =============================================================================
# Define original input data sets, indices and field comparisons for them

# -----------------------------------------------------------------------------
# The publicly available Census data set with synthetic personal names and
# addresses (taken from SecondString data repository).
#
# - The 'entity_id' attribute (2nd attribute) contains entity numbers.
# - No record identifer is available.
#
census_ds_A = dataset.DataSetCSV(description='Census data set A',
                                 access_mode='read',
                                 delimiter='\t',
                                 rec_ident='rec_id',
                                 header_line=False,
                                 field_list=[('relation',0),
                                             ('entity_id',1),
                                             ('surname',2),
                                             ('given_name',3),
                                             ('middle_inital',4),
                                             ('zipcode',5),
                                             ('suburb',6)],
                    file_name = './data/secondstring/censusTextSegmentedA.tab')

census_ds_B = dataset.DataSetCSV(description='Census data set B',
                                 access_mode='read',
                                 delimiter='\t',
                                 rec_ident='rec_id',
                                 header_line=False,
                                 field_list=[('relation',0),
                                             ('entity_id',1),
                                             ('surname',2),
                                             ('given_name',3),
                                             ('middle_inital',4),
                                             ('zipcode',5),
                                             ('suburb',6)],
                    file_name = './data/secondstring/censusTextSegmentedB.tab')

census_index_list = [[['surname', 'surname', False, False, None,
                       [encode.dmetaphone,3]]],
                     [['given_name','given_name', False, False, None,
                       [encode.dmetaphone,3]]],
                     [['surname', 'surname', False, False, 1, []],
                      ['given_name','given_name', False, False, 1, []]],
                     [['zipcode', 'zipcode', False, False, None, []]],
                     [['suburb',  'suburb',   False, False, None,
                       [encode.dmetaphone,3]]]]

# Exact comparison of 'relation' and 'entity_id': If 'relation' is different
# and 'entity_id' is the same it is a match, otherwise not
#
census_entity_id_exact = comparison.FieldComparatorExactString(desc = \
                                                             'entity_id_exact')

census_surname_winkler =    comparison.FieldComparatorWinkler(thres=0,
                                                      desc = 'surname_winkler')
census_given_name_winkler = comparison.FieldComparatorWinkler(thres=0,
                                                   desc = 'given_name_winkler')
census_suburb_winkler =     comparison.FieldComparatorWinkler(thres=0,
                                                       desc = 'suburb_winkler')

census_surname_qgram =    comparison.FieldComparatorQGram(thres=0, q=2,
                                                    common_divisor='average',
                                                        desc = 'surname_qgram')
census_given_name_qgram = comparison.FieldComparatorQGram(thres=0, q=2,
                                                    common_divisor='average',
                                                     desc = 'given_name_qgram')
census_suburb_qgram =     comparison.FieldComparatorQGram(thres=0, q=2,
                                                    common_divisor='average',
                                                    desc = 'suburb_qgram')

census_surname_bagdist =    comparison.FieldComparatorBagDist(thres=0,
                                                      desc = 'surname_bagdist')
census_given_name_bagdist = comparison.FieldComparatorBagDist(thres=0,
                                                   desc = 'given_name_bagdist')
census_suburb_bagdist =     comparison.FieldComparatorBagDist(thres=0,
                                                       desc = 'suburb_bagdist')

census_surname_lcs =    comparison.FieldComparatorLCS(thres=0,
                                                      min_common_len=2,
                                                      common_divisor='average',
                                                      desc = 'surname_lcs')
census_given_name_lcs = comparison.FieldComparatorLCS(thres=0,
                                                      min_common_len=2,
                                                      common_divisor='average',
                                                      desc = 'given_name_lcs')
census_suburb_lcs =     comparison.FieldComparatorLCS(thres=0,
                                                      min_common_len=2,
                                                      common_divisor='average',
                                                      desc = 'suburb_lcs')

census_middle_inital_exact = comparison.FieldComparatorExactString(desc = \
                                                         'middle_inital_exact')
census_zipcode_exact =       comparison.FieldComparatorExactString(desc = \
                                                               'zipcode_exact')

census_zipcode_keydiff = comparison.FieldComparatorKeyDiff(max_key_diff=2,
                                                      desc = 'zipcode_keydiff')

census_fc_list = [(census_entity_id_exact,    'entity_id',    'entity_id'),
                  (census_surname_winkler,    'surname',      'surname'),
                  (census_given_name_winkler, 'given_name',   'given_name'),
                  (census_suburb_winkler,     'suburb',       'suburb'),
                  (census_surname_qgram,      'surname',      'surname'),
                  (census_given_name_qgram,   'given_name',   'given_name'),
                  (census_suburb_qgram,       'suburb',       'suburb'),
                  (census_surname_bagdist,    'surname',      'surname'),
                  (census_given_name_bagdist, 'given_name',   'given_name'),
                  (census_suburb_bagdist,     'suburb',       'suburb'),
                  (census_surname_lcs,        'surname',      'surname'),
                  (census_given_name_lcs,     'given_name',   'given_name'),
                  (census_suburb_lcs,         'suburb',       'suburb'),

                  (census_middle_inital_exact,'middle_inital','middle_inital'),
                  (census_zipcode_exact,      'zipcode',      'zipcode'),
                  (census_zipcode_keydiff,    'zipcode',      'zipcode')]

census_rec_comp = comparison.RecordComparator(census_ds_A, census_ds_B,
                                              census_fc_list,
                                              'Census record comparator')
# List with sub-sets of the field comparisons for experiments
#
census_sel_list = [('Winkler',     [ (1,), (2,), (3,),(13,),(15,)]),
                   ('Q-Gram',      [ (4,), (5,), (6,),(13,),(15,)]),
                   ('Bag-Distance',[ (7,), (8,), (9,),(13,),(15,)]),
                   ('LCS',         [(10,),(11,),(12,),(13,),(15,)])]

cens_bigmatch_index = indexing.BigMatchIndex(descrip = 'Census BigMatch index',
                                       dataset1 = census_ds_A,
                                       dataset2 = census_ds_B,
                                       weight_vec_file = weight_vect_dir + \
                                    'census-bigmatch-index-weight-vectors.csv',
                                       rec_comparator = census_rec_comp,
                                       progress=progress_precentage,
                                       block_method = index_dedup_block_method,
                                       index_def = census_index_list)

cens_full_index = indexing.FullIndex(description = 'Census Full index',
                                     dataset1 = census_ds_A,
                                     dataset2 = census_ds_B,
                                     weight_vec_file = weight_vect_dir + \
                                     'census-full-index-weight-vectors.csv',
                                     rec_comparator = census_rec_comp,
                                     progress=progress_precentage,
                                     index_def = []) # Not needed for full ind.

# Function to be used to check for true matches and non-matches
#
def census_check_funct(rec_id1, rec_id2, weight_vec):
  return (weight_vec[0] == 1.0)

# Function to be used to extract the record identifier from a raw record
#
def census_get_id_funct(rec):
  return rec[1]

experiment_list.append(((census_ds_A,census_ds_B), cens_bigmatch_index,
                        census_sel_list, census_check_funct,
                        census_get_id_funct))
experiment_list.append(((census_ds_A,census_ds_B), cens_full_index,
                        census_sel_list, census_check_funct,
                        census_get_id_funct))

# -----------------------------------------------------------------------------
# The publicly available Restaurant data set (Fodors/ Zagats) restaurant names
# and addresses (taken from SecondString data repository).
#
# - The 'class' attribute contains entity numbers.
# - No record identifer is available.
#
rest_ds = dataset.DataSetCSV(description='Restaurant data set',
                             access_mode='read',
                             rec_ident='rec_id',
                             header_line=True,
                             file_name = './data/secondstring/restaurant.csv')

rest_index_list = [[['phone', 'phone', False, False, None,
                     [encode.get_substring,0,4]]],
                   [['type',  'type', False, False, None, [encode.soundex,4]]],
                   [['city',  'city', False, False, None,
                     [encode.dmetaphone,4]]]]

# Exact comparison of 'class:' If the same it is a match, otherwise not
#
rest_class_exact = comparison.FieldComparatorExactString(desc = 'class_exact')

rest_name_winkler = comparison.FieldComparatorWinkler(thres=0,
                                                      desc = 'name_winkler')
rest_addr_winkler = comparison.FieldComparatorWinkler(thres=0,
                                                      desc = 'addr_winkler')
rest_city_winkler = comparison.FieldComparatorWinkler(thres=0,
                                                      desc = 'city_winkler')

rest_name_qgram =   comparison.FieldComparatorQGram(thres=0, q=2,
                                                    common_divisor='average',
                                                    desc = 'name_qgram')
rest_addr_qgram =   comparison.FieldComparatorQGram(thres=0, q=2,
                                                    common_divisor='average',
                                                    desc = 'addr_qgram')
rest_city_qgram =   comparison.FieldComparatorQGram(thres=0, q=2,
                                                    common_divisor='average',
                                                    desc = 'city_qgram')

rest_name_bagdist = comparison.FieldComparatorBagDist(thres=0,
                                                      desc = 'name_bagdist')
rest_addr_bagdist = comparison.FieldComparatorBagDist(thres=0,
                                                      desc = 'addr_bagdist')
rest_city_bagdist = comparison.FieldComparatorBagDist(thres=0,
                                                      desc = 'city_bagdist')

rest_name_lcs =   comparison.FieldComparatorLCS(thres=0, min_common_len=2,
                                                common_divisor='average',
                                                desc = 'name_lcs')
rest_addr_lcs =   comparison.FieldComparatorLCS(thres=0, min_common_len=2,
                                                common_divisor='average',
                                                desc = 'addr_lcs')
rest_city_lcs =   comparison.FieldComparatorLCS(thres=0, min_common_len=2,
                                                common_divisor='average',
                                                desc = 'city_lcs')

rest_fc_list = [(rest_class_exact,  'class', 'class'),
                (rest_name_winkler, 'name',  'name'),
                (rest_addr_winkler, 'addr',  'addr'),
                (rest_city_winkler, 'city',  'city'),
                (rest_name_qgram,   'name',  'name'),
                (rest_addr_qgram,   'addr',  'addr'),
                (rest_city_qgram,   'city',  'city'),
                (rest_name_bagdist, 'name',  'name'),
                (rest_addr_bagdist, 'addr',  'addr'),
                (rest_city_bagdist, 'city',  'city'),
                (rest_name_lcs,     'name',  'name'),
                (rest_addr_lcs,     'addr',  'addr'),
                (rest_city_lcs,     'city',  'city')]

rest_rec_comp = comparison.RecordComparator(rest_ds, rest_ds, rest_fc_list,
                                            'Restaurant record comparator')

# List with sub-sets of the field comparisons for experiments
#
rest_sel_list = [('Winkler',     [ (1,), (2,), (3,)]),
                 ('Q-Gram',      [ (4,), (5,), (6,)]),
                 ('Bag-Distance',[ (7,), (8,), (9,)]),
                 ('LCS',         [(10,),(11,),(12,)])]

rest_dedup_index = indexing.DedupIndex(description = 'Restaurant Dedup index',
                                       dataset1 = rest_ds,
                                       dataset2 = rest_ds,
                                       weight_vec_file = weight_vect_dir + \
                                   'restaurant-dedup-index-weight-vectors.csv',
                                       rec_comparator = rest_rec_comp,
                                       progress=progress_precentage,
                                       block_method = index_dedup_block_method,
                                       index_def = rest_index_list)

rest_full_index = indexing.FullIndex(description = 'Restaurant Full index',
                                     dataset1 = rest_ds,
                                     dataset2 = rest_ds,
                                     weight_vec_file = weight_vect_dir + \
                                    'restaurant-full-index-weight-vectors.csv',
                                     rec_comparator = rest_rec_comp,
                                     progress=progress_precentage,
                                     index_def = []) # Not needed for full ind.

# Function to be used to check for true matches and non-matches
#
def rest_check_funct(rec_id1, rec_id2, weight_vec):
  return (weight_vec[0] == 1.0)

# Function to be used to extract the record identifier from a raw record
#
def rest_get_id_funct(rec):
  return rec[-1]

experiment_list.append(((rest_ds, rest_ds), rest_dedup_index, rest_sel_list,
                        rest_check_funct, rest_get_id_funct))
experiment_list.append(((rest_ds, rest_ds), rest_full_index, rest_sel_list,
                        rest_check_funct, rest_get_id_funct))

# -----------------------------------------------------------------------------
# The publicly available Cora containing bibliographic citations.

cora_ds = dataset.DataSetCSV(description='Cora data set',
                             access_mode='read',
                             rec_ident='rec_id',
                             delimiter='\t',
                             header_line=False,
                             field_list=[('unknown',0),
                                         ('paper_id',1),
                                         ('author_list',2),
                                         ('pub_details',3),
                                         ('title',4),
                                         ('affiliation',5),
                                         ('conf_journal',6),
                                         ('location',7),
                                         ('publisher',8),
                                         ('year',9),
                                         ('pages',10),
                                         ('editors',11),
                                         ('appear',12),
                                         ('month',13)],
                             file_name = './data/secondstring/cora.tab')

cora_index_list = [[['author_list','author_list', False, False, None,
                     [encode.dmetaphone,3]]],
                   [['title', 'title', False, False, None,
                     [encode.dmetaphone,3]]],
                   [['conf_journal','conf_journal', False, False, None,
                     [encode.dmetaphone,3]]],
                   [['year', 'year', False, False, None,
                     [encode.dmetaphone,3]]]]

# Exact comparison of 'paper_id': If the same it is a match, otherwise not
#
cora_paper_id_exact = comparison.FieldComparatorExactString(desc = \
                                                              'paper_id_exact')

cora_author_list_winkler = comparison.FieldComparatorWinkler(thres=0,
                                                 desc = 'author_list_winkler')
cora_title_winkler = comparison.FieldComparatorWinkler(thres=0,
                                                 desc = 'title_winkler')
cora_conf_journal_winkler = comparison.FieldComparatorWinkler(thres=0,
                                                 desc = 'conf_journal_winkler')

cora_author_list_qgram = comparison.FieldComparatorQGram(thres=0,
                                                   common_divisor='average',
                                                   desc = 'author_list_qgram')
cora_title_qgram = comparison.FieldComparatorQGram(thres=0,
                                                   common_divisor='average',
                                                   desc = 'title_qgram')
cora_conf_journal_qgram = comparison.FieldComparatorQGram(thres=0,
                                                   common_divisor='average',
                                                   desc = 'conf_journal_qgram')

cora_author_list_bagdist = comparison.FieldComparatorBagDist(thres=0,
                                                 desc = 'author_list_bagdist')
cora_title_bagdist = comparison.FieldComparatorBagDist(thres=0,
                                                 desc = 'title_bagdist')
cora_conf_journal_bagdist = comparison.FieldComparatorBagDist(thres=0,
                                                 desc = 'conf_journal_bagdist')

cora_author_list_lcs = comparison.FieldComparatorLCS(thres=0,
                                                   min_common_len=2,
                                                   common_divisor='average',
                                                   desc = 'author_list_qgram')
cora_title_lcs = comparison.FieldComparatorLCS(thres=0,
                                                   min_common_len=2,
                                                   common_divisor='average',
                                                   desc = 'title_qgram')
cora_conf_journal_lcs = comparison.FieldComparatorLCS(thres=0,
                                                   min_common_len=2,
                                                   common_divisor='average',
                                                   desc = 'conf_journal_qgram')

cora_year_exact = comparison.FieldComparatorExactString(desc = 'year_exact')
cora_year_keydiff = comparison.FieldComparatorKeyDiff(max_key_diff=1,
                                                      desc = 'year_keydiff')

cora_pages_exact = comparison.FieldComparatorExactString(desc = 'pages_exact')
cora_pages_keydiff = comparison.FieldComparatorKeyDiff(max_key_diff=2,
                                                      desc = 'pages_keydiff')

cora_fc_list = [(cora_paper_id_exact,       'paper_id',    'paper_id'),
                (cora_author_list_winkler,  'author_list', 'author_list'),
                (cora_title_winkler,        'title',       'title'),
                (cora_conf_journal_winkler, 'conf_journal','conf_journal'),
                (cora_author_list_qgram,    'author_list', 'author_list'),
                (cora_title_qgram,          'title',       'title'),
                (cora_conf_journal_qgram,   'conf_journal','conf_journal'),
                (cora_author_list_bagdist,  'author_list', 'author_list'),
                (cora_title_bagdist,        'title',       'title'),
                (cora_conf_journal_bagdist, 'conf_journal','conf_journal'),
                (cora_author_list_lcs,      'author_list', 'author_list'),
                (cora_title_lcs,            'title',       'title'),
                (cora_conf_journal_lcs,     'conf_journal','conf_journal'),
                (cora_year_exact,           'year',        'year'),
                (cora_year_keydiff,         'year',        'year'),
                (cora_pages_exact,          'pages',       'pages'),
                (cora_pages_keydiff,        'pages',       'pages')]

cora_rec_comp = comparison.RecordComparator(cora_ds, cora_ds, cora_fc_list,
                                            'Cora record comparator')

# List with sub-sets of the field comparisons for experiments
#
cora_sel_list = [('Winkler',     [ (1,), (2,), (3,),(14,),(16,)]),
                 ('Q-Gram',      [ (4,), (5,), (6,),(14,),(16,)]),
                 ('Bag-Distance',[ (7,), (8,), (9,),(14,),(16,)]),
                 ('LCS',         [(10,),(11,),(12,),(14,),(16,)])]

cora_dedup_index = indexing.DedupIndex(description = 'Cora Dedup index',
                                       dataset1 = cora_ds,
                                       dataset2 = cora_ds,
                                       weight_vec_file = weight_vect_dir + \
                                   'cora-dedup-index-weight-vectors.csv',
                                       rec_comparator = cora_rec_comp,
                                       progress=progress_precentage,
                                       block_method = index_dedup_block_method,
                                       index_def = cora_index_list)

# Function to be used to check for true matches and non-matches
#
def cora_check_funct(rec_id1, rec_id2, weight_vec):
  return (weight_vec[0] == 1.0)

# Function to be used to extract the record identifier from a raw record
#
def cora_get_id_funct(rec):
  return rec[1]

experiment_list.append(((cora_ds, cora_ds), cora_dedup_index, cora_sel_list,
                        cora_check_funct, cora_get_id_funct))

# -----------------------------------------------------------------------------
# Synthetic data sets of different sizes generated with Febrl data generator

# Index and field comparisons are the same for all synthetic data sets
#
synth_index_list = [[['surname', 'surname', False, False, None,
                      [encode.soundex]],
                     ['postcode', 'postcode', False, False, None,
                      [encode.get_substring,0,2]]],
                    [['given_name', 'given_name',  False, False, None,
                      [encode.soundex]],
                     ['postcode', 'postcode', False, False, None,
                      [encode.get_substring,1,3]]],
                    [['suburb', 'suburb', False, False, None,
                      [encode.soundex]],
                     ['postcode', 'postcode', False, False, None,
                      [encode.get_substring,2,4]]],
                    [['suburb', 'suburb', False, False, None,
                      [encode.soundex]],
                     ['street_number', 'street_number', False, False, None,
                      []]],
                    [['postcode', 'postcode', False, False, None, []],
                     ['age', 'age', False, False, None, []]],
                    [['address_1', 'address_1', False, False, 3, []],
                     ['age', 'age', False, False, None, []]],
                    [['surname', 'surname', False, False, 3, []],
                     ['state', 'state', False, False, None, []]],
                    [['given_name', 'given_name', False, False, 3, []],
                     ['state', 'state', False, False, None, []]],
                    [['surname', 'surname', False, False, None,
                      [encode.get_substring,0,2]],
                     ['given_name', 'given_name', False, False, None,
                      [encode.get_substring,0,2]]]]

synth_surname_winkler =    comparison.FieldComparatorWinkler(thres=0,
                                                      desc = 'surname_winkler')
synth_given_name_winkler = comparison.FieldComparatorWinkler(thres=0,
                                                   desc = 'given_name_winkler')
synth_address_1_winkler =   comparison.FieldComparatorWinkler(thres=0,
                                                    desc = 'address_1_winkler')
synth_suburb_winkler =     comparison.FieldComparatorWinkler(thres=0,
                                                       desc = 'suburb_winkler')

synth_surname_qgram =    comparison.FieldComparatorQGram(thres=0, q=2,
                                                     common_divisor='average',
                                                     desc = 'surname_qgram')
synth_given_name_qgram = comparison.FieldComparatorQGram(thres=0, q=2,
                                                     common_divisor='average',
                                                     desc = 'given_name_qgram')
synth_address_1_qgram =   comparison.FieldComparatorQGram(thres=0, q=2,
                                                    common_divisor='average',
                                                    desc = 'address_1_qgram')
synth_suburb_qgram =     comparison.FieldComparatorQGram(thres=0, q=2,
                                                     common_divisor='average',
                                                     desc = 'suburb_qgram')

synth_surname_bagdist =    comparison.FieldComparatorBagDist(thres=0,
                                                      desc = 'surname_bagdist')
synth_given_name_bagdist = comparison.FieldComparatorBagDist(thres=0,
                                                   desc = 'given_name_bagdist')
synth_address_1_bagdist =   comparison.FieldComparatorBagDist(thres=0,
                                                    desc = 'address_1_bagdist')
synth_suburb_bagdist =     comparison.FieldComparatorBagDist(thres=0,
                                                       desc = 'suburb_bagdist')

synth_surname_lcs =    comparison.FieldComparatorLCS(thres=0, min_common_len=2,
                                                     common_divisor='average',
                                                     desc = 'surname_lcs')
synth_given_name_lcs = comparison.FieldComparatorLCS(thres=0, min_common_len=2,
                                                     common_divisor='average',
                                                     desc = 'given_name_lcs')
synth_address_1_lcs =  comparison.FieldComparatorLCS(thres=0, min_common_len=2,
                                                    common_divisor='average',
                                                    desc = 'address_1_lcs')
synth_suburb_lcs =     comparison.FieldComparatorLCS(thres=0, min_common_len=2,
                                                     common_divisor='average',
                                                     desc = 'suburb_lcs')

synth_street_num_exact = comparison.FieldComparatorExactString(desc = \
                                                            'street_num_exact')
synth_street_num_keydiff = comparison.FieldComparatorKeyDiff(max_key_diff=1,
                                                   desc = 'street_num_keydiff')

synth_postcode_exact = comparison.FieldComparatorExactString(desc = \
                                                              'postcode_exact')
synth_postcode_keydiff = comparison.FieldComparatorKeyDiff(max_key_diff=1,
                                                     desc = 'postcode_keydiff')

synth_state_exact = comparison.FieldComparatorExactString(desc = \
                                                                 'state_exact')
synth_state_keydiff = comparison.FieldComparatorKeyDiff(max_key_diff=1,
                                                        desc = 'state_keydiff')

synth_fc_list = [(synth_surname_winkler,    'surname',       'surname'),
                 (synth_given_name_winkler, 'given_name',    'given_name'),
                 (synth_address_1_winkler,  'address_1',     'address_1'),
                 (synth_suburb_winkler,     'suburb',        'suburb'),
                 (synth_surname_qgram,      'surname',       'surname'),
                 (synth_given_name_qgram,   'given_name',    'given_name'),
                 (synth_address_1_qgram,    'address_1',     'address_1'),
                 (synth_suburb_qgram,       'suburb',        'suburb'),
                 (synth_surname_bagdist,    'surname',       'surname'),
                 (synth_given_name_bagdist, 'given_name',    'given_name'),
                 (synth_address_1_bagdist,  'address_1',     'address_1'),
                 (synth_suburb_bagdist,     'suburb',        'suburb'),
                 (synth_surname_lcs,        'surname',       'surname'),
                 (synth_given_name_lcs,     'given_name',    'given_name'),
                 (synth_address_1_lcs,      'address_1',     'address_1'),
                 (synth_suburb_lcs,         'suburb',        'suburb'),
                 (synth_street_num_exact,   'street_number', 'street_number'),
                 (synth_street_num_keydiff, 'street_number', 'street_number'),
                 (synth_postcode_exact,     'postcode',      'postcode'),
                 (synth_postcode_keydiff,   'postcode',      'postcode'),
                 (synth_state_exact,        'state',         'state'),
                 (synth_state_keydiff,      'state',         'state')]

# List with sub-sets of the field comparisons for experiments
#
synth_sel_list = [('Winkler',     [ (0,), (1,), (2,), (3,),(17,),(19,),(21,)]),
                  ('Q-Gram',      [ (4,), (5,), (6,), (7,),(17,),(19,),(21,)]),
                  ('Bag-Distance',[ (8,), (9,),(10,),(11,),(17,),(19,),(21,)]),
                  ('LCS',         [(12,),(13,),(14,),(15,),(17,),(19,),(21,)])]

# Function to be used to check for true matches and non-matches
#
def synth_check_funct(rec_id1, rec_id2, weight_vec):
    return (rec_id1[:-1] == rec_id2[:-1])

# Function to be used to extract the record identifier from a raw record
#
def synth_get_id_funct(rec):
  return rec[0][:-1]

# Loop over different data set sizes - - - - - - - - - - - - - - - - - - - - -
#
for size in ['B_1000',  'C_1000',  'B_2500',  'C_2500',
             'B_5000',  'C_5000',  'B_10000', 'C_10000',
             'B_25000', 'C_25000', 'B_50000', 'C_50000']:

  ds_file_name = 'dataset_%s.csv.gz' % (size)

  synth_ds = dataset.DataSetCSV(description='Febrl synthetic data set %s' % \
                                            (size),
                                access_mode='read',
                                rec_ident='rec_id',
                                header_line=True,
                                file_name = './data/dedup-dsgen/%s' % \
                                            (ds_file_name))

  synth_rec_comp = comparison.RecordComparator(synth_ds, synth_ds,
                                               synth_fc_list,
                                               'Febrl synthetic data record' \
                                               + ' comparator')

  synth_dedup_index = indexing.DedupIndex(description = 'Febrl synthetic ' + \
                                                        ' Dedup index',
                                       dataset1 = synth_ds,
                                       dataset2 = synth_ds,
                                       weight_vec_file = weight_vect_dir + \
                            'synth-%s-dedup-index-weight-vectors.csv' % (size),
                                       rec_comparator = synth_rec_comp,
                                       progress=progress_precentage,
                                       block_method = index_dedup_block_method,
                                       index_def = synth_index_list)

  experiment_list.append(((synth_ds, synth_ds), synth_dedup_index,
                          synth_sel_list, synth_check_funct,
                          synth_get_id_funct))

# =============================================================================
# Run all the experiments
#
res_dict = {}  # With data set names as keys and a list with results each

########## Select experiments #######################
#
# Census: 0,1, Restaurant: 2,3, Cora: 4,
# DS-Gen B: 5,7,...15,  DS-Gen C: 6,8,...,16
#

# Started 1 Feb:
#experiment_list =[experiment_list[0],experiment_list[2],experiment_list[4],
#                  experiment_list[6],experiment_list[8],experiment_list[10],
#                  experiment_list[12],experiment_list[14],experiment_list[16]]

# Started 15 Feb:
#experiment_list =[experiment_list[6],experiment_list[8],experiment_list[10],
#                  experiment_list[12],experiment_list[14],experiment_list[16]]

# Started 19 Feb:
#experiment_list =[experiment_list[12],experiment_list[14],experiment_list[16]]

# Started 21 Feb:
#experiment_list =[experiment_list[12]]

# Started 21 Feb (later)
#experiment_list =[experiment_list[4]]

# Started 22 Feb
#experiment_list =[experiment_list[4]]

# Started 25 Feb:
#experiment_list =[experiment_list[12]]

# Started 26 Feb:
#experiment_list =[experiment_list[12]]

# Started 27 Feb:
experiment_list =[experiment_list[4]]

#####################################################

for experiment in experiment_list:
  res_list = []  # List with all results for this data set, a list for each
                 # experiment

  data_set_a =         experiment[0][0]
  data_set_b =         experiment[0][1]
  data_set_index =     experiment[1]
  field_comp_sel =     experiment[2]
  match_check_funct =  experiment[3]
  get_id_funct =       experiment[4]

  res_list_key = (data_set_a.description, data_set_b.description)

  res_file.write(os.linesep + '='*84 + os.linesep + os.linesep)

  if (data_set_a == data_set_b):
    res_file.write('Run deduplication experiments for data set:' + os.linesep)
    res_file.write('===========================================' + os.linesep)
    res_file.write('  %s' % (data_set_a.description) + os.linesep)
    res_file.write(os.linesep)
    res_file.write('  Number of records in data set: %d' % \
                   (data_set_a.num_records) + os.linesep)

  else:  # Linkage of two data sets
    res_file.write('Run linkage experiments for data sets:' + os.linesep)
    res_file.write('======================================' + os.linesep)
    res_file.write('  A: %s' % (data_set_a.description) + os.linesep)
    res_file.write('  B: %s' % (data_set_b.description) + os.linesep)
    res_file.write(os.linesep)
    res_file.write('  Number of records in data sets: %d / %d ' % \
                   (data_set_a.num_records, data_set_b.num_records) + \
                   os.linesep)

  res_file.write(os.linesep)
  res_file.write('  Index: %s' % (data_set_index.description) + os.linesep)
  res_file.write(os.linesep)
  res_file.flush()

  # Check if weight vector file is available for this experiment - - - - - - -
  #
  weight_vec_file = data_set_index.weight_vec_file

  if (weight_vec_file != None):

    # Check if the weight vector file or it's GZipped version is available
    #
    file_avail = (os.access(weight_vec_file, os.R_OK) or \
                  os.access(weight_vec_file+'.gz', os.R_OK) or \
                  os.access(weight_vec_file+'.GZ', os.R_OK))  # Can read file
  else:
    file_avail = False

  if (file_avail == True):  # File is available, so load it - - - - - - - - - -

    res_file.write('  Load weight vector dictionary from file:' + os.linesep)
    res_file.write('    %s' % (weight_vec_file) + os.linesep)
    res_file.flush()

    [field_names_list, w_vec_dict] = \
                                   output.LoadWeightVectorFile(weight_vec_file)
    res_file.write('  Field names from weight vector file:' + os.linesep)
    res_file.write('    %s' % (str(field_names_list)) + os.linesep)
    res_file.flush()

  else:  # Have to calculate weight vectors - - - - - - - - - - - - - - - - - -

    res_file.write('  Build and compact index, then run comparison step' + \
                   os.linesep)
    data_set_index.build()
    data_set_index.compact()

    my_logger.setLevel(logging.INFO)

    w_vec_run_data = data_set_index.run()

    my_logger.setLevel(logging.WARNING)

    if (w_vec_run_data == None):  # Has been written to file - - - - - - - - -

      [field_names_list, w_vec_dict] = \
                                   output.LoadWeightVectorFile(weight_vec_file)
    else:
      [field_names_list, w_vec_dict] = w_vec_run_data

  num_w_vec = len(w_vec_dict)

  res_file.write(os.linesep)
  res_file.write('  Returned dictionary with %d weight vectors' \
                 % (num_w_vec) + os.linesep)
  res_file.flush()

  # Get the true matches and true non-matches in the weight vector dictionary
  #
  true_m_set, true_nm_set = \
                  classification.get_true_matches_nonmatches(w_vec_dict,
                                                             match_check_funct)
  num_true_m =  len(true_m_set)
  num_true_nm = len(true_nm_set)
  assert (num_true_m + num_true_nm) == num_w_vec

  res_file.write('  Number of true matches and non-matches: %d / %d' % \
                 (num_true_m, num_true_nm) + os.linesep)
  res_file.write(os.linesep)
  res_file.flush()

  # Check that true quality is 100% for all measures - - - - - - - - - - - - -
  #
  acc, prec, reca, fmeas = measurements.quality_measures(w_vec_dict,
                                                         true_m_set,
                                                         true_nm_set,
                                                         match_check_funct)
  assert acc ==   1.0, 'Accuracy of true matches and non-matches not 100%'
  assert prec ==  1.0, 'Precision of true matches and non-matches not 100%'
  assert prec ==  1.0, 'Recall of true matches and non-matches not 100%'
  assert fmeas == 1.0, 'F-measure of true matches and non-matches not 100%'

  # Get quality and complexity measures for weight vectors - - - - - - - - -
  #
  if (do_complexity == True):
    rr = measurements.reduction_ratio(w_vec_dict, data_set_a, data_set_b)
    pc = measurements.pairs_completeness(w_vec_dict, data_set_a, data_set_b,
                                         get_id_funct, match_check_funct)
    pq = measurements.pairs_quality(w_vec_dict, match_check_funct)

    res_file.write('  Reduction ratio:    %7.2f%%' % (100.0*rr) + os.linesep)
    res_file.write('  Pairs completeness: %7.2f%%' % (100.0*pc) + os.linesep)
    res_file.write('  Pairs quality:      %7.2f%%' % (100.0*pq) + os.linesep)
    res_file.write(os.linesep)

  res_file.write('-'*84 + os.linesep + os.linesep)

  # ===========================================================================
  # End of initialisation, now perform various classifications

#############
  field_comp_sel = [field_comp_sel[0]] #### ONLY DO WINKLER #################
#############

  # Loop over the field comparison selections ---------------------------------
  #
  for (sel_name, sel_list) in field_comp_sel:

    exp_res_list = [sel_name]  # All results for this experiment

    num_weights = len(sel_list)

    res_file.write('  Field comparison selection: %s' % (sel_name) + \
                   os.linesep)
    res_file.write('  ----------------------------'+'-'*len(sel_name) + \
                   os.linesep + os.linesep)

    res_file.write('    Number of fields comparisons: %d' % (num_weights) + \
                   os.linesep)

    field_names = []
    for tup in sel_list:
      field_names.append(field_names_list[tup[0]])
    res_file.write('    Selected fields:' + os.linesep)
    res_file.write('      %s' % (str(field_names)) + os.linesep + os.linesep)

    sel_w_vec_dict = classification.extract_collapse_weight_vectors(sel_list,
                                                                    w_vec_dict)

    if (do_plotting == True):  # Plot true match and non-match sets
      plot_weight_vectors(sel_w_vec_dict, true_m_set, true_m_set,
                          true_nm_set, true_nm_set, 'True match status')

    # Write header of result tables to file - - - - - - - - - - - - - - - - - -
    #
    res_file.write('    '+'='*80 + os.linesep)
    res_file.write('    Classifier experiment                |  Acc  |  ' + \
                   'Prec |  Reca | F-Meas|  Time (s)' + os.linesep)
    res_file.write('    '+'='*80 + os.linesep)

    # Now conduct the various classification experiments ----------------------

    # -------------------------------------------------------------------------
    # Classify using a one-dimensional classifier that knows true match status
    #
    res_file.write('    ' + \
                   'Optimal threshold classifier (1-dimensional)'.center(79) \
                   + os.linesep)
    res_file.write('    ' + '-'*80 + os.linesep)

    one_dim_sel_list = [tuple(range(len(sel_list)))]
    one_dim_w_vec_dict = \
              classification.extract_collapse_weight_vectors(one_dim_sel_list,
                                                              sel_w_vec_dict)
    assert len(one_dim_w_vec_dict) == num_w_vec

    # Three variations: minimise pos-neg, pos only, neg only - - - - - - - - -
    #
    opt_thres_param_list = [('Minimise false pos-neg ', 'pos-neg'),
                            ('Minimise false pos only', 'pos'),
                            ('Minimise false neg only', 'neg')]

    if (do_opt_thres == False):
      opt_thres_param_list = []  # Don't perform these experiments

    opt_thres_res_list = []

    for (class_name, min_method_name) in opt_thres_param_list:

      opt_thres_classifier = classification.OptimalThreshold(bin_width = 0.01,
                                                  min_method = min_method_name)
      start_time = time.time()
      opt_thres_classifier.train(one_dim_w_vec_dict, true_m_set, true_nm_set)
      res_sets = opt_thres_classifier.classify(one_dim_w_vec_dict)
      opt_time = time.time() - start_time

      test_results = opt_thres_classifier.cross_validate(one_dim_w_vec_dict,
                                                         true_m_set,
                                                         true_nm_set, n)

      opt_acc, opt_prec, opt_reca, opt_fmeas = get_measures(test_results)

      opt_thres_res_list.append((class_name, opt_acc, opt_prec, opt_reca,
                                 opt_fmeas, opt_time))
      exp_res_list.append(opt_thres_res_list[-1])  # Also append to all results

      (class_name, acc, prec, reca, fmeas, rtime) = opt_thres_res_list[-1]
      class_name = class_name.ljust(37)
      res_file.write('    %37s| %5.3f | %5.3f | %5.3f | %5.3f | %9.3f' % \
                     (class_name, acc, prec, reca, fmeas, rtime) + os.linesep)
      res_file.flush()

      del opt_thres_classifier

    res_file.write('    ' + '-'*80 + os.linesep)
    res_file.flush()

    # -------------------------------------------------------------------------
    # Supervised support vector machine (SVM) classifier on all weight vectors
    #
    res_file.write('    ' + 'Support vector machine classifier'.center(79) + \
                   os.linesep)
    res_file.write('    ' + '-'*80 + os.linesep)

    svm_param_list = [('Linear kernel, C=10 ', 'LINEAR',  10),
                      ('Linear kernel, C= 1 ', 'LINEAR',  1),
                      ('Linear kernel, C=0.1', 'LINEAR',  0.1),
                      ('Poly kernel, C=10 ',   'POLY',  10),
                      ('Poly kernel, C= 1 ',   'POLY',  1),
                      ('Poly kernel, C=0.1',   'POLY',  0.1),
                      ('RBF kernel, C=10 ',    'RBF',   10),
                      ('RBF kernel, C= 1 ',    'RBF',   1),
                      ('RBF kernel, C=0.1',    'RBF',   0.1)]
    if (do_svm == False):
      svm_param_list = []  # Don't perform these experiments

    svm_res_list = []

    for (class_name, kernel_name, C_val) in svm_param_list:

      svm_classifier = classification.SuppVecMachine(kernel_type = kernel_name,
                                                     C = C_val)
      start_time = time.time()
      svm_classifier.train(sel_w_vec_dict, true_m_set, true_nm_set)
      res_sets = svm_classifier.classify(sel_w_vec_dict)
      svm_time = time.time() - start_time

      test_results = svm_classifier.cross_validate(sel_w_vec_dict,
                                                   true_m_set, true_nm_set, n)

      svm_acc, svm_prec, svm_reca, svm_fmeas = get_measures(test_results)

      svm_res_list.append((class_name, svm_acc, svm_prec, svm_reca, svm_fmeas,
                           svm_time))
      exp_res_list.append(svm_res_list[-1])  # Also append to all results

      (class_name, acc, prec, reca, fmeas, rtime) = svm_res_list[-1]
      class_name = class_name.ljust(37)
      res_file.write('    %37s| %5.3f | %5.3f | %5.3f | %5.3f | %9.3f' % \
                     (class_name, acc, prec, reca, fmeas, rtime) + os.linesep)
      res_file.flush()

      del svm_classifier

    res_file.write('    ' + '-'*80 + os.linesep)
    res_file.flush()

    # -------------------------------------------------------------------------
    # K-means clustering (with different distance measures)
    #
    res_file.write('    ' + 'K-means clustering'.center(79) + os.linesep)
    res_file.write('    ' + '-'*80 + os.linesep)

    kmeans_param_list = [('Manhatten dist, min/max init', mymath.distL1,
                          'min/max'),
                         ('Euclidean dist, min/max init', mymath.distL2,
                          'min/max'),
#                         ('Cosine dist, min/max init', mymath.distCosine,
#                          'min/max'),
#                         ('L-Infinity dist, min/max init', mymath.distLInf,
#                          'min/max'),
#                         ('Canberra dist, min/max init', mymath.distCanberra,
#                          'min/max'),

                         ('Manhatten dist, random init', mymath.distL1,
                          'random'),
                         ('Euclidean dist, random init', mymath.distL2,
                          'random'),
#                         ('Cosine dist, random init', mymath.distCosine,
#                          'random')
#                         ('L-Infinity dist, random init', mymath.distLInf,
#                          'random')]
#                         ('Canberra dist, random init', mymath.distCanberra,
#                          'random')
                        ]

    if (do_kmeans == False):
      kmeans_param_list = []  # Don't perform these experiments

    kmeans_res_list = []

    for (class_name, dist_meas, centr_init) in kmeans_param_list:

      kmeans_classifier = classification.KMeans(dist_measure = dist_meas,
                                                max_iter_count = 10000,
                                                centroid_init = centr_init)
      start_time = time.time()
      kmeans_classifier.train(sel_w_vec_dict, true_m_set, true_nm_set)
      res_sets = kmeans_classifier.classify(sel_w_vec_dict)
      kmeans_time = time.time() - start_time

      test_results = kmeans_classifier.cross_validate(sel_w_vec_dict,
                                                      true_m_set, true_nm_set,
                                                      n)

      kmeans_acc, kmeans_prec, kmeans_reca, kmeans_fmeas = \
                                                     get_measures(test_results)

      kmeans_res_list.append((class_name, kmeans_acc, kmeans_prec, kmeans_reca,
                              kmeans_fmeas, kmeans_time))
      exp_res_list.append(kmeans_res_list[-1])  # Also append to all results

      (class_name, acc, prec, reca, fmeas, rtime) = kmeans_res_list[-1]
      class_name = class_name.ljust(37)
      res_file.write('    %37s| %5.3f | %5.3f | %5.3f | %5.3f | %9.3f' % \
                     (class_name, acc, prec, reca, fmeas, rtime) + os.linesep)
      res_file.flush()

      del kmeans_classifier

    res_file.write('    ' + '-'*80 + os.linesep)
    res_file.flush()

    # -------------------------------------------------------------------------
    # Farthest-first clustering (with different distance measures)
    #
    res_file.write('    ' + 'Farthest-first clustering'.center(79) + \
                   os.linesep)
    res_file.write('    ' + '-'*80 + os.linesep)

    ffirst_param_list = [('Manhatten dist, traditional init', mymath.distL1,
                          'traditional'),
                         ('Euclidean dist, traditional init', mymath.distL2,
                          'traditional'),
                         ('Cosine dist, traditional init',
                          mymath.distCosine, 'traditional'),
#                        ('L-Infinity dist, traditional init', mymath.distLInf,
#                         'traditional'),
#                        ('Canberra dist, traditional init',
#                         mymath.distCanberra, 'traditional'),

                         ('Manhatten dist, min/max init', mymath.distL1,
                          'min/max'),
                         ('Euclidean dist, min/max init', mymath.distL2,
                          'min/max'),
                         ('Cosine dist, min/max init', mymath.distCosine,
                          'min/max'),
#                         ('L-Infinity dist, min/max init', mymath.distLInf,
#                          'min/max'),
#                         ('Canberra dist, min/max init', mymath.distCanberra,
#                          'min/max'),

                         ('Manhatten dist, mode/max init', mymath.distL1,
                          'mode/max'),
                         ('Euclidean dist, mode/max init', mymath.distL2,
                          'mode/max'),
                         ('Cosine dist, mode/max init', mymath.distCosine,
                          'mode/max')
#                         ('L-Infinity dist, mode/max init', mymath.distLInf,
#                          'mode/max'),
#                         ('Canberra dist, mode/max init', mymath.distCanberra,
#                          'mode/max')
                        ]

    if (do_ffirst == False):
      ffirst_param_list = []  # Don't perform these experiments

    ffirst_res_list = []

    for (class_name, dist_meas, centr_init) in ffirst_param_list:

      ffirst_classifier = classification.FarthestFirst(dist_measu = dist_meas,
                                                       centroid_i = centr_init)
      start_time = time.time()
      ffirst_classifier.train(sel_w_vec_dict, true_m_set, true_nm_set)
      res_sets = ffirst_classifier.classify(sel_w_vec_dict)
      ffirst_time = time.time() - start_time

      test_results = ffirst_classifier.cross_validate(sel_w_vec_dict,
                                                      true_m_set, true_nm_set,
                                                      n)

      ffirst_acc, ffirst_prec, ffirst_reca, ffirst_fmeas = \
                                                     get_measures(test_results)

      ffirst_res_list.append((class_name, ffirst_acc, ffirst_prec, ffirst_reca,
                              ffirst_fmeas, ffirst_time))
      exp_res_list.append(ffirst_res_list[-1])  # Also append to all results

      (class_name, acc, prec, reca, fmeas, rtime) = ffirst_res_list[-1]
      class_name = class_name.ljust(37)
      res_file.write('    %37s| %5.3f | %5.3f | %5.3f | %5.3f | %9.3f' % \
                     (class_name, acc, prec, reca, fmeas, rtime) + os.linesep)
      res_file.flush()

      del ffirst_classifier

    res_file.write('    ' + '-'*80 + os.linesep)
    res_file.flush()

    # -------------------------------------------------------------------------
    # TAILOR unsupervised (hybrid) classifier with different parameters
    #
    res_file.write('    ' + 'TAILOR classifier'.center(79) + os.linesep)
    res_file.write('    ' + '-'*80 + os.linesep)

    tailor_param_list = [
#('Manhatten Linear C=10', mymath.distL1,'LINEAR', 10),
#                         ('Manhatten Linear C= 1', mymath.distL1,'LINEAR',  1),
#                         ('Manhatten Linear C=0.1',mymath.distL1,'LINEAR',0.1),
#                         ('Manhatten Poly C=10', mymath.distL1,  'POLY',   10),
#                         ('Manhatten Poly C= 1', mymath.distL1,  'POLY',    1),
#                         ('Manhatten Poly C=0.1',mymath.distL1,  'POLY',  0.1),
#                         ('Manhatten RBF C=10', mymath.distL1,   'RBF',    10),
#                         ('Manhatten RBF C= 1', mymath.distL1,   'RBF',     1),
#                         ('Manhatten RBF C=0.1',mymath.distL1,   'RBF',   0.1),

                         ('Euclidean Linear C=10', mymath.distL2,'LINEAR', 10),
                         ('Euclidean Linear C= 1', mymath.distL2,'LINEAR',  1),
                         ('Euclidean Linear C=0.1',mymath.distL2,'LINEAR',0.1),
                         ('Euclidean Poly C=10', mymath.distL2,  'POLY',   10),
                         ('Euclidean Poly C= 1', mymath.distL2,  'POLY',    1),
                         ('Euclidean Poly C=0.1',mymath.distL2,  'POLY',  0.1),
                         ('Euclidean RBF C=10', mymath.distL2,   'RBF',    10),
                         ('Euclidean RBF C= 1', mymath.distL2,   'RBF',     1),
                         ('Euclidean RBF C=0.1',mymath.distL2,   'RBF',   0.1),

#                         ('Cosine Linear C=10', mymath.distCosine,'LINEAR',
#                          10),
#                         ('Cosine Linear C= 1', mymath.distCosine,'LINEAR',
#                          1),
#                         ('Cosine Linear C=0.1',mymath.distCosine,'LINEAR',
#                          0.1),
#                         ('Cosine Poly C=10', mymath.distCosine,  'POLY',
#                          10),
#                         ('Cosine Poly C= 1', mymath.distCosine,  'POLY',
#                          1),
#                         ('Cosine Poly C=0.1',mymath.distCosine,  'POLY',
#                          0.1),
#                         ('Cosine RBF C=10', mymath.distCosine,   'RBF',
#                          10),
#                         ('Cosine RBF C= 1', mymath.distCosine,   'RBF',
#                          1),
#                         ('Cosine RBF C=0.1',mymath.distCosine,   'RBF',
#                          0.1),

#                        ('Canberra Linear C=10', mymath.distCanberra,'LINEAR',
#                         10),
#                        ('Canberra Linear C= 1', mymath.distCanberra,'LINEAR',
#                         1),
#                        ('Canberra Linear C=0.1',mymath.distCanberra,'LINEAR',
#                         0.1),
#                        ('Canberra Poly C=10', mymath.distCanberra,  'POLY',
#                         10),
#                        ('Canberra Poly C= 1', mymath.distCanberra,  'POLY',
#                         1),
#                        ('Canberra Poly C=0.1',mymath.distCanberra,  'POLY',
#                         0.1),
#                        ('Canberra RBF C=10', mymath.distCanberra,   'RBF',
#                         10),
#                        ('Canberra RBF C= 1', mymath.distCanberra,   'RBF',
#                         1),
#                        ('Canberra RBF C=0.1',mymath.distCanberra,   'RBF',
#                         0.1),
                        ]

    if (do_tailor == False):
      tailor_param_list = []  # Don't perform these experiments

    tailor_res_list = []

    for (class_name, dist_meas, kernel_name, C_val) in tailor_param_list:

      tailor_classifier = classification.TAILOR(dist_measure = dist_meas,
                                                max_iter_count = 10000,
                                                kernel_type = kernel_name,
                                                C = C_val)
      start_time = time.time()
      tailor_classifier.train(sel_w_vec_dict, true_m_set, true_nm_set)
      res_sets = tailor_classifier.classify(sel_w_vec_dict)
      tailor_time = time.time() - start_time

      test_results = tailor_classifier.cross_validate(sel_w_vec_dict,
                                                      true_m_set, true_nm_set,
                                                      n)

      tailor_acc, tailor_prec, tailor_reca, tailor_fmeas = \
                                                     get_measures(test_results)

      tailor_res_list.append((class_name, tailor_acc, tailor_prec, tailor_reca,
                              tailor_fmeas, tailor_time))
      exp_res_list.append(tailor_res_list[-1])  # Also append to all results

      (class_name, acc, prec, reca, fmeas, rtime) = tailor_res_list[-1]
      class_name = class_name.ljust(37)
      res_file.write('    %37s| %5.3f | %5.3f | %5.3f | %5.3f | %9.3f' % \
                     (class_name, acc, prec, reca, fmeas, rtime) + os.linesep)
      res_file.flush()

      del tailor_classifier

    res_file.write('    ' + '-'*80 + os.linesep)
    res_file.flush()


    # -------------------------------------------------------------------------
    # Two-step unsupervised classifier with thresholds
    #
    res_file.write('    ' + \
                   'Two-step classifier (threshold based)'.center(79) + \
                   os.linesep)
    res_file.write('    ' + '-'*80 + os.linesep)

    two_step_param_list = [
                           ('0.1; NN Manhatten k=1',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('nn', mymath.distL1, 1)),
                           ('0.3; NN Manhatten k=1',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('nn', mymath.distL1, 1)),
                           ('0.5; NN Manhatten k=1',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('nn', mymath.distL1, 1)),
                           ('0.7; NN Manhatten k=1',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('nn', mymath.distL1, 1)),

                           ('0.1; NN Euclidean k=1',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('nn', mymath.distL2, 1)),
                           ('0.3; NN Euclidean k=1',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('nn', mymath.distL2, 1)),
                           ('0.5; NN Euclidean k=1',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('nn', mymath.distL2, 1)),
                           ('0.7; NN Euclidean k=1',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('nn', mymath.distL2, 1)),

#                           ('0.1; NN L-Infinity k=1',
#                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
#                            ('nn', mymath.distLInf, 1)),
#                           ('0.3; NN L-Infinity k=1',
#                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
#                            ('nn', mymath.distLInf, 1)),
#                           ('0.5; NN L-Infinity k=1',
#                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
#                            ('nn', mymath.distLInf, 1)),
#                           ('0.7; NN L-Infinity k=1',
#                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
#                            ('nn', mymath.distLInf, 1)),

                            ('0.1; NN Cosine k=1',
                             (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                             ('nn', mymath.distCosine, 1)),
                            ('0.3; NN Cosine k=1',
                             (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                             ('nn', mymath.distCosine, 1)),
                            ('0.5; NN Cosine k=1',
                             (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                             ('nn', mymath.distCosine, 1)),
                            ('0.7; NN Cosine k=1',
                             (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                             ('nn', mymath.distCosine, 1)),

#                           ('0.1; NN Canberra k=1',
#                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
#                            ('nn', mymath.distCanberra, 1)),
#                           ('0.3; NN Canberra k=1',
#                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
#                            ('nn', mymath.distCanberra, 1)),
#                           ('0.5; NN Canberra k=1',
#                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
#                            ('nn', mymath.distCanberra, 1)),
#                           ('0.7; NN Canberra k=1',
#                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
#                            ('nn', mymath.distCanberra, 1)),


                           ('0.1; NN Manhatten k=3',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('nn', mymath.distL1, 3)),
                           ('0.3; NN Manhatten k=3',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('nn', mymath.distL1, 3)),
                           ('0.5; NN Manhatten k=3',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('nn', mymath.distL1, 3)),
                           ('0.7; NN Manhatten k=3',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('nn', mymath.distL1, 3)),

                           ('0.1; NN Euclidean k=3',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('nn', mymath.distL2, 3)),
                           ('0.3; NN Euclidean k=3',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('nn', mymath.distL2, 3)),
                           ('0.5; NN Euclidean k=3',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('nn', mymath.distL2, 3)),
                           ('0.7; NN Euclidean k=3',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('nn', mymath.distL2, 3)),

#                           ('0.1; NN L-Infinity k=3',
#                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
#                            ('nn', mymath.distLInf, 3)),
#                           ('0.3; NN L-Infinity k=3',
#                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
#                            ('nn', mymath.distLInf, 3)),
#                           ('0.5; NN L-Infinity k=3',
#                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
#                            ('nn', mymath.distLInf, 3)),
#                           ('0.7; NN L-Infinity k=3',
#                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
#                            ('nn', mymath.distLInf, 3)),

                            ('0.1; NN Cosine k=3',
                             (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                             ('nn', mymath.distCosine, 3)),
                            ('0.3; NN Cosine k=3',
                             (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                             ('nn', mymath.distCosine, 3)),
                            ('0.5; NN Cosine k=3',
                             (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                             ('nn', mymath.distCosine, 3)),
                            ('0.7; NN Cosine k=3',
                             (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                             ('nn', mymath.distCosine, 3)),

#                           ('0.1; NN Canberra k=3',
#                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
#                            ('nn', mymath.distCanberra, 3)),
#                           ('0.3; NN Canberra k=3',
#                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
#                            ('nn', mymath.distCanberra, 3)),
#                           ('0.5; NN Canberra k=3',
#                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
#                            ('nn', mymath.distCanberra, 3)),
#                           ('0.7; NN Canberra k=3',
#                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
#                            ('nn', mymath.distCanberra, 3)),


                           ('0.1; NN Manhatten k=9',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('nn', mymath.distL1, 9)),
                           ('0.3; NN Manhatten k=9',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('nn', mymath.distL1, 9)),
                           ('0.5; NN Manhatten k=9',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('nn', mymath.distL1, 9)),
                           ('0.7; NN Manhatten k=9',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('nn', mymath.distL1, 9)),

                           ('0.1; NN Euclidean k=9',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('nn', mymath.distL2, 9)),
                           ('0.3; NN Euclidean k=9',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('nn', mymath.distL2, 9)),
                           ('0.5; NN Euclidean k=9',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('nn', mymath.distL2, 9)),
                           ('0.7; NN Euclidean k=9',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('nn', mymath.distL2, 9)),

#                           ('0.1; NN L-Infinity k=9',
#                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
#                            ('nn', mymath.distLInf, 9)),
#                           ('0.3; NN L-Infinity k=9',
#                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
#                            ('nn', mymath.distLInf, 9)),
#                           ('0.5; NN L-Infinity k=9',
#                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
#                            ('nn', mymath.distLInf, 9)),
#                           ('0.7; NN L-Infinity k=9',
#                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
#                            ('nn', mymath.distLInf, 9)),

                            ('0.1; NN Cosine k=9',
                             (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                             ('nn', mymath.distCosine, 9)),
                            ('0.3; NN Cosine k=9',
                             (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                             ('nn', mymath.distCosine, 9)),
                            ('0.5; NN Cosine k=9',
                             (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                             ('nn', mymath.distCosine, 9)),
                            ('0.7; NN Cosine k=9',
                             (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                             ('nn', mymath.distCosine, 9)),

#                           ('0.1; NN Canberra k=9',
#                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
#                            ('nn', mymath.distCanberra, 9)),
#                           ('0.3; NN Canberra k=9',
#                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
#                            ('nn', mymath.distCanberra, 9)),
#                           ('0.5; NN Canberra k=9',
#                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
#                            ('nn', mymath.distCanberra, 9)),
#                           ('0.7; NN Canberra k=9',
#                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
#                            ('nn', mymath.distCanberra, 9)),


                           ('0.1; K-means Manhatten dist 0',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('kmeans', mymath.distL1, 0)),
                           ('0.3; K-means Manhatten dist 0',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('kmeans', mymath.distL1, 0)),
                           ('0.5; K-means Manhatten dist 0',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('kmeans', mymath.distL1, 0)),
                           ('0.7; K-means Manhatten dist 0',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('kmeans', mymath.distL1, 0)),

                           ('0.1; K-means Euclidean dist 0',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('kmeans', mymath.distL2, 0)),
                           ('0.3; K-means Euclidean dist 0',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('kmeans', mymath.distL2, 0)),
                           ('0.5; K-means Euclidean dist 0',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('kmeans', mymath.distL2, 0)),
                           ('0.7; K-means Euclidean dist 0',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('kmeans', mymath.distL2, 0)),

#                           ('0.1; K-means L-Infinity dist 0',
#                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
#                            ('kmeans', mymath.distLInf, 0)),
#                           ('0.3; K-means L-Infinity dist 0',
#                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
#                            ('kmeans', mymath.distLInf, 0)),
#                           ('0.5; K-means L-Infinity dist 0',
#                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
#                            ('kmeans', mymath.distLInf, 0)),
#                           ('0.7; K-means L-Infinity dist 0',
#                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
#                            ('kmeans', mymath.distLInf, 0)),

                            ('0.1; K-means Cosine dist 0',
                             (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                             ('kmeans', mymath.distCosine, 0)),
                            ('0.3; K-means Cosine dist 0',
                             (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                             ('kmeans', mymath.distCosine, 0)),
                            ('0.5; K-means Cosine dist 0',
                             (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                             ('kmeans', mymath.distCosine, 0)),
                            ('0.7; K-means Cosine dist 0',
                             (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                             ('kmeans', mymath.distCosine, 0)),

#                           ('0.1; K-means Canberra dist 0',
#                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
#                            ('kmeans', mymath.distCanberra, 0)),
#                           ('0.3; K-means Canberra dist 0',
#                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
#                            ('kmeans', mymath.distCanberra, 0)),
#                           ('0.5; K-means Canberra dist 0',
#                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
#                            ('kmeans', mymath.distCanberra, 0)),
#                           ('0.7; K-means Canberra dist 0',
#                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
#                            ('kmeans', mymath.distCanberra, 0)),

                           ('0.1; K-means Manhatten dist',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('kmeans', mymath.distL1, 10000)),
                           ('0.3; K-means Manhatten dist',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('kmeans', mymath.distL1, 10000)),
                           ('0.5; K-means Manhatten dist',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('kmeans', mymath.distL1, 10000)),
                           ('0.7; K-means Manhatten dist',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('kmeans', mymath.distL1, 10000)),

                           ('0.1; K-means Euclidean dist',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('kmeans', mymath.distL2, 10000)),
                           ('0.3; K-means Euclidean dist',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('kmeans', mymath.distL2, 10000)),
                           ('0.5; K-means Euclidean dist',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('kmeans', mymath.distL2, 10000)),
                           ('0.7; K-means Euclidean dist',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('kmeans', mymath.distL2, 10000)),

#                           ('0.1; K-means L-Infinity dist',
#                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
#                            ('kmeans', mymath.distLInf, 10000)),
#                           ('0.3; K-means L-Infinity dist',
#                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
#                            ('kmeans', mymath.distLInf, 10000)),
#                           ('0.5; K-means L-Infinity dist',
#                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
#                            ('kmeans', mymath.distLInf, 10000)),
#                           ('0.7; K-means L-Infinity dist',
#                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
#                            ('kmeans', mymath.distLInf, 10000)),

                            ('0.1; K-means Cosine dist',
                             (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                             ('kmeans', mymath.distCosine, 10000)),
                            ('0.3; K-means Cosine dist',
                             (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                             ('kmeans', mymath.distCosine, 10000)),
                            ('0.5; K-means Cosine dist',
                             (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                             ('kmeans', mymath.distCosine, 10000)),
                            ('0.7; K-means Cosine dist',
                             (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                             ('kmeans', mymath.distCosine, 10000)),

#                           ('0.1; K-means Canberra dist',
#                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
#                            ('kmeans', mymath.distCanberra, 10000)),
#                           ('0.3; K-means Canberra dist',
#                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
#                            ('kmeans', mymath.distCanberra, 10000)),
#                           ('0.5; K-means Canberra dist',
#                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
#                            ('kmeans', mymath.distCanberra, 10000)),
#                           ('0.7; K-means Canberra dist',
#                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
#                            ('kmeans', mymath.distCanberra, 10000)),


                           ('0.1; SVM linear  C=10 0-0',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'LINEAR', 10, 0, 0)),
                           ('0.3; SVM linear  C=10 0-0',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'LINEAR', 10, 0, 0)),
                           ('0.5; SVM linear  C=10 0-0',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'LINEAR', 10, 0, 0)),
                           ('0.7; SVM linear  C=10 0-0',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'LINEAR', 10, 0, 0)),

                           ('0.1; SVM linear  C= 1 0-0',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'LINEAR', 1, 0, 0)),
                           ('0.3; SVM linear  C= 1 0-0',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'LINEAR', 1, 0, 0)),
                           ('0.5; SVM linear  C= 1 0-0',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'LINEAR', 1, 0, 0)),
                           ('0.7; SVM linear  C= 1 0-0',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'LINEAR', 1, 0, 0)),

                           ('0.1; SVM linear  C=0.1 0-0',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'LINEAR', 0.1, 0, 0)),
                           ('0.3; SVM linear  C=0.1 0-0',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'LINEAR', 0.1, 0, 0)),
                           ('0.5; SVM linear  C=0.1 0-0',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'LINEAR', 0.1, 0, 0)),
                           ('0.7; SVM linear  C=0.1 0-0',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'LINEAR', 0.1, 0, 0)),


                           ('0.1; SVM poly  C=10 0-0',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'POLY', 10, 0, 0)),
                           ('0.3; SVM poly  C=10 0-0',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'POLY', 10, 0, 0)),
                           ('0.5; SVM poly  C=10 0-0',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'POLY', 10, 0, 0)),
                           ('0.7; SVM poly  C=10 0-0',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'POLY', 10, 0, 0)),

                           ('0.1; SVM poly  C= 1 0-0',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'POLY', 1, 0, 0)),
                           ('0.3; SVM poly  C= 1 0-0',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'POLY', 1, 0, 0)),
                           ('0.5; SVM poly  C= 1 0-0',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'POLY', 1, 0, 0)),
                           ('0.7; SVM poly  C= 1 0-0',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'POLY', 1, 0, 0)),

                           ('0.1; SVM poly  C=0.1 0-0',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'POLY', 0.1, 0, 0)),
                           ('0.3; SVM poly  C=0.1 0-0',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'POLY', 0.1, 0, 0)),
                           ('0.5; SVM poly  C=0.1 0-0',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'POLY', 0.1, 0, 0)),
                           ('0.7; SVM poly  C=0.1 0-0',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'POLY', 0.1, 0, 0)),


                           ('0.1; SVM RBF  C=10 0-0',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'RBF', 10, 0, 0)),
                           ('0.3; SVM RBF  C=10 0-0',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'RBF', 10, 0, 0)),
                           ('0.5; SVM RBF  C=10 0-0',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'RBF', 10, 0, 0)),
                           ('0.7; SVM RBF  C=10 0-0',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'RBF', 10, 0, 0)),

                           ('0.1; SVM RBF  C= 1 0-0',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'RBF', 1, 0, 0)),
                           ('0.3; SVM RBF  C= 1 0-0',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'RBF', 1, 0, 0)),
                           ('0.5; SVM RBF  C= 1 0-0',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'RBF', 1, 0, 0)),
                           ('0.7; SVM RBF  C= 1 0-0',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'RBF', 1, 0, 0)),

                           ('0.1; SVM RBF  C=0.1 0-0',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'RBF', 0.1, 0, 0)),
                           ('0.3; SVM RBF  C=0.1 0-0',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'RBF', 0.1, 0, 0)),
                           ('0.5; SVM RBF  C=0.1 0-0',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'RBF', 0.1, 0, 0)),
                           ('0.7; SVM RBF  C=0.1 0-0',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'RBF', 0.1, 0, 0)),


                           ('0.1; SVM linear  C=10 25-25',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'LINEAR', 10, 25,25)),
                           ('0.3; SVM linear  C=10 25-25',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'LINEAR', 10, 25,25)),
                           ('0.5; SVM linear  C=10 25-25',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'LINEAR', 10, 25,25)),
                           ('0.7; SVM linear  C=10 25-25',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'LINEAR', 10, 25,25)),

                           ('0.1; SVM linear  C= 1 25-25',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'LINEAR', 1, 25,25)),
                           ('0.3; SVM linear  C= 1 25-25',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'LINEAR', 1, 25,25)),
                           ('0.5; SVM linear  C= 1 25-25',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'LINEAR', 1, 25,25)),
                           ('0.7; SVM linear  C= 1 25-25',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'LINEAR', 1, 25,25)),

                           ('0.1; SVM linear  C=0.1 25-25',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'LINEAR', 0.1, 25,25)),
                           ('0.3; SVM linear  C=0.1 25-25',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'LINEAR', 0.1, 25,25)),
                           ('0.5; SVM linear  C=0.1 25-25',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'LINEAR', 0.1, 25,25)),
                           ('0.7; SVM linear  C=0.1 25-25',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'LINEAR', 0.1, 25,25)),


                           ('0.1; SVM poly  C=10 25-25',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'POLY', 10, 25,25)),
                           ('0.3; SVM poly  C=10 25-25',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'POLY', 10, 25,25)),
                           ('0.5; SVM poly  C=10 25-25',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'POLY', 10, 25,25)),
                           ('0.7; SVM poly  C=10 25-25',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'POLY', 10, 25,25)),

                           ('0.1; SVM poly  C= 1 25-25',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'POLY', 1, 25,25)),
                           ('0.3; SVM poly  C= 1 25-25',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'POLY', 1, 25,25)),
                           ('0.5; SVM poly  C= 1 25-25',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'POLY', 1, 25,25)),
                           ('0.7; SVM poly  C= 1 25-25',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'POLY', 1, 25,25)),

                           ('0.1; SVM poly  C=0.1 25-25',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'POLY', 0.1, 25,25)),
                           ('0.3; SVM poly  C=0.1 25-25',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'POLY', 0.1, 25,25)),
                           ('0.5; SVM poly  C=0.1 25-25',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'POLY', 0.1, 25,25)),
                           ('0.7; SVM poly  C=0.1 25-25',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'POLY', 0.1, 25,25)),


                           ('0.1; SVM RBF  C=10 25-25',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'RBF', 10, 25,25)),
                           ('0.3; SVM RBF  C=10 25-25',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'RBF', 10, 25,25)),
                           ('0.5; SVM RBF  C=10 25-25',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'RBF', 10, 25,25)),
                           ('0.7; SVM RBF  C=10 25-25',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'RBF', 10, 25,25)),

                           ('0.1; SVM RBF  C= 1 25-25',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'RBF', 1, 25,25)),
                           ('0.3; SVM RBF  C= 1 25-25',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'RBF', 1, 25,25)),
                           ('0.5; SVM RBF  C= 1 25-25',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'RBF', 1, 25,25)),
                           ('0.7; SVM RBF  C= 1 25-25',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'RBF', 1, 25,25)),

                           ('0.1; SVM RBF  C=0.1 25-25',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'RBF', 0.1, 25,25)),
                           ('0.3; SVM RBF  C=0.1 25-25',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'RBF', 0.1, 25,25)),
                           ('0.5; SVM RBF  C=0.1 25-25',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'RBF', 0.1, 25,25)),
                           ('0.7; SVM RBF  C=0.1 25-25',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'RBF', 0.1, 25,25)),


                           ('0.1; SVM linear  C=10 25-50',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'LINEAR', 10, 25,50)),
                           ('0.3; SVM linear  C=10 25-50',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'LINEAR', 10, 25,50)),
                           ('0.5; SVM linear  C=10 25-50',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'LINEAR', 10, 25,50)),
                           ('0.7; SVM linear  C=10 25-50',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'LINEAR', 10, 25,50)),

                           ('0.1; SVM linear  C= 1 25-50',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'LINEAR', 1, 25,50)),
                           ('0.3; SVM linear  C= 1 25-50',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'LINEAR', 1, 25,50)),
                           ('0.5; SVM linear  C= 1 25-50',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'LINEAR', 1, 25,50)),
                           ('0.7; SVM linear  C= 1 25-50',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'LINEAR', 1, 25,50)),

                           ('0.1; SVM linear  C=0.1 25-50',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'LINEAR', 0.1, 25,50)),
                           ('0.3; SVM linear  C=0.1 25-50',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'LINEAR', 0.1, 25,50)),
                           ('0.5; SVM linear  C=0.1 25-50',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'LINEAR', 0.1, 25,50)),
                           ('0.7; SVM linear  C=0.1 25-50',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'LINEAR', 0.1, 25,50)),


                           ('0.1; SVM poly  C=10 25-50',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'POLY', 10, 25,50)),
                           ('0.3; SVM poly  C=10 25-50',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'POLY', 10, 25,50)),
                           ('0.5; SVM poly  C=10 25-50',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'POLY', 10, 25,50)),
                           ('0.7; SVM poly  C=10 25-50',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'POLY', 10, 25,50)),

                           ('0.1; SVM poly  C= 1 25-50',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'POLY', 1, 25,50)),
                           ('0.3; SVM poly  C= 1 25-50',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'POLY', 1, 25,50)),
                           ('0.5; SVM poly  C= 1 25-50',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'POLY', 1, 25,50)),
                           ('0.7; SVM poly  C= 1 25-50',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'POLY', 1, 25,50)),

                           ('0.1; SVM poly  C=0.1 25-50',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'POLY', 0.1, 25,50)),
                           ('0.3; SVM poly  C=0.1 25-50',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'POLY', 0.1, 25,50)),
                           ('0.5; SVM poly  C=0.1 25-50',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'POLY', 0.1, 25,50)),
                           ('0.7; SVM poly  C=0.1 25-50',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'POLY', 0.1, 25,50)),


                           ('0.1; SVM RBF  C=10 25-50',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'RBF', 10, 25,50)),
                           ('0.3; SVM RBF  C=10 25-50',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'RBF', 10, 25,50)),
                           ('0.5; SVM RBF  C=10 25-50',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'RBF', 10, 25,50)),
                           ('0.7; SVM RBF  C=10 25-50',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'RBF', 10, 25,50)),

                           ('0.1; SVM RBF  C= 1 25-50',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'RBF', 1, 25,50)),
                           ('0.3; SVM RBF  C= 1 25-50',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'RBF', 1, 25,50)),
                           ('0.5; SVM RBF  C= 1 25-50',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'RBF', 1, 25,50)),
                           ('0.7; SVM RBF  C= 1 25-50',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'RBF', 1, 25,50)),

                           ('0.1; SVM RBF  C=0.1 25-50',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'RBF', 0.1, 25,50)),
                           ('0.3; SVM RBF  C=0.1 25-50',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'RBF', 0.1, 25,50)),
                           ('0.5; SVM RBF  C=0.1 25-50',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'RBF', 0.1, 25,50)),
                           ('0.7; SVM RBF  C=0.1 25-50',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'RBF', 0.1, 25,50)),


                           ('0.1; SVM linear  C=10 50-100',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'LINEAR', 10, 50,100)),
                           ('0.3; SVM linear  C=10 50-100',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'LINEAR', 10, 50,100)),
                           ('0.5; SVM linear  C=10 50-100',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'LINEAR', 10, 50,100)),
                           ('0.7; SVM linear  C=10 50-100',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'LINEAR', 10, 50,100)),

                           ('0.1; SVM linear  C= 1 50-100',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'LINEAR', 1, 50,100)),
                           ('0.3; SVM linear  C= 1 50-100',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'LINEAR', 1, 50,100)),
                           ('0.5; SVM linear  C= 1 50-100',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'LINEAR', 1, 50,100)),
                           ('0.7; SVM linear  C= 1 50-100',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'LINEAR', 1, 50,100)),

                           ('0.1; SVM linear  C=0.1 50-100',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'LINEAR', 0.1, 50,100)),
                           ('0.3; SVM linear  C=0.1 50-100',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'LINEAR', 0.1, 50,100)),
                           ('0.5; SVM linear  C=0.1 50-100',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'LINEAR', 0.1, 50,100)),
                           ('0.7; SVM linear  C=0.1 50-100',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'LINEAR', 0.1, 50,100)),


                           ('0.1; SVM poly  C=10 50-100',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'POLY', 10, 50,100)),
                           ('0.3; SVM poly  C=10 50-100',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'POLY', 10, 50,100)),
                           ('0.5; SVM poly  C=10 50-100',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'POLY', 10, 50,100)),
                           ('0.7; SVM poly  C=10 50-100',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'POLY', 10, 50,100)),

                           ('0.1; SVM poly  C= 1 50-100',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'POLY', 1, 50,100)),
                           ('0.3; SVM poly  C= 1 50-100',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'POLY', 1, 50,100)),
                           ('0.5; SVM poly  C= 1 50-100',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'POLY', 1, 50,100)),
                           ('0.7; SVM poly  C= 1 50-100',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'POLY', 1, 50,100)),

                           ('0.1; SVM poly  C=0.1 50-100',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'POLY', 0.1, 50,100)),
                           ('0.3; SVM poly  C=0.1 50-100',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'POLY', 0.1, 50,100)),
                           ('0.5; SVM poly  C=0.1 50-100',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'POLY', 0.1, 50,100)),
                           ('0.7; SVM poly  C=0.1 50-100',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'POLY', 0.1, 50,100)),


                           ('0.1; SVM RBF  C=10 50-100',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'RBF', 10, 50,100)),
                           ('0.3; SVM RBF  C=10 50-100',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'RBF', 10, 50,100)),
                           ('0.5; SVM RBF  C=10 50-100',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'RBF', 10, 50,100)),
                           ('0.7; SVM RBF  C=10 50-100',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'RBF', 10, 50,100)),

                           ('0.1; SVM RBF  C= 1 50-100',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'RBF', 1, 50,100)),
                           ('0.3; SVM RBF  C= 1 50-100',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'RBF', 1, 50,100)),
                           ('0.5; SVM RBF  C= 1 50-100',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'RBF', 1, 50,100)),
                           ('0.7; SVM RBF  C= 1 50-100',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'RBF', 1, 50,100)),

                           ('0.1; SVM RBF  C=0.1 50-100',
                            (1.0,'threshold',0.1), (0.0,'threshold',0.1),
                            ('svm', 'RBF', 0.1, 50,100)),
                           ('0.3; SVM RBF  C=0.1 50-100',
                            (1.0,'threshold',0.3), (0.0,'threshold',0.3),
                            ('svm', 'RBF', 0.1, 50,100)),
                           ('0.5; SVM RBF  C=0.1 50-100',
                            (1.0,'threshold',0.5), (0.0,'threshold',0.5),
                            ('svm', 'RBF', 0.1, 50,100)),
                           ('0.7; SVM RBF  C=0.1 50-100',
                            (1.0,'threshold',0.7), (0.0,'threshold',0.7),
                            ('svm', 'RBF', 0.1, 50,100))
                           ]

    if (do_two_step_thres == False):
      two_step_param_list = []  # Don't perform these experiments

    two_step_res_list = []

    for (class_name, s1_m, s1_nm, s2_class) in two_step_param_list:

      # Loop over possible random selection methods
      # (for each: string to add to output, random selection tuple, number of
      # iterations to do)
      #
      for rand_sel in [('; ',       None,                     1),
#                       ('; RS-U1',  ('uniform',      1,  1), 10),
#                       ('; RS-U10', ('uniform',     10, 10), 10),
#                       ('; RS-L1',  ('linear',       1,  1), 10),
#                       ('; RS-L10', ('linear',      10, 10), 10),
#                       ('; RS-E1',  ('exponential',  1,  1), 10),
#                       ('; RS-E10', ('exponential', 10, 10), 10)
                       ]:

        (rand_sel_str, rand_sel_tup, num_iter) = rand_sel  # Unpack details

        acc_res_list =   []
        prec_res_list =  []
        reca_res_list =  []
        fmeas_res_list = []
        time_res_list =  []

        for iter in range(num_iter):  # Loop over number of iterations - - - -
          two_step_classifier = classification.TwoStep(s1_match_method = s1_m,
                                                    s1_non_match_meth = s1_nm,
                                                    random_sel = rand_sel_tup,
                                                    s2_classifier = s2_class)
          start_time = time.time()
          two_step_classifier.train(sel_w_vec_dict, true_m_set, true_nm_set)
          res_sets = two_step_classifier.classify(sel_w_vec_dict)
          two_step_time = time.time() - start_time
          time_res_list.append(two_step_time)

          test_results = two_step_classifier.cross_validate(sel_w_vec_dict,
                                                            true_m_set,
                                                            true_nm_set, n)

          two_step_acc, two_step_prec, two_step_reca, two_step_fmeas = \
                                                     get_measures(test_results)
          acc_res_list.append(two_step_acc)
          prec_res_list.append(two_step_prec)
          reca_res_list.append(two_step_reca)
          fmeas_res_list.append(two_step_fmeas)

          del two_step_classifier

        two_step_acc =   sum(acc_res_list) /   num_iter
        two_step_prec =  sum(prec_res_list) /  num_iter
        two_step_reca =  sum(reca_res_list) /  num_iter
        two_step_fmeas = sum(fmeas_res_list) / num_iter
        two_step_time = sum(time_res_list)   / num_iter

        this_class_name = class_name + rand_sel_str

        two_step_res_list.append((this_class_name, two_step_acc, two_step_prec,
                                  two_step_reca, two_step_fmeas,two_step_time))
        exp_res_list.append(two_step_res_list[-1]) # Also append to all results

        (this_class_name, acc, prec, reca, fmeas,rtime) = two_step_res_list[-1]
        this_class_name = this_class_name.ljust(37)
        res_file.write('    %37s| %5.3f | %5.3f | %5.3f | %5.3f | %9.3f' % \
                       (this_class_name, acc, prec, reca, fmeas, rtime) + \
                       os.linesep)
        res_file.flush()

    res_file.write('    ' + '-'*80 + os.linesep)
    res_file.flush()

    # -------------------------------------------------------------------------
    # Two-step unsupervised classifier with nearest
    #
    res_file.write('    ' + \
                   'Two-step classifier (nearest based)'.center(79) + \
                   os.linesep)
    res_file.write('    ' + '-'*80 + os.linesep)

    two_step_param_list = [
##                           ('1% NU; NN Manhatten k=1',
##                            (1.0,'nearest', 1, False),
##                            (0.0,'nearest', 1, False),
##                            ('nn', mymath.distL1, 1)),
##                           ('1%  U; NN Manhatten k=1',
##                            (1.0,'nearest', 1, True),
##                            (0.0,'nearest', 1, True),
##                            ('nn', mymath.distL1, 1)),
#                           ('5% NU; NN Manhatten k=1',
#                            (1.0,'nearest', 5, False),
#                            (0.0,'nearest', 5, False),
#                            ('nn', mymath.distL1, 1)),
#                           ('5%  U; NN Manhatten k=1',
#                            (1.0,'nearest', 5, True),
#                            (0.0,'nearest', 5, True),
#                            ('nn', mymath.distL1, 1)),
#                           ('10% NU; NN Manhatten k=1',
#                            (1.0,'nearest', 10, False),
#                            (0.0,'nearest', 10, False),
#                            ('nn', mymath.distL1, 1)),
#                           ('10%  U; NN Manhatten k=1',
#                            (1.0,'nearest', 10, True),
#                            (0.0,'nearest', 10, True),
#                            ('nn', mymath.distL1, 1)),
#
##                           ('1% NU; NN Euclidean k=1',
##                            (1.0,'nearest', 1, False),
##                            (0.0,'nearest', 1, False),
##                            ('nn', mymath.distL2, 1)),
##                           ('1%  U; NN Euclidean k=1',
##                            (1.0,'nearest', 1, True),
##                            (0.0,'nearest', 1, True),
##                            ('nn', mymath.distL2, 1)),
#                           ('5% NU; NN Euclidean k=1',
#                            (1.0,'nearest', 5, False),
#                            (0.0,'nearest', 5, False),
#                            ('nn', mymath.distL2, 1)),
#                           ('5%  U; NN Euclidean k=1',
#                            (1.0,'nearest', 5, True),
#                            (0.0,'nearest', 5, True),
#                            ('nn', mymath.distL2, 1)),
#                           ('10% NU; NN Euclidean k=1',
#                            (1.0,'nearest', 10, False),
#                            (0.0,'nearest', 10, False),
#                            ('nn', mymath.distL2, 1)),
#                           ('10%  U; NN Euclidean k=1',
#                            (1.0,'nearest', 10, True),
#                            (0.0,'nearest', 10, True),
#                            ('nn', mymath.distL2, 1)),
#
##                           ('1% NU; NN L-Infinity k=1',
##                            (1.0,'nearest', 1, False),
##                            (0.0,'nearest', 1, False),
##                            ('nn', mymath.distLInf, 1)),
##                           ('1%  U; NN L-Infinity k=1',
##                            (1.0,'nearest', 1, True),
##                            (0.0,'nearest', 1, True),
##                            ('nn', mymath.distLInf, 1)),
##                           ('5% NU; NN L-Infinity k=1',
##                            (1.0,'nearest', 5, False),
##                            (0.0,'nearest', 5, False),
##                            ('nn', mymath.distLInf, 1)),
##                           ('5%  U; NN L-Infinity k=1',
##                            (1.0,'nearest', 5, True),
##                            (0.0,'nearest', 5, True),
##                            ('nn', mymath.distLInf, 1)),
##                           ('10% NU; NN L-Infinity k=1',
##                            (1.0,'nearest', 10, False),
##                            (0.0,'nearest', 10, False),
##                            ('nn', mymath.distLInf, 1)),
##                           ('10%  U; NN L-Infinity k=1',
##                            (1.0,'nearest', 10, True),
##                            (0.0,'nearest', 10, True),
##                            ('nn', mymath.distLInf, 1)),
#
##                            ('1% NU; NN Cosine k=1',
##                             (1.0,'nearest', 1, False),
##                             (0.0,'nearest', 1, False),
##                             ('nn', mymath.distCosine, 1)),
##                            ('1%  U; NN Cosine k=1',
##                             (1.0,'nearest', 1, True),
##                             (0.0,'nearest', 1, True),
##                             ('nn', mymath.distCosine, 1)),
##                            ('5% NU; NN Cosine k=1',
##                             (1.0,'nearest', 5, False),
##                             (0.0,'nearest', 5, False),
##                             ('nn', mymath.distCosine, 1)),
##                            ('5%  U; NN Cosine k=1',
##                             (1.0,'nearest', 5, True),
##                             (0.0,'nearest', 5, True),
##                             ('nn', mymath.distCosine, 1)),
##                            ('10% NU; NN Cosine k=1',
##                             (1.0,'nearest', 10, False),
##                             (0.0,'nearest', 10, False),
##                             ('nn', mymath.distCosine, 1)),
##                            ('10%  U; NN Cosine k=1',
##                             (1.0,'nearest', 10, True),
##                             (0.0,'nearest', 10, True),
##                             ('nn', mymath.distCosine, 1)),
#
##                           ('1% NU; NN Canberra k=1',
##                            (1.0,'nearest', 1, False),
##                            (0.0,'nearest', 1, False),
##                            ('nn', mymath.distCanberra, 1)),
##                           ('1%  U; NN Canberra k=1',
##                            (1.0,'nearest', 1, True),
##                            (0.0,'nearest', 1, True),
##                            ('nn', mymath.distCanberra, 1)),
##                           ('5% NU; NN Canberra k=1',
##                            (1.0,'nearest', 5, False),
##                            (0.0,'nearest', 5, False),
##                            ('nn', mymath.distCanberra, 1)),
##                           ('5%  U; NN Canberra k=1',
##                            (1.0,'nearest', 5, True),
##                            (0.0,'nearest', 5, True),
##                            ('nn', mymath.distCanberra, 1)),
##                           ('10% NU; NN Canberra k=1',
##                            (1.0,'nearest', 10, False),
##                            (0.0,'nearest', 10, False),
##                            ('nn', mymath.distCanberra, 1)),
##                           ('10%  U; NN Canberra k=1',
##                            (1.0,'nearest', 10, True),
##                            (0.0,'nearest', 10, True),
##                            ('nn', mymath.distCanberra, 1)),
#
##                           ('1% NU; NN Manhatten k=3',
##                            (1.0,'nearest', 1, False),
##                            (0.0,'nearest', 1, False),
##                            ('nn', mymath.distL1, 3)),
##                           ('1%  U; NN Manhatten k=3',
##                            (1.0,'nearest', 1, True),
##                            (0.0,'nearest', 1, True),
##                            ('nn', mymath.distL1, 3)),
#                           ('5% NU; NN Manhatten k=3',
#                            (1.0,'nearest', 5, False),
#                            (0.0,'nearest', 5, False),
#                            ('nn', mymath.distL1, 3)),
#                           ('5%  U; NN Manhatten k=3',
#                            (1.0,'nearest', 5, True),
#                            (0.0,'nearest', 5, True),
#                            ('nn', mymath.distL1, 3)),
#                           ('10% NU; NN Manhatten k=3',
#                            (1.0,'nearest', 10, False),
#                            (0.0,'nearest', 10, False),
#                            ('nn', mymath.distL1, 3)),
#                           ('10%  U; NN Manhatten k=3',
#                            (1.0,'nearest', 10, True),
#                            (0.0,'nearest', 10, True),
#                            ('nn', mymath.distL1, 3)),
#
##                           ('1% NU; NN Euclidean k=3',
##                            (1.0,'nearest', 1, False),
##                            (0.0,'nearest', 1, False),
##                            ('nn', mymath.distL2, 3)),
##                           ('1%  U; NN Euclidean k=3',
##                            (1.0,'nearest', 1, True),
##                            (0.0,'nearest', 1, True),
##                            ('nn', mymath.distL2, 3)),
#                           ('5% NU; NN Euclidean k=3',
#                            (1.0,'nearest', 5, False),
#                            (0.0,'nearest', 5, False),
#                            ('nn', mymath.distL2, 3)),
#                           ('5%  U; NN Euclidean k=3',
#                            (1.0,'nearest', 5, True),
#                            (0.0,'nearest', 5, True),
#                            ('nn', mymath.distL2, 3)),
#                           ('10% NU; NN Euclidean k=3',
#                            (1.0,'nearest', 10, False),
#                            (0.0,'nearest', 10, False),
#                            ('nn', mymath.distL2, 3)),
#                           ('10%  U; NN Euclidean k=3',
#                            (1.0,'nearest', 10, True),
#                            (0.0,'nearest', 10, True),
#                            ('nn', mymath.distL2, 3)),
#
##                           ('1% NU; NN L-Infinity k=3',
##                            (1.0,'nearest', 1, False),
##                            (0.0,'nearest', 1, False),
##                            ('nn', mymath.distLInf, 3)),
##                           ('1%  U; NN L-Infinity k=3',
##                            (1.0,'nearest', 1, True),
##                            (0.0,'nearest', 1, True),
##                            ('nn', mymath.distLInf, 3)),
##                           ('5% NU; NN L-Infinity k=3',
##                            (1.0,'nearest', 5, False),
##                            (0.0,'nearest', 5, False),
##                            ('nn', mymath.distLInf, 3)),
##                           ('5%  U; NN L-Infinity k=3',
##                            (1.0,'nearest', 5, True),
##                            (0.0,'nearest', 5, True),
##                            ('nn', mymath.distLInf, 3)),
##                           ('10% NU; NN L-Infinity k=3',
##                            (1.0,'nearest', 10, False),
##                            (0.0,'nearest', 10, False),
##                            ('nn', mymath.distLInf, 3)),
##                           ('10%  U; NN L-Infinity k=3',
##                            (1.0,'nearest', 10, True),
##                            (0.0,'nearest', 10, True),
##                            ('nn', mymath.distLInf, 3)),
#
##                            ('1% NU; NN Cosine k=3',
##                             (1.0,'nearest', 1, False),
##                             (0.0,'nearest', 1, False),
##                             ('nn', mymath.distCosine, 3)),
##                            ('1%  U; NN Cosine k=3',
##                             (1.0,'nearest', 1, True),
##                             (0.0,'nearest', 1, True),
##                             ('nn', mymath.distCosine, 3)),
##                            ('5% NU; NN Cosine k=3',
##                             (1.0,'nearest', 5, False),
##                             (0.0,'nearest', 5, False),
##                             ('nn', mymath.distCosine, 3)),
##                            ('5%  U; NN Cosine k=3',
##                             (1.0,'nearest', 5, True),
##                             (0.0,'nearest', 5, True),
##                             ('nn', mymath.distCosine, 3)),
##                            ('10% NU; NN Cosine k=3',
##                             (1.0,'nearest', 10, False),
##                             (0.0,'nearest', 10, False),
##                             ('nn', mymath.distCosine, 3)),
##                            ('10%  U; NN Cosine k=3',
##                             (1.0,'nearest', 10, True),
##                             (0.0,'nearest', 10, True),
##                             ('nn', mymath.distCosine, 3)),
#
##                           ('1% NU; NN Canberra k=3',
##                            (1.0,'nearest', 1, False),
##                            (0.0,'nearest', 1, False),
##                            ('nn', mymath.distCanberra, 3)),
##                           ('1%  U; NN Canberra k=3',
##                            (1.0,'nearest', 1, True),
##                            (0.0,'nearest', 1, True),
##                            ('nn', mymath.distCanberra, 3)),
##                           ('5% NU; NN Canberra k=3',
##                            (1.0,'nearest', 5, False),
##                            (0.0,'nearest', 5, False),
##                            ('nn', mymath.distCanberra, 3)),
##                           ('5%  U; NN Canberra k=3',
##                            (1.0,'nearest', 5, True),
##                            (0.0,'nearest', 5, True),
##                            ('nn', mymath.distCanberra, 3)),
##                           ('10% NU; NN Canberra k=3',
##                            (1.0,'nearest', 10, False),
##                            (0.0,'nearest', 10, False),
##                            ('nn', mymath.distCanberra, 3)),
##                           ('10%  U; NN Canberra k=3',
##                            (1.0,'nearest', 10, True),
##                            (0.0,'nearest', 10, True),
##                            ('nn', mymath.distCanberra, 3)),
#
##                           ('1% NU; NN Manhatten k=9',
##                            (1.0,'nearest', 1, False),
##                            (0.0,'nearest', 1, False),
##                            ('nn', mymath.distL1, 9)),
##                           ('1%  U; NN Manhatten k=9',
##                            (1.0,'nearest', 1, True),
##                            (0.0,'nearest', 1, True),
##                            ('nn', mymath.distL1, 9)),
#                           ('5% NU; NN Manhatten k=9',
#                            (1.0,'nearest', 5, False),
#                            (0.0,'nearest', 5, False),
#                            ('nn', mymath.distL1, 9)),
#                           ('5%  U; NN Manhatten k=9',
#                            (1.0,'nearest', 5, True),
#                            (0.0,'nearest', 5, True),
#                            ('nn', mymath.distL1, 9)),
#                           ('10% NU; NN Manhatten k=9',
#                            (1.0,'nearest', 10, False),
#                            (0.0,'nearest', 10, False),
#                            ('nn', mymath.distL1, 9)),
#                           ('10%  U; NN Manhatten k=9',
#                            (1.0,'nearest', 10, True),
#                            (0.0,'nearest', 10, True),
#                            ('nn', mymath.distL1, 9)),
#
##                           ('1% NU; NN Euclidean k=9',
##                            (1.0,'nearest', 1, False),
##                            (0.0,'nearest', 1, False),
##                            ('nn', mymath.distL2, 9)),
##                           ('1%  U; NN Euclidean k=9',
##                            (1.0,'nearest', 1, True),
##                            (0.0,'nearest', 1, True),
##                            ('nn', mymath.distL2, 9)),
#                           ('5% NU; NN Euclidean k=9',
#                            (1.0,'nearest', 5, False),
#                            (0.0,'nearest', 5, False),
#                            ('nn', mymath.distL2, 9)),
#                           ('5%  U; NN Euclidean k=9',
#                            (1.0,'nearest', 5, True),
#                            (0.0,'nearest', 5, True),
#                            ('nn', mymath.distL2, 9)),
#                           ('10% NU; NN Euclidean k=9',
#                            (1.0,'nearest', 10, False),
#                            (0.0,'nearest', 10, False),
#                            ('nn', mymath.distL2, 9)),
#                           ('10%  U; NN Euclidean k=9',
#                            (1.0,'nearest', 10, True),
#                            (0.0,'nearest', 10, True),
#                            ('nn', mymath.distL2, 9)),
#
##                           ('1% NU; NN L-Infinity k=9',
##                            (1.0,'nearest', 1, False),
##                            (0.0,'nearest', 1, False),
##                            ('nn', mymath.distLInf, 9)),
##                           ('1%  U; NN L-Infinity k=9',
##                            (1.0,'nearest', 1, True),
##                            (0.0,'nearest', 1, True),
##                            ('nn', mymath.distLInf, 9)),
##                           ('5% NU; NN L-Infinity k=9',
##                            (1.0,'nearest', 5, False),
##                            (0.0,'nearest', 5, False),
##                            ('nn', mymath.distLInf, 9)),
##                           ('5%  U; NN L-Infinity k=9',
##                            (1.0,'nearest', 5, True),
##                            (0.0,'nearest', 5, True),
##                            ('nn', mymath.distLInf, 9)),
##                           ('10% NU; NN L-Infinity k=9',
##                            (1.0,'nearest', 10, False),
##                            (0.0,'nearest', 10, False),
##                            ('nn', mymath.distLInf, 9)),
##                           ('10%  U; NN L-Infinity k=9',
##                            (1.0,'nearest', 10, True),
##                            (0.0,'nearest', 10, True),
##                            ('nn', mymath.distLInf, 9)),
#
##                            ('1% NU; NN Cosine k=9',
##                             (1.0,'nearest', 1, False),
##                             (0.0,'nearest', 1, False),
##                             ('nn', mymath.distCosine, 9)),
##                            ('1%  U; NN Cosine k=9',
##                             (1.0,'nearest', 1, True),
##                             (0.0,'nearest', 1, True),
##                             ('nn', mymath.distCosine, 9)),
##                            ('5% NU; NN Cosine k=9',
##                             (1.0,'nearest', 5, False),
##                             (0.0,'nearest', 5, False),
##                             ('nn', mymath.distCosine, 9)),
##                            ('5%  U; NN Cosine k=9',
##                             (1.0,'nearest', 5, True),
##                             (0.0,'nearest', 5, True),
##                             ('nn', mymath.distCosine, 9)),
##                            ('10% NU; NN Cosine k=9',
##                             (1.0,'nearest', 10, False),
##                             (0.0,'nearest', 10, False),
##                             ('nn', mymath.distCosine, 9)),
##                            ('10%  U; NN Cosine k=9',
##                             (1.0,'nearest', 10, True),
##                             (0.0,'nearest', 10, True),
##                             ('nn', mymath.distCosine, 9)),
#
##                           ('1% NU; NN Canberra k=9',
##                            (1.0,'nearest', 1, False),
##                            (0.0,'nearest', 1, False),
##                            ('nn', mymath.distCanberra, 9)),
##                           ('1%  U; NN Canberra k=9',
##                            (1.0,'nearest', 1, True),
##                            (0.0,'nearest', 1, True),
##                            ('nn', mymath.distCanberra, 9)),
##                           ('5% NU; NN Canberra k=9',
##                            (1.0,'nearest', 5, False),
##                            (0.0,'nearest', 5, False),
##                            ('nn', mymath.distCanberra, 9)),
##                           ('5%  U; NN Canberra k=9',
##                            (1.0,'nearest', 5, True),
##                            (0.0,'nearest', 5, True),
##                            ('nn', mymath.distCanberra, 9)),
##                           ('10% NU; NN Canberra k=9',
##                            (1.0,'nearest', 10, False),
##                            (0.0,'nearest', 10, False),
##                            ('nn', mymath.distCanberra, 9)),
##                           ('10%  U; NN Canberra k=9',
##                            (1.0,'nearest', 10, True),
##                            (0.0,'nearest', 10, True),
##                            ('nn', mymath.distCanberra, 9)),


#                           ('1% NU; K-means Manhatten 0',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('kmeans', mymath.distL1, 0)),
#                           ('1%  U; K-means Manhatten 0',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('kmeans', mymath.distL1, 0)),
#                           ('5% NU; K-means Manhatten 0',
#                            (1.0,'nearest', 5, False),
#                            (0.0,'nearest', 5, False),
#                            ('kmeans', mymath.distL1, 0)),
#                           ('5%  U; K-means Manhatten 0',
#                            (1.0,'nearest', 5, True),
#                            (0.0,'nearest', 5, True),
#                            ('kmeans', mymath.distL1, 0)),
#                           ('10% NU; K-means Manhatten 0',
#                            (1.0,'nearest', 10, False),
#                            (0.0,'nearest', 10, False),
#                            ('kmeans', mymath.distL1, 0)),
#                           ('10%  U; K-means Manhatten 0',
#                            (1.0,'nearest', 10, True),
#                            (0.0,'nearest', 10, True),
#                            ('kmeans', mymath.distL1, 0)),

#                           ('1% NU; K-means Euclidean 0',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('kmeans', mymath.distL2, 0)),
#                           ('1%  U; K-means Euclidean 0',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('kmeans', mymath.distL2, 0)),
#                           ('5% NU; K-means Euclidean 0',
#                            (1.0,'nearest', 5, False),
#                            (0.0,'nearest', 5, False),
#                            ('kmeans', mymath.distL2, 0)),
#                           ('5%  U; K-means Euclidean 0',
#                            (1.0,'nearest', 5, True),
#                            (0.0,'nearest', 5, True),
#                            ('kmeans', mymath.distL2, 0)),
#                           ('10% NU; K-means Euclidean 0',
#                            (1.0,'nearest', 10, False),
#                            (0.0,'nearest', 10, False),
#                            ('kmeans', mymath.distL2, 0)),
#                           ('10%  U; K-means Euclidean 0',
#                            (1.0,'nearest', 10, True),
#                            (0.0,'nearest', 10, True),
#                            ('kmeans', mymath.distL2, 0)),

#                           ('1% NU; K-means L-Infinity 0',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('kmeans', mymath.distLInf, 0)),
#                           ('1%  U; K-means L-Infinity 0',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('kmeans', mymath.distLInf, 0)),
#                           ('5% NU; K-means L-Infinity 0',
#                            (1.0,'nearest', 5, False),
#                            (0.0,'nearest', 5, False),
#                            ('kmeans', mymath.distLInf, 0)),
#                           ('5%  U; K-means L-Infinity 0',
#                            (1.0,'nearest', 5, True),
#                            (0.0,'nearest', 5, True),
#                            ('kmeans', mymath.distLInf, 0)),
#                           ('10% NU; K-means L-Infinity 0',
#                            (1.0,'nearest', 10, False),
#                            (0.0,'nearest', 10, False),
#                            ('kmeans', mymath.distLInf, 0)),
#                           ('10%  U; K-means L-Infinity 0',
#                            (1.0,'nearest', 10, True),
#                            (0.0,'nearest', 10, True),
#                            ('kmeans', mymath.distLInf, 0)),

#                            ('1% NU; K-means Cosine 0',
#                             (1.0,'nearest', 1, False),
#                             (0.0,'nearest', 1, False),
#                             ('kmeans', mymath.distCosine, 0)),
#                            ('1%  U; K-means Cosine 0',
#                             (1.0,'nearest', 1, True),
#                             (0.0,'nearest', 1, True),
#                             ('kmeans', mymath.distCosine, 0)),
#                            ('5% NU; K-means Cosine 0',
#                             (1.0,'nearest', 5, False),
#                             (0.0,'nearest', 5, False),
#                             ('kmeans', mymath.distCosine, 0)),
#                            ('5%  U; K-means Cosine 0',
#                             (1.0,'nearest', 5, True),
#                             (0.0,'nearest', 5, True),
#                             ('kmeans', mymath.distCosine, 0)),
#                            ('10% NU; K-means Cosine 0',
#                             (1.0,'nearest', 10, False),
#                             (0.0,'nearest', 10, False),
#                             ('kmeans', mymath.distCosine, 0)),
#                            ('10%  U; K-means Cosine 0',
#                             (1.0,'nearest', 10, True),
#                             (0.0,'nearest', 10, True),
#                             ('kmeans', mymath.distCosine, 0)),

#                           ('1% NU; K-means Canberra 0',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('kmeans', mymath.distCanberra, 0)),
#                           ('1%  U; K-means Canberra 0',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('kmeans', mymath.distCanberra, 0)),
#                           ('5% NU; K-means Canberra 0',
#                            (1.0,'nearest', 5, False),
#                            (0.0,'nearest', 5, False),
#                            ('kmeans', mymath.distCanberra, 0)),
#                           ('5%  U; K-means Canberra 0',
#                            (1.0,'nearest', 5, True),
#                            (0.0,'nearest', 5, True),
#                            ('kmeans', mymath.distCanberra, 0)),
#                           ('10% NU; K-means Canberra 0',
#                            (1.0,'nearest', 10, False),
#                            (0.0,'nearest', 10, False),
#                            ('kmeans', mymath.distCanberra, 0)),
#                           ('10%  U; K-means Canberra 0',
#                            (1.0,'nearest', 10, True),
#                            (0.0,'nearest', 10, True),
#                            ('kmeans', mymath.distCanberra, 0)),

#                           ('1% NU; K-means Manhatten',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('kmeans', mymath.distL1, 10000)),
#                           ('1%  U; K-means Manhatten',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('kmeans', mymath.distL1, 10000)),
#                           ('5% NU; K-means Manhatten',
#                            (1.0,'nearest', 5, False),
#                            (0.0,'nearest', 5, False),
#                            ('kmeans', mymath.distL1, 10000)),
#                           ('5%  U; K-means Manhatten',
#                            (1.0,'nearest', 5, True),
#                            (0.0,'nearest', 5, True),
#                            ('kmeans', mymath.distL1, 10000)),
#                           ('10% NU; K-means Manhatten',
#                            (1.0,'nearest', 10, False),
#                            (0.0,'nearest', 10, False),
#                            ('kmeans', mymath.distL1, 10000)),
#                           ('10%  U; K-means Manhatten',
#                            (1.0,'nearest', 10, True),
#                            (0.0,'nearest', 10, True),
#                            ('kmeans', mymath.distL1, 10000)),

#                           ('1% NU; K-means Euclidean',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('kmeans', mymath.distL2, 10000)),
#                           ('1%  U; K-means Euclidean',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('kmeans', mymath.distL2, 10000)),
#                           ('5% NU; K-means Euclidean',
#                            (1.0,'nearest', 5, False),
#                            (0.0,'nearest', 5, False),
#                            ('kmeans', mymath.distL2, 10000)),
#                           ('5%  U; K-means Euclidean',
#                            (1.0,'nearest', 5, True),
#                            (0.0,'nearest', 5, True),
#                            ('kmeans', mymath.distL2, 10000)),
#                           ('10% NU; K-means Euclidean',
#                            (1.0,'nearest', 10, False),
#                            (0.0,'nearest', 10, False),
#                            ('kmeans', mymath.distL2, 10000)),
#                           ('10%  U; K-means Euclidean',
#                            (1.0,'nearest', 10, True),
#                            (0.0,'nearest', 10, True),
#                            ('kmeans', mymath.distL2, 10000)),

#                           ('1% NU; K-means L-Infinity',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('kmeans', mymath.distLInf, 10000)),
#                           ('1%  U; K-means L-Infinity',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('kmeans', mymath.distLInf, 10000)),
#                           ('5% NU; K-means L-Infinity',
#                            (1.0,'nearest', 5, False),
#                            (0.0,'nearest', 5, False),
#                            ('kmeans', mymath.distLInf, 10000)),
#                           ('5%  U; K-means L-Infinity',
#                            (1.0,'nearest', 5, True),
#                            (0.0,'nearest', 5, True),
#                            ('kmeans', mymath.distLInf, 10000)),
#                           ('10% NU; K-means L-Infinity',
#                            (1.0,'nearest', 10, False),
#                            (0.0,'nearest', 10, False),
#                            ('kmeans', mymath.distLInf, 10000)),
#                           ('10%  U; K-means L-Infinity',
#                            (1.0,'nearest', 10, True),
#                            (0.0,'nearest', 10, True),
#                            ('kmeans', mymath.distLInf, 10000)),

#                            ('1% NU; K-means Cosine',
#                             (1.0,'nearest', 1, False),
#                             (0.0,'nearest', 1, False),
#                             ('kmeans', mymath.distCosine, 10000)),
#                            ('1%  U; K-means Cosine',
#                             (1.0,'nearest', 1, True),
#                             (0.0,'nearest', 1, True),
#                             ('kmeans', mymath.distCosine, 10000)),
#                            ('5% NU; K-means Cosine',
#                             (1.0,'nearest', 5, False),
#                             (0.0,'nearest', 5, False),
#                             ('kmeans', mymath.distCosine, 10000)),
#                            ('5%  U; K-means Cosine',
#                             (1.0,'nearest', 5, True),
#                             (0.0,'nearest', 5, True),
#                             ('kmeans', mymath.distCosine, 10000)),
#                            ('10% NU; K-means Cosine',
#                             (1.0,'nearest', 10, False),
#                             (0.0,'nearest', 10, False),
#                             ('kmeans', mymath.distCosine, 10000)),
#                            ('10%  U; K-means Cosine',
#                             (1.0,'nearest', 10, True),
#                             (0.0,'nearest', 10, True),
#                             ('kmeans', mymath.distCosine, 10000)),

#                           ('1% NU; K-means Canberra',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('kmeans', mymath.distCanberra, 10000)),
#                           ('1%  U; K-means Canberra',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('kmeans', mymath.distCanberra, 10000)),
#                           ('5% NU; K-means Canberra',
#                            (1.0,'nearest', 5, False),
#                            (0.0,'nearest', 5, False),
#                            ('kmeans', mymath.distCanberra, 10000)),
#                           ('5%  U; K-means Canberra',
#                            (1.0,'nearest', 5, True),
#                            (0.0,'nearest', 5, True),
#                            ('kmeans', mymath.distCanberra, 10000)),
#                           ('10% NU; K-means Canberra',
#                            (1.0,'nearest', 10, False),
#                            (0.0,'nearest', 10, False),
#                            ('kmeans', mymath.distCanberra, 10000)),
#                           ('10%  U; K-means Canberra',
#                            (1.0,'nearest', 10, True),
#                            (0.0,'nearest', 10, True),
#                            ('kmeans', mymath.distCanberra, 10000)),


#                           ('1% NU; SVM linear C=10 0-0',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'LINEAR', 10, 0, 0)),
#                           ('1%  U; SVM linear C=10 0-0',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'LINEAR', 10, 0, 0)),
                           ('5% NU; SVM linear C=10 0-0',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'LINEAR', 10, 0, 0)),
                           ('5%  U; SVM linear C=10 0-0',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'LINEAR', 10, 0, 0)),
                           ('10% NU; SVM linear C=10 0-0',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'LINEAR', 10, 0, 0)),
                           ('10%  U; SVM linear C=10 0-0',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'LINEAR', 10, 0, 0)),

#                           ('1% NU; SVM linear C= 1 0-0',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'LINEAR', 1, 0, 0)),
#                           ('1%  U; SVM linear C= 1 0-0',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'LINEAR', 1, 0, 0)),
                           ('5% NU; SVM linear C= 1 0-0',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'LINEAR', 1, 0, 0)),
                           ('5%  U; SVM linear C= 1 0-0',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'LINEAR', 1, 0, 0)),
                           ('10% NU; SVM linear C= 1 0-0',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'LINEAR', 1, 0, 0)),
                           ('10%  U; SVM linear C= 1 0-0',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'LINEAR', 1, 0, 0)),

#                           ('1% NU; SVM linear C=0.1 0-0',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'LINEAR', 0.1, 0, 0)),
#                           ('1%  U; SVM linear C=0.1 0-0',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'LINEAR', 0.1, 0, 0)),
                           ('5% NU; SVM linear C=0.1 0-0',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'LINEAR', 0.1, 0, 0)),
                           ('5%  U; SVM linear C=0.1 0-0',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'LINEAR', 0.1, 0, 0)),
                           ('10% NU; SVM linear C=0.1 0-0',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'LINEAR', 0.1, 0, 0)),
                           ('10%  U; SVM linear C=0.1 0-0',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'LINEAR', 0.1, 0, 0)),


#                           ('1% NU; SVM poly C=10 0-0',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'POLY', 10, 0, 0)),
#                           ('1%  U; SVM poly C=10 0-0',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'POLY', 10, 0, 0)),
                           ('5% NU; SVM poly C=10 0-0',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'POLY', 10, 0, 0)),
                           ('5%  U; SVM poly C=10 0-0',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'POLY', 10, 0, 0)),
                           ('10% NU; SVM poly C=10 0-0',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'POLY', 10, 0, 0)),
                           ('10%  U; SVM poly C=10 0-0',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'POLY', 10, 0, 0)),

#                           ('1% NU; SVM poly C= 1 0-0',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'POLY', 1, 0, 0)),
#                           ('1%  U; SVM poly C= 1 0-0',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'POLY', 1, 0, 0)),
                           ('5% NU; SVM poly C= 1 0-0',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'POLY', 1, 0, 0)),
                           ('5%  U; SVM poly C= 1 0-0',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'POLY', 1, 0, 0)),
                           ('10% NU; SVM poly C= 1 0-0',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'POLY', 1, 0, 0)),
                           ('10%  U; SVM poly C= 1 0-0',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'POLY', 1, 0, 0)),

#                           ('1% NU; SVM poly C=0.1 0-0',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'POLY', 0.1, 0, 0)),
#                           ('1%  U; SVM poly C=0.1 0-0',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'POLY', 0.1, 0, 0)),
                           ('5% NU; SVM poly C=0.1 0-0',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'POLY', 0.1, 0, 0)),
                           ('5%  U; SVM poly C=0.1 0-0',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'POLY', 0.1, 0, 0)),
                           ('10% NU; SVM poly C=0.1 0-0',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'POLY', 0.1, 0, 0)),
                           ('10%  U; SVM poly C=0.1 0-0',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'POLY', 0.1, 0, 0)),


#                           ('1% NU; SVM RBF C=10 0-0',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'RBF', 10, 0, 0)),
#                           ('1%  U; SVM RBF C=10 0-0',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'RBF', 10, 0, 0)),
                           ('5% NU; SVM RBF C=10 0-0',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'RBF', 10, 0, 0)),
                           ('5%  U; SVM RBF C=10 0-0',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'RBF', 10, 0, 0)),
                           ('10% NU; SVM RBF C=10 0-0',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'RBF', 10, 0, 0)),
                           ('10%  U; SVM RBF C=10 0-0',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'RBF', 10, 0, 0)),

#                           ('1% NU; SVM RBF C= 1 0-0',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'RBF', 1, 0, 0)),
#                           ('1%  U; SVM RBF C= 1 0-0',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'RBF', 1, 0, 0)),
                           ('5% NU; SVM RBF C= 1 0-0',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'RBF', 1, 0, 0)),
                           ('5%  U; SVM RBF C= 1 0-0',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'RBF', 1, 0, 0)),
                           ('10% NU; SVM RBF C= 1 0-0',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'RBF', 1, 0, 0)),
                           ('10%  U; SVM RBF C= 1 0-0',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'RBF', 1, 0, 0)),

#                           ('1% NU; SVM RBF C=0.1 0-0',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'RBF', 0.1, 0, 0)),
#                           ('1%  U; SVM RBF C=0.1 0-0',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'RBF', 0.1, 0, 0)),
                           ('5% NU; SVM RBF C=0.1 0-0',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'RBF', 0.1, 0, 0)),
                           ('5%  U; SVM RBF C=0.1 0-0',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'RBF', 0.1, 0, 0)),
                           ('10% NU; SVM RBF C=0.1 0-0',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'RBF', 0.1, 0, 0)),
                           ('10%  U; SVM RBF C=0.1 0-0',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'RBF', 0.1, 0, 0)),


#                           ('1% NU; SVM linear C=10 25-25',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'LINEAR', 10, 25,25)),
#                           ('1%  U; SVM linear C=10 25-25',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'LINEAR', 10, 25,25)),
                           ('5% NU; SVM linear C=10 25-25',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'LINEAR', 10, 25,25)),
                           ('5%  U; SVM linear C=10 25-25',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'LINEAR', 10, 25,25)),
                           ('10% NU; SVM linear C=10 25-25',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'LINEAR', 10, 25,25)),
                           ('10%  U; SVM linear C=10 25-25',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'LINEAR', 10, 25,25)),

#                           ('1% NU; SVM linear C= 1 25-25',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'LINEAR', 1, 25,25)),
#                           ('1%  U; SVM linear C= 1 25-25',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'LINEAR', 1, 25,25)),
                           ('5% NU; SVM linear C= 1 25-25',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'LINEAR', 1, 25,25)),
                           ('5%  U; SVM linear C= 1 25-25',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'LINEAR', 1, 25,25)),
                           ('10% NU; SVM linear C= 1 25-25',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'LINEAR', 1, 25,25)),
                           ('10%  U; SVM linear C= 1 25-25',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'LINEAR', 1, 25,25)),

#                           ('1% NU; SVM linear C=0.1 25-25',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'LINEAR', 0.1, 25,25)),
#                           ('1%  U; SVM linear C=0.1 25-25',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'LINEAR', 0.1, 25,25)),
                           ('5% NU; SVM linear C=0.1 25-25',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'LINEAR', 0.1, 25,25)),
                           ('5%  U; SVM linear C=0.1 25-25',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'LINEAR', 0.1, 25,25)),
                           ('10% NU; SVM linear C=0.1 25-25',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'LINEAR', 0.1, 25,25)),
                           ('10%  U; SVM linear C=0.1 25-25',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'LINEAR', 0.1, 25,25)),


#                           ('1% NU; SVM poly C=10 25-25',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'POLY', 10, 25,25)),
#                           ('1%  U; SVM poly C=10 25-25',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'POLY', 10, 25,25)),
                           ('5% NU; SVM poly C=10 25-25',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'POLY', 10, 25,25)),
                           ('5%  U; SVM poly C=10 25-25',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'POLY', 10, 25,25)),
                           ('10% NU; SVM poly C=10 25-25',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'POLY', 10, 25,25)),
                           ('10%  U; SVM poly C=10 25-25',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'POLY', 10, 25,25)),

#                           ('1% NU; SVM poly C= 1 25-25',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'POLY', 1, 25,25)),
#                           ('1%  U; SVM poly C= 1 25-25',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'POLY', 1, 25,25)),
                           ('5% NU; SVM poly C= 1 25-25',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'POLY', 1, 25,25)),
                           ('5%  U; SVM poly C= 1 25-25',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'POLY', 1, 25,25)),
                           ('10% NU; SVM poly C= 1 25-25',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'POLY', 1, 25,25)),
                           ('10%  U; SVM poly C= 1 25-25',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'POLY', 1, 25,25)),

#                           ('1% NU; SVM poly C=0.1 25-25',
#                            (1.0,'nearest', 1, False),
#                           (0.0,'nearest', 1, False),
#                            ('svm', 'POLY', 0.1, 25,25)),
#                           ('1%  U; SVM poly C=0.1 25-25',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'POLY', 0.1, 25,25)),
                           ('5% NU; SVM poly C=0.1 25-25',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'POLY', 0.1, 25,25)),
                           ('5%  U; SVM poly C=0.1 25-25',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'POLY', 0.1, 25,25)),
                           ('10% NU; SVM poly C=0.1 25-25',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'POLY', 0.1, 25,25)),
                           ('10%  U; SVM poly C=0.1 25-25',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'POLY', 0.1, 25,25)),


#                           ('1% NU; SVM RBF C=10 25-25',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'RBF', 10, 25,25)),
#                           ('1%  U; SVM RBF C=10 25-25',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'RBF', 10, 25,25)),
                           ('5% NU; SVM RBF C=10 25-25',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'RBF', 10, 25,25)),
                           ('5%  U; SVM RBF C=10 25-25',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'RBF', 10, 25,25)),
                           ('10% NU; SVM RBF C=10 25-25',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'RBF', 10, 25,25)),
                           ('10%  U; SVM RBF C=10 25-25',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'RBF', 10, 25,25)),

#                           ('1% NU; SVM RBF C= 1 25-25',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'RBF', 1, 25,25)),
#                           ('1%  U; SVM RBF C= 1 25-25',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'RBF', 1, 25,25)),
                           ('5% NU; SVM RBF C= 1 25-25',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'RBF', 1, 25,25)),
                           ('5%  U; SVM RBF C= 1 25-25',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'RBF', 1, 25,25)),
                           ('10% NU; SVM RBF C= 1 25-25',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'RBF', 1, 25,25)),
                           ('10%  U; SVM RBF C= 1 25-25',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'RBF', 1, 25,25)),

#                           ('1% NU; SVM RBF C=0.1 25-25',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'RBF', 0.1, 25,25)),
#                           ('1%  U; SVM RBF C=0.1 25-25',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'RBF', 0.1, 25,25)),
                           ('5% NU; SVM RBF C=0.1 25-25',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'RBF', 0.1, 25,25)),
                           ('5%  U; SVM RBF C=0.1 25-25',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'RBF', 0.1, 25,25)),
                           ('10% NU; SVM RBF C=0.1 25-25',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'RBF', 0.1, 25,25)),
                           ('10%  U; SVM RBF C=0.1 25-25',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'RBF', 0.1, 25,25)),


#                           ('1% NU; SVM linear C=10 25-50',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'LINEAR', 10, 25,50)),
#                           ('1%  U; SVM linear C=10 25-50',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'LINEAR', 10, 25,50)),
                           ('5% NU; SVM linear C=10 25-50',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'LINEAR', 10, 25,50)),
                           ('5%  U; SVM linear C=10 25-50',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'LINEAR', 10, 25,50)),
                           ('10% NU; SVM linear C=10 25-50',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'LINEAR', 10, 25,50)),
                           ('10%  U; SVM linear C=10 25-50',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'LINEAR', 10, 25,50)),

#                           ('1% NU; SVM linear C= 1 25-50',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'LINEAR', 1, 25,50)),
#                           ('1%  U; SVM linear C= 1 25-50',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'LINEAR', 1, 25,50)),
                           ('5% NU; SVM linear C= 1 25-50',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'LINEAR', 1, 25,50)),
                           ('5%  U; SVM linear C= 1 25-50',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'LINEAR', 1, 25,50)),
                           ('10% NU; SVM linear C= 1 25-50',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'LINEAR', 1, 25,50)),
                           ('10%  U; SVM linear C= 1 25-50',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'LINEAR', 1, 25,50)),

#                           ('1% NU; SVM linear C=0.1 25-50',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'LINEAR', 0.1, 25,50)),
#                           ('1%  U; SVM linear C=0.1 25-50',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'LINEAR', 0.1, 25,50)),
                           ('5% NU; SVM linear C=0.1 25-50',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'LINEAR', 0.1, 25,50)),
                           ('5%  U; SVM linear C=0.1 25-50',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'LINEAR', 0.1, 25,50)),
                           ('10% NU; SVM linear C=0.1 25-50',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'LINEAR', 0.1, 25,50)),
                           ('10%  U; SVM linear C=0.1 25-50',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'LINEAR', 0.1, 25,50)),


#                           ('1% NU; SVM poly C=10 25-50',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'POLY', 10, 25,50)),
#                           ('1%  U; SVM poly C=10 25-50',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'POLY', 10, 25,50)),
                           ('5% NU; SVM poly C=10 25-50',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'POLY', 10, 25,50)),
                           ('5%  U; SVM poly C=10 25-50',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'POLY', 10, 25,50)),
                           ('10% NU; SVM poly C=10 25-50',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'POLY', 10, 25,50)),
                           ('10%  U; SVM poly C=10 25-50',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'POLY', 10, 25,50)),

#                           ('1% NU; SVM poly C= 1 25-50',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'POLY', 1, 25,50)),
#                           ('1%  U; SVM poly C= 1 25-50',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'POLY', 1, 25,50)),
                           ('5% NU; SVM poly C= 1 25-50',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'POLY', 1, 25,50)),
                           ('5%  U; SVM poly C= 1 25-50',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'POLY', 1, 25,50)),
                           ('10% NU; SVM poly C= 1 25-50',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'POLY', 1, 25,50)),
                           ('10%  U; SVM poly C= 1 25-50',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'POLY', 1, 25,50)),

#                           ('1% NU; SVM poly C=0.1 25-50',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'POLY', 0.1, 25,50)),
#                           ('1%  U; SVM poly C=0.1 25-50',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'POLY', 0.1, 25,50)),
                           ('5% NU; SVM poly C=0.1 25-50',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'POLY', 0.1, 25,50)),
                           ('5%  U; SVM poly C=0.1 25-50',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'POLY', 0.1, 25,50)),
                           ('10% NU; SVM poly C=0.1 25-50',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'POLY', 0.1, 25,50)),
                           ('10%  U; SVM poly C=0.1 25-50',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'POLY', 0.1, 25,50)),


#                           ('1% NU; SVM RBF C=10 25-50',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'RBF', 10, 25,50)),
#                           ('1%  U; SVM RBF C=10 25-50',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'RBF', 10, 25,50)),
                           ('5% NU; SVM RBF C=10 25-50',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'RBF', 10, 25,50)),
                           ('5%  U; SVM RBF C=10 25-50',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'RBF', 10, 25,50)),
                           ('10% NU; SVM RBF C=10 25-50',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'RBF', 10, 25,50)),
                           ('10%  U; SVM RBF C=10 25-50',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'RBF', 10, 25,50)),

#                           ('1% NU; SVM RBF C= 1 25-50',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'RBF', 1, 25,50)),
#                           ('1%  U; SVM RBF C= 1 25-50',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'RBF', 1, 25,50)),
                           ('5% NU; SVM RBF C= 1 25-50',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'RBF', 1, 25,50)),
                           ('5%  U; SVM RBF C= 1 25-50',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'RBF', 1, 25,50)),
                           ('10% NU; SVM RBF C= 1 25-50',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'RBF', 1, 25,50)),
                           ('10%  U; SVM RBF C= 1 25-50',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'RBF', 1, 25,50)),

#                           ('1% NU; SVM RBF C=0.1 25-50',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'RBF', 0.1, 25,50)),
#                           ('1%  U; SVM RBF C=0.1 25-50',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'RBF', 0.1, 25,50)),
                           ('5% NU; SVM RBF C=0.1 25-50',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'RBF', 0.1, 25,50)),
                           ('5%  U; SVM RBF C=0.1 25-50',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'RBF', 0.1, 25,50)),
                           ('10% NU; SVM RBF C=0.1 25-50',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'RBF', 0.1, 25,50)),
                           ('10%  U; SVM RBF C=0.1 25-50',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'RBF', 0.1, 25,50)),


#                           ('1% NU; SVM linear C=10 50-100',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'LINEAR', 10, 50,100)),
#                           ('1%  U; SVM linear C=10 50-100',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'LINEAR', 10, 50,100)),
                           ('5% NU; SVM linear C=10 50-100',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'LINEAR', 10, 50,100)),
                           ('5%  U; SVM linear C=10 50-100',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'LINEAR', 10, 50,100)),
                           ('10% NU; SVM linear C=10 50-100',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'LINEAR', 10, 50,100)),
                           ('10%  U; SVM linear C=10 50-100',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'LINEAR', 10, 50,100)),

#                           ('1% NU; SVM linear C= 1 50-100',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'LINEAR', 1, 50,100)),
#                           ('1%  U; SVM linear C= 1 50-100',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'LINEAR', 1, 50,100)),
                           ('5% NU; SVM linear C= 1 50-100',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'LINEAR', 1, 50,100)),
                           ('5%  U; SVM linear C= 1 50-100',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'LINEAR', 1, 50,100)),
                           ('10% NU; SVM linear C= 1 50-100',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'LINEAR', 1, 50,100)),
                           ('10%  U; SVM linear C= 1 50-100',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'LINEAR', 1, 50,100)),

#                           ('1% NU; SVM linear C=0.1 50-100',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'LINEAR', 0.1, 50,100)),
#                           ('1%  U; SVM linear C=0.1 50-100',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'LINEAR', 0.1, 50,100)),
                           ('5% NU; SVM linear C=0.1 50-100',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'LINEAR', 0.1, 50,100)),
                           ('5%  U; SVM linear C=0.1 50-100',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'LINEAR', 0.1, 50,100)),
                           ('10% NU; SVM linear C=0.1 50-100',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'LINEAR', 0.1, 50,100)),
                           ('10%  U; SVM linear C=0.1 50-100',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'LINEAR', 0.1, 50,100)),


#                           ('1% NU; SVM poly C=10 50-100',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'POLY', 10, 50,100)),
#                           ('1%  U; SVM poly C=10 50-100',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'POLY', 10, 50,100)),
                           ('5% NU; SVM poly C=10 50-100',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'POLY', 10, 50,100)),
                           ('5%  U; SVM poly C=10 50-100',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'POLY', 10, 50,100)),
                           ('10% NU; SVM poly C=10 50-100',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'POLY', 10, 50,100)),
                           ('10%  U; SVM poly C=10 50-100',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'POLY', 10, 50,100)),

#                           ('1% NU; SVM poly C= 1 50-100',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'POLY', 1, 50,100)),
#                           ('1%  U; SVM poly C= 1 50-100',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'POLY', 1, 50,100)),
                           ('5% NU; SVM poly C= 1 50-100',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'POLY', 1, 50,100)),
                           ('5%  U; SVM poly C= 1 50-100',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'POLY', 1, 50,100)),
                           ('10% NU; SVM poly C= 1 50-100',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'POLY', 1, 50,100)),
                           ('10%  U; SVM poly C= 1 50-100',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'POLY', 1, 50,100)),

#                           ('1% NU; SVM poly C=0.1 50-100',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'POLY', 0.1, 50,100)),
#                           ('1%  U; SVM poly C=0.1 50-100',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'POLY', 0.1, 50,100)),
                           ('5% NU; SVM poly C=0.1 50-100',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'POLY', 0.1, 50,100)),
                           ('5%  U; SVM poly C=0.1 50-100',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'POLY', 0.1, 50,100)),
                           ('10% NU; SVM poly C=0.1 50-100',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'POLY', 0.1, 50,100)),
                           ('10%  U; SVM poly C=0.1 50-100',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'POLY', 0.1, 50,100)),


#                           ('1% NU; SVM RBF C=10 50-100',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'RBF', 10, 50,100)),
#                           ('1%  U; SVM RBF C=10 50-100',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'RBF', 10, 50,100)),
                           ('5% NU; SVM RBF C=10 50-100',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'RBF', 10, 50,100)),
                           ('5%  U; SVM RBF C=10 50-100',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'RBF', 10, 50,100)),
                           ('10% NU; SVM RBF C=10 50-100',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'RBF', 10, 50,100)),
                           ('10%  U; SVM RBF C=10 50-100',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'RBF', 10, 50,100)),

#                           ('1% NU; SVM RBF C= 1 50-100',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'RBF', 1, 50,100)),
#                           ('1%  U; SVM RBF C= 1 50-100',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'RBF', 1, 50,100)),
                           ('5% NU; SVM RBF C= 1 50-100',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'RBF', 1, 50,100)),
                           ('5%  U; SVM RBF C= 1 50-100',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'RBF', 1, 50,100)),
                           ('10% NU; SVM RBF C= 1 50-100',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'RBF', 1, 50,100)),
                           ('10%  U; SVM RBF C= 1 50-100',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'RBF', 1, 50,100)),

#                           ('1% NU; SVM RBF C=0.1 50-100',
#                            (1.0,'nearest', 1, False),
#                            (0.0,'nearest', 1, False),
#                            ('svm', 'RBF', 0.1, 50,100)),
#                           ('1%  U; SVM RBF C=0.1 50-100',
#                            (1.0,'nearest', 1, True),
#                            (0.0,'nearest', 1, True),
#                            ('svm', 'RBF', 0.1, 50,100)),
                           ('5% NU; SVM RBF C=0.1 50-100',
                            (1.0,'nearest', 5, False),
                            (0.0,'nearest', 5, False),
                            ('svm', 'RBF', 0.1, 50,100)),
                           ('5%  U; SVM RBF C=0.1 50-100',
                            (1.0,'nearest', 5, True),
                            (0.0,'nearest', 5, True),
                            ('svm', 'RBF', 0.1, 50,100)),
                           ('10% NU; SVM RBF C=0.1 50-100',
                            (1.0,'nearest', 10, False),
                            (0.0,'nearest', 10, False),
                            ('svm', 'RBF', 0.1, 50,100)),
                           ('10%  U; SVM RBF C=0.1 50-100',
                            (1.0,'nearest', 10, True),
                            (0.0,'nearest', 10, True),
                            ('svm', 'RBF', 0.1, 50,100))]


    if (do_two_step_near == False):
      two_step_param_list = []  # Don't perform these experiments

    two_step_res_list = []

    for (class_name, s1_m, s1_nm, s2_class) in two_step_param_list:

      # Generate number of examples to select from percentage numbers
      #
      perc_val = s1_m[2]
      sel_num = int(num_w_vec * perc_val/100.0)  # Select same number for
                                                 # matches and non-matches

      # Loop over possible random selection methods
      # (for each: string to add to output, random selection tuple, number of
      # iterations to do)
      #
      for rand_sel in [('; ',       None,                     1),
#                       ('; RS-U1',  ('uniform',      1,  1), 10),
#                       ('; RS-U10', ('uniform',     10, 10), 10),
#                       ('; RS-L1',  ('linear',       1,  1), 10),
#                       ('; RS-L10', ('linear',      10, 10), 10),
#                       ('; RS-E1',  ('exponential',  1,  1), 10),
#                       ('; RS-E10', ('exponential', 10, 10), 10)
                       ]:

        (rand_sel_str, rand_sel_tup, num_iter) = rand_sel  # Unpack details

#        # Balanced version first - - - - - - - - - - - - - - - - - - - - - - -
#        #
#        s1_m =  (s1_m[0],  s1_m[1],  sel_num, s1_m[3])
#        s1_nm = (s1_nm[0], s1_nm[1], sel_num, s1_nm[3])
#        bal_class_name = class_name.replace(';', '  B;')
#
#        acc_res_list =   []
#        prec_res_list =  []
#        reca_res_list =  []
#        fmeas_res_list = []
#        time_res_list =  []
#
#        for iter in range(num_iter):  # Loop over number of iterations - - - -
#
#          two_step_classifier = classification.TwoStep(s1_match_method = s1_m,
#                                                  s1_non_match_method = s1_nm,
#                                                  random_sel = rand_sel_tup,
#                                                  s2_classifier = s2_class)
#          start_time = time.time()
#          two_step_classifier.train(sel_w_vec_dict, true_m_set, true_nm_set)
#          res_sets = two_step_classifier.classify(sel_w_vec_dict)
#          two_step_time = time.time() - start_time
#          time_res_list.append(two_step_time)
#
#          test_results = two_step_classifier.cross_validate(sel_w_vec_dict,
#                                                            true_m_set,
#                                                            true_nm_set, n)
#
#          two_step_acc, two_step_prec, two_step_reca, two_step_fmeas = \
#                                                    get_measures(test_results)
#          acc_res_list.append(two_step_acc)
#          prec_res_list.append(two_step_prec)
#          reca_res_list.append(two_step_reca)
#          fmeas_res_list.append(two_step_fmeas)
#
#          del two_step_classifier
#
#        two_step_acc =   sum(acc_res_list) /   num_iter
#        two_step_prec =  sum(prec_res_list) /  num_iter
#        two_step_reca =  sum(reca_res_list) /  num_iter
#        two_step_fmeas = sum(fmeas_res_list) / num_iter
#        two_step_time = sum(time_res_list)   / num_iter
#
#        this_class_name = bal_class_name + rand_sel_str
#
#        two_step_res_list.append((this_class_name, two_step_acc,two_step_prec,
#                                two_step_reca, two_step_fmeas, two_step_time))
#        exp_res_list.append(two_step_res_list[-1])# Also append to all results
#
#        (this_class_name, acc, prec, reca,fmeas,rtime) = two_step_res_list[-1]
#        this_class_name = this_class_name.ljust(37)
#        res_file.write('    %37s| %5.3f | %5.3f | %5.3f | %5.3f | %9.3f' % \
#                       (this_class_name, acc, prec, reca, fmeas, rtime) + \
#                       os.linesep)
#        res_file.flush()

        # Im-balanced version next - - - - - - - - - - - - - - - - - - - - - -
        # Percentage numbers for matches according to the minimum number of
        # records in data sets (assuming data contains no duplicates, this is
        # the maximum number of matches possible)
        #
        max_num_matches = min(data_set_a.num_records, data_set_b.num_records)

        rec_w_vec_ratio = float(max_num_matches) / \
                          float(num_w_vec-max_num_matches)

        if (rec_w_vec_ratio < 1):  # More non-matches than matches
          m_sel_num = max(1, int(sel_num*rec_w_vec_ratio))
          nm_sel_num = sel_num
        else:  # More matches than non-matches
          m_sel_num = sel_num
          nm_sel_num = max(1, int(sel_num/rec_w_vec_ratio))

        s1_m =  (s1_m[0],  s1_m[1],  m_sel_num,  s1_m[3])
        s1_nm = (s1_nm[0], s1_nm[1], nm_sel_num, s1_nm[3])
        imbal_class_name = class_name.replace(';', ' IB;')

        acc_res_list =   []
        prec_res_list =  []
        reca_res_list =  []
        fmeas_res_list = []
        time_res_list =  []

        for iter in range(num_iter):  # Loop over number of iterations - - - -

          two_step_classifier = classification.TwoStep(s1_match_method = s1_m,
                                                   s1_non_match_method = s1_nm,
                                                   random_sel = rand_sel_tup,
                                                   s2_classifier = s2_class)
          start_time = time.time()
          two_step_classifier.train(sel_w_vec_dict, true_m_set, true_nm_set)
          res_sets = two_step_classifier.classify(sel_w_vec_dict)
          two_step_time = time.time() - start_time
          time_res_list.append(two_step_time)

          test_results = two_step_classifier.cross_validate(sel_w_vec_dict,
                                                            true_m_set,
                                                            true_nm_set, n)

          two_step_acc, two_step_prec, two_step_reca, two_step_fmeas = \
                                                     get_measures(test_results)
          acc_res_list.append(two_step_acc)
          prec_res_list.append(two_step_prec)
          reca_res_list.append(two_step_reca)
          fmeas_res_list.append(two_step_fmeas)

          del two_step_classifier

        two_step_acc =   sum(acc_res_list) /   num_iter
        two_step_prec =  sum(prec_res_list) /  num_iter
        two_step_reca =  sum(reca_res_list) /  num_iter
        two_step_fmeas = sum(fmeas_res_list) / num_iter
        two_step_time = sum(time_res_list)   / num_iter

        this_class_name = imbal_class_name + rand_sel_str

        two_step_res_list.append((this_class_name, two_step_acc, two_step_prec,
                                two_step_reca, two_step_fmeas, two_step_time))
        exp_res_list.append(two_step_res_list[-1]) # Also append to all results

        (this_class_name, acc, prec, reca, fmeas,rtime) = two_step_res_list[-1]
        this_class_name = this_class_name.ljust(37)
        res_file.write('    %37s| %5.3f | %5.3f | %5.3f | %5.3f | %9.3f' % \
                       (this_class_name, acc, prec, reca, fmeas, rtime) + \
                       os.linesep)
        res_file.flush()

    res_file.write('    ' + '-'*80 + os.linesep)
    res_file.flush()

    # -------------------------------------------------------------------------
    # Timing experiments only for selected classifiers
    #
    if (do_timing == True):

      # Write header of result table to file
      #
      res_file.write('    '+'='*80 + os.linesep)
      res_file.write('    Timing experiment' + os.linesep)
      res_file.write('    ' + '-'*80 + os.linesep)

      svm_classifier = classification.SuppVecMachine(kernel_type = 'RBF', C=1)

      start_time = time.time()
      svm_classifier.train(sel_w_vec_dict, true_m_set, true_nm_set)
      res_sets = svm_classifier.classify(sel_w_vec_dict)
      svm_time = time.time() - start_time

      res_file.write('    SVM classifier (RBF kernel, C=1):               ' + \
                     '                         %7.2f' % \
                     (svm_time) + os.linesep)
      res_file.flush()

      kmeans_classifier = classification.KMeans(dist_measure = mymath.distL2,
                                                max_iter_count = 10000,
                                                centroid_init = 'min/max')

      start_time = time.time()
      kmeans_classifier.train(sel_w_vec_dict, true_m_set, true_nm_set)
      res_sets = kmeans_classifier.classify(sel_w_vec_dict)
      kmeans_time = time.time() - start_time

      res_file.write('    K-means clustering (Euclidean distance, ' + \
                     'min/max initialisation):         %7.2f' % \
                     (kmeans_time) + os.linesep)
      res_file.flush()

      ffirst_classifier = classification.FarthestFirst(dist_me = mymath.distL2,
                                                       centroid_i = 'mode/max')
      start_time = time.time()
      ffirst_classifier.train(sel_w_vec_dict, true_m_set, true_nm_set)
      res_sets = ffirst_classifier.classify(sel_w_vec_dict)
      ffirst_time = time.time() - start_time

      res_file.write('    Farthest first clustering (Euclidean distance, ' + \
                     'mode/max initialisation): %7.2f' % (ffirst_time) + \
                     os.linesep)
      res_file.flush()

      two_s_classifier = classification.TwoStep(s1_match=(1.0,'threshold',0.5),
                                     s1_non_match_meth = (0.0,'threshold',0.5),
                                     random_sel = ('exponential', 5, 5),
                                     s2_classifier = ('svm', 'RBF', 1))
      start_time = time.time()
      two_s_classifier.train(sel_w_vec_dict, true_m_set, true_nm_set)
      res_sets = two_s_classifier.classify(sel_w_vec_dict)
      two_s_time = time.time() - start_time

      res_file.write('    Two step classifier (thresholds 0.5/0.5, rand ' + \
                     'expo 5%%/5%%, SVM RBF C=1):  %7.2f' % (two_s_time) + \
                     os.linesep)
      res_file.flush()

      # For nearest, calculate 5% of number of weight vectors - - - - - - - - -
      #
      sel_num = int(num_w_vec * 5.0/100.0)

      two_s_match_method = (1.0,'nearest', sel_num, True)
      two_s_non_match_method = (0.0,'nearest', sel_num, True)

      two_s_classifier = classification.TwoStep(s1_match = two_s_match_method,
                                    s1_non_match_meth = two_s_non_match_method,
                                    random_sel = ('exponential', 5, 5),
                                    s2_classifier = ('svm', 'RBF', 1))
      start_time = time.time()
      two_s_classifier.train(sel_w_vec_dict, true_m_set, true_nm_set)
      res_sets = two_s_classifier.classify(sel_w_vec_dict)
      two_s_time = time.time() - start_time

      res_file.write('    Two step classifier (nearest 5%, B, U, ' + \
                     'rand expo 5%%/5%%, SVM RBF C=1):    %7.2f' % \
                     (two_s_time) + os.linesep)
      res_file.flush()

      # Calculate imbalanced numbers
      #
      max_num_matches = min(data_set_a.num_records, data_set_b.num_records)

      rec_w_vec_ratio = float(max_num_matches) / \
                        float(num_w_vec-max_num_matches)

      if (rec_w_vec_ratio < 1):  # More non-matches than matches
        m_sel_num = max(1, int(sel_num*rec_w_vec_ratio))
        nm_sel_num = sel_num
      else:  # More matches than non-matches
        m_sel_num = sel_num
        nm_sel_num = max(1, int(sel_num/rec_w_vec_ratio))

      two_s_match_method = (1.0,'nearest', m_sel_num, True)
      two_s_non_match_method = (0.0,'nearest', nm_sel_num, True)

      two_s_classifier = classification.TwoStep(s1_match = two_s_match_method,
                                    s1_non_match_meth = two_s_non_match_method,
                                    random_sel = ('exponential', 5, 5),
                                    s2_classifier = ('svm', 'RBF', 1))
      start_time = time.time()
      two_s_classifier.train(sel_w_vec_dict, true_m_set, true_nm_set)
      res_sets = two_s_classifier.classify(sel_w_vec_dict)
      two_s_time = time.time() - start_time

      res_file.write('    Two step classifier (nearest 5%, IB, U,' + \
                     ' rand expo 5%%/5%%, SVM RBF C=1):   %7.2f' % \
                     (two_s_time) + os.linesep)
      res_file.flush()

## PC 25/09 #################################################################

for x in []:
  for y in []:
    for z in []:

      # -----------------------------------------------------------------------
      #
      do_random = True # Set to True or False to do or not do random examples

      # my_logger.setLevel(logging.INFO)  ##################################

      # Now add random weight vectors to hopefully improve classification - -
      #
      num_m_train_data =  len(m_set_tr)
      num_nm_train_data = len(nm_set_tr)

      if ((do_random == True) and \
          (num_m_train_data > 0) and (num_nm_train_data > 0)):

        rand_u_acc =   []  # List of accuracy results for 'uniform' random
        rand_u_prec =  []
        rand_u_reca =  []
        rand_u_fmeas = []
        rand_u_time =  []
        rand_l_acc =   []  # List of accuracy results for 'linear' random
        rand_l_prec =  []
        rand_l_reca =  []
        rand_l_fmeas = []
        rand_l_time =  []
        rand_e_acc =   []  # List of accuracy results for 'exponential' random
        rand_e_prec =  []
        rand_e_reca =  []
        rand_e_fmeas = []
        rand_e_time =  []

        # Add 10% random examples, at least 1 weight vector
        #
        num_m_random =  max(1, int(num_m_train_data*.10))
        num_nm_random = max(1, int(num_nm_train_data*.10))

        # Perform random training example insertion several times - - - - - - -
        #
        for r_iter in range(num_random_select_iterations):

          # Uniform selection of random examples - - - - - - - - - - - - - - -
          #
          start_time = time.time()
          m_setru,nm_setru = classification.get_random_examples(sel_w_vec_dict,
                                                                m_set_tr,
                                                                nm_set_tr,
                                                                'uniform',
                                                                num_m_random,
                                                                num_nm_random)
          r_m_set_tr =  m_set_tr.union(m_setru)
          r_nm_set_tr = nm_set_tr.union(nm_setru)

          conf_matrix = classification.svm_classifier(sel_w_vec_dict,
                                                      r_m_set_tr,
                                                      r_nm_set_tr,
                                                      match_check_funct)
          rand_u_time.append(time.time() - start_time)

          acc =  float(conf_matrix[0] + conf_matrix[3]) / sum(conf_matrix)
          prec = float(conf_matrix[0]) / float(conf_matrix[0] + conf_matrix[2])

          if (float(conf_matrix[0] + conf_matrix[1]) > 0):
            reca = float(conf_matrix[0]) / float(conf_matrix[0]+conf_matrix[1])
          else:
            reca = 0.0

          if ((prec != 0.0) or (reca != 0.0)):
            fmeas = 2*(prec*reca) / (prec+reca)
          else:
            fmeas = 0.0

          rand_u_acc.append(acc)
          rand_u_prec.append(prec)
          rand_u_reca.append(reca)
          rand_u_fmeas.append(fmeas)

          # Linear selection of random examples - - - - - - - - - - - - - - - -
          #
          start_time = time.time()
          m_setrl,nm_setrl = classification.get_random_examples(sel_w_vec_dict,
                                                                m_set_tr,
                                                                nm_set_tr,
                                                                'linear',
                                                                num_m_random,
                                                                num_nm_random)
          r_m_set_tr =  m_set_tr.union(m_setrl)
          r_nm_set_tr = nm_set_tr.union(nm_setrl)

          conf_matrix = classification.svm_classifier(sel_w_vec_dict,
                                                      r_m_set_tr,
                                                      r_nm_set_tr,
                                                      match_check_funct)
          rand_l_time.append(time.time() - start_time)

          acc =  float(conf_matrix[0] + conf_matrix[3]) / sum(conf_matrix)
          prec = float(conf_matrix[0]) / float(conf_matrix[0] + conf_matrix[2])

          if (float(conf_matrix[0] + conf_matrix[1]) > 0):
            reca = float(conf_matrix[0]) / float(conf_matrix[0]+conf_matrix[1])
          else:
            reca = 0.0

          if ((prec != 0.0) or (reca != 0.0)):
            fmeas = 2*(prec*reca) / (prec+reca)
          else:
            fmeas = 0.0

          rand_l_acc.append(acc)
          rand_l_prec.append(prec)
          rand_l_reca.append(reca)
          rand_l_fmeas.append(fmeas)

          # Exponential selection of random examples - - - - - - - - - - - - -
          #
          start_time = time.time()
          m_setrl,nm_setrl = classification.get_random_examples(sel_w_vec_dict,
                                                                m_set_tr,
                                                                nm_set_tr,
                                                                'exponential',
                                                                num_m_random,
                                                                num_nm_random)
          r_m_set_tr =  m_set_tr.union(m_setrl)
          r_nm_set_tr = nm_set_tr.union(nm_setrl)

          conf_matrix = classification.svm_classifier(sel_w_vec_dict,
                                                      r_m_set_tr,
                                                      r_nm_set_tr,
                                                      match_check_funct)
          rand_e_time.append(time.time() - start_time)

          acc =  float(conf_matrix[0] + conf_matrix[3]) / sum(conf_matrix)
          prec = float(conf_matrix[0]) / float(conf_matrix[0] + conf_matrix[2])

          if (float(conf_matrix[0] + conf_matrix[1]) > 0):
            reca = float(conf_matrix[0]) / float(conf_matrix[0]+conf_matrix[1])
          else:
            reca = 0.0

          if ((prec != 0.0) or (reca != 0.0)):
            fmeas = 2*(prec*reca) / (prec+reca)
          else:
            fmeas = 0.0

          rand_e_acc.append(acc)
          rand_e_prec.append(prec)
          rand_e_reca.append(reca)
          rand_e_fmeas.append(fmeas)

        # Get final measurement values
        #
        acc_u =   sum(rand_u_acc) / num_random_select_iterations
        prec_u =  sum(rand_u_prec) / num_random_select_iterations
        reca_u =  sum(rand_u_reca) / num_random_select_iterations
        fmeas_u = sum(rand_u_fmeas) / num_random_select_iterations
        time_u =  sum(rand_u_time) / num_random_select_iterations
        acc_l =   sum(rand_l_acc) / num_random_select_iterations
        prec_l =  sum(rand_l_prec) / num_random_select_iterations
        reca_l =  sum(rand_l_reca) / num_random_select_iterations
        fmeas_l = sum(rand_l_fmeas) / num_random_select_iterations
        time_l =  sum(rand_l_time) / num_random_select_iterations
        acc_e =   sum(rand_e_acc) / num_random_select_iterations
        prec_e =  sum(rand_e_prec) / num_random_select_iterations
        reca_e =  sum(rand_e_reca) / num_random_select_iterations
        fmeas_e = sum(rand_e_fmeas) / num_random_select_iterations
        time_e =  sum(rand_e_time) / num_random_select_iterations

        res_file.write('    %19s |%9f | %9f |%9f | %9f |%9f' % \
                       ('(+ rand uni)'.rjust(19), acc_u, prec_u, reca_u,
                       fmeas_u, time_u) + os.linesep)

        res_file.write('    %19s |%9f | %9f |%9f | %9f |%9f' % \
                       ('(+ rand lin)'.rjust(19), acc_l, prec_l, reca_l,
                       fmeas_l, time_l) + os.linesep)

        res_file.write('    %19s |%9f | %9f |%9f | %9f |%9f' % \
                       ('(+ rand exp)'.rjust(19), acc_e, prec_e, reca_e,
                       fmeas_e, time_e) + os.linesep)
        res_file.flush()

      # my_logger.setLevel(logging.WARNING)  #################################

    res_file.write(os.linesep)

    res_list.append(exp_res_list)  # Add results of this experiment

  res_dict[res_list_key] = res_list  # Save for later analysis

res_file.write('='*80+os.linesep)
res_file.write(os.linesep)
res_file.flush()

#==============================================================================
# Now print summarised results

for data_set_decr_tuple in res_dict:  # Loop over all data set results

  if (data_set_decr_tuple[0] == data_set_decr_tuple[1]):  # A deduplication
    res_file.write('Summary for data set:' + os.linesep)
    res_file.write('---------------------' + os.linesep)
    res_file.write('  %s' % (data_set_decr_tuple[0]) + os.linesep)

  else:  # A linkage
    res_file.write('Summary for data sets:' + os.linesep)
    res_file.write('----------------------' + os.linesep)
    res_file.write('  %s' % (data_set_decr_tuple[0]) + os.linesep)
    res_file.write('  %s' % (data_set_decr_tuple[1]) + os.linesep)

  res_file.write(os.linesep)

  data_set_exp_list = res_dict[data_set_decr_tuple]

  # Average results for each technique for different field comparisons
  #
  avrg_acc_dict =   {}
  avrg_prec_dict =  {}
  avrg_reca_dict =  {}
  avrg_fmeas_dict = {}

  num_exp = 0

  for exp_list in data_set_exp_list:

    for (tech_name_str, acc, prec, reca, fmeas) in exp_list[1:]:

      acc_sum = avrg_acc_dict.get(tech_name_str, 0.0)
      acc_sum += acc
      avrg_acc_dict[tech_name_str] = acc_sum

      prec_sum = avrg_prec_dict.get(tech_name_str, 0.0)
      prec_sum += prec
      avrg_prec_dict[tech_name_str] = prec_sum

      reca_sum = avrg_reca_dict.get(tech_name_str, 0.0)
      reca_sum += reca
      avrg_reca_dict[tech_name_str] = reca_sum

      fmeas_sum = avrg_fmeas_dict.get(tech_name_str, 0.0)
      fmeas_sum += fmeas
      avrg_fmeas_dict[tech_name_str] = fmeas_sum

    num_exp += 1

  tech_name_str_list = avrg_acc_dict.keys()
  tech_name_str_list.sort()

  # Report average values
  #
  res_file.write('  Average measure results:' + os.linesep)
  res_file.write('  ------------------------' + os.linesep)
  res_file.write('                                  Method |  Accuracy |' + \
                 ' Precision |    Recall | F-measure' + os.linesep)
  res_file.write('  '+'-'*86 + os.linesep)

  for tech_name_str in tech_name_str_list:
    avrg_acc =   avrg_acc_dict[tech_name_str] / num_exp
    avrg_prec =  avrg_prec_dict[tech_name_str] / num_exp
    avrg_reca =  avrg_reca_dict[tech_name_str] / num_exp
    avrg_fmeas = avrg_fmeas_dict[tech_name_str] / num_exp

    res_file.write('  %38s | %9f | %9f | %9f | %9f' % \
                   (tech_name_str.rjust(37), avrg_acc, avrg_prec, avrg_reca,
                    avrg_fmeas) + os.linesep)
  res_file.write(os.linesep)

  res_file.write('  '+'='*78+os.linesep)
  res_file.write(os.linesep)

# Each exp_list has:
# 0) name of fields selected
# - Tuples with 5 elements (comment-str, acc, prec, reca, fmeas):
#   1) One-dim optimal classifier, minimise FP and FN
#   2) SVM supervised
#   3) K-means with L1-distance
#   4) K-means with L2-distance
#   5) .. 17 versions of two-step classification


sys.exit()  ####### PC 24/07/07

### OLD syuff below ==============================


# get examples nearest: yop x weight vctors or top x UNIQUE weight vectors!!!


m_sete, nm_sete = classification.get_examples(wv_dict,
                            match_method=(1.0,'threshold',0.3),
                            non_match_method=(0.0,'threshold',0.3))
print 'examples threshold:', len(m_sete),len(nm_sete)
measurements.accuracy(m_sete, nm_sete)

m_sete, nm_sete = classification.get_examples(wv_dict,
                            match_method=(1.0,'nearest',100),
                            non_match_method=(0.0,'nearest',100))
print 'examples nearest 100:', len(m_sete),len(nm_sete)
measurements.accuracy(m_sete, nm_sete)

imbalance = float(rest_ds.num_records) / len(wv_dict)
print 'imblance:', imbalance

m_sete, nm_sete = classification.get_examples(wv_dict,
                            match_method=(1.0,'nearest',100),
                           non_match_method=(0.0,'nearest',int(100/imbalance)))
print 'examples nearest 100 / imbalance:', len(m_sete),len(nm_sete)
measurements.accuracy(m_sete, nm_sete)
