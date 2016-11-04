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
# The Original Software is: "evalIndexing.py"
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

"""Module to evaluate the various index techniques implemented in the Febrl
   module 'indexing.py'.

   Uses Febrl-0.4.

   Peter Christen, 2010/08/27

   Defines a series of data sets and indices, then runs them with a variety of
   parameter settings and evaluates the qualities and complexities of the
   resulting weight vectors generated.

   This version reads the data set and index name from the command line as
   arguments, and only runs the specified tests.

   Usage:  python evalIndexing2.py [data set] [index method]

   If the [index method] is given as 'all' then all index methods will be run.

   This module was used for the experiments presented in the article:

   A Survey of Indexing Techniques for Scalable Record Linkage and
   Deduplication

   Peter Christen

   IEEE Transactions on Knowledge and Data Engineering (TKDE), June 2011.
"""

# =============================================================================
# Imports go here

import logging
import os
import sys
import time

import auxiliary
import classification
import comparison
import dataset
import encode
import indexing
import mymath
import stringcmp

# =============================================================================
# Get command line arguments

arg_data_set_name =     sys.argv[1].lower()
arg_index_method_name = sys.argv[2].lower()

# =============================================================================
# Various settings

progress_precentage = 10

# File name for writing results to (results will be appended to this file)
#
result_file_name = './results/evalIndexing2.res'

# =============================================================================
# Define a project logger

my_logger = logging.getLogger()  # New logger at root level
my_logger.setLevel(logging.WARNING)
#my_logger.setLevel(logging.INFO)

# =============================================================================
# Open the results file
#
res_file = open(result_file_name, 'a')
res_file.write(os.linesep+os.linesep+'Experiment started %s:' % (\
               time.strftime('%Y%m%d-%H%M')) + os.linesep)

# =============================================================================
# Generate four dictionaries with information about:
# 1) data sets, record comparator, and functions to check matches/non-matches
#    and identifiers.
# 2) indexing definitions (using a variety of field value combinations).
# 3) actual index definitions with various parameter settings.
# 4) results achieved from running experiments (RR, PC and PQ).
#
experiment_dict = {}
index_def_dict =  {}
indexing_dict =   {}
result_dict =     {}

# =============================================================================
# Define original input data sets plus quality assessment functions, and index
# definitions

# -----------------------------------------------------------------------------
# The publicly available Census data set with synthetic personal names and
# addresses (taken from SecondString data repository).
#
# - The 'entity_id' attribute (2nd attribute) contains entity numbers.
# - No record identifier is available.
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
                                             ('middle_initial',4),
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
                                             ('middle_initial',4),
                                             ('zipcode',5),
                                             ('suburb',6)],
                    file_name = './data/secondstring/censusTextSegmentedB.tab')

# Define field and record comparators
#
census_entity_id_exact = comparison.FieldComparatorExactString(desc = \
                                                             'entity_id_exact')
census_surname_winkler =    comparison.FieldComparatorWinkler(thres=0,
                                                      desc = 'surname_winkler')
census_given_name_winkler = comparison.FieldComparatorWinkler(thres=0,
                                                   desc = 'given_name_winkler')
census_suburb_winkler =     comparison.FieldComparatorWinkler(thres=0,
                                                       desc = 'suburb_winkler')

census_fc_list = [(census_entity_id_exact,    'entity_id',    'entity_id'),
                  (census_surname_winkler,    'surname',      'surname'),
                  (census_given_name_winkler, 'given_name',   'given_name'),
                  (census_suburb_winkler,     'suburb',       'suburb')]

census_rec_comp = comparison.RecordComparator(census_ds_A, census_ds_B,
                                              census_fc_list,
                                              'Census record comparator')

# Function to be used to check for true matches and non-matches
#
def census_check_funct(rec1, rec2):
  return (rec1[1] == rec2[1])

# Function to be used to extract the record identifier from a raw record
#
def census_get_id_funct(rec):
  return rec[1]

# Insert into data set dictionary
#
experiment_dict['census'] = ['Census', census_ds_A, census_ds_B,
                             census_rec_comp,
                             census_check_funct, census_get_id_funct]

# Set-up index definitions
#
census_index_def1 = \
  [[['surname', 'surname',      False, False, None, [encode.dmetaphone,3]],
    ['given_name','given_name', False, False, None, [encode.dmetaphone,3]]],
   [['suburb', 'suburb',        False, False, None, [encode.dmetaphone,3]],
    ['zipcode', 'zipcode',      False, False, None, [                   ]]]]
census_index_def2 = \
  [[['surname', 'surname', False, False, 3,             [encode.dmetaphone,3]],
    ['middle_initial','middle_initial', False, False,1, []],
    ['zipcode', 'zipcode', False, False, None,          []]],
   [['given_name', 'given_name', False, False, None,    [encode.dmetaphone,3]],
    ['suburb', 'suburb', False, False, None,            [encode.dmetaphone,3]]]]
census_index_def3 = \
  [[['suburb', 'suburb', False, False, None,            [encode.dmetaphone,3]],
    ['surname', 'surname', False, False, 3,             [encode.dmetaphone,3]]],
   [['zipcode', 'zipcode',      False, False, None, []],
    ['given_name', 'given_name', False, False, None,    [encode.dmetaphone,3]]]]

index_def_dict['census'] = [census_index_def1, census_index_def2,
                            census_index_def3]


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

# Define field and record comparators
#
cora_paper_id_exact = comparison.FieldComparatorExactString(desc = \
                                                              'paper_id_exact')
cora_author_list_winkler = comparison.FieldComparatorWinkler(thres=0,
                                                 desc = 'author_list_winkler')
cora_title_winkler = comparison.FieldComparatorWinkler(thres=0,
                                                 desc = 'title_winkler')
cora_conf_journal_winkler = comparison.FieldComparatorWinkler(thres=0,
                                                 desc = 'conf_journal_winkler')

cora_fc_list = [(cora_paper_id_exact,       'paper_id',    'paper_id'),
                (cora_author_list_winkler,  'author_list', 'author_list'),
                (cora_title_winkler,        'title',       'title'),
                (cora_conf_journal_winkler, 'conf_journal','conf_journal')]

cora_rec_comp = comparison.RecordComparator(cora_ds, cora_ds, cora_fc_list,
                                            'Cora record comparator')

# Function to be used to check for true matches and non-matches
#
def cora_check_funct(rec1, rec2):
  return (rec1[1] == rec2[1])

# Function to be used to extract the record identifier from a raw record
#
def cora_get_id_funct(rec):
  return rec[1]

# Insert into data set dictionary
#
experiment_dict['cora'] = ['Cora', cora_ds, cora_ds, cora_rec_comp,
                           cora_check_funct, cora_get_id_funct]

# Set-up index definitions
#
cora_index_def1 = \
  [[['author_list','author_list', False, False, None,   [encode.dmetaphone,3]],
    ['title', 'title', False, False, None,              [encode.dmetaphone,3]]],
   [['conf_journal','conf_journal', False, False, None, [encode.dmetaphone,3]],
    ['year', 'year', False, False, None,                []                  ]]]
cora_index_def2 = \
  [[['affiliation','affiliation', False, False, None, [encode.dmetaphone,3]],
    ['location', 'location', False, False, None,      [encode.dmetaphone,3]]],
   [['year', 'year', False, False, None,              []                   ],
    ['title', 'title', False, False, None,            [encode.dmetaphone,3]]]]
cora_index_def3 = \
  [[['year', 'year', False, False, None,              []],
    ['author_list','author_list', False, False, None,   [encode.dmetaphone,3]]],
   [['conf_journal','conf_journal', False, False, None, [encode.dmetaphone,3]],
    ['title', 'title', False, False, None,            [encode.dmetaphone,3]]]]

# Alternative index definitions
#
#cora_index_def1 = \
#  [[['author_list','author_list', False, False, None, [encode.dmetaphone,3]],
#    ['year', 'year', False, False, None,              [                   ]]],
#   [['year', 'year', False, False, None,              [                   ]],
#    ['title', 'title', False, False, None,            [encode.dmetaphone,3]]]]
#cora_index_def2 = \
#  [[['conf_journal','conf_journal', False, False, None, [encode.dmetaphone,3]],
#    ['year', 'year', False, False, None,                []                  ]],
#   [['title', 'title', False, True, None,               [encode.dmetaphone,3]],
#    ['author_list','author_list', False, True, None,   [encode.dmetaphone,3]]]]
#cora_index_def3 = \
#  [[['author_list','author_list', True, False, None,   [encode.dmetaphone,3]],
#    ['conf_journal','conf_journal', True, False, None, [encode.dmetaphone,3]]],
#   [['author_list','author_list', True, True, None,    [encode.dmetaphone,3]],
#    ['year', 'year', False, False, None,               [                   ]]]]

index_def_dict['cora'] = [cora_index_def1, cora_index_def2, cora_index_def3]


# -----------------------------------------------------------------------------
# The publicly available Restaurant data set (Fodors/ Zagats) restaurant names
# and addresses (taken from SecondString data repository).
#
# - The 'class' attribute contains entity numbers.
# - No record identifier is available.
#
rest_ds = dataset.DataSetCSV(description='Restaurant data set',
                             access_mode='read',
                             rec_ident='rec_id',
                             header_line=True,
                             file_name = './data/secondstring/restaurant.csv')

# Define field and record comparators
#
rest_class_exact = comparison.FieldComparatorExactString(desc = 'class_exact')
rest_name_winkler = comparison.FieldComparatorWinkler(thres=0,
                                                      desc = 'name_winkler')
rest_addr_winkler = comparison.FieldComparatorWinkler(thres=0,
                                                      desc = 'addr_winkler')
rest_city_winkler = comparison.FieldComparatorWinkler(thres=0,
                                                      desc = 'city_winkler')

rest_fc_list = [(rest_class_exact,  'class', 'class'),
                (rest_name_winkler, 'name',  'name'),
                (rest_addr_winkler, 'addr',  'addr'),
                (rest_city_winkler, 'city',  'city')]

rest_rec_comp = comparison.RecordComparator(rest_ds, rest_ds, rest_fc_list,
                                            'Restaurant record comparator')

# Function to be used to check for true matches and non-matches
#
def rest_check_funct(rec1, rec2):
  return (rec1[-1] == rec2[-1])

# Function to be used to extract the record identifier from a raw record
#
def rest_get_id_funct(rec):
  return rec[-1]

# Insert into data set dictionary
#
experiment_dict['rest'] = ['Restaurant', rest_ds, rest_ds, rest_rec_comp,
                           rest_check_funct, rest_get_id_funct]

# Set-up index definitions
#
rest_index_def1 = \
  [[['phone', 'phone', False, False, None, [encode.get_substring,0,4]],
    ['type',  'type',  False, False, None, [encode.dmetaphone,3     ]]],
   [['name', 'name',   False, False, None, [encode.dmetaphone,3     ]],
    ['city', 'city',   False, False, None, [encode.dmetaphone,3     ]]]]
rest_index_def2 = \
  [[['city', 'city',   False, False, None, [encode.dmetaphone,3     ]],
    ['phone', 'phone', False, False, None, [encode.get_substring,0,4]]],
   [['addr', 'addr',   False, False, None, [encode.dmetaphone,3     ]],
    ['city', 'city',   False, False, None, [encode.dmetaphone,3     ]]]]
rest_index_def3 = \
  [[['type', 'type',   False, False, None, [encode.dmetaphone,3     ]],
    ['name', 'name',   False, False, None, [encode.dmetaphone,3     ]]],
   [['addr', 'addr',   False, False, None, [encode.dmetaphone,3     ]],
    ['type',  'type',  False, False, None, [encode.dmetaphone,3     ]]]]

# Alternative index definitions
#
#rest_index_def1 = \
#  [[['phone', 'phone', False, False, None, [encode.get_substring,0,4]],
#    ['type',  'type', False, False, None,  [encode.dmetaphone,3]     ]],
#   [['name', 'name', False, False, None,   [encode.dmetaphone,3]     ],
#    ['city', 'city', False, False, None,   [encode.dmetaphone,3]     ]]]
#rest_index_def2 = \
#  [[['addr', 'addr', False, False, None,   [encode.dmetaphone,3]     ],
#    ['name', 'name', False, False, None,   [encode.dmetaphone,3]     ]],
#   [['type',  'type', False, False, None,  [encode.dmetaphone,3]     ],
#    ['city', 'city', False, False, None,   [encode.dmetaphone,3]     ]]]
#rest_index_def3 = \
#  [[['phone', 'phone', False, False, None, [encode.get_substring,0,4]],
#    ['name', 'name', False, False, None,   [encode.dmetaphone,3]     ]],
#   [['addr', 'addr', False, False, None,   [encode.dmetaphone,3]     ],
#    ['type',  'type', False, False, None,  [encode.dmetaphone,3]     ]]]

index_def_dict['rest'] = [rest_index_def1, rest_index_def2, rest_index_def3]


# -----------------------------------------------------------------------------
# CD data set, as downloaded from:
# http://www.hpi.uni-potsdam.de/naumann/projekte/repeatability/datasets/ \
#   cd_datasets.html
# This data set has been converted from an XML format into a CSV file, and the
# identifiers of know true duplicates have been modified such that groups of
# duplicates have the same identifier values.
#
# - No record identifier is available.
# - The 'identifier' field contains unique values for individual or groups of
#   CDs that are the same entity
#
cddb_ds = dataset.DataSetCSV(description='CDDB data set',
                             access_mode='read',
                             rec_ident='rec_id',
                             header_line=True,
                             file_name = './data/cddb/cddb.csv')

# Define field and record comparators
#
cddb_class_exact = comparison.FieldComparatorExactString(desc = 'class_exact')
cddb_artist_winkler = comparison.FieldComparatorWinkler(thres=0,
                                                        desc = 'artist_winkler')
cddb_title_winkler = comparison.FieldComparatorWinkler(thres=0,
                                                       desc = 'title_winkler')

cddb_fc_list = [(cddb_class_exact,     'identifier', 'identifier'),
                (cddb_artist_winkler,  'artist',     'artist'),
                (cddb_title_winkler,   'title',      'title')]

cddb_rec_comp = comparison.RecordComparator(cddb_ds, cddb_ds, cddb_fc_list,
                                            'CDDB record comparator')

# Function to be used to check for true matches and non-matches
#
def cddb_check_funct(rec1, rec2):
  return (rec1[0] == rec2[0])

# Function to be used to extract the record identifier from a raw record
#
def cddb_get_id_funct(rec):
  return rec[0]

# Insert into data set dictionary
#
experiment_dict['cddb'] = ['CDDB', cddb_ds, cddb_ds, cddb_rec_comp,
                           cddb_check_funct, cddb_get_id_funct]

# Set-up index definitions
#
cddb_index_def1 = \
  [[['artist','artist',      False, False, None, [encode.dmetaphone,3]],
    ['title','title',        False, False, None, [encode.dmetaphone,3]],
    ['genre', 'genre',       False, False, None, [encode.dmetaphone,3]]],
   [['artist','artist',      False, False, None, [encode.dmetaphone,3]],
    ['title','title',        False, False, None, [encode.dmetaphone,3]],
    ['year', 'year',         False, False, None, [                   ]]]]
cddb_index_def2 = \
  [[['artist','artist',      False, False, None, [encode.dmetaphone,3]],
    ['category', 'category', False, False, None, [encode.dmetaphone,3]],
    ['genre', 'genre',       False, False, None, [encode.dmetaphone,3]]],
   [['artist','artist',      False, False, None, [encode.dmetaphone,3]],
    ['category', 'category', False, False, None, [encode.dmetaphone,3]],
    ['year', 'year',         False, False, None, [                   ]]]]
cddb_index_def3 = \
  [[['title','title',        False, False, None, [encode.dmetaphone,3]],
    ['category', 'category', False, False, None, [encode.dmetaphone,3]],
    ['genre', 'genre',       False, False, None, [encode.dmetaphone,3]]],
   [['title','title',        False, False, None, [encode.dmetaphone,3]],
    ['category', 'category', False, False, None, [encode.dmetaphone,3]],
    ['year', 'year',         False, False, None, [                   ]]]]

index_def_dict['cddb'] = [cddb_index_def1, cddb_index_def2, cddb_index_def3]


# -----------------------------------------------------------------------------
# Synthetic data sets of different sizes generated with Febrl data generator

# Define field and record comparators
#
synth_rec_id_exact = comparison.FieldComparatorExactString(desc='rec_id_exact')
synth_surname_winkler =    comparison.FieldComparatorWinkler(thres=0,
                                                      desc = 'surname_winkler')
synth_given_name_winkler = comparison.FieldComparatorWinkler(thres=0,
                                                   desc = 'given_name_winkler')
synth_suburb_winkler =     comparison.FieldComparatorWinkler(thres=0,
                                                       desc = 'suburb_winkler')

synth_fc_list = [(synth_rec_id_exact,       'rec_id',        'rec_id'),
                 (synth_surname_winkler,    'surname',       'surname'),
                 (synth_given_name_winkler, 'given_name',    'given_name'),
                 (synth_suburb_winkler,     'suburb',        'suburb')]

# Function to be used to check for true matches and non-matches
#
def synth_check_funct(rec1, rec2):
    return ((rec1[0][:-1] == rec2[0][:-1]) and (rec1[0][-1] != rec2[0][-1]))

# Function to be used to extract the record identifier from a raw record
#
def synth_get_id_funct(rec):
  return rec[0][:-1]

# Loop over different data set sizes - - - - - - - - - - - - - - - - - - - - -
#
for size in ['B_1000',  'C_1000',  'E_1000',   'F_1000',
             'B_2500',  'C_2500',  'E_2500',   'F_2500',
             'B_5000',  'C_5000',  'E_5000',   'F_5000',
             'B_10000', 'C_10000', 'E_10000',  'F_10000',
             'B_25000', 'C_25000', 'E_25000',  'F_25000',
             'B_50000', 'C_50000', 'E_50000',  'F_50000']:

  if (size[0] in ['A','B','C']):  # A deduplication

    ds_file_name = 'dataset_%s.csv.gz' % (size)

    synth_ds1 = dataset.DataSetCSV(description='Febrl synthetic data set %s' % \
                                                (size),
                                   access_mode='read',
                                   rec_ident='rec_id',
                                   header_line=True,
                                   file_name = './data/dedup-dsgen/%s' % \
                                               (ds_file_name))
    synth_ds2 = synth_ds1

  else:  # A linkage

    ds_file_name1 = 'dataset_%s_org_%s.csv.gz' % (size[0], size[2:])
    ds_file_name2 = 'dataset_%s_dup_%s.csv.gz' % (size[0], size[2:])

    synth_ds1 = dataset.DataSetCSV(description='Febrl synthetic data set %s' % \
                                                (size)+' (org)',
                                   access_mode='read',
                                   rec_ident='rec_id',
                                   header_line=True,
                                   file_name = './data/link-dsgen/%s' % \
                                               (ds_file_name1))

    synth_ds2 = dataset.DataSetCSV(description='Febrl synthetic data set %s' % \
                                                (size)+' (dup)',
                                   access_mode='read',
                                   rec_ident='rec_id',
                                   header_line=True,
                                   file_name = './data/link-dsgen/%s' % \
                                               (ds_file_name2))

  synth_rec_comp = comparison.RecordComparator(synth_ds1, synth_ds2,
                                               synth_fc_list,
                                               'Febrl synthetic data record' \
                                               + ' comparator')

  size = size.lower()  # For inserting into dictionary

  # Insert into data set dictionary
  #
  experiment_dict[size] = ['Synthetic '+size, synth_ds1, synth_ds2,
                           synth_rec_comp,
                           synth_check_funct, synth_get_id_funct]

  # Set-up index definitions
  #
  synth_index_def1 = \
    [[['surname', 'surname', False, False, None,       [encode.dmetaphone,3]],
      ['address_1','address_1',False,False,None,       [encode.dmetaphone,3]],
      ['postcode', 'postcode', False, False, None,     []                   ]],
     [['given_name', 'given_name', False, False, None, [encode.dmetaphone,3]],
      ['address_1','address_1',False,False,None,       [encode.dmetaphone,3]],
      ['postcode', 'postcode', False, False, None,     []                   ]]]
  synth_index_def2 = \
    [[['surname', 'surname', False, False, None,       [encode.dmetaphone,3]],
      ['suburb','suburb',False,False,None,             [encode.dmetaphone,3]],
      ['age', 'age', False, False, None,               []                   ]],
     [['given_name', 'given_name', False, False, None, [encode.dmetaphone,3]],
      ['suburb','suburb',False,False,None,             [encode.dmetaphone,3]],
      ['age', 'age', False, False, None,               []                   ]]]
  synth_index_def3 = \
    [[['surname', 'surname', False, False, None,       [encode.dmetaphone,3]],
      ['address_1','address_1',False,False,None,       [encode.dmetaphone,3]],
      ['age', 'age',           False, False, None,     []                   ]],
     [['given_name', 'given_name', False, False, None, [encode.dmetaphone,3]],
      ['suburb','suburb',False,False,None,       [encode.dmetaphone,3]],
      ['postcode', 'postcode', False, False, None,     []                   ]]]

  index_def_dict[size] = [synth_index_def1, synth_index_def2, synth_index_def3]


for k in experiment_dict.keys():
  assert k in index_def_dict
#  print k
  assert len(experiment_dict[k]) == 6
#  print experiment_dict[k]
  assert len(index_def_dict[k]) == 3
#  print index_def_dict[k]
#  print
assert set(experiment_dict.keys()) == set(index_def_dict.keys())


# =============================================================================
# Define indices for all data sets and their three index definitions each

# Explicitly define experiments order
#
experiment_name_list = ['census', 'cora', 'rest', 'cddb',
                        'b_1000',  'c_1000',  'e_1000',   'f_1000',
                        'b_2500',  'c_2500',  'e_2500',   'f_2500',
                        'b_5000',  'c_5000',  'e_5000',   'f_5000',
                        'b_10000', 'c_10000', 'e_10000',  'f_10000',
                        'b_25000', 'c_25000', 'e_25000',  'f_25000',
                        'b_50000', 'c_50000', 'e_50000',  'f_50000']

# Data set name list is taken from command line argument
#
assert arg_data_set_name in experiment_name_list, \
       ('Unknown data set name given:', arg_data_set_name)

experiment_name_list = [arg_data_set_name]

for ds_name in experiment_name_list:

  data_set_name = experiment_dict[ds_name][0]
  print '='*70
  print
  print 'Data set:', data_set_name, '      ',  time.ctime()
  print '----------'+'-'*len(data_set_name)

  data_set1 = experiment_dict[ds_name][1]
  data_set2 = experiment_dict[ds_name][2]

  if (data_set1 == data_set2):
    task = 'dedup'
    print '  Task: Deduplication'

  else:
    task = 'link'
    print '  Task: Linkage'

  rec_cmp =   experiment_dict[ds_name][3]

  check_match_funct = experiment_dict[ds_name][4]
  get_id_funct =      experiment_dict[ds_name][5]

  #print data_set1.description
  #print data_set2.description
  #print rec_cmp.description
  #print

  index_def_list = index_def_dict[ds_name]
  #print index_def_list

  ds_index_list = []  # All indices defined for this data set (each entry is a
                      # tuple with the index method name and the actual index)

  # Blocking index - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  for this_index_def in index_def_list:

    block_index = indexing.BlockingIndex(description = 'Blocking index',
                                         dataset1 = data_set1,
                                         dataset2 = data_set2,
                                         rec_comparator = rec_cmp,
                                         progress=progress_precentage,
                                         index_def = this_index_def)

    ds_index_list.append(['blocking', block_index])

  # Sorted neighbourhood (inverted index based) index - - - - - - - - - - - - -
  #
  for w in [2,3,5,7,10]:

    for this_index_def in index_def_list:

      sorted_index = indexing.SortingIndex(description = 'Sorting index: ' + \
                                           'w=%d' % (w),
                                           dataset1 = data_set1,
                                           dataset2 = data_set2,
                                           rec_comparator = rec_cmp,
                                           progress=progress_precentage,
                                           index_def = this_index_def,
                                           window_s = w)

      ds_index_list.append(['sorted-inv-index', sorted_index])

  # Sorted neighbourhood (array based) index - - - - - - - - - - - - - - - - -
  #
  for w in [2,3,5,7,10]:

    for this_index_def in index_def_list:

      sorted_index = indexing.SortingArrayIndex(description = 'Sorting ' + \
                                                'array index: w=%d' % (w),
                                                dataset1 = data_set1,
                                                dataset2 = data_set2,
                                                rec_comparator = rec_cmp,
                                                progress=progress_precentage,
                                                index_def = this_index_def,
                                                window_s = w)

      ds_index_list.append(['sorted-array', sorted_index])

  # Adaptive sorted neighbourhood index - - - - - - - - - - - - - - - - - - - -
  #
  for str_cmp_funct in [('Jaro',stringcmp.jaro),('Bigram',stringcmp.bigram),
                        ('ED',stringcmp.editdist),('LCS',stringcmp.lcs)]:
    for str_cmp_thres in [0.8, 0.9]:

      for this_index_def in index_def_list:

        sorted_index = indexing.AdaptSortingIndex(description = 'Adaptive ' + \
                                                  'sorting index: str_cmp=%s' \
                                                  % (str_cmp_funct[0]) + \
                                                  ', thres=%.1f' % \
                                                  (str_cmp_thres),
                                                  dataset1 = data_set1,
                                                  dataset2 = data_set2,
                                                  rec_comparator = rec_cmp,
                                                  progress=progress_precentage,
                                                  index_def = this_index_def,
                                                  str_cmp_funct = \
                                                              str_cmp_funct[1],
                                                  str_cmp_thres = \
                                                              str_cmp_thres)

        ds_index_list.append(['adapt-sorted', sorted_index])

  # Suffix array based indexing - - - - - - - - - - - - - - - - - - - - - - - -
  #
  for min_len in [3,5]:
    for max_block in [5,10,20]:

      for this_index_def in index_def_list:

        sarray_index = indexing.SuffixArrayIndex(desc='Suffix array index:' + \
                                                 ' min_len=%d, max_block=%d' \
                                                 % (min_len, max_block),
                                                 dataset1 = data_set1,
                                                 dataset2 = data_set2,
                                                 rec_comparator = rec_cmp,
                                                 progress=progress_precentage,
                                                 index_def = this_index_def,
                                                 padd = True,
                                                 block_method = (min_len, \
                                                                 max_block),
                                                 suffix_method = 'suffixonly')

        ds_index_list.append(['suffix-array', sarray_index])

  # Suffix array based indexing using all substrings - - - - - - - - - - - - - -
  #
  for min_len in [3,5]:
    for max_block in [5,10,20]:

      for this_index_def in index_def_list:

        sarray_index = indexing.SuffixArrayIndex(desc='Suffix array index' + \
                                                 'sub-string: ' + \
                                                 'min_len=%d, max_block=%d' \
                                                 % (min_len, max_block),
                                                 dataset1 = data_set1,
                                                 dataset2 = data_set2,
                                                 rec_comparator = rec_cmp,
                                                 progress=progress_precentage,
                                                 index_def = this_index_def,
                                                 padd = True,
                                                 block_method = (min_len, \
                                                                 max_block),
                                                 suffix_method = 'allsubstr')

        ds_index_list.append(['suffix-array-substr', sarray_index])


  # Robust suffix array based indexing - - - - - - - - - - - - - - - - - - - -
  #
  for min_len in [3,5]:
    for max_block in [5,10,20]:
      for str_cmp_funct in [('Jaro',stringcmp.jaro),('Bigram',stringcmp.bigram),
                            ('ED',stringcmp.editdist),('LCS',stringcmp.lcs)]:
        for str_cmp_thres in [0.8, 0.9]:

          for this_index_def in index_def_list:

            rsarray_index = indexing.RobustSuffixArrayIndex(desc = \
                                        'Robust suffix array index:' + \
                                        ' min_len=%d, max_block=%d,' \
                                        % (min_len, max_block) + \
                                        'str_cmp=%s, thres=%.1f' % \
                                        (str_cmp_funct[0], str_cmp_thres),
                                        dataset1 = data_set1,
                                        dataset2 = data_set2,
                                        rec_comparator = rec_cmp,
                                        progress=progress_precentage,
                                        index_def = this_index_def,
                                        padd = True,
                                        block_method = (min_len, max_block),
                                        str_cmp_funct = str_cmp_funct[1],
                                        str_cmp_thres = str_cmp_thres)

            ds_index_list.append(['robust-suffix-array', rsarray_index])

  # Q-gram indexing - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  for q in [2,3]:
    for thres in [0.8, 0.9]:

      for this_index_def in index_def_list:

        qgram_index = indexing.QGramIndex(desc = 'Q-gram index: ' + \
                                     'q=%d, thres=%.1f' % (q, thres),
                                     dataset1 = data_set1,
                                     dataset2 = data_set2,
                                     rec_comparator = rec_cmp,
                                     progress=progress_precentage,
                                     index_def = this_index_def,
                                     padd = True,
                                     q = q,
                                     threshold=thres)

        ds_index_list.append(['q-gram', qgram_index])

  # Canopy indexing (threshold based) - - - - - - - - - - - - - - - - - - - - -
  #
  for q in [2,3]:
    for canopy_method in [('tfidf',   'threshold', 0.9, 0.8),
                          ('jaccard', 'threshold', 0.9, 0.8),
                          ('tfidf',   'threshold', 0.8, 0.7),
                          ('jaccard', 'threshold', 0.8, 0.7)]:

      for this_index_def in index_def_list:

        canopy_index = indexing.CanopyIndex(desc='Canopy TH index: q=%d, %s' % \
                                       (q, str(canopy_method)),
                                       dataset1 = data_set1,
                                       dataset2 = data_set2,
                                       rec_comparator = rec_cmp,
                                       progress=progress_precentage,
                                       index_def = this_index_def,
                                       padd = True,
                                       q = q,
                                       canopy_m = canopy_method)

        ds_index_list.append(['canopy-th', canopy_index])

  # Canopy indexing (nearest neighbour based) - - - - - - - - - - - - - - - - -
  #
  for q in [2,3]:
    for canopy_method in [('tfidf',   'nearest',  5, 10),
                          ('jaccard', 'nearest',  5, 10),
                          ('tfidf',   'nearest', 10, 20),
                          ('jaccard', 'nearest', 10, 20)]:

      for this_index_def in index_def_list:

        canopy_index = indexing.CanopyIndex(desc='Canopy NN index: q=%d, %s' % \
                                       (q, str(canopy_method)),
                                       dataset1 = data_set1,
                                       dataset2 = data_set2,
                                       rec_comparator = rec_cmp,
                                       progress=progress_precentage,
                                       index_def = this_index_def,
                                       padd = True,
                                       q = q,
                                       canopy_m = canopy_method)

        ds_index_list.append(['canopy-nn', canopy_index])

  # StringMap based indexing (threshold based) - - - - - - - - - - - - - - - - -
  #
  for (dim, subdim) in [(15,3), (20,5)]:
    for grid in [100,1000]:
      for str_cmp_funct in [('Jaro',stringcmp.jaro),('Bigram',stringcmp.bigram),
                            ('ED',stringcmp.editdist),('LCS',stringcmp.lcs)]:
        for canopy_method in [('threshold', 0.95, 0.85),
                              ('threshold', 0.9, 0.8)]:

          for this_index_def in index_def_list:

            strmap_index = indexing.StringMapIndex(desc='StringMap TH index:' \
                                      + ' index: dim=%d, sub-dim=%d, grid=' % \
                                      (dim, subdim)+'%d, str_cmp=%s, %s' % \
                                      (grid, str_cmp_funct[0], \
                                      str(canopy_method)),
                                      dataset1 = data_set1,
                                      dataset2 = data_set2,
                                      rec_comparator = rec_cmp,
                                      progress=progress_precentage,
                                      index_def = this_index_def,
                                      dim = dim,
                                      sub_dim= subdim,
                                      grid_resol =grid,
                                      sim_fu = str_cmp_funct[1],
                                      canopy_m=canopy_method)

            ds_index_list.append(['string-map-th', strmap_index])

  # StringMap based indexing (nearest neighbour based) - - - - - - - - - - - - -
  #
  for (dim, subdim) in [(15,3), (20,5)]:
    for grid in [100,1000]:
      for str_cmp_funct in [('Jaro',stringcmp.jaro),('Bigram',stringcmp.bigram),
                            ('ED',stringcmp.editdist),('LCS',stringcmp.lcs)]:
        for canopy_method in [('nearest', 5, 10),
                              ('nearest', 10, 20)]:

          for this_index_def in index_def_list:

            strmap_index = indexing.StringMapIndex(desc='StringMap NN index:' \
                                      + ' index: dim=%d, sub-dim=%d, grid=' % \
                                      (dim, subdim)+'%d, str_cmp=%s, %s' % \
                                      (grid, str_cmp_funct[0], \
                                      str(canopy_method)),
                                      dataset1 = data_set1,
                                      dataset2 = data_set2,
                                      rec_comparator = rec_cmp,
                                      progress=progress_precentage,
                                      index_def = this_index_def,
                                      dim = dim,
                                      sub_dim= subdim,
                                      grid_resol =grid,
                                      sim_fu = str_cmp_funct[1],
                                      canopy_m=canopy_method)

            ds_index_list.append(['string-map-nn', strmap_index])

  # ---------------------------------------------------------------------------
  # Run experiments for this data set
  #
  for (index_name, index_method) in ds_index_list:

    # Only do the experiments selected with command line argument, unless
    # argument was a '*' which means do all
    #
    if ((arg_index_method_name != 'all') and \
        (index_name != arg_index_method_name)):
      print
      print '  Skipping index:', index_name

      continue

    print
    print ' ', index_method.description

    time1 = time.time()
    index_method.build()

    time2 = time.time()
    index_method.compact()

    time3 = time.time()

    print '    Time used (in sec): build: %.4f,  compact: %.4f' % \
          (time2-time1, time3-time2)

    memory_usage_str = auxiliary.get_memory_usage()
    if (memory_usage_str != None):
      print '    '+memory_usage_str

    # Calculate complexity and quality measures - - - - - - - - - - - - - - - -
    #
    num_rec_pairs = index_method.num_rec_pairs
    A =  data_set1.num_records

    # Get the number of matches and non-matches in record pair dictionary
    #
    m =  0  # Number of matches
    nm = 0  # Number of non-matches

    for (rec_id1, rec_id2_set) in index_method.rec_pair_dict.iteritems():

      rec1 = index_method.rec_cache1[rec_id1]

      for rec_id2 in rec_id2_set:
        if (task == 'dedup'):
          rec2 = index_method.rec_cache1[rec_id2]  # From same data set
        else:
          rec2 = index_method.rec_cache2[rec_id2]

        if (check_match_funct(rec1, rec2) == True):
          m += 1
        else:
          nm += 1

    assert (m+nm) == num_rec_pairs

    # Get total number of matches possible (without indexing) from data set(s)
    #
    if (task == 'dedup'):

      ent_id_dict = {}  # Dictionary with all unique entity identifers and
                        # counts of how often they appear

      for ent_rec in index_method.rec_cache1.itervalues():
        ent_id = get_id_funct(ent_rec)
        ent_id_count = ent_id_dict.get(ent_id, 0) + 1
        ent_id_dict[ent_id] = ent_id_count

      assert sum(ent_id_dict.values()) == len(index_method.rec_cache1)

      tm = 0  # Total number of true matches (without indexing)

      for (ent_id, ent_count) in ent_id_dict.iteritems():
        tm += ent_count*(ent_count-1)/2

    else:  # A linkage

      ent_id_dict1 = {}  # Dictionaries with all unique entity identifers and
      ent_id_dict2 = {}  # counts of how often they appear

      for ent_rec in index_method.rec_cache1.itervalues():
        ent_id = get_id_funct(ent_rec)
        ent_id_count = ent_id_dict1.get(ent_id, 0) + 1
        ent_id_dict1[ent_id] = ent_id_count

      for ent_rec in index_method.rec_cache2.itervalues():
        ent_id = get_id_funct(ent_rec)
        ent_id_count = ent_id_dict2.get(ent_id, 0) + 1
        ent_id_dict2[ent_id] = ent_id_count

      assert sum(ent_id_dict1.values()) == len(index_method.rec_cache1)
      assert sum(ent_id_dict2.values()) == len(index_method.rec_cache2)

      tm = 0

      if (len(ent_id_dict1) < len(ent_id_dict2)):
        for (ent_id, ent_count) in ent_id_dict1.iteritems():
          if ent_id in ent_id_dict2:
            tm += ent_count*ent_id_dict2[ent_id]

      else:
        for (ent_id, ent_count) in ent_id_dict2.iteritems():
          if ent_id in ent_id_dict1:
            tm += ent_count*ent_id_dict1[ent_id]

    if (tm > 0):
      assert tm >= m
      pc = float(m) / float(tm)
    else:
      pc = -1

    # print '** Total number of true matches:', tm

    print '    Pairs completeness:  %.2f %%' % (pc*100.0)

    pq = float(m) / float(num_rec_pairs)
    print '    Pairs quality:       %.2f %%' % (pq*100.0)


    if (task == 'dedup'):
      rr = 1 - float(num_rec_pairs) / float(0.5*A*(A-1))

    else:  # A linkage
      B = data_set2.num_records
      rr = 1 - float(num_rec_pairs) / (float(A)*float(B))

    print '    Reduction ratio:     %.2f %%' % (rr*100.0)

    memo_use = auxiliary.get_memory_usage_val()
    print '    Memory usage:        %.2f MB' % (memo_use)

    # Write into results file - - - - - - - - - - - - - - - - - - - - - - - - -
    # (data set name, index method name, num_rec_pairs, RR, PC, PQ, time and
    # memory usage for for build and compact)
    #
    res_file_str = '%s, %s, %d, %.4f, %.4f, %.4f, %.4f, %.2f' % \
                   (ds_name, index_name, num_rec_pairs, rr, pc, pq, time3-time1,
                    memo_use)
    print '    Saved into results file: "%s"' % (res_file_str)
    res_file.write(res_file_str+os.linesep)

    del index_method  # Clean up memory

  print
  print

  del ds_index_list
  del data_set1
  del data_set2

print 'Finished experiments at:', time.ctime()

# =============================================================================
