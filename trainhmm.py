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
# The Original Software is: "trainhmm.py"
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

"""Module trainhmm.py - Main module to train a Hidden Markov Model (HMM)

   DESCRIPTION:
     This module can be used to train a Hidden Markov Model (HMM) using tagged
     data as created by a name or address standardiser from standardisation.py

   USAGE:
     python trainhmm.py [hmm_training_file] [hmm_output_file] [hmm_smoothing]

     with the following arguments:

       hmm_training_file  A text file written by the name or address
                          standardisation routine, and modified manually
                          (by adding HMM state sequence details to the provided
                          tag sequences).
       hmm_output_file    The trained hMM in textual representation (so it can
                          be read by the load() routine from simplehmm.py)
       hmm_smoothing      The HMM smoothing to be used, can be one of: 'none',
                          'laplace', or 'absdiscount'.

     The format of the input training file is as follows:
     - Comment lines must start with a hash character (#).
     - Empty lines are possible and are skipped.
     - Each training record consists of two non-commented lines:
       1) A comma separated list with tags (all uppercase and two characters
          long)
       2) The corresponding comma separated list of HMM states for this tag
          sequence. All HMM states must be lowercase.

     If a line starts with '#STOP' then the training process will be ended at
     this position in the file and the remainder will not be processed.

     The output file is a text file with all parameters of the HMM (see
     'simplehmm.py' module for more details).

     For more information on the smoothing methods see the 'simplehmm.py'
     module ('train' routine) or e.g.
       V.Borkar et.al., Automatic Segmentation of Text into Structured Records
       Section 2.2

   EXAMPLE:
     The following lines are two examples of the input training file format
     (given here are names):

     # Input string: peter marco jones
     # Token list:   ['peter', 'mark', 'jones']
     GM:gname1 ,GM:gname2  ,UN:sname1

     # Input string: alison de francesco
     # Token list:   ['alisa', 'de', 'francis']
     GF:gname1,PR:pref1,GM:sname1
"""

# =============================================================================
# Import necessary modules (Python standard modules first, then Febrl modules)

import sys
import time

import simplehmm

# =============================================================================
# Get command line arguments

if (len(sys.argv) != 4):
  print 'USAGE: python trainhmm.py  [hmm_training_file] [hmm_model_file] ' + \
        '[hmm_smoothing]'
  sys.exit()

hmm_training_file = sys.argv[1]
hmm_model_file =    sys.argv[2]
hmm_smoothing =     sys.argv[3].lower()

if (hmm_smoothing not in ['none', 'laplace','absdiscount']):
  print 'Illegal HMM smoothing method: %s ' % (hmm_smoothing) + \
        '(has to be one of: "none", "laplace","absdiscount")'
  sys.exit()

if (hmm_smoothing == 'none'):
  hmm_smoothing = None

# =============================================================================
# Load HMM training file, chek its structure and get all tags and HMM states

try:
  train_file = open(hmm_training_file, 'r')
except:
  logging.exception('Cannot open HMM training file: "%s"' % \
                    (hmm_training_file))
  raise IOError

tag_set =        set()
state_set =      set()
train_rec_list = []

line_cnt = 0

line = train_file.readline()
while (line != ''):
  line = line.strip()

  if (line.startswith('#STOP')):
    break  # End processing training records

  if ((line != '') and (line[0] != '#')):  # Line not empty and not a comment

    if ('#' in line):
      line = line[:line.find('#')].strip()  # Remove comment from end of line

    train_rec = []

    for tag_state_pair in line.split(','):
      tag, state = tag_state_pair.split(':')

      tag =   tag.strip()
      state = state.strip()

      if ((len(tag) != 2) or (tag.isupper() != True)):
        print 'Illegal tag in line %d: %s' % (line_cnt, line)
        sys.exit()
      tag_set.add(tag)

      if ((state == '') or (state.isupper() == True)):
        print 'Empty or illegal state in line %d: %s' % (line_cnt, line)
        sys.exit()
      state_set.add(state)

      train_rec.append((state, tag))

    train_rec_list.append(train_rec)

  line_cnt += 1
  line = train_file.readline()

train_file.close()

tag_list = list(tag_set)
tag_list.sort()
state_list = list(state_set)
state_list.sort()


print 'Set of tags found in HMM training file:'
print '  %s' % (', '.join(tag_list))
print
print 'Set of HMM states found in HMM training file:'
print '  %s' % (', '.join(state_list))
print

print 'Parsed %d training records:' % (len(train_rec_list))
for train_rec in train_rec_list:
  print '  %s' % (train_rec)
print

# Initalise HMM and train it with training data - - - - - - - - - - - - - - - -
#
hmm_name =   'Febrl HMM based on training file "%s"' % (hmm_training_file)
hmm_states = list(state_set)
hmm_observ = list(tag_set)

train_hmm = simplehmm.hmm(hmm_name, hmm_states, hmm_observ)

# Train, print and save the HMM - - - - - - - - - - - - - - - - - - - - - - - -
#
train_hmm.train(train_rec_list, hmm_smoothing)
train_hmm.print_hmm()

# Save trained HMM  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
train_hmm.save_hmm(hmm_model_file)

# =============================================================================
