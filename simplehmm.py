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
# The Original Software is: "simplehmm.py"
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

"""Module simplehmm.py - Routines for simple Hidden Markov Model (HMM)
                         functionality.

   DESCRIPTION:
     This module implements a class 'hmm' for simple Hidden Markov Models
     (HMM).

   PUBLIC FUNCTIONS:
     initialisation (__init__)  Initialise a new Hidden Markov Model
     set_trans_prob             Set transition probability
     set_obser_prob             Set observation probability
     set_init_prob              Set initial state probability
     check_prob                 Check probabilities in HMM for validity
     train                      Train the HMM with annotated training data
     viterbi                    Apply the Viterbi algorithm to get probability
                                of an observation sequence
     save_hmm                   Save the HMM into a text file
     load_hmm                   Load a HMM from a text file
     print_hmm                  Print a HMM

   See doc strings of individual functions for detailed documentation.

   TODO:
   - 6/12/2002: Add start and end states
"""

# =============================================================================
# Import necessary modules (Python standard modules first, then Febrl modules)

import logging
import os
import time

# =============================================================================

class hmm:
  """Routines for simple HMM functionality.
  """

  def __init__(self, description, state_list, obser_list,  inita_prob=None,
               trans_prob=None, obser_prob=None):
    """Initialise a new Hidden Markov Model.

    USAGE:
      myhmm = hmm(description, states, observ)

    ARGUMENTS:
      description  A description for the HMM
      state_list   List of states of the HMM
      obser_list   List of observations
      inita_prob   Initial state probabilities (default: None)
      trans_prob   Transition probabilities (default: None)
      obser_prob   Observation probabilities (default: None)

    DESCRIPTION:
      This routine initialises a new HMM data structure using the given state
      and observation list.

      Optional arguments are initial state probabilities, transition
      probabilities and observation probabilities. If they are not given
      (default) they are set to zero.
    """

    if (not isinstance(description, str)):
      logging.exception('Argument "description" is not a string')
      raise Exception

    if (not isinstance(state_list, list)):
      logging.exception('Argument "state_list" is not a list')
      raise Exception
    if (not isinstance(obser_list, list)):
      logging.exception('Argument "obser_list" is not a list')
      raise Exception

    self.description = description

    self.N = len(state_list)
    self.M = len(obser_list)
    self.S = state_list
    self.O = obser_list

    tmp = range(self.N)  # Temporary array for states
    for i in range(self.N):
      tmp[i] = 0.0

    # Matrix for transition probabilities - - - - - - - - - - - - - - - - - - -
    #
    self.A = []
    for i in range(self.N):  # One line per state
      self.A.append(tmp[:])
    if (trans_prob != None):  # Transition probabilities given as input
      if (not isinstance(trans_prob, list)):
        logging.exception('Argument "trans_prob" is not a list')
        raise Exception
      for i in range(self.N):
        for j in range(self.N):
          if (not isinstance(trans_prob[i][j],float)) or \
             (trans_prob[i][j] < 0.0) or (trans_prob[i][j] > 1.0):
            logging.exception('Argument "trans_prob" at index [%i,%i]' % \
                              (i,j) + ' is not a valid number between 0.0 ' + \
                              'and 1.0')
            raise Exception
          self.A[i][j] = trans_prob[i][j]

    # Vector for initial state probabilities  - - - - - - - - - - - - - - - - -
    #
    self.pi = tmp[:]
    if (inita_prob != None):  # Initial state probabilities given as input
      if (not isinstance(inita_prob, list)):
        logging.exception('Argument "inita_prob" is not a list')
        raise Exception
      for i in range(self.N):
        if (not isinstance(inita_prob[i],float)) or \
           (inita_prob[i] < 0.0) or (inita_prob[i] > 1.0):
          logging.exception('Argument "inita_prob" at index [%i] is not ' % \
                            (i) + 'a valid number between 0.0 and 1.0')
          raise Exception
        self.pi[i] = inita_prob[i]

    tmp = range(self.M)  # Temporary array for observations
    for i in range(self.M):
      tmp[i] = 0.0

    # Matrix for observation probabilities  - - - - - - - - - - - - - - - - - -
    #
    self.B = []
    for i in range(self.N):  # One line per state
      self.B.append(tmp[:])
    if (obser_prob != None):  # Observation probabilities given as input
      if (not isinstance(obser_prob, list)):
        logging.exception('Argument "obser_prob" is not a list')
        raise Exception
      for i in range(self.N):
        for j in range(self.M):
          if (not isinstance(obser_prob[i][j],float)) or \
             (obser_prob[i][j] < 0.0) or (obser_prob[i][j] > 1.0):
            logging.exception('Argument "obser_prob" at index [%i,%i]' % \
                              (i,j) + ' is not a valid number between 0.0 ' + \
                              'and 1.0')
            raise Exception
          self.B[i][j] = obser_prob[i][j]

    # Make dictionaries with indices of state and observation names - - - - - -
    #
    self.S_ind = {}
    for i in range(self.N):
      self.S_ind[self.S[i]] = i
    self.O_ind = {}
    for i in range(self.M):
      self.O_ind[self.O[i]] = i

    # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.info('Initialised HMM:')
    logging.info('  Description:  %s' % (self.description))
    logging.info('  States:       %s' % (str(state_list)))
    logging.info('  Observations: %s' % (str(obser_list)))

  # ---------------------------------------------------------------------------

  def set_trans_prob(self, from_state, to_state, trans_prob):
    """Set transition probability from 'from_state' to 'to_state'

    USAGE:
      myhmm.set_trans_prob(from_state, to_state, trans_prob)

    ARGUMENTS:
      from_state  From state (must be in list of states of the HMM)
      to_state    To state (must be in list of states of the HMM)
      trans_prob  Corresponding transition probability to be set

    DESCRIPTION:
      Sets the transition probability from 'from_state' to 'to_state'
      to the given value.
    """

    if (from_state not in self.S):
      logging.exception('Illegal "from" state: %s' % (str(from_state)))
      raise Exception

    if (to_state not in self.S):
      logging.exception('Illegal "to" state: %s' % (str(to_state)))
      raise Exception

    if (not isinstance(trans_prob,float)) or (trans_prob < 0.0) or \
       (trans_prob > 1.0):
      logging.exception('Argument "trans_prob" is not a valid number ' + \
                        'between 0.0 and 1.0')
      raise Exception

    self.A[self.S_ind[from_state]][self.S_ind[to_state]] = trans_prob

  # ---------------------------------------------------------------------------

  def set_obser_prob(self, state, obser, obser_prob):
    """Set observation probability for 'obser' in 'state'.

    USAGE:
      myhmm.set_obser_prob(state, obser, obser_prob)

    ARGUMENTS:
      state       The state (must be in list of states of the HMM)
      to_state    The observation (must be in list of observations of the HMM)
      obser_prob  Corresponding observation probability to be set (observation
                  must be in the list of observations of the HMM)

    DESCRIPTION:
      For state 'state' sets the probability for observation 'obser' to the
      given value.
     """

    if (state not in self.S):
      logging.exception('Illegal state: %s ' (str(state)))
      raise Exception

    if (obser not in self.O):
      logging.exception('Illegal observation: %s' % (str(obser)))
      raise Exception

    if (not isinstance(obser_prob,float)) or (obser_prob < 0.0) or \
       (obser_prob > 1.0):
      logging.exception('Argument "obser_prob" is not a valid number ' + \
                        'between 0.0 and 1.0')
      raise Exception

    self.B[self.S_ind[state]][self.O_ind[obser]] = obser_prob

  # ---------------------------------------------------------------------------

  def set_init_prob(self, state, init_prob):
    """Set initial state probability for 'state'.

    USAGE:
      myhmm.set_init_prob(state, init_prob)

    ARGUMENTS:
      state      The state (must be in list of states of the HMM)
      init_prob  Corresponding initial probability to be set

    DESCRIPTION:
      For state 'state' sets the initial state probability to the given value.
    """

    if (state not in self.S):
      logging.exception('Illegal state: %s' % (str(state)))
      raise Exception

    if (not isinstance(init_prob,float)) or (init_prob < 0.0) or \
       (init_prob > 1.0):
      logging.exception('Argument "init_prob" is not a valid number ' + \
                        'between 0.0 and 1.0')
      raise Exception

    self.pi[self.S_ind[state]] = init_prob

  # ---------------------------------------------------------------------------

  def check_prob(self):
    """Check probabilities in HMM for validity.

    USAGE:
      myhmm.check_prob()

    ARGUMENTS:
      None

    DESCRIPTION:
      Checks all probabilities in the HMM for validity, i.e. if they sum to
      1.0 in each state (observation and outgoing state probabilities).

      If an error occurs, a negative number is returned, 0 otherwise
    """

    delta = 0.0000000000001  # Account for floating-point rounding errors

    ret = 0

    sum = 0.0
    for i in range(self.N):
      sum = sum+self.pi[i]
    if (abs(sum - 1.0) > delta):
      logging.warn('HMM initial state probabilities sum is not 1: %f' % (sum))
      ret -= 1

    for i in range(self.N):
      sum  = 0.0
      for j in range(self.N):
        sum = sum+self.A[i][j]
      if (abs(sum - 1.0) > delta):
        logging.warn('HMM state "%s" has transition ' % (self.S[i])  + \
                     'probabilities sum not 1.0: %f' % (sum))
        ret -= 1

    for i in range(self.N):
      sum  = 0.0
      for j in range(self.M):
        sum = sum+self.B[i][j]
      if (abs(sum - 1.0) > delta):
        logging.warn('HMM state "%s" has observation ' % (self.S[i]) + \
                     'probabilities sum not 1.0: '+str(sum))
        ret -= 1

    return ret

  # ---------------------------------------------------------------------------

  def train(self, train_data, smoothing=None):
    """Train the HMM with annotated training data (supervised learning).

    USAGE:
      myhmm.train(train_data)

    ARGUMENTS:
      train_data  A set of training data in list form
      smoothing   If smoothing of the observation probabilities is desired (for
                  unknown symbols) then ths argument should be set to either
                  'laplace' (for Laplace smoothing) or 'absdiscount' (for the
                  absolute discounting method).
                  The default is 'None' and no smoothing will be done.

    DESCRIPTION:
      Using the training data, the HMM probabilities are set.

      'train_data' is a Python list with one element per training record.
      Each training record is a list (sequence) with pairs (tuple or list)
      of (state,observation) pairs. These training records can have varying
      length.

      For more information on the smoothing methods, see e.g.
        V.Borkar et.al., Automatic Segmentation of Text into
                         Structured Records
        Section 2.2
    """

    if (smoothing not in [None, 'laplace', 'absdiscount']):
      logging.exception('Illegal value for "smoothing" argument: %s' % \
                        (smoothing)+ ', possible are: None, "laplace" or ' + \
                        '"absdiscount"')
      raise Exception

    # Reset initial state, transition and observation probabilities - - - - - -
    #
    for i in range(self.N):
      self.pi[i] = 0.0
      for j in range(self.N):
        self.A[i][j] = 0.0
      for j in range(self.M):
        self.B[i][j] = 0.0

    # Sum up probabilities from training data - - - - - - - - - - - - - - - - -
    #
    for train_rec in train_data:
      (init_state, init_obser) = train_rec[0]  # Get first tuple
      i = self.S_ind[init_state]
      self.pi[i] = self.pi[i] + 1.0
      prev_i = -1  # No previous index yet

      for pair in train_rec:  # For each pair in this training record
        (state,obser) = pair
        i = self.S_ind[state]
        j = self.O_ind[obser]
        self.B[i][j] = self.B[i][j] + 1.0
        if (prev_i > -1):  # If previous index defined:
          self.A[prev_i][i] = self.A[prev_i][i] + 1.0
        prev_i = i

    # Scale counts into probabilities - - - - - - - - - - - - - - - - - - - - -
    #

    # Scale initial probabilities
    #
    sum = 0.0
    for i in range(self.N):
      sum = sum+self.pi[i]
    if (sum != 0.0):
      for i in range(self.N):
        self.pi[i] = self.pi[i] / sum

    # Scale transition probabilities
    #
    for i in range(self.N):  # For each state
      sum  = 0.0
      for j in range(self.N):
        sum = sum+self.A[i][j]
      if (sum != 0.0):
        for j in range(self.N):
          self.A[i][j] = self.A[i][j] / sum

    # Scale observation probabilities
    #
    if (smoothing == None):  # No smoothing to be done
      for i in range(self.N):  # For each state
        sum  = 0.0
        for j in range(self.M):
          sum = sum+self.B[i][j]
        if (sum != 0.0):
          for j in range(self.M):
            self.B[i][j] = self.B[i][j] / sum

    elif (smoothing == 'laplace'): # Do Laplace smoothing
      for i in range(self.N):  # For each state
        sum  = float(self.M)
        for j in range(self.M):
          sum = sum+self.B[i][j]
        for j in range(self.M):
          self.B[i][j] = (self.B[i][j]+1.0) / sum

    elif (smoothing == 'absdiscount'):  # Do absolute discounting smoothing
      for i in range(self.N):  # For each state
        sum  = 0.0
        mi = 0  # Number of distinct symbols seen in state i

        for j in range(self.M):
          if (self.B[i][j] != 0):  # Symbol has been counted during training
            sum = sum+self.B[i][j]
            mi += 1

        x = 1.0 / float(sum+self.M)

        if (sum != 0.0):
          for j in range(self.M):
            if (self.B[i][j] != 0.0):  # A known observation symbol
              self.B[i][j] = (self.B[i][j] / sum) - x
            else:  # An unknown observation symbol
              self.B[i][j] = float(mi * x) / float(self.M-mi)

    self.check_prob()  # Check if probabilities are OK

    # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.info('Trained HMM with %i training records' % (len(train_data)))
    logging.info('  Smooting technique used: %s' % (str(smoothing)))

  # ---------------------------------------------------------------------------

  def viterbi(self, obser_seq):
    """Apply the Viterbi algorithm.

    USAGE:
      [sequence, seq_prob] = myhmm.viterbi(obser_seq)

    ARGUMENTS:
      obser_seq  An observation sequence (all observations must be in the list
                 of observations of the HMM)

    DESCRIPTION:
      This routine uses the Viterbi algorithm to find the most likely state
      sequence for the given observation sequence. Returns the sequence in a
      list and it's probability.
    """

    tmp = range(self.N)  # Temporary array with zeros
    for i in range(self.N):
      tmp[i] = 0

    obs_len = len(obser_seq)
    obs_ind = []
    for obs in obser_seq:
      obs_ind.append(self.O_ind[obs])

    delta = [tmp[:]]  # Compute initial state probabilities
    for i in range(self.N):
      delta[0][i] = self.pi[i] * self.B[i][obs_ind[0]]

    phi = [tmp[:]]

    for obs in obs_ind[1:]:  # For all observations except the inital one
      delta_t = tmp[:]
      phi_t = tmp[:]
      for j in range(self.N):   # Following formula 33 in Rabiner'89
        tdelta = tmp[:]
        tphimax = -1.0
        for i in range(self.N):
          tphi_tmp = delta[-1][i] * self.A[i][j]
          if (tphi_tmp > tphimax):
            tphimax = tphi_tmp
            phi_t[j] = i
          tdelta[i] = tphi_tmp * self.B[j][obs]
        delta_t[j] = max(tdelta)
      delta.append(delta_t)
      phi.append(phi_t)

    # Backtrack the path through the states  (Formula 34 in Rabiner'89)
    #
    tmax = -1.0
    for i in range(self.N):
      if (delta[-1][i] > tmax):
        tmax = delta[-1][i]
        state_seq = [i]  # Last state with maximum probability

    phi.reverse()  # Because we start from the end of the sequence
    for tphi in phi[:-1]:
      state_seq.append(tphi[state_seq[-1]])

    sequence = []

    for state in state_seq:
      sequence.append(self.S[state])

    sequence.reverse()  # Reverse into correct time direction
    state_seq.reverse()

    # Finally compute probability of this state and observation sequence
    #
    prev_ind = state_seq[0]
    seq_prob = self.pi[prev_ind]
    seq_prob *= self.B[prev_ind][self.O_ind[obser_seq[0]]]

    for i in range(1,len(state_seq)):
      ind = state_seq[i]
      obs = self.O_ind[obser_seq[i]]
      seq_prob *= self.A[prev_ind][ind]
      seq_prob *= self.B[ind][obs]
      prev_ind = ind

    # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.debug('  Viterbi analysis')
    logging.debug('    Input observation sequence: %s' % (str(obser_seq)))
    logging.debug('    Output state sequence:      %s' % (str(sequence)))
    logging.debug('    Output probability: %f' % (seq_prob))

    return [sequence, seq_prob]

  # ---------------------------------------------------------------------------

  def save_hmm(self, file_name):
    """Save the HMM into a text file.

       USAGE:
         myhmm.save_hmm(file_name)

       ARGUMENTS:
         file_name  The name of the text file into which the HMM is written

       DESCRIPTION:
         Writes comments into the text file in Python style, i.e. lines
         beginning with a hash character '#'.

         If the file name does not end with "hmm" and does not contain a "."
         then ".hmm" will be appended to the file name.
    """

    if ((not file_name.endswith('hmm')) and ('.' not in file_name)):
      file_name = file_name+'.hmm'

    # Open file for writing - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    try:
      f = open(file_name, 'w')
    except:
      logging.exception('Cannot write to file: %s' % (file_name))
      raise IOError

    # Write header  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    f.write("# Hidden Markov Model written by 'simplehmm.py'"+ \
            os.linesep)
    f.write("#"+os.linesep)
    f.write("# Created "+time.ctime(time.time())+os.linesep)
    f.write("#"+os.linesep)
    f.write("# File name: "+file_name+os.linesep)
    f.write('#'+'-'*70+os.linesep)
    f.write(os.linesep)

    # Write HMM description - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    f.write("# HMM description"+os.linesep)
    f.write("#"+os.linesep)
    f.write(self.description+os.linesep)
    f.write(os.linesep)

    # Write HMM states  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    f.write("# HMM states"+os.linesep)
    f.write("#"+os.linesep)
    state_str = ', '.join(self.S)
    f.write(state_str+os.linesep)
    f.write(os.linesep)

    # Write HMM observations  - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    f.write("# HMM observations"+os.linesep)
    f.write("#"+os.linesep)
    obser_str = ', '.join(self.O)
    f.write(obser_str+os.linesep)
    f.write(os.linesep)

    # Write HMM initial probabilities - - - - - - - - - - - - - - - - - - - - -
    #
    f.write("# HMM initial probabilities"+os.linesep)
    f.write("#"+os.linesep)
    inital_prob_str = ''
    for i in range(self.N):
      inital_prob_str += '%f' %(self.pi[i])+', '
    f.write(inital_prob_str[:-2]+os.linesep)  # Strip of trailing ', '
    f.write(os.linesep)

    # Write HMM transition probabilities  - - - - - - - - - - - - - - - - - - -
    #
    f.write("# HMM transition probabilities (from state row-wise)"+os.linesep)
    f.write("#"+os.linesep)
    for i in range(self.N):
      trans_prob_str = ''
      for j in range(self.N):
        trans_prob_str +=  '%f' %(self.A[i][j])+', '
      f.write(trans_prob_str[:-2]+os.linesep) # Strip of trailing ', '
    f.write(os.linesep)

    # Write HMM observation probabilities - - - - - - - - - - - - - - - - - - -
    #
    f.write("# HMM observation probabilities (state row-wise)"+os.linesep)
    f.write("#"+os.linesep)
    for i in range(self.N):
      obser_prob_str = ''
      for j in range(self.M):
        obser_prob_str +=  '%f' %(self.B[i][j])+', '
      f.write(obser_prob_str[:-2]+os.linesep) # Strip of trailing ', '
    f.write(os.linesep)

    f.close()

    # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.info('HMM save to file: %s' % (file_name))

  # ---------------------------------------------------------------------------

  def load_hmm(self, file_name):
    """Load a HMM from a text file.

    USAGE:
      myhmm.load_hmm(file_name)

    ARGUMENTS:
      file_name  The name of the text file from which the HMM is read.

    DESCRIPTION:
      This routine reads a HMM from a text file as written by the 'save_hmm'
      routine (see above).

      It is assumed the the HMM data structure is allocated. All elements of
      the HMM are overwritten with the newly loaded values.

      Empty lines and comment lines (starting with a '#') are skipped over.
    """

    # Open file for reading - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    try:
      f = open(file_name, 'r')
    except:
      logging.exception('Cannot read from file: %s' % (file_name))
      raise IOError

    # Skip over header and empty lines  - - - - - - - - - - - - - - - - - - - -
    #
    line = f.readline()
    while (line[0] == '#') or (line.strip() == ''):
      line = f.readline()

    # Read HMM description - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    self.description = line.strip()  # Remove line separators

    # Skip over empty line(s) and comments for states
    #
    line = f.readline()
    while (line[0] == '#') or (line.strip() == ''):
      line = f.readline()

    # Read HMM states - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    line = line.strip()  # Remove line separators
    states = line.split(', ')
    self.N = len(states)
    self.S = states

    # Skip over empty line(s) and comments for observations
    #
    line = f.readline()
    while (line[0] == '#') or (line.strip() == ''):
      line = f.readline()

    # Read HMM observations - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    line = line.strip()  # Remove line separators
    observations = line.split(', ')
    self.M = len(observations)
    self.O = observations

    # Skip over empty line(s) and comments for initial probabilities
    #
    line = f.readline()
    while (line[0] == '#') or (line.strip() == ''):
      line = f.readline()

    # Read HMM initial probabilities  - - - - - - - - - - - - - - - - - - - - -
    #
    line = line.strip()  # Remove line separators
    init_probs = line.split(', ')

    if (len(init_probs) != self.N):
      logging.exception('Illegal number of initial probabilities in file')
      raise Exception

    self.pi = []
    for init_prob in init_probs:
      self.pi.append(float(init_prob))

    # Skip over empty line(s) and comments for transition probabilities
    #
    line = f.readline()
    while (line[0] == '#') or (line.strip() == ''):
      line = f.readline()

    # Read HMM transition probabilities  - - - - - - - - - - - - - - - - - - -
    #
    self.A = []
    for i in range(self.N):  # Read N lines with data
      line = line.strip()  # Remove line separators
      trans_probs = line.split(', ')

      if (len(trans_probs) != self.N):
        logging.exception('Illegal number of transition probabilities in file')
        raise Exception

      tmp_trans_prob = []
      for trans_prob in trans_probs:
        tmp_trans_prob.append(float(trans_prob))
      self.A.append(tmp_trans_prob)
      line = f.readline()  # Read next line with tansition probabilities

    # Skip over empty line(s) and comments for observation probabilities
    #
    line = f.readline()
    while (line[0] == '#') or (line.strip() == ''):
      line = f.readline()

    # Read HMM observation probabilities  - - - - - - - - - - - - - - - - - - -
    #
    self.B = []
    for i in range(self.N):  # Read N lines with data
      line = line.strip()  # Remove line separators
      obser_probs = line.split(', ')

      if (len(obser_probs) != self.M):
        logging.exception('Illegal number of observation probabilities in' + \
                          ' file')
        raise Exception

      tmp_obser_prob = []
      for obser_prob in obser_probs:
        tmp_obser_prob.append(float(obser_prob))
      self.B.append(tmp_obser_prob)
      line = f.readline()  # Read next line with tansition probabilities

    f.close()

    # Make dictionaries with indices of state and observation names - - - - - -
    #
    self.S_ind = {}
    for i in range(self.N):
      self.S_ind[self.S[i]] = i
    self.O_ind = {}
    for i in range(self.M):
      self.O_ind[self.O[i]] = i

    # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.info('HMM loaded from file: %s' % (file_name))

  # ---------------------------------------------------------------------------

  def print_hmm(self):
    """Print a HMM.

    USAGE:
      myhmm.print_hmm()

    ARGUMENTS:
      None

    DESCRIPTION:
      Prints a HMM with all its parameters. Only probabilities with values
      larger than 0.0 are printed.
    """

    msg = []  # Compose a message

    state_list = self.S[:]  # Make a copy
    state_list.sort()
    obser_list = self.O[:]  # Make a copy
    obser_list.sort()

    msg.append('Hidden Markov Model')
    msg.append('  Description:  %s' % (str(self.description)))
    msg.append('  States:       %s' % (str(state_list)))
    msg.append('  Observations: %s' % (str(obser_list)))

    msg.append('')
    msg.append('  Inital state probabilities:')
    for i in range(self.N):
      if (self.pi[i] > 0.0):
        msg.append('    State: '+self.S[i]+': '+str(self.pi[i]))

    msg.append('')
    msg.append('  Transition probabilities:')
    for i in range(self.N):
      msg.append('    From state: '+self.S[i])
      for j in range(self.N):
        if (self.A[i][j] > 0.0):
          msg.append('      to state: '+self.S[j]+': '+str(self.A[i][j]))

    msg.append('')
    msg.append('  Observation probabilities:')
    for i in range(self.N):
      msg.append('    In state: '+self.S[i])
      for j in range(self.M):
        if (self.B[i][j] > 0.0):
          msg.append('      Symbol: '+self.O[j]+': '+str(self.B[i][j]))

    msg.append('')

    for line in msg:
      print line
      logging.info(line)

# =============================================================================
