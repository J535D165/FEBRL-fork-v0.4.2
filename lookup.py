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
# The Original Software is: "lookup.py"
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

"""Module lookup.py - Classes for various types of look-up tables.

   This module contains classes for look-up table and correction lists.
"""

# =============================================================================
# Import necessary modules (Python standard modules first, then Febrl modules)

import logging

import auxiliary

# =============================================================================

class LookupTable(dict):
  """class LookupTable - Based on dictionary type.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor, set general attributes.
    """

    dict.__init__(self)  # Initialise dictionary base type
    self.description =  ''
    self.created =      ''
    self.modified =     ''
    self.file_names =   []
    self.default =      None  # Default return value for non existing keys
    self.length =       None  # Number of entries in the look-up table

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('desc')):
        auxiliary.check_is_string('description', value)
        self.description = value

      elif (keyword.startswith('defau')):
        self.default = value

      elif (keyword.startswith('creat')):
        self.created = value
      elif (keyword.startswith('modif')):
        self.modified = value

      else:
        logging.exception('Illegal constructor argument keyword: "%s"' % \
              (str(keyword)))
        raise Exception

  # ---------------------------------------------------------------------------

  def __getitem__(self, key):
    """Return an item in the look-up table with the given key. If not found,
       return the default value.
    """

    try:
      return dict.__getitem__(self, key)
    except KeyError:
      return self.default

  # ---------------------------------------------------------------------------

  def get(self, key, *args):
    """Return an item in the look-up table with the given key. If not found,
       return the default value.
    """

    if (not args):
      args = (self.default,)
    return dict.get(self, key, *args)

  # ---------------------------------------------------------------------------

  def load(self, file_names):
    """Load one or more files into the look-up table.
       See implementations in derived classes for details.
    """

    logging.exception('Override abstract method in derived class')
    raise Exception

# =============================================================================

class TagLookupTable(LookupTable):
  """A look-up table class for look-up tables with word corrections and tags.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):

    LookupTable.__init__(self, **kwargs)  # Initialise base class

    self.max_key_length = None  # The maximum length of a key in words

    # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.info('Initialised tag look-up table "%s"' %(self.description))
    logging.info('  With default: "%s"' % (str(self.default)))

  # ---------------------------------------------------------------------------

  def load(self, file_names):
    """Load one or more files with word corrections and tags into the look-up
       table.

       See Febrl manual for details on the file format.
    """

    # Check input argument type - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (isinstance(file_names, str)):
      file_names = [file_names]  # Make a list out of a single file name

    auxiliary.check_is_list('file_names', file_names)

    i = 0
    for file_name in file_names:
      auxiliary.check_is_string('file_name[%d]' % (i), file_name[i])
      i += 1

    self.file_names = file_names
    self.clear()  # Remove all items from the look-up table
    self.max_key_length = 0

    # Loop over file names - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    for fn in self.file_names:

      try:  # Open file and read all lines into a list
        f = open(fn,'r')
      except:
        logging.exception('Cannot read from file "%s"' % (fn))
        raise IOError

      file_data = f.readlines()  # Read complete file
      f.close()

      tag = ''  # Start with no tag
      key = ''  # Start with an empty key

      # Now process all lines - - - - - - - - - - - - - - - - - - - - - - - - -
      #
      for line in file_data:
        l = line.strip()  # Remove line separators
        if (len(l) > 0) and (l[0] != '#'):  # Not empty line and not comment

          if (l[:5] == 'tag=<'):  # It's a line with a new tag
            tag = l[5:7]

          else:  # A line with an entry

            # Make sure a tag is set
            #
            if (tag == ''):
              logging.exception('Missing tag specification in file "%s"' % \
                                (fn))
              raise Exception

            line_list = l.split(':')  # Separate key from values

            if (len(line_list) > 2):
              logging.exception('Illegal format in file "%s" in line: %s' % \
                                (fn, l))
              raise Exception

            if (len(line_list) == 2):  # Line contains a key - - - - - - - - -

              key = line_list[0].strip().lower() # Get and clean key

              key_list = key.split(' ')  # Make a list of key words
              if (len(key_list) > self.max_key_length):
                self.max_key_length = len(key_list) # Update maximal key length

              # Insert key itself into lookup table
              #
              dict_val = '_'.join(key_list)
              dict_key = tuple(key_list)
              this_tag = tag

              if (self.__contains__(dict_key)):  # Already in lookup table
                test_item = self.__getitem__(dict_key)
                test_val = test_item[0]  # Value without tag
                test_tag = test_item[1]

                if (dict_val != test_val):
                  logging.warn('Key "%s" already in dictionary with ' % \
                          (str(dict_val)) + 'different value (old value ' + \
                          'will be over written with "%s")' % (str(test_val)))

                if (test_tag.find(this_tag) < 0):  # This tag is new
                  this_tag = test_tag+'/'+this_tag  # Tag for this entry
                else:
                  this_tag = test_tag

              this_val = (dict_val, this_tag)
              self.__setitem__(dict_key,this_val)  # Insert key itself

              vals = line_list[1].lower() # Get values in this line in a string

            elif (len(line_list) == 1):  # Line contains only values - - - - -

              vals = line_list[0].lower() # Get values in this line in a string

            # Porcess all values right of ':' in this line

            val_list = vals.split(',')  # Split values into a list

            for val in val_list:  # Loop over all values  - - - - - - - - - - -

              val_strip = val.strip()

              if (val_strip != ''):  # Only append non-empty values
                key_list = val_strip.split(' ')  # Make a list of key words
                if (len(key_list) >  self.max_key_length):
                  self.max_key_length = len(key_list) # Update maximal key len

                dict_key = tuple(key_list)
                this_tag = tag

                if (self.__contains__(dict_key)):
                  test_item = self.__getitem__(dict_key)
                  test_val = test_item[0]  # Value without tag
                  test_tag = test_item[1]

                  if (dict_val != test_val):
                    logging.warn('Key "%s" already in dictionary with ' % \
                          (str(dict_val)) + 'different value (old value ' + \
                          'will be over written with "%s")' % (str(test_val)))

                  if (test_tag.find(this_tag) < 0):  # This tag is new
                    this_tag = test_tag+'/'+this_tag  # Tag for this entry
                  else:
                    this_tag = test_tag

                this_val = (dict_val, this_tag)
                self.__setitem__(dict_key,this_val)

    self.length = self.__len__()  # Get number of elements in the look-up table

    # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.info('Loaded tag look-up table "%s"' % (self.description))
    logging.info('  From files:         %s' % (str(self.file_names)))
    logging.info('  Number of entries:  %i' % (self.length))
    logging.info('  Maximal key length: %i' % (self.max_key_length))

# =============================================================================

class FrequencyLookupTable(LookupTable):
  """A look-up table class for look-up tables with words and frequencies.
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):

    LookupTable.__init__(self, **kwargs)  # Initialise base class

    self.sum = None  # The sum of all frequency counts
    self.default = 1

    # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.info('Initialised frequency look-up table "%s"' % \
                 (self.description))
    logging.info('  With default: %s' % (str(self.default)))

  # ---------------------------------------------------------------------------

  def load(self, file_names):
    """Load one or more files with words and their frequency counts into the
       look-up table.

       See Febrl manual for details on the file format.
    """

    # Check input argument type - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (isinstance(file_names, str)):
      file_names = [file_names]  # Make a list out of a single file name

    auxiliary.check_is_list('file_names', file_names)

    i = 0
    for file_name in file_names:
      auxiliary.check_is_string('file_name[%d]' % (i), file_name[i])
      i += 1

    self.file_names = file_names
    self.clear()  # Remove all items from the look-up table
    self.sum = 0

    # Loop over file names - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    for fn in self.file_names:

      try:  # Open file and read all lines into a list
        f = open(fn,'r')
      except:
        logging.exception('Cannot read from file "%s"' % (fn))
        raise IOError

      file_data = f.readlines()  # Read complete file
      f.close()

      # Now process all lines - - - - - - - - - - - - - - - - - - - - - - - - -
      #
      for line in file_data:
        l = line.strip()
        if (len(l) > 0) and (l[0] != '#'):  # Not empty line and not comment

          ll = l.split(',')  # Get fields from a line

          # Check for two columns
          #
          if (len(ll) != 2):
            logging.exception('Illegal file format (not 2 columns) in file' + \
                              ': "%s" in line: %s"' % (fn, l))
            raise Exception

          key = ll[0].strip().lower()  # Make sure it's lower case
          val = ll[1].strip().lower()

          try:
            val = int(val)  # Convert the value into an integer
          except:
            logging.exception('Illegal value for frequency count: "%s"' % \
                  (str(val)) + ' in line: "%s"' % (l))
            raise Exception

          if (self.__contains__(key)):
            val += self.__getitem__(key)  # Sum up counts

          self.__setitem__(key, val)
          self.sum += val

    self.length = self.__len__()  # Get number of elements in the look-up table

    # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.info('Loaded frequency look-up table "%s"' % \
                 (self.description))
    logging.info('  From files:        %s' % (str(self.file_names)))
    logging.info('  Number of entries: %i' % (self.length))
    logging.info('  Sum of all value:  %i' % (self.sum))

# =============================================================================

class GeocodeLookupTable(LookupTable):
  """A look-up table class for look-up tables with entires and their locations.

     For each entry in the look-up table, its longitude and latitude are given.
     The file format is three columns comma separated text file (CSV).

     For each entry, the key is the name of the locality, and the value is a
     list with the two entries [longitude,latitude].
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):

    LookupTable.__init__(self, **kwargs)  # Initialise base class

    self.default = []

    # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.info('Initialised geocode look-up table "%s"' % \
                 (self.description))
    logging.info('  With default: %s' % (str(self.default)))

  # ---------------------------------------------------------------------------

  def load(self, file_names):
    """Load one or more files with entries and their localities into the
       table.

       See Febrl manual for details on the file format.
    """

    # Check input argument type - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (isinstance(file_names, str)):
      file_names = [file_names]  # Make a list out of a single file name

    auxiliary.check_is_list('file_names', file_names)

    i = 0
    for file_name in file_names:
      auxiliary.check_is_string('file_name[%d]' % (i), file_name[i])
      i += 1

    self.file_names = file_names
    self.clear()  # Remove all items from the look-up table

    # Loop over file names - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    for fn in self.file_names:

      try:  # Open file and read all lines into a list
        f = open(fn,'r')
      except:
        logging.exception('Cannot read from file "%s"' % (fn))
        raise IOError

      file_data = f.readlines()  # Read complete file
      f.close()

      # Now process all lines - - - - - - - - - - - - - - - - - - - - - - - - -
      #
      for line in file_data:
        l = line.strip()
        if (len(l) > 0) and (l[0] != '#'):  # Not empty line and not comment

          ll = l.split(',')  # Get fields from a line

          # Check for three columns
          #
          if (len(ll) != 3):
            logging.exception('Illegal file format (not 3 columns) in file' + \
                              ': "%s" in line: %s' % (fn, l))
            raise Exception

          key = ll[0].strip().lower()  # Make sure it's lower case
          long = ll[1].strip()
          lati = ll[2].strip()

          # Try to convert into numerical (float) values
          #
          try:
            long = float(long)
          except:
            logging.exception('Longitude: "%s" is not a number in line: "%s"' \
                              % (str(long), l))
            raise Exception
          try:
            lati = float(lati)
          except:
            logging.exception('Lattitude: "%s" is not a number in line: "%s"' \
                              % (str(lati), l))
            raise Exception

          # And check their values
          #
          if (long < -180.0) or (long > 180.0):
            logging.exception('Illegal value for longitude: '+str(long))
            raise Exception
          if (lati < -90.0) or (lati > 90.0):
            logging.exception('Illegal value for latitude: '+str(lati))
            raise Exception

          val = [long,lati]  # Value for dictionary

          if (self.__contains__(key)) and (self.__getitem__(key) != val):
            logging.exception('Key "%s" already in look-up table with ' % \
                  (str(key)) + 'different value')
            raise Exception

          self.__setitem__(key, val)

    self.length = self.__len__()  # Get number of elements in the look-up table

    # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.info('Loaded geocode look-up table "%s"' % (self.description))
    logging.info('  From files:        %s' % (str(self.file_names)))
    logging.info('  Number of entries: %i' % (self.length))

# =============================================================================

class CorrectionList(list):
  """A class for correction lists (containing original and replacement
     strings).
  """

  # ---------------------------------------------------------------------------

  def __init__(self, **kwargs):
    """Constructor, set general attributes.
    """

    list.__init__(self)  # Initialise list base type
    self.description =  ''
    self.created =      ''
    self.modified =     ''
    self.file_name =    ''
    self.length =       None  # Number of entries in the correction list

    for (keyword, value) in kwargs.items():

      if (keyword.startswith('desc')):
        auxiliary.check_is_string('description', value)
        self.description = value

      elif (keyword.startswith('creat')):
        self.created = value
      elif (keyword.startswith('modif')):
        self.modified = value

      else:
        logging.exception('Illegal constructor argument keyword: "%s"' % \
              (str(keyword)))
        raise Exception

    # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.info('Initialised correction list "%s"' % (self.description))

  # ---------------------------------------------------------------------------

  def load(self, file_name):
    """Load one correction list file into a sorted (decreasing length) list.

       See Febrl manual for details on the file format.
    """

    # Check input argument type and open file - - - - - - - - - - - - - - - - -
    #
    auxiliary.check_is_string('file_name', file_name)

    self.file_name = file_name

    # Make sure the list is empty, remove all items from the correction list
    #
    while (self.__len__() > 0):
      self.pop()

    try:  # Open file and read all lines into a list
      f = open(self.file_name, 'r')
    except:
      logging.exception('Cannot read from file "%s"' % (str(self.file_name)))
      raise IOError

    file_data = f.readlines()  # Read complete file
    f.close()

    org_list  = []  # List of original strings (the ones to be replaced)
    repl_list = []  # List of replacement strings
    len_list  = []  # List of original string lengths
    repl = ''       # Set inital replacement to nothing

    # Now process all lines - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    for line in file_data:
      l = line.strip()  # Remove line separators at the end

      if (len(l) > 0) and (l[0] != '#'):  # Not an empty line and not comment
        ll = l.split(':=')  # Separate replacement from values

        if (len(ll) == 2):  # Line contains a replacement - - - - - - - - - - -

          repl = ll[0].strip().lower()  # Make replacement lower and strip

          if (not ((repl[0] == '"') and (repl[-1] == '"') or \
                   (repl[0] == "'") and (repl[-1] == "'"))):
            logging.exception('Replacement string is not properly quoted: '+ \
                  '"%s" in file: "%s"' % (repl, str(self.file_name)))
            raise Exception

          repl = repl[1:-1]  # Remove quotes from replacement string

          v = ll[1].lower() # Get values in a string and make lowercase

        elif (len(ll) == 1):  # Line contains only values - - - - - - - - - - -
          v = ll[0].lower() # Get values in a string and make lowercase

        else:  # More than one ':=' separator in the line - - - - - - - - - - -
          logging.exception('Too many ":=" separators in line: "%s"' % (l))
          raise Exception

        # Now process the values and append them to the list  - - - - - - - - -

        vv = v.split(',')  # Split values into a list

        for v in vv:  # Loop over all values  - - - - - - - - - - - - - - - -
          org = v.strip()  # Get the original string

          if (org != ''):  # Only process non-empty values
            if (not ((org[0] == '"') and (org[-1] == '"') or \
                     (org[0] == "'") and (org[-1] == "'"))):
              logging.exception('Original string is not properly quoted: '+ \
                    '"%s" in file: "%s"' % (org, str(self.file_name)))
              raise Exception

            org = org[1:-1]  # Remove quotes from original string

            if (org != ''):  # Only append non-empty values
              org_list.append(org)
              repl_list.append(repl)
              len_list.append(len(org))

    tmp_list = map(None,len_list,org_list,repl_list)
    tmp_list.sort()
    tmp_list.reverse()

    for (i,org,repl) in tmp_list:
      self.append((org,repl))

    self.length = self.__len__()  # Get number of elements in the look-up table

    # A log message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    logging.info('Loaded correction list "%s"' % (self.description))
    logging.info('  From file:         %s' % (str(self.file_name)))
    logging.info('  Number of entries: %i' % (self.length))

# =============================================================================
