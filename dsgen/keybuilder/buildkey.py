# Script to product keys by concatenating the first field values in a CSV
# file using '-'. This can be used to generate look-up tables that include
# field (or attribute) dependencies.
#
# Agus Pudjijono, 2008.

# =============================================================================

import os
import sys

# =============================================================================
# Start main program

if (len(sys.argv) != 4):  
  print 'Three arguments needed with %s:' % (sys.argv[0])
  print '  - Input File'
  print '  - Output File'
  print '  - Number of keys'
  sys.exit()

file_input  = sys.argv[1]
file_output = sys.argv[2]
num_of_key  = int(sys.argv[3]) # The number of fields to concatenate into a key

if (num_of_key <= 0):
  print 'Error: Number of keys must be positive'
  sys.exit()

try:
   f=open(file_input)
except:
   print 'Error cannot read %s' % (file_input)
   raise IOError

file_data = f.readlines()
f.close()  

out_line = ''

for line in file_data:
  l = line.strip() 
  if (l != ''):

    val_list = l.split(',')
    the_key = ''
    the_val = ''

    for i in range(num_of_key):

      the_key += str(val_list[i]) + '-'

    the_key = the_key[:-1]  # Remove last hyphen

    #print 'key '+the_key 

    for i in range(len(val_list)):
      if (i > num_of_key-1):
        the_val += val_list[i]+','

    the_val = the_val[:-1]  # Remove last comma

    #print 'val '+the_val 

    out_line += the_key + ':' + the_val + os.linesep
 
#Create output file
#
output_file = file_output
try:
  f_out = open(output_file, 'w')
except:
  print 'Error: Can not write to output file "%s"' % (output_file)
  raise IOError
   
f_out.write (out_line+os.linesep)

f_out.close()

# =============================================================================