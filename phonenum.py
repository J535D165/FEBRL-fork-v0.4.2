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
# The Original Software is: "phonenum.py"
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

"""Module phonenum.py - Routines for parsing and standardising US, Australian
                        and international phone numbers.

   Original developer: PV 11/07/2003.

   Extended and improved for Australian phone numbers by Peter Christen,
   01/2005 - 09/2007.

   This module provides one method 'str_to_phonenum' to parse a phone number
   (given as a string) into a list with elements:

     [country_code, country_name, area_code, number, extension]

   If the input phone number cannot be parsed successfully, parts of the
   returned list might be empty strings, with the number being returned in
   the 'number' field.

   If the input number is empty an empty list [] is returned.

   One input argument 'default_country' for country specific parsing needs to
   be given to the 'str_to_phonenum' method. Possible values currently are:
     'australia' (default)
     'canada/usa'
"""

# =============================================================================
# Imports go here

import logging
import string
import re

# =============================================================================

# Define a dictionary of international phone codes
#
phone_code_dict = {'Afghanistan':'93', 'Albania':'355', 'Algeria':'213',
                   'American Samoa':'684', 'Andorra':'376', 'Angola':'244',
                   'Anguilla':'1264', 'Antarctica':'672', 'Antigua':'1268',
                   'Argentina':'54', 'Armenia':'374', 'Aruba':'297',
                   'Ascension Island':'247', 'Australia':'61', 'Austria':'43',
                   'Azerbaijan':'994', 'Bahamas':'1242', 'Bahrain':'973',
                   'Bangladesh':'880', 'Barbados':'1246', 'Barbuda':'1268',
                   'Belarus':'375', 'Belgium':'32', 'Belize':'501',
                   'Benin':'229', 'Bermuda':'1441', 'Bhutan':'975',
                   'Bolivia':'591', 'Bosnia and Herzegovina':'387',
                   'Botswana':'267', 'Brazil':'55',
                   'British Virgin Islands':'1284', 'Brunei':'673',
                   'Bulgaria':'359', 'Burkina Faso':'226', 'Burundi':'257',
                   'Cambodia':'855', 'Cameroon':'237', 'Canada/USA':'1',
                   'Cape Verde Islands':'238', 'Cayman Islands':'1345',
                   'Central African Republic':'236', 'Chad':'235',
                   'Chatham Island (New Zealand)':'64', 'Chile':'56',
                   'China (PRC)':'86', 'Christmas Island':'618',
                   'CocosKeeling Islands':'61', 'Colombia':'57',
                   'Comoros':'269', 'Congo':'242',
                   'Dem. Rep. of Congo (former Zaire) ':'243',
                   'Cook Islands':'682', 'Costa Rica':'506', 'Croatia':'385',
                   'Cuba':'53', 'Cuba (Guantanamo Bay)':'5399',
                   'Curacao':'599', 'Cyprus':'357', 'Czech Republic':'420',
                   'Denmark':'45', 'Diego Garcia':'246', 'Djibouti':'253',
                   'Dominica':'1767', 'Dominican Republic':'1809',
                   'East Timor':'670', 'Easter Island':'56', 'Ecuador':'593',
                   'Egypt':'20', 'El Salvador':'503',
                   'Equatorial Guinea':'240', 'Eritrea':'291', 'Estonia':'372',
                   'Ethiopia':'251', 'Faeroe Islands':'298',
                   'Falkland Islands':'500', 'Fiji Islands':'679',
                   'Finland':'358', 'France':'33', 'French Antilles':'596',
                   'French Guiana':'594', 'French Polynesia':'689',
                   'Gabon':'241', 'Gambia':'220', 'Georgia':'995',
                   'Germany':'49', 'Ghana':'233', 'Gibraltar':'350',
                   'Global Mobile Satellite System (GMSS)':'881',
                   'Greece':'30', 'Greenland':'299', 'Grenada':'1473',
                   'Guadeloupe':'590', 'Guam':'1671', 'Guantanamo Bay':'5399',
                   'Guatemala':'502', 'Guinea Bissau':'245',
                   'Guinea (PRP)':'224', 'Guyana':'592', 'Haiti':'509',
                   'Honduras':'504', 'Hong Kong':'852', 'Hungary':'36',
                   'Iceland':'354', 'India':'91', 'Indonesia':'62',
                   'Inmarsat (Atlantic Ocean  East)':'871',
                   'Inmarsat (Atlantic Ocean  West)':'874',
                   'Inmarsat (Indian Ocean)':'873',
                   'Inmarsat (Pacific Ocean)':'872', 'Inmarsat SNAC':'870',
                   'Iran':'98', 'Iraq':'964', 'Ireland':'353',
                   'Northern Ireland':'48',
                   'Iridium (Mobile Satellite service)':'8816|8817',
                   'Israel':'972', 'Italy':'39', 'Ivory Coast':'225',
                   'Jamaica':'1876', 'Japan':'81', 'Jordan':'962',
                   'Kazakhstan':'7', 'Kenya':'254', 'Kiribati':'686',
                   'Korea (North)':'850', 'Korea (South)':'82', 'Kuwait':'965',
                   'Kyrgyz Republic':'996', 'Laos':'856', 'Latvia':'371',
                   'Lebanon':'961', 'Lesotho':'266', 'Liberia':'231',
                   'Libya':'218', 'Liechtenstein':'423', 'Lithuania':'370',
                   'Luxembourg':'352', 'Macau':'853',
                   'Macedonia (former Yugoslav Rep.)':'389',
                   'Madagascar':'261', 'Malawi':'265', 'Malaysia':'60',
                   'Maldives':'960', 'Mali Republic':'223', 'Malta':'356',
                   'Marshall Islands':'692', 'Martinique':'596',
                   'Mauritania':'222', 'Mauritius':'230',
                   'Mayotte Island':'269', 'Mexico':'52',
                   'Federal States of Micronesia':'691',
                   'Midway Island':'1808', 'Moldova':'373', 'Monaco':'377',
                   'Mongolia':'976', 'Montserrat':'1664', 'Morocco':'212',
                   'Mozambique':'258', 'Myanmar':'95', 'Namibia':'264',
                   'Nauru':'674', 'Nepal':'977', 'Netherlands':'31',
                   'Netherlands Antilles':'599', 'Nevis':'1869',
                   'New Caledonia':'687', 'New Zealand':'64',
                   'Nicaragua':'505', 'Niger':'227', 'Nigeria':'234',
                   'Niue':'683', 'Norfolk Island':'672',
                   'Northern Marianas Islands':'1670', 'Norway':'47',
                   'Oman':'968', 'Pakistan':'92', 'Palau':'680',
                   'Palestine':'970', 'Panama':'507',
                   'Papua New Guinea':'675', 'Paraguay':'595', 'Peru':'51',
                   'Philippines':'63', 'Poland':'48', 'Portugal':'351',
                   'Puerto Rico':'1787|1939', 'Qatar':'974',
                   'Reunion Island':'262', 'Romania':'40', 'Russia':'7',
                   'Rwanda':'250', 'St. Helena':'290',
                   'St. Kitts - Nevis':'1869', 'St. Lucia':'1758',
                   'St. Pierre and Miquelon':'508',
                   'St. Vincent and Grenadines':'1784', 'San Marino':'378',
                   'Sao Tome and Principe':'239', 'Saudi Arabia':'966',
                   'Senegal':'221', 'Serbia':'381', 'Seychelles Islands':'248',
                   'Sierra Leone':'232', 'Singapore':'65',
                   'Slovak Republic':'421', 'Slovenia':'386',
                   'Solomon Islands':'677', 'Somalia':'252',
                   'South Africa':'27', 'Spain':'34', 'Sri Lanka':'94',
                   'Sudan':'249', 'Suriname':'597', 'Swaziland':'268',
                   'Sweden':'46', 'Switzerland':'41', 'Syria':'963',
                   'Taiwan':'886', 'Tajikistan':'992', 'Tanzania':'255',
                   'Thailand':'66', 'Togo':'228', 'Tokelau':'690',
                   'Tonga Islands':'676', 'Trinidad and Tobago':'1868',
                   'Tunisia':'216', 'Turkey':'90', 'Turkmenistan':'993',
                   'Turks and Caicos Islands':'1649', 'Tuvalu':'688',
                   'Uganda':'256', 'Ukraine':'380',
                   'United Arab Emirates':'971', 'United Kingdom':'44',
                   'USA/Canada':'1', 'US Virgin Islands':'1340',
                   'Universal Personal Telecommunications (UPT)':'878',
                   'Uruguay':'598', 'Uzbekistan':'998', 'Vanuatu':'678',
                   'Vatican City':'39', 'Venezuela':'58', 'Vietnam':'84',
                   'Wake Island':'808', 'Wallis and Futuna Islands':'681',
                   'Western Samoa':'685', 'Yemen':'967', 'Yugoslavia':'381',
                   'Zambia':'260', 'Zanzibar':'255', 'Zimbabwe':'263'}

# Define a dictionary of Australian area codes. For more details please see:
# http://www.aca.gov.au/, then go to: "Telecommunications" >
# "Telephone Numbering" > "8 Digit Numbering"
#
australia_area_codes = {'02':'Central east region',
                        '03':'South east region',
                        '04':'Mobile phones',
                        '07':'North east region',
                        '08':'Central and west region'}

# Define a character replace table for data strings - - - - - - - - - - - -
#
string_replace = ["'./,+:-_\\()[]<>{}", \
                  "                 "]

# Characters in the first list are replaced by the corresponding character in
# the second list

replace_table = string.maketrans(string_replace[0], string_replace[1])

# =============================================================================

def parse_australia_phone_number(phonenum_str):
  """A routine to check if the given phone number could be an Australian number
     and if so parse it into area code and number. The routine returns either
     [area_code, number] or if not successful an empty list [].

     It is assumed that the input phone number string only contains digits and
     no spaces or letters.
  """

  if (not phonenum_str.isdigit()):
    logging.exception('Input phone number contains more than just digits: %s' \
                      % (phonenum_str))
    raise Exception

  logging.debug('Attempting to parse Australian phone number: %s' % \
                (phonenum_str))

  valid = False

  # First check if the leading '0' has maybe been removed
  #
  if (len(phonenum_str) == 9) and (phonenum_str[0] != '0'):
    phonenum_str = '0'+phonenum_str

  if (len(phonenum_str) == 8):  # Standard Australian number without area code
    area_code = ''
    number =    phonenum_str[:4]+'-'+phonenum_str[4:]
    valid = True
  elif (len(phonenum_str) == 10):  # Standard Australian number with area code
    area_code = phonenum_str[:2]

    if (area_code in australia_area_codes):
      logging.debug('Valid Australian area code (%s) for: %s' % \
            (area_code, australia_area_codes[area_code]))
      number = phonenum_str[2:6]+'-'+phonenum_str[6:]
      valid =  True
    else:
      logging.warn('Illegal Australian area code: %s' % (area_code))

  else:  # Phone number not 8 or 10 digits is not Australian standard
    logging.warn('Illegal Australian phone number (not 8 or 10 digits): %s' \
          % (phonenum_str))

  if (valid == True):
    return [area_code, number]
  else:
    return []

# =============================================================================

def parse_usa_canada_phone_number(phonenum_str):
  """A routine to check if the given phone number could be a Canadian or US
     number and if so parse it into area code and number. The routine returns
     either [area_code, number] or if not successful an empty list [].

     It is assumed that the input phone number string only contains digits and
     no spaces or letters.

     Original developer: PV, 11/07/2003
  """

  if (not phonenum_str.isdigit()):
    logging.exception('Input phone number contains more than just digits: %s' \
                      % (phonenum_str))
    raise Exception

  logging.debug('Attempting to parse Canadian/US phone number: %s' % \
                (phonenum_str))

  valid = False

  if (len(phonenum_str) == 10):  # Parse as 3+3+4
    area_code = phonenum_str[:3]
    number =    phonenum_str[3:6]+'-'+phonenum_str[6:]
    valid =     True

  elif (len(phonenum_str) == 7):  # Parse as 'UNK'+3+4
    area_code = ''
    number =    phonenum_str[:3]+'-'+phonenum_str[3:]
    valid =     True

  else:  # Phone number not 7 or 10 digits is not Canadian/US standard
    logging.warn('Illegal Canadian/US phone number (not 7 or 10 digits): %s' \
          % (phonenum_str))

  if (valid == True):
    return [area_code, number]
  else:
    return []

# =============================================================================

def str_to_phonenum(phonenum_str, default_country='australia'):
  """A routine that converts a string into an Australian or US centric
     standardized phone number.

  USAGE:
    [country_code, country_name, area_code, number, extension] = \
      str_to_phonenum(phonenum_str, default_country)

  ARGUMENTS:
    phonenum_str     Input phone number (raw) as string
    default_country  For country specific parsing (possible values are
                     'australia' (default) or 'canada/usa')

  DESCRIPTION:
    This routine parses the input raw phonenum string into a list:

      [country_code, country_name, area_code, number, extension]

    The output is a list of five strings, some possibly empty.

    If the input number is an empty string an empty list is returned.
  """

  logging.debug('Input phone number string:   %s' % (phonenum_str))
  logging.debug('Default country for parsing: %s' % (default_country))

  # Copy arg to buffer
  #
  buffer_str = phonenum_str

  # Do some pre-processing  - - - - - - - - - - - - - - - - - - - - - - - - - -

  # Apply character replace table
  #
  buffer_str = buffer_str.translate(replace_table)

  # Remove leading and trailing whitespaces
  #
  buffer_str = buffer_str.strip()

  # Replace triple and double spaces with one space only
  #
  while '  ' in buffer_str:
    buffer_str = buffer_str.replace('  ', ' ')

  # Now check the cleaned phonenum string to make sure it has some content- - -
  #
  if (buffer_str == ''):
    logging.info('Cleaned input string (%s) is empty' % (phonenum_str))
    return []

  logging.debug('Cleaned phone number string:   %s' % (buffer_str))

  international = False
  country_code =  ''
  country_name =  ''
  area_code =     ''
  extension =     ''
  valid =         False

  # Find and remove IDD prefixes (e.g. 0011, 011, 0001, etc.) - - - - - - - - -
  #
  idd_patterns_list = ["^(0+1{2,})", "^(0+)1","^(0+[01]*)"]
  for idd_pattern in idd_patterns_list:
      p = re.compile(idd_pattern)
      m = p.match(buffer_str)

      if m:

        # If a match is found stop looking
        #
        idd_code =   buffer_str[:m.end()]
        buffer_str = buffer_str[m.end():].strip()

        international = True

        logging.debug('Found matching IDD code (%s), new phone number string' \
                      ': %s' % (idd_code, buffer_str))

        break  # Exit the loop

  # Remove and store extension  - - - - - - - - - - - - - - - - - - - - - - - -
  #
  ext_pat = r"[ ]*(x|ext)\.?[ ]*"
  p = re.compile(ext_pat, re.IGNORECASE)
  m = p.search(buffer_str)

  if m:
    extension =  buffer_str[m.end():]    # Extract extension found
    buffer_str = buffer_str[:m.start()]  # Remove extension

    logging.debug('Found an extension (%s), remaining number: %s' % \
          (extension, buffer_str))
  else:
    extension = ''

  # Special case `07' which is a area code in Australia and `7' being Russia's
  # country code
  #
  if (default_country == 'australia'):
    if (buffer_str[0] == '7') and (len(buffer_str) == 10):  # Australian number
      buffer_str = '0'+buffer_str  # So no confusion with a Russian number

  # Match against list of country codes - - - - - - - - - - - - - - - - - - - -
  #
  for name, code in phone_code_dict.items():
      if len(code) == 4:
          match_pat = '1[ ]*' + code[1:]  # Allow spaces with 4-digit codes
      else:
          match_pat = code + ' '

      p = re.compile(match_pat)
      m = p.match(buffer_str)

      if m:

        # If a match is found stop looking and store the normed country code
        #
        buffer_str =   buffer_str[m.end():].strip()  # Remove country code
        country_code = code
        country_name = name
        number =       buffer_str  # Remaining number

        international = True
        valid =         True
        area_code =     ''

        logging.debug('Found country code (%s) of country "%s"' % \
              (country_code, country_name))
        logging.debug('Remaining phone number: %s' % (buffer_str))

        break  # Exit the loop

  # If no match has been found try to strip off leading 1 - - - - - - - - - - -
  #
  if (international == False):

    if buffer_str[0] == '1':
      buffer_str = buffer_str[1:].strip()

      logging.debug('Removed leading "1", remaining phone number: %s' % \
                    (buffer_str))

  # Remove all remaining spaces in the core number - - - - - - - - - - - - - -
  #
  buffer_str = buffer_str.replace(' ','')

  # Now do country specific parsing - - - - - - - - - - - - - - - - - - - - - -
  #
  if (international == False) or \
     (country_name in ['', 'Australia', 'Canada/USA', 'USA/Canada']):

    if (default_country == 'australia'):

      # Try parsing number as Australian first
      #
      parse_res = parse_australia_phone_number(buffer_str)

      if (parse_res != []):  # Successful parsing
        country_code = '61'
        country_name = 'Australia'
        area_code =    parse_res[0]
        number =       parse_res[1]
        valid =        True

      else:
        # Try parsing number as Canadian/US next
        #
        parse_res = parse_usa_canada_phone_number(buffer_str)

        if (parse_res != []):  # Successful parsing
          country_code = '1'
          country_name = 'USA/Canada'
          area_code =    parse_res[0]
          number =       parse_res[1]
          valid =        True

    elif (default_country == 'canada/usa'):

      # Try parsing number as Canadian/US first
      #
      parse_res = parse_usa_canada_phone_number(buffer_str)

      if (parse_res != []):  # Successful parsing
        country_code = '1'
        country_name = 'USA/Canada'
        area_code =    parse_res[0]
        number =       parse_res[1]
        valid =        True

      else:
        # Try parsing number as Australian next
        #
        parse_res = parse_australia_phone_number(buffer_str)

        if (parse_res != []):  # Successful parsing
          country_code = '61'
          country_name = 'Australia'
          area_code =    parse_res[0]
          number =       parse_res[1]
          valid =        True

  # If valid return the found phone number elements - - - - - - - - - - - - - -
  #
  if (valid == True):
    return [country_code, country_name, area_code, number, extension]
  else:
    return [country_code, country_name, area_code, buffer_str, extension]

# =============================================================================
#
# Do some tests if called from command line
#
if (__name__ == '__main__'):

  phone_numbers = ['++61 2 6125 5690',
                   '0061 02 6125 5690',
                   '0061   02    6125-5690',
                   '41 312 17 84',
                   '6125 0010',
                   '1-800-764-0432',
                   '02 6125 0010',
                   '00 1 317-923 4523',
                   '1 317-923 4523',
                   '00111 41 312 17 84',
                   '00001 41 312 17 84',
                   '01 41 312 17 84',
                   '1-541-754-3010',
                   '754-3010',
                   '(541) 754-3010',
                   '+1-541-754-3010',
                   '191 541 754 3010',
                   '001-541-754-3010',
                   '636-48018',
                   '(089) / 636-48018',
                   '+49-89-636-48018',
                   '19-49-89-636-48018',
                   '+61 (02) 6125 0101',
                   '++61 (02) 6125 0101',
                   '++61 (2) 6125 0101',
                   '11 +61 (2) 6125 0101',
                   '0011 ++61 (2) 6125 0101',
                   '0111 ++61 (2) 6125 0101',
                   '0111 61 02 6125 0101',
                   '61 (2) 6125 0101',
                   '07 61316116',
                   '07 6131 6116',
                   '(07) 6131 6116',
                   '(07) 61316116',
                   '09 6131 6116',
                   '09 61316116',
                   '']
  print
  for pn in phone_numbers:
    res = str_to_phonenum(pn)
    print pn, res
    print
