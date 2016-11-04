#!/usr/bin/env python

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
# The Original Software is: "guiFebrl.py"
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

#!/usr/bin/env python

# =============================================================================
# Define some constants
#
FEBRL_GLADE_FILE = "gui/febrl.glade"  # Name of Febrl GUI Glade file

NUM_DATA_ROWS = 15  # Number of rows from data set(s) to show in Data page

# =============================================================================
# Import necessary modules (Python standard modules first, then Febrl modules)

import csv
import gzip
import logging
import math
import os
import sys
import time

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
FEBRL_DIR = sys.path[0]  # Get the directory where Febrl is located
print 'Febrl directory: %s' % (FEBRL_DIR)

# May have to include non system wide site-packages
#
sys.path = ['/home/christen/lib/python2.5/site-packages/']+sys.path

import auxiliary
import classification
import comparison
import dataset
import encode
import indexing
import lookup
import measurements
import mymath
import output
import simplehmm
import stringcmp
import standardisation

# Adjust full file name where Glade file is - - - - - - - - - - - - - - - - - -
#
FEBRL_GLADE_FILE = FEBRL_DIR+os.sep+FEBRL_GLADE_FILE

#sys.path.append("/home/christen/lib/python2.5/site-packages")
#sys.path.append(FEBRL_DIR)  # Make sure Febrl modules can be imported

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Import graphics library for Evaluate page

try:
  import matplotlib
  matplotlib.use('GTK')

  import matplotlib.figure
  import matplotlib.backends.backend_gtk

  imp_matplot = True

except:
  imp_matplot = False
  logging.warn('Matplotlib module not installed.')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Import GTK related modules, make sure correct version is imported

try:
  import pygtk
  pygtk.require("2.0")
except:
  pass

try:
  import gtk
  import gtk.glade
  import pango
except:
  print "GTK and PyGTK not installed."
  sys.exit(1)

# Set the logging level to warnings - - - - - - - - - - - - - - - - - - - - - -
#
my_logger = logging.getLogger()  # New logger at root level
my_logger.setLevel(logging.WARNING)


# =============================================================================
# Main class for the main GUI window, handles all events from the GUI

class MainFebrlWindow:

  def __init__(self):

    # Initialise all Febrl project related variables - - - - - - - - - - - - -
    #
    self.initFebrlProjectVariables()

    #--------------------------------------------------------------------------
    # Initialise GUI related variables

    # Glade file name and top level window names
    #
    self.glade_file_name =         FEBRL_GLADE_FILE
    self.main_window_name =        'febrl_main_window'
    self.about_dialog_name =       'febrl_about_dialog'
    self.license_dialog_name =     'febrl_license_dialog'
    self.new_project_dialog_name = 'febrl_new_project_dialog'
    self.save_file_dialog_name =   'febrl_save_file_dialog'

#    self.current_folder = None  # Currently selected folder

    # Current page and page names in main notebook
    #
    self.main_notebook_curr_page =  0  # Set to Data page at beginning

    self.main_notebook_page_names = ['Data', 'Explore', 'Standardise',
                                     'Index', 'Compare', 'Classify', 'Run',
                                     'Evaluate', 'Review', 'Log']

    # A dictionary, for each notebook page flag if it is active or not
    #
    self.main_notebook_page_active_dict = {'Data':True, 'Explore':False,
                                           'Standardise':False,
                                           'Index':False,
                                           'Compare':False,
                                           'Classify':False, 'Run':False,
                                           'Evaluate':False, 'Review':False,
                                           'Log':True}

    # Mappings for data set types into Data page notebook page numbers
    #
    self.data_set_type_page_dict = {'CSV':0, 'COL':1, 'TAB':2, 'SQL':3}

    # Last value of the save file dialog (file name)
    #
    self.save_dialog_file_name = None

    # Load the Glade file and initialise widgets tree - - - - - - - - - - - - -
    #
    self.mainTree = gtk.glade.XML(self.glade_file_name, self.main_window_name)
    self.mainWin =  self.mainTree.get_widget(self.main_window_name)

    self.setWindowTitle()  # Set window title to project name

    self.setProjectTypeDataButtons()  # Set the project type buttons

    # Set up the log buffer - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    self.log_page_textview = self.mainTree.get_widget('log_text_view')

    # Set font to fixed width
    #
    self.log_page_textview.modify_font(pango.FontDescription("courier 10"))

    self.log_page_buffer =   self.log_page_textview.get_buffer()

    # Build dictionary of signals - - - - - - - - - - - - - - - - - - - - - - -
    #
    dict = {'main_window_quit':self.quitMain,
            'main_window_execute':self.executeMain,
            'main_window_new':self.newMain,
            'main_window_save':self.saveMain,
            'project_type_toggle':self.projectTypeButtonToggle,
            'main_notebook_switch_page':self.mainNotebookSwitchPage,
            'online_help_activate':self.activateOnlineHelp,
            'about_activate':self.activateAbout,
            'license_activate':self.activateLicense,

            # Data page signals
            'data_set_type_toggle':self.dataSetTypeButtonToggle,
            'data_set_header_line_toggle':self.dataSetTypeHeaderLineToggle,
            'data_set_strip_fields_toggle':self.dataSetTypeStripFieldsToggle,
            'file_name_changed':self.dataFileSelect,
            'miss_val_activate':self.dataSetMissValueActivate,
            'delimiter_activate':self.dataSetDelimiterActivate,

            # Explore page signals

            # Indexing page signals
            'index_method_changed':self.indexChangeMethod,

            # Classification page signals
            'classifier_method_changed':self.classifierChangeMethod,

            # Output/run page signals
            'clicked_file_name_button':self.runClickedFileNameButton,
            'save_file_button_toggle':self.runSaveFileButtonToggle,

            }

    self.mainTree.signal_autoconnect(dict)

    # Special handling of click on close window (x) - - - - - - - - - - - - - -
    #
    self.mainWin.connect('delete-event', self.closeMainWindow)

    self.displayCurrentNotebookPage()  # Diplay the current page

  # ===========================================================================

  # Method to initialise all Febrl project related variables
  # (default is a deduplication project)
  #
  def initFebrlProjectVariables(self, project_type='Deduplicate'):

    if (project_type not in ['Standardise', 'Deduplicate', 'Link', 'Geocode']):
      raise Exception, 'Illegal project type provided: %s' % (project_type)

    self.project_type = project_type
    self.project_name = None   # Will be set by New, Load or Save project

    # A dictionary with flags for modified components corresponding to code
    # (set to False initially as no code has been generated)
    #
    self.modified = {'data':False, 'standardise':False, 'indexing':False,
                     'comparison':False, 'classification':False,
                     'output':False} # Set to True when values change

    # A dictionary with flags used when running all code to determine if data
    # set have to be re-initialised or not, and if weight vectors have to be
    # re-calculated or not.
    #
    self.re_run = {'data_init':False, 'w_vec_generate':False}

    # A dictionary with the generated python codes - - - - - - - - - - - - - -
    #
    self.febrl_code = {'data':None, 'standardise':None, 'indexing':None,
                       'comparison':None, 'classification':None, 'output':None}

    # A dictionary with the available string comparison methods (names and
    # corresponding stringcmp methods)
    #
    self.stringcmp_dict = {'Jaro':'jaro', 'Winkler':'winkler',
                           'Q-Gram':'qgram', 'Pos-Q-Gram':'posqgram',
                           'S-Gram':'sgram', 'Edit-Dist':'editdist',
                           'Mod-Edit-Dist':'mod_editdist',
                           'Bag-Dist':'bagdist','Smith-Water-Dist':'swdist',
                           'Syll-Align-Dist':'syllaligndist',
                           'Seq-Match':'seqmatch', 'Long-Common-Seq':'lcs',
                           'Onto-LCS':'ontolcs', 'Editex':'editex',
                           'Perm-Winkler':'permwinkler',
                           'Sort-Winkler':'sortwinkler'}

    # A dictionary with the available string encoding methods (names and
    # corresponding encoding methods)
    #
    self.stringencode_dict = {'None':'', 'Soundex':'soundex',
                              'Mod-Soundex':'mod_soundex', 'Phonex':'phonex',
                              'NYSIIS':'nysiis',
                              'Double-Metaphone':'dmetaphone',
                              'Phonix':'phonix',
                              'Fuzzy-Soundex':'fuzzysoundex',
                              'Substring':'get_substring'}

    # Data (page) related variables -------------------------------------------
    #
    if (self.project_type == 'Link'):  # Initialise data set(s) as of type CSV
      self.data_set_type_list = ['CSV', 'CSV']
    else:
      self.data_set_type_list = ['CSV', None]

    # A list, one entry per data set (both entries only used for a linkage),
    # that will contain default values for data set related variables once they
    # have been initialised (one dictionary per data set)
    #
    self.data_set_default_values = {'header_line':True, 'strip_fields':True,
                                    'header_data':None, 'file_name':None,
                                    'file_data':None,   'miss_val':[''],
                                    'delimiter':',',    'field_names':None,
                                    'rec_id_field':'__rec_id'}
    default_vals = self.data_set_default_values

    # The dictionaries (one per data set) which hold the actual values
    #
    self.data_set_info_list = [{'header_line':default_vals['header_line'],
                                'strip_fields':default_vals['strip_fields'],
                                'header_data':default_vals['header_data'],
                                'file_name':default_vals['file_name'],
                                'file_data':default_vals['file_data'],
                                'miss_val':default_vals['miss_val'],
                                'delimiter':default_vals['delimiter'],
                                'field_names':default_vals['field_names'],
                                'rec_id_field':default_vals['rec_id_field']},
                                {'header_line':default_vals['header_line'],
                                'strip_fields':default_vals['strip_fields'],
                                'header_data':default_vals['header_data'],
                                'file_name':default_vals['file_name'],
                                'file_data':default_vals['file_data'],
                                'miss_val':default_vals['miss_val'],
                                'delimiter':default_vals['delimiter'],
                                'field_names':default_vals['field_names'],
                                'rec_id_field':default_vals['rec_id_field']}]

    self.list_data_store = [None, None]

    # The number of rows (records) to show from data set(s)
    #
    self.num_data_rows = NUM_DATA_ROWS*2  # Double for deduplication

    # Dictionary with file filter data (data set type, name for widget, filter)

    self.file_filters = {'CSV':['CSV files',('*.csv', '*.CSV', '*.csv.gz',
                                             '*.CSV.GZ')],
                         'COL':['COL files',('*.col', '*.COL', '*.col.gz',
                                             '*.COL.GZ')],
                         'TAB':['TAB files',('*.tab', '*.TAB', '*.tab.gz',
                                             '*.TAB.GZ')],
                         'TXT':['Text files',('*.csv', '*.CSV', '*.csv.gz',
                                             '*.CSV.GZ', '*.col', '*.COL',
                                             '*.col.gz', '*.COL.GZ', '*.tab',
                                             '*.TAB', '*.tab.gz', '*.TAB.GZ',
                                             '*.txt', '*.TXT', '*.txt.gz',
                                             '*.TXT.GZ')],
                         'ALL':['All files', ('*',)]}

    # Standardisation (page) related variables --------------------------------

    # A dictionary with component standardiser types names as keys and their
    # corresponding output fields as lists
    #
    self.comp_std_out_fields = {'Date':['day','month','year'],
                                'Phon':['country_code', 'country_name',
                                        'area_code', 'number', 'extension'],
                                'Name':['title', 'gender_guess',
                                        'given_name', 'alt_given_name',
                                        'surname', 'alt_surname'],
                                'Addr':['building_name', 'post_address_type',
                                        'post_address_number',
                                        'lot_number_prefix', 'lot_number',
                                        'lot_number_suffix',
                                        'flat_number_prefix',
                                        'flat_number', 'flat_number_suffix',
                                        'flat_type', 'level_number_prefix',
                                        'level_number', 'level_number_suffix',
                                        'level_type', 'number_first_prefix',
                                        'number_first', 'number_first_suffix',
                                        'number_last_prefix', 'number_last',
                                        'number_last_suffix', 'street_name',
                                        'street_suffix', 'street_type',
                                        'locality_name', 'postcode',
                                        'state_abbrev', 'country']}
    self.comp_std_list = []

    # Date parsing formats for data standardiser
    #
    self.parse_format_str = '%d %m %Y, %d %B %Y, %d %b %Y, %m %d %Y, ' + \
                            '%B %d %Y, %b %d %Y, %Y %m %d, %Y %B %d, ' + \
                            '%Y %b %d, %d %m %y, %d %B %y, %d %b %y, ' + \
                            '%y %m %d, %y %B %d, %y %b %d, %m %d %y, ' + \
                            '%B %d %y, %b %d %y'

    # Indexing (page) related variables ---------------------------------------

    # A list with the names of the possible index method types
    #
    self.index_names = ['FullIndex', 'BlockingIndex', 'SortingIndex',
                        'QGramIndex', 'CanopyIndex', 'StringMapIndex',
                        'SuffixArrayIndex']

    # The following dictionary will include the parameters of the selected
    # indxing method as keys and their values as values. The name will be one
    # of the above names.
    # (default values for general parameters)
    #
    self.index_method = {'name':None, 'index_sep_str':'', 'skip_missing':True}

    # The following list will contain the index definitions, with each being a
    # dictionary with the corresponding index definition details
    #
    self.index_def = []
    self.index_num = 0  # Number of indices defined

    # Comparison (page) related variables -------------------------------------

    # A list with the names of the possible field comparison method types
    #
    self.field_comp_dict = {'Age':'FieldComparatorAge',
                            'Bag-Dist':'FieldComparatorBagDist',
                            'Compression':'FieldComparatorCompress',
                            'Dam-Le-Edit-Dist':'FieldComparatorDaLeDist',
                            'Edit-Dist':'FieldComparatorEditDist',
                            'Editex':'FieldComparatorEditex',
                            'Jaro':'FieldComparatorJaro',
                            'Long-Common-Seq':'FieldComparatorLCS',
                            'Onto-LCS':'FieldComparatorOntoLCS',
                            'Pos-Q-Gram':'FieldComparatorPosQGram',
                            'Q-Gram':'FieldComparatorQGram',
                            'S-Gram':'FieldComparatorSGram',
                            'Smith-Water-Dist':'FieldComparatorSWDist',
                            'Seq-Match':'FieldComparatorSeqMatch',
                            'Syll-Align-Dist':'FieldComparatorSyllAlDist',
                            'Winkler':'FieldComparatorWinkler',
                            'Date':'FieldComparatorDate',
### TODO                            'Distance':'FieldComparatorDistance',
                            'Str-Encode':'FieldComparatorEncodeString',
                            'Str-Exact':'FieldComparatorExactString',
                            'Str-Contains':'FieldComparatorContainsString',
                            'Key-Diff':'FieldComparatorKeyDiff',
                            'Num-Abs':'FieldComparatorNumericAbs',
                            'Num-Perc':'FieldComparatorNumericPerc',
                            'Time':'FieldComparatorTime',
                            'Str-Truncate':'FieldComparatorTruncateString',
                            'Token-Set':'FieldComparatorTokenSet'}

    # The following list will contain the field comparisons, with each being a
    # dictionary with the corresponding details
    #
    self.field_comp_list = []

    # Classification (page) related variables ---------------------------------

    # A list with the names of the possible classification method types
    #
    self.classifier_names = ['FellegiSunter', 'OptimalThreshold', 'KMeans',
                             'FarthestFirst', 'SuppVecMachine', 'TwoStep']

    # The following dictionary will include the parameters of the selected
    # classification method as keys and their values as values. The name will
    # be one of the above names.
    # (default values for general parameters)
    #
    self.classifier_method = {'name':None}

    # Output/Run page related variables ---------------------------------------
    #
    self.output_dict = {'w_vec_file':(False, '(None)'),
                        'histo_file':(False, '(None)', '1.0'),
                        'm_status_file':(False, '(None)'),
                        'm_datasets':(False, ('(None)', 'match_id'),
                                             ('(None)', 'match_id')),
                        'progress_perc':'10',
                        'length_filter_perc':None,
                        'cut_off_threshold':None,
        # Following variables are for a standardisation
                        'std_out_file':'(None)',
                        'pass_field_list':[]}

    self.w_vec_dict =       None
    self.class_w_vec_dict = None

    # Evaluation (page) related variables -------------------------------------
    #
    self.result_sets = None  # Will become a list with three sets (matches,
                             # non-matches, and possible matches)

    # Results will be added into here
    self.evaluation_dict = {'available':False}  # No results available yet

    # Log (page) related variables --------------------------------------------
    #
    self.log_page_text = '# Febrl log.\n'


  # ===========================================================================
  # Methods that handle GUI events
  # ===========================================================================

  # ---------------------------------------------------------------------------
  # Set the project type data buttons
  #
  def setProjectTypeDataButtons(self):

    # Get widgets of the four radion buttons
    #
    standard_project_widget = self.mainTree.get_widget('standard_radio_button')
    dedup_project_widget =    self.mainTree.get_widget('dedup_radio_button')
    link_project_widget =     self.mainTree.get_widget('link_radio_button')
    geocode_project_widget =  self.mainTree.get_widget('geocode_radio_button')

    standard_project_widget.set_active(False)  # De-activate all four first
    dedup_project_widget.set_active(False)
    link_project_widget.set_active(False)
    geocode_project_widget.set_active(False)

    if (self.project_type == 'Standardise'):  # Now activate the correct one
      standard_project_widget.set_active(True)
    elif (self.project_type == 'Deduplicate'):
      dedup_project_widget.set_active(True)
    elif (self.project_type == 'Link'):
      link_project_widget.set_active(True)
    elif (self.project_type == 'Geocode'):
      geocode_project_widget.set_active(True)
    else:
      raise Exception, 'Illegal project type: %s' % (self.project_type)

  # ---------------------------------------------------------------------------
  # Handle toggle of the project type button
  #
  def projectTypeButtonToggle(self, widget):

    if (widget.get_active() != True):  # Do nothing if not the active widget
      return

    widget_name = widget.get_name()
    if (widget_name == 'standard_radio_button'):
      self.project_type = 'Standardise'
    elif (widget_name == 'dedup_radio_button'):
      self.project_type = 'Deduplicate'
    elif (widget_name == 'link_radio_button'):
      self.project_type = 'Link'
    elif (widget_name == 'geocode_radio_button'):
      self.project_type = 'Geocode'
    else:
      raise Exception, 'Illegal project type widget name: %s' % (widget_name)

    print
    print ' *** Changed project type to:', self.project_type,'***'

    # For different project types set settings and visible pages - - - - - - -
    #
    if (self.project_type in ['Standardise','Deduplicate']):
      self.num_data_rows = NUM_DATA_ROWS*2  # Only one data set

    else:  # Linkage
      self.num_data_rows = NUM_DATA_ROWS  # Two data sets

      if (self.data_set_type_list[1] == None):  # If not initialised set to
        self.data_set_type_list[1] = 'CSV'      # default type

    # Make all pages invisible except Data and Log - - - - - - - - - - - - - -
    #
    for page_name in ['Explore', 'Standardise','Index', 'Compare', 'Classify',
                      'Run', 'Evaluate', 'Review']:
      self.main_notebook_page_active_dict[page_name] = False

    self.modified['data'] = True  # Project type has changed

    self.index_def = [] # Delete all previous index definitions
    self.index_num = 0

    self.field_comp_list = []  # Delete previous field comparison functions

    self.setWindowTitle()

    # If current page is already Data only display it, otherwise trigger switch
    #
    if (self.main_notebook_curr_page == 0):
      self.displayCurrentNotebookPage()

    else:  # Set new current page to Data
      notebook_widget = self.mainTree.get_widget('main_notebook')
      notebook_widget.set_current_page(0)

      gtk.main_iteration(False)  # Do a main iteration to update notebook
      gtk.main_iteration(False)  # Do a main iteration to update notebook

      # This will trigger a call to the method mainNotebookSwitchPage()

  # ---------------------------------------------------------------------------
  # Handle switch of main notebook page
  #
  def mainNotebookSwitchPage(self, widget, dummy, curr_page):
    # PyGTK FAQ 17.1: Needs these parameters!

    print
    print '*** Main notebook page switch ***'

    notebook_widget = self.mainTree.get_widget('main_notebook')

    prev_page =  self.main_notebook_curr_page  # Get the previous page number
    prev_page2 = notebook_widget.get_current_page()  # Current before switch

    # Special case where pages are 'hidden' results in a switch-page' signal
    # when an active page is hidden. This happens when the project type is
    # changed. As a consequence, the prev_page and prev_page2 values will be
    # different.
    #
    if (prev_page == prev_page2):
      prev_page_name = self.main_notebook_page_names[prev_page]
      curr_page_name = self.main_notebook_page_names[curr_page]

      curr_page_active = self.main_notebook_page_active_dict[curr_page_name]

      if (curr_page_active == True):
        self.main_notebook_curr_page = curr_page

        print 'Switched from page %s to page %s' % \
              (prev_page_name, curr_page_name)
        print '  Current page number:', curr_page

      else:
        print ' ************ else:', curr_page_name

      self.displayCurrentNotebookPage()  # Re-display current notebook page

    else:  # Must be a page is being hidden, do nothing
      print '  Stored current page: ', self.main_notebook_curr_page
      print '  Widget.get_curr_page:', prev_page2
      print '  Received curr page:  ', curr_page

  # ---------------------------------------------------------------------------
  # Display the current main notebook page
  #
  def displayCurrentNotebookPage(self):

    # Hide all non-active pages first
    #
    for (page_name, active) in self.main_notebook_page_active_dict.items():
      page_box_widget = self.mainTree.get_widget(page_name.lower()+'_page_box')

      if (active == True):
        page_box_widget.show()
      else:
        page_box_widget.hide()

    curr_page =      self.main_notebook_curr_page
    curr_page_name = self.main_notebook_page_names[curr_page]

    print 'Display current page:', curr_page, curr_page_name

    # If switched to Log or Evaluate page disable the Execute button - - - - -
    #
    execute_button_widget = self.mainTree.get_widget('execute_button')

    if (curr_page_name in ['Evaluate','Log']):
      execute_button_widget.set_sensitive(False)
    else:
      execute_button_widget.set_sensitive(True)

    # Handle switches to different pages - - - - - - - - - - - - - - - - - - -
    #
    if (curr_page_name == 'Data'):
      self.dataView()
    elif(curr_page_name == 'Explore'):
      self.exploreView()
    elif(curr_page_name == 'Standardise'):
      self.standardView()
    elif(curr_page_name == 'Index'):
      self.indexView()
    elif(curr_page_name == 'Compare'):
      self.compareView()
    elif(curr_page_name == 'Classify'):
      self.classifyView()
    elif(curr_page_name == 'Run'):
      self.runView()
    elif(curr_page_name == 'Evaluate'):
      self.evaluateView()
    elif(curr_page_name == 'Review'):
      self.reviewView()
    elif(curr_page_name == 'Log'):
      self.logView()
    else:
      raise Exception, 'Illegal notebook page name: %s' % (curr_page_name)

  # ---------------------------------------------------------------------------
  # Handle clicks on Execute button
  #
  def executeMain(self, widget):
    notebook_widget = self.mainTree.get_widget('main_notebook')

    curr_page_name =self.main_notebook_page_names[self.main_notebook_curr_page]

    # Handle executes on different pages - - - - - - - - - - - - - - - - -
    #
    if (curr_page_name == 'Data'):
      self.dataExecute()
    elif(curr_page_name == 'Explore'):
      self.exploreExecute()
    elif(curr_page_name == 'Standardise'):
      self.standardExecute()
    elif(curr_page_name == 'Index'):
      self.indexExecute()
    elif(curr_page_name == 'Compare'):
      self.compareExecute()
    elif(curr_page_name == 'Classify'):
      self.classifyExecute()
    elif(curr_page_name == 'Run'):
      self.runExecute()
    elif(curr_page_name == 'Evaluate'):
      self.evaluateExecute()
    elif(curr_page_name == 'Review'):
      self.reviewExecute()
    else:
      raise Exception

  # ---------------------------------------------------------------------------
  # Handle clicks on Quit button or activate of Quit menu
  #
  def quitMain(self, widget):
    print 'Clicked Quit or activated Quit menu'

    # Check if project is modified or not, if not, ask for save
    #
    if (sum(self.modified.values()) > 0):

      do_save = self.messageDialog('Current project is not saved!\n' + \
                                   'Do you want to save it?', 'c_question')
      if (do_save == True):
        self.saveMain(widget)

        if (sum(self.modified.values()) > 0):  # Saving has been canceled
          do_save = None  # Don't quit

    else:  # Not modified, so no saving needed
      do_save = False  # So quit will be performed

    if (do_save != None):  # If message dialog response was cancel don't quit
      gtk.main_quit()

  # Handle click on close window (x) - - - - - - - - - - - - - - - - - - - - -
  #
  def closeMainWindow(self, widget, event):
    print 'Clicked close window'

    self.quitMain(widget)

    return True

  # ---------------------------------------------------------------------------
  # Handle clicks on New button or activate of New menu
  #
  def newMain(self, widget):
    print 'Clicked New or activated New menu'

    # Check if project is modified or not, if not, ask for save
    #
    if (sum(self.modified.values()) > 0):
      do_save = self.messageDialog('Current project is not saved!\n' + \
                                   'Do you want to save it?', 'question')
      if (do_save == True):
        self.saveMain(widget)

    newProjTree = gtk.glade.XML(self.glade_file_name,
                                self.new_project_dialog_name)
    newProjDialog = newProjTree.get_widget(self.new_project_dialog_name)

    # Add the standard buttons and their responses
    #
    newProjDialog.add_button(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL)
    newProjDialog.add_button(gtk.STOCK_OK, gtk.RESPONSE_OK)

    newProjDialog.show()
    response =  newProjDialog.run()

    # Only get project values if OK was clicked
    #
    if (response == int(gtk.RESPONSE_OK)):

      standard_widget = newProjTree.get_widget('new_standard_button')
      dedup_widget =    newProjTree.get_widget('new_dedup_button')
      link_widget =     newProjTree.get_widget('new_link_button')
      geocode_widget =  newProjTree.get_widget('new_geocode_button')

      # Find the active widget
      #
      if (standard_widget.get_active() == True):
        self.project_type = 'Standardise'
      elif (dedup_widget.get_active() == True):
        self.project_type = 'Deduplicate'
      elif (link_widget.get_active() == True):
        self.project_type = 'Link'
      elif (geocode_widget.get_active() == True):
        self.project_type = 'Geocode'
      else:
        raise Exception, 'No project type activated'

      print
      print ' *** New project type:', self.project_type,'***'

      self.initFebrlProjectVariables(self.project_type)

      # Set the corresponding project type button (first de-activate all)
      #
      standard_proj_widget = self.mainTree.get_widget('standard_radio_button')
      dedup_proj_widget =    self.mainTree.get_widget('dedup_radio_button')
      link_proj_widget =     self.mainTree.get_widget('link_radio_button')
      geocode_proj_widget =  self.mainTree.get_widget('geocode_radio_button')
      standard_proj_widget.set_active(False)  # De-activate all four first
      dedup_proj_widget.set_active(False)
      link_proj_widget.set_active(False)
      geocode_proj_widget.set_active(False)

      if (self.project_type == 'Standardise'):  # Now activate the correct one
        standard_proj_widget.set_active(True)
      elif (self.project_type == 'Deduplicate'):
        dedup_proj_widget.set_active(True)
      elif (self.project_type == 'Link'):
        link_proj_widget.set_active(True)
      elif (self.project_type == 'Geocode'):
        geocode_proj_widget.set_active(True)

    newProjDialog.destroy()

#    self.displayCurrentNotebookPage()  # Re-display current notebook page

  # ---------------------------------------------------------------------------
  # Handle clicks on Save button or activate of Save or Save As menu
  #
  def saveMain(self, widget):
    print 'Clicked Save or activated Save or Save As menu, or Execute on' + \
          'Output/Run page'

    widget_name = widget.get_name()

    # Check if a project file name is set, if not activate save file first, or
    # if Save as was activated, always ask for a file name
    #
    if (self.project_name == None) or (widget_name == 'save_as'):

      saveFileTree = gtk.glade.XML(self.glade_file_name,
                                   self.save_file_dialog_name)
      saveFileDialog = saveFileTree.get_widget(self.save_file_dialog_name)

#      if (self.current_folder != None):
#        saveFileDialog.set_current_folder(self.current_folder)

      # Create a file filter (Python files only)
      #
      file_filter = gtk.FileFilter()  # Add a file filter for python files only
      file_filter.set_name('Python files')
      file_filter.add_pattern('*.py')
      saveFileDialog.add_filter(file_filter)

      saveFileDialog.show()
      save_file_response = saveFileDialog.run()
      self.save_dialog_file_name = saveFileDialog.get_filename()

#      self.current_folder = saveFileDialog.get_current_folder()
#      print 'current folder:', self.current_folder

      saveFileDialog.destroy()

      # Only process if 'OK' was clicked, not 'Cancel' - - - - - - - - - - - -
      #
      if (self.save_dialog_file_name == None):
        return

      print 'Selected file name in save dialog: %s' % \
            (self.save_dialog_file_name)

      # Make sure it is a Python file
      #
      if (self.save_dialog_file_name.endswith('.py') == False):
        self.save_dialog_file_name = self.save_dialog_file_name+'.py'

      self.project_name = self.save_dialog_file_name

      self.setWindowTitle()

    # Generate and save code only if a project file name is given - - - - - - -
    #
    if (self.project_name != None):

      # Check if file exists, if so ask if over-writing is OK
      #
      if (os.access(self.project_name, os.F_OK) == True):

        overwrite = self.messageDialog('File "%s" exists. ' % \
                                       (self.project_name) + \
                                       'Do you want to overwrite it?', 'ques')
        if (overwrite == False):
          return

      print 'Generate code and save it into file:', self.project_name

      febrl_project_code = self.generateAllCode()

      print 'Generated %d lines of code' % (len(febrl_project_code))

      try:
        f = open(self.project_name, 'w')  # Save into a text file
      except:
        self.messageDialog('Cannot write project file: %s' (self.project_name),
                           'error')
        return

      for line in febrl_project_code:

        # Replace 'XXXXXXXXXX' with project file name
        #
        if ('XXXXXXXXXX' in line):
          line = line.replace('XXXXXXXXXX',self.project_name.split(os.sep)[-1])

        f.write(line + os.linesep)
      f.close()

      for component in self.modified:  # Set all components to False
        self.modified[component] = False
      self.setWindowTitle()

  # ---------------------------------------------------------------------------
  # Handle activation of Online Help menu item
  #
  def activateOnlineHelp(self):
    print 'Online help activated - to be implemented...'

  # How to open the Web browser and an URL?

  # from Rattle:
  #on_rattle_menu_activate <- function(action, window)
  #{
  #browseURL("http://rattle.togaware.com")
  #}

  # ---------------------------------------------------------------------------
  # Handle activation of About menu item
  #
  def activateAbout(self, widget):
    print 'About activated'

    self.aboutTree = gtk.glade.XML(self.glade_file_name,
                                   self.about_dialog_name)
    self.aboutDialog = self.aboutTree.get_widget(self.about_dialog_name)
    self.aboutDialog.show()
    self.aboutDialog.run()
    self.aboutDialog.destroy()

    print 'About finished'

  # ---------------------------------------------------------------------------
  # Handle activation of License menu item
  #
  def activateLicense(self, widget):
    print 'License activated'

    self.licenseTree = gtk.glade.XML(self.glade_file_name,
                                     self.license_dialog_name)
    self.licenseDialog = self.licenseTree.get_widget(self.license_dialog_name)
    self.licenseDialog.show()
    self.licenseDialog.run()
    self.licenseDialog.destroy()

    print 'License finished'

  # ---------------------------------------------------------------------------
  # Generate a message dialog for information, warning, error, question or
  # cancel_question (c_question).
  # - Returns None, except for a question where it returns True (if Yes was
  #   answered) or False (No was answered), or None (if cancel was clicked).
  #
  def messageDialog(self, message_str, message_type):

    if (message_type[:4] not in ['info','warn','ques','erro','c_qu']):
      raise Exception, 'Illegal message type: %s' % (message_type)

    gtk_msg_type = {'info':gtk.MESSAGE_INFO,
                    'warn':gtk.MESSAGE_WARNING,
                    'ques':gtk.MESSAGE_QUESTION,
                    'c_qu':gtk.MESSAGE_QUESTION,
                    'erro':gtk.MESSAGE_ERROR}[message_type[:4]]

    # Message dialog buttons depend upon message type
    #
    gtk_msg_buttons = {'info':gtk.BUTTONS_CLOSE,
                       'warn':gtk.BUTTONS_CLOSE,
                       'erro':gtk.BUTTONS_CLOSE,
                       'c_qu':gtk.BUTTONS_CANCEL,
                       'ques':gtk.BUTTONS_YES_NO}[message_type[:4]]

    msgDialog = gtk.MessageDialog(None, gtk.DIALOG_MODAL, gtk_msg_type,
                                  gtk_msg_buttons, message_str)
    if (message_type[:4] == 'c_qu'):
      msgDialog.add_button(gtk.STOCK_NO,  gtk.RESPONSE_NO)
      msgDialog.add_button(gtk.STOCK_YES, gtk.RESPONSE_YES)

    msgDialog.set_transient_for(self.mainWin)
    msgDialog.set_position(gtk.WIN_POS_CENTER_ON_PARENT)

    response = msgDialog.run()
    msgDialog.destroy()

    if (message_type[:4] in ['ques','c_qu']):
      if (response == int(gtk.RESPONSE_YES)):
        ret = True
      elif (response == int(gtk.RESPONSE_NO)):
        ret = False
      else:
        ret = None

    else:
      ret = None
    return ret

  # ---------------------------------------------------------------------------
  # Set the title of the main window (add a '*' if modified)
  # - If title_str is None (default) the value from self.project_name is used
  # - If modified or not is taken from the self.modified flag
  #
  def setWindowTitle(self, title_str=None):

    print 'Set window title'

    if (title_str == None):
      if (self.project_name == None):
        win_project_name = '(None)'
      else:  # Extract file name only for display in main window
        win_project_name = self.project_name.split(os.sep)[-1]
    else:
      win_project_name = title_str

    if (sum(self.modified.values()) > 0):
      win_project_name += '*'

    print '  New project name:', win_project_name

    self.mainWin.set_title('Febrl - '+win_project_name)

  # ---------------------------------------------------------------------------
  # Write a text into the satus bar
  #
  def writeStatusBar(self, text):
    status_bar_widget = self.mainTree.get_widget('main_statusbar')

    status_bar_widget.pop(1)  # Remove previous message

    # Run a single iteration of the main loop to update the status bar
    #
    gtk.main_iteration(False)

    status_bar_widget.push(1, text)  # Push new message

    gtk.main_iteration(False)  # Two iterations to make sure message is shown
    gtk.main_iteration(False)

  # ===========================================================================
  # Methods that handle Data page events
  # ===========================================================================

  # ---------------------------------------------------------------------------
  # Display the Data page according to if data set(s) is/are initialised or not
  #
  def dataView(self):  # A switch to the Data page
    print '  Switched to Data page - Display this page'

    # Get type of project: Linkage (2 data sets) or others (1 data set)
    #
    if (self.project_type == 'Link'):

      data_set_index_suffix_list = [(0,'_a'), (1, '_b')]

      # Make everything for data set B visible
      #
      self.mainTree.get_widget('data_set_type_box_b').show()
      self.mainTree.get_widget('data_set_rec_id_box_b').show()
      self.mainTree.get_widget('data_set_notebook_b').show()
      self.mainTree.get_widget('data_set_scrolled_window_b').show()
      self.mainTree.get_widget("data_set_separator").show()

    else:  # For all other project types there is only one data

      data_set_index_suffix_list = [(0,'_a')]  # Only one data set

      # Make everything for data set B invisible
      #
      self.mainTree.get_widget('data_set_type_box_b').hide()
      self.mainTree.get_widget('data_set_rec_id_box_b').hide()
      self.mainTree.get_widget('data_set_notebook_b').hide()
      self.mainTree.get_widget('data_set_scrolled_window_b').hide()
      self.mainTree.get_widget('data_set_separator').hide()

    # Loop over data set index and suffix name - - - - - - - - - - - - - - - -
    #
    for (data_set_index, data_set_suffix) in data_set_index_suffix_list:

      # Get widget names and then their widgets for this data set
      #
      type_box_name =        'data_set_type_box'+data_set_suffix
      notebook_name =        'data_set_notebook'+data_set_suffix
      scrolled_window_name = 'data_set_scrolled_window'+data_set_suffix

      type_box_widget =        self.mainTree.get_widget(type_box_name)
      notebook_widget =        self.mainTree.get_widget(notebook_name)
      scrolled_window_widget = self.mainTree.get_widget(scrolled_window_name)

      notebook_widget.set_show_tabs(False)  # Hide tabs
      notebook_widget.hide()  # Make notebook invisible while setting it up

      # Get the type of this data set and switch to corresponding notebook page
      #
      data_set_type = self.data_set_type_list[data_set_index]
      page_num =      self.data_set_type_page_dict[data_set_type]
      notebook_widget.set_current_page(page_num)

      # Get the file chooser button widget
      #
      file_chooser_widget = self.mainTree.get_widget(data_set_type.lower() + \
                                       '_file_chooser_button'+ data_set_suffix)

      # For text based data sets set missing values, header line and strip
      # fields
      #
      if (data_set_type in ['CSV', 'COL', 'TAB']):
        miss_val_list = self.data_set_info_list[data_set_index]['miss_val']
        miss_val_widget = self.mainTree.get_widget('miss_val_entry' + \
                                                   data_set_suffix)
        miss_val_widget.set_text(','.join(miss_val_list))

        header_button_name = data_set_type.lower() + '_header_check_button' + \
                             data_set_suffix
        header_flag = self.data_set_info_list[data_set_index]['header_line']
        self.mainTree.get_widget(header_button_name).set_active(header_flag)

        # Set the strip fields check box according to stored setting
        #
        strip_button_name = data_set_type.lower() + '_strip_check_button' + \
                            data_set_suffix
        strip_flag = self.data_set_info_list[data_set_index]['strip_fields']
        self.mainTree.get_widget(strip_button_name).set_active(strip_flag)

      if (data_set_type == 'CSV'):  # Set the CSV delimiter
        delimiter_val = self.data_set_info_list[data_set_index]['delimiter']
        delimiter_widget = self.mainTree.get_widget('csv_delimiter_entry' + \
                                                    data_set_suffix)
        delimiter_widget.set_text(delimiter_val)

      # For text based data sets set the file name and file name filters
      #
      if (data_set_type in ['CSV', 'COL', 'TAB']):

        file_name = self.data_set_info_list[data_set_index]['file_name']

#### PC commented 4/07 ###########################
#        if (file_name != None):  # Set to previously selected file name
#          file_chooser_widget.set_filename(file_name)

#        if (self.current_folder != None):
#          file_chooser_widget.set_current_folder(self.current_folder)

        for filter in file_chooser_widget.list_filters():  # Remove old filters
          file_chooser_widget.remove_filter(filter)

        # Set new file filters appropriate to data set type
        #
        for filter_name in [data_set_type, 'TXT', 'ALL']:
          file_filter = gtk.FileFilter()  # Data set type specific filters
          data_set_filter = self.file_filters[filter_name]
          file_filter.set_name(data_set_filter[0])
          for filter_pattern in data_set_filter[1]:
            file_filter.add_pattern(filter_pattern)
          file_chooser_widget.add_filter(file_filter)

      # Set alignment for missing values label
      #
      miss_val_label_name = 'miss_val_label'+data_set_suffix
      self.mainTree.get_widget(miss_val_label_name).set_alignment(1.0, 0.5)

# commented PC 20/09

#      # Generate the four data set type button names
#      #
#      csv_button_name = 'data_set_csv_radio_button'+data_set_suffix
#      col_button_name = 'data_set_col_radio_button'+data_set_suffix
#      tab_button_name = 'data_set_tab_radio_button'+data_set_suffix
#      sql_button_name = 'data_set_sql_radio_button'+data_set_suffix

      # Deactivate the corresponding widgets first
      #
#      self.mainTree.get_widget(csv_button_name).set_active(False)
#      self.mainTree.get_widget(col_button_name).set_active(False)
#      self.mainTree.get_widget(tab_button_name).set_active(False)
#      self.mainTree.get_widget(sql_button_name).set_active(False)

      # Generate name of the active widget and then activate it
      #
#      active_button_name = 'data_set_' + data_set_type.lower() + \
#                           '_radio_button'+data_set_suffix
#      self.mainTree.get_widget(active_button_name).set_active(True)

      # View the data if it has been initalised - - - - - - - - - - - - - - - -
      #
      self.generateDataView(data_set_index)

      # Generate the combo box for the record identifier field - - - - - - - -
      #
      self.generateRecIDField(data_set_index)

      notebook_widget.show()  # Make widgets visible
      type_box_widget.show()
      scrolled_window_widget.show()

  # ---------------------------------------------------------------------------
  # Handle toggle of data set A and B type radio button
  #
  def dataSetTypeButtonToggle(self, widget):

    if (widget.get_active() == True):  # Check if this is the activated widget
      widget_name = widget.get_name()
      data_set_index = {'_a':0, '_b':1}[widget_name[-2:]]

      data_set_type_name = widget_name[9:12].upper()  # Get data set type name

      if (data_set_type_name not in self.data_set_type_page_dict.keys()):
        raise Exception, 'Illegal data set type widget name: %s' % \
                         (widget_name)

      self.data_set_type_list[data_set_index] = data_set_type_name

      print 'Data set type changed:', data_set_index, self.data_set_type_list

      # Set default information for this data set (keep missing values list)
      #
      data_set_info_dict = self.data_set_info_list[data_set_index]
      default_vals = self.data_set_default_values

      data_set_info_dict['header_line'] =  default_vals['header_line']
      data_set_info_dict['strip_fields'] = default_vals['strip_fields']
      data_set_info_dict['header_data'] =  default_vals['header_data']
      data_set_info_dict['file_name'] =    default_vals['file_name']
      data_set_info_dict['file_data'] =    default_vals['file_data']
      data_set_info_dict['field_names'] =  default_vals['field_names']
      data_set_info_dict['rec_id_field'] = default_vals['rec_id_field']

      if (data_set_type_name == 'CSV'):
        data_set_info_dict['delimiter'] = ','
      elif (data_set_type_name == 'TAB'):
        data_set_info_dict['delimiter'] = '\t'

      self.list_data_store[data_set_index] = None

      # Unselect file name in file chooser of all data set types
      #
      for data_set_type_name in ['csv','col','tab']:
        file_chooser_widget = self.mainTree.get_widget(data_set_type_name + \
                                     '_file_chooser_button' + widget_name[-2:])
        file_chooser_widget.unselect_all()

      # Make all pages invisible except Data and Log - - - - - - - - - - - - -
      #
      for page_name in ['Explore', 'Standardise','Index', 'Compare',
                        'Classify', 'Run', 'Evaluate', 'Review']:
        self.main_notebook_page_active_dict[page_name] = False

      self.modified['data'] = True  # Data set type has changed

      self.index_def = [] # Delete all previous index definitions
      self.index_num = 0

      self.field_comp_list = []  # Delete previous field comparison functions

      self.setWindowTitle()

      self.displayCurrentNotebookPage()  # Re-display current notebook page

  # ---------------------------------------------------------------------------
  # Handle toggle of the header line check buttons
  #
  def dataSetTypeHeaderLineToggle(self, widget):
    print 'Header line toggle'

    widget_name = widget.get_name()
    data_set_ind = {'_a':0, '_b':1}[widget_name[-2:]]

    self.data_set_info_list[data_set_ind]['header_line'] = widget.get_active()

    # Clear the field names and record id field
    #
    self.data_set_info_list[data_set_ind]['field_names'] = \
                                    self.data_set_default_values['field_names']
    self.data_set_info_list[data_set_ind]['rec_id_field'] = \
                                   self.data_set_default_values['rec_id_field']

    self.index_def = []  # Clear all index definitions as they use field names
    self.index_num = 0

    self.comp_std_list = []  # Clear all component standardisers

    self.field_comp_list = []  # And clear field comparisons

    self.modified['data'] = True  # Header line flag has changed
    self.setWindowTitle()

    self.displayCurrentNotebookPage()  # Re-display current notebook page

  # ---------------------------------------------------------------------------
  # Handle toggle of the strip fields check buttons
  #
  def dataSetTypeStripFieldsToggle(self, widget):
    print 'Strip fields toggle'

    widget_name = widget.get_name()
    data_set_ind = {'_a':0, '_b':1}[widget_name[-2:]]

    self.data_set_info_list[data_set_ind]['strip_fields'] = widget.get_active()
    self.modified['data'] = True  # Strip fields flag has changed
    self.setWindowTitle()

    self.displayCurrentNotebookPage()  # Re-display current notebook page

  # ---------------------------------------------------------------------------
  # Handle changes in the missing values text entry field(s)
  #
  def dataSetMissValueActivate(self, widget):
    print 'Missing values activated'

    widget_name = widget.get_name()
    data_set_index = {'_a':0, '_b':1}[widget_name[-2:]]

    miss_val_list = []  # Create a list of values stripped of whitespaces

    for val in widget.get_text().split(','):
      miss_val_list.append(val.strip())

    self.data_set_info_list[data_set_index]['miss_val'] = miss_val_list

    self.modified['data'] = True  # Missing values have changed
    self.setWindowTitle()

    self.displayCurrentNotebookPage()  # Re-display current notebook page

  # ---------------------------------------------------------------------------
  # Handle changes in the delimiter text entry field(s)
  #
  def dataSetDelimiterActivate(self, widget):
    print 'Delimiter activated'

    widget_name = widget.get_name()
    data_set_index = {'_a':0, '_b':1}[widget_name[-2:]]

    delimiter_val = widget.get_text()
    delimiter_val = delimiter_val.strip()  # Remove sourrounding whitespace

    # The delimiter string must not be empty or contain more than one character
    #
    if ((len(delimiter_val) > 1) or (delimiter_val == '')):
      self.messageDialog('Delimiter must be a one-character string.','warn')

      # Set delimiter back to default value
      #
      self.data_set_info_list[data_set_index]['delimiter'] = \
                                      self.data_set_default_values['delimiter']
    else:
      self.data_set_info_list[data_set_index]['delimiter'] = delimiter_val

    # Reset field names back to default
    #
    self.data_set_info_list[data_set_index]['field_names'] = \
                                    self.data_set_default_values['field_names']

    self.modified['data'] = True  # Delimiter value has changed
    self.setWindowTitle()

    self.displayCurrentNotebookPage()  # Re-display current notebook page

  # ---------------------------------------------------------------------------
  # Handle selection of a file name (file name buttons)
  #
  def dataFileSelect(self, widget):
    print 'File selected, name:', widget.get_filename()

    widget_name = widget.get_name()
    data_set_index = {'_a':0, '_b':1}[widget_name[-2:]]

    file_name = widget.get_filename()  # Get the file name returned

    if (file_name != None):  # Do nothing if no file name is returned

      self.current_folder = widget.get_current_folder()
      print 'current folder:', self.current_folder

      data_set_type_name = widget_name[:3].upper()  # First three characters

      # Save new file name and set header data to defaults
      #
      self.data_set_info_list[data_set_index]['file_name'] = file_name
      self.data_set_info_list[data_set_index]['file_data'] = \
                                      self.data_set_default_values['file_data']
      # Commented PC 12/07/07, so to keep previous selections!
      #self.data_set_info_list[data_set_index]['header_line'] = \
      #                             self.data_set_default_values['header_line']

      # Set field names to default to re-generate them for new data set - - - -
      # (PC, 7/08/07)
      #
      self.data_set_info_list[data_set_index]['field_names'] = \
                                   self.data_set_default_values['field_names']

      # Also set record identifier field to default
      #
      self.data_set_info_list[data_set_index]['rec_id_field'] = \
                                   self.data_set_default_values['rec_id_field']

      ##### added PC 30/08/2007
      self.field_comp_list = []  # Remove old field comparisons

      self.comp_std_list = []  # Clear all component standardisers

      self.modified['data'] = True  # Data set file name has changed
      self.setWindowTitle()

      # Reset match file names for output - - - - - - - - - - - - - - - - - - -
      # (keep active flag)
      #
      self.output_dict['m_datasets'] = (self.output_dict['m_datasets'][0],
                                        ('(None)', 'match_id'),
                                        ('(None)', 'match_id'))

      self.getDataSetFromFile(data_set_index)  # Initialise data set

      self.generateDataView(data_set_index)  # Generate the data view

      # Generate the combo box for the record identifier field - - - - - - - -
      #
      self.generateRecIDField(data_set_index)

  # ---------------------------------------------------------------------------
  # Get the first several lines from the data set
  #
  def getDataSetFromFile(self, data_set_index):
    print 'Get data set data'

    file_name = self.data_set_info_list[data_set_index]['file_name']

    # Check if the file is compressed or not
    #
    if (file_name[-3:] in ['.gz', '.GZ']):
      try:
        fp = gzip.open(file_name)
      except:
        fp = 'error'
    else:
      try:
        fp = open(file_name)
      except:
        fp = 'error'

    if (fp == 'error'):
      self.messageDialog('Cannot open file: %s' % (file_name), 'erro')
      return

    file_data = []
    for i in range(self.num_data_rows+1):
      line = fp.readline()
      line = line.rstrip('\n\r')  # Remove line separators ('\n' and/or '\r')
      if (line != ''):  # Not yet reached end of file
        file_data.append(line)
      else:
        break
    fp.close()

    self.data_set_info_list[data_set_index]['file_data'] = file_data

    self.writeStatusBar('Read first %d lines from data set: %s.' % \
                       (self.num_data_rows+1, file_name))

  # ---------------------------------------------------------------------------
  # Generate the combo box for the record identifier field
  #
  def generateRecIDField(self, data_set_index) :
    print 'generate record identifier field combo box'

    data_set_suffix = {0:'_a', 1:'_b'}[data_set_index]

    rec_id_box_name =   'data_set_rec_id_box'+data_set_suffix
    rec_id_box_widget = self.mainTree.get_widget(rec_id_box_name)

    # Generate combo box for record identifier field - - - - - - - - - - - -
    #
    if (len(rec_id_box_widget.get_children()) > 1):  # Remove old combo box
      rec_id_box_widget.remove(rec_id_box_widget.get_children()[1])

    rec_id_field_combo_box = gtk.combo_box_new_text()  # New combo box
    rec_id_field_combo_box.append_text('(None)')

    field_names = self.data_set_info_list[data_set_index]['field_names']
    default_field_names = self.data_set_default_values['field_names']

    if (field_names != default_field_names): # Field names available
      for field_name in field_names:
        rec_id_field_combo_box.append_text(field_name)

    # Set the active field name
    #
    rec_id_field_name = self.data_set_info_list[data_set_index]['rec_id_field']
    if ((rec_id_field_name == self.data_set_default_values['rec_id_field']) or \
        (rec_id_field_name == self.data_set_default_values['rec_id_field'] + \
                              data_set_suffix+'__')):
      rec_id_field_combo_box.set_active(0)  # Default, not a data set field
    else:
      assert rec_id_field_name in field_names  # Check its in field names list

      rec_id_field_index = field_names.index(rec_id_field_name)
      rec_id_field_combo_box.set_active(rec_id_field_index+1)

    rec_id_box_widget.pack_start(rec_id_field_combo_box, False, False, 0)
    rec_id_field_combo_box.show()

  # ---------------------------------------------------------------------------
  # Generate the data view (for either data set A or B)
  #
  def generateDataView(self, data_set_index):
    print 'Generate data set view'

    raw_file_data = self.data_set_info_list[data_set_index]['file_data']

    data_set_suffix = {0:'_a', 1:'_b'}[data_set_index]

    # Get the tree view widget
    #
    tree_view_widget = self.mainTree.get_widget('data_set_treeview' + \
                       data_set_suffix)

    if (raw_file_data == None):  # If no data set initalised, un-set model - -
      tree_view_widget.set_model()
      return  # Nothing more to do

    # Get the data set type, and header line and strip fields flags - - - - - -
    #
    data_set_type = self.data_set_type_list[data_set_index]
    header_line =   self.data_set_info_list[data_set_index]['header_line']
    strip_fields =  self.data_set_info_list[data_set_index]['strip_fields']
    field_names =   self.data_set_info_list[data_set_index]['field_names']

    file_data = []  # List with the cleaned data lines

    # Process the raw file data according to data set type
    #
    if (data_set_type in ['CSV','TAB']):  # - - - - - - - - - - - - - - - - - -

      if (data_set_type == 'TAB'):
        delimiter_val = '\t'
      else:
        delimiter_val = self.data_set_info_list[data_set_index]['delimiter']

      for row in csv.reader(raw_file_data, delimiter = delimiter_val):
        if (strip_fields == True):
          clean_row = []
          for val in row:
            clean_row.append(val.strip())  # Remove sourrounding whitespace
        else:
          clean_row = row
        file_data.append(clean_row)

    elif (data_set_type == 'COL'):  # - - - - - - - - - - - - - - - - - - - -

      # Get width of the column width values from GUI
      #
      col_width_name = 'col_width_entry'+data_set_suffix
      col_width_vals = self.mainTree.get_widget(col_width_name).get_text()
      col_width_vals = col_width_vals.strip()  # Remove sourrounding whitespace

      col_width_list = []  # Convert into a list of integers
      for val in col_width_vals.split(','):
        val = val.strip()

        if ((not val.isdigit()) or (int(val < 1))):
          self.messageDialog('Column widths must be positive integers.',
                             'error')
          return
        col_width_list.append(int(val))

      for row in raw_file_data:
        start = 0
        clean_row = []
        for width in col_width_list:
          val = row[start:start+width]
          if (strip_fields == True):
            val = val.strip()  # Remove sourrounding whitespace
          clean_row.append(val)
          start = start+width
        file_data.append(clean_row)

    elif (data_set_type == 'SQL'):  # - - - - - - - - - - - - - - - - - - -

      print 'SQL data set to be implemented' #############################

    else:
      raise Exception,  'Illegal data set type: %s' % (data_set_type)

    # Now check if any of these values are designated as missing values
    #
    miss_val_list = self.data_set_info_list[data_set_index]['miss_val']

    for i in range(len(file_data)):
      data_list = file_data[i]

      for j in range(len(data_list)):
        data_val = data_list[j]
        for m in miss_val_list:
          if (m == data_val):
            data_val = data_val.replace(m, '')  # Remove that missing value
            data_list[j] = data_val
      file_data[i] = data_list

    # Get the number of fields in the data - - - - - - - - - - - - - - - - -
    #
    num_fields = len(file_data[0])

    # Check if all rows have same number of fields
    #
    for i in range(1,len(file_data)):  ## PC 26/09 self.num_data_rows+1):
      if (len(file_data[i]) != num_fields):
        self.messageDialog('Different number of fields in data in line %d.' \
                           % (i), 'erro')
        return  # Do not initialise data

## changed below, PC 15/07 ########################

    # Generate header line data - - - - - - - - - - - - - - - - - - - - - - -
    #
    default_field_names = self.data_set_default_values['field_names']

    if (field_names != default_field_names):  # Field names have been set
      header_data = field_names  # Use the previously set field names

    else:  # Field names are not set, so generate them

      if (header_line == True):  # Get header line from file data
        header_data = file_data[0]

##        self.data_set_info_list[data_set_index]['header_data'] = header_data
##        self.data_set_info_list[data_set_index]['field_names'] = header_data

      else:  # Get previous header line data or generate field names for it
        header_data = []
        for i in range(num_fields):
          header_data.append('field-%d' % (i))

      ## PC 17/08/07 ############
      self.data_set_info_list[data_set_index]['field_names'] = header_data

    # Generate list store for GUI - - - - - - - - - - - - - - - - - - - - - -

    # Remove all previous columns from tree view
    #
    for col in tree_view_widget.get_columns():
      tree_view_widget.remove_column(col)

    tree_view_widget.set_model(None)  # Remove the previous model

    # One column per data set field plus a boolean column for editing and
    # one column for the text weight (bold or not - 400 is normal weight)
    #
    columns = [str]*num_fields+[bool,int]

    data_store = gtk.ListStore(*columns)  # (see PyGTK FAQ 13.10)

    print 'culumns', columns
    print 'data_store', data_store

    # If header line is ticked (taken from file) don't allow editing of names
    #
    if (header_line == True):
      data_store_extra = [False,600]  # Bold face font for header line
    else:
      data_store_extra = [True,600]  # Allow editing of field names

    data_store.append(header_data+data_store_extra)

    if (header_line == True):
      for i in range(1,len(file_data)): # PC 26/09 self.num_data_rows+1):
        data_store.append(file_data[i]+[False,400])  # No editing, not bold
    else:  # Add PC 21/11/2007
      for i in range(0,len(file_data)-1): # First line is already a data record
        data_store.append(file_data[i]+[False,400])  # No editing, not bold

    tree_view_widget.set_headers_visible(False)  # Don't show headers

    ###########################################################################
    ## PC 7/09/2007: Doesn't work with gtk-2.6.8
    ## tree_view_widget.set_grid_lines(gtk.TREE_VIEW_GRID_LINES_BOTH)
    ###########################################################################

    tree_view_widget.set_rules_hint(True)  # Use this instead of lines

    # Remove all columns from tree view
    #
    for col in tree_view_widget.get_columns():
      tree_view_widget.remove_column(col)

    # Connect to the GUI tree view widget
    #
    tree_view_widget.set_model(data_store)

    self.list_data_store[data_set_index] = data_store  # Save for later

    # Create the tree view columns, cells and renderers - - - - - - - - - - -
    #
    for i in range(num_fields):
      cell_renderer = gtk.CellRendererText()
      cell_renderer.connect('edited', self.headerLineEdit,(data_set_index,i))

      tv_column = gtk.TreeViewColumn('Col %d' % (i), cell_renderer, text=i,
                                   editable=num_fields, weight=num_fields+1)

      tree_view_widget.append_column(tv_column)

  # ---------------------------------------------------------------------------
  # A value of a header line entry has been edited (in either data set A or B)
  #
  def headerLineEdit(self, cell, path, new_value, user_data):

    new_value = new_value.strip()  # Remove whitespaces

    if (new_value == ''):
      self.messageDialog('Field names must not be empty!', 'warning')
      return  # Don't do anything

    data_set_index, col_num = user_data  # Unpack values
    data_store = self.list_data_store[data_set_index]

    # Check if new field names is different from all others
    #
    num_fields = len(self.list_data_store[data_set_index][0])-2

    for c in range(num_fields):
      if ((c != col_num) and (new_value == data_store[path][c])):
        self.messageDialog('All field names must be different!', 'warning')
        return  # Don't do anything


    data_store[path][col_num] = new_value

    print 'Header line column %d edited for data set %d, new value: %s' %\
          (col_num, data_set_index, new_value)

    self.comp_std_list = []  # Clear all component standardisers

    self.field_comp_list = []  # And clear field comparisons

    self.modified['data'] = True  # Header line flag has changed
    self.setWindowTitle()

    # Possibly update record identifier fields???
    # possibly updated field_names??? ##################

  # ---------------------------------------------------------------------------
  # Handle an activate of Execute on the Data page
  #
  def dataExecute(self):
    data_store = self.list_data_store

    self.febrl_code['data'] = None  # Remove all previous data code

    # First check if necessary data sets are initialised - - - - - - - - - - -
    #
    if ((self.project_type == 'Link') and \
        ((data_store[0] == None) or (data_store[1] == None))):
      self.messageDialog('Two data sets must be initalised for a linkage!',
                           'warning')
      return  # Don't do anything

    elif (data_store[0] == None):
      self.messageDialog('No data set initalised!', 'warning')
      return  # Don't do anything

    # The list of strings, each element a line of Python code or comments
    #
    febrl_code = []
    febrl_code.append('# '+'-'*77)
    febrl_code.append('')

    # Get necessary information about data set(s) - - - - - - - - - - - - - - -
    #
    if (self.project_type == 'Link'):
      data_set_index_suffix_list = [(0,'_a'), (1, '_b')]
    else:
      data_set_index_suffix_list = [(0,'_a')]

    # Loop over data set index and suffix name - - - - - - - - - - - - - - - -
    #
    for (data_set_index, data_set_suffix) in data_set_index_suffix_list:

      data_set_info = self.data_set_info_list[data_set_index]
      data_set_type = self.data_set_type_list[data_set_index]

      data_set_type_name = data_set_type
      if (data_set_type == 'TAB'):
        data_set_type_name = 'CSV'  # TAB data set is a CSV data set with \t
                                    # delimiter

      # Get the field names from the data store
      #
      data_set_store = self.list_data_store[data_set_index]
      store_iter = data_set_store.get_iter_first()

      field_names = []  # Get the field names

      # Remember the last two fields are used for editable (True/False) and
      # font weight (for headerline)
      #
      num_fields = data_set_store.get_n_columns() - 2

      for field_num in range(num_fields):
        field_names.append(data_set_store.get_value(store_iter, field_num))

      data_set_info['field_names'] = field_names  # Store field names

      febrl_code.append('# Define input data set %s:' % \
                        data_set_suffix[-1].upper())
      febrl_code.append('#')

      febrl_code.append('data_set_%s = dataset.DataSet' % \
                        (data_set_suffix[-1]) + data_set_type_name + \
                        '(description="Data set generated by Febrl GUI",')

      if (data_set_type == 'TAB'):  # Add a comment for TAB data set
        last_line = febrl_code[-1]
        last_line = last_line[:-2]+' (note TAB data set is implemented as ' + \
                                   'CSV data set with \\t delimiter)",'
        febrl_code[-1] = last_line

      indention_space = ' '*32

      # Add data set type independent arguments - - - - - - - - - - - - - - - -
      #
      febrl_code.append(indention_space+'access_mode="read",')

      strip_fields = data_set_info['strip_fields']
      febrl_code.append(indention_space+'strip_fields=%s,' % \
                        (str(strip_fields)))

      miss_val_list = data_set_info['miss_val']
      febrl_code.append(indention_space+'miss_val=%s,' %
                        (str(miss_val_list)))

      # Get the index for the record identifier field - - - - - - - - - - - - -
      #
      rec_id_box_name =   'data_set_rec_id_box'+data_set_suffix
      rec_id_box_widget = self.mainTree.get_widget(rec_id_box_name)

      # Get and save the selected field index
      #
      rec_id_field_index = rec_id_box_widget.get_children()[1].get_active()

      if (rec_id_field_index <= 0):  # Set to default value
        rec_id_str = self.data_set_default_values['rec_id_field'] + \
                     data_set_suffix+'__'

        # Check identifier field name it is not used
        #
        while (rec_id_str in field_names):
          rec_id_str = '_'+rec_id_str+'_'

        data_set_info['rec_id_field'] = rec_id_str

      else:  # Get name from field names list, adjust for "(None): at index 0
        data_set_info['rec_id_field'] = field_names[rec_id_field_index-1]
        rec_id_str = data_set_info['rec_id_field']

      febrl_code.append(indention_space+'rec_ident="%s",' % (rec_id_str))

      # Now add data set type dependent arguments - - - - - - - - - - - - - - -
      #
      if (data_set_type in ['CSV','COL','TAB']):

        file_name = data_set_info['file_name']
        febrl_code.append(indention_space+'file_name="%s",' % (file_name))

        header_flag = data_set_info['header_line']
        febrl_code.append(indention_space+'header_line=%s,' % \
                          (str(header_flag)))

      if (data_set_type in ['CSV','TAB']):
        delimiter_val = data_set_info['delimiter']
        if (delimiter_val == '\t'):  # Special handling needed for tab
          delimiter_val= '\\t'  # So it is printed correctly

        febrl_code.append(indention_space+'delimiter="%s",' % \
                          (delimiter_val))

      # Add field names - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      #
      if (data_set_type in ['CSV','TAB']):

        for field_num in range(num_fields):
          if (field_num == 0):
            code_str = 'field_list = [("%s",%d),' % \
                       (field_names[field_num], field_num)
          elif (field_num == num_fields-1):  # Last column
            code_str = ' '*14 + '("%s",%d)],' % \
                       (field_names[field_num], field_num)
          else:
            code_str = ' '*14 + '("%s",%d),' % \
                       (field_names[field_num], field_num)
          febrl_code.append(indention_space+code_str)

      elif (data_set_type == 'COL'):

        # Get width of the column width values from GUI
        #
        col_width_name = 'col_width_entry'+data_set_suffix
        col_width_vals = self.mainTree.get_widget(col_width_name).get_text()
        col_width_vals = col_width_vals.strip()
        col_width_list = col_width_vals.split(',')

        for field_num in range(num_fields):
          col_width = col_width_list[field_num]
          col_width = col_width.strip()
          field_name = field_names[field_num]

          if (field_num == 0):
            code_str = 'field_list = [("%s",%s),' % (field_name, col_width)
          elif (field_num == num_fields-1):  # Last column
            code_str = ' '*14 + '("%s",%s)],' % (field_name, col_width)
          else:
            code_str = ' '*14 + '("%s",%s),' % (field_name, col_width)
          febrl_code.append(indention_space+code_str)

      # Add closing bracked to the last generated line - - - - - - - - - - - -
      #
      last_line = febrl_code[-1]
      last_line = last_line[:-1]+')'  # Replace last comma with bracket
      febrl_code[-1] = last_line

      febrl_code.append('')

    febrl_code.append('')

    self.febrl_code['data'] = febrl_code  # Store for later use

    # Finally update the GUI information - - - - - - - - - - - - - - - - - - -
    #
    self.addToLog('')  # Add generated code into log page text
    self.addToLog('='*79)
    self.addToLog('Generated Febrl code for "data" on %s' % (time.asctime()))
    self.addToLog('')
    for line in febrl_code:
      self.addToLog(line)
    self.addToLog('')

    self.modified['data'] = True  # Data set details have been changed
    self.setWindowTitle()

    self.re_run['data_init'] = True  # Need to re-run data set initialisation

    # Update the active and non-active notebook pages
    #
    self.main_notebook_page_active_dict['Explore'] =  True

    self.main_notebook_page_active_dict['Classify'] = False
    self.main_notebook_page_active_dict['Run'] =      False
    self.main_notebook_page_active_dict['Evaluate'] = False

    if (self.project_type in ['Link','Deduplicate']):
      self.main_notebook_page_active_dict['Index'] = True
      self.main_notebook_page_active_dict['Compare'] = True
    elif (self.project_type == 'Standardise'):
      self.main_notebook_page_active_dict['Standardise'] = True
    else:  # Geocoding - TODO
      print 'Geocoding not implemented yet'

    self.index_def = [] # Delete all previous index definitions
    self.index_num = 0

    self.displayCurrentNotebookPage()  # Diplay the current page

    self.writeStatusBar('Generated Febrl Python code for data set ' + \
                        'initialisation (see Log page for generated code).')

  # ===========================================================================
  # Methods that handle Explore page events
  # ===========================================================================

  # ---------------------------------------------------------------------------
  # Display the Explore page
  #
  def exploreView(self):
    print '  Switched to Explore page - Display this page'

    # Set font in text view to fixed width
    #
    explore_textview = self.mainTree.get_widget('explore_text_view')
    explore_textview.modify_font(pango.FontDescription("courier 10"))

    explore_widget_name_b = 'explore_data_set_b_check_button'
    explore_widget_b = self.mainTree.get_widget(explore_widget_name_b)

    # If project is a linkage make second data set tick box sensitive
    #
    if (self.project_type == 'Link'):
      explore_widget_b.set_sensitive(True)
    else:  # Only one data set
      explore_widget_b.set_sensitive(False)

  # ---------------------------------------------------------------------------
  # Handle an activate of Execute on the Explore page
  #
  def exploreExecute(self):
    analyse_data_set = [False, False]  # Which data sets to analyse

    # Get values from select data set tick box(es) - - - - - - - - - - - - - -
    #
    explore_widget_name_a = 'explore_data_set_a_check_button'
    explore_widget_a = self.mainTree.get_widget(explore_widget_name_a)

    if (explore_widget_a.get_active() == True):
       analyse_data_set[0] = True

    if (self.project_type == 'Link'):  # Check second data set
      explore_widget_name_b = 'explore_data_set_b_check_button'
      explore_widget_b = self.mainTree.get_widget(explore_widget_name_b)

      if (explore_widget_b.get_active() == True):
         analyse_data_set[1] = True

    print analyse_data_set

    # Get analysis mode (values or words) - - - - - - - - - - - - - - - - - - -
    #
    values_widget = self.mainTree.get_widget('explore_values_radio_button')
    word_analysis = not values_widget.get_active()

    print 'word analysis:', word_analysis

    # Check if sampling activated and if so get sampling rate - - - - - - - - -
    #
    use_sample_widget = self.mainTree.get_widget('explore_sample_check_button')
    if (use_sample_widget.get_active() == True):
      sample_rate = self.mainTree.get_widget('explore_sample_entry').get_text()

      try:
        sample_rate_val = float(sample_rate)
      except:
        self.messageDialog('Sampling rate must be a positive percentage ' + \
                           'value.', 'error')
        return

      if ((sample_rate_val <= 0.0) or (sample_rate_val > 100.0)):
        self.messageDialog('Sampling rate must be a positive percentage ' + \
                           'value', 'error')
        return

    else:
      sample_rate_val = 100.0

    # Set up the text view and buffer - - - - - - - - - - - - - - - - - - - - -
    #
    explore_textview = self.mainTree.get_widget('explore_text_view')
    explore_buffer =   explore_textview.get_buffer()

    explore_text = ''

    # Initialise data set(s) - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    data_set_index = 0

    for analysis_flag in analyse_data_set:  # Loop over both data sets

      if (analysis_flag == True):  # Analyse this data set

        data_set_info = self.data_set_info_list[data_set_index]
        data_set_type = self.data_set_type_list[data_set_index]

        # Get the data set field names
        #
        field_names = self.data_set_info_list[data_set_index]['field_names']
        num_fields = len(field_names)

#        data_set_store = self.list_data_store[data_set_index]
#        store_iter = data_set_store.get_iter_first()
#
#        field_names = []  # Get the field names
#
#        # Remember the last two fields are used for editable (True/False) and
#        # font weight (for headerline)
#        #
#        num_fields = data_set_store.get_n_columns() - 2
#
#        for field_num in range(num_fields):
#          field_names.append(data_set_store.get_value(store_iter, field_num))

        header_line_flag =  data_set_info['header_line']
        strip_fields_flag = data_set_info['strip_fields']
        miss_val_list =     data_set_info['miss_val']
        rec_id_str =        data_set_info['rec_id_field']

#        if (self.data_set_type_list[data_set_index] in ['CSV','TAB']):
#          if (self.data_set_type_list[data_set_index] == 'TAB'):
#            delimiter_val = data_set_info['delimiter']

        # Create data set object according to data set type - - - - - - - - - -
        #
        if (data_set_type in ['CSV','TAB']):
          if (data_set_type == 'TAB'):
            delimiter_val = '\t'
          else:
            delimiter_val = data_set_info['delimiter']

          file_name_val = data_set_info['file_name']

          self.writeStatusBar('Initialise CSV data set "%s" for reading.' % \
                              (file_name_val))
          field_val_list = []
          for field_num in range(num_fields):
            field_val_list.append((field_names[field_num], field_num))

          data_set = dataset.DataSetCSV(file_name = file_name_val,
                                        rec_ident = rec_id_str,
                                        access_mode = 'read',
                                        header_line = header_line_flag,
                                        strip_fields = strip_fields_flag,
                                        delimiter = delimiter_val,
                                        miss_val = miss_val_list,
                                        field_list = field_val_list)

        elif (data_set_type == 'COL'):  # - - - - - - - - - - - - - - - - - - -

          file_name_val = data_set_info['file_name']
          self.writeStatusBar('Initialise COL data set "%s" for reading.' % \
                              (file_name_val))

          if (data_set_index == 0):
            col_width_name = 'col_width_entry_a'
          else:
            col_width_name = 'col_width_entry_b'

          # Get width of the column width values from GUI
          #
          col_width_vals = self.mainTree.get_widget(col_width_name).get_text()
          col_width_vals = col_width_vals.strip()
          col_width_list = col_width_vals.split(',')

          field_val_list = []
          for field_num in range(num_fields):
            if (header_line_flag == True):
              field_val_list.append(int(col_width_list[field_num]))
            else:
              field_val_list.append((field_names[field_num],
                                     int(col_width_list[field_num])))

          data_set = dataset.DataSetCOL(file_name = file_name_val,
                                        rec_ident = rec_id_str,
                                        access_mode='read',
                                        header_line = header_line_flag,
                                        strip_fields = strip_fields_flag,
                                        miss_val = miss_val_list,
                                        field_list = field_val_list)

        elif (data_set_type == 'SQL'):

          pass  ####### TO DO ##########################

        # Start file analysis - - - - - - - - - - - - - - - - - - - - - - - - -
        # (report progess into status bar every 1000 records)
        #
        analysis_str_list = data_set.analyse(sample_rate_val, word_analysis,
                                             self.writeStatusBar, 1000)

        # Write results into explore text view - - - - - - - - - - - - - - - -
        #
        for line in analysis_str_list:
          explore_text += line
          explore_text += os.linesep

        # Check for valid UTF-8/ASCII unicodes - - - - - - - - - - - - - - - -
        #
        try:
          clean_explore_text = unicode(explore_text)
        except:  # Have to manually remove the illegal characters
          clean_explore_text = '**** Note some non-ASCII characters were ' + \
                               'found in the explore output and replaced ' + \
                               'with "*" ****'+os.linesep
          for c in explore_text:
            if (ord(c) < 128):
              clean_explore_text += c
            else:
              clean_explore_text += '*'

        explore_buffer.set_text(clean_explore_text)

        data_set.finalise()

      data_set_index += 1

  # ===========================================================================
  # Methods that handle Standardise page events
  # ===========================================================================

  # ---------------------------------------------------------------------------
  # Display the Standardise page
  #
  def standardView(self):
    print '  Switched to Standardise page - Display this page'

    std_scrolled_window = self.mainTree.get_widget('standard_scrolled_window')

    # Create box for the component standardiser definitions in scrolled window
    #
    std_viewport = std_scrolled_window.get_child()

    if (std_viewport == None):  # Create new box if none is there yet
      std_box = gtk.VBox()
      std_scrolled_window.add_with_viewport(std_box)
    else:  # Get the box with component std. details (1ast child is viewport)
      std_viewport.show()
      std_box = std_scrolled_window.get_child().get_child()

    for child in std_box.get_children():  # Remove all old children
      std_box.remove(child)

    num_comp_std = len(self.comp_std_list)

    # A dictionary with number of different component standardisers (for adding
    # numbers of output fields)
    #
    self.comp_std_type_count_dict = {'Date':0,'Name':0,'Phon':0,'Addr':0}

    # Now show existing component standardisers - - - - - - - - - - - - - - - -
    #
    for comp_std in self.comp_std_list:
      print 'comp_std:', comp_std

      # A comp_std is a dictionary with keys: type, in_field_list,
      # out_field_list, plus type specific parameters
      #
      comp_std_type = comp_std['type']

      self.comp_std_type_count_dict[comp_std_type]

      comp_std_type_str = {'Date':'Date', 'Phon':'Phone number', 'Name':'Name',
                           'Addr':'Address'}[comp_std_type]
      comp_std_type_box = gtk.HBox()
      comp_std_type_label = gtk.Label('<b>%s standardiser:</b>' % \
                                      (comp_std_type_str))
      comp_std_type_label.set_use_markup(True)
      comp_std_type_box.pack_start(comp_std_type_label, False, False, 0)
      comp_std_type_label.show()
      std_box.pack_start(comp_std_type_box, False, False, 0)
      comp_std_type_box.show()

      comp_std_box = self.new_comp_std(comp_std_type, comp_std)
      std_box.pack_start(comp_std_box, False, False, 0)
      comp_std_box.show()

      horiz_separator = gtk.HSeparator()  # Display a horizontal separator
      std_box.pack_start(horiz_separator, False, False, 5)
      horiz_separator.show()

    # Finally add 'Add' and 'Delete' buttons - - - - - - - - - - - - - - - -
    #
    if (num_comp_std == 0): # No component standardiser so far, only show 'Add'
      del_comp_std = False
    else:
      del_comp_std = True

    comp_std_button_box = self.create_comp_std_button(del_comp_std)

    std_box.pack_start(comp_std_button_box, False, False, 5)
    comp_std_button_box.show()
    std_box.show()

    std_scrolled_window.show()

  # ---------------------------------------------------------------------------
  # Construct a new HBox with buttons to add component standardisers and
  # (if more than one exist, i.e. del_comp_std is set to True, also add a
  # delete button
  # Returns the HBox with the buttons.
  # If the 'del_comp_std' flag is set to False then the corresponding button
  # will not be created
  #
  def create_comp_std_button(self, del_comp_std=True):

    comp_std_button_box = gtk.HBox()

    sep_label = gtk.Label('Add new component standardiser for: ')
    comp_std_button_box.pack_start(sep_label, False, False, 0)
    sep_label.show()

    add_date_button = gtk.Button('Dates')
    add_date_button.connect('clicked', self.clickAddStandardiserButton,
                            comp_std_button_box)
    comp_std_button_box.pack_start(add_date_button, False, False, 5)
    add_date_button.show()

    add_phone_button = gtk.Button('Phone numbers')
    add_phone_button.connect('clicked', self.clickAddStandardiserButton,
                             comp_std_button_box)
    comp_std_button_box.pack_start(add_phone_button, False, False, 5)
    add_phone_button.show()

    add_name_button = gtk.Button('Names')
    add_name_button.connect('clicked', self.clickAddStandardiserButton,
                            comp_std_button_box)
    comp_std_button_box.pack_start(add_name_button, False, False, 5)
    add_name_button.show()

    add_addr_button = gtk.Button('Addresses')
    add_addr_button.connect('clicked', self.clickAddStandardiserButton,
                            comp_std_button_box)
    comp_std_button_box.pack_start(add_addr_button, False, False, 5)
    add_addr_button.show()

    sep_label = gtk.Label('          ')
    comp_std_button_box.pack_start(sep_label, False, False, 0)
    sep_label.show()

    if (del_comp_std == True):
      del_button = gtk.Button('Delete last component standardiser')
      del_button.connect('clicked', self.clickDelCompStdButton,
                         comp_std_button_box)
      comp_std_button_box.pack_start(del_button, False, False, 0)
      del_button.show()

    return comp_std_button_box

  # ---------------------------------------------------------------------------
  #
  def clickAddStandardiserButton(self, widget, add_comp_std_box):
    print 'Clicked "Add new component standardiser"'

    comp_std_type = widget.get_label()[:4]  # The type of comp. std. to add

    self.comp_std_type_count_dict[comp_std_type] += 1

    # A new dictionary for this component standardiser
    #
    self.comp_std_list.append({'type':comp_std_type})

    std_scrolled_window = self.mainTree.get_widget('standard_scrolled_window')

    # Get the box with component standardiser (first child is the viewport)
    #
    std_box = std_scrolled_window.get_child().get_child()

    std_box.remove(add_comp_std_box)  # Remove box with buttons

    # Create a new component standardiser definiton with no values defined
    #
    comp_std_type_str = {'Date':'Date', 'Phon':'Phone number', 'Name':'Name',
                         'Addr':'Address'}[comp_std_type]
    comp_std_type_box = gtk.HBox()
    comp_std_type_label = gtk.Label('<b>%s standardiser:</b>' % \
                                    (comp_std_type_str))
    comp_std_type_label.set_use_markup(True)
    comp_std_type_box.pack_start(comp_std_type_label, False, False, 0)
    comp_std_type_label.show()
    std_box.pack_start(comp_std_type_box, False, False, 0)
    comp_std_type_box.show()

    comp_std_box = self.new_comp_std(comp_std_type)
    std_box.pack_start(comp_std_box, False, False, 0)
    comp_std_box.show()

    horiz_separator = gtk.HSeparator()  # Display a horizontal separator
    std_box.pack_start(horiz_separator, False, False, 5)
    horiz_separator.show()

    # Finally add 'Add' and 'Delete' buttons - - - - - - - - - - - - - - - -
    #
    num_comp_std = len(self.comp_std_list)

    if (num_comp_std == 0): # No component standardiser so far, only show 'Add'
      del_comp_std = False
    else:
      del_comp_std = True

    comp_std_button_box = self.create_comp_std_button(del_comp_std)

    std_box.pack_start(comp_std_button_box, False, False, 5)
    comp_std_button_box.show()
    std_box.show()

  # ---------------------------------------------------------------------------
  #
  def clickDelCompStdButton(self, widget, add_comp_std_box):
    print 'Clicked "Delete last component standardiser"'

    std_scrolled_window = self.mainTree.get_widget('standard_scrolled_window')

    last_comp_std_type = self.comp_std_list[-1]['type']
    self.comp_std_type_count_dict[last_comp_std_type] -= 1

    self.comp_std_list.pop()  # Remove last component standardiser dictionary

    # Get the box with component standardiser definitions
    #
    std_box = std_scrolled_window.get_child().get_child()

    std_box.remove(add_comp_std_box)  # Remove box with buttons
    std_box.remove(std_box.get_children()[-1])  # Remove separator
    std_box.remove(std_box.get_children()[-1])  # Remove last component std.
    std_box.remove(std_box.get_children()[-1])  # Remove its type

    # Check if only one component standardiser left - - - - - - - - - - - - - -
    #
    num_comp_std = len(self.comp_std_list)

    if (num_comp_std == 0): # No component standardiser so far, only show 'Add'
      del_comp_std = False
    else:
      del_comp_std = True

    # Finally add 'Add' and 'Delete' buttons - - - - - - - - - - - - - - - -
    #
    comp_std_button_box = self.create_comp_std_button(del_comp_std)

    std_box.pack_start(comp_std_button_box, False, False, 0)
    comp_std_button_box.show()
    std_box.show()

  # ---------------------------------------------------------------------------
  # Construct a new component standardiser of given type and possibly with the
  # values given in the dictionary.
  # Returns a HBox containing the new component standardiser.
  #
  def new_comp_std(self, comp_std_type, comp_std_dict={}):

    in_field_box = gtk.VBox()   # Left vertical box for input fields
    parameter_box = gtk.VBox()  # Middle vertical box for input fields
    out_field_box = gtk.VBox()  # Right vertical box for input fields

    comp_std_type_str = {'Date':'Date', 'Phon':'Phone number', 'Name':'Name',
                         'Addr':'Address'}[comp_std_type]

    # Left box with input fields (same for all comp. std. types) - - - - - - -
    #
    in_field_label = gtk.Label('Input fields:')
    in_field_box.pack_start(in_field_label, False, False, 5)
    in_field_label.show()

    # Get the list of already defined input fields for this comp. standardiser
    #
    in_field_list = comp_std_dict.get('in_field_list', [])

    # Get the field names from the input data set
    #
    ds_in_field_names = self.data_set_info_list[0]['field_names']

    # Generate and add a combo box for each defined input field
    #
    for in_field in in_field_list:

      field_combo_box = gtk.combo_box_new_text()
      for fn in ds_in_field_names:
        field_combo_box.append_text(fn)
      field_name_index = ds_in_field_names.index(in_field)
      field_combo_box.set_active(field_name_index)
      in_field_box.pack_start(field_combo_box, False, False, 0)
      field_combo_box.show()

    # If no input field defined add one more, empty input field combo box
    #
    if (in_field_list == []):
      field_combo_box = gtk.combo_box_new_text()
      for fn in ds_in_field_names:
        field_combo_box.append_text(fn)
      field_combo_box.set_active(-1)  # Set to no field name
      in_field_box.pack_start(field_combo_box, False, False, 0)
      field_combo_box.show()

    horiz_separator = gtk.HSeparator()  # Display a horizontal separator
    in_field_box.pack_start(horiz_separator, False, False, 5)
    horiz_separator.show()

    # Add button for new (and possibly delete last) input field combo box
    #
    if (len(in_field_list) > 1):
      del_in_field = True
    else:
      del_in_field = False
    in_field_button_box = self.create_in_field_combo_box(del_in_field)
    in_field_box.pack_start(in_field_button_box, False, False, 0)
    in_field_button_box.show()

    # Middle box with parameters - - - - - - - - - - - - - - - - - - - - - - -
    #
    parameter_label = gtk.Label('Parameters:')
    parameter_box.pack_start(parameter_label, False, False, 5)
    parameter_label.show()

    # Build a table, size depending upon comp. standardiser type
    #
    num_table_rows = {'Date':2, 'Phon':2, 'Name':10, 'Addr':7}[comp_std_type]
    parameter_table = gtk.Table(num_table_rows, 2)  # , True)

    # Add the field separator parameter first (not for date component)
    #
    if (comp_std_type != 'Date'):
      field_sep_label = gtk.Label('Field separator: ')
      field_sep_label.set_alignment(xalign=1.0, yalign=0.5)
      field_sep_label.show()
      parameter_table.attach(field_sep_label, 0, 1, 0, 1)
      field_sep_entry = gtk.Entry()
      field_sep_val = comp_std_dict.get('field_separator', '')
      field_sep_entry.set_text(field_sep_val)
      field_sep_entry.show()
      parameter_table.attach(field_sep_entry, 1, 2, 0, 1)

    # Add check word spill only for name and address component
    #
    if (comp_std_type in ['Name','Addr']):
      check_spill_button = gtk.CheckButton('Check word spilling')
      check_spill_val = comp_std_dict.get('check_word_spill', False)
      check_spill_button.set_active(check_spill_val)
      parameter_table.attach(check_spill_button, 1, 2, 1, 2)
      check_spill_button.show()

    if (comp_std_type == 'Date'):  # Date specific parameters

      default_parse_formats = self.parse_format_str

      parse_format_label = gtk.Label('Parse formats: ')
      parse_format_label.set_alignment(xalign=1.0, yalign=0.5)
      parse_format_label.show()
      parameter_table.attach(parse_format_label, 0, 1, 0, 1)
      parse_formats_entry = gtk.Entry()
      # parse_formats_entry.set_width_chars(40)
      parse_formats_val = comp_std_dict.get('parse_formats',
                                            default_parse_formats)
      parse_formats_entry.set_text(parse_formats_val)
      parse_formats_entry.show()
      parameter_table.attach(parse_formats_entry, 1, 2, 0, 1)

      pivot_year_label = gtk.Label('Pivot year: ')
      pivot_year_label.set_alignment(xalign=1.0, yalign=0.5)
      pivot_year_label.show()
      parameter_table.attach(pivot_year_label, 0, 1, 1, 2)
      pivot_year_entry = gtk.Entry()
      # pivot_year_entry.set_width_chars(40)
      pivot_year_val = comp_std_dict.get('pivot_year', '')
      if (pivot_year_val == ''):  # Make default current year + 1
        pivot_year_val = str(time.localtime()[0]+1)[2:]
      pivot_year_entry.set_text(pivot_year_val)
      pivot_year_entry.show()
      parameter_table.attach(pivot_year_entry, 1, 2, 1, 2)

    elif (comp_std_type == 'Phon'):  # Phone number specific parameters

      def_country_label = gtk.Label('Default country: ')
      def_country_label.set_alignment(xalign=1.0, yalign=0.5)
      def_country_label.show()
      parameter_table.attach(def_country_label, 0, 1, 1, 2)
      def_country_combo_box = gtk.combo_box_new_text()  # New combo box
      def_country_combo_box.append_text('Australia')
      def_country_combo_box.append_text('Canada/USA')
      def_country_val = comp_std_dict.get('default_country', 'australia')
      def_country_index = {'australia':0, 'canada/usa':1}[def_country_val]
      def_country_combo_box.set_active(def_country_index)
      def_country_combo_box.show()
      parameter_table.attach(def_country_combo_box, 1, 2, 1, 2)

    else:  # Correction list and tag table used with name and address stand.

      corr_list_label = gtk.Label('Correction list file: ')
      corr_list_label.set_alignment(xalign=1.0, yalign=0.5)
      corr_list_label.show()
      parameter_table.attach(corr_list_label, 0, 1, 2, 3)
      corr_list_file_name = comp_std_dict.get('corr_list', '(None)')
      corr_list_button = gtk.Button()
      corr_list_button.set_name(corr_list_file_name)
      corr_list_button.set_label(corr_list_file_name.split(os.sep)[-1])
      corr_list_button.connect('clicked', self.runClickedStandardFileButton)
      corr_list_button.show()
      parameter_table.attach(corr_list_button, 1, 2, 2, 3)

      tag_table_label = gtk.Label('Tag table file(s): ')
      tag_table_label.set_alignment(xalign=1.0, yalign=0.0)
      tag_table_label.show()
      parameter_table.attach(tag_table_label, 0, 1, 3, 4)

      tag_table_box = gtk.VBox()
      tag_table_file_list = comp_std_dict.get('tag_table', [])
      for tag_table_file_name in tag_table_file_list:
        tag_table_button = gtk.Button()
        tag_table_button.set_name(tag_table_file_name)
        tag_table_button.set_label(tag_table_file_name.split(os.sep)[-1])
        tag_table_button.connect('clicked', self.runClickedStandardFileButton)
        tag_table_box.pack_start(tag_table_button, False, False, 0)
        tag_table_button.show()

      if (tag_table_file_list == []):  # Add new tag table button if none
        tag_table_button = gtk.Button()
        tag_table_button.set_label('(None)')
        tag_table_button.set_name('(None)')
        tag_table_button.connect('clicked', self.runClickedStandardFileButton)
        tag_table_box.pack_start(tag_table_button, False, False, 0)
        tag_table_button.show()

      # Generate buttons to add/delete last tag table
      #
      tag_table_button_box = gtk.HBox()
      add_tag_table_button = gtk.Button('Add tag table')
      add_tag_table_button.connect('clicked', self.clickAddTagTableButton)
      tag_table_button_box.pack_start(add_tag_table_button, False, False, 5)
      add_tag_table_button.show()
      if (len(tag_table_file_list) > 1):  # Delete only if more than one
        del_tag_table_button = gtk.Button('Delete last tag table')
        del_tag_table_button.connect('clicked', self.clickDelTagTableButton)
        tag_table_button_box.pack_start(del_tag_table_button, False, False, 5)
        del_tag_table_button.show()

      tag_table_box.pack_start(tag_table_button_box, False, False, 0)
      tag_table_button_box.show()

      parameter_table.attach(tag_table_box, 1, 2, 3, 4)
      tag_table_box.show()


      hmm_file_label = gtk.Label('%s HMM file: ' % (comp_std_type_str))
      hmm_file_label.set_alignment(xalign=1.0, yalign=0.5)
      hmm_file_label.show()
      parameter_table.attach(hmm_file_label, 0, 1, 4, 5)
      hmm_file_name = comp_std_dict.get('%s_hmm' % (comp_std_type_str.lower()),
                                        '(None)')
      hmm_file_button = gtk.Button()
      hmm_file_button.set_name(hmm_file_name)
      hmm_file_button.set_label(hmm_file_name.split(os.sep)[-1])
      hmm_file_button.connect('clicked', self.runClickedStandardFileButton)
      hmm_file_button.show()
      parameter_table.attach(hmm_file_button, 1, 2, 4, 5)

      hmm_train_label = gtk.Label('HMM training file: ')
      hmm_train_label.set_alignment(xalign=1.0, yalign=0.5)
      hmm_train_label.show()
      parameter_table.attach(hmm_train_label, 0, 1, 5, 6)
      hmm_train_name = comp_std_dict.get('hmm_train_file', '(None)')
      hmm_train_button = gtk.Button()
      hmm_train_button.set_name(hmm_train_name)
      hmm_train_button.set_label(hmm_train_name.split(os.sep)[-1])
      hmm_train_button.connect('clicked', self.runClickedStandardFileButton)
      hmm_train_button.show()
      parameter_table.attach(hmm_train_button, 1, 2, 5, 6)

      hmm_prob_label = gtk.Label('HMM seq. probability file: ')
      hmm_prob_label.set_alignment(xalign=1.0, yalign=0.5)
      hmm_prob_label.show()
      parameter_table.attach(hmm_prob_label, 0, 1, 6, 7)
      hmm_prob_name = comp_std_dict.get('hmm_seq_prob_file', '(None)')
      hmm_prob_button = gtk.Button()
      hmm_prob_button.set_name(hmm_prob_name)
      hmm_prob_button.set_label(hmm_prob_name.split(os.sep)[-1])
      hmm_prob_button.connect('clicked', self.runClickedStandardFileButton)
      hmm_prob_button.show()
      parameter_table.attach(hmm_prob_button, 1, 2, 6, 7)

    if (comp_std_type == 'Name'):  # Name specific parameters

      female_titles_label = gtk.Label('Female titles: ')
      female_titles_label.set_alignment(xalign=1.0, yalign=0.5)
      female_titles_label.show()
      parameter_table.attach(female_titles_label, 0, 1, 7, 8)
      female_titles_entry = gtk.Entry()
      female_titles_val = comp_std_dict.get('female_titles', 'ms, mrs')
      female_titles_entry.set_text(female_titles_val)
      female_titles_entry.show()
      parameter_table.attach(female_titles_entry, 1, 2, 7, 8)

      male_titles_label = gtk.Label('Male titles: ')
      male_titles_label.set_alignment(xalign=1.0, yalign=0.5)
      male_titles_label.show()
      parameter_table.attach(male_titles_label, 0, 1, 8, 9)
      male_titles_entry = gtk.Entry()
      male_titles_val = comp_std_dict.get('male_titles', 'mr')
      male_titles_entry.set_text(male_titles_val)
      male_titles_entry.show()
      parameter_table.attach(male_titles_entry, 1, 2, 8, 9)

      first_name_comp_label = gtk.Label('First name component: ')
      first_name_comp_label.set_alignment(xalign=1.0, yalign=0.5)
      first_name_comp_label.show()
      parameter_table.attach(first_name_comp_label, 0, 1, 9, 10)
      first_name_comp_combo_box = gtk.combo_box_new_text()  # New combo box
      first_name_comp_combo_box.append_text('Given name')
      first_name_comp_combo_box.append_text('Surname')
      first_name_comp_val = comp_std_dict.get('first_name_comp', 'gname')
      first_name_comp_index = {'gname':0, 'sname':1}[first_name_comp_val]
      first_name_comp_combo_box.set_active(first_name_comp_index)
      first_name_comp_combo_box.show()
      parameter_table.attach(first_name_comp_combo_box, 1, 2, 9, 10)

    parameter_box.pack_start(parameter_table, False, False, 0)
    parameter_table.show()

    # Right box with output fields - - - - - - - - - - - - - - - - - - - - - -
    #
    out_field_label = gtk.Label('Output fields:')
    out_field_box.pack_start(out_field_label, False, False, 5)
    out_field_label.show()

    # List of output field names for this component standardiser type
    #
    comp_std_out_fields_list = self.comp_std_out_fields[comp_std_type]
    num_out_fields = len(comp_std_out_fields_list)

    # List of the actual field names for this standardiser
    #
    out_field_list = comp_std_dict.get('out_field_list', [])

    # If no output fields given, generate the list
    #
    if (out_field_list == []):
      comp_std_num_str = str(self.comp_std_type_count_dict[comp_std_type])

      for out_field_name in comp_std_out_fields_list:
        out_field_list.append(out_field_name+comp_std_num_str)

    out_field_table = gtk.Table(num_out_fields, 2, False)

    row_cnt = 0
    for out_field_name in comp_std_out_fields_list:
      out_field_str = out_field_name.replace('_',' ').capitalize()+': '
      out_field_label = gtk.Label(out_field_str)
      out_field_label.set_alignment(xalign=1.0, yalign=0.5)
      out_field_label.show()
      out_field_table.attach(out_field_label, 0, 1, row_cnt, row_cnt+1)
      out_field_entry = gtk.Entry()
      out_field_entry.set_text(out_field_list[row_cnt])
      out_field_entry.show()
      out_field_table.attach(out_field_entry, 1, 2, row_cnt, row_cnt+1)

      row_cnt +=1

    out_field_box.pack_start(out_field_table, False, False, 0)
    out_field_table.show()

    # Finally put all boxes into a HBox, show and return - - - - - - - - - - -
    #
    comp_std_box = gtk.HBox()
    comp_std_box.pack_start(in_field_box, False, False, 0)
    in_field_box.show()

#    sep_label1 = gtk.Label('     ')
#    comp_std_box.pack_start(sep_label1, False, False, 0)
#    sep_label1.show()

    vert_separator = gtk.VSeparator()  # Display a vertical separator
    comp_std_box.pack_start(vert_separator, False, False, 20)
    vert_separator.show()

    comp_std_box.pack_start(parameter_box, False, False, 0)
    parameter_box.show()

    vert_separator = gtk.VSeparator()  # Display a vertical separator
    comp_std_box.pack_start(vert_separator, False, False, 20)
    vert_separator.show()

    comp_std_box.pack_start(out_field_box, False, False, 0)
    out_field_box.show()

    return comp_std_box

  # ---------------------------------------------------------------------------
  # Handle clicks on the add tag table button
  #
  def clickAddTagTableButton(self, widget):
    tag_table_box = widget.get_parent().get_parent()

    # Remove box with Add/Delete button
    #
    tag_table_box.remove(tag_table_box.get_children()[-1])

    tag_table_button = gtk.Button()  # Add new tag table button
    tag_table_button.set_label('(None)')
    tag_table_button.connect('clicked', self.runClickedStandardFileButton)
    tag_table_box.pack_start(tag_table_button, False, False, 0)
    tag_table_button.show()

    # Generate buttons to add/delete last tag table
    #
    tag_table_button_box = gtk.HBox()
    add_tag_table_button = gtk.Button('Add tag table')
    add_tag_table_button.connect('clicked', self.clickAddTagTableButton)
    tag_table_button_box.pack_start(add_tag_table_button, False, False, 5)
    add_tag_table_button.show()

    if (len(tag_table_box.get_children()) > 1):  # Delete only if more than one
      del_tag_table_button = gtk.Button('Delete last tag table')
      del_tag_table_button.connect('clicked', self.clickDelTagTableButton)
      tag_table_button_box.pack_start(del_tag_table_button, False, False, 5)
      del_tag_table_button.show()

    tag_table_box.pack_start(tag_table_button_box, False, False, 0)
    tag_table_button_box.show()

  # ---------------------------------------------------------------------------
  # Handle clicks on the delete tag table button
  #
  def clickDelTagTableButton(self, widget):
    tag_table_box = widget.get_parent().get_parent()

    # Remove box with Add/Delete button and last tag table entry
    #
    tag_table_box.remove(tag_table_box.get_children()[-1])
    tag_table_box.remove(tag_table_box.get_children()[-1])

    # Generate buttons to add/delete last tag table
    #
    tag_table_button_box = gtk.HBox()
    add_tag_table_button = gtk.Button('Add tag table')
    add_tag_table_button.connect('clicked', self.clickAddTagTableButton)
    tag_table_button_box.pack_start(add_tag_table_button, False, False, 5)
    add_tag_table_button.show()

    if (len(tag_table_box.get_children()) > 1):  # Delete only if more than one
      del_tag_table_button = gtk.Button('Delete last tag table')
      del_tag_table_button.connect('clicked', self.clickDelTagTableButton)
      tag_table_button_box.pack_start(del_tag_table_button, False, False, 5)
      del_tag_table_button.show()

    tag_table_box.pack_start(tag_table_button_box, False, False, 0)
    tag_table_button_box.show()

  # ---------------------------------------------------------------------------
  # Handle clicks on the correction list button
  #
  def runClickedStandardFileButton(self, widget):
    widget_name = widget.get_name()
    print 'Clicked on a standardise file button: %s' % (widget_name)

    # Create file dialog
    #
    saveFileTree = gtk.glade.XML(self.glade_file_name,
                                 self.save_file_dialog_name)
    saveFileDialog = saveFileTree.get_widget(self.save_file_dialog_name)

    if (widget_name != '(None)'):
      saveFileDialog.set_current_name(widget_name)

    saveFileDialog.show()
    save_file_response = saveFileDialog.run()
    save_dialog_file_name = saveFileDialog.get_filename()
    saveFileDialog.destroy()

    # Only process if 'OK' was clicked, not 'Cancel'
    #
    if (save_dialog_file_name != None):

      # Special case - if None was typed in, set to (None)
      #
      if (save_dialog_file_name.split(os.sep)[-1].lower() in ['(none)',
                                                              'none']):
        widget.set_name('(None)')
        widget.set_label('(None)')
      else:
        widget.set_name(save_dialog_file_name)
        widget.set_label(save_dialog_file_name.split(os.sep)[-1])

  # ---------------------------------------------------------------------------
  # Construct a new HBox with buttons to add an input field and (if more than
  # one exists) to delete the ast one
  # Returns the HBox with the buttons.
  # If the 'del_in_field' flag is set to False then the corresponding button
  # will not be created
  #
  def create_in_field_combo_box(self, del_in_field=True):

    in_field_button_box = gtk.HBox()

    add_field_button = gtk.Button('Add new input field')
    add_field_button.connect('clicked', self.clickAddStdInFieldButton,
                            in_field_button_box)
    in_field_button_box.pack_start(add_field_button, False, False, 5)
    add_field_button.show()

    if (del_in_field == True):
      del_field_button = gtk.Button('Delete last input field')
      del_field_button.connect('clicked', self.clickDelStdInFieldButton,
                              in_field_button_box)
      in_field_button_box.pack_start(del_field_button, False, False, 5)
      del_field_button.show()

    return in_field_button_box

  # ---------------------------------------------------------------------------
  #
  def clickAddStdInFieldButton(self, widget, add_in_field_box):
    print 'Clicked "Add new comp. std. input field"'

    # The VBox with all input fields
    #
    in_field_box = widget.get_parent().get_parent()

    # Remove button box and horizontal separator first
    #
    in_field_box.remove(in_field_box.get_children()[-1])
    in_field_box.remove(in_field_box.get_children()[-1])

    # Get the field names from the input data set
    #
    ds_in_field_names = self.data_set_info_list[0]['field_names']

    # Add a new, empty input field combo box
    #
    field_combo_box = gtk.combo_box_new_text()
    for fn in ds_in_field_names:
      field_combo_box.append_text(fn)
    field_combo_box.set_active(-1)  # Set to no field name
    in_field_box.pack_start(field_combo_box, False, False, 0)
    field_combo_box.show()

    horiz_separator = gtk.HSeparator()  # Display a horizontal separator
    in_field_box.pack_start(horiz_separator, False, False, 0)
    horiz_separator.show()

    # Number of input fields is box children minusheader, separator and buttons
    #
    num_in_fields = len(in_field_box.get_children())-2

    # Add button for new (and possibly delete last) input field combo box
    #
    if (num_in_fields > 1):
      del_in_field = True
    else:
      del_in_field = False
    in_field_button_box = self.create_in_field_combo_box(del_in_field)
    in_field_box.pack_start(in_field_button_box, False, False, 0)
    in_field_button_box.show()

  # ---------------------------------------------------------------------------
  #
  def clickDelStdInFieldButton(self, widget, add_in_field_box):
    print 'Clicked "Delete last comp. std. input field"'

    # The VBox with all input fields
    #
    in_field_box = widget.get_parent().get_parent()

    # Remove button box, horizontal separator andlast input field
    #
    in_field_box.remove(in_field_box.get_children()[-1])
    in_field_box.remove(in_field_box.get_children()[-1])
    in_field_box.remove(in_field_box.get_children()[-1])

    # Now add new separator and button box
    horiz_separator = gtk.HSeparator()  # Display a horizontal separator
    in_field_box.pack_start(horiz_separator, False, False, 0)
    horiz_separator.show()

    # Number of input fields is box children minusheader, separator and buttons
    #
    num_in_fields = len(in_field_box.get_children())-2

    # Add button for new (and possibly delete last) input field combo box
    #
    if (num_in_fields > 1):
      del_in_field = True
    else:
      del_in_field = False
    in_field_button_box = self.create_in_field_combo_box(del_in_field)
    in_field_box.pack_start(in_field_button_box, False, False, 0)
    in_field_button_box.show()

  # ---------------------------------------------------------------------------
  # Handle an activate of Execute on the Standardise page
  #
  def standardExecute(self):
    self.addToLog('Execute on Standardise page')

    self.febrl_code['standardise'] = None  # Remove all previous std. code

    std_scrolled_window = self.mainTree.get_widget('standard_scrolled_window')
    std_box = std_scrolled_window.get_child().get_child()
    std_box_list = std_box.get_children()

    comp_std_list = self.comp_std_list  # Short cut

    if (comp_std_list == []):  # No component standardiser defined
      self.messageDialog('At least one component standardiser is required.',
                         'error')
      return

    for i in range(len(comp_std_list)):
      this_comp_std_dict = {}  # The details for this component standardiser

      this_comp_std_type = comp_std_list[i]['type']
      this_comp_std_dict['type'] = this_comp_std_type

      type_box_index = i*3  # Hbox with type information
      comp_std_type_label = std_box_list[type_box_index].get_children()[0]

      # Make sure type in list is same as type from GUI
      #
      assert comp_std_type_label.get_text()[:4] == this_comp_std_type

      comp_std_box = std_box_list[type_box_index+1]

      in_field_box =  comp_std_box.get_children()[0]
      parameter_box = comp_std_box.get_children()[2]  # Vertical separators
      out_field_box = comp_std_box.get_children()[4]  # between

      # Extract the input field names and check they are not None - - - - - - -
      #
      in_field_list = []

      # First child is the label, 2nd last the separator, last the 'Add' button
      #
      for in_field_combo_box in in_field_box.get_children()[1:-2]:
        in_field_val = in_field_combo_box.get_active_text()

        if (in_field_val == None):
          self.messageDialog('All input fields must be selected.', 'error')
          return
        in_field_list.append(in_field_val)

      this_comp_std_dict['in_field_list'] = in_field_list
      print 'in field list:', in_field_list

      # Extract the parameters and check their values - - - - - - - - - - - - -
      #
      table_cell_widgets = parameter_box.get_children()[1].get_children()

      # For some reason originally starting from lower right corner
      #
      table_cell_widgets.reverse()

      # Field separator (not for date standardiser)
      #
      if (this_comp_std_type != 'Date'):
        field_sep_val = table_cell_widgets[1].get_text()
        this_comp_std_dict['field_separator'] = field_sep_val

      # Check word spill for name and address only
      #
      if (this_comp_std_type in ['Name','Addr']):
        check_word_spill_val = table_cell_widgets[2].get_active()
        this_comp_std_dict['check_word_spill'] = check_word_spill_val

        # Word spill only possible with no (empty) field separator
        #
        if ((check_word_spill_val == True) and (field_sep_val == '')):
          self.messageDialog('Word spill checking only possible with\nan ' + \
                                 'non-empty field separator.', 'error')
          return

      if (this_comp_std_type == 'Date'):  # Date specific parameters

        parse_formats_val = table_cell_widgets[1].get_text()
        for format_str in parse_formats_val.split(','):  # Check parse formats
          if (len(format_str.strip()) != 8):
            self.messageDialog('Date parse format string has wrong length' + \
                               '\n(should be 8 characters long).', 'error')
            return
          tmp_str = format_str.strip()

          while (tmp_str != ''):  # Check if all directives are valid
            if (tmp_str[:2] not in ['%b','%B','%d','%m','%y','%Y']):
              self.messageDialog('Illegal directive in format string: ' + \
                                 '"%s".' % (format_str), 'error')
              return
            tmp_str = tmp_str[3:]  # Remove directive and following space

        this_comp_std_dict['parse_formats'] = parse_formats_val

        pivot_year_val = table_cell_widgets[3].get_text()
        if ((pivot_year_val.isdigit() != True) or (int(pivot_year_val) < 0) \
            or (int(pivot_year_val) > 99)):
          self.messageDialog('Pivot year must be an integer between 0 and 99.',
                             'error')
          return
        this_comp_std_dict['pivot_year'] = str(int(pivot_year_val))

      elif (this_comp_std_type == 'Phon'):  # Phone number specific parameters

        def_country_val = table_cell_widgets[3].get_active_text().lower()
        this_comp_std_dict['default_country'] = def_country_val

      else:  # Correction list and tag table used name and address stand.

        corr_list_val = table_cell_widgets[4].get_name()
        if (corr_list_val.split(os.sep)[-1].lower() in ['(none)', 'none']):
          corr_list_val = '(None)'
        this_comp_std_dict['corr_list'] = corr_list_val

        tag_table_file_list = []
        for tag_table_button in table_cell_widgets[6].get_children()[:-1]:
          tag_table_val = tag_table_button.get_name()
          if (tag_table_val.split(os.sep)[-1].lower() in ['(none)', 'none']):
            tag_table_val = '(None)'

          if (tag_table_val not in ['(None)','GtkButton']):
            tag_table_file_list.append(tag_table_val)  # Don't append empty

        this_comp_std_dict['tag_table'] = tag_table_file_list

        # Check if either HMM file or training file are not None
        #
        hmm_file_val =       table_cell_widgets[8].get_name()
        hmm_train_file_val = table_cell_widgets[10].get_name()

        if ((hmm_file_val == '(None)') and (hmm_train_file_val == '(None)')):
          self.messageDialog('Both HMM file and HMM training file set to ' + \
                             'None.\nAt least one must be initialised.',
                             'error')
          return

        hmm_file_val = table_cell_widgets[8].get_name()

        if (this_comp_std_type == 'Addr'):
          dict_key_name = 'address_hmm'
        else:
          dict_key_name = 'name_hmm'

        this_comp_std_dict[dict_key_name] = hmm_file_val
        hmm_train_file_val = table_cell_widgets[10].get_name()
        this_comp_std_dict['hmm_train_file'] = hmm_train_file_val

        hmm_seq_prob_file_val = table_cell_widgets[12].get_name()
        this_comp_std_dict['hmm_seq_prob_file'] = hmm_seq_prob_file_val

      if (this_comp_std_type == 'Name'):  # Name specific parameters

        female_titles_val = table_cell_widgets[14].get_text()
        this_comp_std_dict['female_titles'] = female_titles_val

        male_titles_val = table_cell_widgets[16].get_text()
        this_comp_std_dict['male_titles'] = male_titles_val

        first_name_comp_val = table_cell_widgets[18].get_active_text()
        if (first_name_comp_val[0] == 'G'):
          this_comp_std_dict['first_name_comp'] = 'gname'
        else:
          this_comp_std_dict['first_name_comp'] = 'sname'

      # Extract the output field names - - - - - - - - - - - - - - - - - - - -
      #
      out_field_list = []
      num_out_field_none = 0

      out_field_table = out_field_box.get_children()[1]

      table_cell_widgets = out_field_table.get_children()

      # For some reason originally starting from lower right corner
      #
      table_cell_widgets.reverse()

      j = 0
      for widget in table_cell_widgets:
        if (j % 2 == 1):  # Only for text entries (on right side)
          out_field_val = widget.get_text().strip()
          if (out_field_val.lower() in ['', '(none)', 'none']):
            out_field_val = 'None'
            num_out_field_none += 1
          out_field_list.append(out_field_val)
        j += 1

      # Get the number of expected output fields
      #
      num_out_field = len(self.comp_std_out_fields[this_comp_std_type])
      assert len(out_field_list) == num_out_field

      if (num_out_field_none == num_out_field):  # No output field defined
        self.messageDialog('At least one output field must be not None.',
                           'error')
        return

      this_comp_std_dict['out_field_list'] = out_field_list
      print 'out field list:', out_field_list

      comp_std_list[i] = this_comp_std_dict

    print 'comp_std_list:', comp_std_list

    # Now generate the Febrl codes for standardisation - - - - - - - - - - - -
    #
    febrl_code = []
    febrl_code.append('# '+'-'*77)
    febrl_code.append('')
    febrl_code.append('# Define component standardisers')
    febrl_code.append('#')

    corr_list_cnt = 1  # Counters for correction list and tag table definitions
    tag_table_cnt = 1
    comp_std_cnt =  1

    # Loop over all component standardisers - - - - - - - - - - - - - - - - - -
    #
    for comp_std_dict in comp_std_list:
      comp_std_type =  comp_std_dict['type']
      in_field_list =  comp_std_dict['in_field_list']
      out_field_list = comp_std_dict['out_field_list']

      comp_std_type_str = {'Date':'Date', 'Phon':'Phone number', 'Name':'Name',
                           'Addr':'Address'}[comp_std_type]

      corr_list_file_name = comp_std_dict.get('corr_list', '(None)')

      if (corr_list_file_name != '(None)'): # Need to initialise a corr. list
        corr_list_name = '%s_corr_list_%d' % (comp_std_type.lower(),
                                              corr_list_cnt)
        febrl_code.append('%s = lookup.CorrectionList(' % (corr_list_name) + \
                          'descr = "%s correction list")' % \
                          (comp_std_type_str))
        febrl_code.append('%s.load("%s")' % (corr_list_name,
                                             corr_list_file_name))
        febrl_code.append('')

        corr_list_cnt += 1

      tag_table_file_list = comp_std_dict.get('tag_table', [])

      if (tag_table_file_list != []):  # Need to initialise a tag table
        tag_table_name = '%s_tag_table_%d' % (comp_std_type.lower(),
                                              tag_table_cnt)
        febrl_code.append('%s = lookup.TagLookupTable(' % (tag_table_name) + \
                          'descr = "%s tag lookup table")' % \
                          (comp_std_type_str))
        indent_space = ' '*22
        febrl_code.append('%s.load(["%s",' % (tag_table_name,
                                              tag_table_file_list[0]))
        for tag_file_name in tag_table_file_list[1:]:
          febrl_code.append(indent_space+'"%s",' % (tag_file_name))

        # Add closing bracked to the last generated line
        #
        last_line = febrl_code[-1]
        last_line = last_line[:-1]+'])'  # Replace last comma
        febrl_code[-1] = last_line

        febrl_code.append('')

        tag_table_cnt += 1

      # Generate code for hidden Markov models - - - - - - - - - - - - - - - -
      #
      name_hmm_file_name = comp_std_dict.get('name_hmm', '(None)')
      if (name_hmm_file_name != '(None)'):
        febrl_code.append('name_hmm = simplehmm.hmm("Name HMM", [], [])')
        febrl_code.append('name_hmm.load_hmm("%s")' % (name_hmm_file_name))
        febrl_code.append('')

      addr_hmm_file_name = comp_std_dict.get('address_hmm', '(None)')
      if (addr_hmm_file_name != '(None)'):
        febrl_code.append('addr_hmm = simplehmm.hmm("Address HMM", [], [])')
        febrl_code.append('addr_hmm.load_hmm("%s")' % (addr_hmm_file_name))
        febrl_code.append('')

      # Generate the common attributes - - - - - - - - - - - - - - - - - - - -
      #
      in_field_str = '["' + '","'.join(in_field_list) + '"]'
      out_field_str = '['
      for out_field in out_field_list:
        if (out_field == 'None'):
          out_field_str += 'None,'
        else:
          out_field_str += '"%s",' % (out_field)
      out_field_str = out_field_str[:-1]+']'

      # Start generating the component standardiser - - - - - - - - - - - - - -
      #
      if (comp_std_type == 'Date'):  # Date standardiser
        indent_space = ' '*51
        febrl_code.append('date_comp_std_%d = standardisation.' % \
                          (comp_std_cnt) + \
                          'DateStandardiser(input_fields = %s,' % \
                          (in_field_str))

      elif (comp_std_type == 'Phon'):  # Phone number standardiser
        indent_space = ' '*56
        febrl_code.append('phone_comp_std_%d = standardisation.' % \
                          (comp_std_cnt) + \
                          'PhoneNumStandardiser(input_fields = %s,' % \
                          (in_field_str))

      elif (comp_std_type == 'Addr'):  # Address standardiser
        indent_space = ' '*57
        febrl_code.append('address_comp_std_%d = standardisation.' % \
                          (comp_std_cnt) + \
                          'AddressStandardiser(input_fields = %s,' % \
                          (in_field_str))

      else:  # Name standardiser
        indent_space = ' '*51
        febrl_code.append('name_comp_std_%d = standardisation.' % \
                          (comp_std_cnt) + \
                          'NameStandardiser(input_fields = %s,' % \
                          (in_field_str))

      # Some common attributes next
      #
      febrl_code.append(indent_space + 'output_fields = %s,' % \
                        (out_field_str))
      if (comp_std_type != 'Date'):
        febrl_code.append(indent_space + 'field_separator = "%s",' % \
                          (comp_std_dict['field_separator']))
      if (comp_std_type in ['Name','Addr']):
        febrl_code.append(indent_space + 'check_word_spill = %s,' % \
                          (str(comp_std_dict['check_word_spill'])))

      # Some standardiser specific attributes
      #
      if (comp_std_type == 'Date'):
        febrl_code.append(indent_space + 'pivot_year = %s,' % \
                          (comp_std_dict['pivot_year']))
        parse_format_ident_space = indent_space+' '*17

        first_parse_format = True
        for parse_format in comp_std_dict['parse_formats'].split(','):
          if (first_parse_format == True):
            febrl_code.append(indent_space + 'parse_formats = ["%s",' % \
                              (parse_format))
            first_parse_format = False
          else:
            febrl_code.append(indent_space + ' '*17 + '"%s",' % \
                              (parse_format.strip()))
        # Close parse format list and date standardiser
        #
        febrl_code[-1] = febrl_code[-1][:-1]+'])'

      else:  # Correction list and tag tables for all other standardisers
        if (corr_list_file_name != '(None)'):
          febrl_code.append(indent_space+'corr_list = %s,' % (corr_list_name))
        if (tag_table_file_list != []):
          febrl_code.append(indent_space+'tag_table = %s,' % (tag_table_name))

      # Some more standardiser specific attributes
      #
      if (comp_std_type == 'Phon'):  # Phone number standardiser
        febrl_code.append(indent_space + 'default_country = "%s")' % \
                          (str(comp_std_dict['default_country'])))

      elif (comp_std_type == 'Addr'):
        addr_hmm_file_name = comp_std_dict.get('address_hmm', '(None)')
        if (addr_hmm_file_name != '(None)'):
          febrl_code.append(indent_space + 'address_hmm = addr_hmm,')

      elif (comp_std_type == 'Name'):
        female_title_list = comp_std_dict['female_titles'].split(',')
        female_title_str = '['
        for title in female_title_list:
          if (title.strip() != ''):
            female_title_str += '"%s", ' % (title.strip())
        female_title_str = female_title_str[:-2] + ']'
        if (female_title_str == ']'):  # Special case an empty list
          female_title_str = '[]'

        male_title_list = comp_std_dict['male_titles'].split(',')
        male_title_str = '['
        for title in male_title_list:
          if (title.strip() != ''):
            male_title_str += '"%s", ' % (title.strip())
        male_title_str = male_title_str[:-2] + ']'
        if (male_title_str == ']'):  # Special case an empty list
          male_title_str = '[]'

        febrl_code.append(indent_space + 'female_titles = %s,' % \
                          (female_title_str))
        febrl_code.append(indent_space + 'male_titles = %s,' % \
                          (male_title_str))

        febrl_code.append(indent_space + 'first_name_comp = "%s",' % \
                          (comp_std_dict['first_name_comp']))

        name_hmm_file_name = comp_std_dict.get('name_hmm', '(None)')
        if (name_hmm_file_name != '(None)'):
          febrl_code.append(indent_space + 'name_hmm = name_hmm,')

      if (comp_std_type in ['Addr','Name']):
        hmm_train_file_name = comp_std_dict.get('hmm_train_file', '(None')
        if (hmm_train_file_name != '(None)'):
          febrl_code.append(indent_space + 'hmm_train_file = "%s",' % \
                            (hmm_train_file_name))
        hmm_seq_prob_file_name = comp_std_dict.get('hmm_seq_prob_file',
                                                   '(None')
        if (hmm_seq_prob_file_name != '(None)'):
          febrl_code.append(indent_space + 'hmm_seq_prob_file = "%s",' % \
                            (hmm_seq_prob_file_name))

        if (febrl_code[-1][-1] == ','):  # Replace comma with bracket
          last_line = febrl_code[-1]
          last_line = last_line[:-1]+')'  # Replace last comma
          febrl_code[-1] = last_line

      febrl_code.append('')

      comp_std_cnt += 1

    self.febrl_code['standardise'] = febrl_code  # Store for later use

    # Finally update the GUI information - - - - - - - - - - - - - - - - - - -
    #
    self.addToLog('')  # Add generated code into log page text
    self.addToLog('='*79)
    self.addToLog('Generated Febrl code for "standardise" on %s' % \
                  (time.asctime()))
    self.addToLog('')
    for line in febrl_code:
      self.addToLog(line)
    self.addToLog('')

    self.modified['standardise'] = True  # Standardisation details have changed
    self.setWindowTitle()

    # Update the active and non-active notebook pages
    #
    self.main_notebook_page_active_dict['Run'] = True

    self.displayCurrentNotebookPage()  # Diplay the current page

    self.writeStatusBar('Generated Febrl Python code for standardisation ' + \
                        '(see Log page for generated code).')


  # ===========================================================================
  # Methods that handle Index page events
  # ===========================================================================

  # ---------------------------------------------------------------------------
  # Display the Index page
  #
  def indexView(self):
    print '  Switched to Index page - Display this page'

    # Set alignment for separator string label to the right
    #
    self.mainTree.get_widget('index_sep_string_label').set_alignment(1, 0.5)

    indexing_widget = self.mainTree.get_widget('index_page_box')

    # Get the different widgets of the indexing page
    #
    index_widget_children = indexing_widget.get_children()
    index_method_box2 =     index_widget_children[1]  # Is a HBox
    index_method_box3 =     index_widget_children[2]  # Is a HBox
    indexing_separator =    index_widget_children[3]
    index_scrolled_window = index_widget_children[4]

    for child in index_method_box2:   # Remove all old children
      index_method_box2.remove(child)
    for child in index_method_box3:
      index_method_box3.remove(child)

    index_dict = self.index_method  # For quicker access
    print 'index dict:', index_dict

    # Check if an index method has previously been selected
    #
    if (index_dict['name'] == None):
      index_method_box2.hide()  # Hide not used page widgets
      index_method_box3.hide()
      indexing_separator.hide()
      index_scrolled_window.hide()
      return

    # Get the index method name and set general parameters in GUI - - - - - - -
    #
    index_method_name = index_dict['name']

    index_sep_str_widget =self.mainTree.get_widget('index_sep_string_entry')
    index_sep_str_widget.set_text(index_dict['index_sep_str'])

    skip_missing_widget = \
                    self.mainTree.get_widget('index_skip_miss_check_button')
    skip_missing_widget.set_active(index_dict['skip_missing'])

    # Display parameters depending upon index method selected - - - - - - - - -
    #
    if (index_method_name in ['FullIndex', 'BlockingIndex']):
      pass  # No parameters required

    elif (index_method_name == 'SortingIndex'):  # - - - - - - - - - - - - - -

      window_label = gtk.Label('     Window size:')  # Make it indented
      window_text_entry = gtk.Entry(2)
      window_text_entry.set_width_chars(5)
      if ('window_size') in index_dict:
        window_text_entry.set_text(index_dict['window_size'])

      index_method_box2.pack_start(window_label, False, False, 0)
      index_method_box2.pack_start(window_text_entry, False, False, 0)
      window_label.show()
      window_text_entry.show()
      index_method_box2.show()
      index_method_box3.hide()

    elif (index_method_name == 'QGramIndex'):  # - - - - - - - - - - - - - - -

      q_label = gtk.Label('     Length of Q:')  # Make it indented
      q_text_entry = gtk.Entry(2)
      q_text_entry.set_width_chars(5)
      if ('q') in index_dict:
        q_text_entry.set_text(index_dict['q'])

      threshold_label = gtk.Label('  Threshold:')
      threshold_text_entry = gtk.Entry(10)
      threshold_text_entry.set_width_chars(10)
      if ('threshold') in index_dict:
        threshold_text_entry.set_text(index_dict['threshold'])

      padded_check_button = gtk.CheckButton('Padded')
      if ('padded' in index_dict):
        padded_check_button.set_active(index_dict['padded'])
      else:  # Activate per default
        padded_check_button.set_active(True)

      index_method_box2.pack_start(q_label, False, False, 0)
      index_method_box2.pack_start(q_text_entry, False, False, 0)
      index_method_box2.pack_start(threshold_label, False, False, 4)
      index_method_box2.pack_start(threshold_text_entry, False, False, 0)
      index_method_box2.pack_start(padded_check_button, False, False, 4)
      q_label.show()
      q_text_entry.show()
      threshold_label.show()
      threshold_text_entry.show()
      padded_check_button.show()
      index_method_box2.show()

    elif (index_method_name == 'CanopyIndex'):  # - - - - - - - - - - - - - - -

      # First row of parameters is main canopy method and threshold method
      #
      canopy_method_label = gtk.Label('     Canopy method:')  # Indent
      index_method_box2.pack_start(canopy_method_label, False, False, 0)
      canopy_method_label.show()

      canopy_method_combo_box = gtk.combo_box_new_text()
      canopy_method_combo_box.append_text('TF-IDF')
      canopy_method_combo_box.append_text('Jaccard')

      if ('canopy_method' in index_dict):
        if (index_dict['canopy_method'][0] == 'tfidf'):
          canopy_method_combo_box.set_active(0)
        elif (index_dict['canopy_method'][0] == 'jaccard'):
          canopy_method_combo_box.set_active(1)
        else:
          canopy_method_combo_box.set_active(-1)
      else:
        canopy_method_combo_box.set_active(-1)

      index_method_box2.pack_start(canopy_method_combo_box, False, False, 0)
      canopy_method_combo_box.show()

      canopy_using_label = gtk.Label('  Thresholds:')
      index_method_box2.pack_start(canopy_using_label, False, False, 0)
      canopy_using_label.show()

      canopy_thresh_combo_box = gtk.combo_box_new_text()
      canopy_thresh_combo_box.append_text('Nearest')
      canopy_thresh_combo_box.append_text('Global')

      if ('canopy_method' in index_dict):
        if (index_dict['canopy_method'][1] == 'nearest'):
          canopy_thresh_combo_box.set_active(0)
        elif (index_dict['canopy_method'][1] == 'threshold'):
          canopy_thresh_combo_box.set_active(1)
        else:
          canopy_thresh_combo_box.set_active(-1)
      else:
        canopy_thresh_combo_box.set_active(-1)

      index_method_box2.pack_start(canopy_thresh_combo_box, False, False, 0)
      canopy_thresh_combo_box.show()

      canopy_parameters_label = gtk.Label('  Parameters:')
      index_method_box2.pack_start(canopy_parameters_label, False, False, 0)
      canopy_parameters_label.show()

      thresh_text_entry = gtk.Entry()
      thresh_text_entry.set_width_chars(20)
      if ('canopy_method' in index_dict):
        thresh_text_entry.set_text(index_dict['canopy_method'][2])
      index_method_box2.pack_start(thresh_text_entry, False, False, 0)
      thresh_text_entry.show()
      index_method_box2.show()

      # Second row is additional parameters
      #
      q_label = gtk.Label('     Length of Q:')  # Make it indented
      q_text_entry = gtk.Entry(2)
      q_text_entry.set_width_chars(5)
      if ('q') in index_dict:
        q_text_entry.set_text(index_dict['q'])

      del_perc_label = gtk.Label('Delete percentage:')
      del_perc_text_entry = gtk.Entry(3)
      del_perc_text_entry.set_width_chars(5)
      if ('delete_perc') in index_dict:
        del_perc_text_entry.set_text(index_dict['delete_perc'])
      else:
        del_perc_text_entry.set_text('100')

      padded_check_button = gtk.CheckButton('Padded')
      padded_check_button.set_active(True)
      if ('padded' in index_dict):
        padded_check_button.set_active(index_dict['padded'])
      else:  # Activate per default
        padded_check_button.set_active(True)

      index_method_box3.pack_start(q_label, False, False, 0)
      index_method_box3.pack_start(q_text_entry, False, False, 0)
      index_method_box3.pack_start(del_perc_label, False, False, 4)
      index_method_box3.pack_start(del_perc_text_entry, False, False, 0)
      index_method_box3.pack_start(padded_check_button, False, False, 4)
      q_label.show()
      q_text_entry.show()
      del_perc_label.show()
      del_perc_text_entry.show()
      padded_check_button.show()
      index_method_box3.show()

    elif (index_method_name == 'StringMapIndex'):  # - - - - - - - - - - - - -

      # First row of parameters is main canopy method and thresholds
      #
      strmap_method_label = gtk.Label('     Canopy thresholds:')  # Indent
      index_method_box2.pack_start(strmap_method_label, False, False, 0)
      strmap_method_label.show()

      strmap_thresh_combo_box = gtk.combo_box_new_text()
      strmap_thresh_combo_box.append_text('Nearest')
      strmap_thresh_combo_box.append_text('Global')

      if ('canopy_method' in index_dict):
        if (index_dict['canopy_method'][0] == 'nearest'):
          strmap_thresh_combo_box.set_active(0)
        elif (index_dict['canopy_method'][0] == 'threshold'):
          strmap_thresh_combo_box.set_active(1)
        else:
          strmap_thresh_combo_box.set_active(-1)
      else:
        strmap_thresh_combo_box.set_active(-1)

      index_method_box2.pack_start(strmap_thresh_combo_box, False, False, 0)
      strmap_thresh_combo_box.show()

      strmap_parameters_label = gtk.Label('  Parameters:')
      index_method_box2.pack_start(strmap_parameters_label, False, False, 0)
      strmap_parameters_label.show()

      thresh_text_entry = gtk.Entry()
      thresh_text_entry.set_width_chars(20)
      if ('canopy_method' in index_dict):
        thresh_text_entry.set_text(index_dict['canopy_method'][1])
      index_method_box2.pack_start(thresh_text_entry, False, False, 0)
      thresh_text_entry.show()

      strmap_grid_label = gtk.Label('  Grid resolution:')
      index_method_box2.pack_start(strmap_grid_label, False, False, 0)
      strmap_grid_label.show()

      grid_text_entry = gtk.Entry(5)
      grid_text_entry.set_width_chars(5)
      if ('grid_resolution' in index_dict):
        grid_text_entry.set_text(index_dict['grid_resolution'])
      index_method_box2.pack_start(grid_text_entry, False, False, 0)
      grid_text_entry.show()

      index_method_box2.show()

      # Second row is additional parameters
      #
      dim_label = gtk.Label('     Dimension:')  # Make it indented
      dim_text_entry = gtk.Entry(2)
      dim_text_entry.set_width_chars(5)
      if ('dim' in index_dict):
        dim_text_entry.set_text(index_dict['dim'])

      sub_dim_label = gtk.Label('Sub-dimension:')
      sub_dim_text_entry = gtk.Entry(2)
      sub_dim_text_entry.set_width_chars(5)
      if ('sub_dim' in index_dict):
        sub_dim_text_entry.set_text(index_dict['sub_dim'])

      sim_funct_label = gtk.Label('Similarity function:')
      sim_funct_combo_box = gtk.combo_box_new_text()
      sim_funct_name_list = self.stringcmp_dict.keys()
      sim_funct_name_list.sort()
      for f in sim_funct_name_list:
        sim_funct_combo_box.append_text(f)

      if ('sim_funct' in index_dict):
        sim_name_index = sim_funct_name_list.index(index_dict['sim_funct'])
        sim_funct_combo_box.set_active(sim_name_index)
      else:  # No name defined, set to None
        sim_funct_combo_box.set_active(-1)

      cache_check_button = gtk.CheckButton('Cache distance calculations')
      if ('cache_dist' in index_dict):
        cache_check_button.set_active(index_dict['cache_dist'])
      else:  # De-activate per default
        cache_check_button.set_active(False)

      index_method_box3.pack_start(dim_label, False, False, 0)
      index_method_box3.pack_start(dim_text_entry, False, False, 0)
      index_method_box3.pack_start(sub_dim_label, False, False, 4)
      index_method_box3.pack_start(sub_dim_text_entry, False, False, 0)
      index_method_box3.pack_start(sim_funct_label, False, False, 4)
      index_method_box3.pack_start(sim_funct_combo_box, False, False, 0)
      index_method_box3.pack_start(cache_check_button, False, False, 4)
      dim_label.show()
      dim_text_entry.show()
      sub_dim_label.show()
      sub_dim_text_entry.show()
      sim_funct_label.show()
      sim_funct_combo_box.show()
      cache_check_button.show()
      index_method_box3.show()

    elif (index_method_name == 'SuffixArrayIndex'):  # - - - - - - - - - - - -

      min_len_label = gtk.Label('     Minimum length:')  # Make it indented
      min_len_text_entry = gtk.Entry(2)
      min_len_text_entry.set_width_chars(5)
      if ('block_method' in index_dict):
        min_len_text_entry.set_text(index_dict['block_method'][0])

      max_block_label = gtk.Label('  Maximum block size:')
      max_block_text_entry = gtk.Entry(3)
      max_block_text_entry.set_width_chars(5)
      if ('block_method' in index_dict):
        max_block_text_entry.set_text(index_dict['block_method'][1])

      suffix_method_label = gtk.Label('  Suffix method:')
      suffix_method_combo_box = gtk.combo_box_new_text()
      suffix_method_combo_box.append_text('All-Substrings')
      suffix_method_combo_box.append_text('Suffix-Only')
      if ('suffix_method' in index_dict):
        if (index_dict['suffix_method'] == 'allsubstr'):
          suffix_method_combo_box.set_active(0)
        elif (index_dict['suffix_method'] == 'suffixonly'):
          suffix_method_combo_box.set_active(1)
        else:
          suffix_method_combo_box.set_active(-1)
      else:
        suffix_method_combo_box.set_active(-1)

      padded_check_button = gtk.CheckButton('Padded')
      if ('padded' in index_dict):
        padded_check_button.set_active(index_dict['padded'])
      else:  # Activate per default
        padded_check_button.set_active(True)

      index_method_box2.pack_start(min_len_label, False, False, 0)
      index_method_box2.pack_start(min_len_text_entry, False, False, 0)
      index_method_box2.pack_start(max_block_label, False, False, 0)
      index_method_box2.pack_start(max_block_text_entry, False, False, 0)
      index_method_box2.pack_start(suffix_method_label, False, False, 0)
      index_method_box2.pack_start(suffix_method_combo_box, False, False, 0)
      index_method_box2.pack_start(padded_check_button, False, False, 4)
      min_len_label.show()
      min_len_text_entry.show()
      max_block_label.show()
      max_block_text_entry.show()
      suffix_method_label.show()
      suffix_method_combo_box.show()
      padded_check_button.show()
      index_method_box2.show()

    # BigMatch and Dedup index are implementation variations of the Blocking,
    # Sorting and QGram indices
    #
    if (index_method_name in ['BlockingIndex', 'SortingIndex', 'QGramIndex']):

      if (index_method_name == 'BlockingIndex'):
        special_index_box_widget = index_method_box2
      else:
        special_index_box_widget = index_method_box3

      space_label = gtk.Label('     ')  # Make it indented
      special_index_box_widget.pack_start(space_label, False, False, 0)
      space_label.show()

      if (self.project_type == 'Link'):
        special_index_check_button = gtk.CheckButton('Use BigMatch indexing')

        if ('doBigMatch' in index_dict):
          do_special_tick = index_dict['doBigMatch']
        else:
          do_special_tick = True  # Activate per default

      else:
        special_index_check_button = gtk.CheckButton('Use Dedup indexing')

        if ('doDedup' in index_dict):
          do_special_tick = index_dict['doDedup']
        else:
          do_special_tick = True  # Activate per default

      special_index_check_button.set_active(do_special_tick)
      special_index_box_widget.pack_start(special_index_check_button,
                                          False, False,0)
      special_index_check_button.show()
      special_index_box_widget.show()

    # Display the separator between index method and index definitions - - - -
    #
    indexing_separator.show()

    # Create the box for the index definitions in scrolled window - - - - - - -
    #
    index_viewport = index_scrolled_window.get_child()

    if (index_method_name == 'FullIndex'):  # Hide the index definition window

      if (index_viewport != None):
        index_viewport.hide()
      return  # No index defintions needed for full index

    if (index_viewport == None):  # Create new box if none is there yet
      index_def_box = gtk.VBox()
      index_scrolled_window.add_with_viewport(index_def_box)
    else:  # Get the box with index definition details (1st child is viewport)
      index_viewport.show()
      index_def_box = index_scrolled_window.get_child().get_child()

    for child in index_def_box.get_children():  # Remove all old children
      index_def_box.remove(child)

    # Now show existing index definitions or just a 'Add' button - - - - - - -
    #
    index_num = 1
    for index_list in self.index_def:
      print 'index_list:', index_list

      index_num_label = gtk.Label('<b>Index %d:</b>' % (index_num))
      index_num_label.set_use_markup(True)
      index_def_box.pack_start(index_num_label, False, False, 0)
      index_num_label.set_alignment(0.0, 0.5)
      index_num_label.show()

      for index_def_tuple in index_list:  # Each index is one or more tuple
        print 'index_def:', index_def_tuple

        [new_index_def1, new_index_def2] = self.new_index_def(index_def_tuple)
        index_def_box.pack_start(new_index_def1, False, False, 0)
        index_def_box.pack_start(new_index_def2, False, False, 0)

      last_index_def_len = len(index_list)

      # Display separators between indices
      #
      horiz_separator = gtk.HSeparator()
      index_def_box.pack_start(horiz_separator, False, False, 0)
      horiz_separator.show()
      index_num += 1

    index_num -= 1

    # Finally add 'Add' and 'Delete' buttons - - - - - - - - - - - - - - - -
    #
    if (index_num == 0):  # No index defined yet, only show 'Add index' button
      del_index =     False
      add_index_def = False
      del_index_def = False

    else:
      add_index_def = True

      if (index_num == 1):  # Only one index, so don't show 'Delete last index'
        del_index = False
      else:
        del_index = True

      if (last_index_def_len == 1):  # Only one last index definition
        del_index_def = False
      else:
        del_index_def = True

#    if (index_num <= 1):  # Don't allow delete if no or only one index
#      del_index = False
#      add_index_def = False
#    else:
#      del_index = True
#      add_index_def = True

    # Check if only one index definition is left for this index
    #
#    if ((len(index_def_box.get_children()) < 4) or \
#        (type(index_def_box.get_children()[-4]) == type(gtk.Label))):
#      del_index_def = False
#    else:
#      del_index_def = True

    print add_index_def, del_index_def, del_index
###

    button_box = self.create_index_buttons(add_index_def, del_index_def,
                                           del_index)
    index_def_box.pack_start(button_box, False, False, 0)

    index_def_box.show()
    index_scrolled_window.show()

  # ---------------------------------------------------------------------------
  # Handle a click on the 'Add new index' button
  #
  def clickAddIndexButton(self, widget, add_index_box):
    print 'Clicked "Add new index"'

    scrolled_window = self.mainTree.get_widget('index_scrolled_window')

    # Get the box with index definition details (first child is the viewport)
    #
    index_def_box = scrolled_window.get_child().get_child()

    index_def_box.remove(add_index_box)  # Remove box with buttons

    self.index_num = self.index_num + 1  # One more indices

    index_num_label = gtk.Label('<b>Index %d:</b>' % (self.index_num))
    index_num_label.set_use_markup(True)
    index_def_box.pack_start(index_num_label, False, False, 0)
    index_num_label.set_alignment(0.0, 0.5)
    index_num_label.show()

    # Create a new index definiton (two HBoxes) with no values defined
    #
    [index_def_box1, index_def_box2] = self.new_index_def()
    index_def_box.pack_start(index_def_box1, False, False, 0)
    index_def_box.pack_start(index_def_box2, False, False, 0)

    horiz_separator = gtk.HSeparator()  # Display a separator
    index_def_box.pack_start(horiz_separator, False, False, 0)
    horiz_separator.show()

    # Finally add 'Add' and 'Delete' buttons - - - - - - - - - - - - - - - - -
    #
    if (self.index_num == 1):  # Don't allow delete if only one index
      del_index = False
    else:
      del_index = True

    # Create buttons, but not "Delete last index definition" button
    #
    button_box = self.create_index_buttons(True, False, del_index)
    index_def_box.pack_start(button_box, False, False, 0)

  # ---------------------------------------------------------------------------
  # Handle a click on the 'Add new index definition' button
  #
  def clickAddIndexDefButton(self, widget, add_index_box):
    print 'Clicked "Add new index definition"'

    scrolled_window = self.mainTree.get_widget('index_scrolled_window')

    # Get the box with index definition details (first child is the viewport)
    #
    index_def_box = scrolled_window.get_child().get_child()

    index_def_box.remove(add_index_box)  # Remove box with buttons
    index_def_box.remove(index_def_box.get_children()[-1])  # Remove separator

    # Create a new index definiton (two HBoxes) with no values defined
    #
    [index_def_box1, index_def_box2] = self.new_index_def()
    index_def_box.pack_start(index_def_box1, False, False, 0)
    index_def_box.pack_start(index_def_box2, False, False, 0)

    horiz_separator = gtk.HSeparator()  # Display a separator
    index_def_box.pack_start(horiz_separator, False, False, 0)
    horiz_separator.show()

    # Finally add 'Add' and 'Delete' buttons - - - - - - - - - - - - - - - - -
    #
    if (self.index_num == 1):  # Don't allow delete if only one index
      del_index = False
    else:
      del_index = True

    # Check if only one index definition is left for this index
    #
    if (type(index_def_box.get_children()[-4]) == type(gtk.Label)):
      del_index_def = False
    else:
      del_index_def = True

    button_box = self.create_index_buttons(True, del_index_def, del_index)
    index_def_box.pack_start(button_box, False, False, 0)

  # ---------------------------------------------------------------------------
  # Handle a click on the 'Delete last index' button
  #
  def clickDelIndexButton(self, widget, add_index_box):
    print 'Clicked "delete last index"'

    scrolled_window = self.mainTree.get_widget('index_scrolled_window')

    # Get the box with index definition details (first child is the viewport)
    #
    index_def_box = scrolled_window.get_child().get_child()

    index_def_box.remove(add_index_box)  # Remove box with buttons
    index_def_box.remove(index_def_box.get_children()[-1])  # Remove separator

    # Now remove the first two boxes containing the index definition details
    #
    index_def_box.remove(index_def_box.get_children()[-1])  # Remove encode box
    index_def_box.remove(index_def_box.get_children()[-1])  # Remove field box

    previous_widget_type = type(index_def_box.get_children()[-1])
    label_widget_type = type(gtk.Label())

    # Loop over all index definitions in this index
    #
    while (previous_widget_type != label_widget_type):
      index_def_box.remove(index_def_box.get_children()[-1])
      index_def_box.remove(index_def_box.get_children()[-1])
      previous_widget_type = type(index_def_box.get_children()[-1])

    # Remove 'Index %d' label
    #
    index_def_box.remove(index_def_box.get_children()[-1])

    self.index_num -= 1  # Reduce number of indices

    # Finally add 'Add' and 'Delete' buttons - - - - - - - - - - - - - - - - -
    #
    if (self.index_num == 1):  # Don't allow delete if only one index
      del_index = False
    else:
      del_index = True

    # Check if only one index definition is left for this index
    #
    if (type(index_def_box.get_children()[-4]) == type(gtk.Label)):
      del_index_def = False
    else:
      del_index_def = True

    button_box = self.create_index_buttons(True, del_index_def, del_index)
    index_def_box.pack_start(button_box, False, False, 0)

  # ---------------------------------------------------------------------------
  # Handle a click on the 'Delete last index definition' button
  #
  def clickDelIndexDefButton(self, widget, add_index_box):
    print 'Clicked "Delete last index definition"'

    scrolled_window = self.mainTree.get_widget('index_scrolled_window')

    # Get the box with index definition details (first child is the viewport)
    #
    index_def_box = scrolled_window.get_child().get_child()

    index_def_box.remove(add_index_box)  # Remove box with buttons
    index_def_box.remove(index_def_box.get_children()[-1])  # Remove separator

    # Now remove the box containing the index definition details
    #
    index_def_box.remove(index_def_box.get_children()[-1])  # Remove encode box
    index_def_box.remove(index_def_box.get_children()[-1])  # Remove field box

    horiz_separator = gtk.HSeparator()  # Display a separator
    index_def_box.pack_start(horiz_separator, False, False, 0)
    horiz_separator.show()

    # Finally add 'Add' and 'Delete' buttons - - - - - - - - - - - - - - - - -
    #
    if (self.index_num == 1):  # Don't allow delete if only one index
      del_index = False
    else:
      del_index = True

    # Check if only one index definition is left for this index
    #
    if (type(index_def_box.get_children()[-4]) == type(gtk.Label())):
      del_index_def = False
    else:
      del_index_def = True

    button_box = self.create_index_buttons(True, del_index_def, del_index)
    index_def_box.pack_start(button_box, False, False, 0)

  # ---------------------------------------------------------------------------
  # Construct a new HBox with four buttons to add/delete index (definitions)
  # Returns the HBox with the buttons.
  # If the 'add_index_def', 'del_index_def' or 'del_index' flags are set to
  # False then the corresponding buttons will not be created
  #
  def create_index_buttons(self, add_index_def, del_index_def=True,
                           del_index=True):
    index_button_box = gtk.HBox()

    if (add_index_def == True):
      add_index_def_button = gtk.Button('Add new index definition')
      add_index_def_button.connect('clicked', self.clickAddIndexDefButton,
                                   index_button_box)
      index_button_box.pack_start(add_index_def_button, False, False, 0)

    if (del_index_def == True):
      del_index_def_button = gtk.Button('Delete last index definition')
      del_index_def_button.connect('clicked', self.clickDelIndexDefButton,
                                   index_button_box)
      index_button_box.pack_start(del_index_def_button, False, False, 4)

    if ((add_index_def == True) or (del_index_def == True)):
      sep_label = gtk.Label('          ')
      index_button_box.pack_start(sep_label, False, False, 0)

    add_index_button =     gtk.Button('Add new index')
    add_index_button.connect('clicked', self.clickAddIndexButton,
                             index_button_box)
    index_button_box.pack_start(add_index_button, False, False, 4)

    if (del_index == True):
      del_index_button =     gtk.Button('Delete last index')
      del_index_button.connect('clicked', self.clickDelIndexButton,
                               index_button_box)
      index_button_box.pack_start(del_index_button, False, False, 0)

    if (add_index_def == True):
      add_index_def_button.show()
    if (del_index_def == True):
      del_index_def_button.show()
    if ((add_index_def == True) or (del_index_def == True)):
      sep_label.show()
    add_index_button.show()
    if (del_index == True):
      del_index_button.show()
    index_button_box.show()

    return index_button_box

  # ---------------------------------------------------------------------------
  # Construct a new index definition in two HBoxes, with values as given
  # (or empty/None). Returns the two index definition HBoxes
  #
  def new_index_def(self, index_def_tuple=None):

    # Extract index definition details
    #
    if (index_def_tuple != None):
      field_a_name =       index_def_tuple[0]
      field_b_name =       index_def_tuple[1]
      sort_words_flag =    index_def_tuple[2]
      reverse_flag =       index_def_tuple[3]
      max_len =            index_def_tuple[4]
      encode_funct =       index_def_tuple[5]
      encode_funct_param = index_def_tuple[6]
    else:
      field_a_name =       None
      field_b_name =       None
      sort_words_flag =    False
      reverse_flag =       False
      max_len =            None
      encode_funct =       None
      encode_funct_param = None

    index_def_box1 = gtk.HBox()  # A new box for encoding function details

    # Get the field names
    #
    field_names_a = self.data_set_info_list[0]['field_names']
    if (self.project_type == 'Link'):
      field_names_b = self.data_set_info_list[1]['field_names']
    else:
      field_names_b = field_names_a

    if (self.project_type == 'Link'):  # There will be two fields A and B
      field_a_label = gtk.Label('     Field name A:')
    else:  # One field only required, 'A' not needed
      field_a_label = gtk.Label('     Field name:')
    index_def_box1.pack_start(field_a_label, False, False, 0)

    field_a_combo_box = gtk.combo_box_new_text()
    for fn in field_names_a:
      field_a_combo_box.append_text(fn)
    index_def_box1.pack_start(field_a_combo_box, False, False, 0)

    if (field_a_name != None):
      field_name_index = field_names_a.index(field_a_name)
      field_a_combo_box.set_active(field_name_index)
    else:  # No name defined, set to None
      field_a_combo_box.set_active(-1)

    if (self.project_type == 'Link'):  # Only show for a linkage
      field_b_label = gtk.Label('  Field name B:')
      index_def_box1.pack_start(field_b_label, False, False, 4)

      field_b_combo_box = gtk.combo_box_new_text()
      for fn in field_names_b:
        field_b_combo_box.append_text(fn)
      index_def_box1.pack_start(field_b_combo_box, False, False, 0)

      if (field_b_name != None):
        field_name_index = field_names_b.index(field_b_name)
        field_b_combo_box.set_active(field_name_index)
      else:  # No name defined, set to None
        field_b_combo_box.set_active(-1)

    max_len_label = gtk.Label('  Maximum length:')
    index_def_box1.pack_start(max_len_label, False, False, 4)
    max_len_text_entry = gtk.Entry()
    max_len_text_entry.set_width_chars(5)
    index_def_box1.pack_start(max_len_text_entry, False, False, 0)

    if ((max_len != None) and (max_len != 'None')):
      max_len_text_entry.set_text(max_len)

    sort_button = gtk.CheckButton('Sort words')
    sort_button.set_active(sort_words_flag)
    index_def_box1.pack_start(sort_button, False, False, 10)

    rev_button = gtk.CheckButton('Reverse')
    rev_button.set_active(reverse_flag)
    index_def_box1.pack_start(rev_button, False, False, 10)

    field_a_label.show()
    field_a_combo_box.show()
    if (self.project_type == 'Link'):  # Only show for a linkage
      field_b_label.show()
      field_b_combo_box.show()
    max_len_label.show()
    max_len_text_entry.show()
    sort_button.show()
    rev_button.show()
    index_def_box1.show()

    index_def_box2 = gtk.HBox()  # A new box for encoding function details - -

    encode_label = gtk.Label('          Encoding function:')  # Indent
    index_def_box2.pack_start(encode_label, False, False, 0)

    encode_name_list = self.stringencode_dict.keys()  # Encode methods names
    encode_name_list.sort()

    encode_combo_box = gtk.combo_box_new_text()
    for e in encode_name_list:
      encode_combo_box.append_text(e)
    index_def_box2.pack_start(encode_combo_box, False, False, 0)

    if (encode_funct != None):
      encode_name_index = encode_name_list.index(encode_funct)
      encode_combo_box.set_active(encode_name_index)
    else:  # No name defined, set to None
      encode_combo_box.set_active(-1)

    encode_param_label = gtk.Label('  Encoding function parameters:')
    index_def_box2.pack_start(encode_param_label, False, False, 0)
    encode_param_text_entry = gtk.Entry()
    encode_param_text_entry.set_width_chars(10)
    index_def_box2.pack_start(encode_param_text_entry, False, False, 0)

    if (encode_funct_param != None):
      encode_param_text_entry.set_text(encode_funct_param)

    encode_label.show()
    encode_combo_box.show()
    encode_param_label.show()
    encode_param_text_entry.show()
    index_def_box2.show()

    return [index_def_box1, index_def_box2]

  # ---------------------------------------------------------------------------
  # Handle a change of the index method in combo box
  #
  def indexChangeMethod(self, widget):
    print 'Changed index method in combo box to:', widget.get_active()

    new_index_method = self.index_names[widget.get_active()]

    self.index_method['name'] = new_index_method  # Set method
    for k in self.index_method.keys():  # Remove all non general parameters
      if (k not in ['name', 'index_sep_str', 'skip_missing']):
        del self.index_method[k]

    self.indexView()  # Re-display the index page

  # ---------------------------------------------------------------------------
  # Handle an activate of Execute on the Index page
  #
  def indexExecute(self):
    indexing_widget = self.mainTree.get_widget('index_page_box')
    scrolled_window = self.mainTree.get_widget('index_scrolled_window')

    self.febrl_code['indexing'] = None  # Remove all previous indexing code

    # Get the selected index name from GUI
    #
    indexing_met_widget = self.mainTree.get_widget('index_method_combo_box')
    indexing_name_index = indexing_met_widget.get_active()

    if (indexing_name_index >= 0):
      index_method_name = self.index_names[indexing_name_index]

    else:  # No indexing method selected
      self.messageDialog('Please select an indexing method.', 'warn')
      return

    index_dict = self.index_method  # For quicker access

    # Get general parameter values - - - - - - - - - - - - - - - - - - - - - -
    #
    index_sep_str_widget =self.mainTree.get_widget('index_sep_string_entry')
    index_dict['index_sep_str'] = index_sep_str_widget.get_text()

    skip_missing_widget = \
                    self.mainTree.get_widget('index_skip_miss_check_button')
    skip_missing_val = skip_missing_widget.get_active()

    if ((index_method_name == 'CanopyIndex') and (skip_missing_val == False)):
      self.messageDialog('Skip missing must be set (ticked) for canopy index',
                         'warn')
      return
    index_dict['skip_missing'] = skip_missing_val

    # Get the children widgets of first box with parameters
    #
    if (index_method_name != 'FullIndex'):
      parameter_box1_children =indexing_widget.get_children()[1].get_children()

    # Get indexing method specific parameters - - - - - - - - - - - - - - - - -
    #
    if (index_method_name == 'SortingIndex'):  # - - - - - - - - - - - - - -
      win_size = parameter_box1_children[1].get_text().strip()
      if (not self.str_is_pos_int(win_size)):
        self.messageDialog('Window size must be a positive integer.', 'warn')
        return
      index_dict['window_size'] = win_size

    elif (index_method_name == 'QGramIndex'):  # - - - - - - - - - - - - - - -
      q = parameter_box1_children[1].get_text().strip()
      if (not self.str_is_pos_int(q)):
        self.messageDialog('Q must be a positive integer.', 'warn')
        return

      threshold = parameter_box1_children[3].get_text().strip()
      if (not self.str_is_normalised(threshold)):
        self.messageDialog('Threshold value must be larger or equal to 0\n'+ \
                           'and smaller or equal to 1.', 'warn')
        return
      if (threshold == '1'):
        threshold = '1.0'  # Make it a 'float' string
        parameter_box1_children[3].set_text(threshold)
      elif (threshold.startswith('.')):
        threshold = '0'+threshold
        parameter_box1_children[3].set_text(threshold)

      index_dict['q'] = q
      index_dict['threshold'] = threshold
      index_dict['padded'] = parameter_box1_children[4].get_active()

    elif (index_method_name == 'CanopyIndex'):  # - - - - - - - - - - - - - - -
      canopy_method = []

      canopy_method_active = parameter_box1_children[1].get_active()
      if (canopy_method_active == -1):
        self.messageDialog('Please select a canopy method.', 'warn')
        return
      elif (canopy_method_active == 0):
        canopy_method.append('tfidf')
      else:
        canopy_method.append('jaccard')

      canopy_using_active = parameter_box1_children[3].get_active()
      if (canopy_using_active == -1):
        self.messageDialog('Please select a threshold method.', 'warn')
        return
      elif (canopy_using_active == 0):
        canopy_method.append('nearest')
      else:
        canopy_method.append('threshold')
      p = parameter_box1_children[5].get_text().strip()
      if (p.count(',') != 1):  # One comma required to separate two parameters
        self.messageDialog('Two parameter values required, use comma ","\n'+ \
                           'to separate them.', 'warn')
        return
      p1,p2 = p.split(',')
      p1,p2 = p1.strip(), p2.strip()

      if (canopy_method[-1] == 'nearest'):
        if ((not self.str_is_pos_int(p1)) or (not self.str_is_pos_int(p2))):
          self.messageDialog('Parameter values for nearest neighbour ' + \
                             'method must be positive integers.', 'warn')
          return
        if (int(p1) > int(p2)):
          self.messageDialog('First nearest neighbour value must be\n' + \
                             'smaller than second value.', 'warn')
          return
      else:
        if ((not self.str_is_normalised_not_zero(p1)) or \
            (not self.str_is_normalised_not_zero(p2))):
          self.messageDialog('Parameter values for global thresholds ' + \
                             'method must\nbe larger than 0 and smaller or'+ \
                             ' equal to 1.', 'warn')
          return
        if (float(p1) < float(p2)):
          self.messageDialog('First global threshold value must be\n' + \
                             'larger than second value.', 'warn')
          return

      canopy_method.append(p)

      # Get second box with parameters
      #
      parameter_box2_children =indexing_widget.get_children()[2].get_children()

      q = parameter_box2_children[1].get_text().strip()
      if (not self.str_is_pos_int(q)):
        self.messageDialog('Q must be a positive integer.', 'warn')
        return

      del_perc = parameter_box2_children[3].get_text().strip()
      if (not self.str_is_percentage_not_zero(del_perc)):
        self.messageDialog('Delete percentage value must be larger than 0\n'+ \
                           'and smaller or equal to 100.\nSet to 100 if no' + \
                           ' q-grams shall be deleted.', 'warn')
        return

      index_dict['q'] = q
      index_dict['canopy_method'] = canopy_method
      index_dict['delete_perc'] = del_perc
      index_dict['padded'] = parameter_box2_children[4].get_active()

    elif (index_method_name == 'StringMapIndex'):  # - - - - - - - - - - - - -
      canopy_method = []
      canopy_using_active = parameter_box1_children[1].get_active()
      if (canopy_using_active == -1):
        self.messageDialog('Please select a threshold method', 'warn')
        return
      elif (canopy_using_active == 0):
        canopy_method.append('nearest')
      else:
        canopy_method.append('threshold')
      p = parameter_box1_children[3].get_text().strip()
      if (p.count(',') != 1):  # One comma required to separate two parameters
        self.messageDialog('Two parameter values required, use comma ","\n'+ \
                           'to separate them.', 'warn')
        return
      p1,p2 = p.split(',')
      p1,p2 = p1.strip(), p2.strip()

      if (canopy_method[-1] == 'nearest'):
        if ((not self.str_is_pos_int(p1)) or (not self.str_is_pos_int(p2))):
          self.messageDialog('Parameter values for nearest neighbour ' + \
                             'method must be positive integers.', 'warn')
          return
        if (int(p1) > int(p2)):
          self.messageDialog('First nearest neighbour value must be\n' + \
                             'smaller than second value.', 'warn')
          return
      else:
        if ((not self.str_is_normalised_not_zero(p1)) or \
            (not self.str_is_normalised_not_zero(p2))):
          self.messageDialog('Parameter values for global thresholds ' + \
                             'method must\nbe larger than 0 and smaller or'+ \
                             ' equal to 1.', 'warn')
          return
        if (float(p1) < float(p2)):
          self.messageDialog('First global threshold value must be\n' + \
                             'larger than second value.', 'warn')
          return

      canopy_method.append(p)

      grid_res = parameter_box1_children[5].get_text().strip()
      if ((not self.str_is_pos_int(grid_res)) or \
          (grid_res not in ['10','100','1000','10000'])):
        self.messageDialog('Grid resoultion must be a power of 10 value\n' + \
                           '(10, 100, 1000 or 10000).', 'warn')
        return

      # Get second box with parameters
      #
      parameter_box2_children =indexing_widget.get_children()[2].get_children()

      dim = parameter_box2_children[1].get_text().strip()
      if (not self.str_is_pos_int(dim)):
        self.messageDialog('Dimension must be a positive integer.', 'warn')
        return

      sub_dim = parameter_box2_children[3].get_text().strip()
      if (not self.str_is_pos_int(sub_dim)):
        self.messageDialog('Sub-dimension must be a positive integer.', 'warn')
        return
      if (int(sub_dim) > int(dim)):
        self.messageDialog('Sub-dimension must be smaller or equal to ' + \
                           'dimension.', 'warn')
        return

      sim_funct_active = parameter_box2_children[5].get_active()
      if (sim_funct_active == -1):
        self.messageDialog('Please select a similarity function', 'warn')
        return
      sim_funct_name_list = self.stringcmp_dict.keys()
      sim_funct_name_list.sort()

      index_dict['canopy_method'] = canopy_method
      index_dict['grid_resolution'] = grid_res
      index_dict['dim'] = dim
      index_dict['sub_dim'] = sub_dim

      index_dict['sim_funct'] = sim_funct_name_list[sim_funct_active]
      index_dict['cache_dist'] = parameter_box2_children[6].get_active()

    elif (index_method_name == 'SuffixArrayIndex'):  # - - - - - - - - - - - -
      min_len = parameter_box1_children[1].get_text().strip()
      if (not self.str_is_pos_int(min_len)):
        self.messageDialog('Minimum length must be a positive integer.',
                           'warn')
        return
      max_block = parameter_box1_children[3].get_text().strip()
      if (not self.str_is_pos_int(max_block)):
        self.messageDialog('Maximum block size must be a positive integer.',
                           'warn')
        return

      suffix_method_active = parameter_box1_children[5].get_active()
      if (suffix_method_active == -1):
        self.messageDialog('Please select a suffix method', 'warn')
        return
      elif (suffix_method_active == 0):
        index_dict['suffix_method'] = 'allsubstr'
      else:
        index_dict['suffix_method'] = 'suffixonly'

      index_dict['block_method'] = [min_len, max_block]
      index_dict['padded'] = parameter_box1_children[6].get_active()

    # Get check box values for BigMatch or Dedup - - - - - - - - - - - - - - -
    #
    if (index_method_name in ['BlockingIndex', 'SortingIndex', 'QGramIndex']):

      if (index_method_name == 'BlockingIndex'):
        check_box_list = parameter_box1_children
      else:
        check_box_list = indexing_widget.get_children()[2].get_children()

      if (self.project_type == 'Link'):
        index_dict['doBigMatch'] = check_box_list[1].get_active()
      else:
        index_dict['doDedup'] = check_box_list[1].get_active()

    print 'index_dict:', index_dict

    # Process index definition details - - - - - - - - - - - - - - - - - - - -
    #
    scrolled_window = self.mainTree.get_widget('index_scrolled_window')

    if ((scrolled_window.get_child() == None) and
        (index_method_name != 'FullIndex')):
      self.messageDialog('At least one index definition required.', 'warn')
      return

    # No index definitions or parameters for FullIndex
    #
    if (index_method_name != 'FullIndex'):

      # Get the field names
      #
      field_names_a = self.data_set_info_list[0]['field_names']
      if (self.project_type == 'Link'):
        field_names_b = self.data_set_info_list[1]['field_names']
      else:
        field_names_b = field_names_a

      encode_name_list = self.stringencode_dict.keys()  # Encode methods names
      encode_name_list.sort()

      index_def_box =      scrolled_window.get_child().get_child()
      index_def_box_list = index_def_box.get_children()
      num_index_boxes =    len(index_def_box_list)

      if (num_index_boxes < 3):
        self.messageDialog('At least one index definition required.', 'warn')
        return

      self.index_def = []  # Clear old index definitions

      i = 1  # List in index_def_box, item 0 is the first index label
      index_def_list = []  # List with tuples for a new index

      while (i < num_index_boxes):  # Loop over boxes with index and index def.

        # Get first index defintion (2 boxes)
        #
        index_def_box1_children = index_def_box_list[i].get_children()
        index_def_box2_children = index_def_box_list[i+1].get_children()

        field_name_a_ind = index_def_box1_children[1].get_active()
        if (self.project_type == 'Link'):
          field_name_b_ind = index_def_box1_children[3].get_active()
          max_len =          index_def_box1_children[5].get_text().strip()
          sort_words_flag =  index_def_box1_children[6].get_active()
          reverse_flag =     index_def_box1_children[7].get_active()
        else:  # Only one field combo box displayed
          field_name_b_ind = field_name_a_ind
          max_len =          index_def_box1_children[3].get_text().strip()
          sort_words_flag =  index_def_box1_children[4].get_active()
          reverse_flag =     index_def_box1_children[5].get_active()

        if ((field_name_a_ind == -1) or (field_name_b_ind == -1)):
          self.messageDialog('Not all field names are defined!', 'warn')
          return
        field_name_a = field_names_a[field_name_a_ind]
        field_name_b = field_names_b[field_name_b_ind]

        if (max_len == ''):  # No maximum length limit provided
          max_len = 'None'
        else:  # Check max_len is an integer
          if (not self.str_is_pos_int(max_len)):
            self.messageDialog('Maximum length must be empty or a positive' + \
                               ' integer.', 'warn')
            return

        encode_funct_ind =    index_def_box2_children[1].get_active()
        endcode_funct_param = index_def_box2_children[3].get_text().strip()

        if (encode_funct_ind == -1):
          encode_funct_name = 'None'
        else:
          encode_funct_name = encode_name_list[encode_funct_ind]

        new_index_def = (field_name_a, field_name_b, sort_words_flag,
                         reverse_flag, max_len, encode_funct_name,
                         endcode_funct_param)
        if (encode_funct_name == 'Substring'):
          if ((len(endcode_funct_param.split(',')) != 2) or \
              (not endcode_funct_param.split(',')[0].isdigit()) or \
              (not endcode_funct_param.split(',')[1].isdigit())):
            self.messageDialog('Parameters for substring encoding function' + \
                               ' must be two integers (start, end).', 'warn')
            return

        index_def_list.append(new_index_def)

        i += 2

        # Check if next box contains a separator, if so end of this index
        #
        if (type(index_def_box_list[i]) == type(gtk.HSeparator())):
          self.index_def.append(index_def_list)
          index_def_list = []  # Clear for a new index
          i += 2  # Skip separator and 'Index' label

      print 'index_def:', self.index_def

    # Now generate the Febrl codes for indexing - - - - - - - - - - - - - - - -
    #
    febrl_code = []
    febrl_code.append('# '+'-'*77)
    febrl_code.append('')
    febrl_code.append('# Define indices for "blocking"')
    febrl_code.append('#')

    # Process index definitions first
    #
    if (index_method_name != 'FullIndex'):

      index_num = 1
      for index_list in self.index_def:

        index_def_str = 'index_def_%d = [' % (index_num)
        index_def_indent_str = ' '*(14+len('%d' % (index_num)))

        num_index_tuples = len(index_list)
        index_def_num = 0
        for index_def_tuple in index_list:  # Each index is one or more tuple
          field_a_name =       index_def_tuple[0]
          field_b_name =       index_def_tuple[1]
          sort_words_flag =    index_def_tuple[2]
          reverse_flag =       index_def_tuple[3]
          max_len =            index_def_tuple[4]
          encode_funct =       index_def_tuple[5]
          encode_funct_param = index_def_tuple[6]

          this_index_def_str = '["%s", "%s", %s, %s, %s, ' % (field_a_name,
                               field_b_name, str(sort_words_flag),
                               str(reverse_flag), max_len)
          if (encode_funct == 'None'):
            this_index_def_str += '[]'
          else:
            encode_funct_name = self.stringencode_dict[encode_funct]
            this_index_def_str += '[encode.%s' % (encode_funct_name)

            if (encode_funct_param != ''):
              this_index_def_str += ', %s]' % (encode_funct_param)
            else:
              this_index_def_str += ']'

          if (index_def_num == num_index_tuples-1):
            this_index_def_str += ']]'  # Last index definition in this index
          else:
            this_index_def_str += '],'

          if (index_def_num == 0):
            febrl_code.append(index_def_str+this_index_def_str)
          else:
            febrl_code.append(index_def_indent_str+this_index_def_str)
          index_def_num += 1

        febrl_code.append('')
        index_num += 1

      index_num -= 1  # Number of indices

    index_method_name = self.index_method['name']
    sep_str =           self.index_method['index_sep_str']
    skip_missing =      self.index_method['skip_missing']

    if (index_method_name in ['BlockingIndex', 'SortingIndex', 'QGramIndex']):
      if ((self.project_type == 'Link') and \
          (self.index_method['doBigMatch'] == True)):
        index_method_name = 'BigMatchIndex'
      elif ((self.project_type != 'Link') and \
            (self.index_method['doDedup'])):
        index_method_name = 'DedupIndex'

    # Write Febrl python code - - - - - - - - - - - - - - - - - - - - - - - - -

    indention_space = ' '*(18+len(index_method_name))

    # Create code for index - General part first
    #
    febrl_code.append('index = indexing.%s(dataset1 = data_set_a,' % \
                           (index_method_name))
    if (self.project_type == 'Link'):
      febrl_code.append(indention_space+'dataset2 = data_set_b,')
    else:
      febrl_code.append(indention_space+'dataset2 = data_set_a,')
    febrl_code.append(indention_space+'rec_comparator = rec_comp,')
    febrl_code.append(indention_space+'index_sep_str = "%s",' % (sep_str))
    febrl_code.append(indention_space+'skip_missing = %s,' % \
                      (str(skip_missing)))

    if (index_method_name == 'FullIndex'):  # No index definition needed
      febrl_code.append(indention_space+'index_def = [])')

    else:
      index_def_str = 'index_def = [index_def_1'

      if (index_num == 1):  # Only one index
        febrl_code.append(indention_space+index_def_str+'],')
      else:
        febrl_code.append(indention_space+index_def_str+',')

        for i in range(2, index_num):
          febrl_code.append(indention_space+' '*13+'index_def_%d,' % (i))
        febrl_code.append(indention_space+' '*13+'index_def_%d],' % \
                          (index_num))

    # Create code for index - Indexing technique specific parameters
    #
    if (index_method_name == 'BlockingIndex'):  # No parameter, modify code
      last_febrl_code = febrl_code.pop()
      last_febrl_code = last_febrl_code[:-1] + ')'
      febrl_code.append(last_febrl_code)

    elif (index_method_name == 'SortingIndex'):
      win_size = index_dict['window_size']
      febrl_code.append(indention_space+'window_size = %s)' % (win_size))

    elif (index_method_name == 'QGramIndex'):
      q =         index_dict['q']
      threshold = index_dict['threshold']
      padded =    index_dict['padded']
      febrl_code.append(indention_space+'q = %s,' % (q))
      febrl_code.append(indention_space+'threshold = %s,' % (threshold))
      febrl_code.append(indention_space+'padded = %s)' % (str(padded)))

    # Create block_method for BigMatch and Dedup index methods
    #
    elif (index_method_name in ['BigMatchIndex', 'DedupIndex']):
      if (self.index_method['name'] == 'BlockingIndex'):
        block_method_str = '("block",)'
      elif (self.index_method['name'] == 'SortingIndex'):
        win_size = index_dict['window_size']
        block_method_str = '("sort", %s)' % (win_size)
      else:  # QGramIndex
        q =         index_dict['q']
        threshold = index_dict['threshold']
        padded =    index_dict['padded']
        block_method_str = '("qgram", %s, %s, %s)' % (q,padded,str(threshold))
      febrl_code.append(indention_space+'block_method = %s)' % \
                        (block_method_str))

    elif (index_method_name == 'CanopyIndex'):
      canopy_method = index_dict['canopy_method']
      febrl_code.append(indention_space+'canopy_method = ("%s","%s",' % \
                        (canopy_method[0], canopy_method[1]) + '%s),' % \
                        (canopy_method[2]))
      q =        index_dict['q']
      padded =   index_dict['padded']
      del_perc = index_dict['delete_perc']
      febrl_code.append(indention_space+'q = %s,' % (q))
      febrl_code.append(indention_space+'delete_perc = %s,' % (del_perc))
      febrl_code.append(indention_space+'padded = %s)' % (str(padded)))

    elif (index_method_name == 'StringMapIndex'):
      canopy_method = index_dict['canopy_method']
      febrl_code.append(indention_space+'canopy_method = ("%s", %s),' % \
                        (canopy_method[0], canopy_method[1]))
      grid_res =   index_dict['grid_resolution']
      dim =        index_dict['dim']
      sub_dim =    index_dict['sub_dim']
      cache_dist = index_dict['cache_dist']
      sim_funct =  index_dict['sim_funct']
      sim_funct_name = 'stringcmp.'+self.stringcmp_dict[sim_funct]

      febrl_code.append(indention_space+'grid_resolution = %s,' % (grid_res))
      febrl_code.append(indention_space+'dim = %s,' % (dim))
      febrl_code.append(indention_space+'sub_dim = %s,' % (sub_dim))
      febrl_code.append(indention_space+'cache_dist = %s,' % (cache_dist))
      febrl_code.append(indention_space+'sim_funct = %s)' % (sim_funct_name))

    elif (index_method_name == 'SuffixArrayIndex'):
      suffix_method = index_dict['suffix_method']
      block_method =  index_dict['block_method']
      padded =        index_dict['padded']
      febrl_code.append(indention_space+'suffix_method = "%s",' % \
                        (suffix_method))
      febrl_code.append(indention_space+'block_method = (%s, %s),' % \
                        (block_method[0], block_method[1]))
      febrl_code.append(indention_space+'padded = %s)' % (str(padded)))

    febrl_code.append('')

    self.febrl_code['indexing'] = febrl_code  # Store for later use

    # Finally update the GUI information - - - - - - - - - - - - - - - - - - -
    #
    self.addToLog('')  # Add generated code into log page text
    self.addToLog('='*79)
    self.addToLog('Generated Febrl code for "indexing" on %s' % \
                  (time.asctime()))
    self.addToLog('')
    for line in febrl_code:
      self.addToLog(line)
    self.addToLog('')

    self.modified['indexing'] = True  # Indexing details have been changed
    self.setWindowTitle()

    self.re_run['w_vec_generate'] = True  # Need to re-generate weight vectors

    # Update the active and non-active notebook pages
    #
    self.main_notebook_page_active_dict['Evaluate'] = False

    if ((self.febrl_code['indexing'] != None) and \
        (self.febrl_code['comparison'] != None)):
      self.main_notebook_page_active_dict['Classify'] = True

    self.displayCurrentNotebookPage()  # Diplay the current page

    self.writeStatusBar('Generated Febrl Python code for indexing (see Log' + \
                  ' page for generated code).')

  # ===========================================================================
  # Methods that handle Compare page events
  # ===========================================================================

  def compareView(self):  # A switch to the Compare page ----------------------
    print '  Switched to Compare page - Display this page'

    comp_scrolled_window = self.mainTree.get_widget('comp_scrolled_window')

    # Create the box for the comparison definitions in scrolled window - - - -
    #
    comp_viewport = comp_scrolled_window.get_child()

    if (comp_viewport == None):  # Create new box if none is there yet
      comp_box = gtk.VBox()
      comp_scrolled_window.add_with_viewport(comp_box)
    else:  # Get the box with field comparison details (1ast child is viewport)
      comp_viewport.show()
      comp_box = comp_scrolled_window.get_child().get_child()

    for child in comp_box.get_children():  # Remove all old children
      comp_box.remove(child)

    num_comp_funct = len(self.field_comp_list)

    # Now show existing field comparisons - - - - - - - - - - - - - - - - - - -
    #
    for comp_funct in self.field_comp_list:
      print 'comp_funct:', comp_funct

      hbox_list = self.new_field_comp(comp_funct)
      for hbox in hbox_list:
        comp_box.pack_start(hbox, False, False, 0)
        hbox.show()

      horiz_separator = gtk.HSeparator()  # Display a horizontal separator
      comp_box.pack_start(horiz_separator, False, False, 0)
      horiz_separator.show()

    # Finally add 'Add' and 'Delete' buttons - - - - - - - - - - - - - - - -
    #
    if (num_comp_funct <= 1):  # No or only one comparison function defined so
      del_comp_funct = False   # far, only show 'Add...' button
    else:
      del_comp_funct = True

    comp_funct_button_box = self.create_comp_funct_button(del_comp_funct)

    comp_box.pack_start(comp_funct_button_box, False, False, 0)
    comp_funct_button_box.show()
    comp_box.show()

    comp_scrolled_window.show()

  # ---------------------------------------------------------------------------
  # Construct a new HBox with two buttons to add/delete comparison functions
  # Returns the HBox with the buttons.
  # If the 'del_comp_funct' flag is set to False then the corresponding button
  # will not be created
  #
  def create_comp_funct_button(self, del_comp_funct=True):

    comp_funct_button_box = gtk.HBox()

    add_button = gtk.Button('Add new comparison function')
    add_button.connect('clicked', self.clickAddCompFunctButton,
                       comp_funct_button_box)
    comp_funct_button_box.pack_start(add_button, False, False, 0)
    add_button.show()

    sep_label = gtk.Label('          ')
    comp_funct_button_box.pack_start(sep_label, False, False, 0)
    sep_label.show()

    if (del_comp_funct == True):
      del_button = gtk.Button('Delete last comparison function')
      del_button.connect('clicked', self.clickDelCompFunctButton,
                         comp_funct_button_box)
      comp_funct_button_box.pack_start(del_button, False, False, 0)
      del_button.show()

    return comp_funct_button_box

  # ---------------------------------------------------------------------------

  def clickAddCompFunctButton(self, widget, add_comp_funct_box):
    print 'Clicked "Add new comparison function"'

    comp_scrolled_window = self.mainTree.get_widget('comp_scrolled_window')

    # Get the box with field comparison functions (first child is the viewport)
    #
    comp_box = comp_scrolled_window.get_child().get_child()

    comp_box.remove(add_comp_funct_box)  # Remove box with buttons

    # Create a new field comparison function definiton with no values defined
    #
    hbox_list = self.new_field_comp()
    for hbox in hbox_list:
      comp_box.pack_start(hbox, False, False, 0)

    # Display final horizontal separator
    #
    horiz_separator = gtk.HSeparator()
    comp_box.pack_start(horiz_separator, False, False, 0)
    horiz_separator.show()

    # Check only one field comparison left - - - - - - - - - - - - - - - - - -
    #
    if (len(comp_box.get_children()) <= 5):
      del_comp_funct = False
    else:
      del_comp_funct = True

    # Finally add 'Add' and 'Delete' buttons - - - - - - - - - - - - - - - -
    #
    comp_funct_button_box = self.create_comp_funct_button(del_comp_funct)

    comp_box.pack_start(comp_funct_button_box, False, False, 0)
    comp_funct_button_box.show()

    self.field_comp_list.append({})  # A new dictionary for this comparison

  # ---------------------------------------------------------------------------

  def clickDelCompFunctButton(self, widget, add_comp_funct_box):
    print 'Clicked "Delete last comparison function"'

    comp_scrolled_window = self.mainTree.get_widget('comp_scrolled_window')

    # Get the box with field comparison functions (first child is the viewport)
    #
    comp_box = comp_scrolled_window.get_child().get_child()

    comp_box.remove(add_comp_funct_box)  # Remove box with buttons
    comp_box.remove(comp_box.get_children()[-1])  # Remove separator

    previous_widget_type =  type(comp_box.get_children()[-1])
    separator_widget_type = type(gtk.HSeparator())

    # Loop over all boxes for this field comparison (until separator found)
    #
    while (previous_widget_type != separator_widget_type):
      comp_box.remove(comp_box.get_children()[-1])
      previous_widget_type = type(comp_box.get_children()[-1])

    # Check only one field comparison left - - - - - - - - - - - - - - - - - -
    #
    if (len(comp_box.get_children()) <= 5):
      del_comp_funct = False
    else:
      del_comp_funct = True

    # Finally add 'Add' and 'Delete' buttons - - - - - - - - - - - - - - - -
    #
    comp_funct_button_box = self.create_comp_funct_button(del_comp_funct)

    comp_box.pack_start(comp_funct_button_box, False, False, 0)
    comp_funct_button_box.show()

    self.field_comp_list.pop()  # Remove last field comparison dictionary

  # ---------------------------------------------------------------------------
  # An activate of a field comparison function name
  #
  def compChangeFunction(self, widget):
    print 'Changed field comparison function to:', widget.get_active()

    field_comp_name_index = widget.get_active()

    if (field_comp_name_index == -1):  # Nothing to do as no name selected
      return

    comp_scrolled_window = self.mainTree.get_widget('comp_scrolled_window')

    # Get the box with field comparison functions (first child is the viewport)
    #
    comp_box = comp_scrolled_window.get_child().get_child()

    hbox_type = type(gtk.HBox())

    # Get the names of the field comparison functions
    #
    field_comp_names = self.field_comp_dict.keys()
    field_comp_names.sort()
    new_field_comp_name = field_comp_names[field_comp_name_index]
    print '  New name:', new_field_comp_name

    field_comp_index = 0
    child_index = 0
    for child in comp_box.get_children():
      if (type(child) == hbox_type):
        this_hbox_widgets = child.get_children()
        if ((this_hbox_widgets[0].get_text().startswith('Field')) and
            (this_hbox_widgets[1] == widget)):
          self.field_comp_list[field_comp_index]['name'] = new_field_comp_name
          break

      else:  # A separator, so next field comparison follows
        field_comp_index += 1

      child_index += 1

    self.compareView()

  # ---------------------------------------------------------------------------
  # Construct a new field comparison in three or more HBoxes, with values as
  # given in the field comparison dictionary. If an empty dictionary is given
  # (default), an empty new field comparison is created showing default values.
  # Returns a list of HBoxes.
  #
  def new_field_comp(self, field_comp_dict={}):

    # Get the names of the field comparison functions
    #
    field_comp_names = self.field_comp_dict.keys()
    field_comp_names.sort()

    # Get the field names
    #
    field_names_a = self.data_set_info_list[0]['field_names']
    if (self.project_type == 'Link'):
      field_names_b = self.data_set_info_list[1]['field_names']
    else:
      field_names_b = field_names_a

    comp_funct_box1 = gtk.HBox()  # First row with comparison function - - - -

    field_comp_label = gtk.Label('Field comparison function:')
    comp_funct_box1.pack_start(field_comp_label, False, False, 0)

    field_comp_combo_box = gtk.combo_box_new_text()
    for fcn in field_comp_names:
      field_comp_combo_box.append_text(fcn)
    if ('name' in field_comp_dict):
      field_comp_name_index = field_comp_names.index(field_comp_dict['name'])
      field_comp_combo_box.set_active(field_comp_name_index)
    comp_funct_box1.pack_start(field_comp_combo_box, False, False, 0)
    field_comp_combo_box.connect('changed', self.compChangeFunction)

    field_comp_label.show()
    field_comp_combo_box.show()
    comp_funct_box1.show()

    comp_funct_box2 = gtk.HBox()  # First row of parameters - - - - - - - - - -

    # Show the generic parameters
    #
    field_a_label = gtk.Label('     Field name A:')  # Indent
    comp_funct_box2.pack_start(field_a_label, False, False, 0)

    field_a_combo_box = gtk.combo_box_new_text()
    for fn in field_names_a:
      field_a_combo_box.append_text(fn)
    if ('field_a_name' in field_comp_dict):
      field_name_index = field_names_a.index(field_comp_dict['field_a_name'])
    else:
      field_name_index = -1
    field_a_combo_box.set_active(field_name_index)
    comp_funct_box2.pack_start(field_a_combo_box, False, False, 0)

    field_b_label = gtk.Label('  Field name B:')
    comp_funct_box2.pack_start(field_b_label, False, False, 0)

    field_b_combo_box = gtk.combo_box_new_text()
    for fn in field_names_b:
      field_b_combo_box.append_text(fn)
    if ('field_b_name' in field_comp_dict):
      field_name_index = field_names_b.index(field_comp_dict['field_b_name'])
    else:
      field_name_index = -1
    field_b_combo_box.set_active(field_name_index)
    comp_funct_box2.pack_start(field_b_combo_box, False, False, 0)

    cache_field_check_button = gtk.CheckButton('Cache comparisons  ')
    if ('cache_field' in field_comp_dict):
      cache_field_check_button.set_active(field_comp_dict['cache_field'])
    else:
      cache_field_check_button.set_active(False)  # Default don't cache
    comp_funct_box2.pack_start(cache_field_check_button, False, False, 4)

    cache_size_label = gtk.Label('  Maximum cache size:')
    comp_funct_box2.pack_start(cache_size_label, False, False, 0)

    cache_size_text_entry = gtk.Entry(7)
    cache_size_text_entry.set_width_chars(10)
    if ('max_cache_size' in field_comp_dict):
      cache_size_text_entry.set_text(field_comp_dict['max_cache_size'])
    else:
      cache_size_text_entry.set_text('None')
    comp_funct_box2.pack_start(cache_size_text_entry, False, False, 0)

    field_a_label.show()
    field_a_combo_box.show()
    field_b_label.show()
    field_b_combo_box.show()
    cache_field_check_button.show()
    cache_size_label.show()
    cache_size_text_entry.show()
    comp_funct_box2.show()

    comp_funct_box3 = gtk.HBox()  # Second row of parameters - - - - - - - - -

    miss_weight_label = gtk.Label('     Missing value weight:')  # Indent
    comp_funct_box3.pack_start(miss_weight_label, False, False, 0)

    miss_weight_text_entry = gtk.Entry()
    miss_weight_text_entry.set_width_chars(10)
    if ('missing_weight' in field_comp_dict):
      miss_weight_text_entry.set_text(field_comp_dict['missing_weight'])
    else:
      miss_weight_text_entry.set_text('0.0')
    comp_funct_box3.pack_start(miss_weight_text_entry, False, False, 0)

    agree_weight_label = gtk.Label('  Agreeing value weight:')
    comp_funct_box3.pack_start(agree_weight_label, False, False, 0)

    agree_weight_text_entry = gtk.Entry()
    agree_weight_text_entry.set_width_chars(10)
    if ('agree_weight' in field_comp_dict):
      agree_weight_text_entry.set_text(field_comp_dict['agree_weight'])
    else:
      agree_weight_text_entry.set_text('1.0')
    comp_funct_box3.pack_start(agree_weight_text_entry, False, False, 0)

    disagree_weight_label = gtk.Label('  Disagreeing value weight:')
    comp_funct_box3.pack_start(disagree_weight_label, False, False, 0)

    disagree_weight_text_entry = gtk.Entry()
    disagree_weight_text_entry.set_width_chars(10)
    if ('disagree_weight' in field_comp_dict):
      disagree_weight_text_entry.set_text(field_comp_dict['disagree_weight'])
    else:
      disagree_weight_text_entry.set_text('0.0')
    comp_funct_box3.pack_start(disagree_weight_text_entry, False, False, 0)

    miss_weight_label.show()
    miss_weight_text_entry.show()
    agree_weight_label.show()
    agree_weight_text_entry.show()
    disagree_weight_label.show()
    disagree_weight_text_entry.show()
    comp_funct_box3.show()

    field_comp_list = [comp_funct_box1, comp_funct_box2, comp_funct_box3]

    if ('name' not in field_comp_dict):
      return field_comp_list  # No more details to add

    comp_name = self.field_comp_dict[field_comp_dict['name']]
    print 'comp_name:', comp_name

    # Show the field comparison function specific parameters - - - - - - - - -
    #
    comp_funct_box4 = gtk.HBox()
    comp_funct_box5 = gtk.HBox()

    if (comp_name == 'FieldComparatorTruncateString'):
      num_comp_label = gtk.Label('     Number of characters to compare:')
      comp_funct_box4.pack_start(num_comp_label, False, False, 0)
      num_comp_text_entry = gtk.Entry(2)
      num_comp_text_entry.set_width_chars(5)
      if ('num_char_compared' in field_comp_dict):
        num_comp_text_entry.set_text(field_comp_dict['num_char_compared'])
      comp_funct_box4.pack_start(num_comp_text_entry, False, False, 0)
      num_comp_label.show()
      num_comp_text_entry.show()

    elif (comp_name == 'FieldComparatorKeyDiff'):  # - - - - - - - - - - - - -
      max_key_label = gtk.Label('     Maximum key difference:')
      comp_funct_box4.pack_start(max_key_label, False, False, 0)
      max_key_text_entry = gtk.Entry(2)
      max_key_text_entry.set_width_chars(5)
      if ('max_key_diff' in field_comp_dict):
        max_key_text_entry.set_text(field_comp_dict['max_key_diff'])
      comp_funct_box4.pack_start(max_key_text_entry, False, False, 0)
      max_key_label.show()
      max_key_text_entry.show()

    elif (comp_name == 'FieldComparatorNumericPerc'):  # - - - - - - - - - - -
      max_perc_label = gtk.Label('     Maximum percentage difference:')
      comp_funct_box4.pack_start(max_perc_label, False, False, 0)
      max_perc_text_entry = gtk.Entry(4)
      max_perc_text_entry.set_width_chars(5)
      if ('max_perc_diff' in field_comp_dict):
        max_perc_text_entry.set_text(field_comp_dict['max_perc_diff'])
      comp_funct_box4.pack_start(max_perc_text_entry, False, False, 0)
      max_perc_label.show()
      max_perc_text_entry.show()

    elif (comp_name == 'FieldComparatorNumericAbs'):  # - - - - - - - - - - - -
      max_abs_label = gtk.Label('     Maximum absolute difference:')
      comp_funct_box4.pack_start(max_abs_label, False, False, 0)
      max_abs_text_entry = gtk.Entry(4)
      max_abs_text_entry.set_width_chars(5)
      if ('max_abs_diff' in field_comp_dict):
        max_abs_text_entry.set_text(field_comp_dict['max_abs_diff'])
      comp_funct_box4.pack_start(max_abs_text_entry, False, False, 0)
      max_abs_label.show()
      max_abs_text_entry.show()

    elif (comp_name == 'FieldComparatorEncodeString'):  # - - - - - - - - - - -
      encode_label = gtk.Label('     Encoding function:')
      comp_funct_box4.pack_start(encode_label, False, False, 0)

      encode_name_list = self.stringencode_dict.keys()  # Encode methods names
      # Remove 'substring method, as this is not possible here
      #
      encode_name_list.remove('Substring')
      encode_name_list.remove('None')
      encode_name_list.sort()

      encode_combo_box = gtk.combo_box_new_text()
      for e in encode_name_list:
        encode_combo_box.append_text(e)
      comp_funct_box4.pack_start(encode_combo_box, False, False, 0)
      if ('encode_method' in field_comp_dict):
        encode_ind = encode_name_list.index(field_comp_dict['encode_method'])
        encode_combo_box.set_active(encode_ind)
      else:  # No name defined, set to None
        encode_combo_box.set_active(-1)

      reverse_check_button = gtk.CheckButton('Reverse')
      comp_funct_box4.pack_start(reverse_check_button, False, False, 4)
      if ('reverse' in field_comp_dict):
        reverse_check_button.set_active(field_comp_dict['reverse'])
      else:  # De-activate per default
        reverse_check_button.set_active(False)

      max_len_label = gtk.Label('  Maximum code length:')
      comp_funct_box4.pack_start(max_len_label, False, False, 4)
      max_len_text_entry = gtk.Entry()
      max_len_text_entry.set_width_chars(5)
      comp_funct_box4.pack_start(max_len_text_entry, False, False, 0)
      if ('max_code_length' in field_comp_dict):
        max_len_text_entry.set_text(field_comp_dict['max_code_length'])
      else:  # Set to default
        max_len_text_entry.set_text('4')

      encode_label.show()
      encode_combo_box.show()
      reverse_check_button.show()
      max_len_label.show()
      max_len_text_entry.show()

    #elif (comp_name == 'FieldComparatorDistance'):  # TODO ###################
      # geocode_table  A reference to the geocode look-up table (a dictionary)
      # max_distance   A positive number that gives the maximum distance (in
      #                kilometers) tolerated.

    elif (comp_name == 'FieldComparatorDate'):  # - - - - - - - - - - - - - - -
      before_label = gtk.Label('     Maximum day A before day B:')
      comp_funct_box4.pack_start(before_label, False, False, 0)

      before_text_entry = gtk.Entry()
      before_text_entry.set_width_chars(5)
      comp_funct_box4.pack_start(before_text_entry, False, False, 0)
      if ('max_day1_before_day2' in field_comp_dict):
        before_text_entry.set_text(field_comp_dict['max_day1_before_day2'])

      after_label = gtk.Label('  Maximum day A after day B:')
      comp_funct_box4.pack_start(after_label, False, False, 0)

      after_text_entry = gtk.Entry()
      after_text_entry.set_width_chars(5)
      comp_funct_box4.pack_start(after_text_entry, False, False, 0)
      if ('max_day2_before_day1' in field_comp_dict):
        after_text_entry.set_text(field_comp_dict['max_day2_before_day1'])

      date_format_label = gtk.Label('  Date format:')
      comp_funct_box4.pack_start(date_format_label, False, False, 0)
      date_format_text_entry = gtk.Entry(8)
      date_format_text_entry.set_width_chars(10)
      comp_funct_box4.pack_start(date_format_text_entry, False, False, 0)
      if ('date_format' in field_comp_dict):
        date_format_text_entry.set_text(field_comp_dict['date_format'])

      before_label.show()
      before_text_entry.show()
      after_label.show()
      after_text_entry.show()
      date_format_label.show()
      date_format_text_entry.show()

    elif (comp_name == 'FieldComparatorTime'):  # - - - - - - - - - - - - - - -
      before_label = gtk.Label('     Maximum time A before time B:')
      comp_funct_box4.pack_start(before_label, False, False, 0)

      before_text_entry = gtk.Entry()
      before_text_entry.set_width_chars(5)
      comp_funct_box4.pack_start(before_text_entry, False, False, 0)
      if ('max_time1_before_time2' in field_comp_dict):
        before_text_entry.set_text(field_comp_dict['max_time1_before_time2'])

      after_label = gtk.Label('  Maximum time A after time B:')
      comp_funct_box4.pack_start(after_label, False, False, 0)

      after_text_entry = gtk.Entry()
      after_text_entry.set_width_chars(5)
      comp_funct_box4.pack_start(after_text_entry, False, False, 0)
      if ('max_time2_before_time1' in field_comp_dict):
        after_text_entry.set_text(field_comp_dict['max_time2_before_time1'])

      day_start_label = gtk.Label('  Day start:')
      comp_funct_box4.pack_start(day_start_label, False, False, 0)

      day_start_text_entry = gtk.Entry()
      day_start_text_entry.set_width_chars(5)
      comp_funct_box4.pack_start(day_start_text_entry, False, False, 0)
      if ('day_start' in field_comp_dict):
        day_start_text_entry.set_text(field_comp_dict['day_start'])
      else:
        day_start_text_entry.set_text('00:00')

      before_label.show()
      before_text_entry.show()
      after_label.show()
      after_text_entry.show()
      day_start_label.show()
      day_start_text_entry.show()

    elif (comp_name == 'FieldComparatorAge'):  # - - - - - - - - - - - - - - -
      fix_date_label = gtk.Label('     Fix date:')
      comp_funct_box4.pack_start(fix_date_label, False, False, 0)

      fix_date_text_entry = gtk.Entry(12)
      fix_date_text_entry.set_width_chars(10)
      comp_funct_box4.pack_start(fix_date_text_entry, False, False, 0)
      if ('fix_date' in field_comp_dict):
        fix_date_text_entry.set_text(field_comp_dict['fix_date'])
      else:
        fix_date_text_entry.set_text('today')

      max_perc_label = gtk.Label('  Maximum percentage difference:')
      comp_funct_box4.pack_start(max_perc_label, False, False, 0)
      max_perc_text_entry = gtk.Entry(4)
      max_perc_text_entry.set_width_chars(5)
      if ('max_perc_diff' in field_comp_dict):
        max_perc_text_entry.set_text(field_comp_dict['max_perc_diff'])
      comp_funct_box4.pack_start(max_perc_text_entry, False, False, 0)

      date_format_label = gtk.Label('  Date format:')
      comp_funct_box4.pack_start(date_format_label, False, False, 0)
      date_format_text_entry = gtk.Entry(8)
      date_format_text_entry.set_width_chars(10)
      comp_funct_box4.pack_start(date_format_text_entry, False, False, 0)
      if ('date_format' in field_comp_dict):
        date_format_text_entry.set_text(field_comp_dict['date_format'])

      fix_date_label.show()
      fix_date_text_entry.show()
      max_perc_label.show()
      max_perc_text_entry.show()
      date_format_label.show()
      date_format_text_entry.show()

    # All remaining comparators are approximate string comparisons - - - - - -
    #
    elif (comp_name not in ['FieldComparatorExactString',
                            'FieldComparatorContainsString']):

      # All approximate string comparators have a threshold parameter
      #
      threshold_label = gtk.Label('     Threshold:')
      comp_funct_box4.pack_start(threshold_label, False, False, 0)

      threshold_text_entry = gtk.Entry(10)
      threshold_text_entry.set_width_chars(10)
      comp_funct_box4.pack_start(threshold_text_entry, False, False, 0)
      threshold_text_entry.set_text(field_comp_dict.get('threshold', '0.0'))
      threshold_label.show()
      threshold_text_entry.show()

      if (comp_name == 'FieldComparatorWinkler'):  # - - - - - - - - - - - - -

        sim_check_button = gtk.CheckButton('Check similar characters')
        comp_funct_box4.pack_start(sim_check_button, False, False, 4)
        if ('check_sim' in field_comp_dict):
          sim_check_button.set_active(field_comp_dict['check_sim'])
        else:  # Activate per default
          sim_check_button.set_active(True)

        init_check_button = gtk.CheckButton('Check same initial ' + \
                                            'characters')
        comp_funct_box4.pack_start(init_check_button, False, False, 4)
        if ('check_init' in field_comp_dict):
          init_check_button.set_active(field_comp_dict['check_init'])
        else:  # Activate per default
          init_check_button.set_active(True)

        long_check_button = gtk.CheckButton('Check long strings')
        comp_funct_box4.pack_start(long_check_button, False, False, 4)
        if ('check_long' in field_comp_dict):
          long_check_button.set_active(field_comp_dict['check_long'])
        else:  # Activate per default
          long_check_button.set_active(True)

        sim_check_button.show()
        init_check_button.show()
        long_check_button.show()

      elif (comp_name in ['FieldComparatorQGram','FieldComparatorPosQGram']):
        q_label = gtk.Label('  Length of Q:')
        comp_funct_box4.pack_start(q_label, False, False, 0)
        q_text_entry = gtk.Entry(2)
        q_text_entry.set_width_chars(5)
        comp_funct_box4.pack_start(q_text_entry, False, False, 0)
        if ('q' in field_comp_dict):
          q_text_entry.set_text(field_comp_dict['q'])

        common_div_label = gtk.Label('  Common divisor:')
        comp_funct_box4.pack_start(common_div_label, False, False, 0)

        common_div_combo_box = gtk.combo_box_new_text()
        common_div_combo_box.append_text('Shortest')
        common_div_combo_box.append_text('Average')
        common_div_combo_box.append_text('Longest')
        comp_funct_box4.pack_start(common_div_combo_box, False, False, 0)
        if ('common_divisor' in field_comp_dict):
          if (field_comp_dict['common_divisor'] == 'shortest'):
            common_div_combo_box.set_active(0)
          elif (field_comp_dict['common_divisor'] == 'average'):
            common_div_combo_box.set_active(1)
          else:
            common_div_combo_box.set_active(2)
        else:
          common_div_combo_box.set_active(-1)

        padded_check_button = gtk.CheckButton('Padded')
        comp_funct_box4.pack_start(padded_check_button, False, False, 4)
        if ('padded' in field_comp_dict):
          padded_check_button.set_active(field_comp_dict['padded'])
        else:  # Activate per default
          padded_check_button.set_active(True)

        q_label.show()
        q_text_entry.show()
        common_div_label.show()
        common_div_combo_box.show()
        padded_check_button.show()

        if (comp_name == 'FieldComparatorPosQGram'):  # One more parameter - -
          max_dist_label = gtk.Label('   Maximum distance:')
          comp_funct_box4.pack_start(max_dist_label, False, False, 0)
          max_dist_text_entry = gtk.Entry(2)
          max_dist_text_entry.set_width_chars(5)
          comp_funct_box4.pack_start(max_dist_text_entry, False, False, 0)
          if ('max_dist' in field_comp_dict):
            max_dist_text_entry.set_text(field_comp_dict['max_dist'])
          max_dist_label.show()
          max_dist_text_entry.show()

      # Following comparison functions all have common divisor - - - - - - - -
      #
      elif (comp_name in ['FieldComparatorSGram','FieldComparatorSWDist',
                          'FieldComparatorSyllAlDist','FieldComparatorLCS',
                          'FieldComparatorOntoLCS', 'FieldComparatorTokenSet']):
        common_div_label = gtk.Label('  Common divisor:')
        comp_funct_box4.pack_start(common_div_label, False, False, 0)

        common_div_combo_box = gtk.combo_box_new_text()
        common_div_combo_box.append_text('Shortest')
        common_div_combo_box.append_text('Average')
        common_div_combo_box.append_text('Longest')
        comp_funct_box4.pack_start(common_div_combo_box, False, False, 0)
        if ('common_divisor' in field_comp_dict):
          if (field_comp_dict['common_divisor'] == 'shortest'):
            common_div_combo_box.set_active(0)
          elif (field_comp_dict['common_divisor'] == 'average'):
            common_div_combo_box.set_active(1)
          else:
            common_div_combo_box.set_active(2)
        else:
          common_div_combo_box.set_active(-1)

        common_div_label.show()
        common_div_combo_box.show()

        if (comp_name == 'FieldComparatorSGram'):  # - - - - - - - - - - - - -
          padded_check_button = gtk.CheckButton('Padded')
          comp_funct_box4.pack_start(padded_check_button, False, False, 4)
          if ('padded' in field_comp_dict):
            padded_check_button.set_active(field_comp_dict['padded'])
          else:  # Activate per default
            padded_check_button.set_active(True)

          gram_class_label = gtk.Label('  Gram class list:')
          comp_funct_box4.pack_start(gram_class_label, False, False, 0)
          gram_class_text_entry = gtk.Entry()
          #gram_class_text_entry.set_width_chars()
          comp_funct_box4.pack_start(gram_class_text_entry, True, True, 0)
          if ('gram_class_list' in field_comp_dict):
            gram_class_text_entry.set_text(field_comp_dict['gram_class_list'])

          padded_check_button.show()
          gram_class_label.show()
          gram_class_text_entry.show()

        elif (comp_name == 'FieldComparatorSyllAlDist'):  # - - - - - - - - - -
          phonix_check_button = gtk.CheckButton('Do Phonix encoding')
          comp_funct_box4.pack_start(phonix_check_button, False, False, 4)
          if ('do_phonix' in field_comp_dict):
            phonix_check_button.set_active(field_comp_dict['do_phonix'])
          else:  # Activate per default
            phonix_check_button.set_active(True)
          phonix_check_button.show()

        elif (comp_name in ['FieldComparatorLCS', 'FieldComparatorOntoLCS']):
          min_len_label = gtk.Label('  Minimum common length:')
          comp_funct_box4.pack_start(min_len_label, False, False, 0)
          min_len_text_entry = gtk.Entry(2)
          min_len_text_entry.set_width_chars(5)
          comp_funct_box4.pack_start(min_len_text_entry, False, False, 0)
          if ('min_common_len' in field_comp_dict):
            min_len_text_entry.set_text(field_comp_dict['min_common_len'])
          min_len_label.show()
          min_len_text_entry.show()

          if (comp_name == 'FieldComparatorOntoLCS'):
            p_label = gtk.Label('     Constant for Hamacher product ' + \
                                'difference:')
            comp_funct_box5.pack_start(p_label, False, False, 0)
            p_text_entry = gtk.Entry()
            p_text_entry.set_width_chars(5)
            comp_funct_box5.pack_start(p_text_entry, False, False, 0)
            if ('p' in field_comp_dict):
              p_text_entry.set_text(field_comp_dict['p'])
            p_label.show()
            p_text_entry.show()

        elif (comp_name == 'FieldComparatorTokenSet'):  # - - - - - - - - - - -
          stop_word_list_label = gtk.Label('  Stop word list:')
          comp_funct_box4.pack_start(stop_word_list_label, False, False, 0)
          stop_word_list_entry = gtk.Entry()
          comp_funct_box4.pack_start(stop_word_list_entry, True, True, 0)
          if ('stop_word_list' in field_comp_dict):
            stop_word_list_entry.set_text(field_comp_dict['stop_word_list'])
          stop_word_list_label.show()
          stop_word_list_entry.show()

      elif (comp_name == 'FieldComparatorCompress'):  # - - - - - - - - - - - -
        compress_label = gtk.Label('  Compressor:')
        comp_funct_box4.pack_start(compress_label, False, False, 0)

        compress_combo_box = gtk.combo_box_new_text()
        compress_combo_box.append_text('ZIP')
        compress_combo_box.append_text('BZIP2')
        comp_funct_box4.pack_start(compress_combo_box, False, False, 0)
        if ('compressor' in field_comp_dict):
          if (field_comp_dict['compressor'] == 'zlib'):
            compress_combo_box.set_active(0)
          elif (field_comp_dict['compressor'] == 'bz2'):
            compress_combo_box.set_active(1)
          else:
            raise Exception
        else:
          compress_combo_box.set_active(-1)

        compress_label.show()
        compress_combo_box.show()

#      else:
#        raise Exception, 'Illegal funct name: %s' % (comp_name)

    if (len(comp_funct_box4.get_children()) > 0):
      field_comp_list.append(comp_funct_box4)
    if (len(comp_funct_box5.get_children()) > 0):
      field_comp_list.append(comp_funct_box5)

    return field_comp_list

  # ---------------------------------------------------------------------------
  # Handle an activate of Execute on the Compare page
  #
  def compareExecute(self):
    comp_scrolled_window = self.mainTree.get_widget('comp_scrolled_window')
    comp_box = comp_scrolled_window.get_child().get_child()

    self.febrl_code['comparison'] = None  # Remove all previous comparison code

    if (len(comp_box) == 1):  # 'Add' button only
      self.messageDialog('At least one field comparison function is required.',
                         'error')
      return

    # Get the names of the field comparison functions
    #
    field_comp_names = self.field_comp_dict.keys()
    field_comp_names.sort()

    # Get the field names
    #
    field_names_a = self.data_set_info_list[0]['field_names']
    if (self.project_type == 'Link'):
      field_names_b = self.data_set_info_list[1]['field_names']
    else:
      field_names_b = field_names_a

    field_comp_list = self.field_comp_list  # List of field comparison details
    num_comp_funct =  len(field_comp_list)

    print 'Number of comparison functions:', num_comp_funct

    hbox_index = 0  # Counter in the horizontal boxes in the scrolled window

    comp_box_children = comp_box.get_children()

    # Main loop over field comparison functions - - - - - - - - - - - - - - - -
    #
    for fc in range(num_comp_funct):
      field_comp_dict = {}  # PC20/09field_comp_list[fc]  # Short cut

      comp_funct_box = comp_box_children[hbox_index]  # Box with comp. function
      comp_funct_box_children = comp_funct_box.get_children()
      field_comp_funct_index = comp_funct_box_children[1].get_active()

      if (field_comp_funct_index == -1):
        self.messageDialog('All field comparison functions have to be ' + \
                           'selected.', 'error')
        return

      # Name of this comparison function
      #
      field_comp_dict['name'] = field_comp_names[field_comp_funct_index]

      # Actual name of the comparison method in Febrl comparison mdule
      #
      comp_name = self.field_comp_dict[field_comp_dict['name']]

      hbox_index += 1  # Next box with field names - - - - - - - - - - - - - -

      comp_funct_box = comp_box_children[hbox_index]
      comp_funct_box_children = comp_funct_box.get_children()

      field_a_name_index = comp_funct_box_children[1].get_active()
      field_b_name_index = comp_funct_box_children[3].get_active()

      if ((field_a_name_index == -1) or (field_b_name_index == -1)):
        self.messageDialog('All field names have to be selected.', 'error')
        return

      field_comp_dict['field_a_name'] = field_names_a[field_a_name_index]
      field_comp_dict['field_b_name'] = field_names_b[field_b_name_index]
      field_comp_dict['cache_field'] =  comp_funct_box_children[4].get_active()

      # Only check cache size if cache is activated
      #
      if (field_comp_dict['cache_field'] == True):
        max_cache_size = comp_funct_box_children[6].get_text().strip()
        if (max_cache_size.lower() in ['', 'none', '(none)']):
          max_cache_size = 'None'
        if (max_cache_size != 'None'):
          if (not self.str_is_pos_int(max_cache_size)):
            self.messageDialog('Maximum cache size must be "None" or a ' + \
                             'positive integer.', 'warn')
            return
        field_comp_dict['max_cache_size'] = max_cache_size
      else:  # Not active set to None
        field_comp_dict['max_cache_size'] = 'None'

      hbox_index += 1  # Next box with weight parameters - - - - - - - - - - -

      comp_funct_box = comp_box_children[hbox_index]
      comp_funct_box_children = comp_funct_box.get_children()

      miss_weight =     comp_funct_box_children[1].get_text().strip()
      agree_weight =    comp_funct_box_children[3].get_text().strip()
      disagree_weight = comp_funct_box_children[5].get_text().strip()

      if ((not self.str_is_float(miss_weight)) or
          (not self.str_is_float(agree_weight)) or
          (not self.str_is_float(disagree_weight))):
        self.messageDialog('All weights must be set to numbers.', 'warn')
        return

      if (float(agree_weight) <= float(disagree_weight)):
        self.messageDialog('Agreement weights must be larger\nthan' + \
                           ' disagreement weights.', 'warn')
        return

      if (float(agree_weight) < float(miss_weight)):
        self.messageDialog('Agreement weights must be equal to\nor larger ' + \
                           'than missing weights.', 'warn')
        return

      if (float(miss_weight) < float(disagree_weight)):
        self.messageDialog('Missing weights must be equal to\nor larger ' + \
                           'than disagreement weights.', 'warn')
        return

      field_comp_dict['missing_weight'] = miss_weight
      field_comp_dict['agree_weight'] = agree_weight
      field_comp_dict['disagree_weight'] = disagree_weight

      # Exact and 'contains' string comparison do not need more parameters - -
      #
      if (comp_name not in ['FieldComparatorExactString',
                            'FieldComparatorContainsString']): # Process param.
        print 'comp_name:', comp_name

        hbox_index += 1  # Next box with function specific parameters

        comp_funct_box = comp_box_children[hbox_index]
        comp_funct_box_children = comp_funct_box.get_children()

        if (comp_name == 'FieldComparatorTruncateString'):  # - - - - - - - - -
          num_char_compared = comp_funct_box_children[1].get_text().strip()
          if (not self.str_is_pos_int(num_char_compared)):
            self.messageDialog('Number of characters to compare must\nbe a' + \
                               ' positive integer.', 'warn')
            return
          field_comp_dict['num_char_compared'] = num_char_compared

        elif (comp_name == 'FieldComparatorKeyDiff'):  # - - - - - - - - - - -
          max_key_diff = comp_funct_box_children[1].get_text().strip()
          if (not self.str_is_pos_int(max_key_diff)):
            self.messageDialog('Maximum key difference must\nbe a positive' + \
                               ' integer.', 'warn')
            return
          field_comp_dict['max_key_diff'] = max_key_diff

        elif (comp_name == 'FieldComparatorNumericPerc'):  # - - - - - - - - -
          max_perc_diff = comp_funct_box_children[1].get_text().strip()
          if (not self.str_is_percentage(max_perc_diff)):
            self.messageDialog('Maximum percentage difference must\nbe a ' + \
                               'percentage value.', 'warn')
            return
          field_comp_dict['max_perc_diff'] = max_perc_diff

        elif (comp_name == 'FieldComparatorNumericAbs'):  # - - - - - - - - - -
          max_abs_diff = comp_funct_box_children[1].get_text().strip()
          if (not self.str_is_pos_float(max_abs_diff)):
            self.messageDialog('Maximum absolute difference must\nbe zero ' + \
                               'or a positive number.', 'warn')
            return
          field_comp_dict['max_abs_diff'] = max_abs_diff

        elif (comp_name == 'FieldComparatorEncodeString'):  # - - - - - - - - -
          encode_name_list = self.stringencode_dict.keys()  # Encode methods
          # Remove 'substring method, as this is not possible here
          #
          encode_name_list.remove('Substring')
          encode_name_list.remove('None')
          encode_name_list.sort()

          encode_name_ind = comp_funct_box_children[1].get_active()
          if (encode_name_ind == -1):
            self.messageDialog('Encoding method names have to be selected.',
                               'error')
            return
          field_comp_dict['encode_method'] = encode_name_list[encode_name_ind]

          field_comp_dict['reverse'] = comp_funct_box_children[2].get_active()

          max_code_length = comp_funct_box_children[4].get_text().strip()
          if (not self.str_is_pos_int(max_code_length)):
            self.messageDialog('Maximum code length must\nbe a positive ' + \
                               'integer.', 'warn')
            return
          field_comp_dict['max_code_length'] = max_code_length

        #elif (comp_name == 'FieldComparatorDistance'):  # TODO ###############

        elif (comp_name == 'FieldComparatorDate'):  # - - - - - - - - - - - - -
          max_before_day = comp_funct_box_children[1].get_text().strip()
          max_after_day =  comp_funct_box_children[3].get_text().strip()
          date_format =  comp_funct_box_children[5].get_text().strip().lower()

          if (not self.str_is_not_neg_int(max_before_day)):
            self.messageDialog('Maximum day A before day B must\nbe zero ' + \
                               'or a positive number.', 'warn')
            return
          if (not self.str_is_not_neg_int(max_after_day)):
            self.messageDialog('Maximum day A after day B must\nbe zero ' + \
                               'or a positive number.', 'warn')
            return

          if (date_format not in ['ddmmyyyy', 'mmddyyyy', 'yyyymmdd',
                                  'ddmmyy', 'mmddyy']):
            self.messageDialog('Illegal date format value. Possible are:\n' + \
                               '"ddmmyyyy", "mmddyyyy", "yyyymmdd",\n' + \
                               '"ddmmyy", or "mmddyy".', 'warn')
            return

          field_comp_dict['max_day1_before_day2'] = max_before_day
          field_comp_dict['max_day2_before_day1'] = max_after_day
          field_comp_dict['date_format'] = date_format

        elif (comp_name == 'FieldComparatorTime'):  # - - - - - - - - - - - - -
          max_before_time = comp_funct_box_children[1].get_text().strip()
          max_after_time =  comp_funct_box_children[3].get_text().strip()
          if (not self.str_is_not_neg_int(max_before_time)):
            self.messageDialog('Maximum time A before time B must\nbe zero' + \
                               ' or a positive number.', 'warn')
            return
          if (not self.str_is_not_neg_int(max_after_time)):
            self.messageDialog('Maximum time A after time B must\nbe zero ' + \
                               'or a positive number.', 'warn')
            return
          field_comp_dict['max_time1_before_time2'] = max_before_time
          field_comp_dict['max_time2_before_time1'] = max_after_time

          day_start = comp_funct_box_children[5].get_text().strip()
          if (len(day_start) not in [4,5]):
            self.messageDialog('Day start must be of the form HH:MM or HHMM.',
                               'warn')
            return
          hours =   day_start[:2]   # First two digits
          minutes = day_start[-2:]  # Last two digits
          if ((not self.str_is_int(hours)) or (not self.str_is_int(minutes))):
            self.messageDialog('Day start must be of the form HH:MM or HHMM.',
                               'warn')
            return
          hrs = int(hours)
          mins = int(minutes)
          if ((hrs < 0) or (hrs > 23) or (mins < 0) or (mins > 59)):
            self.messageDialog('Day start values out of range: Hours must ' + \
                               'be [0..23]\nand minutes [0..59].', 'warn')
            return
          field_comp_dict['day_start'] = day_start

        elif (comp_name == 'FieldComparatorAge'):  # - - - - - - - - - - - - -
          fix_date = comp_funct_box_children[1].get_text().strip()
          max_perc_diff = comp_funct_box_children[3].get_text().strip()
          date_format =   comp_funct_box_children[5].get_text().strip().lower()

          if (fix_date != 'today'):  # Is not 'today' it must be a tuple
            if ((fix_date[0] != '(') or (fix_date[0] != '(') or
                (len(fix_date.split(',')) != 3)):
              self.messageDialog('Fix date must either be set to "today" ' + \
                                 'or\nbe a tuple of the form (DD,MM,YYY).',
                                 'warn')
              return
            day_str, month_str, year_str = fix_date.split(',')
            day_str =  day_str[1:]  # Remove brackets
            year_str = year_str[:-1]
            if ((len(day_str) != 2) or (len(month_str) != 2) or \
                (len(year_str) != 4)):
              self.messageDialog('Fix date must either be set to "today" ' + \
                                 'or\nbe a tuple of the form (DD,MM,YYYY).',
                                 'warn')
              return
            if ((day_str < '01') or (day_str > '31') or (month_str < '01') or \
                (month_str > '12') or (year_str < '1900') or \
                (year_str > '2030')):
              self.messageDialog('Fix date (DD,MM,YYYY) out of range: Day ' + \
                                 'must be [1..31],\nmonth [1..12] and year' + \
                                 ' [1900..2030].', 'warn')
              return

          if (not self.str_is_percentage(max_perc_diff)):
            self.messageDialog('Maximum percentage difference must\nbe a ' + \
                               'percentage value.', 'warn')
            return

          if (date_format not in ['ddmmyyyy', 'mmddyyyy', 'yyyymmdd',
                                  'ddmmyy', 'mmddyy']):
            self.messageDialog('Illegal date format value. Possible are:\n' + \
                               '"ddmmyyyy", "mmddyyyy", "yyyymmdd",\n' + \
                               '"ddmmyy", or "mmddyy".', 'warn')
            return
          field_comp_dict['fix_date'] = fix_date
          field_comp_dict['max_perc_diff'] = max_perc_diff
          field_comp_dict['date_format'] = date_format

        # All remaining functions are approximate string comparators - - - - -
        #
        else:
          threshold = comp_funct_box_children[1].get_text().strip()
          if ((not self.str_is_normalised(threshold)) or \
              (float(threshold) == 1.0)):
            self.messageDialog('Threshold value must be larger or equal\n'+ \
                               'to 0 and smaller to 1.', 'warn')
            return
          field_comp_dict['threshold'] = threshold

          if (comp_name == 'FieldComparatorWinkler'):  # - - - - - - - - - - -
            field_comp_dict['check_sim'] = \
                                        comp_funct_box_children[2].get_active()
            field_comp_dict['check_init'] = \
                                        comp_funct_box_children[3].get_active()
            field_comp_dict['check_long'] = \
                                        comp_funct_box_children[4].get_active()

          elif (comp_name in ['FieldComparatorQGram',
                              'FieldComparatorPosQGram']):  # - - - - - - - - -
            q = comp_funct_box_children[3].get_text().strip()
            if (not self.str_is_pos_int(q)):
              self.messageDialog('Value of "q" must be a positive integer.', \
                                 'warn')
              return
            field_comp_dict['q'] = q

            common_div_ind = comp_funct_box_children[5].get_active()
            if (common_div_ind == -1):
              self.messageDialog('Common divisor names have to be selected.',
                                 'error')
              return
            common_div_list = ['shortest','average','longest']
            field_comp_dict['common_divisor'] = common_div_list[common_div_ind]

            field_comp_dict['padded'] = comp_funct_box_children[6].get_active()

            if (comp_name == 'FieldComparatorPosQGram'):  # One more parameter
              max_dist = comp_funct_box_children[8].get_text().strip()
              if (not self.str_is_pos_int(max_dist)):
                self.messageDialog('Maximum distance must be a positive ' + \
                                   'integer.', 'warn')
                return
              field_comp_dict['max_dist'] = max_dist

          # Following comparison functions all have common divisor - - - - - -
          #
          elif (comp_name in ['FieldComparatorSGram','FieldComparatorSWDist',
                              'FieldComparatorSyllAlDist','FieldComparatorLCS',
                              'FieldComparatorOntoLCS',
                              'FieldComparatorTokenSet']):
            common_div_ind = comp_funct_box_children[3].get_active()
            if (common_div_ind == -1):
              self.messageDialog('Common divisor names have to be selected.',
                                 'error')
              return
            common_div_list = ['shortest','average','longest']
            field_comp_dict['common_divisor'] = common_div_list[common_div_ind]

            if (comp_name == 'FieldComparatorSGram'):  # - - - - - - - - - - -
              field_comp_dict['padded'] = \
                                        comp_funct_box_children[4].get_active()
              gram_class = comp_funct_box_children[6].get_text().strip()
              if ((gram_class == '') or (gram_class[0] != '[') or
                  (gram_class[-1] != ']') or (len(gram_class) < 4)):
                self.messageDialog('Gram-class value(s) must be a non-' + \
                                   'empty list', 'warn')
                return
              field_comp_dict['gram_class_list'] = gram_class

            elif (comp_name == 'FieldComparatorSyllAlDist'):  # - - - - - - - -
              field_comp_dict['do_phonix'] = \
                                        comp_funct_box_children[4].get_active()

            elif (comp_name in ['FieldComparatorLCS',
                                'FieldComparatorOntoLCS']):  # - - - - - - - -
              min_len = comp_funct_box_children[5].get_text().strip()
              if (not self.str_is_pos_int(min_len)):
                self.messageDialog('Minimum common length must be a ' + \
                                   'positive integer.', 'warn')
                return
              field_comp_dict['min_common_len'] = min_len

              if (comp_name == 'FieldComparatorOntoLCS'):  # - - - - - - - - -
                hbox_index += 1  # Next box with p parameters

                comp_funct_box = comp_box_children[hbox_index]
                comp_funct_box_children = comp_funct_box.get_children()

                p = comp_funct_box_children[1].get_text().strip()
                if (not self.str_is_normalised_not_zero(p)):
                  self.messageDialog('Constant for Hamacher product ' + \
                                     'difference must be\nlarger than zero' + \
                                     ' and smaller or equal to one.', 'warn')
                  return
                field_comp_dict['p'] = p

            elif (comp_name == 'FieldComparatorTokenSet'):  # - - - - - - - - -
              stop_word_list = comp_funct_box_children[5].get_text().strip()
              field_comp_dict['stop_word_list'] = stop_word_list

          elif (comp_name == 'FieldComparatorCompress'):  # - - - - - - - - - -
            compressor_index = comp_funct_box_children[3].get_active()
            if (compressor_index == -1):
              self.messageDialog('Compressor(s) have to be selected.', 'warn')
              return
            compressor_list = ['zlib','bz2']
            field_comp_dict['compressor'] = compressor_list[compressor_index]

      hbox_index += 2  # Jump over box with horizontal separator

      field_comp_list[fc] = field_comp_dict  # Save it into list PC 20/09

      print 'hbox_index:', hbox_index
      print 'FC-dict:', fc, field_comp_dict  #  TEST

    # Now generate the Febrl codes for comparisons - - - - - - - - - - - - - -
    #
    febrl_code = []
    febrl_code.append('# '+'-'*77)
    febrl_code.append('')
    febrl_code.append('# Define field comparison functions')
    febrl_code.append('#')

    field_comp_list_str = 'field_comp_list = ['

    for fc in range(num_comp_funct):
      field_comp_dict = field_comp_list[fc]  # Short cut

      field_a_name = field_comp_dict['field_a_name']
      field_b_name = field_comp_dict['field_b_name']

      comp_name = self.field_comp_dict[field_comp_dict['name']]

      descr_str = field_comp_dict['name']+'-'+field_a_name+'-'+field_b_name

      febrl_code.append('fc_funct_%d = comparison.%s(agree_weight = %s,' % \
                        (fc+1, comp_name, field_comp_dict['agree_weight']))

      indention_space = ' '*(24+len('%s' % (fc))+len(comp_name))

      febrl_code.append(indention_space+'description = "%s",' % (descr_str))

      febrl_code.append(indention_space+'disagree_weight = %s,' % \
                        (field_comp_dict['disagree_weight']))
      febrl_code.append(indention_space+'missing_weight = %s,' % \
                        (field_comp_dict['missing_weight']))
      if (field_comp_dict['cache_field'] == True):
        febrl_code.append(indention_space+'do_caching = True,')
        if (field_comp_dict['max_cache_size'] != 'None'):
          febrl_code.append(indention_space+'max_cache_size = %s,' % \
                            (field_comp_dict['max_cache_size']))

      # Handle field comparison function specific parameters - - - - - - - - -
      #
      if (comp_name == 'FieldComparatorTruncateString'):
        febrl_code.append(indention_space+'num_char_compared = %s,' % \
                          (field_comp_dict['num_char_compared']))

      elif (comp_name == 'FieldComparatorKeyDiff'):
        febrl_code.append(indention_space+'max_key_diff = %s,' % \
                          (field_comp_dict['max_key_diff']))

      elif (comp_name == 'FieldComparatorNumericPerc'):
        febrl_code.append(indention_space+'max_perc_diff = %s,' % \
                          (field_comp_dict['max_perc_diff']))

      elif (comp_name == 'FieldComparatorNumericAbs'):
        febrl_code.append(indention_space+'max_abs_diff = %s,' % \
                          (field_comp_dict['max_abs_diff']))

      elif (comp_name == 'FieldComparatorEncodeString'):
        encode_method = field_comp_dict['encode_method']
        febrl_code.append(indention_space+'encode_method = "%s",' % \
                          (self.stringencode_dict[encode_method]))
        febrl_code.append(indention_space+'reverse = %s,' % \
                          (str(field_comp_dict['reverse'])))
        febrl_code.append(indention_space+'max_code_length = %s,' % \
                          (field_comp_dict['max_code_length']))

      #elif (comp_name == 'FieldComparatorDistance'):  # TODO ###############

      elif (comp_name == 'FieldComparatorDate'):
        febrl_code.append(indention_space+'max_day1_before_day2 = %s,' % \
                          (field_comp_dict['max_day1_before_day2']))
        febrl_code.append(indention_space+'max_day2_before_day1 = %s,' % \
                          (field_comp_dict['max_day2_before_day1']))
        febrl_code.append(indention_space+'date_format = "%s",' % \
                          (field_comp_dict['date_format']))

      elif (comp_name == 'FieldComparatorTime'):
        febrl_code.append(indention_space+'max_time1_before_time2 = %s,' % \
                          (field_comp_dict['max_time1_before_time2']))
        febrl_code.append(indention_space+'max_time2_before_time1 = %s,' % \
                          (field_comp_dict['max_time2_before_time1']))
        febrl_code.append(indention_space+'day_start = "%s",' % \
                          (field_comp_dict['day_start']))

      elif (comp_name == 'FieldComparatorAge'):
        febrl_code.append(indention_space+'fix_date = "%s",' % \
                          (field_comp_dict['fix_date']))
        febrl_code.append(indention_space+'max_perc_diff = %s,' % \
                          (field_comp_dict['max_perc_diff']))
        febrl_code.append(indention_space+'date_format = "%s",' % \
                          (field_comp_dict['date_format']))

      # Hande the various approximate string comparators - - - - - - - - - - -
      #
      elif (comp_name not in ['FieldComparatorExactString',
                              'FieldComparatorContainsString']):
        febrl_code.append(indention_space+'threshold = %s,' % \
                          (field_comp_dict['threshold']))

        if (comp_name == 'FieldComparatorWinkler'):
          febrl_code.append(indention_space+'check_sim = %s,' % \
                            (str(field_comp_dict['check_sim'])))
          febrl_code.append(indention_space+'check_init = %s,' % \
                            (str(field_comp_dict['check_init'])))
          febrl_code.append(indention_space+'check_long = %s,' % \
                            (str(field_comp_dict['check_long'])))

        elif (comp_name in ['FieldComparatorQGram','FieldComparatorPosQGram']):
          febrl_code.append(indention_space+'q = %s,' % \
                            (field_comp_dict['q']))
          febrl_code.append(indention_space+'common_divisor = "%s",' % \
                            (field_comp_dict['common_divisor']))
          febrl_code.append(indention_space+'padded = %s,' % \
                            (str(field_comp_dict['padded'])))
          if (comp_name ==  'FieldComparatorPosQGram'):
            febrl_code.append(indention_space+'max_dist = %s,' % \
                             (field_comp_dict['max_dist']))

        elif (comp_name in ['FieldComparatorSGram','FieldComparatorSWDist',
                            'FieldComparatorSyllAlDist','FieldComparatorLCS',
                            'FieldComparatorOntoLCS',
                            'FieldComparatorTokenSet']):
          febrl_code.append(indention_space+'common_divisor = "%s",' % \
                            (field_comp_dict['common_divisor']))

          if (comp_name ==  'FieldComparatorSGram'):
            febrl_code.append(indention_space+'gram_class_list = %s,' % \
                              (field_comp_dict['gram_class_list']))
            febrl_code.append(indention_space+'padded = %s,' % \
                              (str(field_comp_dict['padded'])))

          elif (comp_name ==  'FieldComparatorSyllAlDist'):
            febrl_code.append(indention_space+'do_phonix = %s,' % \
                              (str(field_comp_dict['do_phonix'])))

          elif (comp_name == 'FieldComparatorLCS'):
            febrl_code.append(indention_space+'min_common_len = %s,' % \
                              (field_comp_dict['min_common_len']))

          elif (comp_name == 'FieldComparatorOntoLCS'):
            febrl_code.append(indention_space+'min_common_len = %s,' % \
                              (field_comp_dict['min_common_len']))
            febrl_code.append(indention_space+'p = %s,' % \
                              (field_comp_dict['p']))

          elif (comp_name == 'FieldComparatorTokenSet'):
            stop_word_list = field_comp_dict['stop_word_list'].split(',')
            stop_word_list_str = ''
            for stop_word in stop_word_list:
              if (stop_word.strip() != ''):
                stop_word_list_str += '"%s",' % (stop_word.strip())
            febrl_code.append(indention_space+'stop_word_list = [%s],' % \
                              (stop_word_list_str[:-1]))  # Remove last comma

        elif (comp_name == 'FieldComparatorCompress'):
          febrl_code.append(indention_space+'compressor = "%s",' % \
                            (field_comp_dict['compressor']))

      # Close the bracket in the last code line
      #
      last_febrl_code = febrl_code.pop()
      last_febrl_code = last_febrl_code[:-1] + ')'
      febrl_code.append(last_febrl_code)

      febrl_code.append('')

      field_comp_list_str += '(fc_funct_%d, "%s", "%s"),\n' % \
                              (fc+1, field_a_name, field_b_name) + ' '*19

    field_comp_list_str = field_comp_list_str[:-21]+']'

    febrl_code.append(field_comp_list_str)

    febrl_code.append('')

    # Add record comparator code - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.project_type == 'Link'):
      febrl_code.append('rec_comp = comparison.RecordComparator(' + \
                        'data_set_a, data_set_b, field_comp_list)')
    else:
      febrl_code.append('rec_comp = comparison.RecordComparator(' + \
                        'data_set_a, data_set_a, field_comp_list)')
    febrl_code.append('')

    self.febrl_code['comparison'] = febrl_code  # Store for later use

    # Finally update the GUI information - - - - - - - - - - - - - - - - - - -
    #
    self.addToLog('')  # Add generated code into log page text
    self.addToLog('='*79)
    self.addToLog('Generated Febrl code for "comparison" on %s' % \
                  (time.asctime()))
    self.addToLog('')
    for line in febrl_code:
      self.addToLog(line)
    self.addToLog('')

    self.modified['comparison'] = True  # Comparison details have been changed
    self.setWindowTitle()

    self.re_run['w_vec_generate'] = True  # Need to re-generate weight vectors

    # Update the active and non-active notebook pages
    #
    self.main_notebook_page_active_dict['Evaluate'] = False

    if ((self.febrl_code['indexing'] != None) and \
        (self.febrl_code['comparison'] != None)):
      self.main_notebook_page_active_dict['Classify'] = True

    self.displayCurrentNotebookPage()  # Diplay the current page

    self.writeStatusBar('Generated Febrl Python code for comparisons (see ' + \
                  'Log page for generated code).')

  # ===========================================================================
  # Methods that handle Classify page events
  # ===========================================================================

  # ---------------------------------------------------------------------------
  # Display the Classify page according to if a classifier has been
  # initialised or not
  #
  def classifyView(self):  # A switch to the Classify page
    print '  Switched to Classify page - Display this page'

    classification_widget = self.mainTree.get_widget('classify_page_box')
    classifier_method_widget = \
                        self.mainTree.get_widget('classifier_method_combo_box')

    # Remove all boxes on this page except the first one
    #
    classification_widget_children = classification_widget.get_children()
    if (len(classification_widget_children) > 1):
      for child in classification_widget_children[1:]:
        classification_widget.remove(child)  # Remove boxes with parameters

    classifier_dict = self.classifier_method  # For quicker access

    # Check if a classification method has previously been selected
    #
    if (classifier_dict['name'] == None):
      return  # Nothing to be done

    # Get the classification method name and set general parameters in GUI - -
    #
    classifier_method_name = classifier_dict['name']

    classifier_param_box1 = gtk.HBox()  # Fourboxes for parameters
    classifier_param_box2 = gtk.HBox()
    classifier_param_box3 = gtk.HBox()
    classifier_param_box4 = gtk.HBox()

    if (classifier_method_name == 'FellegiSunter'):  # - - - - - - - - - - - -

      lower_thres_label = gtk.Label('     Lower threshold:') # Make it indented
      lower_thres_text_entry = gtk.Entry()
      lower_thres_text_entry.set_width_chars(10)
      if ('lower_threshold' in classifier_dict):
        lower_thres_text_entry.set_text(classifier_dict['lower_threshold'])

      upper_thres_label = gtk.Label('  Upper threshold:')
      upper_thres_text_entry = gtk.Entry()
      upper_thres_text_entry.set_width_chars(10)
      if ('upper_threshold' in classifier_dict):
        upper_thres_text_entry.set_text(classifier_dict['upper_threshold'])

      classifier_param_box1.pack_start(lower_thres_label, False, False, 0)
      classifier_param_box1.pack_start(lower_thres_text_entry, False, False, 0)
      classifier_param_box1.pack_start(upper_thres_label, False, False, 0)
      classifier_param_box1.pack_start(upper_thres_text_entry, False, False, 0)
      lower_thres_label.show()
      lower_thres_text_entry.show()
      upper_thres_label.show()
      upper_thres_text_entry.show()

    elif (classifier_method_name == 'OptimalThreshold'):  # - - - - - - - - - -

      min_method_label = gtk.Label('     Minimise false method:')
      min_method_combo_box = gtk.combo_box_new_text()
      min_method_combo_box.append_text('Positives and negatives')
      min_method_combo_box.append_text('Positives')
      min_method_combo_box.append_text('Negatives')
      if ('min_method' in classifier_dict):
        if (classifier_dict['min_method'] == 'pos-neg'):
          min_method_combo_box.set_active(0)
        elif (classifier_dict['min_method'] == 'pos'):
          min_method_combo_box.set_active(1)
        else:
          min_method_combo_box.set_active(2)
      else:  # Default 'pos-neg' mimimise method
        min_method_combo_box.set_active(0)

      bin_width_label = gtk.Label('  Bin width:')
      bin_width_text_entry = gtk.Entry()
      bin_width_text_entry.set_width_chars(10)
      if ('bin_width' in classifier_dict):
        bin_width_text_entry.set_text(classifier_dict['bin_width'])

      classifier_param_box1.pack_start(min_method_label, False, False, 0)
      classifier_param_box1.pack_start(min_method_combo_box, False, False, 0)
      classifier_param_box1.pack_start(bin_width_label, False, False, 0)
      classifier_param_box1.pack_start(bin_width_text_entry, False, False, 0)
      min_method_label.show()
      min_method_combo_box.show()
      bin_width_label.show()
      bin_width_text_entry.show()

    elif (classifier_method_name in ['KMeans', 'FarthestFirst']):  # - - - - -

      dist_measure_label = gtk.Label('     Distance measure:') # Start indented
      dist_measure_combo_box = gtk.combo_box_new_text()
      dist_measure_combo_box.append_text('Manhatten')
      dist_measure_combo_box.append_text('Euclidean')
      dist_measure_combo_box.append_text('L-Infinity')
      dist_measure_combo_box.append_text('Canberra')
      if ('dist_measure' in classifier_dict):
        if (classifier_dict['dist_measure'] == 'Manhatten'):
          dist_measure_combo_box.set_active(0)
        elif (classifier_dict['dist_measure'] == 'Euclidean'):
          dist_measure_combo_box.set_active(1)
        elif (classifier_dict['dist_measure'] == 'L-Infinity'):
          dist_measure_combo_box.set_active(2)
        elif (classifier_dict['dist_measure'] == 'Canberra'):
          dist_measure_combo_box.set_active(3)
      else:
        dist_measure_combo_box.set_active(1)  # Default Euclidean distance

      sample_label = gtk.Label('  Sample:')
      sample_text_entry = gtk.Entry()
      sample_text_entry.set_width_chars(10)
      if ('sample' in classifier_dict):
        sample_text_entry.set_text(classifier_dict['sample'])

      if (classifier_method_name == 'KMeans'):
        max_iter_label = gtk.Label('     Maximum iteration count:')
        max_iter_text_entry = gtk.Entry()
        max_iter_text_entry.set_width_chars(10)
        if ('max_iter_count' in classifier_dict):
          max_iter_text_entry.set_text(classifier_dict['max_iter_count'])

      classifier_param_box1.pack_start(dist_measure_label, False, False, 0)
      classifier_param_box1.pack_start(dist_measure_combo_box, False, False, 0)
      classifier_param_box1.pack_start(sample_label, False, False, 0)
      classifier_param_box1.pack_start(sample_text_entry, False, False, 0)
      dist_measure_label.show()
      dist_measure_combo_box.show()
      sample_label.show()
      sample_text_entry.show()

      if (classifier_method_name == 'KMeans'):
        classifier_param_box1.pack_start(max_iter_label, False, False, 0)
        classifier_param_box1.pack_start(max_iter_text_entry, False, False, 0)
        max_iter_label.show()
        max_iter_text_entry.show()

      centroid_init_label = gtk.Label('     Centroid initialisation:')
      centroid_init_combo_box = gtk.combo_box_new_text()
      if (classifier_method_name == 'KMeans'):
        centroid_init_combo_box.append_text('Min/max')
        centroid_init_combo_box.append_text('Random')
        if ('centroid_init' in classifier_dict):
          if (classifier_dict['centroid_init'] == 'min/max'):
            centroid_init_combo_box.set_active(0)
          else:
            centroid_init_combo_box.set_active(1)
        else:  # Default 'min/max' centroid init method
          centroid_init_combo_box.set_active(0)

      else:  # Different initialisation methods for farthest first
        centroid_init_combo_box.append_text('Traditional')
        centroid_init_combo_box.append_text('Min/max')
        centroid_init_combo_box.append_text('Mode/max')
        if ('centroid_init' in classifier_dict):
          if (classifier_dict['centroid_init'] == 'traditional'):
            centroid_init_combo_box.set_active(0)
          elif (classifier_dict['centroid_init'] == 'min/max'):
            centroid_init_combo_box.set_active(1)
          else:
            centroid_init_combo_box.set_active(2)
        else:  # Default 'traditional' centroid init method
          centroid_init_combo_box.set_active(0)

      fuzzy_thres_label = gtk.Label('  Fuzzy region threshold:')
      fuzzy_thres_text_entry = gtk.Entry()
      fuzzy_thres_text_entry.set_width_chars(10)
      if ('fuzzy_reg_thres' in classifier_dict):
        fuzzy_thres_text_entry.set_text(classifier_dict['fuzzy_reg_thres'])

      classifier_param_box2.pack_start(centroid_init_label, False, False, 0)
      classifier_param_box2.pack_start(centroid_init_combo_box, False, False,0)
      classifier_param_box2.pack_start(fuzzy_thres_label, False, False, 0)
      classifier_param_box2.pack_start(fuzzy_thres_text_entry, False, False, 0)
      centroid_init_label.show()
      centroid_init_combo_box.show()
      fuzzy_thres_label.show()
      fuzzy_thres_text_entry.show()

    elif (classifier_method_name == 'SuppVecMachine'):  # - - - - - - - - - - -

      kernel_type_label = gtk.Label('     SVM kernel type:')
      kernel_type_combo_box = gtk.combo_box_new_text()
      kernel_type_combo_box.append_text('Linear')
      kernel_type_combo_box.append_text('Poly')
      kernel_type_combo_box.append_text('RBF')
      kernel_type_combo_box.append_text('Sigmoid')
      if ('kernel_type' in classifier_dict):
        if (classifier_dict['kernel_type'] == 'LINEAR'):
          kernel_type_combo_box.set_active(0)
        elif (classifier_dict['kernel_type'] == 'POLY'):
          kernel_type_combo_box.set_active(1)
        elif (classifier_dict['kernel_type'] == 'RBF'):
          kernel_type_combo_box.set_active(2)
        else:
          kernel_type_combo_box.set_active(3)
      else:
        kernel_type_combo_box.set_active(0)  # LINEAR is default

      C_label = gtk.Label('  C:')
      C_text_entry = gtk.Entry()
      C_text_entry.set_width_chars(10)
      if ('C' in classifier_dict):
        C_text_entry.set_text(classifier_dict['C'])

      sample_label = gtk.Label('  Sample:')
      sample_text_entry = gtk.Entry()
      sample_text_entry.set_width_chars(10)
      if ('sample' in classifier_dict):
        sample_text_entry.set_text(classifier_dict['sample'])

      classifier_param_box1.pack_start(kernel_type_label, False, False, 0)
      classifier_param_box1.pack_start(kernel_type_combo_box, False, False, 0)
      classifier_param_box1.pack_start(C_label, False, False, 0)
      classifier_param_box1.pack_start(C_text_entry, False, False, 0)
      classifier_param_box1.pack_start(sample_label, False, False, 0)
      classifier_param_box1.pack_start(sample_text_entry, False, False, 0)
      kernel_type_label.show()
      kernel_type_combo_box.show()
      C_label.show()
      C_text_entry.show()
      sample_label.show()
      sample_text_entry.show()

    elif (classifier_method_name == 'TwoStep'):  # - - - - - - - - - - - - - -

      m_value_label = gtk.Label('     Match comparison value:')
      m_value_text_entry = gtk.Entry()
      m_value_text_entry.set_width_chars(5)
      if ('s1_match_method' in classifier_dict):
        m_value_text_entry.set_text(classifier_dict['s1_match_method'][0])
      else:
        m_value_text_entry.set_text('1.0')  # Default value

      m_method_label = gtk.Label('  Match method:')
      m_method_combo_box = gtk.combo_box_new_text()
      m_method_combo_box.append_text('Nearest')
      m_method_combo_box.append_text('Threshold')
      if (('s1_match_method' in classifier_dict) and \
          (classifier_dict['s1_match_method'][1] == 'threshold')):
        m_method_combo_box.set_active(1)
        m_method_threshold = True
      else:
        m_method_combo_box.set_active(0)  # Set as 'Nearest' per default
        m_method_threshold = False

      # Connect to method that will handle changes
      #
      m_method_combo_box.connect('changed',
                                 self.classifyTwoStepChangeMatchMethod)

      classifier_param_box1.pack_start(m_value_label, False, False, 0)
      classifier_param_box1.pack_start(m_value_text_entry, False, False, 0)
      classifier_param_box1.pack_start(m_method_label, False, False, 0)
      classifier_param_box1.pack_start(m_method_combo_box, False, False, 0)
      m_value_label.show()
      m_value_text_entry.show()
      m_method_label.show()
      m_method_combo_box.show()

      # Following parameters depend upon match method
      #
      if (m_method_threshold == True):
        m_thres_label = gtk.Label('  Selection threshold:')
        m_thres_text_entry = gtk.Entry()
        m_thres_text_entry.set_width_chars(5)
        if ('s1_match_method' in classifier_dict):
          m_thres_text_entry.set_text(classifier_dict['s1_match_method'][2])
        classifier_param_box1.pack_start(m_thres_label, False, False, 0)
        classifier_param_box1.pack_start(m_thres_text_entry, False, False, 0)
        m_thres_label.show()
        m_thres_text_entry.show()
      else:  # Nearest
        m_nearest_label = gtk.Label('  Selection nearest:')
        m_nearest_text_entry = gtk.Entry()
        m_nearest_text_entry.set_width_chars(5)
        if ('s1_match_method' in classifier_dict):
          m_nearest_text_entry.set_text(classifier_dict['s1_match_method'][2])
        m_nearest_space = gtk.Label('  ')
        m_unique_check_box = gtk.CheckButton('Select unique weight vectors')
        if ('s1_match_method' in classifier_dict):
          m_unique_check_box.set_active(classifier_dict['s1_match_method'][3])
        classifier_param_box1.pack_start(m_nearest_label, False, False, 0)
        classifier_param_box1.pack_start(m_nearest_text_entry, False, False, 0)
        classifier_param_box1.pack_start(m_nearest_space)
        classifier_param_box1.pack_start(m_unique_check_box, False, False, 0)
        m_nearest_label.show()
        m_nearest_text_entry.show()
        m_nearest_space.show()
        m_unique_check_box.show()

      # Now the same for non-match method - - - - - - - - - - - - - - - - - - -
      #
      nm_value_label = gtk.Label('     Non-match comparison value:')
      nm_value_text_entry = gtk.Entry()
      nm_value_text_entry.set_width_chars(5)
      if ('s1_non_match_method' in classifier_dict):
        nm_value_text_entry.set_text(classifier_dict['s1_non_match_method'][0])
      else:
        nm_value_text_entry.set_text('0.0')  # Default value

      nm_method_label = gtk.Label('  Non-match method:')
      nm_method_combo_box = gtk.combo_box_new_text()
      nm_method_combo_box.append_text('Nearest')
      nm_method_combo_box.append_text('Threshold')
      if (('s1_non_match_method' in classifier_dict) and \
          (classifier_dict['s1_non_match_method'][1] == 'threshold')):
        nm_method_combo_box.set_active(1)
        nm_method_threshold = True
      else:
        nm_method_combo_box.set_active(0)  # Set as 'Nearest' per default
        nm_method_threshold = False

      # Connect to method that will handle changes
      #
      nm_method_combo_box.connect('changed',
                                  self.classifyTwoStepChangeNonMatchMethod)

      classifier_param_box2.pack_start(nm_value_label, False, False, 0)
      classifier_param_box2.pack_start(nm_value_text_entry, False, False, 0)
      classifier_param_box2.pack_start(nm_method_label, False, False, 0)
      classifier_param_box2.pack_start(nm_method_combo_box, False, False, 0)
      nm_value_label.show()
      nm_value_text_entry.show()
      nm_method_label.show()
      nm_method_combo_box.show()

      # Following parameters depend upon match method
      #
      if (nm_method_threshold == True):
        nm_thres_label = gtk.Label('  Selection threshold:')
        nm_thres_text_entry = gtk.Entry()
        nm_thres_text_entry.set_width_chars(5)
        if ('s1_non_match_method' in classifier_dict):
          nm_thres_text_entry.set_text( \
                                     classifier_dict['s1_non_match_method'][2])
        classifier_param_box2.pack_start(nm_thres_label, False, False, 0)
        classifier_param_box2.pack_start(nm_thres_text_entry, False, False, 0)
        nm_thres_label.show()
        nm_thres_text_entry.show()
      else:  # Nearest
        nm_nearest_label = gtk.Label('  Selection nearest:')
        nm_nearest_text_entry = gtk.Entry()
        nm_nearest_text_entry.set_width_chars(5)
        if ('s1_non_match_method' in classifier_dict):
          nm_nearest_text_entry.set_text( \
                                     classifier_dict['s1_non_match_method'][2])
        nm_nearest_space = gtk.Label('  ')
        nm_unique_check_box = gtk.CheckButton('Select unique weight vectors')
        if ('s1_non_match_method' in classifier_dict):
          nm_unique_check_box.set_active( \
                                     classifier_dict['s1_non_match_method'][3])
        classifier_param_box2.pack_start(nm_nearest_label, False, False, 0)
        classifier_param_box2.pack_start(nm_nearest_text_entry, False, False,0)
        classifier_param_box2.pack_start(nm_nearest_space)
        classifier_param_box2.pack_start(nm_unique_check_box, False, False, 0)
        nm_nearest_label.show()
        nm_nearest_text_entry.show()
        nm_nearest_space.show()
        nm_unique_check_box.show()

      # Random selection method - - - - - - - - - - - - - - - - - - - - - - - -
      #
      rand_sel_label = gtk.Label('     Random selection method:')
      rand_sel_combo_box = gtk.combo_box_new_text()
      rand_sel_combo_box.append_text('(None)')
      rand_sel_combo_box.append_text('Uniform')
      rand_sel_combo_box.append_text('Linear')
      rand_sel_combo_box.append_text('Exponential')
      rand_sel_combo_box.set_active(0)  # Default: None
      if (('random_selection' in classifier_dict) and \
          (classifier_dict['random_selection'] != None)):
        if (classifier_dict['random_selection'][0] == 'uniform'):
          rand_sel_combo_box.set_active(1)
        elif (classifier_dict['random_selection'][0] == 'linear'):
          rand_sel_combo_box.set_active(2)
        elif (classifier_dict['random_selection'][0] == 'exponential'):
          rand_sel_combo_box.set_active(3)

      # Connect to method that will handle changes
      #
      rand_sel_combo_box.connect('changed', self.classifyTwoStepChangeRandSel)

      classifier_param_box3.pack_start(rand_sel_label, False, False, 0)
      classifier_param_box3.pack_start(rand_sel_combo_box, False, False, 0)
      rand_sel_label.show()
      rand_sel_combo_box.show()

      # Add parameter boxes if not None
      #
      if (rand_sel_combo_box.get_active() != 0):
        rand_m_perc_label = gtk.Label('  Percentage add to match set:')
        rand_m_perc_text_entry = gtk.Entry()
        rand_m_perc_text_entry.set_width_chars(5)
        rand_m_perc_text_entry.set_text(classifier_dict['random_selection'][1])
        rand_nm_perc_label = gtk.Label('  Percentage add to non-match set:')
        rand_nm_perc_text_entry = gtk.Entry()
        rand_nm_perc_text_entry.set_width_chars(5)
        rand_nm_perc_text_entry.set_text( \
                                        classifier_dict['random_selection'][2])
        classifier_param_box3.pack_start(rand_m_perc_label, False, False, 0)
        classifier_param_box3.pack_start(rand_m_perc_text_entry, False,False,0)
        classifier_param_box3.pack_start(rand_nm_perc_label, False, False, 0)
        classifier_param_box3.pack_start(rand_nm_perc_text_entry,False,False,0)
        rand_m_perc_label.show()
        rand_m_perc_text_entry.show()
        rand_nm_perc_label.show()
        rand_nm_perc_text_entry.show()

      # Step two classifier - - - - - - - - - - - - - - - - - - - - - - - - - -
      #
      classifier_label = gtk.Label('     Step two classification method:')
      classifier_combo_box = gtk.combo_box_new_text()
      classifier_combo_box.append_text('SVM')
      classifier_combo_box.append_text('K-means')
      classifier_combo_box.set_active(0)  # Default: SVM
      if (('s2_classifier' in classifier_dict) and \
          (classifier_dict['s2_classifier'][0] == 'kmeans')):
        classifier_combo_box.set_active(1)

      # Connect to method that will handle changes
      #
      classifier_combo_box.connect('changed', \
                                          self.classifyTwoStepChangeClassifier)

      classifier_param_box4.pack_start(classifier_label, False, False, 0)
      classifier_param_box4.pack_start(classifier_combo_box, False, False, 0)
      classifier_label.show()
      classifier_combo_box.show()

      if (classifier_combo_box.get_active() == 0):  # SVM classifier
        kernel_type_label = gtk.Label('  Kernel type:')
        kernel_type_combo_box = gtk.combo_box_new_text()
        kernel_type_combo_box.append_text('Linear')
        kernel_type_combo_box.append_text('Poly')
        kernel_type_combo_box.append_text('RBF')
        kernel_type_combo_box.append_text('Sigmoid')
        kernel_type_combo_box.set_active(0)  # Default: LINEAR
        if ('s2_classifier' in classifier_dict):
          if (classifier_dict['s2_classifier'][1] == 'POLY'):
            kernel_type_combo_box.set_active(1)
          elif (classifier_dict['s2_classifier'][1] == 'RBF'):
            kernel_type_combo_box.set_active(2)
          elif (classifier_dict['s2_classifier'][1] == 'SIGMOID'):
            kernel_type_combo_box.set_active(3)
        C_label = gtk.Label('  C:')
        C_text_entry = gtk.Entry()
        C_text_entry.set_width_chars(10)
        if ('s2_classifier' in classifier_dict):
          C_text_entry.set_text(classifier_dict['s2_classifier'][2])
        classifier_param_box4.pack_start(kernel_type_label, False, False, 0)
        classifier_param_box4.pack_start(kernel_type_combo_box, False, False,0)
        classifier_param_box4.pack_start(C_label, False, False, 0)
        classifier_param_box4.pack_start(C_text_entry, False, False, 0)
        kernel_type_combo_box.show()
        kernel_type_label.show()
        C_label.show()
        C_text_entry.show()

      elif (classifier_combo_box.get_active() == 1):  # K-means classifier
        dist_measure_label = gtk.Label('  Distance measure:')
        dist_measure_combo_box = gtk.combo_box_new_text()
        dist_measure_combo_box.append_text('Manhatten')
        dist_measure_combo_box.append_text('Euclidean')
        dist_measure_combo_box.append_text('L-Infinity')
        dist_measure_combo_box.append_text('Canberra')
        dist_measure_combo_box.set_active(1)  # Default: Euclidean
        if ('s2_classifier' in classifier_dict):
          if (classifier_dict['s2_classifier'][1] == 'Manhatten'):
            dist_measure_combo_box.set_active(0)
          elif (classifier_dict['s2_classifier'][1] == 'L-Infinity'):
            dist_measure_combo_box.set_active(2)
          elif (classifier_dict['s2_classifier'][1] == 'Canberra'):
            dist_measure_combo_box.set_active(3)
        classifier_param_box4.pack_start(dist_measure_label, False, False, 0)
        classifier_param_box4.pack_start(dist_measure_combo_box, False,False,0)
        dist_measure_label.show()
        dist_measure_combo_box.show()

    # For supervised classification methods match status is needed - - - - - -
    # (assumed to be the exact comparisons of a record attribute done in
    # comparison step)
    #
    if (classifier_method_name in ['OptimalThreshold', 'SuppVecMachine']):

      # Get all exact comparisons, their index in the comparison list (to get
      # index in weight vector), and name of fields involved.
      #
      self.exact_comp_list = []

      i = 0  # Index in field comparison list
      for field_comp_dict in self.field_comp_list:

        if (field_comp_dict['name'] == 'Str-Exact'):  # Exact string comparison
          field_name1 = field_comp_dict['field_a_name']
          field_name2 = field_comp_dict['field_b_name']

          # Only one field name required for deduplication
          #
          if ((field_name1 == field_name2) and (self.project_type != 'Link')):
            match_status_name = 'Exact match on comparison of field "%s"' % \
                                (field_name1)+' with itself'
          else:  # Two different field names or linkage
            match_status_name = 'Exact match on comparison of fields ' + \
                                '"%s" with "%s"' % (field_name1, field_name2)
          self.exact_comp_list.append((match_status_name, i))

        i += 1

      match_status_label = gtk.Label('     Determine match status:')
      match_status_combo_box = gtk.combo_box_new_text()
      for (match_status_name, i) in self.exact_comp_list:
        match_status_combo_box.append_text(match_status_name)
      if ('match_status_ind' in classifier_dict):
        match_status_combo_box.set_active(classifier_dict['match_status_ind'])
      else:  # Set to first
        match_status_combo_box.set_active(0)

      classifier_param_box2.pack_start(match_status_label, False, False, 0)
      classifier_param_box2.pack_start(match_status_combo_box, False, False, 0)
      match_status_label.show()
      match_status_combo_box.show()

      # Warning if no exact comparison available - - - - - - - - - - - - - - -
      #
      if (len(self.exact_comp_list) == 0):
        self.messageDialog('No exact comparisons given, so match status\n' + \
                           'cannot be determined. Please add an exact\n' + \
                           'string comparison of the entity identifiers\n' + \
                           'on "Compare" page.', 'error')

    classification_widget.pack_start(classifier_param_box1, False, False, 0)
    classifier_param_box1.show()

    for hbox in [classifier_param_box2, classifier_param_box3,
                 classifier_param_box4]:
      if (len(hbox.get_children()) > 0):
        classification_widget.pack_start(hbox, False, False, 0)
        hbox.show()

  # ---------------------------------------------------------------------------
  # Handle a change of the classification method in combo box
  #
  def classifierChangeMethod(self, widget):
    print 'Changed classification method in combo box to:', widget.get_active()

    new_classifier_method = self.classifier_names[widget.get_active()]

    # Set method name
    #
    self.classifier_method['name'] = new_classifier_method

    self.classifyView()  # Re-display the classification page

  # ---------------------------------------------------------------------------
  # Handle a change of the match method for two-step classifier
  #
  def classifyTwoStepChangeMatchMethod(self, widget):
    print 'Changed two-step match method'

    classifier_dict = self.classifier_method  # For quicker access

    # Get the horizontal box containing the match parameter widgets
    #
    classification_widget = self.mainTree.get_widget('classify_page_box')
    classification_widget_children = classification_widget.get_children()
    classifier_param_box1 = classification_widget_children[1]

    # Remove all method specific parameters
    #
    for child in classifier_param_box1.get_children()[4:]:
      classifier_param_box1.remove(child)

    # Re-build according to selected method
    #
    if (widget.get_active() == 0):  # nearest
      m_nearest_label = gtk.Label('  Selection nearest:')
      m_nearest_text_entry = gtk.Entry()
      m_nearest_text_entry.set_width_chars(5)
      if ('s1_match_method' in classifier_dict):
        m_nearest_text_entry.set_text(classifier_dict['s1_match_method'][2])
      m_nearest_space = gtk.Label('  ')
      m_unique_check_box = gtk.CheckButton('Select unique weight vectors')
      if ('s1_match_method' in classifier_dict):
        m_unique_check_box.set_active(classifier_dict['s1_match_method'][3])
      classifier_param_box1.pack_start(m_nearest_label, False, False, 0)
      classifier_param_box1.pack_start(m_nearest_text_entry, False, False, 0)
      classifier_param_box1.pack_start(m_nearest_space)
      classifier_param_box1.pack_start(m_unique_check_box, False, False, 0)
      m_nearest_label.show()
      m_nearest_text_entry.show()
      m_nearest_space.show()
      m_unique_check_box.show()

    else:  # Threshold
      m_thres_label = gtk.Label('  Selection threshold:')
      m_thres_text_entry = gtk.Entry()
      m_thres_text_entry.set_width_chars(5)
      if ('s1_match_method' in classifier_dict):
        m_thres_text_entry.set_text(classifier_dict['s1_match_method'][2])
      classifier_param_box1.pack_start(m_thres_label, False, False, 0)
      classifier_param_box1.pack_start(m_thres_text_entry, False, False, 0)
      m_thres_label.show()
      m_thres_text_entry.show()

  # ---------------------------------------------------------------------------
  # Handle a change of the non-match method for two-step classifier
  #
  def classifyTwoStepChangeNonMatchMethod(self, widget):
    print 'Changed two-step non-match method'

    classifier_dict = self.classifier_method  # For quicker access

    # Get the horizontal box containing the non-match parameter widgets
    #
    classification_widget = self.mainTree.get_widget('classify_page_box')
    classification_widget_children = classification_widget.get_children()
    classifier_param_box2 = classification_widget_children[2]

    # Remove all method specific parameters
    #
    for child in classifier_param_box2.get_children()[4:]:
      classifier_param_box2.remove(child)

    # Re-build according to selected method
    #
    if (widget.get_active() == 0):  # nearest
      nm_nearest_label = gtk.Label('  Selection nearest:')
      nm_nearest_text_entry = gtk.Entry()
      nm_nearest_text_entry.set_width_chars(5)
      if ('s1_non_match_method' in classifier_dict):
        nm_nearest_text_entry.set_text( \
                                     classifier_dict['s1_non_match_method'][2])
      nm_nearest_space = gtk.Label('  ')
      nm_unique_check_box = gtk.CheckButton('Select unique weight vectors')
      if ('s1_non_match_method' in classifier_dict):
        nm_unique_check_box.set_active( \
                                     classifier_dict['s1_non_match_method'][3])
      classifier_param_box2.pack_start(nm_nearest_label, False, False, 0)
      classifier_param_box2.pack_start(nm_nearest_text_entry, False, False, 0)
      classifier_param_box2.pack_start(nm_nearest_space)
      classifier_param_box2.pack_start(nm_unique_check_box, False, False, 0)
      nm_nearest_label.show()
      nm_nearest_text_entry.show()
      nm_nearest_space.show()
      nm_unique_check_box.show()

    else:  # Threshold
      nm_thres_label = gtk.Label('  Selection threshold:')
      nm_thres_text_entry = gtk.Entry()
      nm_thres_text_entry.set_width_chars(5)
      if ('s1_non_match_method' in classifier_dict):
        nm_thres_text_entry.set_text(classifier_dict['s1_non_match_method'][2])
      classifier_param_box2.pack_start(nm_thres_label, False, False, 0)
      classifier_param_box2.pack_start(nm_thres_text_entry, False, False, 0)
      nm_thres_label.show()
      nm_thres_text_entry.show()

  # ---------------------------------------------------------------------------
  # Handle a change of the random selection method for two-step classifier
  #
  def classifyTwoStepChangeRandSel(self, widget):
    print 'Changed two-step random selection method'

    classifier_dict = self.classifier_method  # For quicker access

    # Get the horizontal box containing the random selection parameter widgets
    #
    classification_widget = self.mainTree.get_widget('classify_page_box')
    classification_widget_children = classification_widget.get_children()
    classifier_param_box3 = classification_widget_children[3]

    # Check if random selection is None now and was not None before
    #
    if (widget.get_active() == 0) and (len(classifier_param_box3) == 6):

      # Remove parameter widgets
      #
      for child in classifier_param_box3.get_children()[2:]:
        classifier_param_box3.remove(child)  # Remove parameter widgets

    # Check if random selection is not None now but was None before
    #
    elif (widget.get_active() != 0) and (len(classifier_param_box3) < 6):

      # Add parameter widgets
      #
      rand_m_perc_label = gtk.Label('  Percentage add to match set:')
      rand_m_perc_text_entry = gtk.Entry()
      rand_m_perc_text_entry.set_width_chars(5)
      rand_nm_perc_label = gtk.Label('  Percentage add to non-match set:')
      rand_nm_perc_text_entry = gtk.Entry()
      rand_nm_perc_text_entry.set_width_chars(5)
      if (('random_selection' in classifier_dict) and \
          (classifier_dict['random_selection'] != None)):
        rand_m_perc_text_entry.set_text(classifier_dict['random_selection'][1])
        rand_nm_perc_text_entry.set_text( \
                                        classifier_dict['random_selection'][2])
      classifier_param_box3.pack_start(rand_m_perc_label, False, False, 0)
      classifier_param_box3.pack_start(rand_m_perc_text_entry, False,False,0)
      classifier_param_box3.pack_start(rand_nm_perc_label, False, False, 0)
      classifier_param_box3.pack_start(rand_nm_perc_text_entry,False,False,0)
      rand_m_perc_label.show()
      rand_m_perc_text_entry.show()
      rand_nm_perc_label.show()
      rand_nm_perc_text_entry.show()

  # ---------------------------------------------------------------------------
  # Handle a change of the random selection method for two-step classifier
  #
  def classifyTwoStepChangeClassifier(self, widget):
    print 'Changed two-step classifier method'

    classifier_dict = self.classifier_method  # For quicker access

    # Get the horizontal box containing the classifier widgets
    #
    classification_widget = self.mainTree.get_widget('classify_page_box')
    classification_widget_children = classification_widget.get_children()
    classifier_param_box4 = classification_widget_children[4]

    # Remove parameter widgets
    #
    for child in classifier_param_box4.get_children()[2:]:
      classifier_param_box4.remove(child)  # Remove parameter widgets

    if (widget.get_active() == 0):  # SVM classifier
      kernel_type_label = gtk.Label('  Kernel type:')
      kernel_type_combo_box = gtk.combo_box_new_text()
      kernel_type_combo_box.append_text('Linear')
      kernel_type_combo_box.append_text('Poly')
      kernel_type_combo_box.append_text('RBF')
      kernel_type_combo_box.append_text('Sigmoid')
      kernel_type_combo_box.set_active(0)  # Default: LINEAR
      C_label = gtk.Label('  C:')
      C_text_entry = gtk.Entry()
      C_text_entry.set_width_chars(10)
      classifier_param_box4.pack_start(kernel_type_label, False, False, 0)
      classifier_param_box4.pack_start(kernel_type_combo_box, False, False,0)
      classifier_param_box4.pack_start(C_label, False, False, 0)
      classifier_param_box4.pack_start(C_text_entry, False, False, 0)
      kernel_type_combo_box.show()
      kernel_type_label.show()
      C_label.show()
      C_text_entry.show()

    elif (widget.get_active() == 1):  # K-means classifier
      dist_measure_label = gtk.Label('  Distance measure:')
      dist_measure_combo_box = gtk.combo_box_new_text()
      dist_measure_combo_box.append_text('Manhatten')
      dist_measure_combo_box.append_text('Euclidean')
      dist_measure_combo_box.append_text('L-Infinity')
      dist_measure_combo_box.append_text('Canberra')
      dist_measure_combo_box.set_active(1)  # Default: Euclidean
      classifier_param_box4.pack_start(dist_measure_label, False, False, 0)
      classifier_param_box4.pack_start(dist_measure_combo_box, False, False, 0)
      dist_measure_label.show()
      dist_measure_combo_box.show()

  # ---------------------------------------------------------------------------
  # Handle an activate of Execute on the Classify page
  #
  def classifyExecute(self):
    classification_widget = self.mainTree.get_widget('classify_page_box')
    classifier_method_widget = \
                        self.mainTree.get_widget('classifier_method_combo_box')

    self.febrl_code['classification'] = None  # Remove all previous class. code

    # Get all children which will contain classifier parameters
    #
    classification_widget_children = classification_widget.get_children()

    classifier_dict = self.classifier_method  # For quicker access

    classifier_name_box = classification_widget_children[0]
    classifier_name_index = classifier_name_box.get_children()[1].get_active()

    if (classifier_name_index == -1):  # No classifier selected - - - - - - - -
      self.messageDialog('No classifier selected.', 'error')
      return

    classifier_param_list1 = classification_widget_children[1].get_children()
    if (len(classification_widget_children) == 3):
      classifier_param_list2 = classification_widget_children[2].get_children()
    if (len(classification_widget_children) == 5):  # For two-step classifier
      classifier_param_list2 = classification_widget_children[2].get_children()
      classifier_param_list3 = classification_widget_children[3].get_children()
      classifier_param_list4 = classification_widget_children[4].get_children()

    classifier_method_name = self.classifier_names[classifier_name_index]

    # First get the match status information for supervised classifiers - - - -
    #
    if (classifier_method_name in ['OptimalThreshold', 'SuppVecMachine']):

      # Warning if no exact comparison available - - - - - - - - - - - - - - -
      #
      if (len(self.exact_comp_list) == 0):
        self.messageDialog('No exact comparisons given, so match status\n' + \
                           'cannot be determined. Please add an exact\n' + \
                           'string comparison of the entity identifiers\n' + \
                           'on "Compare" page.', 'error')
        return

    if (classifier_method_name == 'FellegiSunter'):  # - - - - - - - - - - - -

      lower_threshold = classifier_param_list1[1].get_text().strip()
      upper_threshold = classifier_param_list1[3].get_text().strip()

      if ((not self.str_is_float(lower_threshold)) or \
          (not self.str_is_float(upper_threshold))):
        self.messageDialog('Both thresholds must be numbers.', 'warn')
        return

      if (float(lower_threshold) > float(upper_threshold)):
        self.messageDialog('Lower threshold must be smaller or\n' + \
                           'equal to upper threshold.', 'warn')
        return

      if (lower_threshold[0] == '.'):
        lower_threshold = '0'+lower_threshold
      if (upper_threshold[0] == '.'):
        upper_threshold = '0'+upper_threshold
      classifier_dict['lower_threshold'] = lower_threshold
      classifier_dict['upper_threshold'] = upper_threshold

    elif (classifier_method_name == 'OptimalThreshold'):  # - - - - - - - - - -

      bin_width_val = classifier_param_list1[3].get_text().strip()
      if ((not self.str_is_float(bin_width_val)) or \
          (float(bin_width_val) == 0)):
        self.messageDialog('Bin width must be a non-zero number.', 'warn')
        return

      if (bin_width_val[0] == '.'):
        bin_width_val = '0'+bin_width_val
      classifier_dict['bin_width'] = bin_width_val

      min_method_index = classifier_param_list1[1].get_active()
      classifier_dict['min_method'] = ['pos-neg', 'pos',
                                       'neg'][min_method_index]

      match_status_index = classifier_param_list2[1].get_active()
      classifier_dict['match_status_ind'] = match_status_index

    elif (classifier_method_name in ['KMeans', 'FarthestFirst']):  # - - - - -

      sample_val = classifier_param_list1[3].get_text().strip()
      if (sample_val != ''):
        if (not self.str_is_percentage_not_zero(sample_val)):
          self.messageDialog('Sample value can be empty (not applied)\n' + \
                             'or it must be a positive percentage value.',
                             'warn')
          return

      fuzzy_thres_val = classifier_param_list2[3].get_text().strip()
      if (fuzzy_thres_val != ''):
        if (not self.str_is_normalised(fuzzy_thres_val)):
          self.messageDialog('Fuzzy threshold can be empty (not applied)\n' + \
                             'or its value must be larger or equal\n' + \
                             'to 0 and smaller or equal to 1.','warn')
          return
        if (fuzzy_thres_val[0] == '.'):
         fuzzy_thres_val = '0'+fuzzy_thres_val

      if (classifier_method_name == 'KMeans'):
        max_iter_count_val = classifier_param_list1[5].get_text().strip()
        if (not self.str_is_not_neg_int(max_iter_count_val)):
          self.messageDialog('Maximum iteration count must be\n' + \
                             'a positive integer number.', 'warn')
          return
        classifier_dict['max_iter_count'] = max_iter_count_val

      dist_measure_index = classifier_param_list1[1].get_active()
      classifier_dict['dist_measure'] = ['Manhatten','Euclidean','L-Infinity',
                                         'Canberra'][dist_measure_index]

      centroid_init_index = classifier_param_list2[1].get_active()
      if (classifier_method_name == 'KMeans'):
        classifier_dict['centroid_init'] = ['min/max',
                                            'random'][centroid_init_index]
      else:  # Farthest first
        classifier_dict['centroid_init'] = ['traditional','min/max',
                                            'mode/max'][centroid_init_index]

      if (sample_val != ''):
        sample_val = str(int(float(sample_val)))
      classifier_dict['sample'] = sample_val
      classifier_dict['fuzzy_reg_thres'] = fuzzy_thres_val

    elif (classifier_method_name == 'SuppVecMachine'):  # - - - - - - - - - - -

      C_val = classifier_param_list1[3].get_text().strip()
      if ((not self.str_is_float(C_val)) or (float(C_val) == 0)):
        self.messageDialog('Value of C must be a non-zero number.', 'warn')
        return

      sample_val = classifier_param_list1[5].get_text().strip()
      if (sample_val != ''):
        if (not self.str_is_percentage_not_zero(sample_val)):
          self.messageDialog('Sample value can be empty (not applied)\n' + \
                             'or it must be a positive percentage value.',
                             'warn')
          return

      kernel_type_index = classifier_param_list1[1].get_active()
      classifier_dict['kernel_type'] = ['LINEAR','POLY','RBF',
                                        'SIGMOID'][kernel_type_index]

      if (C_val[0] == '.'):
        C_val = '0'+C_val
      classifier_dict['C'] = C_val

      if (sample_val != ''):
        sample_val = str(int(float(sample_val)))
      classifier_dict['sample'] = sample_val

      match_status_index = classifier_param_list2[1].get_active()
      classifier_dict['match_status_ind'] = match_status_index

    elif (classifier_method_name == 'TwoStep'):  # - - - - - - - - - - - - - -

      m_method_val = classifier_param_list1[1].get_text().strip()
      if (not self.str_is_float(m_method_val)):
        self.messageDialog('Match comparison value must be a number.', 'warn')
        return

      nm_method_val = classifier_param_list2[1].get_text().strip()
      if (not self.str_is_float(nm_method_val)):
        self.messageDialog('Non-match comparison value must be a number.',
                           'warn')
        return

      if (float(m_method_val) < float(nm_method_val)):
        self.messageDialog('Match comparison value must be larger than\n' + \
                           'non-match comparison value.', 'warn')
        return

      if (m_method_val[0] == '.'):
        m_method_val = '0'+m_method_val
      if (nm_method_val[0] == '.'):
        nm_method_val = '0'+nm_method_val

      m_method_index =  classifier_param_list1[3].get_active()
      m_method_type =   ['nearest','threshold'][m_method_index]
      nm_method_index = classifier_param_list2[3].get_active()
      nm_method_type =  ['nearest','threshold'][nm_method_index]

      if (m_method_type == 'nearest'):
        m_nearest_val = classifier_param_list1[5].get_text().strip()
        if (not self.str_is_pos_int(m_nearest_val)):
          self.messageDialog('Match selection nearest must be a\n' + \
                             'positive integer value.', 'warn')
          return
        m_unique_val = classifier_param_list1[7].get_active()

        match_method_list = [m_method_val, m_method_type, m_nearest_val,
                             m_unique_val]

      else:  # Threshold method
        m_threshold_val = classifier_param_list1[5].get_text().strip()
        if (not self.str_is_normalised_not_zero(m_threshold_val)):
          self.messageDialog('Match selection threshold value must be\n' + \
                             'larger than 0 and smaller or equal to 1.',
                             'warn')
          return
        if (m_threshold_val[0] == '.'):
          m_threshold_val = '0'+m_threshold_val

        match_method_list = [m_method_val, m_method_type, m_threshold_val]

      if (nm_method_type == 'nearest'):
        nm_nearest_val = classifier_param_list2[5].get_text().strip()
        if (not self.str_is_pos_int(nm_nearest_val)):
          self.messageDialog('Non-match selection nearest must be a\n' + \
                             'positive integer value.', 'warn')
          return
        nm_unique_val = classifier_param_list2[7].get_active()

        non_match_method_list = [nm_method_val, nm_method_type, nm_nearest_val,
                                 nm_unique_val]

      else:  # Threshold method
        nm_threshold_val = classifier_param_list2[5].get_text().strip()
        if (not self.str_is_normalised_not_zero(nm_threshold_val)):
          self.messageDialog('Non-match selection threshold value must be' + \
                             '\nlarger than 0 and smaller or equal to 1.',
                             'warn')
          return
        if (nm_threshold_val[0] == '.'):
          nm_threshold_val = '0'+nm_threshold_val

        non_match_method_list = [nm_method_val, nm_method_type,
                                 nm_threshold_val]

      # Random selection parameters
      #
      rand_sel_index = classifier_param_list3[1].get_active()
      rand_sel_type =  [None,'uniform','linear','exponential'][rand_sel_index]

      if (rand_sel_type == None):
        rand_selection_list = None
      else:
        m_perc_val =  classifier_param_list3[3].get_text().strip()
        nm_perc_val = classifier_param_list3[5].get_text().strip()

        if ((not self.str_is_percentage(m_perc_val)) or \
            (not self.str_is_percentage(nm_perc_val))):
          self.messageDialog('Random selection percentage values have\n' + \
                             'to be percentages (0..100).', 'warn')
          return
        if ((float(m_perc_val) + float(nm_perc_val)) > 80):
          self.messageDialog('Sum of random selection percentage values\n' + \
                             'must be less than or equal to 80%.', 'warn')
          return

        rand_selection_list = [rand_sel_type, m_perc_val, nm_perc_val ]

      # Step two classifier parameters
      #
      s2_classifier_index = classifier_param_list4[1].get_active()
      s2_classifier_type =  ['svm','kmeans'][s2_classifier_index]

      if (s2_classifier_type == 'svm'):
        kernel_type_index = classifier_param_list4[3].get_active()
        kernel_type_val = ['LINEAR','POLY','RBF','SIGMOID'][kernel_type_index]
        C_val = classifier_param_list4[5].get_text().strip()
        if ((not self.str_is_float(C_val)) or (float(C_val) == 0)):
          self.messageDialog('Value of C must be a non-zero number.', 'warn')
          return
        if (C_val[0] == '.'):
          C_val = '0'+C_val
        classifier_list = [s2_classifier_type, kernel_type_val, C_val]

      elif (s2_classifier_type == 'kmeans'):
        dist_measure_index = classifier_param_list4[3].get_active()
        dist_measure_val = ['Manhatten','Euclidean','L-Infinity',
                            'Canberra'][dist_measure_index]
        classifier_list = [s2_classifier_type, dist_measure_val]

      classifier_dict['s1_match_method'] =     match_method_list
      classifier_dict['s1_non_match_method'] = non_match_method_list
      classifier_dict['random_selection'] =    rand_selection_list
      classifier_dict['s2_classifier'] =       classifier_list

    # Now generate the Febrl codes for classification - - - - - - - - - - - - -
    #
    febrl_code = []
    febrl_code.append('# '+'-'*77)
    febrl_code.append('')
    febrl_code.append('# Define weight vector (record pair) classifier')
    febrl_code.append('#')

    if (classifier_method_name == 'FellegiSunter'):  # - - - - - - - - - - - -
      indention_space = ' '*42

      febrl_code.append('classifier = classification.FellegiSunter(' + \
                        'lower_threshold = %s,' % \
                        (classifier_dict['lower_threshold']))
      febrl_code.append(indention_space+'upper_threshold = %s)' % \
                        (classifier_dict['upper_threshold']))

    elif (classifier_method_name == 'OptimalThreshold'):  # - - - - - - - - - -
      indention_space = ' '*45

      febrl_code.append('classifier = classification.OptimalThreshold(' \
                        + 'bin_width = %s,' % (classifier_dict['bin_width']))
      febrl_code.append(indention_space+'min_method = "%s")' % \
                        (classifier_dict['min_method']))

    elif (classifier_method_name == 'KMeans'):  # - - - - - - - - - - - - - - -
      indention_space = ' '*35

      if (classifier_dict['dist_measure'] == 'Manhatten'):
        febrl_dist_str = 'mymath.distL1'
      elif (classifier_dict['dist_measure'] == 'Euclidean'):
        febrl_dist_str = 'mymath.distL2'
      elif (classifier_dict['dist_measure'] == 'L-Infinity'):
        febrl_dist_str = 'mymath.distLInf'
      else:
        febrl_dist_str = 'mymath.distCanberra'

      febrl_code.append('classifier = classification.KMeans(' + \
                        'dist_measure = %s,' % (febrl_dist_str))
      febrl_code.append(indention_space+'max_iter_count = %s,' % \
                        (classifier_dict['max_iter_count']))
      if (('sample' in classifier_dict) and (classifier_dict['sample'] != '')):
        febrl_code.append(indention_space+'sample = %s,' % \
                          (classifier_dict['sample']))
      if (('fuzzy_reg_thres' in classifier_dict) and \
          (classifier_dict['fuzzy_reg_thres'] != '')):
        febrl_code.append(indention_space+'fuzz_reg_thres = %s,' % \
                          (classifier_dict['fuzzy_reg_thres']))
      febrl_code.append(indention_space+'centroid_init = "%s")' % \
                        (classifier_dict['centroid_init']))

    elif (classifier_method_name == 'FarthestFirst'):  # - - - - - - - - - - -
      indention_space = ' '*42

      if (classifier_dict['dist_measure'] == 'Manhatten'):
        febrl_dist_str = 'mymath.distL1'
      elif (classifier_dict['dist_measure'] == 'Euclidean'):
        febrl_dist_str = 'mymath.distL2'
      elif (classifier_dict['dist_measure'] == 'L-Infinity'):
        febrl_dist_str = 'mymath.distLInf'
      else:
        febrl_dist_str = 'mymath.distCanberra'

      febrl_code.append('classifier = classification.FarthestFirst(' + \
                        'dist_measure = %s,' % (febrl_dist_str))
      if (('sample' in classifier_dict) and (classifier_dict['sample'] != '')):
        febrl_code.append(indention_space+'sample = %s,' % \
                          (classifier_dict['sample']))
      if (('fuzzy_reg_thres' in classifier_dict) and \
          (classifier_dict['fuzzy_reg_thres'] != '')):
        febrl_code.append(indention_space+'fuzz_reg_thres = %s,' % \
                          (classifier_dict['fuzzy_reg_thres']))
      febrl_code.append(indention_space+'centroid_init = "%s")' % \
                        (classifier_dict['centroid_init']))

    elif (classifier_method_name == 'SuppVecMachine'):  # - - - - - - - - - - -
      indention_space = ' '*43

      febrl_code.append('classifier = classification.SuppVecMachine(' + \
                        'kernel_type = "%s",' % \
                        (classifier_dict['kernel_type']))
      if (('sample' in classifier_dict) and (classifier_dict['sample'] != '')):
        febrl_code.append(indention_space+'sample = %s,' % \
                          (classifier_dict['sample']))
      febrl_code.append(indention_space+'C = %s)' % \
                        (classifier_dict['C']))

    elif (classifier_method_name == 'TwoStep'):  # - - - - - - - - - - - - - -
      indention_space = ' '*36

      if (classifier_dict['s1_match_method'][1] == 'nearest'):
        s1_method_str = '(%s, "nearest", %s, %s)' % \
                        (classifier_dict['s1_match_method'][0],
                         classifier_dict['s1_match_method'][2],
                         classifier_dict['s1_match_method'][3])
      else:  # Threshold
        s1_method_str = '(%s, "threshold", %s)' % \
                        (classifier_dict['s1_match_method'][0],
                         classifier_dict['s1_match_method'][2])

      if (classifier_dict['s1_non_match_method'][1] == 'nearest'):
        s1_non_method_str = '(%s, "nearest", %s, %s)' % \
                            (classifier_dict['s1_non_match_method'][0],
                             classifier_dict['s1_non_match_method'][2],
                             classifier_dict['s1_non_match_method'][3])
      else:  # Threshold
        s1_non_method_str = '(%s, "threshold", %s)' % \
                            (classifier_dict['s1_non_match_method'][0],
                             classifier_dict['s1_non_match_method'][2])

      if (classifier_dict['random_selection'] == None):
        rand_sel_str = None
      else:
        rand_sel_str = '("%s", %s, %s)' % \
                       (classifier_dict['random_selection'][0],
                        classifier_dict['random_selection'][1],
                        classifier_dict['random_selection'][2])

      if (classifier_dict['s2_classifier'][0] == 'svm'):
        s2_classifier_str = '("svm", "%s", %s)' % \
                             (classifier_dict['s2_classifier'][1],
                              classifier_dict['s2_classifier'][2])
      else:  # K-means
        dist_method = classifier_dict['s2_classifier'][1]

        if (dist_method == 'Manhatten'):
          febrl_dist_str = 'mymath.distL1'
        elif (dist_method == 'Euclidean'):
          febrl_dist_str = 'mymath.distL2'
        elif (dist_method == 'L-Infinity'):
          febrl_dist_str = 'mymath.distLInf'
        else:
          febrl_dist_str = 'mymath.distCanberra'

        s2_classifier_str = '("kmeans", %s)' % (febrl_dist_str)

      febrl_code.append('classifier = classification.TwoStep(' + \
                        's1_match_method = %s,' % (s1_method_str))
      febrl_code.append(indention_space+'s1_non_match_method = %s,' % \
                        (s1_non_method_str))

      if (rand_sel_str != None):
        febrl_code.append(indention_space+'random_selection = %s,' % \
                          (rand_sel_str))

      febrl_code.append(indention_space+'s2_classifier = %s)' % \
                        (s2_classifier_str))

    febrl_code.append('')

    self.febrl_code['classification'] = febrl_code  # Store for later use

    # Finally update the GUI information - - - - - - - - - - - - - - - - - - -
    #
    self.addToLog('')  # Add generated code into log page text
    self.addToLog('='*79)
    self.addToLog('Generated Febrl code for "classification" on %s' % \
                  (time.asctime()))
    self.addToLog('')
    for line in febrl_code:
      self.addToLog(line)
    self.addToLog('')

    self.modified['classification'] = True  # Details have been changed
    self.setWindowTitle()

    # Update the active and non-active notebook pages
    #
    self.main_notebook_page_active_dict['Evaluate'] = False
    self.main_notebook_page_active_dict['Run'] = True

    self.displayCurrentNotebookPage()  # Diplay the current page

    self.evaluation_dict['available'] = False  # So evaluation is re-calculated

    self.writeStatusBar('Generated Febrl Python code for classification' + \
                        ' (see Log page for generated code).')

  # ===========================================================================
  # Methods that handle Output/Run page events
  # ===========================================================================

  # ---------------------------------------------------------------------------
  # Display the Output/Run page according to if output files have been
  # initialised or not
  #
  def runView(self):
    print '  Switched to Output/Run page - Display this page'

    output_dict = self.output_dict  # Shortcut
    print output_dict

    # Switch to note book page according to project type
    #
    run_nb_widget = self.mainTree.get_widget('run_page_notebook')

    if (self.project_type == 'Standardise'):  # - - - - - - - - - - - - - - - -
      run_nb_widget.set_current_page(1)

      # Set progress percentage value
      #
      progress_perc_text_entry = \
                            self.mainTree.get_widget('run_std_progress_entry')
      if (output_dict['progress_perc'] != None):
        progress_perc_text_entry.set_text(output_dict['progress_perc'])
      else:
        progress_perc_text_entry.set_text('None')

      # If no output data set file name given generate one from input name - -
      #
      out_file_name = output_dict['std_out_file']

      if (out_file_name == '(None)'):
        in_file_name = self.data_set_info_list[0]['file_name']
        dot_index = in_file_name.rfind('.')
        if (dot_index > 0):
          out_file_name = in_file_name[:dot_index] + '-standardised' + \
                          in_file_name[dot_index:]
        else:
          out_file_name = in_file_name + '-standardised'
        output_dict['std_out_file'] = out_file_name  # Save generated file name

      out_file_button = self.mainTree.get_widget('run_std_out_file_button')
      out_file_button.set_label(out_file_name.split(os.sep)[-1])

      field_name_list = self.data_set_info_list[0]['field_names']

      # Generate table with three columns for pass field check button
      #
      num_rows = int(math.ceil(float(len(field_name_list))/3))

      pass_field_table = gtk.Table(num_rows, 3)

      # List of previously clicked (True) fields (to be passed to output)
      #
      pass_field_clicked_list = output_dict['pass_field_list']

      col_ind = 0
      row_num = 0
      for field_name in field_name_list:
        pass_field_check_button = gtk.CheckButton(field_name, False)
        if (field_name in pass_field_clicked_list):
          pass_field_check_button.set_active(True)
        else:
          pass_field_check_button.set_active(False)

        pass_field_table.attach(pass_field_check_button, col_ind, col_ind+1,
                                row_num, row_num+1)
        pass_field_check_button.show()
        col_ind += 1
        if (col_ind == 3):
          col_ind =  0  # First column again
          row_num += 1  # Next row

      # Append pass field table to notebook page
      #
      run_std_box = self.mainTree.get_widget('run_std_box')
      if (len(run_std_box.get_children()) == 5):  # Remove previous table
        run_std_box.remove(run_std_box.get_children()[-1])
      run_std_box.pack_start(pass_field_table, False, False, 0)
      pass_field_table.show()

    else:  # Deduplication or linkage - - - - - - - - - - - - - - - - - - - - -
      run_nb_widget.set_current_page(0)

      # If no match data set names given, get original data set file names and
      # modify them
      #
      if (output_dict['m_datasets'][1][0] == '(None)'):
        org_file_name = self.data_set_info_list[0]['file_name']
        if (org_file_name[-3:] in ['.gz', '.GZ']):  # Remove suffix
          org_file_name = org_file_name[:-3]
        dot_index = org_file_name.rfind('.')
        if (dot_index > 0):
          match_dataset_a_name = org_file_name[:dot_index] + '-match' + \
                                 org_file_name[dot_index:]
        else:
          match_dataset_a_name = org_file_name
        output_dict['m_datasets'] = (output_dict['m_datasets'][0],
                                     (match_dataset_a_name,
                                      output_dict['m_datasets'][1][1]),
                                      output_dict['m_datasets'][2])

      if (self.project_type == 'Link'):
        if (output_dict['m_datasets'][2][0] == '(None)'):
          org_file_name = self.data_set_info_list[1]['file_name']
          if (org_file_name[-3:] in ['.gz', '.GZ']):  # Remove suffix
            org_file_name = org_file_name[:-3]
          dot_index = org_file_name.rfind('.')
          if (dot_index > 0):
            match_dataset_b_name = org_file_name[:dot_index] + '-match' + \
                                   org_file_name[dot_index:]
          else:
            match_dataset_b_name = org_file_name
          output_dict['m_datasets'] = (output_dict['m_datasets'][0],
                                       output_dict['m_datasets'][1],
                                       (match_dataset_b_name,
                                        output_dict['m_datasets'][2][1]))

      # Set progress percentage value - - - - - - - - - - - - - - - - - - - -
      #
      progress_perc_text_entry = \
                            self.mainTree.get_widget('run_link_progress_entry')
      if (output_dict['progress_perc'] != None):
        progress_perc_text_entry.set_text(output_dict['progress_perc'])
      else:
        progress_perc_text_entry.set_text('None')

      # Set length filtering percentage value - - - - - - - - - - - - - - - - -
      #
      len_filter_text_entry = \
                            self.mainTree.get_widget('run_length_filter_entry')
      if (output_dict['length_filter_perc'] != None):
        len_filter_text_entry.set_text(output_dict['length_filter_perc'])
      else:
        len_filter_text_entry.set_text('None')

      # Set cut-off threshold value - - - - - - - - - - - - - - - - - - - - - -
      #
      coff_th_text_entry = \
                         self.mainTree.get_widget('run_cutoff_threshold_entry')
      if (output_dict['cut_off_threshold'] != None):
        coff_th_text_entry.set_text(output_dict['cut_off_threshold'])
      else:
        coff_th_text_entry.set_text('None')

      # Set ticks and file names for output files - - - - - - - - - - - - - - -
      #
      w_vec_details = output_dict['w_vec_file']
      w_vec_tick_box = self.mainTree.get_widget('w_vec_file_checkbutton')
      w_vec_tick_box.set_active(w_vec_details[0])
      w_vec_file_button = self.mainTree.get_widget('save_w_vec_file_button')
      if (w_vec_details[0] == False):
        w_vec_file_button.set_sensitive(False)
      else:
        w_vec_file_button.set_sensitive(True)
        w_vec_file_button.set_label(w_vec_details[1].split(os.sep)[-1])

      histo_details = output_dict['histo_file']
      histo_tick_box = self.mainTree.get_widget('histogram_file_checkbutton')
      histo_tick_box.set_active(histo_details[0])
      histo_file_button = \
                         self.mainTree.get_widget('save_histogram_file_button')
      histo_bin_width_label = self.mainTree.get_widget('run_bin_width_label')
      histo_bin_width_entry = self.mainTree.get_widget('run_bin_width_entry')
      histo_bin_width_entry.set_text(histo_details[2])
      if (histo_details[0] == False):
        histo_file_button.set_sensitive(False)
        histo_bin_width_label.set_sensitive(False)
        histo_bin_width_entry.set_sensitive(False)
      else:
        histo_file_button.set_sensitive(True)
        histo_file_button.set_label(histo_details[1].split(os.sep)[-1])
        histo_bin_width_label.set_sensitive(True)
        histo_bin_width_entry.set_sensitive(True)

      m_stat_details = output_dict['m_status_file']
      m_stat_tick_box = self.mainTree.get_widget('m_status_file_checkbutton')
      m_stat_tick_box.set_active(m_stat_details[0])
      m_stat_file_button = \
                          self.mainTree.get_widget('save_m_status_file_button')
      if (m_stat_details[0] == False):
        m_stat_file_button.set_sensitive(False)
      else:
        m_stat_file_button.set_sensitive(True)
        m_stat_file_button.set_label(m_stat_details[1].split(os.sep)[-1])

      # If a linkage set second data set widgets to not visible - - - - - - - -
      #
      for wid_name in ['run_data_set_b_label', 'save_data_set_b_file_button',
                       'data_set_b_field_label', 'run_data_set_b_field_entry']:
        if (self.project_type == 'Link'):
          self.mainTree.get_widget(wid_name).show()
        else:  # Hide otherwise
          self.mainTree.get_widget(wid_name).hide()

      m_dataset_details = output_dict['m_datasets']
      m_dataset_tick_box = \
                        self.mainTree.get_widget('m_data_set_file_checkbutton')
      m_dataset_tick_box.set_active(m_dataset_details[0])

      if (m_dataset_details[0] == False):
        ds_a_sens = False
        ds_b_sens = False
      else:
        ds_a_sens = True
        ds_b_sens = (self.project_type == 'Link')  # True for a linkage only

      # Set sensitive for all data set widgets - - - - - - - - - - - - - - - -
      #
      for wid_name in ['run_data_set_a_label', 'save_data_set_a_file_button',
                       'data_set_a_field_label', 'run_data_set_a_field_entry']:
        self.mainTree.get_widget(wid_name).set_sensitive(ds_a_sens)
      for wid_name in ['run_data_set_b_label', 'save_data_set_b_file_button',
                       'data_set_b_field_label', 'run_data_set_b_field_entry']:
        self.mainTree.get_widget(wid_name).set_sensitive(ds_b_sens)

      if (m_dataset_details[0] == True):  # Sensitive
        dsa_details = m_dataset_details[1]
        dsa_f_button = self.mainTree.get_widget('save_data_set_a_file_button')
        dsa_f_button.set_label(dsa_details[0].split(os.sep)[-1])
        dsa_field_ent = self.mainTree.get_widget('run_data_set_a_field_entry')
        dsa_field_ent.set_text(dsa_details[1])

        if (self.project_type == 'Link'):  # Also show second data set
          dsb_details = m_dataset_details[2]
          dsb_f_button = \
                        self.mainTree.get_widget('save_data_set_b_file_button')
          dsb_f_button.set_label(dsb_details[0].split(os.sep)[-1])
          dsb_field_ent = \
                         self.mainTree.get_widget('run_data_set_b_field_entry')
          dsb_field_ent.set_text(dsb_details[1])

  # ---------------------------------------------------------------------------
  # Handle changes of the tick boxes
  #
  def runSaveFileButtonToggle(self, widget):
    widget_name = widget.get_name()
    print 'Toggled one of the save file buttons: %s' % (widget_name)

    if ('w_vec' in widget_name):
      self.output_dict['w_vec_file'] = (not self.output_dict['w_vec_file'][0],
                                        self.output_dict['w_vec_file'][1])
    elif ('histo' in widget_name):
      self.output_dict['histo_file'] = (not self.output_dict['histo_file'][0],
                                        self.output_dict['histo_file'][1],
                                        self.output_dict['histo_file'][2])
    elif ('status' in widget_name):
      self.output_dict['m_status_file'] = \
                                    (not self.output_dict['m_status_file'][0],
                                     self.output_dict['m_status_file'][1])
    elif ('data_set' in widget_name):
      self.output_dict['m_datasets'] = \
                                       (not self.output_dict['m_datasets'][0],
                                        self.output_dict['m_datasets'][1],
                                        self.output_dict['m_datasets'][2])
    else:
      print widget_name
      raise Exception

    self.runView()  # Re-display Output/Run page

  # ---------------------------------------------------------------------------
  # Handle clicks on one of the file name button
  #
  def runClickedFileNameButton(self, widget):
    widget_name = widget.get_name()
    print 'Clicked on of the file name buttons: %s' % (widget_name)

    # Create file dialog
    #
    saveFileTree = gtk.glade.XML(self.glade_file_name,
                                 self.save_file_dialog_name)
    saveFileDialog = saveFileTree.get_widget(self.save_file_dialog_name)

    # Get current name from dictionary
    #
    if ('w_vec') in widget_name:
      curr_file_name = self.output_dict['w_vec_file'][1].split(os.sep)[-1]
    elif ('histo') in widget_name:
      curr_file_name = self.output_dict['histo_file'][1].split(os.sep)[-1]
    elif ('status') in widget_name:
      curr_file_name = self.output_dict['m_status_file'][1].split(os.sep)[-1]
    elif ('data_set_a') in widget_name:
      curr_file_name = self.output_dict['m_datasets'][1][0].split(os.sep)[-1]
    elif ('data_set_b') in widget_name:
      curr_file_name = self.output_dict['m_datasets'][2][0].split(os.sep)[-1]
    elif ('std_out_file') in widget_name:  # For a standardisation
      curr_file_name = self.output_dict['std_out_file'].split(os.sep)[-1]
    else:
      curr_file_name = '(None)'

    if (curr_file_name != '(None)'):
      saveFileDialog.set_current_name(curr_file_name)

    saveFileDialog.show()
    save_file_response = saveFileDialog.run()
    self.save_dialog_file_name = saveFileDialog.get_filename()
    saveFileDialog.destroy()

    # Only process if 'OK' was clicked, not 'Cancel'
    #
    if (self.save_dialog_file_name != None):

      # Special case - if None was typed in, set to (None)
      #
      if (self.save_dialog_file_name.split(os.sep)[-1].lower() in ['(none)',
                                                              'none']):
        self.save_dialog_file_name = '(None)'

      if ('w_vec') in widget_name:
        # Check if file name has changed, if so re-do index.run()
        #
        if (self.output_dict['w_vec_file'][1] != self.save_dialog_file_name):
          self.re_run['w_vec_generate'] = True
          print 'XXXX changed w-vec file name', self.save_dialog_file_name

        self.output_dict['w_vec_file'] = (self.output_dict['w_vec_file'][0],
                                          self.save_dialog_file_name)
      elif ('histo') in widget_name:
        self.output_dict['histo_file'] = (self.output_dict['histo_file'][0],
                                          self.save_dialog_file_name,
                                          self.output_dict['histo_file'][2])
      elif ('status') in widget_name:
        self.output_dict['m_status_file'] = \
             (self.output_dict['m_status_file'][0], self.save_dialog_file_name)

      elif ('data_set_a') in widget_name:
        ds_details = self.output_dict['m_datasets']
        ds_details = (ds_details[0],
                      (self.save_dialog_file_name, ds_details[1][1]),
                      ds_details[2])
        self.output_dict['m_datasets'] = ds_details

      elif ('data_set_b') in widget_name:
        ds_details = self.output_dict['m_datasets']
        ds_details = (ds_details[0],
                      ds_details[1],
                      (self.save_dialog_file_name, ds_details[2][1]))
        self.output_dict['m_datasets'] = ds_details

      elif ('std_out_file') in widget_name:  # For a standardisation
        self.output_dict['std_out_file'] = self.save_dialog_file_name

      else:
        print widget_name, self.save_dialog_file_name
        raise Exception

      self.runView()  # Re-display Output/Run page

  # ---------------------------------------------------------------------------
  # Handle an activate of Execute on the Output/Run page
  #
  def runExecute(self):
    output_dict = self.output_dict  # Shortcut

    print 'output_dict (1)', output_dict

    self.febrl_code['output'] = None  # Remove all previous output code

    # Process page depending upon project type
    #
    if (self.project_type == 'Standardise'):  # - - - - - - - - - - - - - - - -

      # Get and check progress percentage value - - - - - - - - - - - - - - - -
      #
      progr_perc_val = \
                 self.mainTree.get_widget('run_std_progress_entry').get_text()
      progr_perc_val = progr_perc_val.strip()
      if (progr_perc_val.lower() in ['','none']):  # Do not report progress
        output_dict['progress_perc'] = None
      elif ((not self.str_is_percentage_not_zero(progr_perc_val)) or
            (float(progr_perc_val) > 50)):
        self.messageDialog('Progress report percentage must be "None"\n' + \
                           'or a percentage value larger than 0 and \n'+ \
                           'and smaller or equal to 50.', 'warn')
        return
      else:
        output_dict['progress_perc'] = progr_perc_val

      if (output_dict['std_out_file'] == '(None)'):
        self.messageDialog('No output file name given.', 'warn')
        return

      run_std_box = self.mainTree.get_widget('run_std_box')
      pass_field_list = []

      pass_field_table = run_std_box.get_children()[-1]
      for check_button in pass_field_table.get_children():

        if (check_button.get_active() == True):
          pass_field_list.append(check_button.get_label())

      pass_field_list.reverse()  # As children are bottom to top in table
      output_dict['pass_field_list'] = pass_field_list

      print 'output_dict (2)', output_dict

      # Now generate the Febrl codes for output - - - - - - - - - - - - - - - -
      #
      febrl_code = []
      febrl_code.append('# '+'-'*77)
      febrl_code.append('')
      febrl_code.append('# Define standardised output data set and record ' + \
                        'standardiser')
      febrl_code.append('#')

      # Generate code for standardised output data set
      #
      febrl_code.append('data_set_a_std = dataset.DataSetCSV(description=' + \
                        '"Standardised output data set",')
      indent_space = ' '*36
      febrl_code.append(indent_space+'access_mode="write",')
      febrl_code.append(indent_space+'delimiter=",",')

      out_file_name = output_dict['std_out_file']
      febrl_code.append(indent_space+'file_name="%s",' % (out_file_name))

      # Add fields
      #
      field_num = 0  # Field number / column index in CSV output file
      all_field_names = []

      comp_std_name_list = []  # Component standardiser names for record std.

      comp_std_list = self.comp_std_list  # Short cut
      for i in range(len(comp_std_list)):
        cs_out_field_list = comp_std_list[i]['out_field_list']

        comp_std_type = comp_std_list[i]['type']
        if (comp_std_type == 'Date'):
          comp_std_name_list.append('date_comp_std_%d' % (i+1))
        elif (comp_std_type == 'Phon'):
          comp_std_name_list.append('phone_comp_std_%d' % (i+1))
        elif (comp_std_type == 'Addr'):
          comp_std_name_list.append('address_comp_std_%d' % (i+1))
        elif (comp_std_type == 'Name'):
          comp_std_name_list.append('name_comp_std_%d' % (i+1))
        else:
          raise Exception, 'This should not happen: %s' % (comp_std_type)

        for out_field in cs_out_field_list:

          if (out_field != 'None'):
            all_field_names.append(out_field)
            if (field_num == 0):  # First field
              febrl_code.append(indent_space+'field_list = [("%s",%d),' % \
                                (out_field, field_num))
            else:
              febrl_code.append(indent_space + ' '*14 + '("%s",%d),' % \
                                (out_field, field_num))
            field_num += 1

      # Add pass fields
      #
      for pass_field in pass_field_list:
        all_field_names.append(pass_field)
        febrl_code.append(indent_space + ' '*14 + '("%s",%d),' % \
                              (pass_field, field_num))
        field_num += 1

      febrl_code[-1] = febrl_code[-1][:-1]+'],'  # Close field list

      # Generate a record identifier field
      #
      rec_id_str = self.data_set_default_values['rec_id_field']

      # Check identifier field name it is not used
      #
      while (rec_id_str in all_field_names):
        rec_id_str = '_'+rec_id_str+'_'
      febrl_code.append(indent_space+'rec_ident="%s",' % (rec_id_str))

      febrl_code.append(indent_space+'strip_fields=True,')
      febrl_code.append(indent_space+'header_line=True,')
      febrl_code.append(indent_space+'write_header=True)')

      # Generate code for record standardiser
      #
      febrl_code.append('')
      febrl_code.append('# Define record standardiser')
      febrl_code.append('#')
      febrl_code.append('rec_std = standardisation.RecordStandardiser(' + \
                        'descr = "Febrl record standardiser",')
      indent_space = ' '*45
      febrl_code.append(indent_space + 'input_data =  data_set_a,')
      febrl_code.append(indent_space + 'output_data = data_set_a_std,')
      febrl_code.append(indent_space + 'comp_stand_list = [%s],' % \
                        ','.join(comp_std_name_list))
      if (output_dict['progress_perc'] != None):
        febrl_code.append(indent_space + 'progress_report = %s,' %
                          output_dict['progress_perc'])

      if (pass_field_list == []):
        febrl_code.append(indent_space + 'pass_field_list = [])')
      else:
        first_pass_field = True   # Flag for first or not first
        for pass_field in pass_field_list:
          if (first_pass_field == True):  # First field
            febrl_code.append(indent_space+'pass_field_list = ' + \
                              '[("%s","%s"),' % (pass_field, pass_field))
            first_pass_field = False
          else:
            febrl_code.append(indent_space+' '*19 + '("%s","%s"),' % \
                              (pass_field, pass_field))

        # Close pass field list and record standardiser
        #
        febrl_code[-1] = febrl_code[-1][:-1]+'])'

      febrl_code.append('')
      febrl_code.append('# Start standardisation')
      febrl_code.append('#')
      febrl_code.append('rec_std.standardise()')

      febrl_code.append('')

      self.febrl_code['output'] = febrl_code  # Store for later use

    else:  # Deduplication or linkage - - - - - - - - - - - - - - - - - - - - -

      # Get and check progress percentage value - - - - - - - - - - - - - - - -
      #
      progr_perc_val = \
                 self.mainTree.get_widget('run_link_progress_entry').get_text()
      progr_perc_val = progr_perc_val.strip()
      if (progr_perc_val.lower() in ['','none']):  # Do not report progress
        output_dict['progress_perc'] = None
      elif ((not self.str_is_percentage_not_zero(progr_perc_val)) or
            (float(progr_perc_val) > 50)):
        self.messageDialog('Progress report percentage must be "None"\n' + \
                           'or a percentage value larger than 0 and \n'+ \
                           'and smaller or equal to 50.', 'warn')
        return
      else:
        output_dict['progress_perc'] = progr_perc_val

      # Get and check length filtering percentage value - - - - - - - - - - - -
      #
      len_filter_val = \
                 self.mainTree.get_widget('run_length_filter_entry').get_text()
      len_filter_val = len_filter_val.strip()
      if (len_filter_val.lower() in ['','none', '(none)']): # Dont use len-filt
        new_len_filter_val = None
      elif (self.str_is_percentage_not_zero(len_filter_val)):
        new_len_filter_val = len_filter_val
      else:
        self.messageDialog('Length filtering percentage must be "None"\n' + \
                           'or a percentage value larger than 0.', 'warn')
        return

      if (output_dict.get('length_filter_perc', new_len_filter_val) != \
          new_len_filter_val):
        self.re_run['w_vec_generate'] = True
      output_dict['length_filter_perc'] = new_len_filter_val

      # Get and check cut-off threshold value - - - - - - - - - - - - - - - - -
      #
      coff_th_val = \
              self.mainTree.get_widget('run_cutoff_threshold_entry').get_text()
      coff_th_val = coff_th_val.strip()
      if (coff_th_val.lower() in ['','none','(none)']): # Dont use cut-off thr.
        new_coff_th_val = None
      elif (self.str_is_float(coff_th_val)):
        new_coff_th_val = coff_th_val
      else:
        self.messageDialog('Cut-off threshold must be a number.', 'warn')
        return

      if (output_dict.get('cut_off_threshold', new_coff_th_val) != \
          new_coff_th_val):
        self.re_run['w_vec_generate'] = True

      output_dict['cut_off_threshold'] = new_coff_th_val

      # Get output file settings - - - - - - - - - - - - - - - - - - - - - - -
      # (if ticked check if corresponding file names given)
      #
      if (output_dict['w_vec_file'][0] == True):
        if (output_dict['w_vec_file'][1] == '(None)'):
          self.messageDialog('No weight vector file name given.', 'warn')
          return

      if (output_dict['histo_file'][0] == True):
        if (output_dict['histo_file'][1] == '(None)'):
          self.messageDialog('No histogram file name given.', 'warn')
          return

        histo_bin_width_entry = self.mainTree.get_widget('run_bin_width_entry')
        bin_width = histo_bin_width_entry.get_text().strip()
        if ((not self.str_is_float(bin_width)) or (float(bin_width) <= 0.0)):
          self.messageDialog('Bin width must be a positive number.', 'warn')
          return
        output_dict['histo_file'] = \
                                 (True, output_dict['histo_file'][1],bin_width)

      if (output_dict['m_status_file'][0] == True):
        if (output_dict['m_status_file'][1] == '(None)'):
          self.messageDialog('No match status file name given.', 'warn')
          return

      if (output_dict['m_datasets'][0] == True):
        if (self.project_type == 'Link'):
          if ((output_dict['m_datasets'][1][0] == '(None)') or \
              (output_dict['m_datasets'][2][0] == '(None)')):
            self.messageDialog('Two match data set file names\n' + \
                               'need to be given for a linkage.', 'warn')
            return
          if (output_dict['m_datasets'][1][0]== \
              output_dict['m_datasets'][2][0]):
            self.messageDialog('Match data set file names for\n' + \
                               'linkage have to be different.', 'warn')
            return

          # Check output files are different from input files
          #
          in_file_name_a = self.data_set_info_list[0]['file_name']
          in_file_name_b = self.data_set_info_list[1]['file_name']

          if ((in_file_name_a == output_dict['m_datasets'][1][0]) or \
              (in_file_name_b == output_dict['m_datasets'][2][0])):
            self.messageDialog('Match data set file names have to be ' + \
                               'different from input data set file names.',
                               'warn')
            return

          fna = \
              self.mainTree.get_widget('run_data_set_a_field_entry').get_text()
          fna = fna.strip()
          fnb = \
              self.mainTree.get_widget('run_data_set_b_field_entry').get_text()
          fnb = fnb.strip()
          if ((fna == '') or (fnb == '')):
            self.messageDialog('Two match identifier field names\n' + \
                               'need to be given for a linkage.', 'warn')
            return

          # Check if the identifier field names are different from all other
          # field names
          #
          field_name_a_list = self.data_set_info_list[0]['field_names']
          field_name_b_list = self.data_set_info_list[1]['field_names']
          if ((fna in field_name_a_list) or (fnb in field_name_b_list)):
            self.messageDialog('Match identifier field names have to be ' + \
                               'different from all other data set field ' + \
                               'names.', 'warn')
            return

        else:  # Deduplication
          if (output_dict['m_datasets'][1][0] == '(None)'):
            self.messageDialog('No match data set file name given.', 'warn')
            return

          in_file_name_a = self.data_set_info_list[0]['file_name']
          if (in_file_name_a == output_dict['m_datasets'][1][0]):
            self.messageDialog('Match data set file name has to be different' \
                               + ' from input data set file name.', 'warn')
            return

          fna = \
              self.mainTree.get_widget('run_data_set_a_field_entry').get_text()
          fna = fna.strip()
          fnb = None
          if (fna == ''):
            self.messageDialog('A match identifier field name\n' + \
                               'need to be given.', 'warn')
            return

          field_name_a_list = self.data_set_info_list[0]['field_names']
          print field_name_a_list, fna
          if (fna in field_name_a_list):
            self.messageDialog('Match identifier field name has to be ' + \
                               'different from all other data set field ' + \
                               'names.', 'warn')
            return

        self.output_dict['m_datasets'] = (True, \
                                          (output_dict['m_datasets'][1][0],
                                           fna),
                                          (output_dict['m_datasets'][2][0],
                                           fnb))

      print 'output_dict (2)', output_dict

      # Now generate the Febrl codes for output - - - - - - - - - - - - - - - -
      #
      febrl_code = []
      febrl_code.append('# '+'-'*77)
      febrl_code.append('')
      febrl_code.append('# Define output file options')
      febrl_code.append('#')

      # A histogram will be generated in any case - - - - - - - - - - - - - - -
      # (possibly not saved into a file though)
      #
      if (output_dict['histo_file'][0] == False):  # Don't save into file
        febrl_code.append('histo_str_list = output.GenerateHistogram(' + \
                          'class_w_vec_dict, %s)' % \
                          (output_dict['histo_file'][2]))
      else:
        febrl_code.append('histo_str_list = output.GenerateHistogram(' + \
                          'class_w_vec_dict, %s, "%s")' % \
                          (output_dict['histo_file'][2],
                           output_dict['histo_file'][1]))
      febrl_code.append('')

      febrl_code.append('for line in histo_str_list:')  # Print histogram
      febrl_code.append('  print line')

      # Check if match status file shall be written - - - - - - - - - - - - - -
      #
      if (output_dict['m_status_file'][0] == True):
        febrl_code.append('output.SaveMatchStatusFile(class_w_vec_dict, ' + \
                          'm_set, "%s")' % (output_dict['m_status_file'][1]))
        febrl_code.append('')

      # Check if match data set shall be written - - - - - - - - - - - - - - -
      #
      if (output_dict['m_datasets'][0] == True):
        if (self.project_type == 'Link'):
          febrl_code.append('output.SaveMatchDataSet(m_set, data_' + \
                            'set_a, "%s", "%s", data_set_b, "%s", "%s")' % \
                            (output_dict['m_datasets'][1][1],
                             output_dict['m_datasets'][1][0],
                             output_dict['m_datasets'][2][1],
                             output_dict['m_datasets'][2][0]))
        else:  # Deduplication
          febrl_code.append('output.SaveMatchDataSet(m_set, data_set_a, ' + \
                            '"%s", "%s")' % \
                            (output_dict['m_datasets'][1][1],
                             output_dict['m_datasets'][1][0]))
        febrl_code.append('')

      self.febrl_code['output'] = febrl_code  # Store for later use

    # Finally update the GUI information - - - - - - - - - - - - - - - - - - -
    #
    self.addToLog('')  # Add generated code into log page text
    self.addToLog('='*79)
    self.addToLog('Generated Febrl code for "output" on %s' % \
                  (time.asctime()))
    self.addToLog('')
    for line in febrl_code:
      self.addToLog(line)
    self.addToLog('')

    self.modified['classification'] = True  # Details have been changed
    self.setWindowTitle()

    self.writeStatusBar('Generated Febrl Python code for output files ' + \
                        '(see Log page for generated code).')

    # Ask to save project code - - - - - - - - - - - - - - - - - - - - - - - -
    #
    do_save = self.messageDialog('Do you want to save final project code ' + \
                                 'before running it?', 'c_question')
    if (do_save == True):
      self.saveMain(gtk.HSeparator())  # Need to provide a widget

    do_run = self.messageDialog('Do you want to run the project code now?',
                                'c_question')
    if (do_run == True):
      if (self.project_type == 'Standardise'):
        self.runStdCode()
      elif (self.project_type in ['Deduplicate', 'Link']):
        self.runLinkCode()
      else:
        pass  # To do for geocoding

  # ---------------------------------------------------------------------------
  # Method to generate all code for a Febrl project module.
  #
  def generateAllCode(self):
    print 'Generate all code.'

    license_boiler_pate = [
      'AUSTRALIAN NATIONAL UNIVERSITY OPEN SOURCE LICENSE (ANUOS LICENSE)',
      'VERSION 1.3',
      '',
      'The contents of this file are subject to the ANUOS License Version 1.2',
      '(the "License"); you may not use this file except in compliance with',
      'the License. You may obtain a copy of the License at:',
      '',
      '  http://datamining.anu.edu.au/linkage.html',
      '',
      'Software distributed under the License is distributed on an "AS IS"',
      'basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See',
      'the License for the specific language governing rights and limitations',
      'under the License.',
      '',
      'The Original Software is: "XXXXXXXXXX"',  # To be replaced later
      '',
      'The Initial Developers of the Original Software are:',
      '  Peter Christen',
      '',
      'Copyright (C) 2002 - 2011 the Australian National University and',
      'others. All Rights Reserved.',
      '',
      'Contributors:',
      '',
      'Alternatively, the contents of this file may be used under the terms',
      'of the GNU General Public License Version 2 or later (the "GPL"), in',
      'which case the provisions of the GPL are applicable instead of those',
      'above. The GPL is available at the following URL: http://www.gnu.org/',
      'If you wish to allow use of your version of this file only under the',
      'terms of the GPL, and not to allow others to use your version of this',
      'file under the terms of the ANUOS License, indicate your decision by',
      'deleting the provisions above and replace them with the notice and',
      'other provisions required by the GPL. If you do not delete the',
      'provisions above, a recipient may use your version of this file under',
      'the terms of any one of the ANUOS License or the GPL.']

    febrl_code = []

    febrl_code.append('# '+'='*77)  # Start with license boiler plate
    for line in license_boiler_pate:
      febrl_code.append('# %s' % (line))
    febrl_code.append('# '+'='*77)
    febrl_code.append('')

    febrl_code.append('# '+'='*77)
    febrl_code.append('# Start of Febrl project module: "XXXXXXXXXX"')
    febrl_code.append('#')
    febrl_code.append('# Generated using "guiFebrl.py" on %s' % \
                      (time.asctime()))
    febrl_code.append('# '+'='*77)
    febrl_code.append('')

    febrl_code.append('# Import necessary modules (Python standard modules' + \
                      ' first, then Febrl modules)')
    febrl_code.append('')

    febrl_code.append('import logging')  # Standard Python modules
    febrl_code.append('')

    if (self.project_type == 'Standardise'):
      febrl_code.append('import dataset')  # Required Febrl modules
      febrl_code.append('import simplehmm')
      febrl_code.append('import standardisation')

    else:  # Deduplication or linkage
      febrl_code.append('import classification')  # Required Febrl modules
      febrl_code.append('import comparison')
      febrl_code.append('import dataset')
      febrl_code.append('import encode')
      febrl_code.append('import indexing')
      febrl_code.append('import measurements')
      febrl_code.append('import mymath')
      febrl_code.append('import output')
      febrl_code.append('import stringcmp')
    febrl_code.append('')

    febrl_code.append('# '+'-'*77)
    febrl_code.append('# Intialise a logger, set level to info oe warning')
    febrl_code.append('#')
    febrl_code.append('log_level = logging.INFO # logging.WARNING')
    febrl_code.append('')
    febrl_code.append('my_logger = logging.getLogger()')
    febrl_code.append('my_logger.setLevel(log_level)')
    febrl_code.append('')

    # Add a marker that gives the project type - - - - - - - - - - - - - - - -
    #
    febrl_code.append('# '+'-'*77)
    febrl_code.append('# Febrl project type: %s' % (self.project_type))
    febrl_code.append('# '+'-'*77)
    febrl_code.append('')

    # Append data set initialisation code - - - - - - - - - - - - - - - - - - -
    #
    if (self.febrl_code['data'] == None):
      febrl_code.append('# WARNING: No code for data set initialisation ' + \
                        ' provided')
      febrl_code.append('')
    else:
      febrl_code += self.febrl_code['data']

    # Generate code for a standardisation project - - - - - - - - - - - - - - -
    #
    if (self.project_type == 'Standardise'):

      # Append standardisation code
      #
      if (self.febrl_code['standardise'] == None):
        febrl_code.append('# WARNING: No code for standardisation provided')
        febrl_code.append('')
      else:
        febrl_code += self.febrl_code['standardise']

    else:  # Code for a deduplication or linkage - - - - - - - - - - - - - - -

      # Append field and record comparison code - - - - - - - - - - - - - - - -
      #
      if (self.febrl_code['comparison'] == None):
        febrl_code.append('# WARNING: No code for field and record ' + \
                          'comparisons provided')
        febrl_code.append('')
      else:
        febrl_code += self.febrl_code['comparison']

      # Get indexing code and insert progress report and weight vector file - -
      #
      if (self.febrl_code['indexing'] == None):
        febrl_code.append('# WARNING: No code for index definitions provided')
        febrl_code.append('')
      else:
        febrl_indexing_code = self.febrl_code['indexing']

        if ((self.output_dict['w_vec_file'][0] == True) or \
            (self.output_dict['progress_perc'] != None)):

          mod_indexing_code = []  # Modified list of indexing code

          # Find line in indexing code where index is defined
          #
          find_ind = 0
          while ('indexing.' not in febrl_indexing_code[find_ind]):
            mod_indexing_code.append(febrl_indexing_code[find_ind])
            find_ind += 1

          indent_space = ' '*(febrl_indexing_code[find_ind].find('(')+1)

          # Insert after 'dataset' argument(s)
          #
          mod_indexing_code.append(febrl_indexing_code[find_ind])
          mod_indexing_code.append(febrl_indexing_code[find_ind+1])
          find_ind += 2

          if (self.project_type == 'Link'):
            mod_indexing_code.append(febrl_indexing_code[find_ind])
            find_ind += 1

          len_diff = 0  # Number of lines inserted
          if (self.output_dict['w_vec_file'][0] == True):
            mod_indexing_code.append(indent_space+'weight_vec_file = "%s",' \
                                     % (self.output_dict['w_vec_file'][1]))
            len_diff += 1
          if (self.output_dict['progress_perc'] != 'False'):
            mod_indexing_code.append(indent_space+'progress_report = %s,' % \
                                     (self.output_dict['progress_perc']))
            len_diff += 1

          # Append rest of original indexing code
          #
          mod_indexing_code += febrl_indexing_code[find_ind:]

          assert len(mod_indexing_code) == len(febrl_indexing_code)+len_diff

          febrl_code += mod_indexing_code

        else:  # Just add original code

          febrl_code += febrl_indexing_code

      # Add code to build, compact the index and run the comparison step - - -
      #
      if (self.febrl_code['indexing'] != None):

        febrl_code.append('# Build and compact index')
        febrl_code.append('#')
        febrl_code.append('index.build()')
        febrl_code.append('')
        febrl_code.append('index.compact()')
        febrl_code.append('')

        run_parameter_str = ''  # Check if run() will have parameters

        if ((self.output_dict['length_filter_perc'] != None) or \
            (self.output_dict['cut_off_threshold'] != None)):

          if (self.output_dict['length_filter_perc'] != None):
            run_parameter_str += 'length_filter_perc = %s, ' % \
                                 (self.output_dict['length_filter_perc'])
          if (self.output_dict['cut_off_threshold'] != None):
            run_parameter_str += 'cut_off_threshold = %s, ' % \
                                 (self.output_dict['cut_off_threshold'])
          run_parameter_str = run_parameter_str[:-2]  # Remove comma and space

        febrl_code.append('# Do record pair comparisons')
        febrl_code.append('#')

        # run() returns results depending upon weight vector file definition
        #
        if (self.output_dict['w_vec_file'][0] == True):  # Nothing is returned
          if (',' not in run_parameter_str):
            febrl_code.append('index.run(%s)' % (run_parameter_str))
          else:  # Two parameters
            run_parameter_list = run_parameter_str.split(', ')
            febrl_code.append('index.run(%s,' % (run_parameter_list[0]))
            febrl_code.append('          %s)' % (run_parameter_list[1]))

          febrl_code.append('')

          # Now weight vector dictionary has to be loaded from its file
          #
          febrl_code.append('[field_names_list, w_vec_dict] = output.' +\
                            'LoadWeightVectorFile("%s")' % \
                            (self.output_dict['w_vec_file'][1]))

        else:  # run() returns field list and weight vector dictionary
          if (',' not in run_parameter_str):
            febrl_code.append('[field_names_list, w_vec_dict] = ' + \
                              'index.run(%s)' % (run_parameter_str))
          else:  # Two parameters
            run_parameter_list = run_parameter_str.split(', ')
            febrl_code.append('[field_names_list, w_vec_dict] = ' + \
                              'index.run(%s,' % (run_parameter_list[0]))
            febrl_code.append('                                        ' + \
                              '%s)' % (run_parameter_list[1]))
        febrl_code.append('')

      # Append classification code - - - - - - - - - - - - - - - - - - - - - -
      #
      if (self.febrl_code['classification'] == None):
        febrl_code.append('# WARNING: No code for classifiers provided')
        febrl_code.append('')
      else:
        febrl_code += self.febrl_code['classification']  # Initalise classifier

        # For supervised classifiers true matches and non-matches required
        #
        if (self.classifier_method['name'] in ['OptimalThreshold',
                                               'SuppVecMachine']):

          febrl_code.append('# Get true match and non-match status for ' + \
                            'supervised classifier')
          febrl_code.append('')

          match_status_field = self.classifier_method.get('match_status_ind',
                                                          -1)
          if (match_status_field < 0):
            febrl_code.append('# WARNING: Match status field not defined!')
            febrl_code.append('')

          else:  # Get the index of the match status field
            match_status_field_index = \
                                    self.exact_comp_list[match_status_field][1]
            print match_status_field_index

            febrl_code.append('# Define function to get the match status')
            febrl_code.append('#')
            febrl_code.append('def match_status_funct(red_id1,rec_id2,w_vec):')
            febrl_code.append('  return (w_vec[%d] == 1.0)' % \
                              (match_status_field_index))
            febrl_code.append('')

            febrl_code.append('# Get true match and non-match sets')
            febrl_code.append('#')
            febrl_code.append('[tm_set, tnm_set] = classification.' + \
                              'get_true_matches_nonmatches(w_vec_dict, ' + \
                              'match_status_funct)')
            febrl_code.append('')

            # Remove the match status element from the weight vectors - - - - -
            #
            num_comp_funct = len(self.field_comp_list)

            manipulate_list = []
            for i in range(num_comp_funct):
              if (i != match_status_field_index):
                manipulate_list.append((i,))

            febrl_code.append('# Remove the match status element from ' + \
                              'weight vectors')
            febrl_code.append('#')
            febrl_code.append('class_w_vec_dict = classification.' + \
                           'extract_collapse_weight_vectors(%s, w_vec_dict)' \
                           % (str(manipulate_list)))
            febrl_code.append('')

            febrl_code.append('# Supervised training of classifier')  # - - - -
            febrl_code.append('#')

            febrl_code.append('classifier.train(class_w_vec_dict, tm_set,' + \
                              ' tnm_set)')
            febrl_code.append('')

            febrl_code.append('# Test quality of trained classifier')
            febrl_code.append('#')
            febrl_code.append('[tp, fn, fp, tn] = classifier.test(' + \
                              'class_w_vec_dict, tm_set, tnm_set)')
            febrl_code.append('')
            febrl_code.append('print "Number of true positives: ", tp')
            febrl_code.append('print "Number of false positives:", fp')
            febrl_code.append('print "Number of true negatives: ", tn')
            febrl_code.append('print "Number of false negatives:", fn')

        else:  # Unsupervised classifiers

          febrl_code.append('# Unsupervised training of classifier')
          febrl_code.append('#')
          febrl_code.append('class_w_vec_dict = w_vec_dict  # Use orignal ' + \
                            'weight vector dictionary')
          febrl_code.append('')

          febrl_code.append('classifier.train(class_w_vec_dict, ' + \
                            'set(), set())') # Empty training sets

        febrl_code.append('')

        # Classify the weight vectors into matches, non-matches and possible
        # matches
        #
        febrl_code.append('# Classify all weight vectors')
        febrl_code.append('#')
        febrl_code.append('[m_set, nm_set, pm_set] = classifier.classify(' + \
                          'class_w_vec_dict)')
        febrl_code.append('')

        # For supervised classification get linkage quality and complexity - -
        #
        if (self.classifier_method['name'] in ['OptimalThreshold',
                                               'SuppVecMachine']):

          febrl_code.append('# Evaluate linkage quality and complexity')
          febrl_code.append('#')
          febrl_code.append('acc, prec, reca, fmeas = measurements.' + \
                            'quality_measures(w_vec_dict, m_set, nm_set, ' + \
                            'match_status_funct)')
          febrl_code.append('print "Accuracy: ", acc')
          febrl_code.append('print "Precision:", prec')
          febrl_code.append('print "Recall:   ", reca')
          febrl_code.append('print "F-measure:", fmeas')
          febrl_code.append('')

          if (self.project_type == 'Link'):
            febrl_code.append('rr = measurements.reduction_ratio(' + \
                              'w_vec_dict, data_set_a, data_set_b)')
          else:
            febrl_code.append('rr = measurements.reduction_ratio(' + \
                              'w_vec_dict, data_set_a, data_set_a)')
          febrl_code.append('print "Reduction ratio:", rr')
          febrl_code.append('')

          febrl_code.append('pq = measurements.pairs_quality(w_vec_dict, ' + \
                            'match_status_funct)')
          febrl_code.append('print "Pairs quality:  ", pq')
          febrl_code.append('')

          # ADD pairs_completeness(w_vec_dict, csv_ds1, csv_ds2,check_funct) ##

    # Append output files code - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.febrl_code['output'] == None):
      febrl_code.append('# WARNING: No code for output files provided')
      febrl_code.append('')
    else:
      febrl_code += self.febrl_code['output']

    febrl_code.append('')
    febrl_code.append('# '+'='*77)
    febrl_code.append('# End of Febrl project module: "XXXXXXXXXX"')
    febrl_code.append('# '+'='*77)

    return febrl_code

  # ---------------------------------------------------------------------------
  # Method to run code for a Febrl standardisation project
  #
  def runStdCode(self):
    print 'Run standardisation code:'

    # Check if all necessary components have been initialised
    #
    if ((self.febrl_code['data'] ==  None) or
        (self.febrl_code['standardise'] ==  None) or
        (self.febrl_code['output'] ==  None)):
      self.messageDialog('Cannot run standardisation project\n' + \
                         'as not all required code components\n' + \
                         'are initialised.', 'error')
      return

    output_dict =   self.output_dict  # Shortcuts
    comp_std_list = self.comp_std_list

    # Initialise the input data set - - - - - - - - - - - - - - - - - - - - - -
    #
    in_data_set_info = self.data_set_info_list[0]

    if (self.data_set_type_list[0] in ['CSV','TAB']):

      field_val_list = []
      for field_num in range(len(in_data_set_info['field_names'])):
        field_val_list.append((in_data_set_info['field_names'][field_num],
                               field_num))

      ids = dataset.DataSetCSV(access_mode="read",
                               file_name =    in_data_set_info['file_name'],
                               delimiter =    in_data_set_info['delimiter'],
                               rec_ident =    in_data_set_info['rec_id_field'],
                               header_line =  in_data_set_info['header_line'],
                               strip_fields = in_data_set_info['strip_fields'],
                               miss_val =     in_data_set_info['miss_val'],
                               field_list =   field_val_list)

    elif (self.data_set_type_list[0] == 'COL'):

      # Get width of the columns from GUI
      #
      col_width_vals = self.mainTree.get_widget('col_width_entry_a').get_text()
      col_width_vals = col_width_vals.strip()
      col_width_list = col_width_vals.split(',')

      field_val_list = []
      for field_num in range(len(in_data_set_info['field_names'])):
        if (header_line_flag == True):
          field_val_list.append(int(col_width_list[field_num]))
        else:
          field_val_list.append((field_names[field_num],
                                 int(col_width_list[field_num])))

      ids = dataset.DataSetCOL(access_mode="read",
                               file_name =    in_data_set_info['file_name'],
                               rec_ident =    in_data_set_info['rec_id_field'],
                               header_line =  in_data_set_info['header_line'],
                               strip_fields = in_data_set_info['strip_fields'],
                               miss_val =     in_data_set_info['miss_val'],
                               field_list =   field_val_list)

    else:
      pass  # SQL to be implemented ########################################

    self.writeStatusBar('Initalised input data set with %d records.' % \
                        (ids.num_records))

    # Initialise the standardised output data set - - - - - - - - - - - - - - -
    #
    pass_field_list = output_dict['pass_field_list']

    out_field_list =  []
    all_field_names = []
    field_num = 0  # Field number / column index in CSV output file

    for i in range(len(comp_std_list)):
      cs_out_field_list = comp_std_list[i]['out_field_list']

      for out_field in cs_out_field_list:
        if (out_field != 'None'):
          all_field_names.append(out_field)
          out_field_list.append((out_field, field_num))
          field_num += 1

    for pass_field in pass_field_list:  # Add pass fields
      all_field_names.append(pass_field)
      out_field_list.append((pass_field, field_num))
      field_num += 1

    # Generate a record identifier field
    #
    rec_id_str = self.data_set_default_values['rec_id_field']
    while (rec_id_str in all_field_names):
      rec_id_str = '_'+rec_id_str+'_'

    ods = dataset.DataSetCSV(access_mode = "write",
                             file_name = output_dict['std_out_file'],
                             delimiter = ',',
                             rec_ident = rec_id_str,
                             field_list = out_field_list,
                             strip_fields=True,
                             header_line=True,
                             write_header=True)

    # Initialise the component standardisers and if necessary lookup tables,
    # correction lists, etc.
    #
    cs_list = []
    for cs_dict in comp_std_list:
      comp_std_type =  cs_dict['type']
      in_field_list =  cs_dict['in_field_list']
      out_field_list = []

      for out_field in cs_dict['out_field_list']:
        if (out_field != 'None'):
          out_field_list.append(out_field)
        else:
          out_field_list.append(None)

      # Check for correction list initialisation for this comp. standardiser
      #
      corr_list_name = cs_dict.get('corr_list', '(None)')
      if (corr_list_name != '(None)'):
        this_corr_list = lookup.CorrectionList()
        this_corr_list.load(corr_list_name)
      else:
        this_corr_list = []

      # Check for tag lookup table initialisation for this comp. standardiser
      #
      tag_table_name = cs_dict.get('tag_table', '(None)')
      if (tag_table_name != '(None)'):
        this_tag_table = lookup.TagLookupTable()
        this_tag_table.load(tag_table_name)
      else:
        this_tag_table = {}

## Changed PC 25/08/2008
#
#      if (comp_std_type in ['Addr','Name']):
#        field_sep =       cs_dict['field_separator']
#        word_spill_flag = cs_dict['check_word_spill']
      field_sep =       cs_dict.get('field_separator', '')
      word_spill_flag = cs_dict.get('check_word_spill', False)

      # Initialise the component standardisers
      #
      if (comp_std_type == 'Date'):  # - - - - - - - - - - - - - - - - - - - -
        pivot_year = int(cs_dict['pivot_year'])
        parse_f_list = []
        for parse_str in cs_dict['parse_formats'].split(','):
          parse_f_list.append(parse_str.strip())

        this_cs = standardisation.DateStandardiser(input_fiel = in_field_list,
                                                   output_fiel= out_field_list,
                                                   pivot_year = pivot_year,
                                                   parse_formats= parse_f_list)

      elif (comp_std_type == 'Phon'):  # - - - - - - - - - - - - - - - - - - -
        def_country = cs_dict['default_country']

        this_cs = standardisation.PhoneNumStandardiser(input_fi= in_field_list,
                                                   output_fiel= out_field_list,
                                                   corr_list = this_corr_list,
                                                   tag_table = this_tag_table,
                                                   field_separator = field_sep,
                                                   check_wor = word_spill_flag,
                                                   default_co = def_country)

      elif (comp_std_type == 'Addr'):  # - - - - - - - - - - - - - - - - - - -
        hmm_file_name = cs_dict['address_hmm']
        if (hmm_file_name == '(None)'):
          addr_hmm = None
        else:  # Initialise and load the address HMM
          addr_hmm = simplehmm.hmm('Address HMM', [], [])
          addr_hmm.load_hmm(hmm_file_name)
        hmm_train_name = cs_dict['hmm_train_file']
        if (hmm_train_name == '(None)'):
          hmm_train_name = None
        hmm_seq_name = cs_dict['hmm_seq_prob_file']
        if (hmm_seq_name == '(None)'):
          hmm_seq_name = None

        this_cs = standardisation.AddressStandardiser(input_fi = in_field_list,
                                                   output_fiel= out_field_list,
                                                   corr_list = this_corr_list,
                                                   tag_table = this_tag_table,
                                                   field_separator = field_sep,
                                                   check_wor = word_spill_flag,
                                                   address_hmm = addr_hmm,
                                                   hmm_train_f= hmm_train_name,
                                                   hmm_seq_pro = hmm_seq_name)

      elif (comp_std_type == 'Name'):  # - - - - - - - - - - - - - - - - - - -
        hmm_file_name = cs_dict['name_hmm']
        if (hmm_file_name == '(None)'):
          name_hmm = None
        else:  # Initialise and load the name HMM
          name_hmm = simplehmm.hmm('Name HMM', [], [])
          name_hmm.load_hmm(hmm_file_name)
        hmm_train_name = cs_dict['hmm_train_file']
        if (hmm_train_name == '(None)'):
          hmm_train_name = None
        hmm_seq_name = cs_dict['hmm_seq_prob_file']
        if (hmm_seq_name == '(None)'):
          hmm_seq_name = None

        first_name_comp =   cs_dict['first_name_comp']
        female_title_list = cs_dict['female_titles'].split(',')
        male_title_list =   cs_dict['male_titles'].split(',')

        this_cs = standardisation.NameStandardiser(input_fi = in_field_list,
                                                   output_fiel= out_field_list,
                                                   corr_list = this_corr_list,
                                                   tag_table = this_tag_table,
                                                   field_separator = field_sep,
                                                   check_wor = word_spill_flag,
                                                   name_hmm = name_hmm,
                                                   hmm_train_f= hmm_train_name,
                                                   hmm_seq_pro = hmm_seq_name,
                                                   first_nam = first_name_comp,
                                                   male_ti = male_title_list,
                                                   female_t= female_title_list)
      cs_list.append(this_cs)

    # Initialise the record standardiser - - - - - - - - - - - - - - - - - - -
    #
    progress_perc = 1  # Set to 1 for smooth progress bar updating

    pass_tuple_list = []
    for pass_field in pass_field_list:
      pass_tuple_list.append((pass_field, pass_field))

    self.init_progress_bar('Standardise data set:')

    rs = standardisation.RecordStandardiser(input_data = ids,
                                            output_data = ods,
                                            comp_stand_list = cs_list,
                                            progress_report = progress_perc,
                                            pass_field_list = pass_tuple_list,
                                          log_funct = self.update_progress_bar)
    rs.standardise()

    self.close_progress_bar()

    self.writeStatusBar('Standardised data set with %d records.' % \
                        (ids.num_records))

  # ---------------------------------------------------------------------------
  # Method to run code for a Febrl linkage or deduplication project
  #
  def runLinkCode(self):
    print 'Run linkage or deduplication code:', self.re_run

    # Check if all necessary components have been initialised
    #
    if ((self.febrl_code['data'] ==  None) or
        (self.febrl_code['indexing'] ==  None) or
        (self.febrl_code['comparison'] ==  None) or
        (self.febrl_code['classification'] ==  None) or
        (self.febrl_code['output'] ==  None)):
      self.messageDialog('Cannot run linkage or deduplication\n' + \
                         'project as not all required code\n' + \
                         'components are initialised.', 'error')
      return

    # Initialise first data set - - - - - - - - - - - - - - - - - - - - - - - -
    #
    if (self.re_run['data_init'] == True):  # Was modified, so have to re-do
      print 'Re-initialise data sets'

      data_set_info = self.data_set_info_list[0]  # Details for first data set

      if (self.data_set_type_list[0] in ['CSV','TAB']):

        field_val_list = []
        for field_num in range(len(data_set_info['field_names'])):
          field_val_list.append((data_set_info['field_names'][field_num],
                                 field_num))

        ds1 = dataset.DataSetCSV(access_mode="read",
                                 file_name = data_set_info['file_name'],
                                 delimiter = data_set_info['delimiter'],
                                 rec_ident = data_set_info['rec_id_field'],
                                 header_line = data_set_info['header_line'],
                                 strip_fields = data_set_info['strip_fields'],
                                 miss_val=data_set_info['miss_val'],
                                 field_list = field_val_list)

      elif (self.data_set_type_list[0] == 'COL'):

        header_line_flag = data_set_info['header_line']

        # Get width of the columns from GUI
        #
        col_width_vals = \
                       self.mainTree.get_widget('col_width_entry_a').get_text()
        col_width_vals = col_width_vals.strip()
        col_width_list = col_width_vals.split(',')

        field_val_list = []
        for field_num in range(len(data_set_info['field_names'])):
          if (header_line_flag == True):
            field_val_list.append(int(col_width_list[field_num]))
          else:
            field_val_list.append((field_names[field_num],
                                   int(col_width_list[field_num])))

        ds1 = dataset.DataSetCOL(access_mode="read",
                                 file_name = data_set_info['file_name'],
                                 rec_ident = data_set_info['rec_id_field'],
                                 header_line = data_set_info['header_line'],
                                 strip_fields = data_set_info['strip_fields'],
                                 miss_val=data_set_info['miss_val'],
                                 field_list = field_val_list)

      else:
        pass  # SQL to be implemented ########################################

      self.writeStatusBar('Initalised first data set with %d records.' % \
                          (ds1.num_records))

      if (self.project_type != 'Link'):
        ds2 = ds1  # Same data set

      else: # Initialise second data set - - - - - - - - - - - - - - - - - - -

        data_set_info = self.data_set_info_list[1]  # Details for 2nd data set

        if (self.data_set_type_list[1] in ['CSV','TAB']):

          field_val_list = []
          for field_num in range(len(data_set_info['field_names'])):
            field_val_list.append((data_set_info['field_names'][field_num],
                                   field_num))

          ds2 = dataset.DataSetCSV(access_mode="read",
                                   file_name = data_set_info['file_name'],
                                   delimiter = data_set_info['delimiter'],
                                   rec_ident = data_set_info['rec_id_field'],
                                   header_line = data_set_info['header_line'],
                                   strip_field = data_set_info['strip_fields'],
                                   miss_val=data_set_info['miss_val'],
                                   field_list = field_val_list)

        elif (self.data_set_type_list[1] == 'COL'):

          header_line_flag = data_set_info['header_line']

          # Get width of the columns from GUI
          #
          col_width_vals = \
                       self.mainTree.get_widget('col_width_entry_b').get_text()
          col_width_vals = col_width_vals.strip()
          col_width_list = col_width_vals.split(',')

          field_val_list = []
          for field_num in range(len(data_set_info['field_names'])):
            if (header_line_flag == True):
              field_val_list.append(int(col_width_list[field_num]))
            else:
              field_val_list.append((field_names[field_num],
                                     int(col_width_list[field_num])))

          ds2 = dataset.DataSetCOL(access_mode="read",
                                   file_name = data_set_info['file_name'],
                                   rec_ident = data_set_info['rec_id_field'],
                                   header_line = data_set_info['header_line'],
                                   strip_fiel = data_set_info['strip_fields'],
                                   miss_val=data_set_info['miss_val'],
                                   field_list = field_val_list)

        else:
          pass  # SQL to be implemented #######################################

        self.writeStatusBar('Initalised second data set with %d records.' % \
                            (ds2.num_records))

      self.data_sets = [ds1, ds2]  # Save for use with evaluation later

      self.re_run['data_init'] = False  # Can be re-used now

    else:  # Didn't have to re-initialise data sets as they have not changed
      ds1, ds2 = self.data_sets  # Make sure they are available for re-runs

      self.writeStatusBar('No need to re-initialised data set(s).')

    if (self.re_run['w_vec_generate'] == True):  # Re-generate weight vectors
      print 'Re-generate weight vectors'

      # Initialise field and record comparators - - - - - - - - - - - - - - - -
      #
      field_comp_list = []

      # Generate class for existing field comparison functions - - - - - - - -
      #
      for comp_funct_dict in self.field_comp_list:

        # Get common values from dictionary
        #
        comp_funct_name = comp_funct_dict['name']
        field_a_name =    comp_funct_dict['field_a_name']
        field_b_name =    comp_funct_dict['field_b_name']
        cache_field =     comp_funct_dict['cache_field']
        miss_weight =     float(comp_funct_dict['missing_weight'])
        agree_weight =    float(comp_funct_dict['agree_weight'])
        disagree_weight = float(comp_funct_dict['disagree_weight'])
        max_cache_size =  comp_funct_dict['max_cache_size']
        if (max_cache_size == 'None'):
          max_cache_size = None
        else:
          max_cache_size = int(max_cache_size)

        descr_str = comp_funct_dict['name']+'-'+field_a_name+'-'+field_b_name

        if (comp_funct_dict['name'] == 'Str-Exact'):
          cf = comparison.FieldComparatorExactString(do_caching = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight)

        elif (comp_funct_dict['name'] == 'Str-Contains'):
          cf = comparison.FieldComparatorContainsString(do_cachi = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight)

        elif (comp_funct_name == 'Str-Truncate'):
          cf = comparison.FieldComparatorTruncateString(do_cachi = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight,
              num_char_compared = int(comp_funct_dict['num_char_compared']))

        elif (comp_funct_name == 'Key-Diff'):
          cf = comparison.FieldComparatorKeyDiff(do_caching = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight,
              max_key_diff = int(comp_funct_dict['max_key_diff']))

        elif (comp_funct_name == 'Num-Perc'):
          cf = comparison.FieldComparatorNumericPerc(do_caching = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight,
              max_perc_diff = float(comp_funct_dict['max_perc_diff']))

        elif (comp_funct_name == 'Num-Abs'):
          cf = comparison.FieldComparatorNumericAbs(do_caching = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight,
              max_abs_diff = int(comp_funct_dict['max_abs_diff']))

        elif (comp_funct_name == 'Str-Encode'):
          encode_method_name = comp_funct_dict['encode_method']
          encode_method = self.stringencode_dict[encode_method_name]

          cf = comparison.FieldComparatorEncodeString(do_caching = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight,
              encode_method = encode_method,
              reverse = comp_funct_dict['reverse'],
              max_code_length = int(comp_funct_dict['max_code_length']))

        #elif (comp_name == 'FieldComparatorDistance'):  # TODO ###############

        elif (comp_funct_name == 'Date'):
          cf = comparison.FieldComparatorDate(do_caching = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight,
              date_format = comp_funct_dict['date_format'],
             max_day1_before_day2=int(comp_funct_dict['max_day1_before_day2']),
             max_day2_before_day1=int(comp_funct_dict['max_day2_before_day1']))

        elif (comp_funct_name == 'Time'):
          cf = comparison.FieldComparatorTime(do_caching = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight,
              day_start = comp_funct_dict['day_start'],
              max_time1_before_time2 = \
                                int(comp_funct_dict['max_time1_before_time2']),
              max_time2_before_time1 = \
                                int(comp_funct_dict['max_time2_before_time1']))

        elif (comp_funct_name == 'Age'):
          cf = comparison.FieldComparatorAge(do_caching = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight,
              date_format = comp_funct_dict['date_format'],
              fix_date = comp_funct_dict['fix_date'],
              max_perc_diff = float(comp_funct_dict['max_perc_diff']))

        elif (comp_funct_name == 'Jaro'):
          cf = comparison.FieldComparatorJaro(do_caching = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight,
              threshold = float(comp_funct_dict['threshold']))

        elif (comp_funct_name == 'Winkler'):
          cf = comparison.FieldComparatorWinkler(do_caching = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight,
              threshold = float(comp_funct_dict['threshold']),
              check_sim = comp_funct_dict['check_sim'],
              check_init = comp_funct_dict['check_init'],
              check_long = comp_funct_dict['check_long'])

        elif (comp_funct_name == 'Q-Gram'):
          cf = comparison.FieldComparatorQGram(do_caching = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight,
              threshold = float(comp_funct_dict['threshold']),
              q = int(comp_funct_dict['q']),
              common_divisor = comp_funct_dict['common_divisor'],
              padded = comp_funct_dict['padded'])

        elif (comp_funct_name == 'Pos-Q-Gram'):
          cf = comparison.FieldComparatorPosQGram(do_caching = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight,
              threshold = float(comp_funct_dict['threshold']),
              q = int(comp_funct_dict['q']),
              common_divisor = comp_funct_dict['common_divisor'],
              padded = comp_funct_dict['padded'],
              max_dist = int(comp_funct_dict['max_dist']))

        elif (comp_funct_name == 'S-Gram'):
          cf = comparison.FieldComparatorSGram(do_caching = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight,
              threshold = float(comp_funct_dict['threshold']),
              gram_class_list = eval(comp_funct_dict['gram_class_list']),
              common_divisor = comp_funct_dict['common_divisor'],
              padded = comp_funct_dict['padded'])

        elif (comp_funct_name == 'Edit-Dist'):
          cf = comparison.FieldComparatorEditDist(do_caching = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight,
              threshold = float(comp_funct_dict['threshold']))

        elif (comp_funct_name == 'Dam-Le-Edit-Dist'):
          cf = comparison.FieldComparatorDaLeDist(do_caching = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight,
              threshold = float(comp_funct_dict['threshold']))

        elif (comp_funct_name == 'Bag-Dist'):
          cf = comparison.FieldComparatorBagDist(do_caching = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight,
              threshold = float(comp_funct_dict['threshold']))

        elif (comp_funct_name == 'Smith-Water-Dist'):
          cf = comparison.FieldComparatorSWDist(do_caching = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight,
              threshold = float(comp_funct_dict['threshold']),
              common_divisor = comp_funct_dict['common_divisor'])

        elif (comp_funct_name == 'Syll-Align-Dist'):
          cf = comparison.FieldComparatorSyllAlDist(do_caching = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight,
              threshold = float(comp_funct_dict['threshold']),
              common_divisor = comp_funct_dict['common_divisor'],
              do_phonix = comp_funct_dict['do_phonix'])

        elif (comp_funct_name == 'Seq-Match'):
          cf = comparison.FieldComparatorSeqMatch(do_caching = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight,
              threshold = float(comp_funct_dict['threshold']))

        elif (comp_funct_name == 'Editex'):
          cf = comparison.FieldComparatorEditex(do_caching = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight,
              threshold = float(comp_funct_dict['threshold']))

        elif (comp_funct_name == 'Long-Common-Seq'):
          cf = comparison.FieldComparatorLCS(do_caching = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight,
              threshold = float(comp_funct_dict['threshold']),
              common_divisor = comp_funct_dict['common_divisor'],
              min_common_len = int(comp_funct_dict['min_common_len']))

        elif (comp_funct_name == 'Onto-LCS'):
          cf = comparison.FieldComparatorOntoLCS(do_caching = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight,
              threshold = float(comp_funct_dict['threshold']),
              common_divisor = comp_funct_dict['common_divisor'],
              min_common_len = int(comp_funct_dict['min_common_len']),
              p = float(comp_funct_dict['p']))

        elif (comp_funct_name == 'Token-Set'):
          stop_word_list = comp_funct_dict['stop_word_list'].split(',')
          stop_word_list_clean = []
          for stop_word in stop_word_list:
            if (stop_word.strip() != ''):
              stop_word_list_clean.append(stop_word)

          cf = comparison.FieldComparatorTokenSet(do_caching = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight,
              threshold = float(comp_funct_dict['threshold']),
              common_divisor = comp_funct_dict['common_divisor'],
              stop_word_list = stop_word_list_clean)

        elif (comp_funct_name == 'Compression'):
          cf = comparison.FieldComparatorCompress(do_caching = cache_field,
              description = descr_str,
              max_cache_size = max_cache_size,
              missing_weight = miss_weight,
              agree_weight = agree_weight,
              disagree_weight = disagree_weight,
              threshold = float(comp_funct_dict['threshold']),
              compressor = comp_funct_dict['compressor'])

        else:
          raise Exception, 'Illegal field comparison method: %s' % \
                           (comp_funct_name)

        # Append tuple to list of comparison functions
        #
        field_comp_list.append((cf, field_a_name, field_b_name))

      print 'field_comp_list:', field_comp_list ############################

      # Define the record comparator
      #
      rec_comp = comparison.RecordComparator(ds1, ds2, field_comp_list)

      self.writeStatusBar('Initalised field comparisons and record ' + \
                          'comparator.')

      # Initialise index - - - - - - - - - - - - - - - - - - - - - - - - - - -
      #
      index_def_list = []

      for index_list in self.index_def:
        print 'index_list:', index_list

        this_index_list = []

        for index_def_tuple in index_list:
          field_a_name =       index_def_tuple[0]
          field_b_name =       index_def_tuple[1]
          sort_words_flag =    index_def_tuple[2]
          reverse_flag =       index_def_tuple[3]
          max_len =            index_def_tuple[4]
          encode_funct =       index_def_tuple[5]
          encode_funct_param = index_def_tuple[6]

          index_tuple_list = [field_a_name, field_b_name, sort_words_flag,
                              reverse_flag]
          if (max_len in [None, 'None']):
            index_tuple_list.append(None)
          else:
            index_tuple_list.append(int(max_len))

          if (encode_funct in [None,'None']):
            index_tuple_list.append([])
          else:  # Process encoding function value and paramters
            if (encode_funct == 'Soundex'):
              encode_list = [encode.soundex]
            elif (encode_funct == 'Mod-Soundex'):
              encode_list = [encode.mod_soundex]
            elif (encode_funct == 'Phonex'):
              encode_list = [encode.phonex]
            elif (encode_funct == 'NYSIIS'):
              encode_list = [encode.nysiis]
            elif (encode_funct == 'Double-Metaphone'):
              encode_list = [encode.dmetaphone]
            elif (encode_funct == 'Phonix'):
              encode_list = [encode.phonix]
            elif (encode_funct == 'Fuzzy-Soundex'):
              encode_list = [encode.fuzzy_soundex]
            elif (encode_funct == 'Substring'):
              p1,p2 = encode_funct_param.split(',')
              encode_list = [encode.get_substring, int(p1), int(p2)]

            else:
              raise Exception, 'Illegal phonetic encoding function: %s' % \
                               (encode_funct)

            if ((encode_funct != 'Substring') and \
                (encode_funct_param not in ['None', None, ''])):
              encode_list.append(eval(encode_funct_param))
            index_tuple_list.append(encode_list)

          this_index_list.append(index_tuple_list)

        index_def_list.append(this_index_list)

      print 'index_def_list:', index_def_list ############################

      # Initialise the indexing method - - - - - - - - - - - - - - - - - - - -
      #
      progress_perc = 1  # Set to 1 for smooth progress bar updating

      if (self.output_dict['w_vec_file'][0] == False):
        w_vec_file = None
      else:
        w_vec_file = self.output_dict['w_vec_file'][1]

      index_method_name = self.index_method['name']
      index_sep_str =     self.index_method['index_sep_str']
      skip_missing =      self.index_method['skip_missing']

      if (index_method_name == 'FullIndex'):
        index = indexing.FullIndex(dataset1 = ds1,
                                 dataset2 = ds2,
                                 log_funct = self.update_progress_bar,
                                 rec_comparator = rec_comp,
                                 index_sep_str = index_sep_str,
                                 skip_missing = skip_missing,
                                 progress_report = progress_perc,
                                 weight_vec_file = w_vec_file,
                                 index_def = []) # No index definition required

      elif (index_method_name == 'BlockingIndex'):  # - - - - - - - - - - - - -

        if ((self.project_type == 'Link') and \
            (self.index_method['doBigMatch'] == True)):
          index = indexing.BigMatchIndex(dataset1 = ds1,
                                       dataset2 = ds2,
                                       log_funct = self.update_progress_bar,
                                       rec_comparator = rec_comp,
                                       index_sep_str = index_sep_str,
                                       skip_missing = skip_missing,
                                       progress_report = progress_perc,
                                       weight_vec_file = w_vec_file,
                                       index_def = index_def_list,
                                       block_method = ("block",))

        elif ((self.project_type != 'Link') and \
            (self.index_method['doDedup'] == True)):
          index = indexing.DedupIndex(dataset1 = ds1,
                                    dataset2 = ds2,
                                    log_funct = self.update_progress_bar,
                                    rec_comparator = rec_comp,
                                    index_sep_str = index_sep_str,
                                    skip_missing = skip_missing,
                                    progress_report = progress_perc,
                                    weight_vec_file = w_vec_file,
                                    index_def = index_def_list,
                                    block_method = ("block",))

        else:  # Use normal blocking index
          index = indexing.BlockingIndex(dataset1 = ds1,
                                       dataset2 = ds2,
                                       log_funct = self.update_progress_bar,
                                       rec_comparator = rec_comp,
                                       index_sep_str = index_sep_str,
                                       skip_missing = skip_missing,
                                       progress_report = progress_perc,
                                       weight_vec_file = w_vec_file,
                                       index_def = index_def_list)

      elif (index_method_name == 'SortingIndex'):  # - - - - - - - - - - - - -
        window_size = int(self.index_method['window_size'])

        if ((self.project_type == 'Link') and \
            (self.index_method['doBigMatch'] == True)):
          index = indexing.BigMatchIndex(dataset1 = ds1,
                                       dataset2 = ds2,
                                       log_funct = self.update_progress_bar,
                                       rec_comparator = rec_comp,
                                       index_sep_str = index_sep_str,
                                       skip_missing = skip_missing,
                                       progress_report = progress_perc,
                                       weight_vec_file = w_vec_file,
                                       index_def = index_def_list,
                                       block_method = ("sort", window_size))

        elif ((self.project_type != 'Link') and \
            (self.index_method['doDedup'] == True)):
          index = indexing.DedupIndex(dataset1 = ds1,
                                    dataset2 = ds2,
                                    log_funct = self.update_progress_bar,
                                    rec_comparator = rec_comp,
                                    index_sep_str = index_sep_str,
                                    skip_missing = skip_missing,
                                    progress_report = progress_perc,
                                    weight_vec_file = w_vec_file,
                                    index_def = index_def_list,
                                    block_method = ("sort", window_size))

        else:  # Use normal blocking index
          index = indexing.SortingIndex(dataset1 = ds1,
                                      dataset2 = ds2,
                                      log_funct = self.update_progress_bar,
                                      rec_comparator = rec_comp,
                                      index_sep_str = index_sep_str,
                                      skip_missing = skip_missing,
                                      progress_report = progress_perc,
                                      weight_vec_file = w_vec_file,
                                      index_def = index_def_list,
                                      window_size = window_size)

      elif (index_method_name == 'QGramIndex'):  # - - - - - - - - - - - - - -
        q =         int(self.index_method['q'])
        padded =    self.index_method['padded']
        threshold = float(self.index_method['threshold'])

        if ((self.project_type == 'Link') and \
            (self.index_method['doBigMatch'] == True)):
          index = indexing.BigMatchIndex(dataset1 = ds1,
                                       dataset2 = ds2,
                                       log_funct = self.update_progress_bar,
                                       rec_comparator = rec_comp,
                                       index_sep_str = index_sep_str,
                                       skip_missing = skip_missing,
                                       progress_report = progress_perc,
                                       weight_vec_file = w_vec_file,
                                       index_def = index_def_list,
                                       block_method = ("qgram", q, padded,
                                                       threshold))

        elif ((self.project_type != 'Link') and \
            (self.index_method['doDedup'] == True)):
          index = indexing.DedupIndex(dataset1 = ds1,
                                    dataset2 = ds2,
                                    log_funct = self.update_progress_bar,
                                    rec_comparator = rec_comp,
                                    index_sep_str = index_sep_str,
                                    skip_missing = skip_missing,
                                    progress_report = progress_perc,
                                    weight_vec_file = w_vec_file,
                                    index_def = index_def_list,
                                    block_method = ("qgram", q, padded,
                                                    threshold))

        else:  # Use normal blocking index
          index = indexing.QGramIndex(dataset1 = ds1,
                                    dataset2 = ds2,
                                    log_funct = self.update_progress_bar,
                                    rec_comparator = rec_comp,
                                    index_sep_str = index_sep_str,
                                    skip_missing = skip_missing,
                                    progress_report = progress_perc,
                                    weight_vec_file = w_vec_file,
                                    index_def = index_def_list,
                                    q = q,
                                    padded = padded,
                                    threshold = threshold)

      elif (index_method_name == 'CanopyIndex'):  # - - - - - - - - - - - - - -
        q =        int(self.index_method['q'])
        padded =   self.index_method['padded']
        del_perc = float(self.index_method['delete_perc'])

        p1,p2=self.index_method['canopy_method'][2].split(',') # Get parameters
        if (self.index_method['canopy_method'][1] == 'nearest'):
          p1,p2 = int(p1), int(p2)
        else:
          p1,p2 = float(p1), float(p2)
        canopy_method = (self.index_method['canopy_method'][0],
                        self.index_method['canopy_method'][1], p1, p2)

        index = indexing.CanopyIndex(dataset1 = ds1,
                                   dataset2 = ds2,
                                   log_funct = self.update_progress_bar,
                                   rec_comparator = rec_comp,
                                   index_sep_str = index_sep_str,
                                   skip_missing = skip_missing,
                                   progress_report = progress_perc,
                                   weight_vec_file = w_vec_file,
                                   index_def = index_def_list,
                                   q = q,
                                   padded = padded,
                                   delete_perc = del_perc,
                                   canopy_method = canopy_method)

      elif (index_method_name == 'StringMapIndex'):
        p1,p2=self.index_method['canopy_method'][1].split(',') # Get parameters
        if (self.index_method['canopy_method'][0] == 'nearest'):
          p1,p2 = int(p1), int(p2)
        else:
          p1,p2 = float(p1), float(p2)
        canopy_method = (self.index_method['canopy_method'][0], p1, p2)

        dim =             int(self.index_method['dim'])
        sub_dim =         int(self.index_method['sub_dim'])
        grid_resolution = int(self.index_method['grid_resolution'])
        cache_dist =      self.index_method['cache_dist']
        sim_funct_name =  self.index_method['sim_funct']
        sim_funct = {'Jaro':stringcmp.jaro, 'Winkler':stringcmp.winkler,
                     'Q-Gram':stringcmp.qgram, 'Pos-Q-Gram':stringcmp.posqgram,
                     'S-Gram':stringcmp.sgram, 'Edit-Dist':stringcmp.editdist,
                     'Mod-Edit-Dist':stringcmp.mod_editdist,
                     'Bag-Dist':stringcmp.bagdist,
                     'Smith-Water-Dist':stringcmp.swdist,
                     'Syll-Align-Dist':stringcmp.syllaligndist,
                     'Seq-Match':stringcmp.seqmatch,
                     'Long-Common-Seq':stringcmp.lcs,
                     'Onto-LCS':stringcmp.ontolcs, 'Editex':stringcmp.editex,
                     'Perm-Winkler':stringcmp.permwinkler,
                     'Sort-Winkler':stringcmp.sortwinkler}[sim_funct_name]

        index = indexing.StringMapIndex(dataset1 = ds1,
                                      dataset2 = ds2,
                                      log_funct = self.update_progress_bar,
                                      rec_comparator = rec_comp,
                                      index_sep_str = index_sep_str,
                                      skip_missing = skip_missing,
                                      progress_report = progress_perc,
                                      weight_vec_file = w_vec_file,
                                      index_def = index_def_list,
                                      canopy_method = canopy_method,
                                      grid_resolution = grid_resolution,
                                      dim = dim,
                                      sub_dim = sub_dim,
                                      cache_dist = cache_dist,
                                      sim_funct = sim_funct)

      elif (index_method_name == 'SuffixArrayIndex'):  # - - - - - - - - - - -
        suffix_method = self.index_method['suffix_method']
        padded =        self.index_method['padded']
        block_method =  (int(self.index_method['block_method'][0]),
                         int(self.index_method['block_method'][1]))

        index = indexing.SuffixArrayIndex(dataset1 = ds1,
                                        dataset2 = ds2,
                                        log_funct = self.update_progress_bar,
                                        rec_comparator = rec_comp,
                                        index_sep_str = index_sep_str,
                                        skip_missing = skip_missing,
                                        progress_report = progress_perc,
                                        index_def = index_def_list,
                                        padded = padded,
                                        suffix_method = suffix_method,
                                        block_method = block_method)

      else:
        raise Exception, 'Illegal indexing method: %s' % (index_method_name)

      print index

      self.writeStatusBar('Initalised index definitions and indexing method.')

    # Initialise classifier - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    classifier_method_name = self.classifier_method['name']

    if (classifier_method_name == 'FellegiSunter'):
      lower_thres = float(self.classifier_method['lower_threshold'])
      upper_thres = float(self.classifier_method['upper_threshold'])

      classifier = classification.FellegiSunter(lower_threshold = lower_thres,
                                                upper_threshold = upper_thres)

    elif (classifier_method_name == 'OptimalThreshold'):  # - - - - - - - - - -
      min_method = self.classifier_method['min_method']
      bin_width =  float(self.classifier_method['bin_width'])

      classifier = classification.OptimalThreshold(min_method = min_method,
                                                   bin_width = bin_width)

    elif (classifier_method_name == 'KMeans'):  #  - - - - - - - - - - - - - -
      max_iter_count = int(self.classifier_method['max_iter_count'])
      sample = self.classifier_method.get('sample', '')
      if (sample not in ['', None, 'None']):
        sample = int(sample)
      else:
        sample = None

      fuzz_reg_thres = self.classifier_method.get('fuzzy_reg_thres', '')
      if (fuzz_reg_thres not in ['', None, 'None']):
        fuzz_reg_thres = float(fuzz_reg_thres)
      else:
        fuzz_reg_thres = None

      dist_meas_name = self.classifier_method['dist_measure']
      dist_meas = {'Manhatten':mymath.distL1, 'Euclidean':mymath.distL2,
                   'L-Infinity':mymath.distLInf,
                   'Canberra':mymath.distCanberra}[dist_meas_name]
      centroid_init = self.classifier_method['centroid_init']

      classifier = classification.KMeans(max_iter_count=max_iter_count,
                                         dist_measure=dist_meas,
                                         centroid_init=centroid_init,
                                         sample = sample,
                                         fuzz_reg_thres = fuzz_reg_thres)

    elif (classifier_method_name == 'FarthestFirst'):  # - - - - - - - - - - -
      sample = self.classifier_method.get('sample', '')
      if (sample not in ['', None, 'None']):
        sample = int(sample)
      else:
        sample = None

      fuzz_reg_thres = self.classifier_method.get('fuzzy_reg_thres', '')
      if (fuzz_reg_thres not in ['', None, 'None']):
        fuzz_reg_thres = float(fuzz_reg_thres)
      else:
        fuzz_reg_thres = None

      dist_meas_name = self.classifier_method['dist_measure']
      dist_meas = {'Manhatten':mymath.distL1, 'Euclidean':mymath.distL2,
                   'L-Infinity':mymath.distLInf,
                   'Canberra':mymath.distCanberra}[dist_meas_name]
      centroid_init = self.classifier_method['centroid_init']

      classifier = classification.FarthestFirst(dist_measure=dist_meas,
                                                centroid_init=centroid_init,
                                                sample = sample,
                                                fuzz_reg_thres =fuzz_reg_thres)

    elif (classifier_method_name == 'SuppVecMachine'):  # - - - - - - - - - - -
      sample = self.classifier_method.get('sample', '')
      if (sample not in ['', None, 'None']):
        sample = int(sample)
      else:
        sample = None

      kernel_type = self.classifier_method['kernel_type']
      C =           float(self.classifier_method['C'])

      classifier = classification.SuppVecMachine(C = C,
                                                 kernel_type = kernel_type,
                                                 sample = sample)

    elif (classifier_method_name == 'TwoStep'):  # - - - - - - - - - - - - - -
      s1_m_method = (float(self.classifier_method['s1_match_method'][0]),
                     self.classifier_method['s1_match_method'][1])
      if (s1_m_method[-1] == 'nearest'):
        s1_m_method += (int(self.classifier_method['s1_match_method'][2]),
                        self.classifier_method['s1_match_method'][3])
      else:
        s1_m_method += (float(self.classifier_method['s1_match_method'][2]),)

      s1_nm_method = (float(self.classifier_method['s1_non_match_method'][0]),
                      self.classifier_method['s1_non_match_method'][1])
      if (s1_nm_method[-1] == 'nearest'):
        s1_nm_method += (int(self.classifier_method['s1_non_match_method'][2]),
                         self.classifier_method['s1_non_match_method'][3])
      else:
        s1_nm_method += \
                     (float(self.classifier_method['s1_non_match_method'][2]),)

      random_sel = self.classifier_method.get('random_selection', None)
      if (random_sel in ['', None, 'None']):
        random_method = None
      else:
        random_method = (random_sel[0], int(random_sel[1]), int(random_sel[2]))

      if (self.classifier_method['s2_classifier'][0] == 'kmeans'):
        max_iter_count = 1000 # Default value
        dist_meas_name = self.classifier_method['s2_classifier'][1]
        dist_meas = {'Manhatten':mymath.distL1, 'Euclidean':mymath.distL2,
                     'L-Infinity':mymath.distLInf,
                     'Canberra':mymath.distCanberra}[dist_meas_name]
        s2_classifier = (self.classifier_method['s2_classifier'][0],
                         dist_meas, max_iter_count)
      else: # SVM
        increment = 10  # Default value
        train_perc = 75  # Default value
        s2_classifier = (self.classifier_method['s2_classifier'][0],
                         self.classifier_method['s2_classifier'][1],
                         float(self.classifier_method['s2_classifier'][2]),
                         increment,
                         train_perc)

      classifier = classification.TwoStep(s1_match_method = s1_m_method,
                                          s1_non_match_method = s1_nm_method,
                                          random_selection = random_method,
                                          s2_classifier = s2_classifier)

    else:
      raise Exception, 'Illegal classifier method: %s' % \
                       (classifier_method_name)

    print classifier

    self.writeStatusBar('Initialised weight vector classifier.')

    # Now build, compact and run the indexing - - - - - - - - - - - - - - - - -
    #
    if (self.re_run['w_vec_generate'] == True):  # Re-generate weight vectors
      print 'Re-generate weight vectors (part 2: build, compact and run index)'

      self.init_progress_bar('Build index:')
      start_time = time.time()
      index.build()
      end_time_str = auxiliary.time_string(time.time() - start_time)
      self.writeStatusBar('Built index in %s.' % (end_time_str))
      self.close_progress_bar()

      self.init_progress_bar('Compact index:')
      start_time = time.time()
      index.compact()
      end_time_str = auxiliary.time_string(time.time() - start_time)
      self.writeStatusBar('Compacted index in %s.' % (end_time_str))
      self.close_progress_bar()

      self.init_progress_bar('Record pair comparisons:')

      # Check for output options
      #
      output_dict = self.output_dict
      if (output_dict['length_filter_perc'] != None):
        length_filter_perc = int(output_dict['length_filter_perc'])
      else:
        length_filter_perc = None
      if (output_dict['cut_off_threshold'] != None):
        cut_off_threshold = float(output_dict['cut_off_threshold'])
      else:
        cut_off_threshold = None

      if (output_dict['w_vec_file'][0] == True):
        index.run(length_filter_perc, cut_off_threshold)
        w_vec_file_name = output_dict['w_vec_file'][1]

        start_time = time.time()
        [field_names_list, w_vec_dict] = \
                                   output.LoadWeightVectorFile(w_vec_file_name)

      else:
        start_time = time.time()
        [field_names_list, w_vec_dict] = index.run(length_filter_perc,
                                                   cut_off_threshold)
      end_time_str = auxiliary.time_string(time.time() - start_time)
      self.writeStatusBar('Run record pair comparisons in %s.' % \
                          (end_time_str))

      self.close_progress_bar()

      self.addToLog('')  # Add information on record pair comparison to log
      self.addToLog('Compared %d record pairs in %s.' % (len(w_vec_dict),
                    end_time_str))

      self.w_vec_dict = w_vec_dict  # Save for later use in evaluation

      self.re_run['w_vec_generate'] = False  # Can be re-used now

    else:  # Didn't have to re-generate weight vectors
      w_vec_dict = self.w_vec_dict # Make sure they are available for re-runs
      output_dict = self.output_dict

      self.writeStatusBar('No need to re-generate weight vectors.')

    # Classify weight vectors according to classification method - - - - - - -
    #
    if (classifier_method_name in ['OptimalThreshold','SuppVecMachine']):
      print 'Train supervised classifier', classifier_method_name

      # For supervised classifier get true match status and then corresponding
      # true match and non-match sets
      #
      match_status_field = self.classifier_method.get('match_status_ind', -1)
      match_status_field_index = self.exact_comp_list[match_status_field][1]

      global MATCHSTATUSINDEX
      MATCHSTATUSINDEX = match_status_field_index

      def match_status_funct(red_id1, rec_id2, w_vec):
        return (w_vec[MATCHSTATUSINDEX] == 1.0)

      self.match_status_funct = match_status_funct  # Needed for evaluation

      tm_set, tnm_set = classification.get_true_matches_nonmatches(w_vec_dict,
                                                            match_status_funct)

      # Remove the match status fieldfrom weight vector dictionary
      #
      num_comp_funct = len(self.field_comp_list)

      manipulate_list = []
      for i in range(num_comp_funct):
        if (i != match_status_field_index):
          manipulate_list.append((i,))

      self.class_w_vec_dict = \
                classification.extract_collapse_weight_vectors(manipulate_list,
                                                               w_vec_dict)

      self.true_match_set = [tm_set, tnm_set]

    else:  # Unsupervised classification - - - - - - - - - - - - - - - - - - -
      print 'Train un-supervised classifier', classifier_method_name

      self.class_w_vec_dict = w_vec_dict
      self.match_status_funct = None
      tm_set =  set()  # Not needed for unsupervised training
      tnm_set = set()

    start_time = time.time()
    classifier.train(self.class_w_vec_dict, tm_set, tnm_set)
    end_time_str = auxiliary.time_string(time.time() - start_time)
    self.writeStatusBar('Trained classifier in %s.' % (end_time_str))

    # Classify weight vectors into matches, non-matches and possible matches
    #
    start_time = time.time()
    self.result_sets = classifier.classify(self.class_w_vec_dict)
    end_time_str = auxiliary.time_string(time.time() - start_time)

    self.addToLog('')  # Add information on classification to log page
    self.addToLog('Classified %d weight vectors:' % (len(w_vec_dict)))
    self.addToLog('  Number of classified matches:          %d' % \
                  (len(self.result_sets[0])))
    self.addToLog('  Number of classified non-matches:      %d' % \
                  (len(self.result_sets[1])))
    self.addToLog('  Number of classified possible matches: %d' % \
                  (len(self.result_sets[2])))
    self.addToLog('')

    # Generate and possible also save histogram - - - - - - - - - - - - - - - -
    #
    bin_width = float(len(self.field_comp_list)) / 10.0

    if (output_dict['histo_file'][0] == False):  # Don't save into file
      histo_list = output.GenerateHistogram(self.class_w_vec_dict, bin_width)
    else:
      histo_file_name = output_dict['histo_file'][1]
      histo_list = output.GenerateHistogram(self.class_w_vec_dict, bin_width,
                                            histo_file_name)

    self.addToLog('')  # Add histogram to log
    for line in histo_list:
      self.addToLog(line)
    self.addToLog('')

    # If selected, write output files - - - - - - - - - - - - - - - - - - - - -
    #
    match_set = self.result_sets[0]

    if (output_dict['m_status_file'][0] == True):
      match_status_file_name = output_dict['m_status_file'][1]
      output.SaveMatchStatusFile(self.class_w_vec_dict, match_set,
                                 match_status_file_name)
      self.addToLog('Wrote match status into file: %s' % \
                    (match_status_file_name))

    if (output_dict['m_datasets'][0] == True):
      out_ds1_name =  output_dict['m_datasets'][1][0]
      rec_id_field1 = output_dict['m_datasets'][1][1]

      if (self.project_type == 'Link'):
        out_ds2_name =  output_dict['m_datasets'][2][0]
        rec_id_field2 = output_dict['m_datasets'][2][1]

        output.SaveMatchDataSet(match_set, ds1, rec_id_field1, out_ds1_name,
                                ds2, rec_id_field2, out_ds2_name)
      else:
       output.SaveMatchDataSet(match_set, ds1, rec_id_field1, out_ds1_name)

    self.writeStatusBar('Classified weight vectors in %s, see Log page ' % \
                        (end_time_str) + 'for more information and ' + \
                        'matching weight histogram.')

    self.main_notebook_page_active_dict['Evaluate'] = True

    self.evaluation_dict['available'] = False  # So evaluation is re-calculated

    self.displayCurrentNotebookPage()  # Diplay the current page

  # ---------------------------------------------------------------------------
  # Method to initialise and show the progress bar window
  # (change the comment text in the progress window if one is provided)
  #
  def init_progress_bar(self, comment_str = ''):
    progrBarTree = gtk.glade.XML(self.glade_file_name, 'progress_bar_window')
    self.progress_bar_window = progrBarTree.get_widget('progress_bar_window')
    self.progress_bar_window.set_transient_for(self.mainWin)
    self.progress_bar_window.set_position(gtk.WIN_POS_CENTER_ON_PARENT)
    self.progress_bar_window.show()

    self.progress_bar_widget = progrBarTree.get_widget('progress_bar')
    self.progress_bar_widget.set_fraction(0.0)
    self.progress_bar_widget.show()

    if (comment_str != ''):
      self.progress_label_widget = progrBarTree.get_widget('progress_label')
      self.progress_label_widget.set_text('     '+comment_str+'     ')
      self.progress_label_widget.show()

  # ---------------------------------------------------------------------------
  # Method to close the progress bar
  #
  def close_progress_bar(self):
    self.progress_bar_window.hide()
    self.progress_bar_window.destroy()

  # ---------------------------------------------------------------------------
  # Method to update the progress bar, assumed to be called within run()
  # method of indexing, so progress percentage has to be extracted from log
  # text
  #
  def update_progress_bar(self, text):

    for s in text.split(' '):
      if ((len(s) > 2) and (s[0] == '(') and (s[-1] == ')')):
        perc = int(s[1:-2])
        frac = float(perc) / 100.0
        self.progress_bar_widget.set_fraction(frac)
        self.progress_bar_widget.set_text('%d%% done' % (perc))

        gtk.main_iteration(False)  # Two iterations to make sure progress is
        gtk.main_iteration(False)  # shown

        break  # Exit loop

    self.writeStatusBar(text)  # Also write text into status bar

  # ===========================================================================
  # Methods that handle Evaluate page events
  # ===========================================================================

  # ---------------------------------------------------------------------------
  # Display the Evaluate page
  #
  def evaluateView(self):
    print '  Switched to Evaluate page - Display this page'

    # Check that there are weight vectors to classify
    #
    if (len(self.class_w_vec_dict) == 0):
      self.messageDialog('No weight vectors were generated,\nplease check ' + \
                         'Log page.\n(Previous Evaluation page will be ' + \
                         'shown.)', 'error')
      return

    evaluation_dict = self.evaluation_dict  # Shortcut

    # Check if results are available or have to re-calculated - - - - - - - - -
    #
    if ((evaluation_dict['available'] == False) or \
        (sum(self.modified.values()) > 0)):

      [m_set, nm_set, pm_set] = self.result_sets  # Get classified match sets

      if (self.match_status_funct == None):  # Cannot get true match status

        class_m_w_list =  []  # Classified matches list
        class_nm_w_list = []  # Classified non-matches
        class_pm_w_list = []  # Classified possible matches

        min_w = 99999.99  # Global minimum and maximum summed weights needed
        max_w = -9999.99  # for histogram

        for (rec_id_tuple, w_vec) in self.class_w_vec_dict.iteritems():
          w_sum = sum(w_vec)
          min_w = min(min_w, w_sum)
          max_w = max(max_w, w_sum)

          if (rec_id_tuple in m_set):
            class_m_w_list.append(w_sum)
          elif (rec_id_tuple in nm_set):
            class_nm_w_list.append(w_sum)
          else:
            class_pm_w_list.append(w_sum)

        w_m_histo_dict =  {}
        w_nm_histo_dict = {}
        w_pm_histo_dict = {}

        bin_width = (max_w-min_w) / 20  # 20 bars in histogram
        if (bin_width == 0.0):  # Make sure it is not zero
          bin_width = 0.1  # Set to a default value

        half_bin_width = bin_width / 2.0

        for w_sum in class_m_w_list:
          binned_w = w_sum - (w_sum % bin_width) + half_bin_width
          w_count = w_m_histo_dict.get(binned_w, 0) + 1
          w_m_histo_dict[binned_w] = w_count
        for w_sum in class_nm_w_list:
          binned_w = w_sum - (w_sum % bin_width) + half_bin_width
          w_count = w_nm_histo_dict.get(binned_w, 0) + 1
          w_nm_histo_dict[binned_w] = w_count
        for w_sum in class_pm_w_list:
          binned_w = w_sum - (w_sum % bin_width) + half_bin_width
          w_count = w_pm_histo_dict.get(binned_w, 0) + 1
          w_pm_histo_dict[binned_w] = w_count

        evaluation_dict['class_m_weights'] =  w_m_histo_dict
        evaluation_dict['class_nm_weights'] = w_nm_histo_dict
        evaluation_dict['class_pm_weights'] = w_pm_histo_dict

      else:  # True match status is available - - - - - - - - - - - - - - - - -

        [tm_set, tnm_set] = self.true_match_set

        class_tp_w_list = []  # True matches classified matches (TP)
        class_fn_w_list = []  # True matches classified non-matches (FN)
        class_fp_w_list = []  # True non-matches classified matches (FP)
        class_tn_w_list = []  # True non-matches classified non-matches (TN)

        min_w = 99999.99  # Global minimum and maximum summed weights needed
        max_w = -9999.99  # for histogram

        for (rec_id_tuple, w_vec) in self.class_w_vec_dict.iteritems():
          w_sum = sum(w_vec)
          min_w = min(min_w, w_sum)
          max_w = max(max_w, w_sum)

          if ((rec_id_tuple in m_set) and (rec_id_tuple in tm_set)):
            class_tp_w_list.append(w_sum)
          elif ((rec_id_tuple in m_set) and (rec_id_tuple in tnm_set)):
            class_fp_w_list.append(w_sum)
          elif ((rec_id_tuple in nm_set) and (rec_id_tuple in tm_set)):
            class_fn_w_list.append(w_sum)
          elif ((rec_id_tuple in nm_set) and (rec_id_tuple in tnm_set)):
            class_tn_w_list.append(w_sum)

        w_tp_histo_dict = {}
        w_fn_histo_dict = {}
        w_fp_histo_dict = {}
        w_tn_histo_dict = {}

        bin_width = (max_w-min_w) / 20  # 20 bars in histogram
        if (bin_width == 0.0):  # Make sure it is not zero
          bin_width = 0.1  # Set to a default value

        half_bin_width = bin_width / 2.0

        for w_sum in class_tp_w_list:
          binned_w = w_sum - (w_sum % bin_width) + half_bin_width
          w_count = w_tp_histo_dict.get(binned_w, 0) + 1
          w_tp_histo_dict[binned_w] = w_count
        for w_sum in class_fn_w_list:
          binned_w = w_sum - (w_sum % bin_width) + half_bin_width
          w_count = w_fn_histo_dict.get(binned_w, 0) + 1
          w_fn_histo_dict[binned_w] = w_count
        for w_sum in class_fp_w_list:
          binned_w = w_sum - (w_sum % bin_width) + half_bin_width
          w_count = w_fp_histo_dict.get(binned_w, 0) + 1
          w_fp_histo_dict[binned_w] = w_count
        for w_sum in class_tn_w_list:
          binned_w = w_sum - (w_sum % bin_width) + half_bin_width
          w_count = w_tn_histo_dict.get(binned_w, 0) + 1
          w_tn_histo_dict[binned_w] = w_count

        evaluation_dict['class_tp_weights'] = w_tp_histo_dict
        evaluation_dict['class_fn_weights'] = w_fn_histo_dict
        evaluation_dict['class_fp_weights'] = w_fp_histo_dict
        evaluation_dict['class_tn_weights'] = w_tn_histo_dict

        self.addToLog('')
        self.addToLog('Classification results:')
        self.addToLog(' True matches (TP):      %d' % \
                      (sum(w_tp_histo_dict.values())))
        self.addToLog(' False matches (FP):     %d' % \
                      (sum(w_fp_histo_dict.values())))
        self.addToLog(' True non-matches (TN):  %d' % \
                      (sum(w_tn_histo_dict.values())))
        self.addToLog(' False non-matches (FN): %d' % \
                      (sum(w_fn_histo_dict.values())))
        self.addToLog('')

      self.bar_width = (max_w-min_w) / 30  # For histoplots
      if (self.bar_width == 0.0):  #  Make sure its not an empty bar
        self.bar_width = 0.1

      evaluation_dict['available'] = True  # Match results now available

    # Get the page widget and its two filled children - - - - - - - - - - - - -
    #
    evaluate_page_box = self.mainTree.get_widget('evaluate_page_box')
    evaluate_measures_box = self.mainTree.get_widget('evaluate_measures_box')

    for child in evaluate_page_box.get_children():
      evaluate_page_box.remove(child)

    # Check if Matplot lib is available - - - - - - - - - - - - - - - - - - - -
    #
    if (imp_matplot == True):

      # Generate the Matplotlib graph - - - - - - - - - - - - - - - - - - - - -
      #
      histo_figure = matplotlib.figure.Figure()  #figsize=(6,4), dpi=72)
      histo_axis = histo_figure.add_subplot(111)
      histo_axis.set_xlabel('Matching weight')
      histo_axis.set_ylabel('Counts')
      histo_axis.set_title('Summed matching weights histogram')
      histo_axis.grid(True)

      # Check if either only classified or true matches are available
      #
      if (self.match_status_funct == None):

        # Make a list of all x values over all sets
        #
        x_vals = evaluation_dict['class_m_weights'].keys()
        for x in (evaluation_dict['class_nm_weights'].keys() + \
                  evaluation_dict['class_pm_weights'].keys()):
          if x not in x_vals:
            x_vals.append(x)
        x_vals.sort()

        if (len(x_vals) == 1):  # Only one x value
          x_vals = [x_vals[0]-1, x_vals[0], x_vals[0]+1]

        m_y_vals =  []
        nm_y_vals = []
        pm_y_vals = []

        for x in x_vals:
          m_y_vals.append(evaluation_dict['class_m_weights'].get(x,0))
          nm_y_vals.append(evaluation_dict['class_nm_weights'].get(x,0))
          pm_y_vals.append(evaluation_dict['class_pm_weights'].get(x,0))

        p1 = histo_axis.bar(x_vals, m_y_vals, self.bar_width, color='r')
        p2 = histo_axis.bar(x_vals, nm_y_vals, self.bar_width, color='b',
                            bottom = m_y_vals)

        if (len(evaluation_dict['class_pm_weights']) > 0):
          pm_bottom_vals = []
          for i in range(len(x_vals)):
            pm_bottom_vals.append(m_y_vals[i]+nm_y_vals[i])

          p3 = histo_axis.bar(x_vals, pm_y_vals, self.bar_width, color='y',
                              bottom = pm_bottom_vals)
          histo_axis.legend((p1[0], p2[0], p3[0]), ("Matches", "Non-matches",
                             "Possible matches"), shadow = False)
        else:  # No possible matches
          histo_axis.legend((p1[0], p2[0]), ("Matches", "Non-matches"),
                            shadow = False)

      else:  # True match status available - - - - - - - - - - - - - - - - - -

        # Make a list of all x values over all sets
        #
        x_vals = evaluation_dict['class_tp_weights'].keys()
        for x in (evaluation_dict['class_fn_weights'].keys() + \
                  evaluation_dict['class_fp_weights'].keys() + \
                  evaluation_dict['class_tn_weights'].keys()):
          if x not in x_vals:
            x_vals.append(x)
        x_vals.sort()

        tp_y_vals = []
        fn_y_vals = []
        fp_y_vals = []
        tn_y_vals = []

        for x in x_vals:
          tp_y_vals.append(evaluation_dict['class_tp_weights'].get(x,0))
          fp_y_vals.append(evaluation_dict['class_fp_weights'].get(x,0))
          fn_y_vals.append(evaluation_dict['class_fn_weights'].get(x,0))
          tn_y_vals.append(evaluation_dict['class_tn_weights'].get(x,0))

        p1 = histo_axis.bar(x_vals, tp_y_vals, self.bar_width, color='r')
        p2 = histo_axis.bar(x_vals, fn_y_vals, self.bar_width, color='g',
                            bottom = tp_y_vals)

        fp_bottom_vals = []
        tn_bottom_vals = []
        for i in range(len(x_vals)):
          fp_bottom_vals.append(tp_y_vals[i] + fn_y_vals[i])
          tn_bottom_vals.append(tp_y_vals[i] + fn_y_vals[i] + fp_y_vals[i])

        p3 = histo_axis.bar(x_vals, fp_y_vals, self.bar_width, color='y',
                            bottom = fp_bottom_vals)
        p4 = histo_axis.bar(x_vals, tn_y_vals, self.bar_width, color='b',
                            bottom = tn_bottom_vals)

        histo_axis.legend((p1[0], p2[0], p3[0], p4[0]), ("True Matches (TP)",
                           "False Non-matches (FN)", "False matches (FP)",
                           "True non-matches (TN)"), shadow = False)

      # Make it a gtk.DrawingArea
      #
      histo_canvas = \
                  matplotlib.backends.backend_gtk.FigureCanvasGTK(histo_figure)

    # Calculate quality and complexity measure results - - - - - - - - - - - -
    #
    class_w_vec_dict = self.class_w_vec_dict
    [ds1, ds2] =       self.data_sets
    m_status_funct =   self.match_status_funct

    # If no match status available only reduction ratio can be calculated
    #
    evaluation_dict['rr'] = measurements.reduction_ratio(class_w_vec_dict,
                                                         ds1, ds2)

    if (self.match_status_funct != None):
      w_vec_dict = self.w_vec_dict
      evaluation_dict['pq'] = measurements.pairs_quality(w_vec_dict,
                                                         m_status_funct)

      #TODO: #############################
      # evaluation_dict['pq'] = \
      #measurements.pairs_completeness(w_vec_dict, ds1, ds2,
      #  **get_id_funct**, m_check_funct):

      acc, prec, reca, fmeas = measurements.quality_measures(w_vec_dict,
                                                             m_set, nm_set,
                                                             m_status_funct)
      evaluation_dict['acc'] =   acc
      evaluation_dict['prec'] =  prec
      evaluation_dict['reca'] =  reca
      evaluation_dict['fmeas'] = fmeas

    # If there are evaluation values put them into graph and text result fields
    #
    for measure in ['acc', 'prec', 'reca', 'fmeas', 'rr', 'pc', 'pq']:
      measure_label_widget = self.mainTree.get_widget('%s_val_label' % \
                                                     (measure))
      if (measure in evaluation_dict):
        measure_str = '%.3f' % (evaluation_dict[measure])
      else:
        measure_str = '---'
      measure_label_widget.set_text(measure_str)

    # Re-generate the Evaluate box - - - - - - - - - - - - - - - - - - - - - -
    #
    if (imp_matplot == True):
      evaluate_page_box.pack_start(histo_canvas, True, True, 5)
      histo_canvas.show()
    else:
      matplot_label = gtk.Label('Matplotlib not available - sorry no ' + \
                                'graphical histogram output. See Log page ' + \
                                'for textual histogram.')
      evaluate_page_box.pack_start(matplot_label, True, True, 5)
      horiz_separator = gtk.HSeparator()
      evaluate_page_box.pack_start(horiz_separator, False, False, 0)
      matplot_label.show()
      horiz_separator.show()

    evaluate_page_box.pack_start(evaluate_measures_box, False, False, 5)
    evaluate_measures_box.show()
    evaluate_page_box.show()

    self.writeStatusBar('Calculated quality and complexity measures.')

  def evaluateExecute(self):  # An Execute on the Evaluate page ---------------
    self.addToLog('Execute on Evaluate page')


  # ===========================================================================
  # Methods that handle Review page events
  # ===========================================================================

  def reviewView(self):  # A switch to the Review page ------------------------
    print '  Switched to Review page - Display this page'


  def reviewExecute(self):  # An Execute on the Review page -------------------
    self.addToLog('Execute on Review page')


  # ===========================================================================
  # Methods that handle Log page events
  # ===========================================================================

  def logView(self):  # A switch to the Log page ------------------------------
    print '  Switched to Log page - Display this page'

    self.log_page_buffer.set_text(self.log_page_text)  # Print what is in log

    # Create a mark at the end of the text and then scroll down to it
    #
    mark = self.log_page_buffer.create_mark('end',
                                    self.log_page_buffer.get_end_iter(), False)
    self.log_page_textview.scroll_to_mark(mark, 0.05, True, 0.0, 1.0)


  # ===========================================================================
  # Auxilliary methods
  # ===========================================================================

  def addToLog(self, text):
    self.log_page_text += text + os.linesep

  def str_is_int(self, x):
    if (not x.isdigit()):
      return False
    try:
      i = int(x)
    except:
      return False
    return True

  def str_is_pos_int(self, x):
    if (not self.str_is_int(x)):
      return False
    if (int(x) < 1):
      return False
    else:
      return True

  def str_is_not_neg_int(self, x):
    if (not self.str_is_int(x)):
      return False
    if (int(x) < 0):
      return False
    else:
      return True

  def str_is_normalised(self, x):
    try:
      t = float(x)
    except:
      return False
    if ((t < 0) or (t > 1)):
      return False
    else:
      return True

  def str_is_normalised_not_zero(self, x):
    if (not self.str_is_normalised(x)):
      return False
    if (float(x) == 0.0):
      return False
    else:
      return True

  def str_is_percentage(self, x):
    try:
      t = float(x)
    except:
      return False
    if ((t < 0) or (t > 100)):
      return False
    else:
      return True

  def str_is_percentage_not_zero(self, x):
    if (not self.str_is_percentage(x)):
      return False
    if (float(x) == 0.0):
      return False
    else:
      return True

  def str_is_float(self, x):
    try:
      t = float(x)
    except:
      return False
    return True

  def str_is_pos_float(self, x):
    try:
      t = float(x)
    except:
      return False
    if (t < 0.0):
      return False
    return True


# =============================================================================
# Start the GUI

app = MainFebrlWindow()
gtk.main()
