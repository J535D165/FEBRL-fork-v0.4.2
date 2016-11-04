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
# The Original Software is: "evalIndexing.csh"
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

# Shell script to call evalIndexing.py
# Arguments are: [data set] [index method]

python evalIndexing.py census blocking
python evalIndexing.py census sorted-inv-index
python evalIndexing.py census sorted-array
python evalIndexing.py census adapt-sorted
python evalIndexing.py census suffix-array
python evalIndexing.py census suffix-array-substr
python evalIndexing.py census robust-suffix-array
python evalIndexing.py census q-gram
python evalIndexing.py census canopy-th
python evalIndexing.py census canopy-nn
python evalIndexing.py census string-map-th
python evalIndexing.py census string-map-nn


python evalIndexing.py cora blocking
python evalIndexing.py cora sorted-inv-index
python evalIndexing.py cora sorted-array
python evalIndexing.py cora adapt-sorted
python evalIndexing.py cora suffix-array
python evalIndexing.py cora suffix-array-substr
python evalIndexing.py cora robust-suffix-array
python evalIndexing.py cora q-gram
python evalIndexing.py cora canopy-th
python evalIndexing.py cora canopy-nn
python evalIndexing.py cora string-map-th
python evalIndexing.py cora string-map-nn


python evalIndexing.py rest blocking
python evalIndexing.py rest sorted-inv-index
python evalIndexing.py rest sorted-array
python evalIndexing.py rest adapt-sorted
python evalIndexing.py rest suffix-array
python evalIndexing.py rest suffix-array-substr
python evalIndexing.py rest robust-suffix-array
python evalIndexing.py rest q-gram
python evalIndexing.py rest canopy-th
python evalIndexing.py rest canopy-nn
python evalIndexing.py rest string-map-th
python evalIndexing.py rest string-map-nn


python evalIndexing.py cddb blocking
python evalIndexing.py cddb sorted-inv-index
python evalIndexing.py cddb sorted-array
python evalIndexing.py cddb adapt-sorted
python evalIndexing.py cddb suffix-array
python evalIndexing.py cddb suffix-array-substr
python evalIndexing.py cddb robust-suffix-array
python evalIndexing.py cddb q-gram
python evalIndexing.py cddb canopy-th
python evalIndexing.py cddb canopy-nn
python evalIndexing.py cddb string-map-th
python evalIndexing.py cddb string-map-nn

# -------

python evalIndexing.py B_1000 blocking
python evalIndexing.py B_1000 sorted-inv-index
python evalIndexing.py B_1000 sorted-array
python evalIndexing.py B_1000 adapt-sorted
python evalIndexing.py B_1000 suffix-array
python evalIndexing.py B_1000 suffix-array-substr
python evalIndexing.py B_1000 robust-suffix-array
python evalIndexing.py B_1000 q-gram
python evalIndexing.py B_1000 canopy-th
python evalIndexing.py B_1000 canopy-nn
python evalIndexing.py B_1000 string-map-th
python evalIndexing.py B_1000 string-map-nn

python evalIndexing.py B_2500 blocking
python evalIndexing.py B_2500 sorted-inv-index
python evalIndexing.py B_2500 sorted-array
python evalIndexing.py B_2500 adapt-sorted
python evalIndexing.py B_2500 suffix-array
python evalIndexing.py B_2500 suffix-array-substr
python evalIndexing.py B_2500 robust-suffix-array
python evalIndexing.py B_2500 q-gram
python evalIndexing.py B_2500 canopy-th
python evalIndexing.py B_2500 canopy-nn
python evalIndexing.py B_2500 string-map-th
python evalIndexing.py B_2500 string-map-nn

python evalIndexing.py B_5000 blocking
python evalIndexing.py B_5000 sorted-inv-index
python evalIndexing.py B_5000 sorted-array
python evalIndexing.py B_5000 adapt-sorted
python evalIndexing.py B_5000 suffix-array
python evalIndexing.py B_5000 suffix-array-substr
python evalIndexing.py B_5000 robust-suffix-array
python evalIndexing.py B_5000 q-gram
python evalIndexing.py B_5000 canopy-th
python evalIndexing.py B_5000 canopy-nn
python evalIndexing.py B_5000 string-map-th
python evalIndexing.py B_5000 string-map-nn

python evalIndexing.py B_10000 blocking
python evalIndexing.py B_10000 sorted-inv-index
python evalIndexing.py B_10000 sorted-array
python evalIndexing.py B_10000 adapt-sorted
python evalIndexing.py B_10000 suffix-array
python evalIndexing.py B_10000 suffix-array-substr
python evalIndexing.py B_10000 robust-suffix-array
python evalIndexing.py B_10000 q-gram
python evalIndexing.py B_10000 canopy-th
python evalIndexing.py B_10000 canopy-nn
python evalIndexing.py B_10000 string-map-th
python evalIndexing.py B_10000 string-map-nn

python evalIndexing.py B_25000 blocking
python evalIndexing.py B_25000 sorted-inv-index
python evalIndexing.py B_25000 sorted-array
python evalIndexing.py B_25000 adapt-sorted
python evalIndexing.py B_25000 suffix-array
python evalIndexing.py B_25000 suffix-array-substr
python evalIndexing.py B_25000 robust-suffix-array
python evalIndexing.py B_25000 q-gram
python evalIndexing.py B_25000 canopy-th
python evalIndexing.py B_25000 canopy-nn
python evalIndexing.py B_25000 string-map-th
python evalIndexing.py B_25000 string-map-nn

python evalIndexing.py B_50000 blocking
python evalIndexing.py B_50000 sorted-inv-index
python evalIndexing.py B_50000 sorted-array
python evalIndexing.py B_50000 adapt-sorted
python evalIndexing.py B_50000 suffix-array
python evalIndexing.py B_50000 suffix-array-substr
python evalIndexing.py B_50000 robust-suffix-array
python evalIndexing.py B_50000 q-gram
python evalIndexing.py B_50000 canopy-th
python evalIndexing.py B_50000 canopy-nn
python evalIndexing.py B_50000 string-map-th
python evalIndexing.py B_50000 string-map-nn

# -------

python evalIndexing.py C_1000 blocking
python evalIndexing.py C_1000 sorted-inv-index
python evalIndexing.py C_1000 sorted-array
python evalIndexing.py C_1000 adapt-sorted
python evalIndexing.py C_1000 suffix-array
python evalIndexing.py C_1000 suffix-array-substr
python evalIndexing.py C_1000 robust-suffix-array
python evalIndexing.py C_1000 q-gram
python evalIndexing.py C_1000 canopy-th
python evalIndexing.py C_1000 canopy-nn
python evalIndexing.py C_1000 string-map-th
python evalIndexing.py C_1000 string-map-nn

python evalIndexing.py C_2500 blocking
python evalIndexing.py C_2500 sorted-inv-index
python evalIndexing.py C_2500 sorted-array
python evalIndexing.py C_2500 adapt-sorted
python evalIndexing.py C_2500 suffix-array
python evalIndexing.py C_2500 suffix-array-substr
python evalIndexing.py C_2500 robust-suffix-array
python evalIndexing.py C_2500 q-gram
python evalIndexing.py C_2500 canopy-th
python evalIndexing.py C_2500 canopy-nn
python evalIndexing.py C_2500 string-map-th
python evalIndexing.py C_2500 string-map-nn

python evalIndexing.py C_5000 blocking
python evalIndexing.py C_5000 sorted-inv-index
python evalIndexing.py C_5000 sorted-array
python evalIndexing.py C_5000 adapt-sorted
python evalIndexing.py C_5000 suffix-array
python evalIndexing.py C_5000 suffix-array-substr
python evalIndexing.py C_5000 robust-suffix-array
python evalIndexing.py C_5000 q-gram
python evalIndexing.py C_5000 canopy-th
python evalIndexing.py C_5000 canopy-nn
python evalIndexing.py C_5000 string-map-th
python evalIndexing.py C_5000 string-map-nn

python evalIndexing.py C_10000 blocking
python evalIndexing.py C_10000 sorted-inv-index
python evalIndexing.py C_10000 sorted-array
python evalIndexing.py C_10000 adapt-sorted
python evalIndexing.py C_10000 suffix-array
python evalIndexing.py C_10000 suffix-array-substr
python evalIndexing.py C_10000 robust-suffix-array
python evalIndexing.py C_10000 q-gram
python evalIndexing.py C_10000 canopy-th
python evalIndexing.py C_10000 canopy-nn
python evalIndexing.py C_10000 string-map-th
python evalIndexing.py C_10000 string-map-nn

python evalIndexing.py C_25000 blocking
python evalIndexing.py C_25000 sorted-inv-index
python evalIndexing.py C_25000 sorted-array
python evalIndexing.py C_25000 adapt-sorted
python evalIndexing.py C_25000 suffix-array
python evalIndexing.py C_25000 suffix-array-substr
python evalIndexing.py C_25000 robust-suffix-array
python evalIndexing.py C_25000 q-gram
python evalIndexing.py C_25000 canopy-th
python evalIndexing.py C_25000 canopy-nn
python evalIndexing.py C_25000 string-map-th
python evalIndexing.py C_25000 string-map-nn

python evalIndexing.py C_50000 blocking
python evalIndexing.py C_50000 sorted-inv-index
python evalIndexing.py C_50000 sorted-array
python evalIndexing.py C_50000 adapt-sorted
python evalIndexing.py C_50000 suffix-array
python evalIndexing.py C_50000 suffix-array-substr
python evalIndexing.py C_50000 robust-suffix-array
python evalIndexing.py C_50000 q-gram
python evalIndexing.py C_50000 canopy-th
python evalIndexing.py C_50000 canopy-nn
python evalIndexing.py C_50000 string-map-th
python evalIndexing.py C_50000 string-map-nn

# -------

python evalIndexing.py E_1000 blocking
python evalIndexing.py E_1000 sorted-inv-index
python evalIndexing.py E_1000 sorted-array
python evalIndexing.py E_1000 adapt-sorted
python evalIndexing.py E_1000 suffix-array
python evalIndexing.py E_1000 suffix-array-substr
python evalIndexing.py E_1000 robust-suffix-array
python evalIndexing.py E_1000 q-gram
python evalIndexing.py E_1000 canopy-th
python evalIndexing.py E_1000 canopy-nn
python evalIndexing.py E_1000 string-map-th
python evalIndexing.py E_1000 string-map-nn

python evalIndexing.py E_2500 blocking
python evalIndexing.py E_2500 sorted-inv-index
python evalIndexing.py E_2500 sorted-array
python evalIndexing.py E_2500 adapt-sorted
python evalIndexing.py E_2500 suffix-array
python evalIndexing.py E_2500 suffix-array-substr
python evalIndexing.py E_2500 robust-suffix-array
python evalIndexing.py E_2500 q-gram
python evalIndexing.py E_2500 canopy-th
python evalIndexing.py E_2500 canopy-nn
python evalIndexing.py E_2500 string-map-th
python evalIndexing.py E_2500 string-map-nn

python evalIndexing.py E_5000 blocking
python evalIndexing.py E_5000 sorted-inv-index
python evalIndexing.py E_5000 sorted-array
python evalIndexing.py E_5000 adapt-sorted
python evalIndexing.py E_5000 suffix-array
python evalIndexing.py E_5000 suffix-array-substr
python evalIndexing.py E_5000 robust-suffix-array
python evalIndexing.py E_5000 q-gram
python evalIndexing.py E_5000 canopy-th
python evalIndexing.py E_5000 canopy-nn
python evalIndexing.py E_5000 string-map-th
python evalIndexing.py E_5000 string-map-nn

python evalIndexing.py E_10000 blocking
python evalIndexing.py E_10000 sorted-inv-index
python evalIndexing.py E_10000 sorted-array
python evalIndexing.py E_10000 adapt-sorted
python evalIndexing.py E_10000 suffix-array
python evalIndexing.py E_10000 suffix-array-substr
python evalIndexing.py E_10000 robust-suffix-array
python evalIndexing.py E_10000 q-gram
python evalIndexing.py E_10000 canopy-th
python evalIndexing.py E_10000 canopy-nn
python evalIndexing.py E_10000 string-map-th
python evalIndexing.py E_10000 string-map-nn

python evalIndexing.py E_25000 blocking
python evalIndexing.py E_25000 sorted-inv-index
python evalIndexing.py E_25000 sorted-array
python evalIndexing.py E_25000 adapt-sorted
python evalIndexing.py E_25000 suffix-array
python evalIndexing.py E_25000 suffix-array-substr
python evalIndexing.py E_25000 robust-suffix-array
python evalIndexing.py E_25000 q-gram
python evalIndexing.py E_25000 canopy-th
python evalIndexing.py E_25000 canopy-nn
python evalIndexing.py E_25000 string-map-th
python evalIndexing.py E_25000 string-map-nn

python evalIndexing.py E_50000 blocking
python evalIndexing.py E_50000 sorted-inv-index
python evalIndexing.py E_50000 sorted-array
python evalIndexing.py E_50000 adapt-sorted
python evalIndexing.py E_50000 suffix-array
python evalIndexing.py E_50000 suffix-array-substr
python evalIndexing.py E_50000 robust-suffix-array
python evalIndexing.py E_50000 q-gram
python evalIndexing.py E_50000 canopy-th
python evalIndexing.py E_50000 canopy-nn
python evalIndexing.py E_50000 string-map-th
python evalIndexing.py E_50000 string-map-nn

# -------

python evalIndexing.py F_1000 blocking
python evalIndexing.py F_1000 sorted-inv-index
python evalIndexing.py F_1000 sorted-array
python evalIndexing.py F_1000 adapt-sorted
python evalIndexing.py F_1000 suffix-array
python evalIndexing.py F_1000 suffix-array-substr
python evalIndexing.py F_1000 robust-suffix-array
python evalIndexing.py F_1000 q-gram
python evalIndexing.py F_1000 canopy-th
python evalIndexing.py F_1000 canopy-nn
python evalIndexing.py F_1000 string-map-th
python evalIndexing.py F_1000 string-map-nn

python evalIndexing.py F_2500 blocking
python evalIndexing.py F_2500 sorted-inv-index
python evalIndexing.py F_2500 sorted-array
python evalIndexing.py F_2500 adapt-sorted
python evalIndexing.py F_2500 suffix-array
python evalIndexing.py F_2500 suffix-array-substr
python evalIndexing.py F_2500 robust-suffix-array
python evalIndexing.py F_2500 q-gram
python evalIndexing.py F_2500 canopy-th
python evalIndexing.py F_2500 canopy-nn
python evalIndexing.py F_2500 string-map-th
python evalIndexing.py F_2500 string-map-nn

python evalIndexing.py F_5000 blocking
python evalIndexing.py F_5000 sorted-inv-index
python evalIndexing.py F_5000 sorted-array
python evalIndexing.py F_5000 adapt-sorted
python evalIndexing.py F_5000 suffix-array
python evalIndexing.py F_5000 suffix-array-substr
python evalIndexing.py F_5000 robust-suffix-array
python evalIndexing.py F_5000 q-gram
python evalIndexing.py F_5000 canopy-th
python evalIndexing.py F_5000 canopy-nn
python evalIndexing.py F_5000 string-map-th
python evalIndexing.py F_5000 string-map-nn

python evalIndexing.py F_10000 blocking
python evalIndexing.py F_10000 sorted-inv-index
python evalIndexing.py F_10000 sorted-array
python evalIndexing.py F_10000 adapt-sorted
python evalIndexing.py F_10000 suffix-array
python evalIndexing.py F_10000 suffix-array-substr
python evalIndexing.py F_10000 robust-suffix-array
python evalIndexing.py F_10000 q-gram
python evalIndexing.py F_10000 canopy-th
python evalIndexing.py F_10000 canopy-nn
python evalIndexing.py F_10000 string-map-th
python evalIndexing.py F_10000 string-map-nn

python evalIndexing.py F_25000 blocking
python evalIndexing.py F_25000 sorted-inv-index
python evalIndexing.py F_25000 sorted-array
python evalIndexing.py F_25000 adapt-sorted
python evalIndexing.py F_25000 suffix-array
python evalIndexing.py F_25000 suffix-array-substr
python evalIndexing.py F_25000 robust-suffix-array
python evalIndexing.py F_25000 q-gram
python evalIndexing.py F_25000 canopy-th
python evalIndexing.py F_25000 canopy-nn
python evalIndexing.py F_25000 string-map-th
python evalIndexing.py F_25000 string-map-nn

python evalIndexing.py F_50000 blocking
python evalIndexing.py F_50000 sorted-inv-index
python evalIndexing.py F_50000 sorted-array
python evalIndexing.py F_50000 adapt-sorted
python evalIndexing.py F_50000 suffix-array
python evalIndexing.py F_50000 suffix-array-substr
python evalIndexing.py F_50000 robust-suffix-array
python evalIndexing.py F_50000 q-gram
python evalIndexing.py F_50000 canopy-th
python evalIndexing.py F_50000 canopy-nn
python evalIndexing.py F_50000 string-map-th
python evalIndexing.py F_50000 string-map-nn
