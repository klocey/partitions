from __future__ import division

from sage.all import * # Sage is necessary for generating figure 2 of Locey and McGlinn (2013)
# This can be commented out if not generating figure 2 or if the user does not have Sage installed
# If Sage is installed, then from the command line: sage -python Locey_McGlinn_2013.py 

import sys
import os
import partitions as parts
sys.path.append("/home/kenlocey/modules/pymods")
import macroecotools
sys.path.append("/home/kenlocey/modules/FEASIBLE_FUNCTIONS")
import feasible_functions as ff
sys.path.append("/home/kenlocey/partitions/metrics")
import metrics as mt

from os import path, access, R_OK  # W_OK for write permission
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from scipy.stats import gaussian_kde
import  matplotlib.pyplot as plt
from pylab import *
import numpy as np
from scipy import stats
from scipy.stats import kstest
import random
from random import choice
import re
import math
import random, decimal


def test_first_lexical():
    
    partition = [123]
    q = 123
    n = 1
    k = 123
    first_partition = first_lexical(q, n, k)
    if partition is not first_partition:
        print 'first_lexical() is broken. Test 1'
    
    partition = [1,1,1,1,1,1,1,1,1,1]
    q = sum(partition)
    n = len(partition)
    k = max(partition)
    first_partition = first_lexical(q, n, k)
    if partition is not first_partition:
        print 'first_lexical() is broken. Test 2'
    
    partition = [10, 10, 1, 1, 1, 1]
    q = sum(partition)
    n = len(partition)
    k = max(partition)
    first_partition = first_lexical(q, n, k)
    if partition is not first_partition:
        print 'first_lexical() is broken. Test 1'
    
    return
    

def test_last_lexical():
    
    partition = [1]
    q = sum(partition)
    n = len(partition)
    last = last_lexical(q, n)
    if partition is not last:
        print 'last_lexical is broken. Test 1'
    
    partition = [1,1,1,1,1,1,1]
    q = sum(partition)
    n = len(partition)
    last = last_lexical(q, n)
    if partition is not last:
        print 'last_lexical is broken. Test 2'
    
    partition = [11,11,11,10,10,10,10]
    q = sum(partition)
    n = len(partition)
    last = last_lexical(q, n)
    if partition is not last:
        print 'last_lexical is broken. Test 3'
        
    return

def test_next_restricted_part():
    
    partition = [10,9,8,1]
    answer = [10,9,7,2]
    q = sum(partition)
    n = len(partition)
    next = next_restricted_part(partition)
    if next is not answer:
        print 'next_restricted_part is broken. Test 1'
        
    return

def test_min_max():
    
    q = 100
    n = 5
    minmax1 = 20
    minmax2 = min_max(q, n)
    if minmax1 != minmax2:
        print 'min_max() is broken. Test 1'
        
    q = 100
    n = 3
    minmax1 = 34
    minmax2 = min_max(q, n)
    if minmax1 != minmax2:
        print 'min_max() is broken. Test 2'   
   
    return


def test_P():
    
    D= {}
    q = 0
    k = 0
    answer = 1 # by convention
    numparts = parts.P(D, q, k)
    if numparts != answer:
        print 'P() is broken. Test 1'
    
    q = 0
    k = 123
    answer = 1 # by convention
    numparts = parts.P(D, q, k)
    if numparts != answer:
        print 'P() is broken. Test 2'
    
    q = 123
    k = 123
    answer = 2552338241
    numparts = parts.P(D, q, k)
    if numparts != answer:
        print 'P() is broken. Test 3'

    q = 123
    k = 0
    answer = 0
    numparts = parts.P(D, q, k)
    if numparts != answer:
        print 'P() is broken. Test 4'
        
    q = 12345
    k = 123
    answer = 488259162924433580696194373878466788895554319556195978121822180221785381227453675217103501271281020550492
    numparts = parts.P(D, q, k)
    if numparts != answer:
        print 'P() is broken. Test 5'
    
    return
    

def test_NrParts():
    
    q = 0
    numparts = parts.NrParts(q)
    if numparts != 1:
        print 'NrParts is broken. p(0) = 1 !=',numparts
    
    n = 1
    numparts = parts.NrParts(q, n)
    if numparts != 1:
        print 'NrParts is broken. p(0, 1) = 1 !=',numparts
    
    q = 12345
    answer = 69420357953926116819562977205209384460667673094671463620270321700806074195845953959951425306140971942519870679768681736
    numparts = parts.NrParts(q)
    if numparts != answer:
        print 'NrParts is broken. Answer for p(12345) is incorrect.'
    
    n = 876
    answer = 2513021291084594958506275164399019380305807385021575642918576266152785748154845293483022243930843412000601280303427
    numparts = parts.NrParts(q, n)
    if numparts != answer:
        print 'NrParts is broken. Answer for p(12345, 876) is incorrect.'

    return


def test_conjugate():
    
    partition = [1]
    conjugate = [1]
    if conj is not conjugate:
        print 'parts.conjugate() is broken. Test partition is [1]'
        
    partition = [10]
    conjugate = [1,1,1,1,1,1,1,1,1,1]
    if conj is not conjugate:
        print 'parts.conjugate() is broken. Test partition is [10]'
    
    partition = [3,2,1] # symmetrical partition and conjugate
    conjugate = [3,2,1]
    if conj is not conjugate:
        print 'parts.conjugate() is broken. Test partition is [3,2,1]'
    
    partition = [12,9,8,8,7,5,4,2,2,1]
    conjugate = [10, 9, 7, 7, 6, 5, 5, 4, 2, 1, 1, 1]

    conj = parts.conjugate(partition)
    if conj is not conjugate:
        print 'parts.conjugate() is broken'

    return
