#!/usr/bin/env sage -python

from sage.all import *
import sys
import os
import numpy as np
from scipy import stats
import random
from random import choice
import re
import math
import random, decimal
import time
import subprocess
from multiprocessing import Pool, freeze_support

class Timer:    
    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock() 
        self.interval = self.end - self.start
        
        
def rand_part_Sage(Q,N,sample_size):
    parts = []
    while len(parts) < sample_size:
        part = Partitions(Q).random_element()
        if len(part) == N:
            parts.append(part)
    
    return parts

def worker(QN_combo):
    set_random_seed()
    with Timer() as t:
        partitions = subprocess.call(["Rscript", "/home/kenlocey/get_partitions.R", str(QN_combo[0]), str(QN_combo[1]), str(QN_combo[2]), str(QN_combo[3]), str(QN_combo[4]), str(QN_combo[5]), str(QN_combo[6])])
        
    return [QN_combo[2], t.interval]
    

def get_rand_sample_R(NS_combo):
    
    pool = Pool()
    times = pool.map(worker, QN_combos)
    pool.close()
    pool.join()
    return times
    

algorithms = ['multiplicity','top_down','divide_and_conquer','bottom_up']
sample_size = 300 # takes approx 300 random macrostates to capture shape of the feasible set

Qs = range(20,1000,20)
zeros = ['TRUE', 'FALSE']

for zero in zeros:
    for alg in algorithms:
    
        for Q in Qs:
            
            if Q <= 100: step = 1
            elif Q <= 300: step = 5
            elif Q <= 600: step = 10
            else: step = 20
            n1 = step
            nX = 8*step
            
            while nX <= Q: # 
                print zero,alg,Q
                Ns = range(int(n1), int(nX)+step, int(step)) # list of N values
                
                QN_combos = []
                D = {}    
                for N in Ns: QN_combos.append([alg, Q, N, sample_size, zero, "TRUE", "FALSE"])
                
                times = get_rand_sample_R(QN_combos) # list of times to be printed to a file
                
                if zero is 'TRUE':
                    OUT = open('/home/kenlocey/partitions/time_files/timeFiles/R_' + alg + '_zeros_Q=' + str(Q) + '.txt','a+')
                elif zero is 'FALSE':
                    OUT = open('/home/kenlocey/partitions/time_files/timeFiles/R_' + alg + '_Q=' + str(Q) + '.txt','a+')
                
                for t in times:
                    print>>OUT, t[0],t[1]
                OUT.close()
                    
                n1 = nX+step 
                if nX == Q: break
                elif nX + 8*step > Q: nX = Q
                else: nX += 8*step
                    













