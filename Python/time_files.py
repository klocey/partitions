#!/usr/bin/env python

import sys
import subprocess
import os
import numpy as np
from scipy import stats
import random
from random import choice
import re
import math
import random, decimal
import time
from multiprocessing import Pool, freeze_support

""" This file generates .txt files holding times taken for partitioning algorithms to generate random partitions for different combinations of Q and N.
    The time results held in these files were used to generate figure ? of Locey and McGlinn 2013. """

class Timer:    
    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock() 
        self.interval = self.end - self.start
        
        
def worker(QN_combo):
    set_random_seed()
    with Timer() as t:  
        subprocess.call(["nice", "sage", "-python", "get_partitions.py", str(QN_combo[0]), str(QN_combo[1]), str(QN_combo[2]), str(QN_combo[3]), str(QN_combo[4])])
    
    return [QN_combo[2], t.interval]
    

def get_rand_sample(NS_combo):
    
    pool = Pool()
    times = pool.map(worker, QN_combos)
    pool.close()
    pool.join()
    return times
    

algorithms = ['multiplicity','top_down','divide_and_conquer','bottom_up']
sample_size = 900 # takes approx 300 random macrostates to capture shape of the feasible set

Qs = range(20,1000,20)
zeros = [True, False]

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
                #print zero,alg,Q
                Ns = range(int(n1), int(nX)+step, int(step)) # list of N values
                
                QN_combos = []
                D = {}    
                for N in Ns: QN_combos.append([alg, Q, N, sample_size, zero])
                
                times = get_rand_sample(QN_combos) # list of times to be printed to a file
                
                if zero == True:
                    OUT = open('/home/kenlocey/partitions/time_files/timeFiles/Python_' + alg + '_zeros_Q=' + str(Q) + '.txt','a+')
                elif zero == False:
                    OUT = open('/home/kenlocey/partitions/time_files/timeFiles/Python_' + alg + '_Q=' + str(Q) + '.txt','a+')
                
                for t in times:
                    print>>OUT, t[0],t[1]
                OUT.close()
                    
                n1 = nX+step 
                if nX == Q: break
                elif nX + 8*step > Q: nX = Q
                else: nX += 8*step
