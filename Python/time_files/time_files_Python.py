#!/usr/bin/env sage -python

from sage.all import *
import sys
sys.path.append("/home/kenlocey/partitions/Python_script")
import partitions as parts
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
from os import path, access, R_OK  # W_OK for write permission


class Timer:    
    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args):
        self.end = time.time() 
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
        partitions = parts.rand_parts(QN_combo[0], QN_combo[1], QN_combo[2], QN_combo[3], QN_combo[4], QN_combo[5])
        
    return [QN_combo[1], t.interval]
    

def get_rand_sample_Python(NS_combo):
    
    pool = Pool()
    times = pool.map(worker, QN_combos)
    pool.close()
    pool.join()
    return times
    

algorithms = ['multiplicity','top_down','divide_and_conquer','bottom_up']
sample_size = 1 # takes approx 300 random macrostates to capture shape of the feasible set

#Qs1 = range(20,100,1)
#Qs = range(100,1000+10,10)
#Qs = Qs1+Qs2
Qs = [5000,20000,200000,500000] 
Nss = [[10,20,40,80], [20,40,80,160], [40,80,160,320], [80,160,320,640]]

zeros = [False, True]

for zero in zeros:
    for alg in algorithms:
    
        for i, Q in enumerate(Qs):
            
            if Q <= 100: step = 2
            elif Q <= 300: step = 5
            elif Q <= 600: step = 10
            else: step = 20
            n1 = step
            nX = 8*step
            
            while nX <= Q: #
            
                if zero == True: PATH = '/home/kenlocey/partitions/time_files/timeFiles/Python_' + alg + '_zeros_Q=' + str(Q) + 'n=1.txt'
                elif zero == False: PATH = '/home/kenlocey/partitions/time_files/timeFiles/Python_' + alg + '_Q=' + str(Q) + 'n=1.txt'
                
                if path.exists(PATH) and path.isfile(PATH) and access(PATH, R_OK): pass
                
                else:
                    print zero,alg,Q
                    #Ns = range(int(n1), int(nX)+step, int(step)) # list of N values
                    Ns = Nss[i]
                    
                    QN_combos = []
                    D = {}    
                    for N in Ns: QN_combos.append([Q, N, sample_size, alg, D, zero])
                    
                    times = get_rand_sample_Python(QN_combos) # list of times to be printed to a file
                    
                    OUT = open(PATH,'a+')
                    for t in times:
                        print>>OUT, t[0],t[1]
                    
                    OUT.close()
                    
                n1 = nX+step 
                if nX == Q: break
                elif nX + 8*step > Q: nX = Q
                else: nX += 8*step
                    













