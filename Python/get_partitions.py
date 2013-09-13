#!/usr/bin/env sage -python

#from sage.all import *
import sys
import subprocess
import partitions as parts
import os
import numpy as np
from scipy import stats
import random
from random import choice
from os import path, access, R_OK  # W_OK for write permission
import re
import math
import random, decimal
import time
from multiprocessing import Pool, freeze_support


def worker(NS_combo):
    """thread worker function"""
    #set_random_seed()
    random.seed()
    N = NS_combo[0]
    S = NS_combo[1]
    random_macros = parts.rand_parts(N, S, 10, 'divide_and_conquer', zeros=True)
    return random_macros


def get_rand_sample(NS_combo): #choose between worker2 (random partitioning alg derived by KJL (may not return uniform random samples) and worker1 (random partitioning alg provided by SAGE)
    
    unique_SADs = []
    pool = Pool()
    unique_SADs = pool.map(worker, [NS_combo,NS_combo,NS_combo,NS_combo,NS_combo,NS_combo,NS_combo,NS_combo])
    """ worker1 and worker2 call different functions for generating random macrostates. worker1 uses Sage's function (def. worker2 uses the function developed by Ken Locey (faster)."""
    pool.close()
    pool.join()
    return unique_SADs
    
    
def get_SADs(dataset):

    DATA = open('/home/kenlocey/data1/' + dataset + '/' + dataset + '-data.txt','r')
    ct1 = 0
    ct2 = 0
    d = DATA.readline()
    m0 = re.match(r'\A\S*',d).group()
    m2 = int(re.findall(r'\d*\S$',d)[0])
    SAD = [int(m2)]
    SADs = []
        
    for d in DATA:
        ct1+=1
        m1 = re.match(r'\A\S*',d).group()
        if m1 == m0:
            m2 = int(re.findall(r'\d*\S$',d)[0])
            if m2 > 0:
                SAD.append(m2)
                
        else:
            site_name = m0
            m0 = m1
            if len(SAD) > 9 and sum(SAD) <= 1000000:
                SAD.sort()
                SAD.reverse()
                SADs.append(SAD)
                ct2+=1
            SAD = []
            abundance = int(re.findall(r'\d*\S$',d)[0])
            if abundance > 0:SAD.append(abundance)
    DATA.close()
    return(SADs)
    
def get_NS_combos(datasets):
    NS_combos = []
    total_combos = 0
    for dataset in datasets:
        ct=0
        SADs = get_SADs(dataset)
        for SAD in SADs:
            N = sum(SAD)
            S = len(SAD)
            NS_combos.append([N,S])
        print dataset,     
        
    NS_combos = [list(x) for x in set(tuple(x) for x in NS_combos)]
    print len(NS_combos),'unique NS_combos' 
    
    return (NS_combos)
    
        
def get_random_macrostates_for_NScombos(NS_combos):    
    #random.shuffle(NS_combos)
    while NS_combos:
        ct = 0
        for NS_combo in NS_combos:
            
            ct+=1
            N = int(NS_combo[0])
            S = int(NS_combo[1])
            print N, S, len(NS_combos), 'NS combinations left'
            OUT = open('/home/kenlocey/combined1/' + str(N) + '-' + str(S) + '.txt','a+')
            macros = len(set(OUT.readlines()))
            
            if macros < 1000: 
                rand_macros = get_rand_sample(NS_combo) # Use multiprocessing
                for i in rand_macros:
                    for p in i:
                        print>>OUT,p
                OUT.close()
                #NS_combos.remove(NS_combo)
                    
            elif len(NS_combos) == 1:
                NS_combos.remove(NS_combo)
                OUT.close
                break
            else:    
                NS_combos.remove(NS_combo)
                OUT.close()
            if not NS_combos:break
    return


def get_SSADs(dataset):
    
    richness = 301
    abu_class = []
    SSADs = [[] for x in xrange(richness)] # one frequency distribution per species
    
    PATH = '/home/kenlocey/data1/' + dataset + '/' + dataset + '-SSADs.txt'
    if path.exists(PATH) and path.isfile(PATH) and access(PATH, R_OK):
        DATA = open(PATH,'r')
        d = DATA.readline() # first line has no abundance info, so skip it
        
        for d in DATA:
            line = list(str(re.match(r'\A\S*',d).group()).split(',')) # site
            line.pop(0) # first item is the spatial grain and isn't needed, so pop it
            abu_class.append(int(line[0])) # now, the first item is the abundance class, assign it to a variable...
            line.pop(0) # then pop it from the list. Now, only frequencies for species remain in the list
            
            for i, item in enumerate(line):
                SSADs[i].append(int(item))
        
        
        return SSADs 


def hist_to_rank(SSADs):
    
    richness = len(SSADs)
    abu_class = range(0,len(SSADs[0])+1)
    
    rank_vectors = [[] for x in xrange(richness)] # one rank distribution per species
    for i, ssad in enumerate(SSADs):
        for j, ab in enumerate(ssad):
            if ab == 0:
                continue
            else:
                #print abu_class[j], ab
                _list = [abu_class[j]]*ab
                rank_vectors[i].extend(_list)
        
    species_abs = []
    RADs = []
    
    for i, _list in enumerate(rank_vectors):
        if sum(_list) >= 0:
            _list.reverse() 
            species_abs.append(sum(_list))
            RADs.append(_list)
    # uncommenting these lines reveals that the abundances for species on BCI
    # matches the BCI RAD, meaning nothing has gone wrong.
    #species_abs.sort()
    #species_abs.reverse()
    #for i in species_abs: print i
    #print species_abs
    #print len(rank_vectors), max(species_abs)
        
    # So, at this point we have a rank distribution for each species.        
    
    return RADs


def NScombos_SSAD(datasets):
    
    NScombos = []
    
    SSADs = get_SSADs(datasets[0])
    RADs = hist_to_rank(SSADs)
    
    for rad in RADs:
        N = sum(rad)
        S = len(rad)
        if N >= 2000 and N < 10000: 
            NScombos.append([N,S])
        
    return NScombos
    


datasets = ['BCI']

NScombos = NScombos_SSAD(datasets)
print len(NScombos)
#NScombos = get_NS_combos(datasets)

get_random_macrostates_for_NScombos(NScombos)







#def worker(QN_combo):
#    set_random_seed()
#        subprocess.call(["nice", "sage", "-python", "get_partitions.py", str(QN_combo[0]), str(QN_combo[1]), str(QN_combo[2]), str(QN_combo[3]), str(QN_combo[4])])
#    return [QN_combo[2], t.interval]
