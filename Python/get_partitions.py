#! /usr/bin/env python

from __future__ import division
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


""" Code to generate random integer partitions for a total q having n parts. 
    combinations of q and n are obtained from data files formatted in columns.
    For example, with site as the first column, species as the second, and the
    species abundance as the third. This code uses the multiprocessing module.
    Randon integer partitions are saved to a file. """

def worker(NS_combo):
    """thread worker function"""
    #set_random_seed()
    random.seed()
    q = qn_combo[0]
    n = qn_combo[1]
    random_macros = parts.rand_partitions(q, n, 30, 'divide_and_conquer', zeros=True)
    # partitions including zeros = True; excluding zeros = False
    return random_macros


def get_rand_sample(qn_combo):
    
    unique_partitions = []
    pool = Pool()
    unique_partitions = pool.map(worker, [qn_combo, qn_combo, qn_combo, qn_combo, qn_combo, qn_combo, qn_combo, qn_combo])
    pool.close()
    pool.join()
    
    return unique_partitions
    
    
def get_partitions(dataset):

    DATA = open('/home/kenlocey/data1/' + dataset + '/' + dataset + '-data.txt','r')
    ct1 = 0
    ct2 = 0
    d = DATA.readline()
    label = re.match(r'\A\S*',d).group()
    val = int(re.findall(r'\d*\S$',d)[0])
    partition = [int(val)]
    partitions = []
        
    for d in DATA:
        ct1+=1
        label1 = re.match(r'\A\S*',d).group()
        if label1 == label:
            val = int(re.findall(r'\d*\S$',d)[0])
            if val > 0:
                partition.append(val)
                
        else:
            label = label1
            if len(partition) > 9 and sum(partition) <= 1000000:
                partition.sort()
                partition.reverse()
                partitions.append(partition)
                ct2+=1
            partition = []
            val = int(re.findall(r'\d*\S$',d)[0])
            if val > 0:partition.append(val)
    partitions.append(partition)
    DATA.close()
    return(partitions)
    
def get_qn_combos(datasets):
    qn_combos = []
    total_combos = 0
    for dataset in datasets:
        ct=0
        partitions = get_partitions(dataset)
        for partition in partitions:
            q = sum(partition)
            n = len(partition)
            qn_combos.append([q, n])
        #print dataset,     
        
    qn_combos = [list(x) for x in set(tuple(x) for x in qn_combos)]
    print len(qn_combos),'unique qn combos' 
    
    return (qn_combos)
    
        
def get_random_partitions_for_qn_combos(qn_combos):    
    #random.shuffle(qn_combos)
    qn_combos.sort(key=lambda x: float(x[0]))
    while qn_combos:
        ct = 0
        for qn_combo in qn_combos:
            
            ct+=1
            q = int(qn_combo[0])
            n = int(qn_combo[1])
            print q, n, len(qn_combos), 'NS combinations left'
            OUT = open('/home/kenlocey/partition_files/' + str(q) + '-' + str(n) + '.txt','a+')
            partitions = len(set(OUT.readlines()))
            
            if partitions < 500: 
                rand_partitions = get_rand_sample(qn_combo) # Use multiprocessing
                for i in rand_macros:
                    for p in i:
                        print>>OUT,p
                OUT.close()
                #NS_combos.remove(qn_combo)
                    
            elif len(qn_combos) == 1:
                qn_combos.remove(qn_combo)
                OUT.close
                break
            else:    
                qn_combos.remove(qn_combo)
                OUT.close()
            if not qn_combos:break
    return


datasets = ['BBS']
qncombos = get_qn_combos(datasets)
get_random_partitions_for_qn_combos(qncombos)
