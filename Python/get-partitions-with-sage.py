from __future__ import division
from sage.all import *
import sys
import os
from os import path, access, R_OK  # W_OK for write permission
import random
from random import choice
import re


def rand_parts_Sage(q, n, sample_size, zeros):
    """" A Sage function for generating random partitons of an integer q
    having exactly n parts (i.e. Partitions(q).random_element())."""
    
    if zeros == True:
        q += n
    
    partitions = []
    while len(partitions) < sample_size:
        random.seed()
        partition = Partitions(q).random_element() # may have to be commented out
        if len(partition) == n:
            partitions.append(partition)
          
        else:
            partition = Partition(partition).conjugate()
            if len(partition) == n and sum(partition) == q:
                partitions.append(partition)
            
    if zeros == False:
        return partitions

    elif zeros == True:
        weakpartitions = []
        for partition in partitions:
            weakpartition = []
            for part in partition:
                weakpartition.append(part-1) 
            if sum(weakpartition) == q-n and len(weakpartition) == n:
                weakpartitions.append(weakpartition)
            else:
                print 'weakpartitions is broke'    
        return weakpartitions


def random_partitions_with_Sage(): 
    
    qn_combos = [[50,10],[100,20],[200,40]]
    sample_size = 300 
    
    for combo in qn_combos:
        q = combo[0]
        n = combo[1]
        
        OUT = open('/home/kenlocey/partitions/new/testfiles/sage_zeros_q=' + str(q) + '_n='+str(n)+'.txt','a+')
        partitions = rand_parts_Sage(q, n, sample_size, True)
        
        for partition in partitions:
            print>>OUT, partition
        OUT.close()
        
        OUT = open('/home/kenlocey/partitions/new/testfiles/sage_q=' + str(q) + '_n='+str(n)+'.txt','a+')
        partitions = rand_parts_Sage(q, n, sample_size, False)
    
        for partition in partitions:
            print>>OUT, partition
        OUT.close() 


      
random_partitions_with_Sage()

















