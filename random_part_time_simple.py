#!/usr/bin/python

##!/usr/bin/env sage -python # uncomment if using the Sage environment

import sys
import partitions as parts
import re
import time

class Timer:    
    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock()
        self.interval = self.end - self.start

def rand_part_Sage(N,S,sample_size):
    parts = []
    while len(parts) < sample_size:
        part = Partitions(N).random_element()
        if len(part) == S:
            parts.append(part)
    
    return parts
    
Q = 10000
N = 100
zeros = 'no'
sample_size = 1
print '\nQ =',Q,' N =',N,' Sample size =',sample_size

witches = ['multiplicity','top_down','divide_and_conquer','bottom_up']
print_statements = ['\nMultiplicity:','Top down:','Divide and conquer:','Bottom up:']

for i, x in enumerate(witches):
    with Timer() as t:
        which = witches[i]
        partitions = parts.rand_parts(Q,N,sample_size,which,zeros)
    t1 = t.interval   
    print print_statements[i],t1,'seconds'
