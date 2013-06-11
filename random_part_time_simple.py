#!/usr/bin/env sage -python

import sys
sys.path.append("/home/kenlocey/modules/partitions")
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
    
N = 1000
S = 100
sample_size = 1

with Timer() as t:
    x = parts.rand_parts2(N,S,sample_size)
t1 = t.interval    
print 'N:',N,'S:',S,'sample size:',sample_size,'time:',t1,'seconds'
