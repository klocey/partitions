#!/usr/bin/env sage -python

from sage.all import *
import sys
sys.path.append("/home/kenlocey/modules/partitions")
import partitions as parts
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import os
import  matplotlib.pyplot as plt
from pylab import *
#import numpy as np
from scipy import stats
import random
from random import choice
import re
import math
import random, decimal
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

""" The code below compares the speed of the random partitioning function
    of Sage to the functions developed in Locey and McGlinn (20??), for
    cases when zero values are and are not allowed. The code also shows 
    a tradeoff in speed between the two methods (i.e. 'divide and conquer',
    'bottum-up') developed in L&M (20??). 
    
    The code generates figure 2 of Locey and McGlinn (20??) """


fig = plt.figure()

ii = 1
a = 0.9

Ns = [100,500,1000] # N values
sample_size = 100 # number of partitions to generate for each N-S combo

while ii <= 3: # for the first 3 subplots (i.e. 1st row)
    for N in Ns:
        print N
        ax = fig.add_subplot(2,3,ii)
    
        R1_times = [] # using 'bottom-up' method
        Ss = range(int(N/20.0),int(0.8*N),int(N/20.0)) # S values
        for S in Ss:    
            with Timer() as t:
                x = parts.rand_parts1(N,S,sample_size)
            for i in range(1):
                R1_times.append(round(t.interval,2))
        print 'R1',N,R1_times
        plt.plot(Ss,R1_times,lw=3,color='r',label='R1, N='+str(N),alpha=a)
 
        R2_times = [] # using 'divide-and-conquer' method
        Ss = range(int(N/20.0),int(0.8*N),int(N/20.0))
        for S in Ss:    
            with Timer() as t:
                x = parts.rand_parts2(N,S,sample_size)
            for i in range(1):
                R2_times.append(round(t.interval,2))
        print 'R2',N,R2_times
        plt.plot(Ss,R2_times,lw=3,color='b',label='R2, N='+str(N),alpha=a)
    
        R3_times = [] # using either 'divide-and-conquer' or 'bottom up'
        Ss = range(int(N/20.0),int(0.8*N),int(N/20.0))
        for S in Ss:    
            with Timer() as t:
                if N < 250 or S >= int(N/2.0): x = parts.rand_parts1(N,S,sample_size)
                else: x = parts.rand_parts2(N,S,sample_size)
            for i in range(1):
                R3_times.append(round(t.interval,2))
        print 'R3',N,R3_times      
        plt.plot(Ss,R3_times,lw=3,color='m',label='R3, N='+str(N),alpha=a)
    
        plt.tick_params(axis='both', which='major', labelsize=8)
        plt.xlabel("S",fontsize=8)
        if ii == 1:
            plt.ylim(0.0,0.1)
            plt.ylabel("Seconds",fontsize=8)
        #plt.setp(axins, xticks=Ns,yticks=[0.2,0.4,0.6,0.8,1.0])
        leg = plt.legend(loc=1,prop={'size':8})
        leg.draw_frame(False)
        ii+=1   


Ns = [50,100,200]
sample_size = 100
    
while ii <= 6: # for the last three subplots (i.e. 2nd row)
    for N in Ns:
        ax = fig.add_subplot(2,3,ii)
    
        R3_times = [] # using either 'divide-and-conquer' or 'bottom up'
        Ss = range(int(N/20.0),int(0.5*N),int(N/20.0))
        for S in Ss:    
            print 'KJL',S
            with Timer() as t:
                if N < 250 or S >= int(N/1.5): x = parts.rand_parts1(N,S,sample_size)
                else: x = parts.rand_parts2(N,S,sample_size)
            for i in range(1):
                R3_times.append(round(t.interval,2))
        print 'R3',N,R3_times

        Sage_times = []
        for S in Ss:    
            print 'Sage',S
            with Timer() as t:
                x = rand_part_Sage(N,S,sample_size)
            for i in range(1):
                Sage_times.append(round(t.interval,2))
        print 'Sage',N,Sage_times,'\n'

        Y = []
        for i, t in enumerate(R3_times):
            if Sage_times[i] > 0.0:
                Y.append(t/float(Sage_times[i]))
            else:
                Y.append(1.0)
        plt.plot(Ss,Y,lw=3,color='m',label='R3, N='+str(N),alpha=a)
    
        if ii == 4:
            plt.ylabel("time(algorithm)/time(sage)",fontsize=8)
        
        ii+=1
        plt.tick_params(axis='both', which='major', labelsize=8)
        plt.xlabel("S",fontsize=8)
        plt.ylim(0.0,0.2)
        #plt.setp(axins, xticks=Ns,yticks=[0.2,0.4,0.6,0.8,1.0])
        leg = plt.legend(loc=1,prop={'size':8})
        leg.draw_frame(False)
        
plt.savefig('/home/kenlocey/Fig2-'+str(sample_size)+'.png', dpi=500, pad_inches=0)
