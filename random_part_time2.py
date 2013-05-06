#!/usr/bin/env sage -python

# when calling from the command line enter: sage -python random_part_time2.py

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
##################################################################################
##################### Compare functions of Sage and K. Locey #########

fig = plt.figure()

Ns = [100,150,200]  
sample_sizes = [10,50]
ii = 1
a = 0.9

for sample_size in sample_sizes:
    """ code in this loop generates figures that compare Sage to rand_parts1
    and rand_parts2.
    
    It uses repeated chunks of code, so clearly, I wasn't going for
    efficiency or elegance when I was writing this -KL  """
    
    for N in Ns:
    
        ax = fig.add_subplot(3,3,ii)
    
        R1_times = [] # times for rand_parts1
        Ss = range(int(N/20.0),int(0.5*N),int(N/20.0))
        for S in Ss:    
            with Timer() as t:
                x = parts.rand_parts1(N,S,sample_size)  # get me some random partitions
            for i in range(1):
                R1_times.append(round(t.interval,2))
        print 'R1',N,R1_times
 
        R2_times = [] # times for rand_parts2
        Ss = range(int(N/20.0),int(0.5*N),int(N/20.0))
        for S in Ss:    
            with Timer() as t:
                x = parts.rand_parts2(N,S,sample_size) # get me some random partitions
            for i in range(1):
                R2_times.append(round(t.interval,2))
        print 'R2',N,R2_times


        R3_times = []
        Ss = range(int(N/20.0),int(0.5*N),int(N/20.0))
        for S in Ss:    
            with Timer() as t:
                if N < 250 or S >= int(N/1.5): x = parts.rand_parts1(N,S,sample_size)
                else: x = parts.rand_parts2(N,S,sample_size) # get me some random partitions ...depending on N and N/S
            for i in range(1):
                R3_times.append(round(t.interval,2))
        print 'R3',N,R3_times

    
        Sage_times = []
        for S in Ss:    
            with Timer() as t:
                x = rand_part_Sage(N,S,sample_size)
            for i in range(1):
                Sage_times.append(round(t.interval,2))
        print 'Sage',N,Sage_times,'\n'


        Y = []
        for i, t in enumerate(R1_times):
            if Sage_times[i] > 0.0:
                Y.append(t/float(Sage_times[i]))
            else:
                Y.append(1.0)
       
        plt.plot(Ss,Y,lw=3,color='r',label='R1, N='+str(N),alpha=a)
      
        Y = []
        for i, t in enumerate(R2_times):
            if Sage_times[i] > 0.0:
                Y.append(t/float(Sage_times[i]))
            else:
                Y.append(1.0)
       
        plt.plot(Ss,Y,lw=3,color='b',label='R2, N='+str(N),alpha=a)
    
        Y = []
        for i, t in enumerate(R3_times):
            if Sage_times[i] > 0.0:
                Y.append(t/float(Sage_times[i]))
            else:
                Y.append(1.0)
       
        plt.plot(Ss,Y,lw=3,color='m',label='R3, N='+str(N),alpha=a)
    
        if ii == 1 or ii == 4 or ii == 7:
            plt.ylabel("time(algorithm)/time(sage)",fontsize=8)
        
        ii+=1
        
        plt.tick_params(axis='both', which='major', labelsize=8)
        plt.xlabel("S",fontsize=8)
        plt.ylim(0.0,0.2)
        #plt.setp(axins, xticks=Ns,yticks=[0.2,0.4,0.6,0.8,1.0])
        leg = plt.legend(loc=1,prop={'size':8})
        leg.draw_frame(False)
        

Ns = [100,200,300]
sample_size = 20

for N in Ns:
    """ code in this loop compares the speed of different algorithms 
    derived by K. Locey across ratios of the total N to the number 
    of parts S, for which Sage would be impractical.
    
    It uses repeated chunks of code, so clearly, I wasn't going for
    efficiency or elegance when I was writing this -KL
     """

    ax = fig.add_subplot(3,3,ii)
    
    R1_times = []
    Ss = range(int(N/20.0),int(0.8*N),int(N/20.0))
    for S in Ss:    
        with Timer() as t:
            x = parts.rand_parts1(N,S,sample_size) # get me some random partitions
        for i in range(1):
            R1_times.append(round(t.interval,2))
    print 'R1',N,R1_times
    plt.plot(Ss,R1_times,lw=3,color='r',label='R1, N='+str(N),alpha=a)
 
    R2_times = []
    Ss = range(int(N/20.0),int(0.8*N),int(N/20.0))
    for S in Ss:    
        with Timer() as t:
            x = parts.rand_parts2(N,S,sample_size) # get me some random partitions
        for i in range(1):
            R2_times.append(round(t.interval,2))
    print 'R2',N,R2_times
    plt.plot(Ss,R2_times,lw=3,color='b',label='R2, N='+str(N),alpha=a)
    
    R3_times = []
    Ss = range(int(N/20.0),int(0.8*N),int(N/20.0))
    for S in Ss:    
        with Timer() as t:
            if N < 250 or S >= int(N/1.5): x = parts.rand_parts1(N,S,sample_size)
            else: x = parts.rand_parts2(N,S,sample_size) # get me some random partitions ...depending on N and N/S
        for i in range(1):
            R3_times.append(round(t.interval,2))
    print 'R3',N,R3_times      
    plt.plot(Ss,R3_times,lw=3,color='m',label='R3, N='+str(N),alpha=a)
    
    plt.tick_params(axis='both', which='major', labelsize=8)
    plt.xlabel("S",fontsize=8)
    if ii == 7:
        #plt.ylim(0.0,0.015)
        plt.ylabel("time(algorithm)/time(sage)",fontsize=8)
    #plt.setp(axins, xticks=Ns,yticks=[0.2,0.4,0.6,0.8,1.0])
    leg = plt.legend(loc=1,prop={'size':8})
    leg.draw_frame(False)
    
    ii+=1   
plt.savefig('/home/kenlocey/time3-'+str(sample_size)+'.png', dpi=400, pad_inches=0)

