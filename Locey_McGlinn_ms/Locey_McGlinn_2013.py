from __future__ import division

from sage.all import * # Sage is necessary for generating figure 2 and 3 of Locey and McGlinn (2013)
# This can be commented out if not generating figure 2 or 3 or if the user does not have Sage installed
# If Sage is installed, then from the command line: sage -python Locey_McGlinn_2013.py 

import sys
import os
from os import path, access, R_OK  # W_OK for write permission
import re
import pypartitions as parts

sys.path.append("metrics")
import metrics as mt

from mpl_toolkits.axes_grid.inset_locator import inset_axes
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
from matplotlib.pylab import *
import numpy as np
from scipy import stats
import random
from random import choice
import math
import decimal
import time


def get_kdens(summands):
    """ Finds the kernel density function across a sample of parts
    of partitions for a given total (N) and number of parts (S) """
    
    density = gaussian_kde(summands)
    n = len(summands)
    xs = np.linspace(float(min(summands)),float(max(summands)),n)
    density.covariance_factor = lambda : .5
    density._compute_covariance()
    D = [xs,density(xs)]
    return D
    
    
class Timer:    
    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock()
        self.interval = self.end - self.start
        
        
def rand_part_Sage(q, n, sample_size):
    """" A Sage function for generating random partitons of an integer q
    having exactly n parts (i.e. Partitions(q).random_element())."""
    
    partitions = []
    while len(partitions) < sample_size:
        partition = Partitions(q).random_element() # may have to be commented out
        # if Sage is not installed and imported
        if len(partition) == n:
            partitions.append(partition)
    
    return partitions


def obs_pred_rsquare(obs, pred):
    """Determines the prop of variability in a data set accounted for by a model.
    In other words, this determines the proportion of variation explained by
    the 1:1 line in an observed-predicted plot. """
    return 1 - sum((obs - pred) ** 2) / sum((obs - np.mean(obs)) ** 2)


def time_trials_sage():
    """ Compare the speed of the random partitioning function of Sage to the
    functions developed in Locey and McGlinn (2013), for cases when zero
    values are and are not allowed."""

    fig = plt.figure()
    algorithms = ['multiplicity', 'top_down', 'divide_and_conquer', 'bottom_up']
    LABELS = ['Multiplicity:', 'Top down:', 'Divide and conquer:', 'Bottom up:']
    COLORS = ['#00CED1', '#FF1493', 'k', 'gray']
    qs = [50, 100, 150, 200]
    sample_size = 300
    
    for i, q in enumerate(qs):
    
        ax = fig.add_subplot(2,2,i+1)
        print '\n',q
        
        step = int(q/10)
        ns = range(step, int(0.5*q)+step,step)
        
        Sage_times = []
        for n in ns:    
            print 'Sage', n
            
            with Timer() as t:
                x = rand_part_Sage(q, n, sample_size)
            Sage_times.append(round(t.interval,2))
        
        print 'Sage', q, Sage_times, '\n'
        
        for ii, name in enumerate(algorithms):
            times = []
            for n in ns:
                zeros=False
                D = {}    
                with Timer() as t:
                    x = parts.rand_partitions(q, n, sample_size, name, D, zeros)
                times.append(round(t.interval,2))
            
            Y = []
            for iii, t in enumerate(times):
                if Sage_times[ii] > 0.0:
                    Y.append(float(Sage_times[iii])/t)
                else:
                    Y.append(1.0)
            
            plt.plot(ns, Y, lw=3, color=COLORS[ii], label=LABELS[ii])
            plt.yscale('log')
            
        if i == 0:
            plt.text(1.5,11,'Sage/algorithm',fontsize=10,rotation='90')
            legend(bbox_to_anchor=(-0.03, 1.1, 2.25, .2), loc=10, ncol=4, mode="expand",prop={'size':9})#, borderaxespad=0.)
        
        plt.tick_params(axis='both', which='major', labelsize=8)
        plt.xlabel('number of parts',fontsize=10)
        
    plt.savefig('time_trials_sage.png', dpi=500) #, pad_inches=0)
    return


def time_trials_bigQ(  ):
    """ Compare the speed of the random partitioning functions developed
    in Locey and McGlinn (2013), for cases when zero values are and
    are not allowed. The code generates figure ? of the Appendix from
    Locey and McGlinn (2013) """
    
    fig = plt.figure()
    
    sample_size=1
    zeros = False
    qs = [5000, 20000, 200000, 500000] # at these values, this script will take >1 day to run
    
    for i, q in enumerate(qs):
        ax = fig.add_subplot(2,2,i+1)
        if i==0:
            algorithms = ['multiplicity','top_down','divide_and_conquer','bottom_up']
            LABELS = ['Multiplicity:','Top down:','Divide and conquer:','Bottom up:']
            COLORS = ['#00CED1','#FF1493','k','gray']
            ns = [10,20,40,80]
           
        if i==1:
            algorithms = ['multiplicity','top_down','divide_and_conquer']
            LABELS = ['Multiplicity:','Top down:','Divide and conquer:']
            COLORS = ['#00CED1','#FF1493','k']
            ns = [20,40,80,160]
        
        elif i==2:
            algorithms = ['multiplicity','top_down']
            LABELS = ['Multiplicity:','Top down:']
            COLORS = ['#00CED1','#FF1493']
            ns = [40,80,160,320]
            
        elif i==3:
            algorithms = ['multiplicity']
            LABELS = ['Multiplicity:']
            COLORS = ['#00CED1']
            ns = [80,160,320,640]
            
        for ii, name in enumerate(algorithms):
            times = []
            for n in ns:
                D = {}
                with Timer() as t:
                    x = parts.rand_partitions(q, n, sample_size, name, D, zeros=False)
                times.append(round(t.interval,2))
            
            plt.plot(ns,times,color=COLORS[ii],lw=3,label=LABELS[ii])
            
        plt.xlabel('number of parts, sample size='+str(sample_size)+')',fontsize=8)    
        plt.tick_params(axis='both', which='major', labelsize=8)
        
        if i==0 or i==2:
                plt.ylabel("Seconds",fontsize=10)
                if i==0:
                    legend(bbox_to_anchor=(0.0, 1.1, 2.2, .2), loc=10, ncol=4, mode="expand",prop={'size':9})#, borderaxespad=0.)
    
    plt.savefig('time_trials_bigQ.png', dpi=500)#, pad_inches=0)
    return


def winner():
    """ Reveal how the random partitioning functions developed in Locey and McGlinn (2013),
    compare across many combinations of Q and N for cases when zero values are and are not
    allowed. The code generates figure 3 of the Appendix from Locey and McGlinn (2013) """
    
    algorithms = [['multiplicity','#00CED1'], ['top_down','#FF1493'], ['divide_and_conquer','k'], ['bottom_up','gray']]
    qs1 = range(20,100,5)
    qs2 = range(100,1000+20,20)
    qs = qs1+qs2
    zeros = [True, False]

    BigList = [[],[]]
    for zero in zeros:
    
        for q in qs:
        
            bestTime = 10**10 # starting with a ridiculously high time that is easy to beat
            bestAlg = None 
            
                        
            for alg in algorithms:
                if zero == True:
                    DATA = open('timeFiles/Python_' + alg[0] + '_zeros_Q=' + str(q) + '.txt','r')
                    for d in DATA:
                        n = int(re.findall(r'\d*\S\d*',d)[0])
                        t = float(re.findall(r'\d*\S\d*',d)[1])
                        BigList[0].append([alg[1], q, n, t]) # assigns a color value representing the algorithm
            
                elif zero == False:
                    DATA = open('timeFiles/Python_' + alg[0] + '_Q=' + str(q) + '.txt','r')
                    for d in DATA:
                        n = int(re.findall(r'\d*\S\d*',d)[0])
                        t = float(re.findall(r'\d*\S\d*',d)[1])
                        BigList[1].append([alg[1], q, n, t]) # assigns a color value representing the algorithm
            
                        
                DATA.close()

    fig = plt.figure()
    x = 1
    for _list in BigList:
        print 'plot:',x
        ax = fig.add_subplot(2,2,x)
        D = {}
        for i in _list:
            colorVal = i[0]
            q = i[1]
            n = i[2]
            t = i[3]

            if (q, n) not in D:
                D[(q, n)] = [colorVal, t]
         
            elif (q, n) in D:
                if t < D[(q, n)][1]:
                    D[(q, n)] = [colorVal, t]
                elif t == D[(q, n)][1]:
                    D[(q, n)] = ['y', t]
    
        for key, value in D.items():
        
            if len(value) > 2: # this indicates there's been a tie between algorithms
                print 'There has been a tie between algorithms. Good-bye.'
                sys.exit()
            #    value.sort(key=lambda x: float(x[1])) # sort list of lists by times
            #    value.reverse()
            #    value = value[0]
            
            qn_combo = list(key)
            q = int(qn_combo[0])
            n = int(qn_combo[1])
            c = value[0]
            t = value[1]
        
            sz=0
            if q <= 100: sz = 0.75
            elif q <= 300: sz = 1
            elif q <= 600: sz = 2
            else: sz = 3
            plt.scatter(q, n, marker='s', s=sz, color=c, edgecolor=c,alpha=0.9)
    
        plt.axis([20,1000,2,1000])
        plt.tick_params(axis='both', which='major', labelsize=8)
    
        if x == 1:
            plt.text(-170,800,'Number of parts, n',fontsize=10,rotation='90')
            plt.text(1000,-170,'Total, q',fontsize=10)
            for alg in algorithms:
                plt.plot(-1, -1, lw=3, color=alg[1], label=alg[0])
            plt.legend(bbox_to_anchor=(-0.02, 1.00, 2.24, .2), loc=10, ncol=4, mode="expand",prop={'size':9})#, borderaxespad=0.)
        
        # insert inset of low q values
        if x == 1: a = plt.axes([0.16, 0.75, .13, .13], axisbg='w')
        elif x == 2: a = plt.axes([0.58,0.75, .13, .13], axisbg='w')
        for key, value in D.items():
            qn_combo = list(key)
            q = int(qn_combo[0])
            if q <= 100:
                n = int(qn_combo[1])
                c = value[0]
                t = value[1]
                plt.scatter(q, n, marker='s', s=3, color=c, edgecolor=c,alpha=0.9)
    
        plt.setp(a, xlim=(20,101), ylim=(2,101), xticks=[20,40,60,80,100], yticks=[20,40,60,80,100])
        plt.tick_params(axis='both', which='major', labelsize=7)
    
        x+=1
       
    plt.savefig('algs_winner.png', dpi=500, bbox_inches = 'tight', pad_inches=0.4)
    plt.close()


def kdens_unbias(): 
    """ The code below compares random partitioning nplottions of Sage and Locey and McGlinn (2013)
    to full feasible sets. These analyses confirm that the algorithms are unbiased. The code
    uses full feasible sets, the random partitioning function in Sage, and the random partitioning
    for cases when 0' are or are not allowed."""

    algs = ['multiplicity','top_down','divide_and_conquer','bottom_up']
    colors = ['#00CED1','#FF1493','k','gray']

    fig = plt.figure()
    nplot = 1 # a variable used to designate subplots
    sample_size = 10000 # min number of macrostates needed to safely capture distributional
                      # features across the feasible set
    
    metrics = ['gini', 'variance', 'median', 'skewness', 'evar']
    metric = metrics[2]
        
    while nplot <= 4:
        ax =fig.add_subplot(2,2,nplot)

        if nplot < 3:
            q = 50 # values of q and n small enough to generate 
            n = 10 # the entire feasible set
        else:
            q = 100 # values of q and n requiring random samples
            n = 20  # of feasible sets
        
        partitions = []
        for i, alg in enumerate(algs):
            if nplot == 1 or nplot == 3:
                zeros = False
                D = {}
                partitions = parts.rand_partitions(q, n, sample_size, alg, D, zeros)
            else:
                D = {}
                zeros = True
                partitions = parts.rand_partitions(q, n, sample_size, alg, D, zeros)
                
            kdens = mt.get_kdens_obs(partitions, metric)
            plt.xlim(min(kdens[0]), max(kdens[0]))
            plt.plot(kdens[0], kdens[1], color=colors[i], lw=0.7)
            
        if nplot == 1: # using the full feasible set, no zero values (i.e. proper integer partitions)
            
            partitions = []
            numparts = parts.NrParts(q, n)    
            partition = parts.first_lexical(q, n, None)
            partitions.append(partition)
            ct2 = 0
            while len(partitions) < numparts:
                    
                partition = parts.next_restricted_part(partition)
                if len(partition) == n: partitions.append(partition) 
                else:
                    print 'bug in next_restricted_part()'
                    sys.exit()    
                
            kdens = mt.get_kdens_obs(partitions, metric)
            plt.xlim(min(kdens[0]), max(kdens[0]))
            plt.plot(kdens[0], kdens[1], color='r', lw=3.0, alpha=0.5)
                
        elif nplot == 2: # using the full feasible set, zero values included
            partitions = []    
                
            for p in Partitions(q):
                partition = list(p)
                
                if len(partition) == n:
                    partitions.append(partition) 
                    
                elif len(partition) < n:
                    zeros = [0]*(n-len(partition))
                    partition.extend(zeros)
                    partitions.append(partition)
            
            kdens = mt.get_kdens_obs(partitions, metric)
            plt.xlim(min(kdens[0]), max(kdens[0]))
            plt.plot(kdens[0], kdens[1], color='r', lw=3.0, alpha=0.5)
    
        elif nplot == 3: 
            partitions = []
            while len(partitions) < sample_size: # Use the random partition nplottion in Sage to generate partitions for q and n
                partition = Partitions(q).random_element()
                if len(partition) == n:
                    partitions.append(partition)
                     
                else:
                    partition = parts.conjugate(partition)
                    if len(partition) == n:
                        partitions.append(partition)
                             
            kdens = mt.get_kdens_obs(partitions, metric)
            plt.xlim(min(kdens[0]), max(kdens[0]))
            plt.plot(kdens[0], kdens[1], color='r', lw=3.0, alpha=0.5)
        
        elif nplot == 4:
            partitions = []
            while len(partitions) < sample_size: # Use the random partition nplottion in Sage to generate partitions for q and n
                part = list(Partitions(q).random_element())
                if len(part) == n:
                    partitions.append(part)
                
                elif len(part) < n:
                    zeros = [0]*(n - len(part))
                    part.extend(zeros)
                    partitions.append(part)
                
            kdens = mt.get_kdens_obs(partitions, metric)
            plt.xlim(min(kdens[0]), max(kdens[0]))
            plt.plot(kdens[0], kdens[1], color='r', lw=3.0, alpha=0.5)
                     
        if nplot == 1:
            plt.plot([0],[0], color='#00CED1', lw=2, label = 'Multiplicity')
            plt.plot([0],[0], color='#FF1493',lw=2, label='Top-down')    
            plt.plot([0],[0], color='k',lw=2, label='Divide & Conquer')
            plt.plot([0],[0], color='gray',lw=2, label='Bottom-up')
            plt.plot([0],[0], color='r',lw=2, label='FS q='+str(q)+', n='+str(n),alpha=0.5)
            plt.legend(bbox_to_anchor=(-0.02, 1.00, 2.24, .2), loc=10, ncol=5, mode="expand",prop={'size':8})#, borderaxespad=0.)
            
        if nplot == 1 or nplot == 3:
            plt.ylabel("density", fontsize=12)    
        
        if nplot == 3 or nplot == 4:
            plt.xlabel(metric, fontsize=12)
        
        print nplot
        nplot+=1
        
        plt.tick_params(axis='both', which='major', labelsize=8)
        
    plt.savefig('kdens_'+metric+'_'+str(sample_size)+'.png', dpi=500, pad_inches=0)


def get_skews(_list):

    skews = []
    for i in _list:
        skews.append(stats.skew(i))
    
    return skews
    

def vectorTohist(vector, zeros):
    
    bins = range(0,40)#[0,1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768]
    n = len(vector)
    SSAD = [0]*len(bins)
    for p in vector:
        if p < 40: SSAD[p] += 1/n
        
    while SSAD[-1] == 0:
        SSAD.pop()   
    
    return SSAD    


def get_all_SSADs(qnlist): # Figure X Locey and McGlinn (2013)      
    
    s_size = 300 
    i = 1
    fig = plt.figure()
    q = qnlist[0]
    nlist = qnlist[1]
    
    colors = ['#00CED1', '#FF1493', 'gray']
    
    for ind, num in enumerate(nlist):
        ax = fig.add_subplot(3,3,i)
        
        n = num  
        print i
        ct = 0
        clr = colors[ind]
        
        sample_size = s_size
        zeros = True    
        D = {}
        name = 'divide_and_conquer'
        SSADs = parts.rand_partitions(q, n, sample_size, name, D, zeros)
        
        for SSAD in SSADs:
            Y = vectorTohist(SSAD, zeros)
            
            while len(Y) > 40:
                Y.pop()   
            while len(Y) < 40:
                Y.append(0)
                
            x = range(0,40)
            plt.bar(x,Y, color=clr, linewidth=0, align='center', alpha = 0.015)
                
        plt.bar([0],[0], color=clr, linewidth=0, align='center', label= 'q='+str(q)+', n='+str(n))      
        ct+=1           
        plt.xlim(-1,8)
        plt.xlabel("Abundance class",fontsize=8)
        plt.ylabel("Frequency",fontsize=8)
        plt.tick_params(axis='both', which='major', labelsize=5)
        plt.setp(ax, xticks=[0,1,2,3,4,5,6,7,8])
        leg = plt.legend(loc=1,prop={'size':7})
        leg.draw_frame(False)
          
        i+=1
        
        ct = 0
        ax = fig.add_subplot(3,3,i)
        print i
        clr = colors[ind]
        sample_size = s_size
        zeros = True    
        D = {}
        name = 'divide_and_conquer'
        partitions = parts.rand_partitions(q, n, sample_size, name, D, zeros)
            
        skews = []
        for partition in partitions:
            skews.append(stats.skew(partition))
        
        D = get_kdens(skews)
        plt.plot(D[0],D[1],color = clr,lw=3, alpha = 0.99,label= 'n='+str(n))
        plt.xlabel("Skewnness",fontsize=8)
        
        plt.tick_params(axis='both', which='major', labelsize=5)
        plt.ylabel("Density",fontsize=8)
        i+=1
                
        ct = 0
        ax = fig.add_subplot(3,3,i)
        print i
        clr = colors[ind]
        sample_size = s_size
        zeros = False    
        D = {}
        name = 'divide_and_conquer'
        
        RADs = parts.rand_partitions(q, n, sample_size, name, D, zeros)
            
        ranks = range(1,n+1)
        max_ab = 0
        varlist = []
        for rad in RADs:
            log_rad = list(np.log(rad))
            
            variance = np.var(rad, ddof=1)
            varlist.append(variance)
            
            if max(rad) > max_ab: max_ab = max(rad)
            plt.plot(ranks,rad, color=colors[ind], lw=1.0,alpha=0.04)    
        print ' log(mean) vs. log(variance):', np.log(q/n), np.log(np.mean(varlist))
        plt.tick_params(axis='both', which='major', labelsize=5)
        plt.yscale('log')
        plt.xlabel("Rank",fontsize=8)
        plt.ylabel("Abundance",fontsize=8)
        i+=1
        
        
    plt.subplots_adjust(wspace=0.4, hspace=0.12)    
    plt.savefig('SSADfig-'+str(q)+'-'+str(n)+'.png', dpi=600, bbox_inches = 'tight', pad_inches=0.01)
    print 'done'
    
    
#kdens_unbias() # figures 2 Locey & McGlinn (2013), as well as figs 1 and 2 of the appendix 
#print 'Fig 2: finished'            
#time_trials_sage() # figure 3 Locey & McGlinn (2013)
#print 'Fig 3: finished'
winner() # figure 4 Locey & McGlinn (2013)
print 'Fig 4: finished'
#get_all_SSADs([1000,[100,500]])
#time_trials_bigQ() # figure 3 in the appendix of Locey & McGlinn (2013)





