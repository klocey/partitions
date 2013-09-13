from sage.all import * # Sage is necessary for generating figure 2 of Locey and McGlinn (2013)
# Can be commented out if not generating figure 2 or if the user does not have Sage installed
# If Sage is installed, then from the command line: sage -python Locey_McGlinn_2013.py 
import sys

sys.path.append("/home/kenlocey/modules/pymods")
import macroecotools
sys.path.append("/home/kenlocey/modules/FEASIBLE_FUNCTIONS")
import feasible_functions as ff
#sys.path.append("/home/kenlocey/modules/partitions")
import partitions as parts
sys.path.append("/home/kenlocey/partitions/metrics")
import metrics

from os import path, access, R_OK  # W_OK for write permission
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from scipy.stats import gaussian_kde
import os
import  matplotlib.pyplot as plt
from pylab import *
import numpy as np
from scipy import stats
from scipy.stats import kstest
import random
from random import choice
import re
import math
import random, decimal
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
        
        
def rand_part_Sage(Q,N,sample_size):
    parts = []
    while len(parts) < sample_size:
        part = Partitions(Q).random_element()
        if len(part) == N:
            parts.append(part)
    
    return parts


def time_trials(  ):

    """ Compare the speed of the random partitioning functions developed
    in Locey and McGlinn (2013), for cases when zero values are and
    are not allowed.The code generates figure ? of the Appendix from
    Locey and McGlinn (2013) """

    fig = plt.figure()
    
    algorithms = ['multiplicity','top_down','divide_and_conquer','bottom_up']
    LABELS = ['Multiplicity:','Top down:','Divide and conquer:','Bottom up:']
    COLORS = ['#00CED1','#FF1493','k','gray']
    Qs = [100,200,300,400]
    sample_size = 10
    
    for i, Q in enumerate(Qs):
        
        ax = fig.add_subplot(2,2,i+1)
        print '\n',Q
            
        for ii, name in enumerate(algorithms):
            times = [] 
            Ns = range(int(Q/20.0),int(0.8*Q),int(Q/20.0)) # S values
                
            for N in Ns:
                zeros = False    
                D = {}
                with Timer() as t:
                    x = parts.rand_parts(Q,int(N),sample_size,name,D,zeros)
                times.append(round(t.interval,2))
                
            print LABELS[0],Q,times
            plt.plot(Ns,times,lw=3,color=COLORS[ii],label=LABELS[ii])
            
        plt.tick_params(axis='both', which='major', labelsize=8)
        plt.xlabel('N, (Q='+str(Q)+', n='+str(sample_size)+')',fontsize=8)
            
        if i == 0:
            plt.ylabel("Seconds",fontsize=8)
            legend(bbox_to_anchor=(0.0, 1.1, 2.25, .2), loc=10, ncol=4, mode="expand",prop={'size':9})#, borderaxespad=0.)
        
    plt.savefig('/home/kenlocey/time_trials.png', dpi=500)#, pad_inches=0)
    return


def time_trials_sage(  ):

    """ Compare the speed of the random partitioning function of Sage to the
    functions developed in Locey and McGlinn (20??), for cases when zero
    values are and are not allowed. The code also shows a tradeoff in speed
    between the two methods (i.e. 'divide and conquer', 'bottum-up')
    developed in L&M (20??).   
    The code generates figure 2 of Locey and McGlinn (2013) """


    fig = plt.figure()
    algorithms = ['multiplicity','top_down','divide_and_conquer','bottom_up']
    LABELS = ['Multiplicity:','Top down:','Divide and conquer:','Bottom up:']
    COLORS = ['#00CED1','#FF1493','k','gray']
    Qs = [50,100,150,200]
    sample_size = 300
    
    for i, Q in enumerate(Qs):
    
        ax = fig.add_subplot(2,2,i+1)
        print '\n',Q
        
        Q = Qs[i]
        step = int(Q/10.0)
        Ns = range(step,int(0.5*Q)+step,step)
        
        Sage_times = []
        for N in Ns:    
            print 'Sage',N
            with Timer() as t:
                x = rand_part_Sage(Q,N,sample_size)
            Sage_times.append(round(t.interval,2))
        print 'Sage',Q,Sage_times,'\n'
        
        for ii, name in enumerate(algorithms):
            times = []
            for N in Ns:
                zeros=False
                D = {}    
                with Timer() as t:
                    x = parts.rand_parts(Q, int(N), sample_size, name, D, zeros)
                times.append(round(t.interval,2))
            
            Y = []
            for iii, t in enumerate(times):
                if Sage_times[ii] > 0.0:
                    Y.append(float(Sage_times[iii])/t)
                else:
                    Y.append(1.0)
            
            plt.plot(Ns,Y,lw=3,color=COLORS[ii],label=LABELS[ii])
            plt.yscale('log')
            
        if i == 0:
            plt.text(1.5,11,'Sage/algorithm',fontsize=10,rotation='90')
            legend(bbox_to_anchor=(-0.03, 1.1, 2.25, .2), loc=10, ncol=4, mode="expand",prop={'size':9})#, borderaxespad=0.)
        
        plt.tick_params(axis='both', which='major', labelsize=8)
        plt.xlabel('N, (Q='+str(Q)+')',fontsize=10)
        
    plt.savefig('/home/kenlocey/time_trials_sage.png', dpi=500)#, pad_inches=0)
    return



def time_trials_bigQ(  ):
    
    """ Compare the speed of the random partitioning functions developed
    in Locey and McGlinn (2013), for cases when zero values are and
    are not allowed. The code generates figure ? of the Appendix from
    Locey and McGlinn (2013) """
    
    fig = plt.figure()
    
    sample_size=1
    zeros = False
    Qs = [5000,20000,200000,500000]
    
    for i, Q in enumerate(Qs):
        #print '\n',Q
        ax = fig.add_subplot(2,2,i+1)
        if i==0:
            algorithms = ['multiplicity','top_down','divide_and_conquer','bottom_up']
            LABELS = ['Multiplicity:','Top down:','Divide and conquer:','Bottom up:']
            COLORS = ['#00CED1','#FF1493','k','gray']
            Ns = [10,20,40,80]
           
        if i==1:
            algorithms = ['multiplicity','top_down','divide_and_conquer']
            LABELS = ['Multiplicity:','Top down:','Divide and conquer:']
            COLORS = ['#00CED1','#FF1493','k']
            Ns = [20,40,80,160]
        elif i==2:
            algorithms = ['multiplicity','top_down']
            LABELS = ['Multiplicity:','Top down:']
            COLORS = ['#00CED1','#FF1493']
            Ns = [40,80,160,320]
            
        elif i==3:
            algorithms = ['multiplicity']
            LABELS = ['Multiplicity:']
            COLORS = ['#00CED1']
            Ns = [80,160,320,640]
            
        for ii, name in enumerate(algorithms):
            times = []
            for N in Ns:
                D = {}
                with Timer() as t:
                    x = parts.rand_parts(Q, int(N), sample_size, name, D, zeros=False)
                times.append(round(t.interval,2))
            
            plt.plot(Ns,times,color=COLORS[ii],lw=3,label=LABELS[ii])
            #print name,times
            
        plt.xlabel('N, (Q='+str(Q)+', n='+str(sample_size)+')',fontsize=8)    
        plt.tick_params(axis='both', which='major', labelsize=8)
        
        if i==0 or i==2:
                plt.ylabel("Seconds",fontsize=10)
                if i==0:
                    legend(bbox_to_anchor=(0.0, 1.1, 2.2, .2), loc=10, ncol=4, mode="expand",prop={'size':9})#, borderaxespad=0.)
    
    plt.savefig('/home/kenlocey/time_trials_bigQ.png', dpi=500)#, pad_inches=0)
    return


def winner():
    
    """ Reveal how the random partitioning functions developed in Locey and McGlinn (2013),
    compare across many combinations of Q and N for cases when zero values are and are not
    allowed. The code generates figure 3 of the Appendix from Locey and McGlinn (2013) """
    
    algorithms = [['multiplicity','#00CED1'], ['top_down','#FF1493'], ['divide_and_conquer','k'], ['bottom_up','gray']]
    Qs1 = range(20,100,5)
    Qs2 = range(100,1000+20,20)
    Qs = Qs1+Qs2
    zeros = [True, False]

    BigList = [[],[]]
    for zero in zeros:
    
        for Q in Qs:
        
            bestTime = 10**10 # starting with a ridiculously high time that is easy to beat
            bestAlg = None 
            
                        
            for alg in algorithms:
                if zero == True:
                    DATA = open('/home/kenlocey/partitions/time_files/timeFiles/Python_' + alg[0] + '_zeros_Q=' + str(Q) + '.txt','r')
                    for d in DATA:
                        n = int(re.findall(r'\d*\S\d*',d)[0])
                        t = float(re.findall(r'\d*\S\d*',d)[1])
                        BigList[0].append([alg[1], Q, n, t]) # assigns a color value representing the algorithm
            
                elif zero == False:
                    DATA = open('/home/kenlocey/partitions/time_files/timeFiles/Python_' + alg[0] + '_Q=' + str(Q) + '.txt','r')
                    for d in DATA:
                        n = int(re.findall(r'\d*\S\d*',d)[0])
                        t = float(re.findall(r'\d*\S\d*',d)[1])
                        BigList[1].append([alg[1], Q, n, t]) # assigns a color value representing the algorithm
            
                        
                DATA.close()

    fig = plt.figure()
    x = 1
    for _list in BigList:
        print 'plot:',x
        ax = fig.add_subplot(2,2,x)
        D = {}
        for i in _list:
            colorVal = i[0]
            Q = i[1]
            N = i[2]
            t = i[3]

            if (Q, N) not in D:
                D[(Q, N)] = [colorVal, t]
         
            elif (Q, N) in D:
                if t < D[(Q, N)][1]:
                    D[(Q, N)] = [colorVal, t]
                elif t == D[(Q, N)][1]:
                    D[(Q, N)] = ['y', t]
    
        for key, value in D.items():
        
            if len(value) > 2: # this indicates there's been a tie between algorithms
                print 'There has been a tie between algorithms. Good-bye.'
                sys.exit()
            #    value.sort(key=lambda x: float(x[1])) # sort list of lists by times
            #    value.reverse()
            #    value = value[0]
            
            QN_combo = list(key)
            Q = int(QN_combo[0])
            N = int(QN_combo[1])
            c = value[0]
            t = value[1]
        
            sz=0
            if Q <= 100: sz = 0.75
            elif Q <= 300: sz = 1
            elif Q <= 600: sz = 2
            else: sz = 3
            plt.scatter(Q, N, marker='s', s=sz, color=c, edgecolor=c,alpha=0.9)
    
        plt.axis([20,1000,2,1000])
        plt.tick_params(axis='both', which='major', labelsize=8)
    
        if x == 1:
            plt.text(-170,800,'Number of parts, N',fontsize=10,rotation='90')
            plt.text(1000,-170,'Total, Q',fontsize=10)
            for alg in algorithms:
                plt.plot(-1, -1, lw=3, color=alg[1], label=alg[0])
            plt.legend(bbox_to_anchor=(-0.02, 1.00, 2.24, .2), loc=10, ncol=4, mode="expand",prop={'size':9})#, borderaxespad=0.)
        
        # insert inset of low Q values
        if x == 1: a = plt.axes([0.16, 0.75, .13, .13], axisbg='w')
        elif x == 2: a = plt.axes([0.58,0.75, .13, .13], axisbg='w')
        for key, value in D.items():
            QN_combo = list(key)
            Q = int(QN_combo[0])
            if Q <= 100:
                N = int(QN_combo[1])
                c = value[0]
                t = value[1]
                plt.scatter(Q, N, marker='s', s=3, color=c, edgecolor=c,alpha=0.9)
    
        plt.setp(a, xlim=(20,101), ylim=(2,101), xticks=[20,40,60,80,100], yticks=[20,40,60,80,100])
        plt.tick_params(axis='both', which='major', labelsize=7)
    
        x+=1
       
    plt.savefig('/home/kenlocey/algs_winner.png', dpi=500, bbox_inches = 'tight', pad_inches=0.4)
    plt.close()


def kdens(): 

    """ The code below compares random partitioning nplottions of Sage and Locey and McGlinn (2013)
        to full feasible sets. These analyses confirm that the algorithms are unbiased. The code
        uses full feasible sets, the random partitioning function in Sage, and the random partitioning
        for cases when 0' are or are not allowed.
            
        The code generates figures 1-4 of the Appendix of Locey and McGlinn (2013) """

    algs = ['multiplicity','top_down','divide_and_conquer','bottom_up']
    colors = ['#00CED1','#FF1493','k','gray']

    fig = plt.figure()
    nplot = 1 # a variable used to designate which nplottions and analyses are used for particular subplots
    sample_size = 500 # min number of macrostates needed to safely capture distributional features
                      # across the feasible set

    while nplot <= 4:
        ax =fig.add_subplot(2,2,nplot)

        if nplot < 3:
            Q = 40
            N = 10
        else:
            Q = 200 # The full feasible set can't be generated for N = 500 & S = 50
            N = 30
    
        for i in range(1,6):
            partitions = []
            for i, alg in enumerate(algs):
                if nplot == 1 or nplot == 3:
                    zeros = False
                    D = {}
                    partitions = parts.rand_parts(Q, N, sample_size, alg, D, zeros)
                else:
                    D = {}
                    zeros = True
                    partitions = parts.rand_parts(Q, N, sample_size, alg, D, zeros)
                #D = metrics.get_kdens_obs_Evar(partitions) # evenness
                #D = metrics.get_kdens_obs_gini(partitions) # inequality
                #D = metrics.get_kdens_obs_var(partitions) # variance
                #D = metrics.get_kdens_obs_skew(partitions) # skewness
                D = metrics.get_kdens_obs_MD(partitions)  # median summand
                plt.xlim(min(D[0]),max(D[0]))
                plt.plot(D[0],D[1],color=colors[i],lw=0.7)
        
        if nplot == 1: # using the full feasible set, no zero values (i.e. proper integer partitions)
            for i in range(1,2):
                partitions = []
                for p in Partitions(Q,length=N):
                    partitions.append(p)
        
            #D = metrics.get_kdens_obs_Evar(partitions) # evenness
            #D = metrics.get_kdens_obs_gini(partitions) # inequality
            #D = metrics.get_kdens_obs_var(partitions) # variance
            #D = metrics.get_kdens_obs_skew(partitions) # skewness
            D = metrics.get_kdens_obs_MD(partitions)  # median summand
            plt.xlim(min(D[0]),max(D[0]))
            plt.plot(D[0],D[1],color='r',lw=3,alpha=0.5)
        
        elif nplot == 2: # using the full feasible set, zero values included
            for i in range(1,2):
                partitions = []    
                n = 1
                while n <= N:
                
                    numparts = parts.NrParts(Q,n)    
                    part = parts.firstpart(Q,n,None)
                    ct2 = 0
                    while ct2 < numparts:
                    
                        part = parts.next_restricted_part(part)
                        if len(part) == N: partitions.append(part) 
                        else:
                            part2 = list(part)
                            zeros = [0]*(N-len(part))
                            part2.extend(zeros)
                            partitions.append(part2)
                        
                        #print nplot,numparts-ct2
                        ct2+=1
                    n+=1
            #D = metrics.get_kdens_obs_Evar(partitions) # evenness
            #D = metrics.get_kdens_obs_gini(partitions) # inequality
            #D = metrics.get_kdens_obs_var(partitions) # variance
            #D = metrics.get_kdens_obs_skew(partitions) # skewness
            D = metrics.get_kdens_obs_MD(partitions)  # median summand
            plt.xlim(min(D[0]),max(D[0]))
            plt.plot(D[0],D[1],color='r',lw=3,alpha=0.5)
    
        elif nplot == 3: 
            for i in range(1,6):
                partitions = []
                while len(partitions) < sample_size: # Use the random partition nplottion in Sage to generate a sample of partitions for N and S
                    part = Partitions(Q).random_element()
                    if len(part) == N:
                        partitions.append(part)
                       #print nplot,i, sample_size - len(partitions)   
                    else:
                        part = list(Partition(part).conjugate())
                        if len(part) == N:
                            partitions.append(part)
                            #print nplot,i,sample_size - len(partitions)   
                #D = metrics.get_kdens_obs_Evar(partitions) # evenness
                #D = metrics.get_kdens_obs_gini(partitions) # inequality
                #D = metrics.get_kdens_obs_var(partitions) # variance
                #D = metrics.get_kdens_obs_skew(partitions) # skewness
                D = metrics.get_kdens_obs_MD(partitions)  # median summand
                plt.xlim(min(D[0]),max(D[0]))
                plt.plot(D[0],D[1],color='r',lw=0.7,alpha=0.5)
        
        
        elif nplot == 4:
            for i in range(1,6):
                partitions = []
                while len(partitions) < sample_size: # Use the random partition nplottion in Sage to generate a sample of partitions for N and S
                    part = list(Partitions(Q).random_element())
                    if len(part) == N:
                        partitions.append(part)
                    elif len(part) < N:
                        zeros = [0]*(N-len(part))
                        part.extend(zeros)
                        partitions.append(part)
                
                    #print nplot,i,sample_size - len(partitions)   
                #D = metrics.get_kdens_obs_Evar(partitions) # evenness
                #D = metrics.get_kdens_obs_gini(partitions) # inequality
                #D = metrics.get_kdens_obs_var(partitions) # variance
                #D = metrics.get_kdens_obs_skew(partitions) # skewness
                D = metrics.get_kdens_obs_MD(partitions)  # median summand
                plt.xlim(min(D[0]),max(D[0]))
                #if nplot == 4:  plt.ylim(0,0.01)
                #elif nplot == 3: plt.ylim(0,0.01)
                #else: plt.ylim(0,max(D[1])+0.02)
                plt.plot(D[0],D[1],color='r',lw=0.7,alpha=0.5)
                
                     
        if nplot == 1:
            plt.plot([0],[0], color='#00CED1', lw=2, label = 'Multiplicity')
            plt.plot([0],[0], color='#FF1493',lw=2, label='Top-down')    
            plt.plot([0],[0], color='k',lw=2, label='Divide & Conquer')
            plt.plot([0],[0], color='gray',lw=2, label='Bottom-up')
            plt.plot([0],[0], color='r',lw=2, label='FS Q='+str(Q)+',N='+str(N),alpha=0.5)
            plt.legend(bbox_to_anchor=(-0.02, 1.00, 2.24, .2), loc=10, ncol=5, mode="expand",prop={'size':8})#, borderaxespad=0.)
            
        if nplot == 1 or nplot == 3:
            plt.ylabel("pdf",fontsize=12)    
        
        if nplot == 3 or nplot == 4:
            plt.xlabel("Median",fontsize=12)
        
        print nplot
        nplot+=1
        
        plt.tick_params(axis='both', which='major', labelsize=8)
        
    plt.savefig('/home/kenlocey/partitions/kdens_Median'+str(sample_size)+'.png', dpi=500, pad_inches=0)



def get_obs_pred(dataset,pattern):
    
    PATH = '/home/kenlocey/data1/' + dataset + '/' + dataset + '_' + pattern + '_obs_pred.txt'
    if path.exists(PATH) and path.isfile(PATH) and access(PATH, R_OK):
        DATA = open(PATH,'r')
        ct1 = 0
        ct2 = 0
        d = DATA.readline()
        m0 = re.match(r'\A\S*',d).group() # site
        m2 = list(str(re.findall(r'\d*\s*\d*\S$',d)[0]).split())
        obsSAD = [int(m2[0])]
        predSAD = [int(m2[1])]
        
        for d in DATA:
            ct1+=1
            m1 = re.match(r'\A\S*',d).group()
            if m1 == m0:
                m2 = list(str(re.findall(r'\d*\s*\d*\S$',d)[0]).split())
                obsSAD.append(int(m2[0]))
                predSAD.append(int(m2[1]))
            
            else:
                site_name = m0
                m0 = m1
                ct2+=1
                SAD = []
                m2 = list(str(re.findall(r'\d*\s*\d*\S$',d)[0]).split())
                obsSAD =  [int(m2[0])]
                predSAD = [int(m2[1])]
        
        DATA.close()
        return([obsSAD, predSAD])


def get_pattern(dataset, pattern, obs):
    
    if obs == True:
        ind = 0
    else: ind = 1
    
    PATH = '/home/kenlocey/data1/' + dataset + '/' + dataset + '_' + pattern + '_obs_pred.txt'
    if path.exists(PATH) and path.isfile(PATH) and access(PATH, R_OK):
        DATA = open(PATH,'r')
        ct1 = 0
        ct2 = 0
        d = DATA.readline()
        m0 = re.match(r'\A\S*',d).group() # site
        m2 = list(str(re.findall(r'\d*\s*\d*\S$',d)[0]).split())
        SAD = [int(m2[ind])]
        
        for d in DATA:
            ct1+=1
            m1 = re.match(r'\A\S*',d).group()
            if m1 == m0:
                m2 = list(str(re.findall(r'\d*\s*\d*\S$',d)[0]).split())
                SAD.append(int(m2[ind]))
                    
            else:
                site_name = m0
                m0 = m1
                ct2+=1
                SAD = []
                m2 = list(str(re.findall(r'\d*\s*\d*\S$',d)[0]).split())
                SAD = [int(m2[ind])]
                    
        DATA.close()
        return(SAD)


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
    #print len(rank_vectors)
    for i, _list in enumerate(rank_vectors):
        if sum(_list) >= 100:
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

    
def obs_pred_rsquare(obs, pred):
    """Determines the prop of variability in a data set accounted for by a model
    
    In other words, this determines the proportion of variation explained by
    the 1:1 line in an observed-predicted plot.
    
    """
    return 1 - sum((obs - pred) ** 2) / sum((obs - np.mean(obs)) ** 2)


def histOutline(dataIn, *args, **kwargs): # found at http://www.scipy.org/Cookbook/Matplotlib/UnfilledHistograms
    (histIn, binsIn) = np.histogram(dataIn, *args, **kwargs)
    
    stepSize = binsIn[1] - binsIn[0]
    bins = np.zeros(len(binsIn)*2 + 2, dtype=np.float)
    data = np.zeros(len(binsIn)*2 + 2, dtype=np.float)
    for bb in range(len(binsIn)):
        bins[2*bb + 1] = binsIn[bb]
        bins[2*bb + 2] = binsIn[bb] + stepSize
        if bb < len(histIn):
            data[2*bb + 1] = histIn[bb]
            data[2*bb + 2] = histIn[bb]
    bins[0] = bins[1]
    bins[-1] = bins[-2]
    data[0] = 0
    data[-1] = 0
    
    return (bins, data)


def rank_to_PrestonPlot(RAD):

    bins = [0,1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768]
    SAD = [0]*len(bins)
    for p in RAD:
        if p > 16384: SAD[16]+=1
        elif p > 8192: SAD[15]+=1
        elif p > 4096: SAD[14]+=1
        elif p > 2048: SAD[13]+=1
        elif p > 1024: SAD[12]+=1
        elif p > 512: SAD[11]+=1
        elif p > 256: SAD[10]+=1
        elif p > 128: SAD[9]+=1
        elif p > 64: SAD[8]+=1
        elif p > 32: SAD[7]+=1
        elif p > 16: SAD[6]+=1
        elif p > 8: SAD[5]+=1
        elif p > 4: SAD[4]+=1
        elif p > 2: SAD[3]+=1
        elif p > 1: SAD[2]+=1
        elif p == 1: SAD[1]+=1
        elif p == 0: SAD[0]+=1
        
    SAD2 = [SAD[0],SAD[1],SAD[2]]
    for i, freq in enumerate(SAD):
        if i > 2:
            x = freq/(bins[i]/2.0) 
            SAD2.append(x)
    return SAD2



def lose_trailing_zeros(_list):
    for dx in (0, -1):
        while _list and _list[dx] == 0:
            _list.pop(dx)
            
    return _list


            
def SAD_SSAD(fs):
    switch = 0
    dataset = 'BCI'
    fig = plt.figure()
        
    """ This does this ... figure 3 of Locey and McGilnn 2013"""
    
    ax =fig.add_subplot(2,2,1)
    pattern = 'SAD'
    RADs = get_obs_pred(dataset, pattern)
    obsRAD = RADs[0]
    predRAD = RADs[1]
    print obs_pred_rsquare(np.log(obsRAD), np.log(predRAD))
    N = sum(obsRAD)
    S = len(obsRAD)
    
    DATA = open('/home/kenlocey/combined1/' + str(N) + '-' + str(S) + '.txt','r')
    RADs = []
    for d in DATA:
        rad = eval(d)
        RADs.append(rad)
    DATA.close()
    
    x = []
    y = []
    for rad in RADs:
        rank = 1
        for j in rad:
            x.append(rank)
            y.append(np.log(j))
            rank += 1

    plt.hexbin(x,y,bins='log',gridsize=100,mincnt=1,cmap=cm.Reds) # bins can be 'log' or None
    plt.tick_params(axis='both', which='major', labelsize=fs-3)
    ranks = range(1,301+1)
    plt.scatter(ranks,np.log(obsRAD), color='0.3', marker = '.')
    plt.text(225, 8, '0.28',fontsize=fs+2)
    plt.xlabel("Rank in abundance",fontsize=fs)
    plt.ylabel("ln(abundance)",fontsize=fs)
    #plt.ylim(-0.1,19.5)
    plt.xlim(0.01,310)
    
    
    ax =fig.add_subplot(2,2,2)
    
    SSADs_freqs = get_SSADs(dataset)
    SSADs = hist_to_rank(SSADs_freqs)
    
    Pvals = []
    for i, obsSSAD in enumerate(SSADs):
        
        Q = sum(obsSSAD)
        N = len(obsSSAD)
        
        partitions = []
        PATH = '/home/kenlocey/combined1/' + str(Q) + '-' + str(N) + '.txt'
        if path.exists(PATH) and path.isfile(PATH) and access(PATH, R_OK):
            DATA = open(PATH,'r')
            DATA.readline() # first line has no abundance info, so skip it
            for d in DATA:
                partitions.append(eval(d))
        else: continue
        
        if len(partitions) == 0:
            print Q,N,'something is wrong'
            sys.exit()  
              
        predSSAD = parts.central_tendency(partitions)
        print sum(obsSSAD),len(obsSSAD),' ',sum(predSSAD),len(predSSAD),
        
        pval = stats.ks_2samp(obsSSAD, predSSAD)[1]
        Pvals.append(pval)
        print pval
        
        if pval > 0.05 and Q > 1000 and switch == 0:
            
            SSADs = []
            for part in partitions:
                ssad = rank_to_PrestonPlot(part) 
                ssad2 = lose_trailing_zeros(ssad)
                SSADs.append(ssad2)
            
            #print SSADs[0]
            #sys.exit()
            
            x = []
            y = []
            for ssad in SSADs:
                ab_class = 0
                for j in ssad:
                    x.append(ab_class)
                    y.append(j)
                    ab_class += 1
            
            plt.hexbin(x,y,bins=None,gridsize=25,mincnt=1,cmap=cm.Reds) 
            obsSSAD2 = rank_to_PrestonPlot(obsSSAD)
            obsSSAD3 = lose_trailing_zeros(obsSSAD2)
            bins = range(0,len(obsSSAD2))
            #plt.plot(obsSSAD, color='0.3', lw=2)
            plt.scatter(bins, obsSSAD3, color='0.3', marker = 'o')
            #plt.xlim(-0.5,max(bins))
            #plt.ylim(-0.5,max(max(predSSAD),max(obsSSAD))+10)
            plt.xlabel("Abundance class",fontsize=fs)
            plt.ylabel("frequency",fontsize=fs)
            plt.tick_params(axis='both', which='major', labelsize=fs-3)
            
            switch += 1
            
        if len(Pvals) > 50 and switch == 1: break
    # Create inset for histogram of site level r^2 values
    axins = inset_axes(ax, width="30%", height="30%", loc=1)
    D = get_kdens(Pvals)
    plt.plot(D[0],D[1],color='0.2',lw=2.0)
    plt.setp(axins, xticks=[0.0, 0.5, 1.0])
    plt.tick_params(axis='both', which='major', labelsize=fs-4)
    plt.xlim(0.0,1.0)
    #plt.ylim(0.0, max(D[1])+0.003)
    plt.subplots_adjust(wspace=0.3)
    plt.savefig('/home/kenlocey/partitions/SAD-SSAD.png', dpi=400, bbox_inches = 'tight', pad_inches=0.03)

#fs = 10 # fontsize    
#SAD_SSAD(fs)

#time_trials_sage() # figure 2 Locey & McGlinn (2013)

#winner() # figure 3 Locey & McGlinn (2013)

#time_trials() # figure ? Appendix of Locey & McGlinn (2013)

#time_trials_bigQ() # figure ? Appendix of Locey & McGlinn (2013)

#kdens() # figures 1-4 in the Appendix of Locey & McGlinn (2013)







