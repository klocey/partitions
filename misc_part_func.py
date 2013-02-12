# place to store misc. functions, little organization as of yet


Q = 10 # the total
N = 2  # number of parts
K = 8  # largest part
i = 0
#print Partitions(Q,length=N,max_part=K).cardinality() # Sage is very slow at this

Q = 4 # the total
N = 2 # number of parts
K = 2 # largest part

i = 0
#print Partitions(Q,length=N,max_part=K).cardinality() # Sage is very slow at this

print Q,N,K,'\n'
print 'i Q  N K  _sum'

def parts_nsx(i,Q,N,K):
    print i,Q,N,K
    if Q==0 and N==0:
        print 'found one, K=',K
        return 1
    if Q<=0 or N<=0 or K<=0:
        print 'change i'
        return 0
    if Q>0 and N>0 and K>0:
        _sum = 0
        for i in range(0,N+1):
            print 'K=',K,'decrease K'
            _sum += parts_nsx(i,Q-i*K, N-i, K-1)
            
        return _sum
            
print parts_nsx(i,Q,N,K)



######################### following all pertains to accomplishing one task

def most_even_partition(n,s):
    most_even = [int(floor(float(n)/float(s)))]*s
    _remainder = int(n%s)
    
    j = 0
    while _remainder > 0:
        most_even[j] += 1
        _remainder -= 1
        j += 1
    return most_even


def portion(alist, indices):

    return [alist[i:j] for i, j in zip([0]+indices, indices+[None])]

def next_restricted_part(p,n,s):
    
    if p == most_even_partition(n,s):return Partitions(n,length=s).first()
    
    for i in enumerate(reversed(p)):
        if i[1] - p[-1] > 1:
            if i[0] == (s-1):
                return Partitions(n,length=s,max_part=(i[1]-1)).first()
            else:
                parts = portion(p,[s-i[0]-1]) # split p into the part that won't change and the part that will
                h1 = parts[0]
                h2 = parts[1]
                next = list(Partitions(sum(h2),length=len(h2),max_part=(h2[0]-1)).first())
                return h1+next
                
N = 10
S = 2
p = list(Partitions(N,length=S).first())

for i in range(1,6):
    p = next_restricted_part(p,sum(p),len(p))         
    print p


