Code that is currently underdevelopment for...you'll find out later


def last_part(Q,N):
    _last = [int(floor(float(Q)/float(N)))]*N
    _remainder = int(Q%N)
    
    j = 0
    while _remainder > 0:
        _last[j] += 1
        _remainder -= 1
        j += 1
    return _last


def last_part_minmax(Q,N,K):
    part = [K]
    Q -= K
    N -= 1
    
    _last = [int(floor(float(Q)/float(N)))]*N
    _remainder = int(Q%N)
    
    j = 0
    while _remainder > 0:
        _last[j] += 1
        _remainder -= 1
        j += 1
    
    part.extend(_last)
    
    return part
   

def int_to_list(i):
    return [int(x) for x in str(i)]

def int_to_list_fill(i,fill):
    return [x for x in str(i).zfill(fill)]

def list_to_int(l):
    return "".join(str(x) for x in l)
  
def part_to_num(part,fill):
    p_list = []
    for p in part:
        if len(int_to_list(p)) != fill:
            l = int_to_list_fill(p,fill)
            p = list_to_int(l)
        p_list.append(p)
    _list = list_to_int(p_list)
    return _list  
    
def num_to_part(num,fill,S):
    
    _list = int_to_list(num)
    if len(_list) != fill*S:
        ct = fill*S - len(_list) 
        while ct > 0:    
            _list.insert(0,0)
            ct -= 1    
    new_list1 = []
    new_list2 = []
    for i in _list:
        new_list1.append(i)
        if len(new_list1) == fill:
            new_list2.append(new_list1)
            new_list1 = []
    
    part = []
    for i in new_list2:
        j = int(list_to_int(i))
        part.append(j)
        
    return part
    
    
    
Q = int(25)
N = int(4)
print 'P(Q,N) = ',number_of_partitions(Q,N)

first_part = list(Partitions(Q,length=N).first())
fill = len(int_to_list(max(first_part)))
first_num = int(part_to_num(first_part,fill))

last_part = last_part(Q,N)
last_num = part_to_num(last_part,fill)
big_diff = int(first_num - int(last_num))

min_max = int(last_part[0])
print first_num,'\n',last_num,'\n',int(last_num),'\n',big_diff,'\n'

big_list = []

for i in range(min_max,Q-N):
    part = list(Partitions(Q, length=N, max_part=i).first())
    _num = int(part_to_num(part,fill))
    previous_part = last_part(Q,N,i)
    
    _diff = first_num - _num
    print _diff
    big_list.append(_diff)
    

_avg = round(mean(big_list))    
print 'avg among largest:',_avg

_len = len(int_to_list(first_num))
nines = [9,9]
zeros = [0]*(_len - 4)
nines.extend(zeros)
nines = float(list_to_int(nines))

print big_diff, nines

_sum = 0
while nines >= 99:
    
    _sum += (big_diff/nines) - ((big_diff%nines)/nines)
    print _sum
    
    remainder = (big_diff%nines)/nines
    big_diff 
    
    nines = int_to_list(int(nines))
    if len(nines) > 2:
        nines.remove(0)
        nines.remove(0)
    else: nines = 0.0
    nines = float(list_to_int(nines))
    
print _sum
