
setwd('~/partitions')

source('bci_functions.R')

## clean up raw bci data
dat = read.table('./data/bci_census7.txt', sep='\t', header=TRUE)
  
goodData = dat$Status=='alive' & !is.na(dat$DBH) & 
           !is.na(dat$gx) & !is.na(dat$gy) & 
           dat$Latin != 'Unidentified species' &
           dat$Stem != 'secondary'

dat = dat[goodData ,]
  
write.csv(dat, file='./data/bci_census7_filtered.csv', row.names=F)
#dat = read.csv('./data/bci_census7_filtered.csv')

## generate a site x species matrix for each spatial scale
uniSpeciesNames = as.character(sort(unique(dat$Latin)))
dat$spnum = match(dat$Latin,uniSpeciesNames)
S = max(dat$spnum) 

## generate the empirical SAD files
sad = sort(table(dat$spnum), dec=T)

write.table(matrix(sad, nrow=1), file='./data/bci_sad.csv',
            sep=',', row.names=FALSE, col.names=FALSE)


## generate bci community matrix at several different grains
i_bisections = seq(5, 13, 2)
n_quadrats = 2^i_bisections
domain = c(0,1000,0,500) # spatial domain in meters defined here
comm = make_comm_matrix(dat$spnum, S, cbind(dat$gx, dat$gy), n_quadrats, domain)
write.csv(comm, file='./data/bci_comm.csv', row.names=FALSE)

## compute SSADs
## find which species has the max abu
for (i in i_bisections){
  true = comm[ , 1] == round((1e3 * 5e2) / 2^i, 2)
  ssad = get_SSAD(comm[true, ])
  #max_abu = max(colSums(comm[ , -(1:3)]))
  #sp_max = which(colSums(comm[ , -(1:3)]) == max_abu)
  #ssad = ssad[ , c(1:2, sp_max + 2)]
  write.csv(ssad, file=paste('./data/bci_ssad_', i,'.csv', sep=''),
            row.names=FALSE)
}


## compute OFD
ofd = get_OFD(comm)
write.csv(ofd, file='./data/bci_ofd.csv', row.names=FALSE)

## choose 100 random coordinates for the 64 m2 plots to compute SAR on
n = 100
len = 8
domain = c(0, 1000, 0, 500) # spatial domain in meters defined here
x_rand = runif(n, domain[1], domain[2] - len)
y_rand = runif(n, domain[3], domain[4] - len)

## for each coordinate sample an 8 x 8 m square at 6 grains,
## 1, 2, 4, 8, 16, 64 m2
comms = vector('list', n)
for(i in seq_along(x_rand)){
  i_bisections = log2(len^2):1
  n_quadrats = 2^i_bisections
  true = dat$gx >= x_rand[i] & dat$gx < x_rand[i] + len &
         dat$gy >= y_rand[i] & dat$gy < y_rand[i] + len
  dat_tmp = dat[true, ]
  uni_ids = unique(dat_tmp$spnum)
  S_tmp = length(uni_ids)
  spnum_tmp = match(dat_tmp$spnum, uni_ids)
  domain = c(x_rand[i], x_rand[i] + len, y_rand[i], y_rand[i] + len)
  comms_tmp = make_comm_matrix(spnum_tmp, S_tmp, cbind(dat_tmp$gx, dat_tmp$gy),
                               n_quadrats, domain)
  comms[[i]] = comms_tmp
}

save(comms, file='./data/bci_comms.Rdata')
#load('./data/bci_comms.Rdata')

## compute SAR
Amin = 1
for (i in seq_along(comms))
  comms[[i]] = (comms[[i]][comms[[i]][ , 1] == Amin, ])

bisect_fine = 6
Ns = n_pixels_long(bisect_fine)
Ms = n_pixels_wide(bisect_fine)
psp = vector('list',length(comms))
for(i in seq_along(comms))
  psp[[i]] = mat2psp(comms[[i]][ , -(1:3)], comms[[i]][ , 2:3], Ns, Ms)

## compute SARs, non-movinging window 
grains = 2^(0:bisect_fine)

## the finest grain is at 1 m2
grain_fine = 1

sar = vector('list',length(psp))
for (i in seq_along(psp)) {
  sar[[i]] = get_SAR(psp[[i]], grains)
  ## add area m2 column
  sar[[i]] = data.frame(sar[[i]], area = sar[[i]][ , 1] * grain_fine)
}

sar_avg = sar[[1]]
for (i in 2:length(sar)){
  sar_avg = sar_avg + sar[[i]]
}
sar_avg = sar_avg / length(sar)
sar_avg = sar_avg[ , c('area', 'richness')]

write.csv(sar_avg, file='./data/bci_sar.csv', row.names=FALSE)



  