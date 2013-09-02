## functions for generating results for the BCI dataset

make_comm_matrix = function(spnum, S, coords, n_quadrats, domain, abu = NULL,
                            grainSuffix=NULL)
{ 
  ## Output: 
  ## A community matrix where each row is a differnet pixel on a grid.  
  ## Arguments:
  ## spnum : an integer specifying species identities
  ## S : the size of the species pool may be larger than the number of unique 
  ##     spnum
  ## coords : two column matrix (x,y) specifying the spatial coordinates of each stem
  ## n_quadrats : the number of quadrats at each spatial grain
  ## domain : specifies the spatial domain of the area:  (xmin, xmax, ymin, ymax)
  ## abu: abundance associated with each record, if NULL then it is set to 1
  ##      individual per record
  ## grainSuffix : if supplied the grain column will have this appended to it
  ##               so that it is clear what community this corresponds with
  xdiff = abs(domain[2] - domain[1])
  ydiff = abs(domain[4] - domain[3])
  if (xdiff > ydiff) {
    xlengths = xdiff / n_pixels_long(log2(n_quadrats))
    ylengths = ydiff / n_pixels_wide(log2(n_quadrats))
  }  
  else if (xdiff < ydiff) {
    xlengths = xdiff / n_pixels_wide(log2(n_quadrats))
    ylengths = ydiff / n_pixels_long(log2(n_quadrats))
  }
  else if (xdiff == ydiff) {
    xlengths = ylengths = rep(NA, length(n_quadrats))
    for (i in seq_along(n_quadrats)) {
      if (log2(n_quadrats[i]) %% 2 == 1) {
        ## if # of bisections is odd then arbitrarily 
        ## make x dimension have more cells
        xlengths[i] = xdiff / n_pixels_long(log2(n_quadrats[i]))
        ylengths[i] = ydiff / n_pixels_wide(log2(n_quadrats[i]))
      }  
      else {
        xlengths[i] = xdiff / sqrt(n_quadrats[i])
        ylengths[i] = ydiff / sqrt(n_quadrats[i])
      } 
    }  
  }
  else
    stop('Function cannot figure out how to split up the area')
  comms = matrix(NA, nrow=sum(n_quadrats), ncol=S + 3)
  colnames(comms) = c('grain', 'x', 'y', paste('sp', 1:S, sep=''))
  irow = 1
  for (i in seq_along(n_quadrats)) {
    xbreaks = seq(domain[1], domain[2], xlengths[i])
    ybreaks = seq(domain[3], domain[4], ylengths[i]) 
    for (x in 1:(length(xbreaks) - 1)) {
      for (y in 1:(length(ybreaks) - 1)) {
        inQuad =  xbreaks[x] <= coords[ , 1] & coords[ , 1] < xbreaks[x + 1] & 
          ybreaks[y] <= coords[ , 2] & coords[ , 2] < ybreaks[y + 1]
        if (is.null(grainSuffix)) {
          comms[irow, c(1:3)] = c(as.numeric(paste(round(xlengths[i] * ylengths[i], 2), sep='')),
                                  x, y)
        }
        else {
          comms[irow, c(1:3)] = c(paste(round(xlengths[i] * ylengths[i], 2),
                                        grainSuffix, sep=''), x, y)
        }
        if (is.null(abu) ){
          comms[irow, -c(1:3)] = as.integer(table(c(spnum[inQuad],1:S)) - 1)
        }
        else {
          comms[irow, -c(1:3)] =  as.integer(table(c(unlist(mapply(
            rep, spnum[inQuad], abu[inQuad])), 1:S)) - 1)
        }  
        irow = irow + 1 
      }
    }
  }
  return(comms)
}

mat2psp = function(sp_mat, xy_coord, N=NULL, M=NULL)
{
  ##place site by species matrix (x) into an S x N x M array where N >= M
  ##a multidimensional array that eases computing 
  ##replaces the old function 'grid.pres'
  ##Note: if area rectangular the first dimension of xy_coord does NOT
  ## have to be larger than the second dimension
  if (is.null(N) )
    N = sqrt(nrow(sp_mat))
  if (is.null(M) )
    M = N
  if (nrow(sp_mat) != nrow(xy_coord))
    stop('Number of samples in species matrix must match number of samples in xy-coordinates')
  if (N < M)
    stop('N should be >= M')
  if (N * M != nrow(sp_mat))
    stop('Number of specified samples (N X M) must be equal to the number of samples in sp_mat')
  S = ncol(sp_mat)
  ## order rows of sp_mat based on xy_coords so that they are placed
  ## in the multidimensional array in the correct arrangement
  proper_order = order(xy_coord[ , 2], xy_coord[ , 1])
  sp_mat = sp_mat[proper_order, ]
  ## check that its numberic
  if (is.character(sp_mat[1,1]))
    sp_mat = matrix(as.numeric(sp_mat), nrow=nrow(sp_mat), ncol=ncol(sp_mat))
  psp = array(sp_mat, dim=c(N, M, S))
  psp = aperm(psp, c(3, 1, 2))
  spSums = apply(psp, 1, sum)
  ## drop species that never occur
  if(any(spSums %in% 0))
    psp = psp[spSums > 0, , ]
  return(psp)
}

n_pixels_long = function(i_bisect)
{
  ## returns the number of pixels on the side of a grid with more or equal pixels
  ## after i bisection events
  ## old function name: len
  2^floor((i_bisect + 1) / 2) 
}

n_pixels_wide = function(i_bisect)
{
  ## returns the number of pixels on the side of a grid with less or equal pixels
  ## after i bisecttion events
  ## old function name: wid
  2^floor(i_bisect / 2) 
} 

get_SAR = function(psp, grains, mv_window=FALSE)
{
  ## Purpose: to construct spatially explict SAR based upon a
  ## mapped grid of occurances
  ## this function replaces the older function 'grid.SAR'
  ## Arguments:
  ## psp: community array (i.e., S x N x M abundance array where N >= M)
  ## grains: the areas in pixels for which to compute the SAR
  ##         only grains that have integer log base 2 are considered
  ## mv_window: FALSE indicates that a non-moving window SAR will be calculated
  if (class(psp) != 'array')
    stop('psp must be a community array (S X N X M)')
  grains = grains[log2(grains) == round(log2(grains))]
  ## define the size of sampling units on each side
  lenN = n_pixels_long(log2(grains))
  lenM = n_pixels_wide(log2(grains))
  sr = rep(0, length(grains))
  ind = rep(0, length(grains))
  std = rep(0, length(grains))
  cs = rep(0, length(grains))
  S = dim(psp)[1]
  N = dim(psp)[2]
  M = dim(psp)[3]
  if (M > N) {
    stop('The first spatial dimension of psp must be larger than or equal to the second 
         (i.e. psp[S,N,M] where N >= M)')
  }       
  for (l in seq_along(grains)) {
    if (grains[l] == 1) {  # if area=1
      sr[l] = sum(psp > 0)
      std[l] = sd(as.vector(psp > 0))
      ind[l] = sum(psp)
      cs[l] = N * M
    }
    else{
      if (mv_window) {
        brksN = 1:(N - lenN[l] + 1)
        brksM = 1:(M - lenM[l] + 1)
      }
      else{
        brksN = seq(1, N, lenN[l])
        brksM = seq(1, M, lenM[l])
      }  
      sr_vec = NULL
      for (n in brksN) {
        for (m in brksM) {
          psp_tmp = psp[ , n:(n + (lenN[l] - 1)),
                        m:(m + (lenM[l] - 1))]
          if (S == 1) {
            sr_vec = c(sr_vec, any(psp_tmp > 0) * 1)
            ind[l] = ind[l] + sum(psp_tmp)
          }
          else {
            sr_vec = c(sr_vec, sum(apply(psp_tmp > 0, 1, sum) > 0))
            ind[l] = ind[l] + sum(apply(psp_tmp, 1, sum))
          }
          cs[l] = cs[l] + 1
        }
      }
      sr[l] = sum(sr_vec)
      std[l] = sd(sr_vec)
    }
  }
  out = cbind(grains, sr / cs, ind / cs, cs, std)  
  colnames(out) = c('grains', 'richness', 'indiv', 'count', 'sr_std')
  return(out)
} 

get_SSAD = function(comms)
{
  if (is.character(comms[1,1])) {
    tmp_comms = matrix(as.numeric(comms), nrow=nrow(comms), ncol=ncol(comms))
    colnames(tmp_comms) = colnames(comms)
    comms = tmp_comms
  }  
  grains = unique(comms[ , 1])
  for (g in seq_along(grains)) {
    tmp_comm = comms[comms[ , 1] == grains[g], -(1:3)]
    uni_abu = sort(unique(as.vector(tmp_comm)))
    ssad = apply(tmp_comm, 2, function(x) table(c(x, uni_abu)) - 1)
    ssad = cbind(grains[g], uni_abu, ssad)
    if (g == 1)
      out = ssad
    else
      out = rbind(out, ssad)
  }
  colnames(out)[1:2] = c('grain', 'abu')
  return(out)
}

get_OFD = function(comms)
{
  if (is.character(comms[1,1])) {
    tmp_comms = matrix(as.numeric(comms), nrow=nrow(comms), ncol=ncol(comms))
    colnames(tmp_comms) = colnames(comms)
    comms = tmp_comms
  }  
  grains = unique(comms[ , 1])
  for (g in seq_along(grains)) {
    tmp_comm = (comms[comms[ , 1] == grains[g], -(1:3)] > 0) * 1
    occ = colSums(tmp_comm)
    uni_occ = min(occ):max(occ)
    ofd = as.vector(table(c(occ, uni_occ)) - 1)
    ofd = cbind(grains[g], uni_occ, ofd)
    if (g == 1)
      out = ofd
    else
      out = rbind(out, ofd)
  }
  colnames(out) = c('grain', 'occ', 'num_sp')
  return(out)
}
