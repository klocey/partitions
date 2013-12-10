
void NrParts(int *Q, int *N, double *p) {
  int i, m ;
  for (i = 2 ; i < *N + 1 ; i++) {
    for (m = (i + 1) ; m < (*Q - i + 1 + 1) ; m++) {
      p[m] = p[m] + p[m - i] ;
    }
  }
}


void conjugate(int *l, int *j, int *partition, int *conj) {
  int i, k, times ;
  for (i = (*l - 1) ; i > 0; i--) {
    times = partition[i - 1] - partition[i] ;
    if (times > 0) {
      for(k = *j ; k < *j + times ; k++) {
        conj[k - 1] = i ;
      }
      *j += times ;
    }
  }
}
