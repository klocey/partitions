
void NrParts(int *n, int *k, double *p) {
  int i, m ;
  for (i = 2 ; i < (*k + 1) ; i++) {
    if ((i + 1) <= (*n - i + 1 + 1)) {
      for (m = (i + 1) ; m < (*n - i + 1 + 1) ; m++) {
        p[m] = p[m] + p[m - i] ;
      }
    }
  }
}

void conjugate(int *l, int *j, int *part, int *conj) {
  int i, k, times ;
  for (i = (*l - 1) ; i > 0; i--) {
    times = part[i - 1] - part[i] ;
    if (times > 0) {
      for(k = *j ; k < *j + times ; k++) {
        conj[k - 1] = i ;
      }
      *j += times ;
    }
  }
}
