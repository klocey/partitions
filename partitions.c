
void NrParts(int *n, int *k, int *p) {
  int i, m ; 
  for (i = 2 ; i < (*k + 1) ; i++) {
    if ((i + 1) <= (*n - i + 1 + 1)) {
      for (m = (i + 1) ; m < (*n - i + 1 + 1) ; m++) {
        p[m] = p[m] + p[m - i] ; 
      }
    }
  }  
}