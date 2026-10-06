#ifndef INCL_ANALYSIS_WRAPPER
#define INCL_ANALYSIS_WRAPPER 1

#include "G4INCLHistogram1D.hh"
#include <iostream>
#include <vector>

extern "C" {
  int histoID;
  std::vector<G4INCL::Analysis::Histogram1D*> *histos;
  void analysisinit_();
  void h1d_(int *nbins, double *xmin, double *xmax, int *id);
  void h1dfill_(int *id, double *x);
  void h1dprintdata_(int *id);
}

#endif
