#ifndef G4INCLHistogram1D_hh
#define G4INCLHistogram1D_hh 1

#include <string>

namespace G4INCL {
  namespace Analysis {
    
    class Histogram1D {
    public:
      Histogram1D(int nbins, double xmin, double xmax);
      ~Histogram1D();

      void fill(double x, double weight=1.0);

      int getNbins() { return bins; };
      double getBinContent(int i) { return binContent[i]; };

      void normalize(double norm=1.0);

      std::string print();
      std::string printHistogramSteps();

    private:
      int bins;
      double min, max;
      double underflow, overflow;
      double *binEdges;
      double *binContent;
    };
  }
}

#endif
