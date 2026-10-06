#ifndef INCLVectorBench_hh
#define INCLVectorBench_hh 1

/**
 * Benchmark different three vector implementations.
 */
template< class VectorT >
class INCLVectorBench {

public:
  INCLVectorBench()
  { };
  ~INCLVectorBench()
  { };

  bool test() {
    std::cout <<"INCLVectorBench< T > : Running the tests" << std::endl;
    VectorT v1(1.0, 2.0, 3.0);
    VectorT v2(2.0, -54.0, -443.0);
    VectorT v3(-22.0, 5.0, 6663.0);
    VectorT v4(6.777, 2.14, 0.0056);
    for(long i = 0; i != 1000000; i++) {
      VectorT r1 = v1 + v2;
    }

    for(long i = 0; i != 1000000; i++) {
      VectorT r1 = v1 + v2 + v3 + v4 + (v1 + v2)*(2.0 - 1.0 * i);
      VectorT r2 = r1 + v4 - v2; 
    }

    for(long i = 0; i != 1000000; i++) {
      VectorT r1 = v1 + v2 + v3 + v4 + (v1 + v2)*(2.0 * i);
    }


    for(long i = 0; i != 1000000; i++) {
      VectorT r1 = v1 + v2 + v3 + v4 + (v1 + v2)*(2.0 * i);
    }


    for(long i = 0; i != 1000000; i++) {
      VectorT r1 = v1 + v2 + v3 + v4 + (v1 + v2)*(2.0 * i);
    }

    std::cout <<"INCLVectorBench< T > : Tests complete" << std::endl;
  }
};

#endif
