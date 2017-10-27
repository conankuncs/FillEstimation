// On Linux and MacOS, you can compile and run this program like so:
//   g++ -std=c++11 -O3 -DNDEBUG -DTACO -I../taco/include/ spmv.cc -o spmv -L../taco/build/lib -ltaco -lm -ldl
//   LD_LIBRARY_PATH=../../build/lib ./spmv

#include <random>
#include <iostream>
#include <string>
#include "taco.h"

using namespace std;
using namespace taco;
using namespace storage;
int main(int argc, char* argv[]) {
  if(argc != 2) {
    cout<<"Usage: ./spmv matrix/[matrix name].mtx"<<endl;
    return 0;
  }
  Format csr({Dense,Sparse});
  // Declare a new tensor "A" of double-precision floats with dimensions
  // Read Compressed Sparse Fiber Tensor from argv[1] (input file)
  Tensor<double> A ("A", {3,3}, csr);
  A.insert({1,1},1);
  A.insert({2,1},-2);
  A.insert({2,2},1);
  A.insert({3,1},0.5);
  A.insert({3,3},1);

  /*

  Expected)

  Aptr[] = { 0, 1, 3, 5 };
  Aind[] = { 0, 0, 1, 0, 2 };
  Aval[] = { 1, -2, 1, 0.5, 1 };

  Result)

  Aptr[] = {3}
  Aind[] = {0,0,1,3}
  Aval[] = {1,1,2}

  It also has a bug on read function since reading same ddta from file produces different results.
  */
  A.pack();
  // Get Compressed Tensor
  Storage s = A.getStorage();

  // Extract nnz
  Array values = s.getValues();
  int size = values.getSize();

  cout<<s<<endl;

  double *vals = (double *)values.getData();

  int i,j;

  cout<<"nnz: "<<size<<endl;

  cout<<"Values: ";
  for(i=0;i<size;i++) cout<<vals[i]<<" ";
  cout<<endl;

  // Extract Index arrays
  Index ind = s.getIndex();

  int size2 = ind.getSize();
  write("test2.mtx", A);
}
