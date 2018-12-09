#ifndef SIMU_ARRAY_H
#define SIMU_ARRAY_H

struct Array1
{
  double* ptr;
  int n;

  Array1() : ptr(nullptr), n(0) {} // constructor

  void allocate(int an)
  {
    free();
    n = an;
    ptr = new double[an];
  }

  void free()
  {
    delete[] ptr;
    ptr = nullptr;
    n = 0;
  }

  void fill(double value)
  {
    for(int i=0; i<n; ++i){
      ptr[i] = value;
    }
  }

  // access element in array by providing reference to it
  double& operator()(int i) { return *(ptr + i); }
};

#endif
