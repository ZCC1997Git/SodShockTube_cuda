#ifndef solver_h
#define solver_h
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <fstream>
#include <iostream>
using namespace std;
/*the number of nodes*/
inline constexpr int NUM = 2048;
/*the information of GPU*/
inline int SM;
inline int Block;

/*define some variables which respond to the physical status */
struct status {
  float x;   /*position*/
  float rou; /*density;*/
  float u;   /*velocity;*/
  float p;   /*pressure;*/
  float E;   /*total energy;*/
  float T;   /*tmpeture*/
  float c;   /*the speed of sound*/
};
struct gpuStatus {
  float *x = nullptr;
  float *rou = nullptr;
  float *u = nullptr;
  float *p = nullptr;
  float *E = nullptr;
  float *T = nullptr;
  float *c = nullptr;
};
struct gpuU {
  float *u1 = nullptr;
  float *u2 = nullptr;
  float *u3 = nullptr;
};
struct Vector {
  float rou_r;
  float rou_u;
  float rou_e;
};
/*help to specific the physical meaning of U and F for different index*/
typedef enum {
  rou_r, /*density*/
  rou_u, /*density x velocity*/
  rou_e, /*density x energy*/
} U;
//------------------------------------------------------//
class Euler {
private:
  float R0 = 287.0;
  float gamma = 1.4;
  int size;
  float dt, dx, t;
  status *data = nullptr;
  Vector *U = nullptr;
  /*gpu*/
  gpuStatus gpu_data;
  gpuU gpu_U;
  cudaDeviceProp device_prop;
  cudaEvent_t kernel_start, kernel_end;
  float runtime = 0.0;
  /*update data*/
  void refresh();

public:
  Euler(int size, float t);
  ~Euler();
  void Runge_Kutta();
  void solver();
  void memCpy();
  void output(string name) const;
};
//------------------------------------------------------//
#endif
