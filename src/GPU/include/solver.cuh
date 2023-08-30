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
  double x;   /*position*/
  double rou; /*density;*/
  double u;   /*velocity;*/
  double p;   /*pressure;*/
  double E;   /*total energy;*/
  double T;   /*tmpeture*/
  double c;   /*the speed of sound*/
};
struct gpuStatus {
  double *x = nullptr;
  double *rou = nullptr;
  double *u = nullptr;
  double *p = nullptr;
  double *E = nullptr;
  double *T = nullptr;
  double *c = nullptr;
};
struct gpuU {
  double *u1 = nullptr;
  double *u2 = nullptr;
  double *u3 = nullptr;
};
struct Vector {
  double rou_r;
  double rou_u;
  double rou_e;
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
  double R0 = 287.0;
  double gamma = 1.4;
  int size;
  double dt, dx, t;
  status *data = nullptr;
  Vector *U = nullptr;
  /*gpu*/
  gpuStatus gpu_data;
  gpuU gpu_U;
  cudaDeviceProp device_prop;
  cudaEvent_t kernel_start, kernel_end;
  double runtime = 0.0;
  /*update data*/
  void refresh();

public:
  Euler(int size, double t);
  ~Euler();
  void Runge_Kutta();
  void solver();
  void memCpy();
  void output(string name) const;
};
//------------------------------------------------------//
#endif
