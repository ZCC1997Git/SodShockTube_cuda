#include "solver.cuh"
#include <math.h>
#include <stdio.h>
using namespace std;
/*constan value stored in GPU*/
inline __constant__ double gpu_R;
inline __constant__ double gpu_gamma;
inline __constant__ double gpu_dx;
inline __constant__ double gpu_dt;

__device__ double F21[NUM];
__device__ double F22[NUM];
__device__ double F23[NUM];
/*initializa the fild*/
__global__ void initia(gpuStatus gpu_data, gpuU gpu_U) {
  /*one dimension grid and block*/
  unsigned int id = blockDim.x * blockIdx.x + threadIdx.x;

  gpu_data.x[id] = (id + 1) * gpu_dx;
  if (gpu_data.x[id] < 0.5) {
    gpu_data.u[id] = 0;
    gpu_data.rou[id] = 1;
    gpu_data.p[id] = 1;
  } else {
    gpu_data.u[id] = 0;
    gpu_data.rou[id] = 0.125;
    gpu_data.p[id] = 0.1;
  }
  gpu_data.E[id] = gpu_data.p[id] / (gpu_gamma - 1) / gpu_data.rou[id] +
                   0.5 * pow(gpu_data.u[id], 2);
  gpu_data.c[id] = pow(gpu_gamma * gpu_data.p[id] / gpu_data.rou[id], 0.5);
  gpu_data.T[id] = gpu_data.p[id] / gpu_data.rou[id] * gpu_gamma *
                   pow(gpu_data.u[id] / gpu_data.c[id], 2);
  gpu_U.u1[id] = gpu_data.rou[id];
  gpu_U.u2[id] = gpu_data.rou[id] * gpu_data.u[id];
  gpu_U.u3[id] = gpu_data.rou[id] * gpu_data.E[id];
}
/*rou average process*/
__device__ __forceinline__ void rouAverage(double *L, double *R, double &uu,
                                           double &hh, double &cc) {
  double roul = L[rou_r], rour = R[rou_r];
  double ul = L[rou_u] / roul, ur = R[rou_u] / rour;
  double el = L[rou_e] / roul, er = R[rou_e] / rour;
  double hl, hr, pl, pr;
  pl = (el * roul - 0.5 * roul * pow(ul, 2)) * (gpu_gamma - 1);
  pr = (er * rour - 0.5 * rour * pow(ur, 2)) * (gpu_gamma - 1);
  hl = el + pl / roul;
  hr = er + pr / rour;
  uu = (pow(roul, 0.5) * ul + pow(rour, 0.5) * ur) /
       (pow(roul, 0.5) + pow(rour, 0.5));
  hh = (pow(roul, 0.5) * hl + pow(rour, 0.5) * hr) /
       (pow(roul, 0.5) + pow(rour, 0.5));
  cc = (gpu_gamma - 1) * (hh - 0.5 * pow(uu, 2));
  cc = pow(cc, 0.5);
}
/*A*/
__device__ __forceinline__ void A_cal(double uu, double cc, double *result) {
  for (int i = 0; i <= 8; i++)
    result[i] = 0;
  result[0] = abs(uu - cc);
  result[4] = abs(uu);
  result[8] = abs(uu + cc);
}
/*R*/
__device__ __forceinline__ void R_cal(double uu, double cc, double H,
                                      double *result) {
  result[0] = result[1] = result[2] = 1;
  result[3] = uu - cc;
  result[4] = uu;
  result[5] = uu + cc;
  result[6] = H - uu * cc;
  result[7] = pow(uu, 2) / 2;
  result[8] = H + uu * cc;
}
/*L*/
__device__ __forceinline__ void L_cal(double uu, double cc, double *result) {
  result[0] = 0.5 * (0.5 * (gpu_gamma - 1) * pow(uu, 2) / pow(cc, 2) + uu / cc);
  result[1] = -0.5 * ((gpu_gamma - 1) * uu / pow(cc, 2) + 1 / cc);
  result[2] = 0.5 * (gpu_gamma - 1) / pow(cc, 2);
  result[3] = 1 - 0.5 * (gpu_gamma - 1) * pow(uu, 2) / pow(cc, 2);
  result[4] = (gpu_gamma - 1) * uu / pow(cc, 2);
  result[5] = -(gpu_gamma - 1) / pow(cc, 2);
  result[6] = 0.5 * (0.5 * (gpu_gamma - 1) * pow(uu, 2) / pow(cc, 2) - uu / cc);
  result[7] = -0.5 * ((gpu_gamma - 1) * uu / pow(cc, 2) - 1 / cc);
  result[8] = 0.5 * (gpu_gamma - 1) / pow(cc, 2);
}
/*calculate F using U*/
__device__ __forceinline__ void U2F(double *U, double *F) {
  double density, velocity, energy, pressure;
  density = U[rou_r];
  velocity = U[rou_u] / density;
  energy = U[rou_e] / density;
  pressure =
      (energy * density - 0.5 * density * pow(velocity, 2)) * (gpu_gamma - 1);
  F[0] = density * velocity;
  F[1] = density * pow(velocity, 2) + pressure;
  F[2] = (density * energy + pressure) * velocity;
}
/*calculate the mulplication of Matrix and vector*/
__device__ __forceinline__ void MMV(double *M, double *V) {
  double tmp[3];
  tmp[0] = V[0];
  tmp[1] = V[1];
  tmp[2] = V[2];
  for (int i = 0; i < 3; i++) {
    V[i] = 0;
    for (int j = i * 3; j < (i + 1) * 3; j++)
      V[i] += M[j] * tmp[j - i * 3];
  }
}
/*rou formate*/
__device__ __forceinline__ void Rou(double *R, gpuU U, int size, int id) {
  double Ul[3], Ur[3], F1[3], F2[3];
  double uu, hh, cc; // average value;
  double mA[9], mR[9], mL[9];
  if (id < size - 1) {
    Ur[0] = U.u1[id + 1];
    Ur[1] = U.u2[id + 1];
    Ur[2] = U.u3[id + 1];
  } else {
    Ur[0] = U.u1[id];
    Ur[1] = U.u2[id];
    Ur[2] = U.u3[id];
  }
  Ul[0] = U.u1[id];
  Ul[1] = U.u2[id];
  Ul[2] = U.u3[id];
  rouAverage(Ul, Ur, uu, hh, cc);
  A_cal(uu, cc, mA);
  R_cal(uu, cc, hh, mR);
  L_cal(uu, cc, mL);

  U2F(Ul, F1);
  U2F(Ur, F2);
  F21[id] = 0.5 * (F1[0] + F2[0]);
  F22[id] = 0.5 * (F1[1] + F2[1]);
  F23[id] = 0.5 * (F1[2] + F2[2]);

  Ur[0] = Ur[0] - Ul[0];
  Ur[1] = Ur[1] - Ul[1];
  Ur[2] = Ur[2] - Ul[2];
  MMV(mL, Ur);
  MMV(mA, Ur);
  MMV(mR, Ur);
  F21[id] -= 0.5 * Ur[0];
  F22[id] -= 0.5 * Ur[1];
  F23[id] -= 0.5 * Ur[2];
  /*waiting for all threads to update shared memory*/
  __syncthreads();
  if (id > 0) {
    R[0] = (-1.0) * (F21[id] - F21[id - 1]) / gpu_dx;
    R[1] = (-1.0) * (F22[id] - F22[id - 1]) / gpu_dx;
    R[2] = (-1.0) * (F23[id] - F23[id - 1]) / gpu_dx;
  } else
    R[0] = R[1] = R[2] = 0;
}
/*calculate the filed*/
__global__ void GPU_RungeKutta(gpuU gpu_U, int size) {
  double U0[3], R[3];
  /*one dimension grid and block*/
  int id = blockDim.x * blockIdx.x + threadIdx.x;
  /*copy gpu_U*/
  U0[0] = gpu_U.u1[id];
  U0[1] = gpu_U.u2[id];
  U0[2] = gpu_U.u3[id];
  /*rou*/
  Rou(R, gpu_U, size, id);
  /*update U using R*/
  gpu_U.u1[id] += gpu_dt * R[0];
  gpu_U.u2[id] += gpu_dt * R[1];
  gpu_U.u3[id] += gpu_dt * R[2];
  __syncthreads();
  /*rou*/
  Rou(R, gpu_U, size, id);
  /*update U using R*/
  double ratio = 1.0 / 4;
  gpu_U.u1[id] =
      (1 - ratio) * U0[0] + ratio * gpu_U.u1[id] + (ratio * gpu_dt) * R[0];
  gpu_U.u2[id] =
      (1 - ratio) * U0[1] + ratio * gpu_U.u2[id] + (ratio * gpu_dt) * R[1];
  gpu_U.u3[id] =
      (1 - ratio) * U0[2] + ratio * gpu_U.u3[id] + (ratio * gpu_dt) * R[2];
  __syncthreads();
  /*rou*/
  Rou(R, gpu_U, size, id);
  /*update U using R*/
  ratio = 2.0 / 3;
  gpu_U.u1[id] =
      (1 - ratio) * U0[0] + ratio * gpu_U.u1[id] + (ratio * gpu_dt) * R[0];
  gpu_U.u2[id] =
      (1 - ratio) * U0[1] + ratio * gpu_U.u2[id] + (ratio * gpu_dt) * R[1];
  gpu_U.u3[id] =
      (1 - ratio) * U0[2] + ratio * gpu_U.u3[id] + (ratio * gpu_dt) * R[2];
}
