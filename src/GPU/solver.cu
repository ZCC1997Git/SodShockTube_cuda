#include "solver.cuh"
#include <math.h>
#include <stdio.h>
using namespace std;
void displayProgressBarGPU(int progress, int total, int barWidth = 60);
__global__ void initia(gpuStatus gpu_data, gpuU gpu_U);
__global__ void GPU_RungeKutta(gpuU gpu_U, int size);

/*constan value stored in GPU*/
extern __constant__ float gpu_R;
extern __constant__ float gpu_gamma;
extern __constant__ float gpu_dx;
extern __constant__ float gpu_dt;

/*update data*/
void Euler::refresh() {
  for (int i = 0; i <= size - 1; i++) {
    data[i].rou = U[i].rou_r;
    data[i].u = U[i].rou_u / U[i].rou_r;
    data[i].E = U[i].rou_e / U[i].rou_r;
    data[i].p =
        (data[i].E * data[i].rou - 0.5 * data[i].rou * pow(data[i].u, 2)) *
        (gamma - 1);
    data[i].c = pow(gamma * data[i].p / data[i].rou, 0.5);
    data[i].T = data[i].p / data[i].rou * gamma * pow(data[i].u / data[i].c, 2);
  }
}

/*construction function*/
Euler::Euler(int size, float t) {
  this->size = size;
  this->t = t;
  dx = 1.0 / size;
  dt = 0.0001;
  data = new status[size];
  U = new Vector[size];
  /*gpu_data*/
  cudaMalloc((void **)&(gpu_data.x), size * sizeof(float));
  cudaMalloc((void **)&(gpu_data.rou), size * sizeof(float));
  cudaMalloc((void **)&(gpu_data.u), size * sizeof(float));
  cudaMalloc((void **)&(gpu_data.p), size * sizeof(float));
  cudaMalloc((void **)&(gpu_data.E), size * sizeof(float));
  cudaMalloc((void **)&(gpu_data.T), size * sizeof(float));
  cudaMalloc((void **)&(gpu_data.c), size * sizeof(float));
  /*gpu_u*/
  cudaMalloc((void **)&(gpu_U.u1), size * sizeof(float));
  cudaMalloc((void **)&(gpu_U.u2), size * sizeof(float));
  cudaMalloc((void **)&(gpu_U.u3), size * sizeof(float));
  /*set constant memory*/
  cudaMemcpyToSymbol(gpu_R, &R0, sizeof(float));
  cudaMemcpyToSymbol(gpu_gamma, &gamma, sizeof(float));
  cudaMemcpyToSymbol(gpu_dx, &dx, sizeof(float));
  cudaMemcpyToSymbol(gpu_dt, &dt, sizeof(float));
  /*event created*/
  cudaEventCreate(&kernel_start);
  cudaEventCreate(&kernel_end);
  cudaEventCreateWithFlags(&kernel_end, cudaEventBlockingSync);
  /*get property*/
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  if (deviceCount == 0) {
    cerr << "No CUDA-capable devices found." << endl;
    exit(0);
  } else
    cout << deviceCount << " GPUs has been found!" << endl;

  cudaGetDeviceProperties(&device_prop, 0);
  cout << "running on " << device_prop.name << endl;
  SM = device_prop.multiProcessorCount;
  cout << "Block number in  grid : Thread in each Block (" << SM << " : "
       << size / SM << ")" << endl;
  /*call gpu_kernel*/
  initia<<<SM, size / SM>>>(gpu_data, gpu_U);
}

/*the process of runge-kutta*/
void Euler::Runge_Kutta() {
  float time;
  cudaEventRecord(kernel_start, 0);
  GPU_RungeKutta<<<SM, size / SM>>>(gpu_U, size);
  cudaEventRecord(kernel_end, 0);
  cudaEventSynchronize(kernel_end);
  cudaEventElapsedTime(&time, kernel_start, kernel_end);
  runtime += time;
}
/*solve the equation*/
void Euler::solver() {
  float tt = 0;
  int count = 0;
  while (tt < t) {
    count++;
    displayProgressBarGPU(count, t / dt);
    Runge_Kutta();
    tt = tt + dt;
  }
  cout << endl;
}

/*delete the memory*/
Euler::~Euler() {
  if (data)
    delete[] data;
  if (U)
    delete[] U;
  cudaEventDestroy(kernel_start);
  cudaEventDestroy(kernel_end);
  cudaFree(gpu_data.x);
  cudaFree(gpu_data.rou);
  cudaFree(gpu_data.u);
  cudaFree(gpu_data.p);
  cudaFree(gpu_data.E);
  cudaFree(gpu_data.T);
  cudaFree(gpu_data.c);
  cudaFree(gpu_U.u1);
  cudaFree(gpu_U.u2);
  cudaFree(gpu_U.u3);
}

/*copy data back to cpu*/
void Euler::memCpy() {
  float *tmp = (float *)malloc(size * sizeof(float));
  cudaMemcpy(tmp, gpu_data.x, size * sizeof(float), cudaMemcpyDeviceToHost);
  for (int i = 0; i < size; i++)
    data[i].x = tmp[i];
  cudaMemcpy(tmp, gpu_U.u1, size * sizeof(float), cudaMemcpyDeviceToHost);
  for (int i = 0; i < size; i++)
    U[i].rou_r = tmp[i];
  cudaMemcpy(tmp, gpu_U.u2, size * sizeof(float), cudaMemcpyDeviceToHost);
  for (int i = 0; i < size; i++)
    U[i].rou_u = tmp[i];
  cudaMemcpy(tmp, gpu_U.u3, size * sizeof(float), cudaMemcpyDeviceToHost);
  for (int i = 0; i < size; i++)
    U[i].rou_e = tmp[i];
  free(tmp);
  refresh();
}

/*output result*/
void Euler::output(string name) const {
  cout << "runtime on gpu is:" << runtime << "ms" << endl;
  name.append(".dat");
  ofstream out(name);
  for (int i = 0; i <= size - 1; i++)
    // out << data[i].x << "\t" << data[i].p << endl;
    out << data[i].x << '\t' << data[i].u << '\t' << data[i].p << '\t'
        << data[i].rou << '\t' << data[i].T << endl;
  out.close();
}
