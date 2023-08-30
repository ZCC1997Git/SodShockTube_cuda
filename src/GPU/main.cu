#include "solver.cuh"
using namespace std;

int main() {
  Euler test(NUM, 0.2);
  test.solver();
  test.memCpy();
  test.output("./output/gpu");
  return 0;
}
