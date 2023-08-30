#include <iostream>
using namespace std;

void displayProgressBarGPU(int progress, int total, int barWidth = 60) {
  float fraction = static_cast<float>(progress) / total;
  int width = static_cast<int>(barWidth * fraction);

  std::cout << "[";
  for (int i = 0; i < barWidth; ++i) {
    if (i < width)
      std::cout << "=";
    else
      std::cout << " ";
  }
  std::cout << "] " << int(fraction * 100.0) << "%"
            << "\r";
  std::cout.flush();
}
