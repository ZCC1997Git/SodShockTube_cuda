#include <iostream>
#include "solver.h"
#include <chrono>
using namespace std;

int main()
{
	{
		std::ifstream cpuinfo("/proc/cpuinfo");
		std::string line;
		std::string cpuModel;

		while (std::getline(cpuinfo, line))
		{
			if (line.substr(0, 10) == "model name")
			{
				cpuModel = line.substr(line.find(":") + 2);
				break;
			}
		}
		std::cout << "running on " << cpuModel << std::endl;
	}
	Euler test(2048, 0.2);
	auto start_cpu = std::chrono::high_resolution_clock::now();
	test.solver();
	auto stop_cpu = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop_cpu - start_cpu);
	auto CPU_time = (double)duration.count() / 1000;
	std::cout << "CPU runing time: " << CPU_time << " ms" << std::endl;
	std::cout << std::endl;
	test.output("./output/cpu");
	return 0;
}
