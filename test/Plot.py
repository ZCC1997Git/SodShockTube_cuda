import numpy as np
import matplotlib.pyplot as plt

print("\nThe results are plotting through Python")
# x= np.loadtxt('x.dat')
cpu= np.loadtxt('./output/cpu.dat')
gpu= np.loadtxt('./output/gpu.dat')
cpu=cpu.T 
gpu=gpu.T 

plt.plot(cpu[0],cpu[1],color='k',label='CPU_u')
plt.plot(gpu[0],gpu[1],color='g',label='GPU_u')
plt.plot(cpu[0],cpu[2],color='k',linestyle='--',label='CPU_p')
plt.plot(gpu[0],gpu[2],color='g',linestyle='--',label='GPU_p')
plt.plot(cpu[0],cpu[4],color='k',linestyle=':',label='CPU_T')
plt.plot(gpu[0],gpu[4],color='g',linestyle=':',label='GPU_T')

plt.xlabel("X")
plt.ylabel("U/P/T")
plt.title("SodShockWave code simluation result")
plt.legend()
plt.grid(True, linestyle='--',color='gray',alpha=0.5)

plt.savefig("SimulationResult.png")

print("Plotting is over!")
