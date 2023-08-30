import numpy as np
import matplotlib.pyplot as plt

print("\nThe results are plotting through Python")
# x= np.loadtxt('x.dat')
cpu= np.loadtxt('./output/cpu.dat')
gpu= np.loadtxt('./output/gpu.dat')
cpu=cpu.T 
gpu=gpu.T 

plt.plot(cpu[0],cpu[1],color='k',label='CPU')
plt.plot(gpu[0],gpu[1],color='r',label='GPU')
plt.xlabel("X")
plt.ylabel("P")
plt.title("SodShockWave code simluation result")
plt.legend()

plt.savefig("SimulationResult.png")

print("Plotting is over!")
