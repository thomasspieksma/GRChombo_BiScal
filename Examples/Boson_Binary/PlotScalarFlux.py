import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("/lustre/astro/spieksma/Sims_Scalar/Boson_Cloud_LR/data/ScalarFlux_mode_11.dat")
data2 = np.genfromtxt("/lustre/astro/spieksma/Sims_Scalar/Boson_Cloud_LR/data/ScalarFlux_mode_22.dat")
data3 = np.genfromtxt("/lustre/astro/spieksma/Sims_Scalar/Boson_Cloud_LR/data/ScalarFlux_mode_20.dat")

fig,(ax1,ax2) = plt.subplots(2,1,sharex=True,figsize = (9,8))
ax1.plot(data[:,0], data[:,1], label = r'$(1,1)$')
ax1.plot(data2[:,0], data2[:,1], label = r'$(2,2)$')
ax1.plot(data3[:,0], data3[:,1], label = r'$(2,0)$')

ax1.set_ylabel(r"$\mathrm{Re}[\Phi_{\ell m}]$", fontsize = 18)
ax2.set_xlabel(r"$t/M$", fontsize = 18)
ax2.set_ylabel(r"$\mathrm{Re}[\Phi_{\ell m}]$", fontsize = 18)

ax1.set_title(r"$r = 10M$", fontsize = 18)

ax2.plot(data[:,0], abs(data[:,1]))
ax2.plot(data2[:,0], abs(data2[:,1]))
ax2.plot(data3[:,0], abs(data3[:,1]))
ax2.set_yscale('log')
ax2.set_ylim(1e-10, 1e-1)
ax1.legend(fontsize = 15)
ax1.set_xlim(0,100),ax2.set_xlim(0,100)

plt.savefig("phiFlux_11.pdf")
plt.show()
