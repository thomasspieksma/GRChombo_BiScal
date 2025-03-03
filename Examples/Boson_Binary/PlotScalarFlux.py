import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("/lustre/astro/spieksma/Sims_Scalar/Boson_Cloud_Binary3/data/ScalarFlux_mode_11.dat")

fig,(ax1) = plt.subplots(1,1,sharex=True,figsize = (9,8))
ax1.plot(data[:,0], data[:,1], label = r'$x = 10M$')

ax1.set_ylabel(r"$\Phi_{11}$", fontsize = 12)


ax1.legend(fontsize = 12, ncol = 2)

plt.savefig("phiFlux_11.pdf")
plt.show()
