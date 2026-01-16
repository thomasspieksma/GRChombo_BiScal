import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("/lustre/astro/spieksma/Sims_Scalar/Sims_Scalar/data_line_ex_phi/phi_lineout.dat")


fig,ax1 = plt.subplots(1,1,sharex=True,figsize = (9,8))
ax1.plot(data[0][1:], color = 'k')

ax1.set_ylabel(r"$\Phi$", fontsize = 12)
ax1.set_xlabel(r"$x/M$", fontsize = 12)
ax1.set_yscale('log')
#plt.legend(fontsize = 12)
plt.savefig("phi_plot_space.pdf")
plt.show()
