import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("/lustre/astro/spieksma/Sims_Scalar/Boson_Cloud_Binary/data/data_out.dat")
# data2 = np.genfromtxt("/lustre/astro/spieksma/Sims_Scalar/Boson_Cloud_Binary/data/phi_lineout.dat")
# data_lineout = np.genfromtxt("/lustre/astro/spieksma/Sims_Scalar/Sims_Scalar/data_line_ex_phi/phi_lineout.dat")

fig,(ax1,ax2) = plt.subplots(2,1,sharex=True,figsize = (9,8))
ax1.plot(data[:,0], data[:,-2], label = r'$x = 0M$')
# ax1.plot(data_lineout[:,0], data_lineout[:,2], label = r'$x = 1.6M$')
# ax1.plot(data_lineout[:,0], data_lineout[:,3], label = r'$x = 3.2M$')
# ax1.plot(data_lineout[:,0], data_lineout[:,4], label = r'$x = 4.8M$')
# ax1.plot(data_lineout[:,0], data_lineout[:,5], label = r'$x = 6.4M$')
# ax1.plot(data_lineout[:,0], data_lineout[:,6], label = r'$x = 8M$')
# ax1.plot(data_lineout[:,0], data_lineout[:,7], label = r'$x = 9.6M$')
# ax1.plot(data_lineout[:,0], data_lineout[:,8], label = r'$x = 11.2M$')
# ax1.plot(data_lineout[:,0], data_lineout[:,9], label = r'$x = 12.8M$')
# ax1.plot(data_lineout[:,0], data_lineout[:,10], l4abel = r'$x = 14.4M$')

# ax2.plot(data2[:,0], data2[:,-1], label = r'$x = 14.40M$')
# ax2.plot(data[:,0], data[:,2], label = r'$||\mathcal{M}||$')

# ax2.plot(data_hr2[:,0], data_hr2[:,1], label = r'$||\mathcal{H}||_{2}$')
# ax2.plot(data_hr2[:,0], data_hr2[:,2], label = r'$||\mathcal{M}||_{2}$')

ax1.set_ylabel(r"$\Phi$", fontsize = 12)

ax2.set_xlabel(r"$t/M$", fontsize = 12)

ax1.legend(fontsize = 12, ncol = 2)
ax2.legend(fontsize = 12)

plt.savefig("phi_plot_time.pdf")
plt.show()
