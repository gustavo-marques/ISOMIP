import numpy as np
import matplotlib.pyplot as plt

#gamma_com = np.arange(0.01,0.21,0.01)
#melt_com = [3.20, 7.15, 10.96, 14.24, 17.09, 19.52, 21.57, 23.34, 24.74, 26.03, 27.33, 28.11, 28.94, 29.74, 30.30, 30.84, 31.36, 31.85, 32.32, 32.57]

#gamma_tpy = [0.01, 0.05, 0.1, 0.15, 0.2, 0.5, 1., 1.5, 2., 2.5, 3.]
#melt_tpy = [9.56, 14.4, 15.45, 16.24, 19.96, 22.76, 30.05, 34.47, 37.71, 41.68, 45.17]

gamma_h10 = [0.01, 0.02, 0.03, 0.04, 0.05, 0.075,0.1,0.12,0.14,0.15,0.16,0.18,0.2]
melt_h10 = [3.19,7.13,10.94,14.21,16.98,22.48,26.20,28.24,29.90,30.45,30.98,32.,32.75]

#coeffs = np.polyfit(gamma_com,melt_com,3) # third-order poly
coeffs = np.polyfit(gamma_h10,melt_h10,4) # third-order poly
new_gamma = np.arange(0.01,0.21,0.0001)
new_melt = coeffs[0]*new_gamma**4 + coeffs[1]*new_gamma**3 + coeffs[2]*new_gamma**2 + coeffs[3]*new_gamma + coeffs[4]
# find values that yields melt=30
tmp = np.nonzero(new_melt<=30.)[-1][-1]
print('Gamma is: '+str(new_gamma[tmp]))

plt.figure()
#plt.plot(gamma_com,melt_com,'ko')
plt.plot(gamma_h10,melt_h10,'ko')
plt.plot(new_gamma,new_melt,'k-',lw=1.0)
plt.plot(new_gamma[tmp],new_melt[tmp],'b*',markersize=12)
plt.xlabel(r'$\Gamma_T$')
plt.ylabel('mean melt rate (m a$^{-1}$)')
plt.xlim(0,0.2001)
plt.savefig('melt_gamma_COM.png')

#plt.figure()
#plt.plot(gamma_com,melt_com,'k-o',label='COM',lw=1)
#plt.plot(gamma_tpy,melt_tpy,'r-o',label='TPY',lw=1)
#plt.legend(loc='upper left', shadow=True)
#plt.xlabel(r'$\Gamma_T$')
#plt.ylabel('mean melt rate (m a$^{-1}$)')
#plt.xlim(0,3.01)
#plt.savefig('melt_gamma_COMvsTPY.png')

