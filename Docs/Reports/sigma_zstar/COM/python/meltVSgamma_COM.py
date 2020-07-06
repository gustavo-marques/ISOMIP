import numpy as np
import matplotlib.pyplot as plt

gamma_h10 = [0.01, 0.02, 0.03, 0.04, 0.05, 0.07,0.09,0.1, 0.12]
melt_h10 = [4.6, 12.27, 20.0, 27.3, 33.19, 41.3, 46.0, 47.4, 49.4]

#coeffs = np.polyfit(gamma_com,melt_com,3) # third-order poly
coeffs = np.polyfit(gamma_h10,melt_h10,4) # third-order poly
new_gamma = np.arange(0.01,0.12,0.0001)
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
plt.show()
#plt.figure()
#plt.plot(gamma_com,melt_com,'k-o',label='COM',lw=1)
#plt.plot(gamma_tpy,melt_tpy,'r-o',label='TPY',lw=1)
#plt.legend(loc='upper left', shadow=True)
#plt.xlabel(r'$\Gamma_T$')
#plt.ylabel('mean melt rate (m a$^{-1}$)')
#plt.xlim(0,3.01)
#plt.savefig('melt_gamma_COMvsTPY.png')

