import math
from math import pi as pi
from math import sqrt as sqrt
import array
import numpy as np
from scipy import stats
from scipy import special
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import binom
from scipy.stats import poisson
from tabulate import tabulate


prev = np.array([780, 500, 200, 20]) #previous data
prvang = np.array([0, 7.5, 15, 22.5]) #previous data angles

m3, b3 = np.polyfit(prev, prvang, 1)

plt.plot(prvang, prev, 'o', c='orange', label='Ground State (Previous)')
#plt.plot(prev*m3 + b3, prev, c='orange')

ang = np.array([0, 3.093, 6.185, 9.277, 12.368, 15.458, 18.547]) #present angles


#online
online = np.array([779, 813, 572, 491, 351, 190, 124])

online_err = np.array([20, 19, 17, 20, 13, 10, 14])

m2, b2 = np.polyfit(online, ang, 1)

#plt.errorbar(ang, online, c='pink', yerr=online_err, fmt='none')
#plt.plot(ang, online, 'o', c='pink', label='Online')
#plt.plot(online*m2 + b2, online, c='pink')

#ground
ground = np.array([792, 629, 604, 474, 380, 187, 98]) #present data

ground_err = np.array([22, 19, 20, 21, 17, 12, 16])

m, b = np.polyfit(ground, ang, 1)

plt.errorbar(ang, ground, c='blue', yerr=ground_err,fmt='none')
plt.plot(ang, ground, 'o', c='blue', label='Ground State (Present)')
#plt.plot(ground*m + b, ground, c='blue')

#plt.title('Cross Section vs. Scattering Angle; 11 MeV', fontsize=20)
#plt.xlabel('Center of Mass Angle (Degrees)', fontsize=15)
#plt.ylabel('Cross Section (\u03BC/sr)', fontsize=15)
#plt.grid()
#legend = plt.legend(loc='upper right', fontsize=13)


#plt.show()

####################################################################################


#2+ first excited state (1.5 MeV)
tfes = np.array([63, 80, 78, 11, 37, 53, 18])

tfes_err = np.array([20, 16, 18, 19, 18, 9, 14])

#plt.errorbar(ang, tfes, c='red', yerr=tfes_err,fmt='none')
#plt.plot(ang, tfes, 'o', c='red', label='2+ First Excited State')

m4, b4 = np.polyfit(tfes, ang, 1)
#plt.plot(tfes*m4 + b4, tfes, c='red')


#0+ first excited state (1.8 MeV)
zfes = np.array([35, 10, 19, 57, 25, 0, 5])

zfes_err = np.array([22, 19, 19, 21, 21, 0, 14])

#plt.errorbar(ang, zfes, c='orange', yerr=zfes_err,fmt='none')
#plt.plot(ang, zfes, 'o', c='orange', label='0+ First Excited State')

#m5, b5 = np.polyfit(zfes, ang, 1)
#plt.plot(zfes*m5 + b5, zfes, c='orange')


#2+ second excited state (2.4 MeV)
tses = np.array([13, 21, 30, 0, 17, 28, 0])

tses_err = np.array([11, 8, 11, 0, 10, 7, 0])

#plt.errorbar(ang, tses, c='yellow', yerr=tses_err,fmt='none')
#plt.plot(ang, tses, 'o', c='yellow', label='2+ Second Excited State')

m6, b6 = np.polyfit(tses, ang, 1)
#plt.plot(tses*m6 + b6, tses, c='yellow')


#0+ second excited state (3.3 MeV)
zses = np.array([69, 43, 55, 43, 50, 48, 40])

zses_err = np.array([12, 10, 10, 12, 11, 10, 14])

#plt.errorbar(ang, zses, c='violet', yerr=zses_err,fmt='none')
#plt.plot(ang, zses, 'o', c='violet', label='0+/2+ Second Excited State')

m7, b7 = np.polyfit(zses, ang, 1)
#plt.plot(zses*m7 + b7, zses, c='violet')


#plt.xlim(-2,23)
#plt.ylim(-15, 250)


plt.title('Cross Section vs. Scattering Angle; 11 MeV', fontsize=20)
plt.xlabel('Center of Mass Angle (Degrees)', fontsize=15)
plt.ylabel('Cross Section (\u03BC/sr)', fontsize=15)
plt.grid()
legend = plt.legend(loc='upper right', fontsize=10)

plt.show()

################################################################################










