from scipy.optimize import leastsq
import numpy as np


print "\nThis script works on the assumption that data is in format created by running dataAnalysis_volScan.sh"
#vols = np.array([13.71, 14.82, 16.0, 17.23, 18.52])
#energies = np.array([-56.29, -56.41, -56.46, -56.463, -56.41])
print "\nGrabbing data ..."
with open('Energy.dat', 'rb') as f:
    content = f.readlines()
vols = np.array([float(x.split()[0]) for x in content])
energies = np.array([float(x.split()[3]) for x in content]) 
scalefactors = np.array([float(x.strip('\n').split()[-1]) for x in content])

def Murnaghan(parameters, vol):
    'From Phys. Rev. B 28, 5480 (1983)'
    E0, B0, BP, V0 = parameters

    E = E0 + B0 * vol / BP * (((V0 / vol)**BP) / (BP - 1) + 1) - V0 * B0 / (BP - 1.0)

    return E

def objective(pars, y, x):
    #we will minimize this function
    err =  y - Murnaghan(pars, x)
    return err

print "\n!!!!!!!!!!!!!!!!!!!!!!!!!!"
print "Fitting to volume now"
print "!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
#x0 = [ -56.0, 0.54, 2.0, 16.5] #initial guess of parameters
# The initial guess for volume needs to be far larger when working with
# the volumes in angstroms cubed. Compare with x0 below.
x0 = [-56, 0.54, 2.0, 1000]

plsq = leastsq(objective, x0, args=(energies, vols), maxfev=10000)

print 'Fitted parameters = {0}'.format(plsq[0])
print "\nEquilibrium volume = {0:.2f} $\AA^3$".format(plsq[0][3])

import matplotlib.pyplot as plt
import matplotlib
font = {'size'   : 16}

matplotlib.rc('font', **font)

plt.figure(1)
plt.plot(vols,energies, 'ro')

#plot the fitted curve on top
x = np.linspace(min(vols), max(vols), 50)
y = Murnaghan(plsq[0], x)
plt.plot(x, y, 'k-')
plt.annotate(r'V$_0$ = {0:.2f} $\AA^3$'.format(plsq[0][3]), xy=(0.6, 0.2),  xycoords='axes fraction')
plt.xlabel(r'Volume ($\AA^3$)')
plt.ylabel('Energy (eV)')
plt.savefig('equation_of_state.pdf')
plt.close(1)

print "\n!!!!!!!!!!!!!!!!!!!!!!!!!!"
print "Fitting to scalefactor now"
print "!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
x0 = [ -56.0, 0.54, 2.0, 1.0] #initial guess of parameters

plsq = leastsq(objective, x0, args=(energies, scalefactors), maxfev=10000)

print 'Fitted parameters = {0}'.format(plsq[0])
print "\nEquilibrium scalefactor = {0:.3f}".format(plsq[0][3])
print ""
plt.figure(2)
plt.plot(scalefactors,energies, 'ro')

#plot the fitted curve on top
x = np.linspace(min(scalefactors), max(scalefactors), 50) 
y = Murnaghan(plsq[0], x)
plt.plot(x, y, 'k-')
plt.annotate(r'Scale$_0$ = {0:.3f}'.format(plsq[0][3]), xy=(0.6, 0.2),  xycoords='axes fraction')
plt.xlabel(r'Scale')
plt.ylabel('Energy (eV)')
plt.savefig('equation_of_state_scale.pdf')
plt.close(2)
