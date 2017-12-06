import numpy as np
from StringIO import StringIO
import math
#function to calculate rot partition function
#at temp T (K) with rotational symmetry number
#sigma_rot (see http://link.springer.com/article/10.1007%2Fs00214-007-0328-0 
#- Symmetry numbers and reaction rates, Ramos et al).
#IA,IB,IC should be in units of amu Angstrom^2, they will be converted
#inside the function.
hbar=6.582E-16 #eV s
conversion_I=(10E-20*931.494E6/9E16) #convert IA etc from amu Angstrom^2 to eV s^2
			      #this is 931.494 MeV/c^2 x 10E-20 m^2
			
def Z_rot(T,IA,IB,IC, sigma_rot):
	k_BT=(float(T)/300)*0.025 #eV
	Z_rot=math.sqrt(k_BT*conversion_I*IA/(hbar*hbar) )*math.sqrt(k_BT*conversion_I*IB/(hbar*hbar) )*math.sqrt(k_BT*conversion_I*IC/(hbar*hbar) )
	Z_rot=math.sqrt(math.pi)*Z_rot
	Z_rot=Z_rot/sigma_rot
	return Z_rot

#function to return the configurational/rotational free energy for the molecule as a function of temp
def F_rot(Z_rot,T):
	if (Z_rot==0):
		print "Error z_rot==0"	
		return 0
	else:
		F_rot=-(float(T)/300)*0.025*math.log(Z_rot)
		return F_rot
	 	

#Read in bare bones xyz file (first two lines deleted) and calculate moment of intertia of the molecule.

data_raw=np.genfromtxt('CH3NH3.xyz',dtype=("|S10", float, float, float),skip_header=2)

data_with_mass = np.zeros(len(data_raw),dtype=[('Atom','|S10'),
						 ('x','float'),
						 ('y','float'),	
						 ('z','float'),
						 ('Mass','float')] )
M=0.0 #total mass

for i in range(len(data_raw)):
	data_with_mass[i][0]=data_raw[i][0]
	data_with_mass[i][1]=data_raw[i][1]
	data_with_mass[i][2]=data_raw[i][2]
	data_with_mass[i][3]=data_raw[i][3]
	if data_raw[i][0]=='C':
		data_with_mass[i][4]=12
		M+=12
	if data_raw[i][0]=='H':
		data_with_mass[i][4]=1
		M+=1
	if data_raw[i][0]=='N':
		data_with_mass[i][4]=14
		M+=14

print data_with_mass, M, "\n"

I=np.zeros( (3,3) )

mixi=0.0 #sum of m_i x_i
miyi=0.0
mizi=0.0

mixiyi = 0.0
mixizi = 0.0 
miyizi = 0.0

for i in range(len(data_with_mass)):
	I[0][0] += data_with_mass[i][4]*(data_with_mass[i][2]*data_with_mass[i][2] + 	
			data_with_mass[i][3]*data_with_mass[i][3]) 
	
	I[1][1] += data_with_mass[i][4]*(data_with_mass[i][1]*data_with_mass[i][1] +	
			data_with_mass[i][3]*data_with_mass[i][3]) 
	
	I[2][2] += data_with_mass[i][4]*(data_with_mass[i][1]*data_with_mass[i][1] +	
			data_with_mass[i][2]*data_with_mass[i][2]) 

	mixi += data_with_mass[i][4]*data_with_mass[i][1]
	mixiyi += data_with_mass[i][4]*data_with_mass[i][1]*data_with_mass[i][2]
 
	miyi += data_with_mass[i][4]*data_with_mass[i][2]
	mixizi += data_with_mass[i][4]*data_with_mass[i][1]*data_with_mass[i][3]
	
	mizi += data_with_mass[i][4]*data_with_mass[i][3]
	miyizi += data_with_mass[i][4]*data_with_mass[i][2]*data_with_mass[i][3]


I[0][0] = I[0][0] - miyi*miyi/M - mizi*mizi/M
I[1][1] = I[1][1] - mixi*mixi/M - mizi*mizi/M
I[2][2] = I[2][2] - mixi*mixi/M - miyi*miyi/M

I[1][0] = -mixiyi + mixi*miyi/M
I[0][1] = I[1][0]

I[2][0] = -mixizi + mixi*mizi/M
I[0][2] = I[2][0]

I[1][2] = -miyizi + miyi*mizi/M
I[2][1] = I[1][2]

print "\n", I, "\n"
print "Principal moments of inertia = ", np.linalg.eig(I)[0][0], ",", np.linalg.eig(I)[0][1], ",",  np.linalg.eig(I)[0][2], "amu Angstrom^2"

print "Product I_A*I_B*I_C = ", np.linalg.eig(I)[0][0]*np.linalg.eig(I)[0][1]*np.linalg.eig(I)[0][2], "amu^3 Angstrom^6\n"

Tmin=50
Tmax=1000
dT=10

F=np.zeros( ( (Tmax-Tmin)/dT,2) )
Z=np.zeros( ( (Tmax-Tmin)/dT,2) )
print F

for i in range(0, (Tmax-Tmin)/dT):
	T=(i+1)*dT
	Z[i][0]=T
	Z[i][1]=Z_rot(T,np.linalg.eig(I)[0][0], np.linalg.eig(I)[0][1], np.linalg.eig(I)[0][2], 3)
	F[i][0]=T
	F[i][1]=F_rot(Z[i][1], T)

print "F=\n", F, "\n"
print "Z=\n", Z, "\n"

import matplotlib.pyplot as plt
plt.xlim(250,400)
plt.ylim(-0.3,0.0)
plt.plot(F[:,0],F[:,1])
plt.show()
np.savetxt('F_conf.dat',np.c_[F[:,0],F[:,1]])
