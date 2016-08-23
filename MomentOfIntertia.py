import numpy as np
from StringIO import StringIO

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
print "Principal moments of inertia = ", np.linalg.eig(I)[0], "amu Angstrom^2"

print "Product I_A*I_B*I_C = ", np.linalg.eig(I)[0][0]*np.linalg.eig(I)[0][1]*np.linalg.eig(I)[0][2], "amu^3 Angstrom^6\n"
print "                    = ", 1.1087*1E-84*np.linalg.eig(I)[0][0]*np.linalg.eig(I)[0][1]*np.linalg.eig(I)[0][2], "eV s^2"

print "\n " 
