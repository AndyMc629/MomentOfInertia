import numpy as np
from StringIO import StringIO

"""
with open('CH3NH3.xyz') as data:
	for line in data:
		if line(0)='C':
			print "Carbon"

"""
data_raw=np.genfromtxt('CH3NH3.xyz',dtype=("|S10", float, float, float),skip_header=2)
print data_raw.shape
print data_raw
data_with_mass = np.zeros(len(data_raw),dtype=[('Atom','S'),
						 ('x','float'),
						 ('y','float'),	
						 ('z','float'),
						 ('Mass','float')] )


