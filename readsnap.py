"""
Example to read snapshot from FIRE simulation
Author: Pouria A. Mistani
email: p.a.mistani@gmail.com
"""
import numpy as np
import matplotlib.pyplot as plt
import snapHDF5 as snap
import pdb
import PySmoother as PyS
import matplotlib.pyplot as plt









filename = '../../../snapshot_440'
filename = '../m12qq_hr_Dec16_2013/snapshot_440'
## -- read stellar components ---
pos = snap.read_block(filename, "POS ", parttype=4) #Notice that to read POS and VEL of gas and DM, only change parttype=0 or 1, respectively
print 'pase pos'
vel = snap.read_block(filename, "VEL ", parttype=4)
print 'pase vel'
mass = snap.read_block(filename, "MASS", parttype=4)
print 'pase mass'
'''
age = snap.read_block(filename, "AGE ", parttype=4)
print 'pase age'
metal = snap.read_block(filename, "Z   ", parttype=4)
print 'pase metal'
## -- read some gas-exclusive info ----- ##
sfr = snap.read_block(filename, "SFR ", parttype=0)
print 'pase sfr'
'''
print 'done reading',pos.shape




X = pos[:,0]
Y = pos[:,1]
Z = pos[:,2]
Vx = vel[:,0]
Vy = vel[:,1]
Vz = vel[:,2]
mass = np.ones(len(X))

#---
Animation = True
	
if Animation:
		image = PyS.PyS(X, Y, Z, mass, Vx, Vy,Vz, r=2, w=500, h=500, offset=10)
		
		image.align()	
		
		img = image.smoothStars(plane='xz')
		extent = image.extent()
			
		plt.figure()			
		ax = plt.subplot(111)
		ax.imshow(img, extent = extent)
			
		
		plt.savefig('test.png')
		


#r, g, b = img.split()

pdb.set_trace()
