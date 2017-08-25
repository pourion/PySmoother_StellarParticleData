#convert -delay 100 -loop 0 *.png reference.g
"""
    PySmoother version 2.0
    Utilizing Stellar Particles to make RGB color animations.
    Author: Pouria A. Mistani
    email: p.a.mistani@gmail.com
"""
import numpy as np
import scipy.ndimage as ndi
from PIL import Image
import colorsys
import pdb
import matplotlib.pyplot as plt
from scipy import spatial
import scipy
#-----------------------------------------------------------
class PyS:
	def __init__(self, xpos, ypos, zpos, mass, Velx, Vely, Velz, r=2, w=1000, h=1000, offset=5):
		self.Vx = Velx
		self.Vy = Vely
		self.Vz = Velz
		self.Xpos = xpos
		self.Ypos = ypos
		self.Zpos = zpos		
		self.mass = mass
		self.w = w
		self.h = h
		self.r = r
		self.offset = offset
		self.centering()

	def centering(self):
		xcen = np.median(self.Xpos)
		ycen = np.median(self.Ypos)
		zcen = np.median(self.Zpos)		
		for i in range(20):
				xcen = np.median(self.Xpos[abs(self.Xpos-xcen)<0.3*self.offset])
				ycen = np.median(self.Ypos[abs(self.Ypos-ycen)<0.3*self.offset])
				zcen = np.median(self.Zpos[abs(self.Zpos-zcen)<0.3*self.offset])
				self.xcen = xcen
				self.ycen = ycen
				self.zcen = zcen
		self.Xpos -= self.xcen
		self.Ypos -= self.ycen
		self.Zpos -= self.zcen

	def centeringXY(self):
		xcen = np.median(self.Xpos)
		ycen = np.median(self.Ypos)
		for i in range(10):
				xcen = np.median(self.Xpos[abs(self.Xpos-xcen)<0.3*self.offset])
				ycen = np.median(self.Ypos[abs(self.Ypos-ycen)<0.3*self.offset])
				self.xcen = xcen
				self.ycen = ycen

		self.limits = [self.xcen - self.offset,self.xcen + self.offset, self.ycen - self.offset, self.ycen + self.offset]

	def centeringXZ(self):
		xcen = np.median(self.Xpos)
		ycen = np.median(self.Zpos)
		for i in range(10):
				xcen = np.median(self.Xpos[abs(self.Xpos-xcen)<0.3*self.offset])
				ycen = np.median(self.Zpos[abs(self.Zpos-ycen)<0.3*self.offset])
				self.xcen = xcen
				self.ycen = ycen

		self.limits = [self.xcen - self.offset,self.xcen + self.offset, self.ycen - self.offset, self.ycen + self.offset]

	def centeringYZ(self):
		xcen = np.median(self.Zpos)
		ycen = np.median(self.Ypos)
		for i in range(10):
				xcen = np.median(self.Zpos[abs(self.Zpos-xcen)<0.3*self.offset])
				ycen = np.median(self.Ypos[abs(self.Ypos-ycen)<0.3*self.offset])
				self.xcen = xcen
				self.ycen = ycen

		self.limits = [self.xcen - self.offset,self.xcen + self.offset, self.ycen - self.offset, self.ycen + self.offset]


	def extent(self):
		return self.limits


	def smooth(self, plane = 'xy'):
	    print 'smoothing ...'
	    if plane == 'xy':
	    	X = np.array([x for x, y, z in zip(self.Xpos, self.Ypos, self.Zpos)])
	    	Y = np.array([y for x, y, z in zip(self.Xpos, self.Ypos, self.Zpos)])
		Z = np.array([z for x, y, z in zip(self.Xpos, self.Ypos, self.Zpos)])
		self.centeringXY()
	    if plane == 'xz':
	    	X = np.array([x for x, y in zip(self.Xpos, self.Zpos)])
	    	Y = np.array([y for x, y in zip(self.Xpos, self.Zpos)])
		self.centeringXZ() 
	    if plane == 'yz':
	    	X = np.array([x for x, y in zip(self.Ypos, self.Zpos)])
	    	Y = np.array([y for x, y in zip(self.Ypos, self.Zpos)])
	        self.centeringYZ()

	    MASS = self.mass
	    velX = self.Vx
	    velY = self.Vy
	    velZ = self.Vz
	    x0, x1, y0, y1 = self.limits
	    idx = (X < x1) * (X > x0)*(Y < y1)*(Y>y0)
	    X = X[idx]
	    Y = Y[idx]
	    Z = Z[idx]
            MASS = MASS[idx]
            velX = velX[idx]
	    velY = velY[idx]
	    velZ = velZ[idx]
	    
	    
	    
	    
	    kx = (self.w - 1) / float(x1 - x0)
	    ky = (self.h - 1) / float(y1 - y0)
            self.kx = kx
	    self.ky = ky


	    print "start building tree ..."
	    tree = scipy.spatial.cKDTree(zip(X,Y,Z),leafsize=100)
	    print "tree built!"
	    neighbors = []
	    density = []
	    sigma = []
	    print "start sph measurements"
	    for item in zip(X,Y,Z):
	    	distance, idx = tree.query(item, k=64, distance_upper_bound=10)
	    	perm_idx = distance < np.inf
	    	ngb_idx = idx[perm_idx]
		neighbors.append([ngb_idx])
		density.append(np.sum(MASS[ngb_idx])/np.max(distance[perm_idx])**3.0)   
	
		sigma.append((np.std(velX[ngb_idx])**2 + np.std(velY[ngb_idx])**2 + np.std(velZ[ngb_idx])**2)**0.5/np.sqrt(3))
	    print "done sph measurements!"
	    density = np.array(density)
	    density[density == np.inf] = np.max(density[density < np.inf])
	
	    sigma = np.array(sigma)
	    sigma[sigma == np.inf] = np.max(sigma[sigma < np.inf])

	    img = [1] * (self.w * self.h)
	    imgSigma = [0 for i in range(self.w * self.h)]
	    print "building image values/hues"
	    for x, y, p, s in zip(X,Y, density, sigma):
	                ix = int((x - x0) * kx)
		        iy = int((y - y0) * ky)
	                img[iy * self.w + ix] += p**2   # density
			
	                imgSigma[iy * self.w + ix] += s*p**2   #velocity dispersion
	    print "image values/hues built!"




	    img = np.array(img)
   	    imgSigma = np.array(imgSigma)
   	    
	    
	    

	    img2 = imgSigma/img
	    img3 = (img2 - np.min(img2))/(np.max(img2) - np.min(img2))    

	    '''
	    img3 = 0.7*img3 + 0.48   ### arbitrary shift/scaling of color hue to go from blue to red
	    img3[img3 > 1.] -= 1.0
	    '''


	    img = np.log10(img)
	    img = (img - np.min(img))/(np.max(img) - np.min(img))


	    R = []#np.zeros(self.w*self.h, 'uint8')
	    G = []#np.zeros(self.w*self.h, 'uint8')
	    B = []#np.zeros(self.w*self.h, 'uint8')
	    print "start getting RGB colors"
	    for h, v in zip(img3, img):
		#(r, g, b) = HSV_2_RGB(int(h*255), 125, int(v*255))
	        #saturation = (90 + np.random.rand() * 10)/100.
                (r, g, b) = colorsys.hls_to_rgb(h, v, 0.9)

	    	#(r, g, b) = colorsys.hsv_to_rgb(h, .5, v)	
	    	R.append(r)
	    	G.append(g)
	    	B.append(b)
	 
	    R = np.array(R).reshape(self.w,self.h)
	    G = np.array(G).reshape(self.w,self.h)
	    B = np.array(B).reshape(self.w,self.h)
	    print "done RGB colors!"
	    
	    
 	    R = ndi.gaussian_filter(R, (self.r,self.r))
	    R = np.flipud(R)
 	    print 'R smoothed!'
 	    G = ndi.gaussian_filter(G, (self.r,self.r))
	    G = np.flipud(G)
 	    print 'G smoothed!'
 	    B = ndi.gaussian_filter(B, (self.r,self.r))
	    B = np.flipud(B)
	    print 'B smoothed!'
	    
	    

	    rgbArray = np.zeros((self.w,self.h,3), 'uint8')
            rgbArray[..., 0] = R*256
            rgbArray[..., 1] = G*256
            rgbArray[..., 2] = B*256

            img = Image.fromarray(rgbArray)

	    self.image = img
	    return img

	
	def rotateX(self, theta):
		x_new = []
		y_new = []
		z_new = []
		theta *= np.pi/180
		for x,y,z in zip(self.Xpos, self.Ypos, self.Zpos):
			x_new.append(x)
			y_new.append(y*np.cos(theta) - z*np.sin(theta))
			z_new.append(y*np.sin(theta) + z*np.cos(theta))		
		self.Xpos = np.array(x_new)
		self.Ypos = np.array(y_new)
		self.Zpos = np.array(z_new)


	def rotateY(self, theta):
		x_new = []
		y_new = []
		z_new = []
		theta *= np.pi/180
		for x,y,z in zip(self.Xpos, self.Ypos, self.Zpos):
			x_new.append(x*np.cos(theta) + z*np.sin(theta))
			y_new.append(y)
			z_new.append(-x*np.sin(theta) + z*np.cos(theta))		
		self.Xpos = np.array(x_new)
		self.Ypos = np.array(y_new)
		self.Zpos = np.array(z_new)

	def rotateZ(self, theta):
		x_new = []
		y_new = []
		z_new = []
		theta *= np.pi/180
		for x,y,z in zip(self.Xpos, self.Ypos, self.Zpos):
			x_new.append(x*np.cos(theta) - y*np.sin(theta))
			y_new.append(x*np.sin(theta) + y*np.cos(theta))
			z_new.append(z)		
		self.Xpos = np.array(x_new)
		self.Ypos = np.array(y_new)
		self.Zpos = np.array(z_new)


	def Virial_Cut(self, r_vir):
		mass = self.mass
		X = self.Xpos
		Y = self.Ypos
		Z = self.Zpos
		velX = self.Vx
	    	velY = self.Vy
	    	velZ = self.Vz

		xcen = np.median(self.Xpos)
		ycen = np.median(self.Ypos)
		zcen = np.median(self.Zpos)		
		for i in range(10):
				xcen = np.median(self.Xpos[abs(self.Xpos-xcen)<0.3*self.offset])
				ycen = np.median(self.Ypos[abs(self.Ypos-ycen)<0.3*self.offset])
				zcen = np.median(self.Zpos[abs(self.Zpos-zcen)<0.3*self.offset])
		rad = ((X-xcen)**2+(Y-ycen)**2+(Z-zcen)**2)**0.5
		idx = rad < 1.1*r_vir
		M_vir = np.sum(self.mass[rad <= r_vir])
		self.Xpos = X[idx]
		self.Ypos = Y[idx]
		self.Zpos = Z[idx]
		self.mass = mass[idx]
		self.Vx = velX[idx]
		self.Vy = velY[idx]
		self.Vz = velZ[idx]


		return M_vir

	def align(self):
		print "start alignment ..."
		X = self.Xpos
		Y = self.Ypos
		Z = self.Zpos
		vx = self.Vx
		vy = self.Vy
		vz = self.Vz
		mass = self.mass

		L = mass*np.cross([X, Y, Z],[vx, vy, vz],axis=0)
		Lx = np.sum(L[0])
		Ly = np.sum(L[1])
		Lz = np.sum(L[2])
		lx = float(Lx)/np.linalg.norm([Lx, Ly, Lz])
		ly = float(Ly)/np.linalg.norm([Lx, Ly, Lz])
		lz = float(Lz)/np.linalg.norm([Lx, Ly, Lz])

		v = np.cross([lx,ly,lz], [0,0,1], axis=0)
		s = np.linalg.norm(v)
		c = np.dot([0,0,1], [lx,ly,lz])
		a = (1- c)/(s*s)


		R11 = 1 - a*(v[1]**2 + v[2]**2)
		R22 = 1 - a*(v[0]**2 + v[2]**2)
		R33 = 1 - a*(v[0]**2 + v[1]**2)

		R12 = -v[2] + a*v[0]*v[1]
		R13 = v[1] + a*v[0]*v[2]
		
		R21 = v[2] + a*v[0]*v[1]		
		R23 = -v[0] + a*v[1]*v[2]

		R31 = -v[1] + a*v[0]*v[2]
		R32 = v[0] + a*v[1]*v[2]

		Xp = []
		Yp = []
		Zp = []

		VXp = []
		VYp = []
		VZp = []

		for x, y, z, v_x, v_y, v_z in zip(X, Y, Z, vx, vy, vz):
			Xp.append(R11*x + R12*y + R13*z)
			Yp.append(R21*x + R22*y + R23*z)
			Zp.append(R31*x + R32*y + R33*z)
			VXp.append(R11*v_x + R12*v_y + R13*v_z)
			VYp.append(R21*v_x + R22*v_y + R23*v_z)
			VZp.append(R31*v_x + R32*v_y + R33*v_z)
		self.Xpos = np.array(Xp)
		self.Ypos = np.array(Yp)
		self.Zpos = np.array(Zp)
		
		self.Vx = np.array(VXp)
		self.Vy = np.array(VYp)
		self.Vz = np.array(VZp)
		
		print "done alignment!"	



	def smoothStars(self, plane = 'xy'):
	    print 'smoothing ...'
	    if plane == 'xy':
	    	X = np.array([x for x, y, z in zip(self.Xpos, self.Ypos, self.Zpos)])
	    	Y = np.array([y for x, y, z in zip(self.Xpos, self.Ypos, self.Zpos)])
		Z = np.array([z for x, y, z in zip(self.Xpos, self.Ypos, self.Zpos)])
		self.centeringXY()
	    if plane == 'xz':
	    	X = np.array([x for x, y, z in zip(self.Xpos, self.Zpos, self.Ypos)])
	    	Y = np.array([y for x, y, z in zip(self.Xpos, self.Zpos, self.Ypos)])
	    	Z = np.array([z for x, y, z in zip(self.Xpos, self.Zpos, self.Ypos)])
		self.centeringXZ() 
	    if plane == 'yz':
	    	X = np.array([x for x, y, z in zip(self.Ypos, self.Zpos, self.Xpos)])
	    	Y = np.array([y for x, y, z in zip(self.Ypos, self.Zpos, self.Xpos)])
	    	Z = np.array([z for x, y, z in zip(self.Ypos, self.Zpos, self.Xpos)])
	        self.centeringYZ()

	    MASS = self.mass
	    x0, x1, y0, y1 = self.limits
	    idx = (X < x1) * (X > x0)*(Y < y1)*(Y>y0)
	    X = X[idx]
	    Y = Y[idx]
	    Z = Z[idx]
            MASS = MASS[idx]
	    
	    
	    
	    
	    kx = (self.w - 1) / float(x1 - x0)
	    ky = (self.h - 1) / float(y1 - y0)
            self.kx = kx
	    self.ky = ky


	    print "start building tree ..."
	    tree = scipy.spatial.cKDTree(zip(X,Y,Z),leafsize=100)
	    print "tree built!"
	    neighbors = []
	    density = []
	    sigma = []
	    print "start sph measurements"
	    for item in zip(X,Y,Z):
	    	distance, idx = tree.query(item, k=64, distance_upper_bound=10)
	    	perm_idx = distance < np.inf
	    	ngb_idx = idx[perm_idx]
		neighbors.append([ngb_idx])
		density.append(np.sum(MASS[ngb_idx])/np.max(distance[perm_idx])**3.0)   
	
	    print "done sph measurements!"
	    density = np.array(density)
	    density[density == np.inf] = np.max(density[density < np.inf])
	

	    img = [1] * (self.w * self.h)
	    print "building image values"
	    for x, y, p in zip(X,Y, density):
	                ix = int((x - x0) * kx)
		        iy = int((y - y0) * ky)
	                img[iy * self.w + ix] += p**2   # density
			
	    print "image values built!"

	    img = np.array(img)
	    img = np.log10(img)
	    img = (img - np.min(img))/(np.max(img) - np.min(img))

	    img *= 255.0
	    img = np.array(img).reshape(self.w,self.h)
	    
	    
 	    img = ndi.gaussian_filter(img, (self.r,self.r))
	    img = np.flipud(img)
 	    print 'image smoothed!'
	    self.image = img
	    return img









def HSV_2_RGB(H,S,V):
    # Converts an integer HSV tuple (value range from 0 to 255) to an RGB tuple 

    # Check if the color is Grayscale
    if S == 0:
        R = V
        G = V
        B = V
        return (R, G, B)

    # Make hue 0-5
    region = H // 43;

    # Find remainder part, make it from 0-255
    remainder = (H - (region * 43)) * 6; 

    # Calculate temp vars, doing integer multiplication
    P = (V * (255 - S)) >> 8;
    Q = (V * (255 - ((S * remainder) >> 8))) >> 8;
    T = (V * (255 - ((S * (255 - remainder)) >> 8))) >> 8;


    # Assign temp vars based on color cone region
    if region == 0:
        R = V
        G = T
        B = P

    elif region == 1:
        R = Q; 
        G = V; 
        B = P;

    elif region == 2:
        R = P; 
        G = V; 
        B = T;

    elif region == 3:
        R = P; 
        G = Q; 
        B = V;

    elif region == 4:
        R = T; 
        G = P; 
        B = V;

    else: 
        R = V; 
        G = P; 
        B = Q;


    return (R, G, B)

