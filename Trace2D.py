import numpy as np
from scipy.optimize import fsolve
from scipy.misc import derivative
import pylab as pl
import matplotlib

numerror = 0.00001
maxsteps = 10


def distance(pos1,pos2):
	if pos1 == None or pos2 == None:
		return 999999999
	return np.sqrt((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2)


class World:
	
	def __init__(self, box, n):
		self.size = box
		self.refractiveindex = n
		
	def inside(self, pos):
		if pos[0] < 0 or pos[0] > self.size[0] or pos[1] < 0 or pos[1] > self.size[1]:
			return False
		else:
			return True

		
class Visualizer:
	
	def __init__(self, world, source, objList, detector):
		self.world = world
		self.source = source
		self.detector = detector
		self.objects = objList
		self.drawstatics()
		
	def drawstatics(self):
		pl.ion()
		self.fig = pl.figure()
		self.ax = self.fig.add_subplot(111)
		source = matplotlib.patches.Rectangle(self.source.position, self.source.width, self.source.width*0.01, color='orange')
		self.ax.add_patch(source)
		det = matplotlib.patches.Rectangle(self.detector.position, self.detector.size, self.detector.size*0.01, color='grey')
		self.ax.add_patch(det)
		for obj in self.objects:
			obj.visualize(self.world.size[0])
		pl.xlim(0,self.world.size[0])
		pl.ylim(0,self.world.size[1])
		pl.axes().set_aspect('equal')
		pl.draw()
		

class Source:
	
	def __init__(self, pos, width, angdist, intdist):
		self.position = pos
		self.width = width
		self.angles = angdist
		self.intensities = intdist
		self.start = 0
		self.count = 0
		
	def emit(self):
		#pos = self.position + np.array([self.width/1000*(1000-self.count),0]) 
		#pos = np.array([np.random.normal(self.position[0], self.width), 0])
		pos = self.position + np.array([np.random.uniform() * self.width, 0])
		self.count += 1
		self.start = pos
		vec = np.array([0,1])
		return (pos,vec)


class Detector:
	
	def __init__(self, pos, width, qe=1):
		self.position = pos
		self.size = width
		self.efficiency = qe
		self.lostrays = 0
		self.detectedrays = 0
		self.detectedvec = []
		self.detectedpos = []
		self.startedpos = []
		
	def intersect(self, pos, vec):
		intersectionx = vec[0]/vec[1]*(self.position[1]-pos[1])+pos[0]
		if self.position[0] <= intersectionx <= self.position[0]+self.size:
			return intersectionx
		else:
			return False
			
	def showlog(self):
		print(self.lostrays)
		print(self.detectedrays)
		print(self.detectedvec)
		print(self.detectedpos)
		
	def statistics(self):
		return [ self.detectedvec, self.detectedpos, self.startedpos, self.detectedrays, self.lostrays ]
		

class Object:
	
	def __init__(self, fsurf, bsurf, dfsurf=None, dbsurf=None, n=1.5):
		self.frontsurface = fsurf
		self.dfrontsurface = dfsurf
		self.backsurface = bsurf
		self.dbacksurface = dbsurf
		self.refractiveindex = n
		
	def intersection(self,pos,vec):
		if vec[0] == 0:
			secpos0_front = pos[0]
			secpos0_back = pos[0]
		else:
			try:
				secpos0_front = fsolve(lambda x : self.frontsurface(x) - (vec[1]/vec[0]*(x-pos[0])+pos[1]),pos[0])[0]
				secpos0_back = fsolve(lambda x : self.backsurface(x) - (vec[1]/vec[0]*(x-pos[0])+pos[1]),[pos[0],pos[0],pos[0]])[0]
			except RuntimeWarning:
				return None
			if self.frontsurface(secpos0_front) - (vec[1]/vec[0]*(secpos0_front-pos[0])+pos[1]) > numerror:
				return None
			if self.backsurface(secpos0_back) - (vec[1]/vec[0]*(secpos0_back-pos[0])+pos[1]) > numerror:
				return None
		secpos_front = np.array([secpos0_front, self.frontsurface(secpos0_front)])
		secpos_back = np.array([secpos0_back, self.backsurface(secpos0_back)])
		# Check which section point lies "in front" of the ray
		k_front = vec[1]/(secpos_front-pos)[1]
		k_back = vec[1]/(secpos_back-pos)[1]
		#print k_front
		#print k_back
		#pl.plot
		if k_front > 0 and k_back < 0:
			secpos = secpos_front
		elif k_front < 0 and k_back > 0:
			secpos = secpos_back
		elif k_front < 0 and k_back < 0:
			return None
		else:
			# If both are possible, take the nearest
			if distance(secpos_front,pos) < distance(secpos_back,pos):
				secpos = secpos_front
			else:
				secpos = secpos_back
		return secpos
	
	def inside(self, pos):
		if self.frontsurface(pos[0]) < pos[1] and self.backsurface(pos[0]) > pos[1]:
			return True
		else:
			return False
			
	def normal(self, pos):
		if np.abs(self.frontsurface(pos[0]) - pos[1]) < np.abs(self.backsurface(pos[0]) - pos[1]):
			if self.dfrontsurface == None:
				df = derivative(self.frontsurface, pos[0], dx=1e-10)
			else:
				df = self.dfrontsurface(pos[0])
			if df == 0:
				return np.array([0,1])
			m = -1/df
		else:
			if self.dbacksurface == None:
				df = derivative(self.backsurface, pos[0], dx=1e-10)
			else:
				df = self.dbacksurface(pos[0])
			if df == 0:
				return np.array([0,1])
			m = -1/df
		normalvec = np.array([1,m])
		normalvec = normalvec / np.linalg.norm(normalvec)
		#z = pos + normalvec*0.01
		#x = [pos[0],z[0]]
		#y = [pos[1],z[1]]
		#pl.plot(x,y,"r-")
		return normalvec
		
	def visualize(self, width):
		x = np.linspace(0,width,1000)
		pl.plot(x,self.frontsurface(x),"b-")
		pl.plot(x,self.backsurface(x),"b-")
		pl.fill_between(x, self.frontsurface(x), self.backsurface(x), facecolor='blue', alpha=0.15)
		

class Ray:
	
	def __init__(self, world, source, objectList, detector, visualize=False):
		self.pos = None
		self.vec = None
		self.active = False
		self.world = world
		self.currentn = world.refractiveindex
		self.source = source
		self.detector = detector
		self.objectList = objectList
		self.visualize = visualize
		if self.visualize:
			self.pathx = []
			self.pathy = []
			self.line, = pl.plot(self.pathx,self.pathy,"k-", alpha=0.25)
			#self.dots, = pl.plot(self.pathx,self.pathy,"ko", alpha=0.25)
			pl.draw()
		
	def generate(self):
		self.active = True
		self.pos, self.vec = self.source.emit()
		if self.visualize:
			self.pathx.append(self.pos[0])
			self.pathy.append(self.pos[1])
			self.line.set_xdata(self.pathx)
			self.line.set_ydata(self.pathy)
			#self.dots.set_xdata(self.pathx)
			#self.dots.set_ydata(self.pathy)
			pl.draw()
		
	def nextobjectintersection(self):
		intersectionpoint = None
		intersectionobject = None
		for i,obj in enumerate(self.objectList):
			point = obj.intersection(self.pos,self.vec)
			if point == None:
				continue
			else:
				if distance(self.pos,point) < distance(self.pos,intersectionpoint):
					intersectionpoint = point
					intersectionobject = i
		if intersectionpoint == None:
			return [None, None]
		else:
			return [intersectionpoint, intersectionobject]
		
	def refract(self,n1,n2,normalvec):
		u = np.array([self.vec[0],self.vec[1]])
		n = np.array([normalvec[0],normalvec[1]])
		w = (np.dot(u,n))**2 + (n2/n1)**2 - 1
		#total refelection
		if w < 0:
			#print 'tot ref.'
			v = u - 2*np.dot(n,u)*n
			self.currentn = n1
		#refraction
		else:
			k = np.sign(np.dot(u,n)) * np.sqrt(w) - np.dot(u,n) 
			v = (u + k*n) / np.linalg.norm(u + k*n)
			self.currentn = n2
		if self.visualize:
			self.pathx.append(self.pos[0])
			self.pathy.append(self.pos[1])
			self.line.set_xdata(self.pathx)
			self.line.set_ydata(self.pathy)
			#self.dots.set_xdata(self.pathx)
			#self.dots.set_ydata(self.pathy)
			pl.draw()
		#z = self.pos + v*0.01
		#x = [self.pos[0],z[0]]
		#y = [self.pos[1],z[1]]
		#pl.plot(x,y,"r-")
		#set ray's new vec and move it a tiny bit into the new direction
		self.vec = v
		self.pos = self.pos + self.vec*1e-3
	
	def checkdetector(self):
		dectpos = self.detector.intersect(self.pos,self.vec)
		if dectpos:
			self.active = False
			self.detector.detectedrays += 1
			self.detector.detectedvec.append(self.vec)
			self.detector.detectedpos.append(dectpos)
			self.detector.startedpos.append(self.source.start)
			if self.visualize:
				self.pathx.append(dectpos)
				self.pathy.append(self.detector.position[1])
				self.line.set_xdata(self.pathx)
				self.line.set_ydata(self.pathy)
				#self.dots.set_xdata(self.pathx)
				#self.dots.set_ydata(self.pathy)
				pl.draw()
		
	def checkoutbound(self):
		#if not self.world.inside(self.pos,self.vec):
		if self.active:
			self.active = False
			self.detector.lostrays += 1
				
	def trace(self):
		self.generate()
		step = 0
		while self.active and step < maxsteps:
			step += 1
			# Get next intersection
			nextpos, nextobj = self.nextobjectintersection()
			#print nextpos
			# If there is no intersection, check for detection or leaving the world
			if nextpos == None or self.world.inside(nextpos) == False:
				self.checkdetector()
				self.checkoutbound()
			#Change pos and then vec according to snell's law
			else:
				obj = self.objectList[nextobj]
				if self.currentn == self.world.refractiveindex:
					newn = obj.refractiveindex
				else:
					newn = self.world.refractiveindex
				self.pos = nextpos
				normalvec = obj.normal(nextpos)
				self.refract(self.currentn,newn,normalvec)
		
