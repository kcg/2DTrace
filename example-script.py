import numpy as np
from Trace2D import *
import pylab as pl

def fsurf(x):
	return 0.1*np.cos(4*x)+0.5

def bsurf(x):
	return -0.1*np.cos(4*x)+1
	
def df(x):
	return -0.1*np.sin(4*x)*4

def db(x):
	return 0.1*np.sin(4*x)*4
	

world = World(box=(2,1.6), n=1.0)
source = Source(pos=(0.5,0.01), width=1, angdist=None, intdist=None)
detector = Detector(pos=(0,1.5), width=2)
objectList = [Object(fsurf, bsurf, dfsurf=df, dbsurf=db, n=1.3)]

visualizer = Visualizer(world, source, objectList, detector)

runs = 100
for i in range(runs):
	# Set visualize to "False" to speed up simulation for high numbers of runs
	r = Ray(world, source, objectList, detector, visualize=True)
	r.trace()
	del r

pl.ioff()
pl.show()

vec, pos, startpos, detected, lost = detector.statistics()

angle = [np.arccos(v[0])/np.pi*180.0-90.0 for v in vec]
pl.hist(pos, 50)
pl.xlabel("Position")
pl.ylabel("Counts")
pl.show()
