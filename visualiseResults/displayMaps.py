# Display the SOM or heatmaps using these functions.
from heatmapVisualiser import *
import os
import glob
import numpy as np

def viewHM(path):
	m = MAPVisualizer(path, path)
	print(path)
	m.readMap()
	m.showMap()
	raw_input("Press ENTER to continue.")
	
def viewMultipleHM(paths):
	for p in paths:
		viewHM(p)
		
def viewAllInFolder(path):
	#os.chdir(path)
	f = []
	for e in glob.glob( path + "*.bin"): f.append(e)
	f = np.array(f)
	
	f = f[np.argsort(f)]
	
	viewMultipleHM(f)
