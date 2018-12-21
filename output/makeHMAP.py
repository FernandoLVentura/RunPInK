
import numpy as np
import struct
import os
import os.path


# This script contains functions used to generate the hexagonal heatmap 
# in the same format as the SOM is saved in and saves it to disk in binary format


#Make a heat map element from a value and its standar deviation
#Assume the value is in [0,1] and the standard deviation in [0,1]
def createCell(v, s, d=64):
	
    val  = v
    sd = s

	# Cell size is 64x64 by default
    e = d
	
    img = np.zeros( (e,e,4) )
	#make the element:
    el = (1, 1, 1, 0)
	
    if( not np.isnan(val) ):
		#Convert Std. Dev. to range [0.5, 1]
        r = val
        b = (1-val)
        el = (r, 0, b, sd)
	
    img[:,:] = el
    
    return img

#Make a heatmap from an array of values and their standard deviation
#Assume input in any range including NaNs.
def generateMap( cellValues, stdValues, dimensions = 64):
    
    mn = np.nanmin(cellValues)
    smn = np.nanmin(stdValues)
    mx = np.nanmax(cellValues)
    smx = np.nanmax(stdValues)
	
    rng = mx - mn
    srng = smx - smn
	
    cellValues = np.array(cellValues)
    stdValues = np.array(stdValues)
    
    if (rng > 0):
        vals = (cellValues-mn)/rng;
    else:
        vals = np.ones( len(cellValues) );
        
    # St. Dev. is inversed at this point so that a higher number is better.
    if(srng > 0):
        stds = (stdValues-smn)/srng
        stds = (1 - stds)/2  + 0.5
    else:
        stds = np.ones( len(stdValues) )
	
    cells = []
    for i in range(len(cellValues)):
        #Make the cell
        cell = createCell( vals[i], stds[i], dimensions)
        cells.append(cell)
		
    cells = np.array(cells)
    return cells

#Save a map passed in.
#Takes in an array of cells in order
def saveMap( cellArr, targetPath, dimc = 64, dimm = 11):
	
    numCells = len(cellArr)
    m = open( targetPath, 'wb')		#Open: Write bytes mode.
	
    m.write(struct.pack('i', 1))	#No. Channels
    m.write(struct.pack('i', dimm))	#Map width
    m.write(struct.pack('i', dimm))	#Map height
    m.write(struct.pack('i', 1))	#Map depth
    
    m.write(struct.pack('i', dimc))	#Cell width
    m.write(struct.pack('i', dimc))	#Cell height
    
    for c in cellArr:
        c.astype('f').tofile(m)
    


#Generate and save a heatmap
def writeHeatmap( cVals, stdVals, path, dmap=11, dcell=64 ):
	
    map = generateMap( cVals, stdVals, dcell)
    saveMap(map, path, dcell, dmap)

    

