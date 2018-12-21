import numpy as np
import struct

# This is a short script to display the indices of those images which are the greatest outliers.
# These images match their prototypical image the least.

# Important variables
numProto = 91               # Number of prototypes in the SOM
numOutliers = 3             # Number of outlier indices to display
mapF = "stampmap150k.bin"   # Mapping file. Must be in this directory.

# Find the indices of the least fitting images

with open(mapF, 'rb') as inputStream:
    
    numImages = struct.unpack("i", inputStream.read(4))[0]
    
    width = struct.unpack("i", inputStream.read(4))[0]
    height = struct.unpack("i", inputStream.read(4))[0]
    depth = struct.unpack("i", inputStream.read(4))[0]

    mapping = np.ones( (numImages, numProto) )
    for i in range(numImages):
        for t in range(numProto):
            mapping[i,t] = struct.unpack_from("f", inputStream.read(4))[0]

# Find and priont outlier indices.
distances = np.min(mapping, axis=1)         # Array of the distance to the closest prototype for each image.
worstIndices = np.argsort(distances)[::-1]  # Array of the indices, sorted by distance descending.

for ind in worstIndices[:numOutliers]: print(ind)
