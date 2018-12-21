#from future import division, print_function
import struct
import numpy
from matplotlib import pyplot
import matplotlib
import math
import sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--function", dest="outputFunction",
        help="Choose between lin (linear) or log (logarithmic), default = lin", metavar="FUN")
parser.add_option("-q", "--quiet",
        action="store_false", dest="verbose", default=True,
        help="don't print status messages to stdout")

outputFunction = ''
(options, args) = parser.parse_args()
print(outputFunction)

class MAPVisualizer():
    def __init__(self, filePathName, fileName):
        self.__filePathName = filePathName
        self.__fileName = fileName
        self.__numberOfChannels = 0
        self.__somWidth = 0
        self.__somHeight = 0
        self.__somDepth = 0
        self.__neuronWidth = 0
        self.__neuronHeight = 0
        self.__neurons = []

    def readMap(self):
        inputStream = open(self.__filePathName, 'rb')
        self.__numberOfChannels = struct.unpack("i", inputStream.read(4))[0]
        self.__somWidth = struct.unpack("i", inputStream.read(4))[0]
        self.__somHeight = struct.unpack("i", inputStream.read(4))[0]
        self.__somDepth = struct.unpack("i", inputStream.read(4))[0]
        self.__neuronWidth = struct.unpack("i", inputStream.read(4))[0]
        self.__neuronHeight = struct.unpack("i", inputStream.read(4))[0]

        print ("width: " + str(self.__somWidth))
        print ("height: " + str(self.__somHeight))
        print ("depth: " + str(self.__somDepth))
        print ("neurons: " + str(self.__neuronWidth) +"x" + str(self.__neuronHeight))
        try:
            while True:
                data = numpy.ones(self.__neuronWidth * self.__neuronHeight)
                for i in range(self.__neuronWidth * self.__neuronHeight):
                    data[i] = struct.unpack_from("f", inputStream.read(4))[0]
                self.__neurons.append(data)

        except:
            inputStream.close()
        self.__neurons = numpy.array(self.__neurons)
        print (str(len(self.__neurons)) + " neurons loaded")

    def isHexMap(self):
        return len(self.__neurons) < self.__somHeight * self.__somWidth

    def calculateMap(self, neurons, shareIntensity = True, border = 0, shape="box"):
        if shape == "box":
            mapSize = numpy.array([self.__somWidth, self.__somHeight])
            neuronSize = numpy.array([self.__neuronWidth, self.__neuronHeight])
            size = numpy.multiply(mapSize,numpy.array(neuronSize) + border) + border
            image = numpy.empty(size)
            image[:] = numpy.NAN
            for x in range(self.__somWidth):
                for y in range(self.__somHeight):
                    data = neurons[x + y*self.__somWidth].reshape(self.__neuronWidth, self.__neuronHeight)
                    if not shareIntensity:
                        if numpy.max(data)-numpy.min(data) != 0:
                            data = 1.0 * (data - numpy.min(data)) / (numpy.max(data) - numpy.min(data))
                    image[x*(self.__neuronWidth + border): (x+1) * (self.__neuronWidth + border) -border, y * (self.__neuronHeight + border): (y+1) * (self.__neuronHeight + border) - border] = data 
#             for position in indices:
#                 data = neurons[position]
#                 print data.shape


#                 for pos in numpy.ndindex(data.shape):
#                     image[(numpy.multiply(position, numpy.array(neuronSize) + border) + pos)[0]+border][(numpy.multiply(position, numpy.array(neuronSize) + border) + pos)[1]+border] = data[pos]
            return image
        if shape == "hex":
            size = numpy.multiply((self.__somWidth+0.5,self.__somHeight),numpy.array((self.__neuronWidth, self.__neuronHeight)) + border) + border
            size[1] = math.ceil(size[1] - (self.__somHeight-1.0) * self.__neuronHeight / 4.0)
            size[0] = math.ceil(size[0])
            image = numpy.empty(list(map(int,size)))
            image[:] = numpy.NAN
            mapY = 0
            mapX = abs((self.__somHeight-1)/2 - mapY)
            mapX = mapX / 2 + mapX % 2 - 1
            for neuron in self.__neurons:
                mapX = mapX + 1
                off = abs((self.__somHeight-1)/2 - mapY)
                if mapX >= self.__somWidth - math.floor(off / 2) - off % 2 * (mapY) % 2:
                    mapY = mapY + 1
                    mapX = abs((self.__somHeight-1)/2 - mapY)
                    mapX = math.floor(mapX / 2) + mapX % 2 * (1-mapY) % 2
                if mapY >= self.__somHeight:
                    print("abort")
                    return image
                if not shareIntensity:
                    neuron = 1.0 * (neuron - numpy.min(neuron)) / (numpy.max(neuron) - numpy.min(neuron))
                for xPos in range(self.__neuronWidth):
                    for yPos in range(self.__neuronHeight):
                        if math.floor(xPos/2.0 + yPos) < self.__neuronHeight / 4.0 or \
                           math.floor(xPos/2.0 + yPos + 1) > self.__neuronWidth + self.__neuronHeight / 4.0 or \
                           math.floor(xPos/2.0 - yPos + 1) > self.__neuronHeight / 4.0 or \
                           math.floor(xPos/2.0 - yPos) < -self.__neuronWidth + self.__neuronHeight / 4.0:
                            continue
                        else:
                            x = int(math.floor((1.0 * mapX + (mapY%2) / 2.0) * (self.__neuronWidth + border) + xPos + border))
                            y = int(math.floor(1.0 * mapY * (self.__neuronHeight * 3.0 / 4.0 + border) + yPos + border))
                            image[x][y] = neuron[xPos + yPos * self.__neuronWidth]
            return image

    def showMap(self):
        figure = pyplot.figure("SOM", figsize=(16, 16))
        figure.patch.set_alpha(0.0)
        figure.patch.set_facecolor('#ffaadd') 
        pyplot.subplots_adjust(left = 0.01, right = 0.99, bottom = 0.01, top = 0.99, wspace = 0.1, hspace = 0.1)
        if self.isHexMap():
            print ("hexagonal map")
            image = self.calculateMap(self.__neurons, border=2, shareIntensity=False, shape="hex")
        else:
            print ("quadratic map")
            image = self.calculateMap(self.__neurons, border=2, shareIntensity=False, shape="box")
        ax = pyplot.subplot()
        cmap = matplotlib.cm.jet
        cmap.set_bad('#ffaadd',1.)
        ax.imshow(image.T, aspect='auto', interpolation="nearest", cmap=cmap)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis("off")
        pyplot.show()

    def printMap(self):
        figure = pyplot.figure("SOM", figsize=(16, 16))
        figure.patch.set_alpha(0.0)
        figure.patch.set_facecolor('#ffaadd') 
        pyplot.subplots_adjust(left = 0.01, right = 0.99, bottom = 0.01, top = 0.99, wspace = 0.1, hspace = 0.1)
        if self.isHexMap():
            image = self.calculateMap(self.__neurons, border=2, shareIntensity=False, shape="hex")
        else:
            image = self.calculateMap(self.__neurons, border=2, shareIntensity=False, shape="box")
        ax = pyplot.subplot()
        cmap = matplotlib.cm.jet
        cmap.set_bad('#ffaadd',1.)
        ax.imshow(image.T, aspect='auto', interpolation="nearest", cmap=cmap)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis("off")
        pyplot.savefig(self.__fileName + '.pdf', dpi=300)
        pyplot.savefig(self.__fileName + '.png', dpi=300)
        print('Figure saved to pdf and png formats at 300 dpi.')

#if len(sys.argv) < 2:
#    raise Exception('needs a <result-image>.bin as input')
#filename = sys.argv[1]

#myVisualizer = MAPVisualizer(filename) #"resultBootes_10x1x128_1.4_0.05_max.bin")
#myVisualizer = MAPVisualizer("../Results/resultHetdex_10x10x128_1.4_0.05_max.bin")
#myVisualizer = MAPVisualizer("../Results/resultHetdex_hex21x21x128_1.4_0.05_max.bin")
# myVisualizer = MAPVisualizer("../Results/resultHetdex_10x10x128_1.4_0.05_std.bin")
# myVisualizer = MAPVisualizer("../Results/resultBootes_hex11x11x128_1.4_0.05_max.bin")
#myVisualizer.readMap()
#myVisualizer.showMap()
