#
#   Philatelist (Stamp Collector)
#   ______________________________________________________________
#
#   -Go through a catalogue and seek useful sources.
#   -Collect any filenames for matching images and add to a list.
#   -Generate treated stamps from the image files of certain size.
#   -Take the stamps and collate them into a binary file for PInK. 
#

#
#   Imports
#   _______
#

from __future__ import print_function
import numpy as np
import astropy
from astropy.io import fits as fits
from astropy import units as u
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.visualization as vis 
from astLib.astCoords import hms2decimal,dms2decimal
import os
import os.path
import warnings
import struct
import glob
import csv




#
#   Important Variables
#   ___________________
#

minPeak = 0.0                                       #Exclude lower peaks
maxPeak = 1000000                                   #Exclude higher peaks
clip_threshold = 3.0                                # No. times signal must exc. std. dev. for admission.

binFile = 'output2/stamps2.bin'                     #Binary filename
catFile = 'catalog_13jun05.bin'                     #File location of catalogue
FIRST_dirs = '/data/users/deane/vla_first/data/'    #File location of FIRST directories
album = 'stamps/'                                   #File location of stamps
resultingCat = 'output2/stampCat.csv'               #Catalogue of training sources
errList = 'output2/files_with_problems.txt'         #List of corrupt files

numStamps = 150000                                  #Maximum number of training stamps allowed
stampSize = 128                                     #Size of the stamp in pixels




#
#   Go through a catalogue and seek useful sources.
#   _______________________________________________
#

print('Obtaining useful sources from catalogue...')

#Get relevant data from the catalogue:
cat = np.loadtxt(catFile, usecols=[0,1,2,3,4,5,7],skiprows=2)

#Get the R.A. in degrees for each source
h,m,s=cat[:,0],cat[:,1],cat[:,2]
RAdeg=[]
for i in range(len(h)):
    RAdeg.append( hms2decimal("%s"%(h[i])+" "+"%s"%(m[i])+" "+"%s"%(s[i]),' ' ))

#Get the Dec. in degrees for each source
d,m,s=cat[:,3],cat[:,4],cat[:,5]
DECdeg=[]
for i in range(len(h)):
    DECdeg.append( dms2decimal("%s"%(d[i])+" "+"%s"%(m[i])+" "+"%s"%(s[i]),' ' ))

#Get the the peak and coordinates for each source as np arrays
Fpeak = cat[:,6]
RAdeg = np.array(RAdeg)
DECdeg = np.array(DECdeg)

#Get initial list of possible filenames
print('Getting initial list of possible associated filenames...')

#The below process should be handled more elegantly, as above.
#Might rewrite this.
filenames = []
f = open(catFile,'r')
for line in f:
    if (line[0] == '#'): continue
    else:
        splitline = line.split()
        filenames.append(FIRST_dirs + splitline[16][:5] + '/' + splitline[16]+'.fits')

#Filter suitable sources
ind = np.where((Fpeak > minPeak) & (Fpeak < maxPeak) == True)[0]
filenames_subset = [filenames[i] for i in ind]

subcat = np.array([RAdeg[ind],DECdeg[ind],Fpeak[ind],ind])
RAdeg = RAdeg[ind]
DECdeg = DECdeg[ind]

#Remove the below after testing
#print('\nREQUIRED FILES:\n')
#for fn in range(len(filenames_subset)):
#    print(filenames_subset[fn])
#print('\nEND OF FILES\n')

print('Ready to collect files...')





#
#   Collect any filenames for matching images and add to a list.
#   ____________________________________________________________
#

files_on_hand = []
ra_foh = []
dec_foh = []

ns = 0

haveFile = 0
altFile = 0
noFile = 0

numFilesWithBuffer = 1.2*numStamps

for obj in range(len(ind)):

    
    #Check if it exists.
    try:
        if os.path.isfile(filenames_subset[obj]):
            
            files_on_hand.append(filenames_subset[obj])
            ra_foh.append(RAdeg[obj])
            dec_foh.append(DECdeg[obj])
            
            ns += 1
            haveFile += 1
            #If it does not exist, seek an alternative:
        else:
            alternate = glob.glob(filenames_subset[obj][:-6] + '*.fits')
            if (len(alternate) > 0 and os.path.isfile(alternate[-1])):
                
                files_on_hand.append(alternate[-1])
                ra_foh.append(RAdeg[obj])
                dec_foh.append(DECdeg[obj])
            
                altFile += 1
                ns += 1
            else:
                noFile += 1
    except IOError:
        print(' FILE ACCESS ERROR: ' + filenames_subset[obj])
	
    if ( ns == numFilesWithBuffer ): break   #Stop when we reach the max required.
    #If no suitable file was found or an error occurred, move to the next image.

files_on_hand = np.array(files_on_hand)
ra_foh = np.array(ra_foh)
dec_foh = np.array(dec_foh)

print('Of ' + str(ns+noFile) + " files processed: ")
print('ON HAND: ' + str(haveFile) + "\t" + str(100.0*haveFile/(ns+noFile)) + '%')
print('ALTERNATE: ' + str(altFile) + "\t" + str(100.0*altFile/(ns+noFile)) + '%')
print('NONE: ' + str(noFile) + "\t" + str(100.0*noFile/(ns+noFile)) + '%')
    
    
#
#   Generate treated stamps from the image files of certain size.
#   Take the stamps and collate them into a binary file for PInK. 
#   _____________________________________________________________
#

bin_file_path = binFile
size = (stampSize,stampSize)                # Pixel dimensions of each cutout
np.random.seed(7)

#catalog = np.genfromtxt(resultingCat, delimiter = ',')
warnings.filterwarnings('ignore')
#myHiPSfs = HiPSfs(hips_file, noCache=True)

output = open(bin_file_path, 'wb')                  # output file opened for byte writing
output.write(struct.pack('i', numStamps))           # number of objects
output.write(struct.pack('i', 1))                   # number of channels
output.write(struct.pack('i', stampSize) )          # width
output.write(struct.pack('i', stampSize) )          # height


len_catalog = len(files_on_hand)
count = 0
info = []

numErrs = 0
errFiles = []

print('Number of sources: ' + str(len_catalog))

for i in range(len(files_on_hand)):
    
    if count == numStamps: break #Stop when we have enough
    
    # Notify when every N objects has been completed
    if i % 1000 == 0:
        print( "\n" + str(i) + '/' + str(len_catalog))
        
    #Generate cutout
    
    
    try: 
        ra, dec = (ra_foh[i], dec_foh[i]) #Coordinates of this source
    
        fileData = fits.open(files_on_hand[i])
        hdr = fileData[0].header
    
        hdr.remove('CRPIX3')
        hdr.remove('CRVAL3')
        hdr.remove('CDELT3')
        hdr.remove('CTYPE3')
        hdr.remove('CROTA3')
        hdr.remove('NAXIS3')
    
        hdr.remove('CRPIX4')
        hdr.remove('CRVAL4')
        hdr.remove('CDELT4')
        hdr.remove('CTYPE4')
        hdr.remove('CROTA4')
        hdr.remove('NAXIS4')
        
        hdr['NAXIS'] = 2
        
        w = WCS( header = hdr ) #WCS object for this source
        
        x, y = w.all_world2pix(ra, dec, 1, adaptive=True) #Pixel coordinates of target
        
        position = (x, y)
        imData = fileData[0].data[0][0] #Shape it correctly. Assumes 4D. Might rewrite to generalise.
        cutout = Cutout2D(imData, position, size, mode='partial') #Anything past the border filled with NaN
        
        image = np.array(cutout.data) #Get the image data as a np array.
        
        empty = np.isnan(image) #Get all NaN pixels and replace with normalised noise.
        image[empty] = np.random.normal(loc=np.mean(image[~empty]), scale=np.std(image[~empty]), size=image.shape)[empty] # replace missing values with noise  
        image = image - np.min(image) #Set base intensity to 0
        
        # clip background
        image = np.clip(image, clip_threshold*np.std(image), 1e10)
        
        # Normalize using minmax
        mminterval = vis.MinMaxInterval()
        image = mminterval(image)
        
        image.astype('f').tofile(output) #Write to the Binary file
        
        #Append the info to our catalogue
        entry = []
        entry.append(str(ra))
        entry.append(str(dec))
        entry.append(files_on_hand[i])
        file = album + str(count) + '-' + str(int(round(ra))) + ':' + str(int(round(dec))) + '.stamp.fits'
        entry.append(file)
        info.append(entry)
        
        #Create a FITS image of the stamp
        hdu = fits.PrimaryHDU()
        hdu.data = image
        hdu.header['CRVAL1'] = ra
        hdu.header['CRVAL2'] = dec
        
        hdu.writeto(file, clobber=True)
        
        count += 1
        
    except:
        numErrs += 1
        errFiles.append(files_on_hand[i])
        
output.close()

print('Number of images in binary:', count, ' of ', len_catalog)
info = np.array(info)

print("Number of files with errors: " + str(numErrs))
errInfo = np.array(errFiles)
np.savetxt(errList, errInfo, fmt="%s", delimiter='\n', newline='\n')

np.savetxt(resultingCat, info, fmt="%s", delimiter=',', newline='\n', header='RA,DEC,Origin File,Stamp File', comments='# ')

