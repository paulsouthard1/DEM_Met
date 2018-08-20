# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 16:45:49 2018

@author: ps29626
"""

import os
import numpy as np
import shutil
import pickle
import argparse
import fileinput

#ArgParse stuff
parser = argparse.ArgumentParser(description='Extract longitudinal profile and list of points on flow accumulation line')
parser.add_argument('folder', help='Folder to store resulting files in A.K.A. name describing analysis')
parser.add_argument('demin', help='DEM of Region')
parser.add_argument('reachin', help='X,Y Coordinates of reach')
parser.add_argument('dataset', help='dataset to extract swaths from')
parser.add_argument('--spacing',  type=int, nargs='?', default=50, help='Spacing between lines, length of sampling window')
# parser.add_argument('clip', help='Anuga result used to clip lines')
# parser.add_argument('lines', help='Do you need to make line shapefiles?  "Yes" or "No".')
#parser.add_argument('swaths', type=bool, help='Do you need to make swath shapefiles?')
args = parser.parse_args()



#Define neighbor function
def neighbors(im, i, j, d=1):
	b = im[i-d:i+d+1, j-d:j+d+1].flatten()
	# remove the element (i,j)
	n = np.hstack((b[:len(b)//2],b[len(b)//2+1:] ))
	return n

#Define loadascii function
def loadascii(DEMA,FAA):
	print("Loading ASCII file " + DEMA)
	#Load DEM
	DArray=np.loadtxt(DEMA,skiprows=6)
	#Load FAA
	FAArray=np.loadtxt(FAA,skiprows=6)
	return DArray,FAArray

#Define read header function
def readheader(DEMA):
	#Create blank dictionary to store header info
	Header = {}
	with open(DEMA,'r') as f:
		print("Reading header info from " + DEMA)
		#Read in number of columns
		ncols_f = str(f.readlines(1))
		ncols = ''.join(filter(str.isdigit, ncols_f))
		#Read in number of rows
		nrows_f = str(f.readlines(2))
		nrows = ''.join(filter(str.isdigit, nrows_f))
		#Read in X coordinate of lower left corner
		xllcorner_f = str(f.readlines(3))
		xllcorner = ''.join(filter(str.isdigit, xllcorner_f))
		xll = float(xllcorner[:6]+"."+xllcorner[6:])
		#Read in y coordinate of lower left corner
		yllcorner_f = str(f.readlines(4))
		yllcorner = ''.join(filter(str.isdigit, yllcorner_f))
		yll = float(yllcorner[:7] + "." + yllcorner[7:])
		#Read in cellsize
		cellsize_f = str(f.readlines(5))
		cellsize = int(''.join(filter(str.isdigit, cellsize_f)))
		#Read in cell value for NoData cells
		NODATA_f = str(f.readlines(6))
		NODATA = ''.join(filter(str.isdigit, NODATA_f))
		#Find remainders
		xllr = xll%1
		yllr = yll%1
		#Find y coordinate lower left corner of upper left cell for indexing cells
		yul = yll + (float(nrows)-1)
	#Store header lines in dictionary
	Header[1] = ncols_f
	Header[2] = nrows_f
	Header[3] = xllcorner_f
	Header[4] = yllcorner_f
	Header[5] = cellsize_f
	Header[6] = NODATA_f
	return ncols,nrows,xll,yll,cellsize,NODATA,xllr,yllr,yul,Header


#Define read reach coordinates function
def readreach(reachin,storepath,folder):
	print("Reading reach beginning and terminus from "+reachin)
	#Read in starting point for reach
	with open(reachin,'r') as f:
		#Variable reachlist contains coordinates for each reach to be analyzed
		reachlist = f.readlines()
		reachlist = reachlist[0]
		shutil.copyfile(reachin,storepath+"reach_"+folder+".txt")
	return reachlist

def pulldata(reachlist,FAArray,DArray,xll,yll,yul,xllr,yllr,nrows,ncols,cellsize,spacing):
	# Create blank dictionary to store profile points along flow accumulation line
	points = {}
	# Create blank dictionary to store profile indices along flow accumulation line
	inds = {}
	# Create blank vectors to store dataset values 
	Elevs = np.zeros(1000000)
	Longs = np.zeros(1000000)
	streamline = np.zeros((int(nrows),int(ncols)),dtype=np.int)
	#Define exact coordinates of start and end cell corners
	x1 = float(reachlist[0:6])+xllr
	y1 = float(reachlist[7:14])+yllr
	x2 = float(reachlist[15:21])+xllr
	y2 = float(reachlist[22:29])+yllr
	print("Pulling dataset values for reach ("+str(x1)+", "+str(y1)+") to ("+str(x2)+", "+str(y2)+")")
	#Find starting indices
	l1 = round((x1-xll)/cellsize)
	k1 = round((yul-y1)/cellsize)
	#Define initial k,l
	l = l1
	k = k1
	#Find ending indices
	l2 = round((x2-xll)/cellsize)
	k2 = round((yul-y2)/cellsize)
	#Create zero-value variable to test whether spacing distance has been achieved on each iteration
	longthresh = 0
	#Create zero-value variable to order profile points
	seg = 0
	# Create blank variables to store k and l coordinates for each reach
	klist = np.zeros(1000000)
	llist = np.zeros(1000000)
	ma0 = 10**20
	m = 0
	while ((k!=k2) == True) & ((l != l2) == True):
		#Store k's and l's in variable
		klist[m] = k
		llist[m] = l
		#Store Elevation value
		Elevs[m] = DArray[k,l]
		#Mark cell k,l in streamline for later flowline raster
		streamline[k,l] = m
		#Create a new cross-section point if spacing has been exceeded
		if longthresh >= spacing:
			#Define point number
			seg = seg + 1
			#Compute and store coordinates of center of cell that is 
			tempx = xll + (l * cellsize) + (cellsize/2)
			tempy = yul - (k * cellsize) + (cellsize/2)
			#Store indices where points were selected
			inds["seg_{}".format(seg)] = [m]
			#Store coordinates to base cross-section line
			points["seg_{}".format(seg)] = [tempx,tempy]
			#Reset threshold counter
			longthresh = 0
			#Run Neighbor function again to get last point
			Temp = neighbors(FAArray,k,l,d=1)
			j = np.argmax(Temp)
			ma1 = Temp[j]
			Temp2 = np.delete(Temp,j)
			ma2 = np.amax(Temp2)
			j2 = np.where(Temp==ma2)
			j2 = j2[0][0]
			if (ma2>ma0) & (ma1>ma2):
				j = j2
				ma1 = ma2
			if j == 0:
				k = k-1
				l = l-1
				Long = np.sqrt(2)*cellsize
			elif j == 1:
				k = k-1
				l = l
				Long = cellsize
			elif j == 2:
				k = k-1
				l = l+1
				Long = np.sqrt(2)*cellsize
			elif j == 3:
				k = k
				l = l-1
				Long = cellsize
			elif j == 4:
				k = k
				l = l+1
				Long = cellsize
			elif j == 5:
				k = k+1
				l = l-1
				Long = np.sqrt(2)*cellsize
			elif j == 6:
				k = k+1
				l = l
				Long = cellsize
			else:
				k = k+1
				l = l+1
				Long = np.sqrt(2)*cellsize
			m = m+1
			ma0 = ma1
			longthresh = longthresh + Long
		else:
			#Run Neighbor function again to get last point
			Temp = neighbors(FAArray,k,l,d=1)
			j = np.argmax(Temp)
			ma1 = Temp[j]
			Temp2 = np.delete(Temp,j)
			ma2 = np.amax(Temp2)
			j2 = np.where(Temp==ma2)
			j2 = j2[0][0]
			if (ma2>ma0) & (ma1>ma2):
				j = j2
				ma1 = ma2
			if j == 0:
				k = k-1
				l = l-1
				Long = np.sqrt(2)*cellsize
			elif j == 1:
				k = k-1
				l = l
				Long = cellsize
			elif j == 2:
				k = k-1
				l = l+1
				Long = np.sqrt(2)*cellsize
			elif j == 3:
				k = k
				l = l-1
				Long = cellsize
			elif j == 4:
				k = k
				l = l+1
				Long = cellsize
			elif j == 5:
				k = k+1
				l = l-1
				Long = np.sqrt(2)*cellsize
			elif j == 6:
				k = k+1
				l = l
				Long = cellsize
			else:
				k = k+1
				l = l+1
				Long = np.sqrt(2)*cellsize
			# Keep track of how far away the last cross-section is
			m = m+1
			Longs[m] = Longs[m-1] + Long
			ma0 = ma1
			longthresh = longthresh + Long
	Elevs = Elevs[0:m]
	Longs = Longs[0:m]
	streamline.astype(int)
	return points,inds,Elevs,Longs,streamline,seg

def StreamRast(dir_path,storepath,streamline,Header):
	# Create text file of header for streamline raster with write function
	rastsave = "Header.txt"
	f = open(rastsave, 'w')
	for j in range(6):
		htext = Header[j+1]
		htext = htext[2:-4]
		f.write(htext+"\n")
	f.close()
	# Create text file with streamline array using np.savetxt
	rastsave = "Raster.txt"
	np.savetxt(rastsave,streamline,fmt='%.1i' ,delimiter=' ', newline='\n', header="")
	#Concatenate header and array files
	outfilename = "Streamline.txt"
	headername = "Header.txt"
	rastername = "Raster.txt" 
	with open(outfilename, 'w') as fout:
		fin = fileinput.input(files=(headername, rastername))
		for line in fin:
			fout.write(line)
		fin.close()
	os.remove(headername)
	os.remove(rastername)
	shutil.move(dir_path+"Streamline.txt",storepath+"streamline.txt")

def storemeta(dir_path,storepath,folder,demin,spacing,seg):
	print("Storing Metadata")
	#Create blank processnotes dictionary
	sect = {}
	#Store number of files for export as pickle file
	sect["NumSections"] = seg
	#Save spacing in pickle file
	sect["Spacing"] = spacing
	sect["Folder"] = folder
	sect["Reach"] = "reach_"+folder+".txt"
	file = open(dir_path + "sectnum", "wb")
	pickle.dump(sect,file,protocol=2)
	file.close()
	shutil.move(dir_path+"sectnum",storepath+"sectnum")

def storedata(dir_path,storepath,points,inds,Elevs,Longs,streamline):
	All = {}
	All["points"] = points
	All["inds"] = inds
	All["Elevs"] = Elevs
	All["Longs"] = Longs
	file = open(dir_path + "slope", "wb")
	pickle.dump(points,file,protocol=2)
	file.close()
	shutil.move(dir_path+"slope",storepath+"slope")

def Main(folder,demin,reachin,dataset,spacing):
	dir_path = os.path.dirname(os.path.realpath(__file__)) + "\\"
	arcpath = dir_path + "ArcFiles\\"
	crosspath = dir_path + "CrossFiles\\"
	if not os.path.exists(crosspath):
		os.makedirs(crosspath)
	#Define directory for results storage
	storepath = crosspath + folder + "\\"
	if not os.path.exists(storepath):
		os.makedirs(storepath)
	#Pull out filename
	filein = os.path.basename(demin)
	filein = filein[:filein.find(".")]
	#Define names for ASCII files based on DEM file name
	DEMA = arcpath + filein + "_demasc.txt"
	FAA = arcpath + filein + "_faasc.txt"
	#Load in ASCII's as Arrays
	DArray,FAArray = loadascii(DEMA,FAA)
	#Load in Header information from DEM ASCII
	ncols,nrows,xll,yll,cellsize,NODATA,xllr,yllr,yul,Header = readheader(DEMA)
	#Load in reach endpoints for analysis
	reachlist = readreach(reachin,storepath,folder)
	points,inds,Elevs,Longs,streamline,seg = pulldata(reachlist,FAArray,DArray,xll,yll,yul,xllr,yllr,nrows,ncols,cellsize,spacing)
	storedata(dir_path,storepath,points,inds,Elevs,Longs,streamline)
	storemeta(dir_path,storepath,folder,demin,spacing,seg)
	StreamRast(dir_path,storepath,streamline,Header)

Main(args.folder,args.demin,args.reachin,args.dataset,args.spacing)