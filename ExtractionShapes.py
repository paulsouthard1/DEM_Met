# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 12:58:02 2018

@author: ps29626
"""

import numpy as np
import os
import pickle
import shutil
import argparse
from cmath import rect, phase
from math import radians, degrees

parser = argparse.ArgumentParser(description='Produce some vertical lines along a flowpath to sample channel characteristics longitudinally')
parser.add_argument('folder', help='Folder to store resulting files in A.K.A. name describing analysis')
parser.add_argument('demin', help='DEM of Region')
parser.add_argument('reachin', help='X,Y Coordinates of reach')
parser.add_argument('--spacing',  type=int, nargs='?', default=20, help='Spacing between lines, length of sampling window')
parser.add_argument('--width',  type=int, nargs='?', default=150, help='Width of lines')
args = parser.parse_args()


#Define neighbor function
def neighbors(im,i,j,d=1):
	b = im[i-d:i+d+1, j-d:j+d+1].flatten()
	# remove the element (i,j)
	n = np.hstack((b[:len(b)//2],b[len(b)//2+1:] ))
	return n

#Define mean angles function
from math import sin,cos,atan2,pi
def meanangle(angles,weights=0,setting='degrees'):
	'''computes the mean angle'''
	if weights==0:
		 weights=np.ones(len(angles))
	sumsin=0
	sumcos=0
	if setting=='degrees':
		angles=np.array(angles)*pi/180
	for i in range(len(angles)):
		sumsin+=weights[i]/sum(weights)*sin(angles[i])
		sumcos+=weights[i]/sum(weights)*cos(angles[i])
	average=atan2(sumsin,sumcos)
	if setting=='degrees':
		average=average*180/pi
	return average


def mean_angle(deg):
	return degrees(phase(sum(rect(1, radians(d)) for d in deg)/len(deg)))

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
	with open(DEMA,'r') as f:
		print("Reading header info from " + DEMA)
		#Read in number of columns
		ncols = str(f.readlines(1))
		ncols = ''.join(filter(str.isdigit, ncols))
		#Read in number of rows
		nrows = str(f.readlines(2))
		nrows = ''.join(filter(str.isdigit, nrows))
		#Read in X coordinate of lower left corner
		xllcorner = str(f.readlines(3))
		xllcorner = ''.join(filter(str.isdigit, xllcorner))
		xll = float(xllcorner[:6]+"."+xllcorner[6:])
		#Read in y coordinate of lower left corner
		yllcorner = str(f.readlines(4))
		yllcorner = ''.join(filter(str.isdigit, yllcorner))
		yll = float(yllcorner[:7] + "." + yllcorner[7:])
		#Read in cellsize
		cellsize = str(f.readlines(5))
		cellsize = int(''.join(filter(str.isdigit, cellsize)))
		#Read in cell value for NoData cells
		NODATA = str(f.readlines(6))
		NODATA = ''.join(filter(str.isdigit, NODATA))
		#Find remainders
		xllr = xll%1
		yllr = yll%1
		#Find y coordinate lower left corner of upper left cell for indexing cells
		yul = yll + (float(nrows)-1)
	return ncols,nrows,xll,yll,cellsize,NODATA,xllr,yllr,yul


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

def pointlist(dir_path,reachlist,FAArray,xll,yll,yul,xllr,yllr,cellsize,spacing,width,storepath):
	# Create blank dictionary to store profile points along flow accumulation line
	points = {}
	# Create blank dictionary to store profile indices along flow accumulation line
	inds = {}
	#Define exact coordinates of start and end cell corners
	x1 = float(reachlist[0:6])+xllr
	y1 = float(reachlist[7:14])+yllr
	x2 = float(reachlist[15:21])+xllr
	y2 = float(reachlist[22:29])+yllr
	print("Creating cross-sections for reach ("+str(x1)+", "+str(y1)+") to ("+str(x2)+", "+str(y2)+")")
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
	#Create max value variable so that FA stepthrough works
	ma0 = 10**20
	# Create blank variables to store k and l coordinates for each reach
	klist = np.zeros(1000000)
	llist = np.zeros(1000000)
	#Run loop that works through reach using flow accumulation array and selects cell with maximum accumulation
	m = 0
	while ((k!=k2) == True) & ((l != l2) == True):
		#Store k's and l's in variable
		klist[m] = k
		llist[m] = l
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
			#Run Neighbor function again to get last point
			Temp = neighbors(FAArray,k,l,d=1)
			j = np.argmax(Temp)
			ma1 = Temp[j]
			Temp2 = np.delete(Temp,j)
			ma2 = np.amax(Temp2)
			j2 = np.where(Temp==ma2)
			j2 = j2[0][0]
			print(str(ma0)+" "+str(ma1)+" "+str(ma2))
			if (ma2>ma0) & (ma1>ma2):
				print("a")
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
			longthresh = longthresh + Long
			ma0 = ma1
			m = m+1
		else:
			Temp = neighbors(FAArray,k,l,d=1)
			j = np.argmax(Temp)
			ma1 = Temp[j]
			Temp2 = np.delete(Temp,j)
			ma2 = np.amax(Temp2)
			j2 = np.where(Temp==ma2)
			j2 = j2[0][0]
			print(str(ma0)+" "+str(ma1)+" "+str(ma2))
			if (ma2>ma0) & (ma1>ma2):
				print("a")
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
			ma0 = ma1
			longthresh = longthresh + Long
	return seg,points,klist,llist,inds
			
def storemeta(dir_path,storepath,folder,demin,spacing,width,seg):
	print("Storing Metadata")
	#Create blank processnotes dictionary
	sect = {}
	#Store number of files for export as pickle file
	sect["NumSections"] = seg
	#Save spacing in pickle file
	sect["Spacing"] = spacing
	sect["Width"] = width
	sect["Folder"] = folder
	sect["Reach"] = "reach_"+folder+".txt"
	file = open(dir_path + "sectnum", "wb")
	pickle.dump(sect,file,protocol=2)
	file.close()
	shutil.move(dir_path+"sectnum",storepath+"sectnum")
	
def makesect(storepath,seg,inds,klist,llist,points,width):
	angles = {}
	for m in np.arange(2, seg):
		# Call index along flow accumulation line where a given cross-section was stored
		step1 = inds["seg_{}".format(m-1)]
		step1 = [llist[step1],klist[step1]]
		step2 = inds["seg_{}".format(m)]
		step2 = [llist[step2],klist[step2]]
		step3 = inds["seg_{}".format(m+1)]
		step3 = [llist[step3],klist[step3]]
		theta1 = (np.pi/2)-np.arctan2(step2[1]-step1[1], step2[0]-step1[0])
		theta2 = (np.pi/2)-np.arctan2(step3[1]-step2[1], step3[0]-step2[0])
		theta = [theta1,theta2]
		theta = meanangle(theta)
		print(str(theta1)+", "+str(theta2)+" -> "+str(theta))
		angles["sct_{}".format(m)] = [theta]
		temppoints = points["seg_{}".format(m)]
		#Calculate endpoints of line centered on the flow accumulation line
		RBX = temppoints[0]+(np.cos(theta)*(width/2))
		RBY = temppoints[1]+(np.sin(theta)*(width/2))
		LBX = temppoints[0]-(np.cos(theta)*(width/2))
		LBY = temppoints[1]-(np.sin(theta)*(width/2))
		# Save endpoints to a textfile
		savename = storepath + "CS_{}".format(m)+".txt"
		print("Theta of " + str(theta) + "for " + savename)
		f = open(savename, 'w')
		f.write("RBX, RBY, LBX, LBY" + "\n")
		f.write(str(RBX) + ", " + str(RBY) + ", " + str(LBX) + ", " + str(LBY))
		f.close()
	return angles
	
def Main(folder,demin,reachin,spacing,width):
	#Define root directory
	dir_path = os.path.dirname(os.path.realpath(__file__))+"\\"
	#Define Arc File Directory
	arcpath = dir_path + "ArcFiles\\"
	#Define Cross File Directory
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
	ncols,nrows,xll,yll,cellsize,NODATA,xllr,yllr,yul = readheader(DEMA)
	#Load in reach endpoints for analysis
	reachlist = readreach(reachin,storepath,folder)
	#Create list of points on thalweg to sample cross-sections at
	seg,points,klist,llist,inds = pointlist(dir_path,reachlist,FAArray,xll,yll,yul,xllr,yllr,cellsize,spacing,width,storepath)
	storemeta(dir_path,storepath,folder,demin,spacing,width,seg)
	angles = makesect(storepath,seg,inds,klist,llist,points,width)
Main(args.folder,args.demin,args.reachin,args.spacing,args.width)

