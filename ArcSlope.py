
# coding: utf-8

# In[58]:

import numpy as np
import matplotlib.pyplot as plt
import os
#Define neighbor function
def neighbors(im, i, j, d=1):
    b = im[i-d:i+d+1, j-d:j+d+1].flatten()
    # remove the element (i,j)
    n = np.hstack((b[:len(b)//2],b[len(b)//2+1:] ))
    return n

direc = os.path.dirname(os.path.realpath(__file__)) + "\\"


# In[17]:

#Define file names
filein = "BE_08"
DEMA = direc + filein + "_demasc.txt"
FAA = direc + filein + "_faasc.txt"


# In[18]:

#Load DEM
DArray=np.loadtxt(DEMA,skiprows=6)
#Load FAA
FAArray=np.loadtxt(FAA,skiprows=6)


# In[19]:

#Create blank dictionary to store header info
Header = {}
#Read in pertinent header info
with open(DEMA,'r') as f:
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


# In[20]:

#Read in endpoints for reaches
ReachFile = direc + "reach.txt"
with open(ReachFile,'r') as f:
#Variable reachlist contains coordinate pairs for each reach
    reachlist = f.readlines()
#Create an array of zeroes based on dimensions of input raster
dim = (int(nrows),int(ncols))


# In[21]:

# Create blank dictionaries
elevations = {}
longs = {}
title = {}
reachrast = {}
#Define points along reach for analysis, store them
for i in range(len(reachlist)):
    tempreach = reachlist[i]
    # Define exact coordinates of start and end cell corners, and then create title for plots
    x1 = float(tempreach[0:6])+xllr
    y1 = float(tempreach[7:14])+yllr
    x2 = float(tempreach[15:21])+xllr
    y2 = float(tempreach[22:29])+yllr
    title[i] = "Slope: Reach " + str(i+1)
    #Find starting indices
    l = round((x1-xll)/cellsize)
    k = round((yul-y1)/cellsize)
    #Find ending indices
    l2 = round((x2-xll)/cellsize)
    k2 = round((yul-y2)/cellsize)
    #Define blank gradient vectors
    Elev = np.zeros(1000000)
    Long = np.zeros(1000000)
    streamline = np.zeros(dim,dtype=np.int)
    #Run loop that works through reach using flow accumulation array and selects cell with maximum accumulation
    for m in range(0,999999):
        #Store elevation at k,l in dictionary
        Elev[m] = DArray[k,l]
        #Mark cell k,l in streamline for later flowline raster
        streamline[k,l] = 1
        streamline.astype(int)
        #Break loop if at end of reach
        if (k == k2) & (l == l2):
            break
        #If not at end, continue to work through neighborhoods looking for max accumulation
        else:
            Temp = neighbors(FAArray,k,l, d=1)
            j = np.argmax(Temp)
            if j == 0:
                k = k-1
                l = l-1
                Long[m+1] = Long[m]+np.sqrt(2)*cellsize
            elif j == 1:
                k = k-1
                l = l
                Long[m+1] = Long[m]+cellsize
            elif j == 2:
                k = k-1
                l = l+1
                Long[m+1] = Long[m]+np.sqrt(2)*cellsize
            elif j == 3:
                k = k
                l = l-1
                Long[m+1] = Long[m]+cellsize
            elif j == 4:
                k = k
                l = l+1
                Long[m+1] = Long[m]+cellsize
            elif j == 5:
                k = k+1
                l = l-1
                Long[m+1] = Long[m]+np.sqrt(2)*cellsize
            elif j == 6:
                k = k+1
                l = l
                Long[m+1] = Long[m]+cellsize
            else:
                k = k+1
                l = l+1
                Long[m+1] = Long[m]+np.sqrt(2)*cellsize
    #Store elevation, longitudinal distance and flowline arrays for each reach as individual arrays in dictionary
    elevations["elev_{}".format(i)] = Elev[0:m]
    longs["long_{}".format(i)] = Long[0:m]
    reachrast["reachrast_{}".format(i)] = streamline


# In[60]:

# Plot slope of each reach
for i in range(len(reachlist)):
    x = "long_{}".format(i)
    y = "elev_{}".format(i)
    t = title[i]
    savename = t[7:-2]+"_"+t[-1]+".png"
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(longs[x], elevations[y])
    plt.title(t)
    plt.ylabel("Elevation (m)")
    plt.xlabel("Longitudinal Distance (m)")
    plt.savefig(savename,bbox_inches = 'tight')


# In[61]:

# Create text file of header for streamline raster with .write function
for i in range(len(reachlist)):
    rastsave = "Header_" + str(i) + ".txt"
    f = open(rastsave, 'w')
    for j in range(6):
        htext = Header[j+1]
        htext = htext[2:-4]
        f.write(htext+"\n")
    f.close()


# In[62]:

# Create text file with streamline array using np.savetxt
for i in range(len(reachlist)):
    rastsave = "Raster_" + str(i) + ".txt"
    rastcall = "reachrast_{}".format(i)
    np.savetxt(rastsave, reachrast[rastcall],fmt='%.1i' ,delimiter=' ', newline='\n', header="")


# In[63]:

#Concatenate header and array files
import fileinput
for i in range(len(reachlist)):
    outfilename = "Streamline_" + str(i) + ".txt"
    headername = "Header_" + str(i) + ".txt"
    rastername = "Raster_" + str(i) + ".txt"
    with open(outfilename, 'w') as fout:
        fin = fileinput.input(files=(headername, rastername))
        for line in fin:
            fout.write(line)
        fin.close()
    os.remove(headername)
    os.remove(rastername)

