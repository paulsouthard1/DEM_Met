import numpy as np
import os
import pickle
import fiona
import argparse
from bokeh.palettes import RdGy as palette
import matplotlib.pyplot as plt
from matplotlib import rcParams
from shapely import geometry
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma']
rcParams['font.size'] = 20
rcParams['font.style'] = 'italic'
rcParams['axes.labelsize'] = 10
rcParams['axes.labelweight'] = 'medium'
rcParams['xtick.labelsize'] = 10
rcParams['ytick.labelsize'] = 10


def ReadData(lookuppath):
    #Read in pickle file with number of sections
    print("Loading Metadata Pickle File")
    sect = pickle.load(open(lookuppath + "sectnum","rb"))
    #Assign number of sections to variable
    NumSect = sect["NumSections"]
    return NumSect

def GeomCalc(NumSect,lookuppath,dataname,sectionfolder):
    #Make blank dictionaries
    CS = {}
    CSA = {}
    WP = {}
    HR = {}
    LW = {}
    #Perform calculations
    for i in np.arange(2,NumSect):
        #Define dictionary name, look and store locations
        keyname = dataname + "_{}".format(i)
        savename = keyname + ".txt"
        lookupname = lookuppath + savename
        length = 0.0
        sectname = "cl_CS_{}".format(i) + ".shp"
        sectionlookup = sectionfolder + sectname
        #Use fiona module to get shapefile length
        with fiona.open(sectionlookup, 'r') as shp:
            for line in shp:
                length += geometry.shape(line['geometry']).length
        LW[keyname] = length
        CSTemp = np.loadtxt(lookupname,delimiter = ',',skiprows = 1,usecols = (1,2))
        CS[keyname] = CSTemp
        area = 0
        perimeter = 0
        if len(CSTemp) > 0:
            #Calculate CSA
            dist = CSTemp[:,0]
            depth = CSTemp[:,1]
            area = np.zeros(len(dist)-1)
            perimeter = np.zeros(len(dist)-1)
            lb = depth[0]
            rb = depth[-1]
            for j in range(len(dist)-1):
                length = (dist[j+1]-dist[j])
                depthavg = (depth[j+1]+depth[j])/2
                area[j] = depthavg * length
                lb = depth[0]
                rb = depth[-1]
                perimeter[j] = np.sqrt(((depth[j+1]-depth[j])**2+(dist[j+1]-dist[j])**2))
            CSA[keyname] = np.sum(area)
            WP[keyname] = lb + rb + np.sum(perimeter)
            HR[keyname] = (CSA[keyname]/WP[keyname])
        else:
            CSA[keyname] = 0
            WP[keyname] = 0
            HR[keyname] = 0
            print(keyname + " does not intersect target dataset")
    return CS,CSA,WP,HR,LW

def GeomArrange(NumSect,CS,CSA,WP,HR,LW,dataname):
    CSA_a = np.zeros(len(np.arange(2,NumSect)))
    WP_a = np.zeros(len(np.arange(2,NumSect)))
    HR_a = np.zeros(len(np.arange(2,NumSect)))
    LW_a = np.zeros(len(np.arange(2,NumSect)))
    keynum = np.zeros(len(np.arange(2,NumSect)))
    sects = np.arange(2,NumSect)
    for i in np.arange(2,NumSect):
        keyname = dataname + "_{}".format(i)
        CSA_a[i-2] = CSA[keyname]
        WP_a[i-2] = WP[keyname]
        HR_a[i-2] = HR[keyname]
        LW_a[i-2] = LW[keyname]
    return CSA_a,WP_a,HR_a,LW_a,sects

def MovingAvg(CSA_a,WP_a,HR_a,LW_a,wind):
    look = (wind-1)/2
    datalist = [CSA_a,WP_a,HR_a,LW_a]
    CSA_S = np.zeros(len(CSA_a))
    WP_S = np.zeros(len(WP_a))
    HR_S = np.zeros(len(HR_a))
    LW_S = np.zeros(len(LW_a))
    newdata = [CSA_S,WP_a,HR_a,LW_a]
    for i in range(len(datalist)):
        for j in range(len(datalist[i])):
            back = int(j-look)
            if back < 0:
                back = 0
            front = int(j+look)
            if front > len(datalist[i]):
                front = len(datalist[i])
            reach = datalist[i][back:front]
            val = np.average(reach)
            newdata[i][j] = val
    CSA_S = newdata[0]
    WP_S = newdata[1]
    HR_S = newdata[2]
    LW_S = newdata[3]
    return CSA_S,WP_S,HR_S,LW_S


def PlotGeom(CSA_a,WP_a,HR_a,LW_a,sects,NumSect,show,fileout):
    #Plot Results
    #Define colorlist
    colorlist = palette[10]
    #Create fig and define axes
    fig = plt.subplots(nrows = 2, ncols = 2,figsize = (20,10))
    ax1 = plt.subplot(221)
    ax2 = plt.subplot(222)
    ax3 = plt.subplot(223)
    ax4 = plt.subplot(224)
    #Set listnum to 0 for colorlist
    listnum = 0
    #Put Section markers on graphs
    for i in np.arange(0,len(CSA_a),show):
        if listnum == 10:
            listnum = 0
        else:
            listnum = listnum
        col = colorlist[listnum]
        ax1.plot(sects[i], HR_a[i], marker = 's',color=col, mec = 'k', markersize = 10)
        ax2.plot(sects[i], CSA_a[i], marker = 's',color=col, mec = 'k', markersize = 10)
        ax3.plot(sects[i], WP_a[i], marker = 's',color=col, mec = 'k', markersize = 10)
        ax4.plot(sects[i], LW_a[i], marker = 's',color=col, mec = 'k', markersize = 10)
        listnum = listnum + 1
    #Plot all data
    ax1.plot(sects, HR_a, color = '#206f89',linewidth = 2, alpha = 0.9)
    ax4.plot(sects, LW_a, color = '#893a20',linewidth = 2, alpha = 0.9)
    ax2.plot(sects, CSA_a, color = '#206f89',linewidth = 2, alpha = 0.9)
    ax3.plot(sects, WP_a, color = '#893a20',linewidth = 2, alpha = 0.9)
    #Set axis limits
    ax1.set_ylim(0,np.max(HR_a))
    ax2.set_ylim(0,np.max(CSA_a))
    ax3.set_ylim(0,np.max(WP_a))
    ax4.set_ylim(0,np.max(LW_a))
    ax1.set_xlim(0,NumSect-2)
    ax2.set_xlim(0,NumSect-2)
    ax3.set_xlim(0,NumSect-2)
    ax4.set_xlim(0,NumSect-2)
    #Make grid
    ax1.grid()
    ax2.grid()
    ax3.grid()
    ax4.grid()
    #Label graph and save
    ax1.set_ylabel("Hydraulic Radius $(m)$", size = 20, color = '#206f89', alpha = 0.9)
    ax2.set_ylabel("Cross-Sectional Area $(m^2)$", size = 20, color = '#206f89', alpha = 0.9)
    ax3.set_ylabel("Wetted Perimeter $(m)$", size = 20, color = '#893a20', alpha = 0.9)
    ax4.set_ylabel("Cross-section Length $(m)$", size = 20, color = '#893a20', alpha = 0.9)
    ax4.set_xlabel("Cross-Section Number", size = 20, color = 'k')
    ax2.set_xlabel("Cross-Section Number", size = 20, color = 'k')
    ax1.set_xlabel("Cross-Section Number", size = 20, color = 'k')
    ax3.set_xlabel("Cross-Section Number", size = 20, color = 'k')
    plt.suptitle("Hydraulic Measurements for Modeled 20 CMS discharge in Woodruff Canyon",fontsize = 25)
    plt.savefig(fileout,bbox_inches = 'tight')

#ArgParse stuff
parser = argparse.ArgumentParser(description='Calculate and plot metrics of hydraulic geometry along the reach')
parser.add_argument('folder', help='Folder to store resulting files in A.K.A. name describing analysis')
parser.add_argument('dataset', help='depth dataset')
parser.add_argument('show', help='Show every _th cross-section', type=int)
parser.add_argument('fileout', help='File path to store plot')
parser.add_argument('smooth', help='Smooth data? Yes or No')
parser.add_argument('wind', help='Smoothing window to work over',type = int)
args = parser.parse_args()


def Main(folder,dataset,show,fileout,smooth,wind):
    #Define necessary pathnames
    dir_path = os.path.dirname(os.path.realpath(__file__))+"\\"
    arcpath = dir_path + "ArcFiles\\"
    crosspath = dir_path + "CrossFiles\\"
    lookuppath = crosspath + folder + "\\"
    sectionfolder = arcpath + folder + "\\"
    dataname = os.path.basename(dataset)
    dataname = dataname[:dataname.find(".")]
    #Run function to find number of sections
    NumSect = ReadData(lookuppath)
    #Run function to calculate geometries
    CS,CSA,WP,HR,LW = GeomCalc(NumSect,lookuppath,dataname,sectionfolder)
    CSA_a,WP_a,HR_a,LW_a,sects = GeomArrange(NumSect,CS,CSA,WP,HR,LW,dataname)
    if smooth == "No":
        #Run function to make 4X4 plot of geometry
        PlotGeom(CSA_a,WP_a,HR_a,LW_a,sects,NumSect,show,fileout)
    if smooth == "Yes":
        CSA_S,WP_S,HR_S,LW_S = MovingAvg(CSA_a,WP_a,HR_a,LW_a,wind)
        PlotGeom(CSA_S,WP_S,HR_S,LW_S,sects,NumSect,show,fileout)


Main(args.folder,args.dataset,args.show,args.fileout,args.smooth,args.wind)