import pickle
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mpl
from bokeh.palettes import RdGy as palette
from matplotlib import rcParams
import os
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma']
rcParams['font.size'] = 20
rcParams['font.style'] = 'italic'
rcParams['axes.labelsize'] = 10
rcParams['axes.labelweight'] = 'medium'
rcParams['xtick.labelsize'] = 10
rcParams['ytick.labelsize'] = 10


def ReadPickle(lookuppath):
	#Load Pickle File w/ slope data
	print("Loading Metadata Pickle File")
	sect = pickle.load(open(lookuppath + "MetaSlope ","rb"))
	print("Loading Slope Data Pickle File")
	slope = pickle.load(open(lookuppath + "slope","rb"))
	Elevs = slope["Elevs"]
	Longs = slope["Longs"]
	inds = slope["inds"]
	points = slope["points"]
	return sect, slope, Elevs, Longs, inds, points

def SegmentMatch(sect, Elevs, Longs, inds):
	print("Matching segments to line")
	segments = sect["NumSections"]
	ind = np.zeros((segments,3))
	for i in np.arange(segments):
		lookup = "seg_{}".format(i+1)
		ind[i,0] = inds[lookup][0]
		ind[i,1] = Longs[inds[lookup]]
		ind[i,2] = Elevs[inds[lookup]]
	return ind, segments

def MakePlt(fileout,SP,show,segments,Longs,Elevs,Surfs,ind):
	#Make Slope plot with squares at cross-section locations, optionally plot modeled water surface also
	print("Making Plot")
	fig = plt.subplots(nrows = 1, ncols = 1,figsize = (15,10))
	ax1 = plt.subplot(111)
	colorlist = palette[10]
	listnum = 0
	for i in np.arange(0,int(segments),int(show)):
	    if listnum == 10:
	        listnum = 0
	    else:
	        listnum = listnum
	    col = colorlist[listnum]
	    ax1.plot(ind[i,1], ind[i,2], marker = 's',color=col, mec = 'k', markersize = 10)
	    listnum = listnum + 1
	ax1.plot(Longs,Elevs,color = '#893a20',linewidth = 2, alpha = 0.9)
	if SP == True:
	    ax2 = ax1.twinx()
	    ax2.plot(Longs,Surfs,color = '#206f89',linewidth = 2, alpha = 0.9)
	    ax1.set_ylabel("Elevation (m)", size = 20, color = '#893a20')
	    ax2.set_ylabel("Modeled River Stage (m)", size = 20, color = '#893a20')
	else:
	    ax1.set_ylabel("Elevation (m)", size = 20, color = 'k')
	ax1.set_xlim(-20,np.max(Longs)+20)
	ax1.set_xlabel("Longitudinal Distance (m)", size = 20, color = 'k')
	ax1.set_ylabel("Elevation (m)", size = 20, color = 'k')
	plt.title("Longitudinal Profile",size = 25)
	plt.savefig(fileout,dpi=500,bbox_inches = 'tight')
	plt.close()

#ArgParse stuff
parser = argparse.ArgumentParser(description='Plot longitudinal profile with optional water surface, including cross-sections')
parser.add_argument('folder', help='Folder for analysis')
parser.add_argument('fileout', help='Output plot filename')
parser.add_argument('SP', help='Are we plotting a surface profile? "Yes" or "No".')
parser.add_argument('show', help='Show every _th cross section point')
args = parser.parse_args()



def Main(folder,fileout,SP,show):
	dir_path = os.path.dirname(os.path.realpath(__file__)) + "\\"
	arcpath = dir_path + "ArcFiles\\"
	crosspath = dir_path + "CrossFiles\\"
	if not os.path.exists(crosspath):
		os.makedirs(crosspath)
	#Define directory for results storage
	lookuppath = crosspath + folder + "\\"
	if not os.path.exists(lookuppath):
		os.makedirs(lookuppath)
	Surfs = 0
	sect, slope, Elevs, Longs, inds, points = ReadPickle(lookuppath)
	ind, segments = SegmentMatch(sect, Elevs, Longs, inds)
	MakePlt(fileout,SP,show,segments,Longs,Elevs,Surfs,ind)

Main(args.folder,args.fileout,args.SP,args.show)
