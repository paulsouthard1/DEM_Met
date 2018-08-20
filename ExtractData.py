# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 16:45:49 2018

@author: ps29626
"""

import os
import arcpy
import numpy as np
import pickle
import argparse

arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")
arcpy.CheckOutExtension("3d")

#ArgParse stuff
parser = argparse.ArgumentParser(description='Extract lines and swaths of data along a reach and make shapefiles for visualization')
parser.add_argument('folder', help='Folder to store resulting files in A.K.A. name describing analysis')
parser.add_argument('demin', help='DEM of Region')
parser.add_argument('dataset', help='dataset to extract swaths from')
parser.add_argument('clip', help='Anuga result used to clip lines')
parser.add_argument('lines', help='Do you need to make line shapefiles?  "Yes" or "No".')
#parser.add_argument('swaths', type=bool, help='Do you need to make swath shapefiles?')
args = parser.parse_args()

#Define function to read MetaData
def ReadData(lookuppath):
    #Read in pickle file with number of sections
    print("Loading Metadata Pickle File")
    sect = pickle.load(open(lookuppath + "sectnum","rb"))
    #Assign number of sections to variable
    NumSect = sect["NumSections"]
    return NumSect

#Define function to produce shapefile of flow extent
def RastDomain(storepath,clip,clipname):
    print("Creating Raster Domain")
    clipshape = storepath + clipname + ".shp"
    arcpy.RasterDomain_3d(clip,clipshape,"POLYGON")
    print("Raster Domain created")
    return clipshape

#Define function to Load Cross Sections as shapefiles
def LoadSect(lookuppath,storepath,arcpath,demname,NumSect):
    #Iterate through number of segments generated
    startx_field = "RBX"
    starty_field = "RBY"
    endx_field = "LBX"
    endy_field = "LBY"
    prj = arcpath + demname + "_fdemasc.prj"
    for j in np.arange(2,NumSect):
        #Create key to find text file
        lookup = "CS2_{}".format(j) + ".txt"
        print("Creating Line from "+ lookup)
        #Assign parameters for ArcTool
        in_table = lookuppath + "CS2_{}".format(j)+".txt"
        out_featureclass = storepath+"CS2_{}".format(j)+".shp"
        arcpy.XYToLine_management(in_table, out_featureclass, startx_field, starty_field, endx_field, endy_field, "GEODESIC", "", prj)
        
#Clip lines to modeled flow width
def cliplines(storepath,clipshape,NumSect):
    for j in np.arange(2,NumSect):
        #Create key to find text file
        lookup = "CS2_{}".format(j) + ".shp"
        inshape = storepath + lookup
        outshape = storepath + "cl_"+ lookup
        print("Clipping Line "+ lookup)
        arcpy.Clip_analysis(inshape,clipshape,outshape)

#Pull data from given dataset along cross-section lines
def LinePull(storepath,lookuppath,dataset,dataname,NumSect):
    for j in np.arange(2,NumSect):
        #Create key to find text file
        lookup = "CS2_{}".format(j) + ".shp"
        print("Pulling data from " + lookup)
        in_line_features = storepath + "cl_"+ lookup
        out_table = lookuppath+dataname+"_"+str(j)+".txt"
        arcpy.StackProfile_3d(in_line_features,dataset,out_table)

#Pull data from given dataset for given swath
def SwathPull(lookuppath,storepath,dataset,dataname,NumSect):
    #Read in pickle file with points
    print("Loading points Pickle File")
    points = pickle.load(open(lookuppath + "points","rb"))
    print("Successfully loaded points")
    for j in np.arange(2,NumSect):
        key2 = "seg_{}".format(j)
        key1 = "seg_{}".format(j-1)
        print("Extracting swath between "+key1+" and "+key2)
        line1 = points[key1]
        line2 = points[key2]
        maxx = np.max([line1[0],line2[0]])
        minx = np.min([line1[0],line2[0]])
        maxy = np.max([line1[1],line2[1]])
        miny = np.min([line1[2],line2[2]])
        if maxx == minx or maxy == miny:
            print("Invalid swath "+key1)
        else:
            rectangle = str(minx)+" "+str(miny)+" "+str(maxx)+" "+str(maxy)
            outraster = storepath+dataname+"_"+str(j)+".tif"
            arcpy.Clip_management(dataset,rectangle,outraster,"#","-9999","NONE","MAINTAIN_EXTENT")
            out_ascii_file = lookuppath+dataname+"_"+key1+".txt"
            arcpy.RasterToASCII_conversion(outraster,out_ascii_file)
            
def Main(folder,demin,dataset,clip,lines):
    dir_path = os.path.dirname(os.path.realpath(__file__))+"\\"
    arcpath = dir_path + "ArcFiles\\"
    crosspath = dir_path + "CrossFiles\\"
    if not os.path.exists(crosspath):
        os.makedirs(crosspath)
    storepath = arcpath + folder + "\\"
    if not os.path.exists(storepath):
        os.makedirs(storepath)
    lookuppath = crosspath + folder + "\\"
    demname = os.path.basename(demin)
    demname = demname[:demname.find(".")]
    NumSect = ReadData(lookuppath)
    dataname = os.path.basename(dataset)
    dataname = dataname[:dataname.find(".")]
    clipname = os.path.basename(clip)
    clipname = clipname[:clipname.find(".")]
    clipshape = RastDomain(storepath,clip,clipname)
    if lines == "Yes":
        print("lines = " + str(lines))
        print(NumSect)
        LoadSect (lookuppath,storepath,arcpath,demname,NumSect)
        cliplines(storepath,clipshape,NumSect)
    else:
        print("Did not make new cross-section lines")
    LinePull(storepath,lookuppath,dataset,dataname,NumSect)

#    if swaths = True:
#        
#        SwathPull(lookuppath,storepath,dataset,dataname,NumSect)
    
parser = argparse.ArgumentParser(description='Extract lines and swaths of data along a reach and make shapefiles for visualization')
parser.add_argument('folder', help='Folder to store resulting files in A.K.A. name describing analysis')
parser.add_argument('demin', help='DEM of Region')
parser.add_argument('dataset', help='dataset to extract swaths from')
parser.add_argument('clip', help='Anuga result used to clip lines')
parser.add_argument('lines', help='Do you need to make line shapefiles?  "Yes" or "No".')
#parser.add_argument('swaths', type=bool, help='Do you need to make swath shapefiles?')
args = parser.parse_args()

Main(args.folder,args.demin,args.dataset,args.clip,args.lines)