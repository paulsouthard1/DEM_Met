import arcpy
from arcpy.sa import Hillshade
import os
import argparse

arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")

def Filename(filein):
	demname = os.path.basename(filein)
	demname = demname[:demname.find(".")]
	dirname = os.path.dirname(filein)
	return demname,dirname

def HS(filein,demname,dirname):
	fileout = dirname +"\\" + demname + "_HS.tif"
	outHS = Hillshade(filein)
	outHS.save(fileout)
	print(fileout)

def Main(filein):
	demname,dirname = Filename(filein)
	HS(filein,demname,dirname)
	print("Hillshade created")


parser = argparse.ArgumentParser(description='Make a Hillshade of DEM usiung ArcGIS')
parser.add_argument('filein', help='Name of input DEM')
args = parser.parse_args()

Main(args.filein)