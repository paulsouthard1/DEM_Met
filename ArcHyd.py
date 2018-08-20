import arcpy as ap
from arcpy import sa
from arcpy import env
import os
import shutil
import argparse

dir_path = os.path.dirname(os.path.realpath(__file__))+"\\"
#env.workspace = dir_path
#env.scratchWorkspace = dir_path
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

parser = argparse.ArgumentParser(description='Run ArcHydro tools DEM to produce flow accumulation and accompanying ASCII files, as well as ASCII files of Anuga Results')
parser.add_argument('demin', help='DEM of Region')
parser.add_argument('momin', help='Anuga Momentum results raster')
parser.add_argument('depin', help='Anuga Depth results raster')
parser.add_argument('stgin', help='Anuga Stage results raster')
args = parser.parse_args()

newpath = dir_path + "ArcFiles\\" 
if not os.path.exists(newpath):
    os.makedirs(newpath)
env.workspace = newpath
env.scratchWorkspace = newpath
#USER SET FILENAME HERE - text before extension
demin = args.demin
filein = os.path.basename(demin)
filein = filein[:filein.find(".")]

#Create ASCII of DEM
asciiout = newpath + filein + "_demasc.txt"
ap.RasterToASCII_conversion(demin,asciiout)

#Name Output ASCII File, then run ArcHydro tools and export ASCII of filled DEM and flow accumulation
asciiout = newpath + filein + "_fdemasc.txt"
outfill = sa.Fill(demin)
ap.RasterToASCII_conversion(outfill,asciiout)
asciiout = newpath + filein + "_faasc.txt"
outfd = sa.FlowDirection(outfill)
outfa = sa.FlowAccumulation(outfd,"","INTEGER")
ap.RasterToASCII_conversion(outfa,asciiout)

#Name Anuga Result ASCII files and export as ASCII
momin=args.momin
filein = os.path.basename(momin)
filein = filein[:filein.find(".")]
asciiout = newpath + filein + "_asc.txt"
ap.RasterToASCII_conversion(momin,asciiout)

depin=args.depin
filein = os.path.basename(depin)
filein = filein[:filein.find(".")]
asciiout = newpath + filein + "_asc.txt"
ap.RasterToASCII_conversion(depin,asciiout)

stgin=args.stgin
filein = os.path.basename(stgin)
filein = filein[:filein.find(".")]
asciiout = newpath + filein + "_asc.txt"
ap.RasterToASCII_conversion(stgin,asciiout)