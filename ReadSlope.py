import pickle


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
	ind = np.zeros((segments,3))
	segments = sect["NumSections"]
	for i in np.arange(segments):
    	lookup = "seg_{}".format(i+1)
    	ind[i,0] = inds[lookup][0]
    	ind[i,1] = Longs[inds[lookup]]
    	ind[i,2] = Elevs[inds[lookup]]
    return ind


def Main(folder):
	dir_path = os.path.dirname(os.path.realpath(__file__)) + "\\"
	arcpath = dir_path + "ArcFiles\\"
	crosspath = dir_path + "CrossFiles\\"
	if not os.path.exists(crosspath):
		os.makedirs(crosspath)
	#Define directory for results storage
	lookuppath = crosspath + folder + "\\"
	if not os.path.exists(lookuppath):
		os.makedirs(lookuppath)
	sect, slope, Elevs, Longs, inds, points = ReadPickle(lookuppath)
