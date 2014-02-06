#!/usr/bin/python

from os.path import isfile
from os import listdir
from os import sys
from optparse import OptionParser # Command line parsing
usage = "usage: %prog summary files"
version = "%prog."
parser = OptionParser(usage=usage,version=version)
parser.add_option("-x", "--xsec", action="store", dest="xsec", type="float", default=0.001, help="set xsec scale [default 0.001]")

(options, args) = parser.parse_args()

in_file = args[0] #"/local/cms/user/pastika/heavyNuAnalysis_2012/skims/WW_TuneZ2star_8TeV_pythia6_tauola_START53_V7A-v1/res/"


logFiles = []
for logFile in listdir(in_file):
	fullPath = in_file + '/' + logFile
	if logFile.endswith("-summary.log") and isfile(fullPath):
		logFiles.append(fullPath)

of = open(in_file + "/limit_summary.txt", "w")

output = []
for logFile in logFiles:
	nums = []
	openFile = open(logFile)
	for line in openFile:
		sLine = line.split("=")
		nums.append(float(sLine[1]))
	if len(nums) > 7:
		#of.write("%i %i %e %e %e %e %e %e\n"%(nums[0], nums[1], nums[2]/100, nums[3]/100, nums[7]/100, nums[5]/100, nums[4]/100, nums[6]/100))
		output.append([nums[0], nums[1], nums[2]*options.xsec, nums[3]*options.xsec, nums[7]*options.xsec, nums[5]*options.xsec, nums[4]*options.xsec, nums[6]*options.xsec])

output.sort()
for i in output:
	of.write("%i %i %e %e %e %e %e %e\n"%(i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7]))
