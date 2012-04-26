#!/usr/bin/python

from beadlib import *
import sys


def printHelp():
	print "This program is designed to optimized sphere radius based on a given molecule and given molecule center, Proposed Usage is optimizeRadius.py [molecule file in .xyz format],[Atom Index of sphere 1],[Atom Index of sphere 2]...."


if(len(sys.argv) <= 2):
    printHelp()
    exit()

molefile = sys.argv[1]
sphereCenter = range(len(sys.argv)-2)
for x in range(2,len(sys.argv)):
	sphereCenter[x-2] = int(sys.argv[x])


mol=readMolecule(molefile)
radius = optRadius(mol,sphereCenter)
print(radius)

