#! /usr/bin/python

from psurf import *
import sys, math


def readChains(chainFile):
    chains = []
    with open(chainFile, 'r') as f:
        for line in f.readlines():
            chains.append(line.split())
    return chains



def printHelp():
    print "getInteriorResidues.py [surf file] [chain file] [surface cutoff] [occlusion margin] [center axis unit vector] [minMax] [rmax]"
    print "The chain file should have a list of chains which are identical on each line. Only the first chain of each line will be calculated and the results will be replicated"

def testIntersectAlgorithm(origin, points, centralAxis):
    for p in points:
        #Get the point on the central axis that is closest to the origin atom
        #cv = s * axis, s = P_o dot axis. Assuming axis is a unit vector
        tempDP = sum([x * y for x,y in zip(origin, centralAxis)])
        centerVector = [x * tempDP for x in centralAxis]
        
        avector = [x - y for x,y in zip(p, origin)]
        bvector = [x - y for x,y in zip(centerVector, origin)]
        
            #Check to see if we need to consider occlusion (because the other atom is too far away)
        if(sum([x**2 for x in avector]) > sum([x**2 for x in bvector])):
            print "Bad distance"
            continue
        
            #|A|^2 - (A * B)^2 / |B|^2
        
            #(A * B)
        tempDP = sum([x * y for x,y in zip(avector, bvector)])
        
            #Check if other atom is "behind" the first atom. As in the angle is obtuse
        if(tempDP < 0):
            print "Bad angle"
            continue
        
        radiusSqr = sum([x**2 for x in avector]) - tempDP**2 / sum([x**2 for x in bvector])
        
        print(radiusSqr)

def getInterior(chain, protein, centralAxis, minMax, margin, rmax, heavyMargin = 3):
    residueIndices = []
    count = 0
    rlength = len(protein.getResidues(chain))
    for rorigin in protein.getResidues(chain):
        if(rorigin.isSurf()):
            for atom in rorigin.getAtoms():
                if(atom.isHeavy()):
                    atom.occluded = False
            for rother in protein.getResidues():
                if(rother != rorigin and rother.isSurf()):
                    centerIntersect(rorigin, rother, centralAxis, minMax, margin, rmax)

            clearCount = 0
            heavyCount = 0
            for atom in rorigin.getAtoms():
                if(atom.isHeavy()):
                    heavyCount += 1
                    if(not atom.occluded):
                        clearCount += 1
            if(clearCount >= heavyMargin):#float(clearCount) / heavyCount >= heavyMargin):
                residueIndices.append(rorigin.getIndex())

        count += 1
    return residueIndices
                    
            
def centerIntersect(rorigin, rother, centralAxis, minMax, margin, rmax):
    for a1 in rorigin.getAtoms():
        if(not a1.isHeavy()):
            continue
        for a2 in rother.getAtoms():
            if(not a2.isHeavy()):
                continue

            #Get the point on the central axis that is closest to the origin atom
            #cv = s * axis, s = P_o dot axis. Assuming axis is a unit vector
            tempDP = sum([x * y for x,y in zip(a1.getCoord(), centralAxis)])
            if(tempDP < minMax[0] or tempDP > minMax[1]):
                a1.occluded = True
                continue

            centerVector = [x * tempDP for x in centralAxis]

            avector = [x - y for x,y in zip(a2.getCoord(), a1.getCoord())]
            bvector = [x - y for x,y in zip(centerVector, a1.getCoord())]

            #Check if the atom is too far away to be considered in the interior
            if(sum([x**2 for x in bvector]) > rmax*rmax):
                a1.occluded = True
                continue
            
            #Check to see if we need to consider occlusion (because the other atom is too far away)
            if(sum([x**2 for x in avector]) > sum([x**2 for x in bvector])):
                continue

            #|A|^2 - (A * B)^2 / |B|^2
            
            #(A * B)
            tempDP = sum([x * y for x,y in zip(avector, bvector)])
            
            #Check if other atom is "behind" the first atom. As in the angle is obtuse
            if(tempDP < 0):
                continue
            
            radiusSqr = sum([x**2 for x in avector]) - tempDP**2 / sum([x**2 for x in bvector])

            if(radiusSqr < margin**2):
                a1.occluded = True


def main():
    
    if(len(sys.argv) != 11):
        printHelp()
        exit()
    
    inputFile = sys.argv[1]
    chainFile = sys.argv[2]
    surfCutoff = float(sys.argv[3])
    margin = float(sys.argv[4])
    centralAxis = [float(x) for x in sys.argv[5:8]]
    minMax = [float(x) for x in sys.argv[8:10]]
    rmax = float(sys.argv[10])
 
    chains = readChains(chainFile)

    p = readProteinSA(inputFile)

    #center the protein
    p.centerProtein()
    
    setSurfCutoff(surfCutoff)

    for r in p.getResidues():
        for a in r.getAtoms():
            if(r.isSurf()):
                a.setBeta(1.0)
            else:
                a.setBeta(0.0)

    for chainList in chains:
        residueIndices = getInterior(chainList[0], p, centralAxis, minMax, margin, rmax)
        for c in chainList:
            print c
            for ri in residueIndices:
                print ri
                res = p.getResidue(ri, c)
                if(res):
                    for a in res.getAtoms():
                        a.setBeta(-1.0)
                else:
                    print "Could not find residue index %d in chain %s" % (ri, c)
    p.writePDB("%s_beta.pdb" % p.getName())

    #write an XYZ file containing the interior CA's
    p.writeXYZAtoms("%s_beta.xyz" % p.getName(), beta=-1)

main()

#points = [[2,3,0], [3,2,0], [1.95,1.75, 0], [-1.95,1.75, 0], [-2,1, 0], [-1.95,-1.75, 0], [2.1, 1, 0]]
#origin = [2,1,0]
#ca = [0,1,0]
#testIntersectAlgorithm(origin, points, ca)
