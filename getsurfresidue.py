#!/usr/bin/python

import psurf
import sys


def printHelp():   #system help for the program
    print "getsurfresidue [input pdb file] [output residue list file] [output atom list file] [output distance file]"
    exit(0)

if( len(sys.argv) != 5 ):
    printHelp()

inf = sys.argv[1]
outf = sys.argv[2]
outf2 = sys.argv[3]
outf3 = sys.argv[4]

f = psurf.readProteinSA(inf)

data = []
center = [0,0,0]
atomnum = [0]

with open(outf, 'w') as g:
    with open(outf2, 'w') as k:
        g.write("Chain Rnum Res SA\n")
        k.write("Chain Rnum Atom x y z\n")
        for n in range(0,len(f)):
            r1 = f.getResidue(n)
            if r1.isSurf() == True :
               g.write("%s %s %s %s\n" %(r1.getChain(), r1.getresNum(), r1.getType(), r1.getSA()))
               atomnum[0] += r1.getAtomNum()
               for j in range(0,r1.getAtomNum()):
                   a1 = r1.getAtom(j)
                   k.write("%s %s %s %s %s %s\n" %(r1.getChain(), r1.getresNum(), a1.getType(), a1.getCoord()[0], a1.getCoord()[1], a1.getCoord()[2]))
                   center[0] += a1.getCoord()[0]
                   center[1] += a1.getCoord()[1]
                   center[2] += a1.getCoord()[2]
                   

center[0] = center[0] / atomnum[0]
center[1] = center[1] / atomnum[0]
center[2] = center[2] / atomnum[0]
print("Center Point Coordinate is %s" %(center))

with open(outf3, 'w') as d:
    d.write("Chain Rnum Atom dist\n")
    for n in range(0,len(f)):
    	for j in range(0,f.getResidue(n).getAtomNum()):
            a1 = f.getResidue(n).getAtom(j)
            d.write("%s %s %s %s\n" %(f.getResidue(n).getChain(), f.getResidue(n).getresNum(), a1.getType(), a1.centerdistance(center)))

with open('OutSurRes.txt', 'w') as g:
    with open('InSurRes.txt', 'w') as k:
        g.write("Chain Rnum Res\n")
        k.write("Chain Rnum Res\n")
        for n in range(0,len(f)):
            r1 = f.getResidue(n)
            if r1.isSurf() == True :
                 for j in range(0,len(f)):
                      r2 = f.getResidue(j)
                      if r2.getChain() == r1.getChain():
                          if j > n+1 or j < n-1:
                              d = r1.ResIntersect(r2, center, 1)
                              #if d != None:
                              #   g.write("%s\n" %(d))
                              control = []
                              if d == "Outside" :
                                 control = [r1.getChain(), r1.getresNum(), r1.getType()]
                                 g.write("%s %s %s\n" %(r1.getChain(), r1.getresNum(), r1.getType()))
                                 break
                 r = [r1.getChain(), r1.getresNum(), r1.getType()]
                 if r != control :
                    k.write("%s %s %s\n" %(r1.getChain(), r1.getresNum(), r1.getType()))         
