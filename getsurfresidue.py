#!/usr/bin/python

from psurf import *
import sys


def printHelp():   #system help for the program
    print "getsurfresidue [input pdb file] [output file name]"
    exit(0)

if( len(sys.argv) != 3 ):
    printHelp()

inf = sys.argv[1]
outf = sys.argv[2]

f = readProteinSA(inf)

data = []
center = [0,0,0]
atomnum = [0]

#with open(outf, 'w') as g:
#    with open(outf2, 'w') as k:
#        g.write("Chain Rnum Res SA\n")
#        k.write("Chain Rnum Atom x y z\n")

#Calculate the center of the molecule
for n in range(0,len(f)):
    r1 = f.getResidues(n)
    if r1.isSurf() == True :
               #g.write("%s %s %s %s\n" %(r1.getChain(), r1.getIndex(), r1.getType(), r1.getSA()))
        atomnum[0] += r1.getAtomNum()
        for j in range(0,r1.getAtomNum()):
            a1 = r1.getSingleAtom(j)
                   #k.write("%s %s %s %s %s %s\n" %(r1.getChain(), r1.getIndex(), a1.getType(), a1.getCoord()[0], a1.getCoord()[1], a1.getCoord()[2]))
            center[0] += a1.getCoord()[0]
            center[1] += a1.getCoord()[1]
            center[2] += a1.getCoord()[2]
                   
center[0] = center[0] / atomnum[0]
center[1] = center[1] / atomnum[0]
center[2] = center[2] / atomnum[0]
print("Center Point Coordinate is %s" %(center))

#with open(outf3, 'w') as d:
#    d.write("Chain Rnum Atom dist\n")

#Not sure what this part does
for n in range(0,len(f)):
    for j in range(0,f.getResidues(n).getAtomNum()):
        a1 = f.getResidues(n).getSingleAtom(j)
            #d.write("%s %s %s %s\n" %(f.getResidues(n).getChain(), f.getResidues(n).getIndex(), a1.getType(), a1.centerdistance(center)))

with open(outf, 'w') as g:
    g.write("PDBID, Chain, Rnum, Res, Location, SAfrac\n")
    for n in range(0,len(f)):
        r1 = f.getResidues(n)
        if r1.isSurf() == True :
            for j in range(0,len(f)):
                r2 = f.getResidues(j)
                if r2.getChain() == r1.getChain():
                    if j > n+1 or j < n-1:
                        d = r1.ResIntersect(r2, center, 1)
                        control = []
                        if d == "Outside" :
                            control = [r1.getChain(), r1.getIndex(), r1.getType()]
                            g.write("%s, %s, %s, %s, %s, " %(inf[0:4], r1.getChain(), r1.getIndex(), r1.getType(),"Out"))
                            try:
                                g.write("%6.4f\n" % (r1.getSA() / aaSA[r1.getType()]))
                            except KeyError:
                                pass
                            break
            r = [r1.getChain(), r1.getIndex(), r1.getType()]
            if r != control :
                g.write("%s, %s, %s, %s, %s, " %(inf[0:4], r1.getChain(), r1.getIndex(), r1.getType(), "In"))  
                try:
                    g.write("%6.4f\n" % (r1.getSA() / aaSA[r1.getType()]))
                except KeyError:
                    pass
        else:
           g.write("%s, %s, %s, %s, %s, " %(inf[0:4], r1.getChain(), r1.getIndex(), r1.getType(), "Buried"))  
           try:
               g.write("%6.4f\n" % (r1.getSA() / aaSA[r1.getType()]))
           except KeyError:
               pass 
