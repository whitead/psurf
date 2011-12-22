#!/usr/bin/python

from psurf import *
import sys


def printHelp():   #system help for the program
    print "getsurfresidue [input pdb file] [output file name] [optional XYZ coordinates for interior]"
    exit(0)

if( len(sys.argv) < 3 ):
    printHelp()

inf = sys.argv[1]
outf = sys.argv[2]

xyzfile = ""
if(len(sys.argv) == 4):
    xyzfile = sys.argv[3]

protein = readProteinSA(inf)

data = []
center = [0,0,0]
atomnum = [0]

#Calculate the center of the molecule
for n in range(0,len(protein)):
    r1 = protein.getResidues()[n]
    if r1.isSurf() == True :

        atomnum[0] += r1.getAtomNum()
        for j in range(0,r1.getAtomNum()):
            a1 = r1.getSingleAtom(j)

            center[0] += a1.getCoord()[0]
            center[1] += a1.getCoord()[1]
            center[2] += a1.getCoord()[2]
                   
center[0] = center[0] / atomnum[0]
center[1] = center[1] / atomnum[0]
center[2] = center[2] / atomnum[0]
print("Center Point Coordinate is %s" %(center))


#Not sure what this part does
for n in range(0,len(protein)):
    for j in range(0,protein.getResidues()[n].getAtomNum()):
        a1 = protein.getResidues()[n].getSingleAtom(j)

xyzOut = False
if(xyzfile != ""):
    xyzOut = open(xyzfile, "w")

with open(outf, 'w') as g:
    g.write("pdb_id, chain, res_index, res_type, location, res_surface_area_ratio\n")
    for n in range(0,len(protein)):
        r1 = protein.getResidues()[n]
        if r1.isSurf() == True :
            for j in range(0,len(protein)):
                r2 = protein.getResidues()[j]
                if r2.getChain() == r1.getChain():
                    if abs(i-j) > 1:
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
                #write XYZ coordinates of this interior residue
                if(xyzOut):
                    for atom in r.getAtoms():
                        g.write("%s\n" % atom.printCoord())
                    
        else:
           g.write("%s, %s, %s, %s, %s, " %(inf[0:4], r1.getChain(), r1.getIndex(), r1.getType(), "Buried"))  
           try:
               g.write("%6.4f\n" % (r1.getSA() / aaSA[r1.getType()]))
           except KeyError:
               pass 

if(xyzOut):
    xyzOut.close()
