#!/usr/bin/python2.7
import sys


def printHelp():
    print "analyzeSA [file] [SA column number] [output file]"
    exit(0)

if( len(sys.argv) != 4 ):
    printHelp()

inf = sys.argv[1]
cnumber = int(sys.argv[2])
outf = sys.argv[3]

#Load file and sort into list of residues and SAs
masses = {'C':12.011, 'N':14.007, 'H':1.0079, 'S':32.065, 'O':15.999}
mx = 0
my = 0
mz = 0
vx = 0
vy = 0
vz = 0
residues = []
lresID = 0
lres = ""
sa = 0
netSA = 0
heavyAtomN = 0
masstot = 0
with open(inf, 'r') as f:
    for l in f.readlines():
        if(len(l) < cnumber):
            continue
        resID = l[22:27]
        res = l[17:20]
        if(resID != lresID):
            if(lresID != 0):
                residues.append((lres,  heavyAtomN, mx / masstot, my / masstot, mz / masstot, sa, vx, vy, vz))
            netSA += sa
            sa = 0
            heavyAtomN = 0
            lres = res
            lresID = resID
            masstot = 0
            mx = 0
            my = 0
            mz = 0
            vx = 0
            vy = 0
            vz = 0
        sa += float(l[cnumber:-2]) #skip new line
        heavyAtomN += 1
        #calculate center of mass
        if( not masses.has_key(l[13]) ):
            mass = 12
        else:
            mass = masses[l[13]]
        masstot += mass
        mx += float(l[30:37]) * mass
        my += float(l[38:45]) * mass
        mz += float(l[46:54]) * mass
        if(l[13:15] == "CB"):
            vx -= float(l[30:37])
            vy -= float(l[38:45])
            vz -= float(l[46:54])
        elif(l[13:15] == "CA" and l[17:20] != "GLY"):
            vx += float(l[30:37])
            vy += float(l[38:45])
            vz += float(l[46:54])
            

residues.append((lres,  heavyAtomN, mx / masstot, my / masstot, mz / masstot, sa, vx, vy, vz))
netSA += sa
print netSA
with open(outf, 'w') as o:
    o.write("res num x y z sa vx vy vz\n")
    for p in residues:
        o.write("%s %d %g %g %g %g %g %g %g\n" % (p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8]))
