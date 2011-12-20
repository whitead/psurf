#!/usr/bin/python
import sys
from psurf import *

if(len(sys.argv) != 6):
    print "Usage: [countSurfacecontacts.py] [surf_id list] [output] [cutoff] [ionic only] [include backbone]"
    exit()

setSurfCutoff(float(sys.argv[3]))
ionicOnly = False
if(sys.argv[4] in ['T', 't', 'True', 'true']):
    ionicOnly = True

backbone = False
if(sys.argv[5] in ['T', 't', 'True', 'true']):
    backbone = True

if(ionicOnly):
    print "Disabling non-ionic radii"
    OPLSSigmas['C'] = 0.0
    OPLSSigmas['S'] = 0.0

listFile = sys.argv[1]
pdbIDs = []

with open(listFile, 'r') as f:
    for line in f.readlines():
        if(len(line) > 1):
            if(line[-1] == "\n"):
                pdbIDs.append(line[:-1])
            else:
                pdbIDs.append(line)

counts = [[0 for x in AAs] for x in AAs]
counts.append([0 for x in AAs]) #water
counts.append([0 for x in AAs]) #hydrated
counts.append([0 for x in AAs]) #non-contact counts
counts.append([0 for x in AAs]) # total count    

def clearCounts():
    for i in range(len(counts)):
        for j in range(len(counts[i])):
            counts[i][j] = 0
        

def incrementCount(t1, t2):
    try:
        if(t1 == "WATER"):
            counts[len(counts) - 4][AAs.index(t2)] += 1
        elif(t1 == "HYDRATED"):
            counts[len(counts) - 3][AAs.index(t2)] += 1
        elif(t1 == "FREE"):
            counts[len(counts) - 2][AAs.index(t2)] += 1
        elif(t1 == "TOTAL"):
            counts[len(counts) - 1][AAs.index(t2)] += 1
        else:
            counts[AAs.index(t1)][AAs.index(t2)] += 1
    except ValueError:
        print "Unrecognizable type pair (%s, %s)" % (t1,t2)


output = []

output.append("%-5s, %-8s, " % ("pdb_id", "res_type"))
for i in range(len(counts[0]) - 1):
    output.append("%6s, " % AAs[i])
output.append("%6s\n" % AAs[-1])

with open(sys.argv[2], 'w') as f:
    f.write("".join(output))

for i in pdbIDs:
    output = []
    clearCounts()
    p = readProteinSA(i)
    print i,
    for r1 in range(len(p)):	
        res1 = p.residues[r1]
        if(res1.isSurf() or res1.isHydrated()):
            contact = False
            hydrationNumber = res1.hydrationNumber(backbone)
            if(hydrationNumber > 0):
                incrementCount("HYDRATED", res1.getType())
                for h in range(hydrationNumber):
                    incrementCount("WATER", res1.getType())
            for r2 in range(len(p)):
                res2 = p.residues[r2]
                if(not res1.isNeighbor(res2) and res1.inContact(res2, backbone)):
                    incrementCount(res1.getType(), res2.getType())
                    contact = True
            if(not contact):
                incrementCount("FREE", res1.getType())
            incrementCount("TOTAL", res1.getType())
        print "\r%d/%d" % (r1+1, len(p)),
    print ""

    for j in range(len(counts)):
        output.append("%-5s, " % p.getName())
        if(j < len(AAs)):
            output.append("%-8s, " % AAs[j])
        elif(j == len(AAs)):
            output.append("%-8s, " % "WATER")
        elif(j == len(AAs) + 1):
            output.append("%-8s, " % "HYDRATED")
        elif(j == len(AAs) + 2):
            output.append("%-8s, " % "FREE")
        else:
            output.append("%-8s, " % "TOTAL")
        for k in range(len(counts[j])- 1):
            output.append("%6g, " % counts[j][k])
        output.append("%6g\n" % counts[j][-1])
    with open(sys.argv[2], 'a') as f:
        f.write("".join(output))
