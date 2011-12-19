#! /usr/bin/python

import psurf, shlex, subprocess

pcaProgName = "/home/whitead/Documents/ProteinSurfaces/psurf/pca.R"

ids = []
charges = []
SAs = []
scounts = []
lengths = []
gaps = []
lambdaMatrix = []

def getLambdas(filename):
    args = shlex.split("%s %s" % (pcaProgName, filename))
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    p.wait()
    output = p.stdout.readline()
    lambdas = [float(x) for x in output.split()]
    return lambdas
    
 
with open("chargeList", 'r') as f:
    for line in f.readlines():
        pdbid = None
        try:
            pdbid = line.split()[0]
            if(len(line.split()) == 3):
                charge = line.split()[1]
                SA = float(line.split()[2])
            else:
                SA = float(line.split()[1])
                charge = "NULL"
        except:
            print "Skipping line \"%s\"" % line[:-1]
            continue
        print "Reading %s" % pdbid
        p = psurf.readProteinSA("%s.surf" % pdbid, "%s.dssp" % pdbid)
        ids.append(pdbid)
        charges.append(charge)
        SAs.append(SA)
        scounts.append(p.getSurfaceCounts())
        lengths.append(len(p))
        gaps.append(p.getGapCount())
        p.writeXYZAtoms("%s.xyz" % pdbid)
        lambdas = getLambdas("%s.xyz" % pdbid)
        lambdaMatrix.append(lambdas)
        print "Writing %s" % pdbid
        p.writeCSVAtoms("%s_3.csv" % pdbid)
        p.writeCSVResidues("%s_2.csv" % pdbid)

with open("chargeList.csv", 'w') as f:
    f.write("pdb_id,charge,surface_area, res_num, gaps, lambda_x, lambda_y, lambda_z")
    for k in psurf.AAs:
        f.write(",%s" % k)
    f.write("\n")
    for (i,c,s,l,g,lm, sf) in zip(ids, charges, SAs, lengths, gaps, lambdaMatrix, scounts):
        f.write("%s,%s,%f,%d,%d" % (i,c,s,l,g))
        for lv in lm:
            f.write(",%g" % lv)
        for k in psurf.AAs:
            f.write(",%d" % sf[k])
        f.write("\n")
            
