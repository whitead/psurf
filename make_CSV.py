#! /usr/bin/python

import psurf

ids = []
charges = []
SAs = []
scounts = []
lengths = []
gaps = []

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
        p = psurf.readProteinSA("%s.surf" % pdbid)
        ids.append(pdbid)
        charges.append(charge)
        SAs.append(SA)
        scounts.append(p.getSurfaceCounts())
        lengths.append(len(p))
        gaps.append(p.getGapCount())
        print "Writing %s" % pdbid
        p.writeCSVAtoms("%s_3.csv" % pdbid)
        p.writeCSVResidues("%s_2.csv" % pdbid)

with open("chargeList.csv", 'w') as f:
    f.write("pdb_id,charge,surface_area, res_num, gaps")
    for k in psurf.AAs:
        f.write(",%s" % k)
    f.write("\n")
    for (i,c,s,l,g,sf) in zip(ids, charges, SAs, lengths, gaps, scounts):
        f.write("%s,%s,%f,%d,%d" % (i,c,s,l,g))
        for k in psurf.AAs:
            f.write(",%d" % sf[k])
        f.write("\n")
            
