import math, re

#Created by Andrew White, 2011

#0 atom number
#1 atom type
#2 residue type
#3 chain
#4 residue number
#5 x-cord
#6 y-cord
#7 z-cord
#8 occupancy
#9 beta factor
                     #ATOM    atomnum     atom type             residue type     chain           residue number  x                  y                    z          occupancy   beta factor
regexp = re.compile("ATOM[\s]*([\d]*)[\s]*([\w\d]{1,4})[\s]{0,3}([\w]{1,3})[\s]*([\w]*)[\s]{0,5}([\d]{1,5})[\s]*([-\d\.]{1,8})[\s]*([-\d\.]{1,8})[\s]*([-\d\.]{1,8})[\s]*([\d\.]{1,6})[\s]*([-\d\.]{1,6})[\s]*(.*)")

#added SA column
regexpSA = re.compile("ATOM[\s]*([\d]*)[\s]*([\w\d]{1,4})[\s]{0,3}([\w]{1,3})[\s]*([\w]*)[\s]{0,5}([\d]{1,5})[\s]*(-{0,1}[\d\.]{3,7})[\s]*(-{0,1}[\d\.]{3,7})[\s]*(-{0,1}[\d\.]{3,7})[\s]*(-{0,1}[\d\.]{1,5})[\s]*(-{0,1}[\d\.]{1,5})[\s]*([-\d\.]*)")

regexpWater = re.compile(
"HETATM([\d\s]{5})  O  .HOH [\w][\d\s]{4}.\s{3}([\d\.\s]{8})([-\d\.\s]{8})([-\d\.\s]{8}).*")

pdbFormat = "%7s%5d%4s %3s%1s%4d %8.3f%8.3f%8.3f%6.2f%6.2f    "

#Sequence Reader
regexpSeq = re.compile("SEQRES[\s]*[\d]*[\s]*[\w]*[\s]*[\d]*[\s]*(.*)")

conversion = {'ALA':'A', 'ARG':'R', 'ASP':'D', 'ASN':'N', 'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}

MW = {'A':'89.09', 'R':'174.20', 'D':'133.10', 'N':'132.12', 'C':'121.16', 'Q':'146.14', 'E':'147.13', 'G':'75.07', 'H':'155.15', 'I':'131.17', 'L':'131.17', 'K':'146.19', 'M':'149.21', 'F':'165.19', 'P':'115.13', 'S':'105.09', 'T':'119.12', 'W':'204.23', 'Y':'181.19', 'V':'117.15'}

chis = {'GLU':[0,1], 'LYS':[0,1,2,3], 'SER':[], 'GLY':[], 'PRO':[]}
chiTypes = [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD'], ['CB', 'CG', 'CD' , 'CE'], ['CB', 'CG', 'CD' , 'NZ']]
OPLSSigmas = {'N':3.25, 'O':2.96, 'C':3.75, 'S':3.55, 'H':0.00, 'HOH':3.1507}
backbone = ['N', 'O', 'CA', 'C']

aaSA = {'ALA':113, 'ARG':241, 'ASP':151, 'ASN':158, 'CYS':140, 'GLN':189, 'GLU':183, 'GLY':85, 'HIS':194, 'ILE':182, 'LEU':180, 'LYS':211, 'MET':204, 'PHE':218, 'PRO':143, 'SER':122, 'THR':146, 'TRP':259, 'TYR':229, 'VAL':160}
surfCutoff = 0.4

class Atom: 
   def __init__(self):
      self.coord = (0,0,0)
      self.type = "H"
      self.index = 0
      self.beta = 0
      self.sa = 0
      self.occ = 1
      self.sig = OPLSSigmas["H"]
      self.backbone = False

   def fromData(self, coordinates, atomType, atomIndex, occupancy, beta):
      self.coord = coordinates
      self.type = atomType
      if(self.type in backbone):
         self.backbone = True
      self.index = atomIndex
      self.ba = beta
      self.occ = occupancy
      try:
         if(atomType == "HOH"):
            self.sig = OPLSSigmas[atomType]
         else:
            self.sig = OPLSSigmas[atomType[0]]
      except KeyError:
         self.sig = 0
         print "Warning: could not find atom type for %s" % atomType

   def __str__(self):
      return("Type: \"%s\"; Index: %03d; Coordinates: (%05.3f, %05.3f, %05.3f)" % (self.type, self.index, self.coord[0], self.coord[1], self.coord[2]))
   
   def setSA(self, sa):
      self.sa = sa

   def getSA(self):
      return self.sa

   def distance(self, otherAtom):
      return math.sqrt(self.distanceSqr(otherAtom))

   def distanceSqr(self, otherAtom):
      dist = 0
      for (x1,x2) in zip(self.coord, otherAtom.coord):
         dist += (x1 - x2)**2
      return dist

   def centerDistance(self, center):
      dist = 0
      for (x1,x2) in zip(self.coord, center):
         dist += (x1 - x2)**2
      return math.sqrt(dist)

   def inContact(self, otherAtom):
      if(self.sig + otherAtom.sig > 0.00001):
         if(self.distanceSqr(otherAtom) < self.sig * otherAtom.sig * 2.**(1./6)):
            return True
      return False

   def getBeta(self):
      return self.beta

   def getOccupancy(self):
      return self.occ

   def getType(self):
      return self.type

   def getIndex(self):
      return self.index

   def getCoord(self):
      return self.coord

   def centerdistance(self, center):
      dist = 0
      for (x1,x2) in zip(self.coord, (center[0],center[1],center[2])):
         dist += (x1 - x2)**2
      return math.sqrt(dist)

class Residue:
   def __init__(self):
       self.atoms = []
       self.type = "X"
       self.parent = None
       self.hydrated = None

   def fromData(self, atoms, rType, chain, resNum):
       self.type = rType
       self.atoms = atoms
       self.chain = chain
       self.index = resNum

   def __len__(self):
      return len(self.atoms)

   def setParent(self, p):
      self.parent = p
          
   def inContact(self, other):

      if(self != other):
         if(self.getCA() != None and other.getCA() != None):
            if(self.getCA().distanceSqr(other.getCA()) < 10 ** 2):
               for a in self.atoms:
                  for b in other.atoms:
                     if(a.getType()[0] != "H" and b.getType() != "H"):
                        if(not a.backbone and not b.backbone):
                           if(a.inContact(b)):
                              return True
      return False

   def inWaterContact(self, waterAtom):
      if(self.getCA() != None and  self.getCA().distanceSqr(waterAtom) < 10 ** 2):
         for a in self.atoms:            
            if(a.getType()[0] != "H" and (not a.backbone)):
               if(a.inContact(waterAtom)):
                  return True
      return False

   def isHydrated(self):
      if(self.hydrated == True):
         return True
      if(self.hydrated == None):
         if(self.parent != None):
            for water in self.parent.waters:
               if(self.inWaterContact(water)):
                  self.hydrated = True
                  return True
      self.hydrated = False
      return False

   def hydrationNumber(self):
      count = 0
      if(self.parent != None):
         for water in self.parent.waters:
            if(self.inWaterContact(water)):
               count += 1
      return count

   def isNeighbor(self, other):
      if(self.getChain() == other.getChain()):
         if(abs(self.getIndex() - other.getIndex()) == 1):
            return True
      return False

   def getAtoms(self):
      return self.atoms

   def getChain(self):
      return self.chain

   def __str__(self):
       return("Residue of type %s with %d atoms" % (self.type,
                                                    len(self.atoms)))

   def getAtom(self, index):
      for a in self.atoms:
         if(a.getIndex() == index):
            return a
      return None
  
   def getSingleAtom(self, index):
      return self.atoms[index]

   def getAtomByType(self, atype):
      for a in self.atoms:
           if(a.getType() == atype):
               return a
      return None

   def getAtomNum(self):
       return len(self.atoms)

   def getType(self):
       return self.type

   def assignAtom(self, atom):
       self.atoms.append(atom)

   def getCA(self):
      return self.getAtomByType("CA")

   def getSA(self):
      sa = 0
      for a in self.atoms:
         sa += a.getSA()
      return sa

   def isSurf(self):
      if(not (self.getType() in aaSA.keys())):
         return False
      return isSurf(self.getSA() / aaSA[self.getType()])

   def getIndex(self):
       return self.index

   def ResIntersect(self, otherRes, centerPt, cutoffdist):
       for j in range(0,self.getAtomNum()):
             coord1 = self.getSingleAtom(j).getCoord()
             for n in range(0, otherRes.getAtomNum()):
                   coord2 = otherRes.getSingleAtom(n).getCoord()
                   t = ((coord1[0]-coord2[0]) * (coord1[0]-centerPt[0]) + (coord1[1]-coord2[1]) * (coord1[1]-centerPt[1]) + (coord1[2]-coord2[2]) * (coord1[2]-centerPt[2])) / ((coord1[0]-centerPt[0])**2 + (coord1[1]-centerPt[1])**2 + (coord1[2]-centerPt[2])**2)
                   if 0 < t < 1:
                   #print (t)
                       coord3 = [(coord1[0] * (1-t) + t * centerPt[0]), (coord1[1] * (1-t) + t * centerPt[1]), (coord1[2] * (1-t) + t * centerPt[2])]
                       d = 0
                       for (x1,x2) in zip(coord2, coord3):
                           d += (x1 - x2)**2
                       dist = math.sqrt(d)
                       print("Chain %s Res %s Atom %s with Chain %s Res %s Atom %s, distance is %s" %(self.getChain(), self.getIndex(), self.getSingleAtom(j), otherRes.getChain(), otherRes.getIndex(), otherRes.getSingleAtom(n), dist))
                       #print(dist)
                       #return dist
                       
                       if dist <= cutoffdist:
                          #print('Outside')
                          return("Outside")#, "Residue on Chain %s Number %s Type %s" % (self.getChain(), self.getresNum(), self.getType()))
                          break
             #if dist <= cutoffdist:
                #break
       #return("Inside")#, "Residue on Chain %s Number %s Type %s" % (self.getChain(), self.getresNum(), self.getType()))
                                  

class Protein:
    def __init__(self):
        self.residues = []
        self.waters = []
        self.id = ""
        self.chains = None
        self.seq = []
        
    def __len__(self):
        return len(self.residues)

    def setType(self, pid):
        self.id = pid

    def getGapCount(self):
       gaps = 0
       index = 0

       while(index < len(self.residues)):
          start = self.residues[index].getIndex()
          chain = self.residues[index].getChain()
          index_s = index
          while(index < len(self.residues) and chain == self.residues[index].getChain()):
             index += 1
          gaps += (self.residues[index - 1].getIndex() - start) - (index - index_s - 1)

       return gaps
       
    def getName(self):
       return self.id

    def getChains(self):
       if(self.chains == None):
          self.chains = []
          for r in self.residues:
             if(not r.getChain() in self.chains):
                self.chains.append(r.getChain())

       return self.chains

    def getSequence(self, chain):
       seq = []
       for r in self.residues:
          if(r.getChain() == chain):
             seq.append(conversion[r.getType()])
       seq.append("X")
       return seq
             
    def calMolWeight(self):
       sum = 0
       for i in range(len(self.seq)):
          sum = sum + float(MW[str(self.seq[1])])
       return sum

    def getSurfaceCounts(self):
       scounts = aaSA.copy()
       for k in scounts.keys():
          scounts[k] = 0
       for r in self.residues:
          if(r.isSurf()):
             if(r.getType() in scounts.keys()):
                scounts[r.getType()] += 1
             else:
                print "Skipping type %s" % r.getType()

       return scounts

    def removeBackboneAttribute(self):
       for r in self.residues:
          for a in r.getAtoms():
             a.backbone = False

    def fromData(self, residues, pid):
        self.residues = residues
        self.id = pid

    def addResidue(self, r):
        self.residues.append(r)
        r.setParent(self)

    def addWater(self, w):
       self.waters.append(w)

    def getWaterNumber(self):
       return len(self.waters)

    def getResidues(self,r):
        return self.residues[r]

    def getResidue(self, rindex, chain):
       rindex = self.getResidueIndex(rindex, chain)
       if(rindex != None):
          return self.residues[rindex]

       return None

    def getResidueIndex(self, rindex, chain):
       for i,r in zip(range(len(self.residues)),self.residues):
          if(r.getIndex() == rindex and r.getChain() == chain):
             return i
       return None


    def getLocalAlignment(self, resindex, dist):
        ca = self.residues[resindex].getCA()
        chain = self.residues[resindex].getChain()
        la = []
        for ri in range(len(self.residues)):
            oca = self.residues[ri].getCA()
            if(ri != resindex and chain == self.residues[ri].getChain()):
                if(ca.distance(oca) < dist):
                    la.append(ri)
        return la
              

    def writeCSVAtoms(self, filename, header=False):
       lines = []
       if(header):
          lines.append("pdb_id, atom_index, atom_type, res_type, res_type_sh, chain, res_index, x, y, z, occupancy, beta, atom_surface_area, res_surface_area, res_surface_area_ratio\n")

       for r in self.residues:
          for a in r.getAtoms():
             lines.append(self.id)
             lines[-1] += "," + str(a.getIndex())
             lines[-1] += "," + a.getType()
             lines[-1] += "," + r.getType()
             lines[-1] += "," + convRID(r.getType())
             lines[-1] += "," + r.getChain()
             lines[-1] += "," + str(r.getIndex())
             lines[-1] += "," + str(a.getCoord()[0])
             lines[-1] += "," + str(a.getCoord()[1])
             lines[-1] += "," + str(a.getCoord()[2])
             lines[-1] += "," + str(a.getOccupancy())
             lines[-1] += "," + str(a.getBeta())
             lines[-1] += "," + str(a.getSA())
             lines[-1] += "," + str(r.getSA())
             try:
                lines[-1] += ",%6.4f" % (r.getSA() / aaSA[r.getType()])
             except KeyError:
                print "Could not find residue %s in SA dataset" % r.getType()
                lines[-1] += ",NULL"
             lines[-1] += "\n"
             
       with open(filename, 'w') as f:
          f.writelines(lines)

    def writeCSVResidues(self, filename, header=False):
       lines = []
       if(header):
          lines.append("pdb_id, res_index, res_type, res_type_sh, chain, res_surface_area, res_surface_area_ratio, phi, psi\n")

       for i,r in zip(range(len(self.residues)), self.residues):
          lines.append(self.id)
          lines[-1] += ","  + str(r.getIndex())
          lines[-1] += "," + r.getType()
          lines[-1] += "," + convRID(r.getType())
          lines[-1] += "," + r.getChain()
          lines[-1] += "," + str(r.getSA())
          try:
             lines[-1] += ",%6.4f" % (r.getSA() / aaSA[r.getType()])
          except KeyError:
             print "Could not find residue %s in SA dataset" % r.getType()
             lines[-1] += ",NULL"
          if(i == 0):
             lines[-1] += ",NULL"
          else:
             try:
                lines[-1] += ",%6.4f" % (self._getPhi(i))
             except:
                lines[-1] += ",NULL"
          if(i + 1 == len(self.residues)):
             lines[-1] += ",NULL"
          else:
             try:
                lines[-1] += ",%6.4f" % (self._getPsi(i))
             except:
                lines[-1] += ",NULL"
          lines[-1] += "\n"

       with open(filename, 'w') as f:
          f.writelines(lines)
          
          
          

    def getTorsion(self, indices, types):
       v1 = self._makeVectorFromTypes(indices[0], indices[1], types[0], types[1])
       v2 = self._makeVectorFromTypes(indices[1], indices[2], types[1], types[2])
       v3 = self._makeVectorFromTypes(indices[2], indices[3], types[2], types[3])
       return self._dihedralAngle(v1,v2, v3)
    
    def _getChis(self, rindex):
       result = []
       for i in chis[self.getResidue(rindex).getType()]:
          result.append(self.getTorsion([rindex, rindex, rindex, rindex], chiTypes[i]))
       return result

    def _getPhi(self, rindex):
       if(rindex == 0):
          raise IndexError('cannot use N-terminus residue for Phi angle')       
       return self.getTorsion([rindex - 1, rindex, rindex, rindex], ["C", "N", "CA", "C"])
         
    def _getPsi(self, rindex):
       if(rindex + 1 == len(self.residues)):
          raise IndexError('cannot use C-terminus residue for Psi angle')
       return self.getTorsion([rindex, rindex, rindex, rindex + 1], ["N", "CA", "C","N"])
         
    def getChis(self, rindex, chain):
       self._getChis(self.getResidueIndex(rindex, chain))
       
    
    def getPhi(self, rindex, chain):
       return self_getPhi(self.getResidueIndex(rindex, chain))

    def getPsi(self, rindex, chain):
       return self._getPsi(self.getResidueIndex(rindex, chain))

    def _makeVectorFromTypes(self, r1, r2, t1, t2):
       try:
          v = (self.residues[r2].getAtomByType(t2).getCoord()[0] - self.residues[r1].getAtomByType(t1).getCoord()[0],
             self.residues[r2].getAtomByType(t2).getCoord()[1] - self.residues[r1].getAtomByType(t1).getCoord()[1],
             self.residues[r2].getAtomByType(t2).getCoord()[2] - self.residues[r1].getAtomByType(t1).getCoord()[2])
       except:
          if(self.residues[r1].getAtomByType(t1) == None):
             raise RuntimeError("Missing Atom type %s in residue %s" %( t1, self.residues[r1].__str__()))
          if(self.residues[r2].getAtomByType(t2) == None):
             raise RuntimeError("Missing Atom type %s in residue %s" %( t2, self.residues[r2].__str__()))
          

       return v

    def _dot(self, v1,v2):
       return (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])

    def _cross(self, v1, v2):
       return (v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0])

    def _dihedralAngle(self, v1, v2, v3):

       x = math.sqrt(self._dot(v2,v2)) * (self._dot(v1, self._cross(v2,v3)))
       y = self._dot(self._cross(v1,v2), self._cross(v2,v3))

       altAngle = math.atan2(x,y)
       
       return altAngle
       

class sqMatrix:

    def __init__(self, dim):
        self.dim = dim
        self.matrix = [None] * self.dim * self.dim
    
    def __len__(self):
        return self.dim

    def __str__(self):
        result = []
        for i in range(self.dim):
            result.append("|  ")
            for j in range(self.dim): 
                result.append("%05.3e " % self.get(i,j))
            result.append(" |\n")
        return "".join(result)

    def get(self, row, col):
        return self.matrix[row * self.dim + col]

    def set(self, row, col, x):
        self.matrix[row * self.dim + col] = x
    
    def setRow(self, row, array):
        self.matrix[row:(row + dim)] = array


def readProtein(pdbfile):

    prot = Protein()
    prot.setType(pdbfile.split("/")[-1].split(".")[0])
    currentR = Residue()
    proteinnum = [0,-1]
    atomn = 0

    with open(pdbfile, "r") as f:  #open the file
        for line in f.readlines():       #read file line by line
            m = regexp.match(line)       #matches the regular expression line by line

            if(m): 
                proteinnum[0] = int(m.group(5))  #assign current residue number
                atomn = atomn+1
                if(proteinnum[0] != proteinnum[1]):  #check residue number              
                    currentR = Residue()               #initialize a new residue
                    proteinnum[1] = proteinnum[0]        #assign new residue number
                    currentR.fromData([], m.group(3), m.group(4), int(m.group(5)))  #initialization
                    prot.addResidue(currentR)              #add new residue to protein

            
                atoml = Atom()          #import the atom
                atoml.fromData(((float(m.group(6)),float(m.group(7)),float(m.group(8)))), m.group(2), int(m.group(1)), float(m.group(9)), float(m.group(10)))
                currentR.assignAtom(atoml)   #put current atom into current residue                
            
            else:
               #check if it is a water
               m = regexpWater.match(line)
               if(m):
                  atomw = Atom()
                  atomw.fromData((float(m.group(2)), float(m.group(3)), float(m.group(4))), 'HOH', -1, 0, 0)
                  prot.addWater(atomw)

    return prot

def readProteinSeq(pdbfile):
    prot = Protein()
    with open(pdbfile, "r") as f:
        for line in f.readlines():
            m = regexpSeq.match(line)

            if(m):
               for i in range(len(m.group(1).split())):
                    if m.group(1).split()[i] in conversion:
                       prot.seq.append(conversion[m.group(1).split()[i]])   
    return prot           


def readProteinSA(pdbfile):

    prot = Protein()
    prot.setType(pdbfile.split("/")[-1].split(".")[0])
    currentR = Residue()
    proteinnum = [0,-1]
    atomn = 0

    with open(pdbfile, "r") as f:  #open the file
        for line in f.readlines():       #read file line by line
            m = regexpSA.match(line)       #matches the regular expression line by line

            if(m): 
                proteinnum[0] = int(m.group(5))  #assign current residue number
                atomn = atomn+1
                if(proteinnum[0] != proteinnum[1]):  #check residue number              
                    currentR = Residue()               #initialize a new residue
                    proteinnum[1] = proteinnum[0]        #assign new residue number
                    currentR.fromData([], m.group(3), m.group(4), int(m.group(5)))  #initialization
                    prot.addResidue(currentR)              #add new residue to protein

            
                atoml = Atom()          #import the atom
                atoml.fromData(((float(m.group(6)),float(m.group(7)),float(m.group(8)))), m.group(2), int(m.group(1)), float(m.group(9)), float(m.group(10)))
                try:
                   atoml.setSA(float(m.group(11)))
                except:
                   print "Failure in reading"
                   print m.groups()
                   print line
                   exit()
                currentR.assignAtom(atoml)   #put current atom into current residue

            else:
               #check if it is a water
               m = regexpWater.match(line)
               if(m):
                  atomw = Atom()
                  atomw.fromData((float(m.group(2)), float(m.group(3)), float(m.group(4))), 'HOH', -1, 0, 0)
                  prot.addWater(atomw)

    return prot


def makePDB(var1, var2, pdbfile, outfile, chains):

     #extract pdb name
   name = pdbfile.split(".")[0]
   index = 0
   res = ""

   with open(pdbfile, 'r') as p:
      with open(outfile, 'w') as o:
         for line in p.readlines():
            m = regexp.match(line)
            if(m and int(m.group(5)) < len(var1) and m.group(4) in chains):
               index = int(m.group(5))
               o.write(line[:54])
               o.write("%6.2f" % var1[index])
               o.write("%6.2f" % var2[index])
               o.write(line[66:])
            else:
               o.write(line)

def convRID(r):
   if(len(r) == 3):
      if r in conversion:
         return conversion[r]
      else:
         return "Z"
   else:
      for k,v in conversion.items():
         if(v == r):
            return k
   return None

def chiNumber(resType):
   return len(chis[resType])

def setCutoff(cut):
   global surfCutoff
   surfCutoff = cut

def isSurf(ratio):
   if(surfCutoff == 0):
      return True
   elif(surfCutoff < 0):
      if(ratio < abs(surfCutoff)):
         return True
   elif(ratio > surfCutoff):
      return True

   return False
