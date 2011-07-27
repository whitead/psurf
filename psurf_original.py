import math, re



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
regexp = re.compile("ATOM[\s]*([\d]*)[\s]*([\w]{1,3}[\d]{0,2}[\']{0,2})[\s]{0,3}([\w]{1,4})[\s]*([\w]*)[\s]{0,5}([\d]{1,5})[\s]*([-\d\.]{3,7})[\s]*([-\d\.]*)[\s]*([-\d\.]*)[\s]*([\d\.]*)[\s]*([-\d\.]*)[\s]*(.*)")

#added SA column
regexpSA = re.compile("ATOM[\s]*([\d]*)[\s]*([\w]{1,3}[\d]{0,2}[\']{0,2})[\s]{0,3}([\w]{1,4})[\s]*([\w]*)[\s]{0,5}([\d]{1,5})[\s]*(-{0,1}[\d\.]{3,7})[\s]*(-{0,1}[\d\.]{3,7})[\s]*([-{0,1}\d\.]{3,7})[\s]*([\d\.]*)[\s]*([-\d\.]*)[\s]*([-\d\.]*)")

conversion = {'ALA':'A', 'ARG':'R', 'ASP':'D', 'ASN':'N', 'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}

chis = {'GLU':[0,1], 'LYS':[0,1,2,3], 'SER':[], 'GLY':[], 'PRO':[]}
chiTypes = [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD'], ['CB', 'CG', 'CD' , 'CE'], ['CB', 'CG', 'CD' , 'NZ']]

peptideSurfData = "/home/wenjunh/Documents/ProteinSurfaces/aalist.txt"
aaSA = {}
surfCutoff = 0.3

with open(peptideSurfData, 'r') as f:
   for line in f.readlines():
      sline = line.split()
      try:
         aaSA[sline[0]] = float(sline[1])
      except:
         pass

class Atom: 
   def __init__(self):
      self.coord = (0,0,0)
      self.type = "CA"
      self.index = 0

   def fromData(self, coordinates, atomType, atomIndex):
      self.coord = coordinates
      self.type = atomType
      self.index = atomIndex

   def __str__(self):
      return("Type: \"%s\"; Index: %03d; Coordinates: (%05.3f, %05.3f, %05.3f)" % (self.type, self.index, self.coord[0], self.coord[1], self.coord[2]))
   
   def setSA(self, sa):
      self.sa = sa

   def getSA(self):
      return self.sa

   def distance(self, otherAtom):
      dist = 0
      for (x1,x2) in zip(self.coord, otherAtom.coord):
         dist += (x1 - x2)**2
      return math.sqrt(dist)

   def centerdistance(self, center):
      dist = 0
      for (x1,x2) in zip(self.coord, (center[0],center[1],center[2])):
         dist += (x1 - x2)**2
      return math.sqrt(dist)

   def getType(self):
      return self.type

   def getIndex(self):
      return self.index

   def getCoord(self):
      return self.coord

class Residue:
   def __init__(self):
       self.atoms = []
       self.type = "X"

   def fromData(self, atoms, rType, chain, resNum):
       self.type = rType
       self.atoms = atoms
       self.chain = chain
       self.resNum = resNum

   def getChain(self):
      return self.chain

   def __str__(self):
       return("Residue of type %s with %d atoms" % (self.type,
len(self.atoms)))

   def getAtom(self, index):
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
      return isSurf(self.getSA() / aaSA[self.getType()])

   def getresNum(self):
       return self.resNum

   def ResIntersect(self, otherRes, centerPt, cutoffdist):
       for j in range(0,self.getAtomNum()):
             coord1 = self.getAtom(j).getCoord()
             for n in range(0, otherRes.getAtomNum()):
                   coord2 = otherRes.getAtom(n).getCoord()
                   t = ((coord1[0]-coord2[0]) * (coord1[0]-centerPt[0]) + (coord1[1]-coord2[1]) * (coord1[1]-centerPt[1]) + (coord1[2]-coord2[2]) * (coord1[2]-centerPt[2])) / ((coord1[0]-centerPt[0])**2 + (coord1[1]-centerPt[1])**2 + (coord1[2]-centerPt[2])**2)
                   if 0 < t < 1:
                   #print (t)
                       coord3 = [(coord1[0] * (1-t) + t * centerPt[0]), (coord1[1] * (1-t) + t * centerPt[1]), (coord1[2] * (1-t) + t * centerPt[2])]
                       d = 0
                       for (x1,x2) in zip(coord2, coord3):
                           d += (x1 - x2)**2
                       dist = math.sqrt(d)
                       print("Chain %s Res %s Atom %s with Chain %s Res %s Atom %s, distance is %s" %(self.getChain(), self.getresNum(), self.getAtom(j), otherRes.getChain(), otherRes.getresNum(), otherRes.getAtom(n), dist))
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
        self.id = ""
        self.chains = None
        
    def __len__(self):
        return len(self.residues)

    def setType(self, pid):
        self.id = pid

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
             
             

    def fromData(self, residues, pid):
        self.residues = residues
        self.id = pid

    def addResidue(self, r):
        self.residues.append(r)

    def getResidues(self):
        return self.residuex

    def getResidue(self, r):
        return self.residues[r]

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
              
    def getTorsion(self, indices, types):
       v1 = self._makeVectorFromTypes(indices[0], indices[1], types[0], types[1])
       v2 = self._makeVectorFromTypes(indices[1], indices[2], types[1], types[2])
       v3 = self._makeVectorFromTypes(indices[2], indices[3], types[2], types[3])
       return self._dihedralAngle(v1,v2, v3)
                      
    def getChis(self, rindex):
       result = []
       for i in chis[self.getResidue(rindex).getType()]:
          result.append(self.getTorsion([rindex, rindex, rindex, rindex], chiTypes[i]))
       return result
       
  
    def getPhi(self, rindex):
       if(rindex == 0):
          raise IndexError('cannot use N-terminus residue for Phi angle')       
       return self.getTorsion([rindex - 1, rindex, rindex, rindex], ["C", "N", "CA", "C"])

    def getPsi(self, rindex):
       if(rindex == len(self.residues)):
          raise IndexError('cannot use C-terminus residue for Psi angle')
       return self.getTorsion([rindex, rindex, rindex, rindex + 1], ["N", "CA", "C","N"])

    def _makeVectorFromTypes(self, r1, r2, t1, t2):
       try:
          v = (self.residues[r2].getAtomByType(t2).getCoord()[0] - self.residues[r1].getAtomByType(t1).getCoord()[0],
             self.residues[r2].getAtomByType(t2).getCoord()[1] - self.residues[r1].getAtomByType(t1).getCoord()[1],
             self.residues[r2].getAtomByType(t2).getCoord()[2] - self.residues[r1].getAtomByType(t1).getCoord()[2])
       except:
          if(self.residues[r1].getAtomByType(t1) == None):
             raise RuntimeError("Missing Atom type %s in residue %s", t1, self.residues[r1].__str__())
          if(self.residues[r2].getAtomByType(t2) == None):
             raise RuntimeError("Missing Atom type %s in residue %s", t2, self.residues[r2].__str__())
          

       return v

    def _dot(self, v1,v2):
       return (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])

    def _cross(self, v1, v2):
       return (v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0])

    def _dihedralAngle(self, v1, v2, v3):

       #u = self._cross(v1,v2)
       #v = self._cross(v2,v3)
       #w = self._cross(u,v)
       #angle = math.acos(self._dot(u,v) / (math.sqrt(self._dot(u,u)) *  math.sqrt(self._dot(v,v))))
       #if(math.acos(self._dot(v2,w) / (math.sqrt(self._dot(v2,v2)) *  math.sqrt(self._dot(w,w)))) > 0.001):
          #angle = -angle
       
          

        

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
    prot.setType(pdbfile)
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
                    currentR.fromData([], m.group(3), m.group(4), m.group(5))  #initialization
                    prot.addResidue(currentR)              #add new residue to protein

            
                atoml = Atom()          #import the atom
                atoml.fromData(((float(m.group(6)),float(m.group(7)),float(m.group(8)))), m.group(2), int(m.group(1)))
                currentR.assignAtom(atoml)   #put current atom into current residue                
    return prot

def readProteinSA(pdbfile):

    prot = Protein()
    prot.setType(pdbfile)
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
                    currentR.fromData([], m.group(3), m.group(4), m.group(5))  #initialization
                    prot.addResidue(currentR)              #add new residue to protein

            
                atoml = Atom()          #import the atom
                atoml.fromData(((float(m.group(6)),float(m.group(7)),float(m.group(8)))), m.group(2), int(m.group(1)))
                try:
                   atoml.setSA(float(m.group(11)))
                except:
                   print "Failure in reading"
                   print m.groups()
                   print line
                   exit()
                currentR.assignAtom(atoml)   #put current atom into current residue
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
      return conversion[r]
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
   if(surfCutoff < 0):
      if(ratio < abs(surfCutoff)):
         return True
   elif(ratio > surfCutoff):
      return True

   return False
