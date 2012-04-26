

import re, math, scipy.optimize, random

#xyz Format
xyzSliceIndices = [1,10,18,25,33,40,48]
AtomName = 0
AtomCoords = (2, 4, 6)

backbone = ['N', 'O', 'C']
atomlist = ['N', 'O', 'C', 'H']

class Atom:
    def __init__(self):
        self.coord = (0,0,0)
        self.type = "H"
        self.index = -1
        self.backbone = False
    
    def getIndex(self):
        return self.Index

    def fromData(self, coordinates, atomType, atomIndex):
        self.coord = coordinates
        self.type = atomType
        if(self.type in backbone):
            self.backbone = True
        self.index = atomIndex
        
    def distance(self, otherAtom):
        return math.sqrt(self.distanceSqr(otherAtom))

    def distanceSqr(self, otherAtom):
        dist = 0
        for (x1,x2) in zip(self.coord, otherAtom.coord):
            dist += (x1 - x2)**2
        return dist

    def setParent(self, s):
        self.partent = s

    def __str__(self):
        return("Type: \"%s\"; Index: %03d; Coordinates: (%05.3f, %05.3f, %05.3f)" % (self.type, self.index, self.coord[0], self.coord[1], self.coord[2]))

    def getType(self):
        return self.type

    def getIndex(self):
        return self.index
    
    def getCoord(self):
        return self.coord

class Molecule:
    def __init__(self):
        self.atoms = []
        self.spheres = []

    def setName(self,name):
        self.name = name

    def getAtoms(self):
        return self.atoms

    def getSpheres(self):
        return self.spheres

    def getSingleShpere(self, spherenum):
        return self.spheres[spherenum]

    def getName(self):
        return self.name
    
    def getAtomNums(self):
        return len(self.atoms)

    def addAtom(self, a):
        self.atoms.append(a)

    def addSphere(self, s, spherenum):
        self.spheres.append(s)
        self.spheres.index = spherenum
        

def _slices(s, indices):
    position = 0
    for next in indices:
        yield s[position:next]
        position = next

def _readXYZLine(line):
    sliced = list(_slices(line, xyzSliceIndices))
    result = {'atomName':sliced[AtomName].strip(), 'atomCoords':[float(x) for x in [sliced[x] for x in AtomCoords]]}
    return result

def readMolecule(xyzfile):
    
    mole = Molecule()
    mole.setName(xyzfile.split("/")[-1].split(".")[0])
    atomn = 0

    with open(xyzfile, "r") as f:
        for line in f.readlines():
            if line[0] in atomlist:
                lineData = _readXYZLine(line)
                atomn = atomn+1

                atoml = Atom()
                atoml.fromData(tuple(lineData['atomCoords']), lineData['atomName'],atomn)
                mole.addAtom(atoml)
    
    return mole


def genRadius(num):
    sphereRadius = range(num)
    
    for x in sphereRadius:
        sphereRadius[x] = abs(10 * random.random())
#        sphereRadius[x] = 10
    
    return sphereRadius
    

def sphereVol(sphereRadius, molecule, sphereCenter):
    
    vol = 0
    for x in sphereRadius:
#        vol = vol + 4 * math.pi * math.pow(float(x),3) / 3
        vol += x
    
    penalty1 = 0 #keep the atoms in only one sphere
    
    for x in range(len(molecule.getAtoms())):
        if molecule.getAtoms()[x].getIndex() in sphereCenter:
            continue
        else:
            l = 0
            include = []
            exclude = []
            for y in range(len(sphereCenter)):
                exclude.append(molecule.getAtoms()[x].distance(molecule.getAtoms()[int(sphereCenter[y])-1]) - sphereRadius[y-1])
                if molecule.getAtoms()[x].distance(molecule.getAtoms()[int(sphereCenter[y])-1]) <= sphereRadius[y-1] and molecule.getAtoms()[int(sphereCenter[y])-1].getIndex() != molecule.getAtoms()[x].getIndex():
                    #print(x)
                    include.append(sphereRadius[y-1] - molecule.getAtoms()[x].distance(molecule.getAtoms()[int(sphereCenter[y])-1]))
                    l = l + 1
        
            include = sorted(include)
            #print(x+1)
            #print(include)
            #penalty1 = penalty1 + math.log(math.pow(sum(include[1:]),3)+1)
            penalty1 += 10 * sum([abs(x) for x in include[1:]])
            
            if l == 0:
                #print(exclude)
                #penalty1 = penalty1 + math.log(math.pow(sum(exclude),10)+1)
                penalty1 += 25 * sum([abs(x) for x in exclude])
    
    
#    print(vol)
        

    penalty2 = 0 #keep radius positive

    for x in sphereRadius:
        if x <= 0 : 
            penalty2 = penalty2 + math.pow(10,10)
        
    
    sphereVol = vol + penalty1 + penalty2
#    print(sphereVol)

    return sphereVol

def test(x):
    y = math.pow(x,2) - 2*x + 1
    return(y)

def optRadius(molecule, sphereCenter):
    
    sample = 25

    optRadius = []
    optValue = []
    
    for x in range(sample):
        sphereRadius = genRadius(len(sphereCenter))
        result = scipy.optimize.fmin_cg(sphereVol, sphereRadius, args=(molecule, sphereCenter), full_output = True)
        optRadius.append(result[0])
        optValue.append(result[1])

    optRadius = optRadius[optValue.index(min(optValue))]

    return optRadius
                
            


