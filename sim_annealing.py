#!/usr/bin/python

import numpy.random
import math


a = [1]*20
R = numpy.random.mtrand.dirichlet(a, size=1)


countMatrix = open("countMatrix.txt",'r')
countMatrix = str(countMatrix.read()).split()
groEL_Close = open("GroCl.txt",'r')
groEL_Open = open("GroOp.txt",'r')
surFrac = open("surfrac.txt",'r')
surFrac = float(str(surFrac.read()))
pSurf = open("psurf.txt",'r')
pSurf = str(pSurf.read()).split()
pBuried = open("pburied.txt",'r')
pBuried = str(pBuried.read()).split()
pUnfold = open("punfold.txt",'r')
pUnfold = str(pUnfold.read()).split()
contact = open("contact.txt",'r')
contact = str(contact.read()).split()
hydration = open("hydration.txt",'r')
hydration = str(hydration.read()).split()
hydrationsum = 0
for i in range(20) :
  hydrationsum = hydrationsum + float(hydration[i])



def potential(ini):
  sigmaXln = [0,0,0]
  sigmaY = [0]

  X = [0] * 20
  for i in range(20) :
     Y = [0] *  len(ini)
     for j in range(20) :
       t1 = float(countMatrix[i * 20 + j]) * ini[j] * float(contact[j])
       t2 = float(hydration[i]) / hydrationsum * float(hydration[j]) * ini[j]
       Y[j] = t1 + t2      
     sigmaY = math.log(sum(Y))
     temp = 0
     for k in range(903) :
       temp = temp + float(pSurf[20 * k + i])
     X[i] = temp / 903 * sigmaY    
  sigmaXln[0] = sum(X) * surFrac

  X = [0] * 20
  for i in range(20) :
     Y = [0] *  len(ini)
     for j in range(20) :
       t1 = float(countMatrix[i * 20 + j]) * float(pBuried[]) * float(contact[j])  #need the index for pBuried here
       t2 = float(hydration[i]) / hydrationsum * float(hydration[j]) * ini[j]
       Y[j] = t1 + t2      
     sigmaY = math.log(sum(Y))
     temp = 0
     for k in range(903) :
       temp = temp + float(pBuried[20 * k + i])
     X[i] = temp / 903 * sigmaY    
  sigmaXln[1] = sum(X) * surFrac

  X = [0] * 20
  for i in range(20) :
     Y = [0] *  len(ini)
     for j in range(20) :
       t1 = float(countMatrix[i * 20 + j]) * ini[j] * float(contact[j])
       t2 = float(hydration[i]) / hydrationsum * float(hydration[j]) * ini[j]
       Y[j] = t1 + t2      
     sigmaY = math.log(sum(Y))
     temp = 0
     for k in range(903) :
       temp = temp + float(pUnfold[20 * k + i])
     X[i] = temp / 903 * sigmaY    
  sigmaXln[0] = sum(X) * (1-surFrac)

  return(sigmaXln[0] + sigmaXln[1] - sigmaXln[2])


def MH(T,ini) :



#for i in range(1000000) :
ini = R[0]
potential(ini)
    










