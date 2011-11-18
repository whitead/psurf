#!/usr/bin/python

import numpy.random, math, sys

def printHelp():   #system help for the program
    print "sto_annealing.py [total step size] [annealing step size] [annealing cycle]"
    exit(0)

if( len(sys.argv) != 4 ):
  printHelp()

a = [1]*20                              # Initialization of the parameter for dirichlet distribution and randomness
alpha = [0.8]

ss = int(sys.argv[1])                   # Take in step information from the user input
anss = int(sys.argv[2])
cycle = int(sys.argv[3])

countMatrix = open("countMatrix.txt",'r')          # Read in all the required information
countMatrix = str(countMatrix.read()).split()
groEL_Close = open("GroCl.txt",'r')
groEL_Close = str(groEL_Close.read()).split()
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

# Begining of the MC program
with open("distRec_test.txt","w") as f:
  with open("potRec_test.txt","w") as g:
    g.write("StepNum Potential Temp\n")
    f.write("StepNum CYS ASP SER GLN LYS ILE PRO THR PHE ALA GLY HIS GLU LEU ARG TRP VAL ASN TYR MET\n")

    # The probablity Model
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
          t1 = float(countMatrix[i * 20 + j]) * float(pBuried[i * 20 + j]) * float(contact[j])  #need the index for pBuried here
          t2 = float(hydration[i]) / hydrationsum * float(hydration[j]) * ini[j]
          Y[j] = t1 + t2
        sigmaY = math.log(sum(Y))
        temp = 0
        for k in range(903) :
          temp = temp + float(pBuried[20 * k + i])
        X[i] = temp / 903 * sigmaY    
      sigmaXln[1] = sum(X) * (1-surFrac)

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
      sigmaXln[2] = sum(X) * (surFrac)

      return(sigmaXln[0] + sigmaXln[1] - sigmaXln[2])

    
    # The annealing alogorithm
    def MH(step,cycle,pot,potP,R):     #pot is current potential, potP is Previous potential
                                       #Temp = [0.2,0.1,0.05,0.02,0.01,0.005,0.001,0.0005,0.001,0.005,0.01,0.02,0.05,0.1]
#      Tscale = 10
      Tscale = 1  
      Thot = 1 * Tscale
#      Tcold = 0.5*10**-3 * Tscale
      Tcold = 1 * Tscale
      StepSize = (math.log(Thot)-math.log(Tcold)) / (step / (cycle * 2))   # define stepsize for temperature, and all the temperature values
      Temp1 = [0] * (step / (cycle * 2))
      Temp2 = [0] * (step / (cycle * 2))
      for n in range(step / (cycle * 2)) :
          Temp1[n] = math.exp(math.log(Thot) - StepSize * n)
          Temp2[n] = math.exp(math.log(Tcold) + StepSize * n)
      Temp = Temp1 + Temp2
      #print(Temp)

      #Randsize = (math.log(0.9) - math.log(0.01)) / (step / (cycle * 2))  # define stepsize for randomness, and all the randomness values corresponding to different temperature
      #Rand1 = [0] * (step / (cycle * 2))
      #Rand2 = [0] * (step / (cycle * 2))
      #for n in range(step / (cycle * 2)) :
      #    Rand1[n] = math.exp(math.log(0.9) - Randsize * n)
      #    Rand2[n] = math.exp(math.log(0.01) + Randsize * n)
      #Rand = Rand1 + Rand2
      
      tP = pot                      #Initialize the annealing alogorithm
      T = Temp * cycle
      #Rand = Rand * cycle
      MH = [0] * len(T)
      accept = 0
      total = 0
      totalcyc = 0
      frac = 0.3
      Rs = R

      for j in range(0,len(T)):
          MH[j] = -50*(pot - potP) / T[j]
          MH[j] = math.exp(MH[j])            #evaluate the MH value, pot(new) - potP(old)
          r = numpy.random.random()
          #print(pot,potP)
          #print(r)
          total = total + 1
          #print(T[j])
          if(r > MH[j]):
              tR = []
              for x in R:
                  tR.append(str(x))
                  tR.append(" ")
              tR = "".join(tR)
              f.write("%s %s\n" %(j+1,tR))
              g.write("%s %s %s\n" %(j+1,str(potP),T[j]))
              Rnew = numpy.random.mtrand.dirichlet(a,size=1).tolist()[0]
              for k in range(20):
                  R[k] = Rs[k] * (1-float(alpha[0])) + Rnew[k] * float(alpha[0])
              pot = potential(R)
              
          else:
              accept = accept + 1
              tR = []
              for x in R:
                  tR.append(str(x))
                  tR.append(" ")
              tR = "".join(tR)
              f.write("%s %s\n" %(j+1,tR))
              g.write("%s %s %s\n" %(j+1,str(potP),T[j]))
              Rs = R
              Rnew = numpy.random.mtrand.dirichlet(a,size=1).tolist()[0]
              for k in range(20):
                  R[k] = R[k] * (1-float(alpha[0])) + Rnew[k] * float(alpha[0])
              potP = pot
              pot = potential(R)
          if total == (step / cycle) :
              currenttotal = totalcyc * step / cycle
              frac = (frac * currenttotal +  float(accept)) / (currenttotal + float(total))
              accept = 0
              total = 0
              totalcyc = totalcyc+1
              #print(frac)
      print(frac)

#for i in range(20) :
#  if float(groEL_Close[i]) == 0 :
#    groEL_Close[i] = 0.00001
#
#print(potential(groEL_Close))


    pot = [0] * ss                       #big steps, ss value is small for now
    R = numpy.random.mtrand.dirichlet(a, size=1).tolist()[0]     #draw initial guess from dirichlet distribution

##    the following three lines are using the GroEL inside surface distribution as the input (for obtaining reference values)
##    c = []
##    for x in range(len(groEL_Close)):
##        c.append(float(groEL_Close[x]))

    for i in range(ss) :     
        pot[i] = potential(R)
        print("\r%s / %s" %(i+1, ss))
        if i > 0 and pot[i] > pot[i-1]:     # for minimization, plug in current Potential(pot) and previous potential(potP) in to MH program
            MH(anss,cycle,pot[i],pot[i-1],R)
            a = [1] * 20
        else:
            tR = []
            Rnew = numpy.random.mtrand.dirichlet(a,size=1).tolist()[0]
            for x in R:
                tR.append(str(x))
                tR.append(" ")
            tR = "".join(tR)
            #f.write("%s\n" %(tR))
            #g.write("%s\n" %(str(pot[i])))
            for j in range(20):
                R[j] = R[j] * (1-float(alpha[0])) + Rnew[j] * float(alpha[0])
          
        

  
    










