movie <- function() {
  path = ''
  ICs = makeICs(shape='cube',stretch=c(1,2.1,1.5),v=c(-0.5,0,5),r=c(0,0,5),L=c(6,3.09,2),u=0.1)
  x = simulation(ICs,0.01,nIterationsMax=1e4)
  summarizeSimulation(x)
  showSimulation(x,showVectors=T,path=path)
}

checkEnergy <- function() {
  ICs = makeICs(stretch=c(2,4,6),k=0.7,u=0.3,g=c(0,0,-9.81),r=c(0,0,10),v=rnorm(3),A=diag(3),L=c(10,10,10))
  x = simulation(ICs,0.1)
  Eold = 1e99
  E = array()
  t = array()
  for (i in seq(length(x$dyn))) {
    y = x$dyn[[i]]
    Etra = t(y$v)%*%y$v/2
    Epot = -t(x$init$g)%*%y$r
    Erot = t(y$w)%*%y$L/2
    #print(c(x$t,Epot,Etra,Erot,Epot+Etra,Epot+Etra+Erot))
    E[i] = Etra+Epot+Erot
    t[i] = y$t
    if ((E[i]-Eold)/Eold>1e-5) {
      print("error") 
    }
    Eold = E[i]
  }
  plot(t,E,col='blue',log='y',type='l')
}

checkTime <- function() {
  
  ICs = makeICs(stretch=c(2,4,6),k=0.7,u=0.3,g=c(0,0,-9.81),r=c(0,0,10),v=c(0,0,0),A=diag(3),L=c(10,10,10))
  x = simulation(ICs,dtOutput=1e99,tmax=5,maxRotationAngle=0.01)
  summarizeSimulation(x)
  showSimulation(x,length(x$dyn))
  
  # load reference calculation and/or save more accurate reference calculation
  filename = '/Users/do/Dropbox/Code/R/cuboid/reference.img'
  fileExists = !is.na(file.info(filename)[[1]])
  if (fileExists) {
    load(file=filename)
    if (x$stats$numberOfIterations>ref$stats$numberOfIterations) {
      ref = x
      save(ref,file=filename)
    }
  } else {
    ref = x
    save(ref,file=filename)
  }
  i = 2
  dr = fct.vectNorm(x$dyn[[i]]$r-ref$dyn[[i]]$r)/fct.vectNorm(ref$dyn[[i]]$r)
  dA = fct.vectNorm(x$dyn[[i]]$A-ref$dyn[[i]]$A)/fct.vectNorm(ref$dyn[[i]]$A)
  dv = fct.vectNorm(x$dyn[[i]]$v-ref$dyn[[i]]$v)/fct.vectNorm(ref$dyn[[i]]$v)
  dL = fct.vectNorm(x$dyn[[i]]$L-ref$dyn[[i]]$L)/fct.vectNorm(ref$dyn[[i]]$L)
  cat('r-vector: ',dr,'\n')
  cat('A-vector: ',dA,'\n')
  cat('v-vector: ',dv,'\n')
  cat('L-vector: ',dL,'\n')
  cat('TOTAL:    ',fct.vectNorm(c(dr,dA,dv,dL)),'\n')
}