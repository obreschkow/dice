#########################################################################################
### Description #########################################################################
#########################################################################################

# HOW TO USE THIS CODE?
#
# This codes runs on the programming language R, although the core routines are written
# in C and precompiled. Follow these steps to get started.
#
# 1) Install the free programming language R.
# 2) Optionally install the free programming environment Rstudio.
# 3) Start Rstudio, or, if only R itself is installed, run R in a terminal
# 4) Before running the dice simulation for the first time, install some R-packages
#    by running this command in R:
#    > install.packages(c("compiler","inline","rgl"))
# 5) Compile the dice simluation by typing
#    > source('[path]/diceSimulation.R')
# 6) Run the default example by typing
#    > toss()
#    or
#    > endlessToss()
# 7) Run other tosses, using different settings, for example:
#
#    bounce on a face
#    > toss(stretch=c(2,2,2),k=0.7,u=0.5,r=c(0,0,6),v=c(0,0,0),L=c(0,0,2),sim_dtOutput=0.05)
#    
#    bounce on an edge
#    > toss(stretch=c(2,2,2),k=0.7,u=0.5,r=c(0,0,6),v=c(0,0,0),L=c(0,0,2),A=fct.rotationMatrix(c(pi/4,0,0))%*%diag(3),sim_dtOutput=0.05)
#    
#    bounce on a corner
#    > toss(stretch=c(1,1,1),k=0.7,u=0.5,r=c(0,0,2),v=c(0,0,0),L=c(0,0,1),A=fct.rotationMatrix(c(0,atan(1/sqrt(2)),0))%*%fct.rotationMatrix(c(45/180*pi,0,0))%*%diag(3),sim_dtOutput=0.02)
#
#    show 'precession'
#    > toss(stretch=c(3,3,5),k=0,u=0,g=c(0,0,0),r=c(0,0,10),v=c(0,0,0),L=c(5,10,5),animation_showVectors=TRUE,sim_nIterationsMax=1000)
#    
#    show 'nutation'
#    > toss(stretch=c(3,4,5),k=0,u=0,g=c(0,0,0),r=c(0,0,10),v=c(0,0,0),L=c(10,10,10),animation_showVectors=TRUE,sim_nIterationsMax=1000)
#
#    stability of shortest and longest axis, instability of intermediate axis
#    > toss(stretch=c(1,3,6),k=0,u=0,g=c(0,0,0),r=c(0,0,10),v=c(0,0,0),L=c(10,3e-2,3e-2),sim_dtOutput=0.05,sim_nIterationsMax=20000,animation_showVectors=TRUE)
#    > toss(stretch=c(1,3,6),k=0,u=0,g=c(0,0,0),r=c(0,0,10),v=c(0,0,0),L=c(3e-2,10,3e-2),sim_dtOutput=0.05,sim_nIterationsMax=20000,animation_showVectors=TRUE)
#    > toss(stretch=c(1,3,6),k=0,u=0,g=c(0,0,0),r=c(0,0,10),v=c(0,0,0),L=c(3e-2,3e-2,10),sim_dtOutput=0.05,sim_nIterationsMax=20000,animation_showVectors=TRUE)
#
#    2D die
#    > toss(stretch=c(0,2,4),k=0.7,u=0.5,r=c(0,0,10),v=c(0,0,0),A=fct.rotationMatrix(c(atan(0.5),0,0))%*%diag(3),L=c(-3,0,0),sim_dtOutput=0.05)
#
#    tetrahedron
#    > toss('tetrahedron',stretch=c(3,3,5))

# AUTHORSHIP
# Algorithm and code developped by Danail Obreschkow, 1998-2015.
# Last modification on 20 March 2015
# Code can be shared under the Apache License Version 2.0, January 2004
# (see http://www.apache.org/licenses/)

# CODE INFORMATION

# Each tossing event is stored in a structured variable x with the following components.

# Required as simulation input
# x$init$k - bouncing coefficient, scalar [0,1]
# x$init$u - friction coefficient, scalar [0,1]
# x$init$g - gravitational acceleration 3-vector [m/s]
# x$init$c - position vectors of vertices (corners) [3,]-matrix
# x$init$I - tensor of inertia in comoving coordinates
# x$init$r - initial position vector of the center of gravity [m]
# x$init$v - initial velocity fector of the center of gravity [m/s]
# x$init$A - initial rotation matrix discribing the orientation
# x$init$L - initial angular momentum vector [kg m^2/s]

# Triangle list x$triangle required only for graphics and face statistics
# x$triangle$vertex - [n,3]-matrix with the indices of the vertices (corners)
# x$triangle$color - n-vector with colors
# x$triangle$face - n-vector with face index (e.g. 1-6 in case of a cuboid)

# Changes made by simulation: dyn converted a list dyn[[]]
# x$dyn[[]]$c - position vectors of vertices (corners) [3,]-matrix
# x$dyn[[]]$r - position vector of the center of gravity [m]
# x$dyn[[]]$v - velocity fector of the center of gravity [m/s]
# x$dyn[[]]$A - rotation matrix discribing the orientation
# x$dyn[[]]$L - angular momentum vector [kg m^2/s]
# x$dyn[[]]$w - angular velocity vector
# x$dyn[[]]$n - number of bounces so far
# x$dyn[[]]$t - current time [s]
# x$stats$numberOfIterations - number of main iterations performed
# x$stats$computationTime - time required to run the simulation [s]
# x$stats$physicalTime - simulated physical time [s]
# x$stats$iterationsPerSecond = number of main iterations per second
# x$stats$success - boolean indicating if the simulation was successful
# x$stats$nbounces - total number of bounces
# x$stats$touchingVertex - vector with the indices of the vertices that land on the floor
# x$stats$faceDown - face value of the lowest triangle 

#########################################################################################
### Load libraries  #####################################################################
#########################################################################################

library(compiler)
library(inline)
library(rgl)

#########################################################################################
### Example #############################################################################
#########################################################################################

toss <- function(shape='cube',
                 stretch=c(1,2,3),              # stretch factors of the shape in x,y,z
                 k=0.5,                         # bouncing coefficient (between 0 and 1)
                 u=0.1,                         # friction coefficient (between 0 and 1)
                 g=c(0,0,-9.81),                # [m/s^2] gravitational acceleration
                 r=c(0,0,3),                    # [m] initial center of mass position
                 v=c(0,0,5),                    # [m/s] initial center of mass velocity
                 A=diag(3),                     # rotation matrix describing the initial orientation
                 L=c(2,3,3),                    # [kg m^2/s] # initial centre of mass angular momentum
                 text_out=T,                    # defines whether a text summary is given
                 animation_out=T,               # defines whether an animation is displayed
                 animation_frames=0,            # if 0, show movie, if >0, show a single frame
                 animation_pixels=600,          # size of the animation window in pixels
                 animation_path='',             # if path given, frames are saved as images
                 animation_newWindow=TRUE,      # set whether to open a new animation window
                 animation_showVectors=FALSE,   # set whether to show vectors
                 sim_dtOutput=0.02,             # [s] time step between discrete outputs; if 0, no intermediate timesteps are produced
                 sim_maxRotationAngle=0.01,     # [rad] maximal rotation angle; smaller values give higher accuracy
                 sim_tmax=1e99,                 # [s] maximal simulation time before abortion
                 sim_nIterationsMax=1e4,        # maximal number of iterations before abortion
                 sim_stopASAP=FALSE)            # if set to TRUE, the simulation is stopped as soon as the landing face is fixed 
{
  ICs = makeICs(shape=shape,stretch=stretch,k=k,u=u,g=g,r=r,v=v,A=A,L=L)
  x = simulation(ICs,
                 dtOutput=sim_dtOutput,
                 maxRotationAngle=sim_maxRotationAngle,
                 tmax=sim_tmax,
                 nIterationsMax=sim_nIterationsMax,
                 stopASAP=sim_stopASAP)
  if (text_out) summarizeSimulation(x)
  if (animation_out) showSimulation(x,
                                    frames=animation_frames,
                                    pixels=animation_pixels,
                                    path=animation_path,
                                    newWindow=animation_newWindow,
                                    showVectors=animation_showVectors)
}

endlessToss <- function() {
  set.seed(1)
  i = 0
  repeat{
    i = i+1
    ICs = randomICs(e=0.4)
    ICs$init$r[3] = ICs$init$r[3]*5
    x = simulation(ICs,0.02,nIterationsMax=2e4)
    showSimulation(x,newWindow=(i==1))
  }
}

#########################################################################################
### Output Functins #####################################################################
#########################################################################################

summarizeSimulation <- function(x) {
  cat('Computation time (s):  ',x$stats$computationTime,'\n')
  cat('Physical time (s):     ',x$stats$physicalTime,'\n')
  cat('Number of iterations:  ',x$stats$numberOfIterations,'\n')
  cat('Number of it per sec:  ',x$stats$iterationsPerSecond,'\n')
  cat('Face touching ground:  ',x$stats$faceDown,'\n')
  cat('Number of bounces:     ',x$stats$nbounces,'\n')
  cat('Simulation successful: ',x$stats$success,'\n')
  cat('=========================================\n')
}

showSimulation <- function(x,frames=0,pixels=600,path='',newWindow=TRUE,showVectors=FALSE) {
  
  lightVec = c(0,1,-1.5)
  numberOfTriangles = dim(x$triangle$vertex)[1]
  
  # define x,y,z limits of 3D plot
  if (frames==0) {frames=seq(length(x$dyn))}
  xmin = ymin = 1e99
  xmax = ymax = zmax = -1e99
  for (frame in frames) {
    xmin = min(xmin,x$dyn[[frame]]$c[1,])
    xmax = max(xmax,x$dyn[[frame]]$c[1,])
    ymin = min(ymin,x$dyn[[frame]]$c[2,])
    ymax = max(ymax,x$dyn[[frame]]$c[2,])
    zmax = max(zmax,x$dyn[[frame]]$c[3,])
  }
  h = sqrt(max(apply(x$init$c^2,2,sum)))
  d = max(xmax-xmin,ymax-ymin,zmax)+2*h
  mx = (xmax+xmin)/2
  my = (ymax+ymin)/2
  xmax = mx+d/2
  xmin = mx-d/2
  ymax = my+d/2
  ymin = my-d/2
  zmax = d
  s = max(xmax-xmin,ymax-ymin,zmax)
  zmin = -3e-3*s
  
  # delete old plots
  if (newWindow) {
    while (rgl.cur()[[1]]>0) {rgl.close()}
    open3d(windowRect = c(50,50,pixels,pixels))
  } else {
    clear3d()
  }
  
  # initiate plt
  rot = list()
  angle = c(0,-90,-130)/180*pi
  for (i in seq(3)) {
    j = i+1
    if (j==4) {j=1}
    ci = cos(angle[i])
    si = sin(angle[i])
    rot[[i]] = diag(3)
    rot[[i]][c(i,j),c(i,j)] = rot[[i]][c(i,j),c(i,j)]%*%matrix(c(ci,si,-si,ci),2,2)
  }
  rotation = diag(4)
  rotation[1:3,1:3] = rot[[1]]%*%rot[[3]]%*%rot[[2]]
  rgl.viewpoint(userMatrix=rotation,fov=60,zoom=0.88)
  rgl.clear(type="lights")
  rgl.light(viewpoint.rel=FALSE,ambient='black',diffuse='white',specular='black',x=mx,y=ymax+d*0.5,z=zmax*0.2)
  rgl.light(viewpoint.rel=FALSE,ambient='black',diffuse='white',specular='black',x=xmax+d*0.5,y=my,z=zmax*0.5)
  decorate3d(xlim=c(xmin,xmax),ylim=c(ymin,ymax),zlim=c(zmin,zmax),axes=FALSE,xlab='',ylab='',zlab='')
  
  # draw floor
  planes3d(0, 0, -1, -2e-3*s, alpha=1,col=rgb(0.9,0.9,0.9))
  
  # draw cuboid
  for (frame in frames) {
    
    # stop drawing (i.e. draw frame in memory and only update figure when done)
    par3d(skipRedraw=TRUE)
    
    # delete previous frame
    if (frame!=frames[1]) {
      for (i in seq(numberOfTriangles*2)) {rgl.pop()}
      if (showVectors) {
        for (i in seq(2)) {rgl.pop()}
      }
    }
    
    # draw poition and velocity vectors
    if (showVectors) {
      
      r = matrix(x$dyn[[frame]]$r,3,1)
      points3d(r[1],r[2],r[3],size=2,col='black')
      
      v = matrix(x$dyn[[frame]]$w,3,1)
      v = v/sqrt(sum(v^2))*h*2
      points3d(r[1]+c(0,v[1]),r[2]+c(0,v[2]),r[3]+c(0,v[3]),size=2,col='blue')
      lines3d(r[1]+c(0,v[1]),r[2]+c(0,v[2]),r[3]+c(0,v[3]),lwd=3)
      
      v = matrix(x$dyn[[frame]]$L,3,1)
      v = v/sqrt(sum(v^2))*h*2
      lines3d(r[1]+c(0,v[1]),r[2]+c(0,v[2]),r[3]+c(0,v[3]))
    }
    
    # draw triangles with shadows
    corner = x$dyn[[frame]]$c
    corner[3,] = pmax(0,corner[3,])
    for (i in seq(numberOfTriangles)) {
      
      # triangles
      c = col2rgb(x$triangle$color[i])/255
      triangle = x$triangle$vertex[i,]
      rgl.triangles(corner[1,triangle],corner[2,triangle],corner[3,triangle],col=rgb(c[1],c[2],c[3]),alpha=0.6)
      
      # shadow
      shadow = matrix(0,3,3)
      for (j in seq(3)) {
        shadow[1,j] = min(xmax,max(xmin,corner[1,triangle[j]]-corner[3,triangle[j]]*lightVec[1]/lightVec[3]))
        shadow[2,j] = min(ymax,max(ymin,corner[2,triangle[j]]-corner[3,triangle[j]]*lightVec[2]/lightVec[3]))
      }
      rgl.triangles(shadow[1,],shadow[2,],shadow[3,]-2e-4*s,col=rgb(0.5,0.5,0.5))
    }
    
    # update figure
    par3d(skipRedraw=FALSE)
    
    # save image to file
    if (nchar(path)>0) {
      file = sprintf(paste0(path,"frame_%04d.png"),frame)
      rgl.snapshot(filename=file,fmt="png",top=TRUE )
    }
    
  }
}

#########################################################################################
### Make initial conditions #############################################################
#########################################################################################

makeICs <- function(shape='cube',stretch=c(1,1,1),k=0.5,u=0.1,g=c(0,0,-9.81),r=c(0,0,3),v=c(1,0,0),A=diag(3),L=c(0,1,1),nPoints=50) {
  
  x = list()
  
  # shape independent initial variables
  x$init$k = k
  x$init$u = u
  x$init$g = g
  x$init$r = r
  x$init$v = v
  x$init$A = A
  x$init$L = L
  
  # definition of the shape
  if (shape=='cube') {
    
    # cube of sidelength 1 and mass 1
    x$init$c = matrix(c(1,1,1,1,1,-1,1,-1,1,1,-1,-1,-1,1,1,-1,1,-1,-1,-1,1,-1,-1,-1),3,8)/2
    x$init$I = diag(3)/6 # correct, verified numerically
    x$triangle$vertex = rbind(c(1,2,3),c(2,3,4),c(1,2,5),c(2,5,6),c(1,5,3),c(5,3,7),c(2,6,4),c(6,4,8),c(3,4,7),c(4,7,8),c(5,6,7),c(6,7,8))
    x$triangle$color = c('red','red','green','green','blue','blue','yellow','yellow','magenta','magenta','cyan','cyan')
    x$triangle$face = c(1,1,2,2,3,3,4,4,5,5,6,6)
    
  } else if (shape=='tetrahedron') {
    
    # tetrahedron of sidelength 1 and mass 1
    x$init$c = cbind(c(1,0,-1/sqrt(2)),c(-1,0,-1/sqrt(2)),c(0,1,1/sqrt(2)),c(0,-1,1/sqrt(2)))/2
    print(fct.vectNorm(x$init$c[,1]-x$init$c[,4]))
    x$init$I = diag(3)/20 # correct, verified numerically (wiki was wrong)
    x$triangle$vertex = rbind(c(1,2,3),c(1,2,4),c(1,3,4),c(2,3,4))
    x$triangle$color = c('red','green','blue','yellow')
    x$triangle$face = seq(4)
    
  } else if (shape=='cone') {
    
    # cone of height 1, diameter 1, and mass 1
    nphi = max(3,nPoints-2)
    h = 1
    r = 0.5
    dphi = 2*pi/nphi
    phi = seq(0,2*pi-dphi/2,by=dphi)
    vertex = rbind(r*cos(phi),r*sin(phi),rep(-h/4,nphi))
    vertex = cbind(vertex,c(0,0,-h/4),c(0,0,3/4*h))
    triangle = cbind(rep(seq(nphi),2),rep(seq(nphi)%%nphi+1,2),c(rep(nphi+1,nphi),rep(nphi+2,nphi)))
    x$init$c = vertex
    x$init$I = diag(3)*3/40 # correct, verified numerically
    x$triangle$vertex = triangle
    x$triangle$color = c(rep('red',nphi),rep('yellow',nphi/2),rep('blue',nphi/2))
    x$triangle$face = seq(2*nphi)
    
  } else {
    
    stop('Shape name unknown.')
    
  }
  
  # stretch shape
  x$init$c = stretch*x$init$c
  S = diag(c(sum(stretch[2:3]^2),sum(stretch[c(1,3)]^2),sum(stretch[1:2]^2)))/2
  x$init$I = S%*%x$init$I
  
  return(x)
  
}

randomICs <- function(e=1,k=0.5,u=0.7) {
  
  # fixed variables
  g = c(0,0,-9.81)
  s = (c(runif(3))*0.9+0.1)*2
  
  # dynamical variables
  if (e>=1) { # equipartition
    Epot = e*fct.vectNorm(s)*fct.vectNorm(g)/2
    Etra = Epot*3
    Erot = Epot*3
  } else {
    if (e>=1/7) {
      Epot = fct.vectNorm(s)*fct.vectNorm(g)/2
      Etot = 7*e*fct.vectNorm(s)*fct.vectNorm(g)/2
      Etra = (Etot-Epot)/2
      Erot = (Etot-Epot)/2
    } else {
      stop('normalized energy must be larger than 1/7')
    }
  }
  r = c(0,0,Epot/fct.vectNorm(g))
  v = fct.randomUnitVector()*sqrt(2*Etra)
  A = fct.randomRotationMatrix()
  I = A%*%diag(c(sum(s[2:3]^2),sum(s[c(1,3)]^2),sum(s[1:2]^2))/12)%*%t(A)
  Ltmp = fct.randomUnitVector()
  wtmp = solve(I)%*%Ltmp
  Erottmp = t(wtmp)%*%Ltmp/2
  L = Ltmp*sqrt(Erot/Erottmp)
  
  # finalize
  return(makeICs(shape='cube',k=k,u=u,g=g,s=s,r=r,v=r,A=A,L=L))
  
}

#########################################################################################
### Vector & Matrix Operations ##########################################################
#########################################################################################

vectProd_slow <- function(x,y){
  # cross-product of two 3-vectors
  return(c(x[2]*y[3]-x[3]*y[2],x[3]*y[1]-x[1]*y[3],x[1]*y[2]-x[2]*y[1]))
}
fct.vectProd <- cmpfun(vectProd_slow)

vectNorm_slow <- function(x){
  # norm of an n-vector
  return(sqrt(sum(x*x)))
}
fct.vectNorm <- cmpfun(vectNorm_slow)

rotationMatrix_slow <- function(x){
  # 3D-rotation matrix for a rotation around the axis parallel to x and an angle=|x|
  theta = fct.vectNorm(x)
  x = x/theta
  xsi = x*sin(theta)
  ci = cos(theta)
  xcim = x*(1-ci)
  return(matrix(c(ci+x[1]*xcim[1],x[1]*xcim[2]+xsi[3],x[1]*xcim[3]-xsi[2],
                  x[1]*xcim[2]-xsi[3],ci+x[2]*xcim[2],x[2]*xcim[3]+xsi[1],
                  x[1]*xcim[3]+xsi[2],x[2]*xcim[3]-xsi[1],ci+x[3]*xcim[3]),3))
}
fct.rotationMatrix <- cmpfun(rotationMatrix_slow)

fct.randomUnitVector <- function(){
  # draw random 3D unit vector from an isotropic distribution
  repeat {
    x = c(runif(3)-0.5)
    if (fct.vectNorm(x)<0.5) {break}
  }
  return(x/fct.vectNorm(x))
}

fct.randomRotationMatrix <- function(){
  # draw random 3D rotation matrix from and isotropic distribution
  x1 = fct.randomUnitVector()
  x2 = fct.randomUnitVector()
  x3 = fct.vectProd(x1,x2)
  x1 = fct.vectProd(x2,x3)
  return(cbind(x1/fct.vectNorm(x1),x2/fct.vectNorm(x2),x3/fct.vectNorm(x3)))
}

antisymmetricTensor_slow <- function(x){
  return(matrix(c(0,x[3],-x[2],-x[3],0,x[1],x[2],-x[1],0),3))
}
fct.antisymmetricTensor <- cmpfun(antisymmetricTensor_slow)

#########################################################################################
### Free Fall Motion ###################3################################################
#########################################################################################

{code <- "
 // translation
 for (int i=0; i < 3; i++) {
 r[i] = r[i]+v[i]*(*dt)+g[i]*(*dt)*(*dt)*0.5;
 v[i] = v[i]+g[i]*(*dt);
 }
 
 // rotation
 if (L[0]*L[0]+L[1]*L[1]+L[2]*L[2]>0) {
 double tsub = 0;
 while(tsub<*dt) {
 double d = *accuracy/sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
 if (d+tsub>*dt) {d=*dt-tsub;}
 tsub = tsub+d;
 double dw[] = {d*w[0],d*w[1],d*w[2]};
 double theta = sqrt(dw[0]*dw[0]+dw[1]*dw[1]+dw[2]*dw[2]);
 double x[] = {dw[0]/theta,dw[1]/theta,dw[2]/theta};
 double s = sin(theta);
 double c = cos(theta);
 double xs[] = {x[0]*s,x[1]*s,x[2]*s};
 double xcm[] = {x[0]*(1-c),x[1]*(1-c),x[2]*(1-c)};
 double rot[] = {c+x[0]*xcm[0],x[0]*xcm[1]-xs[2],x[0]*xcm[2]+xs[1],
 x[0]*xcm[1]+xs[2],c+x[1]*xcm[1],x[1]*xcm[2]-xs[0],
 x[0]*xcm[2]-xs[1],x[1]*xcm[2]+xs[0],c+x[2]*xcm[2]};
 double Atmp[9];
 for (int i=0;i<3;i++) {
 for (int j=0;j<3;j++) {
 Atmp[i*3+j] = rot[i*3]*A[j*3]+rot[i*3+1]*A[j*3+1]+rot[i*3+2]*A[j*3+2];
 }
 }
 double M1[9];
 for (int i=0;i<3;i++) {
 for (int j=0;j<3;j++) {
 M1[i*3+j] = Atmp[i*3]*Iinv_[j]+Atmp[i*3+1]*Iinv_[j+3]+Atmp[i*3+2]*Iinv_[j+6];
 }
 }
 double M[9];
 for (int i=0;i<3;i++) {
 for (int j=0;j<3;j++) {
 M[i*3+j] = M1[i*3]*Atmp[j*3]+M1[i*3+1]*Atmp[j*3+1]+M1[i*3+2]*Atmp[j*3+2];
 }
 }
 double wtmp[3];
 wtmp[0] = M[0]*L[0]+M[1]*L[1]+M[2]*L[2];
 wtmp[1] = M[3]*L[0]+M[4]*L[1]+M[5]*L[2];
 wtmp[2] = M[6]*L[0]+M[7]*L[1]+M[8]*L[2];
 dw[0] = 0.5*d*(w[0]+wtmp[0]);
 dw[1] = 0.5*d*(w[1]+wtmp[1]);
 dw[2] = 0.5*d*(w[2]+wtmp[2]);
 theta = sqrt(dw[0]*dw[0]+dw[1]*dw[1]+dw[2]*dw[2]);
 x[0] = dw[0]/theta;
 x[1] = dw[1]/theta;
 x[2] = dw[2]/theta;
 s = sin(theta);
 c = cos(theta);
 xs[0] = x[0]*s;
 xs[1] = x[1]*s;
 xs[2] = x[2]*s;
 xcm[0] = x[0]*(1-c);
 xcm[1] = x[1]*(1-c);
 xcm[2] = x[2]*(1-c);
 rot[0] = c+x[0]*xcm[0];
 rot[1] = x[0]*xcm[1]-xs[2];
 rot[2] = x[0]*xcm[2]+xs[1];
 rot[3] = x[0]*xcm[1]+xs[2];
 rot[4] = c+x[1]*xcm[1];
 rot[5] = x[1]*xcm[2]-xs[0];
 rot[6] = x[0]*xcm[2]-xs[1];
 rot[7] = x[1]*xcm[2]+xs[0];
 rot[8] = c+x[2]*xcm[2];
 double Anew[9];
 for (int i=0;i<3;i++) {
 for (int j=0;j<3;j++) {
 Anew[j*3+i] = rot[i*3]*A[j*3]+rot[i*3+1]*A[j*3+1]+rot[i*3+2]*A[j*3+2];
 }
 }
 for (int i=0;i<9;i++) {
 A[i] = Anew[i];
 }
 for (int i=0;i<3;i++) {
 for (int j=0;j<3;j++) {
 M[i*3+j] = A[i]*Iinv_[j*3]+A[i+3]*Iinv_[j*3+1]+A[i+6]*Iinv_[j*3+2];
 }
 }
 double Iinv[9];
 for (int i=0;i<3;i++) {
 for (int j=0;j<3;j++) {
 Iinv[j*3+i] = M[i*3]*A[j]+M[i*3+1]*A[j+3]+M[i*3+2]*A[j+6];
 }
 }
 w[0] = Iinv[0]*L[0]+Iinv[3]*L[1]+Iinv[6]*L[2];
 w[1] = Iinv[1]*L[0]+Iinv[4]*L[1]+Iinv[7]*L[2];
 w[2] = Iinv[2]*L[0]+Iinv[5]*L[1]+Iinv[8]*L[2];
 }
 }
 *t = *t+(*dt);
 "}
fct.freeFall.c = cfunction(signature(r='n',v='n',L='n',w='n',A='n',Iinv_='n',g='n',t='n',dt='n',accuracy='n'),code,convention=".C")
freeFall_slow <- function(dyn,g,Iinv,dt,accuracy) {
  
  # accuracy is in radians
  
  CcodeAcceleration = T
  
  if (CcodeAcceleration) {
    
    ff = fct.freeFall.c(dyn$r,dyn$v,dyn$L,dyn$w,dyn$A,Iinv,g,dyn$t,dt,accuracy)
    dyn$r = ff$r
    dyn$v = ff$v
    dyn$A = matrix(ff$A,3)
    dyn$w = ff$w
    dyn$t = ff$t
    
  } else {
    
    # translation
    dyn$r = dyn$r+dyn$v*dt+g*dt*dt/2
    dyn$v = dyn$v+g*dt
    
    # rotation
    if (fct.vectNorm(dyn$L)>0) {
      tsub = 0
      while(tsub<dt) {
        dtsub = min(accuracy/fct.vectNorm(dyn$w),dt-tsub)
        tsub = tsub+dtsub
        Atmp = fct.rotationMatrix(dyn$w*dtsub)%*%dyn$A
        wtmp = Atmp%*%Iinv%*%t(Atmp)%*%dyn$L
        dyn$A = fct.rotationMatrix((dyn$w+wtmp)/2*dtsub)%*%dyn$A
        dyn$w = dyn$A%*%Iinv%*%t(dyn$A)%*%dyn$L
      }
    }
    
    # time step
    dyn$t = dyn$t+dt
    
  }
  
  # time step
  return(dyn)
}
freeFall <- cmpfun(freeFall_slow)

#########################################################################################
### Main simulation #####################################################################
#########################################################################################

simulation_slow <- function(x,dtOutput=0,maxRotationAngle=0.01,tmax=1e99,nIterationsMax=1e6,stopASAP=FALSE) {
  
  # complete initial conditions and initialize variables
  Iinv = solve(x$init$I)
  dyn = list(r=x$init$r,v=x$init$v,A=x$init$A,L=x$init$L)
  dyn$w = dyn$A%*%Iinv%*%t(dyn$A)%*%dyn$L
  dyn$n = 0
  dyn$t = 0
  nc = dim(x$init$c)[2]
  cavg = mean(sqrt(colSums(x$init$c^2)))
  eps = 1e-13*cavg
  threshold = 1e-3*cavg
  radius = max(sqrt(colSums(x$init$c^2)))
  last_down = rep(FALSE,dim(x$init$c)[2])
  nIterations = 0
  timeResolution = 0.1
  bouncingCorners = matrix(FALSE,nc)
  
  # save initial conditions
  if (dtOutput>0) {
    x$dyn[[1]] = dyn
    tOutput=dtOutput
  } else {tOutput=1e99}
  
  # time-integration
  tStart = Sys.time()
  while(nIterations<nIterationsMax && dyn$t<tmax) {
    
    # compute adaptive time step
    vmax = fct.vectNorm(dyn$v)+fct.vectNorm(dyn$w)*radius
    dt = min(timeResolution*cavg/vmax,timeResolution*sqrt(2*cavg/fct.vectNorm(x$init$g)),tmax-dyn$t,tOutput-dyn$t)
    nIterations = nIterations+1
    
    # integrate free-fall motion
    dynOld = dyn
    dyn = freeFall(dyn,x$init$g,Iinv,dt,maxRotationAngle)
    
    # handle contact
    if (dyn$r[3]<radius) {
      pc3 = dyn$A[3,]%*%x$init$c+dyn$r[3]
      
      if (any(pc3<=0)) {
        
        # move back to precise contact time
        gc = dyn$A%*%x$init$c
        vc3 = dyn$v[3]+dyn$w[1]*gc[2,]-dyn$w[2]*gc[1,]
        touch = (gc[3,]<=+eps-dyn$r[3]) & (vc3<0)
        if (!any(touch & bouncingCorners)) {
          vmax = fct.vectNorm(dyn$v)+fct.vectNorm(dyn$w)*radius
          pc3old = dynOld$A[3,]%*%x$init$c+dynOld$r[3]
          for (i in seq(nc)) {
            gc = dyn$A%*%x$init$c[,i]
            vc3 = dyn$v[3]+dyn$w[1]*gc[2]-dyn$w[2]*gc[1]
            touch = (gc[3]<=-eps-dyn$r[3]) & (vc3<0) & (pc3old[i]>=0)
            if (touch) {
              pc3 = gc[3]+dyn$r[3]
              repeat {
                pc3before = pc3
                factor = min(2,pc3old[i]/(abs(pc3old[i]-pc3)+eps))
                dtold = dt
                dt = dt*factor
                dyn = freeFall(dynOld,x$init$g,Iinv,dt,maxRotationAngle)
                pc3 = dyn$A[3,]%*%x$init$c[,i]+dyn$r[3]
                if (abs(pc3)+2*eps>=abs(pc3before)) {
                  dyn = freeFall(dynOld,x$init$g,Iinv,dtold,maxRotationAngle)
                  break
                }
                if (abs(1-factor)<1e-8) {break}
              }
            }
          }
        } else {
          timeResolution = 0.01
        }
        bouncingCorners = matrix(FALSE,nc)
        
        # check end
        pc3 = dyn$A[3,]%*%x$init$c+dyn$r[3]
        down = pc3<threshold
        if (length(which(down&last_down))>=3) {
          Etra = t(dyn$v)%*%dyn$v/2
          Epot = -t(x$init$g)%*%dyn$r
          Erot = t(dyn$w[1:2])%*%dyn$L[1:2]/2
          if (Etra+Erot<1e-2*Epot | stopASAP) {break}
        }
        last_down = down
        
        # compute collision
        gc = dyn$A%*%x$init$c
        vc3 = dyn$v[3]+dyn$w[1]*gc[2,]-dyn$w[2]*gc[1,]
        touch = (gc[3,]<=+eps-dyn$r[3]) & (vc3<0)
        if (any(touch)) {
          id = which(touch)
          gcMean = rowMeans(matrix(gc[,id],3))
          vc = dyn$v+fct.vectProd(dyn$w,gcMean)
          dvc = c(-vc[1],-vc[2],-{1+x$init$k}*vc[3])
          C = fct.antisymmetricTensor(gcMean)
          H = diag(3)-C%*%dyn$A%*%Iinv%*%t(dyn$A)%*%C
          Q = solve(H)%*%dvc
          Qh = sqrt(Q[1]^2+Q[2]^2)
          if (Qh>x$init$u*Q[3]) {
            q = c(x$init$u*Q[1]/Qh,x$init$u*Q[2]/Qh,1)
            Q3 = -{1+x$init$k}*vc[3]/sum(H[3,]*q)
            Q = q*Q3
          }
          dyn$v = dyn$v+Q
          dyn$L = dyn$L+fct.vectProd(gcMean,Q)
          dyn$w = dyn$A%*%Iinv%*%t(dyn$A)%*%dyn$L
          dyn$n = dyn$n+1
          bouncingCorners = touch
        } else {bouncingCorners = matrix(FALSE,nc)}
      } else {bouncingCorners = matrix(FALSE,nc)}
    } else {bouncingCorners = matrix(FALSE,nc)}
    
    # save data
    if ((dtOutput>0) && (dyn$t>=tOutput)) {
      x$dyn[[length(x$dyn)+1]] = dyn
      tOutput = tOutput+dtOutput
    }
    
  }
  tEnd = Sys.time()
  
  # finalize output
  if (dtOutput>0) {
    x$dyn[[length(x$dyn)+1]] = dyn
    nout = length(x$dyn)
    for (i in seq(nout)) {
      x$dyn[[i]]$c = x$dyn[[i]]$A%*%x$init$c+matrix(x$dyn[[i]]$r,3,nc)
    }
  }
  x$stats$numberOfIterations = nIterations
  x$stats$computationTime = as.double(tEnd)-as.double(tStart)
  x$stats$physicalTime = dyn$t
  x$stats$iterationsPerSecond = round(nIterations/x$stats$computationTime)
  x$stats$success=nIterations<nIterationsMax
  x$stats$nbounces=dyn$n
  pc3 = dyn$A[3,]%*%x$init$c+dyn$r[3]
  x$stats$touchingVertex = which(pc3<threshold)
  
  # find triangle facing down
  if (!is.null(x$triangle)) {
    zmean = (pc3[x$triangle$vertex[,1]]+pc3[x$triangle$vertex[,2]]+pc3[x$triangle$vertex[,3]])/3
    x$stats$faceDown = x$triangle$face[which.min(zmean)]
  } else {
    x$stats$faceDown = NA
  }

  return(x)
  
}
simulation <- cmpfun(simulation_slow)