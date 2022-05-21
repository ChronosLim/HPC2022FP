Programm main
  implicit none

!Define what parameter will be used.

  real*8 ::pi=4.0de*atan(1.0d0),L=1.0d0,rho=1.0d0,c=1.0d0,kappa=1.0d0,alpha=0.5d0
  real*8 ::dx,dt,tend,alpha,sum,umax
  real*8,dimension(N+2) :: uold,u,x,uexact
  integer::j,it
  integer,parameter::nsteps=1000
  real*8,dimension(nsteps)::k,rate

!Select one preseted grid densities or set your onw one as your like to do experiment.
  
  integer,paramater ::N=32
!  integer,paramater ::N=64
!  integer,paramater ::N=128
!  integer,paramater ::N=256

!Initializing of parameter.

  sum=0.0d0
  k=0.0d0
  rate=0.0d0
  dx=L/N
  dt=0.001d0
!  dt=alpha*rho*c*dx**2/kappa
  write(*,*)"dt is",dt
  u=0.0d0
  uold=0.0d0
  
  do j=1,N+2
    x(j)=(j-1.5d0)*dx
    uold(j)=DEXP(x(j))
    
!Modified the temperature of vitual points to meet the B.C.     
    
    u(1)=-u(2)
    u(N+2)=-u(N+1)
    
  end do
  
  do j=2,N+1
    sum=sum+uold(j)
  end do
  
  sum=sum/N
  write(*,*)"Initial average temperature is",sum
  
!To check whether the I.C. has been deployed successfully.

!  write(*,*)uold

!The main iterative process with explicit Euler FDM scheme. 
  
  do it=1,nsteps
    
    do j=2,N+1

!Update the temperature of vitual points to meet the B.C.
      
      u(1)=-u(2)
      u(N+2)=-u(N+1)
      
      k(it)=k(it)+u(j)
      rate(it)=
      umax=maxval(u)
    end do
    
    if(mod(it,1000).eq.0.) then
    write(*,*)umax
    endif
    
    k(it)=k(it)/N
    rate(it)=rate(it)/N
    
    do j=1,N+2
      uold=u
    end do
    
  end do
  
  tend=dt*nsteps
  write(*,*)"Programm stop at t =",tend
  
  write(54,*)x
  write(55,*)u(j)
  write(56,*)k
  write(57,*)vis

end Programm main
