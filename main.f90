Programm main
  implicit none

!Define what parameter will be used.

  real*8 ::pi=4.0de*atan(1.0d0),L=1.0d0,u0=1.0d0dx
  real*8 ::dt,tend,alpha,sum,umax
  real*8,dimension(N+2) :: uold,u,x
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
  write(*,*)"dt is",dt
  u=0.0d0
  uold=0.0d0
  
  do j=1,N+2
    x(j)=(j-1.5d0)*dx
  end do
  
  do j=2,N+1
    sum=sum+uold(j)
  end do
  
  sum=sum/N
  write(*,*)"Initial average temperature is",sum
  
!To check whether the I.C. has been deployed successfully.

!  write(*,*)uold

!The main iterative process with FDM scheme1. 
  
  do it=1,nsteps
    
    do j=2,N+1
      
      u(1)=u(N+1)
      u(N+2)=u(2)
      
      k(it)=
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
  write(*,*)"Programm stop at t=",tend
  
  write(54,*)x
  write(55,*)u(j)
  write(56,*)k
  write(57,*)vis

end Programm main
