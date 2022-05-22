Programm main
  implicit none

!Define what parameter will be used.

  real*8 ::pi=4.0de*atan(1.0d0),L=1.0d0,rho=1.0d0,c=1.0d0,kappa=1.0d0,alpha=0.5d0
  real*8 ::dx,dt,tend,alpha,sum,umax,e
  real*8,dimension(N+2) :: uold,u,x,uexact,h,f,error
  integer::j,it
  integer,parameter::nsteps=1000
  real*8,dimension(nsteps)::k,rate

!Select one preseted grid densities or set your onw one as your like to do experiment.
  
  integer,paramater ::N=32
!  integer,paramater ::N=64
!  integer,paramater ::N=128
!  integer,paramater ::N=256

!Initializing of parameter.

  e=0.0d0
  error=0.0d0
  sum=0.0d0
  k=0.0d0
  rate=0.0d0
  h=0.0d0
  f=0.0f0
  dx=L/N
  dt=0.001d0
!  dt=alpha*rho*c*dx**2/kappa
  write(*,*)"dt is",dt
  u=0.0d0
  uold=0.0d0
  
!Calculate the excact solution under I.C.:f=sin(pi*l*x),u0=e^x, u(0,t)=u(1,t)=0
!Consider that when t→∞, the temperature distrubition of our 1D domain should be stable, so the spatial integral of heat-source term should be 0.
!Thus, parameter l should be an odd number.
!The analytical solution should be T(x,t)=-rho*c*sin(pi*l*x)/(kappa*(pi*l)**2)+C1x+C2 when t→∞
  
  real*8 ::l,C1,C2

  do j=1,N+2
    x(j)=(j-1.5d0)*dx
  end do
  
  do j=2,N+1
    f(j)=dsin(pi*l*x(j))
  end do
  
  l=2.0d0
  C1=rho*c*dsin(pi*l)/(kappa*(pi*l)**2)
  C2=0.0d0
  
  do j=2,N+1
    uexact(j)=-rho*c*f(j)/(kappa*(pi*l)**2)+C1*x(j)+C2
  end do


!Go to explicit Euler FDM scheme to solve our PDE
!The numerical scheme is written as: rho*c/dt*(u^n+1_j-u^n_j)-kappa/dx**2*(u^n_j+1-2*u^n_j+u^n_j-1)=f(x_j)
  
  do j=1,N+2
    uold(j)=dexp(x(j))
    
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

!The main iterative process under explicit Euler FDM scheme. 
  
  do it=1,nsteps
    
    do j=2,N+1
      u(j)=uold(j)+dt*(f(j)+kappa*(uold(j+1)-2*uold(j)+uold(j-1))/dx**2)/(rho*c)

!Update the temperature of vitual points to meet the B.C.
      
      u(1)=-u(2)
      u(N+2)=-u(N+1)
      
      k(it)=k(it)+u(j)
      h(j)=kappa*(u(j+1)-u(j-1))/(2*dx)
      rate(it)=rate(it)+h(j)
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
  
  do j=2,N+1
    error(j)=dabs(u(j)-uexact(j))
  end do
  
  e=maxval(error)
  write(*,*)"The error of explicit Euler FDM scheme is",e
  
  write(54,*)x
  write(55,*)u(j)
  write(56,*)h(j)
  write(57,*)k
  write(58,*)rate



!Go to implicit Euler FDM scheme to solve our PDE
!The numerical scheme is written as: rho*c/dt*(u^n+1_j-u^n_j)-kappa/dx**2*(u^n+1_j+1-2*u^n+1_j+u^n+1_j-1)=f(x_j)
!Initializing of parameter.

  e=0.0d0
  error=0.0d0
  sum=0.0d0
  k=0.0d0
  rate=0.0d0
  h=0.0d0
  f=0.0f0
  dx=L/N
  dt=0.001d0
!  dt=alpha*rho*c*dx**2/kappa
  write(*,*)"dt is",dt
  u=0.0d0
  uold=0.0d0
  
  do j=1,N+2
    uold(j)=dexp(x(j))
    
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

!The main iterative process under implicit Euler FDM scheme. 
  
  do it=1,nsteps
    
    do j=2,N+1
      u(j)=uold(j)+dt*(f(j)+kappa*(uold(j+1)-2*uold(j)+uold(j-1))/dx**2)/(rho*c)

!Update the temperature of vitual points to meet the B.C.
      
      u(1)=-u(2)
      u(N+2)=-u(N+1)
      
      k(it)=k(it)+u(j)
      h(j)=kappa*(u(j+1)-u(j-1))/(2*dx)
      rate(it)=rate(it)+h(j)
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
  
  do j=2,N+1
    error(j)=dabs(u(j)-uexact(j))
  end do
  
  e=maxval(error)
  write(*,*)"The error of explicit Euler FDM scheme is",e
  
  write(54,*)x
  write(55,*)u(j)
  write(56,*)h(j)
  write(57,*)k
  write(58,*)rate
  
end Programm main
