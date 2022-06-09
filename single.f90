Program main
    implicit none

!Select one preseted grid densities or set your onw one as your like to do experiment.
  
!  integer,parameter ::N=32
!  integer,parameter ::N=64
!  integer,parameter ::N=128
  integer,parameter ::N=100

!Define what parameter will be used.

  real*8 :: pi=4.0d0*atan(1.0d0),Length=1.0d0,rho=1.0d0,c=1.0d0,kappa=1.0d0,l=2.0d0
  real*8 :: dx,dt,tend,alpha,sum,sb,umax,e,tol,fo
  real*8,dimension(N+1) :: uold,u,x,uexact
  real*8,dimension(N+1) :: h,f,error
  real*8,dimension(N+1) :: ahe,beh,con,dia
  integer :: j,it,tex_start,tex_end,tim_start,tim_end,content
  integer,parameter :: nsteps=5000000
  real*8,dimension(nsteps)::k,rate
  real*8 ::C_1,C_2
  character(len=32) ::arg

!Initializing of parameter.

  e=0.0d0
  error=0.0d0
  sum=0.0d0
  k=0.0d0
  rate=0.0d0
  h=0.0d0
  f=0.0d0
  dx=Length/N
  dt=0.00004d0
  tol=1.0e-16
  sb=0.0d0
  fo=0.0d0
  
!Create files to store datas.

  open(53,file='uexact.dat',form='formatted',status='unknown')
  open(54,file='x.dat',form='formatted',status='unknown')
  open(55,file='u_EEFDMS.dat',form='formatted',status='unknown')
  open(56,file='h_EEFMDS.dat',form='formatted',status='unknown')
  open(57,file='k_EEFMDS.dat',form='formatted',status='unknown')
  open(58,file='rate_EEFMDS.dat',form='formatted',status='unknown')
  open(59,file='u_IEFMDS.dat',form='formatted',status='unknown')
  open(60,file='h_IEFMDS.dat',form='formatted',status='unknown')
  open(61,file='k_IEFMDS.dat',form='formatted',status='unknown')
  open(62,file='rate_IEFMDS.dat',form='formatted',status='unknown')
  
!Reset the initial value of N and dt from command line.
!The first argument should be an integer for N, and the second argument should be a &
!double precision number for dt.
  
  call get_command_argument(1,arg)
    if(len_trim(arg)==0) goto 996
      read(arg,'(F5.0)') content
      dt=content
      
  996 write(*,*)"Without out-source input"
  continue

!  alpha=0.5d0
!  dt=alpha*rho*c*dx**2/kappa
  alpha=dt*kappa/(rho*c*dx**2)
  write(*,*)"cfl number is",alpha
  write(*,*)"dt for EEFDMS is",dt
  
!Assertion for explicit Euler scheme.  
  
  if(alpha .gt. 0.5 .or. alpha .lt. 0) then
    write(*,*)"Stability criterion is out of range."
    write(*,*)"Program EEFDMS stopped."
    goto 999
  end if
  
  u=0.0d0
  uold=0.0d0
  tex_start=0
  tex_end=0
  


!Calculate the excact solution under I.C.:f=sin(pi*l*x),u0=e^x, u(0,t)=u(1,t)=0
!Consider that when t→∞, the temperature distrubition of our 1D domain should be stable, so the spatial integral of heat-source term should be 0.
!Thus, parameter l should be an odd number.
!In factor, due to the B.C.,this system is not a insulation system, however, the &
!source term f is not a function of t, so always can the temperature distribution be &
!stable even when the spatial integral of heat-source term is not 0.
!The analytical solution should be T(x,t)=rho*c*sin(pi*l*x)/(kappa*(pi*l)**2)+C1x+C2 when t→∞

  do j=1,N+1
    x(j)=(j-1.0d0)*dx
  end do
  
  do j=2,N
    f(j)=dsin(pi*l*x(j))
  end do
  
  C_1=-rho*c*dsin(pi*l*Length)/(kappa*(pi*l)**2)
  C_2=0.0d0
  
  do j=2,N
    uexact(j)=rho*c*f(j)/(kappa*(pi*l)**2)+C_1*x(j)+C_2
    
!Check the value of uexact.    

!    write(*,*)uexact(j)
    
  end do
  
  do j=2,N
    sum=sum+uexact(j)
  end do
  
  sum=sum/(N-1)
  write(*,*)"The exact average temperature is",sum
  write(53,*)uexact


!Go to explicit Euler FDM scheme to solve our PDE
!The numerical scheme is written as: rho*c/dt*(u^n+1_j-u^n_j)-kappa-dx**2 &
!*(u^n_j+1-2*u^n_j+u^n_j-1)=f(x_j)
  
  call system_clock(tex_start)
  
!Modified the temperature of vitual points to meet the B.C.     
    
  uold(1)=0
  uold(N+1)=0
  
  do j=2,N
    uold(j)=dexp(x(j))
  end do
  
  do j=2,N
    sum=sum+uold(j)
  end do
  
  sum=sum/(N-1)
  write(*,*)"Initial average temperature is",sum
  
!To check whether the I.C. has been deployed successfully.

!  do j=2,N
!    write(*,*)uold(j)
!  end do

!The main iterative process under explicit Euler FDM scheme. 
  
  do it=1,nsteps
    
    sum=0.0d0
    sb=0.0d0
    
!Update the temperature of vitual points to meet the B.C.
      
    u(1)=0
    u(N+1)=0
    
    do j=2,N
      u(j)=uold(j)+dt*(f(j)+kappa*(uold(j+1)-2*uold(j)+uold(j-1))/dx**2)/(rho*c)
    
      k(it)=k(it)+u(j)
      h(j)=kappa*(u(j+1)-u(j-1))/(2*dx)
      rate(it)=rate(it)+h(j)
      umax=maxval(u)
    end do
    
!    if(mod(it,1000).eq.0.) then
!    write(*,*)umax
!    endif
    
    k(it)=k(it)/N
    rate(it)=rate(it)/N
    
    do j=2,N
      sum=sum+dabs(u(j)-uold(j))
    end do
    sum=sum/(N-1)
    
    do j=1,N+1
      uold=u
    end do
    
    do j=2,N
      error(j)=dabs(u(j)-uexact(j))
    end do
  
    e=maxval(error)
    tend=dt*it
    sum=sum/(N-1)
    fo=dabs(sum-sb)
    
    if(mod(it,10000).eq.0.) then
!    write(*,*)fo
    endif
    
    if(fo .lt. tol) write(*,*)"Calculation reaches stable."
    if(fo .lt. tol) exit
    
    sb=sum
  
  end do

  
  
  write(*,*)"Programm stop at t =",tend
  
  write(*,*)"The error of explicit Euler FDM scheme is",e
  
  call system_clock(tex_end)
  
  write(*,*)"The time cost of singlethreading program to solve EEFDMS:",tex_end&
  -tex_start,"ms"
  
  do j=2,N
    sum=sum+u(j)
  end do
  write(*,*)"The latest average temperature is",sum
  


!  do j=2,N
!    write(*,*)u(j)
!  end do  

  write(54,*)x
  write(55,*)u
  write(56,*)h
  write(57,*)k
  write(58,*)rate



!Go to implicit Euler FDM scheme to solve our PDE
!The numerical scheme is written as: rho*c/dt*(u^n+1_j-u^n_j)-kappa &
!/dx**2*(u^n+1_j+1-2*u^n+1_j+u^n+1_j-1)=f(x_j)
!Initializing of parameter.

  e=0.0d0
  error=0.0d0
  sum=0.0d0
  k=0.0d0
  rate=0.0d0
  h=0.0d0
  f=0.0d0
  sb=0.0d0
  fo=0.0d0
  
  999 continue
  
!  alpha=0.5d0
!  dt=alpha*rho*c*dx**2/kappa
  alpha=dt*kappa/(rho*c*dx**2)
  write(*,*)"dt for IEFDMS is",dt
  
!Assertion for implicit Euler scheme.  
  
  if(alpha .lt. 0) then
    write(*,*)"Stability criterion is out of range."
    write(*,*)"Program IEFDMS stopped."
    stop 2333
  end if

  uold=0.0d0
  u=0.0d0
  tim_start=0
  tim_end=0
  
  call system_clock(tim_start)
  
  do j=2,N
    uold(j)=dexp(x(j))
    
!Modified the temperature of vitual points to meet the B.C.     
    
    uold(1)=0
    uold(N+1)=0
    
    
    
  end do
  
  do j=2,N
    sum=sum+uold(j)
  end do
  
  sum=sum/(N-1)
  write(*,*)"Initial average temperature is",sum
  
!To check whether the I.C. has been deployed successfully.

!  do j=2,N
!    write(*,*)uold(j)
!  end do

!The main iterative process under implicit Euler FDM scheme. 
  
  ahe=0.0d0
  beh=0.0d0
  dia=0.0d0
  ahe(3:N)=-alpha
  beh(2:N-1)=-alpha
  dia(2:N)=1+2*alpha
  
  do it=1,nsteps
    
    sum=0.0d0
    sb=0.0d0
    
    con(1)=0.0d0
    con(N+1)=0.0d0
    
    do j=2,N
      con(1)=0.0d0
      con(N+1)=0.0d0
      con(j)=uold(j)+f(j)
    end do
    
    call sy(2,N,ahe,dia,beh,con)
    
!Update the temperature of vitual points to meet the B.C.

    u(1)=0
    u(N+1)=0

    do j=2,N
      u(j)=con(j)
      
      k(it)=k(it)+u(j)
      h(j)=kappa*(u(j+1)-u(j-1))/(2*dx)
      rate(it)=rate(it)+h(j)
      umax=maxval(u)
    end do  
    
!    if(mod(it,1000).eq.0.) then
!    write(*,*)umax
!    endif
    
    k(it)=k(it)/N
    rate(it)=rate(it)/N
    
    do j=2,N
      sum=sum+dabs(u(j)-uold(j))
    end do
    sum=sum/(N-1)
    
    do j=1,N+1
      uold=u
    end do
    
    do j=2,N
      error(j)=dabs(u(j)-uexact(j))
    end do
  
    e=maxval(error)
    tend=dt*it
    sum=sum/(N-1)
    fo=dabs(sum-sb)
    
    if(mod(it,10000).eq.0.) then
!    write(*,*)fo
    endif
    
    if(fo .lt. tol) write(*,*)"Calculation reaches stable."
    if(fo .lt. tol) exit
    
    sb=sum
  
  end do

  
  
  write(*,*)"Programm stop at t =",tend
  
  write(*,*)"The error of implicit Euler FDM scheme is",e
  
  call system_clock(tim_end)
  
  write(*,*)"The time cost of singlethreading program to solve IEFDMS:",tim_end&
  -tim_start,"ms"
  
  do j=2,N
    sum=sum+u(j)
  end do
  write(*,*)"The latest average temperature is",sum
  
!  do j=2,N
!    write(*,*)u(j)
!  end do  

  write(59,*)u
  write(60,*)h
  write(61,*)k
  write(62,*)rate
  
  
  
end Program main



!******************************
! Subroutine SY, for solving a tridiagonal system of equations.  Based
! on the Thomas algorithm, from "Computational Fluid Mechanics and
! Heat Transfer," by Anderson, Tannehill, and Pletcher.
! 27 May 1992   Jim DeSpirito
!
      subroutine sy(il,iu,bb,dd,aa,cc)
!
!      double precision aa(iu),bb(iu),cc(iu),dd(iu)
      real*8 aa(iu),bb(iu),cc(iu),dd(iu)
!
!.....il = subscript of first equation
!.....iu = subscript of last equation
!.....bb = coefficient behind diagonal
!.....dd = coefficient on diagonal
!.....aa = coefficient ahead of diagonal
!.....cc = element of constant vector
!
!.....establish upper triangular matrix
!
      lp = il + 1
      do 10 i = lp,iu
       r = bb(i)/dd(i-1)
       dd(i) = dd(i) - r * aa(i-1)
       cc(i) = cc(i) - r * cc(i-1)
   10 continue
!
!.....back substitution
!
      cc(iu) = cc(iu) / dd(iu)
      do 20 i = lp,iu
       j = iu - i + il
       cc(j) = (cc(j) - aa(j) * cc(j+1)) / dd(j)
   20 continue
!
!.....solution stored in cc
!
      return
      end
!
!********************************************************
