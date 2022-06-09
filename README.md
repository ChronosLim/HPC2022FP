# HPC2022FP
This is the final project for HPC2022
In this package, two numerical scheme of finite differential method based algorithm are included to solve the transient heat equation.
The temporal differential scheme is always downwind scheme, the difference between those two included schemes is as follow:
Explicit Euler finite differential scheme:rho*c*(u^(n+1)_j-u^n_j)/dt-kappa*(u^n_(j+1)-2*u^n_j+u^n_(j-1))/(dx^2)=f(x_j)
Implicit Euler finite differential scheme:rho*c*(u^(n+1)_j-u^n_j)/dt-kappa*(u^(n+1)_(j+1)-2*u^(n+1)_j+u^(n+1)_(j-1))/(dx^2)=f(x_j)

The program was completed by fortran 90, and the compiler I used is gnu compiler.
The multithreading version is written for petcs which can run on a cluster.

Please make appropriate modifications to the CMakelist.txt or makefile to install the package.
If you want to run those codes on cluster, the jobscript file may be helpful.
