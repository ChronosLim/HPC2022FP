static char help[] = "Solve IEFDMS \n"

#include <petscmat.h>
#include <stdio.h>
#include <math.h>
  
#define PI acos(-1)

int main(int argc,char **args)
{
  Vec            zk0,zk,yk0,yk,Az;
  Mat            A;
  PetscReal      err,tol=1e-8,dt,dx,alpha,rho=1.0,c=1.0,kappa=1.0,l=2.0;
  PetscErrorCode ierr;
  PetscInt       i,n = 100,col[3],its,rstart,rend,nlocal;
  PetscScalar    zero = 0.0, one = 1.0, value[3],uold,u,f;
  
  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-dt",&dt,NULL);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"\nmatrix size n=%d \n",n);
  
  dt = 1e-6;
  dx = 1/n;
  alpha = kappa * dt / (rho * c * dx *dx);
  
  ierr = VecCreate(PETSC_COMM_WORLD,&zk0);CHKERRQ(ierr);
  ierr = VecSetSizes(zk0,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(zk0);CHKERRQ(ierr);
  ierr = VecDuplicate(zk0,&zk);CHKERRQ(ierr);
  ierr = VecDuplicate(zk0,&yk0);CHKERRQ(ierr);
  ierr = VecDuplicate(zk0,&yk);CHKERRQ(ierr);
  ierr = VecDuplicate(zk0,&Az);CHKERRQ(ierr);

  ierr = VecGetOwnershipRange(zk0,&rstart,&rend);CHKERRQ(ierr);
  ierr = VecGetLocalSize(zk0,&nlocal);CHKERRQ(ierr);



  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,nlocal,nlocal,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);



  if (!rstart) 
  {
    rstart = 1;
    i      = 0; col[0] = 0; col[1] = 1; value[0] = 1.0 + 2.0 * alpha; value[1] = -alpha;
    ierr   = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);

  }  

  if (rend == n) 
  {
    rend = n-1;
    i    = n-1; col[0] = n-2; col[1] = n-1; value[0] = -alpha; value[1] = 1.0 + 2.0 * alpha;
    ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }

  value[0] = -alpha; value[1] = 1.0 + 2.0 * alpha; value[2] = -alpha;
  for (i=rstart; i<rend; i++) 
  {
    col[0] = i-1; col[1] = i; col[2] = i+1;
    ierr   = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  //ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);



  ierr = VecSet(zk0,zero);CHKERRQ(ierr);
  ierr = VecSetValue(zk0,0,one,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(zk0);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(zk0);CHKERRQ(ierr);
  //ierr = VecView(z,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  

  ierr = MatMult(A,zk0,yk);CHKERRQ(ierr);
  //ierr = VecView(yk0,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = VecNorm(yk,NORM_2,&u);CHKERRQ(ierr); 
  ierr = VecScale(yk,1.0/u);CHKERRQ(ierr); 
  //ierr = VecView(yk0,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = VecAYPX(zk,0,yk);CHKERRQ(ierr);
  ierr = MatMult(A,zk,yk0);CHKERRQ(ierr);
  ierr = VecNorm(yk0,NORM_2,&uold);CHKERRQ(ierr); 
  err=PetscAbsReal((uold-u)/u);
  //printf("error=%f \n",err);

  its = 0;

  while(err>tol)
  {     
    ierr = VecScale(yk0,1.0/uold);CHKERRQ(ierr); 
    ierr = VecAYPX(zk0,0,yk0);CHKERRQ(ierr);
    ierr = MatMult(A,zk0,yk);CHKERRQ(ierr);
    ierr = VecNorm(yk,NORM_2,&u);CHKERRQ(ierr); 
    ierr=PetscAbsReal((u-uold)/uold);
    ierr = VecAYPX(yk0,0,yk);CHKERRQ(ierr);
    uold = u;
    its += 1;   
  }

  

  PetscPrintf(PETSC_COMM_WORLD,"its=%d, error=%.16g \n",its,err);
  //ierr = VecView(z,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  


  ierr = MatMult(A,zk0,Az);CHKERRQ(ierr);
  //ierr = VecView(Az,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = VecTDot(zk0,Az,&f);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"temperature=%.6g \n",u);



  ierr = VecDestroy(&zk0);CHKERRQ(ierr);  
  ierr = VecDestroy(&zk);CHKERRQ(ierr);
  ierr = VecDestroy(&Az);CHKERRQ(ierr);  
  ierr = VecDestroy(&yk0);CHKERRQ(ierr);
  ierr = VecDestroy(&yk);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  //ierr = KSPDestroyk0(&ksp);CHKERRQ(ierr);



  ierr = PetscFinalize();
  return ierr;
}
