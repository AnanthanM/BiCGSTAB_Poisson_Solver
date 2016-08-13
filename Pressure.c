//Every function for solving Pressure Poisson 
//Equation is in this file

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>

#include "fvm.h"

void Get_RHS_of_Pressure_Poisson_Eqn(Domain domain,Field * RHS,Constant constant)
{
  Field * p   = domain.p;
   
  double dx = constant.dx;
  double dy = constant.dy;
  double dt = constant.dt;
  
  int N_cells_C = p->N;
  int Nx_C      = p->N_x;
  int Ny_C      = p->N_y;

  int l,i,j;

  for(l = 0;l<N_cells_C;l++)
  {
     RHS->val[l] = 0.0;
  }

  return;
}

void Pressure_Poisson(Field * p, Constant constant,double * Ap_vector,Domain domain)
{
  double dx = constant.dx;
  double dy = constant.dy;

  int N_cells_C    = p->N;
  int Nx_C         = p->N_x;
  int Ny_C         = p->N_y;

  int l,i,j;

  int      N,S,E,W;                       //for rho values
  
  double   pP,
           pN,pS,pE,pW;                  
  

  double   aP,
           aN,aS,aE,aW;                  

  int       row;
  
  double  * val_pP,* val_pE,* val_pW,* val_pN,* val_pS;

  for(l=0;l<N_cells_C;l++)
    Ap_vector[l] = 0.0;

  //BC Setup
  
  //West Side Neumann 
  
  i = 0;
  for(j=0;j<Ny_C;j++)
  {
    l = j*Nx_C + i;
    p->val[l] = p->val[l+1];
  }

  //East Side Neumann
  
  i = Nx_C - 1;
  for(j=0;j<Ny_C;j++)
  {
    l = j*Nx_C + i;
    p->val[l] = p->val[l-1];
  }
  
  //South side Dirchlet

  j = 0;
  for(i=0;i<Nx_C;i++)
  {
    l = j*Nx_C + i;
    p->val[l] = 2*(p->val[l]) - p->val[l+Nx_C];
  }
  
  //North side Dirchlet

  j = Ny_C - 1;
  for(i=0;i<Nx_C;i++)
  {
    l = j*Nx_C + i;
    p->val[l] = 2*(p->val[l]) - p->val[l-Nx_C];
  }
  
  // Loops for all the inner cells 
  
  for(j=1;j<(Ny_C-1);j++)
  {
    row = j*Nx_C;

    val_pP = &p->val[row];
    val_pE = &p->val[row+1];
    val_pW = &p->val[row-1];
    val_pN = &p->val[row+Nx_C];
    val_pS = &p->val[row-Nx_C];

    for(i=1;i<(Nx_C-1);i++)
    {
      l = j*Nx_C + i;

      pP = val_pP[i];
      pE = val_pE[i];
      pW = val_pW[i];
      pN = val_pN[i];
      pS = val_pS[i];

      aE = ( 1 );  
      aW = ( 1 );  
      aN = ( 1 );  
      aS = ( 1 );
      aP = aE + aW + aN + aS ;

      Ap_vector[l] = -pP*aP + pE*aE + pW*aW + pN*aN + pS*aS; 
    }
  }

  return;
}  
