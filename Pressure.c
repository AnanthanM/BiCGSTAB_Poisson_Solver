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
    if(p->bc_type[l] == NONE)
    {
        i = l%Nx_C;
        j = (int) l/Nx_C;

        RHS->val[l] = 0.0;
    }
    else
    {
        RHS->val[l] = 0.0;
    }
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
  
  for(l = 0;l<N_cells_C;l++)
  {
    if(p->bc_type[l] == NONE)
    {
        i = l%Nx_C;
        j = (int) l/Nx_C;
  
        S = (j-1)*Nx_C + i;
        N = (j+1)*Nx_C + i;
        W = j*Nx_C + (i-1);      // or l-1
        E = j*Nx_C + (i+1);      // or l+1
         
        aE = ( 1 );  
        aW = ( 1 );  
        aN = ( 1 );  
        aS = ( 1 );
        aP = aE + aW + aN + aS ;
 
        pP =p->val[l];
        pE =p->val[E];
        pW =p->val[W];
        pN =p->val[N];
        pS =p->val[S];

        if(p->bc_type[E] != NONE)
          pE = pP;

        if(p->bc_type[W] != NONE)
          pW = pP;

        if(p->bc_type[N] != NONE)
          pN = 2*p->val[N] - pP;        //Dirchelet

        if(p->bc_type[S] != NONE)
          pS = 2*p->val[S] - pP;        //Dirchelet
        
        Ap_vector[l] = -pP*aP + pE*aE + pW*aW + pN*aN + pS*aS; 

    }
    else
    {
        Ap_vector[l] = 0.0;
    }
    
  }

  return;

}  
