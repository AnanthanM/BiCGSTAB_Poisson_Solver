/****************************************************************************
 *   
 *   Boundary related functions where we set which type of ghost cell 
 *   is there for Different Fields
 *   
 *****************************************************************************/   

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>

#include "fvm.h"

// It takes input as the initialissed domain struct then assigns what kind of
// Boundary type  is there in each ghost cell

void set_ghost_cells_type(Domain domain)
{
  int i,j,l;

  Field * p   = domain.p;

  //*******GHOSTS CELLS FOR COLOCATED VARIABLES***********
  
  int N_Cells_x_p = p->N_x;
  int N_Cells_y_p = p->N_y;
  int N_Cells_p   = N_Cells_x_p*N_Cells_y_p;
  
  for(l=0;l<N_Cells_p;l++)
  {
    i = l%N_Cells_x_p;
    j = (int) l/N_Cells_x_p;

    if( j == 0 || j == (N_Cells_y_p-1) )
    {
      p->bc_type[l]   = DIRCHLET;       // N and S values given
    }
    else if( i == 0 || i == (N_Cells_x_p-1) )
    {
      p->bc_type[l]   = NEUMANN;       // E and W side Neumann
    }
    else
    {
      p->bc_type[l]   = NONE; 
    }
  }
    
  return;

}

// We call this function to set the Boundary cell
// values for ghost cells for any input Field having
// Dirchelet Boundary Condition

void set_ghost_cells_value(Field * phi)
{
  int i,j,l;

  int N_Cells_x = phi->N_x;
  int N_Cells_y = phi->N_y;
  int N_Cells   = N_Cells_x*N_Cells_y;
  
  //BC Value Setup
  
  //West Side Neumann 
  
  i = 0;
  for(j=0;j<N_Cells_y;j++)
  {
    l = j*N_Cells_x + i;
    phi->val[l] = 0.0;
  }

  //East Side Neumann
  
  i = N_Cells_x - 1;
  for(j=0;j<N_Cells_y;j++)
  {
    l = j*N_Cells_x + i;
    phi->val[l] = 0.0;
  }
  
  //South side Dirchlet

  j = 0;
  for(i=0;i<N_Cells_x;i++)
  {
    l = j*N_Cells_x + i;
    phi->val[l] = phi->BC_Value[YMIN];
  }
  
  //North side Dirchlet

  j = N_Cells_y - 1;
  for(i=0;i<N_Cells_x;i++)
  {
    l = j*N_Cells_x + i;
    phi->val[l] = phi->BC_Value[YMAX];
  }

  return;
}

