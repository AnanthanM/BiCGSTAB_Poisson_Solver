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

    if( i == 0 || i == (N_Cells_x_p-1) || j == 0 || j == (N_Cells_y_p-1) )
    {
      p->bc_type[l]   = INSIDE_GC;       // means inside Ghost cell
    }
    else
    {
      p->bc_type[l]   = NONE; 
    }
  }
  
  
  return;

}

// We call this function to set the Boundary cell
// values for ghost cells for any input Field

void set_ghost_cells_value(Field * phi)
{
  int i,j,l;

  int N_Cells_x = phi->N_x;
  int N_Cells_y = phi->N_y;
  int N_Cells   = N_Cells_x*N_Cells_y;
  
  for(l=0;l<N_Cells;l++)
  {
    i = l%N_Cells_x;
    j = (int) l/N_Cells_x;
    
    if(phi->bc_type[l] != NONE)
    {
      if(i==0)
        phi->val[l] = phi->BC_Value[XMIN];
      else if(i==(N_Cells_x-1))
        phi->val[l] = phi->BC_Value[XMAX];
      if(j==0)
        phi->val[l] = phi->BC_Value[YMIN];
      else if(j==(N_Cells_y-1))
        phi->val[l] = phi->BC_Value[YMAX];
    }
  }

  return;
}
