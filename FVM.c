// A Poisson Equation Solver

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>

#include <omp.h>

#include "fvm.h"


int main(int argc, char *argv[])
{
  double start_wall_time = omp_get_wtime();

  Constant constant;
  Domain   domain; 

  double L_x = 1.0;                                                    // Domain Dimensions
  double L_y = 1.0;
  
  int N_cells_x   = 100+2;                                            //For the entire domain for colocated variables
  int N_cells_y   = 100+2;                                            //For the entire domain for colocated variablesy
  int N_cells_z   = 1;
  int N_cells     = N_cells_x * N_cells_y * N_cells_z ;               //Total Number of Colocated cells

  /******Initialising the Constants **************/
  
  constant.L_x = L_x;
  constant.L_y = L_y; 
  constant.dx  = L_x/(N_cells_x-2);
  constant.dy  = L_y/(N_cells_y-2);
  constant.dz  = 1.0;

  /*****Initialising the Field Structs**********/

  // p   -> Pressure
  

  Field * p    = Allocate_Field( N_cells_x, N_cells_y, N_cells_z );
  Field * RHS  = Allocate_Field( N_cells_x, N_cells_y, N_cells_z );
  
  Bicgstab * PP = Allocate_Bicgstab( N_cells ); 

  /*****Initialising the Domain struct***************/

  domain.p   = p;
  domain.PP  = PP;

  /****Ghost cell types are set for each Field*********/

  set_ghost_cells_type(domain);

  /****Boundary Cell values are set for each Field******/

  int i;
// i = 0 -> XMIN -> west
// i = 1 -> XMAX -> east
// i = 2 -> YMIN -> south
// i = 3 -> YMAX -> north
  for(i=0;i<4;i++)
  {
    p->BC_Value[i] = 0.0;
  }
  
  p->BC_Value[YMIN] = 1.0;
  p->BC_Value[YMAX] = 0.0;

  for(i =0;i < N_cells;i++)
  {
    p->val[i]   = 0.0; 
  }
  
  set_ghost_cells_value(p);

  /******STARTING*******/  
  Write_VTK(0,domain,constant);                    //Writing Output at the time t=0
  
  int si_no     = 1;
  int test;      

  /******Solving Pressure Poisson ********/
  Get_RHS_of_Pressure_Poisson_Eqn(domain,RHS,constant);
  test = solve_Pressure_Poisson_BiCGSTAB(p,RHS,constant,domain);
  if(test)
  {
    printf("\n Problem Solving Pressure Poisson Equation\n");
  }
  
  set_ghost_cells_value(p);  // This call is only needed if Dirchlet BCs are there
                             // or When writing VTK output
  
  /*********Writing Output*************************/
  Write_VTK(si_no,domain,constant); 
  printf("VTK file is written \n");

  printf("Time taken for the simulation: %lf", omp_get_wtime() - start_wall_time);

  return 0;
}



















