/*************** Helper Functions *****************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>

#include "fvm.h"

// Allocate Field takes 3 integers as input and initialise 
// a Field with those values as Nx, Ny and Nz points
// in x,y and z direction.It returs the adress of the field.

Field * Allocate_Field(int N_x,int N_y,int N_z)
{
  Field * phi;

  phi           = malloc(sizeof(Field));

  phi->N_x      = N_x;
  phi->N_y      = N_y;
  phi->N_z      = N_z;
  phi->N        = N_x*N_y*N_z;
  phi->bc_type  = malloc(phi->N * sizeof(BC_Type)); 
  phi->val      = malloc(phi->N * sizeof(double));

  return phi;
}

// Allocates struct Bicgstab which contains arrays necessary
// for BiCGSTAB Solver, Takes input as an integer

Bicgstab * Allocate_Bicgstab(int N)
{
  Bicgstab * phi;

  phi             = malloc(sizeof(Bicgstab));

  phi->r0_cap     = malloc(N*sizeof(double));
  phi->sj         = malloc(N*sizeof(double));
  phi->rj         = malloc(N*sizeof(double));
  phi->pj         = malloc(N*sizeof(double));
  phi->pcap       = malloc(N*sizeof(double));
  phi->scap       = malloc(N*sizeof(double));
  phi->Ax_vector  = malloc(N*sizeof(double));
  phi->h          = malloc(N*sizeof(double));
  phi->vj         = malloc(N*sizeof(double));

  return phi;
}
