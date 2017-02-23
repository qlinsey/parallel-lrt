//
//  prefix.h
//  LRT
//
//  Copyright (c) 2013 The University of Sydney. All rights reserved.
//

#ifndef __PREFIX_H__
#define __PREFIX_H__

#include "lrt.h"

////////////////////////////////////////////////////////////////////////////
//  type definitions
//

//
// prefix sums for a data grid
typedef struct prefix_sums {
  grid *gridref;    // grid
  cell *A,       // all four corners 
       *B, 
       *X, 
       *Y;  
} prefix_sums;

////////////////////////////////////////////////////////////////////////////
// function prototypes
//

//
// compute prefix sum of a grid
//
prefix_sums *compute_psums(grid *);   

//
// compute prefix sum of a grid
//
void free_psums(prefix_sums *);   

#endif

