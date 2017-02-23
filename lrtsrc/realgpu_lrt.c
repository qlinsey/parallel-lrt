//
//  fgpu_lrt.c
//  LRT
//
//  Copyright (c) 2013 The University of Sydney. All rights reserved.
//
//  gpu lrt computation with high complexity and high memory usage
//

#include "lrt.h"
#include "heap.h"
#include "prefix.h"
#include <cuda.h>
// prototype wrapper function

extern void wrapper_lrt(prefix_sums *, rectangle *, int *, int , int, int, int ,int , int, int ,float, float);

//
// fast implementation of LRT using prefix sums
//
rectangle *gpu_lrt(grid *g,int kbest)
{
   // set width and height
   int width = g->width;
   int height = g->height;

   // compute prefix sums
   prefix_sums *p = compute_psums(g);
   // compute grid totals, ratio, and likelihood
   int g_n = p->A[width * height - 1].n;
   int g_k = p->A[width * height - 1].k;
   float g_q = (float)g_k / (float)g_n;
   float g_l = g_k * log(g_q);

   // define block size for grid
   int bx_size = 2;
   int by_size = 2;


   // define tile size for a thread
   int tx_size = 2;
   int ty_size = 2;

   // rectangular length
   int rx_size = (width + 1)/2;
   int ry_size = (height + 1)/2;

   // A [(N+1],m]^2 grid needs to be mappeds to a two-dimensional
   // grid [(N+1)*m]^2 and divided by the block size. Ensure that if the
   // block size does not divide grid length, one extra block is added
   // (NB: This is done by adding bx_size-1 before division).

   int gx_size=((width+1)*rx_size + bx_size*tx_size-1)/(bx_size*tx_size);
   int gy_size=((height+1)*ry_size + by_size*ty_size-1)/(by_size*ty_size);

   int gx,gy,    // index variables for grid and block
       bx,by;    // loops to simulate GPGPU

   // allocate a heap for each thread with size k
   rectangle *r_global = (rectangle *)ALLOCV(rectangle,kbest*bx_size*by_size*gx_size*gy_size);
   int *heapsize_global = (int *)ALLOCV(int,bx_size*by_size*gx_size*gy_size);
   wrapper_lrt(p,r_global,heapsize_global,tx_size,ty_size,width,height,kbest,g_n,g_k,g_q,g_l);

}
