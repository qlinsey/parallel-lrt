//
//  naive_lrt.c
//  LRT
//
//  Copyright (c) 2013 The University of Sydney. All rights reserved.
//
//  naive lrt computation with high complexity and high memory usage 
//

#include "lrt.h"

////////////////////////////////////////////////////////////////////////////
//  functions
//

//
// naive implementation of LRT
// 
// the complexity of the naive implementation is O(n^6) and 
// hence not practical for larger grids. 
// 
// Steps:
//   1. allocate memory for resulting rectangles and their scores
//   2. compute grid totals
//   2. compute lrt for each rectangle
//     2.a) compute rectangle total
//     2.b) compute score
//   3. sort all rectangles according to their score 
// 
rectangle *naive_lrt(grid *g,int kbest)
{
   int idx=0;           // output index 
   int i,j;             // index for regions in grid 

   // set width and height
   int width=g->width;
   int height=g->height;
   
   // compute number of rectangles in the grid and allocate 
   // memory to store result 
   size_t size = (size_t)(width + 1) *
          (size_t)width  * 
          (size_t)(height + 1) * 
          (size_t)height / 4;
   rectangle *r = (rectangle *) ALLOCV(rectangle,size);
   
   // compute grid totals, ratio, and likelihood
   int g_n = 0;
   int g_k = 0;
   for (i=0;i<height;i++) {
      for (j=0;j<width;j++) {
         g_n += g->cells[I(i,j)].n;
         g_k += g->cells[I(i,j)].k;
      }
   }
   float g_q = (float)g_k / (float)g_n;
   float g_l = g_k * log(g_q);
   
   // compute LRT for all rectangle in grid
   int j1,i1,j2,i2;  // coordinates of a rectangle
   for (i1=0;i1<height;i1++){ 
       for (j1=0;j1<width;j1++){ 
          for (i2=i1;i2<height;i2++){ 
             for (j2=j1;j2<width;j2++){ 
               
               // set rectangle width and height 
               int r_width = j2 - j1 + 1; 
               int r_height = i2 - i1 + 1; 
               
               // compute rectangle totals, ratio and likelihood
               int r_n=0; 
               int r_k=0; 
               cell *c = &g->cells[i1*width+j1]; 
               for (i=0;i<r_height;i++) { 
                  for (j=0;j<r_width;j++) { 
                     r_n += c[I(i,j)].n; 
                     r_k += c[I(i,j)].k; 
                  }  
               } 
               float r_q = (float)r_k / (float)r_n;
               float r_l = r_k * log(r_q) - r_k; 
               
               // compute rectangle's complement totals, ratio, and likelihood
               int   c_n = g_n - r_n;
               int   c_k = g_k - r_k;
               float c_q = (float)c_k / (float)c_n;
               float c_l = c_k * log(c_q) - c_k;
               
               // compute score 
               float score = r_l + c_l - g_l; 

               // store result in array r
               r[idx].score = score; 
               r[idx].i1 = i1;
               r[idx].j1 = j1;
               r[idx].i2 = i2;
               r[idx].j2 = j2;
               idx = idx + 1;
            }
         } 
      }
   }
   
   // sort result and return it
   qsort(r,size,sizeof(rectangle),rectangle_compare);  
   return r; 
}
