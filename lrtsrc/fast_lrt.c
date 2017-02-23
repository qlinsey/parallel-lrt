//
//  fast_lrt.c
//  LRT
//
//  Copyright (c) 2013 The University of Sydney. All rights reserved.
//
//  naive lrt computation with high complexity and high memory usage 
//

#include "lrt.h"
#include "heap.h"
#include "prefix.h"

//
// fast implementation of LRT using prefix sums
//
rectangle *fast_lrt(grid *g,int kbest)
{
   // heap size (initialized with zero and grows up to kbest)
   int heap_size=0;  

   // set width and height
   int width = g->width;
   int height = g->height;

   // allocate heap with size k
   rectangle *r = (rectangle *) ALLOCV(rectangle,kbest);
   
   // compute prefix sums
   prefix_sums *p = compute_psums(g);

   // compute grid totals, ratio, and likelihood
   int g_n = p->A[width * height - 1].n;
   int g_k = p->A[width * height - 1].k;
   float g_q = (float)g_k / (float)g_n;
   float g_l = g_k * log(g_q);
   
   // compute LRT for each rectangle
   int j1,i1,j2,i2;                       // coordinates of a rectangle
   for (i1=0;i1<height;i1++){ 
      for (j1=0;j1<width;j1++){
         for (i2=i1;i2<height;i2++){ 
            for (j2=j1;j2<width;j2++){ 
               
               // compute rectangle totals, ratio and likelihood
               int   a_n = p->A[I(i2,j2)].n;
               int   a_k = p->A[I(i2,j2)].k;
               int   b_n = p->B[I(i1,j1)].n;
               int   b_k = p->B[I(i1,j1)].k;
               int   y_n = p->Y[I(i1,j2)].n;
               int   y_k = p->Y[I(i1,j2)].k;
               int   x_n = p->X[I(i2,j1)].n;
               int   x_k = p->X[I(i2,j1)].k;
               int   r_n = a_n + b_n + x_n + y_n - g_n;
               int   r_k = a_k + b_k + x_k + y_k - g_k;
               float r_q = (float)r_k / (float)r_n;
               float r_l = r_k * log(r_q) - r_k; 
               
               // check whether prefix sums are correct
               #if DEBUG_LEVEL > 0
               int t_n=0; 
               int t_k=0; 
               int i,j;
               int r_width = j2 - j1 + 1; 
               int r_height = i2 - i1 + 1; 
               cell *c = &g->cells[i1*width+j1]; 
               for (i=0;i<r_height;i++) { 
                  for (j=0;j<r_width;j++) { 
                     t_n += c[I(i,j)].n; 
                     t_k += c[I(i,j)].k; 
                  }  
               } 
               assert((t_n == r_n && t_k == r_k) && "problems with prefix");
               #endif

               // compute rectangle's complement totals, ratio, and likelihood
               int   c_n = g_n - r_n;
               int   c_k = g_k - r_k;
               float c_q = (float)c_k / (float)c_n;
               float c_l = c_k * log(c_q) - c_k;
               
               // compute score 
               float score = r_l + c_l - g_l; 

               // populate current rectangle
               rectangle current;
               current.score = score; 
               current.i1 = i1;
               current.j1 = j1;
               current.i2 = i2;
               current.j2 = j2;
                  
               // store result in heap
               // if the heap size is still smaller than kbest, add rectangle to the end of the heap
               // and rise the last element until the heap condition holds 
               if (heap_size < kbest) { 
                  r[heap_size].score = score; 
                  r[heap_size].j1 = j1;
                  r[heap_size].i1 = i1;
                  r[heap_size].j2 = j2;
                  r[heap_size].i2 = i2;
                  heap_size = heap_size + 1;
                  minheap_rise(r,heap_size,sizeof(rectangle),heap_rcmp);
               // otherwise check whether the current score is greater than the smallest element (=first element) 
               // if so, replace the first element with current rectangle and sink the current element until 
               // the heap condition holds. 
               } else if (heap_rcmp(&current,r) > 0) {
                  r[0] = current;
                  minheap_sink(r,kbest,sizeof(rectangle),heap_rcmp);
               } 
            }
         } 
      }
   }

   // sort result and return it
   qsort(r,kbest,sizeof(rectangle),rectangle_compare);

   // free prefix sums
   free_psums(p); 

   return r; 
}
