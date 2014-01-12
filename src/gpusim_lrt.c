//
//  gpusim_lrt.c
//  LRT
//
//  Copyright (c) 2013 The University of Sydney. All rights reserved.
//
//  naive lrt computation with high complexity and high memory usage 
//

#include "lrt.h"
#include "heap.h"
#include "prefix.h"

////////////////////////////////////////////////////////////////////////////
//  functions
//

//
// comparison function for qsort that compares the scores of two rectangles.
//
static int rectangle_cmp(rectangle *r_a, rectangle *r_b)
{
   return  - rectangle_compare(r_a,r_b);
}

//
// sinks the first element of the min-heap until heap condition holds
//
static void rheap_sink(rectangle *heap, size_t n)
{
   // allocate memory for swap on heap
   rectangle tmp;
   int k;

   for (k=0;2*k<n;) {

      // set j to first child of k
      int j = 2 * k + 1;

      // check which of the children (if there are two) is smaller
      // move j to the smaller child so that the smaller child is 
      // swapped with a[k]
      if (j < n-1) {
         if (rectangle_cmp(&heap[j+1],&heap[j]) < 0) {
            j++;
         }
      }

      // check whether min-mheap condition holds, i.e a[k] is smaller 
      // than smallest child. If so, stop.
      if (rectangle_cmp(&heap[k],&heap[j]) < 0)
         break;

      // swap elements k and j 
      tmp = heap[k]; 
      heap[k] = heap[j];
      heap[j] = tmp;

      // set new n
      k = j;
   }
}

//
// raises the last element of the min-heap until heap condition holds
//
static void rheap_rise(rectangle *heap, size_t n)
{
   // allocate memory for swap on heap
   rectangle tmp;
   int k;

   for (k=n-1;k > 0; ) { 
       // set j to parent
       int j = k/2;

       // check whether a[k] is smaller than parent a[j] 
       if (rectangle_cmp(&heap[k],&heap[j]) < 0) { 
          // swap elements k and j 
          tmp = heap[k];
          heap[k] = heap[j];
          heap[j] = tmp; 

          // set k to current node 
          k = j; 
       } else break;
   } 
}

//
// gpu simulation of LRT computation
// 
rectangle *gpusim_lrt(grid *g,int kbest)
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

   // define grid size for a thread
   int gx_size = 4; 
   int gy_size = 4; 

   // rectangular length
   int rx_size = (width + 1)/2;
   int ry_size = (height + 1)/2;

   // compute the number of iterations necessary for the loops of a thread
   // (for each dimension)
   int tx_size=((width+1)*rx_size + bx_size*gx_size-1)/(bx_size*gx_size);
   int ty_size=((height+1)*ry_size + by_size*gy_size-1)/(by_size*gy_size);
 

   // allocate a heap for each thread with size k 
   rectangle *r_global = ALLOCV(rectangle,kbest*bx_size*by_size*gx_size*gy_size);
   int *heapsize_global = ALLOCV(int,bx_size*by_size*gx_size*gy_size);

   int gx,gy,    // index variables for grid and block 
       bx,by;    // loops to simulate GPGPU

   // enumerate all blocks 
   for (gx=0;gx<gx_size;gx++){
      for (gy=0;gy<gy_size;gy++){

         // enumerate all threads in a block. 
         for (bx=0;bx<bx_size;bx++){
            for (by=0;by<by_size;by++){

               //////////////////////////////////////////////////
               // begin of kernel 
               //////////////////////////////////////////////////

               // allocate heap with size k for each thread  locally
               // in cuda replace with "rectangle r[kbest];"
               rectangle *r = ALLOCV(rectangle,kbest);
   
               // heap size (initialized with zero and grows up to kbest)
               int heap_size=0;  

               // thread loops 
               int x,y;
               for (x=0;x<tx_size;x++) { 
                  for (y=0;y<ty_size;y++) { 
                  
                     // convert grid and block coordinates back to 
                     // a four dimensional grid.
                     int grid_x=(gx*bx_size+bx)*tx_size+x;
                     int grid_y=(gy*by_size+by)*ty_size+y;
 
                     // check whether thread is out of range 
                     // (if block sizes don't divide grid length) 
                     if(grid_x>=(width+1)*rx_size) continue;
                     if(grid_y>=(height+1)*ry_size) continue;
              
                     // get the four dimensional rectangular 
                     // coordinates using inverse Horner scheme. 
                     int ti1=(grid_y)%ry_size;
                     int ti2=(grid_y)/ry_size;
                     int tj1=(grid_x)%rx_size;
                     int tj2=(grid_x)/rx_size;

                     assert(0 <= ti1 && ti1 < ry_size && 0 <= ti2 && ti2 < height+1 && "Rectangular coordinate for rows out of range");
                     assert(0 <= tj1 && tj1 < rx_size && 0 <= tj2 && tj2 < width+1 && "Rectangular coordinate for columns out of range");

                     // transform rectangular coordinates to 
                     // triangular coordinates. 
                     int i1,j1,i2,j2;
                     if ((tj2 < width - tj1)) {
                        j1=tj1;
                        j2=tj2+tj1;
                     } else if( (tj1+1)*2 < width+1){
                        j1=width-tj1-1;
                        j2=width-tj2+j1;
                     } else continue;
                     if ((ti2 < height - ti1)) {
                         i1=ti1;
                         i2=ti2+ti1;
                     } else if( (ti1+1)*2 < height+1){
                         i1=height-ti1-1;
                         i2=height-ti2+i1;
                     } else continue;

                     assert(0 <= i1 && i1 <= i2 && i2 < height && "Triangular coordinate for rows out of range");
                     assert(0 <= j1 && j1 <= j2 && j2 < width && "Triangular coordinate for columns out of range");
 
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
                        r[heap_size].i1 = i1;
                        r[heap_size].j1 = j1;
                        r[heap_size].i2 = i2;
                        r[heap_size].j2 = j2;
                        heap_size = heap_size + 1;
                        rheap_rise(r,heap_size);
                     // otherwise check whether the current score is greater than the smallest element (=first element) 
                     // if so, replace the first element with current rectangle and sink the current element until 
                     // the heap condition holds. 
                     } else if (heap_rcmp(&current,r) > 0) {
                        r[0] = current;
                        rheap_sink(r,kbest);
                     } 
                  }
               }

               // store heap_size
               heapsize_global[((gy*by_size+by)*gx_size*bx_size+gx*bx_size+bx)] = heap_size;

               // Write heap back to global memory
               memcpy(&r_global[((gy*by_size+by)*gx_size*bx_size+gx*bx_size+bx)*kbest],r,heap_size*sizeof(rectangle));

               //////////////////////////////////////////////////
               // end of kernel 
               //////////////////////////////////////////////////
            }
         }
      }
   }

   ////////////////////////////////////////////////////////////////////////////////////
   // This part needs to be replaced by a second kernel / lookup Mark Harris slides 
   // about reduction.
   rectangle *r = ALLOCV(rectangle,kbest);
   int heap_size = 0;

   // enumerate all blocks 
   for (gx=0;gx<gx_size;gx++){
      for (gy=0;gy<gy_size;gy++){


         // enumerate all threads in a block. 
         for (bx=0;bx<bx_size;bx++){
            for (by=0;by<by_size;by++){

               int hs = heapsize_global[((gy*by_size+by)*gx_size*bx_size+gx*bx_size+bx)];
               rectangle *rs = &r_global[((gy*by_size+by)*gx_size*bx_size+gx*bx_size+bx)*kbest];
               int i;
               for (i=0;i<hs;i++) {
                  if (heap_size < kbest) { 
                     r[heap_size] = rs[i];
                     heap_size = heap_size + 1;
                     rheap_rise(r,heap_size);
                  } else if (heap_rcmp(&rs[i],r) > 0) {
                     r[0] = rs[i];
                     rheap_sink(r,kbest);
                  } 
               }
            }
         }
      } 
   } 

   // sort result and return it
   qsort(r,heap_size,sizeof(rectangle),rectangle_compare);  

   return r; 
}
