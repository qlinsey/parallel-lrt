//
//  gpu_lrt.cu
//  LRT
//
//  Copyright (c) 2013 The University of Sydney. All rights reserved.
//
//  GPU implementation of LRT in CUDA
//


#include <cuda.h> 
#include <assert.h> 
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string>
#include <sstream>
#include <stdexcept>
#include <iostream> 

extern "C" {
#include "lrt.h"
#include "prefix.h"
#include "gpu.h"
#include "heap.h"
}

//
// compare two rectangles. If the scores are equal, the coordinates
// are used to create a total order 
// 
__device__ __forceinline__ int rectangle_compare(rectangle *r_a, rectangle *r_b)
{

   if(isnan(r_a->score)) return -1;
   if(isnan(r_b->score)) return 1;
   
   float score_a = r_a->score;
   float score_b = r_b->score;

   if (score_a < score_b) return -1;
   if (score_a > score_b) return 1;
   
   int   area_a=(r_a->j2-r_a->j1+1)*(r_a->i2-r_a->i1+1);
   int   area_b=(r_b->j2-r_b->j1+1)*(r_b->i2-r_b->i1+1);

   if(area_a < area_b) return 1;
   if(area_a > area_b) return -1;
 
   if(r_a->i1 < r_b->i1) return 1;
   if(r_a->i1 > r_b->i1) return -1;

   if(r_a->i2 < r_b->i2) return 1;
   if(r_a->i2 > r_b->i2) return -1;

   if(r_a->j1 < r_b->j1) return 1;
   if(r_a->j1 > r_b->j1) return -1;

   if(r_a->j2 < r_b->j2) return 1;
   if(r_a->j2 > r_b->j2) return -1;

   return 0;
}

//
// raises the last element of the min-heap until heap condition holds
//
__device__  __forceinline__ void rheaprise(rectangle *heap, int n)
{
   // allocate memory for swap on heap
   rectangle tmp;
   int k;

   for (k=n-1;k > 0; ) { 
       // set j to parent
       int j = k/2;
       // check whether a[k] is smaller than parent a[j] 
       if (rectangle_compare(&heap[k],&heap[j]) < 0) { 
          // swap elements k and j 
          tmp=heap[k];
          heap[k]=heap[j];
          heap[j]=tmp;

          // set k to current node 
          k = j; 
       } else  break;
   }   
}

//
// sinks the first element of the min-heap until heap condition holds
//
__device__  __forceinline__ void rheapsink(rectangle *heap, int n)
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
         if (rectangle_compare(&heap [j+1], &heap[j]) < 0) {
            j++;
         }
      }
      // check whether min-mheap condition holds, i.e a[k] is smaller 
      // than smallest child. If so, stop.
      if (rectangle_compare(&heap[k], &heap[j]) <0)
         break;
      // swap elements k and j 
       tmp=heap[k];
       heap[k]=heap[j];
       heap[j]=tmp;
       // set new n
       k = j;
   }
}


//
// execute LRT (fully enumeration mapping+heap) on GPGPU
//
__global__ void lrt_kernel(int width, 
                           int height,
                           int kbest, 
                           cell *pA, 
                           cell *pB, 
                           cell *pX,  
                           cell *pY, 
                           rectangle *r_global, 
                           int *heapsize_global
                          )
{
   // compute grid totals, ratio, and likelihood
   int g_n = pA[width * height - 1].n;
   int g_k = pA[width * height - 1].k;
   float g_q = (float)g_k / (float)g_n;
   float g_l = g_k * log(g_q);

   // rectangular length
   int rx_size = (width + 1)/2;
   int ry_size = (height + 1)/2;

   // compute the number of iterations necessary for the loops of a thread
   // (for each dimension)
   int tx_size=((width+1)*rx_size + blockDim.x*gridDim.x-1)/(blockDim.x*gridDim.x);
   int ty_size=((height+1)*ry_size + blockDim.y*gridDim.y-1)/(blockDim.y*gridDim.y);
 
   // allocate heap with size k for each thread  locally
   extern __shared__ rectangle sr[];
   rectangle *r = &sr[(threadIdx.y * blockDim.x + threadIdx.x)*kbest]; 
   
   // heap size (initialized with zero and grows up to kbest)
   int heap_size=0;  

   // thread loops 
   int x,y;
   for (x=0;x<tx_size;x++) { 
      for (y=0;y<ty_size;y++) { 
      
         // convert grid and block coordinates back to 
         // a four dimensional grid.
         int grid_x=(blockIdx.x*blockDim.x+threadIdx.x)*tx_size+x;
         int grid_y=(blockIdx.y*blockDim.y+threadIdx.y)*ty_size+y;

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

         // compute rectangle totals, ratio and likelihood
         int   a_n = pA[I(i2,j2)].n;
         int   a_k = pA[I(i2,j2)].k;
         int   b_n = pB[I(i1,j1)].n;
         int   b_k = pB[I(i1,j1)].k;
         int   y_n = pY[I(i1,j2)].n;
         int   y_k = pY[I(i1,j2)].k;
         int   x_n = pX[I(i2,j1)].n;
         int   x_k = pX[I(i2,j1)].k;
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
            rheaprise(r,heap_size);
         // otherwise check whether the current score is greater than the smallest element (=first element) 
         // if so, replace the first element with current rectangle and sink the current element until 
         // the heap condition holds. 
         } else if (rectangle_compare(&current,r) > 0) {
            r[0] = current;
            rheapsink(r,kbest);
         } 
      }
   }

   // store heap_size
   heapsize_global[((blockIdx.y*blockDim.y+threadIdx.y)*gridDim.x*blockDim.x+blockIdx.x*blockDim.x+threadIdx.x)] = heap_size;

   // Write heap back to global memory
   memcpy(&r_global[((blockIdx.y*blockDim.y+threadIdx.y)*gridDim.x*blockDim.x+blockIdx.x*blockDim.x+threadIdx.x)*kbest],r,heap_size*sizeof(rectangle));
}


//
// gpu simulation of LRT computation
// 
__host__ rectangle *gpu_lrt(grid *g,int kbest)
{

   // set width and height
   int width = g->width;
   int height = g->height;
   
   // compute prefix sums
   prefix_sums *p = compute_psums(g);

   // allocate memory for prefix sums 
   // and transfer them to the CUDA device
   cell *pa_d,
        *pb_d, 
        *px_d, 
        *py_d;
   size_t num_cells=sizeof(cell)*width*height;
   cudaMalloc(&pa_d,num_cells);
   cudaMalloc(&pb_d,num_cells);
   cudaMalloc(&px_d,num_cells);
   cudaMalloc(&py_d,num_cells);
   cudaMemcpy(pa_d,p->A,num_cells,cudaMemcpyHostToDevice);
   cudaMemcpy(pb_d,p->B,num_cells,cudaMemcpyHostToDevice);
   cudaMemcpy(px_d,p->X,num_cells,cudaMemcpyHostToDevice);
   cudaMemcpy(py_d,p->Y,num_cells,cudaMemcpyHostToDevice);
  
   // allocate memory for resulting heaps for each thread 
   // and heap sizes.  
   rectangle *r_global_d; 
   cudaEvent_t start,stop;
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   cudaMalloc(&r_global_d,sizeof(rectangle)*kbest*THREADS_IN_GRID);
   int *heapsize_global_d; 
   cudaMalloc(&heapsize_global_d,sizeof(int)*THREADS_IN_GRID);
   cudaEventRecord(start);
   // run kernel
   dim3 dimg(GRIDSIZE_X,GRIDSIZE_Y);
   dim3 dimb(BLOCKSIZE_X,BLOCKSIZE_Y);
   lrt_kernel<<<dimg,dimb,sizeof(rectangle)*BLOCKSIZE_X*BLOCKSIZE_Y*kbest>>>(
      width,
      height,
      kbest,
      pa_d, 
      pb_d,
      px_d,
      py_d,  
      r_global_d,
      heapsize_global_d);
   cudaEventRecord(stop);
   // retrieve result from GPGPU
   rectangle *r_global = ALLOCV(rectangle, sizeof(rectangle)*kbest*THREADS_IN_GRID);
   int *heapsize_global = ALLOCV(int, sizeof(int) * THREADS_IN_GRID); 
   cudaMemcpy(r_global,r_global_d,sizeof(rectangle)*kbest*THREADS_IN_GRID,cudaMemcpyHostToDevice);
   cudaMemcpy(heapsize_global,heapsize_global_d,sizeof(int)*THREADS_IN_GRID,cudaMemcpyHostToDevice);
   cudaEventSynchronize(stop);
   float mseconds=0;
   cudaEventElapsedTime(&mseconds,start,stop);
   printf("lrt kernel computation time=%f (s)\n",((double)mseconds/1000));
  
   ////////////////////////////////////////////////////////////////////////////////////
   // join heaps of all threads to a single heap 
   rectangle *r = ALLOCV(rectangle,kbest);
   int heap_size = 0;

   // enumerate all blocks 
   int i;
   for (i=0;i<THREADS_IN_GRID;i++) { 
      int hs = heapsize_global[i];
      rectangle *rs = &r_global[i*kbest];
      int j;
      for (j=0;j<hs;j++) {
         if (heap_size < kbest) { 
            r[heap_size] = rs[j];
            heap_size = heap_size + 1;
            minheap_rise(r,heap_size,sizeof(rectangle),heap_rcmp);
         } else if (heap_rcmp(&rs[j],r) > 0) {
            r[0] = rs[j];
            minheap_sink(r,kbest,sizeof(rectangle),heap_rcmp);
         } 
      } 
   } 

   // sort result and return it
   qsort(r,heap_size,sizeof(rectangle),rectangle_compare);  

   return r; 
}
