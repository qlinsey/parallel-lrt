//
//  real gpu_lrt.cu
//  LRT
//
//  Copyright (c) 2013 The University of Sydney. All rights reserved.
//
//  gpu lrt computation with range mapping scheme and heap reduction 
//


#include <cuda.h> 
#include <assert.h> 
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
//#include "util/cuPrintf.cu" 
#include <string>
#include <sstream>
#include <stdexcept>
#include <iostream> 


using namespace std;

#ifndef __Kernel_CU__
#define __Kernel_CU__
__global__void kernel(int a ){
  int x=threadIdx.x+a;
  int y=threadIdx.y+a;

}
extern "C" void wrapper_lrt(){

   kernel<<<2,2>>>(2);

}
#endif
////////////////////////////////////////////////////////////////////////////
//  device functions
//

//
// sinks the first element of the min-heap until heap condition holds
//

/*__device__ void minheap_sink(int t){


}*/

//
// fast GPGPU implementation of LRT using prefix sums 
//
/*__global__ void gpu_lrt(prefix_sums *p, rectangle * r_global, int * heapsize_global,  int tx_size, int ty_size, int width, int height, int kbest, int g_n, int g_k, float g_q, float g_l) {
   // rectangular length
   int rx_size = (width + 1)/2;
   int ry_size = (height + 1)/2;
   int heap_size=0;
   int x,y;
   // allocate heap with size k for each thread  locally
   // in cuda replace with "rectangle r[kbest];"
   //rectangle* r =(rectangle *)malloc(sizeof(rectangle)*kbest);
   rectangle r[10];
   for(x=0;x<tx_size;x++) { 
      for (y=0;y<ty_size;y++) { 
         // convert grid and block coordinates back to 
         // a four dimensional grid.
         int grid_x=blockIdx.x*blockDim.x*tx_size+threadIdx.x*tx_size+x;
         int grid_y=blockIdx.y*blockDim.y*ty_size+threadIdx.y*ty_size+y;
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
//         assert(0 <= ti1 && ti1 < ry_size && 0 <= ti2 && ti2 < height+1 && "Rectangular coordinate for rows out of range");
//         assert(0 <= tj1 && tj1 < rx_size && 0 <= tj2 && tj2 < width+1 && "Rectangular coordinate for columns out of range");

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

//                     assert(0 <= i1 && i1 <= i2 && i2 < height && "Triangular coordinate for rows out of range");
//                     assert(0 <= j1 && j1 <= j2 && j2 < width && "Triangular coordinate for columns out of range");
 
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
                        //minheap_rise(r,heap_size,sizeof(rectangle),heap_rcmp);
                     // otherwise check whether the current score is greater than the smallest element (=first element) 
                     // if so, replace the first element with current rectangle and sink the current element until 
                     // the heap condition holds. 
                     } /*else if (heap_rcmp(&current,r) > 0) {
                        r[0] = current;
                        //minheap_sink(r,kbest,sizeof(rectangle),heap_rcmp);
                     } */
                  }
               }

               // store heap_size
               //blockIdx.y*blockDim.y*gridDim.x*blockDim.x+blockIdx.x*blockDim.x*blockDim.y+threadIdx.y*blockDim.x+threadIdx.x
               heapsize_global[((blockIdx.y*blockDim.y+threadIdx.y)*gridDim.x*blockDim.x+blockIdx.x*blockDim.x+threadIdx.x)] = heap_size;
               // Write heap back to global memory
               memcpy(&r_global[((blockIdx.y*blockDim.y+threadIdx.y)*gridDim.x*blockDim.x+blockIdx.x*blockDim.x+threadIdx.x)*kbest],r,heap_size*sizeof(rectangle));
}*/


// test gpu lrt code
/*int main(int argc, char* argv[]){
   // number of the best scoring rectangles printed by the output
   int kbest=DEFAULT_KBEST; 
   size_t size;         // number of rectangles in the grid
   //generate grid 
   grid      *g;           // input grid
   FILE *input = stdin;
    if (optind < argc){
        if (argc-optind > 1){
            error("only one input file permitted."); 
        }
        if ((input=fopen(argv[optind],"rb"))==NULL){
            error("cannot open input file %s.",optarg);
        }
    }
   g = read_grid(input);
   int width=g->width;
   int height=g->height;
   // compute prefix sums
   prefix_sums *p = compute_psums(g);
   // compute grid totals, ratio, and likelihood
   int g_n = p->A[width * height - 1].n;
   int g_k = p->A[width * height - 1].k;
   float g_q = (float)g_k / (float)g_n;
   float g_l = g_k * log(g_q);

   //block size ,means the number of theards in each block 
   int bx_size = 2;
   int by_size = 2;
   // define tile size for a thread
   int tx_size = 4; 
   int ty_size = 4; 
   // rectangular length
   int rx_size = (width + 1)/2;
   int ry_size = (height + 1)/2;
   int gx_size=((width+1)*rx_size + bx_size*tx_size-1)/(bx_size*tx_size);
   int gy_size=((height+1)*ry_size + by_size*ty_size-1)/(by_size*ty_size);
   
   // allocate a heap for each thread with size k 
   rectangle *r_global = (rectangle *)ALLOCV(rectangle,kbest*bx_size*by_size*gx_size*gy_size);
   int *heapsize_global = (int *)ALLOCV(int,bx_size*by_size*gx_size*gy_size);

   // A [(N+1],m]^2 grid needs to be mappeds to a two-dimensional
   // grid [(N+1)*m]^2 and divided by the block size. Ensure that if the
   // block size does not divide grid length, one extra block is added 
   // (NB: This is done by adding bx_size-1 before division). 

    size = (size_t)(g->width + 1) *
           (size_t)g->width  *
           (size_t)(g->height + 1) *
           (size_t)g->height / 4;
    // adjust kbest
    kbest = min(size, kbest);
    
   //dimg encodes the dimension of grid in terms of the number of
   //blocks in each dimension
   dim3 dimg(gx_size, gy_size);
   //dimb encodes the dimension of each block in terms of the number of threads per block 
   dim3 dimb(bx_size, by_size);
   gpu_lrt<<<dimg, dimb>>>(p, r_global,heapsize_global,tx_size,ty_size,width, height, kbest,g_n,g_k,g_q,g_l);
   
}*/
