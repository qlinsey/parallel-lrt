//
//  multicore_lrt.c
//  LRT
//
//  Copyright (c) 2013 The University of Sydney. All rights reserved.
//
//  naive lrt computation with high complexity and high memory usage 
//

#include <pthread.h>
#include <time.h>
#include <sys/time.h>
#include "lrt.h"
#include "heap.h"
#include "prefix.h"
////////////////////////////////////////////////////////////////////////////
//  type definition
//

// parameter for a single thread
typedef struct thread_param {
   long long start,              // rectangle range
             end; 
   prefix_sums *p;               // prefix sums
   int kbest;                    // kbest
} thread_param;
   

////////////////////////////////////////////////////////////////////////////
//  global variables
//

extern int num_threads;
extern int profiling_flag;

////////////////////////////////////////////////////////////////////////////
//  functions
//


//
// work function for a single thread
//
// a rectangle is encoded as a single 
// number. First the rectanlge coordinates
// need to be retrieved from the number
// and then the LRT is to be computed. 
// 
static void *work_function(void *arg)
{
   // get info from parameters of the thread
   thread_param *param = (thread_param *) arg;
   prefix_sums *p = param->p;
   grid *g = p->gridref; 
   int kbest=param->kbest;
   
   // current size of the heap that stores the kbest
   int heap_size=0;

   // set size of grid, and rectangular compressed space
   int width=g->width;
   int height=g->height;
   long long m_x=(width+1)/2;
   long long m_y=(height+1)/2;

   // compute grid totals, ratio, and likelihood
   int   g_n=p->A[width * height -1].n ;
   int   g_k=p->A[width * height -1].k ;
   float g_q=(float)g_k/(float)g_n;
   float g_l=g_k*log(g_q);
  
   // heap to store kbest
   rectangle *r = (rectangle *) ALLOCV(rectangle, kbest);
   
   long long idx;
   for(idx=param->start;idx<param->end;idx++){

      // convert index to rectangular coordinates
      long long val=idx; 
      int tj1 = val % m_x;
      val /= m_x;
      int tj2 = val % (width + 1);
      val /= width+1;
      int ti1= val % m_y;
      val /= m_y;
      int ti2 = val;
      assert(0 <= ti1 && ti1 < m_y && 0 <= ti2 && ti2 < height+1 && "rectangular coordinate for row out of range");
      assert(0 <= tj1 && tj1 < m_x && 0 <= tj2 && tj2 < width+1 && "rectangular coordinate for column out of range");

      // 
      // convert rectangular coordinates to triangular coordinates.
      int j1,j2,i1,i2;
      if (tj2 < width- tj1) {
           j1=tj1;
           j2=tj2+tj1;
      } else if( (tj1+1)*2 < width+1){
           j1=width-tj1-1;
           j2=width-tj2+j1;
      } else continue;
      if (ti2 < height - ti1) {
           i1=ti1;
           i2=ti2+ti1;
      } else if( (ti1+1)*2 < height+1){
           i1=height-ti1-1;
           i2=height-ti2+i1;
      } else continue;
      assert(0 <= i1 && i1 <= i2 && i2 < height && "triangular coordinate for row out of range");
      assert(0 <= j1 && j1 <= j2 && j2 < width && "triangular coordinate for column out of range");

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

   return r;
}

//
// parallel multicore implementation of LRT using pthreads
// 
// the complexity of the naive implementation is O(n^6) and 
// hence not practical for larger grids. 
//
rectangle *multicore_lrt(grid *g,int kbest)
{
   // compute number of rectangles and amount of work
   // a thread needs to do
   pthread_t thread[num_threads];
   thread_param param[num_threads];
   int i;

   assert(g != NULL && "no grid.");
   assert(kbest > 0 && "kbest is zero or smaller");
   assert(num_threads > 0 && "number of threads is zero or smaller");
   //time measurement using wallclock time
   struct timespec begin, end;
   double elapsed;
   clock_gettime(CLOCK_MONOTONIC, &begin);
   
   // allocate memory for result  
   rectangle *r = (rectangle *) ALLOCV(rectangle, num_threads * kbest);

   // compute prefix sums 
   prefix_sums *p = compute_psums(g);

   // compute number of rectangles and amount of work 
   // a thread needs to do 
   long long start = 0; 
   long long size = (long long)(g->height + 1) * 
                    (long long)((g->height + 1) / 2) *
                    (long long)((g->width + 1) *
                    (long long)((g->width + 1) / 2));

   long long stride = size / num_threads;
   // initialize parameter structure and spawn thread
   for (i=0;i<num_threads;i++){

       // set parameters  for threads
       param[i].start = start;
       if (i == num_threads-1) {
          param[i].end = size;
       } else {
          param[i].end = start + stride;
       }
       param[i].p = p;
       param[i].kbest = kbest; 

       // create thread 
       pthread_create(&thread[i], NULL, &work_function, &param[i]);

       // update start for next thread
       start = start + stride;
   }
  
   // wait for termination and copy result
   for (i=0;i<num_threads;i++){
       // rectangles of current thread
       rectangle *thread_r;
       // wait for termination of thread[i]
       pthread_join(thread[i], (void **)&thread_r);

       // copy result 
       memcpy(r+i*kbest,thread_r,sizeof(rectangle)*kbest);  
   }

   // sort all result and return it
   qsort(r,num_threads*kbest,sizeof(rectangle),rectangle_compare);
   clock_gettime(CLOCK_MONOTONIC, &end);
   elapsed = end.tv_sec - begin.tv_sec;
   elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;
   if(profiling_flag)
   fprintf(stderr,"time to compute LRT on multicore: %f (s)\n",elapsed);
   return r;
}
