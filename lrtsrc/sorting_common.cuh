//
//  sorting_common.cuh
//  LRT
//
//  Copyright (c) 2013 The University of Sydney. All rights reserved.
//
//  define comparators 
//  

#include "lrt.h"

#ifndef SORTING_COMMON_CUH
#define SORTING_COMMON_CUH

__device__ inline int Comparator(rectangle r_a, rectangle r_b){


   if(isnan(r_a.score)) return 1;
   if(isnan(r_b.score)) return -1;
   
   float score_a = r_a.score;
   float score_b = r_b.score;

   if (score_a < score_b) return -1;

   if (score_a > score_b) return 1;



   int   area_a=(r_a.j2-r_a.j1+1)*(r_a.i2-r_a.i1+1);
   int   area_b=(r_b.j2-r_b.j1+1)*(r_b.i2-r_b.i1+1);

   if(area_a < area_b) return 1;
   if(area_a > area_b) return -1;
 
   if(r_a.i1 < r_b.i1) return 1;
   if(r_a.i1 > r_b.i1) return -1;

   if(r_a.i2 < r_b.i2) return 1;
   if(r_a.i2 > r_b.i2) return -1;

   if(r_a.j1 < r_b.j1) return 1;
   if(r_a.j1 > r_b.j1) return -1;

   if(r_a.j2 < r_b.j2) return 1;
   if(r_a.j2 > r_b.j2) return -1;

   return 0;

}
#endif