//
//  misc.c
//  LRT
//
//  Copyright (c) 2013 The University of Sydney. All rights reserved.
//
//  command line tool for computing likilhood ratio tests of 2d grids
//  

#include "lrt.h"

////////////////////////////////////////////////////////////////////////////
//  global variables
//

//
// number of bytes used dynamically
static size_t mem_usage;

//
// program name of executable 
extern char *progname;

////////////////////////////////////////////////////////////////////////////
//  functions
//

//
// allocate memory
//
void *mymalloc(size_t n)
{
    void *ptr;

    // allocate and check memory    
    if ((ptr= (void *)malloc(n)) == NULL) {
        error("memory exhausted");
    }
    mem_usage += n;
    return ptr;
}

//
// get memory usage
//
size_t get_memusage(void)
{
   return mem_usage;
}

//
// error routine
// 
void error(char *fmt, ...)
{
    va_list va;
    va_start(va, fmt);
    fprintf(stderr,"%s:",progname); 
    vfprintf(stderr,fmt,va);
    fprintf(stderr,"\n"); 
    exit(1);
}

//
// comparison function for qsort that compares the scores of two rectangles.
//
int rectangle_compare(const void * a, const void * b)
{
   rectangle *r_a = (rectangle *)a; 
   rectangle *r_b = (rectangle *)b;

   float score_a = r_a->score;
   float score_b = r_b->score;

   if (score_a < score_b) return 1;

   if (score_a > score_b) return -1;

   if(isnan(score_a)) return 1;
   if(isnan(score_b)) return -1;

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
// compare two rectangles a and b: if the score of a is smaller than score of b, then 
// return a number smaller than zero. If they are the same, return zero. Otherwise
// return a number greater than zero. 
// 
int heap_rcmp(const void * a, const void * b)
{
   return  - rectangle_compare(a,b);
}


