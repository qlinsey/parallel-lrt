//
//  lrt.h
//  LRT
//
//  Copyright (c) 2013 The University of Sydney. All rights reserved.
//

#ifndef __LRT_H__
#define __LRT_H__

////////////////////////////////////////////////////////////////////////////
// includes
//

#include <stdio.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <getopt.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>

// disable assert statements to speedup 
// program execution if definition FAST 
// is set. 
#ifndef FAST
#include <assert.h>
#else
#define assert(a)
#endif

////////////////////////////////////////////////////////////////////////////
//  defines
//

//
// constants  
#define VERSION        ("0.1")
#define MAGIC_NUMBER   (27182818)
#define DEFAULT_KBEST  (10)
#define DEFAULT_METHOD (0) 
#define DEFAULT_WIDTH  (10)
#define DEFAULT_HEIGHT (10)
#define DEFAULT_MAX_POPULATION (10)
#define DEFAULT_LAMBDA (5)

//
//define the outlier region 
#define OUTLIER_LAMBDA (70)
#define DEFAULT_OUTLIER_ROWS (1)
#define DEFAULT_OUTLIER_COLS (1)
#define DEFAULT_OUTLIER_RIDX (1)
#define DEFAULT_OUTLIER_CIDX (1)

//
// memory allocation macros
#define ALLOC(t)       (t *)(mymalloc(sizeof(t)))
#define ALLOCV(t,num)  (t *)(mymalloc(sizeof(t)*num)) 

//
// debug macros allows several debug levels
#ifndef DEBUG_LEVEL 
#define DEBUG(l,fmt,...)
#else 
#define DEBUG(l,fmt,...) { if(l <= DEBUG_LEVEL ) { fprintf(stderr,"DEBUG:");fprintf(stderr,fmt,##__VA_ARGS__);fprintf(stderr,"\n");}}
#endif

//
// min & max macros
#define min(a,b) (((a)<(b))?(a):(b))

// 
// index for two-dimensional array assuming width is defined
#define I(row,col) ((row)*width+(col))

////////////////////////////////////////////////////////////////////////////
//  type definitions
//

//
// bolean data type
typedef enum {FALSE=0, TRUE=1} Bool; 

//
// define a cell
typedef struct cell {
   int n; // total population count in a cell
   int k; // total number of cases reported in a cell
} cell;

//
// defines a grid 
typedef struct grid {
 int width,   // width and height of grid 
     height; 
 cell *cells; // cells of the grid
} grid;

//
// define a rectangle with its score
typedef struct rectangle {
 float score; // score of rectangle 
              // coordinates of rectangle 
 short i1,      //  row of upper-left corner
     j1,      //  column of upper-left corner
     i2,      //  row of lower-right corner
     j2;      //  column of lower-right corner
} rectangle;

////////////////////////////////////////////////////////////////////////////
// function prototypes
//

void *mymalloc(size_t);                                // allocate memory
void error(char *, ...);                               // error function 
grid* read_grid(FILE *);                               // read grid 
void write_grid(FILE *, grid *);                       // write grid to a file
void print_rectangles(FILE *, rectangle *, int);       // print rectangles 
size_t get_memusage(void);                             // query memory usage
rectangle *naive_lrt(grid *,int);                      // naive lrt computation
rectangle *fast_lrt(grid *,int);                       // fast lrt computation
rectangle *multicore_lrt(grid *,int);                  // multicore lrt computation
rectangle *gpu_lrt(grid *,int);                        // gpu simulation lrt computation
rectangle *gpusim_lrt(grid *,int);                     // gpu lrt computation
int rectangle_compare(const void *, const void *);     // compare two rectangles for qsort
int heap_rcmp(const void * a, const void * b);         // compare two rectangles for heapsort
void free_grid(grid *);                                // free grid

#endif

