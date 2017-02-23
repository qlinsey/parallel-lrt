//
//  io.c
//  LRT
//
//  Copyright (c) 2013 The University of Sydney. All rights reserved.
//
//  input/output functions for lrt computation
//

#include "lrt.h"

////////////////////////////////////////////////////////////////////////////
// simple grid file format 
// 
//   - header consisting of 3 integers 
//      * magic number to identifier grid file format  
//      * width of grid 
//      * heigh of grid 
//   - cells of the grid
//      * cell data structure directly written width*heigh times

////////////////////////////////////////////////////////////////////////////
//  functions
//

//
// read grid
// 
grid* read_grid(FILE *in)
{
    int   width, 
          height;
    int   size;
    grid *g;
    int header[3];
    cell *c;

    // read and check header 
    if (fread(header,sizeof(int),3,in)!=  3){
       error("cannot read grid file header");  
    }
    if (header[0] != MAGIC_NUMBER) { 
       error("wrong grid file format"); 
    } 
    width = header[1];
    height = header[2]; 
    size = width * height; 
    assert(width >0 && height > 0 && "wrong grid dimensions whilst reading grid file"); 

    DEBUG(1,"read grid size: %d x %d",width, height); 
   
    // read grid data  
    c = (cell *) ALLOCV(cell, size);
    if (fread(c,sizeof(cell),size,in)!=  size){
       error("cannot read grid file data");  
    }

    // allocate and initialize grid data structure
    g = (grid *) ALLOC(grid);
    g->width = width; 
    g->height = height; 
    g->cells = c;

    return g; 
}

//
// write grid to a file
// 
void write_grid(FILE *out, grid *g)
{
    int header[3];
    int size;

    // write header 
    header[0]=MAGIC_NUMBER;
    header[1]=g->width;
    header[2]=g->height; 
    size = g->width * g->height; 
    if (fwrite(header,sizeof(int),3,out) !=  3) { 
       error("cannot write grid file header"); 
    }

    // write grid data 
    if (fwrite(g->cells,sizeof(cell),size,out)!= size){
       error("cannot write grid file data");  
    }
}

//
// print rectangles 
// 
void print_rectangles(FILE *out, rectangle *r, int n)
{
    int i; 
    
    assert(out != NULL && "no file descriptor");
    assert(r != NULL && "no rectangles"); 
    
    for (i=0;i<n;i++) {
        fprintf(out,"(%d,%d)-(%d,%d) %f\n",r[i].i1,r[i].j1,r[i].i2,r[i].j2,r[i].score);    
    }
}
