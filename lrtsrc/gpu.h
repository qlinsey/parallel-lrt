//
//  gpu.h
//  LRT
//
//  Copyright (c) 2013 The University of Sydney. All rights reserved.
//

#ifndef __GPU_H__
#define __GPU_H__

// threads in a block 
#define BLOCKSIZE_X (16) 
#define BLOCKSIZE_Y (8) 

// blocks in a grid
#define GRIDSIZE_X (8) 
#define GRIDSIZE_Y (32) 

// Total number of threds in a grid
#define THREADS_IN_GRID (BLOCKSIZE_X * BLOCKSIZE_Y * GRIDSIZE_X * GRIDSIZE_Y) 

#endif

