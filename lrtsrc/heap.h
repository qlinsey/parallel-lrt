//
//  heap.h
//  LRT
//
//  Copyright (c) 2013 The University of Sydney. All rights reserved.
//

#ifndef __HEAP_H__
#define __HEAP_H__

//
// rises an element in the heap
//
void minheap_rise(void *ptr, size_t num, size_t size, int (*) (const void *, const void * )); 

//
// sinks an element in the heap
//
void minheap_sink(void *ptr, size_t num, size_t size, int (*) (const void *, const void * )); 

#endif

