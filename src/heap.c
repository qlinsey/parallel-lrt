//
//  heap.c
//  LRT
//
//  Copyright (c) 2013 The University of Sydney. All rights reserved.
//
//  command line tool for computing likilhood ratio tests of 2d grids
//  

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "heap.h"

////////////////////////////////////////////////////////////////////////////
//  functions
//

//
// sinks the first element of the min-heap until heap condition holds
//
void minheap_sink(void *heap, size_t n, size_t size, int (*cmp) (const void *,const  void * ))
{
   // allocate memory for swap on heap
   unsigned char tmp[size];
   int k;

   for (k=0;2*k<n;) {

      // set j to first child of k
      int j = 2 * k + 1;

      // check which of the children (if there are two) is smaller
      // move j to the smaller child so that the smaller child is 
      // swapped with a[k]
      if (j < n-1) {
         if ((*cmp)(heap + (j+1)*size, heap + j*size) < 0) {
            j++;
         }
      }

      // check whether min-mheap condition holds, i.e a[k] is smaller 
      // than smallest child. If so, stop.
      if ((*cmp)(heap + k*size, heap + j*size) < 0)
         break;

      // swap elements k and j 
      memcpy(tmp,heap + k*size,size);
      memcpy(heap+k*size,heap + j * size,size);
      memcpy(heap+j*size,tmp,size);

      // set new n
      k = j;
   }
}

//
// raises the last element of the min-heap until heap condition holds
//
void minheap_rise(void *heap, size_t n, size_t size, int (*cmp) (const void *, const void * ))
{
   // allocate memory for swap on heap
   unsigned char tmp[size];
   int k;

   for (k=n-1;k > 0; ) { 
       // set j to parent
       int j = k/2;

       // check whether a[k] is smaller than parent a[j] 
       if ((*cmp)(heap + k*size,heap + j*size) < 0) { 
          // swap elements k and j 
          memcpy(tmp,heap + k*size,size);
          memcpy(heap+k*size,heap + j * size,size);
          memcpy(heap+j*size,tmp,size);

          // set k to current node 
          k = j; 
       } else  break;
   } 
}

//
// simple unit test for heap functions
//
// int int_cmp(void *x, void *y)
// {
//    return *((int *)x) - *((int *)y);
// }
// 
// int main() {
//    int x[]={10,5,4,11,12,8,9};
//    int y[]={1,5,4,11,12,8,9,0};
//    int i;
// 
//    minheap_sink(x,sizeof(x)/sizeof(x[0]),sizeof(x[0]),int_cmp);
// 
//    for(i=0;i<sizeof(x)/sizeof(x[0]);i++){ 
//      printf("%d. %d\n",i,x[i]);
//    }
// 
//    minheap_rise(y,sizeof(y)/sizeof(y[0]),sizeof(y[0]),int_cmp);
//    for(i=0;i<sizeof(y)/sizeof(y[0]);i++){ 
//      printf("%d. %d\n",i,y[i]);
//    }
// }
