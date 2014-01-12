//
//  prefix.c
//  LRT
//
//  Copyright (c) 2013 The University of Sydney. All rights reserved.
//
//  command line tool for computing likilhood ratio tests of 2d grids
//  

#include "prefix.h"

////////////////////////////////////////////////////////////////////////////
// Functions

////////////////////////////////////////////////////////////////////////////
// Recursions for prefix sums A,B,X, and Y 

//
// compute prefix sum A (upper left corner), i.e.,
//  A[i,j] = sum{1<=l<=i,1<=k<=j} G[l,k]
//
void compute_sumA(prefix_sums *p)
{
   int i,j,k; 
   
   assert(p != NULL && "invalid prefix sum pointer");
   assert(p->gridref != NULL && "no associated grid");
  
   cell *A = p->A;
   cell *G = p->gridref->cells;

   // set height and width of grid
   int width = p->gridref->width; 
   int height = p->gridref->height; 
   
   // set first element
   A[0].n = G[0].n;
   A[0].k = G[0].k;
   
   // first column
   for (i=1;i<height;i++) {
      A[I(i,0)].n = G[I(i,0)].n + A[I(i-1,0)].n; 
      A[I(i,0)].k = G[I(i,0)].k + A[I(i-1,0)].k;
   }

   // first row 
   for (j=1;j<width;j++) {
      A[I(0,j)].n = G[I(0,j)].n + A[I(0,j-1)].n; 
      A[I(0,j)].k = G[I(0,j)].k + A[I(0,j-1)].k;
   }

   // iterate over all diagonal elements
   for (k=1;k<min(width,height);k++) {
      // associate column of diagonal element k
      for (i=k;i<height;i++) {
         A[I(i,k)].n = G[I(i,k)].n + A[I(i-1,k)].n + A[I(i,k-1)].n - A[I(i-1,k-1)].n;
         A[I(i,k)].k = G[I(i,k)].k + A[I(i-1,k)].k + A[I(i,k-1)].k - A[I(i-1,k-1)].k;
      }
      // associate rows of diagonal element k
      for(j=k;j<width;j++) {
         A[I(k,j)].n = G[I(k,j)].n + A[I(k,j-1)].n + A[I(k-1,j)].n - A[I(k-1,j-1)].n;
         A[I(k,j)].k = G[I(k,j)].k + A[I(k,j-1)].k + A[I(k-1,j)].k - A[I(k-1,j-1)].k;
      }
   }
}

//
// compute prefix sum B (lower right corner), i.e.,
//  B[i,j] = sum{i<=l<=n,j<=k<=m} G[l,k]
//
void compute_sumB(prefix_sums *p)
{
   int i,j,k; 

   assert(p != NULL && "invalid prefix sum pointer");
   assert(p->gridref != NULL && "no associated grid");
   
   cell *B = p->B;
   cell *G = p->gridref->cells;

   // set height and width of grid
   int width = p->gridref->width; 
   int height = p->gridref->height; 

   // set last element
   B[I(height-1,width-1)].n =G[I(height-1,width-1)].n;
   B[I(height-1,width-1)].k =G[I(height-1,width-1)].k;

   // last column
   for(i=height-1;i>=1;i--) {
      B[I(i-1,width-1)].n =G[I(i-1,width-1)].n + B[I(i,width-1)].n;
      B[I(i-1,width-1)].k =G[I(i-1,width-1)].k + B[I(i,width-1)].k;
   }

   // last row
   for(j=width-1;j>=1;j--){
      B[I(height-1,j-1)].n =G[I(height-1,j-1)].n + B[I(height-1,j)].n;
      B[I(height-1,j-1)].k =G[I(height-1,j-1)].k + B[I(height-1,j)].k;
   }

   // iterate over all diagonal elements
   for(k=1;k<min(width,height);k++){
      // associate column of diagonal element k
      for(i=height-1;i>=k;i--){
         B[I(i-k,width-1-k)].n = G[I(i-k,width-1-k)].n + B[I(i-k+1,width-1-k)].n+B[I(i-k,width-k)].n-B[I(i-k+1,width-k)].n;
         B[I(i-k,width-1-k)].k = G[I(i-k,width-1-k)].k + B[I(i-k+1,width-1-k)].k+B[I(i-k,width-k)].k-B[I(i-k+1,width-k)].k;
      }
      // associate row of diagonal element k
      for(j=width-1;j>=k;j--) {
         B[I(height-1-k,j-k)].n = G[I(height-1-k,j-k)].n + B[I(height-1-k,j-k+1)].n+B[I(height-k,j-k)].n-B[I(height-k,j-k+1)].n;
         B[I(height-1-k,j-k)].k = G[I(height-1-k,j-k)].k + B[I(height-1-k,j-k+1)].k+B[I(height-k,j-k)].k-B[I(height-k,j-k+1)].k;
      }
   }
}

//
// compute prefix sum X (lower left corner ), i.e.,
// X[i,j] = sum {1<=l<=i-1, j+1 <= k <= m } G[l,k]
//
void compute_sumX(prefix_sums *p)
{
   int i,j,k; 
   assert(p != NULL && "invalid prefix sum pointer");
   assert(p->gridref != NULL && "no associated grid");
   
   cell *X = p->X;
   cell *G = p->gridref->cells;

   // set height and width of grid
   int width = p->gridref->width; 
   int height = p->gridref->height; 

   // first column
   for(i=0;i<height;i++) {
      X[I(i,0)].n = 0;
      X[I(i,0)].k = 0;
   }

   // last row
   for(j=0;j<width;j++) {
      X[I(height-1,j)].n = 0;
      X[I(height-1,j)].k = 0;
   }

   // iterate over all diagonal elements
   for(k=1;k<min(width,height);k++) {
      // associate column of diagonal element k
      for(i=height-1;i>=k;i--) {
         X[I(i-k,k)].n = G[I(i-k+1,k-1)].n + X[I(i-k+1,k)].n+X[I(i-k,k-1)].n-X[I(i-k+1,k-1)].n;
         X[I(i-k,k)].k= G[I(i-k+1,k-1)].k + X[I(i-k+1,k)].k+X[I(i-k,k-1)].k-X[I(i-k+1,k-1)].k;
      }
      // associate row of diagonal element k
      for(j=k;j<width;j++) {
         X[I(height-1-k,j)].n = G[I(height-k,j-1)].n + X[I(height-1-k,j-1)].n+X[I(height-k,j)].n-X[I(height-k,j-1)].n;
         X[I(height-1-k,j)].k = G[I(height-k,j-1)].k + X[I(height-1-k,j-1)].k+X[I(height-k,j)].k-X[I(height-k,j-1)].k;
      }
   }
}

//
// compute prefix sum Y (upper right corner ), i.e.,
//  Y[i,j] = sum{i<=l<=n,j<=k<=m} G[l,k]
//
void compute_sumY(prefix_sums *p)
{
   int i,j,k; 
   assert(p != NULL && "invalid prefix sum pointer");
   assert(p->gridref != NULL && "no associated grid");
   
   cell *Y = p->Y;
   cell *G = p->gridref->cells;

   // set height and width of grid
   int width = p->gridref->width; 
   int height = p->gridref->height; 

   // last column
   for(i=1;i<height;i++) {
      Y[I(i,width-1)].n =0;
      Y[I(i,width-1)].k =0;
   }

   // first row
   for(j=0;j<width;j++) {
      Y[j].n =0;
      Y[j].k =0;
   }

   // iterate over all diagonal elements
   for(k=1;k<min(width,height);k++) {
      // associate column of diagonal element k
      for(i=k;i<height;i++) {
         Y[I(i,width-1-k)].n = G[I(i-1,width-k)].n + Y[I(i-1,width-1-k)].n + Y[I(i,width-k)].n-Y[I(i-1,width-k)].n;
         Y[I(i,width-1-k)].k = G[I(i-1,width-k)].k + Y[I(i-1,width-1-k)].k + Y[I(i,width-k)].k-Y[I(i-1,width-k)].k;
      }
      // associate row of diagonal element k
      for(j=width-1;j>=k;j--) {
         Y[I(k,j-k)].n = G[I(k-1,j-k+1)].n + Y[I(k,j-k+1)].n+Y[I(k-1,j-k)].n-Y[I(k-1,j-k+1)].n;
         Y[I(k,j-k)].k = G[I(k-1,j-k+1)].k + Y[I(k,j-k+1)].k+Y[I(k-1,j-k)].k-Y[I(k-1,j-k+1)].k;
      }
   }
}



////////////////////////////////////////////////////////////////////////////
// Unit tests for prefix sums A,B,X, and Y 

//
// unit test for prefix sum A (upper left corner)
// brute-force prefix sum A vs. compute_sumA()
//
Bool check_sumA(prefix_sums *p)
{
   cell *A = p->A;
   cell *G = p->gridref->cells;

   // set height and width of grid
   int width = p->gridref->width;
   int height = p->gridref->height;

   // check for all indices of the prefix sum
   int i, j;
   for (j=0;j<height;j++ ){
      for (i=0;i<width;i++ ){
         int l,k;
         // compute prefix sum for (i,j)
         int a_n=0;
         int a_k=0;
         for(l=0;l<=i;l++){
            for(k=0;k<=j;k++){
               a_n+=G[I(k,l)].n;
               a_k+=G[I(k,l)].k;
            }
         }
         DEBUG(2,"A: (%d,%d)=(%d,%d),A=(%d,%d)\n",i,j,a_n,a_k,A[I(j,i)].n,A[I(j,i)].k);
         // check equivalence
         if (A[I(j,i)].n != a_n || A[I(j,i)].k != a_k)  {
            return FALSE;
         }
      }
   }
   return TRUE;
}

//
// unit test for prefix sum B (lower right corner)
// brute-force prefix sum B vs. compute_sumB()
//
Bool check_sumB(prefix_sums *p)
{
   cell *B = p->B;
   cell *G = p->gridref->cells;
   int width = p->gridref->width;
   int height = p->gridref->height;
   int b_n=0;
   int b_k=0;
   Bool flag=TRUE;
   int i,j,l,k;
   // check for all indices of the prefix sum
   for(i=width-1;i>=0;i-- ){
      for( j=height-1;j>=0;j-- ){
         b_n=0;
         b_k=0;
         // compute prefix sum for (i,j)
         for(l=i;l<width;l++){
            for(k=j;k<height;k++){
               b_n+=G[I(k,l)].n;
               b_k+=G[I(k,l)].k;
            }
         }
         DEBUG(2,"b_(%d,%d)=(%d,%d),B=(%d,%d)\n",i,j,b_n,b_k,B[I(j,i)].n,B[I(j,i)].k);
         // check equivalence
         if((b_n!=B[I(j,i)].n)||(b_k!=B[I(j,i)].k)){
            flag=FALSE;
         }
      }
   }

   return flag;
}
//
// unit test for prefix sum X (lower left corner)
// brute-force prefix sum X vs. compute_sumX()
//
Bool check_sumX(prefix_sums *p)
{
   cell *X = p->X;
   cell *G = p->gridref->cells;
   int width = p->gridref->width;
   int height = p->gridref->height;
   int x_n=0;
   int x_k=0;
   Bool flag=TRUE;
   int i,j,l,k;

   // check for all indices of the prefix sum except last row, first column
   for(i=0;i<width-1;i++){
      for(j=height-1;j>0;j--){
         x_n=0;
         x_k=0;
         // compute prefix sum for (i+1,j-1)
         for(l=j;l<height;l++){
            for(k=0;k<=i;k++){
               x_n+=G[I(l,k)].n;
               x_k+=G[I(l,k)].k;
            }
         }
         DEBUG(2,"x_(%d,%d)=(%d,%d),X_(%d,%d)=(%d,%d)\n",i,j,x_n,x_k,i,j,X[I(j-1,i+1)].n,X[I(j-1,i+1)].k);
         // check equivalence
         if((x_n!=X[I(j-1,i+1)].n)||(x_k!=X[I(j-1,i+1)].k)){
            flag=FALSE;
         }
      }
   }
   return flag;
}

//
// unit test for prefix sum Y (upper right corner)
// brute-force prefix sum vs. compute_sumY()
//
Bool check_sumY(prefix_sums *p)
{

   cell *Y = p->Y;
   cell *G = p->gridref->cells;
   int width = p->gridref->width;
   int height = p->gridref->height;
   int y_n=0;
   int y_k=0;
   Bool flag=TRUE;
   int i,j,l,k;
   // check for all indices of the prefix sum except first row, last column
   for(i=width-1;i>0;i--){
      for(j=0;j<height-1;j++){
         y_n=0;
         y_k=0;
         // compute prefix sum for (i-1,j+1)
         for(k=width-1;k>=i;k--){
            for(l=0;l<=j;l++){
               y_n+=G[I(l,k)].n;
               y_k+=G[I(l,k)].k;
            }
         }
         DEBUG(2,"y_(%d,%d)=(%d,%d),Y_(%d,%d)=(%d,%d)\n",i,j,y_n,y_k,i,j,Y[I(j+1,i-1)].n,Y[I(j+1,i-1)].k);
         // check equivalence
         if((y_n!=Y[I(j+1,i-1)].n)||(y_k!=Y[I(j+1,i-1)].k)){
            flag=FALSE;
         }
      }
   }
   return flag;
}


////////////////////////////////////////////////////////////////////////////
// main function for prefix sums A,B,X, and Y 

//
// compute prefix sums of a grid
// 
prefix_sums *compute_psums(grid *g)
{
   // allocate structure 
   prefix_sums *p = (prefix_sums *) ALLOC(prefix_sums);
   p->gridref = g;
   
   // allocate memory for prefix sums
   int size = g->width * g->height; 
   p->A = (cell *) ALLOCV(cell,size);
   p->B = (cell *) ALLOCV(cell,size);
   p->X = (cell *) ALLOCV(cell,size);
   p->Y = (cell *) ALLOCV(cell,size);

   // fill in computation
   compute_sumA(p);
   compute_sumB(p);
   compute_sumX(p);
   compute_sumY(p);

   // check whether results are correct
   assert (check_sumA(p) && "computation of prefix A failed");
   assert (check_sumB(p) && "computation of prefix B failed");
   assert (check_sumX(p) && "computation of prefix X failed");
   assert (check_sumY(p) && "computation of prefix Y failed");

   return p;
}
