//
//  multicore_lrt.c
//  LRT
//
//  Copyright (c) 2013 The University of Sydney. All rights reserved.
//
//  naive lrt computation with high complexity and high memory usage 
//

#include <mpi.h>
#include "lrt.h"
#include "heap.h"
#include "prefix.h"

////////////////////////////////////////////////////////////////////////////
//  global variables
//

//
// program name (name of executable)
char *progname;

//
// file descriptor and name of input file. Default is stdin. 
FILE *input; 

//
// file descriptor and name of output file. Default is stdout. 
FILE *output; 

//
// number of the best scoring rectangles printed by the output
int kbest=DEFAULT_KBEST; 


////////////////////////////////////////////////////////////////////////////
//  functions
//


// 
// user defined reduce operation as MPI that joins two heaps
//
void heap_reduce(void *in, void *inout, int *len, MPI_Datatype *dptr)
{ 
   int i; 

   rectangle *in_heap = in;
   rectangle *inout_heap = inout; 

   // and rise the last element until the heap condition holds 
   for(i=0; i < *len; i++) { 
      if (heap_rcmp(&in_heap[i],inout_heap) > 0) {
         inout_heap[0] = in_heap[i];
         minheap_sink(&inout_heap[i],*len,sizeof(rectangle),heap_rcmp);
      } 
   }
}


//
// work function for a single thread
//
// a rectangle is encoded as a single 
// number. First the rectanlge coordinates
// need to be retrieved from the number
// and then the LRT is to be computed. 
// 
static rectangle *compute(grid *g, int kbest, long long start, long long end)
{
   // compute prefix sums 
   prefix_sums *p = compute_psums(g);
   
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
   for(idx=start;idx<end;idx++){

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
rectangle *mpi_lrt(grid *g,int kbest, int myid, int numprocs)
{
   // compute number of rectangles and amount of work
   // a thread needs to do
   int i;

   assert(g != NULL && "no grid.");
   assert(kbest > 0 && "kbest is zero or smaller");


   // compute number of rectangles and amount of work 
   // a thread needs to do 
   long long size = (long long)(g->height + 1) * 
                    (long long)((g->height + 1) / 2) *
                    (long long)((g->width + 1) *
                    (long long)((g->width + 1) / 2));
   // adjust kbest
   kbest = min(size, kbest);

   long long stride = (size + numprocs - 1 ) / numprocs;
   long long start = myid*stride;
   long long end = (myid+1)*stride;
   if (end > size) { 
     end = size;
   }

   // do the computation for rectangles in the range [start,end]
   return compute(g,kbest,start,end); 
}

//
// parse command line arguments
// 
void parse_arguments(int argc, char ** argv)
{
    int c;
    
    opterr = 0;
    
    input = stdin;
    output = NULL;
    while ((c = getopt (argc, argv, "o:k:h")) != -1) {
        switch (c) {
            case 'k':
                kbest = atoi(optarg); 
                if (kbest <= 0 ) {
                   error("k-best parameter must greater than zero.");
                }
                break;

            case 'o':
                if ((output=fopen(optarg,"wb"))==NULL){
                    error("cannot open output file %s.",optarg);
                }
                break;

            case 'h': 
                fprintf(stderr,"lrt: likely ratio test computation for 2-d grids\n Version %s\n",VERSION);
                fprintf(stderr,"Common Options:\n");
                fprintf(stderr,"  -o <filename> location of output file\n");
                fprintf(stderr,"  -k <#> kbest\n");                
                fprintf(stderr,"  -h help\n");                
                break;

            case '?':
                if (optopt == 'o')
                    error("Option -%c requires an argument.\n", optopt);
                else if (isprint (optopt))
                    error("Unknown option `-%c'.\n", optopt);
                else
                    error("Unknown option character `\\x%x'.\n",optopt);

            default:
                abort ();
        }
    }
    if (optind < argc){
        if (argc-optind > 1){
            error("only one input file permitted."); 
        }
        if ((input=fopen(argv[optind],"rb"))==NULL){
            error("cannot open input file %s.",optarg);
        }
    }
}

//
// declare MPI datatypes for struct cell and struct rectangle 
// 
void declare_types(MPI_Datatype *c_type, MPI_Datatype *r_type)
{
   // define MPI datastructure for cell
   cell c;
   MPI_Datatype c_elements[2]={MPI_INT, MPI_INT}; 
   MPI_Aint c_disp[2]; 
   int c_blocklen[2]={1,1}; 
   c_disp[0] = (void *)&c.n - (void *)&c; 
   c_disp[1] = (void *)&c.k - (void *)&c.n; 
   MPI_Type_create_struct(2,c_blocklen, c_disp, c_elements, c_type);
   MPI_Type_commit(c_type); 
  
   // define MPI datastructure for rectangle 
   rectangle r; 
   MPI_Datatype r_elements[5]={MPI_FLOAT, MPI_INT, MPI_INT, MPI_INT, MPI_INT}; 
   MPI_Aint r_disp[5]; 
   int r_blocklen[5]={1,1,1,1,1}; 
   r_disp[0] = (void *)&r.score - (void *)&r; 
   r_disp[1] = (void *)&r.i1 - (void *)&r.score;
   r_disp[2] = (void *)&r.j1 - (void *)&r.i1;
   r_disp[3] = (void *)&r.i2 - (void *)&r.j1;
   r_disp[4] = (void *)&r.j2 - (void *)&r.i2;
   MPI_Type_create_struct(5,r_blocklen, r_disp, r_elements, r_type);
   MPI_Type_commit(r_type); 
}

//
// main program
// 
// 1. parses input parameters
// 2. read an LRT problem
// 3. compute LRT 
// 4. print result
//
int main(int argc, char **argv)
{   
   grid      *g;           // input grid
   progname = argv[0];     // set name of executable

   int myid,               // MPI process id
       numprocs;           // MPI number of processors
   
   int ret;

   MPI_Datatype c_type;     // MPI datatype for a cell
   MPI_Datatype r_type;     // MPI datatype for a cell
   MPI_Op op;              // MPI operation for reduction

   // initialize MPI
   MPI_Init(&argc,&argv); 

   // declare MPI types 
   declare_types(&c_type,&r_type); 

   // define reduction operation 
   MPI_Op_create(heap_reduce, TRUE, &op); 

   // get size & rank
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);

   // if current node is the root node, parse args and read grid
   if (myid == 0) {     
      // parse arguments
      parse_arguments(argc,argv);      
      // only root process reads grid otherwise allocate data structure
      g = read_grid(input);
   } else {
      g = ALLOC(grid); 
   } 
       
   // send grid dimensions to all nodes  
   ret=MPI_Bcast(&g->height, 1, MPI_INT, 0, MPI_COMM_WORLD); 
   assert(ret == MPI_SUCCESS && "MPI: 1. BC failed");
   
   ret=MPI_Bcast(&g->width,  1, MPI_INT, 0, MPI_COMM_WORLD); 
   assert(ret == MPI_SUCCESS && "MPI: 2. BC failed");

   // if not root node allocate cells for grid
   if (myid != 0) { 
      g->cells = ALLOCV(cell, g->width*g->height); 
   } 

   // send grid data to all nodes
   ret=MPI_Bcast(g->cells, 2 * g->width * g->height, MPI_INT, 0, MPI_COMM_WORLD); 
   assert(ret == MPI_SUCCESS && "MPI: 3. BC failed");


   // compute number of rectangles and amount of work 
   // a thread needs to do 
   long long size = (long long)(g->height + 1) * 
                    (long long)((g->height + 1) / 2) *
                    (long long)((g->width + 1) *
                    (long long)((g->width + 1) / 2));
   // adjust kbest
   kbest = min(size, kbest);

   long long stride = (size + numprocs - 1 ) / numprocs;
   long long start = myid*stride;
   long long end = (myid+1)*stride;
   if (end > size) { 
     end = size;
   }

   // assume grid is large enough 
   assert(end-start > kbest && "grid to small"); 

   // do computation 
   rectangle *r = compute(g,kbest,start,end); 

   // allocate answer for result 
   rectangle *result = ALLOCV(rectangle, kbest); 
   memset(result, 0, sizeof(rectangle)*kbest); 

   // this might be an issue
   ret=MPI_Reduce(r, result, kbest, r_type, op, 0, MPI_COMM_WORLD); 
   assert(ret == MPI_SUCCESS && "MPI: reduction failed");

   // root node writes results 
   if (myid == 0)  {
      // sort result and return it
      qsort(result,kbest,sizeof(rectangle),rectangle_compare);
      if (output != NULL) {
         print_rectangles(output, result, kbest);
      }
   }

   MPI_Finalize();
    
   return 0;
}
