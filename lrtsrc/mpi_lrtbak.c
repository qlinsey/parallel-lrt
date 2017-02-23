//
//  mpi_lrt.c
//  LRT
//
//  Copyright (c) 2013 The University of Sydney. All rights reserved.
//
//  MPI implementation of LRT
//

#include <mpi.h>
#include "lrt.h"
#include "heap.h"
#include "prefix.h"

#define send_data_tag 2013
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

// profiling flag: if switch on
Bool profiling_flag=FALSE;
//
// sinks the first element of the min-heap until heap condition holds
//
 __inline__ void mpiheapsink(rectangle *heap, int n)
{
   // allocate memory for swap on heap
   rectangle tmp;
   int k;
   for (k=0;2*k<n;) {

      // set j to first child of k
      int j = 2 * k + 1;

      // check which of the children (if there are two) is smaller
      // move j to the smaller child so that the smaller child is 
      // swapped with a[k]
      if (j < n-1) {
         if (heap_rcmp(&heap [j+1], &heap[j]) < 0) {
            j++;
         }
      }
      // check whether min-mheap condition holds, i.e a[k] is smaller 
      // than smallest child. If so, stop.
      if (heap_rcmp(&heap[k], &heap[j]) <0)
         break;
      // swap elements k and j 
       tmp=heap[k];
       heap[k]=heap[j];
       heap[j]=tmp;
       // set new n
       k = j;
   }
}

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
         mpiheapsink(inout_heap ,*len);
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
   init_rectangles(r,kbest);
   
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
      //since mpi process do initialization of kbest rectangles, so just do heap sink
      //heap rise couldn't rise up the -FLT_MAX value
      if (heap_rcmp(&current,r) > 0) {
         r[0] = current;
         mpiheapsink(r ,kbest); 
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
    while ((c = getopt (argc, argv, "po:k:h")) != -1) {
        switch (c) {
            case 'k':
                kbest = atoi(optarg); 
                if (kbest <= 0 ) {
                   error("k-best parameter must greater than zero.");
                }
                break;

            case 'p':
                profiling_flag=TRUE;
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
   clock_t   t_read,       // time for reading LRT
             t_compute,    //time for computing LRT
             t_send,       //time for measuring the send data copy overhead
             t_reduce;     // time for measuring the reduce overhead
   double    p_start,      //measure each process compute time 
             p_end,
             p_last;
   double    p_cstart,
             p_cend,
             p_clast;
   struct    timespec begin, end;
   double    elapsed;
   int myid,               // MPI process id
       numprocs,           // MPI number of processors
       sender;             // sender id, MPI resource 
   int ret;

   MPI_Datatype r_type;     // MPI datatype for a cell
   MPI_Op op;              // MPI operation for reduction
   MPI_Status status;      // MPI_STATUS 
   // initialize MPI
   MPI_Init(&argc,&argv); 

   // declare MPI types 
   rectangle rs; 
   MPI_Datatype r_elements[2]={MPI_FLOAT, MPI_SHORT}; 
   MPI_Aint r_disp[2]; 
   int r_blocklen[2]={1,4}; 
   r_disp[0] = (void *)&rs.score - (void *)&rs; 
   r_disp[1] = (void *)&rs.i1 - (void *)&rs.score;
   MPI_Type_create_struct(2,r_blocklen, r_disp, r_elements, &r_type);
   MPI_Type_commit(&r_type); 

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
      t_read = clock();                 // read the grid (and measure time) 
      g = read_grid(input);
      t_read = clock() - t_read; 
   } else {
      g = ALLOC(grid); 
   } 
      
 
   // start barrier for measuring broadcast overhead
   MPI_Barrier(MPI_COMM_WORLD);
   if (myid == 0) {
      t_send=clock();                   // measure for sending overhead
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

   // end barrier for measuring broadcast overhead
   MPI_Barrier(MPI_COMM_WORLD);
   // start compute time 
   p_start = MPI_Wtime(); 
   
   if (myid == 0) {
      t_send=clock()-t_send;
      t_compute = clock();  
      p_cstart=MPI_Wtime();
      clock_gettime(CLOCK_MONOTONIC, &begin);        
   }

   // compute number of rectangles and amount of work 
   // a thread needs to do 
   long long size = (long long)(g->height + 1) * 
                    (long long)((g->height + 1) / 2) *
                    (long long)((g->width + 1) *
                    (long long)((g->width + 1) / 2));
                    
   size_t    org_size = (size_t)(g->width + 1) *
                        (size_t)g->width  *
                        (size_t)(g->height + 1) *
                        (size_t)g->height / 4;

   // adjust kbest
   kbest = min(org_size, kbest);
   long long stride = (size + numprocs - 1 ) / numprocs;
   long long start = myid*stride;
   long long end = (myid+1)*stride;
   if (end > size) { 
     end = size;
   }

   // do computation 
   rectangle *r = compute(g,kbest,start,end); 

   // start barrier for measuring reduce overhead
   MPI_Barrier(MPI_COMM_WORLD);
   //end of compute time 
   p_end = MPI_Wtime();
   p_last=p_end-p_start;
   if(myid!=0) 
      MPI_Send(&p_last, 1, MPI_DOUBLE, 0, send_data_tag, MPI_COMM_WORLD);
   if (myid == 0) {
      t_compute = clock() - t_compute;
      p_cend=MPI_Wtime();
      p_clast=p_cend-p_cstart;
      t_reduce=clock();                 // measure the reduce overhead
   }

   // allocate answer for result 
   rectangle *result = ALLOCV(rectangle, kbest); 
   init_rectangles(result,kbest);

   ret=MPI_Reduce(r, result, kbest, r_type, op, 0, MPI_COMM_WORLD); 
   assert(ret == MPI_SUCCESS && "MPI: reduction failed");
      
   // end barrier for measuring reduce overhead
   MPI_Barrier(MPI_COMM_WORLD);
   // root node writes results 
   if (myid == 0)  {
      t_reduce=clock()-t_reduce;      // measure the reduce overhead

      // sort result and return it
      qsort(result,kbest,sizeof(rectangle),rectangle_compare);
      if (output != NULL) {
         print_rectangles(output, result, kbest);
      }

      if (profiling_flag) {             // output time measurements and memory consumption
         clock_gettime(CLOCK_MONOTONIC, &end);
         elapsed = end.tv_sec - begin.tv_sec;
         elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;
         assert(output != NULL && "no file descriptor");
         fprintf(stderr,"read-input-time-s,%f\n",((double)t_read)/CLOCKS_PER_SEC);
         fprintf(stderr,"total-compute-time-s,%f\n",((double)t_compute)/CLOCKS_PER_SEC);
         fprintf(stderr,"memory-usage-kb,%f\n",((double)get_memusage()) / 1024);
         fprintf(output,"time-to-read-input-file-s,%f\n",((double)t_read)/CLOCKS_PER_SEC);
         fprintf(output,"time-to-broadcast-grid-s,%f\n",((double)t_send)/CLOCKS_PER_SEC);
         fprintf(output,"time-to-reduce-kbest-s,%f\n",((double)t_reduce)/CLOCKS_PER_SEC);
         fprintf(output,"time-to-compute-LRT-s,%f\n",((double)t_compute)/CLOCKS_PER_SEC);
         fprintf(output,"process-root-compute-s,%f\n",p_clast);
         fprintf(output,"process-elapsed-root-compute-s,%f\n",elapsed);
         int tp;
         for(tp=1;tp<numprocs;tp++){
            MPI_Recv(&p_last, 1, MPI_DOUBLE, MPI_ANY_SOURCE, send_data_tag, MPI_COMM_WORLD,&status);
            sender = status.MPI_SOURCE;
            fprintf(output,"process-%i-compute-s,%f\n",sender,p_last);
        }
      }
   }
   MPI_Finalize();
    
   return 0;
}
