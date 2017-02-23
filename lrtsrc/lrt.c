//
//  lrt.c
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

//
// profiling flag: if switch on
Bool profiling_flag=FALSE;

//
// number of processors
int num_threads=1;

//
// array of function pointers pointing to lrt computation methods
rectangle* (*compute[])(grid *,int) = { 
   naive_lrt
   ,fast_lrt
   ,multicore_lrt
   ,gpusim_lrt
#ifdef CUDA
   ,gpu_lrt
#endif
};


//
// selector for lrt  method (range between 0 and # of 
// elements of array above)
int lrt_method=DEFAULT_METHOD;

////////////////////////////////////////////////////////////////////////////
//  functions
//

//
// parse command line arguments
// 
void parse_arguments(int argc, char ** argv)
{
    int c;
    
    opterr = 0;
    
    input = stdin;
    output = NULL;
    while ((c = getopt (argc, argv, "pk:m:t:ho:")) != -1) {
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

            case 'm':
                lrt_method = atoi(optarg); 
                if (lrt_method < 0 || lrt_method >= sizeof(compute) / sizeof(compute[0])) {
                   error("unknown lrt method. number ranges from 0 to %d",sizeof(compute)/sizeof(compute[0])-1);
                }
                break;

            case 't':
                num_threads = atoi(optarg); 
                if (num_threads <= 0 ) {
                   error("number of processors must greater than zero.");
                }
                break;

            case 'h': 
                fprintf(stderr,"lrt: likely ratio test computation for 2-d grids\n Version %s\n",VERSION);
                fprintf(stderr,"Common Options:\n");
                fprintf(stderr,"  -o <filename> location of output file\n");
                fprintf(stderr,"  -p enable profiling\n");
                fprintf(stderr,"  -m <#> sets lrt method\n");                
                fprintf(stderr,"     0 naive\n");                
                fprintf(stderr,"     1 fast\n");                
                fprintf(stderr,"     2 multi-core (parallelized fast lrt)\n");                
                fprintf(stderr,"     3 gpu simulator\n");                
                fprintf(stderr,"     4 gpu\n");                
                fprintf(stderr,"  -t <#> number of threads for parallel execution of lrt\n");                
                fprintf(stderr,"  -h help\n");                
                break;

            case 'p':
                profiling_flag = TRUE;
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
    rectangle *r;           // output rectangles 
    clock_t   t_read,       // time for reading LRT
              t_compute;    // time for computing LRT
    size_t size;         // number of rectangles in the grid
    progname = argv[0];               // set name of executable
    
    parse_arguments(argc,argv);       // parse arguments
    
    t_read = clock();                 // read the grid (and measure time) 
    g = read_grid(input);
    t_read = clock() - t_read; 
    size = (size_t)(g->width + 1) *
           (size_t)g->width  *
           (size_t)(g->height + 1) *
           (size_t)g->height / 4;
    // adjust kbest
    kbest = min(size, kbest);

    t_compute = clock();              // compute LRT
    r = (*compute[lrt_method])(g,kbest);
    if (output != NULL) {
        print_rectangles(output,r,kbest);
    }
    t_compute = clock() - t_compute;
    
    if (profiling_flag) {             // output time measurements and memory consumption
        fprintf(stderr,"time to read input file: %f (s)\n",((double)t_read)/CLOCKS_PER_SEC);
        fprintf(stderr,"time to compute LRT: %f (s)\n",((double)t_compute)/CLOCKS_PER_SEC);
        fprintf(stderr,"memory usage: %f (KB)\n",((double)get_memusage()) / 1024);
    }

    // release memory for rectangles and grid
    free_grid(g);
    free(r); 
    
    return 0;
}
