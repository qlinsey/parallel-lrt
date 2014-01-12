//
//  mk_grid.c
//  LRT
//
//  Copyright (c) 2013 The University of Sydney. All rights reserved.
//
//  command line tool for making a 2d data grid with outliers. 
//  

#include "lrt.h"

////////////////////////////////////////////////////////////////////////////
//  global variables
//

//
// program name (name of executable)
char *progname;

//
// file descriptor and name of output file. Default is stdout. 
FILE *output; 

//
// width of grid
int width=DEFAULT_WIDTH;

//
// height of grid
int height=DEFAULT_HEIGHT;

//
// maximum population
int max_population=DEFAULT_MAX_POPULATION;

//
// lambda
float lambda=DEFAULT_LAMBDA;

//
// outlier lambd 
float outlier_lambda=OUTLIER_LAMBDA;

//
// outlier rectangle size in rows and columns
int outlier_rows=DEFAULT_OUTLIER_ROWS;
int outlier_cols=DEFAULT_OUTLIER_COLS;

//
// outlier rectangle index 
int outlier_rowidx=DEFAULT_OUTLIER_RIDX;
int outlier_colidx=DEFAULT_OUTLIER_CIDX;

////////////////////////////////////////////////////////////////////////////
//  functions
//

//
// Compute poisson 
// from Knuth, via wikipedia complexity is linear in lambda
//
int simple_poisson(float lam) 
{
   float L = exp(-lam);
   int k = 0;
   float p = 1;
   do  {
      k++;
      p *= (float) rand() / (float)RAND_MAX;
   } while (p > L);
   return k - 1;
}

//
// generate grid and choose a region as outlier
//
grid *make_grid(int width, int height, int srows, int scols, int srowIdx,int scolIdx)
{
   grid *g = (grid *) ALLOC(grid);
   int i,j;
   
   g->width=width;
   g->height = height; 
   g->cells = (cell *) ALLOCV(cell,width*height);
   
   for (i = 0; i < width; i++) {
      for (j = 0; j < height; j++) {
         g->cells[j*width + i].n = max_population;
         int kc=0;
         while(kc==0){
           kc = simple_poisson(lambda);
           g->cells[j*width + i].k = kc;
         }
      }
      
   }
   //fill in small region 
   for(j=scolIdx;j<scolIdx+scols;j++){
      for(i=srowIdx;i<srowIdx+srows;i++){ 
         g->cells[j*width + i].n = max_population;
         int kc=0;
         while(kc==0){
           kc = simple_poisson(outlier_lambda);
           g->cells[j*width + i].k = kc;
         }
         DEBUG(2, "S(%d,%d) ",g->cells[j*width+i].n,g->cells[j*width+i].k);
      }
   }
   
   for (i = 0; i < width; i++) {
      for (j = 0; j < height; j++) {
         DEBUG(2,"(%d,%d) ",g->cells[j*width+i].n,g->cells[j*width+i].k);
      }
   }
   return g;
}

//
// parse command line arguments
// 
void parse_arguments(int argc, char ** argv)
{
   int c;
   
   opterr = 0;
   
   while ((c = getopt (argc, argv, "o:l:n:m:h:r:c:d:e:q")) != -1) {
      switch (c) {
         case 'n':
            width = atoi(optarg); 
            if (width <= 0 ) {
               error("width must greater than zero.");
            }
            break;
            
         case 'm':
            height = atoi(optarg);
            if (height <= 0 ) {
               error("height must greater than zero.");
            }
            break;
            
         case 'l':
            lambda = atof(optarg);
            break;
            
         case 'o': //specify outlier lambda 
            outlier_lambda = atof(optarg);
            if(outlier_lambda<=lambda){
               error("outlier region lambda must greater than normal lambda");
            }
            break;
            
         case 'r':
            outlier_rows = atoi(optarg);
            break;
            
         case 'c':
            outlier_cols = atoi(optarg);
            break;
            
         case 'd':
            outlier_rowidx = atoi(optarg);
            break;
            
         case 'e':
            outlier_colidx = atoi(optarg);
            break;
            
         case 'p':
            max_population = atoi(optarg);
            break;
            
         case 'h': 
            fprintf(stderr,"mk_grid: generate a random 2-d grids\n Version %s\n",VERSION);
            fprintf(stderr,"Common Options:\n");
            fprintf(stderr,"  -n <width>\n");
            fprintf(stderr,"  -m <height>\n");
            fprintf(stderr,"  -l <lambda>\n"); 
            fprintf(stderr,"  -o <outlierlambda>\n");
            fprintf(stderr,"  -r <outlierwidth>\n");
            fprintf(stderr,"  -c <outlierheight>\n");
            fprintf(stderr,"  -d <outlierrowidx>\n");
            fprintf(stderr,"  -e <outliercolumnidx>\n");
            fprintf(stderr,"  -p <max population>\n");
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
   if (optind >= argc) { 
      error("no output file provided."); 
   } 
   if (argc-optind > 1){
      error("only one output file permitted."); 
   }
   if ((output=fopen(argv[optind],"wb"))==NULL){
      error("cannot open output file %s.",optarg);
   }
}

//
// main program
// 
// 1. parses input parameters
// 2. generate grid
// 3. write result to file
//
int main(int argc, char **argv)
{   
   grid *g;                    
   
   progname = argv[0];               // set name of executable for error function
   
   parse_arguments(argc,argv);       // parse arguments
   
   g=make_grid(width,height,outlier_rows,outlier_cols,outlier_rowidx,outlier_colidx);        // make grid
   write_grid(output,g);             // write output 
   
   DEBUG(1,"gridsize %d x %d with lambda %f written.\n",width,height,lambda) 
   
   return 0;
}
