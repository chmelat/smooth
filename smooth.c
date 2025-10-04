/*  Balance of table points by polynom p-th degree at sliding window size
 *  n points. Least square principe.
 *  V5.1/2025-10-04/ Simplified Tikhonov method (removed adaptive_weights option)
 *  V5.0/2025-05-27/ Modularized version with separate method implementations
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <libgen.h>
#include <unistd.h>

#include "decomment.h"
#include "revision.h"
#include "tikhonov.h"
#include "polyfit.h"
#include "savgol.h"

#define BUF 512
#define N 5
#ifndef DPMAX
#define DPMAX 12
#endif
#define DP 2
#define LAMBDA_DEFAULT 0.1

/* Method flags */
#define METHOD_POLYFIT 0
#define METHOD_SAVGOL 1
#define METHOD_TIKHONOV 2

/* Local declare functions */
static void usage(void);
static void help(void);

/* Global variables */
char *progname;

int main(int argc, char **argv)
{
  char *filename;
  FILE *fp;
  double *x = NULL;
  double *y = NULL;
  int n=0;
  int i;
  int ch;
  int sp = N;
  int dp = DP;
  int abuf = 0;
  int method = METHOD_POLYFIT;
  double lambda = LAMBDA_DEFAULT;
  int auto_lambda = 0;
  int show_derivative = 0;

  progname = basename(argv[0]);

/* None arguments */  
  if (argc == 1) {
    usage();
    exit (EXIT_FAILURE);
  }

/* Options command line */
  while ( (ch = getopt(argc, argv, "n:p:m:l:dh?")) != -1 ) {
    switch (ch) {
      case 'n':
        sp = atoi(optarg);
        break;
      case 'p':
        dp = atoi(optarg);
        break;
      case 'm':
        if (optarg[0] >= '0' && optarg[0] <= '9') {
          int method_num = atoi(optarg);
          if (method_num >= METHOD_POLYFIT && method_num <= METHOD_TIKHONOV) {
            method = method_num;
          } else {
            fprintf(stderr, "Unknown method number: %d\n", method_num);
            fprintf(stderr, "Valid methods: 0 (polyfit), 1 (savgol), 2 (tikhonov)\n");
            exit(EXIT_FAILURE);
          }
        } else {
          if (strcmp(optarg, "polyfit") == 0) {
            method = METHOD_POLYFIT;
          } else if (strcmp(optarg, "savgol") == 0) {
            method = METHOD_SAVGOL;
          } else if (strcmp(optarg, "tikhonov") == 0) {
            method = METHOD_TIKHONOV;
          } else {
            fprintf(stderr, "Unknown method: %s\n", optarg);
            fprintf(stderr, "Valid methods: polyfit, savgol, tikhonov, 0, 1, 2\n");
            exit(EXIT_FAILURE);
          }
        }
        break;
      case 'l':
        if (strcmp(optarg, "auto") == 0) {
          auto_lambda = 1;
        } else {
          lambda = atof(optarg);
          if (lambda < 0) {
            fprintf(stderr, "Lambda must be non-negative!\n");
            exit(EXIT_FAILURE);
          }
        }
        break;
      case 'd':
        show_derivative = 1;
        break;
      case 'h':
        help();
        exit(EXIT_SUCCESS);
      case '?':
        help();
        exit(EXIT_SUCCESS);
      default:
        fprintf(stderr,"Unknown option!\n");
        help();
        exit(EXIT_FAILURE);
    }
  }

  if ( sp<3 || !(sp%2) ) {
    fprintf(stderr,"Incorrect number points in moving windows (odd >= 3)!\n");
    exit (EXIT_FAILURE);
  }

  if (dp < 0 || dp > DPMAX) {
    fprintf(stderr,"Incorrect degree of approx. polynom, (0<=p<=%d)!\n",DPMAX);
    exit (EXIT_FAILURE);
  }

  if (dp > 6) {
    fprintf(stderr, "Warning: High polynomial degree (%d) may cause numerical instability\n", dp);
  }

/* Argument - filename */  
  if (argv[optind] == NULL) {
    fprintf(stderr,"Data file missing!\n");
    exit (EXIT_FAILURE);
  }
  else 
    filename = argv[optind];

  fp = decomment(filename);

/* Read data table from file */
  {
    double rx, ry;
    while (fscanf(fp,"%lf %lf",&rx,&ry) != EOF) {
      if (n == abuf) {
        abuf += BUF;
        
        /* Safe realloc for x */
        double *temp_x = (double *)realloc(x, abuf*sizeof(double));
        if (temp_x == NULL) {
          free(x);
          free(y);
          fprintf(stderr,"No memory for data table!\n");
          exit(EXIT_FAILURE);
        }
        x = temp_x;
        
        /* Safe realloc for y */
        double *temp_y = (double *)realloc(y, abuf*sizeof(double));
        if (temp_y == NULL) {
          free(x);
          free(y);
          fprintf(stderr,"No memory for data table!\n");
          exit(EXIT_FAILURE);
        }
        y = temp_y;
      }
      x[n] = rx;
      y[n] = ry; 
      n++;
    }
    fclose(fp);
  }

  if (n < sp && method != METHOD_TIKHONOV) {
    fprintf(stderr,"Need more data (n < %d)!\n",sp);
    exit (EXIT_FAILURE);
  }

  /* Process data according to selected method */
  switch (method) {
    case METHOD_TIKHONOV:
      {
        TikhonovResult *result;
        
        /* Find optimal lambda if requested */
        if (auto_lambda) {
          lambda = find_optimal_lambda_gcv(x, y, n);
          printf("# Automatic lambda selection using GCV: lambda = %.6e\n", lambda);
        }

        /* Apply Tikhonov smoothing */
        result = tikhonov_smooth(x, y, n, lambda);

        if (result == NULL) {
          fprintf(stderr, "Tikhonov smoothing failed!\n");
          exit(EXIT_FAILURE);
        }

        /* Output header with functional information */
        printf("# Data smooth - Tikhonov regularization with lambda = %.6lG\n", lambda);
        printf("# Functional J = %.6lG (Data: %.6lG + Regularization: %.6lG)\n", 
               result->total_functional, result->data_term, result->regularization_term);
        printf("# Data/Total ratio = %.3f, Regularization/Total ratio = %.3f\n",
               result->data_term / result->total_functional,
               result->regularization_term / result->total_functional);
        if (show_derivative) {
          printf("#    x          y          y'\n");
        } else {
          printf("#    x          y\n");
        }

        /* Output results */
        for (i = 0; i < n; i++) {
          if (show_derivative) {
            printf("%12.8lG %10.6lG %10.6lG\n", x[i], result->y_smooth[i], result->y_deriv[i]);
          } else {
            printf("%12.8lG %10.6lG\n", x[i], result->y_smooth[i]);
          }
        }

        /* Clean up */
        free_tikhonov_result(result);
      }
      break;
      
    case METHOD_SAVGOL:
      {
        SavgolResult *result;
        
        result = savgol_smooth(x, y, n, sp, dp);
        
        if (result == NULL) {
          fprintf(stderr, "Savitzky-Golay smoothing failed!\n");
          exit(EXIT_FAILURE);
        }
        
        printf("# Data smooth - Savitzky-Golay filter, poly deg %d from %d points of moving window\n", dp, sp);
        if (show_derivative) {
          printf("#    x          y          y'\n");
        } else {
          printf("#    x          y\n");
        }
        
        for (i = 0; i < n; i++) {
          if (show_derivative) {
            printf("%12.8lG %10.6lG %10.6lG\n", x[i], result->y_smooth[i], result->y_deriv[i]);
          } else {
            printf("%12.8lG %10.6lG\n", x[i], result->y_smooth[i]);
          }
        }
        
        free_savgol_result(result);
      }
      break;
      
    case METHOD_POLYFIT:
    default:
      {
        PolyfitResult *result;
        
        result = polyfit_smooth(x, y, n, sp, dp);
        
        if (result == NULL) {
          fprintf(stderr, "Polynomial fitting failed!\n");
          exit(EXIT_FAILURE);
        }
        
        printf("# Data smooth - aprox. pol. %ddg from %d points of moving window (least square)\n", dp, sp);
        if (show_derivative) {
          printf("#    x          y          y'\n");
        } else {
          printf("#    x          y\n");
        }
        
        for (i = 0; i < n; i++) {
          if (show_derivative) {
            printf("%12.8lG %10.6lG %10.6lG\n", x[i], result->y_smooth[i], result->y_deriv[i]);
          } else {
            printf("%12.8lG %10.6lG\n", x[i], result->y_smooth[i]);
          }
        }
        
        free_polyfit_result(result);
      }
      break;
  }

  /* Clean up */
  free(x);
  free(y);

  return EXIT_SUCCESS;
}

/* Usage */
static void usage(void)
{
  fprintf(stderr,"Data smooth by approximation polynom from moving window of data (least square)\n");
  fprintf(stderr, "Usage: %s [options] data_file\n" ,progname);
  fprintf(stderr,"-h or -? for help\n");
}

/* Help */
static void help(void)
{
  static char *msg[] = {
    "-h, -?\tHelp",
    "-n\tPoints in moving window, default 5 (not used for Tikhonov)",  
    "-p\tDegree of approx. polynom, default 2",
    "-m\tMethod: 0 (polyfit, default), 1 (savgol), or 2 (tikhonov)",
    "-l\tLambda regularization parameter for Tikhonov method, default 0.1",
    "\tUse '-l auto' for automatic lambda selection using GCV",
    "-d\tShow first derivative in output",
    "\nMethods:",
    "  polyfit  - Local polynomial fitting (least squares)",
    "  savgol   - Savitzky-Golay filter for better derivatives", 
    "  tikhonov - Global smoothing with regularization (correct for non-uniform grids)",
    "\nExamples:",
    "  smooth -m 0 -n 7 -p 3 data.txt     # Polyfit with 7-point window, 3rd degree",
    "  smooth -m 1 -n 5 -p 2 -d data.txt  # Savitzky-Golay with derivatives",
    "  smooth -m 2 -l 0.01 data.txt       # Tikhonov with manual lambda",
    "  smooth -m 2 -l auto -d data.txt    # Tikhonov with auto lambda and derivatives",
    0
  };
  char **p = msg;
  usage();
  fprintf(stderr,"%s, Version %s, %s\n\n", progname, VERSION, REVDATE);
  while (*p)
    fprintf(stderr, "%s\n", *p++);
}
