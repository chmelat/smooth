/*  Balance of table points by polynom p-th degree at sliding window size
 *  n points. Least square principe.
 *  V5.0/2025-05-27/ Modularized version with separate method implementations
 *  V4.2/2025-05-26/ Added adaptive weights option for Tikhonov regularization
 *  V4.1/2025-05-24/ Tikhonov regularization moved to separate module with LAPACK
 *  V3.0/2025-03-08/ add Savitzky-Golay filter for better derivative estimation
 *  V2.1/2020-01-10/ add BUF and abuf as group allocation */ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <libgen.h>  /* char *basename(char *), char *dirname(char *) */
#include <unistd.h>  /* UNIX standard function definitions including getopt */

#include "decomment.h"
#include "revision.h"
#include "tikhonov.h"
#include "polyfit.h"
#include "savgol.h"

#define BUF 512  /* Block data for memory allocation */
#define N 5  /* Odd number >= 3, data for least square */
#ifndef DPMAX
#define DPMAX 12  /* Max degree of approximation polynom - configurable */
#endif
#define DP 2  /* Default degree of approx. polynom */
#define LAMBDA_DEFAULT 0.1  /* Default regularization parameter for Tikhonov */

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
  double *x = NULL;  /* Abscissa */
  double *y = NULL;  /* Ordinate */
  int n=0;  /* Number of table points */
  int i;
  int ch;
  int sp = N;  /* Sliding points for least square (default) */
  int dp = DP;  /* Degree of approximation polynom */
  int abuf = 0; /* Dynamic data allocation */
  int method = METHOD_POLYFIT; /* Default method is POLYFIT filter */
  double lambda = LAMBDA_DEFAULT; /* Regularization parameter for Tikhonov */
  int auto_lambda = 0; /* Flag for automatic lambda selection */
  int adaptive_weights = 0; /* Flag for adaptive weights in Tikhonov method */
  int show_derivative = 0; /* Flag for showing derivative in output */

  progname = basename(argv[0]);

/* None arguments */  
  if (argc == 1) {
    usage();
    exit (EXIT_FAILURE);
  }

/* Options command line */
  while ( (ch = getopt(argc, argv, "n:p:m:l:adh?")) != -1 ) {
    switch (ch) {
      case 'n':
        sp = atoi(optarg);
        break;
      case 'p':
        dp = atoi(optarg);
        break;
      case 'm':
        /* Check if input is a number */
        if (optarg[0] >= '0' && optarg[0] <= '9') {
          /* Parse as number */
          int method_num = atoi(optarg);
          if (method_num >= METHOD_POLYFIT && method_num <= METHOD_TIKHONOV) {
            method = method_num;
          } else {
            fprintf(stderr, "Unknown method number: %d\n", method_num);
            fprintf(stderr, "Valid methods: 0 (polyfit), 1 (savgol), 2 (tikhonov)\n");
            exit(EXIT_FAILURE);
          }
        } else {
          /* Traditional string method names for backward compatibility */
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
      case 'a':
        adaptive_weights = 1;
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
        x = (double *)realloc((void *)x,abuf*sizeof(double));
        if (x == NULL) {
          fprintf(stderr,"No memory for data table!\n");
          exit(EXIT_FAILURE);
        }
        y = (double *)realloc((void *)y,abuf*sizeof(double));
        if (y == NULL) {
          fprintf(stderr,"No memory for data table!\n");
          exit(EXIT_FAILURE);
        }
      }
      x[n] = rx, y[n] = ry; 
      n++;
    }
    fclose(fp);
  }

  if (n < sp && method != METHOD_TIKHONOV) {
    fprintf(stderr,"Need moore data (n < %d)!\n",sp);
    exit (EXIT_FAILURE);
  }

  /* Process data according to selected method */
  switch (method) {
    case METHOD_TIKHONOV:
      {
        /* Tikhonov regularization - global smoothing */
        TikhonovResult *result;
        
        /* Validate adaptive weights flag */
        if (adaptive_weights && method != METHOD_TIKHONOV) {
          fprintf(stderr, "Warning: -a flag only applies to Tikhonov method, ignoring.\n");
          adaptive_weights = 0;
        }
        
        /* Find optimal lambda if requested */
        if (auto_lambda) {
          lambda = find_optimal_lambda_gcv(x, y, n);
          printf("# Automatic lambda selection using GCV: lambda = %.6e\n", lambda);
        }

        /* Apply Tikhonov smoothing with adaptive weights option */
        result = tikhonov_smooth_adaptive(x, y, n, lambda, adaptive_weights);

        if (result == NULL) {
          fprintf(stderr, "Tikhonov smoothing failed!\n");
          exit(EXIT_FAILURE);
        }

        /* Output header with functional information */
        printf("# Data smooth - Tikhonov regularization with lambda = %.6lG\n", lambda);
        printf("# Adaptive weights: %s\n", adaptive_weights ? "ENABLED" : "DISABLED");
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
        /* Savitzky-Golay filtering method */
        SavgolResult *result;
        
        result = savgol_smooth(x, y, n, sp, dp);
        
        if (result == NULL) {
          fprintf(stderr, "Savitzky-Golay smoothing failed!\n");
          exit(EXIT_FAILURE);
        }
        
        /* Output header */
        printf("# Data smooth - Savitzky-Golay filter, poly deg %d from %d points of moving window\n", dp, sp);
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
        free_savgol_result(result);
      }
      break;
      
    case METHOD_POLYFIT:
    default:
      {
        /* Polynomial fitting method */
        PolyfitResult *result;
        
        result = polyfit_smooth(x, y, n, sp, dp);
        
        if (result == NULL) {
          fprintf(stderr, "Polynomial fitting failed!\n");
          exit(EXIT_FAILURE);
        }
        
        /* Output header */
        printf("# Data smooth - aprox. pol. %ddg from %d points of moving window (least square)\n", dp, sp);
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
  static char *msg[] = {"-h, -?\tHelp",
                 "-n\tPoints in moving window, default 5 (not used for Tikhonov)",  
                 "-p\tDegree of approx. polynom, default 2",
                 "-m\tMethod: 0 (polyfit, default), 1 (savgol), or 2 (tikhonov)",
                 "-l\tLambda regularization parameter for Tikhonov method, default 0.1",
                 "\tUse '-l auto' for automatic lambda selection using GCV",
                 "-a\tEnable adaptive local weights for non-uniform grids (Tikhonov only)",
                 "-d\tShow first derivative in output",
                 "\nMethods:",
                 "  polyfit  - Local polynomial fitting (least squares)",
                 "  savgol   - Savitzky-Golay filter for better derivatives", 
                 "  tikhonov - Global smoothing with regularization",
                 "\nTikhonov options:",
                 "  Default: Uses average coefficient method for non-uniform grids",
                 "  With -a: Uses local interval-based weights (better for highly non-uniform grids)",
                 "\nExamples:",
                 "  smooth -m 0 -n 7 -p 3 data.txt     # Polyfit with 7-point window, 3rd degree",
                 "  smooth -m 1 -n 5 -p 2 -d data.txt  # Savitzky-Golay with derivatives",
                 "  smooth -m 2 -l 0.01 data.txt       # Tikhonov with manual lambda",
                 "  smooth -m 2 -l auto -d data.txt    # Tikhonov with auto lambda and derivatives",
                 "  smooth -m 2 -l 0.01 -a data.txt    # Tikhonov with adaptive weights",
	           0};
  char **p = msg;
  usage();
  fprintf(stderr,"%s, Version %s, %s\n\n", progname, VERSION, REVDATE);
  while (*p)
    fprintf(stderr, "%s\n", *p++);
}
