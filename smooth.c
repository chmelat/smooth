/*  Balance of table points by polynom p-th degree at sliding window size
 *  n points. Least square principe.
 *  V5.5/2025-11-03/ Added Butterworth filtfilt method
 *  V5.2/2025-10-13/ Added -g option for detailed grid uniformity analysis
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
#include "butterworth.h"
#include "grid_analysis.h"
#include "timestamp.h"

#define BUF 512
#define N 5
#ifndef DPMAX
#define DPMAX 12
#endif
#define DP 2
#define LAMBDA_DEFAULT 0.1
#define CUTOFF_DEFAULT 0.2

/* Method flags */
#define METHOD_POLYFIT 0
#define METHOD_SAVGOL 1
#define METHOD_TIKHONOV 2
#define METHOD_BUTTERWORTH 3

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
  double cutoff_freq = CUTOFF_DEFAULT;
  int auto_cutoff = 0;
  int show_derivative = 0;
  int show_grid_analysis = 0;
  int y_column = 2;  /* Default: second column (1=first, 2=second, etc.) */
  int timestamp_mode = 0;  /* Flag for timestamp input mode */
  char **timestamp_strings = NULL;  /* Array of timestamp strings */
  TimestampContext *ts_ctx = NULL;  /* Timestamp conversion context */

  progname = basename(argv[0]);

/* None arguments - default to stdin */
  if (argc == 1) {
    /* Will read from stdin below */
  }

/* Options command line */
  while ( (ch = getopt(argc, argv, "n:p:m:l:f:k:dgTh?")) != -1 ) {
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
          if (method_num >= METHOD_POLYFIT && method_num <= METHOD_BUTTERWORTH) {
            method = method_num;
          } else {
            fprintf(stderr, "Unknown method number: %d\n", method_num);
            fprintf(stderr, "Valid methods: 0 (polyfit), 1 (savgol), 2 (tikhonov), 3 (butterworth)\n");
            exit(EXIT_FAILURE);
          }
        } else {
          if (strcmp(optarg, "polyfit") == 0) {
            method = METHOD_POLYFIT;
          } else if (strcmp(optarg, "savgol") == 0) {
            method = METHOD_SAVGOL;
          } else if (strcmp(optarg, "tikhonov") == 0) {
            method = METHOD_TIKHONOV;
          } else if (strcmp(optarg, "butterworth") == 0) {
            method = METHOD_BUTTERWORTH;
          } else {
            fprintf(stderr, "Unknown method: %s\n", optarg);
            fprintf(stderr, "Valid methods: polyfit, savgol, tikhonov, butterworth, 0, 1, 2, 3\n");
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
      case 'f':
        if (strcmp(optarg, "auto") == 0) {
          auto_cutoff = 1;
        } else {
          cutoff_freq = atof(optarg);
          if (cutoff_freq <= 0.0 || cutoff_freq >= 1.0) {
            fprintf(stderr, "Cutoff frequency must be in range (0, 1)!\n");
            fprintf(stderr, "where fc = 1 corresponds to Nyquist frequency (fs/2)\n");
            exit(EXIT_FAILURE);
          }
        }
        break;
      case 'k':
        y_column = atoi(optarg);
        if (y_column < 1) {
          fprintf(stderr, "Column number must be >= 1!\n");
          exit(EXIT_FAILURE);
        }
        break;
      case 'd':
        show_derivative = 1;
        break;
      case 'g':
        show_grid_analysis = 1;
        break;
      case 'T':
        timestamp_mode = 1;
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

/* Argument - filename or stdin */
  if (argv[optind] == NULL || strcmp(argv[optind], "-") == 0) {
    /* Read from stdin */
    fp = stdin;
    filename = "stdin";
  }
  else {
    filename = argv[optind];
    fp = decomment(filename);
  }

/* Read data table from file */
  {
    char line[4096];
    int line_number = 0;

    /* Allocate arrays based on mode */
    if (timestamp_mode) {
      timestamp_strings = malloc(BUF * sizeof(char*));
      if (!timestamp_strings) {
        fprintf(stderr, "No memory for timestamp strings!\n");
        exit(EXIT_FAILURE);
      }
      /* Allocate y array for timestamp mode */
      y = malloc(BUF * sizeof(double));
      if (!y) {
        fprintf(stderr, "No memory for data table!\n");
        free(timestamp_strings);
        exit(EXIT_FAILURE);
      }
      abuf = BUF;  /* Track allocated size */
    }

    while (fgets(line, sizeof(line), fp) != NULL) {
      line_number++;

      /* Parse based on mode */
      if (timestamp_mode) {
        /* Timestamp mode: format is "YYYY-MM-DD HH:MM:SS.fff value" (date and time separated by space) */
        char date_str[20], time_str[20];
        char timestamp_str[100];
        double y_value;

        /* Skip leading whitespace */
        char *ptr = line;
        while (*ptr == ' ' || *ptr == '\t') ptr++;

        /* Skip empty lines */
        if (*ptr == '\n' || *ptr == '\0') continue;

        /* Parse date, time, and y-value
         * Format with space separator: "YYYY-MM-DD HH:MM:SS.fff value"
         * Format with T separator: "YYYY-MM-DDTHH:MM:SS.fff value"
         * sscanf will parse the timestamp correctly if it uses T, or as two tokens if it uses space */
        int parsed = sscanf(ptr, "%19s %19s %lf", date_str, time_str, &y_value);

        if (parsed == 2) {
          /* Only got date and value - maybe format is "YYYY-MM-DDTHH:MM:SS.fff value" (T separator) */
          y_value = atof(time_str);  /* time_str actually contains the y-value */
          snprintf(timestamp_str, sizeof(timestamp_str), "%s", date_str);
        } else if (parsed == 3) {
          /* Check if time_str starts with '.' - this means T format split the subseconds */
          if (time_str[0] == '.') {
            /* Format was "YYYY-MM-DDTHH:MM:SS.fff" - combine date_str (which has everything up to '.') and time_str (which has '.fff') */
            snprintf(timestamp_str, sizeof(timestamp_str), "%s%s", date_str, time_str);
          } else {
            /* Got date, time, and value - format is "YYYY-MM-DD HH:MM:SS.fff value" (space separator) */
            snprintf(timestamp_str, sizeof(timestamp_str), "%s %s", date_str, time_str);
          }
        } else {
          /* Could not parse - skip this line silently */
          continue;
        }

        /* Reallocate arrays if needed */
        if (n == abuf) {
          abuf += BUF;

          /* Realloc timestamp_strings */
          char **temp_ts = (char**)realloc(timestamp_strings, abuf * sizeof(char*));
          if (!temp_ts) {
            fprintf(stderr, "No memory for timestamp strings!\n");
            for (int i = 0; i < n; i++) free(timestamp_strings[i]);
            free(timestamp_strings);
            free(y);
            exit(EXIT_FAILURE);
          }
          timestamp_strings = temp_ts;

          /* Realloc y array */
          double *temp_y = (double*)realloc(y, abuf * sizeof(double));
          if (!temp_y) {
            fprintf(stderr, "No memory for data table!\n");
            for (int i = 0; i < n; i++) free(timestamp_strings[i]);
            free(timestamp_strings);
            free(y);
            exit(EXIT_FAILURE);
          }
          y = temp_y;
        }

        /* Store timestamp string and y-value */
        timestamp_strings[n] = strdup(timestamp_str);
        if (!timestamp_strings[n]) {
          fprintf(stderr, "No memory for timestamp string!\n");
          for (int i = 0; i < n; i++) free(timestamp_strings[i]);
          free(timestamp_strings);
          free(y);
          exit(EXIT_FAILURE);
        }
        y[n] = y_value;
        n++;

      } else {
        /* Normal mode: parse numeric columns */
        double values[100];  /* Support up to 100 columns */
        int col_count = 0;
        char *ptr = line;
        char *endptr;

        /* Skip leading whitespace */
        while (*ptr == ' ' || *ptr == '\t') ptr++;

        /* Skip empty lines */
        if (*ptr == '\n' || *ptr == '\0') continue;

        /* Parse all numeric values on the line */
        while (*ptr != '\n' && *ptr != '\0' && col_count < 100) {
          errno = 0;
          values[col_count] = strtod(ptr, &endptr);

          /* Check if parsing was successful */
          if (ptr == endptr || errno != 0) {
            break;  /* Not a number, stop parsing this line */
          }

          col_count++;
          ptr = endptr;

          /* Skip whitespace between columns */
          while (*ptr == ' ' || *ptr == '\t') ptr++;
        }

        /* Check if we have enough columns */
        if (col_count < 1) {
          /* No valid numbers found, skip this line */
          continue;
        }

        if (col_count < y_column) {
          fprintf(stderr, "Error: Line %d has only %d column(s), but column %d was requested for y-data\n",
                  line_number, col_count, y_column);
          free(x);
          free(y);
          if (fp != stdin) fclose(fp);
          exit(EXIT_FAILURE);
        }

        /* Reallocate arrays if needed */
        if (n == abuf) {
          abuf += BUF;

          /* Safe realloc for x */
          double *temp_x = (double *)realloc(x, abuf*sizeof(double));
          if (temp_x == NULL) {
            free(x);
            free(y);
            fprintf(stderr,"No memory for data table!\n");
            if (fp != stdin) fclose(fp);
            exit(EXIT_FAILURE);
          }
          x = temp_x;

          /* Safe realloc for y */
          double *temp_y = (double *)realloc(y, abuf*sizeof(double));
          if (temp_y == NULL) {
            free(x);
            free(y);
            fprintf(stderr,"No memory for data table!\n");
            if (fp != stdin) fclose(fp);
            exit(EXIT_FAILURE);
          }
          y = temp_y;
        }

        /* Store x from first column and y from specified column */
        x[n] = values[0];           /* Always first column for x */
        y[n] = values[y_column - 1]; /* User specifies 1-indexed, array is 0-indexed */
        n++;
      }  /* end else (normal mode) */
    }  /* end while (fgets) */

    if (fp != stdin) {
      fclose(fp);
    }

    /* Convert timestamps to relative time if in timestamp mode */
    if (timestamp_mode) {
      if (n == 0) {
        fprintf(stderr, "Error: No valid data points found\n");
        for (int i = 0; i < n; i++) free(timestamp_strings[i]);
        free(timestamp_strings);
        free(y);
        exit(EXIT_FAILURE);
      }

      int first_error_line = -1;
      ts_ctx = convert_timestamps_to_relative(timestamp_strings, n, &x, &first_error_line);

      if (ts_ctx == NULL) {
        fprintf(stderr, "Error: No valid timestamps found in input\n");
        if (first_error_line > 0) {
          fprintf(stderr, "First invalid timestamp at line %d\n", first_error_line);
        }
        for (int i = 0; i < n; i++) free(timestamp_strings[i]);
        free(timestamp_strings);
        free(y);
        exit(EXIT_FAILURE);
      }

      /* Print warning if some timestamps were invalid */
      if (ts_ctx->errors_encountered > 0) {
        fprintf(stderr, "Warning: Skipped %d line(s) with invalid timestamps (first error at line %d)\n",
                ts_ctx->errors_encountered, first_error_line);
      }

      /* Update n to reflect actual number of valid points */
      n = ts_ctx->n;

      /* Free original timestamp_strings array (ts_ctx has its own copy) */
      for (int i = 0; i < n + ts_ctx->errors_encountered; i++) {
        free(timestamp_strings[i]);
      }
      free(timestamp_strings);
      timestamp_strings = NULL;
    }
  }

  if (n < sp && method != METHOD_TIKHONOV) {
    fprintf(stderr,"Need more data (n < %d)!\n",sp);
    exit (EXIT_FAILURE);
  }

  /* Perform grid uniformity analysis (always done before smoothing) */
  GridAnalysis *grid_info = analyze_grid(x, n, 0);  /* store_spacings=0 - no histogram needed */
  if (grid_info == NULL) {
    fprintf(stderr, "Error: Grid analysis failed\n");
    free(x);
    free(y);
    exit(EXIT_FAILURE);
  }

  /* If -g flag: show detailed analysis and exit (no smoothing) */
  if (show_grid_analysis) {
    print_grid_analysis(grid_info, 1, "# ");  /* verbose=1 - basic stats + recommendations only */
    free_grid_analysis(grid_info);
    free(x);
    free(y);
    return EXIT_SUCCESS;
  }

  /* Show warnings if grid has reliability concerns */
  if (grid_info->reliability_warning) {
    printf("# Grid analysis warnings:\n");
    print_grid_analysis(grid_info, 0, "# ");  /* verbose=0 - basic statistics + warnings */
  }

  /* Process data according to selected method */
  switch (method) {
    case METHOD_TIKHONOV:
      {
        TikhonovResult *result;
        
        /* Find optimal lambda if requested */
        if (auto_lambda) {
          lambda = find_optimal_lambda_gcv(x, y, n, grid_info);
          printf("# Automatic lambda selection using GCV: lambda = %.6e\n", lambda);
        }

        /* Apply Tikhonov smoothing */
        result = tikhonov_smooth(x, y, n, lambda, grid_info);

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
        if (timestamp_mode) {
          if (show_derivative) {
            printf("# Derivative units: dy/dt (t in seconds)\n");
            printf("#    timestamp          y          y'\n");
          } else {
            printf("#    timestamp          y\n");
          }
        } else {
          if (show_derivative) {
            printf("#    x          y          y'\n");
          } else {
            printf("#    x          y\n");
          }
        }

        /* Output results */
        for (i = 0; i < n; i++) {
          if (timestamp_mode) {
            if (show_derivative) {
              printf("%s %10.6lG %10.6lG\n", ts_ctx->original_timestamps[i], result->y_smooth[i], result->y_deriv[i]);
            } else {
              printf("%s %10.6lG\n", ts_ctx->original_timestamps[i], result->y_smooth[i]);
            }
          } else {
            if (show_derivative) {
              printf("%12.8lG %10.6lG %10.6lG\n", x[i], result->y_smooth[i], result->y_deriv[i]);
            } else {
              printf("%12.8lG %10.6lG\n", x[i], result->y_smooth[i]);
            }
          }
        }

        /* Clean up */
        free_tikhonov_result(result);
      }
      break;
      
    case METHOD_SAVGOL:
      {
        SavgolResult *result;
        
        result = savgol_smooth(x, y, n, sp, dp, grid_info);
        
        if (result == NULL) {
          fprintf(stderr, "Savitzky-Golay smoothing failed!\n");
          exit(EXIT_FAILURE);
        }
        
        printf("# Data smooth - Savitzky-Golay filter, poly deg %d from %d points of moving window\n", dp, sp);
        if (timestamp_mode) {
          if (show_derivative) {
            printf("# Derivative units: dy/dt (t in seconds)\n");
            printf("#    timestamp          y          y'\n");
          } else {
            printf("#    timestamp          y\n");
          }
        } else {
          if (show_derivative) {
            printf("#    x          y          y'\n");
          } else {
            printf("#    x          y\n");
          }
        }

        for (i = 0; i < n; i++) {
          if (timestamp_mode) {
            if (show_derivative) {
              printf("%s %10.6lG %10.6lG\n", ts_ctx->original_timestamps[i], result->y_smooth[i], result->y_deriv[i]);
            } else {
              printf("%s %10.6lG\n", ts_ctx->original_timestamps[i], result->y_smooth[i]);
            }
          } else {
            if (show_derivative) {
              printf("%12.8lG %10.6lG %10.6lG\n", x[i], result->y_smooth[i], result->y_deriv[i]);
            } else {
              printf("%12.8lG %10.6lG\n", x[i], result->y_smooth[i]);
            }
          }
        }
        
        free_savgol_result(result);
      }
      break;
      
    case METHOD_BUTTERWORTH:
      {
        ButterworthResult *result;

        /* Warn if derivative output was requested */
        if (show_derivative) {
          fprintf(stderr, "# WARNING: Derivative output (-d) not supported for Butterworth method\n");
          fprintf(stderr, "#          Outputting smoothed values only\n");
        }

        /* Apply Butterworth filtfilt */
        result = butterworth_filtfilt(x, y, n, cutoff_freq, auto_cutoff, grid_info);

        if (result == NULL) {
          fprintf(stderr, "Butterworth filtering failed!\n");
          exit(EXIT_FAILURE);
        }

        /* Output header */
        printf("# Data smooth - Butterworth filter (order %d, filtfilt)\n", result->order);
        printf("# Sample rate: fs = %.6lG\n", result->sample_rate);
        printf("# Nyquist frequency: f_Nyquist = %.6lG (= fs/2)\n", result->sample_rate / 2.0);
        printf("# Normalized cutoff frequency: fc = %.6lG (where 1.0 = f_Nyquist)\n", result->cutoff_freq);
        printf("# Actual cutoff frequency: f_cutoff = %.6lG (= fc × f_Nyquist)\n",
               result->cutoff_freq * result->sample_rate / 2.0);
        printf("# Effective order after filtfilt: %d\n", 2 * result->order);
        if (timestamp_mode) {
          printf("#    timestamp          y\n");
        } else {
          printf("#    x          y\n");
        }

        /* Output results - Butterworth doesn't have derivatives */
        for (i = 0; i < n; i++) {
          if (timestamp_mode) {
            printf("%s %10.6lG\n", ts_ctx->original_timestamps[i], result->y_smooth[i]);
          } else {
            printf("%12.8lG %10.6lG\n", x[i], result->y_smooth[i]);
          }
        }

        /* Clean up */
        free_butterworth_result(result);
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
        if (timestamp_mode) {
          if (show_derivative) {
            printf("# Derivative units: dy/dt (t in seconds)\n");
            printf("#    timestamp          y          y'\n");
          } else {
            printf("#    timestamp          y\n");
          }
        } else {
          if (show_derivative) {
            printf("#    x          y          y'\n");
          } else {
            printf("#    x          y\n");
          }
        }

        for (i = 0; i < n; i++) {
          if (timestamp_mode) {
            if (show_derivative) {
              printf("%s %10.6lG %10.6lG\n", ts_ctx->original_timestamps[i], result->y_smooth[i], result->y_deriv[i]);
            } else {
              printf("%s %10.6lG\n", ts_ctx->original_timestamps[i], result->y_smooth[i]);
            }
          } else {
            if (show_derivative) {
              printf("%12.8lG %10.6lG %10.6lG\n", x[i], result->y_smooth[i], result->y_deriv[i]);
            } else {
              printf("%12.8lG %10.6lG\n", x[i], result->y_smooth[i]);
            }
          }
        }

        free_polyfit_result(result);
      }
      break;
  }

  /* Clean up */
  free_grid_analysis(grid_info);
  free(x);
  free(y);
  if (timestamp_mode && ts_ctx) {
    free_timestamp_context(ts_ctx);
  }

  return EXIT_SUCCESS;
}

/* Usage */
static void usage(void)
{
  fprintf(stderr,"Data smooth by approximation polynom from moving window of data (least square)\n");
  fprintf(stderr, "Usage: %s [options] [data_file|-]\n" ,progname);
  fprintf(stderr, "If data_file is omitted or '-', reads from stdin\n");
  fprintf(stderr,"-h or -? for help\n");
}

/* Help */
static void help(void)
{
  static char *msg[] = {
    "-h, -?\tHelp",
    "-n\tPoints in moving window, default 5 (not used for Tikhonov or Butterworth)",
    "-p\tDegree of approx. polynom, default 2 (not used for Butterworth)",
    "-m\tMethod: 0 (polyfit, default), 1 (savgol), 2 (tikhonov), or 3 (butterworth)",
    "-l\tLambda regularization parameter for Tikhonov method, default 0.1",
    "\tUse '-l auto' for automatic lambda selection using GCV",
    "-f\tNormalized cutoff frequency for Butterworth filter, default 0.2",
    "\tRange: 0 < fc < 1, where fc = f_cutoff / f_Nyquist",
    "\t(fc = 1 corresponds to Nyquist frequency = f_sample/2)",
    "\tUse '-f auto' for automatic cutoff selection",
    "-k\tColumn number for y-data (default 2), x-data always from column 1",
    "\tColumns are numbered starting from 1",
    "-T\tTimestamp mode: first column is RFC3339-style timestamp, second is y-value",
    "\tSupports formats: YYYY-MM-DD HH:MM:SS[.fff] or YYYY-MM-DDTHH:MM:SS[.fff]",
    "\tTimestamps are converted to seconds for smoothing, original format preserved in output",
    "-d\tShow first derivative in output (not available for Butterworth)",
    "\tWith -T flag, derivatives are in units dy/dt where t is in seconds",
    "-g\tShow detailed grid uniformity analysis",
    "\nMethods:",
    "  polyfit     - Local polynomial fitting (least squares)",
    "  savgol      - Savitzky-Golay filter for better performance at uniform grids",
    "  tikhonov    - Global smoothing with regularization (correct for non-uniform grids)",
    "  butterworth - Digital Butterworth low-pass filter with filtfilt (zero-phase)",
    "\nExamples:",
    "  smooth -m 0 -n 7 -p 3 data.txt        # Polyfit with 7-point window, 3rd degree",
    "  smooth -m 1 -n 5 -p 2 -d data.txt     # Savitzky-Golay with derivatives",
    "  smooth -m 2 -l 0.01 data.txt          # Tikhonov with manual lambda",
    "  smooth -m 2 -l auto -d data.txt       # Tikhonov with auto lambda and derivatives",
    "  smooth -m 3 -f 0.1 data.txt           # Butterworth with fc=0.1",
    "  smooth -m 3 -f auto data.txt          # Butterworth with auto cutoff",
    "  smooth -k 3 -m 2 -l auto data.txt     # Use 3rd column for y-data",
    "  smooth -T -m 2 -l auto data.txt       # Timestamp mode with Tikhonov",
    "  smooth -T -m 1 -n 5 -p 2 -d data.txt  # Timestamps with derivatives (dy/dt in sec)",
    "  smooth -g data.txt                    # Detailed grid analysis",
    "  cat data.txt | smooth -m 1 -n 5       # Use as Unix filter (stdin -> stdout)",
    "  smooth -m 2 -l 0.01 < input.txt       # Read from stdin with redirection",
    0
  };
  char **p = msg;
  usage();
  fprintf(stderr,"%s, Version %s, %s\n\n", progname, VERSION, REVDATE);
  while (*p)
    fprintf(stderr, "%s\n", *p++);
}
