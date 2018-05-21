/*  Vyrovnani tabulkovych bodu polynomem 2 stupně, pocitany nejmensimi
 *  ctverci z 5-ti okolnich bodu  */

#include <stdio.h>
#include <stdlib.h>
//#include <string.h>
//#include <errno.h>
#include <libgen.h>  /* char *basename(char *), char *dirname(char *) */
#include <unistd.h>  /* UNIX standard function definitions */

#include "decomment.h"
#include "revision.h"
#include "matrix.h"
#include "create_matrix.h"
#include "gauss_elimination.h"
#include "free_matrix.h"

#define N 5  /* Odd number > 3, data for least square */
#define DPMAX 6  /* Max degree of approximation polynom */
#define DP 2  /* Default degree of approx. polynom */

/* Lokal declare functions */
static void usage(void);
static void help(void);
static double power(double x, int n);

/* Global variables */
char *progname;

int main(int argc, char **argv)
{
  char *filename;
  FILE *fp;
  double *x = NULL;  /* Abscisa */
  double *y = NULL;  /* Ordinata */
  int n=0; /* Pocet tabulkovych bodu */
  int i,j,k;
  int ch;
  double sx[DPMAX+3], sy[DPMAX+1];  /* Sums */
  MATRIX C; /* Matrix of coefficients */
  MATRIX B; /* Right side vector */
  int sp = N; /* Sliding points for least square (default) */
  int o;      /* Offset for central of moving windows */
  int dp = DP; /* Degree of approximation polynom */
  double fi, dfi;

  progname = basename(argv[0]);

/* None arguments */  
  if (argc == 1) {
    usage();
    exit (EXIT_FAILURE);
  }

/* Options command line */
  while ( (ch = getopt(argc, argv, "n:p:h?")) != -1 ) {
    switch (ch) {
      case 'n':
        sp = atoi(optarg);
        break;
      case 'p':
        dp = atoi(optarg);
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

  if ( sp<5 || !(sp%2) ) {
    fprintf(stderr,"Incorrect number points in moving windows (odd > 3)!\n");
    exit (EXIT_FAILURE);
  }

    o = sp/2;  /* Offset */

  if (dp < 0 || dp > DPMAX) {
    fprintf(stderr,"Incorrect degree of approx. polynom, (0<=p<=%d)!\n",DPMAX);
    exit (EXIT_FAILURE);
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
      x = (double *)realloc((void *)x,(n+1)*sizeof(double));
      y = (double *)realloc((void *)y,(n+1)*sizeof(double));
      x[n] = rx, y[n] = ry; 
      n++;
      if (x == NULL || y == NULL) {
        fprintf(stderr,"No memory for data table!\n");
        exit(EXIT_FAILURE);
      }
    }
    fclose(fp);
  }

  if (n < sp) {
    fprintf(stderr,"Need moore data (n < %d)!\n",sp);
    exit (EXIT_FAILURE);
  }

  C = create_matrix(dp+1,dp+1);
  B = create_matrix(dp+1,1);

/* Approximation points (x[i]) */ 
  printf("# Data smooth - aprox. pol. %ddg from %d points of moving windows\
 (least square)\n",dp,sp);
  printf("#    x          y          y'\n");

  for (i=o; i<n-o; i++) {

    for (k=0; k<dp+3; k++)
      sx[k] = 0;
    for (k=0; k<dp+1; k++)
      sy[k] = 0;

  /* Sums for matrix off coefficients */  
    for (j=i-o; j<=i+o; j++) {
      double dx = x[j]-x[i];
     
      for (k=0; k<dp+3; k++)
        sx[k] += power(dx,k);

      for (k=0; k<dp+1; k++)
        sy[k] += y[j]*power(dx,k);
    }


    {
      int i,j;
      for (i=0; i<=dp; i++)
        for (j=0; j<=dp; j++)
          C.M[i][j] = sx[i+j];  /* Matrix of coefficients */

      for (i=0; i<=dp; i++)
        B.M[i][0] = sy[i];  /* Right side vector */
    }

    gauss_elimination(C,B);  /* Solve linear equations */ 
        
/* Left boundary */
    if (i==o) {
      double dx[o];
      int m;
      for (k=0;k<o;k++) {
        dx[k] = x[k]-x[o];
        fi = dfi = 0;
        for (m=0; m<=dp; m++)
          fi += B.M[m][0]*power(dx[k],m);
        for (m=1; m<=dp; m++)
          dfi += m*B.M[m][0]*power(dx[k],m-1);

        printf("%10.6lG %10.6lG %10.6lG\n",x[k],fi,dfi);
      }
    }

/* Print table */
    fi = B.M[0][0];
    if (dp>0)
      dfi = B.M[1][0];
    else
      dfi = 0;
    printf("%10.6lG %10.6lG %10.6lG\n", x[i],fi,dfi); 

/* Right boundary */
    if (i == n-o-1) {
      double dx[o];
      int m;
      for (k=0;k<o;k++) {
        dx[k] = x[n-o+k]-x[n-o-1];
        fi = dfi = 0;
        for (m=0; m<=dp; m++)
          fi += B.M[m][0]*power(dx[k],m);
        for (m=1; m<=dp; m++)
          dfi += m*B.M[m][0]*power(dx[k],m-1);

        printf("%10.6lG %10.6lG %10.6lG\n",x[n-o+k],fi,dfi);
      }
    }
  }

  free_matrix(C);
  free_matrix(B);

  return EXIT_SUCCESS;
}

/* Usage */
static void usage(void)
{
  fprintf(stderr,"Data smooth by approximation polynom from moving windows of data\
 (least square)\n");
  fprintf(stderr, "Usage: %s [options] data_file\n" ,progname);
  fprintf(stderr,"-h or -? for help\n");
}

/* Help */
static void help(void)
{
  static char *msg[] = {"-h, -?\tHelp",
                 "-n\tPoints in moving windows, default 5",  
                 "-p\tDegree of approx. polynom, default 2",    
	           0};
  char **p = msg;
  usage();
  fprintf(stderr,"%s, Version %s, %s\n\n", progname, VERSION, REVDATE);
  while (*p)
    fprintf(stderr, "%s\n", *p++);
}

/* Power x^n */
static double power(double x, int n)
{
  double p = 1;
  
  while (n-- > 0)
    p *= x;

  return p;
}

