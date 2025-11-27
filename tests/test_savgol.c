/*  Unit testy pro savgol.c
 *  Používá Unity testing framework
 *  V1.0/2025-11-27/ První verze testů - pouze uniformní mřížky
 */

#include "unity.h"           // Unity testing framework
#include "../savgol.h"       // Modul který testujeme
#include "../grid_analysis.h" // Pro analýzu mřížky
#include <math.h>            // Pro fabs() - absolutní hodnota
#include <stdlib.h>          // Pro malloc(), free(), rand(), srand()
#include <time.h>            // Pro time() - seed pro rand()


/* ============================================================================
 * TEST 1: Konstantní fce
 * ============================================================================
 */
void test_savgol_smooth_constant_data(void) {

  /* ARRANGE (příprava testovacích dat) */

  #define N 100 // Pocet bodu
  SavgolResult *result;
  GridAnalysis *grid_info;

  double x[N];
  double y[N];
  double d[N];
  for ( int i=0; i<N; i++ ) {
    x[i] = 0.1*i;
    y[i] = 3.14;
    d[i] = 0.0;
  }

  int P = 0; /* Stupen polynomu */
  int W = 5; /* Window size */

  /* Analyzuj mřížku (nutné pro savgol) */
  grid_info = analyze_grid(x, N, 0);
  TEST_ASSERT_NOT_NULL(grid_info);

  for (P=0; P<3; P++) {
    /* ACT (provedení testované funkce)
     * Parametry:
     *   x - pole x-ových souřadnic
     *   y - pole y-novych souradnic
     *   N - počet bodů
     *   W - window size
     *   P - stupen polynomu
     *   grid_info - informace o mřížce
     */
    result =  savgol_smooth(x, y, N, W, P, grid_info);

    // 1. Zkontroluj že funkce vrátila validní pointer (ne NULL)
    TEST_ASSERT_NOT_NULL(result);

    // 2. Number of points
    TEST_ASSERT_EQUAL_INT(N, result->n);

    // 3. Poly degree
    TEST_ASSERT_EQUAL_INT(P, result->poly_degree);

    // 4. Window size
    TEST_ASSERT_EQUAL_INT(W, result->window_size);

    //5. Test smoothed data
    TEST_ASSERT_DOUBLE_ARRAY_WITHIN(0.001, y, result->y_smooth, N);

    //6. Test first derivatives
    TEST_ASSERT_DOUBLE_ARRAY_WITHIN(0.001, d, result->y_deriv, N);

    /* CLEANUP */
    free_savgol_result(result);
  }

  free_grid_analysis(grid_info);
  #undef N
}

/* ============================================================================
 * TEST 2: Linearni fce
 * ============================================================================
 */
void test_savgol_smooth_linear_data(void) {

  /* ARRANGE (příprava testovacích dat) */

  #define N 100 // Pocet bodu
  SavgolResult *result;
  GridAnalysis *grid_info;

  double x[N];
  double y[N];
  double d[N];
  for ( int i=0; i<N; i++ ) {
    x[i] = 0.1*i;
    y[i] = 3.14 + 0.2*x[i];
    d[i] = 0.2;
  }

  int P = 0; /* Stupen polynomu */
  int W = 5; /* Window size */

  /* Analyzuj mřížku */
  grid_info = analyze_grid(x, N, 0);
  TEST_ASSERT_NOT_NULL(grid_info);

  for (P=1; P<4; P++) {
    /* ACT */
    result =  savgol_smooth(x, y, N, W, P, grid_info);

    // 1. Zkontroluj že funkce vrátila validní pointer (ne NULL)
    TEST_ASSERT_NOT_NULL(result);

    // 2. Number of points
    TEST_ASSERT_EQUAL_INT(N, result->n);

    // 3. Poly degree
    TEST_ASSERT_EQUAL_INT(P, result->poly_degree);

    // 4. Window size
    TEST_ASSERT_EQUAL_INT(W, result->window_size);

    //5. Test smoothed data
    TEST_ASSERT_DOUBLE_ARRAY_WITHIN(0.001, y, result->y_smooth, N);

    //6. Test first derivatives
    TEST_ASSERT_DOUBLE_ARRAY_WITHIN(0.001, d, result->y_deriv, N);

    /* CLEANUP */
    free_savgol_result(result);
  }

  free_grid_analysis(grid_info);
  #undef N
}

/* ============================================================================
 * TEST 3: Quadratic fce
 * ============================================================================
 */
void test_savgol_smooth_quadratic_data(void) {

  /* ARRANGE (příprava testovacích dat) */

  #define N 100 // Pocet bodu
  SavgolResult *result;
  GridAnalysis *grid_info;

  double x[N];
  double y[N];
  double d[N];
  for ( int i=0; i<N; i++ ) {
    x[i] = 0.1*i;
    y[i] = 3.14 + 0.2*x[i] + 0.7*x[i]*x[i];
    d[i] = 0.2 + 1.4*x[i];
  }

  int P = 0; /* Stupen polynomu */
  int W = 5; /* Window size */

  /* Analyzuj mřížku */
  grid_info = analyze_grid(x, N, 0);
  TEST_ASSERT_NOT_NULL(grid_info);

  for (P=2; P<5; P++) {
    /* ACT */
    result =  savgol_smooth(x, y, N, W, P, grid_info);

    // 1. Zkontroluj že funkce vrátila validní pointer (ne NULL)
    TEST_ASSERT_NOT_NULL(result);

    // 2. Number of points
    TEST_ASSERT_EQUAL_INT(N, result->n);

    // 3. Poly degree
    TEST_ASSERT_EQUAL_INT(P, result->poly_degree);

    // 4. Window size
    TEST_ASSERT_EQUAL_INT(W, result->window_size);

    //5. Test smoothed data
    TEST_ASSERT_DOUBLE_ARRAY_WITHIN(0.001, y, result->y_smooth, N);

    //6. Test first derivatives
    TEST_ASSERT_DOUBLE_ARRAY_WITHIN(0.001, d, result->y_deriv, N);

    /* CLEANUP */
    free_savgol_result(result);
  }

  free_grid_analysis(grid_info);
  #undef N
}


/* ============================================================================
 * TEST 4: Edge case - Minimální počet bodů (N == W)
 * ============================================================================
 */
void test_savgol_edge_case_n_equals_w(void) {
  /* ARRANGE */
  #define N_MIN 5
  SavgolResult *result;
  GridAnalysis *grid_info;

  double x[N_MIN];
  double y[N_MIN];
  for (int i = 0; i < N_MIN; i++) {
    x[i] = i;
    y[i] = 2.0 + 0.5 * i;  // Lineární funkce
  }

  int W = N_MIN;  // Okno stejně velké jako počet bodů
  int P = 1;      // Lineární fit

  grid_info = analyze_grid(x, N_MIN, 0);
  TEST_ASSERT_NOT_NULL(grid_info);

  /* ACT */
  result = savgol_smooth(x, y, N_MIN, W, P, grid_info);

  /* ASSERT */
  TEST_ASSERT_NOT_NULL(result);
  TEST_ASSERT_EQUAL_INT(N_MIN, result->n);
  TEST_ASSERT_EQUAL_INT(P, result->poly_degree);
  TEST_ASSERT_EQUAL_INT(W, result->window_size);

  /* CLEANUP */
  free_savgol_result(result);
  free_grid_analysis(grid_info);
  #undef N_MIN
}

/* ============================================================================
 * TEST 5: Edge case - N < W (mělo by vrátit NULL nebo upravit W)
 * ============================================================================
 */
void test_savgol_edge_case_n_less_than_w(void) {
  /* ARRANGE */
  #define N_SMALL 5
  SavgolResult *result;
  GridAnalysis *grid_info;

  double x[N_SMALL];
  double y[N_SMALL];
  for (int i = 0; i < N_SMALL; i++) {
    x[i] = i;
    y[i] = 2.0;
  }

  int W = 7;  // Okno větší než počet bodů
  int P = 0;

  grid_info = analyze_grid(x, N_SMALL, 0);
  TEST_ASSERT_NOT_NULL(grid_info);

  /* ACT */
  result = savgol_smooth(x, y, N_SMALL, W, P, grid_info);

  /* ASSERT */
  // Funkce by měla vrátit NULL
  TEST_ASSERT_NULL(result);

  /* CLEANUP */
  free_grid_analysis(grid_info);
  #undef N_SMALL
}

/* ============================================================================
 * TEST 6: Edge case - Velký počet bodů
 * ============================================================================
 */
void test_savgol_edge_case_large_n(void) {
  /* ARRANGE */
  #define N_LARGE 1000
  SavgolResult *result;
  GridAnalysis *grid_info;

  double *x = malloc(N_LARGE * sizeof(double));
  double *y = malloc(N_LARGE * sizeof(double));

  TEST_ASSERT_NOT_NULL(x);
  TEST_ASSERT_NOT_NULL(y);

  for (int i = 0; i < N_LARGE; i++) {
    x[i] = i * 0.01;
    y[i] = 5.0 + 3.0 * x[i];  // Lineární funkce
  }

  int W = 11;
  int P = 1;

  grid_info = analyze_grid(x, N_LARGE, 0);
  TEST_ASSERT_NOT_NULL(grid_info);

  /* ACT */
  result = savgol_smooth(x, y, N_LARGE, W, P, grid_info);

  /* ASSERT */
  TEST_ASSERT_NOT_NULL(result);
  TEST_ASSERT_EQUAL_INT(N_LARGE, result->n);
  TEST_ASSERT_EQUAL_INT(P, result->poly_degree);
  TEST_ASSERT_EQUAL_INT(W, result->window_size);

  /* CLEANUP */
  free_savgol_result(result);
  free_grid_analysis(grid_info);
  free(x);
  free(y);
  #undef N_LARGE
}

/* ============================================================================
 * TEST 7: Edge case - Vysoký stupeň polynomu (P = 6)
 * ============================================================================
 */
void test_savgol_edge_case_high_poly_degree(void) {
  /* ARRANGE */
  #define N_TEST 20
  SavgolResult *result;
  GridAnalysis *grid_info;

  double x[N_TEST];
  double y[N_TEST];
  for (int i = 0; i < N_TEST; i++) {
    x[i] = i;
    y[i] = 1.0 + 0.5*i;  // Jednoduchá lineární funkce
  }

  int W = 15;  // Musí být > P
  int P = 6;   // Maximální povolený stupeň

  grid_info = analyze_grid(x, N_TEST, 0);
  TEST_ASSERT_NOT_NULL(grid_info);

  /* ACT */
  result = savgol_smooth(x, y, N_TEST, W, P, grid_info);

  /* ASSERT */
  TEST_ASSERT_NOT_NULL(result);
  TEST_ASSERT_EQUAL_INT(N_TEST, result->n);
  TEST_ASSERT_EQUAL_INT(P, result->poly_degree);
  TEST_ASSERT_EQUAL_INT(W, result->window_size);

  /* CLEANUP */
  free_savgol_result(result);
  free_grid_analysis(grid_info);
  #undef N_TEST
}

/* ============================================================================
 * TEST 8: Edge case - Velké okno (W blízké N)
 * ============================================================================
 */
void test_savgol_edge_case_large_window(void) {
  /* ARRANGE */
  #define N_TEST 21
  SavgolResult *result;
  GridAnalysis *grid_info;

  double x[N_TEST];
  double y[N_TEST];
  for (int i = 0; i < N_TEST; i++) {
    x[i] = i;
    y[i] = 10.0 + 2.0*i + 0.3*i*i;  // Kvadratická funkce
  }

  int W = 19;  // Velké okno (musí být liché, blízké N)
  int P = 2;

  grid_info = analyze_grid(x, N_TEST, 0);
  TEST_ASSERT_NOT_NULL(grid_info);

  /* ACT */
  result = savgol_smooth(x, y, N_TEST, W, P, grid_info);

  /* ASSERT */
  TEST_ASSERT_NOT_NULL(result);
  TEST_ASSERT_EQUAL_INT(N_TEST, result->n);
  TEST_ASSERT_EQUAL_INT(W, result->window_size);

  // Pro velké okno a kvadratickou funkci by měl být výsledek velmi blízký vstupu
  TEST_ASSERT_DOUBLE_ARRAY_WITHIN(0.01, y, result->y_smooth, N_TEST);

  /* CLEANUP */
  free_savgol_result(result);
  free_grid_analysis(grid_info);
  #undef N_TEST
}

/* ============================================================================
 * TEST 9: Edge case - Minimální velikost okna (W = 3)
 * ============================================================================
 */
void test_savgol_edge_case_min_window(void) {
  /* ARRANGE */
  #define N_TEST 16
  SavgolResult *result;
  GridAnalysis *grid_info;

  double x[N_TEST];
  double y[N_TEST];
  for (int i = 0; i < N_TEST; i++) {
    x[i] = i;
    y[i] = 5.0 + 1.5*i;  // Lineární funkce
  }

  int W = 3;   // Minimální povolené okno
  int P = 1;   // Lineární fit

  grid_info = analyze_grid(x, N_TEST, 0);
  TEST_ASSERT_NOT_NULL(grid_info);

  /* ACT */
  result = savgol_smooth(x, y, N_TEST, W, P, grid_info);

  /* ASSERT */
  TEST_ASSERT_NOT_NULL(result);
  TEST_ASSERT_EQUAL_INT(W, result->window_size);
  TEST_ASSERT_EQUAL_INT(P, result->poly_degree);

  /* CLEANUP */
  free_savgol_result(result);
  free_grid_analysis(grid_info);
  #undef N_TEST
}

/* ============================================================================
 * TEST 10: Edge case - Nevalidní parametry (NULL pointery)
 * ============================================================================
 */
void test_savgol_edge_case_null_pointers(void) {
  /* ARRANGE */
  #define N_TEST 10
  double x[N_TEST];
  double y[N_TEST];
  GridAnalysis *grid_info;

  for (int i = 0; i < N_TEST; i++) {
    x[i] = i;
    y[i] = i;
  }

  grid_info = analyze_grid(x, N_TEST, 0);
  TEST_ASSERT_NOT_NULL(grid_info);

  SavgolResult *result;

  /* ACT & ASSERT - NULL pointer pro x */
  result = savgol_smooth(NULL, y, N_TEST, 5, 1, grid_info);
  TEST_ASSERT_NULL(result);

  /* ACT & ASSERT - NULL pointer pro y */
  result = savgol_smooth(x, NULL, N_TEST, 5, 1, grid_info);
  TEST_ASSERT_NULL(result);

  /* ACT & ASSERT - NULL pointer pro grid_info */
  result = savgol_smooth(x, y, N_TEST, 5, 1, NULL);
  TEST_ASSERT_NULL(result);

  /* CLEANUP */
  free_grid_analysis(grid_info);
  #undef N_TEST
}

/* ============================================================================
 * TEST 11: Edge case - Sudá velikost okna (mělo by vrátit NULL)
 * ============================================================================
 */
void test_savgol_edge_case_even_window(void) {
  /* ARRANGE */
  #define N_TEST 16
  double x[N_TEST];
  double y[N_TEST];
  GridAnalysis *grid_info;

  for (int i = 0; i < N_TEST; i++) {
    x[i] = i;
    y[i] = 3.0;
  }

  SavgolResult *result;
  int W = 6;  // Sudé číslo (nepovolené)
  int P = 1;

  grid_info = analyze_grid(x, N_TEST, 0);
  TEST_ASSERT_NOT_NULL(grid_info);

  /* ACT */
  result = savgol_smooth(x, y, N_TEST, W, P, grid_info);

  /* ASSERT */
  // Funkce by měla vrátit NULL pro sudé W
  TEST_ASSERT_NULL(result);

  /* CLEANUP */
  free_grid_analysis(grid_info);
  #undef N_TEST
}

/* ============================================================================
 * TEST 12: Edge case - P >= W (stupeň polynomu >= velikost okna)
 * ============================================================================
 */
void test_savgol_edge_case_p_greater_or_equal_w(void) {
  /* ARRANGE */
  #define N_TEST 20
  double x[N_TEST];
  double y[N_TEST];
  GridAnalysis *grid_info;

  for (int i = 0; i < N_TEST; i++) {
    x[i] = i;
    y[i] = 2.0 + 0.3*i;
  }

  SavgolResult *result;
  int W = 5;
  int P = 5;  // P >= W (nedostatek bodů pro fit)

  grid_info = analyze_grid(x, N_TEST, 0);
  TEST_ASSERT_NOT_NULL(grid_info);

  /* ACT */
  result = savgol_smooth(x, y, N_TEST, W, P, grid_info);

  /* ASSERT */
  // Toto by mělo vrátit NULL
  TEST_ASSERT_NULL(result);

  /* CLEANUP */
  free_grid_analysis(grid_info);
  #undef N_TEST
}

/* ============================================================================
 * TEST 13: Konstantní funkce + Gaussovský šum
 * ============================================================================
 */
void test_savgol_constant_with_noise(void) {
  /* ARRANGE */
  #define N_NOISE 100
  SavgolResult *result;
  GridAnalysis *grid_info;

  double x[N_NOISE];
  double y[N_NOISE];
  double y_noisy[N_NOISE];

  // Konstantní funkce y = 5.0
  for (int i = 0; i < N_NOISE; i++) {
    x[i] = i * 0.1;
    y[i] = 5.0;
  }

  // Přidej šum
  srand(12345);  // Fixní seed
  double noise_amplitude = 0.3;
  for (int i = 0; i < N_NOISE; i++) {
    double noise = noise_amplitude * ((double)rand() / RAND_MAX * 2.0 - 1.0);
    y_noisy[i] = y[i] + noise;
  }

  int W = 11;
  int P = 0;

  grid_info = analyze_grid(x, N_NOISE, 0);
  TEST_ASSERT_NOT_NULL(grid_info);

  /* ACT */
  result = savgol_smooth(x, y_noisy, N_NOISE, W, P, grid_info);

  /* ASSERT */
  TEST_ASSERT_NOT_NULL(result);
  TEST_ASSERT_EQUAL_INT(N_NOISE, result->n);
  TEST_ASSERT_EQUAL_INT(P, result->poly_degree);

  // Vyhlazená data by měla být blíže původní konstantě
  double error_noisy = 0.0;
  double error_smooth = 0.0;
  for (int i = W/2; i < N_NOISE - W/2; i++) {
    error_noisy += fabs(y_noisy[i] - y[i]);
    error_smooth += fabs(result->y_smooth[i] - y[i]);
  }

  TEST_ASSERT_LESS_THAN_DOUBLE(error_noisy / 2.0, error_smooth);

  // Vyhlazená data by měla být blízko 5.0
  for (int i = W/2; i < N_NOISE - W/2; i++) {
    TEST_ASSERT_DOUBLE_WITHIN(0.2, 5.0, result->y_smooth[i]);
  }

  /* CLEANUP */
  free_savgol_result(result);
  free_grid_analysis(grid_info);
  #undef N_NOISE
}

/* ============================================================================
 * TEST 14: Lineární funkce + Gaussovský šum
 * ============================================================================
 */
void test_savgol_linear_with_noise(void) {
  /* ARRANGE */
  #define N_NOISE 100
  SavgolResult *result;
  GridAnalysis *grid_info;

  double x[N_NOISE];
  double y[N_NOISE];
  double y_noisy[N_NOISE];

  // Lineární funkce y = 2.0 + 0.5*x
  for (int i = 0; i < N_NOISE; i++) {
    x[i] = i * 0.1;
    y[i] = 2.0 + 0.5 * x[i];
  }

  // Přidej šum
  srand(23456);
  double noise_amplitude = 0.4;
  for (int i = 0; i < N_NOISE; i++) {
    double noise = noise_amplitude * ((double)rand() / RAND_MAX * 2.0 - 1.0);
    y_noisy[i] = y[i] + noise;
  }

  int W = 11;
  int P = 1;

  grid_info = analyze_grid(x, N_NOISE, 0);
  TEST_ASSERT_NOT_NULL(grid_info);

  /* ACT */
  result = savgol_smooth(x, y_noisy, N_NOISE, W, P, grid_info);

  /* ASSERT */
  TEST_ASSERT_NOT_NULL(result);
  TEST_ASSERT_EQUAL_INT(N_NOISE, result->n);
  TEST_ASSERT_EQUAL_INT(P, result->poly_degree);

  // Vyhlazená data by měla být blíže původní přímce
  double error_noisy = 0.0;
  double error_smooth = 0.0;
  for (int i = W/2; i < N_NOISE - W/2; i++) {
    error_noisy += fabs(y_noisy[i] - y[i]);
    error_smooth += fabs(result->y_smooth[i] - y[i]);
  }

  TEST_ASSERT_LESS_THAN_DOUBLE(error_noisy / 2.0, error_smooth);

  // Vyhlazená data by měla být blízko původní funkce
  for (int i = W/2; i < N_NOISE - W/2; i++) {
    TEST_ASSERT_DOUBLE_WITHIN(0.3, y[i], result->y_smooth[i]);
  }

  /* CLEANUP */
  free_savgol_result(result);
  free_grid_analysis(grid_info);
  #undef N_NOISE
}

/* ============================================================================
 * TEST 15: Kvadratická funkce + šum
 * ============================================================================
 */
void test_savgol_quadratic_with_noise(void) {
  /* ARRANGE */
  #define N_NOISE 100
  SavgolResult *result;
  GridAnalysis *grid_info;

  double x[N_NOISE];
  double y[N_NOISE];
  double y_noisy[N_NOISE];
  double dy[N_NOISE];

  // Kvadratická funkce y = 1.0 + 0.2*x + 0.3*x²
  for (int i = 0; i < N_NOISE; i++) {
    x[i] = i * 0.1;
    y[i] = 1.0 + 0.2 * x[i] + 0.3 * x[i] * x[i];
    dy[i] = 0.2 + 0.6 * x[i];
  }

  // Přidej šum
  srand(34567);
  double noise_amplitude = 0.5;
  for (int i = 0; i < N_NOISE; i++) {
    double noise = noise_amplitude * ((double)rand() / RAND_MAX * 2.0 - 1.0);
    y_noisy[i] = y[i] + noise;
  }

  int W = 13;
  int P = 2;

  grid_info = analyze_grid(x, N_NOISE, 0);
  TEST_ASSERT_NOT_NULL(grid_info);

  /* ACT */
  result = savgol_smooth(x, y_noisy, N_NOISE, W, P, grid_info);

  /* ASSERT */
  TEST_ASSERT_NOT_NULL(result);
  TEST_ASSERT_EQUAL_INT(N_NOISE, result->n);
  TEST_ASSERT_EQUAL_INT(P, result->poly_degree);

  // Vyhlazená data by měla být blíže původní funkce
  double error_noisy = 0.0;
  double error_smooth = 0.0;
  for (int i = W/2; i < N_NOISE - W/2; i++) {
    error_noisy += fabs(y_noisy[i] - y[i]);
    error_smooth += fabs(result->y_smooth[i] - y[i]);
  }

  TEST_ASSERT_LESS_THAN(error_noisy, error_smooth * 2.0);

  // Vyhlazená data by měla být blízko původní funkce
  for (int i = W/2; i < N_NOISE - W/2; i++) {
    TEST_ASSERT_DOUBLE_WITHIN(0.7*noise_amplitude, y[i], result->y_smooth[i]);
  }

  /* CLEANUP */
  free_savgol_result(result);
  free_grid_analysis(grid_info);
  #undef N_NOISE
}

/* ============================================================================
 * TEST 16: Test odmítnutí neuniformní mřížky
 * ============================================================================
 * Savgol MUSÍ odmítnout neuniformní mřížku (CV > 0.05)
 */
void test_savgol_rejects_nonuniform_grid(void) {
  /* ARRANGE */
  #define N_NONUNIF 50
  SavgolResult *result;
  GridAnalysis *grid_info;

  double x[N_NONUNIF];
  double y[N_NONUNIF];

  // Vytvoř neuniformní mřížku s výraznou variabilitou
  x[0] = 0.0;
  for (int i = 1; i < N_NONUNIF; i++) {
    // Střídání hustých a řídkých oblastí
    if (i % 2 == 0) {
      x[i] = x[i-1] + 0.05;  // Husté
    } else {
      x[i] = x[i-1] + 0.15;  // Řídké
    }
  }

  // Konstantní y hodnoty
  for (int i = 0; i < N_NONUNIF; i++) {
    y[i] = 5.0;
  }

  int W = 7;
  int P = 1;

  grid_info = analyze_grid(x, N_NONUNIF, 0);
  TEST_ASSERT_NOT_NULL(grid_info);

  // Ověř že mřížka je opravdu neuniformní (CV > 0.05)
  TEST_ASSERT_GREATER_THAN_DOUBLE(0.05, grid_info->cv);

  /* ACT */
  result = savgol_smooth(x, y, N_NONUNIF, W, P, grid_info);

  /* ASSERT */
  // Savgol MUSÍ odmítnout neuniformní mřížku - vrací NULL
  TEST_ASSERT_NULL(result);

  /* CLEANUP */
  free_grid_analysis(grid_info);
  #undef N_NONUNIF
}

/* ============================================================================
 * POZNÁMKY K TESTŮM:
 * ============================================================================
 *
 * 1. STRUKTURA KAŽDÉHO TESTU (AAA pattern):
 *    - ARRANGE: Připrav testovací data
 *    - ACT:     Zavolej testovanou funkci
 *    - ASSERT:  Zkontroluj výsledky
 *    - CLEANUP: Uvolni paměť
 *
 * 2. DŮLEŽITÉ PRO SAVGOL:
 *    - VŽDY je nutné zavolat analyze_grid() před savgol_smooth()
 *    - Savgol vyžaduje uniformní mřížku (CV < 0.05)
 *    - Musíme uvolnit jak SavgolResult, tak GridAnalysis
 *
 * 3. ROZDÍL OD POLYFIT TESTŮ:
 *    - Přidán parametr grid_info do všech volání
 *    - Přidán test na odmítnutí neuniformní mřížky (test 16)
 *    - Vynechány testy s neuniformní mřížkou (savgol je nepodporuje)
 *
 * 4. MEMORY MANAGEMENT:
 *    Je KRITICKÉ volat:
 *    - free_savgol_result(result)
 *    - free_grid_analysis(grid_info)
 *    Jinak dostaneme memory leak.
 *
 * 5. NUMERICKÁ TOLERANCE:
 *    Pro double používáme vždy WITHIN s malou tolerancí (0.001).
 *    Nikdy neporovnáváme floaty pomocí == !
 *
 */
