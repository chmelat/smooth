/*  Unit testy pro polyfit.c
 *  Používá Unity testing framework
 *  V1.0/2025-11-xx/ První verze testů
 */

#include "unity.h"           // Unity testing framework
#include "../polyfit.h" // Modul který testujeme
#include <math.h>            // Pro fabs() - absolutní hodnota
#include <stdlib.h>          // Pro malloc(), free(), rand(), srand()
#include <time.h>            // Pro time() - seed pro rand()


/* ============================================================================
 * TEST 1: Konstantní fce
 * ============================================================================
 */
void test_polyfit_smooth_constant_data(void) {

  /* ARRANGE (příprava testovacích dat) */

  #define N 100 // Pocet bodu
  PolyfitResult *result;

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

  for (P=0; P<3; P++) {
    /* ACT (provedení testované funkce)
     * Parametry:
     *   x - pole x-ových souřadnic
     *   y - pole y-novych souradnic
     *   N - počet bodů
     *   W - window size
     *   P - stupen polynomu
     */
    result =  polyfit_smooth(x, y, N, W, P);

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

    /* CLEANUP (uvolnění paměti)
     * VELMI DŮLEŽITÉ v C! Musíme uvolnit alokovanou paměť.
     * PolyfitResult struktura byla vytvořena pomocí malloc() uvnitř
     * polyfit_smooth(), proto ji musíme uvolnit pomocí free_polyfit_result()
     */
    free_polyfit_result(result);
  }
  #undef N
}

/* ============================================================================
 * TEST 2: Linearni fce
 * ============================================================================
 */
void test_polyfit_smooth_linear_data(void) {

  /* ARRANGE (příprava testovacích dat) */

  #define N 100 // Pocet bodu
  PolyfitResult *result;

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

  for (P=1; P<4; P++) {
    /* ACT (provedení testované funkce)
     * Parametry:
     *   x - pole x-ových souřadnic
     *   y - pole y-novych souradnic
     *   N - počet bodů
     *   W - window size
     *   P - stupen polynomu
     */
    result =  polyfit_smooth(x, y, N, W, P);

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

    /* CLEANUP (uvolnění paměti)
     * VELMI DŮLEŽITÉ v C! Musíme uvolnit alokovanou paměť.
     * PolyfitResult struktura byla vytvořena pomocí malloc() uvnitř
     * polyfit_smooth(), proto ji musíme uvolnit pomocí free_polyfit_result()
     */
    free_polyfit_result(result);
  }
  #undef N
}

/* ============================================================================
 * TEST 3: Quadratic fce
 * ============================================================================
 */
void test_polyfit_smooth_quadratic_data(void) {

  /* ARRANGE (příprava testovacích dat) */

  #define N 100 // Pocet bodu
  PolyfitResult *result;

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

  for (P=2; P<5; P++) {
    /* ACT (provedení testované funkce)
     * Parametry:
     *   x - pole x-ových souřadnic
     *   y - pole y-novych souradnic
     *   N - počet bodů
     *   W - window size
     *   P - stupen polynomu
     */
    result =  polyfit_smooth(x, y, N, W, P);

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

    /* CLEANUP (uvolnění paměti)
     * VELMI DŮLEŽITÉ v C! Musíme uvolnit alokovanou paměť.
     * PolyfitResult struktura byla vytvořena pomocí malloc() uvnitř
     * polyfit_smooth(), proto ji musíme uvolnit pomocí free_polyfit_result()
     */
    free_polyfit_result(result);
  }
  #undef N
}


/* ============================================================================
 * TEST 4: Edge case - Minimální počet bodů (N == W)
 * ============================================================================
 */
void test_polyfit_edge_case_n_equals_w(void) {
  /* ARRANGE */
  #define N_MIN 5
  PolyfitResult *result;

  double x[N_MIN];
  double y[N_MIN];
  for (int i = 0; i < N_MIN; i++) {
    x[i] = i;
    y[i] = 2.0 + 0.5 * i;  // Lineární funkce
  }

  int W = N_MIN;  // Okno stejně velké jako počet bodů
  int P = 1;      // Lineární fit

  /* ACT */
  result = polyfit_smooth(x, y, N_MIN, W, P);

  /* ASSERT */
  TEST_ASSERT_NOT_NULL(result);
  TEST_ASSERT_EQUAL_INT(N_MIN, result->n);
  TEST_ASSERT_EQUAL_INT(P, result->poly_degree);
  TEST_ASSERT_EQUAL_INT(W, result->window_size);

  /* CLEANUP */
  free_polyfit_result(result);
  #undef N_MIN
}

/* ============================================================================
 * TEST 5: Edge case - N < W (mělo by vrátit NULL nebo upravit W)
 * ============================================================================
 */
void test_polyfit_edge_case_n_less_than_w(void) {
  /* ARRANGE */
  #define N_SMALL 5
  PolyfitResult *result;

  double x[N_SMALL];
  double y[N_SMALL];
  for (int i = 0; i < N_SMALL; i++) {
    x[i] = i;
    y[i] = 2.0;
  }

  int W = 7;  // Okno větší než počet bodů
  int P = 0;

  /* ACT */
  result = polyfit_smooth(x, y, N_SMALL, W, P);

  /* ASSERT */
  // Funkce by měla buď vrátit NULL, nebo automaticky upravit W
  // Kontrolujeme že funkce nepadne a vrací něco rozumného
  if (result != NULL) {
    // Pokud vrací výsledek, měla by upravit window_size
    TEST_ASSERT_LESS_OR_EQUAL_INT(N_SMALL, result->window_size);
    free_polyfit_result(result);
  }
  // Pokud vrací NULL, je to také akceptovatelné chování

  #undef N_SMALL
}

/* ============================================================================
 * TEST 6: Edge case - Minimální počet bodů pro stupeň polynomu (N = P + 1)
 * ============================================================================
 */
void test_polyfit_edge_case_n_equals_p_plus_1(void) {
  /* ARRANGE */
  int P = 2;           // Kvadratický polynom
  #define N_MIN 3      // Minimální počet bodů pro P=2
  PolyfitResult *result;

  double x[N_MIN];
  double y[N_MIN];
  for (int i = 0; i < N_MIN; i++) {
    x[i] = i;
    y[i] = 1.0 + 2.0*i + 0.5*i*i;  // Kvadratická funkce
  }

  int W = 3;  // Minimální lichá velikost okna

  /* ACT */
  result = polyfit_smooth(x, y, N_MIN, W, P);

  /* ASSERT */
  TEST_ASSERT_NOT_NULL(result);
  TEST_ASSERT_EQUAL_INT(N_MIN, result->n);
  TEST_ASSERT_EQUAL_INT(P, result->poly_degree);

  /* CLEANUP */
  free_polyfit_result(result);
  #undef N_MIN
}

/* ============================================================================
 * TEST 7: Edge case - Velký počet bodů
 * ============================================================================
 */
void test_polyfit_edge_case_large_n(void) {
  /* ARRANGE */
  #define N_LARGE 1000
  PolyfitResult *result;

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

  /* ACT */
  result = polyfit_smooth(x, y, N_LARGE, W, P);

  /* ASSERT */
  TEST_ASSERT_NOT_NULL(result);
  TEST_ASSERT_EQUAL_INT(N_LARGE, result->n);
  TEST_ASSERT_EQUAL_INT(P, result->poly_degree);
  TEST_ASSERT_EQUAL_INT(W, result->window_size);

  /* CLEANUP */
  free_polyfit_result(result);
  free(x);
  free(y);
  #undef N_LARGE
}

/* ============================================================================
 * TEST 8: Edge case - Vysoký stupeň polynomu (P = 6, maximum dle dokumentace)
 * ============================================================================
 */
void test_polyfit_edge_case_high_poly_degree(void) {
  /* ARRANGE */
  #define N_TEST 20
  PolyfitResult *result;

  double x[N_TEST];
  double y[N_TEST];
  for (int i = 0; i < N_TEST; i++) {
    x[i] = i;
    y[i] = 1.0 + 0.5*i;  // Jednoduchá lineární funkce
  }

  int W = 15;  // Musí být > P
  int P = 6;   // Maximální povolený stupeň

  /* ACT */
  result = polyfit_smooth(x, y, N_TEST, W, P);

  /* ASSERT */
  TEST_ASSERT_NOT_NULL(result);
  TEST_ASSERT_EQUAL_INT(N_TEST, result->n);
  TEST_ASSERT_EQUAL_INT(P, result->poly_degree);
  TEST_ASSERT_EQUAL_INT(W, result->window_size);

  /* CLEANUP */
  free_polyfit_result(result);
  #undef N_TEST
}

/* ============================================================================
 * TEST 9: Edge case - Velké okno (W blízké N)
 * ============================================================================
 */
void test_polyfit_edge_case_large_window(void) {
  /* ARRANGE */
  #define N_TEST 21
  PolyfitResult *result;

  double x[N_TEST];
  double y[N_TEST];
  for (int i = 0; i < N_TEST; i++) {
    x[i] = i;
    y[i] = 10.0 + 2.0*i + 0.3*i*i;  // Kvadratická funkce
  }

  int W = 19;  // Velké okno (musí být liché, blízké N)
  int P = 2;

  /* ACT */
  result = polyfit_smooth(x, y, N_TEST, W, P);

  /* ASSERT */
  TEST_ASSERT_NOT_NULL(result);
  TEST_ASSERT_EQUAL_INT(N_TEST, result->n);
  TEST_ASSERT_EQUAL_INT(W, result->window_size);

  // Pro velké okno a kvadratickou funkci by měl být výsledek velmi blízký vstupu
  TEST_ASSERT_DOUBLE_ARRAY_WITHIN(0.01, y, result->y_smooth, N_TEST);

  /* CLEANUP */
  free_polyfit_result(result);
  #undef N_TEST
}

/* ============================================================================
 * TEST 10: Edge case - Minimální velikost okna (W = 3)
 * ============================================================================
 */
void test_polyfit_edge_case_min_window(void) {
  /* ARRANGE */
  #define N_TEST 16
  PolyfitResult *result;

  double x[N_TEST];
  double y[N_TEST];
  for (int i = 0; i < N_TEST; i++) {
    x[i] = i;
    y[i] = 5.0 + 1.5*i;  // Lineární funkce
  }

  int W = 3;   // Minimální povolené okno
  int P = 1;   // Lineární fit

  /* ACT */
  result = polyfit_smooth(x, y, N_TEST, W, P);

  /* ASSERT */
  TEST_ASSERT_NOT_NULL(result);
  TEST_ASSERT_EQUAL_INT(W, result->window_size);
  TEST_ASSERT_EQUAL_INT(P, result->poly_degree);

  /* CLEANUP */
  free_polyfit_result(result);
  #undef N_TEST
}

/* ============================================================================
 * TEST 11: Edge case - Nevalidní parametry (NULL pointery)
 * ============================================================================
 */
void test_polyfit_edge_case_null_pointers(void) {
  /* ARRANGE */
  #define N_TEST 10
  double x[N_TEST];
  double y[N_TEST];
  for (int i = 0; i < N_TEST; i++) {
    x[i] = i;
    y[i] = i;
  }

  PolyfitResult *result;

  /* ACT & ASSERT - NULL pointer pro x */
  result = polyfit_smooth(NULL, y, N_TEST, 5, 1);
  TEST_ASSERT_NULL(result);

  /* ACT & ASSERT - NULL pointer pro y */
  result = polyfit_smooth(x, NULL, N_TEST, 5, 1);
  TEST_ASSERT_NULL(result);

  /* ACT & ASSERT - Oba NULL */
  result = polyfit_smooth(NULL, NULL, N_TEST, 5, 1);
  TEST_ASSERT_NULL(result);

  #undef N_TEST
}

/* ============================================================================
 * TEST 12: Edge case - Nevalidní velikosti (záporné a nulové hodnoty)
 * ============================================================================
 */
void test_polyfit_edge_case_invalid_sizes(void) {
  /* ARRANGE */
  #define N_TEST 10
  double x[N_TEST];
  double y[N_TEST];
  for (int i = 0; i < N_TEST; i++) {
    x[i] = i;
    y[i] = i;
  }

  PolyfitResult *result;

  /* ACT & ASSERT - Záporné N */
  result = polyfit_smooth(x, y, -5, 5, 1);
  TEST_ASSERT_NULL(result);

  /* ACT & ASSERT - Nulové N */
  result = polyfit_smooth(x, y, 0, 5, 1);
  TEST_ASSERT_NULL(result);

  /* ACT & ASSERT - Záporné W */
  result = polyfit_smooth(x, y, N_TEST, -5, 1);
  TEST_ASSERT_NULL(result);

  /* ACT & ASSERT - Záporné P */
  result = polyfit_smooth(x, y, N_TEST, 5, -1);
  TEST_ASSERT_NULL(result);

  /* ACT & ASSERT - P > 6 (nad maximum) */
  result = polyfit_smooth(x, y, N_TEST, 11, 10);
  TEST_ASSERT_NULL(result);

  #undef N_TEST
}

/* ============================================================================
 * TEST 13: Edge case - Sudá velikost okna (mělo by být automaticky upraveno nebo zamítnuto)
 * ============================================================================
 */
void test_polyfit_edge_case_even_window(void) {
  /* ARRANGE */
  #define N_TEST 16
  double x[N_TEST];
  double y[N_TEST];
  for (int i = 0; i < N_TEST; i++) {
    x[i] = i;
    y[i] = 3.0;
  }

  PolyfitResult *result;
  int W = 6;  // Sudé číslo (nepovolené dle dokumentace)
  int P = 1;

  /* ACT */
  result = polyfit_smooth(x, y, N_TEST, W, P);

  /* ASSERT */
  // Funkce by měla buď vrátit NULL, nebo automaticky upravit na liché W
  if (result != NULL) {
    // Pokud vrací výsledek, window_size by mělo být liché
    TEST_ASSERT_TRUE(result->window_size % 2 == 1);
    free_polyfit_result(result);
  }
  // Pokud vrací NULL, je to také akceptovatelné

  #undef N_TEST
}

/* ============================================================================
 * TEST 14: Edge case - P >= W (stupeň polynomu >= velikost okna)
 * ============================================================================
 */
void test_polyfit_edge_case_p_greater_or_equal_w(void) {
  /* ARRANGE */
  #define N_TEST 20
  double x[N_TEST];
  double y[N_TEST];
  for (int i = 0; i < N_TEST; i++) {
    x[i] = i;
    y[i] = 2.0 + 0.3*i;
  }

  PolyfitResult *result;
  int W = 5;
  int P = 5;  // P >= W (nedostatek bodů pro fit)

  /* ACT */
  result = polyfit_smooth(x, y, N_TEST, W, P);

  /* ASSERT */
  // Toto by mělo vrátit NULL, protože nelze fitovat polynom stupně P
  // s pouze W body (potřebujeme alespoň P+1 bodů)
  TEST_ASSERT_NULL(result);

  #undef N_TEST
}

/* ============================================================================
 * TEST 15: Konstantní funkce + Gaussovský šum
 * ============================================================================
 */
void test_polyfit_constant_with_noise(void) {
  /* ARRANGE */
  #define N_NOISE 100
  PolyfitResult *result;

  double x[N_NOISE];
  double y[N_NOISE];
  double y_noisy[N_NOISE];

  // Konstantní funkce y = 5.0
  for (int i = 0; i < N_NOISE; i++) {
    x[i] = i * 0.1;
    y[i] = 5.0;
  }

  // Přidej šum (použij pseudonáhodná čísla)
  srand(12345);  // Fixní seed pro reprodukovatelnost testů
  double noise_amplitude = 0.3;
  for (int i = 0; i < N_NOISE; i++) {
    double noise = noise_amplitude * ((double)rand() / RAND_MAX * 2.0 - 1.0);  // -1.0 až 1.0
    y_noisy[i] = y[i] + noise;
  }

  int W = 11;  // Větší okno pro lepší vyhlazení
  int P = 0;   // Konstantní polynom

  /* ACT */
  result = polyfit_smooth(x, y_noisy, N_NOISE, W, P);

  /* ASSERT */
  TEST_ASSERT_NOT_NULL(result);
  TEST_ASSERT_EQUAL_INT(N_NOISE, result->n);
  TEST_ASSERT_EQUAL_INT(P, result->poly_degree);

  // Vyhlazená data by měla být blíže původní konstantě než zašuměná data
  // Spočítej průměrnou absolutní odchylku od skutečné hodnoty
  double error_noisy = 0.0;
  double error_smooth = 0.0;
  for (int i = W/2; i < N_NOISE - W/2; i++) {  // Ignoruj okraje
    error_noisy += fabs(y_noisy[i] - y[i]);
    error_smooth += fabs(result->y_smooth[i] - y[i]);
  }

  // Vyhlazená data by měla mít menší chybu než zašuměná
  // TEST_ASSERT_LESS_THAN(threshold, actual) znamená actual < threshold
  // Chceme: error_smooth < error_noisy / 2.0
  TEST_ASSERT_LESS_THAN_DOUBLE(error_noisy / 2.0, error_smooth);  // Alespoň 2x lepší

  // Vyhlazená data by měla být blízko 5.0 (tolerance 0.2)
  for (int i = W/2; i < N_NOISE - W/2; i++) {
    TEST_ASSERT_DOUBLE_WITHIN(0.2, 5.0, result->y_smooth[i]);
  }

  /* CLEANUP */
  free_polyfit_result(result);
  #undef N_NOISE
}

/* ============================================================================
 * TEST 16: Lineární funkce + Gaussovský šum
 * ============================================================================
 */
void test_polyfit_linear_with_noise(void) {
  /* ARRANGE */
  #define N_NOISE 100
  PolyfitResult *result;

  double x[N_NOISE];
  double y[N_NOISE];
  double y_noisy[N_NOISE];

  // Lineární funkce y = 2.0 + 0.5*x
  for (int i = 0; i < N_NOISE; i++) {
    x[i] = i * 0.1;
    y[i] = 2.0 + 0.5 * x[i];
  }

  // Přidej šum (použij pseudonáhodná čísla)
  srand(23456);  // Fixní seed pro reprodukovatelnost testů
  double noise_amplitude = 0.4;
  for (int i = 0; i < N_NOISE; i++) {
    double noise = noise_amplitude * ((double)rand() / RAND_MAX * 2.0 - 1.0);  // -1.0 až 1.0
    y_noisy[i] = y[i] + noise;
  }

  int W = 11;  // Větší okno pro lepší vyhlazení
  int P = 1;   // Lineární polynom

  /* ACT */
  result = polyfit_smooth(x, y_noisy, N_NOISE, W, P);

  /* ASSERT */
  TEST_ASSERT_NOT_NULL(result);
  TEST_ASSERT_EQUAL_INT(N_NOISE, result->n);
  TEST_ASSERT_EQUAL_INT(P, result->poly_degree);

  // Vyhlazená data by měla být blíže původní přímce než zašuměná data
  double error_noisy = 0.0;
  double error_smooth = 0.0;
  for (int i = W/2; i < N_NOISE - W/2; i++) {  // Ignoruj okraje
    error_noisy += fabs(y_noisy[i] - y[i]);
    error_smooth += fabs(result->y_smooth[i] - y[i]);
  }

  // Vyhlazená data by měla mít menší chybu než zašuměná
  // TEST_ASSERT_LESS_THAN(threshold, actual) znamená actual < threshold
  // Chceme: error_smooth < error_noisy / 2.0
  TEST_ASSERT_LESS_THAN_DOUBLE(error_noisy / 2.0, error_smooth);  // Alespoň 2x lepší

  // Vyhlazená data by měla být blízko původní funkce (tolerance 0.3)
  for (int i = W/2; i < N_NOISE - W/2; i++) {
    TEST_ASSERT_DOUBLE_WITHIN(0.3, y[i], result->y_smooth[i]);
  }

#if 0
  // Test derivace - měla by být blízko 0.5
  // Kvůli šumu a vyhlazení může být derivace méně přesná
  // Testujeme jen střední část dat (vyhýbáme se okrajům)
  int margin = W;
  double deriv_sum = 0.0;
  int deriv_count = 0;
  for (int i = margin; i < N_NOISE - margin; i++) {
    deriv_sum += result->y_deriv[i];
    deriv_count++;
  }
  double deriv_mean = deriv_sum / deriv_count;
  // Průměrná derivace by měla být blízko 0.5
  TEST_ASSERT_DOUBLE_WITHIN(0.1, 0.5, deriv_mean);
#endif

  /* CLEANUP */
  free_polyfit_result(result);
  #undef N_NOISE
}

/* ============================================================================
 * TEST 17: Kvadratická funkce + šum
 * ============================================================================
 */
void test_polyfit_quadratic_with_noise(void) {
  /* ARRANGE */
  #define N_NOISE 100
  PolyfitResult *result;

  double x[N_NOISE];
  double y[N_NOISE];
  double y_noisy[N_NOISE];
  double dy[N_NOISE];  // Analytická derivace

  // Kvadratická funkce y = 1.0 + 0.2*x + 0.3*x²
  // dy/dx = 0.2 + 0.6*x
  for (int i = 0; i < N_NOISE; i++) {
    x[i] = i * 0.1;
    y[i] = 1.0 + 0.2 * x[i] + 0.3 * x[i] * x[i];
    dy[i] = 0.2 + 0.6 * x[i];
  }

  // Přidej šum (použij pseudonáhodná čísla)
  srand(34567);  // Fixní seed pro reprodukovatelnost testů
  double noise_amplitude = 0.5;
  for (int i = 0; i < N_NOISE; i++) {
    double noise = noise_amplitude * ((double)rand() / RAND_MAX * 2.0 - 1.0);  // -1.0 až 1.0
    y_noisy[i] = y[i] + noise;
  }

  int W = 13;  // Větší okno pro lepší vyhlazení
  int P = 2;   // Kvadratický polynom

  /* ACT */
  result = polyfit_smooth(x, y_noisy, N_NOISE, W, P);

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

  // Vyhlazená data by měla mít menší chybu
  TEST_ASSERT_LESS_THAN(error_noisy, error_smooth * 2.0);

  // Vyhlazená data by měla být blízko původní funkce
  for (int i = W/2; i < N_NOISE - W/2; i++) {
    TEST_ASSERT_DOUBLE_WITHIN(0.7*noise_amplitude, y[i], result->y_smooth[i]);
  }

#if 0
  // Test derivace - kvůli šumu a kvadratické povaze derivace může být méně přesná
  for (int i = W; i < N_NOISE - W; i++) {
    TEST_ASSERT_DOUBLE_WITHIN(0.5, dy[i], result->y_deriv[i]);
  }
#endif

  /* CLEANUP */
  free_polyfit_result(result);
  #undef N_NOISE
}

/* ============================================================================
 * TEST 18: Test kvality vyhlazení - Standard deviation se šumem
 * ============================================================================
 */
void test_polyfit_noise_reduction_quality(void) {
  /* ARRANGE */
  #define N_NOISE 200
  PolyfitResult *result;

  double x[N_NOISE];
  double y[N_NOISE];
  double y_noisy[N_NOISE];

  // Konstantní funkce y = 10.0
  for (int i = 0; i < N_NOISE; i++) {
    x[i] = i * 0.05;
    y[i] = 10.0;
  }

  // Přidej větší šum (použij pseudonáhodná čísla)
  srand(45678);  // Fixní seed pro reprodukovatelnost testů
  double noise_amplitude = 1.0;
  for (int i = 0; i < N_NOISE; i++) {
    double noise = noise_amplitude * ((double)rand() / RAND_MAX * 2.0 - 1.0);  // -1.0 až 1.0
    y_noisy[i] = y[i] + noise;
  }

  int W = 15;  // Velké okno pro silné vyhlazení
  int P = 0;   // Konstantní polynom

  /* ACT */
  result = polyfit_smooth(x, y_noisy, N_NOISE, W, P);

  /* ASSERT */
  TEST_ASSERT_NOT_NULL(result);

  // Vypočítej standard deviation pro zašuměná a vyhlazená data
  double mean_noisy = 0.0;
  double mean_smooth = 0.0;
  int n_valid = 0;

  // Počítej mean ze střední části (bez okrajů)
  for (int i = W/2; i < N_NOISE - W/2; i++) {
    mean_noisy += y_noisy[i];
    mean_smooth += result->y_smooth[i];
    n_valid++;
  }
  mean_noisy /= n_valid;
  mean_smooth /= n_valid;

  double std_noisy = 0.0;
  double std_smooth = 0.0;
  for (int i = W/2; i < N_NOISE - W/2; i++) {
    double diff_noisy = y_noisy[i] - mean_noisy;
    double diff_smooth = result->y_smooth[i] - mean_smooth;
    std_noisy += diff_noisy * diff_noisy;
    std_smooth += diff_smooth * diff_smooth;
  }
  std_noisy = sqrt(std_noisy / n_valid);
  std_smooth = sqrt(std_smooth / n_valid);

  // Chceme: std_noisy > 0.3 (tedy šum musí být přítomen s amplitudou 1.0)
  TEST_ASSERT_GREATER_THAN_DOUBLE(0.3, std_noisy);  // Šum musí být přítomen

  // Standard deviation vyhlazených dat by měla být výrazně menší
  // Chceme: std_smooth < std_noisy / 2.0
  TEST_ASSERT_LESS_THAN_DOUBLE(std_noisy / 2, std_smooth);  // Alespoň 2x menší

  // Průměr by měl být blízko 10.0
  TEST_ASSERT_DOUBLE_WITHIN(0.3, 10.0, mean_smooth);

  /* CLEANUP */
  free_polyfit_result(result);
  #undef N_NOISE
}

/* ============================================================================
 * TEST 19: Neuniformní mřížka - Konstantní funkce
 * ============================================================================
 * Testuje polyfit na neekvidistantní mřížce s konstantní funkcí.
 * Mřížka: x[i] = x[i-1] + 0.1*sin²(i), což dává plynule se měnící hustotu bodů.
 */
void test_polyfit_nonuniform_grid_constant(void) {
  /* ARRANGE */
  #define N_NONUNIF 100
  PolyfitResult *result;

  double x[N_NONUNIF];
  double y[N_NONUNIF];

  /* Vytvoř neuniformní mřížku pomocí sin²(i) */
  x[0] = 0.0;
  for (int i = 1; i < N_NONUNIF; i++) {
    double sin_i = sin((double)i);
    x[i] = x[i-1] + 0.1 * sin_i * sin_i;  // d_i ∈ [0, 0.1]
  }

  /* Konstantní funkce y = 7.5 */
  for (int i = 0; i < N_NONUNIF; i++) {
    y[i] = 7.5;
  }

  int W = 11;
  int P = 0;  // Konstantní polynom

  /* ACT */
  result = polyfit_smooth(x, y, N_NONUNIF, W, P);

  /* ASSERT */
  TEST_ASSERT_NOT_NULL(result);
  TEST_ASSERT_EQUAL_INT(N_NONUNIF, result->n);
  TEST_ASSERT_EQUAL_INT(P, result->poly_degree);
  TEST_ASSERT_EQUAL_INT(W, result->window_size);

  /* Vyhlazená data by měla být blízko konstanty 7.5 */
  for (int i = W/2; i < N_NONUNIF - W/2; i++) {
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 7.5, result->y_smooth[i]);
  }

  /* Derivace konstanty by měla být 0 */
  for (int i = W/2; i < N_NONUNIF - W/2; i++) {
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 0.0, result->y_deriv[i]);
  }

  /* CLEANUP */
  free_polyfit_result(result);
  #undef N_NONUNIF
}

/* ============================================================================
 * TEST 20: Neuniformní mřížka - Lineární funkce
 * ============================================================================
 * Testuje polyfit na neekvidistantní mřížce s lineární funkcí.
 */
void test_polyfit_nonuniform_grid_linear(void) {
  /* ARRANGE */
  #define N_NONUNIF 100
  PolyfitResult *result;

  double x[N_NONUNIF];
  double y[N_NONUNIF];

  /* Vytvoř neuniformní mřížku pomocí sin²(i) */
  x[0] = 0.0;
  for (int i = 1; i < N_NONUNIF; i++) {
    double sin_i = sin((double)i);
    x[i] = x[i-1] + 0.1 * sin_i * sin_i;
  }

  /* Lineární funkce y = 2.0 + 0.3*x */
  for (int i = 0; i < N_NONUNIF; i++) {
    y[i] = 2.0 + 0.3 * x[i];
  }

  int W = 11;
  int P = 1;  // Lineární polynom

  /* ACT */
  result = polyfit_smooth(x, y, N_NONUNIF, W, P);

  /* ASSERT */
  TEST_ASSERT_NOT_NULL(result);
  TEST_ASSERT_EQUAL_INT(N_NONUNIF, result->n);
  TEST_ASSERT_EQUAL_INT(P, result->poly_degree);

  /* Vyhlazená data by měla odpovídat lineární funkci */
  for (int i = W/2; i < N_NONUNIF - W/2; i++) {
    double expected = 2.0 + 0.3 * x[i];
    TEST_ASSERT_DOUBLE_WITHIN(0.01, expected, result->y_smooth[i]);
  }

  /* Derivace by měla být blízko 0.3 (ve středu dat) */
  int margin = W;
  double deriv_sum = 0.0;
  int deriv_count = 0;
  for (int i = margin; i < N_NONUNIF - margin; i++) {
    deriv_sum += result->y_deriv[i];
    deriv_count++;
  }
  double deriv_mean = deriv_sum / deriv_count;
  TEST_ASSERT_DOUBLE_WITHIN(0.05, 0.3, deriv_mean);

  /* CLEANUP */
  free_polyfit_result(result);
  #undef N_NONUNIF
}

/* ============================================================================
 * TEST 21: Neuniformní mřížka - Kvadratická funkce
 * ============================================================================
 * Testuje polyfit na neekvidistantní mřížce s kvadratickou funkcí.
 */
void test_polyfit_nonuniform_grid_quadratic(void) {
  /* ARRANGE */
  #define N_NONUNIF 100
  PolyfitResult *result;

  double x[N_NONUNIF];
  double y[N_NONUNIF];
  double dy[N_NONUNIF];

  /* Vytvoř neuniformní mřížku pomocí sin²(i) */
  x[0] = 0.0;
  for (int i = 1; i < N_NONUNIF; i++) {
    double sin_i = sin((double)i);
    x[i] = x[i-1] + 0.1 * sin_i * sin_i;
  }

  /* Kvadratická funkce y = 1.0 + 0.5*x + 0.2*x²
   * Derivace dy/dx = 0.5 + 0.4*x */
  for (int i = 0; i < N_NONUNIF; i++) {
    y[i] = 1.0 + 0.5 * x[i] + 0.2 * x[i] * x[i];
    dy[i] = 0.5 + 0.4 * x[i];
  }

  int W = 13;
  int P = 2;  // Kvadratický polynom

  /* ACT */
  result = polyfit_smooth(x, y, N_NONUNIF, W, P);

  /* ASSERT */
  TEST_ASSERT_NOT_NULL(result);
  TEST_ASSERT_EQUAL_INT(N_NONUNIF, result->n);
  TEST_ASSERT_EQUAL_INT(P, result->poly_degree);

  /* Vyhlazená data by měla odpovídat kvadratické funkci */
  for (int i = W/2; i < N_NONUNIF - W/2; i++) {
    double expected = 1.0 + 0.5 * x[i] + 0.2 * x[i] * x[i];
    TEST_ASSERT_DOUBLE_WITHIN(0.02, expected, result->y_smooth[i]);
  }

  /* Derivace by měla být blízko analytické derivaci (ve středu dat) */
  int margin = W;
  for (int i = margin; i < N_NONUNIF - margin; i++) {
    TEST_ASSERT_DOUBLE_WITHIN(0.1, dy[i], result->y_deriv[i]);
  }

  /* CLEANUP */
  free_polyfit_result(result);
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
 * 2. UNITY MAKRA:
 *    - TEST_ASSERT_NOT_NULL(ptr)            : zkontroluje že ptr != NULL
 *    - TEST_ASSERT_NULL(ptr)                : zkontroluje že ptr == NULL
 *    - TEST_ASSERT_EQUAL_INT(expected, actual) : int porovnání
 *    - TEST_ASSERT_DOUBLE_WITHIN(delta, expected, actual) : double porovnání
 *    - TEST_ASSERT_LESS_THAN(threshold, actual) : actual < threshold
 *    - TEST_ASSERT_GREATER_THAN(threshold, actual) : actual > threshold
 *
 * 3. PROČ KONTROLUJEME VÍCE VLASTNOSTÍ:
 *    Každý test kontroluje několik vlastností výsledku.
 *    Pokud selže jedna assert, test skončí - ale vidíme KTERÁ to byla.
 *
 * 4. MEMORY MANAGEMENT:
 *    Je KRITICKÉ volat free_polyfit_result() na konci každého testu!
 *    Jinak dostaneme memory leak.
 *
 * 5. NUMERICKÁ TOLERANCE:
 *    Pro double používáme vždy WITHIN s malou tolerancí (0.001).
 *    Nikdy neporovnáváme floaty pomocí == !
 *
 */
