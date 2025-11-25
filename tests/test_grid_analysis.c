/*  Unit testy pro grid_analysis.c
 *  Používá Unity testing framework
 *  V1.0/2025-11-25/ První verze testů
 */

#include "unity.h"           // Unity testing framework
#include "../grid_analysis.h" // Modul který testujeme
#include <math.h>            // Pro fabs() - absolutní hodnota

/* ============================================================================
 * SETUP A TEARDOWN FUNKCE
 * ============================================================================
 * Tyto funkce se volají automaticky před a po každém testu.
 * setUp()    - příprava před testem (alokace paměti, inicializace)
 * tearDown() - úklid po testu (uvolnění paměti, zavření souborů)
 */

void setUp(void) {
    /* Tato funkce se spustí PŘED KAŽDÝM testem
     * Zde bychom mohli:
     * - alokovat paměť pro testovací data
     * - otevřít soubory
     * - inicializovat globální proměnné
     *
     * Pro naše jednoduché testy není potřeba nic dělat.
     */
}

void tearDown(void) {
    /* Tato funkce se spustí PO KAŽDÉM testu
     * Zde bychom měli:
     * - uvolnit alokovanou paměť (free)
     * - zavřít otevřené soubory
     * - vyčistit globální stav
     *
     * Pro naše jednoduché testy není potřeba nic dělat.
     */
}

/* ============================================================================
 * TEST 1: Uniformní grid (ideální případ)
 * ============================================================================
 * Testujeme, zda funkce správně rozpozná perfektně uniformní grid.
 */
void test_grid_perfectly_uniform(void) {
    /* ARRANGE (příprava testovacích dat)
     * Vytvoříme perfektně uniformní grid s konstantním spacingem = 1.0
     */
    double x[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    int n = 6;

    /* ACT (provedení testované funkce)
     * Voláme analyze_grid() která analyzuje uniformitu gridu.
     * Parametry:
     *   x - pole x-ových souřadnic
     *   n - počet bodů
     *   0 - neukládáme spacing array (šetříme paměť)
     */
    GridAnalysis *result = analyze_grid(x, n, 0);

    /* ASSERT (kontrola výsledků)
     * Používáme Unity makra pro kontrolu:
     * TEST_ASSERT_NOT_NULL - kontrola že pointer není NULL
     * TEST_ASSERT_EQUAL_INT - kontrola že int hodnoty jsou stejné
     * TEST_ASSERT_DOUBLE_WITHIN - kontrola že double hodnoty jsou blízko
     */

    // 1. Zkontroluj že funkce vrátila validní pointer (ne NULL)
    TEST_ASSERT_NOT_NULL(result);

    // 2. Grid by měl být označen jako uniformní
    TEST_ASSERT_EQUAL_INT(1, result->is_uniform);

    // 3. Koeficient variace (CV) by měl být téměř 0
    //    CV = std_dev / avg - pro perfektní uniformitu je CV = 0
    //    TEST_ASSERT_DOUBLE_WITHIN(tolerance, expected, actual)
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 0.0, result->cv);

    // 4. Průměrný spacing by měl být 1.0
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 1.0, result->h_avg);

    // 5. Min spacing = max spacing = 1.0
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 1.0, result->h_min);
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 1.0, result->h_max);

    // 6. Ratio max/min by měl být 1.0 (všechny spacingy stejné)
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 1.0, result->ratio_max_min);

    /* CLEANUP (uvolnění paměti)
     * VELMI DŮLEŽITÉ v C! Musíme uvolnit alokovanou paměť.
     * GridAnalysis struktura byla vytvořena pomocí malloc() uvnitř
     * analyze_grid(), proto ji musíme uvolnit pomocí free_grid_analysis()
     */
    free_grid_analysis(result);
}

/* ============================================================================
 * TEST 2: Neuniformní grid (velká variace)
 * ============================================================================
 * Testujeme detekci výrazně neuniformního gridu.
 */
void test_grid_nonuniform(void) {
    /* ARRANGE
     * Vytvoříme grid s různými spacingy:
     * 0.0 -> 1.0 : spacing = 1.0
     * 1.0 -> 1.5 : spacing = 0.5
     * 1.5 -> 4.0 : spacing = 2.5  <- výrazně větší!
     * 4.0 -> 5.0 : spacing = 1.0
     * 5.0 -> 10.0: spacing = 5.0  <- velmi velký!
     */
    double x[] = {0.0, 1.0, 1.5, 4.0, 5.0, 10.0};
    int n = 6;

    /* ACT */
    GridAnalysis *result = analyze_grid(x, n, 0);

    /* ASSERT */
    // 1. Funkce by měla úspěšně proběhnout
    TEST_ASSERT_NOT_NULL(result);

    // 2. Grid by NEMĚL být označen jako uniformní
    TEST_ASSERT_EQUAL_INT(0, result->is_uniform);

    // 3. Koeficient variace by měl být VYSOKÝ (>> 0.05)
    //    Pro neuniformní grid očekáváme CV > 0.05
    //    Unity syntax: TEST_ASSERT_GREATER_THAN(threshold, actual) = actual > threshold
    TEST_ASSERT_GREATER_THAN_DOUBLE(0.05, result->cv);

    // 4. Min spacing by měl být 0.5
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 0.5, result->h_min);

    // 5. Max spacing by měl být 5.0
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 5.0, result->h_max);

    // 6. Ratio max/min by měl být 10.0 (5.0 / 0.5 = 10)
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 10.0, result->ratio_max_min);

    /* CLEANUP */
    free_grid_analysis(result);
}

/* ============================================================================
 * TEST 3: Minimální počet bodů (edge case)
 * ============================================================================
 * Testujeme chování s minimálním počtem bodů (n=2).
 * To je nejmenší možný dataset - máme pouze 1 interval.
 */
void test_grid_minimum_points(void) {
    /* ARRANGE
     * Pouze 2 body -> 1 interval -> nelze počítat variaci
     */
    double x[] = {0.0, 1.0};
    int n = 2;

    /* ACT */
    GridAnalysis *result = analyze_grid(x, n, 0);

    /* ASSERT */
    // Funkce by měla vždy vrátit validní výsledek
    TEST_ASSERT_NOT_NULL(result);

    // S pouze 1 intervalem je grid technicky "uniformní"
    // (nemáme s čím porovnávat)
    TEST_ASSERT_EQUAL_INT(1, result->is_uniform);

    // Spacing by měl být 1.0
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 1.0, result->h_avg);
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 1.0, result->h_min);
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 1.0, result->h_max);

    // CV by měl být 0 (žádná variace)
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 0.0, result->cv);

    /* CLEANUP */
    free_grid_analysis(result);
}

/* ============================================================================
 * TEST 4: Téměř uniformní grid (malá variace)
 * ============================================================================
 * Testujeme grid který je "téměř" uniformní - malé odchylky.
 * Tohle je realistický případ z měření.
 */
void test_grid_nearly_uniform(void) {
    /* ARRANGE
     * Grid s malými odchylkami od uniformity
     * Ideální spacing by byl 1.0, máme malé odchylky ±0.02
     */
    double x[] = {0.0, 0.98, 2.01, 2.99, 4.02, 5.0};
    int n = 6;

    /* ACT */
    GridAnalysis *result = analyze_grid(x, n, 0);

    /* ASSERT */
    TEST_ASSERT_NOT_NULL(result);

    // CV by měl být velmi malý (< 0.05)
    // Ale ne přesně 0 jako u perfektního gridu
    //    Unity syntax pro double: použij _DOUBLE verzi maker
    TEST_ASSERT_LESS_THAN_DOUBLE(0.05, result->cv);
    TEST_ASSERT_GREATER_THAN_DOUBLE(0.0, result->cv);

    // Průměrný spacing by měl být blízko 1.0
    TEST_ASSERT_DOUBLE_WITHIN(0.05, 1.0, result->h_avg);

    // Ratio by měl být blízko 1.0 (ale ne přesně 1.0)
    TEST_ASSERT_DOUBLE_WITHIN(0.1, 1.0, result->ratio_max_min);
    TEST_ASSERT_GREATER_THAN_DOUBLE(1.0, result->ratio_max_min);

    /* CLEANUP */
    free_grid_analysis(result);
}

/* ============================================================================
 * TEST 5: NULL pointer handling (error case)
 * ============================================================================
 * Testujeme, jak funkce reaguje na špatný vstup (NULL pointer).
 * Dobrý kód by měl gracefully zvládnout chybový vstup.
 */
void test_grid_null_pointer(void) {
    /* ARRANGE
     * Předáme NULL pointer místo pole dat
     */
    double *x = NULL;
    int n = 5;

    /* ACT */
    GridAnalysis *result = analyze_grid(x, n, 0);

    /* ASSERT
     * Funkce by měla vrátit NULL když dostane nevalidní vstup
     * (Poznámka: Tohle závisí na implementaci analyze_grid.
     *  Pokud funkce neřeší NULL check, test selže a upozorní nás!)
     */
    TEST_ASSERT_NULL(result);

    /* CLEANUP
     * Není potřeba - result je NULL
     */
}

/* ============================================================================
 * TEST 6: Velký dataset (performance check)
 * ============================================================================
 * Testujeme funkci s větším datasetem - kontrola že funguje i pro n=1000
 */
void test_grid_large_dataset(void) {
    /* ARRANGE
     * Vytvoříme velký uniformní grid
     */
    int n = 1000;
    double x[1000];

    // Generujeme uniformní grid: x[i] = i * 0.1
    for (int i = 0; i < n; i++) {
        x[i] = i * 0.1;
    }

    /* ACT */
    GridAnalysis *result = analyze_grid(x, n, 0);

    /* ASSERT */
    TEST_ASSERT_NOT_NULL(result);

    // Velký uniformní grid by měl být rozpoznán jako uniformní
    TEST_ASSERT_EQUAL_INT(1, result->is_uniform);

    // CV by měl být téměř 0
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 0.0, result->cv);

    // Spacing by měl být 0.1
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 0.1, result->h_avg);

    // Počet intervalů by měl být n-1 = 999
    TEST_ASSERT_EQUAL_INT(999, result->n_intervals);

    /* CLEANUP */
    free_grid_analysis(result);
}

/* ============================================================================
 * TEST 7: Grid s outlierem (realistický problém)
 * ============================================================================
 * Testujeme grid kde většina spacingů je uniformní, ale 1 je výrazně odlišný.
 * Tohle může nastat při chybě měření nebo missing data point.
 */
void test_grid_with_outlier(void) {
    /* ARRANGE
     * Grid kde většina spacingů je 1.0, ale jeden je 10.0
     */
    double x[] = {0.0, 1.0, 2.0, 3.0, 13.0, 14.0};  // Gap mezi 3.0 a 13.0!
    int n = 6;

    /* ACT */
    GridAnalysis *result = analyze_grid(x, n, 0);

    /* ASSERT */
    TEST_ASSERT_NOT_NULL(result);

    // Grid by NEMĚL být označen jako uniformní kvůli outlieru
    TEST_ASSERT_EQUAL_INT(0, result->is_uniform);

    // Max spacing by měl být 10.0 (outlier)
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 10.0, result->h_max);

    // Min spacing by měl být 1.0
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 1.0, result->h_min);

    // Ratio by měl být 10.0
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 10.0, result->ratio_max_min);

    // CV by měl být vysoký kvůli velké variaci
    TEST_ASSERT_GREATER_THAN_DOUBLE(0.5, result->cv);

    /* CLEANUP */
    free_grid_analysis(result);
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
 *    Je KRITICKÉ volat free_grid_analysis() na konci každého testu!
 *    Jinak dostaneme memory leak.
 *
 * 5. NUMERICKÁ TOLERANCE:
 *    Pro double používáme vždy WITHIN s malou tolerancí (0.001).
 *    Nikdy neporovnáváme floaty pomocí == !
 *
 * 6. CO TESTUJEME:
 *    ✓ Základní funkcionalita (uniformní grid)
 *    ✓ Edge cases (minimum bodů)
 *    ✓ Error handling (NULL pointer)
 *    ✓ Realistické případy (téměř uniformní, outlier)
 *    ✓ Performance (velký dataset)
 *
 * 7. CO NETESTUJEME (zatím):
 *    - Memory leaks (na to použijeme Valgrind)
 *    - Thread safety (single-threaded kód)
 *    - Performance benchmarking (jiný typ testů)
 */
