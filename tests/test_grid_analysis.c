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
    GridAnalysis *result = analyze_grid(x, n);

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
    GridAnalysis *result = analyze_grid(x, n);

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
    GridAnalysis *result = analyze_grid(x, n);

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
    GridAnalysis *result = analyze_grid(x, n);

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
    GridAnalysis *result = analyze_grid(x, n);

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
    GridAnalysis *result = analyze_grid(x, n);

    /* ASSERT */
    TEST_ASSERT_NOT_NULL(result);

    // Velký uniformní grid by měl být rozpoznán jako uniformní
    TEST_ASSERT_EQUAL_INT(1, result->is_uniform);

    // CV by měl být téměř 0
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 0.0, result->cv);

    // Spacing by měl být 0.1
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 0.1, result->h_avg);

    // Počet bodů by měl odpovídat vstupu
    TEST_ASSERT_EQUAL_INT(1000, result->n_points);

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
    GridAnalysis *result = analyze_grid(x, n);

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

/* TEST 8: Detekce chybějících vzorků v pravidelné mřížce
 *
 * A lost sample leaves a gap of k*h_base for integer k >= 2. Dropping index 10
 * leaves one gap of 2h (1 missing); dropping 20, 21, 22 leaves one gap of 4h
 * (3 missing). Counts are exact integers, so no tolerance is needed on them.
 */
void test_grid_dropouts_detected(void) {
    /* ARRANGE: uniform base 1.0, 40 nominal points, drop 10 and 20-22 */
    double x[36];
    int n = 0;
    for (int i = 0; i < 40; i++) {
        if (i == 10 || i == 20 || i == 21 || i == 22) continue;
        x[n++] = i * 1.0;
    }

    /* ACT */
    GridAnalysis *result = analyze_grid(x, n);

    /* ASSERT */
    TEST_ASSERT_NOT_NULL(result);
    TEST_ASSERT_EQUAL_INT(36, n);
    TEST_ASSERT_EQUAL_INT(1, result->has_dropouts);
    TEST_ASSERT_EQUAL_INT(4, result->n_missing);
    TEST_ASSERT_EQUAL_INT(2, result->n_gaps);
    TEST_ASSERT_EQUAL_INT(3, result->max_run);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 1.0, result->h_base);

    /* CLEANUP */
    free_grid_analysis(result);
}

/* TEST 9: Kompletní mřížka nesmí hlásit žádné výpadky */
void test_grid_dropouts_none_when_complete(void) {
    /* ARRANGE */
    double x[40];
    for (int i = 0; i < 40; i++) x[i] = i * 1.0;

    /* ACT */
    GridAnalysis *result = analyze_grid(x, 40);

    /* ASSERT */
    TEST_ASSERT_NOT_NULL(result);
    TEST_ASSERT_EQUAL_INT(0, result->has_dropouts);
    TEST_ASSERT_EQUAL_INT(0, result->n_missing);
    TEST_ASSERT_EQUAL_INT(0, result->n_gaps);
    /* Every spacing is exactly 1*h_base */
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 1.0, result->integer_fit);

    /* CLEANUP */
    free_grid_analysis(result);
}

/* TEST 10: Stupňovaná mřížka NESMÍ být hlášena jako mřížka s výpadky
 *
 * This is the discriminator test. On a geometrically graded mesh the spacings
 * vary continuously, so few of them land near an integer multiple of the median
 * and integer_fit falls below DROPOUT_FIT_MIN. Without this test the gate is
 * untested and the detector would happily invent missing samples in data that
 * is merely non-uniform.
 */
void test_grid_dropouts_rejects_graded_mesh(void) {
    /* ARRANGE: geometric grading r = 1.05 */
    double x[60];
    double h = 0.1;
    x[0] = 0.0;
    for (int i = 1; i < 60; i++) { x[i] = x[i-1] + h; h *= 1.05; }

    /* ACT */
    GridAnalysis *result = analyze_grid(x, 60);

    /* ASSERT */
    TEST_ASSERT_NOT_NULL(result);
    TEST_ASSERT_EQUAL_INT(0, result->has_dropouts);
    TEST_ASSERT_LESS_THAN_DOUBLE(0.90, result->integer_fit);

    /* CLEANUP */
    free_grid_analysis(result);
}

/* TEST 11: Medián musí přežít těžkou ztrátu vzorků
 *
 * Guards the median-vs-mean choice. Keeping only i % 5 >= 2 drops 40% of the
 * samples and leaves spacings of 1, 1, 3 repeating: the median is 1.0 (correct)
 * while the MEAN is ~1.67. A mean-based h_base gives integer_fit ~= 0.33, below
 * the gate, so has_dropouts would be 0 and this test fails.
 */
void test_grid_dropouts_median_survives_heavy_loss(void) {
    /* ARRANGE */
    double x[200];
    int n = 0;
    for (int i = 0; i < 200; i++) {
        if (i % 5 < 2) continue;
        x[n++] = i * 1.0;
    }

    /* ACT */
    GridAnalysis *result = analyze_grid(x, n);

    /* ASSERT */
    TEST_ASSERT_NOT_NULL(result);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 1.0, result->h_base);
    TEST_ASSERT_EQUAL_INT(1, result->has_dropouts);
    /* Mean spacing here is ~1.67 — well outside the tolerance above */
    TEST_ASSERT_GREATER_THAN_DOUBLE(1.5, result->h_avg);

    /* CLEANUP */
    free_grid_analysis(result);
}

/* TEST 12: Krátká mřížka — detekce se přeskočí, analýza běží dál
 *
 * Below DROPOUT_MIN_SPACES the median is not meaningful, so the dropout fields
 * stay at their calloc zeros. analyze_grid must still return a valid analysis
 * with cv computed — a diagnostic extra must never turn a working run into a
 * failure.
 */
void test_grid_dropouts_short_grid_guard(void) {
    /* ARRANGE: 6 points -> 5 spacings, below the 10-spacing floor */
    double x[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};

    /* ACT */
    GridAnalysis *result = analyze_grid(x, 6);

    /* ASSERT */
    TEST_ASSERT_NOT_NULL(result);
    TEST_ASSERT_EQUAL_INT(0, result->has_dropouts);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.0, result->h_base);
    TEST_ASSERT_EQUAL_INT(0, result->n_missing);
    /* The rest of the analysis is unaffected */
    TEST_ASSERT_EQUAL_INT(1, result->is_uniform);
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 1.0, result->h_avg);

    /* CLEANUP */
    free_grid_analysis(result);
}

/* TEST 13: Změna vzorkovací frekvence NENÍ ztráta vzorků
 *
 * This is the corrective test. Both halves of the record have spacings that are
 * integer multiples of the median, so integer_fit reaches 1.00 and the v5.11.51
 * gate passed: it reported "2241 missing samples, coverage 18.2%" for a record
 * where nothing is missing at all. regime_shift separates the two cases: it is
 * 1.0 for dropouts of any severity and 10.0 here.
 */
void test_grid_regime_rate_change_detected(void) {
    /* ARRANGE: 250 points at spacing 1.0, then 250 at 0.1 */
    double x[500];
    for (int i = 0; i < 250; i++) x[i] = i * 1.0;
    for (int i = 0; i < 250; i++) x[250 + i] = x[249] + (i + 1) * 0.1;

    /* ACT */
    GridAnalysis *result = analyze_grid(x, 500);

    /* ASSERT */
    TEST_ASSERT_NOT_NULL(result);
    TEST_ASSERT_EQUAL_INT(1, result->multi_regime);
    TEST_ASSERT_EQUAL_INT(0, result->has_dropouts);
    TEST_ASSERT_DOUBLE_WITHIN(0.1, 10.0, result->regime_shift);
    TEST_ASSERT_DOUBLE_WITHIN(0.1, 10.0, result->max_jump);

    /* CLEANUP */
    free_grid_analysis(result);
}

/* TEST 14: Stupňovaná mřížka NENÍ smíšený režim
 *
 * Pins the two-part multi_regime condition. A graded mesh has a low mode
 * fraction (spacings vary continuously, so few sit at the median) but no jump
 * at all. Testing regime_shift alone would let an `||` slip in and every
 * graded mesh would be misreported -- the exact mistake made while prototyping
 * this design.
 */
void test_grid_regime_graded_mesh_is_not_mixed(void) {
    /* ARRANGE: geometric grading r = 1.10 */
    double x[60];
    double h = 0.1;
    x[0] = 0.0;
    for (int i = 1; i < 60; i++) { x[i] = x[i-1] + h; h *= 1.10; }

    /* ACT */
    GridAnalysis *result = analyze_grid(x, 60);

    /* ASSERT */
    TEST_ASSERT_NOT_NULL(result);
    TEST_ASSERT_EQUAL_INT(0, result->multi_regime);
    /* The halves differ enormously -- a graded mesh grows without bound ... */
    TEST_ASSERT_GREATER_THAN_DOUBLE(1.5, result->regime_shift);
    /* ... yet no two NEIGHBOURING spacings differ by more than 1.1x, which is
     * why this is a gradual trend and not a regime change. Assert both: this
     * fixture is the reason multi_regime needs two conditions. */
    TEST_ASSERT_LESS_THAN_DOUBLE(2.0, result->max_jump);

    /* CLEANUP */
    free_grid_analysis(result);
}

/* TEST 15: Uniformní mřížka nemá žádný skok */
void test_grid_regime_uniform_has_no_jump(void) {
    /* ARRANGE */
    double x[40];
    for (int i = 0; i < 40; i++) x[i] = i * 1.0;

    /* ACT */
    GridAnalysis *result = analyze_grid(x, 40);

    /* ASSERT */
    TEST_ASSERT_NOT_NULL(result);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 1.0, result->max_jump);
    TEST_ASSERT_EQUAL_INT(0, result->n_jumps);
    TEST_ASSERT_EQUAL_INT(0, result->multi_regime);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 1.0, result->regime_shift);

    /* CLEANUP */
    free_grid_analysis(result);
}

/* TEST 16: Izolovaná mezera se najde A lokalizuje
 *
 * The case the removed cluster detector scored 0 on: it required
 * h_prev < 0.1*h_avg while h_prev here is 1.0 and 0.1*h_avg is ~0.55.
 *
 * max_jump_x is x[i+1] -- the point BETWEEN the two spacings compared, i.e.
 * the last sample before the gap. Here the compared spacings are h[18] = 1.0
 * and h[19] = 200.0, so the location is x[19] = 19.0.
 *
 * This fixture also sets has_dropouts (a 200x gap is an integer multiple, so
 * 199 "missing"). That is the instrument-restart case in miniature and is
 * correct; the assertions here are about the jump, not the classification.
 */
void test_grid_regime_isolated_gap_located(void) {
    /* ARRANGE: 20 points at 1.0, a 200.0 gap, then 20 more at 1.0 */
    double x[40];
    for (int i = 0; i < 20; i++) x[i] = i * 1.0;
    for (int i = 0; i < 20; i++) x[20 + i] = x[19] + 200.0 + i * 1.0;

    /* ACT */
    GridAnalysis *result = analyze_grid(x, 40);

    /* ASSERT */
    TEST_ASSERT_NOT_NULL(result);
    TEST_ASSERT_DOUBLE_WITHIN(1.0, 200.0, result->max_jump);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 19.0, result->max_jump_x);

    /* CLEANUP */
    free_grid_analysis(result);
}

/* TEST 17: Zpřísněná brána nesmí rozbít detekci výpadků z plánu 005
 *
 * Same fixture as test_grid_dropouts_detected. A regular grid with scattered
 * dropouts keeps a high mode fraction, so multi_regime must stay 0 and the
 * dropout report must survive unchanged.
 */
void test_grid_regime_dropouts_still_reported(void) {
    /* ARRANGE: uniform base 1.0, 40 nominal points, drop 10 and 20-22 */
    double x[36];
    int n = 0;
    for (int i = 0; i < 40; i++) {
        if (i == 10 || i == 20 || i == 21 || i == 22) continue;
        x[n++] = i * 1.0;
    }

    /* ACT */
    GridAnalysis *result = analyze_grid(x, n);

    /* ASSERT */
    TEST_ASSERT_NOT_NULL(result);
    TEST_ASSERT_EQUAL_INT(1, result->has_dropouts);
    TEST_ASSERT_EQUAL_INT(0, result->multi_regime);
    TEST_ASSERT_EQUAL_INT(4, result->n_missing);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 1.0, result->regime_shift);

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
