/*  Test runner pro smooth projekt
 *  Spouští všechny unit testy pomocí Unity frameworku
 *  V1.0/2025-11-25/ První verze
 */

#include "unity.h"
#include <stdio.h>

/* ============================================================================
 * DEKLARACE TESTOVACÍCH FUNKCÍ
 * ============================================================================
 * Musíme deklarovat všechny testovací funkce které chceme spustit.
 * Tyto funkce jsou definované v test_grid_analysis.c
 */

// Testy pro grid_analysis.c
void test_grid_perfectly_uniform(void);
void test_grid_nonuniform(void);
void test_grid_minimum_points(void);
void test_grid_nearly_uniform(void);
void test_grid_null_pointer(void);
void test_grid_large_dataset(void);
void test_grid_with_outlier(void);

// Testy pro polyfit.c
void test_polyfit_smooth_constant_data(void);
void test_polyfit_smooth_linear_data(void);
void test_polyfit_smooth_quadratic_data(void);

// Testy pro polyfit.c - Edge cases
void test_polyfit_edge_case_n_equals_w(void);
void test_polyfit_edge_case_n_less_than_w(void);
void test_polyfit_edge_case_n_equals_p_plus_1(void);
void test_polyfit_edge_case_large_n(void);
void test_polyfit_edge_case_high_poly_degree(void);
void test_polyfit_edge_case_large_window(void);
void test_polyfit_edge_case_min_window(void);
void test_polyfit_edge_case_null_pointers(void);
void test_polyfit_edge_case_invalid_sizes(void);
void test_polyfit_edge_case_even_window(void);
void test_polyfit_edge_case_p_greater_or_equal_w(void);

// Testy pro polyfit.c - Noise handling
void test_polyfit_constant_with_noise(void);
void test_polyfit_linear_with_noise(void);
void test_polyfit_quadratic_with_noise(void);
void test_polyfit_noise_reduction_quality(void);

// Testy pro polyfit.c - Non-uniform grid
void test_polyfit_nonuniform_grid_constant(void);
void test_polyfit_nonuniform_grid_linear(void);
void test_polyfit_nonuniform_grid_quadratic(void);

// Testy pro savgol.c - Basic functionality
void test_savgol_smooth_constant_data(void);
void test_savgol_smooth_linear_data(void);
void test_savgol_smooth_quadratic_data(void);

// Testy pro savgol.c - Edge cases
void test_savgol_edge_case_n_equals_w(void);
void test_savgol_edge_case_n_less_than_w(void);
void test_savgol_edge_case_large_n(void);
void test_savgol_edge_case_high_poly_degree(void);
void test_savgol_edge_case_large_window(void);
void test_savgol_edge_case_min_window(void);
void test_savgol_edge_case_null_pointers(void);
void test_savgol_edge_case_even_window(void);
void test_savgol_edge_case_p_greater_or_equal_w(void);

// Testy pro savgol.c - Noise handling
void test_savgol_constant_with_noise(void);
void test_savgol_linear_with_noise(void);
void test_savgol_quadratic_with_noise(void);

// Testy pro savgol.c - Grid uniformity
void test_savgol_rejects_nonuniform_grid(void);

// Testy pro tikhonov.c - Mathematical correctness
void test_tikhonov_constant_function(void);
void test_tikhonov_linear_function(void);
void test_tikhonov_quadratic_function(void);
void test_tikhonov_sine_function(void);
void test_tikhonov_with_noise_constant(void);
void test_tikhonov_with_noise_linear(void);
void test_tikhonov_lambda_effect(void);

// Testy pro tikhonov.c - Discretization methods
void test_tikhonov_uniform_grid_average_method(void);
void test_tikhonov_nonuniform_grid_local_method(void);
void test_tikhonov_discretization_threshold_cv015(void);
void test_tikhonov_grid_info_null_fallback(void);

// Testy pro tikhonov.c - GCV optimization
void test_gcv_optimal_lambda_constant_with_noise(void);
void test_gcv_optimal_lambda_quadratic_with_noise(void);
void test_gcv_trace_penalty_overfitting(void);

// Testy pro tikhonov.c - Edge cases
void test_tikhonov_minimum_points_n3(void);
void test_tikhonov_minimum_points_n2(void);
void test_tikhonov_lambda_zero(void);
void test_tikhonov_lambda_very_large(void);
void test_tikhonov_invalid_inputs(void);
void test_tikhonov_large_dataset(void);

// Testy pro tikhonov.c - Functional diagnostics
void test_tikhonov_functional_computation_consistency(void);
void test_tikhonov_functional_discretization_consistency(void);
void test_tikhonov_boundary_conditions_natural(void);

// Testy pro tikhonov.c - Memory management
void test_tikhonov_memory_allocation_success(void);
void test_tikhonov_memory_error_handling(void);

// Testy pro butterworth.c - Basic functionality
void test_butterworth_returns_valid_result(void);
void test_butterworth_smooths_noisy_signal(void);
void test_butterworth_preserves_constant_signal(void);
void test_butterworth_preserves_linear_trend(void);

// Testy pro butterworth.c - Cutoff frequency
void test_butterworth_higher_cutoff_less_smoothing(void);
void test_butterworth_invalid_cutoff_frequency(void);

// Testy pro butterworth.c - Zero-phase filtering
void test_butterworth_zero_phase_no_delay(void);

// Testy pro butterworth.c - Edge behavior
void test_butterworth_edge_points_no_excessive_oscillation(void);

// Testy pro butterworth.c - Grid uniformity
void test_butterworth_works_on_uniform_grid(void);
void test_butterworth_works_on_nearly_uniform_grid(void);
void test_butterworth_handles_moderately_nonuniform_grid(void);

// Testy pro butterworth.c - Edge cases and robustness
void test_butterworth_minimum_points(void);
void test_butterworth_too_few_points(void);
void test_butterworth_null_inputs(void);
void test_butterworth_large_dataset(void);

// Testy pro butterworth.c - Memory management
void test_butterworth_free_null_safe(void);
void test_butterworth_no_memory_leaks(void);

/* ============================================================================
 * MAIN FUNKCE - SPOUŠTÍ VŠECHNY TESTY
 * ============================================================================
 * Unity framework vyžaduje následující strukturu:
 * 1. UNITY_BEGIN()     - inicializace Unity
 * 2. RUN_TEST(...)     - spuštění každého testu
 * 3. UNITY_END()       - ukončení Unity a vrácení výsledku
 */

int main(void) {
    /* Vytiskni úvodní banner */
    printf("\n");
    printf("========================================\n");
    printf("  SMOOTH PROJECT - UNIT TESTS\n");
    printf("========================================\n");
    printf("Running tests for grid_analysis module\n");
    printf("========================================\n\n");

    /* UNITY_BEGIN() inicializuje Unity framework
     * Musí být volán PŘED jakýmkoliv RUN_TEST()
     */
    UNITY_BEGIN();

    /* ========================================
     * SPUŠTĚNÍ TESTŮ PRO GRID_ANALYSIS
     * ========================================
     * RUN_TEST(function_name) spustí jeden test.
     * Unity automaticky:
     * - Zavolá setUp() před testem
     * - Spustí test
     * - Zavolá tearDown() po testu
     * - Zaznamená výsledek (pass/fail)
     */

    printf("--- Basic functionality tests ---\n");
    RUN_TEST(test_grid_perfectly_uniform);
    RUN_TEST(test_grid_nonuniform);

    printf("\n--- Edge case tests ---\n");
    RUN_TEST(test_grid_minimum_points);
    RUN_TEST(test_grid_null_pointer);

    printf("\n--- Realistic scenario tests ---\n");
    RUN_TEST(test_grid_nearly_uniform);
    RUN_TEST(test_grid_with_outlier);

    printf("\n--- Performance tests ---\n");
    RUN_TEST(test_grid_large_dataset);

    printf("\n");
    printf("========================================\n");
    printf("Running tests for polyfit module\n");
    printf("========================================\n\n");
    printf("--- Functionality tests for constant data - polynoms dg. 0..2 ---\n");
    RUN_TEST(test_polyfit_smooth_constant_data);

    printf("\n--- Functionality tests for linear data - polynoms dg. 1..3 ---\n");
    RUN_TEST(test_polyfit_smooth_linear_data);

    printf("\n--- Functionality tests for quadratic data - polynoms dg. 2..4 ---\n");
    RUN_TEST(test_polyfit_smooth_quadratic_data);

    printf("\n--- Edge case tests for polyfit module ---\n");
    RUN_TEST(test_polyfit_edge_case_n_equals_w);
    RUN_TEST(test_polyfit_edge_case_n_less_than_w);
    RUN_TEST(test_polyfit_edge_case_n_equals_p_plus_1);
    RUN_TEST(test_polyfit_edge_case_large_n);
    RUN_TEST(test_polyfit_edge_case_high_poly_degree);
    RUN_TEST(test_polyfit_edge_case_large_window);
    RUN_TEST(test_polyfit_edge_case_min_window);
    RUN_TEST(test_polyfit_edge_case_null_pointers);
    RUN_TEST(test_polyfit_edge_case_invalid_sizes);
    RUN_TEST(test_polyfit_edge_case_even_window);
    RUN_TEST(test_polyfit_edge_case_p_greater_or_equal_w);

    printf("\n--- Noise handling tests for polyfit module ---\n");
    RUN_TEST(test_polyfit_constant_with_noise);
    RUN_TEST(test_polyfit_linear_with_noise);
    RUN_TEST(test_polyfit_quadratic_with_noise);
    RUN_TEST(test_polyfit_noise_reduction_quality);

    printf("\n--- Non-uniform grid tests for polyfit module ---\n");
    RUN_TEST(test_polyfit_nonuniform_grid_constant);
    RUN_TEST(test_polyfit_nonuniform_grid_linear);
    RUN_TEST(test_polyfit_nonuniform_grid_quadratic);

    printf("\n");
    printf("========================================\n");
    printf("Running tests for savgol module\n");
    printf("========================================\n\n");

    printf("--- Functionality tests for constant data - polynoms dg. 0..2 ---\n");
    RUN_TEST(test_savgol_smooth_constant_data);

    printf("\n--- Functionality tests for linear data - polynoms dg. 1..3 ---\n");
    RUN_TEST(test_savgol_smooth_linear_data);

    printf("\n--- Functionality tests for quadratic data - polynoms dg. 2..4 ---\n");
    RUN_TEST(test_savgol_smooth_quadratic_data);

    printf("\n--- Edge case tests for savgol module ---\n");
    RUN_TEST(test_savgol_edge_case_n_equals_w);
    RUN_TEST(test_savgol_edge_case_n_less_than_w);
    RUN_TEST(test_savgol_edge_case_large_n);
    RUN_TEST(test_savgol_edge_case_high_poly_degree);
    RUN_TEST(test_savgol_edge_case_large_window);
    RUN_TEST(test_savgol_edge_case_min_window);
    RUN_TEST(test_savgol_edge_case_null_pointers);
    RUN_TEST(test_savgol_edge_case_even_window);
    RUN_TEST(test_savgol_edge_case_p_greater_or_equal_w);

    printf("\n--- Noise handling tests for savgol module ---\n");
    RUN_TEST(test_savgol_constant_with_noise);
    RUN_TEST(test_savgol_linear_with_noise);
    RUN_TEST(test_savgol_quadratic_with_noise);

    printf("\n--- Grid uniformity tests for savgol module ---\n");
    RUN_TEST(test_savgol_rejects_nonuniform_grid);

    printf("\n");
    printf("========================================\n");
    printf("Running tests for tikhonov module\n");
    printf("========================================\n\n");

    printf("--- Mathematical correctness tests ---\n");
    RUN_TEST(test_tikhonov_constant_function);
    RUN_TEST(test_tikhonov_linear_function);
    RUN_TEST(test_tikhonov_quadratic_function);
    RUN_TEST(test_tikhonov_sine_function);
    RUN_TEST(test_tikhonov_with_noise_constant);
    RUN_TEST(test_tikhonov_with_noise_linear);
    RUN_TEST(test_tikhonov_lambda_effect);

    printf("\n--- Discretization method tests ---\n");
    RUN_TEST(test_tikhonov_uniform_grid_average_method);
    RUN_TEST(test_tikhonov_nonuniform_grid_local_method);
    RUN_TEST(test_tikhonov_discretization_threshold_cv015);
    RUN_TEST(test_tikhonov_grid_info_null_fallback);

    printf("\n--- GCV optimization tests ---\n");
    RUN_TEST(test_gcv_optimal_lambda_constant_with_noise);
    RUN_TEST(test_gcv_optimal_lambda_quadratic_with_noise);
    RUN_TEST(test_gcv_trace_penalty_overfitting);

    printf("\n--- Edge case tests ---\n");
    RUN_TEST(test_tikhonov_minimum_points_n3);
    RUN_TEST(test_tikhonov_minimum_points_n2);
    RUN_TEST(test_tikhonov_lambda_zero);
    RUN_TEST(test_tikhonov_lambda_very_large);
    RUN_TEST(test_tikhonov_invalid_inputs);
    RUN_TEST(test_tikhonov_large_dataset);

    printf("\n--- Functional diagnostics tests ---\n");
    RUN_TEST(test_tikhonov_functional_computation_consistency);
    RUN_TEST(test_tikhonov_functional_discretization_consistency);
    RUN_TEST(test_tikhonov_boundary_conditions_natural);

    printf("\n--- Memory management tests ---\n");
    RUN_TEST(test_tikhonov_memory_allocation_success);
    RUN_TEST(test_tikhonov_memory_error_handling);

    printf("\n");
    printf("========================================\n");
    printf("Running tests for butterworth module\n");
    printf("========================================\n\n");

    printf("--- Basic functionality tests ---\n");
    RUN_TEST(test_butterworth_returns_valid_result);
    RUN_TEST(test_butterworth_smooths_noisy_signal);
    RUN_TEST(test_butterworth_preserves_constant_signal);
    RUN_TEST(test_butterworth_preserves_linear_trend);

    printf("\n--- Cutoff frequency tests ---\n");
    RUN_TEST(test_butterworth_higher_cutoff_less_smoothing);
    RUN_TEST(test_butterworth_invalid_cutoff_frequency);

    printf("\n--- Zero-phase filtering tests ---\n");
    RUN_TEST(test_butterworth_zero_phase_no_delay);

    printf("\n--- Edge behavior tests ---\n");
    RUN_TEST(test_butterworth_edge_points_no_excessive_oscillation);

    printf("\n--- Grid uniformity tests ---\n");
    RUN_TEST(test_butterworth_works_on_uniform_grid);
    RUN_TEST(test_butterworth_works_on_nearly_uniform_grid);
    RUN_TEST(test_butterworth_handles_moderately_nonuniform_grid);

    printf("\n--- Edge cases and robustness tests ---\n");
    RUN_TEST(test_butterworth_minimum_points);
    RUN_TEST(test_butterworth_too_few_points);
    RUN_TEST(test_butterworth_null_inputs);
    RUN_TEST(test_butterworth_large_dataset);

    printf("\n--- Memory management tests ---\n");
    RUN_TEST(test_butterworth_free_null_safe);
    RUN_TEST(test_butterworth_no_memory_leaks);

    /* UNITY_END() ukončí Unity framework a vrátí výsledek
     * Návratová hodnota:
     *   0 = všechny testy prošly (success)
     *   1 = některý test selhal (failure)
     *
     * To je standardní Unix konvence: 0 = success, non-zero = error
     */
    return UNITY_END();
}

/* ============================================================================
 * JAK TO FUNGUJE:
 * ============================================================================
 *
 * 1. Když spustíme program, main() zavolá UNITY_BEGIN()
 *
 * 2. Pro každý RUN_TEST():
 *    a) Unity zavolá setUp() (z test_grid_analysis.c)
 *    b) Unity spustí testovací funkci (např. test_grid_perfectly_uniform)
 *    c) Pokud test projde, Unity zaznamená "PASS"
 *    d) Pokud test selže (některý assert), Unity zaznamená "FAIL"
 *    e) Unity zavolá tearDown() (z test_grid_analysis.c)
 *
 * 3. Po všech testech UNITY_END() vytiskne shrnutí:
 *    - Kolik testů proběhlo
 *    - Kolik prošlo
 *    - Kolik selhalo
 *    - Které konkrétně selhaly
 *
 * 4. Program vrátí exit code:
 *    - 0 pokud všechny testy prošly
 *    - 1 pokud některý test selhal
 *
 * PŘÍKLAD VÝSTUPU:
 * ========================================
 *   SMOOTH PROJECT - UNIT TESTS
 * ========================================
 * Running tests for grid_analysis module
 * ========================================
 *
 * --- Basic functionality tests ---
 * test_grid_analysis.c:44:test_grid_perfectly_uniform:PASS
 * test_grid_analysis.c:78:test_grid_nonuniform:PASS
 *
 * --- Edge case tests ---
 * test_grid_analysis.c:108:test_grid_minimum_points:PASS
 * test_grid_analysis.c:178:test_grid_null_pointer:FAIL: Expected NULL Was 0x12345
 *                                                         ^^^^^^^^^^^^^^^^^^^^^^^^^^^
 *                                                         Tento test selhal!
 *
 * -----------------------
 * 4 Tests 1 Failures 0 Ignored
 * FAIL
 *
 * INTERPRETACE:
 * - Celkem 4 testy
 * - 3 prošly (PASS)
 * - 1 selhal (FAIL)
 * - Vidíme přesně který test a na kterém řádku (řádek 178)
 * - Exit code = 1 (failure)
 */

/* ============================================================================
 * PŘIDÁVÁNÍ NOVÝCH TESTŮ:
 * ============================================================================
 *
 * Když napíšeš nový test (např. v test_savgol.c):
 *
 * 1. Napiš testovací funkci v test_savgol.c:
 *    void test_savgol_constant(void) { ... }
 *
 * 2. Deklaruj ji zde nahoře:
 *    void test_savgol_constant(void);
 *
 * 3. Přidej RUN_TEST() do main():
 *    RUN_TEST(test_savgol_constant);
 *
 * 4. Rekompiluj a spusť:
 *    make test
 */
