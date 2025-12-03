/* test_timestamp.c - Unit tests for timestamp module */

#include "unity.h"
#include "../timestamp.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Test: parse_timestamp with space separator */
void test_parse_timestamp_space_separator(void) {
    double epoch;
    int result = parse_timestamp("2025-09-25 14:06:06.390", &epoch);

    TEST_ASSERT_EQUAL(0, result);
    TEST_ASSERT_TRUE(epoch > 0.0);

    /* Check subseconds are preserved (should have .390) */
    double fractional = epoch - floor(epoch);
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 0.390, fractional);
}

/* Test: parse_timestamp with T separator (RFC3339) */
void test_parse_timestamp_T_separator(void) {
    double epoch;
    int result = parse_timestamp("2025-09-25T14:06:06.390", &epoch);

    TEST_ASSERT_EQUAL(0, result);
    TEST_ASSERT_TRUE(epoch > 0.0);
}

/* Test: parse_timestamp without subseconds */
void test_parse_timestamp_no_subseconds(void) {
    double epoch;
    int result = parse_timestamp("2025-09-25 14:06:06", &epoch);

    TEST_ASSERT_EQUAL(0, result);
    TEST_ASSERT_TRUE(epoch > 0.0);

    /* Check no fractional part */
    double fractional = epoch - floor(epoch);
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 0.0, fractional);
}

/* Test: parse_timestamp with various subsecond precisions */
void test_parse_timestamp_subsecond_precision(void) {
    double epoch1, epoch2, epoch3;

    /* 1 digit (0.1 seconds) */
    parse_timestamp("2025-09-25 14:06:06.1", &epoch1);
    double frac1 = epoch1 - floor(epoch1);
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 0.1, frac1);

    /* 2 digits (0.12 seconds) */
    parse_timestamp("2025-09-25 14:06:06.12", &epoch2);
    double frac2 = epoch2 - floor(epoch2);
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 0.12, frac2);

    /* 3 digits (0.123 seconds) */
    parse_timestamp("2025-09-25 14:06:06.123", &epoch3);
    double frac3 = epoch3 - floor(epoch3);
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 0.123, frac3);
}

/* Test: parse_timestamp rejects invalid separator */
void test_parse_timestamp_invalid_separator(void) {
    double epoch;
    int result = parse_timestamp("2025-09-25X14:06:06.390", &epoch);

    TEST_ASSERT_EQUAL(-1, result);
}

/* Test: parse_timestamp rejects invalid date */
void test_parse_timestamp_invalid_date(void) {
    double epoch;

    /* Invalid month */
    TEST_ASSERT_EQUAL(-1, parse_timestamp("2025-13-25 14:06:06", &epoch));

    /* Invalid day */
    TEST_ASSERT_EQUAL(-1, parse_timestamp("2025-09-32 14:06:06", &epoch));

    /* Invalid hour */
    TEST_ASSERT_EQUAL(-1, parse_timestamp("2025-09-25 25:06:06", &epoch));
}

/* Test: parse_timestamp rejects NULL inputs */
void test_parse_timestamp_null_inputs(void) {
    double epoch;

    TEST_ASSERT_EQUAL(-1, parse_timestamp(NULL, &epoch));
    TEST_ASSERT_EQUAL(-1, parse_timestamp("2025-09-25 14:06:06", NULL));
}

/* Test: parse_timestamp rejects malformed strings */
void test_parse_timestamp_malformed(void) {
    double epoch;

    TEST_ASSERT_EQUAL(-1, parse_timestamp("not a timestamp", &epoch));
    TEST_ASSERT_EQUAL(-1, parse_timestamp("2025-09-25", &epoch));  /* Missing time */
    TEST_ASSERT_EQUAL(-1, parse_timestamp("14:06:06", &epoch));    /* Missing date */
}

/* Test: convert_timestamps_to_relative with valid data */
void test_convert_timestamps_basic(void) {
    char *timestamps[] = {
        "2025-09-25 14:06:06.000",
        "2025-09-25 14:06:07.000",
        "2025-09-25 14:06:08.500"
    };
    int n = 3;
    double *x_out = NULL;
    int first_error = -1;

    TimestampContext *ctx = convert_timestamps_to_relative(timestamps, n, &x_out, &first_error);

    TEST_ASSERT_NOT_NULL(ctx);
    TEST_ASSERT_NOT_NULL(x_out);
    TEST_ASSERT_EQUAL(3, ctx->n);
    TEST_ASSERT_EQUAL(0, ctx->errors_encountered);
    TEST_ASSERT_EQUAL(-1, first_error);

    /* First timestamp should be t=0 */
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 0.0, x_out[0]);

    /* Second timestamp: +1 second */
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 1.0, x_out[1]);

    /* Third timestamp: +2.5 seconds */
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 2.5, x_out[2]);

    /* Check original timestamps preserved */
    TEST_ASSERT_EQUAL_STRING("2025-09-25 14:06:06.000", ctx->original_timestamps[0]);
    TEST_ASSERT_EQUAL_STRING("2025-09-25 14:06:07.000", ctx->original_timestamps[1]);
    TEST_ASSERT_EQUAL_STRING("2025-09-25 14:06:08.500", ctx->original_timestamps[2]);

    free(x_out);
    free_timestamp_context(ctx);
}

/* Test: convert_timestamps_to_relative with mixed valid/invalid */
void test_convert_timestamps_with_errors(void) {
    char *timestamps[] = {
        "2025-09-25 14:06:06.000",
        "invalid timestamp",           /* Error on line 2 */
        "2025-09-25 14:06:08.000",
        "also invalid",                /* Error on line 4 */
        "2025-09-25 14:06:09.000"
    };
    int n = 5;
    double *x_out = NULL;
    int first_error = -1;

    TimestampContext *ctx = convert_timestamps_to_relative(timestamps, n, &x_out, &first_error);

    TEST_ASSERT_NOT_NULL(ctx);
    TEST_ASSERT_NOT_NULL(x_out);
    TEST_ASSERT_EQUAL(3, ctx->n);  /* Only 3 valid timestamps */
    TEST_ASSERT_EQUAL(2, ctx->errors_encountered);
    TEST_ASSERT_EQUAL(2, first_error);  /* First error on line 2 */

    /* Valid timestamps should have relative times */
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 0.0, x_out[0]);
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 2.0, x_out[1]);
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 3.0, x_out[2]);

    free(x_out);
    free_timestamp_context(ctx);
}

/* Test: convert_timestamps_to_relative with all invalid */
void test_convert_timestamps_all_invalid(void) {
    char *timestamps[] = {
        "invalid1",
        "invalid2",
        "invalid3"
    };
    int n = 3;
    double *x_out = NULL;
    int first_error = -1;

    TimestampContext *ctx = convert_timestamps_to_relative(timestamps, n, &x_out, &first_error);

    TEST_ASSERT_NULL(ctx);  /* Should return NULL when no valid timestamps */
    TEST_ASSERT_EQUAL(1, first_error);  /* First error on line 1 */
}

/* Test: convert_timestamps_to_relative preserves format */
void test_convert_timestamps_preserves_format(void) {
    char *timestamps[] = {
        "2025-09-25 14:06:06.1",      /* 1 decimal place */
        "2025-09-25T14:06:07.12",     /* T separator, 2 decimals */
        "2025-09-25 14:06:08.123"     /* 3 decimal places */
    };
    int n = 3;
    double *x_out = NULL;
    int first_error = -1;

    TimestampContext *ctx = convert_timestamps_to_relative(timestamps, n, &x_out, &first_error);

    TEST_ASSERT_NOT_NULL(ctx);

    /* Original strings should be preserved exactly */
    TEST_ASSERT_EQUAL_STRING("2025-09-25 14:06:06.1", ctx->original_timestamps[0]);
    TEST_ASSERT_EQUAL_STRING("2025-09-25T14:06:07.12", ctx->original_timestamps[1]);
    TEST_ASSERT_EQUAL_STRING("2025-09-25 14:06:08.123", ctx->original_timestamps[2]);

    free(x_out);
    free_timestamp_context(ctx);
}

/* Test: convert_timestamps_to_relative with NULL inputs */
void test_convert_timestamps_null_inputs(void) {
    char *timestamps[] = {"2025-09-25 14:06:06"};
    double *x_out = NULL;
    int first_error = -1;

    TEST_ASSERT_NULL(convert_timestamps_to_relative(NULL, 1, &x_out, &first_error));
    TEST_ASSERT_NULL(convert_timestamps_to_relative(timestamps, 0, &x_out, &first_error));
    TEST_ASSERT_NULL(convert_timestamps_to_relative(timestamps, 1, NULL, &first_error));
    TEST_ASSERT_NULL(convert_timestamps_to_relative(timestamps, 1, &x_out, NULL));
}

/* Test: free_timestamp_context with NULL */
void test_free_timestamp_context_null(void) {
    /* Should not crash */
    free_timestamp_context(NULL);
    TEST_ASSERT_TRUE(1);  /* If we get here, test passed */
}

/* Test: Time difference calculation accuracy */
void test_convert_timestamps_subsecond_accuracy(void) {
    /* Test data from example.dat (first 3 lines) */
    char *timestamps[] = {
        "2025-09-25 14:06:06.390",
        "2025-09-25 14:06:06.391",
        "2025-09-25 14:06:06.763"
    };
    int n = 3;
    double *x_out = NULL;
    int first_error = -1;

    TimestampContext *ctx = convert_timestamps_to_relative(timestamps, n, &x_out, &first_error);

    TEST_ASSERT_NOT_NULL(ctx);
    TEST_ASSERT_EQUAL(0, ctx->errors_encountered);

    /* First point: t=0 */
    TEST_ASSERT_DOUBLE_WITHIN(0.0001, 0.0, x_out[0]);

    /* Second point: +0.001 seconds (1 millisecond) */
    TEST_ASSERT_DOUBLE_WITHIN(0.0001, 0.001, x_out[1]);

    /* Third point: +0.373 seconds */
    TEST_ASSERT_DOUBLE_WITHIN(0.0001, 0.373, x_out[2]);

    free(x_out);
    free_timestamp_context(ctx);
}
