/* test_timestamp.c - Unit tests for timestamp module */

#include "unity.h"
#include "../timestamp.h"
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

/* Test: parse_timestamp rejects calendar-impossible dates (audit A2) */
void test_parse_timestamp_nonexistent_date(void) {
    double epoch;

    /* Days that pass the 1-31 range check but do not exist; timegm() would
     * silently normalize these forward without the post-conversion check. */
    TEST_ASSERT_EQUAL(-1, parse_timestamp("2025-02-31 00:00:00", &epoch));
    TEST_ASSERT_EQUAL(-1, parse_timestamp("2025-04-31 00:00:00", &epoch));
    TEST_ASSERT_EQUAL(-1, parse_timestamp("2025-02-29 00:00:00", &epoch));  /* 2025 not leap */

    /* Valid leap day must still be accepted (regression guard) */
    TEST_ASSERT_EQUAL(0, parse_timestamp("2024-02-29 00:00:00", &epoch));   /* 2024 is leap */
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

/* Test: free_timestamp_context with NULL */
void test_free_timestamp_context_null(void) {
    /* Should not crash */
    free_timestamp_context(NULL);
    TEST_ASSERT_TRUE(1);  /* If we get here, test passed */
}

/* Test: DST transition does not corrupt relative timestamps.
 * 2025-03-30 is CET→CEST spring-forward in Europe/Prague.
 * 01:59 CET + 62 min = 03:01 CEST.  With timegm (UTC) the
 * difference must be exactly 3720s regardless of local TZ. */
void test_parse_timestamp_dst_invariant(void) {
    double epoch_before, epoch_after;

    int r1 = parse_timestamp("2025-03-30 01:59:00", &epoch_before);
    int r2 = parse_timestamp("2025-03-30 03:01:00", &epoch_after);

    TEST_ASSERT_EQUAL(0, r1);
    TEST_ASSERT_EQUAL(0, r2);

    double diff = epoch_after - epoch_before;
    TEST_ASSERT_DOUBLE_WITHIN(0.001, 3720.0, diff);  /* exactly 62 minutes */
}

/* Test: sub-millisecond differences survive the epoch arithmetic */
void test_parse_timestamp_subsecond_differences(void) {
    /* Test data from example.dat (first 3 lines) */
    double t0, t1, t2;
    TEST_ASSERT_EQUAL(0, parse_timestamp("2025-09-25 14:06:06.390", &t0));
    TEST_ASSERT_EQUAL(0, parse_timestamp("2025-09-25 14:06:06.391", &t1));
    TEST_ASSERT_EQUAL(0, parse_timestamp("2025-09-25 14:06:06.763", &t2));
    TEST_ASSERT_DOUBLE_WITHIN(0.0001, 0.001, t1 - t0);
    TEST_ASSERT_DOUBLE_WITHIN(0.0001, 0.373, t2 - t0);
}
