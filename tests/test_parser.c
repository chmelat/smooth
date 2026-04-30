/* test_parser.c - End-to-end tests for the smooth input parser.
 *
 * Drives the ./smooth binary via popen() on small generated fixtures
 * under /tmp, parses data rows from its stdout, and verifies column
 * counts, values, and the "# Skipped N ..." summary message.
 *
 * Requires that ./smooth is built (Makefile target `test` adds the
 * dependency).
 */

#include "unity.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SMOOTH_BIN "./smooth"

typedef struct {
    int data_rows;
    int has_skip_msg;
    int skip_count;
    double first_x, first_y;
    double last_x,  last_y;
} SmoothRun;

static void write_fixture(const char *path, const char *content) {
    FILE *f = fopen(path, "w");
    TEST_ASSERT_NOT_NULL_MESSAGE(f, "could not open fixture file for writing");
    fputs(content, f);
    fclose(f);
}

static SmoothRun run_smooth(const char *args, const char *fixture) {
    char cmd[1024];
    snprintf(cmd, sizeof(cmd), SMOOTH_BIN " %s %s 2>&1", args, fixture);
    FILE *p = popen(cmd, "r");
    TEST_ASSERT_NOT_NULL_MESSAGE(p, "popen failed");

    SmoothRun r = {0};
    r.skip_count = -1;

    char ln[1024];
    while (fgets(ln, sizeof(ln), p)) {
        if (ln[0] == '#') {
            int n;
            if (sscanf(ln, "# Skipped %d data row", &n) == 1) {
                r.has_skip_msg = 1;
                r.skip_count = n;
            }
            continue;
        }
        char *q = ln;
        while (*q == ' ' || *q == '\t') q++;
        if (*q == '\n' || *q == '\r' || *q == '\0') continue;

        double xv, yv;
        if (sscanf(ln, "%lf %lf", &xv, &yv) == 2) {
            if (r.data_rows == 0) { r.first_x = xv; r.first_y = yv; }
            r.last_x = xv;
            r.last_y = yv;
            r.data_rows++;
        }
    }
    pclose(p);
    return r;
}

/* T-separated ISO timestamp counts as one whitespace token, so a row
 * "2026-04-29T11:40:00 0 1 2 3" has 5 logical columns. -k 2:5 selects
 * the seconds counter (col 2) as x and the last value (col 3) as y. */
void test_parser_iso_t_timestamp_one_column(void) {
    const char *path = "/tmp/test_parser_iso_t.dat";
    write_fixture(path,
        "2026-04-29T11:40:00 0 1.0 2.0 3.0\n"
        "2026-04-29T11:40:01 1 1.1 2.1 3.1\n"
        "2026-04-29T11:40:02 2 1.2 2.2 3.2\n"
        "2026-04-29T11:40:03 3 1.3 2.3 3.3\n"
        "2026-04-29T11:40:04 4 1.4 2.4 3.4\n"
        "2026-04-29T11:40:05 5 1.5 2.5 3.5\n");
    SmoothRun r = run_smooth("-m0 -n3 -p1 -k 2:5", path);
    TEST_ASSERT_EQUAL_INT(6, r.data_rows);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 0.0, r.first_x);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 5.0, r.last_x);
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 3.5, r.last_y);
    TEST_ASSERT_FALSE(r.has_skip_msg);
    remove(path);
}

/* Space-separated timestamp is two whitespace tokens (date + time), so the
 * same data has 6 logical columns; -k 3:6 selects the seconds counter and
 * the last value. */
void test_parser_iso_space_timestamp_two_columns(void) {
    const char *path = "/tmp/test_parser_iso_space.dat";
    write_fixture(path,
        "2026-04-29 11:40:00 0 1.0 2.0 3.0\n"
        "2026-04-29 11:40:01 1 1.1 2.1 3.1\n"
        "2026-04-29 11:40:02 2 1.2 2.2 3.2\n"
        "2026-04-29 11:40:03 3 1.3 2.3 3.3\n"
        "2026-04-29 11:40:04 4 1.4 2.4 3.4\n"
        "2026-04-29 11:40:05 5 1.5 2.5 3.5\n");
    SmoothRun r = run_smooth("-m0 -n3 -p1 -k 3:6", path);
    TEST_ASSERT_EQUAL_INT(6, r.data_rows);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 0.0, r.first_x);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 5.0, r.last_x);
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 3.5, r.last_y);
    TEST_ASSERT_FALSE(r.has_skip_msg);
    remove(path);
}

/* "3.5e2x" — strtod consumes 3.5e2 but leaves trailing "x" — placeholder. */
void test_parser_partial_numeric_token_is_placeholder(void) {
    const char *path = "/tmp/test_parser_partial.dat";
    write_fixture(path,
        "1      10\n"
        "2      11\n"
        "3.5e2x 12\n"
        "4      13\n"
        "5      14\n"
        "6      15\n");
    SmoothRun r = run_smooth("-m0 -n3 -p1", path);
    TEST_ASSERT_EQUAL_INT(5, r.data_rows);
    TEST_ASSERT_TRUE(r.has_skip_msg);
    TEST_ASSERT_EQUAL_INT(1, r.skip_count);
    remove(path);
}

/* NaN literal in y must not propagate into the output. */
void test_parser_nan_in_y_skips_row(void) {
    const char *path = "/tmp/test_parser_nan_y.dat";
    write_fixture(path,
        "1 10\n"
        "2 11\n"
        "3 NaN\n"
        "4 13\n"
        "5 14\n"
        "6 15\n");
    SmoothRun r = run_smooth("-m0 -n3 -p1", path);
    TEST_ASSERT_EQUAL_INT(5, r.data_rows);
    TEST_ASSERT_TRUE(r.has_skip_msg);
    TEST_ASSERT_EQUAL_INT(1, r.skip_count);
    TEST_ASSERT_FALSE(isnan(r.last_y));
    TEST_ASSERT_FALSE(isinf(r.last_y));
    remove(path);
}

/* Inf literal in x column also triggers skip. */
void test_parser_inf_in_x_skips_row(void) {
    const char *path = "/tmp/test_parser_inf_x.dat";
    write_fixture(path,
        "1   10\n"
        "2   11\n"
        "Inf 12\n"
        "4   13\n"
        "5   14\n"
        "6   15\n");
    SmoothRun r = run_smooth("-m0 -n3 -p1", path);
    TEST_ASSERT_EQUAL_INT(5, r.data_rows);
    TEST_ASSERT_TRUE(r.has_skip_msg);
    TEST_ASSERT_EQUAL_INT(1, r.skip_count);
    remove(path);
}

/* A non-numeric label in a column outside the selected x/y is harmless. */
void test_parser_label_outside_xy_is_harmless(void) {
    const char *path = "/tmp/test_parser_label.dat";
    write_fixture(path,
        "expA 1 10\n"
        "expA 2 11\n"
        "expA 3 12\n"
        "expA 4 13\n"
        "expA 5 14\n"
        "expA 6 15\n");
    SmoothRun r = run_smooth("-m0 -n3 -p1 -k 2:3", path);
    TEST_ASSERT_EQUAL_INT(6, r.data_rows);
    TEST_ASSERT_FALSE(r.has_skip_msg);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 1.0, r.first_x);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 6.0, r.last_x);
    remove(path);
}
