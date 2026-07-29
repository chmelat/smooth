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

/* Timestamp-mode counterpart of run_smooth().
 *
 * A -T output row is "<timestamp tokens...> <y>", so the x field cannot be read
 * with sscanf("%lf %lf") — "2026-01-01" would parse as 2026 followed by -1.
 * Instead the last whitespace token is parsed as y, and a row is only accepted
 * when it starts with a digit (timestamps do; "Warning:"/"ERROR:" lines do not).
 *
 * The two "# Skipped ..." messages share a prefix, so they are told apart by
 * substring, not by sscanf format.
 */
typedef struct {
    int    data_rows;
    int    skip_nonnumeric;   /* -1 if the message was absent */
    int    skip_malformed;    /* -1 if the message was absent */
    double first_y, last_y;
} SmoothRunTs;

static SmoothRunTs run_smooth_ts(const char *args, const char *fixture) {
    char cmd[1024];
    snprintf(cmd, sizeof(cmd), SMOOTH_BIN " %s %s 2>&1", args, fixture);
    FILE *p = popen(cmd, "r");
    TEST_ASSERT_NOT_NULL_MESSAGE(p, "popen failed");

    SmoothRunTs r = {0};
    r.skip_nonnumeric = -1;
    r.skip_malformed  = -1;

    char ln[1024];
    while (fgets(ln, sizeof(ln), p)) {
        if (ln[0] == '#') {
            int k;
            if (sscanf(ln, "# Skipped %d", &k) == 1) {
                if (strstr(ln, "malformed timestamp") != NULL) {
                    r.skip_malformed = k;
                } else {
                    r.skip_nonnumeric = k;
                }
            }
            continue;
        }

        char *q = ln;
        while (*q == ' ' || *q == '\t') q++;
        if (*q < '0' || *q > '9') continue;  /* not a data row */

        /* Find the last whitespace-separated token and require it to be a
         * complete number. */
        size_t len = strlen(ln);
        while (len > 0 && (ln[len-1] == '\n' || ln[len-1] == '\r' ||
                           ln[len-1] == ' '  || ln[len-1] == '\t')) {
            ln[--len] = '\0';
        }
        if (len == 0) continue;
        char *tok = ln + len;
        while (tok > ln && tok[-1] != ' ' && tok[-1] != '\t') tok--;

        char *endptr;
        double yv = strtod(tok, &endptr);
        if (*endptr != '\0') continue;  /* trailing junk: not a data row */

        if (r.data_rows == 0) r.first_y = yv;
        r.last_y = yv;
        r.data_rows++;
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

/* '#' comments (full-line and inline) and blank lines are stripped by the
 * parser itself. A comment must not be counted as a skipped data row. */
void test_parser_comments_are_stripped(void) {
    const char *path = "/tmp/test_parser_comments.dat";
    write_fixture(path,
        "# header comment\n"
        "1 10\n"
        "2 11   # inline note\n"
        "\n"
        "3 12\n"
        "   # indented comment\n"
        "4 13\n"
        "5 14\n"
        "6 15\n");
    SmoothRun r = run_smooth("-m0 -n3 -p1", path);
    TEST_ASSERT_EQUAL_INT(6, r.data_rows);
    TEST_ASSERT_FALSE(r.has_skip_msg);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 1.0, r.first_x);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 6.0, r.last_x);
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 10.0, r.first_y);
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

/* -T mode, T-separated timestamp (ts_token_count = 1), default columns
 * (x = timestamp, y = last column). Data is exactly linear so -m0 -n3 -p1
 * reproduces y exactly. */
void test_parser_ts_t_format_default_columns(void) {
    const char *path = "/tmp/test_parser_ts_t.dat";
    write_fixture(path,
        "2026-01-01T00:00:00 10.0\n"
        "2026-01-01T00:00:01 11.0\n"
        "2026-01-01T00:00:02 12.0\n"
        "2026-01-01T00:00:03 13.0\n"
        "2026-01-01T00:00:04 14.0\n");
    SmoothRunTs r = run_smooth_ts("-T -m0 -n3 -p1", path);
    TEST_ASSERT_EQUAL_INT(5, r.data_rows);
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 10.0, r.first_y);
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 14.0, r.last_y);
    TEST_ASSERT_EQUAL_INT(-1, r.skip_nonnumeric);
    TEST_ASSERT_EQUAL_INT(-1, r.skip_malformed);
    remove(path);
}

/* -T mode, space-separated timestamp (ts_token_count = 2), default columns.
 * Same data as the previous test, just with the timestamp split across two
 * whitespace tokens. */
void test_parser_ts_space_format_default_columns(void) {
    const char *path = "/tmp/test_parser_ts_space.dat";
    write_fixture(path,
        "2026-01-01 00:00:00 10.0\n"
        "2026-01-01 00:00:01 11.0\n"
        "2026-01-01 00:00:02 12.0\n"
        "2026-01-01 00:00:03 13.0\n"
        "2026-01-01 00:00:04 14.0\n");
    SmoothRunTs r = run_smooth_ts("-T -m0 -n3 -p1", path);
    TEST_ASSERT_EQUAL_INT(5, r.data_rows);
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 10.0, r.first_y);
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 14.0, r.last_y);
    TEST_ASSERT_EQUAL_INT(-1, r.skip_nonnumeric);
    TEST_ASSERT_EQUAL_INT(-1, r.skip_malformed);
    remove(path);
}

/* -T mode, T-separated timestamp, y_column (3) > x_column (1): shift is 0
 * because the timestamp is already a single whitespace token. */
void test_parser_ts_k_maps_y_after_timestamp_t_format(void) {
    const char *path = "/tmp/test_parser_ts_k_t.dat";
    write_fixture(path,
        "2026-01-01T00:00:00 expA 10.0\n"
        "2026-01-01T00:00:01 expA 11.0\n"
        "2026-01-01T00:00:02 expA 12.0\n"
        "2026-01-01T00:00:03 expA 13.0\n"
        "2026-01-01T00:00:04 expA 14.0\n");
    SmoothRunTs r = run_smooth_ts("-T -k 1:3 -m0 -n3 -p1", path);
    TEST_ASSERT_EQUAL_INT(5, r.data_rows);
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 10.0, r.first_y);
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 14.0, r.last_y);
    TEST_ASSERT_EQUAL_INT(-1, r.skip_nonnumeric);
    TEST_ASSERT_EQUAL_INT(-1, r.skip_malformed);
    remove(path);
}

/* -T mode, space-separated timestamp, y_column (3) > x_column (1): shift is
 * +1 because the timestamp occupies two whitespace tokens. This pins the
 * logical-column abstraction: even though the space-separated row has one
 * more whitespace token than the T-separated row above, the SAME -k 1:3
 * flag is correct for both, because the timestamp counts as a single
 * logical column regardless of its token width. (-k 1:4 would instead fail
 * with "insufficient columns for y column 4".) */
void test_parser_ts_k_maps_y_after_timestamp_space_format(void) {
    const char *path = "/tmp/test_parser_ts_k_space.dat";
    write_fixture(path,
        "2026-01-01 00:00:00 expA 10.0\n"
        "2026-01-01 00:00:01 expA 11.0\n"
        "2026-01-01 00:00:02 expA 12.0\n"
        "2026-01-01 00:00:03 expA 13.0\n"
        "2026-01-01 00:00:04 expA 14.0\n");
    SmoothRunTs r = run_smooth_ts("-T -k 1:3 -m0 -n3 -p1", path);
    TEST_ASSERT_EQUAL_INT(5, r.data_rows);
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 10.0, r.first_y);
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 14.0, r.last_y);
    TEST_ASSERT_EQUAL_INT(-1, r.skip_nonnumeric);
    TEST_ASSERT_EQUAL_INT(-1, r.skip_malformed);
    remove(path);
}

/* -T mode, y_column (1) < x_column (2): logical columns before the timestamp
 * are unaffected by its token width, so no shift is applied. */
void test_parser_ts_y_before_timestamp(void) {
    const char *path = "/tmp/test_parser_ts_y_first.dat";
    write_fixture(path,
        "10.0 2026-01-01 00:00:00\n"
        "11.0 2026-01-01 00:00:01\n"
        "12.0 2026-01-01 00:00:02\n"
        "13.0 2026-01-01 00:00:03\n"
        "14.0 2026-01-01 00:00:04\n");
    SmoothRunTs r = run_smooth_ts("-T -k 2:1 -m0 -n3 -p1", path);
    TEST_ASSERT_EQUAL_INT(5, r.data_rows);
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 10.0, r.first_y);
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 14.0, r.last_y);
    TEST_ASSERT_EQUAL_INT(-1, r.skip_nonnumeric);
    TEST_ASSERT_EQUAL_INT(-1, r.skip_malformed);
    remove(path);
}

/* Regression guard for plan 001's silent-data-loss fix: a row whose timestamp
 * field fails to parse must be dropped *and reported* via "# Skipped N ...
 * malformed timestamp". Before 001 it was dropped silently, so the count is
 * what this guards -- the drop itself already worked. */
void test_parser_ts_malformed_row_is_reported(void) {
    const char *path = "/tmp/test_parser_ts_malformed.dat";
    write_fixture(path,
        "2026-01-01 00:00:00 1\n"
        "2026-01-01 00:00:01 2\n"
        "BADROW\n"
        "2026-01-01 00:00:03 4\n"
        "alsobad\n"
        "2026-01-01 00:00:05 6\n"
        "2026-01-01 00:00:06 7\n");
    SmoothRunTs r = run_smooth_ts("-T -m0 -n3 -p1", path);
    TEST_ASSERT_EQUAL_INT(5, r.data_rows);
    TEST_ASSERT_EQUAL_INT(2, r.skip_malformed);
    remove(path);
}

/* The skipped_nonnumeric path inside the -T branch (a syntactically valid
 * timestamp but a non-numeric y value), as distinct from a malformed
 * timestamp. */
void test_parser_ts_nonnumeric_y_is_reported(void) {
    const char *path = "/tmp/test_parser_ts_nan.dat";
    write_fixture(path,
        "2026-01-01 00:00:00 10.0\n"
        "2026-01-01 00:00:01 11.0\n"
        "2026-01-01 00:00:02 NaN\n"
        "2026-01-01 00:00:03 13.0\n"
        "2026-01-01 00:00:04 14.0\n"
        "2026-01-01 00:00:05 15.0\n");
    SmoothRunTs r = run_smooth_ts("-T -m0 -n3 -p1", path);
    TEST_ASSERT_EQUAL_INT(5, r.data_rows);
    TEST_ASSERT_EQUAL_INT(1, r.skip_nonnumeric);
    TEST_ASSERT_EQUAL_INT(-1, r.skip_malformed);
    remove(path);
}

/* CRLF input, -T with a space-separated timestamp.
 *
 * The timestamp tokenizer treated ' ', '\t' and '\n' as separators but not
 * '\r', so on a CRLF file the CR stayed attached to the last token and strtod
 * rejected it: every row was counted as non-numeric and the run died with
 * "ERROR: No valid data points found". The numeric branch never had the bug —
 * it lists '\r' explicitly — which is why this went unnoticed. Fixed in
 * v5.11.50 by stripping the terminator once after fgets.
 *
 * A regression here reports 0 data rows and skip_nonnumeric == 5. */
void test_parser_crlf_timestamp_space_format(void) {
    const char *path = "/tmp/test_parser_crlf_ts_space.dat";
    write_fixture(path,
        "2026-01-01 00:00:00 10.0\r\n"
        "2026-01-01 00:00:01 11.0\r\n"
        "2026-01-01 00:00:02 12.0\r\n"
        "2026-01-01 00:00:03 13.0\r\n"
        "2026-01-01 00:00:04 14.0\r\n");
    SmoothRunTs r = run_smooth_ts("-T -m0 -n3 -p1", path);
    TEST_ASSERT_EQUAL_INT(5, r.data_rows);
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 10.0, r.first_y);
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 14.0, r.last_y);
    TEST_ASSERT_EQUAL_INT(-1, r.skip_nonnumeric);
    TEST_ASSERT_EQUAL_INT(-1, r.skip_malformed);
    remove(path);
}

/* Same bug via the T-separated timestamp form (ts_token_count = 1), which
 * takes a different assembly path at parser.c:135-137. */
void test_parser_crlf_timestamp_t_format(void) {
    const char *path = "/tmp/test_parser_crlf_ts_t.dat";
    write_fixture(path,
        "2026-01-01T00:00:00 10.0\r\n"
        "2026-01-01T00:00:01 11.0\r\n"
        "2026-01-01T00:00:02 12.0\r\n"
        "2026-01-01T00:00:03 13.0\r\n"
        "2026-01-01T00:00:04 14.0\r\n");
    SmoothRunTs r = run_smooth_ts("-T -m0 -n3 -p1", path);
    TEST_ASSERT_EQUAL_INT(5, r.data_rows);
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 10.0, r.first_y);
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 14.0, r.last_y);
    TEST_ASSERT_EQUAL_INT(-1, r.skip_nonnumeric);
    TEST_ASSERT_EQUAL_INT(-1, r.skip_malformed);
    remove(path);
}

/* CRLF on the plain numeric path, with a comment line and a blank line mixed
 * in. This path already handled '\r'; the test pins that the shared strip did
 * not change it, and that comment stripping still works when '#' is followed
 * by CRLF rather than LF. */
void test_parser_crlf_numeric_with_comments(void) {
    const char *path = "/tmp/test_parser_crlf_num.dat";
    write_fixture(path,
        "# header comment\r\n"
        "\r\n"
        "1.0 10.0\r\n"
        "2.0 20.0  # inline comment\r\n"
        "3.0 30.0\r\n"
        "4.0 40.0\r\n"
        "5.0 50.0\r\n");
    SmoothRun r = run_smooth("-m0 -n3 -p1", path);
    TEST_ASSERT_EQUAL_INT(5, r.data_rows);
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 1.0,  r.first_x);
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 10.0, r.first_y);
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 5.0,  r.last_x);
    TEST_ASSERT_DOUBLE_WITHIN(0.01, 50.0, r.last_y);
    TEST_ASSERT_EQUAL_INT(0, r.has_skip_msg);
    remove(path);
}
