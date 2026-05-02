/* parser.h - Input parser for the smooth program.
 *
 * Reads a whitespace-separated table of numeric values (normal mode) or
 * timestamp+value pairs (timestamp mode) from an open stream, and produces
 * parallel x[]/y[] arrays plus, in timestamp mode, a TimestampContext
 * holding the original timestamp strings and the relative-time conversion.
 */

#ifndef PARSER_H
#define PARSER_H

#include <stdio.h>
#include "timestamp.h"

typedef struct {
  double *x;                 /* allocated array of length n */
  double *y;                 /* allocated array of length n */
  int n;                     /* number of valid data points */
  TimestampContext *ts_ctx;  /* NULL unless timestamp_mode */
} ParseResult;

/* Parse a smooth-format data table from `fp`. The caller owns `fp` and
 * is responsible for closing it. On success the function returns 0 and
 * populates `result`; ownership of `result->x`, `result->y`, and
 * `result->ts_ctx` transfers to the caller. On failure it returns a
 * non-zero status, prints `ERROR:` to stderr, and frees any partial
 * allocations (the fields of `result` are NULL/0 on return).
 *
 * `timestamp_mode`: 0 for two numeric columns, 1 for timestamp + value.
 * `x_column`, `y_column`: 1-indexed logical columns (in timestamp mode
 *   `x_column` is the position of the timestamp itself).
 *
 * Diagnostic output written by the parser:
 *   stdout `# Skipped %d ...`  rows whose selected x/y column was
 *                              non-numeric or NaN/Inf
 *   stderr `Warning: Skipped %d line(s) with invalid timestamps ...`
 *                              timestamp_mode only
 */
int parse_input(FILE *fp,
                int timestamp_mode,
                int x_column,
                int y_column,
                ParseResult *result);

/* Free the contents of `result` (does not free the struct itself). */
void free_parse_result(ParseResult *result);

#endif
