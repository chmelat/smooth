/* parser.c - Input parser for the smooth program. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "parser.h"
#include "timestamp.h"

#define BUF      512
#define MAX_LINE 4096
#define MAX_COLS 100

int parse_input(FILE *fp,
                int timestamp_mode,
                int x_column,
                int y_column,
                ParseResult *result)
{
  char line[MAX_LINE];
  int line_number = 0;
  int skipped_nonnumeric = 0;
  int n = 0;
  int abuf = 0;
  double *x = NULL;
  double *y = NULL;
  char **timestamp_strings = NULL;
  TimestampContext *ts_ctx = NULL;

  result->x = NULL;
  result->y = NULL;
  result->n = 0;
  result->ts_ctx = NULL;

  if (timestamp_mode) {
    timestamp_strings = malloc(BUF * sizeof(char*));
    if (!timestamp_strings) {
      fprintf(stderr, "ERROR: No memory for timestamp strings\n");
      goto fail;
    }
    y = malloc(BUF * sizeof(double));
    if (!y) {
      fprintf(stderr, "ERROR: No memory for data table\n");
      goto fail;
    }
    abuf = BUF;
  }

  while (fgets(line, sizeof(line), fp) != NULL) {
    line_number++;

    /* Detect line overflow: buffer filled without trailing newline AND
     * stream still has data — line was truncated mid-content (audit B9). */
    {
      size_t llen = strlen(line);
      if (llen == sizeof(line) - 1 && line[llen-1] != '\n') {
        int c = fgetc(fp);
        if (c != EOF) {
          ungetc(c, fp);
          fprintf(stderr,
                  "ERROR: Line %d exceeds %zu-byte read buffer (MAX_LINE). "
                  "Increase MAX_LINE in parser.c or shorten the input line.\n",
                  line_number, sizeof(line));
          goto fail;
        }
      }
    }

    if (timestamp_mode) {
      /* Timestamp mode with logical-column model: timestamp lives at logical
       * column x_column (default 1), y at logical column y_column (default 2).
       * Timestamp itself spans 1 (T-separator) or 2 (space-separator) whitespace
       * tokens; the logical-column abstraction hides that from the user. */
      char timestamp_str[100];
      double y_value;

      /* Tokenize line on whitespace (destructive on local MAX_LINE buffer) */
      char *tokens[MAX_COLS];
      int ntok = 0;
      char *p = line;
      while (*p && ntok < MAX_COLS) {
        while (*p == ' ' || *p == '\t' || *p == '\n') *p++ = '\0';
        if (*p == '\0') break;
        tokens[ntok++] = p;
        while (*p && *p != ' ' && *p != '\t' && *p != '\n') p++;
      }
      if (ntok == 0) continue;  /* blank or whitespace-only line */

      /* Detect token overflow: hit cap with more non-whitespace data on line (audit B9) */
      if (ntok == MAX_COLS) {
        char *check = p;
        while (*check == ' ' || *check == '\t' || *check == '\n') check++;
        if (*check != '\0') {
          fprintf(stderr,
                  "ERROR: Line %d has more than %d tokens (MAX_COLS). "
                  "Increase MAX_COLS in parser.c.\n",
                  line_number, MAX_COLS);
          goto fail;
        }
      }

      /* Validate timestamp position (x_column is 1-indexed logical column) */
      int ts_tok_start = x_column - 1;
      if (ts_tok_start >= ntok) {
        fprintf(stderr, "ERROR: Line %d has %d token(s), but timestamp column %d was requested\n",
                line_number, ntok, x_column);
        goto fail;
      }

      /* Detect timestamp format and assemble timestamp_str */
      int ts_token_count;
      if (strchr(tokens[ts_tok_start], 'T') != NULL) {
        ts_token_count = 1;
        snprintf(timestamp_str, sizeof(timestamp_str), "%s", tokens[ts_tok_start]);
      } else if (ts_tok_start + 1 < ntok) {
        ts_token_count = 2;
        snprintf(timestamp_str, sizeof(timestamp_str), "%s %s",
                 tokens[ts_tok_start], tokens[ts_tok_start + 1]);
      } else {
        continue;  /* malformed: space-format expects two tokens, only one present */
      }

      /* Map logical y_column to whitespace-token index. Logical columns before
       * the timestamp are unaffected by its width; columns after shift by
       * ts_token_count - 1. y_column != x_column is enforced at -k parse time. */
      int y_token_idx = (y_column < x_column)
                        ? y_column - 1
                        : y_column - 1 + (ts_token_count - 1);
      if (y_token_idx >= ntok) {
        fprintf(stderr, "ERROR: Line %d has insufficient columns for y column %d\n",
                line_number, y_column);
        goto fail;
      }

      char *endptr;
      errno = 0;
      y_value = strtod(tokens[y_token_idx], &endptr);
      char *y_tok_end = tokens[y_token_idx] + strlen(tokens[y_token_idx]);
      if (endptr != y_tok_end || errno != 0 || isnan(y_value) || isinf(y_value)) {
        skipped_nonnumeric++;
        continue;  /* token is not a fully numeric finite value */
      }

      /* Reallocate arrays if needed */
      if (n == abuf) {
        abuf = abuf ? abuf * 2 : BUF;

        char **temp_ts = (char**)realloc(timestamp_strings, abuf * sizeof(char*));
        if (!temp_ts) {
          fprintf(stderr, "ERROR: No memory for timestamp strings\n");
          goto fail;
        }
        timestamp_strings = temp_ts;

        double *temp_y = (double*)realloc(y, abuf * sizeof(double));
        if (!temp_y) {
          fprintf(stderr, "ERROR: No memory for data table\n");
          goto fail;
        }
        y = temp_y;
      }

      timestamp_strings[n] = strdup(timestamp_str);
      if (!timestamp_strings[n]) {
        fprintf(stderr, "ERROR: No memory for timestamp string\n");
        goto fail;
      }
      y[n] = y_value;
      n++;

    } else {
      /* Normal mode: parse whitespace-separated tokens.
       * Each token = one logical column. A token that strtod cannot fully
       * consume (e.g. ISO timestamp 2026-04-29T11:40:00, label "abc",
       * partially numeric "1.5e2x") is a placeholder: the column position
       * is preserved, but no numeric value is available. Rows where the
       * selected x_column or y_column lands on a placeholder (or NaN/Inf)
       * are skipped with a per-file summary on stdout. */
      double values[MAX_COLS];
      int placeholder[MAX_COLS] = {0};
      int col_count = 0;
      char *ptr = line;

      while (*ptr == ' ' || *ptr == '\t') ptr++;
      if (*ptr == '\n' || *ptr == '\r' || *ptr == '\0') continue;

      while (*ptr && *ptr != '\n' && *ptr != '\r' && col_count < MAX_COLS) {
        char *tok = ptr;
        while (*ptr && *ptr != ' ' && *ptr != '\t' && *ptr != '\n' && *ptr != '\r')
          ptr++;
        char *tok_end = ptr;

        char *endptr;
        errno = 0;
        double v = strtod(tok, &endptr);
        if (endptr == tok_end && errno == 0 && !isnan(v) && !isinf(v)) {
          values[col_count] = v;
        } else {
          values[col_count] = 0.0;
          placeholder[col_count] = 1;
        }
        col_count++;

        while (*ptr == ' ' || *ptr == '\t') ptr++;
      }

      /* Detect column overflow: hit cap with more tokens still on line (audit B9) */
      if (col_count == MAX_COLS && *ptr != '\n' && *ptr != '\r' && *ptr != '\0') {
        fprintf(stderr,
                "ERROR: Line %d has more than %d columns (MAX_COLS). "
                "Increase MAX_COLS in parser.c.\n",
                line_number, MAX_COLS);
        goto fail;
      }

      if (col_count < 1) continue;

      int max_col = (x_column > y_column) ? x_column : y_column;
      if (col_count < max_col) {
        fprintf(stderr, "ERROR: Line %d has only %d column(s), but columns %d (x) and %d (y) were requested\n",
                line_number, col_count, x_column, y_column);
        goto fail;
      }

      if (placeholder[x_column - 1] || placeholder[y_column - 1]) {
        skipped_nonnumeric++;
        continue;
      }

      if (n == abuf) {
        abuf = abuf ? abuf * 2 : BUF;

        double *temp_x = (double *)realloc(x, abuf * sizeof(double));
        if (temp_x == NULL) {
          fprintf(stderr, "ERROR: No memory for data table\n");
          goto fail;
        }
        x = temp_x;

        double *temp_y = (double *)realloc(y, abuf * sizeof(double));
        if (temp_y == NULL) {
          fprintf(stderr, "ERROR: No memory for data table\n");
          goto fail;
        }
        y = temp_y;
      }

      x[n] = values[x_column - 1];
      y[n] = values[y_column - 1];
      n++;
    }
  }

  if (skipped_nonnumeric > 0) {
    if (timestamp_mode) {
      printf("# Skipped %d data row(s) with non-numeric or NaN/Inf value in column %d (y)\n",
             skipped_nonnumeric, y_column);
    } else {
      printf("# Skipped %d data row(s) with non-numeric or NaN/Inf value in column %d (x) or %d (y)\n",
             skipped_nonnumeric, x_column, y_column);
    }
  }

  if (timestamp_mode) {
    if (n == 0) {
      fprintf(stderr, "ERROR: No valid data points found\n");
      goto fail;
    }

    int first_error_line = -1;
    ts_ctx = convert_timestamps_to_relative(timestamp_strings, n, &x, &first_error_line);
    if (ts_ctx == NULL) {
      fprintf(stderr, "ERROR: No valid timestamps found in input\n");
      if (first_error_line > 0) {
        fprintf(stderr, "First invalid timestamp at line %d\n", first_error_line);
      }
      goto fail;
    }

    if (ts_ctx->errors_encountered > 0) {
      fprintf(stderr, "Warning: Skipped %d line(s) with invalid timestamps (first error at line %d)\n",
              ts_ctx->errors_encountered, first_error_line);
    }

    int n_parsed = n;
    n = ts_ctx->n;

    for (int i = 0; i < n_parsed; i++) free(timestamp_strings[i]);
    free(timestamp_strings);
    timestamp_strings = NULL;
  }

  result->x = x;
  result->y = y;
  result->n = n;
  result->ts_ctx = ts_ctx;
  return 0;

fail:
  if (timestamp_strings) {
    for (int i = 0; i < n; i++) free(timestamp_strings[i]);
    free(timestamp_strings);
  }
  free(x);
  free(y);
  if (ts_ctx) free_timestamp_context(ts_ctx);
  return 1;
}

void free_parse_result(ParseResult *result)
{
  if (!result) return;
  free(result->x);
  free(result->y);
  if (result->ts_ctx) free_timestamp_context(result->ts_ctx);
  result->x = NULL;
  result->y = NULL;
  result->n = 0;
  result->ts_ctx = NULL;
}
