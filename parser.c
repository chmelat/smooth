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
  int skipped_malformed_ts = 0;
  int first_malformed_ts_line = 0;
  int n = 0;
  int abuf = 0;
  double *x = NULL;
  double *y = NULL;
  int *line_no = NULL;   /* input-file line of each accepted row, for messages */
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
    line_no = malloc(BUF * sizeof(int));
    if (!y || !line_no) {
      fprintf(stderr, "ERROR: No memory for data table\n");
      goto fail;
    }
    abuf = BUF;
  }

  while (fgets(line, sizeof(line), fp) != NULL) {
    line_number++;

    /* Detect line overflow: buffer filled without trailing newline AND the
     * next character is not the line's terminator — line was truncated
     * mid-content (audit B9). A line of exactly sizeof(line)-1 bytes leaves
     * only its LF (or CR LF) unread; that is consumed here, not an overflow.
     * Must run before the comment strip below, which changes strlen(line). */
    int truncated = 0;
    {
      size_t llen = strlen(line);
      if (llen == sizeof(line) - 1 && line[llen-1] != '\n') {
        int c = fgetc(fp);
        if (c == '\r')
          c = fgetc(fp);
        if (c != '\n' && c != EOF) {
          ungetc(c, fp);
          truncated = 1;
        }
      }
    }

    /* Strip the line terminator, CR included, so CRLF input parses. The two
     * tokenizers below disagree about '\r': the numeric branch treats it as a
     * terminator, the timestamp branch does not, so on a CRLF file the CR stuck
     * to the last token and strtod rejected it — every row of the file was
     * dropped as non-numeric. Stripping once here fixes both branches and any
     * future one. Must run after the truncation check above, which needs the
     * raw '\n', and before the comment strip below. */
    {
      size_t llen = strlen(line);
      while (llen > 0 && (line[llen-1] == '\n' || line[llen-1] == '\r'))
        line[--llen] = '\0';
    }

    /* Strip '#' comment (full-line or inline) through end of line. A line that
     * is nothing but a comment becomes empty and is dropped by the blank-line
     * checks below, so it never counts as a skipped data row. */
    char *hash = strchr(line, '#');
    if (hash) *hash = '\0';

    if (truncated) {
      if (hash) {
        /* Truncation fell inside the comment — everything before '#' was read
         * intact. Discard the rest of the physical line and parse what we have. */
        int c;
        while ((c = fgetc(fp)) != '\n' && c != EOF)
          ;
      } else {
        fprintf(stderr,
                "ERROR: Line %d exceeds %zu-byte read buffer (MAX_LINE). "
                "Increase MAX_LINE in parser.c or shorten the input line.\n",
                line_number, sizeof(line));
        goto fail;
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

      /* Assemble the timestamp (x_column is its 1-indexed logical column):
       * one token if it contains 'T', else two (date, time). A row whose
       * timestamp is missing or does not parse -- a header such as
       * "date value" (audit A6), a damaged or short row -- is skipped; the
       * summary names the first such line.
       * ponytail: convert_timestamps_to_relative() parses it again; keep the
       * epoch from here if that ever shows up in a profile. */
      int ts_tok_start = x_column - 1;
      int ts_token_count = 0;
      if (ts_tok_start < ntok && strchr(tokens[ts_tok_start], 'T') != NULL) {
        ts_token_count = 1;
        snprintf(timestamp_str, sizeof(timestamp_str), "%s", tokens[ts_tok_start]);
      } else if (ts_tok_start + 1 < ntok) {
        ts_token_count = 2;
        snprintf(timestamp_str, sizeof(timestamp_str), "%s %s",
                 tokens[ts_tok_start], tokens[ts_tok_start + 1]);
      }
      double epoch;
      if (ts_token_count == 0 || parse_timestamp(timestamp_str, &epoch) != 0) {
        if (skipped_malformed_ts++ == 0) first_malformed_ts_line = line_number;
        continue;
      }

      /* Map logical y_column to whitespace-token index. Logical columns before
       * the timestamp are unaffected by its width; columns after shift by
       * ts_token_count - 1. y_column != x_column is enforced at -k parse time. */
      int y_token_idx = (y_column < x_column)
                        ? y_column - 1
                        : y_column - 1 + (ts_token_count - 1);
      if (y_token_idx >= ntok) {  /* valid timestamp, no y: broken data row */
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

        int *temp_line = (int*)realloc(line_no, abuf * sizeof(int));
        if (!temp_line) {
          fprintf(stderr, "ERROR: No memory for data table\n");
          goto fail;
        }
        line_no = temp_line;
      }

      timestamp_strings[n] = strdup(timestamp_str);
      if (!timestamp_strings[n]) {
        fprintf(stderr, "ERROR: No memory for timestamp string\n");
        goto fail;
      }
      y[n] = y_value;
      line_no[n] = line_number;
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

        int *temp_line = (int *)realloc(line_no, abuf * sizeof(int));
        if (temp_line == NULL) {
          fprintf(stderr, "ERROR: No memory for data table\n");
          goto fail;
        }
        line_no = temp_line;
      }

      x[n] = values[x_column - 1];
      y[n] = values[y_column - 1];
      line_no[n] = line_number;
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

  if (skipped_malformed_ts > 0) {
    printf("# Skipped %d data row(s) with malformed timestamp in column %d "
           "(first at line %d)\n",
           skipped_malformed_ts, x_column, first_malformed_ts_line);
  }

  if (timestamp_mode) {
    if (n == 0) {
      fprintf(stderr, "ERROR: No valid data points found\n");
      goto fail;
    }

    /* Every timestamp was validated above, so this fails only on memory. */
    int first_error_line;
    ts_ctx = convert_timestamps_to_relative(timestamp_strings, n, y, line_no, &x, &first_error_line);
    if (ts_ctx == NULL) {
      fprintf(stderr, "ERROR: Timestamp conversion failed\n");
      goto fail;
    }

    for (int i = 0; i < n; i++) free(timestamp_strings[i]);
    free(timestamp_strings);
    timestamp_strings = NULL;
  }

  /* Checked here, not only in analyze_grid(), because only the parser knows
   * which file line a row came from (audit A5). */
  for (int i = 1; i < n; i++) {
    if (x[i] <= x[i-1]) {
      fprintf(stderr, "ERROR: %s not strictly increasing at line %d "
              "(previous data row: line %d)\n",
              timestamp_mode ? "Timestamps" : "x data", line_no[i], line_no[i-1]);
      goto fail;
    }
  }
  free(line_no);

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
  free(line_no);
  if (ts_ctx) free_timestamp_context(ts_ctx);
  return 1;
}
