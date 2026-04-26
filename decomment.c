/*
 *  Decomment, remove all comments starting with '#' character until end of line
 *  Optimized version with improved error handling
 *  Original: V4/2016-12-24/TCh
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <stdbool.h>

#include "decomment.h"

/**
 * Processes an already-open stream by removing all text from '#' character
 * to the end of line and removing empty lines.
 *
 * Does not close the input stream — caller retains ownership. This lets the
 * same stripping logic run on stdin (which must not be closed) and on regular
 * files opened by decomment() below.
 */
FILE *decomment_stream(FILE *fr, const char *source_name)
{
  FILE *fw_tmp = NULL;
  int c;
  int p = 0;
  bool success = false;

  if (fr == NULL) {
    fprintf(stderr, "Error: NULL input stream provided\n");
    return NULL;
  }

  fw_tmp = tmpfile();
  if (fw_tmp == NULL) {
    fprintf(stderr, "Failed to open tmpfile (%d: %s)\n",
            errno, strerror(errno));
    return NULL;
  }

  while ((c = getc(fr)) != EOF) {
    if (c == '#') {
      while ((c = getc(fr)) != '\n' && c != EOF)
        ;
      if (c == EOF)
        break;
    }

    if (c == '\n' && p == 0)
      continue;

    if (c == '\n')
      p = 0;
    else
      p++;

    if (putc(c, fw_tmp) == EOF) {
      fprintf(stderr, "Error writing to temporary file (%d: %s)\n",
              errno, strerror(errno));
      goto cleanup;
    }
  }

  if (ferror(fr)) {
    fprintf(stderr, "Error reading from %s (%d: %s)\n",
            source_name ? source_name : "input stream",
            errno, strerror(errno));
    goto cleanup;
  }

  rewind(fw_tmp);
  success = true;

cleanup:
  if (!success && fw_tmp != NULL) {
    fclose(fw_tmp);
    fw_tmp = NULL;
  }
  return fw_tmp;
}

/**
 * Processes a file by removing all text from '#' character to the end of line
 * and removing empty lines. Thin wrapper around decomment_stream().
 */
FILE *decomment(const char *name)
{
  FILE *fr;
  FILE *out;

  if (name == NULL) {
    fprintf(stderr, "Error: NULL filename provided\n");
    return NULL;
  }

  fr = fopen(name, "r");
  if (fr == NULL) {
    fprintf(stderr, "Failed to open file '%s' (%d: %s)\n",
            name, errno, strerror(errno));
    return NULL;
  }

  out = decomment_stream(fr, name);

  if (fclose(fr) != 0) {
    fprintf(stderr, "Error closing input file (%d: %s)\n",
            errno, strerror(errno));
    if (out != NULL) {
      fclose(out);
      out = NULL;
    }
  }

  return out;
}
