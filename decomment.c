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

/**
 * Processes a file by removing all text from '#' character to the end of line
 * and removing empty lines.
 *
 * @param name Path to the input file
 * @return FILE pointer to a temporary file containing processed content or NULL on error
 */
FILE *decomment(const char *name)
{
  FILE *fr = NULL;      /* Input file descriptor */
  FILE *fw_tmp = NULL;  /* Output (temporary) file descriptor */
  int c;                /* Current character */
  int p = 0;            /* Position in current line */
  bool success = false; /* Flag to track successful completion */

  /* Validate input */
  if (name == NULL) {
    fprintf(stderr, "Error: NULL filename provided\n");
    return NULL;
  }

  /* Open input file */
  fr = fopen(name, "r");
  if (fr == NULL) {
    fprintf(stderr, "Failed to open file '%s' (%d: %s)\n",
            name, errno, strerror(errno));
    return NULL;
  }

  /* Create output (tmpfile) */
  fw_tmp = tmpfile();
  if (fw_tmp == NULL) {
    fprintf(stderr, "Failed to open tmpfile (%d: %s)\n",
            errno, strerror(errno));
    fclose(fr);
    return NULL;
  }

  /* Process file contents */
  while ((c = getc(fr)) != EOF) {
    /* Handle comment - skip all characters until newline or EOF */
    if (c == '#') {
      while ((c = getc(fr)) != '\n' && c != EOF)
        ; /* Skip until newline or EOF */
      
      /* If we hit EOF during comment, exit loop */
      if (c == EOF)
        break;
    }

    /* Skip empty lines (newline at position 0) */
    if (c == '\n' && p == 0)
      continue;

    /* Update position counter */
    if (c == '\n')
      p = 0;  /* Reset position at newline */
    else
      p++;    /* Increment position for other characters */

    /* Write character to output file */
    if (putc(c, fw_tmp) == EOF) {
      fprintf(stderr, "Error writing to temporary file (%d: %s)\n",
              errno, strerror(errno));
      goto cleanup;
    }
  }

  /* Check if we stopped due to an error rather than EOF */
  if (ferror(fr)) {
    fprintf(stderr, "Error reading from file '%s' (%d: %s)\n",
            name, errno, strerror(errno));
    goto cleanup;
  }

  /* Close input file */
  if (fclose(fr) != 0) {
    fprintf(stderr, "Error closing input file (%d: %s)\n",
            errno, strerror(errno));
    fr = NULL;  /* Avoid double-close in cleanup */
    goto cleanup;
  }
  fr = NULL;  /* Set to NULL to avoid double-close in cleanup */

  /* Reset file pointer to beginning of temporary file */
  rewind(fw_tmp);
  success = true;

cleanup:
  /* Clean up in case of error */
  if (!success) {
    if (fr != NULL)
      fclose(fr);
    if (fw_tmp != NULL) {
      fclose(fw_tmp);
      fw_tmp = NULL;
    }
  }

  return fw_tmp;
}
