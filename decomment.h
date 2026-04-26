/**
 * @file decomment.h
 * @brief Header file for decomment functionality
 * 
 * Declares the decomment function which removes comments and empty lines from a file
 */

#ifndef DECOMMENT_H
#define DECOMMENT_H

#include <stdio.h>

/**
 * Processes a file by removing all text from '#' character to the end of line
 * and removing empty lines.
 *
 * @param name Path to the input file
 * @return FILE pointer to a temporary file containing processed content or NULL on error
 */
extern FILE *decomment(const char *name);

/**
 * Processes an already-open stream by removing all text from '#' character
 * to the end of line and removing empty lines. Does not close the input stream
 * (caller retains ownership).
 *
 * @param fr           Open input stream
 * @param source_name  Name used in error messages (e.g. "stdin", filename)
 * @return FILE pointer to a temporary file containing processed content or NULL on error
 */
extern FILE *decomment_stream(FILE *fr, const char *source_name);

#endif /* DECOMMENT_H */
