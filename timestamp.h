/* timestamp.h - RFC3339-style timestamp parsing
 * Supports formats: YYYY-MM-DD HH:MM:SS[.fff] or YYYY-MM-DDTHH:MM:SS[.fff]
 * No timezone support (assumes all timestamps in same timezone)
 */

#ifndef TIMESTAMP_H
#define TIMESTAMP_H

/* Original timestamp strings of the accepted rows, kept for output.
 * Built by parse_input(); free with free_timestamp_context(). */
typedef struct {
    char **original_timestamps;  /* Array of original timestamp strings */
    int n;                       /* Number of timestamps */
} TimestampContext;

/* Parse timestamp string to Unix epoch (seconds since 1970-01-01 00:00:00 UTC)
 * Accepts both space and 'T' as date/time separator
 * Format: YYYY-MM-DD HH:MM:SS[.fff] or YYYY-MM-DDTHH:MM:SS[.fff]
 *
 * Parameters:
 *   str: Input timestamp string
 *   epoch_seconds: Output pointer for epoch seconds (double for subsecond precision)
 *
 * Returns:
 *   0 on success
 *   -1 on parse error
 */
int parse_timestamp(const char *str, double *epoch_seconds);

/* Free timestamp context and all allocated memory */
void free_timestamp_context(TimestampContext *ctx);

#endif /* TIMESTAMP_H */
