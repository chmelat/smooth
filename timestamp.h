/* timestamp.h - RFC3339-style timestamp parsing and conversion
 * Supports formats: YYYY-MM-DD HH:MM:SS[.fff] or YYYY-MM-DDTHH:MM:SS[.fff]
 * No timezone support (assumes all timestamps in same timezone)
 */

#ifndef TIMESTAMP_H
#define TIMESTAMP_H

/* Context for timestamp conversion - preserves original strings */
typedef struct {
    char **original_timestamps;  /* Array of original timestamp strings */
    double reference_epoch;      /* Unix epoch of first timestamp (seconds) */
    int n;                       /* Number of timestamps */
    int errors_encountered;      /* Count of skipped invalid timestamps */
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

/* Convert array of timestamp strings to relative time in seconds
 * First valid timestamp becomes reference (t=0)
 * Subsequent timestamps are converted to seconds since reference
 * Invalid timestamps are skipped with warning on first occurrence
 *
 * Parameters:
 *   timestamp_strings: Array of timestamp strings from input
 *   n: Number of timestamps
 *   x_out: Output pointer for relative time array (allocated by function)
 *   first_error_line: Output pointer for line number of first error (or -1 if none)
 *
 * Returns:
 *   TimestampContext pointer on success (must be freed with free_timestamp_context)
 *   NULL on error (no valid timestamps found)
 */
TimestampContext* convert_timestamps_to_relative(
    char **timestamp_strings,
    int n,
    double **x_out,
    int *first_error_line
);

/* Free timestamp context and all allocated memory */
void free_timestamp_context(TimestampContext *ctx);

#endif /* TIMESTAMP_H */
