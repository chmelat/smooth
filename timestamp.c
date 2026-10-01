/* timestamp.c - RFC3339-style timestamp parsing */

#define _DEFAULT_SOURCE  /* for timegm() */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "timestamp.h"

/* Parse timestamp string to Unix epoch seconds
 * Accepts: YYYY-MM-DD HH:MM:SS[.fff] or YYYY-MM-DDTHH:MM:SS[.fff]
 */
int parse_timestamp(const char *str, double *epoch_seconds)
{
    if (!str || !epoch_seconds) {
        return -1;
    }

    struct tm tm_time = {0};
    double subseconds = 0.0;
    int year, month, day, hour, min, sec;
    char separator;
    int pos = 0;

    /* Parse date and time without subseconds */
    int matched = sscanf(str, "%d-%d-%d%c%d:%d:%d%n",
                        &year, &month, &day, &separator,
                        &hour, &min, &sec, &pos);

    if (matched != 7) {
        return -1;  /* Parse failed */
    }

    /* Validate separator (must be space or 'T') */
    if (separator != ' ' && separator != 'T') {
        return -1;
    }

    /* Parse optional subseconds if present */
    if (str[pos] == '.') {
        const char *frac_start = &str[pos + 1];
        int digits = 0;
        long frac_value = 0;

        /* Accumulate up to 9 digits (nanosecond resolution). Beyond that the
         * value would overflow `long` on 32-bit platforms and the precision
         * is meaningless for double-precision epoch arithmetic anyway. */
        while (digits < 9 && frac_start[digits] >= '0' && frac_start[digits] <= '9') {
            frac_value = frac_value * 10 + (frac_start[digits] - '0');
            digits++;
        }

        if (digits > 0) {
            subseconds = frac_value / pow(10.0, digits);
        }
    }

    /* Validate ranges */
    if (year < 1970 || year > 2100 ||
        month < 1 || month > 12 ||
        day < 1 || day > 31 ||
        hour < 0 || hour > 23 ||
        min < 0 || min > 59 ||
        sec < 0 || sec > 59 ||
        subseconds < 0.0 || subseconds >= 1.0) {
        return -1;
    }

    /* Fill tm structure */
    tm_time.tm_year = year - 1900;  /* Years since 1900 */
    tm_time.tm_mon = month - 1;     /* Months since January (0-11) */
    tm_time.tm_mday = day;
    tm_time.tm_hour = hour;
    tm_time.tm_min = min;
    tm_time.tm_sec = sec;
    tm_time.tm_isdst = 0;           /* UTC has no DST */

    /* Convert to Unix epoch (UTC, no DST issues) */
    time_t epoch = timegm(&tm_time);
    if (epoch == -1) {
        return -1;  /* mktime failed (invalid date) */
    }

    /* timegm() normalizes overflowed fields in place (e.g. Feb 31 -> Mar 3),
     * so a calendar-impossible date would be silently accepted. Reject it by
     * checking the normalized struct still matches the requested date. */
    if (tm_time.tm_year != year - 1900 ||
        tm_time.tm_mon  != month - 1 ||
        tm_time.tm_mday != day) {
        return -1;
    }

    /* Combine integer seconds and subseconds */
    *epoch_seconds = (double)epoch + subseconds;

    return 0;
}

/* Free timestamp context */
void free_timestamp_context(TimestampContext *ctx)
{
    if (!ctx) {
        return;
    }

    if (ctx->original_timestamps) {
        for (int i = 0; i < ctx->n; i++) {
            free(ctx->original_timestamps[i]);
        }
        free(ctx->original_timestamps);
    }

    free(ctx);
}
