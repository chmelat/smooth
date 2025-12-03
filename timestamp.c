/* timestamp.c - RFC3339-style timestamp parsing and conversion */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

    /* Try parsing with subseconds first */
    int matched = sscanf(str, "%d-%d-%d%c%d:%d:%d.%lf",
                        &year, &month, &day, &separator,
                        &hour, &min, &sec, &subseconds);

    if (matched >= 7) {
        /* Got at least date and time, subseconds optional */
        if (matched == 7) {
            subseconds = 0.0;  /* No subseconds provided */
        } else {
            /* Subseconds parsed as decimal part (e.g., 390 instead of 0.390) */
            /* Count digits to normalize */
            const char *dot_pos = strchr(str, '.');
            if (dot_pos) {
                int digits = 0;
                dot_pos++;
                while (*dot_pos >= '0' && *dot_pos <= '9') {
                    digits++;
                    dot_pos++;
                }
                /* Convert to fractional seconds (e.g., 390 with 3 digits -> 0.390) */
                subseconds /= pow(10.0, digits);
            }
        }

        /* Validate separator (must be space or 'T') */
        if (separator != ' ' && separator != 'T') {
            return -1;
        }
    } else {
        return -1;  /* Parse failed */
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
    tm_time.tm_isdst = -1;          /* Let mktime determine DST */

    /* Convert to Unix epoch */
    time_t epoch = mktime(&tm_time);
    if (epoch == -1) {
        return -1;  /* mktime failed (invalid date) */
    }

    /* Combine integer seconds and subseconds */
    *epoch_seconds = (double)epoch + subseconds;

    return 0;
}

/* Convert timestamp strings to relative time array */
TimestampContext* convert_timestamps_to_relative(
    char **timestamp_strings,
    int n,
    double **x_out,
    int *first_error_line)
{
    if (!timestamp_strings || n <= 0 || !x_out || !first_error_line) {
        return NULL;
    }

    *first_error_line = -1;

    /* Allocate context */
    TimestampContext *ctx = malloc(sizeof(TimestampContext));
    if (!ctx) {
        return NULL;
    }

    ctx->original_timestamps = NULL;
    ctx->reference_epoch = 0.0;
    ctx->n = 0;
    ctx->errors_encountered = 0;

    /* Allocate arrays for valid timestamps */
    ctx->original_timestamps = malloc(n * sizeof(char*));
    double *x_temp = malloc(n * sizeof(double));

    if (!ctx->original_timestamps || !x_temp) {
        free(ctx->original_timestamps);
        free(x_temp);
        free(ctx);
        return NULL;
    }

    /* Initialize all pointers to NULL for safe cleanup */
    for (int i = 0; i < n; i++) {
        ctx->original_timestamps[i] = NULL;
    }

    /* Parse timestamps and build arrays */
    int valid_count = 0;
    int reference_set = 0;

    for (int i = 0; i < n; i++) {
        double epoch;

        if (parse_timestamp(timestamp_strings[i], &epoch) != 0) {
            /* Invalid timestamp */
            ctx->errors_encountered++;
            if (*first_error_line == -1) {
                *first_error_line = i + 1;  /* Line numbers start at 1 */
            }
            continue;  /* Skip this timestamp */
        }

        /* Valid timestamp - store it */
        ctx->original_timestamps[valid_count] = strdup(timestamp_strings[i]);
        if (!ctx->original_timestamps[valid_count]) {
            /* Memory allocation failed - cleanup and abort */
            for (int j = 0; j < valid_count; j++) {
                free(ctx->original_timestamps[j]);
            }
            free(ctx->original_timestamps);
            free(x_temp);
            free(ctx);
            return NULL;
        }

        /* Set reference epoch from first valid timestamp */
        if (!reference_set) {
            ctx->reference_epoch = epoch;
            reference_set = 1;
        }

        /* Calculate relative time in seconds */
        x_temp[valid_count] = epoch - ctx->reference_epoch;
        valid_count++;
    }

    /* Check if we got any valid timestamps */
    if (valid_count == 0) {
        for (int i = 0; i < n; i++) {
            free(ctx->original_timestamps[i]);
        }
        free(ctx->original_timestamps);
        free(x_temp);
        free(ctx);
        return NULL;
    }

    /* Store final count */
    ctx->n = valid_count;

    /* Allocate final output array */
    *x_out = malloc(valid_count * sizeof(double));
    if (!*x_out) {
        for (int i = 0; i < valid_count; i++) {
            free(ctx->original_timestamps[i]);
        }
        free(ctx->original_timestamps);
        free(x_temp);
        free(ctx);
        return NULL;
    }

    /* Copy valid x values to output */
    memcpy(*x_out, x_temp, valid_count * sizeof(double));
    free(x_temp);

    return ctx;
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
