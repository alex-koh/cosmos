#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>

#include "gravity.h"

enum {HEADER, BODY, ERROR };

#define BUFFER 100

typedef struct {
    double a2, b2, f;
} geoid_t;

int error(const char * msg, ...) {
    va_list argptr;
    va_start(argptr, msg);
    vfprintf(stderr, msg, argptr);
    va_end(argptr);
    return ERROR;
}

static int header_set(geoid_t * geoid, const char * line) {
    while (*line && isspace(*line)) line++;
    if (!*line) return HEADER;

#define set_field(field, field_name, len) \
    if (strncmp(line, #field_name, len) == 0) { \
        return (geoid->field = atof(line+len)) == 0.0 ? error("can not parse %s\n", line) : HEADER; \
    }

    set_field(f , flatrefpot, 10)
    set_field(a2, radiusrefpot, 12)

    fprintf(stderr, "can not find field %s \n", line);
    return HEADER;
}

static int geoid_calc(geoid_t * geoid) {
    geoid->b2 = geoid->a2 * (1 - geoid->f);
    fprintf(stderr, "calc geoid a = %.12e b = %.12e f = %.12e\n", geoid->a2, geoid->b2, geoid->f);
    geoid->a2 *= geoid->a2;
    geoid->b2 *= geoid->b2;
    return BODY;
}

static int test_case_potential (const geoid_t * geoid, const char * line) {
    while (*line && isspace(*line)) line++;
    double lmb = atof(line);
    while (*line && !isspace(*line)) line++;
    while (*line && isspace(*line)) line++;
    double phi = atof(line);
    while (*line && !isspace(*line)) line++;
    double expected = atof(line);

    double cphi = cos(phi*M_PI/180.), sphi = sin(phi*M_PI/180.);
    double clmb = cos(lmb*M_PI/180.), slmb = sin(lmb*M_PI/180.);

    double alp = sqrt(geoid->a2*cphi*cphi + geoid->b2*sphi*sphi);
    double x = geoid->a2 * clmb * cphi / alp;
    double y = geoid->a2 * slmb * cphi / alp;
    double z = geoid->b2 * sphi / alp;

    double actual = potential(&gravity, x, y, z);
    double delta = 0;
    if ((delta  = 2 * fabs((actual - expected) / (actual + expected))) > 1e-14) { // 1.6e-10
        fprintf(stderr, " angles : %.12e, %.12e \n", lmb, phi);
        fprintf(stderr, " r : %.12e, %.12e, %.12e \n", x, y, z);
        fprintf(stderr, "delta = %+.12e actual = %+.12e expected = %+.12e\n", delta, actual, expected);
        return ERROR;
    }
    return BODY;
}

int main(size_t argc, const char ** argv) {
    if (argc != 2) return error("usage: %s <resource>.gfc\n", argv[0]);
    geoid_t geoid;
    char line[BUFFER];
    int state = HEADER;
    fprintf (stderr, "open %s\n", argv[1]);
    FILE * file = fopen (argv[1], "r");
    while (NULL != fgets(line, BUFFER, file)) {
        switch(state) {
        case HEADER:
            state = strncmp (line, "end_of_head", 11) ? header_set (&geoid, line) : geoid_calc(&geoid);
            break;
        case BODY:
            state = test_case_potential (&geoid, line);
            break;
        default: 
            return error("unknown state\n");
        }
        if (state == ERROR) return error ("stop with error\n");
    }
    fclose (file);
    fprintf (stderr, "close %s\n", argv[1]);
    return 0;
}
