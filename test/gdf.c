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
        return (sscanf(line+len, "%le", &geoid->field)) ? HEADER : \
            error("can not parse %s\n", line); \
    }

    set_field(f , flatrefpot, 10)
    set_field(a2, radiusrefpot, 12)

    fprintf(stderr, "can not find field %s \n", line);
    return HEADER;
}

static int geoid_calc(geoid_t * geoid) {
    geoid->b2 = geoid->a2 * (1 - geoid->f);
    fprintf(stderr, "calc geoid a = %.16e b = %.16e f = %.16e\n", geoid->a2, geoid->b2, geoid->f);
    geoid->a2 *= geoid->a2;
    geoid->b2 *= geoid->b2;
    return BODY;
}

static int test_case_potential (const geoid_t * geoid, const char * line) {
    double lmb, phi, expected, s;
    if (sscanf (line, "%le %le %le", &lmb, &phi, &expected) != 3) {
        fprintf(stderr, "can not parse %s\n", line);
        return ERROR;
    }
    double cphi = cos(phi*M_PI/180.), sphi = sin(phi*M_PI/180.);
    double clmb = cos(lmb*M_PI/180.), slmb = sin(lmb*M_PI/180.);

    double alp = sqrt(geoid->a2*cphi*cphi + geoid->b2*sphi*sphi);
    vector_t r = { geoid->a2 * clmb * cphi, geoid->a2 * slmb * cphi, geoid->b2 * sphi };
    vector_scale(r, s, 1. / alp)

    double actual = potential(&gravity, &r);
    double delta = 0;
    if ((delta  = fabs((actual - expected) / (actual + expected))) > 1e-14) { 
        fprintf(stderr, " angles : %.16e, %.16e \n", lmb, phi);
        fprintf(stderr, " r : %.16e, %.16e, %.16e \n", r.x, r.y, r.z);
        fprintf(stderr, "delta = %+.16e actual = %+.16e expected = %+.16e\n", delta, actual, expected);
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
