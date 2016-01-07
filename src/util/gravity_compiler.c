#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>

enum { BEGIN, HEADER, BODY_START, BODY, ERROR };

#define BUFFER 100

#define print_space(out) fprintf(out, "    ");

typedef struct {double c,s;} complex_t;

typedef struct {
    char modelname[BUFFER];
    int max_degree;
    double earth_gravity_constant;
    double radius;
    complex_t ** cs;
} model_t;

int error(const char * msg, ...) {
    va_list argptr;
    va_start(argptr, msg);
    vfprintf(stderr, msg, argptr);
    va_end(argptr);
    return ERROR;
}

static const char * next_token(const char * line) {
    while (*line && isspace(*line)) line++; 
    while (*line && !isspace(*line)) line++; 
    return line;
}

static int header_set(model_t * model, const char * line) {
    if (!*line) return HEADER;
    int r;

#define set_field(field, len, fmt) \
    if (strncmp(line, #field, len) == 0) { \
        return sscanf(line + len, fmt, &model->field) ? HEADER : \
            error("can not parse %s\n"); \
    }

    if (strncmp(line, "modelname", 9) == 0) {
        line = next_token(line);
        while (*line && isspace(*line)) line++; 
        const char * start = line;
        line = next_token(line);
        if (!*line) return error("can not find modelname\n");
        strncpy(model->modelname, start, line-start);
        model->modelname[line-start] = 0;
        return HEADER;
    }

    set_field(radius, 6, "%lf")
    set_field(max_degree, 10, "%d")
    set_field(earth_gravity_constant, 22, "%lf")
    fprintf(stderr, "can not find field %s \n", line);
    return HEADER;
}

static int body_start(model_t * model) {
    fprintf (stderr, "body_start\n");
    model->cs = calloc((model->max_degree + 1), sizeof(void*));
    if (model->cs == NULL) return error("can not alloc memory for cs\n");
    int i;
    for (i=0; i<=model->max_degree; i++) {
        model->cs[i] = calloc((model->max_degree + 1), 2*sizeof(double));
        if (model->cs[i] == NULL) return error("can not alloc memory for cs[i]\n");
    }
    return BODY;
}

static void model_free(model_t * model) {
    fprintf (stderr, "model_free\n");
    int i;
    for (i=0; i<=model->max_degree; i++) free(model->cs[i]);
    free(model->cs);
}

static int body_add(model_t * model, const char * line) {
    int n,m;
    complex_t cs;
    if (sscanf (line, "%d %d %le %le", &n, &m, &cs.c, &cs.s) != 4) 
        return error("can not parse %s\n", line);
    if (m*m < 0 || model->max_degree < m || model->max_degree < n) 
        return error("wrong n|m %d %s\n", m, line);
    if (model->cs[n][m].c != 0 || model->cs[n][m].s != 0) 
        return error("already set %d %d\n", n, m);
    model->cs[n][m] = cs;
    return BODY;
}

static int model_print(const model_t * model, const char * filename) {
    fprintf (stderr, "model_print\n");
    int i, n, m;

    fprintf (stderr, "open %s\n", filename);
    FILE * out = fopen(filename, "w");
    fprintf(out, "#include \"gravity.h\"\n\n");

    i = 0;
    fprintf(out,"static const complex_t gravity_%s_CS[] = {\n", model->modelname);
    for (m=0; m <= model->max_degree; m++) {
        fprintf(out,"\n");
        for (n=m; n <= model->max_degree; n++) {
            print_space(out);
            fprintf(out,"/*%7d %4d %4d */ {%+.12e, %+.12e},\n", i++, n, m, 
                    model->cs[n][m].c, model->cs[n][m].s);
        }
    }
    fprintf(out,"};\n\n");
    fflush(out);

    i = 0;
    fprintf(out,"static const double gravity_%s_K[] = {\n", model->modelname);
    for (m = 0; m <= model->max_degree; m++) {
        fprintf(out,"\n");
        print_space(out);
        if (m==0) {
            fprintf(out,"/*%7d %4d %4d */ %+.12e,\n", i++, 0, 0, 1.);
        } else {
            fprintf(out,"/*%7d %4d %4d */ %+.12e,\n", i++, m, m, 
                    sqrt((m + .5) / m) / (2*m - 1));
        }
        for (n = m + 1; n <= model->max_degree; n++) {
            print_space(out);
            fprintf(out,"/*%7d %4d %4d */ %+.12e,\n", i++, n, m, 
                    sqrt((2*n + 1.) * (n - m) / (2*n - 1.) / (n + m)));
        }
    }
    fprintf(out,"};\n\n");
    fflush(out);

    fprintf(out,"const gravity_t gravity = {\n");
    print_space(out); fprintf(out,"%d,\n",    model->max_degree);
    print_space(out); fprintf(out,"%.12le,\n", model->radius);
    print_space(out); fprintf(out,"%.12le,\n", model->earth_gravity_constant);
    print_space(out); fprintf(out,"gravity_%s_CS,\n", model->modelname);
    print_space(out); fprintf(out,"gravity_%s_K\n",   model->modelname);
    fprintf(out,"};\n\n");
    fclose(out);
    fprintf (stderr, "close %s\n", filename);
    return 0;
}

/**
 *  Handle {@code line} withing to specifed state and return next state
 * */
static int handle(model_t * model, int state, char * line) {
//    fprintf (stderr, "handler %d\n", state);
    while (*line && isspace(*line)) line++;
    switch(state) {
    case BEGIN:
        return *line ? BEGIN : HEADER;
    case HEADER:
        return strncmp(line, "key", 3) ? header_set(model, line) : BODY_START;
    case BODY_START:
        return strncmp(line, "end_of_head", 11) ? error("token \"end_of_head\" missed\n") : body_start(model);
    case BODY:
        return strncmp(line, "gfc", 3) ? error("token \"gfc\" missed\n") : body_add(model, line+3);
    default:
        return error("unknown error\n");
    }
}

int main(size_t argc, const char ** argv) {
    if (argc != 3) return error("usage: %s <output.c> <resource>\n", argv[0]);
    model_t model;
    char buffer[BUFFER];
    int state = BEGIN;
    fprintf (stderr, "open %s\n", argv[2]);
    FILE * file = fopen (argv[2], "r");
    model.cs = NULL;
    while ((NULL != fgets(buffer, BUFFER, file)) &&
        ((state = handle (&model, state, buffer)) != ERROR)) {
    }
    if (state != ERROR) 
        return model_print (&model, argv[1]);
    else 
        return error ("stop with error\n");
    model_free (&model);
    fclose (file);
    fprintf (stderr, "close %s\n", argv[2]);
    return 0;
}
