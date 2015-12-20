#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>

enum { BEGIN, HEADER, BODY_START, BODY, ERROR };

#define BUFFER 100

#define print_space(out) fprintf(out, "    ");

typedef struct {
    char * modelname;
    int max_degree;
    float earth_gravity_constant;
    float radius;
    struct {double c, s;} ** cs;
} model_t;

int error(const char * msg, ...) {
    va_list argptr;
    va_start(argptr, msg);
    vfprintf(stderr, msg, argptr);
    va_end(argptr);
    return ERROR;
}

const char * next_token(const char * line) {
    while (*line && isspace(*line)) line++; 
    while (*line && !isspace(*line)) line++; 
    return line;
}

static int header_set(model_t * model, const char * line) {
    if (!*line) return HEADER;

#define set_field(field, len, func) \
    if (strncmp(line, #field, len) == 0) { \
        line = next_token(line); \
        model->field = func(line); \
        return model->field == 0.0 ? error("can not parse %s %s\n", #field, line) : HEADER;\
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

    set_field(radius, 6, atof)
    set_field(max_degree, 10, atof)
    set_field(earth_gravity_constant, 22, atof)
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
    int n = atoi(line = next_token(line));
    if (n < 0 || model->max_degree < n) return error("wrong n %d %s\n", n, line);
    int m = atoi(line = next_token(line));
    if (m < 0 || model->max_degree < m) return error("wrong m %d %s\n", m, line);
    if (model->cs[n][m].c != 0 || model->cs[n][m].s != 0) return error("already set %d %d\n", n, m);
    model->cs[n][m].c = atof(line = next_token(line));
    model->cs[n][m].s = atof(line = next_token(line));
    return BODY;
}

static int model_print(const model_t * model, FILE * out) {
    fprintf (stderr, "model_print\n");
    int i, n, m;

    i = 0;
    fprintf(out,"static const complex_t %s_CS[] = {\n", model->modelname);
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
    fprintf(out,"static const double %s_K[] = {\n", model->modelname);
    for (m = 0; m <= model->max_degree; m++) {
        fprintf(out,"\n");
        print_space(out);
        if (m==0) {
            fprintf(out,"/*%7d %4d %4d */ %+.12e,\n", i++, 0, 0, 1.);
        } else {
            fprintf(out,"/*%7d %4d %4d */ %+.12e,\n", i++, m, m, 
                    sqrt((m + 1.5) / m) / (2*m - 1));
        }
        for (n = m + 1; n <= model->max_degree; n++) {
            print_space(out);
            fprintf(out,"/*%7d %4d %4d */ %+.12e,\n", i++, n, m, 
                    sqrt((2*n + 1.) * (n - m) / (2*n - 1.) / (n + m)));
        }
    }
    fprintf(out,"};\n\n");
    fflush(out);

    fprintf(out,"static const gravity_t MODEL_%s = {\n", model->modelname);
    print_space(out); fprintf(out,"%d,\n", model->max_degree);
    print_space(out); fprintf(out,"%.12e,\n", model->radius);
    print_space(out); fprintf(out,"%.12e,\n", model->earth_gravity_constant);
    print_space(out); fprintf(out,"%s_CS,\n", model->modelname);
    print_space(out); fprintf(out,"%s_K\n", model->modelname);
    fprintf(out,"};\n\n");
    fflush(out);
    return 0;
}

static int main_print(char ** names, int size, FILE * out) {
    fprintf(out, "const gravity_t * gravity_get(const char * name) {\n");
    int i;
    for (i=0; i<size; i++) {
        print_space(out);
        fprintf(out, "if (strcmp(\"%s\", name) == 0) return &MODEL_%s;\n", names[i], names[i]);
    }
    print_space(out);
    fprintf(out, "return NULL;\n");
    fprintf(out, "}\n\n");
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
        return strncmp(line, "gfc", 3) ? error("token \"gfc\" missed\n") : body_add(model, line);
    default:
        return error("unknown error\n");
    }
}

int main(size_t argc, const char ** argv) {
    if (argc < 3) return error("usage: %s <output.c> <resource> [<resource> ...]]\n", argv[0]);
    fprintf (stderr, "open %s\n", argv[1]);
    FILE * out = fopen(argv[1], "w");
    fprintf(out, "#include \"gravity.h\"\n\n");
    model_t model;
    char names[BUFFER][argc-2];
    char *pointers[argc-2];
    char buffer[BUFFER];
    int state = BEGIN;
    int i;
    for (i=2; i<argc; i++) {
        pointers[i-2] = names[i-2];
        fprintf (stderr, "open %s\n", argv[i]);
        FILE * file = fopen (argv[i], "r");
        model.modelname = names[i-2];
        model.cs = NULL;
        while ((NULL != fgets(buffer, BUFFER, file)) &&
            ((state = handle (&model, state, buffer)) != ERROR)) {
        }
        if (state != ERROR) model_print(&model, out);
        model_free(&model);
        fclose(file);
        fprintf (stderr, "close %s\n", argv[i]);
        if (state == ERROR) return error("stop with error\n");
    }
    main_print(pointers, argc-2, out);
    fclose(out);
    fprintf (stderr, "close %s\n", argv[1]);
    return 0;
}
