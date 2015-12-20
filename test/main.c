#include "test_gravity.h"

#include <stdio.h>

#include <stdarg.h>

#define TEST(name, arg, params) \
    if (strcmp(#name, arg) == 0) { \
        fprintf(stderr, "start test %s\n", #name); \
        if (res = test_ ## name params) { \
            fprintf(stderr, "test %s fail with code %d\n", #name, res); \
            return res; \
        } \
        return 0; \
    }
int error(const char * msg, ...) {
    va_list argptr;
    va_start(argptr, msg);
    vfprintf(stderr, msg, argptr);
    va_end(argptr);
    return 1;
}

int main(int argc, const char ** argv) {
    fprintf(stderr, "START TESTS\n");
    if (argc == 1) return error("%s potential\n", argv[0]);

    int res = 0;
    
    TEST (potential, argv[1], ())

    const char * resources = argv[1];
    printf("resources : %s\n", resources);
    
    return 0;
}


