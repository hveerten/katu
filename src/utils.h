#ifndef UTILS_H
#define UTILS_H
#pragma once

#include <string.h>

#define MEM_PREPARE(X,S,type) \
    do \
    {\
        posix_memalign((void **) &(X), 32, sizeof(type) * (S));\
        memset(X, '\0', sizeof(type) * (S));\
        \
    } while(0);

#endif /* end of include guard: UTILS_H */
