#include "list.h"

extern "C" int test_listsum(List<int> a) {
    int sum = 0;
    for (int i = 0; i < listLength(a); ++i) {
        sum += listGet(a, i);
    }
    return sum;
}

extern "C" int test_listlen(List<int> a) {
    int len = 0;
    for (int i = 0; i < listLength(a); ++i) {
        len++;
    }
    return len;
}

// Expected target code - bzero(a, sizeof(int) * n)
extern "C" List<int> test_bzero(List<int> a) {
    List<int> b = newList<int>();
    for (int i = 0; i < listLength(a); i++)
        b = listAppend(b, listGet(a, i));
    return b;
}

// cblas_saxpy computes y = a*x + y
//
// below we return a new list z to map into pure function semantics
// z = a*x + y
extern "C" List<int> test_cblas1_saxpy(int a, List<int> x, List<int> y) {
    List<int> z = newList<int>();

    if (listLength(x) != listLength(y)) return z;

    for (int i = 0; i < listLength(x); i++)
        z = listAppend(z, a * listGet(x, i) + listGet(y, i));
    return z;
}