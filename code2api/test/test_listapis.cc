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