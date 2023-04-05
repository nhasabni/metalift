#include "list.h"

// Expected target code - bzero(a, sizeof(int) * n)
extern "C" List<int> test_bzero(List<int> a) {
    List<int> b = newList<int>();
    for (int i = 0; i < listLength(a); i++)
        b = listAppend(b, 0);
    return b;
}