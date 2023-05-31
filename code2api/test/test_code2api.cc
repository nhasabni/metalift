#include "list.h"
#include "matrix.h"

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
//
// TODO: Currently not modelling incx, incy yet.
// https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-1/cblas-axpy.html
extern "C" List<int> test_cblas_saxpy(int a, List<int> x, List<int> y) {
    List<int> z = newList<int>();

    if (listLength(x) != listLength(y)) return z;

    for (int i = 0; i < listLength(x); i++)
        z = listAppend(z, a * listGet(x, i) + listGet(y, i));
    return z;
}

// cblas_sdot computes res = sum(x * y)
//
// TODO: Currently not modelling incx, incy yet.
// https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-1/cblas-dot.html#GUID-93DA36DC-40CA-4C01-B883-DABAB0D378D4
extern "C" int test_cblas_sdot(List<int> x, List<int> y) {
    int res = 0;

    if (listLength(x) != listLength(y)) return res;

    for (int i = 0; i < listLength(x); i++)
        res += listGet(x, i) * listGet(y, i);
    return res;
}

// cblas_sdot computes z = alpha * a * x + beta * y
//
// TODO: Currently not modelling transpose, row/col major yet.
// https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-1/cblas-gemv.html#GUID-25178576-05F1-4A33-8A0E-3694F0CCD242
extern "C" List<int> test_cblas_sgemv(int alpha, nestedList<int> a,
                                List<int>x, int beta, List<int> y) {
    List<int> z = newList<int>();

    if (nestedLength(a) <= 0 || listLength(nestedGet(a, 0)) <= 0) return z;

    // value of n
    if (listLength(nestedGet(a, 0)) != listLength(x)) return z;

    // value of m
    if (nestedLength(a) != listLength(y)) return z;

    int m = listLength(y);
    int n = listLength(x);

    // a is of size m * n, x is of size n x 1, y is of size m x 1
    for (int i = 0; i < m; i++) {
        int res = 0;
        for (int j = 0; j < n; j++) {
            // summation(a[i][j] * x[j]) over i=1 to m, j=1 to n
            res += listGet(nestedGet(a, i), j) * listGet(x, j);
        }
        z = listAppend(z, alpha * res + beta * listGet(y, i));
    }

    return z;
}