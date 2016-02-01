
#include <stdio.h>
#include <math.h>

int main() {
    double x;
    int i = 3;
    x = 1.0 / (double) i;  /* cast i to double */
    double y = sin(x);
    printf("i = %d, x = %f, y = %f\n", i, x, y);
    return 0;
}
