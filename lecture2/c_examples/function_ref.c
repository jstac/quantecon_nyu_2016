#include <stdio.h>

void f(double x, double *y_pointer)
{
    *y_pointer = 1.0 / x;
}

int main(void)
{
    double x = 4.0;
    double y;
    f(x, &y);
    printf("y = %f\n", y);
    return 0;
}
