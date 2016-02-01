#include <stdio.h>

double f(double x)
{
    return 1.0 / x;
}

int main(void)
{
    double x = 4.0;
    double y = f(x);
    printf("y = %f\n", y);
    return 0;
}
