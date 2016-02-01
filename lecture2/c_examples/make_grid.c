#include <stdio.h>

/* Creates a linear array.  */
void linspace(double *ls, double a, double b, int n)
{
    double step = (b - a) / (n - 1);
    int i;
    for (i = 0; i < n; i++) 
    {
        ls[i] = a;
        a += step;
    }
}

int main(void) 
{
    int n = 10;
    double grid[n];
    linspace(grid, 1, 2, n);
    int i;
    for (i = 0; i < n; i++) 
    {
        printf("%g\n", grid[i]);
    }
    return 0;
}


