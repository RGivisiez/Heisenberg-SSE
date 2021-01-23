/* Computing log (1 + x) in C. */
#include <math.h>
#include <quadmath.h>
double log1px_ (double *x)
{
  return log1p (*x);
}
__float128 log1pxq_ (__float128 *x)
{
  return log1pq (*x);
}
