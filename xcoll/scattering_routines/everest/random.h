#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* generate a random value weighted within the normal (gaussian) distribution */
static double get_random_gauss(void)
{
  double x = (double)random() / RAND_MAX,
         y = (double)random() / RAND_MAX,
         z = sqrt(-2 * log(x)) * cos(2 * M_PI * y);
  return z;
}

/* aggregate 100k cycles and display */
// main(void) {
//   static long g[11], i = 0;
//   srandom(time(NULL));
//   while (i++ < 100000)
//     g[(long)floor(gauss()+0.5) + 5]++;
//   for (i = 0; i < 11; i++)
//     printf("%2ld: %ld\n", i, g[i]);
//   return 0;
// }


// #include <stdlib.h>
// #include <math.h>
// #ifndef M_PI
// #define M_PI 3.14159265358979323846
// #endif

// double random_normal();