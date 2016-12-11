#ifndef UTILS_H
#define UTILS_H
#include <time.h>
typedef struct{
	time_t tv_sec; /* seconds */
	long tv_nsec; /* nanoseconds */
}timespect;
typedef enum
{
  true = 1, false = 0
}bool;
bool checkerror(const double *resp, const double *ress, int dim);
void getmul(const double *val, const double *vec, const int *rIndex,
	const int *cIndex, int nz, double *res);

#endif
