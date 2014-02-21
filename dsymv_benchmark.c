#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <cblas.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>

typedef unsigned int uint;
const int steps = 5;
const double eps = 0.0001;

extern void my_dsymv(const enum CBLAS_ORDER , const enum CBLAS_UPLO ,
	      const int , const double , const double*,
	      const int , const double*, const int ,
	      const double , double*, const int );
extern void mypp_dsymv(const enum CBLAS_ORDER , const enum CBLAS_UPLO ,
	      const int , const double , const double*,
	      const int , const double*, const int ,
	      const double , double*, const int );
extern void mymm_dsymv(const enum CBLAS_ORDER , const enum CBLAS_UPLO ,
	      const int , const double , const double*,
	      const int , const double*, const int ,
	      const double , double*, const int );

int compare_results(const double* Y1, const double* Y2, uint size)
{
	int i;
	for(i=0;i<size;i++)
	{
		if(abs((*Y1)-(*Y2))>eps)
			return 0;
		Y1++;Y2++;
	}

	return 1;
}

void flush_cache()
{
	char* input = malloc(1<<26);
	char* output = malloc(1<<26);
	memcpy(output, input, 1<<26);
	free(input);
	free(output);
}

int main(int argc, char* argv[]) {
	const double alpha = 1.33;
	const double beta = 0.77;
	uint size = 256;
	struct timeval start, end;
	double *A, *X, *Y, *p, *Y2, *Y3;
	long int total_time;
	int i,j;
	if(argc==2)
		size = abs(atoi(argv[1]));
	srand(time(NULL));

	posix_memalign( (void**)&A, 16, size*size*sizeof(double));
	if(!A)
		exit(-2);
	memset(A, 0, size*size*sizeof(double));

	posix_memalign( (void**)&X, 16, size*sizeof(double));
	if(!X)
	{
		free(A);exit(-2);
	}
	memset(X, 0, size*sizeof(double));

	posix_memalign( (void**)&Y, 16, size*sizeof(double));
	if(!Y)
	{
		free(X);free(A);exit(-2);
	}
	posix_memalign( (void**)&Y2, 16, size*sizeof(double));
	if(!Y2)
	{
		free(Y);free(X);free(A);exit(-2);
	}
	posix_memalign( (void**)&Y3, 16, size*sizeof(double));
	if(!Y3)
	{
		free(Y2);free(Y);free(X);free(A);exit(-2);
	}

	// random data initialization
	for(i=0;i<size;i++)
	{
		p = A + i*size;
		for(j=i;j<size;j++)
			*(p+j) = ((double)(rand()%8192))/17.33;
	}
	for(i=0;i<size;i++)
		*(X+i) = ((double)(rand()%8192))/17.33;
	for(i=0;i<size;i++)
		*(Y+i) = ((double)(rand()%8192))/17.33;
	memcpy(Y2, Y, size*sizeof(double));
	memcpy(Y3, Y, size*sizeof(double));

	// BLAS library function
	printf("size %d: ", size);
	total_time = 0;
	for(i=0;i<steps;i++)
	{
		gettimeofday(&start, NULL);
		cblas_dsymv(CblasRowMajor, CblasUpper, size, 
			alpha, A, size, X, 1, beta, Y2, 1);
		gettimeofday(&end, NULL);
		total_time += (end.tv_sec * 1000000 + end.tv_usec) - 
		 		(start.tv_sec * 1000000 + start.tv_usec);
		flush_cache();
	}
	printf("%ld ", total_time/steps);
	total_time = 0;

	flush_cache();
	// my unoptimized version for cblas_dsymv
	for(i=0;i<steps;i++)
	{
		gettimeofday(&start, NULL);
		my_dsymv(CblasRowMajor, CblasUpper, size, 
			alpha, A, size, X, 1, beta, Y3, 1);
		gettimeofday(&end, NULL);
		total_time += (end.tv_sec * 1000000 + end.tv_usec) - 
		 		(start.tv_sec * 1000000 + start.tv_usec);
		flush_cache();
	}
	printf("%ld ", total_time/steps);
	total_time = 0;

#ifndef NDEBUG
	printf("compare Y2 & my_dsymv Y3 = %d\n", compare_results( Y2, Y3, size));
#endif

	memcpy(Y3, Y, size*sizeof(double));
	flush_cache();
	
	// my optimized version of cblas_dsymv
	for(i=0;i<steps;i++)
	{
		gettimeofday(&start, NULL);
		mypp_dsymv(CblasRowMajor, CblasUpper, size, 
			  alpha, A, size, X, 1, beta, Y3, 1);
		gettimeofday(&end, NULL);
		total_time += (end.tv_sec * 1000000 + end.tv_usec) - 
		 		(start.tv_sec * 1000000 + start.tv_usec);
		flush_cache();
	}
	printf("%ld ", total_time/steps);
	total_time = 0;

#ifndef NDEBUG
	printf("compare Y2 & mypp_dsymv Y3 = %d\n", compare_results( Y2, Y3, size));
#endif

	memcpy(Y3, Y, size*sizeof(double));
	flush_cache();

	// my sse version of cblas_dsymv
	for(i=0;i<steps;i++)
	{
		gettimeofday(&start, NULL);
		mymm_dsymv(CblasRowMajor, CblasUpper, size, 
			  alpha, A, size, X, 1, beta, Y3, 1);
		gettimeofday(&end, NULL);
		total_time += (end.tv_sec * 1000000 + end.tv_usec) - 
		 		(start.tv_sec * 1000000 + start.tv_usec);
		flush_cache();
	}
	printf("%ld ", total_time/steps);
	total_time = 0;

#ifndef NDEBUG
	printf("compare Y2 & mypp_dsymv Y3 = %d\n", compare_results( Y2, Y3, size));
#endif

//	memcpy(Y3, Y, size*sizeof(double));
	flush_cache();
	printf("\n");

	free(Y3);
	free(Y2);
	free(Y);
	free(X);
	free(A);
	return 0;
}
