#include<cblas.h>
#include<assert.h>

void my_dsymv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
	      const int N, const double alpha, const double *A,
	      const int lda, const double *X, const int incX,
	      const double beta, double *Y, const int incY)
{
	// limited implementation 
	assert(Order==CblasRowMajor);
	assert(Uplo==CblasUpper);
	assert(lda==N);

	int i,j;
	double temp;
	for(i=0;i<lda;i++)
	{			
		temp=0;
		for(j=0;j<N;j++)
		{
			if(j<i)
				temp += A[j*N+i]*X[j*incX];
			else
				temp += A[i*N+j]*X[j*incX];	
		}
		Y[i*incY] = alpha*temp + beta*Y[i*incY];
	}
}
