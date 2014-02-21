#include <cblas.h>
#include <assert.h>

void mypp_dsymv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
	      const int N, const double alpha, const double *A,
	      const int lda, const double *X, const int incX,
	      const double beta, double *Y, const int incY)
{
	// limited implementation 
	assert(Order==CblasRowMajor);
	assert(Uplo==CblasUpper);
	assert(N==lda);

	__builtin_prefetch (Y, 1, 3);
	__builtin_prefetch (X, 1, 3);
	int i,j;
	double temp, reg1, reg2;
	const double *pA, *pX;
	double* pY = Y;
	pA = A;
	pX = X;

	// y = beta*y
	for(i=0;i<lda;i++,pY+=incY)
		(*pY) = beta * (*pY);

	// reset pointers
	pY = Y;

	for(i=0;i<lda;i++,pA+=i,pY+=incY)
	{
		pX = X + i*incX;
		reg1 = (*pX++);
		(*pY) += alpha * (*pA++) * reg1;
		temp = 0.0;
		for(j=i+1;j<N;j++,pA++,pX+=incX)
		{
			reg2 = alpha * (*pA);
			temp += reg2 * (*pX);
			Y[j*incY] += reg2 * reg1;
		}
		(*pY) += temp;
	}
}
	
