#include <cblas.h>
#include <assert.h>
#include <emmintrin.h>

void mymm_dsymv(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
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
	double temp , reg1;
	double* unpack;
	const double *pA=A, *pX=X;
	double* pY = Y;
	__m128d mm_beta = _mm_set1_pd(beta);
	__m128d mm_alpha = _mm_set1_pd(alpha);
	__m128d mm_reg1 __attribute__((aligned (16)));
	__m128d mm_reg2 __attribute__((aligned (16)));
	__m128d mm_temp __attribute__((aligned (16)));
	posix_memalign((void**)&unpack, 16, 2*sizeof(double));

	// y = beta*y
	for(i=0;i<lda-1;i+=2,pY+=(2*incY))
	{
		mm_reg1 = _mm_loadu_pd(pY);
		mm_reg1 = _mm_mul_pd(mm_reg1, mm_beta);
		_mm_storeu_pd( pY, mm_reg1);
	}
	for(;i<lda;i++,pY++)
		(*pY) = beta * (*pY);

	pY = Y;
	for(i=0;i<lda;i++,pA+=i,pY+=incY)
	{
		pX = X + i*incX;
		temp = 0.0;
		mm_temp = _mm_set1_pd(0.0);
		j=i;
		for(;j<N-1;j+=2,pA+=2,pX+=(2*incX))
		{
			mm_reg1 = _mm_loadu_pd(pA);
			mm_reg2 = _mm_loadu_pd(pX);
			mm_reg1 = _mm_mul_pd(mm_reg1, mm_reg2);
			mm_reg2 = _mm_mul_pd(mm_alpha, mm_reg1);
			mm_temp = _mm_add_pd(mm_temp, mm_reg2);
		}
		_mm_store_pd(unpack, mm_temp);
		temp+=unpack[0];
		temp+=unpack[1];
		for(;j<N;j++,pA++,pX+=incX)
			temp += alpha * (*pA) * (*pX);
		(*pY) += temp;
	}

	pA = A;
	pX = X;
	for(i=0;i<lda;i++,pA+=i,pX+=incX)
	{
		reg1 = (*pX);
		mm_reg1 = _mm_set1_pd(reg1);
		pA++;
		pY=Y+(i+1)*incY;
		j=i+1;
		for(;j<N-1;j+=2,pA+=2,pY+=(2*incY))
		{
			mm_reg2 = _mm_loadu_pd(pA);
			mm_reg2 = _mm_mul_pd(mm_reg2, mm_reg1);
			mm_reg2 = _mm_mul_pd(mm_alpha, mm_reg2);
			_mm_storeu_pd( pY, _mm_add_pd(_mm_loadu_pd(pY),mm_reg2));
		}
		for(;j<N;j++,pA++,pY+=incY)
			*pY += alpha * (*pA) * reg1;
	}
}
