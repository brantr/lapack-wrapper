#include"lapack_wrapper.h"
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<clapack/clapack.h>
#include<clapack/f2c.h>
/*! \fn void eigenvalue_decomposition(double **matrix, int nmatrix, double *&lambda, double **e_vectors)
 *  \brief Compute the eigenvalues of a symmetric. 
 *  real matrix via lapack.
 * 
 *  The eigenvalues are calculated and put in lambda. 
 *  The eigenvectors are calculated and put in the columns of e_vectors. */
void eigenvalue_decomposition(double **matrix, int nmatrix, double *lambda, double **e_vectors)
{
	integer n = nmatrix;
	integer lda = nmatrix;
	integer ldz = nmatrix;
	char uplo  = 'U';
	char compz = 'V';
	doublereal *A;		//symmetric matrix A
	doublereal *D;		//diagonal elements of tridiagonal matrix T
	doublereal *E;		//off-diagonal elemetns of tridiagonal matrix T
	doublereal *tau;	//scalar factors of elementary reflectors
	doublereal *ptr;
	doublereal *work;
	doublereal *z;

	integer  info;
	integer lgN = ((int)  log(n)/log(2)) + 1;
	integer  lwork  = 1 + 3*n + 2*n*lgN + 3*n*n;
	integer  liwork = 6 + 6*n + 5*n*lgN;
	integer  *iwork;

	int ntest = n*n;

	lwork  = n*ntest;

	A      = (doublereal *) malloc( lda*n*sizeof(doublereal) );
	D      = (doublereal *) malloc( n*sizeof(doublereal) );
	E      = (doublereal  *) malloc( (n-1)*sizeof(doublereal ) );
	tau    = (doublereal  *) malloc( n*sizeof(doublereal ) );
	work   = (doublereal *) malloc( lwork*sizeof(doublereal) );

	//first, store the symmetrical matrix into A
	ptr = A;
	for(int j=0; j<n; j++)
	{
		for(int i=0;i<ldz; i++)
		{
			*ptr++ = matrix[i][j];
		}
	}

	//perform the tridiagonal reduction
	dsytrd_(&uplo, &n, A, &lda, D, E, tau, work, &lwork, &info);


	//find the orthogonal matrix used in the reduction

	dorgtr_(&uplo, &n, A, &lda, tau, work, &lwork, &info);
	
	free(work);

	//allocate for eigenvalue decomposition

	lgN = ((int)  log(n)/log(2)) + 1;
	lwork  = 1 + 3*n + 2*n*lgN + 3*n*n;
	liwork = 6 + 6*n + 5*n*lgN;

	z      = A;
	work   = (doublereal *) malloc( lwork*sizeof(doublereal) );
	iwork  = (integer *)    malloc( liwork*sizeof(integer) );

	//perform eigenvalue decomposition
	dstedc_(&compz, &n, D, E, z, &ldz, work, &lwork, iwork, &liwork, &info);

	ptr = z;
	for(int j=0; j<n; j++)
	{
		lambda[j] = D[j];
		for(int i=0;i<ldz; i++)
		{
			e_vectors[i][j] = *ptr++;
		}
	}



	free(A);
	free(D);
	free(E);
	free(tau);
	free(work);
	free(iwork);
}
