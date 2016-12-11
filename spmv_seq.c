/*
*********************************************
*  314 Principles of Programming Languages  *
*  Fall 2016                                *
*********************************************
*
*   Read a real (non-complex) sparse matrix from a Matrix Market (v. 2.0) file
*   and a vector from a txt file, perform matrix multiplication and store the
*   result to output.txt. This is the sequential version of sparse matrix vector
*   multiplication.
*
*
*
*   NOTES:
*
*   1) Matrix Market files are always 1-based, i.e. the index of the first
*      element of a matrix is (1,1), not (0,0) as in C.  ADJUST THESE
*      OFFSETS ACCORDINGLY offsets accordingly when reading and writing
*      to files.
*
*   2) ANSI C requires one to use the "l" format modifier when reading
*      double precision floating point numbers in scanf() and
*      its variants.  For example, use "%lf", "%lg", or "%le"
*      when reading doubles, otherwise errors will occur.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mmio.h"
#include <sys/time.h>
#include "utils.h"

//sorting according to the index
void quicksort(double* a, double* vindex, int* rindex, int* cindex, int n)
{
	int i, j, m;

	double p, t, s;
	if (n < 2)
		return;
	p = vindex[n / 2];

	for (i = 0, j = n - 1;; i++, j--) {
		while (vindex[i]<p)
			i++;
		while (p<vindex[j])
			j--;
		if (i >= j)
			break;
		t = a[i];
		a[i] = a[j];
		a[j] = t;

		s = vindex[i];
		vindex[i] = vindex[j];
		vindex[j] = s;

		m = rindex[i];
		rindex[i] = rindex[j];
		rindex[j] = m;

		m = cindex[i];
		cindex[i] = cindex[j];
		cindex[j] = m;
	}
	quicksort(a, vindex, rindex, cindex, i);
	quicksort(a + i, vindex + i, rindex + i, cindex + i, n - i);
}


int main(int argc, char *argv[])
{
	int ret_code;
	MM_typecode matcode;
	FILE *f;
	int M, N, nz;   //M is row number, N is column number and nz is the number of non-zero entries
	int i, j, vecdim, *rIndex, *cIndex, *rsIndex, *reIndex;
	double *val, *res, *vec, *vIndex;
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s [martix-market-filename] [input-vector-filename] \n", argv[0]);
		exit(1);
	}

	printf("\nOpening input matrix file: %s\n", argv[1]);
	if ((f = fopen(argv[1], "r")) == NULL)
	{
		printf("Fail to open the input matrix file!\n");
		exit(1);
	}
	if (mm_read_banner(f, &matcode) != 0)
	{
		printf("Could not process Matrix Market banner.\n");
		exit(1);
	}

	/*  This is how one can screen matrix types if their application */
	/*  only supports a subset of the Matrix Market data types.      */
	if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
		mm_is_sparse(matcode))
	{
		printf("Sorry, this application does not support ");
		printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
		exit(1);
	}

	/* find out size of sparse matrix .... */
	if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0)
		exit(1);

	/* reserve memory for matrices */
	rIndex = (int *)malloc(nz * sizeof(int));
	cIndex = (int *)malloc(nz * sizeof(int));
	val = (double *)malloc(nz * sizeof(double));

	/* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
	/*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
	/*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */
	for (i = 0; i<nz; i++)
	{
		fscanf(f, "%d %d %lg\n", &rIndex[i], &cIndex[i], &val[i]);
		rIndex[i]--;  /* adjust from 1-based to 0-based */
		cIndex[i]--;
	}

	if (f != stdin) fclose(f);

	//Open and load the input vector  
	printf("Opening input vector file: %s\n", argv[2]);
	if ((f = fopen(argv[2], "r")) == NULL)
	{
		printf("Fail to open the input vector file!\n");
		exit(1);
	}
	fscanf(f, "%d\n", &vecdim);
	if (vecdim != M)
	{
		printf("dimension mismatch!\n");
		exit(1);
	}
	vec = (double*)malloc(vecdim * sizeof(double));
	for (i = 0; i<vecdim; i++)
	{
		fscanf(f, "%lg\n", &vec[i]);
	}
	if (f != stdin) fclose(f);

	//the original calculation result
	double* res_seq = (double*)malloc(M*sizeof(double));
	memset(res_seq, 0, M*sizeof(double));

	getmul(val, vec, rIndex, cIndex, nz, res_seq);

	//preprocess the dataset to make the calculation can be parallelized
	vIndex = (double*)malloc(nz*sizeof(double));
	memset(vIndex, 0, nz*sizeof(double));
	for (i = 0; i < nz; i++)
	{
		vIndex[i] = (double)rIndex[i] * N + cIndex[i];
		if (vIndex[i] < 0)
		{	
	           printf("Error!\n");
	           exit(1);
	        }
	}

	quicksort(val, vIndex, rIndex, cIndex, nz);

	//we use rsIndex/reIndex to keep the start/end position of each row. The intial values are 
	//-1 for all entries.  rsIndex[i] indicates the start poistion of the i-th row. Hence 
	//the position index of the i-th row is from rsIndex[i] to reIndex[i]
	rsIndex = (int*)malloc(M*sizeof(int)); //start/end position of each row
	memset(rsIndex, -1, M*sizeof(int));

	reIndex = (int*)malloc(M*sizeof(int));
	memset(reIndex, -2, M*sizeof(int));

	for (i = 0; i<nz; i++)
	{
		int tmp = (int)(vIndex[i] / N);
		if (rsIndex[tmp] == -1)
		{
			rsIndex[tmp] = i;
			reIndex[tmp] = i;
		}
		else
			reIndex[tmp] = i;
	}

  printf("\n Start computation ... \n");
	struct timeval start, end;
	
  /**********************************************************/
	/* Start the spmarse matrix vector multiplication         */
	/**********************************************************/
	res = (double*)malloc(M*sizeof(double));
	memset(res, 0, M*sizeof(double));
	
	gettimeofday(&start, NULL);
	for (i = 0; i<M; i++)
	{
		for (j = rsIndex[i]; j <= reIndex[i]; j++)
		{
		  int tmp = cIndex[j];
			res[i] += val[j] * vec[tmp];
		}
	}
	gettimeofday(&end, NULL);

	
  printf(" End of computation ... \n \n");

	long elapsed_time = ((end.tv_sec * 1000000 + end.tv_usec)
		  - (start.tv_sec * 1000000 + start.tv_usec));
  
	
	if (!checkerror(res, res_seq, M))
	{
		printf("Calculation Error!\n");
		exit(1);
	}
	else {
		printf(" Test Result Passed ... \n");
	}


	printf(" Total time: %ld micro-seconds\n\n",  elapsed_time);

	// save the result
	if ((f = fopen("output.txt", "w")) == NULL)
	{
		printf("Fail to open the output file!\n");
		exit(1);
	}
	for (i = 0; i<vecdim; i++)
	{
		fprintf(f, "%lg\n", res[i]);
	}
	fclose(f);

	free(res_seq);
	free(vIndex);
	free(res);
	free(vec);
	free(rIndex);
	free(cIndex);
	free(val);
	free(rsIndex);
	free(reIndex);
	return 0;
}
