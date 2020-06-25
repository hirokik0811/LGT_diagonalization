#pragma once
#ifndef TRACE_INCLUDED
#define TRACE_INCLUDED
#include <stdio.h>
#include <mkl.h>
#include <string.h>


// Recover Hermitian from triangular matrix. 

// For sparse Hermitian matrices in CSR format. The input contains only upper triangle part. 
MKL_Complex16 triangle_to_hermitian_sparse_z_csr(sparse_matrix_t A) {

	/* To avoid constantly repeating the part of code that checks inbound SparseBLAS functions' status,
   use macro CALL_AND_CHECK_STATUS */
#define CALL_AND_CHECK_STATUS(function, error_message) do { \
          if(function != SPARSE_STATUS_SUCCESS)             \
          {                                                 \
          printf(error_message); fflush(0);                 \
          status = function;                                \
          goto memory_free;                                 \
          }                                                 \
    } while(0)

	sparse_status_t status = SPARSE_STATUS_SUCCESS; // stores the status of MKL function evaluations.  
	int rowA, colA, i; // indices to specify the diagonal values

	// output locations for mkl_sparse_c_export_csr
	MKL_INT n_rowsA, n_colsA;
	MKL_INT* rows_startA = NULL, * rows_endA = NULL, * col_indxA = NULL;
	MKL_Complex16* temp_valuesA = NULL; // temp locations to store matrix values
	MKL_Complex16* valuesA = NULL;

	// number of non-zero elements of matrices A, B, and C
	MKL_INT nnzA;

	sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

	// export A
	CALL_AND_CHECK_STATUS(mkl_sparse_z_export_csr(A, &indexing, &n_rowsA, &n_colsA, &rows_startA, &rows_endA, &col_indxA, &valuesA), "Error occurs when exporting A\n");

	// count the number of non-zero elements in each matrix
	nnzA = rows_endA[n_rowsA - 1];

	// Matrix B to store the hermitian of the off-diagonal elements of A (i.e. B is lower triangle without diagonal elements)
	MKL_INT nnzB = nnzA; // number of non-zero elements of matrix B
	MKL_INT n_rowsC = K * K, n_colsC = K * K; // shape of the output matrix C
	MKL_Complex16* valuesC = (MKL_Complex16*)mkl_malloc(sizeof(MKL_Complex16) * n_all, 64); // memory to store the values of C
	MKL_INT* rows_indxC = (MKL_INT*)mkl_malloc(sizeof(MKL_INT) * (n_rowsC + 1), 64),
		* cols_indxC = (MKL_INT*)mkl_malloc(sizeof(MKL_INT) * n_all, 64); // coordinates of non-zero elements of C (in coo format)
	sparse_matrix_t cooC = NULL; // storage of C in coo format
	int cnt = 0;

	// compute the trace
	for (rowA = 0; rowA < n_rowsA; ++rowA) {
		for (i = rows_startA[rowA]; i < rows_endA[rowA]; ++i) {
			if (col_indxA[i] != rowA) { // if the value A_{ij} is off-diagonal, get its value and store its conjugate to A_{ji}
				MKL_Complex16 val = valuesA[i];
				
				continue;
			}
		}
	}
	MKL_Complex16 tr = { tr_re, tr_im };
memory_free:

	return tr;
}
#endif#pragma once
