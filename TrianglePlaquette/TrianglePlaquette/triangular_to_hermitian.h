#pragma once
#ifndef TRACE_INCLUDED
#define TRACE_INCLUDED
#include <stdio.h>
#include <mkl.h>
#include <string.h>


// Recover the  Hermitian matrix from the upper triangular form in CSR format. 

// For sparse Hermitian matrices in CSR format. 
sparse_status_t triangular_to_hermitian_sparse_z_csr(sparse_matrix_t A) {

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
	double tr_re = 0, tr_im = 0; // real and imaginary parts of the trace

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


	// reallocate the arrays to have enough memory
	valuesA = (MKL_Complex16*)mkl_realloc(valuesA, sizeof(MKL_Complex16) * nnzA*2);
	col_indxA = (MKL_INT*)mkl_realloc(col_indxA, sizeof(MKL_INT) * nnzA * 2);

	int cnt = 0;
	// compute the trace
	for (rowA = 0; rowA < n_rowsA; ++rowA) {
		for (i = rows_startA[rowA]; i < rows_endA[rowA]; ++i) {
			if (col_indxA[i] != rowA) { // get the off-diagonal element A_{ij} and store its conjutate to A_{ji}
				MKL_Complex16 val = valuesA[i];
				MKL_Complex16 conj_val = { val.real, -val.imag };
				valuesA[cnt] = conj_val;
				cols_indxB[cnt] = rowA;
				cnt++;
			}
		}
	}
	
	

memory_free:
	
	return status;
}
#endif#pragma once
