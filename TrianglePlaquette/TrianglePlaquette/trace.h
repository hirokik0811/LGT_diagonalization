#pragma once
#ifndef TRACE_INCLUDED
#define TRACE_INCLUDED
#include <stdio.h>
#include <mkl.h>
#include <string.h>


// Compute a trace of a matrix A.

// For sparse Hermitian matrices in CSR format. The output contains only upper triangle part. Destroys the two input matrices afterwards. 
MKL_Complex16 trace_sparse_z_csr(sparse_matrix_t A) {

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

	// compute the trace
	for (rowA = 0; rowA < n_rowsA; ++rowA) {
		for (i = rows_startA[rowA]; i < rows_endA[rowA]; i++) {
			int colA = col_indxA[i];
			if (colA == rowA) { // if the value is diagonal, get its value and add it to tr_re and tr_im
				tr_re += valuesA[i].real;
				tr_im += valuesA[i].imag;
			}
		}
	}
	MKL_Complex16 tr = { tr_re, tr_im };
memory_free:
	
	return tr;
}
#endif#pragma once
