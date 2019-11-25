
#ifndef KRONECKER_INCLUDED
#define KRONECKER_INCLUDED
#pragma once
#include <stdio.h>
#include <mkl.h>


// Compute a Kronecker product of two matrices A \otimes B = C.
// c_{(i-1)*nrowB+k, (j-1)ncolB+l} = a_{i,j} * b_{k, l}

// For sparse matrices in CSR format 
sparse_status_t kronecker_sparse_c_csr(sparse_matrix_t *C, sparse_matrix_t const A, sparse_matrix_t const B) {
	
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
	int cnt = 0; // count from 0 to nnzC to specify the location of a value to be stored. Used as computing C matrix. 
	// output locations for mkl_sparse_c_export_csr

	MKL_INT n_rowsA, n_rowsB, n_colsA, n_colsB,
		*rows_startA = NULL, *rows_startB = NULL, *rows_endA = NULL, *rows_endB = NULL, *col_indxA = NULL, *col_indxB = NULL;
	MKL_Complex8 *valuesA = NULL, *valuesB = NULL;

	// number of non-zero elements of matrices A, B, and C
	MKL_INT nnzA, nnzB, nnzC;

	// shape of the output matrix C
	MKL_INT n_rowsC, n_colsC;

	// memory to store the values of C
	MKL_Complex8* valuesC = NULL;

	// coordinates of non-zero elements of C (in coo format)
	MKL_INT *rows_indxC = NULL, *cols_indxC = NULL;

	// storage of C in coo format
	sparse_matrix_t cooC = NULL;

	sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;
	CALL_AND_CHECK_STATUS(mkl_sparse_c_export_csr(A, &indexing, &n_rowsA, &n_colsA, &rows_startA, &rows_endA, &col_indxA, &valuesA), "Error occurs in mkl_sparse_c_export_csr(A) with ID: %d\n", status);
	CALL_AND_CHECK_STATUS(mkl_sparse_c_export_csr(B, &indexing, &n_rowsB, &n_colsB, &rows_startB, &rows_endB, &col_indxB, &valuesB), "Error occurs in mkl_sparse_c_export_csr(B) with ID: %d\n", status);

	// count the number of non-zero elements in each matrix
	nnzA = rows_endA[n_rowsA - 1];
	nnzB = rows_endB[n_rowsB - 1];

	// allocate memories for the output matrix C. C will be stored in coo format first and converted to csr afterwards. 
	nnzC = nnzA * nnzB;
	n_rowsC = n_rowsA * n_rowsB;
	n_colsC = n_colsA * n_colsB;
	valuesC = (MKL_Complex8*)mkl_malloc(sizeof(MKL_Complex8) * nnzC, 64);
	rows_indxC = (MKL_INT*)mkl_malloc(sizeof(MKL_INT) * nnzC, 64);
	cols_indxC = (MKL_INT*)mkl_malloc(sizeof(MKL_INT) * nnzC, 64);

	// compute C
	for (int rowA = 0; rowA < n_rowsA; ++rowA) {
		for (int i = rows_startA[rowA]; i < rows_endA[rowA]; i++) {
			int colA = col_indxA[i];
			for (int rowB = 0; rowB < n_rowsB; ++rowB) {
				for (int j = rows_startB[rowB]; j < rows_endB[rowB]; j++) {
					int colB = col_indxB[j];
					//printf("a_%d%d x b_%d%d \n", rowA, colA, rowB, colB);
					valuesC[cnt] = { valuesA[i].real * valuesB[j].real - valuesA[i].imag * valuesB[j].imag, valuesA[i].real * valuesB[j].imag + valuesA[i].imag * valuesB[j].real };
					rows_indxC[cnt] = rowA * n_rowsB + rowB;
					cols_indxC[cnt] = colA * n_colsB + colB;
					//printf("%5.3f + %5.3f I at (%d %d) \n", valuesC[cnt].real, valuesC[cnt].imag, rows_indxC[cnt], cols_indxC[cnt]);
					cnt++;
				}
			}
		}
	}
	// create C in coo format
	CALL_AND_CHECK_STATUS(mkl_sparse_c_create_coo(&cooC, indexing, n_rowsC, n_colsC, nnzC, rows_indxC, cols_indxC, valuesC), "Error occurs in mkl_sparse_c_create_coo(C) with ID: %d\n", status);
	// convert C to csr format
	CALL_AND_CHECK_STATUS(mkl_sparse_convert_csr(cooC, SPARSE_OPERATION_NON_TRANSPOSE, C), "Error occurs in mkl_sparse_convert_csr(C) with ID: %d\n", status);

/* Deallocate memory */
memory_free :

	// delocate memories used for extracting data of matrices A and B.
	mkl_free(valuesC); 
	mkl_free(rows_indxC); 
	mkl_free(cols_indxC);

	//Release matrix handle and deallocate arrays for which we allocate memory ourselves.
	status = mkl_sparse_destroy(cooC);
	if (status != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error after MKL_SPARSE_DESTROY(cooC) \n"); fflush(0);
	}

    // return status
	return status;

}
#endif