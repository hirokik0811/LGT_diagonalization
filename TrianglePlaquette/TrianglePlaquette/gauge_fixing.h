#pragma once
#ifndef GAUGE_FIXING
#define GAUGE_FIXING

#include <stdio.h>
#include <math.h>
#include <mkl.h>
#include "pauli_hamiltonian.h"

#define CALL_AND_CHECK_STATUS(function, error_message) do { \
          if(function != SPARSE_STATUS_SUCCESS)             \
          {                                                 \
          printf(error_message); fflush(0);                 \
          status = function;                                       \
          goto memory_free;                                 \
          }                                                 \
} while(0)


typedef struct indexed_double {
	int ind;
	double val;
} indexed_double;
/*
void quick_sort(indexed_double* const low, indexed_double* const high){
	if (low < high)
	{
		// pi is partitioning index, arr[pi] is now
		// at right place 
		indexed_double* pi = partition(low, high);

		quickSort(low, pi - 1);  // Before pi
		quickSort(pi + 1, high); // After pi
	}
}


// This function takes last element as pivot, places
// the pivot element at its correct position in sorted
// array, and places all smaller (smaller than pivot)
// to left of pivot and all greater elements to right
// of pivot 
indexed_double* partition(indexed_double* const low, indexed_double* const high) {
	// pivot (Element to be placed at right position)
	double pivot_val = (*high).val;

	indexed_double* smPtr = low - 1;  // Index of smaller element
	indexed_double* curPtr;
	for (curPtr = low; curPtr < high; curPtr++)
	{
		// If current element is smaller than the pivot
		if ((*curPtr).val < pivot_val)
		{
			smPtr++;    // increment index of smaller element

			// swap values at smPtr and curPtr
			indexed_double val_sm = *smPtr;
			*smPtr = *curPtr;
			*curPtr = val_sm;
		}
	}
	// swap values at i+1 and high
	indexed_double val_i1 = *(smPtr +1);
	*smPtr = *high;
	*high = val_i1;
	return smPtr + 1;
}
*/
sparse_status_t zero_gauge_block(sparse_matrix_t* const dest, sparse_matrix_t const H, int const nLayers) {
	/*
	Compute diagonal elements (i.e. eigenvalues) of the gauge operators J12, J23, J31. 
	Return the block Hamiltonian with J12^2+J23^2+J31^2 = 0
	*/
	const int nQubits = nLayers * 3;
	const int matDim = (int)pow(2, nQubits);
	const int pauliLength = nLayers * 2;
	sparse_matrix_t J12 = NULL, J23 = NULL, J31 = NULL; // gauge operators

	int i, j;
	sparse_status_t status = SPARSE_STATUS_SUCCESS; // stores the status of MKL function evaluations. 

	// To store the data of the gauge operators
	sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

	indexed_double* gauge_square_sum = NULL; // array of squared sums of the diagonal elements of the gauge operators
	int* indZeros = NULL; // indices of H corresponding to zero flux

	// list of coefficients (\pm 1)
	double* listCoefJ12 = malloc(pauliLength * sizeof(double*));
	double* listCoefJ23 = (double*)malloc(pauliLength * sizeof(double*));
	double* listCoefJ31 = (double*)malloc(pauliLength * sizeof(double*));

    // list of pauli matrices 0: I, 1: X, 2: Y, 3: Z
	int** listPauliJ12 = (int**)malloc(pauliLength * sizeof(int*));
	int** listPauliJ23 = (int**)malloc(pauliLength * sizeof(int*));
	int** listPauliJ31 = (int**)malloc(pauliLength * sizeof(int*));

	for (i = 0; i < 2*nLayers; ++i) {
		listPauliJ12[i] = (int*)calloc(nQubits, sizeof(int));
		listPauliJ23[i] = (int*)calloc(nQubits, sizeof(int));
		listPauliJ31[i] = (int*)calloc(nQubits, sizeof(int));
	}
	
	for (i = 0; i < nLayers; ++i) {
		// set coefficients
		listCoefJ12[i * 2] = 1.0; listCoefJ23[i * 2] = 1.0; listCoefJ31[i * 2] = 1.0;
		listCoefJ12[i * 2 + 1] = -1.0; listCoefJ23[i * 2 + 1] = -1.0; listCoefJ31[i * 2 + 1] = -1.0;

		// set pauli matrices
		for (j = 0; j < nQubits; ++j) {
			if (j == 0 + 3 * i) {
				listPauliJ12[i * 2][j] = 3; listPauliJ31[i * 2 + 1][j] = 3;
			}
			else {
				listPauliJ12[i * 2][j] = 0; listPauliJ31[i * 2 + 1][j] = 0;
			}
			if (j == 1 + 3 * i) {
				listPauliJ23[i * 2][j] = 3; listPauliJ12[i * 2 + 1][j] = 3;
			}
			else {
				listPauliJ23[i * 2][j] = 0; listPauliJ12[i * 2 + 1][j] = 0;
			}
			if (j == 2 + 3 * i) {
				listPauliJ31[i * 2][j] = 3; listPauliJ23[i * 2 + 1][j] = 3;
			}
			else {
				listPauliJ31[i * 2][j] = 0; listPauliJ23[i * 2 + 1][j] = 0;
			}
		}
	}

	// Compute the gauge operators
	CALL_AND_CHECK_STATUS(pauli_hamiltonian_matrix(&J12, nQubits, pauliLength, listPauliJ12, listCoefJ12), "Error during computing J12\n");
	CALL_AND_CHECK_STATUS(pauli_hamiltonian_matrix(&J23, nQubits, pauliLength, listPauliJ23, listCoefJ23), "Error during computing J23\n");
	CALL_AND_CHECK_STATUS(pauli_hamiltonian_matrix(&J31, nQubits, pauliLength, listPauliJ31, listCoefJ31), "Error during computing J31\n");

	// Export the gauge operators
	MKL_Complex16* valuesJ12 = NULL;
	MKL_INT* pointerB_J12 = NULL;
	MKL_INT* pointerE_J12 = NULL;
	MKL_INT* columns_J12 = NULL;
	MKL_INT n_rowsJ12, n_colsJ12;

	MKL_Complex16* valuesJ23 = NULL;
	MKL_INT* pointerB_J23 = NULL;
	MKL_INT* pointerE_J23 = NULL;
	MKL_INT* columns_J23 = NULL;
	MKL_INT n_rowsJ23, n_colsJ23;

	MKL_Complex16* valuesJ31 = NULL;
	MKL_INT* pointerB_J31 = NULL;
	MKL_INT* pointerE_J31 = NULL;
	MKL_INT* columns_J31 = NULL;
	MKL_INT n_rowsJ31, n_colsJ31;

	// Export J12
	CALL_AND_CHECK_STATUS(mkl_sparse_z_export_csr(J12, &indexing, &n_rowsJ12, &n_colsJ12, &pointerB_J12, &pointerE_J12, &columns_J12, &valuesJ12),
		"Error after MKL_SPARSE_Z_EXPORT_CSR  J12\n");
	// Export J23
	CALL_AND_CHECK_STATUS(mkl_sparse_z_export_csr(J23, &indexing, &n_rowsJ23, &n_colsJ23, &pointerB_J23, &pointerE_J23, &columns_J23, &valuesJ23),
		"Error after MKL_SPARSE_Z_EXPORT_CSR  J23\n");
	// Export J31
	CALL_AND_CHECK_STATUS(mkl_sparse_z_export_csr(J31, &indexing, &n_rowsJ31, &n_colsJ31, &pointerB_J31, &pointerE_J31, &columns_J31, &valuesJ31),
		"Error after MKL_SPARSE_Z_EXPORT_CSR  J31\n");
	
	// compute the squared sums
	gauge_square_sum = (indexed_double*)malloc(matDim * sizeof(indexed_double));
	indZeros = (int*)malloc(matDim * sizeof(int));
	int nZeros = 0;
	for (i = 0; i < matDim; ++i) {
		double val = pow(valuesJ12[i].real, 2.0) + pow(valuesJ23[i].real, 2.0) + pow(valuesJ31[i].real, 2.0);
		gauge_square_sum[i].val = i;
		gauge_square_sum[i].ind = i;
		if (val == 0.0) {
			indZeros[nZeros] = i;
			nZeros++;
		}
	}
	indZeros = (int*)realloc(indZeros, nZeros* sizeof(int));
	

	/*
	// sort the squared sums and get the permutation order
	quick_sort(gauge_square_sum, gauge_square_sum + matDim - 1);
	int* perm = (int*)malloc(matDim * sizeof(int));
	for (i = 0; i < matDim; ++i) {
		perm[i] = gauge_square_sum[i].ind;
	}
	*/

	// export H
	MKL_Complex16* valuesH = NULL;
	MKL_INT* pointerB_H = NULL;
	MKL_INT* pointerE_H = NULL;
	MKL_INT* columns_H = NULL;
	MKL_INT n_rowsH, n_colsH;
	CALL_AND_CHECK_STATUS(mkl_sparse_z_export_csr(H, &indexing, &n_rowsH, &n_colsH, &pointerB_H, &pointerE_H, &columns_H, &valuesH),
		"Error during exporting H");

	// construct the block Hamiltonian
	
	MKL_Complex16* values = (MKL_Complex16*)mkl_malloc(sizeof(MKL_Complex16) * (nZeros*nZeros), 64);
	MKL_INT* columns = (MKL_INT*)mkl_malloc(sizeof(MKL_INT) * (nZeros * nZeros), 64);
	MKL_INT* rowIndex = (MKL_INT*)mkl_malloc(sizeof(MKL_INT) * (nZeros+1), 64);

	int cnt = 0;
	for (i = 0; i < nZeros; ++i) {
		rowIndex[i] = cnt;
		MKL_INT* rowB = columns_H + pointerB_H[indZeros[i]], *rowE = columns_H + pointerE_H[indZeros[i]];
		MKL_INT* curPtr;
		// check if the row contains column indices in indZeros
		for (curPtr = rowB; curPtr < rowE; ++curPtr) {
			for (j = 0; j < nZeros; ++j) {
				if (*curPtr == indZeros[j]) {
					columns[cnt] = j;
					values[cnt] = *(valuesH + (curPtr - columns_H));
					cnt++;
				}
			}
		}
	}
	rowIndex[nZeros] = cnt;
	values = (MKL_Complex16*)mkl_realloc(values, sizeof(MKL_Complex16) * cnt);
	columns = (MKL_INT*)mkl_realloc(columns, sizeof(MKL_INT) * cnt);

	CALL_AND_CHECK_STATUS(mkl_sparse_z_create_csr(dest, SPARSE_INDEX_BASE_ZERO, nZeros, nZeros, rowIndex, rowIndex + 1, columns, values),
		"Error in MKL_SPARSE_Z_CREATE_CSR zeroH\n");

memory_free:
	free(listCoefJ12);
	free(listCoefJ23);
	free(listCoefJ31);
	for (i = 0; i < pauliLength; ++i) {
		free(listPauliJ12[i]); free(listPauliJ23[i]); free(listPauliJ31[i]);
	}
	free(listPauliJ12); free(listPauliJ23); free(listPauliJ31);
	free(gauge_square_sum); free(indZeros);
	status = mkl_sparse_destroy(J12);
	if (status != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error during MKL_SPARSE_DESTROY(J12) \n"); fflush(0);
	}
	status = mkl_sparse_destroy(J23);
	if (status != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error during MKL_SPARSE_DESTROY(J23) \n"); fflush(0);
	}
	status = mkl_sparse_destroy(J31);
	if (status != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error during MKL_SPARSE_DESTROY(J31) \n"); fflush(0);
	}

	return status;

}

#endif