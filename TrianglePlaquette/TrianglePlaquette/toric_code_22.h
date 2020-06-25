#pragma once
#pragma once
#ifndef TORIC_CODE_22_INCLUDED
#define TORIC_CODE_22_INCLUDED
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

sparse_status_t toric_code_22_hamiltonian_matrix(sparse_matrix_t* const dest, int const nLayers, double const Je, double const Jm) {
	//
	// One layer contains 8 links
	// Structure: 
	// +---+---+
	// |   |   |
	// +---+---+
	// |   |   |
	// +---+---+
	
	//
	// Indexing:
	// A 0 B 4 A
	// 1   3   1
	// C 2 D 6 C
	// 5   7   5
	// A 0 B 4 A
	//
	//
	// E = -Je \sum_a \prod_s \prod_{j \in +_{a, s}} Z_j
	// B = -Jm \sum_s (\prod_{j \in \square}S+_j S+_{3s+1} + h.c.)
	// 

	sparse_status_t status = SPARSE_STATUS_SUCCESS; // stores the status of MKL function evaluations. 
	int nQubits = nLayers * 8;
	int pauliLength = 4+ nLayers * 4;
	double* coefs = (double*)calloc(pauliLength, sizeof(double));
	int** listOfPauliList = (int**)malloc(sizeof(int*) * pauliLength);
	int i, j;

	// initialize all values to zero with calloc
	for (i = 0; i < pauliLength; ++i) {
		listOfPauliList[i] = (int*)calloc(nQubits, sizeof(int));
	}
	int cnt = 0, s = 0;

	// Electric terms
	coefs[0] = -Je; coefs[1] = -Je; coefs[2] = -Je; coefs[3] = -Je;
	for (s = 0; s < nLayers; ++s) {
		listOfPauliList[0][0 + 8 * s] = 3; listOfPauliList[0][1 + 8 * s] = 3; listOfPauliList[0][4 + 8 * s] = 3; listOfPauliList[0][5 + 8 * s] = 3;
		listOfPauliList[1][0 + 8 * s] = 3; listOfPauliList[1][3 + 8 * s] = 3; listOfPauliList[1][4 + 8 * s] = 3; listOfPauliList[1][7 + 8 * s] = 3;
		listOfPauliList[2][1 + 8 * s] = 3; listOfPauliList[2][2 + 8 * s] = 3; listOfPauliList[2][5 + 8 * s] = 3; listOfPauliList[2][6 + 8 * s] = 3;
		listOfPauliList[3][2 + 8 * s] = 3; listOfPauliList[3][3 + 8 * s] = 3; listOfPauliList[3][6 + 8 * s] = 3; listOfPauliList[3][7 + 8 * s] = 3;
	}

	// Plaqutte terms
	for (s = 0; s < nLayers; ++s) {
		coefs[4 + 4 * s] = -Jm; coefs[5 + 4 * s] = -Jm; coefs[6 + 4 * s] = -Jm; coefs[7 + 4 * s] = -Jm;
		listOfPauliList[4 + 4 * s][0 + 8 * s] = 1; listOfPauliList[4 + 4 * s][1 + 8 * s] = 1; listOfPauliList[4 + 4 * s][2 + 8 * s] = 1; listOfPauliList[4 + 4 * s][3 + 8 * s] = 1;
		listOfPauliList[5 + 4 * s][1 + 8 * s] = 1; listOfPauliList[5 + 4 * s][3 + 8 * s] = 1; listOfPauliList[5 + 4 * s][4 + 8 * s] = 1; listOfPauliList[5 + 4 * s][6 + 8 * s] = 1;
		listOfPauliList[6 + 4 * s][0 + 8 * s] = 1; listOfPauliList[6 + 4 * s][2 + 8 * s] = 1; listOfPauliList[6 + 4 * s][5 + 8 * s] = 1; listOfPauliList[6 + 4 * s][7 + 8 * s] = 1;
		listOfPauliList[7 + 4 * s][4 + 8 * s] = 1; listOfPauliList[7 + 4 * s][5 + 8 * s] = 1; listOfPauliList[7 + 4 * s][6 + 8 * s] = 1; listOfPauliList[7 + 4 * s][7 + 8 * s] = 1;
	}

	// Compute the Hamiltonian of the model in Pauli basis
	CALL_AND_CHECK_STATUS(pauli_hamiltonian_matrix(dest, nQubits, pauliLength, listOfPauliList, coefs), "Error during computing the pauli hamiltonian\n");

memory_free:
	free(coefs);
	for (i = 0; i < 1; ++i) {
		free(listOfPauliList[i]);
	}
	free(listOfPauliList);

	return status;
}


typedef struct indexed_double {
	int ind;
	double val;
} indexed_double;

/*
sparse_status_t zero_gauge_block(sparse_matrix_t* const dest, sparse_matrix_t const H, int const nLayers) {
	
	//Compute diagonal elements (i.e. eigenvalues) of the gauge operators J12, J23, J31.
	//Return the block Hamiltonian with J12^2+J23^2+J31^2 = 0
	
	const int nQubits = nLayers * 12;
	const int matDim = (int)pow(2, nQubits);
	const int pauliLength = nLayers * 6;
	sparse_matrix_t JA = NULL, JB = NULL, JC = NULL, JD = NULL; // gauge operators

	int i, j, l;
	sparse_status_t status = SPARSE_STATUS_SUCCESS; // stores the status of MKL function evaluations. 

	// To store the data of the gauge operators
	sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

	indexed_double* gauge_square_sum = NULL; // array of squared sums of the diagonal elements of the gauge operators
	int* indZeros = NULL; // indices of H corresponding to zero flux

	// list of coefficients (\pm 1)
	double* listCoefJA = (double*)malloc(pauliLength * sizeof(double*));
	double* listCoefJB = (double*)malloc(pauliLength * sizeof(double*));
	double* listCoefJC = (double*)malloc(pauliLength * sizeof(double*));
	double* listCoefJD = (double*)malloc(pauliLength * sizeof(double*));

	// list of pauli matrices 0: I, 1: X, 2: Y, 3: Z
	int** listPauliJA = (int**)malloc(pauliLength * sizeof(int*));
	int** listPauliJB = (int**)malloc(pauliLength * sizeof(int*));
	int** listPauliJC = (int**)malloc(pauliLength * sizeof(int*));
	int** listPauliJD = (int**)malloc(pauliLength * sizeof(int*));

	// initialize all values to 0 (I). 
	for (i = 0; i < pauliLength; ++i) {
		listPauliJA[i] = (int*)calloc(nQubits, sizeof(int));
		listPauliJB[i] = (int*)calloc(nQubits, sizeof(int));
		listPauliJC[i] = (int*)calloc(nQubits, sizeof(int));
		listPauliJD[i] = (int*)calloc(nQubits, sizeof(int));
	}

	for (l = 0; l < nLayers; ++l) {
		// set coefficients
		for (i = 0; i < 3; ++i) {
			listCoefJA[l * 6 + i * 2] = 1.0; listCoefJB[l * 6 + i * 2] = 1.0; listCoefJC[l * 6 + i * 2] = 1.0; listCoefJD[l * 6 + i * 2] = 1.0;
			listCoefJA[l * 6 + i * 2 + 1] = -1.0; listCoefJB[l * 6 + i * 2 + 1] = -1.0; listCoefJC[l * 6 + i * 2 + 1] = -1.0; listCoefJD[l * 6 + i * 2 + 1] = -1.0;
		}

		// set pauli matrices
		listPauliJA[0][l * 12 + 1] = 3; listPauliJA[1][l * 12 + 2] = 3;
		listPauliJA[2][l * 12 + 5] = 3;
		listPauliJA[3][l * 12 + 3] = 3; listPauliJA[4][l * 12 + 9] = 3; listPauliJA[5][l * 12 + 10] = 3;
		listPauliJB[0][l * 12 + 3] = 3; listPauliJB[1][l * 12 + 4] = 3; listPauliJB[2][l * 12 + 7] = 3; listPauliJB[3][l * 12 + 8] = 3; listPauliJB[4][l * 12 + 11] = 3; listPauliJB[5][l * 12 + 9] = 3;
		listPauliJC[0][l * 12 + 2] = 3; listPauliJC[1][l * 12 + 0] = 3; listPauliJC[2][l * 12 + 6] = 3; listPauliJC[3][l * 12 + 7] = 3; listPauliJC[4][l * 12 + 10] = 3; listPauliJC[5][l * 12 + 11] = 3;
		listPauliJD[0][l * 12 + 0] = 3; listPauliJD[1][l * 12 + 1] = 3; listPauliJD[2][l * 12 + 4] = 3; listPauliJD[3][l * 12 + 5] = 3; listPauliJD[4][l * 12 + 8] = 3; listPauliJD[5][l * 12 + 6] = 3;

	}

	// Compute the gauge operators
	CALL_AND_CHECK_STATUS(pauli_hamiltonian_matrix(&JA, nQubits, pauliLength, listPauliJA, listCoefJA), "Error during computing JA\n");
	CALL_AND_CHECK_STATUS(pauli_hamiltonian_matrix(&JB, nQubits, pauliLength, listPauliJB, listCoefJB), "Error during computing JB\n");
	CALL_AND_CHECK_STATUS(pauli_hamiltonian_matrix(&JC, nQubits, pauliLength, listPauliJC, listCoefJC), "Error during computing JC\n"); \
		CALL_AND_CHECK_STATUS(pauli_hamiltonian_matrix(&JD, nQubits, pauliLength, listPauliJC, listCoefJD), "Error during computing JD\n");

	// Export the gauge operators
	MKL_Complex16* valuesJA = NULL;
	MKL_INT* pointerB_JA = NULL;
	MKL_INT* pointerE_JA = NULL;
	MKL_INT* columns_JA = NULL;
	MKL_INT n_rowsJA, n_colsJA;

	MKL_Complex16* valuesJB = NULL;
	MKL_INT* pointerB_JB = NULL;
	MKL_INT* pointerE_JB = NULL;
	MKL_INT* columns_JB = NULL;
	MKL_INT n_rowsJB, n_colsJB;

	MKL_Complex16* valuesJC = NULL;
	MKL_INT* pointerB_JC = NULL;
	MKL_INT* pointerE_JC = NULL;
	MKL_INT* columns_JC = NULL;
	MKL_INT n_rowsJC, n_colsJC;

	MKL_Complex16* valuesJD = NULL;
	MKL_INT* pointerB_JD = NULL;
	MKL_INT* pointerE_JD = NULL;
	MKL_INT* columns_JD = NULL;
	MKL_INT n_rowsJD, n_colsJD;

	// Export JA
	CALL_AND_CHECK_STATUS(mkl_sparse_z_export_csr(JA, &indexing, &n_rowsJA, &n_colsJA, &pointerB_JA, &pointerE_JA, &columns_JA, &valuesJA),
		"Error after MKL_SPARSE_Z_EXPORT_CSR  JA\n");
	// Export JB
	CALL_AND_CHECK_STATUS(mkl_sparse_z_export_csr(JB, &indexing, &n_rowsJB, &n_colsJB, &pointerB_JB, &pointerE_JB, &columns_JB, &valuesJB),
		"Error after MKL_SPARSE_Z_EXPORT_CSR  JB\n");
	// Export JC
	CALL_AND_CHECK_STATUS(mkl_sparse_z_export_csr(JC, &indexing, &n_rowsJC, &n_colsJC, &pointerB_JC, &pointerE_JC, &columns_JC, &valuesJC),
		"Error after MKL_SPARSE_Z_EXPORT_CSR  JC\n");
	// Export JD
	CALL_AND_CHECK_STATUS(mkl_sparse_z_export_csr(JD, &indexing, &n_rowsJD, &n_colsJD, &pointerB_JD, &pointerE_JD, &columns_JD, &valuesJD),
		"Error after MKL_SPARSE_Z_EXPORT_CSR  JD\n");

	// compute the squared sums
	gauge_square_sum = (indexed_double*)malloc(matDim * sizeof(indexed_double));
	indZeros = (int*)malloc(matDim * sizeof(int));
	int nZeros = 0;
	for (i = 0; i < matDim; ++i) {
		double val = pow(valuesJA[i].real, 2.0) + pow(valuesJB[i].real, 2.0) + pow(valuesJC[i].real, 2.0) + pow(valuesJD[i].real, 2.0);
		gauge_square_sum[i].val = val;
		gauge_square_sum[i].ind = i;
		if (val == 0.0) {
			indZeros[nZeros] = i;
			nZeros++;
		}
	}
	indZeros = (int*)realloc(indZeros, nZeros * sizeof(int));

	// export H
	MKL_Complex16* valuesH = NULL;
	MKL_INT* pointerB_H = NULL;
	MKL_INT* pointerE_H = NULL;
	MKL_INT* columns_H = NULL;
	MKL_INT n_rowsH, n_colsH;
	CALL_AND_CHECK_STATUS(mkl_sparse_z_export_csr(H, &indexing, &n_rowsH, &n_colsH, &pointerB_H, &pointerE_H, &columns_H, &valuesH),
		"Error during exporting H");

	// construct the block Hamiltonian

	MKL_Complex16* values = (MKL_Complex16*)mkl_malloc(sizeof(MKL_Complex16) * (nZeros * nZeros), 64);
	MKL_INT* columns = (MKL_INT*)mkl_malloc(sizeof(MKL_INT) * (nZeros * nZeros), 64);
	MKL_INT* rowIndex = (MKL_INT*)mkl_malloc(sizeof(MKL_INT) * (nZeros + 1), 64);

	int cnt = 0;
	for (i = 0; i < nZeros; ++i) {
		rowIndex[i] = cnt;
		MKL_INT* rowB = columns_H + pointerB_H[indZeros[i]], * rowE = columns_H + pointerE_H[indZeros[i]];
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
	free(listCoefJA);
	free(listCoefJB);
	free(listCoefJC);
	free(listCoefJD);
	for (i = 0; i < pauliLength; ++i) {
		free(listPauliJA[i]); free(listPauliJB[i]); free(listPauliJC[i]); free(listPauliJD[i]);
	}
	free(listPauliJA); free(listPauliJB); free(listPauliJC); free(listPauliJD);
	free(gauge_square_sum); free(indZeros);
	status = mkl_sparse_destroy(JA);
	if (status != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error during MKL_SPARSE_DESTROY(JA) \n"); fflush(0);
	}
	status = mkl_sparse_destroy(JB);
	if (status != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error during MKL_SPARSE_DESTROY(JB) \n"); fflush(0);
	}
	status = mkl_sparse_destroy(JC);
	if (status != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error during MKL_SPARSE_DESTROY(JC) \n"); fflush(0);
	}
	status = mkl_sparse_destroy(JD);
	if (status != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error during MKL_SPARSE_DESTROY(JD) \n"); fflush(0);
	}

	return status;

}

*/

#endif
