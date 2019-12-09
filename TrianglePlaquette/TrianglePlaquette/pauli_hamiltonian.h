#pragma once
#ifndef PAULIHAMILTONIAN_INCLUDED
#define PAULIHAMILTONIAN_INCLUDED
#include <mkl.h>
#include "pauli_operator.h"

#define CALL_AND_CHECK_STATUS(function, error_message) do { \
          if(function != SPARSE_STATUS_SUCCESS)             \
          {                                                 \
          printf(error_message); fflush(0);                 \
          status = function;                                       \
          goto memory_free;                                 \
          }                                                 \
} while(0)

sparse_status_t pauli_hamiltonian_matrix(sparse_matrix_t* const dest, int const nQubits, int const PauliLength, int** const listOfPauliList, double* const coefs) {

	sparse_status_t status = SPARSE_STATUS_SUCCESS; // stores the status of MKL function evaluations. 
	sparse_matrix_t P = NULL; // a matrix to store the next pauli to be added

	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	descr.mode = SPARSE_FILL_MODE_UPPER;
	descr.diag = SPARSE_DIAG_NON_UNIT;

	int* pauliList = (int *)malloc(nQubits * sizeof(int));
	int i, j;
	for (j = 0; j < nQubits; ++j) {
		pauliList[j] = listOfPauliList[0][j];
	}
	printf("0 th Pauli Operator\n");
	CALL_AND_CHECK_STATUS(pauli_operator_matrix(dest, nQubits, pauliList, coefs[0]), "Error during computing 0th pauli term\n"); // compute the first pauli term

	for (i = 1; i < PauliLength; ++i) {
		printf("%d th Pauli Operator\n", i);
		for (j = 0; j < nQubits; ++j) {
			pauliList[j] = listOfPauliList[i][j];
		}
		CALL_AND_CHECK_STATUS(pauli_operator_matrix(&P, nQubits, pauliList, coefs[i]), "Error when computing a pauli term\n"); // define the ith pauli term
		MKL_Complex16 one = { 1.0, 0.0 };
		CALL_AND_CHECK_STATUS(mkl_sparse_z_add(SPARSE_OPERATION_NON_TRANSPOSE, *dest, one, P, dest),
			"Error when addition of two pauli terms\n");
	}	

memory_free:

	free(pauliList);
	status = mkl_sparse_destroy(P);
	if (status != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error during MKL_SPARSE_DESTROY(P) \n"); fflush(0);
	}

	return status;
}

#endif