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
	sparse_matrix_t Sum = NULL; // a matrix to store the sum
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
	printf("\n");
	CALL_AND_CHECK_STATUS(pauli_operator_matrix(&Sum, nQubits, pauliList, coefs[0]), "Error during computing 0th pauli term\n"); // compute the first pauli term
	//CALL_AND_CHECK_STATUS(mkl_sparse_copy(P, descr, &Sum), "Error during copying a pauli operator at the 0th site\n");

	for (i = 1; i < PauliLength; ++i) {
		for (j = 0; j < nQubits; ++j) {
			pauliList[j] = listOfPauliList[i][j];
		}
		CALL_AND_CHECK_STATUS(pauli_operator_matrix(&P, nQubits, pauliList, coefs[i]), "Error during computing %d th pauli term\n", i); // define the ith pauli term
		MKL_Complex16 one = { 1.0, 0.0 };
		CALL_AND_CHECK_STATUS(mkl_sparse_z_add(SPARSE_OPERATION_NON_TRANSPOSE, Sum, one, P, &Sum),
			"Error during addition of two pauli terms at %d th site\n", i);
	}	
	CALL_AND_CHECK_STATUS(mkl_sparse_copy(Sum, descr, dest), "Error during copying a pauli operator\n");

memory_free:

	free(pauliList);
	status = mkl_sparse_destroy(Sum);
	if (status != SPARSE_STATUS_SUCCESS)
	{
		printf(" Error during MKL_SPARSE_DESTROY(Sum) \n"); fflush(0);
	}

	return status;
}

#endif