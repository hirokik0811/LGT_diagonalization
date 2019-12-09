#pragma once
#ifndef PAULIHAMILTONIAN_INCLUDED
#define PAULIHAMILTONIAN_INCLUDED
#include <mkl.h>
#include "pauli_operator.h"
#include <stdbool.h> 

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
	sparse_matrix_t temp = NULL; // a matrix to store the computation result temporarily
	bool is_P_destroyed = false, is_temp_destroyed = false;

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
		// Compute the Pauli term
		CALL_AND_CHECK_STATUS(pauli_operator_matrix(&P, nQubits, pauliList, coefs[i]), "Error when computing a pauli term\n"); // define the ith pauli term
		is_P_destroyed = false;

		// Add P to the sum
		MKL_Complex16 one = { 1.0, 0.0 };
		CALL_AND_CHECK_STATUS(mkl_sparse_z_add(SPARSE_OPERATION_NON_TRANSPOSE, *dest, one, P, &temp),
			"Error when addition of two matrices\n");
		is_temp_destroyed = false;

		// Destroy old dest
		status = mkl_sparse_destroy(*dest);
		if (status != SPARSE_STATUS_SUCCESS)
		{
			printf(" Error during MKL_SPARSE_DESTROY(*dest) \n"); fflush(0);
		}
		
		// Destroy P
		status = mkl_sparse_destroy(P);
		if (status != SPARSE_STATUS_SUCCESS)
		{
			printf(" Error during MKL_SPARSE_DESTROY(P) \n"); fflush(0);
		}
		else {
			is_P_destroyed = true;
		}

		// Copy the sum to dest
		CALL_AND_CHECK_STATUS(mkl_sparse_copy(temp, descr, dest),
			"Error when copying temp to dest\n");

		// Destroy temp
		status = mkl_sparse_destroy(temp);
		if (status != SPARSE_STATUS_SUCCESS)
		{
			printf(" Error during MKL_SPARSE_DESTROY(temp) \n"); fflush(0);
		}
		else {
			is_temp_destroyed = true;
		}
	}	

memory_free:
	free(pauliList);
	if (!is_P_destroyed) {
		status = mkl_sparse_destroy(P);
		if (status != SPARSE_STATUS_SUCCESS)
		{
			printf(" Error during MKL_SPARSE_DESTROY(P) \n"); fflush(0);
		}
	}
	if (!is_temp_destroyed) {
		status = mkl_sparse_destroy(temp);
		if (status != SPARSE_STATUS_SUCCESS)
		{
			printf(" Error during MKL_SPARSE_DESTROY(temp) \n"); fflush(0);
		}
	}
	return status;
}

#endif