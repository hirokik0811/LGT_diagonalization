#pragma once
#ifndef CONVERSION_TENSOR_INCLUDED
#define CONVERSION_TENSOR_INCLUDED
#include <stdio.h>
#include <math.h>
#include <mkl.h>
#include <stdbool.h>
#include "pauli_operator.h"
#include "trace.h"

#define CALL_AND_CHECK_STATUS(function, error_message) do { \
          if(function != SPARSE_STATUS_SUCCESS)             \
          {                                                 \
          printf(error_message); fflush(0);                 \
          status = function;                                       \
          goto memory_free;                                 \
          }                                                 \
} while(0)

int** cartesian_int(int const n, int const repeat) {
	int** prod = malloc(pow(n, repeat) * sizeof(int*));
	for (int i = 0; i < pow(n, repeat); ++i) {
		prod[i] = malloc(repeat * sizeof(int));
	}
	for (int i = 0; i < pow(n, repeat); ++i) {
		for (int j = 0; j < repeat; ++j) {
			int val = (int)fmod(i / pow(n, j), n);
			prod[i][repeat - j - 1] = val;
		}
	}
	return prod;
}
sparse_status_t conversion_tensor(sparse_matrix_t* const dest, int const nQubits) {
	// Compute the tensor Tr[P_k P_i P_l P_j]/d required to calculate the PTM from the QPM 

	sparse_matrix_t Pi = NULL, Pj = NULL, Pk = NULL, Pl = NULL; // four sigma matrices 
	sparse_matrix_t subProd1 = NULL, subProd2 = NULL; // temporary matrices to store the result of matrix multiplication
	sparse_matrix_t prod = NULL; // a matrix to store the matrix product

	sparse_status_t status = SPARSE_STATUS_SUCCESS; // stores the status of MKL function evaluations. 
	bool is_Pi_destroyed = true, is_Pj_destroyed = true, is_Pk_destroyed = true, is_Pl_destroyed = true, 
		is_subProd1_destroyed = true, is_subProd2_destroyed = true, is_prod_destroyed = true;

	int i, j, k, l;
	int K = (int)pow(4, nQubits); // number of sigma matrices to consider for each index of i, j, k, l


	MKL_INT nnzC = 0; // number of non-zero elements of matrix C
	MKL_INT n_all = K * K * K * K;
	MKL_INT n_rowsC = K * K, n_colsC = K * K; // shape of the output matrix C
	MKL_Complex16* valuesC = (MKL_Complex16*)mkl_malloc(sizeof(MKL_Complex16) * n_all, 64); // memory to store the values of C
	MKL_INT* rows_indxC = (MKL_INT*)mkl_malloc(sizeof(MKL_INT) * (n_rowsC+1), 64),
		   * cols_indxC = (MKL_INT*)mkl_malloc(sizeof(MKL_INT) * n_all, 64); // coordinates of non-zero elements of C (in coo format)
	sparse_matrix_t cooC = NULL; // storage of C in coo format

	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
	descr.mode = SPARSE_FILL_MODE_UPPER;
	descr.diag = SPARSE_DIAG_NON_UNIT;

	int** inds = cartesian_int(4, nQubits);
	for (i = 0; i < K; ++i) {
		// Compute P_i
		CALL_AND_CHECK_STATUS(pauli_operator_matrix(&Pi, nQubits, inds[i], 1.0), "Error constructing P_i\n");
		is_Pi_destroyed = false;
		for (j = 0; j < K; ++j) {
			// Compute P_j
			CALL_AND_CHECK_STATUS(pauli_operator_matrix(&Pj, nQubits, inds[j], 1.0), "Error constructing P_j\n");
			is_Pj_destroyed = false;
			rows_indxC[i*K + j] = nnzC;
			for (k = 0; k < K; ++k) {
				// Compute P_k
				CALL_AND_CHECK_STATUS(pauli_operator_matrix(&Pk, nQubits, inds[k], 1.0), "Error constructing P_k\n");
				is_Pk_destroyed = false;
				for (l = 0; l < K; ++l) {
					// Compute P_l
					CALL_AND_CHECK_STATUS(pauli_operator_matrix(&Pl, nQubits, inds[l], 1.0), "Error constructing P_l\n");
					is_Pl_destroyed = false;
					CALL_AND_CHECK_STATUS(mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, Pl, Pj, &subProd1), "Error P_l*P_j \n");
					is_subProd1_destroyed = false;

					CALL_AND_CHECK_STATUS(mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, Pi, subProd1, &subProd2), "Error P_i*(P_l*P_j)\n");
					is_subProd2_destroyed = false;
					CALL_AND_CHECK_STATUS(mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, Pk, subProd2, &prod), "Error P_k*(P_i*P_l*P_j)\n");
					is_prod_destroyed = false;
					MKL_Complex16 tr = trace_sparse_z_csr(prod);
					if (tr.real != 0. || tr.imag != 0.) {
						cols_indxC[nnzC] = k * K + l;
						MKL_Complex16 tr_normalized = { tr.real / pow(2, nQubits), tr.imag / pow(2, nQubits) };
						valuesC[nnzC] = tr_normalized;
						nnzC++;
					}
					
					
					// Destroy the prod
					status = mkl_sparse_destroy(prod);
					if (status != SPARSE_STATUS_SUCCESS)
					{
						printf(" Error during MKL_SPARSE_DESTROY(prod) \n"); fflush(0);
					}
					else {
						is_prod_destroyed = true;
					}
					// Destroy subProd1 and subProd2
					status = mkl_sparse_destroy(subProd2);
					if (status != SPARSE_STATUS_SUCCESS)
					{
						printf(" Error during MKL_SPARSE_DESTROY(subProd2) \n"); fflush(0);
					}
					else {
						is_subProd2_destroyed = true;
					}

					status = mkl_sparse_destroy(subProd1);
					if (status != SPARSE_STATUS_SUCCESS)
					{
						printf(" Error during MKL_SPARSE_DESTROY(subProd1) \n"); fflush(0);
					}
					else {
						is_subProd1_destroyed = true;
					}
					
					// Destroy P_l
					status = mkl_sparse_destroy(Pl);
					if (status != SPARSE_STATUS_SUCCESS)
					{
						printf(" Error during MKL_SPARSE_DESTROY(Pl) \n"); fflush(0);
					}
					else {
						is_Pl_destroyed = true;
					}
				}
				// Destroy P_k
				status = mkl_sparse_destroy(Pk);
				if (status != SPARSE_STATUS_SUCCESS)
				{
					printf(" Error during MKL_SPARSE_DESTROY(Pk) \n"); fflush(0);
				}
				else {
					is_Pk_destroyed = true;
				}
			}
			rows_indxC[i * K + j + 1] = nnzC;
			// Destroy P_j
			status = mkl_sparse_destroy(Pj);
			if (status != SPARSE_STATUS_SUCCESS)
			{
				printf(" Error during MKL_SPARSE_DESTROY(Pj) \n"); fflush(0);
			}
			else {
				is_Pj_destroyed = true;
			}
		}
		// Destroy P_i
		status = mkl_sparse_destroy(Pi);
		if (status != SPARSE_STATUS_SUCCESS)
		{
			printf(" Error during MKL_SPARSE_DESTROY(Pi) \n"); fflush(0);
		}
		else {
			is_Pi_destroyed = true;
		}
	}

	valuesC = (MKL_Complex16*)mkl_realloc(valuesC, sizeof(MKL_Complex16) * nnzC);
	cols_indxC = (MKL_INT*)mkl_realloc(cols_indxC, sizeof(MKL_INT) * nnzC);
	// create the tensor in the csr format
	CALL_AND_CHECK_STATUS(mkl_sparse_z_create_csr(dest, SPARSE_INDEX_BASE_ZERO, n_rowsC, n_colsC, rows_indxC, rows_indxC + 1, cols_indxC, valuesC), 
		"Error when creating the conversion tensor in csr format\n");


memory_free:
	if (!is_Pi_destroyed) {
		status = mkl_sparse_destroy(Pi);
		if (status != SPARSE_STATUS_SUCCESS)
		{
			printf(" Error during MKL_SPARSE_DESTROY(Pi) \n"); fflush(0);
		}
	}
	if (!is_Pj_destroyed) {
		status = mkl_sparse_destroy(Pj);
		if (status != SPARSE_STATUS_SUCCESS)
		{
			printf(" Error during MKL_SPARSE_DESTROY(Pj) \n"); fflush(0);
		}
	}
	if (!is_Pk_destroyed) {
		status = mkl_sparse_destroy(Pk);
		if (status != SPARSE_STATUS_SUCCESS)
		{
			printf(" Error during MKL_SPARSE_DESTROY(Pk) \n"); fflush(0);
		}
	}
	if (!is_Pl_destroyed) {
		status = mkl_sparse_destroy(Pl);
		if (status != SPARSE_STATUS_SUCCESS)
		{
			printf(" Error during MKL_SPARSE_DESTROY(Pl) \n"); fflush(0);
		}
	}
	if (!is_subProd1_destroyed) {
		status = mkl_sparse_destroy(subProd1);
		if (status != SPARSE_STATUS_SUCCESS)
		{
			printf(" Error during MKL_SPARSE_DESTROY(subProd1) \n"); fflush(0);
		}
	}
	if (!is_subProd2_destroyed) {
		status = mkl_sparse_destroy(subProd2);
		if (status != SPARSE_STATUS_SUCCESS)
		{
			printf(" Error during MKL_SPARSE_DESTROY(subProd2) \n"); fflush(0);
		}
	}

	if (!is_prod_destroyed) {
		status = mkl_sparse_destroy(prod);
		if (status != SPARSE_STATUS_SUCCESS)
		{
			printf(" Error during MKL_SPARSE_DESTROY(prod) \n"); fflush(0);
		}
	}

	for (i = 0; i < K; ++i) {
		free(inds[i]);
	}
	free(inds);
	return status;
}

#endif