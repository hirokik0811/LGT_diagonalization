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

	sparse_matrix_t subProd1 = NULL, subProd2 = NULL; // temporary matrices to store the result of matrix multiplication
	sparse_matrix_t prod = NULL; // a matrix to store the matrix product

	sparse_status_t status = SPARSE_STATUS_SUCCESS; // stores the status of MKL function evaluations. 
	bool is_Pi_destroyed = true, is_Pj_destroyed = true, is_Pk_destroyed = true, is_Pl_destroyed = true, 
		is_subProd1_destroyed = true, is_subProd2_destroyed = true, is_prod_destroyed = true;

	int i, j, k, l;
	const int K = (int)pow(4, nQubits); // number of sigma matrices to consider for each index of i, j, k, l


	MKL_INT nnzC = 0; // number of non-zero elements of matrix C
	const MKL_INT n_all = K * K * K * K;
	MKL_INT n_rowsC = K * K, n_colsC = K * K; // shape of the output matrix C
	MKL_Complex16* valuesC = (MKL_Complex16*)mkl_malloc(sizeof(MKL_Complex16) * n_all, 64); // memory to store the values of C
	MKL_INT* rows_indxC = (MKL_INT*)mkl_malloc(sizeof(MKL_INT) * (n_rowsC+1), 64),
		   * cols_indxC = (MKL_INT*)mkl_malloc(sizeof(MKL_INT) * n_all, 64); // coordinates of non-zero elements of C

	int** inds = cartesian_int(4, nQubits);
	sparse_matrix_t* PList = mkl_malloc(sizeof(sparse_matrix_t) * K, 64);
	for (i = 0; i < K; ++i) {
		sparse_matrix_t S = NULL; // temporary matrix to store P_i
		CALL_AND_CHECK_STATUS(pauli_operator_matrix(&S, nQubits, inds[i], 1.0), "Error constructing list of P_k\n");

		
		// Recover the Hermitian matrix from the upper triangle matrix S
		struct matrix_descr descr;
		descr.type = SPARSE_MATRIX_TYPE_HERMITIAN;
		descr.mode = SPARSE_FILL_MODE_UPPER;
		descr.diag = SPARSE_DIAG_NON_UNIT;
		CALL_AND_CHECK_STATUS(mkl_sparse_copy(S, descr, PList + i), "Error recovering the Hermitian from S\n");

		
		// Destroy  S
		status = mkl_sparse_destroy(S);
		if (status != SPARSE_STATUS_SUCCESS)
		{
			printf(" Error during MKL_SPARSE_DESTROY(S) \n"); fflush(0);
		}
	}
	
	for (i = 0; i < K; ++i) {
		//printf("i: %d\n", i);
		for (j = 0; j < K; ++j) {
			//printf(" j: %d\n", j);
			rows_indxC[i*K + j] = nnzC;
			for (k = 0; k < K; ++k) {
				//printf("  k: %d\n", k);
				for (l = 0; l < K; ++l) {
					//printf("   l: %d\n", l);
					CALL_AND_CHECK_STATUS(mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, PList[l], PList[j], &subProd1), "Error P_l*P_j \n");
					is_subProd1_destroyed = false;

					CALL_AND_CHECK_STATUS(mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, PList[i], subProd1, &subProd2), "Error P_i*(P_l*P_j)\n");
					is_subProd2_destroyed = false;
					CALL_AND_CHECK_STATUS(mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, PList[k], subProd2, &prod), "Error P_k*(P_i*P_l*P_j)\n");
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
				}
			}
			rows_indxC[i * K + j + 1] = nnzC;
		}
	}

	valuesC = (MKL_Complex16*)mkl_realloc(valuesC, sizeof(MKL_Complex16) * nnzC);
	cols_indxC = (MKL_INT*)mkl_realloc(cols_indxC, sizeof(MKL_INT) * nnzC);
	// create the tensor in the csr format
	CALL_AND_CHECK_STATUS(mkl_sparse_z_create_csr(dest, SPARSE_INDEX_BASE_ZERO, n_rowsC, n_colsC, rows_indxC, rows_indxC + 1, cols_indxC, valuesC), 
		"Error when creating the conversion tensor in csr format\n");


memory_free:
	for (i = 0; i < K; ++i) {
		status = mkl_sparse_destroy(*(PList+i));
		if (status != SPARSE_STATUS_SUCCESS)
		{
			printf(" Error during MKL_SPARSE_DESTROY(*(PList+i)) \n"); fflush(0);
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