#pragma once
#ifndef TRIANGLE_PLAQUETTE_HAMILTONIAN_INCLUDED
#define TRIANGLE_PLAQUETTE_HAMILTONIAN_INCLUDED
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

sparse_status_t triangle_plaquette_hamiltonian_matrix(sparse_matrix_t* const dest, int const nLayers, double const g, double const alpha) {
	sparse_status_t status = SPARSE_STATUS_SUCCESS; // stores the status of MKL function evaluations. 
	int nQubits = nLayers * 3;
	int pauliLength = nQubits * 3 + nLayers * 4 + 3;
	double* coefs = (double *)malloc(sizeof(double) * pauliLength);
	int** listOfPauliList = (int **)malloc(sizeof(int*) * pauliLength);
	int i, j;
	for (i = 0; i < pauliLength; ++i) {
		listOfPauliList[i] = (int *)malloc(sizeof(int) * nQubits);
	}
	int cnt = 0, s= 0;
	for (s = 0; s < nLayers; ++s) {
		// Electric terms
		for (i = 0; i < 3; ++i) {
			coefs[cnt] = g * g / 2;
			for (j = 0; j < nQubits; ++j) {
				if (j == 3 * s + i || j == (3 * (s + 1) + i) % nQubits)
					listOfPauliList[cnt][j] = 3;
				else
					listOfPauliList[cnt][j] = 0;
			}
			cnt++;
		}
		// XX terms
		for (i = 0; i < 3; ++i) {
			coefs[cnt] = -alpha / (2 * g * g);
			for (j = 0; j < nQubits; ++j) {
				if (j == 3 * s + i || j == (3 * (s + 1) + i) % nQubits)
					listOfPauliList[cnt][j] = 1;
				else
					listOfPauliList[cnt][j] = 0;
			}
			cnt++;
		}
		// YY terms
		for (i = 0; i < 3; ++i) {
			coefs[cnt] = -alpha / (2 * g * g);
			for (j = 0; j < nQubits; ++j) {
				if (j == 3 * s + i || j == (3 * (s + 1) + i) % nQubits)
					listOfPauliList[cnt][j] = 2;
				else
					listOfPauliList[cnt][j] = 0;
			}
			cnt++;
		}

		// Plaqutte terms
		// +XXX
		coefs[cnt] = -1 / (8 * g * g);
		for (j = 0; j < nQubits; ++j) {
			if (j / 3 == s)
				listOfPauliList[cnt][j] = 1;
			else
				listOfPauliList[cnt][j] = 0;
		}
		cnt++;

		// -XYY
		coefs[cnt] = 1 / (8 * g * g);
		for (j = 0; j < nQubits; ++j) {
			if (j / 3 == s) {
				if (j % 3 == 0)
					listOfPauliList[cnt][j] = 1;
				else
					listOfPauliList[cnt][j] = 2;
			}
			else {
				listOfPauliList[cnt][j] = 0;
			}
		}
		cnt++;

		// -YXY
		coefs[cnt] = 1 / (8 * g * g);
		for (j = 0; j < nQubits; ++j) {
			if (j / 3 == s) {
				if (j % 3 == 1)
					listOfPauliList[cnt][j] = 1;
				else
					listOfPauliList[cnt][j] = 2;
			}
			else {
				listOfPauliList[cnt][j] = 0;
			}
		}
		cnt++;

		// -XYY
		coefs[cnt] = 1 / (8 * g * g);
		for (j = 0; j < nQubits; ++j) {
			if (j / 3 == s) {
				if (j % 3 == 2)
					listOfPauliList[cnt][j] = 1;
				else
					listOfPauliList[cnt][j] = 2;
			}
			else {
				listOfPauliList[cnt][j] = 0;
			}
		}
		cnt++;
	}
	// shift the matrix to PSD by adding identities
	coefs[cnt] = (g * g / 2) * (nQubits - 6 * (nLayers % 2));
	for (j = 0; j < nQubits; ++j)
		listOfPauliList[cnt][j] = 0;
	cnt++;
	coefs[cnt] = (alpha / (2 * g * g)) * (2 * nQubits);
	for (j = 0; j < nQubits; ++j)
		listOfPauliList[cnt][j] = 0;
	cnt++;
	coefs[cnt] = (1 / (2 * g * g)) * (nLayers);
	for (j = 0; j < nQubits; ++j)
		listOfPauliList[cnt][j] = 0;
	cnt++;


	// load the constructor of the parent class
	CALL_AND_CHECK_STATUS(pauli_hamiltonian_matrix(dest, nQubits, pauliLength, listOfPauliList, coefs), "Error during computing the pauli hamiltonian\n");

memory_free:
	free(coefs);
	for (int i = 0; i < 1; ++i) {
		free(listOfPauliList[i]);
	}
	free(listOfPauliList);
}

#endif