#pragma once
#include "Data_Types.h"
#include <cmath>
// Could be micro-optimised, but should really look into the assembly code for this to be done properly.

void calculate_Function_Values(double* input_values, numerical_Parameters &parameters, double* output) {
	// Calculation of the alpha4 value.
	double alpha4 = 1;
	for (int i = 0; i < 2; i++) {
		for (int j = i + 1; j < 3; j++) {
			alpha4 -= 2 * input_values[i] * input_values[j] * parameters.matrices[0][i][j];
		}
		alpha4 -= input_values[i] * input_values[i] * parameters.matrices[0][i][i];
	}
	alpha4 = sqrt(alpha4);
	parameters.alpha4 = alpha4;
	double alpha4_derivative;
	// Calculation of the reusable sum value.
	double sum1 = 0, sum2 = 0;
	for (int i = 0; i < 4; i++) {
		for (int n = 1; n < 4; n++) {
			sum1 += (parameters.matrices[n][3][i] + parameters.matrices[n][i][3]);
		}
		sum1 *= input_values[i];
		sum2 += input_values[i] * (parameters.matrices[0][3][i]);
	}
	parameters.multi_use_sum1 = sum1;
	parameters.multi_use_sum2 = sum2;
	// Calculation of the 3 function values.
	for (int k = 0; k < 3; k++) {
		// Calcualtion of the value of the alpha4 derivative.
		alpha4_derivative = 0;
		for (int i = 0; i < 3; i++) {
			alpha4_derivative += input_values[i] * parameters.matrices[0][k][i];
		}
		alpha4_derivative /= -alpha4;
		parameters.alpha4_derivatives[k] = alpha4_derivative;
		// Calcualtion of the value of the f_k.
		double summand1 = 0, summand2 = 0;
		for (int i = 0; i < 4; i++) {
			for (int n = 1; n < 4; n++) {
				summand1 = parameters.matrices[n][k][i] + parameters.matrices[n][i][k];
			}
			summand1 *= input_values[i];
			summand2 += input_values[i] * parameters.matrices[0][k][i];
		}
		output[k] = summand1 + sum1 * alpha4_derivative + 2 * input_values[3] * (summand2  + sum2 * alpha4_derivative);
	}
	output[3] = 1;
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 4; j++) {
			output[3] -= 2 * input_values[i] * input_values[j] * parameters.matrices[0][i][j];
		}
		output[3] -= input_values[i] * input_values[i] * parameters.matrices[0][i][i];
	}
}

void calculate_Jacobi_Matrix(double* input_values, numerical_Parameters& parameters, double** output) {
	for (int k = 0; k < 3; k++) {
		for (int m = k + 1; m < 3; m++) {
			output[k][m] = (parameters.matrices[0][m][k] + parameters.alpha4_derivatives[m] * parameters.matrices[0][4][k] +
				(-1 / parameters.alpha4) * (parameters.multi_use_sum1 + parameters.multi_use_sum2) * (parameters.alpha4_derivatives[k] * parameters.alpha4_derivatives[m] + parameters.matrices[0][m][k]) +
				parameters.alpha4_derivatives[k] * (parameters.matrices[0][m][4] + parameters.alpha4_derivatives[m] * parameters.matrices[0][k][4])) * 2 * input_values[3];
			for (int n = 1; n < 4; n++) {
				output[k][m] += parameters.matrices[n][k][m] + parameters.matrices[n][m][k] + parameters.alpha4_derivatives[m] * (parameters.matrices[n][k][3] + parameters.matrices[n][3][k]) +
					parameters.alpha4_derivatives[k] * (parameters.matrices[n][3][m] + parameters.matrices[n][m][3] + 2 * parameters.alpha4_derivatives[m] * parameters.matrices[n][3][3]);
			}
			output[m][k] = output[k][m];
		}
		output[k][k] = (parameters.matrices[0][k][k] + parameters.alpha4_derivatives[k] * parameters.matrices[0][4][k] +
				(-1 / parameters.alpha4) * (parameters.multi_use_sum1 + parameters.multi_use_sum2) * (parameters.alpha4_derivatives[k] * parameters.alpha4_derivatives[k] + parameters.matrices[0][k][k]) +
				parameters.alpha4_derivatives[k] * (parameters.matrices[0][k][4] + parameters.alpha4_derivatives[k] * parameters.matrices[0][k][4])) * 2 * input_values[3];
		for (int n = 1; n < 4; n++) {
			output[k][k] += 2 * parameters.matrices[n][k][k] + parameters.alpha4_derivatives[k] * (parameters.matrices[n][k][3] + parameters.matrices[n][3][k]) +
				parameters.alpha4_derivatives[k] * (parameters.matrices[n][3][k] + parameters.matrices[n][k][3] + 2 * parameters.alpha4_derivatives[k] * parameters.matrices[n][3][3]);
		}
		output[3][k] = parameters.alpha4_derivatives[k] * parameters.multi_use_sum2;
		for (int i = 0; i < 4; i++) {
			output[3][k] += input_values[i] * parameters.matrices[0][i][k];
		}
		output[3][k] *= 2;
		output[k][3] = output[3][k];
	}
	output[3][3] = 0;
}

