#pragma once
#include "Data_Types.h"
// Theoretically this header file could be micro-optimised to cut execution time, but ideally this requires an analysis of the compiled .exe file.

// Function to calculate the values of the derivatives of the Lagrange function with given variables/coordinates as input (passed in as the input_values array).
// The calculated values are saved to the output array. The numerical parameters are passed in the numerical_Parameters struct. Some of the reusable parameters are also written into this struct.
void calculate_Function_Values(double* input_values, numerical_Parameters &parameters, double* output) {
	// Calculation of the alpha4 parameter. Represents the coordinate before the orthogonal quaternion. Depends fully on the current variables' values.
	// The alpha4 function value and its 3 derivative values (at the current coordinates/variables) are stroed for future reuse.
	double alpha4 = 1;
	double alpha4_derivative;
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 3; j++) {
			alpha4 -= 2 * input_values[i] * input_values[j] * parameters.matrices[0][i][j];
		}
		alpha4 -= input_values[i] * input_values[i] * parameters.matrices[0][i][i];
	}
	alpha4 = sqrt(abs(alpha4));
	parameters.alpha4 = alpha4;
	// Calculation of the values of the 2 reusable sums.
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
	// Calculation of the values of the first 3 function values: differentiation with respect to the variables/coordinates.
	for (int k = 0; k < 3; k++) {
		// Calculation of the kth (k=1,2,3) derivative of the alpha4 function at the current coordinates/varibles.
		alpha4_derivative = 0;
		for (int i = 0; i < 3; i++) {
			alpha4_derivative += input_values[i] * parameters.matrices[0][k][i];
		}
		alpha4_derivative /= -alpha4;
		parameters.alpha4_derivatives[k] = alpha4_derivative;
		// Calculation of the kth derivative function values at the current coordinates/varibles.
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
	// Calculation of the 4th function value: differentiation with respect to Lagrange's coefficient.
	output[3] = 1;
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 4; j++) {
			output[3] -= 2 * input_values[i] * input_values[j] * parameters.matrices[0][i][j];
		}
		output[3] -= input_values[i] * input_values[i] * parameters.matrices[0][i][i];
	}
}

// Function to calculate the Jacobi matrix of the transformation (R^4 -> R^4) defined by the 4 derivative functions at the current coordinates/variables (passed in as the input_values array).
// The calculate values are saved to the output 2D array. The numerical parameters and the reusable coefficients are passed in the numerical_Parameters struct.
void calculate_Jacobi_Matrix(double* input_values, numerical_Parameters parameters, double (*output)[4]) {
	for (int k = 0; k < 3; k++) {
		// The Jacobi matrix is symmetic, which is used to cut execution time.
		// Calculation of the left upper corner 3x3 block.
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
		// Calculation of the right upper 3x1 block.
		output[3][k] = parameters.alpha4_derivatives[k] * parameters.multi_use_sum2;
		for (int i = 0; i < 4; i++) {
			output[3][k] += input_values[i] * parameters.matrices[0][i][k];
		}
		output[3][k] *= 2;
		output[k][3] = output[3][k];
	}
	output[3][3] = 0;
}

// Function to calculate the inverse of a symmetrical matrix using the adjugate matrix. The fact that the inverse matrix also is symmetric is used to cut execution time.
// The input matrix is passed as 2D array input; the output matrix is the 2D array output.
// This algorithm can be replaced with a more efficient algorithm based upon factorising the input matrix.
void calculate_Matrix_Inverse(double (*input)[4], double (*output)[4]) {
	// First the first row of the adjugate matrix is calculated, which is then used to calculate the determinant of the input matrix.
	// Once the value of the determinant is calculated, the first row of the conjugate is modified to become the first row of the inverse matrix.
	//  All the consequent elements of the matrix are calculated using the traditional formula.
	double determinant = 1;
	for (int i = 0; i < 4; i++) {
		for (int j = i + 1; j < 4; j++) {
			int i_ref[3], j_ref[3];
			int i_counter = 0, j_counter = 0;
			for (int k = 0; k < 4; k++) {
				if (k != i) {
					i_ref[i_counter] = k;
					i_counter += 1;
				}
				if (k != j) {
					j_ref[j_counter] = k;
					j_counter++;
				}
			}
			output[i][j] = (input[i_ref[0]][j_ref[0]] * (input[i_ref[1]][j_ref[1]] * input[i_ref[2]][j_ref[2]] - input[i_ref[1]][j_ref[2]] * input[i_ref[2]][j_ref[1]])
				- input[i_ref[0]][j_ref[1]] * (input[i_ref[1]][j_ref[0]] * input[i_ref[2]][j_ref[2]] - input[i_ref[1]][j_ref[2]] * input[i_ref[2]][j_ref[0]])
				+ input[i_ref[0]][j_ref[2]] * (input[i_ref[1]][j_ref[0]] * input[i_ref[2]][j_ref[1]] - input[i_ref[1]][j_ref[1]] * input[i_ref[2]][j_ref[0]])) * (pow(-1.0, i+j)) / determinant;
			output[j][i] = output[i][j];
		}
		int ref[3];
		int counter = 0;
		for (int k = 0; k < 4; k++) {
			if (k != i) {
				ref[counter] = k;
				counter += 1;
			}
		}
		output[i][i] = (input[ref[0]][ref[0]] * (input[ref[1]][ref[1]] * input[ref[2]][ref[2]] - input[ref[1]][ref[2]] * input[ref[2]][ref[1]])
			- input[ref[0]][ref[1]] * (input[ref[1]][ref[0]] * input[ref[2]][ref[2]] - input[ref[1]][ref[2]] * input[ref[2]][ref[0]])
			+ input[ref[0]][ref[2]] * (input[ref[1]][ref[0]] * input[ref[2]][ref[1]] - input[ref[1]][ref[1]] * input[ref[2]][ref[0]])) / determinant;
		// Calculatio of the determinant and modification of the first row & column of the adjugate matrix.
		if (i == 0) {
			determinant = input[0][0] * output[0][0] + input[0][1] * output[0][1] + input[0][2] * output[0][2] + input[0][3] * output[0][3];
			for (int k = 1; k < 4; k++) {
				output[0][k] /= determinant;
				output[k][0] /= determinant;
			}
			output[0][0] /= determinant;
		}
	}
}

// Whaba Loss Function is used as a criterion for the optimality of a given quaternion. This function returns a positive loss value.
double Wahba_Loss_Function(numerical_Parameters parameters, double* coordinates) {
	double loss_Value = 1;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			loss_Value -= coordinates[i] * coordinates[j] * (parameters.matrices[1][i][j] + parameters.matrices[2][i][j] + parameters.matrices[3][i][j]);
		}
	}
	return loss_Value;
}

// Newton's numerical method is used to solve a system of non-linear equations, consequently allowing to find the conditional minimum of Wahba's Loss Function.
// The numerical_Parameters struct contains the numerical data used to calculate the values of the derivative functions and the Jacobi matrix. The final variables/coorrdinates are stored in the output_coordinates array.
void Newtons_Algorithm(numerical_Parameters& parameters, double* output_coordinates) {
	// Calculation of the initial coordinates/variables' values. These are calculated using the loss values for each of the 3 approximative quaternions.
	// To do: this block of code should be replaced with a smarter function to calculate the starting parameters. For now the initial value of Lagrange's coefficient is set to 1.
	double loss_Array[3];
	for (int i = 0; i < 3; i++) {
		double coordinates[4] = { 0, 0, 0, 0 };
		coordinates[i] = 1;
		loss_Array[i] = Wahba_Loss_Function(parameters, coordinates);
	}
	double variables[4] = { 0, 0, 0, 1 };
	double sum = 1 / loss_Array[0] + 1 / loss_Array[1] + 1 / loss_Array[2];
	for (int i = 0; i < 3; i++) {
		variables[i] = 1 / (loss_Array[i] * sum);
	}
	// Main loop to realise the iterative process. First the derivative function values, the Jacobi matrix and its inverse are calculated. Then the variables/coordinates are updated.
	// To do: should be replaced with a while loop which is set to exit once the desired accuracy is achieved.
	for (int n = 0; n < 100; n++) {
		double function_values[4];
		double Jacobi_Matrix[4][4], Jacobi_Inverse[4][4];
		calculate_Function_Values(variables, parameters, function_values);
		calculate_Jacobi_Matrix(variables, parameters, Jacobi_Matrix);
		calculate_Matrix_Inverse(Jacobi_Matrix, Jacobi_Inverse);
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				variables[i] -= Jacobi_Inverse[i][j] * function_values[j];
			}
		}
	}
	// Calcualtion of the alpha4 value and the formation of the final result of the algorithm.
	for (int i = 0; i < 3; i++) {
		output_coordinates[i] = variables[i];
	}
	double alpha4 = 1;
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 3; j++) {
			alpha4 -= 2 * variables[i] * variables[j] * parameters.matrices[0][i][j];
		}
		alpha4 -= variables[i] * variables[i] * parameters.matrices[0][i][i];
	}
	output_coordinates[3] = sqrt(abs(alpha4));
}