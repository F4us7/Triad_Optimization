#pragma once
#include "Data_Types.h"
#include <iostream>

// Function to calculate the values of the derivatives of the Lagrange function with given variables/coordinates as input (passed in as the input_values array).
// The calculated values are saved to the output array. The numerical parameters are passed in the numerical_Parameters struct. Some of the reusable parameters are also written into this struct.
void calculate_Function_Values(optimisation_Output& variables, numerical_Parameters& parameters, double* output) {
	// Calculation of the alpha4 parameter. Represents the coordinate before the orthogonal quaternion. Depends fully on the current variables' values.
	// The alpha4 function value and its 3 derivative values (at the current coordinates/variables) are stroed for future reuse.
	variables.coordinates[3] = 1;
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 3; j++) {
			variables.coordinates[3] -= 2 * variables.coordinates[i] * variables.coordinates[j] * parameters.matrices[0][i][j];
		}
		variables.coordinates[3] -= variables.coordinates[i] * variables.coordinates[i] * parameters.matrices[0][i][i];
	}
	variables.coordinates[3] = sqrt(abs(variables.coordinates[3]));
	// Calculation of the values of the 2 reusable sums.
	double sum = 0;
	for (int i = 0; i < 4; i++) {
		for (int n = 1; n < 4; n++) {
			sum += (parameters.matrices[n][3][i] + parameters.matrices[n][i][3]);
		}
		sum *= variables.coordinates[i];
	}
	parameters.multi_use_sum = sum;
	// Calculation of the values of the first 3 function values: differentiation with respect to the variables/coordinates.
	for (int k = 0; k < 3; k++) {
		// Calculation of the kth (k=1,2,3) derivative of the alpha4 function at the current coordinates/varibles.
		double alpha4_derivative = 0;
		for (int i = 0; i < 3; i++) {
			alpha4_derivative += variables.coordinates[i] * parameters.matrices[0][k][i];
		}
		alpha4_derivative /= -variables.coordinates[3];
		parameters.alpha4_derivatives[k] = alpha4_derivative;
		// Calculation of the kth derivative function values at the current coordinates/varibles.
		double summand1 = 0, summand2 = 0;
		for (int i = 0; i < 4; i++) {
			for (int n = 1; n < 4; n++) {
				summand1 = parameters.matrices[n][k][i] + parameters.matrices[n][i][k];
			}
			summand1 *= variables.coordinates[i];
			summand2 += variables.coordinates[i] * parameters.matrices[0][k][i];
		}
		output[k] = summand1 + sum * alpha4_derivative + 2 * variables.mu * (summand2 + alpha4_derivative);
	}
	// Calculation of the 4th function value: differentiation with respect to Lagrange's coefficient.
	output[3] = 1;
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 4; j++) {
			output[3] -= 2 * variables.coordinates[i] * variables.coordinates[j] * parameters.matrices[0][i][j];
		}
		output[3] -= variables.coordinates[i] * variables.coordinates[i] * parameters.matrices[0][i][i];
	}
}

// Function to calculate the Jacobi matrix of the transformation (R^4 -> R^4) defined by the 4 derivative functions at the current coordinates/variables (passed in as the input_values array).
// The calculate values are saved to the output 2D array. The numerical parameters and the reusable coefficients are passed in the numerical_Parameters struct.
void calculate_Jacobi_Matrix(optimisation_Output variables, numerical_Parameters parameters, double(*output)[4]) {
	for (int k = 0; k < 3; k++) {
		// The Jacobi matrix is symmetic, which is used to cut execution time.
		// Calculation of the left upper corner 3x3 block.
		double sum = 0;
		for (int m = k + 1; m < 3; m++) {
			output[k][m] = (parameters.matrices[0][m][k] + parameters.alpha4_derivatives[m] * parameters.matrices[0][3][k] +
				(-1 / variables.coordinates[3]) * (parameters.multi_use_sum + 1) * (parameters.alpha4_derivatives[k] * parameters.alpha4_derivatives[m] + parameters.matrices[0][m][k]) +
				parameters.alpha4_derivatives[k] * (parameters.matrices[0][m][3] + parameters.alpha4_derivatives[m] * parameters.matrices[0][k][3])) * 2 * variables.mu;
			/*sum += (parameters.matrices[0][m][k] + parameters.alpha4_derivatives[m] * parameters.matrices[0][3][k]) * 2 * variables.mu;
			sum += (-1 / variables.coordinates[3]);
			sum += parameters.alpha4_derivatives[k] * (parameters.matrices[0][m][3] + parameters.alpha4_derivatives[m] * parameters.matrices[0][k][3]) * 2 * variables.mu;
			output[k][m] = sum;*/
			for (int n = 1; n < 4; n++) {
				output[k][m] += parameters.matrices[n][k][m] + parameters.matrices[n][m][k] + parameters.alpha4_derivatives[m] * (parameters.matrices[n][k][3] + parameters.matrices[n][3][k]) +
					parameters.alpha4_derivatives[k] * (parameters.matrices[n][3][m] + parameters.matrices[n][m][3] + 2 * parameters.alpha4_derivatives[m] * parameters.matrices[n][3][3]);
			}
			output[m][k] = output[k][m];
		}
		output[k][k] = (parameters.matrices[0][k][k] + parameters.alpha4_derivatives[k] * parameters.matrices[0][3][k] +
			(-1 / variables.coordinates[3]) * (parameters.multi_use_sum + 1) * (parameters.alpha4_derivatives[k] * parameters.alpha4_derivatives[k] + parameters.matrices[0][k][k]) +
			parameters.alpha4_derivatives[k] * (parameters.matrices[0][k][3] + parameters.alpha4_derivatives[k] * parameters.matrices[0][k][3])) * 2 * variables.mu;
		for (int n = 1; n < 4; n++) {
			output[k][k] += 2 * parameters.matrices[n][k][k] + parameters.alpha4_derivatives[k] * (parameters.matrices[n][k][3] + parameters.matrices[n][3][k]) +
				parameters.alpha4_derivatives[k] * (parameters.matrices[n][3][k] + parameters.matrices[n][k][3] + 2 * parameters.alpha4_derivatives[k] * parameters.matrices[n][3][3]);
		}
		// Calculation of the right upper 3x1 block.
		output[3][k] = parameters.alpha4_derivatives[k];
		for (int i = 0; i < 4; i++) {
			output[3][k] += variables.coordinates[i] * parameters.matrices[0][i][k];
		}
		output[3][k] *= 2;
		output[k][3] = output[3][k];
	}
	output[3][3] = 0;
}

// Function to calculate the values of the derivatives of the Lagrange function with given variables/coordinates as input (passed in as the input_values array).
// The calculated values are saved to the output array. The numerical parameters are passed in the numerical_Parameters struct. Some of the reusable parameters are also written into this struct.
void calculate_ALT_Function_Values(optimisation_Output variables, numerical_Parameters& parameters, double* output) {
	// Calculation of the values of the first 3 function values: differentiation with respect to the variables/coordinates.
	for (int k = 0; k < 3; k++) {
		double summand1 = 0, summand2 = 0;
		for (int i = 0; i < 3; i++) {
			for (int n = 1; n < 4; n++) {
				summand1 = parameters.matrices[n][k][i] + parameters.matrices[n][i][k];
			}
			summand1 *= variables.coordinates[i];
			summand2 += variables.coordinates[i] * parameters.matrices[0][k][i];
		}
		output[k] = summand1 + 2 * variables.mu * summand2;
	}
	// Calculation of the 4th function value: differentiation with respect to Lagrange's coefficient.
	output[3] = 1;
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 3; j++) {
			output[3] -= 2 * variables.coordinates[i] * variables.coordinates[j] * parameters.matrices[0][i][j];
		}
		output[3] -= variables.coordinates[i] * variables.coordinates[i] * parameters.matrices[0][i][i];
	}
}

// Function to calculate the Jacobi matrix of the transformation (R^4 -> R^4) defined by the 4 derivative functions at the current coordinates/variables (passed in as the input_values array).
// The calculate values are saved to the output 2D array. The numerical parameters and the reusable coefficients are passed in the numerical_Parameters struct.
void calculate_ALT_Jacobi_Matrix(optimisation_Output variables, numerical_Parameters parameters, double(*output)[4]) {
	for (int k = 0; k < 3; k++) {
		// The Jacobi matrix is symmetic, which is used to cut execution time.
		// Calculation of the left upper corner 3x3 block.
		for (int m = k + 1; m < 3; m++) {
			output[k][m] = parameters.matrices[0][m][k] * 2 * variables.mu;
			for (int n = 1; n < 4; n++) {
				output[k][m] += parameters.matrices[n][k][m] + parameters.matrices[n][m][k];
			}
			output[m][k] = output[k][m];
		}
		output[k][k] = parameters.matrices[0][k][k] * 2 * variables.mu;
		for (int n = 1; n < 4; n++) {
			output[k][k] += 2 * parameters.matrices[n][k][k];
		}
		// Calculation of the right upper 3x1 block.
		output[3][k] = 0;
		for (int i = 0; i < 3; i++) {
			output[3][k] += variables.coordinates[i] * parameters.matrices[0][i][k];
		}
		output[3][k] *= 2;
		output[k][3] = output[3][k];
	}
	output[3][3] = 0;
}

// Function to calculate the inverse of a symmetrical matrix using the adjugate matrix. The fact that the inverse matrix also is symmetric is used to cut execution time.
// The input matrix is passed as 2D array input; the output matrix is the 2D array output.
// This algorithm can be replaced with a more efficient algorithm based upon factorising the input matrix.
void calculate_Matrix_Inverse(double(*input)[4], double(*output)[4]) {
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
				+ input[i_ref[0]][j_ref[2]] * (input[i_ref[1]][j_ref[0]] * input[i_ref[2]][j_ref[1]] - input[i_ref[1]][j_ref[1]] * input[i_ref[2]][j_ref[0]])) * (pow(-1.0, i + j)) / determinant;
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
		// Calculation of the determinant and modification of the first row & column of the adjugate matrix.
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

// Experimental function to sample different quaternions and select the "best" starting values for Newton's iteration.
int calculate_Initial_Values(numerical_Parameters parameters, TRIAD_Output quaternions, double* output_coordinates, double& loss) {
	// Calculation of the initial loss values and the search for the approximative quaternion with a minimum loss value.
	double min_initial_loss, loss_Array[9], coordinates_matrix[9][4] = {};
	int min_index = 0;
	coordinates_matrix[0][0] = 1;
	min_initial_loss = Wahba_Loss_Function(parameters, coordinates_matrix[0]);
	loss_Array[0] = min_initial_loss;
	for (int i = 1; i < 3; i++) {
		coordinates_matrix[i][i] = 1;
		loss_Array[i] = Wahba_Loss_Function(parameters, coordinates_matrix[i]);
		if (loss_Array[i] < min_initial_loss) {
			min_initial_loss = loss_Array[i];
			min_index = i;
		}
	}
	int order_array[3] = { 0, 1, 2 };
	order_array[min_index] = 0;
	order_array[0] = min_index;
	// Calculation of the mirrored quaternions' coordinates and loss values.
	coordinates_matrix[3][order_array[1]] = -1;
	coordinates_matrix[4][order_array[2]] = -1;
	for (int i = 0; i < 4; i++) {
		coordinates_matrix[3][order_array[0]] += quaternions.quaternion_array[order_array[0]][i] * quaternions.quaternion_array[order_array[1]][i];
		coordinates_matrix[4][order_array[0]] += quaternions.quaternion_array[order_array[0]][i] * quaternions.quaternion_array[order_array[2]][i];
	}
	coordinates_matrix[3][order_array[0]] *= 2;
	coordinates_matrix[4][order_array[0]] *= 2;
	loss_Array[3] = Wahba_Loss_Function(parameters, coordinates_matrix[3]);
	loss_Array[4] = Wahba_Loss_Function(parameters, coordinates_matrix[4]);
	// Calculation of the first triangulated quaternion.
	double sum = 1 / loss_Array[0] + 1 / loss_Array[1] + 1 / loss_Array[2];
	for (int i = 0; i < 3; i++) {
		coordinates_matrix[5][i] = 1 / (loss_Array[i] * sum);
	}
	double norm = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 3; j++) {
			norm += 2 * coordinates_matrix[5][i] * coordinates_matrix[5][j] * parameters.matrices[0][i][j];
		}
		norm += pow(coordinates_matrix[5][i], 2) * parameters.matrices[0][i][i];
	}
	norm = sqrt(norm);
	for (int i = 0; i < 4; i++) { coordinates_matrix[5][i] /= norm; }
	loss_Array[5] = Wahba_Loss_Function(parameters, coordinates_matrix[5]);
	// Calculation of the second triangulated quaternion.
	sum = 1 / loss_Array[order_array[0]] + 1 / loss_Array[order_array[1]] + 1 / loss_Array[4];
	coordinates_matrix[6][order_array[0]] = 1 / (loss_Array[order_array[0]] * sum) + coordinates_matrix[1][order_array[0]] / (loss_Array[4] * sum);
	coordinates_matrix[6][order_array[1]] = 1 / (loss_Array[order_array[1]] * sum);
	coordinates_matrix[6][order_array[2]] = -1 / (loss_Array[4] * sum);
	norm = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 3; j++) {
			norm += 2 * coordinates_matrix[6][i] * coordinates_matrix[6][j] * parameters.matrices[0][i][j];
		}
		norm += pow(coordinates_matrix[6][i], 2) * parameters.matrices[0][i][i];
	}
	norm = sqrt(norm);
	for (int i = 0; i < 4; i++) { coordinates_matrix[6][i] /= norm; }
	loss_Array[6] = Wahba_Loss_Function(parameters, coordinates_matrix[6]);
	// Calculation of the third triangulated quaternion.
	sum = 1 / loss_Array[order_array[0]] + 1 / loss_Array[order_array[2]] + 1 / loss_Array[3];
	coordinates_matrix[7][order_array[0]] = 1 / (loss_Array[order_array[0]] * sum) + coordinates_matrix[0][order_array[0]] / (loss_Array[3] * sum);
	coordinates_matrix[7][order_array[2]] = 1 / (loss_Array[order_array[2]] * sum);
	coordinates_matrix[7][order_array[1]] = -1 / (loss_Array[3] * sum);
	norm = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 3; j++) {
			norm += 2 * coordinates_matrix[7][i] * coordinates_matrix[7][j] * parameters.matrices[0][i][j];
		}
		norm += pow(coordinates_matrix[7][i], 2) * parameters.matrices[0][i][i];
	}
	norm = sqrt(norm);
	for (int i = 0; i < 4; i++) { coordinates_matrix[7][i] /= norm; }
	loss_Array[7] = Wahba_Loss_Function(parameters, coordinates_matrix[7]);
	// Calculation of the fourth triangulated quaternion. 
	sum = 1 / loss_Array[order_array[0]] + 1 / loss_Array[3] + 1 / loss_Array[4];
	coordinates_matrix[8][order_array[0]] = 1 / (loss_Array[order_array[0]] * sum) + coordinates_matrix[0][order_array[0]] / (loss_Array[3] * sum) + coordinates_matrix[1][order_array[0]] / (loss_Array[4] * sum);
	coordinates_matrix[8][order_array[1]] = -1 / (loss_Array[3] * sum);
	coordinates_matrix[8][order_array[2]] = -1 / (loss_Array[4] * sum);
	norm = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < 3; j++) {
			norm += 2 * coordinates_matrix[8][i] * coordinates_matrix[8][j] * parameters.matrices[0][i][j];
		}
		norm += pow(coordinates_matrix[8][i], 2) * parameters.matrices[0][i][i];
	}
	norm = sqrt(norm);
	for (int i = 0; i < 4; i++) { coordinates_matrix[8][i] /= norm; }
	loss_Array[8] = Wahba_Loss_Function(parameters, coordinates_matrix[8]);
	// Selection of the "best" starting values out of all the samples.
	for (int i = 3; i < 9; i++) {
		if (loss_Array[i] < min_initial_loss) {
			min_initial_loss = loss_Array[i];
			min_index = i;
		}
	}
	for (int i = 0; i < 4; i++) {
		output_coordinates[i] = coordinates_matrix[min_index][i];
	}
	//std::cout << min_initial_loss << "\n";
	loss = min_initial_loss;
	return min_index;
}

optimisation_Output Newtons_Algorithm(numerical_Parameters& parameters, TRIAD_Output quaternions) {
	// Calculation of the initial coordinates/variables' values. These are calculated using the loss values for each of the 3 approximative quaternions.
	optimisation_Output result = {};
	result.initial_values_case = calculate_Initial_Values(parameters, quaternions, result.coordinates, result.loss_Value);
	result.mu = 0;
	// Main loop to realise the iterative process. First the derivative function values, the Jacobi matrix and its inverse are calculated. Then the variables/coordinates are updated.
	while (result.iteration_count < 1 || result.delta > 1.0e-6) {
		result.iteration_count++;
		double function_values[4];
		double Jacobi_Matrix[4][4], Jacobi_Inverse[4][4];
		/*if (result.iteration_count == 1 || result.coordinates[3] < 1.0e-9) {
			calculate_Alt_Function_Values(result, parameters, function_values);
			calculate_Alt_Jacobi_Matrix(result, parameters, Jacobi_Matrix);
		}
		else {
			calculate_Function_Values(result, parameters, function_values);
			calculate_Jacobi_Matrix(result, parameters, Jacobi_Matrix);
		}*/
		//if (result.initial_values_case < 3 && result.iteration_count == 1) {
		//	calculate_Alt_Function_Values(result, parameters, function_values);
		//	calculate_Alt_Jacobi_Matrix(result, parameters, Jacobi_Matrix);
		//}
		calculate_Function_Values(result, parameters, function_values);
		calculate_Jacobi_Matrix(result, parameters, Jacobi_Matrix);
		calculate_Matrix_Inverse(Jacobi_Matrix, Jacobi_Inverse);
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 4; j++) {
				result.coordinates[i] -= Jacobi_Inverse[i][j] * function_values[j];
			}
		}
		for (int j = 0; j < 4; j++) {
			result.mu -= Jacobi_Inverse[3][j] * function_values[j];
		}
		// Calcualtion of the alpha4 value and the formation of the final result of the algorithm.
		result.coordinates[3] = 1;
		for (int i = 0; i < 3; i++) {
			for (int j = i + 1; j < 3; j++) {
				result.coordinates[3] -= 2 * result.coordinates[i] * result.coordinates[j] * parameters.matrices[0][i][j];
			}
			result.coordinates[3] -= result.coordinates[i] * result.coordinates[i] * parameters.matrices[0][i][i];
		}
		result.coordinates[3] = sqrt(abs(result.coordinates[3]));
		double previous_Loss_Value = result.loss_Value;
		result.loss_Value = Wahba_Loss_Function(parameters, result.coordinates);
		result.delta = abs(previous_Loss_Value - result.loss_Value);
		//std::cout << result.loss_Value << "\n";
	}
	for (int i = 0; i < 4; i++) {
		result.optimal_quaternion[i] = result.coordinates[0] * quaternions.quaternion_array[0][i] + result.coordinates[1] * quaternions.quaternion_array[1][i]
			+ result.coordinates[0] * quaternions.quaternion_array[2][i] + result.coordinates[3] * quaternions.quaternion_array[3][i];
	}
	return result;
}