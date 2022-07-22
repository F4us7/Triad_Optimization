#pragma once
#include "Triad.h"
#include <cmath>

// Function to calculate the quaternion orthogonal to approximation hyperplane in the R^4 space.
// This function may be rewritten more compactly using a determinant function, yet this may result in the loss of speed. Deserves to be looked into.
void calculate_Orthogonal_Quaternion(TRIAD_Output input, double* output) {
	output[0] = input.quaternion1[1] * (input.quaternion2[2] * input.quaternion3[3] - input.quaternion2[3] * input.quaternion3[2])
		- input.quaternion2[1] * (input.quaternion1[2] * input.quaternion3[3] - input.quaternion1[3] * input.quaternion3[2])
		+ input.quaternion3[1] * (input.quaternion1[2] * input.quaternion2[3] - input.quaternion1[3] * input.quaternion2[2]);
	output[1] = - input.quaternion1[0] * (input.quaternion2[2] * input.quaternion3[3] - input.quaternion2[3] * input.quaternion3[2])
		+ input.quaternion2[0] * (input.quaternion1[2] * input.quaternion3[3] - input.quaternion1[3] * input.quaternion3[2])
		- input.quaternion3[0] * (input.quaternion1[2] * input.quaternion2[3] - input.quaternion1[3] * input.quaternion2[2]);
	output[2] = input.quaternion1[0] * (input.quaternion2[1] * input.quaternion3[3] - input.quaternion2[3] * input.quaternion3[1])
		- input.quaternion2[0] * (input.quaternion1[1] * input.quaternion3[3] - input.quaternion1[3] * input.quaternion3[1])
		+ input.quaternion3[0] * (input.quaternion1[1] * input.quaternion2[3] - input.quaternion1[3] * input.quaternion2[1]);
	output[3] = - input.quaternion1[0] * (input.quaternion2[1] * input.quaternion3[2] - input.quaternion2[2] * input.quaternion3[1])
		+ input.quaternion2[0] * (input.quaternion1[1] * input.quaternion3[2] - input.quaternion1[2] * input.quaternion3[1])
		- input.quaternion3[0] * (input.quaternion1[1] * input.quaternion2[2] - input.quaternion1[2] * input.quaternion2[1]);
	// Scaling the calculated orthogonal quaternion such that it has a norm = 1.
	double norm = sqrt(pow(output[0], 2) + pow(output[1], 2) + pow(output[2], 2) + pow(output[3], 2));
	for (int i = 0; i < 4; i++) {
		output[i] /= norm;
	}
};

// Function to calculate the 64 parameters (stored in the 4x4x4 3D array output.matrices) which fully define the function systems used to calculate the optimal quaternion.
// The intial 3D vector projections used to calculate the matrices output.matrices[n] (n = 1,2,3) are passed in the input_Data struct.
// The 3 quaternions calculated from the intial 3D vector projections are stored in a TRIAD_output struct. 
// Theoretically this function could be micro-optimised to cut execution time, but ideally this requires an analysis of the compiled .exe file.
void calculate_Numerical_Parameters(input_Data input, TRIAD_Output &quaternions, numerical_Parameters &output) {
	input_Package vector_Pair;
	// Calculation of the 3 approximative quaternions from the intial 3D vector projections.
	// To do (before working on data): Add a control sequence to verify that the initial 3D vector projections are non-collinear and are not zero-vectors, before calling the TRIAD algorithm.
	vector_Pair = { 1, input.model_Vectors[0], input.model_Vectors[1], input.real_Vectors[0], input.real_Vectors[1] };
	TRIAD(vector_Pair, quaternions);
	vector_Pair = { 2, input.model_Vectors[1], input.model_Vectors[2], input.real_Vectors[1], input.real_Vectors[2] };
	TRIAD(vector_Pair, quaternions);
	vector_Pair = { 3, input.model_Vectors[2], input.model_Vectors[0], input.real_Vectors[2], input.real_Vectors[0] };
	TRIAD(vector_Pair, quaternions);
	// Calculation of the quaternion orthogonal to the hyperplane spaned by the 3 approximative quaternions.
	double quaternion4[4];
	double* quaternion_Array[4] = { quaternions.quaternion1, quaternions.quaternion2, quaternions.quaternion3, quaternion4 };
	calculate_Orthogonal_Quaternion(quaternions, quaternion4);
	// Calculation of output.matrices[0]. It is a symmetrical matrix, which is used to cut execution time.
	double running_sum;
	for (int i = 0; i < 4; i++) {
		for (int j = i + 1; j < 4; j++) {
			running_sum = 0;
			for (int k = 0; k < 4; k++) { running_sum += quaternion_Array[i][k] * quaternion_Array[j][k]; };
			output.matrices[0][i][j] = running_sum;
			output.matrices[0][j][i] = running_sum;
		}
		running_sum = 0;
		for (int k = 0; k < 4; k++) { running_sum += quaternion_Array[i][k] * quaternion_Array[i][k]; };
		output.matrices[0][i][i] = running_sum;
	}
	// Calculation of matrices output.matrices[n] (n = 1,2,3). It should be noted that these matrices theoretically mustn't be symmetrical.
	double scalar_component;
	double vector_component[3];
	double vector_summand[3];
	// Iterator n defines which matrix is calculated; i and j define the current destination in the selected matrix.
	// The accuracy coefficients are accounted for at this stage in the algorithm.
	for (int n = 0; n < 3; n++) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				scalar_component = -scalar_Product((quaternion_Array[i] + 1), input.model_Vectors[n]);
				vector_Product((quaternion_Array[i] + 1), input.model_Vectors[n], vector_component);
				for (int k = 0; k < 3; k++) {
					vector_component[k] += quaternion_Array[i][0] * input.model_Vectors[n][k];
				}
				vector_Product((quaternion_Array[j] + 1), vector_component, vector_summand);
				output.matrices[n + 1][i][j] = 0;
				for (int k = 0; k < 3; k++) {
					output.matrices[n + 1][i][j] += input.accuracies[n] * input.real_Vectors[n][k] * (vector_summand[k] - scalar_component * quaternion_Array[j][k + 1] + quaternion_Array[j][0] * vector_component[k]);
				}
			}
		}
	}
}