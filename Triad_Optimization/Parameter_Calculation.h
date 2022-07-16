#pragma once
#include "Triad.h"

// See if this can be simplified without much loss of speed.
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
};

// Could be micro-optimised, but should really look into the assembly code for this to be done properly.
void calculate_Numerical_Parameters(input_Data input, TRIAD_Output &quaternions, numerical_Parameters &output) {
	input_Package vector_Pair;
	// Should add a control element to test whether the Triad algorithm will work (ie. 0 norm or collinear vectors).
	vector_Pair = { 1, input.model_Vectors[0], input.model_Vectors[1], input.real_Vectors[0], input.real_Vectors[1] };
	TRIAD(vector_Pair, quaternions);
	vector_Pair = { 2, input.model_Vectors[1], input.model_Vectors[2], input.real_Vectors[1], input.real_Vectors[2] };
	TRIAD(vector_Pair, quaternions);
	vector_Pair = { 3, input.model_Vectors[2], input.model_Vectors[0], input.real_Vectors[2], input.real_Vectors[0] };
	TRIAD(vector_Pair, quaternions);
	double quaternion4[4];
	double* quaternion_Array[4] = { quaternions.quaternion1, quaternions.quaternion2, quaternions.quaternion3, quaternion4 };
	calculate_Orthogonal_Quaternion(quaternions, quaternion4);
	// The Matrix0 should be symmetrical, thus it is possible to save some time its calculation.
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
	// Micro-management theory: Using 3 loops might save time, by eliminating the need to for the switch statement.
	double scalar_component;
	double vector_component[3];
	double vector_summand[3];
	double (*current_Matrix)[4];
	double *current_model_Vector, *current_real_Vector;
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
					output.matrices[n + 1][i][j] += input.real_Vectors[n][k] * (vector_summand[k] - scalar_component * quaternion_Array[j][k + 1] + quaternion_Array[j][0] * vector_component[k]);
				}
			}
		}
	}
}