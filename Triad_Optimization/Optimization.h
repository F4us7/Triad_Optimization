#pragma once
#include "Triad.h"

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

void calculate_Numerical_Parameters(input_Data input, TRIAD_Output* quaternions, numerical_Parameters* output) {
	input_Package vector_Pair;
	vector_Pair = { 1, input.model_vector1, input.model_vector2, input.real_vector1, input.real_vector2 };
	TRIAD(vector_Pair, quaternions);
	vector_Pair = { 2, input.model_vector2, input.model_vector3, input.real_vector2, input.real_vector3 };
	TRIAD(vector_Pair, quaternions);
	vector_Pair = { 3, input.model_vector3, input.model_vector1, input.real_vector3, input.real_vector1 };
	TRIAD(vector_Pair, quaternions);
	double quaternion4[4];
	double* quaternion_Array[4] = { quaternions->quaternion1, quaternions->quaternion2, quaternions->quaternion3, quaternion4 };
	calculate_Orthogonal_Quaternion((*quaternions), quaternion4);
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			double running_sum = 0;
			for (int k = 0; k < 4; k++) { running_sum += quaternion_Array[i][k] * quaternion_Array[j][k]; };
			(output->matrix0)[i][j] = running_sum;
		}
	}
	// Theory: Using 3 loops might save time, by eliminating the need to look up the matrix & vector pointers each time from the following reference arrays.
	//double* matrix_Array[3] = { output->matrix1[0], output->matrix2[0], output->matrix3[0] };
	//double* model_Vector_Array[3] = { input.model_vector1, input.model_vector2, input.model_vector3 };
	//double* real_Vector_Array[3] = { input.real_vector1, input.real_vector2, input.real_vector3 };
	//for (int n = 1; n < 4; n++) {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			double scalar_component = -scalar_Product((quaternion_Array[i] + 1), input.model_vector1);
			double vector_component[3];
			vector_Product((quaternion_Array[i] + 1), input.model_vector1, vector_component);
			for (int k = 0; k < 3; k++) {
				vector_component[k] += quaternion_Array[i][0] * input.model_vector1[k];
			}
			double vector_summand[3];
			vector_Product((quaternion_Array[j] + 1), vector_component, vector_summand);
			output->matrix1[i][j] = 0;
			for (int k = 0; k < 3; k++) {
				output->matrix1[i][j] += input.real_vector1[k] * (vector_summand[k] - scalar_component * quaternion_Array[j][k + 1] + quaternion_Array[j][0] * vector_component[k]);
			}
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			double scalar_component = -scalar_Product((quaternion_Array[i] + 1), input.model_vector2);
			double vector_component[3];
			vector_Product((quaternion_Array[i] + 1), input.model_vector2, vector_component);
			for (int k = 0; k < 3; k++) {
				vector_component[k] += quaternion_Array[i][0] * input.model_vector2[k];
			}
			double vector_summand[3];
			vector_Product((quaternion_Array[j] + 1), vector_component, vector_summand);
			output->matrix2[i][j] = 0;
			for (int k = 0; k < 3; k++) {
				output->matrix2[i][j] += input.real_vector2[k] * (vector_summand[k] - scalar_component * quaternion_Array[j][k + 1] + quaternion_Array[j][0] * vector_component[k]);
			}
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			double scalar_component = -scalar_Product((quaternion_Array[i] + 1), input.model_vector3);
			double vector_component[3];
			vector_Product((quaternion_Array[i] + 1), input.model_vector3, vector_component);
			for (int k = 0; k < 3; k++) {
				vector_component[k] += quaternion_Array[i][0] * input.model_vector3[k];
			}
			double vector_summand[3];
			vector_Product((quaternion_Array[j] + 1), vector_component, vector_summand);
			output->matrix3[i][j] = 0;
			for (int k = 0; k < 3; k++) {
				output->matrix3[i][j] += input.real_vector3[k] * (vector_summand[k] - scalar_component * quaternion_Array[j][k + 1] + quaternion_Array[j][0] * vector_component[k]);
			}
		}
	}
}