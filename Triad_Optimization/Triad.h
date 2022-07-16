#pragma once
#include "Data_Types.h"

double scalar_Product(double* vector1, double* vector2) {
	double running_sum = 0;
	for (int i = 0; i < 3; i++) { running_sum += vector1[i] * vector2[i]; };
	return running_sum;
};

void vector_Product(double* vector1, double* vector2, double* output) {
	output[0] = (vector1[1] * vector2[2]) - (vector1[2] * vector2[1]);
	output[1] = (vector1[2] * vector2[0]) - (vector1[0] * vector2[2]);
	output[2] = (vector1[0] * vector2[1]) - (vector1[1] * vector2[0]);
}

void TRIAD(input_Package input, TRIAD_Output &output) {
	//Auxiliary matrix from model data (M), auxiliary matrix from real data (R).
	double *matrix_Model[3], *matrix_Real[3]; 
	// Calculation of the first vector pair of the 2 auxiliary Triads. 
	double model_Unit1[3], real_Unit1[3];
	double model_Norm1 = sqrt(scalar_Product(input.model_vector1, input.model_vector1));
	double real_Norm1 = sqrt(scalar_Product(input.real_vector1, input.real_vector1));
	for (int i = 0; i < 3; i++) {
		model_Unit1[i] = input.model_vector1[i] / model_Norm1;
		real_Unit1[i] = input.real_vector1[i] / real_Norm1;
	}
	matrix_Model[0] = model_Unit1;
	matrix_Real[0] = real_Unit1;
	// Calculation of the third vector pair of the 2 auxiliary Triads.
	double model_Unit3[3], real_Unit3[3];
	vector_Product(input.model_vector1, input.model_vector2, model_Unit3);
	vector_Product(input.real_vector1, input.real_vector2, real_Unit3);
	double model_Norm3 = sqrt(scalar_Product(model_Unit3, model_Unit3));
	double real_Norm3 = sqrt(scalar_Product(real_Unit3, real_Unit3));
	for (int i = 0; i < 3; i++) {
		model_Unit3[i] /= model_Norm3;
		real_Unit3[i] /= real_Norm3;
	}
	matrix_Model[2] = model_Unit3;
	matrix_Real[2] = real_Unit3;
	// Calculation of the second vector pair of the 2 auxiliary Triads.
	double model_Unit2[3], real_Unit2[3];
	vector_Product(model_Unit3, model_Unit1, model_Unit2);
	vector_Product(real_Unit3, real_Unit1, real_Unit2);
	matrix_Model[1] = model_Unit2;
	matrix_Real[1] = real_Unit2;
	// Calculation of the transition matrix C=((R^T)*(M^T))^T.
	double matrix_Transition[3][3] = { {}, {}, {} };
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				matrix_Transition[j][i] += matrix_Real[i][k] * matrix_Model[k][j];
			}
		}
	}
	//Choosing the output variable and calculating the Quaternion.
	double* output_Quaternion;
	switch (input.counter) {
	case 1:
		output_Quaternion = output.quaternion1;
		break;
	case 2:
		output_Quaternion = output.quaternion2;
		break;
	case 3:
		output_Quaternion = output.quaternion3;
	}
	output_Quaternion[0] = 0.5 * sqrt(1 + matrix_Transition[0][0] + matrix_Transition[1][1] + matrix_Transition[2][2]);
	output_Quaternion[1] = (matrix_Transition[1][2] - matrix_Transition[2][1]) / (4 * output_Quaternion[0]);
	output_Quaternion[2] = (matrix_Transition[2][0] - matrix_Transition[0][2]) / (4 * output_Quaternion[0]);
	output_Quaternion[3] = (matrix_Transition[0][1] - matrix_Transition[1][0]) / (4 * output_Quaternion[0]);
};