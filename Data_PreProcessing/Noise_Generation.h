#pragma once
#include <cmath>
#include <iostream>
using namespace std;

// Function to calculate the scalar product of 2 3D vectors. 
double scalar_Product(double* vector1, double* vector2) {
	double running_sum = 0;
	for (int i = 0; i < 3; i++) { running_sum += vector1[i] * vector2[i]; };
	return running_sum;
};

// Function to calculate the vector product of 2 3D vectors.
void vector_Product(double* vector1, double* vector2, double* output) {
	output[0] = (vector1[1] * vector2[2]) - (vector1[2] * vector2[1]);
	output[1] = (vector1[2] * vector2[0]) - (vector1[0] * vector2[2]);
	output[2] = (vector1[0] * vector2[1]) - (vector1[1] * vector2[0]);
}

// Function to apply a random measurement error to a vector, which is passed in as a pointer input_vector.
// The random angles used to apply the measurement error are random_theta (size of measurement deviation from the real value) and random_phi (direction of the diviation).
// The resulting vector with noise is saved in the output array.
void apply_Random_Noise(double* input_vector, double random_theta, double random_phi, double* output) {
	// Calculation of the parameters of the transform which maps the vector (0, 0, 1) onto the input vector.
	double z_axis_direction[3] = { 0, 0, 1 }, x_axis_direction[3] = { 1, 0 ,0 }, K_axis_direction[3];
	vector_Product(z_axis_direction, input_vector, K_axis_direction);
	double K_axis_norm = sqrt(scalar_Product(K_axis_direction, K_axis_direction));
	for (int i = 0; i < 3; i++) {
		K_axis_direction[i] /= K_axis_norm;
	}
	double transition_psi = acos(scalar_Product(K_axis_direction, x_axis_direction));
	if (K_axis_direction[1] < 0) { transition_psi *= -1; }
	double transition_theta = acos(scalar_Product(input_vector, z_axis_direction));
	// Calculation of the matrix defining the aforementioned transform.
	double transition_matrix[3][3] = { {cos(transition_psi), -cos(transition_theta) * sin(transition_psi), sin(transition_theta) * sin(transition_psi)},
	{sin(transition_psi), cos(transition_theta) * cos(transition_psi), -sin(transition_theta) * cos(transition_psi)},
	{0, +sin(transition_theta), cos(transition_theta)} };
	// Calculation of the vector with measurement noise.
	double random_vector[3] = { sin(random_theta) * cos(random_phi), sin(random_theta) * sin(random_phi), cos(random_theta) };
	for (int i = 0; i < 3; i++) {
		output[i] = 0;
		for (int k = 0; k < 3; k++) {
			output[i] += transition_matrix[i][k] * random_vector[k];
		}
	}
}