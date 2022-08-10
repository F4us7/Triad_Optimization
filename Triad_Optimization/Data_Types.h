#pragma once

// Structure to store the fundamental input data, including the projections of the 3 3D vectors and the sensor accuracy coefficients.
struct input_Data {
	double time;
	double model_Vectors[3][3];
	double real_Vectors[3][3];
};

// Structure used to pass data into the TRIAD algorithm.
struct input_Package {
	int counter;
	double *model_vector1, *model_vector2;
	double *real_vector1, *real_vector2;
};

// Structure used to store the output form the TRIAD algorithm: the 3 approximative quaternions.
struct TRIAD_Output {
	double quaternion_array[4][4];
};

// Structure used to store all the numerical data required for the efficient execution of Newton's method.
struct numerical_Parameters {
	double matrices[4][4][4];
	double multi_use_sum, alpha4;
	double alpha4_derivatives[3];
};

// Structure to contain the logged output from the TRIAD_optimisation algorithm.
struct optimisation_Output {
	double mu, coordinates[4];
	double loss_Value, delta, optimal_quaternion[4];
	int initial_values_case, iteration_count;
};