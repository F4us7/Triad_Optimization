#pragma once

// Structure to store the fundamental input data, including the projections of the 3 3D vectors and the sensor accuracy coefficients.
struct input_Data {
	double time;
	double model_Vectors[3][3];
	double real_Vectors[3][3];
	double accuracies[3];
};

// Structure used to pass data into the TRIAD algorithm.
struct input_Package {
	int counter;
	double *model_vector1, *model_vector2;
	double *real_vector1, *real_vector2;
};

// Structure used to store the output form the TRIAD algorithm: the 3 approximative quaternions.
struct TRIAD_Output {
	double quaternion1[4], quaternion2[4], quaternion3[4];
};

// Structure used to store all the numerical data required for the efficient execution of Newton's method.
struct numerical_Parameters {
	double matrices[4][4][4];
	double multi_use_sum1, multi_use_sum2, alpha4;
	double alpha4_derivatives[3];
};