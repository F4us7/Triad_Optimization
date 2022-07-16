#pragma once

struct input_Data {
	double time;
	double model_Vectors[3][3];
	double real_Vectors[3][3];
};

struct input_Package {
	int counter;
	double *model_vector1, *model_vector2;
	double *real_vector1, *real_vector2;
};

// The 3 quaternion outputs of TRIAD are not essential for the algorithm and do not require external storage, but can be used to evaluate the initial loss values.
struct TRIAD_Output {
	double quaternion1[4], quaternion2[4], quaternion3[4];
};

struct numerical_Parameters {
	double matrices[4][4][4];
	double multi_use_sum1, multi_use_sum2, alpha4;
	double alpha4_derivatives[3];
};