#pragma once

struct input_Data {
	double time;
	double model_vector1[3], model_vector2[3], model_vector3[3];
	double real_vector1[3], real_vector2[3], real_vector3[3];
};

struct input_Package {
	int counter;
	double model_vector1[3], model_vector2[3];
	double real_vector1[3], real_vector2[3];
};

struct TRIAD_Output {
	double time;
	double quaternion1[4], quaternion2[4], quaternion3[4];
};