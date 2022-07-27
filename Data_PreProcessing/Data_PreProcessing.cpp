#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <string>
#include "Noise_Generation.h"
#include "..\\Triad_Optimization\\Data_Types.h"
using namespace std;

int main()
{
	// Set up of pseudo-random generator and of the needed distributions.
	default_random_engine generator;
	normal_distribution<double> normal_Distribution(0.5, 0.166667);
	uniform_real_distribution<double> uniform_Distribuiton(0, 6.283185);
	// Initialisations of the required i/o streams.
	ifstream a1, a2, a3, a7, a8, a9;
	ofstream vector_data, text_data;
	a1.open("..\\Data\\a1.txt");
	a2.open("..\\Data\\a2.txt");
	a3.open("..\\Data\\a3.txt");
	a7.open("..\\Data\\a7.txt");
	a8.open("..\\Data\\a8.txt");
	a9.open("..\\Data\\a9.txt");
	vector_data.open("..\\Data\\Vector_Data.dat", ios::out | ios::binary | ios::trunc);
	text_data.open("..\\Data\\Vector_Data.txt", ios::out | ios::trunc);
	// The first 499 (~= 2 minutes of data) should be skipped, as they represent a transition period in the model.
	a1.seekg(32934);
	a2.seekg(32934);
	a3.seekg(32934);
	a7.seekg(32934);
	a8.seekg(32934);
	a9.seekg(32934);
	// The precisions of the sensors, ie. the maximum angle between the measured vector and its real value.
	double precisions[3] = { 0.0349065850399, 0.0349065850399, 0.0698131700798 };
	// Set up of the i/o stream for the readable text output.
	text_data << scientific;
	text_data.precision(7);
	while (!a1.eof()) {
		input_Data data;
		// Reading the input data for the real vector projections without noise.
		double time, real_Vectors[3][3];
		a1 >> time >> real_Vectors[0][0] >> real_Vectors[0][1] >> real_Vectors[0][2];
		a2 >> time >> real_Vectors[1][0] >> real_Vectors[1][1] >> real_Vectors[1][2];
		a3 >> time >> real_Vectors[2][0] >> real_Vectors[2][1] >> real_Vectors[2][2];
		// Check to skip the periods of time when the measured direction to the sun becomes a 0-vector.
		// Once the TRIAD-optimisation algorithm is expanded to operate on such data, this fragment of code should be removed.
		double a1_norm = scalar_Product(real_Vectors[0], real_Vectors[0]);
		if (a1_norm == 0) { continue; }
		// Reading the input data for the model vector projections.
		a7 >> time >> data.model_Vectors[0][0] >> data.model_Vectors[0][1] >> data.model_Vectors[0][2];
		a8 >> time >> data.model_Vectors[1][0] >> data.model_Vectors[1][1] >> data.model_Vectors[1][2];
		a9 >> time >> data.model_Vectors[2][0] >> data.model_Vectors[2][1] >> data.model_Vectors[2][2];
		// Generation of the random parameters and the application of the measurement noise to the real vectors.
		double random_theta_component, random_phi;
		for (int i = 0; i < 3; i++) {
			random_phi = uniform_Distribuiton(generator);
			random_theta_component = normal_Distribution(generator);
			while (random_theta_component < 0 || random_theta_component > 1)
			{
				random_theta_component = normal_Distribution(generator);
			}
			random_theta_component *= precisions[i];
			apply_Random_Noise(real_Vectors[i], random_theta_component, random_phi, data.real_Vectors[i]);
		}
		// Writing of the data to the .dat and .txt files.
		data.time = time;
		vector_data.write((char*)& data, sizeof(input_Data));
		text_data << time << "\n";
		for (int i = 0; i < 3; i++) {
			text_data << "a" << i+1 << ":   " << setw(15) << data.real_Vectors[i][0] << "  " << setw(15) << data.real_Vectors[i][1] << "  " << setw(15) << data.real_Vectors[i][2] << "\n";
		}
		for (int i = 0; i < 3; i++) {
			text_data << "a" << i + 7 << ":   " << setw(15) << data.model_Vectors[i][0] << "  " << setw(15) << data.model_Vectors[i][1] << "  " << setw(15) << data.model_Vectors[i][2] << "\n";
		}
		text_data << "\n";
	}
	a1.close();
	a2.close();
	a3.close();
	a7.close();
	a8.close();
	a9.close();
	vector_data.close();
	text_data.close();
}