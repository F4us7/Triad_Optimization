#include <iostream>
#include <fstream>
#include <iomanip>
#include "Parameter_Calculation.h"
#include "Newton_Method.h"
using namespace std;


// First complete and stable version. Currently set to test the algorithm on single data entries.
// To do: update the code operate over the cpmplete dataset.
int main()
{
	// Current test data:
	//input_Data test = { 131.43419, {{9.9667318e-01, -6.7044180e-17, 8.1502031e-02}, {-2.4004261e-02, -8.8525934e-01, 4.6447787e-01}, {-9.7717547e-01, 1.6826065e-01, -1.2967827e-01}},
	//{{-9.9983266e-01, 1.7191604e-02, 6.2537704e-03}, {3.2220124e-03, 9.9997612e-01, -6.1137945e-03}, {9.8159624e-01, -1.8871919e-01, 2.9221237e-02}}, {0.33333, 0.33333, 0.33333} };
	//input_Data test = { 131.43419, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, {{0.9845, 0.01254, 0.0058}, {0, -1, 0}, {0.1254, 0.3345, -1.0141}}, {0.2, 0.7, 0.1} }

	ifstream vector_data;
	ofstream text_results;
	vector_data.open("..\\Data\\Vector_Data.dat", ios::in | ios::binary);
	//vector_data.seekg(11896 * sizeof(input_Data));
	//vector_data.seekg(23 * sizeof(input_Data));
	text_results.open("..\\Data\\Results.txt", ios::in | ios::trunc);
	int counter = 0;
	while (!vector_data.eof()) {
		counter++;
		input_Data data;
		vector_data.read((char*)& data, sizeof(input_Data));
		TRIAD_Output TRIAD_results;
		numerical_Parameters parameters;
		calculate_Numerical_Parameters(data, TRIAD_results, parameters);
		optimisation_Output results = Newtons_Algorithm(parameters, TRIAD_results);
		text_results << scientific;
		text_results.precision(7);
		text_results << setw(15) << data.time << "  ->  Loss Value: " << setw(15) << results.loss_Value << ";  Delta: " << setw(15) << results.delta
			<< ";  Iterations: " << setw(6) << results.iteration_count << ";  Initial case: " << results.initial_values_case << ";\n";
		text_results << setw(40) << "    Approximative basis coordinates: ";
		for (int i = 0; i < 4; i++) {
			text_results << setw(15) << results.coordinates[i];
		}
		text_results << "\n" << setw(40) << "    Optimal quaternion: ";
		for (int i = 0; i < 4; i++) {
			text_results << setw(15) << results.optimal_quaternion[i];
		}
		text_results << "\n\n";
		cout << counter << "\n";
	}
	vector_data.close();
	text_results.close();
	return 0;
}