#include <iostream>
#include <fstream>
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
	//input_Data test = { 131.43419, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, {{0.9845, 0.01254, 0.0058}, {0, -1, 0}, {0.1254, 0.3345, -1.0141}}, {0.2, 0.7, 0.1} };
	
	input_Data test = { 1, { {0, 1, 0}, {-1, 0, 0}, {0, 0, 1} }, { {0.999035926328132, -0.043474931522077, 0.006094935179475}, {-0.047106940640310, 0.998877184237770, 0.005030601627420}, {0.006798923593416, 0.013350913621106, 0.999887757572545} }, { 0.230769230769231, 0.461538461538461, 0.307692307692308 } };
	TRIAD_Output test_res = { {}, {}, {} };
	numerical_Parameters param;
	calculate_Numerical_Parameters(test, test_res, param);
	double res_coordinates[4];
	Newtons_Algorithm(param, res_coordinates);
	cout << Wahba_Loss_Function(param, res_coordinates);
}