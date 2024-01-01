
/**
 ******************************************************************************
 * file           : main.c
 * author         : Aly_Mustafa_Enaya
 * brief          : Get the angels using Newton-Raphson Method
 * Date			  : 1/1/2024
 ******************************************************************************
 */


#include <stdio.h>
#include <math.h>

//Put The Desired voltage you want
double Vdesired = 200.0;

//Put The Peak voltage you want
double Vpeak = 400.0;

//Multiply the two matrices
void multiplyMatrices(double x[6][6], double y[6][1], double result[6][1]) {
	// Perform matrix multiplication
	for (int i = 0; i < 6; ++i) {
		result[i][0] = 0;
		for (int j = 0; j < 6; ++j) {
			result[i][0] += x[i][j] * y[j][0];
		}
	}
}

//Subtract the two matrices
void subtractMatrices(double matrix1[][1], double matrix2[][1], double result[][1], int rows) {
	for (int i = 0; i < rows; ++i) {
		result[i][0] = matrix1[i][0] - matrix2[i][0];
	}
}
// Functions definition: f(alpha_1, alpha_2)
double function_1(double alpha_1,double alpha_2,double alpha_3,double alpha_4,double alpha_5,double alpha_6 ) {
	return ( (2*Vpeak/M_PI)*(1 - 2*cos(alpha_1) + 2*cos(alpha_2) - 2*cos(alpha_3) + 2*cos(alpha_4) - 2*cos(alpha_5) + 2*cos(alpha_6))   -   Vdesired  );}
double function_3(double alpha_1,double alpha_2,double alpha_3,double alpha_4,double alpha_5,double alpha_6 ) {
	return ( (2*Vpeak/M_PI)*(1 - 2*cos(3*alpha_1) + 2*cos(3*alpha_2) - 2*cos(3*alpha_3) + 2*cos(3*alpha_4) - 2*cos(3*alpha_5) + 2*cos(3*alpha_6)));}
double function_5(double alpha_1,double alpha_2,double alpha_3,double alpha_4,double alpha_5,double alpha_6 ) {
	return ( (2*Vpeak/M_PI)*(1 - 2*cos(5*alpha_1) + 2*cos(5*alpha_2) - 2*cos(5*alpha_3) + 2*cos(5*alpha_4) - 2*cos(5*alpha_5) + 2*cos(5*alpha_6)));}
double function_7(double alpha_1,double alpha_2,double alpha_3,double alpha_4,double alpha_5,double alpha_6 ) {
	return ( (2*Vpeak/M_PI)*(1 - 2*cos(7*alpha_1) + 2*cos(7*alpha_2) - 2*cos(7*alpha_3) + 2*cos(7*alpha_4) - 2*cos(7*alpha_5) + 2*cos(7*alpha_6)));}
double function_9(double alpha_1,double alpha_2,double alpha_3,double alpha_4,double alpha_5,double alpha_6 ) {
	return ( (2*Vpeak/M_PI)*(1 - 2*cos(9*alpha_1) + 2*cos(9*alpha_2) - 2*cos(9*alpha_3) + 2*cos(9*alpha_4) - 2*cos(9*alpha_5) + 2*cos(9*alpha_6)));}
double function_11(double alpha_1,double alpha_2,double alpha_3,double alpha_4,double alpha_5,double alpha_6 ) {
	return ( (2*Vpeak/M_PI)*(1 - 2*cos(11*alpha_1) + 2*cos(11*alpha_2) - 2*cos(11*alpha_3) + 2*cos(11*alpha_4) - 2*cos(11*alpha_5) + 2*cos(11*alpha_6)));}

//*****************************************
// Function to calculate the inverse of a 6x6 matrix
//*****************************************
#define N 6

void getCofactor(double mat[N][N], double temp[N][N], int p, int q, int n) {
	int i = 0, j = 0;

	for (int row = 0; row < n; row++) {
		for (int col = 0; col < n; col++) {
			if (row != p && col != q) {
				temp[i][j++] = mat[row][col];

				if (j == n - 1) {
					j = 0;
					i++;
				}
			}
		}
	}
}

double determinant(double mat[N][N], int n) {
	double D = 0.0; // Initialize result

	if (n == 1)
		return mat[0][0];

	double temp[N][N]; // To store cofactors

	int sign = 1; // To store sign multiplier

	for (int f = 0; f < n; f++) {
		getCofactor(mat, temp, 0, f, n);
		D += sign * mat[0][f] * determinant(temp, n - 1);
		// terms are to be added with alternate sign
		sign = -sign;
	}

	return D;
}

void adjoint(double A[N][N], double adj[N][N]) {
	if (N == 1) {
		adj[0][0] = 1.0;
		return;
	}
	int sign = 1;
	double temp[N][N];

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			// Get cofactor of A[i][j]
			getCofactor(A, temp, i, j, N);

			// sign of adj[j][i] positive if sum of row
			// and column indexes is even.
			sign = ((i + j) % 2 == 0) ? 1 : -1;

			// Interchanging rows and columns to get the
			// transpose of the cofactor matrix
			adj[j][i] = sign * determinant(temp, N - 1);
		}
	}
}

// Function to calculate and store inverse, returns false if
// matrix is singular
int inverse(double A[N][N], double inv[N][N]) {
	// Find determinant of A[][]
	double det = determinant(A, N);
	if (det == 0.0) {
		printf("Singular matrix, can't find its inverse");
		return 0;
	}

	// Find adjoint
	double adj[N][N];
	adjoint(A, adj);

	// Find Inverse using formula "inverse(A) = adj(A)/det(A)"
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			inv[i][j] = adj[i][j] / det;

	return 1;
}

//##########################################################################

int main() {
	double InvMatrix[6][6];
	double result1[6][1];
	double final_matrix[6][1];
	//arrays
	double Init_alpha_arr[6][1];
	double Init_Fn_arr[6][1];

	double alpha_1 = 20 *M_PI/180.0;
	double alpha_2 = 30 *M_PI/180.0;
	double alpha_3 = 40 *M_PI/180.0;
	double alpha_4 = 50 *M_PI/180.0;
	double alpha_5 = 60 *M_PI/180.0;
	double alpha_6 = 70 *M_PI/180.0;

	double h = 0.0001;  // Small value for numerical approximation
	Init_alpha_arr[0][0] = alpha_1;
	Init_alpha_arr[1][0] = alpha_2;
	Init_alpha_arr[2][0] = alpha_3;
	Init_alpha_arr[3][0] = alpha_4;
	Init_alpha_arr[4][0] = alpha_5;
	Init_alpha_arr[5][0] = alpha_6;

	Init_Fn_arr[0][0] = function_1(alpha_1 , alpha_2 , alpha_3 , alpha_4 , alpha_5 , alpha_6);
	Init_Fn_arr[1][0] = function_3(alpha_1 , alpha_2 , alpha_3 , alpha_4 , alpha_5 , alpha_6);
	Init_Fn_arr[2][0] = function_5(alpha_1 , alpha_2 , alpha_3 , alpha_4 , alpha_5 , alpha_6);
	Init_Fn_arr[3][0] = function_7(alpha_1 , alpha_2 , alpha_3 , alpha_4 , alpha_5 , alpha_6);
	Init_Fn_arr[4][0] = function_9(alpha_1 , alpha_2 , alpha_3 , alpha_4 , alpha_5 , alpha_6);
	Init_Fn_arr[5][0] = function_11(alpha_1 , alpha_2 , alpha_3 , alpha_4 , alpha_5 , alpha_6);
	// Calculate partial derivatives
	//Jacobian_array
	double Jacobian_arr[6][6] =
	{
			//F1
			{(function_1(alpha_1+h,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6) 	- function_1(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_1(alpha_1 ,alpha_2+h,alpha_3,alpha_4,alpha_5,alpha_6)  	- function_1(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_1(alpha_1 ,alpha_2,alpha_3+h,alpha_4,alpha_5,alpha_6) 	- function_1(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_1(alpha_1 ,alpha_2,alpha_3,alpha_4+h,alpha_5,alpha_6) 	- function_1(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_1(alpha_1 ,alpha_2,alpha_3,alpha_4,alpha_5+h,alpha_6) 	- function_1(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_1(alpha_1 ,alpha_2,alpha_3,alpha_4+h,alpha_5,alpha_6+h) 	- function_1(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
			},

			//F3
			{(function_3(alpha_1+h,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6) 	- function_3(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_3(alpha_1 ,alpha_2+h,alpha_3,alpha_4,alpha_5,alpha_6)  	- function_3(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_3(alpha_1 ,alpha_2,alpha_3+h,alpha_4,alpha_5,alpha_6) 	- function_3(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_3(alpha_1 ,alpha_2,alpha_3,alpha_4+h,alpha_5,alpha_6) 	- function_3(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_3(alpha_1 ,alpha_2,alpha_3,alpha_4,alpha_5+h,alpha_6) 	- function_3(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_3(alpha_1 ,alpha_2,alpha_3,alpha_4+h,alpha_5,alpha_6+h) 	- function_3(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
			},

			//F5
			{(function_5(alpha_1+h,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6) 	- function_5(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_5(alpha_1 ,alpha_2+h,alpha_3,alpha_4,alpha_5,alpha_6)  	- function_5(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_5(alpha_1 ,alpha_2,alpha_3+h,alpha_4,alpha_5,alpha_6) 	- function_5(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_5(alpha_1 ,alpha_2,alpha_3,alpha_4+h,alpha_5,alpha_6) 	- function_5(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_5(alpha_1 ,alpha_2,alpha_3,alpha_4,alpha_5+h,alpha_6) 	- function_5(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_5(alpha_1 ,alpha_2,alpha_3,alpha_4+h,alpha_5,alpha_6+h) 	- function_5(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
			},

			//F7
			{(function_7(alpha_1+h,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6) 	- function_7(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_7(alpha_1 ,alpha_2+h,alpha_3,alpha_4,alpha_5,alpha_6)  	- function_7(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_7(alpha_1 ,alpha_2,alpha_3+h,alpha_4,alpha_5,alpha_6) 	- function_7(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_7(alpha_1 ,alpha_2,alpha_3,alpha_4+h,alpha_5,alpha_6) 	- function_7(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_7(alpha_1 ,alpha_2,alpha_3,alpha_4,alpha_5+h,alpha_6) 	- function_7(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_7(alpha_1 ,alpha_2,alpha_3,alpha_4+h,alpha_5,alpha_6+h) 	- function_7(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
			},

			//F9
			{(function_9(alpha_1+h,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6) 	- function_9(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_9(alpha_1 ,alpha_2+h,alpha_3,alpha_4,alpha_5,alpha_6)  	- function_9(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_9(alpha_1 ,alpha_2,alpha_3+h,alpha_4,alpha_5,alpha_6) 	- function_9(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_9(alpha_1 ,alpha_2,alpha_3,alpha_4+h,alpha_5,alpha_6) 	- function_9(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_9(alpha_1 ,alpha_2,alpha_3,alpha_4,alpha_5+h,alpha_6) 	- function_9(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_9(alpha_1 ,alpha_2,alpha_3,alpha_4+h,alpha_5,alpha_6+h) 	- function_9(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
			},
			//F11
			{(function_11(alpha_1+h,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6) 	- function_11(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_11(alpha_1 ,alpha_2+h,alpha_3,alpha_4,alpha_5,alpha_6)  	- function_11(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_11(alpha_1 ,alpha_2,alpha_3+h,alpha_4,alpha_5,alpha_6) 	- function_11(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_11(alpha_1 ,alpha_2,alpha_3,alpha_4+h,alpha_5,alpha_6) 	- function_11(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_11(alpha_1 ,alpha_2,alpha_3,alpha_4,alpha_5+h,alpha_6) 	- function_11(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
					,(function_11(alpha_1 ,alpha_2,alpha_3,alpha_4+h,alpha_5,alpha_6+h) 	- function_11(alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6)) /h
			}
	};

	inverse(Jacobian_arr, InvMatrix);
	multiplyMatrices(InvMatrix, Init_Fn_arr,result1);
	subtractMatrices(Init_alpha_arr, result1, final_matrix, 6);
	int i=0;
	setvbuf(stdout,NULL,_IONBF,0); setvbuf(stderr,NULL,_IONBF,0);
	printf("\n<========= For Desired Voltage = %0.1lfv =========>\n", Vdesired );

	for(i=0 ; i<6 ;i++)
	{
		setvbuf(stdout,NULL,_IONBF,0); setvbuf(stderr,NULL,_IONBF,0);
		printf("\t\tangle %d = %lf\n", i+1,final_matrix[i][0]*180/M_PI );
	}

	return 0;
}
