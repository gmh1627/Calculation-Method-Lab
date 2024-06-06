#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>


#define FILENAME "M_8.txt"

int M = 20;
double eps = 1e-6;
int Satisfied_Count;

double ax(double t) {
	return sin(t) / (sqrt(t) + 1);
}

double ay(double t) {
	return log(t + 1) / (t + 1);
}

double Romberg(double a, double b, double eps, int M, double f(double t), int is_dis);

double vx(double t) {
	return Romberg(0.0, t, eps, M, ax, 0);
}
double vy(double t) {
	return Romberg(0.0, t, eps, M, ay, 0);
}


int main() {
	double Succ_rate;

	printf("%d\n", M);
	Satisfied_Count = 0;

	FILE *fp = fopen(FILENAME, "w");
	if (fp == NULL) {
		printf("Error opening file.\n");
		return 1;
	}
	for (double t = 0.1; t <= 10; t += 0.1) {
		double v_x = vx(t);
		double v_y = vy(t);

		double x = Romberg(0.0, t, eps, M, vx, 1);
		double y = Romberg(0.0, t, eps, M, vy, 1);
		printf("t=%f,vx=%f, vy=%f, (%f, %f)\n", t, v_x, v_y, x, y);

		//fprintf(fp, "%.6f %.6f\n", x, y);

	}

	fclose(fp);
	printf("Trajectory data saved to %s\n", FILENAME);

	Succ_rate = Satisfied_Count / 200.0;
	printf("%f%\n", Succ_rate * 100.0);


}

double Romberg(double a, double b, double eps, int M, double f(double t), int is_dis) {
	double R[M][M];
	double h[M];

	h[0] = b - a;
	R[0][0] = (f(a) + f(b)) * h[0] / 2;
	for (int k = 1; k < M; k++) {
		h[k] = h[k - 1] / 2;
		double sum = 0.0;
		for (int i = 1; i <= pow(2, k - 1); i++) {
			sum += f(a + (2 * i - 1) * h[k]);
		}
		R[k][0] = (R[k - 1][0] + h[k - 1] * sum) / 2;


		for (int j = 1; j <= k; j++) {
			R[k][j] = R[k][j - 1] + (R[k][j - 1] - R[k - 1][j - 1]) / (pow(4, j) - 1);
		}
		if (k > 1 && fabs(R[k][k] - R[k - 1][k - 1]) < eps) {
			if (is_dis) {
				Satisfied_Count++;
				return R[k][k];
			}
		}
	}
	return R[M - 1][M - 1];
}
