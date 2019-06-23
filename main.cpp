#include <iostream>
#include <math.h>
#include "lib/rkf45.h"
#include "lib/rkf45.cpp"

double Exact(double t){
  return 1.0/t;
}

// double eulerFun(double x, double y, double t){
//   double expT = exp(t);
//   return (2.0 * x + expT * y)/(expT + 1.0);
// }

void F(double t, double *y, double *dy){
  dy[0] = y[1];
  dy[1] = -(t*t*t*y[1] + (t*t - 2.0)*y[0])/t*t;
}

double PreciseY0(double t)
{
	return Exact(t);
}


double PreciseY1(double t) {
	return -1 * Exact(t);
}


void adams (void (*func) (double t, double *y, double *dy), int NEQN, double y_last[], double
yp_last[], double y_curr[], double yp_curr[], double t, double h) {

	func(t - h, y_last, yp_last);
	func(t, y_curr, yp_curr);
	for (int i = 0; i < NEQN; i++) {

		y_last[i] = y_curr[i];
		y_curr[i] = y_curr[i] + h / 2 * (3 * yp_curr[i] - yp_last[i]);
	}
}

void computeWithAdams(void (*func)(double t, double *y, double *dy), int NEQN, double h, double
t0, double t1) {

	double lastY[NEQN], lastYP[NEQN], currY[NEQN], currYP[NEQN], lastY1[NEQN], lastYP1[NEQN],
	currY1[NEQN], currYP1[NEQN];
	double lclErrSum = 0;
	double t = t0;
	double tt = t;
	lastY[0] = PreciseY1(t - h);
	lastY[1] = PreciseY0(t - h);
	currY[0] = PreciseY1(t);
	currY[1] = PreciseY0(t);

	printf(" %.2f % 15.8f % 15.8f % 15.8f % 15.8f\n",
	t, currY[1], currY[0], fabs(currY[1]-PreciseY0(t)), fabs(currY[0]-PreciseY1(t)));
	while (t < t1) {

		adams(func, NEQN, lastY, lastYP, currY, currYP, t, h);

		printf(" %.2f % 15.8f % 15.8f % 15.8f % 15.8f ",
		t+h, currY[1], currY[0], fabs(currY[1]-PreciseY0(t+h)), fabs(currY[0]-
		PreciseY1(t+h)));

		lastY1[0] = PreciseY1(t - h);
		lastY1[1] = PreciseY0(t - h);
		currY1[0] = PreciseY1(t);
		currY1[1] = PreciseY0(t);

		adams(func, NEQN, lastY1, lastYP1, currY1, currYP1, tt, h);

		printf(" y_lcl=%11.8f 	err_lcl=%11.8f\n", currY1[1], currY1[1] - PreciseY0(t + h));
		lclErrSum += currY1[1] - PreciseY0(t + h);

		t += h; tt = t;
	}
	printf("\n 		 lclErrSum=%11.8f\n", lclErrSum);
}


int main(){
  double t = 1.0, tmax = 2.0, h = 0.1, x = -1.0, y = 1.0;
  int arraySize = 1 + (tmax-t)/h;
  
  double *exact = new double[arraySize],
	*euler = new double[arraySize],
	*rkf45 = new double[arraySize],
	*eulerEps = new double[arraySize],
	*rkf45Eps = new double[arraySize],
	*res = new double[2],
	globalEulerEps = 0.0,
	globalRkf45Eps = 0.0,
	
	
	*Y = new double[2],
	*WORK = new double[15],
	T, TOUT,
	RELERR = 1.0e-8, ABSERR = 1.0e-8;
  
  int *IWORK = new int[5],
	IFLAG = 1,
	N = 2;

	exact[0] = y;
	// euler[0] = y;
	rkf45[0] = y;

	Y[0] = y;
	Y[1] = x;

	printf(" %u ", IFLAG);

	for (int index = 1; index < arraySize; index++){
		exact[index] = Exact(t+h*index);

		// Euler(t+h*index, h, x, euler[index - 1], res);
		// euler[index] = res[0];
		// x = res[1];

		T = t+h*(index - 1);
		TOUT = T+h;
		
		RKF45(F, N, Y, &T, &TOUT, &RELERR, &ABSERR, &IFLAG, WORK, IWORK);
		rkf45[index] = Y[0];
		printf(" %u ", IFLAG);

		// eulerEps[index] = abs(exact[index] - euler[index]);
		// globalEulerEps += eulerEps[index];
		// rkf45Eps[index] = abs(exact[index] - rkf45[index]);
		// globalRkf45Eps += rkf45Eps[index];
	}

	printf("\nT\tExact\t\tRKF45\n");

	for (int index = 0; index < arraySize; index++){
		printf("%.2f\t%.7f\t%.7f\n",t+h*index, exact[index], rkf45[index]);
	}
  
	// printf("\nT\tEXACT\t\tEULER\t\tEULER_EPS\tRKF45\t\tRKF45_EPS\n");

	// for (int index = 0; index < arraySize; index++){
	// 	printf("%.2f\t%.7f\t%.7f\t%.7f\t%.7f\t%.20f\n", t+h*index,
	// 	 exact[index],
	// 	 euler[index], eulerEps[index],
	// 	 rkf45[index], rkf45Eps[index]);
	// }

	// printf("\nGLOBAL EPS\': EULER = %.7f,\t RKF45 = %.15f", globalEulerEps, globalRkf45Eps);
	printf("\n 		Adams\n");

	printf(" t 		y 		y' 		eps1		eps2 \n");


	int Negn = 2;
	computeWithAdams(F, Negn, h, 1, 2);

  return 0;
}