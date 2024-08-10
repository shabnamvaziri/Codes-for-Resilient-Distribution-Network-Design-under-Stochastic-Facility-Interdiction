#include "def.h"

extern int         S, S2, iter, J, K, I, Iperim, *keepx, max_iter, iter1, sGlobal, *keepx_Keep, *y, *tempx, *H, *tempH, *set, filled, fillcounter, *selectedscenario, tempi, *tempx_keep;
extern double     *x1, theta, theta0, B, *b, *c, *cperim, Q, theta1, *lambda, *temps, Qcurrent;
extern double     MAX_DOUBLE;
extern double     MAX_DOUBLE2;
extern double     *xi, *xi2, *p, *d, *alpha, *alpha1, *delta, *epsilon, *beta, *gamma1, *Delta1, obj, *delta_Keep, *epsilon_Keep, *beta_Keep, *Delta1_Keep, obj1;
extern double     *vPerim, *v, rho, *f;
extern int        *jPerim, jCounter;
extern int        *supplier, *openFacility;
extern double UB;
extern double objective;
int counter;

void createScenario(void)
{
	S = pow((double)2, (double)J);
	//xi = create_double_matrix(jPerim[iter], S[iter]);
	//tempxi = create_double_matrix(jPerim, pow(2, jPerim));
	counter = 0;
	rec(0, J);
	rec(1, J);
}

void rec2(int val, int count, int b) {
	if (count <= 1) {
		int i;
		for (i = b - 1; i >= 0; i--) {
//			printf("%d", (val >> i) & 1);
			xi[(b - i - 1) * S + counter] = (val >> i) & 1;
		}
		counter++;
	//	printf("\n");
	}
	else {
		rec2(val * 2, count - 1, b);
		rec2(val * 2 + 1, count - 1, b);
	}
}

void rec(int val, int count) {
	rec2(val, count, count);
}

	