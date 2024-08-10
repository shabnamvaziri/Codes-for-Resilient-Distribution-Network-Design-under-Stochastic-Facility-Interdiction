#include "def.h"////

extern int         S, S2, iter, J, K, I, Iperim, *keepx, max_iter, iter1, sGlobal, *keepx_Keep, *y, *tempx, H, *tempH, *set, filled, fillcounter, *selectedscenario, tempi, *tempx_keep;
extern double     *x1, theta, theta0, B, *b, *c, *cperim, Q, theta1, *lambda, *temps, Qcurrent, * correlation, * tau, * newbeta, **newbeta_Keep;
extern double     MAX_DOUBLE;
extern double     MAX_DOUBLE2;
extern double     *xi, *xi2, *p, *d, *alpha, *alpha1, *delta, *epsilon, *beta, *gamma1, *Delta1, obj, *delta_Keep, *epsilon_Keep, *beta_Keep, *Delta1_Keep, obj1, *prob;
extern double     *vPerim, *v, rho, *f;
extern int        *jPerim, jCounter;
extern int        *supplier, *openFacility, *Hconstraint, sumy;
extern double UB;
extern double objective;
extern double		MasterTime, SubTime;
extern clock_t		startTemp, endTemp;

double solve_SubProblem(void)
{
	int    max_iter1 = 100;
	int    interv_dec_lambda = 10;// number of iterations for epsilon/2, shabnam
	int    last = 1;
	double epsilonPrim = 2.25;
	double StepLength, LB;
	int i, j, k, s, ii, jj;
	double  best_upper_bound = MAX_DOUBLE;
	theta1 = MAX_DOUBLE;
	double best_lower_bound = 0;
	//int *Hconstraint;
	int inside;
	//Hconstraint = create_int_vector((J + 1)* S);

	for (s = 0; s < S; s++){
		for (k = 0; k < K; k++){
			beta[(k*S + s)] = 0;
		}
		for (i = 0; i < I; i++){
			delta[(i*S + s)] = 0;
		}
		for (j = 0; j < J + 1; j++){
			for (jj = 0; jj < J + 1; jj++) {
				tau[(j * (J + 1) + jj) + (s * (J + 1) * (J + 1))] = 0;
			}
		}
		for (j = 0; j < J + 1; j++){
			for (jj = 0; jj < J + 1; jj++) {
				newbeta[(j * (J + 1) + jj) + (s * (J + 1) * (J + 1))] = 0;
			}
		}

	}
	for (j = 0; j < J + 1; j++){
		lambda[j] = 0;
		temps[j] = 0;
	}

	iter1 = 0;
	y[J] = 1;
	alpha1[J] = 0;
	for (s = 0; s < S; s++){
		xi[J * S + s] = 0;
	}
	b[J] = 100;
	for (k = 0; k < K; k++){
		c[J*K + k] = 150000;
		alpha1[J] += d[k];
	}
	for (i = 0; i < Iperim; i++){
	cperim[i*(J + 1) + J] = 150000;
}
for (j = 0; j < J + 1; j++) {
	correlation[J * (J + 1) + j] = 0;
}
	if (sumy == J){


		gettimeofday(&start, NULL);
		//startTemp = clock();
		for (iter1 = 0; iter1 < max_iter1; iter1++){
			Q = 0;
			for (sGlobal = 0; sGlobal < S; sGlobal++){  ////adding on march 17
				solveCPLEX1();
			}


			double *xx;
			double *xxx;
			double *xxxx;
			xx = create_double_vector(S);
			xxx = create_double_vector(S);
			xxxx = create_double_vector(S);
			for (sGlobal = 0; sGlobal < S; sGlobal++){
				xx[tempi] = 0;
				xxx[tempi] = 0;
				xxxx[tempi] = 0;
			}
			for (sGlobal = 0; sGlobal < S; sGlobal++){
				for (k = 0; k < K; k++){
					xx[tempi] += p[sGlobal] * d[k] * beta[(k*S + sGlobal)];
				}
				for (i = 0; i < I; i++){
					xxx[tempi] += p[sGlobal]* alpha[i] * delta[(i*S + sGlobal)];
				}
				for (j = 0; j < J + 1; j++){
					for (jj = 0; jj < J + 1; jj++) {
						xxxx[tempi] += p[sGlobal] * (alpha1[j] * y[j] * newbeta[(j * (J + 1) + jj) + (sGlobal * (J + 1) * (J + 1))]);
					}
				}
				Q += (xx[tempi] - xxx[tempi] - xxxx[tempi]);
			}


			/////updating lower and upper bound
			if (Q < best_upper_bound){
				best_upper_bound = Q;
				LB = LowerBound1();
				if (LB > best_lower_bound){                          //Update best lower bound
					best_lower_bound = LB;
				}
			}
			///////updating lagrangean multiplier/////
			double SqrNorm = 0;
			for (j = 0; j < J + 1; j++){
				temps[j] = 0;
				for (sGlobal = 0; sGlobal < S; sGlobal++){
					temps[j] += p[sGlobal] * (Hconstraint[sGlobal] * tempx[j*S + sGlobal]);
				}
				SqrNorm += pow((double)temps[j], 2);
			}

			//Check stoping criteria
			if ((best_upper_bound - best_lower_bound) / best_upper_bound * 100 < 0.5)
				break;

			//Compute steplength
			if (iter1 - last > interv_dec_lambda) {
				last = iter1;
				epsilonPrim = epsilonPrim / 2;
			}
			if (epsilonPrim < 0.1)
				epsilonPrim = 2;

			StepLength = epsilonPrim* (best_upper_bound - best_lower_bound) / SqrNorm;

			//Update Lagrange multipliers
			for (j = 0; j < J + 1; j++){
				lambda[j] = lambda[j] + StepLength*temps[j];
			}
		}

		int counter = 0;
		for (j = 0; j < J; j++){
			for (sGlobal = 0; sGlobal < S; sGlobal++){
				if (tempx[j*S + sGlobal] == tempx[j * S + 0]){
					counter++;
				}
				else
					break;
			}
		}

		filled = 0;
		int count2;
		int identical;

		if (counter != J*S){
			for (sGlobal = 0; sGlobal < S; sGlobal++){
				identical = 0;
				for (fillcounter = 0; fillcounter < filled; fillcounter++){
					count2 = 0;
					for (j = 0; j < J; j++){
						if (tempx[j*S + sGlobal] == set[(fillcounter*(J + 1) + j)])
							count2++;
						else
							break;
					}
					if (count2 == J)
						identical = 1;
					if (identical == 1)
						break;
				}
				if (identical == 0){
					filled++;
					for (j = 0; j < J + 1; j++)
						set[(fillcounter*(J + 1) + j)] = tempx[j*S + sGlobal];
				}
			}
			double currentLB = 0;
			for (fillcounter = 0; fillcounter < filled; fillcounter++){

				Qcurrent = 0;
				for (sGlobal = 0; sGlobal < S; sGlobal++){
					XCPLEX1();
				}


				double *xx;
				double *xxx;
				double *xxxx;
				xx = create_double_vector(S);
				xxx = create_double_vector(S);
				xxxx = create_double_vector(S);
				for (sGlobal = 0; sGlobal < S; sGlobal++){
					xx[sGlobal] = 0;
					xxx[sGlobal] = 0;
					xxxx[sGlobal] = 0;
				}
				for (sGlobal = 0; sGlobal < S; sGlobal++){
					for (k = 0; k < K; k++){
						xx[sGlobal] += p[sGlobal] * d[k] * beta[(k*S + sGlobal)];
					}
					for (i = 0; i < I; i++){
						xxx[sGlobal] += p[sGlobal] * alpha[i] * delta[(i*S + sGlobal)];
					}
					for (j = 0; j < J + 1; j++){
						for (jj = 0; jj < J + 1; jj++) {
							xxxx[sGlobal] += p[sGlobal] * (alpha1[j] * y[j] * newbeta[(j * (J + 1) + jj) + (sGlobal * (J + 1) * (J + 1))]);
						}
					}
					Qcurrent += (xx[sGlobal] - xxx[sGlobal] - xxxx[sGlobal]);
				}
				if (Qcurrent > currentLB){
					currentLB = Qcurrent;
					for (sGlobal = 0; sGlobal < S; sGlobal++){
						for (k = 0; k < K; k++){
							beta_Keep[(k*S + sGlobal) + ((iter)*S*K)] = beta[(k*S + sGlobal)];
						}
						for (i = 0; i < I; i++){
							delta_Keep[(i*S + sGlobal) + ((iter)*I*S)] = delta[(i*S + sGlobal)];
						}
						for (j = 0; j < J; j++){
							for (jj = 0; jj < J + 1; jj++) {

								newbeta_Keep[(j * (J + 1) + jj) + (sGlobal * (J + 1) * (J + 1))][iter] = newbeta[(j * (J + 1) + jj) + (sGlobal * (J + 1) * (J + 1))];
							}
							tempx_keep[j*S + sGlobal] = set[(fillcounter*(J + 1) + j)];
						}
					}
				}
			}
			gettimeofday(&stop, NULL);
			//endTemp = clock();
			//SubTime = (double)(endTemp - startTemp) / (double)(CLOCKS_PER_SEC);
			SubTime = ((double)(stop.tv_sec - start.tv_sec) * 1000 + (double)(stop.tv_usec - start.tv_usec) / 1000) / 1000;
			return currentLB;
		}
		else{
			for (sGlobal = 0; sGlobal < S; sGlobal++){
				for (k = 0; k < K; k++){
					beta_Keep[(k*S + sGlobal) + ((iter)*S*K)] = beta[(k*S + sGlobal)];
				}
				for (i = 0; i < I; i++){
					delta_Keep[(i*S + sGlobal) + ((iter)*I*S)] = delta[(i*S + sGlobal)];
				}
				for (j = 0; j < J; j++){
					for (jj = 0; jj < J + 1; jj++) {
						newbeta_Keep[(j * (J + 1) + jj) + (sGlobal * (J + 1) * (J + 1))][iter] = newbeta[(j * (J + 1) + jj) + (sGlobal * (J + 1) * (J + 1))];
					}
				}
			}
			
			gettimeofday(&stop, NULL);
			//endTemp = clock();
			//SubTime = (double)(endTemp - startTemp) / (double)(CLOCKS_PER_SEC);
			SubTime = ((double)(stop.tv_sec - start.tv_sec) * 1000 + (double)(stop.tv_usec - start.tv_usec) / 1000) / 1000;
			return Q;
		}

	}
	else{
printf("sub4\n");
		prob = create_double_vector(pow(2, jPerim[iter])* max_iter);
		xi2 = create_double_vector((jPerim[iter] + 1)* pow(2, jPerim[iter]));//adding on March 17
		createScenario2();
		for (tempi = 0; tempi < S2; tempi++){
			prob[tempi*max_iter + iter] = 1;
		}
		for (j = 0; j < jPerim[iter]; j++){
			for (tempi = 0; tempi < S2; tempi++){
				if (xi2[j*S2 + tempi] == 1)
					prob[tempi*max_iter + iter] *= 0.75;
				else
					prob[tempi*max_iter + iter] *= 0.25;
			}
		}
		selectedscenario = create_int_vector(S2);////adding on march 17
		for (tempi = 0; tempi < S2; tempi++){  ////adding on march 17
			for (sGlobal = 0; sGlobal < S; sGlobal++){
				////adding on march 17
				int count = 0;
				for (j = 0; j < jPerim[iter]; j++){
					if (xi2[j*S2 + tempi] == xi[openFacility[j*max_iter + iter] * S + sGlobal]){
						count++;
					}
					else
						break;
				}
				if (count == jPerim[iter]){
					selectedscenario[tempi] = sGlobal;
					break;
				}
			}
		}
		gettimeofday(&start, NULL);
		//startTemp = clock();
		for (iter1 = 0; iter1 < max_iter1; iter1++){
			Q = 0;
			for (tempi = 0; tempi < S2; tempi++){  ////adding on march 17
				solveCPLEX();
			}

			double *xx;
			double *xxx;
			double *xxxx;
			xx = create_double_vector(S2);
			xxx = create_double_vector(S2);
			xxxx = create_double_vector(S2);
			for (tempi = 0; tempi < S2; tempi++){
				xx[tempi] = 0;
				xxx[tempi] = 0;
				xxxx[tempi] = 0;
			}
			for (tempi = 0; tempi < S2; tempi++){
				for (k = 0; k < K; k++){
					xx[tempi] += prob[tempi*max_iter + iter] * d[k] * beta[(k*S + selectedscenario[tempi])];
				}
				for (i = 0; i < I; i++){
					xxx[tempi] += prob[tempi*max_iter + iter] * alpha[i] * delta[(i*S + selectedscenario[tempi])];
				}
				for (j = 0; j < J + 1; j++){
					for (jj = 0; jj < J + 1; jj++) {
						xxxx[tempi] += prob[tempi*max_iter + iter] * (alpha1[j] * y[j] * newbeta[(j * (J + 1) + jj) + (selectedscenario[tempi]* (J + 1) * (J + 1))]);
					}
				}
				Q += (xx[tempi] - xxx[tempi] - xxxx[tempi]);
			}
printf("sub3\n");
			/////updating lower and upper bound
			if (Q < best_upper_bound){
				best_upper_bound = Q;
				LB = LowerBound();
				if (LB > best_lower_bound){                          //Update best lower bound
					best_lower_bound = LB;
				}
			}
			///////updating lagrangean multiplier/////
			double SqrNorm = 0;
			for (j = 0; j < J + 1; j++){
				temps[j] = 0;
				for (tempi = 0; tempi < S2; tempi++){
					temps[j] += prob[tempi*max_iter + iter] * (Hconstraint[selectedscenario[tempi]] * tempx[j*S + selectedscenario[tempi]]);
				}
				SqrNorm += pow((double)temps[j], 2);
			}

			//Check stoping criteria
			if ((best_upper_bound - best_lower_bound) / best_upper_bound * 100 < 0.5)
				break;

			//Compute steplength
			if (iter1 - last > interv_dec_lambda) {
				last = iter1;
				epsilonPrim = epsilonPrim / 2;
			}
			if (epsilonPrim < 0.1)
				epsilonPrim = 2;

			StepLength = epsilonPrim* (best_upper_bound - best_lower_bound) / SqrNorm;

			//Update Lagrange multipliers
			for (j = 0; j < J + 1; j++){
				lambda[j] = lambda[j] + StepLength*temps[j];
			}
		}

		int counter = 0;
		for (j = 0; j < J; j++){
			for (tempi = 0; tempi < S2; tempi++){
				if (tempx[j*S + selectedscenario[tempi]] == tempx[j * S + 0]){
					counter++;
				}
				else
					break;
			}
		}

		filled = 0;
		int count2;
		int identical;
		//set = create_int_vector(100*J* S);

		if (counter != J*S2){
			for (tempi = 0; tempi < S2; tempi++){
				identical = 0;
				for (fillcounter = 0; fillcounter < filled; fillcounter++){
					count2 = 0;
					for (j = 0; j < J; j++){
						if (tempx[j*S + selectedscenario[tempi]] == set[(fillcounter*(J + 1) + j)])
							count2++;
						else
							break;
					}
					if (count2 == J)
						identical = 1;
					if (identical == 1)
						break;
				}
				if (identical == 0){
					filled++;
					for (j = 0; j < J + 1; j++)
						set[(fillcounter*(J + 1) + j)] = tempx[j*S + selectedscenario[tempi]];
				}
			}
			double currentLB = 0;
			for (fillcounter = 0; fillcounter < filled; fillcounter++){

				Qcurrent = 0;
				for (tempi = 0; tempi < S2; tempi++){
					XCPLEX();
				}
printf("sub1\n");
				double *xx;
				double *xxx;
				double *xxxx;
				xx = create_double_vector(S);
				xxx = create_double_vector(S);
				xxxx = create_double_vector(S);
				for (tempi = 0; tempi < S2; tempi++){
					xx[tempi] = 0;
					xxx[tempi] = 0;
					xxxx[tempi] = 0;
				}
				for (tempi = 0; tempi < S2; tempi++){
					for (k = 0; k < K; k++){
						xx[tempi] += prob[tempi*max_iter + iter] * d[k] * beta[(k*S + selectedscenario[tempi])];
					}
					for (i = 0; i < I; i++){
						xxx[tempi] += prob[tempi*max_iter + iter] * alpha[i] * delta[(i*S + selectedscenario[tempi])];
					}
					for (j = 0; j < J + 1; j++){
						for (jj = 0; jj < J + 1; jj++) {
							xxxx[tempi] += prob[tempi*max_iter + iter] * (alpha1[j] * y[j] * newbeta[(j * (J + 1) + jj) + (selectedscenario[tempi]* (J + 1) * (J + 1))]);
						}
					}
					Qcurrent += (xx[tempi] - xxx[tempi] - xxxx[tempi]);
				}
printf("sub2\n");

				if (Qcurrent > currentLB){

printf("sub5\n");
					currentLB = Qcurrent;
					for (tempi = 0; tempi < S2; tempi++){
						for (sGlobal = 0; sGlobal < S; sGlobal++){
							int count = 0;
							for (j = 0; j < jPerim[iter]; j++){
								if (xi2[j*S2 + tempi] == xi[openFacility[j*max_iter + iter] * S + sGlobal]){
									count++;
								}
								else
									break;
							}
printf("sub6\n");
printf("count = %d\n",count);

printf("jperim = %d\n", jPerim[iter]);



							if (count == jPerim[iter]){
printf("sub7\n");
								for (k = 0; k < K; k++){
									beta_Keep[(k*S + sGlobal) + ((iter)*S*K)] = beta[(k*S + selectedscenario[tempi])];
								}
printf("sub8\n");
								for (i = 0; i < I; i++){
									delta_Keep[(i*S + sGlobal) + ((iter)*I*S)] = delta[(i*S + selectedscenario[tempi])];
								}
printf("sub9\n");
								for (j = 0; j < J; j++){
									for (jj = 0; jj < J + 1; jj++) {
										newbeta_Keep[(j * (J + 1) + jj) + (sGlobal * (J + 1) * (J + 1))][iter] = newbeta[(j * (J + 1) + jj) + (selectedscenario[tempi] * (J + 1) * (J + 1))];
									}
printf("sub10\n");
									tempx_keep[j*S + sGlobal] = set[(fillcounter*(J + 1) + j)];
printf("sub11\n");
								}
printf("sub12\n");
							}
printf("sub13\n");
						}
printf("sub14\n");
					}
				}
printf("sub15\n");
			}
printf("sub16\n");
			gettimeofday(&stop, NULL);
			//endTemp = clock();
			//SubTime = (double)(endTemp - startTemp) / (double)(CLOCKS_PER_SEC);
			SubTime = ((double)(stop.tv_sec - start.tv_sec) * 1000 + (double)(stop.tv_usec - start.tv_usec) / 1000) / 1000;
			return currentLB;
		}

		else{
			//for (j = 0; j < J + 1; j++){
			//	for (sGlobal = 0; sGlobal < S; sGlobal++){
			//		tempx_keep[j*S + sGlobal] = tempx[j*S + sGlobal];
			//	}
			//}
			for (tempi = 0; tempi < S2; tempi++){
				for (sGlobal = 0; sGlobal < S; sGlobal++){
					int count = 0;
					for (j = 0; j < jPerim[iter]; j++){
						if (xi2[j*S2 + tempi] == xi[openFacility[j*max_iter + iter] * S + sGlobal] && sGlobal != selectedscenario[tempi]){
							count++;
						}
						else
							break;
					}
					if (count == jPerim[iter]){
						for (k = 0; k < K; k++){
							beta_Keep[(k*S + sGlobal) + ((iter)*S*K)] = beta[(k*S + selectedscenario[tempi])];
						}
						for (i = 0; i < I; i++){
							delta_Keep[(i*S + sGlobal) + ((iter)*I*S)] = delta[(i*S + selectedscenario[tempi])];
						}
						for (j = 0; j < J; j++){
							for (jj = 0; jj < J + 1; jj++) {
										newbeta_Keep[(j * (J + 1) + jj) + (sGlobal * (J + 1) * (J + 1))][iter] = newbeta[(j * (J + 1) + jj) + (selectedscenario[tempi] * (J + 1) * (J + 1))];
									}
							tempx_keep[j*S + sGlobal] = tempx[j*S + selectedscenario[tempi]];
						}
					}
				}
			}
			gettimeofday(&stop, NULL);
			//endTemp = clock();
			//SubTime = (double)(endTemp - startTemp) / (double)(CLOCKS_PER_SEC);
			SubTime = ((double)(stop.tv_sec - start.tv_sec) * 1000 + (double)(stop.tv_usec - start.tv_usec) / 1000) / 1000;

			return Q;
		}
	}
}


/* This simple routine frees up the pointer *ptr, and sets *ptr
to NULL */

static void free_and_null(char **ptr)
{
	if (*ptr != NULL) {
		free(*ptr);
		*ptr = NULL;
	}
} /* END free_and_null */