#include "def.h"

int         S, S2, iter, J, K, I, Iperim, *keepx, max_iter, iter1, sGlobal, *keepx_Keep, *y, *tempx, H, *tempH, *set, filled, fillcounter, *selectedfacility, *selectedscenario, tempi, *tempx_keep;
double     *x1, theta, theta0, B, *b, *c, *cperim, Q, theta1, tetha00, *lambda, *temps, Qcurrent;
double     MAX_DOUBLE = 10000000000;
double     MAX_DOUBLE2 = 10000000;
double     *xi, *xi2, *p, *d, *alpha, *alpha1, *delta, *epsilon, *beta, *gamma1, *Delta1, obj, *delta_Keep, *epsilon_Keep, *beta_Keep, *Delta1_Keep, obj1, *prob;
double     *vPerim, *v, rho, *f;
int        *jPerim, jCounter, *Hconstraint, sumy;
int        *supplier, *openFacility;
double objective;
double UB;
double timeTotal = 0;
double		MasterTime, SubTime;
clock_t		startTemp, endTemp;

char OutName[100];
FILE *Out = NULL;
char outfile[20];

int main(int argc, char *argv[])
{
	char	instance[20];
	char	path[50];
	FILE		*ini;
	clock_t  start, end;
	double obj_value_ch, obj_value_ls, opt_value;
	double cputime;
	int i, j, k, s, MaxNumInst, numSample;
	start = clock();

	if (argc == 1) {
		printf("Error: Input file not specified \n");
		exit(8);
	}
	ini = Open_File(argv[1], "r");
	fscanf(ini, "%d", &MaxNumInst);
	fscanf(ini, "%s", &outfile);

	ttt = time(NULL);
	tm = localtime(&ttt);
	Out = Open_File(outfile, "a+");
	fprintf(Out, "\n %s\n", asctime(tm));
	fclose(Out);

	for (numSample = 1; numSample < MaxNumInst; numSample++){
		fscanf(ini, "%s", &instance);
		sprintf(path, "./Data/");
		strcat(path, instance);

		//start = clock();
		gettimeofday(&startTotal, NULL);
		read_instance(path);

		Out = Open_File(outfile, "a+");					//Writing what instance we are solving
		fprintf(Out, "\n%s | %s ;\t", argv[1], instance);
		fprintf(Out, "%d\t%d\t%d\t%f\t|\t", I, J, K, B);
		fclose(Out);

		theta = 0;
		theta0 = -1 * MAX_DOUBLE;
		UB = 10000;
		createScenario();
		for (s = 0; s < S; s++){
			p[s] = 1;
		}
		for (j = 0; j < J; j++){
			for (s = 0; s < S; s++){
				if (xi[j*S + s] == 1)
					p[s] *= 0.75;
				else
					p[s] *= 0.25;
			}
		}
                Out = Open_File(outfile, "a+");
		fprintf(Out, "%f\t%f\t|\t", 0.75, rho);
		fclose(Out);

		iter = 0;
		for (s = 0; s < S; s++){
			if (s == 0){
				Hconstraint[s] = S - 1;
			}
			else
				Hconstraint[s] = -1;
		}
		while (fabs((UB - theta)) >= 0.1){
			obj = solve_MasterProblem();
			UB = solve_SubProblem();
			Out = Open_File(outfile, "a+");
			fprintf(Out, "%d\t%f\t%f\t%f\t%f\t%f\t|\t", iter + 1, obj, theta, UB, MasterTime, SubTime);
			fclose(Out);
			iter++;
		}
		//end = clock();
		gettimeofday(&stopTotal, NULL);
		//timeTotal = (double)(end - start) / (double)(CLOCKS_PER_SEC);
		timeTotal = ((double)(stopTotal.tv_sec - startTotal.tv_sec) * 1000 + (double)(stopTotal.tv_usec - startTotal.tv_usec) / 1000) / 1000;


		Out = Open_File(outfile, "a+");
		fprintf(Out, "%f\t", timeTotal);
		fprintf(Out, "|\tY:\t");
		for (j = 0; j < J; j++){
			if (y[j] > 0)
			{
				fprintf(Out, "%d\t,\t", j);
			}
		}
		fprintf(Out, "|\tX:\t");
		for (j = 0; j < J + 1; j++){
			if (tempx_keep[j*S + 0]>0){
				fprintf(Out, "%d\t,\t", j);
			}
		}

		fclose(Out);

		free_memory();

	}
	fclose(ini);
}


void Print_solution(void)
{
}
