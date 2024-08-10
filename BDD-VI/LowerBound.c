#include "def.h"

extern int         S, S2, iter, J, K, I, Iperim, *keepx, max_iter, iter1, sGlobal, *keepx_Keep, *y, *tempx, H, *tempH, *set, filled, fillcounter, *selectedscenario, tempi, *tempx_keep;
extern double     *x1, theta, theta0, B, *b, *c, *cperim, Q, theta1, *lambda, *temps, Qcurrent;
extern double     MAX_DOUBLE;
extern double     MAX_DOUBLE2;
extern double     *xi, *xi2, *p, *d, *alpha, *alpha1, *delta, *epsilon, *beta, *gamma1, *Delta1, obj, *delta_Keep, *epsilon_Keep, *beta_Keep, *Delta1_Keep, obj1, *prob;
extern double     *vPerim, *v, rho, *f;
extern int        *jPerim, jCounter, *Hconstraint;
extern int        *supplier, *openFacility;
extern double UB;
extern double objective;
extern double		MasterTime, SubTime;
extern clock_t		startTemp, endTemp;


double LowerBound(void){
	int i, j, k, s;
	int index, index1;  // auxiliar indices to fill in the constraint matrix
	double best_upper_bound3, best_lower_bound;
	int nodecount;     //Variables to call cplex
	CPXLPptr  lp4;      // data strucutre to store a problem in cplex ...................
	CPXENVptr env4;     // cplex environment.............................................
	int       numcols; // number of variables ..........................................
	int       numrows; // number of constraints.........................................
	int       numnz;   // number of non-zero elements in the matrix ....................
	int       objsen;  // optimization sense (min:1, max:-1 ) ..........................
	double    *obj;    // objective function coefficients ..............................
	double    *rhs;    // right and side of constraints ................................
	char      *sense;  // constraints sense (<=: 'L', =:'E', >=:'G') ...................
	int       *matbeg; // index of first non-zero element in each row...................
	int       *matind; // associated column of each non-zelo element ...................
	double    *matval; // coefficient values fo the non-zero elements of constraints....
	double    *lb;     // lower bounds of variables.....................................
	double    *ub;     // upper bounds of variables.....................................
	int       status;  // optimization status......................... .................
	double    *x4;      // solution vector (double, even if the problem is integer) .....
	char		probname[16]; // problem name for cplex .......................................
	char      *ctype;  // variable type ('C', 'I', 'B') only if integer.................
	double    value;   // objevtive value of solution ..................................
	int      num_x_var;
	int      *pos_x;
	double		temp1 = 0, temp2 = 0, temp3 = 0, *temp4;
	int     num_Delta1_var;
	int     *pos_Delta1;
	double	*M;
	int sperim;
	//int *Hstar;

	pos_x = create_int_vector((J + 1)* S2);
	M = create_double_vector(J + 1);
	pos_Delta1 = create_int_vector((J + 1)* S2);
	//Hstar = create_int_vector((J + 1)* S);

	/*for (j = 0; j < J + 1; j++){
	M[j] = MAX_DOUBLE2;
	}*/
	//for (s = 0; s < S; s++){
	//	for (j = 0; j < J + 1; j++){
	//		Hstar[j*S + s] = 0;
	//	}
	//}
	//for (s = 0; s < S; s++){
	//	for (j = 0; j < J + 1; j++){
	//		for (sperim = 0; sperim < S - 1; sperim++){
	//			if (s == 0){
	//				Hstar[j*S + s] += 1;
	//			}
	//			else if (sperim == s - 1 & s >0)
	//				Hstar[j*S + s] += -1;
	//			else
	//				Hstar[j*S + s] += 0;
	//			//Hstar[j*S+s] += tempH[(j*(S - 1) + sperim) + (s*(J + 1)*(S - 1))];
	//		}
	//	}
	//}


	//Initialize CPLEX environment
	env4 = CPXopenCPLEX(&status);
	if (env4 == NULL) {
		char  errmsg[1024];
		printf("Could not open CPLEX. \n");
		CPXgeterrorstring(env4, status, errmsg);
		printf("%s", errmsg);
	}

	// Create the problem in CPLEX 
	strcpy(probname, "UFLP");
	lp4 = CPXcreateprob(env4, &status, probname);
	if (env4 == NULL) {
		char  errmsg[1024];
		printf("Could not create LP. \n");
		CPXgeterrorstring(env4, status, errmsg);
		printf("%s", errmsg);
	}

	CPXchgobjsen(env4, lp4, CPX_MAX);

	//Define x variables
	index1 = 0;  // index of columns
	numcols = S2 * (J + 1);
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");


	for (j = 0; j < J + 1; j++){
		for (tempi = 0; tempi < S2; tempi++){
			//pos_x[j*S2 + selectedscenario[tempi]] = index1;
			pos_x[j*S2 + tempi] = index1;
			obj[index1] = 0;
			ctype[index1] = 'B';
			lb[index1] = 0;
			ub[index1] = 1;
			index1++;
		}
	}
	status = CPXnewcols(env4, lp4, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_x_var = index1;


	//Define Delta variables
	index1 = 0;  // index of columns
	numcols = (J + 1) * S2;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (j = 0; j < J + 1; j++){
		for (tempi = 0; tempi < S2; tempi++){
			//sprintf(colname[(int)(index1 + num_x_var + num_beta_var + num_delta_var + num_epsilon_var + num_gamma_var)], "Delta1%3d_%3d", j, s);
			//pos_Delta1[j*S + selectedscenario[tempi]] = index1 + num_x_var;
			pos_Delta1[j*S2 +tempi] = index1 + num_x_var;
			obj[index1] = prob[tempi*max_iter + iter] * alpha1[j] * y[j] * xi[j*S + selectedscenario[tempi]];
			ctype[index1] = 'C';
			lb[index1] = 0;
			ub[index1] = CPX_INFBOUND;
			index1++;
		}
	}

	status = CPXnewcols(env4, lp4, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_Delta1_var = index1;

	//Add budget constraint 
	numrows = S2;
	numnz = S2 * (J + 1);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (tempi = 0; tempi < S2; tempi++){
		sense[index1] = 'L';
		rhs[index1] = B;
		matbeg[index1++] = index;
		for (j = 0; j < J + 1; j++){
			matind[index] = pos_x[j * S2 + tempi];
			matval[index++] = b[j];
		}

	}

	status = CPXaddrows(env4, lp4, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//Add constraint 4 Delta-Mx  
	numrows = S2 * (J + 1);
	numnz = 2 * S2 * (J + 1);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (j = 0; j < J + 1; j++){
		for (tempi = 0; tempi < S2; tempi++){
			sense[index1] = 'L';
			rhs[index1] = 0;
			matbeg[index1++] = index;
			matind[index] = pos_Delta1[j*S2 + tempi];
			matval[index++] = 1;
			matind[index] = pos_x[j*S2 + tempi];
			matval[index++] = -MAX_DOUBLE2;
		}
	}

	status = CPXaddrows(env4, lp4, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//Add constraint 5 Delta-epsilon  
	numrows = S2 * (J + 1);
	numnz = S2 * (J + 1);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (j = 0; j < J + 1; j++){
		for (tempi = 0; tempi < S2; tempi++){
			sense[index1] = 'L';
			rhs[index1] = epsilon[(j * S + selectedscenario[tempi])];
			matbeg[index1++] = index;
			matind[index] = pos_Delta1[j * S2 + tempi];
			matval[index++] = 1;
		}
	}

	status = CPXaddrows(env4, lp4, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);


	//Add constraint 6 Delta - epsilon - Mx >= - M  
	numrows = S2 * (J + 1);
	numnz = 2 * S2 * (J + 1);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (j = 0; j < J + 1; j++){
		for (tempi = 0; tempi < S2; tempi++){
			sense[index1] = 'G';
			rhs[index1] = -MAX_DOUBLE2 + epsilon[(j*S + selectedscenario[tempi])];
			matbeg[index1++] = index;
			matind[index] = pos_Delta1[j*S2 + tempi];
			matval[index++] = 1;
			matind[index] = pos_x[j*S2 + tempi];
			matval[index++] = -MAX_DOUBLE2;
		}
	}
	status = CPXaddrows(env4, lp4, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);


	//Add constraint 7 nonanticipity constraint 
	numrows = (J + 1);
	numnz = S2 * (J + 1);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (j = 0; j < J + 1; j++){
		sense[index1] = 'E';
		rhs[index1] = 0;
		matbeg[index1++] = index;
		for (tempi = 0; tempi < S2; tempi++){
			matind[index] = pos_x[j*S2 + tempi];
			//matval[index++] = Hconstraint[selectedscenario[tempi]];
			if (tempi == 0){
				matval[index++] = S2-1;
			}else
				matval[index++] = - 1;
		}
	}
	status = CPXaddrows(env4, lp4, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);
	//CPXwriteprob(env4, lp4, "model4.lp", NULL);                          //write the model in .lp format if needed (to debug)

	CPXsetintparam(env4, CPX_PARAM_SCRIND, CPX_ON); //output display
	CPXsetintparam(env4, CPX_PARAM_THREADS, 4);
	//CPXsetintparam(env,CPX_PARAM_INTSOLLIM,1);    //stops after finding first integer sol.
	CPXsetintparam(env4, CPX_PARAM_MIPDISPLAY, 3); //different levels of output display
	//CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,1);//0:balanced; 1:feasibility; 2:optimality,3:bestbound, 4:hiddenfeas
	CPXsetdblparam(env4, CPX_PARAM_TILIM, 86400); // time limit
	status = CPXsetintparam(env4, CPX_PARAM_MEMORYEMPHASIS, 1);	//conserve memory where possible
	CPXsetintparam(env4, CPX_PARAM_NODEFILEIND, 0);
	CPXsetdblparam(env4, CPX_PARAM_EPGAP, 0.0001); // e-optimal solution (%gap)
	//CPXsetdblparam(env,CPX_PARAM_TRELIM, 14000); // B&B memory limit
	//CPXsetdblparam(env,CPX_PARAM_EPAGAP, 0.0000000001); // e-optimal solution (absolute value)
	//CPXsetdblparam(env,CPX_PARAM_EPINT, 0.0000000001); // integer precision
	//CPXsetintparam(env,CPX_PARAM_THREADS, 1); // Number of threads to use
	//CPXsetdblparam(env,CPX_PARAM_EPRHS, 0.0000001);
	//CPXsetintparam(env,CPX_PARAM_REDUCE, 0);  // only needed when adding lazy constraints
	//CPXsetintparam(env,CPX_PARAM_HEURFREQ, -1); //heuristic frequency and intensisty 
	//CPXsetdblparam(env, CPX_PARAM_CUTSFACTOR, 1.0);  //limit the number of cuts added by cplex 1.0002
	//CPXsetdblparam(env,CPX_PARAM_CUTUP,UpperBound+.01); // provide an initial upper bound
	//CPXsetintparam(env4, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_OPTIMALITY);  // MIP emphasis: optimality, feasibility, moving best bound
	//CPXsetintparam(env,CPX_PARAM_PARALLELMODE, 1); 
	//CPXsetintparam(env,CPX_PARAM_PREIND,0);
	//CPXsetintparam(env,CPX_PARAM_MIPORDIND,CPX_ON); // Turn on or off the use of priorities on bracnhing variables
	//CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,CPX_MIPEMPHASIS_BESTBOUND);  // MIP emphasis: optimality, feasibility, moving best bound

	CPXmipopt(env4, lp4);  //solve the integer program

	i = CPXgetstat(env4, lp4);
	if (i == 101)
		printf("Optimal solution found\n");
	else if (i == 102)
		printf("e-optimal solution found\n");
	else if (i == 103)
		printf(" infeasible solution\n");
	else if (i == 107)
		printf("Time limit reached\n");
	else
		printf("Unknown stopping criterion (%d)\n", i);

	// retrive solution values
	CPXgetmipobjval(env4, lp4, &value);
	printf("Upper bound: %.2f   ", value);
	best_upper_bound3 = value;
	// If CPLEX was able to find the optimal solution, the previous function provides the optimal solution value
	//if not, it provides the best upper bound
	CPXgetbestobjval(env4, lp4, &value);  //best lower bound in case the problem was not solved to optimality
	best_lower_bound = value;
	printf("Lower bound: %.2f  \n", value);

	nodecount = CPXgetnodecnt(env4, lp4);
	printf("Number of BB nodes : %ld  \n", nodecount);

	numcols = CPXgetnumcols(env4, lp4);
	d_vector(&x4, numcols, "open_cplex:0");
	CPXgetmipx(env4, lp4, x4, 0, numcols - 1);  // obtain the values of the decision variables

	if (lp4 != NULL) {
		status = CPXfreeprob(env4, &lp4);
		if (status) {
			fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
		}
	}
	if (env4 != NULL) {
		status = CPXcloseCPLEX(&env4);
		if (status) {
			char  errmsg[1024];
			fprintf(stderr, "Could not close CPLEX environment.\n");
			CPXgeterrorstring(env4, status, errmsg);
			fprintf(stderr, "%s", errmsg);
		}
	}

	free(pos_x);
	free(pos_Delta1);
	free(M);
	free(x4);


	double *calculatebeta;
	double *calculatedelta;
	double *calculatepsilon;
	double Qcalculate;

	calculatebeta = create_double_vector(S);
	calculatedelta = create_double_vector(S);
	calculatepsilon = create_double_vector(S);
	Qcalculate = 0;

	for (tempi = 0; tempi < S2; tempi++){
		for (k = 0; k < K; k++){
			calculatebeta[tempi] += prob[tempi*max_iter + iter] * d[k] * beta[(k*S + selectedscenario[tempi])];
		}
		for (i = 0; i < I; i++){
			calculatedelta[tempi] += prob[tempi*max_iter + iter] * alpha[i] * delta[(i*S + selectedscenario[tempi])];
		}
		for (j = 0; j < J + 1; j++){
			calculatepsilon[tempi] += prob[tempi*max_iter + iter] * alpha1[j] * y[j] * epsilon[(j*S + selectedscenario[tempi])];
		}
		Qcalculate += (calculatebeta[tempi] - calculatedelta[tempi] - calculatepsilon[tempi]);
	}
printf("end of lowerbound\n");
	//return best_upper_bound ;
	return best_upper_bound3 + Qcalculate;
}
static void free_and_null(char **ptr)
{
	if (*ptr != NULL) {
		free(*ptr);
		*ptr = NULL;
	}
} /* END free_and_null */
