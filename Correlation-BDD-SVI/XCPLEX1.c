#include "def.h"

extern int         S, S2, iter, J, K, I, Iperim, *keepx, max_iter, iter1, sGlobal, *keepx_Keep, *y, *tempx, H, *tempH, *set, filled, fillcounter, *selectedscenario, tempi, *tempx_keep;
extern double     *x1, theta, theta0, B, *b, *c, *cperim, Q, theta1, *lambda, *temps, Qcurrent, * correlation, * tau, ** newbeta_Keep, *newbeta;
extern double     MAX_DOUBLE;
extern double     MAX_DOUBLE2;
extern double     *xi, *xi2, *p, *d, *alpha, *alpha1, *delta, *epsilon, *beta, *gamma1, *Delta1, obj, *delta_Keep, *epsilon_Keep, *beta_Keep, *Delta1_Keep, obj1, *prob;
extern double     *vPerim, *v, rho, *f;
extern int        *jPerim, jCounter;
extern int        *supplier, *openFacility;
extern double UB;
extern double objective;
extern double		MasterTime, SubTime;
extern clock_t		startTemp, endTemp;

double XCPLEX1(void)
{
	int i, j, k, s,jj;
	int index, index1;  // auxiliar indices to fill in the constraint matrix
	double best_upper_bound2, best_lower_bound;
	int nodecount;     //Variables to call cplex
	CPXLPptr  lp7;      // data strucutre to store a problem in cplex ...................
	CPXENVptr env7;     // cplex environment.............................................
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
	double    *x5;      // solution vector (double, even if the problem is integer) .....
	char		probname[16]; // problem name for cplex .......................................
	char      *ctype;  // variable type ('C', 'I', 'B') only if integer.................
	double    value;   // objevtive value of solution ..................................
	int      num_beta_var, num_delta_var, num_tau_var, num_gamma_var, num_newbeta_var;
	int		    *pos_beta, *pos_delta, *pos_tau, *pos_gamma, *pos_newbeta;
	double		*M;
	double		lagrangeanResult;
	int sperim;

	d_vector(&M, J + 1, "open_cplex:1");
	i_vector(&pos_beta, K, "open_cplex:1");
	i_vector(&pos_delta, I, "open_cplex:1");
	i_vector(&pos_tau, (J + 1) * (J + 1), "open_cplex:1");
	i_vector(&pos_gamma, J + 1, "open_cplex:1");
	i_vector(&pos_newbeta, (J + 1) * (J + 1), "open_cplex:1");


	for (j = 0; j < J + 1; j++){
	M[j] = MAX_DOUBLE2;
	}
	
	//Initialize CPLEX environment
	env7 = CPXopenCPLEX(&status);
	if (env7 == NULL) {
		char  errmsg[1024];
		printf("Could not open CPLEX. \n");
		CPXgeterrorstring(env7, status, errmsg);
		printf("%s", errmsg);
	}

	// Create the problem in CPLEX 
	strcpy(probname, "UFLP");
	lp7 = CPXcreateprob(env7, &status, probname);
	if (env7 == NULL) {
		char  errmsg[1024];
		printf("Could not create LP. \n");
		CPXgeterrorstring(env7, status, errmsg);
		printf("%s", errmsg);
	}

	CPXchgobjsen(env7, lp7, CPX_MAX);

	//Define beta variables
	index1 = 0;  // index of columns
	numcols = K;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (k = 0; k < K; k++){
		pos_beta[k] = index1;
		obj[index1] = d[k];
		ctype[index1] = 'C';
		lb[index1] = 0;
		ub[index1] = CPX_INFBOUND;
		index1++;
	}

	status = CPXnewcols(env7, lp7, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_beta_var = index1;

	//Define delta variables
	index1 = 0;  // index of columns
	numcols = I;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");
	for (i = 0; i < I; i++){
		pos_delta[i] = index1 + num_beta_var;
		obj[index1] = -1 * (alpha[i]);
		ctype[index1] = 'C';
		lb[index1] = 0;
		ub[index1] = CPX_INFBOUND;
		index1++;
	}
	status = CPXnewcols(env7, lp7, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_delta_var = index1;

	////tau
	index1 = 0;  // index of columns
	numcols = (J + 1) * (J + 1);
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (j = 0; j < J + 1; j++) {
		for (jj = 0; jj < J + 1; jj++) {
			pos_tau[(j * (J + 1) + jj)] = index1 + num_beta_var + num_delta_var;
			obj[index1] = 0;
			ctype[index1] = 'C';
			lb[index1] = 0;
			ub[index1] = CPX_INFBOUND;
			index1++;
		}
	}

	status = CPXnewcols(env7, lp7, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_tau_var = index1;

	//Define gamma variables
	index1 = 0;  // index of columns
	numcols = (J + 1);
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (j = 0; j < J + 1; j++){
		pos_gamma[j] = index1 + num_beta_var + num_delta_var + num_tau_var;
		obj[index1] = 0;
		ctype[index1] = 'C';
		lb[index1] = 0;
		ub[index1] = CPX_INFBOUND;
		index1++;
	}

	status = CPXnewcols(env7, lp7, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_gamma_var = index1;

	//Define newbeta variables
	index1 = 0;  // index of columns
	numcols = (J + 1) * (J+1);
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (j = 0; j < J + 1; j++) {
		for (jj = 0; jj < J + 1; jj++) {
				pos_newbeta[(j * (J + 1) + jj)] = index1 + num_beta_var + num_delta_var + num_tau_var + num_gamma_var;
				obj[index1] = alpha1[j] * y[j] * -1;
				ctype[index1] = 'C';
				lb[index1] = 0;
				ub[index1] = CPX_INFBOUND;
				index1++;
		}
	}

	status = CPXnewcols(env7, lp7, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_newbeta_var = index1;


	//Add constraint beta-gamma
	numrows = (J + 1) * K;
	numnz = 2 * (J + 1) * K;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (j = 0; j < J + 1; j++){
		for (k = 0; k < K; k++){
			sense[index1] = 'L';
			rhs[index1] = c[j*K + k];
			matbeg[index1++] = index;
			matind[index] = pos_beta[k];
			matval[index++] = 1;
			matind[index] = pos_gamma[j];
			matval[index++] = -1;
		}
	}

	status = CPXaddrows(env7, lp7, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//Add constraint 2 gamma-delta-epsilon  
	numrows = (J + 1) * I;
	numnz = 3 * (J + 1) * I  * (J + 1);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (i = 0; i < I; i++){
		for (j = 0; j < J + 1; j++){
			sense[index1] = 'L';
			rhs[index1] = cperim[i * (J + 1) + j];
			matbeg[index1++] = index;
			matind[index] = pos_gamma[j];
			matval[index++] = 1;
			matind[index] = pos_delta[i];
			matval[index++] = -1;
			for (jj = 0; jj < J + 1; jj++) {
				matind[index] = pos_tau[(j * (J + 1) + jj)];
				matval[index++] = -1;
				}
		}
	}

	status = CPXaddrows(env7, lp7, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//Add constraint 35  
	numrows =(J + 1) * (J + 1);
	numnz = 2 * (J + 1) * (J + 1);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (j = 0; j < J + 1; j++) {
		for (jj = 0; jj < J + 1; jj++) {
			sense[index1] = 'G';
			rhs[index1] = 0;
			matbeg[index1++] = index;
			matind[index] = pos_newbeta[(j * (J + 1) + jj)];
			matval[index++] = 1;
			matind[index] = pos_tau[(j * (J + 1) + jj)];
			matval[index++] = (1 - correlation[j * (J + 1) + jj]) * -1;
		}
	}

	status = CPXaddrows(env7, lp7, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//Add constraint 36  
	numrows = (J + 1) * (J + 1);
	numnz = (2 * (J + 1) * (J + 1) + (J + 1)) * 2;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (j = 0; j < J + 1; j++) {
		for (jj = 0; jj < J + 1; jj++) {
			sense[index1] = 'G';
			rhs[index1] = -set[(fillcounter*(J + 1) + j)] * M[j] *  xi[jj * S + sGlobal];
			matbeg[index1++] = index;
			matind[index] = pos_newbeta[(j * (J + 1) + jj)];
			matval[index++] = 1;
			matind[index] = pos_tau[(j * (J + 1) + jj)];
			matval[index++] = -1;
		}
	}

	status = CPXaddrows(env7, lp7, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);
	


	//CPXwriteprob(env7, lp7, "model5.lp", NULL);                          //write the model in .lp format if needed (to debug)

	CPXsetintparam(env7, CPX_PARAM_SCRIND, CPX_ON); //output display
	CPXsetintparam(env7, CPX_PARAM_THREADS, 4);
	//CPXsetintparam(env,CPX_PARAM_INTSOLLIM,1);    //stops after finding first integer sol.
	CPXsetintparam(env7, CPX_PARAM_MIPDISPLAY, 3); //different levels of output display
	//CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,1);//0:balanced; 1:feasibility; 2:optimality,3:bestbound, 4:hiddenfeas
	CPXsetdblparam(env7, CPX_PARAM_TILIM, 86400); // time limit
	status = CPXsetintparam(env7, CPX_PARAM_MEMORYEMPHASIS, 1);	//conserve memory where possible
	CPXsetintparam(env7, CPX_PARAM_NODEFILEIND, 0);
	CPXsetdblparam(env7, CPX_PARAM_EPGAP, 0.0001); // e-optimal solution (%gap)
	//CPXsetdblparam(env,CPX_PARAM_TRELIM, 14000); // B&B memory limit
	//CPXsetdblparam(env,CPX_PARAM_EPGAP, 0.0000000001); // e-optimal solution (%gap)
	//CPXsetdblparam(env,CPX_PARAM_EPAGAP, 0.0000000001); // e-optimal solution (absolute value)
	//CPXsetdblparam(env,CPX_PARAM_EPINT, 0.0000000001); // integer precision
	//CPXsetintparam(env,CPX_PARAM_THREADS, 1); // Number of threads to use
	//CPXsetdblparam(env,CPX_PARAM_EPRHS, 0.0000001);
	//CPXsetintparam(env,CPX_PARAM_REDUCE, 0);  // only needed when adding lazy constraints
	//CPXsetintparam(env,CPX_PARAM_HEURFREQ, -1); //heuristic frequency and intensisty 
	//CPXsetdblparam(env, CPX_PARAM_CUTSFACTOR, 1.0);  //limit the number of cuts added by cplex 1.0002
	//CPXsetdblparam(env,CPX_PARAM_CUTUP,UpperBound+.01); // provide an initial upper bound
	//CPXsetintparam(env7, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_OPTIMALITY);  // MIP emphasis: optimality, feasibility, moving best bound
	//CPXsetintparam(env,CPX_PARAM_PARALLELMODE, 1); 
	//CPXsetintparam(env,CPX_PARAM_PREIND,0);
	//CPXsetintparam(env,CPX_PARAM_MIPORDIND,CPX_ON); // Turn on or off the use of priorities on bracnhing variables
	//CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,CPX_MIPEMPHASIS_BESTBOUND);  // MIP emphasis: optimality, feasibility, moving best bound

	CPXmipopt(env7, lp7);  //solve the integer program
	i = CPXgetstat(env7, lp7);
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
	CPXgetmipobjval(env7, lp7, &value);
	printf("Upper bound: %.2f   ", value);
	best_upper_bound2 = value;
	// If CPLEX was able to find the optimal solution, the previous function provides the optimal solution value
	//if not, it provides the best upper bound
	CPXgetbestobjval(env7, lp7, &value);  //best lower bound in case the problem was not solved to optimality
	best_lower_bound = value;
	printf("Lower bound: %.2f  \n", value);
	nodecount = CPXgetnodecnt(env7, lp7);
	printf("Number of BB nodes : %ld  \n", nodecount);

	numcols = CPXgetnumcols(env7, lp7);
	d_vector(&x5, numcols, "open_cplex:0");
	CPXgetmipx(env7, lp7, x5, 0, numcols - 1);  // obtain the values of the decision variables
	if (lp7 != NULL) {
		status = CPXfreeprob(env7, &lp7);
		if (status) {
			fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
		}
	}
	if (env7 != NULL) {
		status = CPXcloseCPLEX(&env7);
		if (status) {
			char  errmsg[1024];
			fprintf(stderr, "Could not close CPLEX environment.\n");
			CPXgeterrorstring(env7, status, errmsg);
			fprintf(stderr, "%s", errmsg);
		}
	}

	index = 0;
	for (k = 0; k < K; k++){
		beta[(k*S + sGlobal)] = x5[index];
		index++;
	}
	for (i = 0; i < I; i++){
		delta[(i*S + sGlobal)] = x5[index];
		index++;
	}
	double sumepsi = 0;
	for (j = 0; j < J + 1; j++) {
		for (jj = 0; jj < J + 1; jj++) {
				tau[(j * (J + 1) + jj) + (sGlobal * (J + 1) * (J + 1))] = x5[index];
				index++;
		}
	}

	for (j = 0; j < J + 1; j++){
		gamma1[(j*S + sGlobal)] = x5[index];
		index++;
	}
	for (j = 0; j < J + 1; j++) {
		for (jj = 0; jj < J + 1; jj++) {
				newbeta[(j * (J + 1) + jj) + (sGlobal * (J + 1) * (J + 1))]= x5[index];
				/*if (newbeta[(j * (J + 1) + jj) + (selectedscenario[tempi] * (J + 1) * (J + 1))][iter] > 0) {
					printf("hi");
				}*/
				index++;
		}
	}

	printf("best_upper_bound = %f\n", best_upper_bound2 + 0.0000001*sumepsi);

	free(pos_beta);
	free(pos_tau);
	free(pos_gamma);
	free(pos_newbeta);
	free(x5);
	free(M);
printf("end of XCPLEX1\n");
	//return best_upper_bound2 ;
	return best_upper_bound2 + 0.0000001*sumepsi;
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

