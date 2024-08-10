LIB = ../../../../../encs/pkg/cplex-12.10.0/root/cplex
DualDecomposition : main.o createScenario.o createScenario2.o ReadData.o LowerBound.o solve_MasterProblem_model.o solve_SubProblem_model.o sloveCPLEX.o XCPLEX.o solveCPLEX1.o XCPLEX1.o LowerBound1.o
	gcc -O2 -I$(LIB)/include main.o createScenario.o createScenario2.o ReadData.o LowerBound.o solve_MasterProblem_model.o solve_SubProblem_model.o sloveCPLEX.o XCPLEX.o solveCPLEX1.o XCPLEX1.o LowerBound1.o -L$(LIB)/lib/x86-64_linux/static_pic  -lcplex -lm -lpthread -ldl -o DualDecomposition

main.o: main.c def.h
	gcc -O2 -c main.c
createScenario.o: createScenario.c def.h
	gcc -O2 -c createScenario.c
createScenario2.o: createScenario2.c def.h
	gcc -O2 -c createScenario2.c
ReadData.o: ReadData.c def.h
	gcc -O2 -c ReadData.c
LowerBound.o: LowerBound.c def.h
	gcc -O2 -c LowerBound.c
solve_MasterProblem_model.o: solve_MasterProblem_model.c def.h
	gcc -O2 -c solve_MasterProblem_model.c
solve_SubProblem_model.o: solve_SubProblem_model.c def.h
	gcc -O2 -c solve_SubProblem_model.c
sloveCPLEX.o: sloveCPLEX.c def.h
	gcc -O2 -c sloveCPLEX.c
XCPLEX.o: XCPLEX.c def.h
	gcc -O2 -c XCPLEX.c
solveCPLEX1.o: solveCPLEX1.c def.h
	gcc -O2 -c solveCPLEX1.c
XCPLEX1.o: XCPLEX1.c def.h
	gcc -O2 -c XCPLEX1.c
LowerBound1.o: LowerBound1.c def.h
	gcc -O2 -c LowerBound1.c
clean :
	rm DualDecomposition main.o createScenario.o createScenario2.o ReadData.o LowerBound.o solve_MasterProblem_model.o solve_SubProblem_model.o sloveCPLEX.o XCPLEX.o solveCPLEX1.o XCPLEX1.o LowerBound1.o