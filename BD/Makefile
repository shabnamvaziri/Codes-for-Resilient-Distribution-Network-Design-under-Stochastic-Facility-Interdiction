LIB = ../../../../../../encs/pkg/cplex-12.10.0/root/cplex
Benders : main.o createScenario.o ReadData.o solve_MasterProblem_model.o solve_SubProblem_model.o 
	gcc -O2 -I$(LIB)/include main.o createScenario.o ReadData.o solve_MasterProblem_model.o solve_SubProblem_model.o -L$(LIB)/lib/x86-64_linux/static_pic  -lcplex -lm -lpthread -ldl -o Benders

main.o: main.c def.h
	gcc -O2 -c main.c
createScenario.o: createScenario.c def.h
	gcc -O2 -c createScenario.c
ReadData.o: ReadData.c def.h
	gcc -O2 -c ReadData.c
solve_MasterProblem_model.o: solve_MasterProblem_model.c def.h
	gcc -O2 -c solve_MasterProblem_model.c
solve_SubProblem_model.o: solve_SubProblem_model.c def.h
	gcc -O2 -c solve_SubProblem_model.c
clean :
	rm Benders main.o createScenario.o ReadData.o solve_MasterProblem_model.o solve_SubProblem_model.o 
