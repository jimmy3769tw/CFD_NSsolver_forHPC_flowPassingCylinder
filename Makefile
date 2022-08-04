# Nom du compilateur
# CC = icpc
# CC = g++



CC = mpicxx 

STD_FLAG = -std=c++17


OCL_FLAG = -lOpenCL

# !==================== OPT

OPT = $(OMP_CC_ICPC)


AMGCL_OPT = -O3 -std=c++17 -fopenmp -march=native

OMP_CC_ICPC = -O3 -fopenmp -std=c++17 -xHost -qopt-report=5 -qopt-report-phase=vec #mpicpc


OPT = $(STD_FLAG) -fopenmp -O3 -march=native

# OPT = -O3 -g -fopenmp -std=c++17 -march=native -mcmodel=large# for g++(debug)
# OPT = -O3 -fopenmp -std=c++17 -march=native -mcmodel=large# for g++

# !==================== OPT



# !====================INC

Linking =

INC = -I inc/ \
	$(INCBOOST) \
	-I inc/ocl \
	-I inc/matrix \
	-I inc/Dfib \
	-I inc/Debug \
	-I /mnt/c/Users/user/includes/amgcl

INCBOOST = -I /home/Jimmy/installBoost/boost_1_76_0

# !====================INC


# -ipo -xCORE-AVX2 -align array32byte 
# -mcmodel=large
# -heap-arrays 64 
# -mcmodel=medium 


DIR_OBJ = ./obj
DIR_BIN = ./bin






# Defining the objects (OBJS) variables
OBJS = \
	CFD_MX.o \
	RUN_CPU.o \
	1_2_Omp.o \
	1_3_Mpi.o \
	# OutputPLOT3D.o \
	# gridder.o \
	# readingData.o \
	# InitialConditions.o \
	# filer.o \
	# ConvectionScheme.o \
	# PressureMatrix.o \
	# PressureSolvers.o \
	# calNewVelocity.o \
	# Ocl.o \



# 2_0_2_0_GridGenerator

# Linking object files
EXE = mx

exe :  $(OBJS) 
	$(CC) -o $(EXE) $(OBJS) $(OPT) $(Linking)

# echo something
	@echo "    |  Author:  | Chen Jing Ming                  "
	@echo "    |  Version: | 0.1                             "
	@echo "    |  Web:     | http://smetana.me.ntust.edu.tw/ "
	@echo "    |  Editor:  | JM                              "
	@echo "    |  run:     | ./mx                            "

# Defining the flags of objects
%.o: src/%.cpp
	@$(CC) $< $(OPT) $(INC) -c 


%.o: src/Debug/%.cpp
	@$(CC) $< $(OPT) $(INC) -c 



%.o: src/%.c
	@$(C) $< $(OPT) $(INC) -c 

# Removeing data file 
.PHONY: cleanData
clean :
	@/bin/rm -f mx_out/*.x
	@/bin/rm -f mx_out/*.q
	@/bin/rm -f mx_Diffrent/*.x
	@/bin/rm -f mx_Diffrent/*.q
	@/bin/rm -f Information/*.dat
	@/bin/rm -f Information/Uprofile/*.dat
	@/bin/rm -f *.optrpt
	@echo "clean"


cO :
	@/bin/rm -f *.optrpt
	@echo "clean opt"



# Removing object files
.PHONY: cleanall
cleanall : 
	@/bin/rm -f mx_out/*.x
	@/bin/rm -f mx_out/*.q
	@/bin/rm -f mx_Diffrent/*.x
	@/bin/rm -f mx_Diffrent/*.q
	@/bin/rm -f Information/*.csv
	@/bin/rm -f Information/*.dat
	@/bin/rm -f Information/Vprofile/*.dat
	@/bin/rm -f Information/Uprofile/*.dat
	@/bin/rm -f *.optrpt
	@/bin/rm -f Information/Chronograph/*.dat
	@/bin/rm -f Information/Chronograph/*.csv
	@/bin/rm -f $(OBJS) $(EXE)  *.mod
	@echo "cleanall"

config :
	if [ ! -d obj ] ; then mkdir obj ; fi
	if [ ! -d run ] ; then mkdir bin ; fi



nohup :
	nohup ./mx &>> nohup.out &



.PHONY: mv
mv :
	mv *x *q ../mx_in/



.PHONY: perf-test
perf-test :
	

.PHONY: test
test : 
	valgrind ./mx

# .PHONY:plot 
# plot :
# 	./
