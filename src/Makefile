# Makefile
#
# author      : Matt Bellis
# date        : Jan., 2012
#

################################################################################
# This section contains some environment variables that you may want to edit, or
# perhaps need to edit, depending on your CUDA installation.
################################################################################

# INSTALL_DIR is where the CUDA and C executables will be installed. You may find it 
# useful to install in some system-wide directory (e.g. /usr/local/bin, 
# $(HOME)/bin, etc.).
INSTALL_DIR       := ../bin/

# SDK_INSTALL_PATH and CUDA_INSTALL_PATH are probably the more significant
# environment variables and will depend on your system and where the CUDA
# libraries and header files are installed.
#SDK_INSTALL_PATH  := ${HOME}/CUDA/SDK_4.0/
#CUDA_INSTALL_PATH := /opt/cuda-5.5/
#SDK_INSTALL_PATH  := ${CUDA_INSTALL_PATH}/samples/
#CUDA_INSTALL_PATH := /usr/local/cuda/
#SDK_INSTALL_PATH  := /usr/include
#CUDA_INSTALL_PATH := /usr/local/cuda-5.0/
#CUDA_INSTALL_PATH := /usr/local/cuda-6.5/
#CUDA_INSTALL_PATH := /usr/local/cuda-5.5/
CUDA_INSTALL_PATH := /opt/apps/cuda/5.5/ # TACC
SDK_INSTALL_PATH  := ${CUDA_INSTALL_PATH}/samples/

# While the INCLUDES and LIBS flags are dependent on the install path, there may 
# still be some subdirectories that are machine-dependent. 
INCLUDES          += -I. -I$(SDK_INSTALL_PATH)/C/common/inc/ -I$(CUDA_INSTALL_PATH)/include/
LIBS              += -L. -L$(SDK_INSTALL_PATH)/C/lib/ -L$(CUDA_INSTALL_PATH)//usr/lib/x86_64-linux-gnu/
LIBS              += -L$(CUDA_INSTALL_PATH)/lib64

LDFLAGS := -lcudart

COMMONFLAGS       := -m64  
OBJ_DIR           := ./objects/
CXXFLAGS          := $(COMMONFLAGS)

# compilers
CC              := g++
NVCC              := $(CUDA_INSTALL_PATH)/bin/nvcc -arch sm_20 $(COMMONFLAGS)

################################################################################

################################################################################
# Rules for Makefile
################################################################################
all: 3pcf_Bellis_totes_different 3pcf_Bellis 3pcf_Bellis_naive 3pcf_Bellis_basically_1D 3pcf_C_version 3pcf_Bellis_2D_histo_on_GPU
#all: 3pcf_C_version

3pcf_Bellis : $(OBJ_DIR)/3pcf_Bellis.cu_o
	mkdir -p $(INSTALL_DIR)
	$(CXX) -fPIC -o $(INSTALL_DIR)/3pcf_Bellis $(OBJ_DIR)//3pcf_Bellis.cu_o $(LDFLAGS) $(INCLUDES) $(LIBS)

$(OBJ_DIR)/3pcf_Bellis.cu_o : 3pcf_Bellis.cu
	mkdir -p $(OBJ_DIR)
	$(NVCC) -g -c $(INCLUDES) -o $(OBJ_DIR)/3pcf_Bellis.cu_o $<

3pcf_Bellis_totes_different : $(OBJ_DIR)/3pcf_Bellis_totes_different.cu_o
	mkdir -p $(INSTALL_DIR)
	$(CXX) -fPIC -o $(INSTALL_DIR)/3pcf_Bellis_totes_different $(OBJ_DIR)//3pcf_Bellis_totes_different.cu_o $(LDFLAGS) $(INCLUDES) $(LIBS)

$(OBJ_DIR)/3pcf_Bellis_totes_different.cu_o : 3pcf_Bellis_totes_different.cu
	mkdir -p $(OBJ_DIR)
	$(NVCC) -g -c $(INCLUDES) -o $(OBJ_DIR)/3pcf_Bellis_totes_different.cu_o $<

3pcf_Bellis_naive : $(OBJ_DIR)/3pcf_Bellis_naive.cu_o
	mkdir -p $(INSTALL_DIR)
	$(CXX) -fPIC -o $(INSTALL_DIR)/3pcf_Bellis_naive $(OBJ_DIR)//3pcf_Bellis_naive.cu_o $(LDFLAGS) $(INCLUDES) $(LIBS)

$(OBJ_DIR)/3pcf_Bellis_naive.cu_o : 3pcf_Bellis_naive.cu
	mkdir -p $(OBJ_DIR)
	$(NVCC) -g -c $(INCLUDES) -o $(OBJ_DIR)/3pcf_Bellis_naive.cu_o $<

3pcf_Bellis_basically_1D : $(OBJ_DIR)/3pcf_Bellis_basically_1D.cu_o
	mkdir -p $(INSTALL_DIR)
	$(CXX) -fPIC -o $(INSTALL_DIR)/3pcf_Bellis_basically_1D $(OBJ_DIR)//3pcf_Bellis_basically_1D.cu_o $(LDFLAGS) $(INCLUDES) $(LIBS)

$(OBJ_DIR)/3pcf_Bellis_basically_1D.cu_o : 3pcf_Bellis_basically_1D.cu
	mkdir -p $(OBJ_DIR)
	$(NVCC) -g -c $(INCLUDES) -o $(OBJ_DIR)/3pcf_Bellis_basically_1D.cu_o $<

3pcf_Bellis_2D_histo_on_GPU : $(OBJ_DIR)/3pcf_Bellis_2D_histo_on_GPU.cu_o
	mkdir -p $(INSTALL_DIR)
	$(CXX) -fPIC -o $(INSTALL_DIR)/3pcf_Bellis_2D_histo_on_GPU $(OBJ_DIR)//3pcf_Bellis_2D_histo_on_GPU.cu_o $(LDFLAGS) $(INCLUDES) $(LIBS)

$(OBJ_DIR)/3pcf_Bellis_2D_histo_on_GPU.cu_o : 3pcf_Bellis_2D_histo_on_GPU.cu
	mkdir -p $(OBJ_DIR)
	$(NVCC) -g -c $(INCLUDES) -o $(OBJ_DIR)/3pcf_Bellis_2D_histo_on_GPU.cu_o $<


3pcf_CUDA_version_triangle_parametrization : $(OBJ_DIR)/3pcf_CUDA_version_triangle_parametrization.cu_o
	mkdir -p $(INSTALL_DIR)
	$(CXX) -fPIC -o $(INSTALL_DIR)/3pcf_CUDA_version_triangle_parametrization $(OBJ_DIR)//3pcf_CUDA_version_triangle_parametrization.cu_o $(LDFLAGS) $(INCLUDES) $(LIBS)

$(OBJ_DIR)/3pcf_CUDA_version_triangle_parametrization.cu_o : 3pcf_CUDA_version_triangle_parametrization.cu
	mkdir -p $(OBJ_DIR)
	$(NVCC) -g -c $(INCLUDES) -o $(OBJ_DIR)/3pcf_CUDA_version_triangle_parametrization.cu_o $<


3pcf_C_version : 3pcf_C_version.cc
	mkdir -p $(INSTALL_DIR)
	$(CC) -O1 -o $(INSTALL_DIR)/$@ $@.cc 

3pcf_C_version_triangle_parametrization : 3pcf_C_version_triangle_parametrization.cc
	mkdir -p $(INSTALL_DIR)
	$(CC) -lm -O1 -g -o $(INSTALL_DIR)/$@ $@.cc 

2pcf_C_version : 2pcf_C_version.cc
	mkdir -p $(INSTALL_DIR)
	$(CC) -O1 -o $(INSTALL_DIR)/$@ $@.cc 

clean:
	rm -f $(INSTALL_DIR)/3pcf_C_version 
	rm -f $(INSTALL_DIR)/3pcf_C_version_triangle_parametrization 
	rm -f $(INSTALL_DIR)/3pcf_Bellis
	rm -f $(INSTALL_DIR)/3pcf_Bellis_totes_different
	rm -f $(INSTALL_DIR)/3pcf_Bellis_naive
	rm -f $(INSTALL_DIR)/3pcf_Bellis_basically_1D
	rm -f $(INSTALL_DIR)/3pcf_Bellis_2D_histo_on_GPU
	rm -f $(OBJ_DIR)/*.o 
	rm -f $(OBJ_DIR)/*.cu_o 

################################################################################
