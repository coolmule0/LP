# 'make'		build executable file $(TARGET) within $(TARGETDIR)
# 'make clean`		removes all .o and executable files


#Cuda compiler nvcc
NVCC = /usr/local/cuda-8.0/bin/nvcc
#nvcc compiler arguments
NVFLAGS = -arch=sm_52 -O3

#Include directories
INCLUDES = -I$(BOOSTDIR) -I$(GLPKDIR) -I$(LOCALDIR)
#Library paths
LFLAGS = -L$(BOOSTLINKEDDIR) -L$(GLPKLIBDIR)
#Library files
LIBS = -lboost_timer -lboost_chrono -lboost_system -lglpk

##########
#Includes
#boost library include location
BOOSTDIR=/home/john/Documents/codeLibraries/boost_1_65_1/
#GLPK library include location
GLPKDIR=/home/john/Documents/codeLibraries/glpk-4.65/src
#Local library files
LOCALDIR=/home/john/Documents/RGBLP/include

#########
#Libraries
#boost linked libraries directory
BOOSTLINKEDDIR=/home/john/Documents/codeLibraries/boost_1_65_1/stage/lib
#GLPK library directory
GLPKLIBDIR=/home/john/Documents/codeLibraries/glpk-4.65/lib
#boost library names


#build target executable
MAIN = RGBSolve
#executable directory
TARGETDIR = /home/john/Documents/RGBLP

#########
#########
#Make functions

.PHONY: clean

all: $(MAIN)

$(MAIN): main.o FileIO.o
	$(NVCC) $(NVFLAGS) $(INCLUDES) -o $(TARGETDIR)/$(MAIN) main.o FileIO.o $(LFLAGS) $(LIBS)

main.o: main.cu
	$(NVCC) $(NVFLAGS) $(INCLUDES) -o $@ -c $< $(LFLAGS) $(LIBS)

FileIO.o: FileIO.cpp FileIO.hpp
	$(NVCC) $(NVFLAGS) $(INCLUDES) -o $@ -c $< $(LFLAGS) $(LIBS)

clean:
	rm -r main.o $(TARGETDIR)/$(MAIN)
