# 'make'		build executable file $(TARGET) within $(TARGETDIR)
# 'make clean`		removes all .o and executable files


#Cuda compiler nvcc
NVCC = /usr/local/cuda-8.0/bin/nvcc
#nvcc compiler arguments
NVFLAGS = -O3 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_52,code=sm_52

#build target executable
MAIN = RGBLP
#executable directory
TARGETDIR = ../


#Include directories
INCLUDES = -I$(LOCALDIR)

##########
#Includes
#Local library files
LOCALDIR=../include

#########
#########
#Make functions

.PHONY: clean

all: $(MAIN)

$(MAIN): main.o FileIO.o
	$(NVCC) $(NVFLAGS) $(INCLUDES) -o $(TARGETDIR)/$(MAIN) main.o FileIO.o

main.o: main.cu
	$(NVCC) $(NVFLAGS) $(INCLUDES) -o $@ -c $<

FileIO.o: FileIO.cpp FileIO.h
	$(NVCC) $(NVFLAGS) $(INCLUDES) -o $@ -c $<

clean:
	rm main.o FileIO.o $(TARGETDIR)/$(MAIN)
