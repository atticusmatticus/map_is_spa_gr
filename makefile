
flags = -ftree-vectorize -O3 -fopenmp -ffree-line-length-none

all: map_3D_gr.f08 stringlib.f90
	gfortran map_3D_gr.f08 stringlib.f90 $(flags) -o map_3D_gr.x

clean : 
	rm -f *.o *.mod *.x
