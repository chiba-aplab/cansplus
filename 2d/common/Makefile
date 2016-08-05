CANS_LIB_DIR = ..
include $(CANS_LIB_DIR)/Makefile_inc
TARGET = $(CANS_LIB_DIR)/libcans3d.a
OBJS =  getNewdt.o lr_state.o flux_calc.o boundary.o convert.o dac_header.o mpi_setup.o

'.update' : $(OBJS)
	$ ar -r $(TARGET) $(OBJS)
	touch .update

.PHONY : all 
.PHONY : clean
.SUFFIXES :
.SUFFIXES : .o .f90

.f90.o:
	$(FC) $(FFLAGS) -c $<

#Dependencies
boundary.o : mpi_setup.o
getNewdt.o : mpi_setup.o

clean :
	rm -f $(OBJS) *.mod
