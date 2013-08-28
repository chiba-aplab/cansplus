CANS_LIB_DIR = ..
include $(CANS_LIB_DIR)/Makefile_inc
TARGET = $(CANS_LIB_DIR)/libcans2_3d.a
OBJS =  boundary.o

'.update' : $(OBJS)
	$ ar -r $(TARGET) $(OBJS)
	touch .update

.PHONY : all 
.PHONY : clean
.SUFFIXES : .o .f90

.f90.o:
	$(FC) $(FFLAGS) -c $<

clean :
	rm -f $(OBJS) *.mod
