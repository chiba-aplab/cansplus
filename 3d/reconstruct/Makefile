CANS_LIB_DIR = ..
include $(CANS_LIB_DIR)/Makefile_inc
TARGET = $(CANS_LIB_DIR)/libcans2_3d.a
OBJS =  esystem_glmmhd.o reconstructionConstant.o \
	lr_state_MP5_2.o lr_state_MP5_3.o

'.update' : $(OBJS)
	$ ar -r $(TARGET) $(OBJS)
	touch .update

.PHONY : all 
.PHONY : clean
.SUFFIXES : .o .f90

.f90.o:
	$(FC) $(FFLAGS) -c $<

clean :
	rm -f $(OBJS)
