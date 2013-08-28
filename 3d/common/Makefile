CANS_LIB_DIR = ..
include $(CANS_LIB_DIR)/Makefile_inc
TARGET = $(CANS_LIB_DIR)/libcans2_3d.a
OBJS =  getcurrent.o getcurrent_cyl.o getEta.o cal_resflux.o cal_resflux_fee.o \
	getNewdt_glm.o getNewdt_glmcyl.o  

'.update' : $(OBJS)
	$ ar -r $(TARGET) $(OBJS)
	touch .update

.PHONY : all 
.PHONY : clean
.SUFFIXES :
.SUFFIXES : .o .f90

.f90.o:
	$(FC) $(FFLAGS) -c $<

clean :
	rm -f $(OBJS)
