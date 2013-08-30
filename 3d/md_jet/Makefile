CANS_LIB_DIR = ..
include $(CANS_LIB_DIR)/Makefile_inc
TARGET = a.out
OBJS = mpi_domain_xz.o exchangeMpixz_res.o exchangeMpixz.o \
	const.o init.o perturb.o openfile.o convert.o \
	integrate_cyl.o\
	switch_MP5toMC2.o switch_MP5toMC2_2.o MP5toMC2.o MP5toMC2_2.o \
	model.o bnd.o main.o 

.PHONY : all 
.PHONY : clean
.SUFFIXES :
.SUFFIXES : .o .f90

all: $(TARGET)

$(TARGET) : $(OBJS)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJS) -L$(CANS_LIB_DIR) -lcans2_3d -lcansnc

.f90.o:
	$(FC) $(FFLAGS) -L$(CANS_LIB_DIR) -I$(CANS_LIB_DIR)/bc -c $< -lcans2_3d -lcansnc

$(CANS_LIB_DIR)/libcans2_3d.a : 
	cd $(CANS_LIB_DIR); make

$(CANS_LIB_DIR)/libcansnc.a : 
	cd $(CANS_LIB_DIR)/cansnc; make

#Dependencies
params.txt : a.out
model.o : $(CANS_LIB_DIR)/libcans2_3d.a mpi_domain_xz.o init.o const.o
bnd.o : $(CANS_LIB_DIR)/libcans2_3d.a mpi_domain_xz.o const.o init.o
exchangeMpixz.o : mpi_domain_xz.o
exchangeMpixz_res.o : mpi_domain_xz.o
integrate_cyl.o : convert.o mpi_domain_xz.o
main.o : mpi_domain_xz.o openfile.o const.o init.o
openfile.o : $(CANS_LIB_DIR)/libcansnc.a mpi_domain_xz.o 
perturb.o : mpi_domain_xz.o
init.o : mpi_domain_xz.o const.o

clean:
	rm -f params.txt* *.o *.mod *.oo *.dac.* a.out $(OBJS)
