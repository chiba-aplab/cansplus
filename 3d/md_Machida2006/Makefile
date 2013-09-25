CANS_LIB_DIR = ..
include $(CANS_LIB_DIR)/Makefile_inc
TARGET = a.out
OBJS = mpi_setup.o exchangeMpixz.o \
	getEta2.o openfile.o const.o model.o init.o \
	integrate_cyl.o\
	switch_MP5toMC2.o switch_MP5toMC2_2.o MP5toMC2.o MP5toMC2_2.o \
	bnd.o main.o 

.PHONY : all 
.PHONY : clean
.SUFFIXES :
.SUFFIXES : .o .f90

all: $(TARGET)

$(TARGET) : $(OBJS)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJS) -L$(CANS_LIB_DIR) -lcans3d -lcansnc

.f90.o:
	$(FC) $(FFLAGS) -L$(CANS_LIB_DIR) -I$(CANS_LIB_DIR)/common -c $< -lcans3d -lcansnc

$(CANS_LIB_DIR)/libcans3d.a : 
	cd $(CANS_LIB_DIR); make

$(CANS_LIB_DIR)/libcansnc.a : 
	cd $(CANS_LIB_DIR)/cansnc; make

#Dependencies
params.txt : $(TARGET)
init.o : $(CANS_LIB_DIR)/libcans3d.a const.o mpi_setup.o model.o
model.o : const.o mpi_setup.o
bnd.o : $(CANS_LIB_DIR)/libcans3d.a mpi_setup.o
exchangeMpixz.o : mpi_setup.o
integrate_cyl.o : $(CANS_LIB_DIR)/libcans3d.a bnd.o exchangeMpixz.o
main.o : $(CANS_LIB_DIR)/libcans3d.a init.o const.o mpi_setup.o openfile.o
openfile.o : $(CANS_LIB_DIR)/libcansnc.a mpi_setup.o 

clean:
	rm -f params.txt* *.o *.mod *.oo *.dac.* $(TARGET) $(OBJS)
