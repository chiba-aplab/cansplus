CANS_LIB_DIR = ..
include $(CANS_LIB_DIR)/Makefile_inc
TARGET = a.out
OBJS =  getEta2.o openfile.o const.o model.o init.o \
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
	cd $(CANS_LIB_DIR)/common; make

$(CANS_LIB_DIR)/libcansnc.a : 
	cd $(CANS_LIB_DIR)/cansnc; make

#Dependencies
params.txt : $(TARGET)
init.o : $(CANS_LIB_DIR)/libcans3d.a const.o model.o bnd.o
model.o : $(CANS_LIB_DIR)/libcans3d.a const.o
bnd.o : $(CANS_LIB_DIR)/libcans3d.a
integrate_cyl.o : $(CANS_LIB_DIR)/libcans3d.a bnd.o
main.o : $(CANS_LIB_DIR)/libcans3d.a init.o const.o openfile.o
openfile.o : $(CANS_LIB_DIR)/libcansnc.a 

clean:
	rm -f params.txt* *.o *.mod *.oo *.dac.* *.dat $(TARGET) $(OBJS)
