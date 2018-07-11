CANS_LIB_DIR = ..
include $(CANS_LIB_DIR)/Makefile_inc
TARGET = a.out
OBJS = 	openfile.o const.o model.o init.o integrate.o bnd.o main.o 

.PHONY : all 
.PHONY : clean
.SUFFIXES :
.SUFFIXES : .o .f90

all: $(TARGET)

$(TARGET) : $(OBJS)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJS) -L$(CANS_LIB_DIR) -lcans3d

.f90.o:
	$(FC) $(FFLAGS) -L$(CANS_LIB_DIR) -I$(CANS_LIB_DIR)/common -c $< -lcans3d

$(CANS_LIB_DIR)/libcans3d.a : 
	cd $(CANS_LIB_DIR)/common; make

#Dependencies
init.o : $(CANS_LIB_DIR)/libcans3d.a const.o bnd.o model.o
model.o : $(CANS_LIB_DIR)/libcans3d.a const.o
bnd.o : $(CANS_LIB_DIR)/libcans3d.a const.o
integrate.o : $(CANS_LIB_DIR)/libcans3d.a bnd.o
main.o : $(CANS_LIB_DIR)/libcans3d.a init.o const.o openfile.o integrate.o
openfile.o : $(CANS_LIB_DIR)/libcans3d.a const.o

clean:
	rm -f *.mod data/* $(TARGET) $(OBJS)
