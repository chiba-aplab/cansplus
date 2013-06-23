#FC = ftn -O2 -Minline
#FC = ftn -fastsse -Mipa=fast -O2 -Minline -Minfo=all
FC = mpifrtpx -Kfast
#FC = ifort -parallel -inline
#FC = gfortran
CANS_LIB_DIR = ..

TARGET = $(CANS_LIB_DIR)/libcans2_3d.a
OBJS =  esystem_glmmhd.o getMP5value.o reconstructionConstant.o \
	MP5_reconstruction_charGlmMhd.o MP5_reconstruction_charGlmMhd2.o \
	MP5_reconstruction_charGlmMhd3.o

'.update' : $(OBJS)
	$ ar -r $(TARGET) $(OBJS)
	touch .update

.SUFFIXES:.f90

.f90.o:
	$(FC) -c $<

clean :
	rm -f $(OBJS)
