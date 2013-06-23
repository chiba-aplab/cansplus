#FC = mpifrtpx -fastsse -Mipa=fast -O2 -Minline -Minfo=all
FC = mpifrtpx -Kfast
#FC = ifort -parallel -inline
#FC = gfortran
#FC = ifort -parallel -inline
CANS_LIB_DIR = ..

TARGET = $(CANS_LIB_DIR)/libcans2_3d.a
OBJS =  bd_frex.o bd_frez.o bd_frey.o \
	bd_perx.o bd_pery.o bd_perz.o \
	bd_consx.o bd_consy.o bd_consz.o \
	bd_openbx_cyl.o bd_openbz_cyl.o \
	bd_synpx.o bd_synnx.o \
	bd_synpx_car.o bd_synnx_car.o \
	bd_synpy_car.o bd_synny_car.o \
	bd_synpz_car.o bd_synnz_car.o

'.update' : $(OBJS)
	$ ar -r $(TARGET) $(OBJS)
	touch .update

.SUFFIXES:.f90

.f90.o:
	$(FC) -c $<

clean :
	rm -f $(OBJS)
