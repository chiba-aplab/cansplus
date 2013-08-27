#FC = mpifrtpx -fastsse -Mipa=fast -O2 -Minline -Minfo=all
#FC = mpifrtpx -Kfast
#FC = ifort -parallel -inline
#FC = gfortran
#FC = ifort -parallel -inline
CANS_LIB_DIR = ..

TARGET = $(CANS_LIB_DIR)/libcans2_3d.a
OBJS =  boundary.o

'.update' : $(OBJS)
	$ ar -r $(TARGET) $(OBJS)
	touch .update

.SUFFIXES:.f90

.f90.o:
	$(FC) -c $<

clean :
	rm -f $(OBJS) *.mod
