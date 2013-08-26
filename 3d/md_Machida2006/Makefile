#FC = ftn -fastsse -Mipa=fast -O2 -Minline -Minfo=all
#FC = mpifrtpx -Kfast
#FC = mpifrtpx -Kfast,optmsg=2,simd=2 -Qt -x-
FC = mpif90
LIB_DIR1 = ..
OBJS =  getEta2.o perturb.o openfile.o convert.o \
	mpi_domain_xz.o exchangeMpixz_res.o exchangeMpixz.o \
	integrate_cyl.o\
	switch_MP5toMC2.o switch_MP5toMC2_2.o MP5toMC2.o MP5toMC2_2.o \
	model_machida.o bnd.o main.o 

params.txt : a.out
#	./a.out

a.out : $(OBJS) $(LIB_DIR1)/libcans2_3d.a
	$(FC) $(FFLAGS) -o a.out $(OBJS) \
	-L$(LIB_DIR1) -lcans2_3d -lcansnc

.SUFFIXES:.f90

.f90.o:
	$(FC) -c $<

clean:
	rm -f params.txt* *.o *.mod *.oo *.dac.* a.out $(OBJS)
