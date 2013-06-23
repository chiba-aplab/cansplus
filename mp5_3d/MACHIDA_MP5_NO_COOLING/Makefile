#FC = ftn -fastsse -Mipa=fast -O2 -Minline -Minfo=all
FC = mpifrtpx -Kfast
FC = mpifrtpx -Kfast,optmsg=2,simd=2 -Qt -x-
LIB_DIR1 = ..
OBJS =  convert.o openfile.o \
	mpi_domain_xz.o exchangeMpixz_res.o exchangeMpixz2.o \
	integrate_cyl_mp5_1.o integrate_cyl_mp5_2.o integrate_cyl_mp5_3.o \
	integrate_cyl_mp5_1st.o integrate_cyl_mp5_2nd.o \
	perturb.o absorb.o perturb2.o \
	model_machida.o bnd.o bnd_absoub.o \
	MC2.o MC2_dxyz.o switch_MP5toMC2.o switch_MP5toMC2_2.o MP5toMC2.o MP5toMC2_2.o \
	get_temperature.o getNewdt_glmcyl_res_mpi.o main.o 

params.txt : a.out
#	./a.out

a.out : $(OBJS) $(LIB_DIR1)/libcans2_3d.a
	$(FC) $(FFLAGS) -o a.out $(OBJS) \
	-L$(LIB_DIR1) -lcans2_3d -lcansnc

.SUFFIXES:.f90

.f90.o:
	$(FC) -c $<

clean:
	rm -f params.txt *.o *.mod *.oo *.dac.* a.out $(OBJS)
