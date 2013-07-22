FC = ftn -fastsse -Mipa=fast -O2 -Minline -Minfo=all
LIB_DIR1 = ..
OBJS =  convert.o openfile.o \
	mpi_domain.o exchangeMpi.o \
	integrate_mp5_1st.o integrate_mp5_2nd.o \
	model_asahina.o bnd_asahina.o \
	main.o

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
