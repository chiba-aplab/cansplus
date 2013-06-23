#FC = ftn -O2 -Minline
#FC = ftn -fastsse -Mipa=fast -O2 -Minline -Minfo=all
FC = mpifrtpx -Kfast
#FC = ifort -parallel -inline
#FC = gfortran
CANS_LIB_DIR = ..

TARGET = $(CANS_LIB_DIR)/libcans2_3d.a
OBJS =  saveInit_3d.o cal_surface.o chkNaN.o \
	getcurrent.o getcurrent_cyl.o getEta.o getEta2.o cal_resflux.o cal_resflux_fee.o \
	getNewdt_glmcyl_res.o getNewdt_glmcyl_res_mpi.o  \
	getNewdt_glm_res2.o getNewdt_glm.o \
	getNewdt_glm_res_mpi.o

'.update' : $(OBJS)
	$ ar -r $(TARGET) $(OBJS)
	touch .update

.SUFFIXES:.f90

.f90.o:
	$(FC) -c $<

clean :
	rm -f $(OBJS)
