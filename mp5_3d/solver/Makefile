#FC = ftn -O2 -Minline
#FC = ftn -fastsse -Mipa=fast -O2 -Minline -Minfo=all
FC = mpifrtpx -Kfast
#FC = ifort -parallel -inline
#FC = gfortran
CANS_LIB_DIR = ..

TARGET = $(CANS_LIB_DIR)/libcans2_3d.a
OBJS =  hlld_flux.o hll_flux.o hlld_flux2.o hlld_flux3.o hlld_flux_org.o\
	glm_flux.o cal_interface_BP.o hancock_diff.o hancock_Predictor.o\
	hancock_diff_e.o lr_state_minmod.o lr_state_MC.o lr_state_1D.o \
	lr_state_MC2.o lr_state_MC2_dxyz.o \


'.update' : $(OBJS)
	$ ar -r $(TARGET) $(OBJS)
	touch .update

.SUFFIXES:.f90

.f90.o:
	$(FC) -c $<

clean :
	rm -f $(OBJS)
