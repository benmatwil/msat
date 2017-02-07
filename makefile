FC = gfortran
FFLAGS = -O3
MODULES = -Jmod

ifneq ($(strip $(coord))),)
	DEFINE = -D$(coord)
else
	DEFINE = 
endif
ifeq ($(strip $(mode)),debug)
	FFLAGS = -O0 -g -fbounds-check
endif

# DEFINE = -D$(COORD)
# DEBUG = -g -fbounds-check

all: nf sf ssf readin #nfnew #sf_gordon

nf : params.f90 nf_mod.f90 nf.f90
	$(FC) $(FFLAGS) $(MODULES) $^ -o $@

sf : params.f90 sf_mod.f90 sf.f90
	$(FC) $(FFLAGS) $(MODULES) $^ -o $@

ssf : params.f90 common.F90 trace.F90 ring.f90 ssf.F90
	$(FC) $(FFLAGS) $(MODULES) $(DEFINE) -fopenmp $^ -o $@
	
readin : params.f90 readin.f90
	$(FC) $(FFLAGS) $(MODULES) $^ -o $@

###########################################################

clean:
	@rm mod/*.mod nf sf ssf readin

tidy:
	@rm output/*.dat

###########################################################

writedata : params.f90 writedata.f90
	$(FC) $(FFLAGS) $(MODULES) $^ -o $@

nfnew : params.f90 nfnew_mod.f90 nfnew.f90
	$(FC) $(FFLAGS) $(MODULES) -g $^ -o $@

writedata_ws : params.f90 writedata_ws.f90
	$(FC) $(FFLAGS) $(MODULES) $^ -o $@

# writedata_ws_em : params.f90 writedata_ws_em.f90
# 	$(FC) $(FFLAGS) $(MODULES) params.f90 writedata_ws_em.f90 -o writedata_ws_em

# writedata_ws_gen : params.f90 writedata_ws_gen.f90
# 	$(FC) $(FFLAGS) $(MODULES) params.f90 writedata_ws_gen.f90 -o writedata_ws_gen

# sf_gordon : params.f90 gordon/sf_gordon_mod.f90 gordon/sf_gordon.f90
# 	$(FC) $(FFLAGS) $(MODULES) params.f90 gordon/sf_gordon_mod.f90 gordon/sf_gordon.f90 -o gordon/sf_gordon
