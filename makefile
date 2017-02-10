FC = gfortran

ifneq ($(coord),)
	DEFINECOORD = -D$(coord)
endif
ifeq ($(mode),debug)
	FLAGS = -O0 -g -fbounds-check
	DEFINEMODE += -D$(mode)
else
	FLAGS = -O3
endif
FLAGS += -Jmod

all: nf sf ssf readin #nfnew #sf_gordon
	@echo "Current number of OpenMP threads: $(OMP_NUM_THREADS)"

nf : params.f90 nf_mod.f90 nf.f90
	$(FC) $(FLAGS) $^ -o $@

sf : params.f90 sf_mod.f90 sf.F90
	$(FC) $(FLAGS) $(DEFINEMODE) -fopenmp $^ -o $@

ssf : params.f90 common.F90 trace.F90 ring.f90 ssf.F90
	$(FC) $(FLAGS) $(DEFINECOORD) $(DEFINEMODE) -fopenmp $^ -o $@
	
readin : params.f90 readin.f90
	$(FC) $(FLAGS) $^ -o $@

###########################################################

clean:
	@rm mod/*.mod nf sf ssf readin

tidy:
	@rm output/*.dat

###########################################################

writedata : params.f90 writedata.f90
	$(FC) $(FLAGS) $^ -o $@

nfnew : params.f90 nfnew_mod.f90 nfnew.f90
	$(FC) $(FLAGS) -g $^ -o $@

writedata_ws : params.f90 writedata_ws.f90
	$(FC) $(FLAGS) $^ -o $@

# writedata_ws_em : params.f90 writedata_ws_em.f90
# 	$(FC) $(FLAGS) params.f90 writedata_ws_em.f90 -o writedata_ws_em

# writedata_ws_gen : params.f90 writedata_ws_gen.f90
# 	$(FC) $(FLAGS) params.f90 writedata_ws_gen.f90 -o writedata_ws_gen

# sf_gordon : params.f90 gordon/sf_gordon_mod.f90 gordon/sf_gordon.f90
# 	$(FC) $(FLAGS) params.f90 gordon/sf_gordon_mod.f90 gordon/sf_gordon.f90 -o gordon/sf_gordon
