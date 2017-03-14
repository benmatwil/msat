FC = gfortran

ifeq ($(coord),)
	coord = cartesian
endif
DEFINECOORD = -D$(coord)

ifeq ($(debug),on)
	FLAGS = -O0 -g -fbounds-check
	DEFINEDEBUG = -Ddebug
else
	FLAGS = -O3
endif
FLAGS += -Jmod

all: nf sf ssf hcs bp make_cut #nfnew #sf_gordon
	@echo "Current number of OpenMP threads: $(OMP_NUM_THREADS)"
	@echo "Using coordinate system: $(coord)"

nf : params.f90 nf_mod.F90 nf.F90
	$(FC) $(FLAGS) $(DEFINEDEBUG) $^ -o $@

sf : params.f90 sf_mod.f90 sf.F90
	$(FC) $(FLAGS) $(DEFINEDEBUG) -fopenmp $^ -o $@

ssf : params.f90 common.F90 trace.F90 ring.F90 ssf.F90
	$(FC) $(FLAGS) $(DEFINECOORD) $(DEFINEDEBUG) -fopenmp $^ -o $@

hcs : params.f90 common.F90 trace.F90 ring.F90 hcs.F90
	$(FC) $(FLAGS) $(DEFINECOORD) $(DEFINEDEBUG) -fopenmp $^ -o $@

bp : params.f90 common.F90 trace.F90 ring.F90 bp.F90
	$(FC) $(FLAGS) $(DEFINECOORD) $(DEFINEDEBUG) -fopenmp $^ -o $@

make_cut : common.F90 params.f90 make_cut.f90
	$(FC) $(FLAGS) $(DEFINEDEBUG) $^ -o $@

###########################################################

clean:
	@rm mod/*.mod nf sf ssf hcs bp make_cut

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
