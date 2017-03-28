FC = gfortran

ifeq ($(openmp),off)
else
	FOPENMP = -fopenmp
endif

ifeq ($(coord),)
	coord = cartesian
endif
DEFINECOORD = -D$(coord)

ifeq ($(debug),on)
	FLAGS = -O0 -g -fbounds-check -Ddebug
else
	FLAGS = -O3
endif
FLAGS += -Jmod

ifeq ($(data),)
else
	DATADIR = $(data)
endif

ifeq ($(output),)
else
	OUTPUTDIR = $(output)
endif

ALLEXE = writedata nf sf ssf hcs bp make_cut

all: $(ALLEXE) check
	
writedata : params.f90 src/writedata.f90
	$(FC) $(FLAGS) $^ -o $@

nf : params.f90 src/nf_mod.F90 src/nf.F90
	$(FC) $(FLAGS) $^ -o $@

sf : params.f90 src/sf_mod.f90 src/sf.F90
	$(FC) $(FLAGS) $(FOPENMP) $^ -o $@

ssf : params.f90 src/common.F90 src/trace.F90 src/ring.F90 src/ssf.F90
	$(FC) $(FLAGS) $(DEFINECOORD) $(FOPENMP) $^ -o $@

hcs : params.f90 src/common.F90 src/trace.F90 src/ring.F90 src/hcs.F90
	$(FC) $(FLAGS) $(DEFINECOORD) $(FOPENMP) $^ -o $@

bp : params.f90 src/common.F90 src/trace.F90 src/ring.F90 src/bp.F90
	$(FC) $(FLAGS) $(DEFINECOORD) $(FOPENMP) $^ -o $@

make_cut : params.f90 src/common.F90 src/make_cut.f90
	$(FC) $(FLAGS) $^ -o $@

###########################################################

check:
#	@echo "Current number of OpenMP threads: $(OMP_NUM_THREADS)"
#	@echo "Using coordinate system: $(coord)"

doc: doc/manual.tex
	pdflatex doc/manual
	@pdflatex doc/manual
	@rm -f *.aux *.log *.out

setup:
	@rm -f data output
	ln -s -f $(DATADIR) data
	ln -s -f $(OUTPUTDIR) output
	@mkdir -p mod figures

clean:
	@rm -rf mod/*.mod $(ALLEXE)

tidy:
	@rm -rf output/*.dat

.PHONY: check

###########################################################

nfnew : params.f90 nfnew_mod.f90 nfnew.f90
	$(FC) $(FLAGS) -g $^ -o $@

# writedata_ws_em : params.f90 writedata_ws_em.f90
# 	$(FC) $(FLAGS) params.f90 writedata_ws_em.f90 -o writedata_ws_em

# writedata_ws_gen : params.f90 writedata_ws_gen.f90
# 	$(FC) $(FLAGS) params.f90 writedata_ws_gen.f90 -o writedata_ws_gen

# sf_gordon : params.f90 gordon/sf_gordon_mod.f90 gordon/sf_gordon.f90
# 	$(FC) $(FLAGS) params.f90 gordon/sf_gordon_mod.f90 gordon/sf_gordon.f90 -o gordon/sf_gordon
