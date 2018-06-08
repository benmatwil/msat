FC = gfortran

ifeq ($(openmp),off)
else
	ifneq ($(debug),on)
		FOPENMP = -fopenmp
	endif
endif

ifeq ($(debug),on)
	FLAGS = -O0 -g -fbounds-check -Ddebug
else ifeq ($(debug),basic)
	FLAGS = -O0 -g
else
	FLAGS = -O3
endif
FLAGS += -Jmod -Wno-zerotrip

ifeq ($(data),)
else
	DATADIR = $(data)
endif

ifeq ($(output),)
else
	OUTPUTDIR = $(output)
endif

ALLEXE = writedata nf sfxyz sfrpz sfrtp ssfxyz ssfrpz ssfrtp hcs make_cut

SSFFILES = params.f90 src/common.F90 src/trace.F90 src/ring.F90 src/ssf.F90

all: $(ALLEXE) check

sf: sfxyz sfrpz sfrtp

ssf: ssfxyz ssfrpz ssfrtp

cartesian: nf sfxyz ssfxyz

cylindrical: nf sfrpz ssfrpz

spherical: nf sfrtp ssfrtp hcs

xyz: cartesian

rpz: cylindrical

rtp: spherical
	
writedata : params.f90 src/writedata.f90
	$(FC) $(FLAGS) $^ -o $@

nf : params.f90 src/common.F90 src/nf_mod.F90 src/nf.F90
	$(FC) $(FLAGS) $^ -o $@

sfxyz : params.f90 src/common.F90 src/trace.F90 src/sf_mod.F90 src/sf.F90
	$(FC) $(FLAGS) -Dcartesian $(FOPENMP) $^ -o $@

sfrpz : params.f90 src/common.F90 src/trace.F90 src/sf_mod.F90 src/sf.F90
	$(FC) $(FLAGS) -Dcylindrical $(FOPENMP) $^ -o $@

sfrtp : params.f90 src/common.F90 src/trace.F90 src/sf_mod.F90 src/sf.F90
	$(FC) $(FLAGS) -Dspherical $(FOPENMP) $^ -o $@

ssfxyz : $(SSFFILES)
	$(FC) $(FLAGS) -Dcartesian -Dssf_code $(FOPENMP) $^ -o $@

ssfrpz : $(SSFFILES)
	$(FC) $(FLAGS) -Dcylindrical -Dssf_code $(FOPENMP) $^ -o $@

ssfrtp : $(SSFFILES)
	$(FC) $(FLAGS) -Dspherical -Dssf_code $(FOPENMP) $^ -o $@

hcs : params.f90 src/common.F90 src/trace.F90 src/ring.F90 src/hcs.F90
	$(FC) $(FLAGS) -Dspherical -Dssf_code $(FOPENMP) $^ -o $@

bp : params.f90 src/common.F90 src/trace.F90 src/ring.F90 src/bp.F90
	$(FC) $(FLAGS) -Dspherical -Dssf_code $(FOPENMP) $^ -o $@

make_cut : params.f90 src/common.F90 src/make_cut.f90
	$(FC) $(FLAGS) $^ -o $@

###########################################################

check:
#	@echo "Current number of OpenMP threads: $(OMP_NUM_THREADS)"
#	@echo "Using coordinate system: $(coord)"

doc:
	pdflatex --output-directory=doc doc/manual
	@biber doc/manual
	@pdflatex --output-directory=doc doc/manual
	@pdflatex --output-directory=doc doc/manual
	@rm -f doc/*.aux doc/*.log doc/*.out doc/*.xml doc/*.blg doc/*.toc doc/*.bcf doc/*.bbl

setup:
	@rm -f data output
	ln -s -f $(DATADIR) data
	ln -s -f $(OUTPUTDIR) output
	@mkdir -p mod figures

clean:
	@rm -rf mod/*.mod $(ALLEXE)

tidy:
	@rm -rf output/*.dat

.PHONY: check doc

###########################################################

# nfnew : params.f90 nfnew_mod.f90 nfnew.f90
# 	$(FC) $(FLAGS) -g $^ -o $@

# writedata_ws_em : params.f90 writedata_ws_em.f90
# 	$(FC) $(FLAGS) params.f90 writedata_ws_em.f90 -o writedata_ws_em

# writedata_ws_gen : params.f90 writedata_ws_gen.f90
# 	$(FC) $(FLAGS) params.f90 writedata_ws_gen.f90 -o writedata_ws_gen

# sf_gordon : params.f90 gordon/sf_gordon_mod.f90 gordon/sf_gordon.f90
# 	$(FC) $(FLAGS) params.f90 gordon/sf_gordon_mod.f90 gordon/sf_gordon.f90 -o gordon/sf_gordon
