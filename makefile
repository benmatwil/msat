FC = gfortran
FFLAGS = -O3
MODULES = -Jmod
DEFINE = -D$(COORD)
DEBUG = -g -fbounds-check

NF = params.f90 nf_mod.f90 nf.f90
SF = params.f90 sf_mod.f90 sf.f90
SSF = params.f90 common.F90 trace.F90 ring.f90 ssf.F90
RI = params.f90 readin.f90

all: nf sf ssf readin #nfnew #sf_gordon

nf : $(NF)
	$(FC) $(FFLAGS) $(MODULES) $(NF) -o nf

sf : $(SF)
	$(FC) $(FFLAGS) $(MODULES) -g $(SF) -o sf

ssf : $(SSF)
	$(FC) $(FFLAGS) $(MODULES) $(DEFINE) $(SSF) -fopenmp -o ssf
	
readin : $(RI)
	$(FC) $(FFLAGS) $(MODULES) $(RI) -o readin

writedata : params.f90 writedata.f90
	$(FC) $(FFLAGS) $(MODULES) params.f90 writedata.f90 -o writedata

nfnew : params.f90 nfnew_mod.f90 nfnew.f90
	$(FC) $(FFLAGS) $(MODULES) -g params.f90 nfnew_mod.f90 nfnew.f90 -o nfnew

clean:
	@rm mod/*.mod nf sf ssf readin

tidy:
	@rm output/*.dat

writedata_ws : params.f90 writedata_ws.f90
	$(FC) $(FFLAGS) $(MODULES) params.f90 writedata_ws.f90 -o writedata_ws

# writedata_ws_em : params.f90 writedata_ws_em.f90
# 	$(FC) $(FFLAGS) $(MODULES) params.f90 writedata_ws_em.f90 -o writedata_ws_em

# writedata_ws_gen : params.f90 writedata_ws_gen.f90
# 	$(FC) $(FFLAGS) $(MODULES) params.f90 writedata_ws_gen.f90 -o writedata_ws_gen

# sf_gordon : params.f90 gordon/sf_gordon_mod.f90 gordon/sf_gordon.f90
# 	$(FC) $(FFLAGS) $(MODULES) params.f90 gordon/sf_gordon_mod.f90 gordon/sf_gordon.f90 -o gordon/sf_gordon
