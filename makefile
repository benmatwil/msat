FC = gfortran
FFLAGS = -O3
MODULES = -Jmod

all: writedata writedata_ws writedata_ws_em writedata_ws_gen nf sf ssf readin #sf_gordon

nf : params.f90 nf_mod.f90 nf.f90
	$(FC) $(FFLAGS) $(MODULES) params.f90 nf_mod.f90 nf.f90 -o nf

sf : params.f90 sf_mod.f90 sf.f90
	$(FC) $(FFLAGS) $(MODULES) params.f90 sf_mod.f90 sf.f90 -o sf

ssf : params.f90 common.f90 trace.f90 ring.f90 ssf.f90
	$(FC) $(FFLAGS) $(MODULES) -fopenmp params.f90 common.f90 trace.f90 ring.f90 ssf.f90 -o ssf
#	$(FC) $(FFLAGS) -fopenmp -g -fcheck=all -fbounds-check -Wall params.f90 common.f90 trace.f90 ring.f90 ssf.f90 -o ssf
	
readin : params.f90 common.f90 readin.f90
	$(FC) $(FFLAGS) $(MODULES) params.f90 common.f90 readin.f90 -o readin

writedata : params.f90 writedata.f90
	$(FC) $(FFLAGS) $(MODULES) params.f90 writedata.f90 -o writedata

writedata_ws : params.f90 writedata_ws.f90
	$(FC) $(FFLAGS) $(MODULES) params.f90 writedata_ws.f90 -o writedata_ws

writedata_ws_em : params.f90 writedata_ws_em.f90
	$(FC) $(FFLAGS) $(MODULES) params.f90 writedata_ws_em.f90 -o writedata_ws_em

writedata_ws_gen : params.f90 writedata_ws_gen.f90
	$(FC) $(FFLAGS) $(MODULES) params.f90 writedata_ws_gen.f90 -o writedata_ws_gen

#sf_gordon : params.f90 gordon/sf_gordon_mod.f90 gordon/sf_gordon.f90
#	$(FC) $(FFLAGS) $(MODULES) params.f90 gordon/sf_gordon_mod.f90 gordon/sf_gordon.f90 -o gordon/sf_gordon
