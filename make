FC= gfortran
FFLAGS = -O3

all: writedata nf sf ssf readin
	
writedata : params.f90 writedata.f90
	$(FC) $(FFLAGS) params.f90 writedata.f90 -o writedata

nf : params.f90 nf_mod.f90 nf.f90
	$(FC) $(FFLAGS) params.f90 nf_mod.f90 nf.f90 -o nf

sf : params.f90 sf_mod.f90 sf.f90
	$(FC) $(FFLAGS) params.f90 sf_mod.f90 sf.f90 -o sf

ssf : params.f90 common.f90 trace.f90 ring.f90 ssf.f90
	$(FC) $(FFLAGS) -fopenmp -fbacktrace params.f90 common.f90 trace.f90 ring.f90 ssf.f90 -o ssf
	
readin : readin.f90
	$(FC) $(FFLAGS) readin.f90 -o readin

