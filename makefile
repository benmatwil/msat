FC = gfortran
FFLAGS = -O3

all: writedata nf sf ssf readin sf_converge
	
writedata : params.f90 writedata.f90
	$(FC) $(FFLAGS) params.f90 writedata.f90 -o writedata

nf : params.f90 nf_mod.f90 nf.f90
	$(FC) $(FFLAGS) params.f90 nf_mod.f90 nf.f90 -o nf

sf : params.f90 sf_mod.f90 sf.f90
	$(FC) $(FFLAGS) params.f90 sf_mod.f90 sf.f90 -o sf
    
sf_converge : params.f90 sf_converge_mod.f90 sf_converge.f90
	$(FC) $(FFLAGS) params.f90 sf_converge_mod.f90 sf_converge.f90 -o sf_converge

ssf : params.f90 common.f90 trace.f90 ring.f90 ssf.f90
	$(FC) $(FFLAGS) -fopenmp params.f90 common.f90 trace.f90 ring.f90 ssf.f90 -o ssf
#	$(FC) $(FFLAGS) -fopenmp -g -fcheck=all -fbounds-check -Wall params.f90 common.f90 trace.f90 ring.f90 ssf.f90 -o ssf
	
readin : params.f90 common.f90 readin.f90
	$(FC) $(FFLAGS) params.f90 common.f90 readin.f90 -o readin

