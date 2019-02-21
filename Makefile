 #F95		= f95
# need f90 instead of f95 to get around bug in GNU f95 that produces
# bad initializers in d1mach()
#F95		= f90


#FC =gfortran
FC= ifort

#FCFLAGS = -check bounds -g -O0 -traceback #Debug flags

FCFLAGS = -g  -xHost  -r8 #-fopenmp #Release flags for ifort -fast -O3
#FCFLAGS = -g  -fno-range-check -fdefault-real-8   -fdefault-double-8 #Release flags for gfortran-fast -O3


PROGRAMS	= DropletCode

RM		= /bin/rm -f

SHELL		= /bin/sh

all:	$(PROGRAMS)

clean:
	-$(RM) *.i
	-$(RM) *.o
	-$(RM) *~
	-$(RM) \#*
	-$(RM) a.out
	-$(RM) core core.*
	-$(RM) *.dat
	-$(RM) *.png

clobber:	distclean

distclean:	mostlyclean
	-$(RM) $(PROGRAMS)

DropletCode: main.o Besseli.o Besselk.o gegenC.o cylindervel.o ellipke2.o fft.o interp.o  intlib3.o pertvelT.o rk4.o stress.o stressG2.o stressI2.o stressA.o trapz1.o trapzN.o wstresscoeffcalculator.o StressCalculatorNew.o NewtonInterp.o #quartic.o
	#$(F95) $(F95FLAGS) $(LDFLAGS) -o $@ $^
	$(FC) $(FCFLAGS) $(LDFFLAGS) -o $@ $^
Besseli.o:Besseli.f90
	#$(F95) $(F95FLAGS) -c $<
	$(FC) $(FCFLAGS) -c  $<
Besselk.o:Besselk.f90
	#$(F95) $(F95FLAGS) -c $<
	$(FC) $(FCFLAGS) -c  $<
gegenC.o:gegenC.f90
	#$(F95)$(F95FLAGS) -c $<
	$(FC) $(FCFLAGS) -c  $<
cylindervel.o:cylindervel.f90
	#$(F95) $(F95FLAGS) -c $<
	$(FC) $(FCFLAGS) -c  $<
ellipke2.o:ellipke2.f90
	#$(F95) $(F95FLAGS) -c $<
	$(FC) $(FCFLAGS) -c  $<
fft.o:fft.f90
	#$(F95) $(F95FLAGS) -c $<
	$(FC) $(FCFLAGS) -c  $<
main.o:main.f90
	#$(F95) $(F95FLAGS) -c $<
	$(FC) $(FCFLAGS) -c  $<
interp.o:interp.f90
	#$(F95) $(F95FLAGS) -c $<
	$(FC) $(FCFLAGS) -c  $<
intlib3.o:intlib3.f90
	#$(F95) $(F95FLAGS) -c $<
	$(FC) $(FCFLAGS) -c  $<
trapz1.o:trapz1.f90
	$(FC) $(FCFLAGS) -c  $<
trapzN.o:trapzN.f90
	$(FC) $(FCFLAGS) -c  $<
pertvelT.o:pertvelT.f90
	#$(F95) $(F95FLAGS) -c $<
	$(FC) $(FCFLAGS) -c  $<
rk4.o:rk4.f90
	#$(F95) $(F95FLAGS) -c $<
	$(FC) $(FCFLAGS) -c  $<
stress.o:stress.f90
	#$(F95) $(F95FLAGS) -c $<
	$(FC) $(FCFLAGS) -c  $<
stressG2.o:stressG2.f90
	#$(F95) $(F95FLAGS) -c $<
	$(FC) $(FCFLAGS) -c  $<
stressI2.o:stressI2.f90
	#$(F95) $(F95FLAGS) -c $<
	$(FC) $(FCFLAGS) -c  $<
stressA.o:stressA.f90
	#$(F95) $(F95FLAGS) -c $<
	$(FC) $(FCFLAGS) -c  $<
wstresscoeffcalculator.o:wstresscoeffcalculator.f90 
	#$(F95) $(F95FLAGS) -c $<
	$(FC) $(FCFLAGS) -c  $<

StressCalculatorNew.o:StressCalculatorNew.f90 
	#$(F95) $(F95FLAGS) -c $<
	$(FC) $(FCFLAGS) -c  $<

NewtonInterp.o:NewtonInterp.f90 
	#$(F95) $(F95FLAGS) -c $<
	$(FC) $(FCFLAGS) -c  $<

#quartic.o:quartic.f90
#	$(FC) $(FCFLAGS) -c  $<
	
maintainer-clean:	distclean
	@echo "This command is intended for maintainers to use;"
	@echo "it deletes files that may require special tools to rebuild."

mostlyclean: clean


