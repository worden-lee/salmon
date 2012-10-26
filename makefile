####################################################################
################ makefile for salmon modeling project ##############
####################################################################

############### first mechanics  ###################################

# where is the vxl library directory?  Default value defined here,
# if it's not provided in an environment variable.
VXL ?= ../../vxl

CXX = g++

CFLAGS += -g -DVNL -I$(VXL)/core -I$(VXL)/core/vnl -I$(VXL)/core/vnl/algo -I$(VXL)/vcl -I$(VXL)/lib

# different versions of g++ need different flags
# this test will probably only work on lw's and amh's computers
GXX_VERSION = $(shell g++ -dumpversion)
ifeq ($(GXX_VERSION),2.96)
CFLAGS += -DNEED_STL_FUNCTION -DGNU_OSTREAM=1
else
CFLAGS += -DGNU_OSTREAM=0
endif

# this one isn't so robust either
# what program to view eps files with
GV = $(shell if [ -e /usr/X11R6/bin/ghostview ]; then echo ghostview; \
	elif [ -e /usr/bin/kghostview ]; then echo kghostview; fi)

LDFLAGS += -L$(VXL)/lib -lvnl_algo -lvnl -lvcl -lnetlib -lv3p_netlib
#LDFLAGS += -L$(VXL)/lib -lvcl -lvnl -lvnl_algo

LD_LIBRARY_PATH = $(VXL)/lib
#LD_LIBRARY_PATH += $(VXL)/lib
export LD_LIBRARY_PATH

LDPATH += $(VXL)/lib
export LDPATH

# SHELL is important for the _sweep target and some others
SHELL = /bin/bash

######## now compiling the program and related rules ##################

default: salmon

SALMON_OBJS = salmon.o random.o population.o simulation.o rand.o genrand2.o
.SECONDARY: salmon #don't delete it if it's made incidentally
salmon: $(SALMON_OBJS)
	$(CXX) $(CFLAGS) -o $@ $(SALMON_OBJS) $(LDFLAGS)

test-rand: test-rand.o rand.o genrand2.o
	$(CXX) $^ -o $@

TARGETS = salmon test-rand

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $<

salmon.o: population.h simulation.h

random.o: random.h

population.o: population.h

simulation.o: simulation.h

.PHONY: clean
clean: clear-out clear-figures
	$(RM) -rf *.o *~ core.* $(TARGETS)

clear clear-out:
	$(RM) -r out

clear-outdir:
	$(RM) -r $(OUTDIR)/*

clear-figures:
	$(RM) -r figures

# % is for instance, eps
clear-%:
	find out -name '*.$*' -exec $(RM) "{}" ";"

#################### parameters to the program ######################

# max age per salmon:
ifeq ($(EXPT),coho)
N = 4
else
N = 5
endif

# nest site competition intensity factor
BETA = 1.7E-4

# maximum number of offspring for a given salmon:
ALPHA = 60.0

# time (in years) to run the simulation:
#T = 10000
#T = 16384
T=1024

# output quantity
Q = t

# varying quantity
VARYING = ac

# threshold for extinction time measurements
#  this is only used if EXTINCTION=yes
#EXTINCTION_THRESHOLD = 1080000

# per-cohort threshold for extinction time measurements
#  this is only used if EXTINCTION=yes and COHORT_EXTINCTION=yes
#COHORT_EXTINCTION_THRESHOLD = 150000

# different defaults in case of extinction experiments
ifeq ($(EXTINCTION),yes)
ifeq ($(VARYING),ac)
#VARIANCE=0.5
#EXTINCTION_THRESHOLD=500000
EXTINCTION_THRESHOLD=1000000
COHORT_EXTINCTION_THRESHOLD = 150000
T=100000
endif
ifeq ($(VARYING),s)
#VARIANCE=0.1
#EXTINCTION_THRESHOLD=50000
EXTINCTION_THRESHOLD=900000
COHORT_EXTINCTION_THRESHOLD = 150000
T=100000
endif
endif

# option to set variance magnitude
# you can set VARIANCE_EXP
# or set VARIANCE
# or set either S_VARIANCE or AC_VARIANCE, depending
ifneq ($(VARIANCE_EXP),)
VARIANCE = 1e$(VARIANCE_EXP)
endif

# mean survivorship proportion:
S_MEAN = 0.85

# variance of survivorship proportion:
ifeq ($(VARYING),s)
ifeq ($(VARIANCE),)
ifeq ($(RECOMMENDED_S_VARIANCE),)
#S_VARIANCE = 0.0025
S_VARIANCE = 0.007225# this is from actual data
VARIANCE = $(S_VARIANCE)
else
S_VARIANCE = $(RECOMMENDED_S_VARIANCE)
VARIANCE = $(S_VARIANCE)
endif
else
S_VARIANCE = $(VARIANCE)
endif
else
S_VARIANCE = 0
endif

# mean ac:
ifeq ($(EXPT),coho)
#AC_MEAN=3
AC_MEAN=2.75
else
AC_MEAN=4
endif

# variance of ac:
ifeq ($(VARYING),ac)
ifeq ($(VARIANCE),)
#AC_VARIANCE = 0.0025
AC_VARIANCE = 0.04# estimated from actual data
VARIANCE = $(AC_VARIANCE)
else
AC_VARIANCE = $(VARIANCE)
endif
else
AC_VARIANCE = 0
endif

# spread of spawning age
ifeq ($(EXPT),coho)
#SIGMA = 0.3
SIGMA = 0.2# makes delta_e~=0.1, delta_l small
else
SIGMA = 0.4
endif

# real eigenvalue for colored noise
# used when $(COLOR)=yes and $(COMPLEX_COLOR)=no
# positive for red, negative for blue, zero for white
A = 0

# complex eigenvalue for colored noise
# used when $(COLOR)=yes and $(COMPLEX_COLOR)=yes
R = 0
THETA = 0

# special case, no delta_l
THREE_D = no

# special case, impulsive noise
IMPULSE = no

# special case, only early ocean mortality varies
EARLY_OCEAN = no

############## internal rules for invoking the program ##################

# construct a long directory for holding output files
# that depends on parameter values
ifeq ($(IMPULSE),yes)
IMP_SUBDIR = /impulse
endif
ifeq ($(THREE_D),yes)
THREE_D_SUBDIR = /3d
endif
ifeq ($(EARLY_OCEAN),yes)
EARLY_SUBDIR = /early-ocean
endif
ifeq ($(EXTINCTION),yes)
ifeq ($(COHORT_EXTINCTION),yes)
EXT_SUBDIR = /COHORT_EXTINCTION_THRESHOLD_$(COHORT_EXTINCTION_THRESHOLD)
else
EXT_SUBDIR = /EXTINCTION_THRESHOLD_$(EXTINCTION_THRESHOLD)
endif
endif
ifeq ($(COLOR),yes)
ifeq ($(COMPLEX_COLOR),yes)
COLOR_SUBDIR = /R_$(R)/THETA_$(THETA)
else
COLOR_SUBDIR = /A_$(A)
endif
endif
OUTDIR = out$(IMP_SUBDIR)$(THREE_D_SUBDIR)$(EARLY_SUBDIR)/N_$(N)/BETA_$(BETA)/ALPHA_$(ALPHA)/S_MEAN_$(S_MEAN)/S_VARIANCE_$(S_VARIANCE)/AC_MEAN_$(AC_MEAN)/AC_VARIANCE_$(AC_VARIANCE)/SIGMA_$(SIGMA)/Q_$(Q)/T_$(T)$(EXT_SUBDIR)$(COLOR_SUBDIR)
#OUTDIR = out$(IMP_SUBDIR)$(THREE_D_SUBDIR)/N_$(N)/BETA_$(BETA)/ALPHA_$(ALPHA)/S_MEAN_$(S_MEAN)/S_VARIANCE_$(S_VARIANCE)/AC_MEAN_$(AC_MEAN)/AC_VARIANCE_$(AC_VARIANCE)/SIGMA_$(SIGMA)/Q_$(Q)/VARYING_$(VARYING)/T_$(T)$(EXT_SUBDIR)$(COLOR_SUBDIR)

# construct command line argument list
# holding the parameter values
SALMON_ARGS = $(N) $(BETA) $(ALPHA) $(T) $(S_MEAN) $(S_VARIANCE) $(AC_MEAN) $(AC_VARIANCE) $(SIGMA) $(Q) $(OUTDIR)
ifeq ($(EXTINCTION),yes)
ifeq ($(COHORT_EXTINCTION),yes)
SALMON_ARGS += -ec $(COHORT_EXTINCTION_THRESHOLD)
else
SALMON_ARGS += -e $(EXTINCTION_THRESHOLD)
endif
endif
ifeq ($(COLOR),yes)
ifeq ($(COMPLEX_COLOR),yes)
SALMON_ARGS += -cp $(R) $(THETA)
else
SALMON_ARGS += -c $(A)
endif
endif
ifeq ($(IMPULSE),yes)
SALMON_ARGS += -imp
endif
ifeq ($(THREE_D),yes)
SALMON_ARGS += -3
endif
ifeq ($(EARLY_OCEAN),yes)
SALMON_ARGS += -se
endif

# call the program and get its output in the right place
run: salmon
	mkdir -p $(OUTDIR)
	./salmon $(SALMON_ARGS)

# rules expressing how the output files come from running the program
MODE_NAMES = mode0 mode1 mode3
ifneq ($(N),4)
MODE_NAMES += mode4
endif
# this list can easily get out of date, may cause some confusion
RUN_PRODUCTS_FN = analytic-extinction-time.out analytic-variance-components.out\
 analytic-variance-components-pairwise.out analytic-variance-contributions.out\
 analytic-variance.out deltas.out eigensystem.txt eigenvalues-mag.out\
 eigenvalues.out eigenvectors.tex forcing.out mifi.out mfr.out output-quantity.out\
 population.out resonances.out transformed-forcing.out transformed-weights.out\
 u.out variance-components.tex variance-objects.txt\
 $(addsuffix .out,$(MODE_NAMES))
# censored products
#  mifi.mag.out resonances.mag.out transformed-forcing.mag.out transformed-weights.mag.out u%.mag.out v%.mag.out
ifeq ($(EXTINCTION),yes)
RUN_PRODUCTS_FN += extinction-objs.txt extinction-time.out
else
RUN_PRODUCTS_FN += analytic-transfer-components.out\
 analytic-transfer-components-parametric.out analytic-transfer.out\
 observed-transfer.out population.fft.mag.out transfer-peaks.out u.fft.mag.out
endif
RUN_PRODUCTS = $(addprefix $(OUTDIR)/, $(RUN_PRODUCTS_FN))

#$(OUTDIR)/%.out $(OUTDIR)/%.tex: salmon
$(RUN_PRODUCTS) : salmon
	$(MAKE) run

# sometimes useful for re-plotting things, or something
touch-products:
	touch $(RUN_PRODUCTS)

run-if-needed: $(OUTDIR)/analytic-variance.out

############## various particular settings ##################

TODO = run

# setting one of these can select some nondefault parameters
coho:
	@$(MAKE) EXPT=coho $(TODO)

chinook:
	@$(MAKE) EXPT=chinook $(TODO)

chinook-small-s:
	$(MAKE) EXPT=chinook S_MEAN=0.2828427 RECOMMENDED_S_VARIANCE=0.01 XTAG=-small-s $(TODO)

coho-small-s:
	$(MAKE) EXPT=coho S_MEAN=0.2828427 RECOMMENDED_S_VARIANCE=0.01 XTAG=-small-s $(TODO)

coho-very-small-s:
	$(MAKE) EXPT=coho S_MEAN=0.2 RECOMMENDED_S_VARIANCE=0.01 XTAG=-very-small-s $(TODO)

check-transfers:
	$(MAKE) LD_LIBRARY_PATH=$(LD_LIBRARY_PATH) LDPATH=$(LDPATH) VARYING=ac AC_VARIANCE=1E-12 run plot-transfers
	$(MAKE) LD_LIBRARY_PATH=$(LD_LIBRARY_PATH) LDPATH=$(LDPATH) VARYING=s S_VARIANCE=1E-12 run plot-transfers

check-cohort-extinction:
	$(MAKE) run save-mode-figures EXTINCTION=yes COHORT_EXTINCTION=yes

3d:
	$(MAKE) THREE_D=yes TAG=-3d coho

ac:
	$(MAKE) LD_LIBRARY_PATH=$(LD_LIBRARY_PATH) LDPATH=$(LDPATH) VARYING=ac $(TODO)

s:
	$(MAKE) LD_LIBRARY_PATH=$(LD_LIBRARY_PATH) LDPATH=$(LDPATH) VARYING=s $(TODO)

extinction: 
	$(MAKE) $(TODO) EXTINCTION=yes

ac-extinction-path-experiment:
#	$(MAKE) extinction epath-plots EXTINCTION=yes VARYING=ac VARIANCE=0.2 EXTINCTION_THRESHOLD=650000 T=1000000
	$(MAKE) extinction epath-plots EXTINCTION=yes VARYING=ac VARIANCE=0.5 EXTINCTION_THRESHOLD=500000 T=100000

s-extinction-path-experiment:
#	$(MAKE) extinction epath-plots EXTINCTION=yes VARYING=s VARIANCE=0.01 EXTINCTION_THRESHOLD=650000 T=100000
	$(MAKE) extinction epath-plots EXTINCTION=yes VARYING=s VARIANCE=0.1 EXTINCTION_THRESHOLD=50000 T=100000

epath-plots: $(OUTDIR)/population.tail.cum.cohorts.eps.gv $(OUTDIR)/population.tail.eps.gv $(OUTDIR)/u.tail.eps.gv


output-variance-gains:
	@echo $(VARYING) $(Q) \
	  `$(MAKE) --no-print-directory print-analytic-variance-gain` \
	  `$(MAKE) --no-print-directory print-observed-variance-gain`

print-observed-variance: ./scripts/variance
	@./scripts/variance $(OUTDIR)/output-quantity.out

print-observed-variance-gain:
	@perl -e "print `./scripts/variance $(OUTDIR)/output-quantity.out` / $(VARIANCE), \"\\n\""

print-analytic-variance: 
	@cat $(OUTDIR)/analytic-variance.out

print-analytic-variance-gain:
	@perl -e "print `cat $(OUTDIR)/analytic-variance.out` / $(VARIANCE), \"\\n\""

print-outdir:
	@echo $(OUTDIR)

6variances:
	$(MAKE) --no-print-directory VARYING=ac Q=t VARIANCE=$(VARIANCE) run-if-needed
	$(MAKE) --no-print-directory VARYING=ac Q=r VARIANCE=$(VARIANCE) run-if-needed
	$(MAKE) --no-print-directory VARYING=ac Q=d VARIANCE=$(VARIANCE) run-if-needed
	$(MAKE) --no-print-directory VARYING=s  Q=t VARIANCE=$(VARIANCE) run-if-needed
	$(MAKE) --no-print-directory VARYING=s  Q=r VARIANCE=$(VARIANCE) run-if-needed
	$(MAKE) --no-print-directory VARYING=s  Q=d VARIANCE=$(VARIANCE) run-if-needed
	$(MAKE) --no-print-directory print-6variances

print-6variances:
	@echo VARYING Q analytic observed
	@$(MAKE) --no-print-directory VARYING=ac Q=t VARIANCE=$(VARIANCE) output-variance-gains
	@$(MAKE) --no-print-directory VARYING=ac Q=r VARIANCE=$(VARIANCE) output-variance-gains
	@$(MAKE) --no-print-directory VARYING=ac Q=d VARIANCE=$(VARIANCE) output-variance-gains
	@$(MAKE) --no-print-directory VARYING=s  Q=t VARIANCE=$(VARIANCE) output-variance-gains
	@$(MAKE) --no-print-directory VARYING=s  Q=r VARIANCE=$(VARIANCE) output-variance-gains
	@$(MAKE) --no-print-directory VARYING=s  Q=d VARIANCE=$(VARIANCE) output-variance-gains

emphases:
	@$(MAKE) coho Q=r
	@$(MAKE) coho Q=t
	@$(MAKE) chinook Q=r
	@$(MAKE) chinook Q=t
	@$(MAKE) print-emphases

print-emphases:
	@echo coho
	@$(MAKE) --no-print-directory coho Q=r TODO=output-emphases
	@$(MAKE) --no-print-directory coho Q=t TODO=output-emphases
	@echo chinook
	@$(MAKE) --no-print-directory chinook Q=r TODO=output-emphases
	@$(MAKE) --no-print-directory chinook Q=t TODO=output-emphases

output-emphases:
	@echo "Q =" $(Q)
	@cat $(OUTDIR)/transformed-weights.mag.out

forcing:
	@$(MAKE) coho VARYING=ac
	@$(MAKE) coho VARYING=s
	@$(MAKE) chinook VARYING=ac
	@$(MAKE) chinook VARYING=s
	@$(MAKE) print-forcing

print-forcing:
	@echo coho
	@$(MAKE) --no-print-directory coho VARYING=ac TODO=output-forcing
	@$(MAKE) --no-print-directory coho VARYING=s TODO=output-forcing
	@echo chinook
	@$(MAKE) --no-print-directory chinook VARYING=ac TODO=output-forcing
	@$(MAKE) --no-print-directory chinook VARYING=s TODO=output-forcing

output-forcing:
	@echo "VARYING =" $(VARYING)
	@cat $(OUTDIR)/transformed-forcing.mag.out

############### shared code for generating figures ###############

# whether figures should be black+white
BW = no
ifeq ($(BW),yes)
BW_ARG = -bw
OC_ARG = -oc
else
BW_ARG = 
endif

ifeq ($(Q),t)
QNAME = total
endif
ifeq ($(Q),r)
QNAME = recruitment
endif
ifeq ($(Q),c)
QNAME = catch
endif
ifeq ($(Q),d)
QNAME = difference
endif

ifeq ($(VARYING),ac)
VARYNAME = age
else
ifeq ($(EARLY_OCEAN),yes)
VARYNAME = early-survivorship
else
VARYNAME = survivorship
endif
endif
# set only if not set
#EXPT ?= "AC_MEAN=$(AC_MEAN)"

########### targets for generating figures for paper 1 ############

## these make figures and table in figures/
##  copying into the latex project directory is done by hand

p01-figures:
	$(MAKE) _p01-figures BW=yes

_p01-figures: 3d-figures fft-figure eigenvalue-figures emphasis-figures \
	forcing-figures mifi-figures mfr-figures transfer-figures noisy-transfer-figures \
	small-s-figures very-small-s-figures no-sigma-figures \
	p01-impulse-figures sweep-eigenvalue-figures tables
# note sweep-eigenvalue-figures should possibly be made separately because it takes
#  a little while
sweep-eigenvalue-figures:
	$(MAKE) _sweep-eigenvalue-figures BW=yes SHOW_CIRCLE=no

_sweep-eigenvalue-figures: \
 sweep-eigenvalue-chinook-figure sweep-eigenvalue-coho-figure

sweep-eigenvalue-chinook-figure:
	$(MAKE) chinook T=1 TODO=save-sweep-eigenvalues

sweep-eigenvalue-coho-figure:
	$(MAKE) coho T=1 TODO=save-sweep-eigenvalues

save-sweep-eigenvalues:
	$(MAKE) sweep-sigma
	@mkdir -p figures
	cp out/sweep-SIGMA/eigenvalues.eps figures/eigenvalues-$(EXPT)-sweep.eps

3d-figures: 
	$(MAKE) 3d TODO=save-3d-figures

save-3d-figures: save-eigenvalues save-eigenvectors

fft-figure:
#	$(MAKE) chinook Q=d VARYING=ac T=1024 TODO=save-fft-figure
	$(MAKE) chinook Q=t VARYING=s S=0.2828427 RECOMMENDED_S_VARIANCE=0.01 STAG=-small-s T=1024 TODO=save-fft-figure

save-fft-figure: $(OUTDIR)/population-fft.eps
	@mkdir -p figures
	cp $< figures/fft-$(EXPT)-$(QNAME)-$(VARYNAME)$(XTAG)$(STAG).eps

xeigenvalue-figures:
	$(MAKE) chinook TODO=save-eigenvalues
	$(MAKE) coho TODO=save-eigenvalues
	$(MAKE) chinook SIGMA=0 XTAG=-no-sigma TODO=save-eigenvalues
	$(MAKE) coho SIGMA=0 XTAG=-no-sigma TODO=save-eigenvalues

eigenvalue-figures:
	$(MAKE) grid GRID_TARGET=save-eigenvalues

# set GRID_TARGET, this will do it for all 5 cases
grid:
	$(MAKE) chinook TODO=$(GRID_TARGET)
	$(MAKE) chinook TODO=$(GRID_TARGET) S_MEAN=0.2828427 RECOMMENDED_S_VARIANCE=0.01 STAG=-small-s
	$(MAKE) coho TODO=$(GRID_TARGET)
	$(MAKE) coho TODO=$(GRID_TARGET) S_MEAN=0.2828427 RECOMMENDED_S_VARIANCE=0.01 STAG=-small-s
	$(MAKE) coho TODO=$(GRID_TARGET) S_MEAN=0.2 RECOMMENDED_S_VARIANCE=0.01 STAG=-very-small-s

save-eigenvalues: $(OUTDIR)/eigenvalues.eps
	@mkdir -p figures
	cp $< figures/eigenvalues-$(EXPT)$(XTAG)$(STAG).eps

emphasis-figures:
	$(MAKE) chinook Q=t TODO=save-emphasis
	$(MAKE) chinook Q=r TODO=save-emphasis
	$(MAKE) chinook Q=c TODO=save-emphasis
	$(MAKE) chinook Q=d TODO=save-emphasis
	$(MAKE) coho Q=t TODO=save-emphasis
	$(MAKE) coho Q=r TODO=save-emphasis
	$(MAKE) coho Q=c TODO=save-emphasis
	$(MAKE) coho Q=d TODO=save-emphasis

xsave-emphasis: $(OUTDIR)/transformed-weights.eps
	@mkdir -p figures
	cp $< figures/emphasis-$(EXPT)-$(QNAME).eps

save-emphasis: $(OUTDIR)/transformed-weights.mag.eps
	@mkdir -p figures
	cp $< figures/emphasis-$(EXPT)-$(QNAME).eps

forcing-figures:
	@mkdir -p figures
	$(MAKE) chinook VARYING=ac TODO=save-forcing
	$(MAKE) chinook VARYING=s TODO=save-forcing
	$(MAKE) chinook VARYING=s EARLY_OCEAN=yes TODO=save-forcing
	$(MAKE) coho VARYING=ac TODO=save-forcing
	$(MAKE) coho VARYING=s TODO=save-forcing
	$(MAKE) coho VARYING=s EARLY_OCEAN=yes TODO=save-forcing

xsave-forcing: $(OUTDIR)/transformed-forcing.eps
	@mkdir -p figures
	cp $< figures/forcing-$(EXPT)-$(VARYNAME).eps

save-forcing: $(OUTDIR)/transformed-forcing.mag.eps
	@mkdir -p figures
	cp $< figures/forcing-$(EXPT)-$(VARYNAME).eps

mifi-figures:
	@mkdir -p figures
	$(MAKE) chinook VARYING=s  Q=t TODO=save-mifi
	$(MAKE) chinook VARYING=s  Q=r TODO=save-mifi
	$(MAKE) chinook VARYING=s  Q=c TODO=save-mifi
	$(MAKE) chinook VARYING=s  Q=d TODO=save-mifi
	$(MAKE) coho VARYING=s  Q=t TODO=save-mifi
	$(MAKE) coho VARYING=s  Q=r TODO=save-mifi
	$(MAKE) coho VARYING=s  Q=c TODO=save-mifi
	$(MAKE) coho VARYING=s  Q=d TODO=save-mifi
	$(MAKE) chinook VARYING=ac Q=t TODO=save-mifi
	$(MAKE) chinook VARYING=ac Q=r TODO=save-mifi
	$(MAKE) chinook VARYING=ac Q=c TODO=save-mifi
	$(MAKE) chinook VARYING=ac Q=d TODO=save-mifi
	$(MAKE) coho VARYING=ac Q=t TODO=save-mifi
	$(MAKE) coho VARYING=ac Q=r TODO=save-mifi
	$(MAKE) coho VARYING=ac Q=c TODO=save-mifi
	$(MAKE) coho VARYING=ac Q=d TODO=save-mifi
	$(MAKE) chinook VARYING=s EARLY_OCEAN=yes Q=t TODO=save-mifi
	$(MAKE) chinook VARYING=s EARLY_OCEAN=yes Q=r TODO=save-mifi
	$(MAKE) chinook VARYING=s EARLY_OCEAN=yes Q=c TODO=save-mifi
	$(MAKE) chinook VARYING=s EARLY_OCEAN=yes Q=d TODO=save-mifi
	$(MAKE) coho VARYING=s EARLY_OCEAN=yes Q=t TODO=save-mifi
	$(MAKE) coho VARYING=s EARLY_OCEAN=yes Q=r TODO=save-mifi
	$(MAKE) coho VARYING=s EARLY_OCEAN=yes Q=c TODO=save-mifi
	$(MAKE) coho VARYING=s EARLY_OCEAN=yes Q=d TODO=save-mifi

xsave-mifi: $(OUTDIR)/mifi.eps
	@mkdir -p figures
	cp $< figures/mifi-$(EXPT)-$(QNAME)-$(VARYNAME).eps

save-mifi: $(OUTDIR)/mifi.mag.eps
	@mkdir -p figures
	cp $< figures/mifi-$(EXPT)-$(QNAME)-$(VARYNAME)$(TRTAG).eps

mfr-figures:
	@mkdir -p figures
	$(MAKE) chinook VARYING=s  Q=t TODO=save-mfr
	$(MAKE) chinook VARYING=s  Q=r TODO=save-mfr
	$(MAKE) chinook VARYING=s  Q=c TODO=save-mfr
	$(MAKE) chinook VARYING=s  Q=d TODO=save-mfr
	$(MAKE) coho VARYING=s  Q=t TODO=save-mfr
	$(MAKE) coho VARYING=s  Q=r TODO=save-mfr
	$(MAKE) coho VARYING=s  Q=c TODO=save-mfr
	$(MAKE) coho VARYING=s  Q=d TODO=save-mfr
	$(MAKE) chinook VARYING=ac Q=t TODO=save-mfr
	$(MAKE) chinook VARYING=ac Q=r TODO=save-mfr
	$(MAKE) chinook VARYING=ac Q=c TODO=save-mfr
	$(MAKE) chinook VARYING=ac Q=d TODO=save-mfr
	$(MAKE) coho VARYING=ac Q=t TODO=save-mfr
	$(MAKE) coho VARYING=ac Q=r TODO=save-mfr
	$(MAKE) coho VARYING=ac Q=c TODO=save-mfr
	$(MAKE) coho VARYING=ac Q=d TODO=save-mfr
	$(MAKE) chinook VARYING=s EARLY_OCEAN=yes Q=t TODO=save-mfr
	$(MAKE) chinook VARYING=s EARLY_OCEAN=yes Q=r TODO=save-mfr
	$(MAKE) chinook VARYING=s EARLY_OCEAN=yes Q=c TODO=save-mfr
	$(MAKE) chinook VARYING=s EARLY_OCEAN=yes Q=d TODO=save-mfr
	$(MAKE) coho VARYING=s EARLY_OCEAN=yes Q=t TODO=save-mfr
	$(MAKE) coho VARYING=s EARLY_OCEAN=yes Q=r TODO=save-mfr
	$(MAKE) coho VARYING=s EARLY_OCEAN=yes Q=c TODO=save-mfr
	$(MAKE) coho VARYING=s EARLY_OCEAN=yes Q=d TODO=save-mfr

xsave-mfr: $(OUTDIR)/mfr.eps
	@mkdir -p figures
	cp $< figures/mfr-$(EXPT)-$(QNAME)-$(VARYNAME).eps

save-mfr: $(OUTDIR)/mfr.mag.eps
	@mkdir -p figures
	cp $< figures/mfr-$(EXPT)-$(QNAME)-$(VARYNAME)$(TRTAG).eps

transfer-figures:
	@mkdir -p figures
	$(MAKE) chinook VARYING=ac Q=t TODO=save-transfer CLIP_FACTOR=2
	$(MAKE) chinook VARYING=ac Q=r TODO=save-transfer CLIP_FACTOR=2.5
	$(MAKE) chinook VARYING=ac Q=c TODO=save-transfer CLIP_FACTOR=2
	$(MAKE) chinook VARYING=ac Q=d TODO=save-transfer CLIP_FACTOR=2.5
	$(MAKE) chinook VARYING=s  Q=t TODO=save-transfer
	$(MAKE) chinook VARYING=s  Q=r TODO=save-transfer
	$(MAKE) chinook VARYING=s  Q=c TODO=save-transfer
	$(MAKE) chinook VARYING=s  Q=d TODO=save-transfer
	$(MAKE) chinook VARYING=s  Q=t EARLY_OCEAN=yes TODO=save-transfer
	$(MAKE) chinook VARYING=s  Q=r EARLY_OCEAN=yes TODO=save-transfer
	$(MAKE) chinook VARYING=s  Q=c EARLY_OCEAN=yes TODO=save-transfer
	$(MAKE) chinook VARYING=s  Q=d EARLY_OCEAN=yes TODO=save-transfer
	$(MAKE) coho VARYING=ac Q=t TODO=save-transfer CLIP_FACTOR=3
	$(MAKE) coho VARYING=ac Q=r TODO=save-transfer CLIP_FACTOR=5
	$(MAKE) coho VARYING=ac Q=c TODO=save-transfer CLIP_FACTOR=2.5
	$(MAKE) coho VARYING=ac Q=d TODO=save-transfer CLIP_FACTOR=4
	$(MAKE) coho VARYING=s  Q=t TODO=save-transfer
	$(MAKE) coho VARYING=s  Q=r TODO=save-transfer
	$(MAKE) coho VARYING=s  Q=c TODO=save-transfer
	$(MAKE) coho VARYING=s  Q=d TODO=save-transfer
	$(MAKE) coho VARYING=s  Q=t EARLY_OCEAN=yes TODO=save-transfer
	$(MAKE) coho VARYING=s  Q=r EARLY_OCEAN=yes TODO=save-transfer
	$(MAKE) coho VARYING=s  Q=c EARLY_OCEAN=yes TODO=save-transfer
	$(MAKE) coho VARYING=s  Q=d EARLY_OCEAN=yes TODO=save-transfer

NOISY_AC_VARIANCE_CHINOOK=0.04# 4 * 5 percent squared
NOISY_AC_VARIANCE_COHO=0.0189# 2.75 * 5 percent, squared
NOISY_S_VARIANCE=0.007225# 0.085 squared
NOISY_EARLY_S_VARIANCE=$(NOISY_S_VARIANCE)
noisy-transfer-figures:
	@mkdir -p figures
	$(MAKE) chinook VARYING=ac AC_VARIANCE=$(NOISY_AC_VARIANCE_CHINOOK) Q=t TODO=save-noisy-transfer
	$(MAKE) chinook VARYING=ac AC_VARIANCE=$(NOISY_AC_VARIANCE_CHINOOK) Q=r TODO=save-noisy-transfer
	$(MAKE) chinook VARYING=ac AC_VARIANCE=$(NOISY_AC_VARIANCE_CHINOOK) Q=c TODO=save-noisy-transfer
	$(MAKE) chinook VARYING=ac AC_VARIANCE=$(NOISY_AC_VARIANCE_CHINOOK) Q=d TODO=save-noisy-transfer
	$(MAKE) chinook VARYING=s S_VARIANCE=$(NOISY_S_VARIANCE) Q=t TODO=save-noisy-transfer
	$(MAKE) chinook VARYING=s S_VARIANCE=$(NOISY_S_VARIANCE) Q=r TODO=save-noisy-transfer
	$(MAKE) chinook VARYING=s S_VARIANCE=$(NOISY_S_VARIANCE) Q=c TODO=save-noisy-transfer
	$(MAKE) chinook VARYING=s S_VARIANCE=$(NOISY_S_VARIANCE) Q=d TODO=save-noisy-transfer
	$(MAKE) chinook VARYING=s S_VARIANCE=$(NOISY_EARLY_S_VARIANCE) Q=t EARLY_OCEAN=yes TODO=save-noisy-transfer
	$(MAKE) chinook VARYING=s S_VARIANCE=$(NOISY_EARLY_S_VARIANCE) Q=r EARLY_OCEAN=yes TODO=save-noisy-transfer
	$(MAKE) chinook VARYING=s S_VARIANCE=$(NOISY_EARLY_S_VARIANCE) Q=c EARLY_OCEAN=yes TODO=save-noisy-transfer
	$(MAKE) chinook VARYING=s S_VARIANCE=$(NOISY_EARLY_S_VARIANCE) Q=d EARLY_OCEAN=yes TODO=save-noisy-transfer
	$(MAKE) coho VARYING=ac AC_VARIANCE=$(NOISY_AC_VARIANCE_COHO) Q=t TODO=save-noisy-transfer
	$(MAKE) coho VARYING=ac AC_VARIANCE=$(NOISY_AC_VARIANCE_COHO) Q=r TODO=save-noisy-transfer
	$(MAKE) coho VARYING=ac AC_VARIANCE=$(NOISY_AC_VARIANCE_COHO) Q=c TODO=save-noisy-transfer
	$(MAKE) coho VARYING=ac AC_VARIANCE=$(NOISY_AC_VARIANCE_COHO) Q=d TODO=save-noisy-transfer
	$(MAKE) coho VARYING=s S_VARIANCE=$(NOISY_S_VARIANCE) Q=t TODO=save-noisy-transfer
	$(MAKE) coho VARYING=s S_VARIANCE=$(NOISY_S_VARIANCE) Q=r TODO=save-noisy-transfer
	$(MAKE) coho VARYING=s S_VARIANCE=$(NOISY_S_VARIANCE) Q=c TODO=save-noisy-transfer
	$(MAKE) coho VARYING=s S_VARIANCE=$(NOISY_S_VARIANCE) Q=d TODO=save-noisy-transfer
	$(MAKE) coho VARYING=s S_VARIANCE=$(NOISY_EARLY_S_VARIANCE) Q=t EARLY_OCEAN=yes TODO=save-noisy-transfer
	$(MAKE) coho VARYING=s S_VARIANCE=$(NOISY_EARLY_S_VARIANCE) Q=r EARLY_OCEAN=yes TODO=save-noisy-transfer
	$(MAKE) coho VARYING=s S_VARIANCE=$(NOISY_EARLY_S_VARIANCE) Q=c EARLY_OCEAN=yes TODO=save-noisy-transfer
	$(MAKE) coho VARYING=s S_VARIANCE=$(NOISY_EARLY_S_VARIANCE) Q=d EARLY_OCEAN=yes TODO=save-noisy-transfer

small-s-figures:
	$(MAKE) transfer-figures mifi-figures mfr-figures S_MEAN=0.2828427 RECOMMENDED_S_VARIANCE=0.01 TRTAG=-small-s

very-small-s-figures:
	$(MAKE) transfer-figures mifi-figures mfr-figures S_MEAN=0.2 RECOMMENDED_S_VARIANCE=0.01 TRTAG=-very-small-s

no-sigma-figures:
	$(MAKE)	transfer-figures mifi-figures SIGMA=0 TRTAG=-no-sigma

save-transfer:
	$(MAKE) _save-transfer T=1024

_save-transfer: $(OUTDIR)/transfers.eps# $(OUTDIR)/transfers.eps.gv
	@mkdir -p figures
	cp $< figures/transfer-$(EXPT)-$(QNAME)-$(VARYNAME)$(TRTAG).eps

xsave-transfer: $(OUTDIR)/transfers.eps# $(OUTDIR)/transfers.eps.gv
	@mkdir -p figures
	cp $< figures/transfer-$(EXPT)-$(QNAME)-$(VARYNAME)$(XTAG)$(STAG).eps

save-noisy-transfer:
	$(MAKE) _save-noisy-transfer T=1024

_save-noisy-transfer: $(OUTDIR)/transfers.eps
	@mkdir -p figures
	cp $< figures/transfer-$(EXPT)-$(QNAME)-noisy-$(VARYNAME).eps

p01-impulse-figures:
	$(MAKE) chinook TODO="save-u save-pop save-q" VARYING=s Q=r VARIANCE=0.25 IMPULSE=yes T=16 TAG_MODE=impulse-s-
	$(MAKE) chinook TODO="save-q" VARYING=s Q=t VARIANCE=0.25 IMPULSE=yes T=16 TAG_MODE=impulse-s-
	$(MAKE) chinook TODO="save-u save-pop save-q" VARYING=s EARLY_OCEAN=yes Q=r VARIANCE=0.25 IMPULSE=yes T=16 TAG_MODE=impulse-early-s-
	$(MAKE) chinook TODO="save-q" VARYING=s EARLY_OCEAN=yes Q=t VARIANCE=0.25 IMPULSE=yes T=16 TAG_MODE=impulse-early-s-
	$(MAKE) chinook TODO="save-u save-pop save-q" VARYING=ac Q=r VARIANCE=0.25 IMPULSE=yes T=16 TAG_MODE=impulse-ac-
	$(MAKE) chinook TODO="save-q" VARYING=ac Q=t VARIANCE=0.25 IMPULSE=yes T=16 TAG_MODE=impulse-ac-

tables: eigenvector-tables #components-tables

xeigenvector-tables:
	@mkdir -p figures
	$(MAKE) chinook TODO=save-eigenvectors
	$(MAKE) coho TODO=save-eigenvectors
#	$(MAKE) chinook SIGMA=0 TODO=save-eigenvectors
#	$(MAKE) coho SIGMA=0 TODO=save-eigenvectors

eigenvector-tables:
	$(MAKE) grid GRID_TARGET=save-eigenvectors T=1

save-eigenvectors: $(OUTDIR)/eigenvectors.tex
	@mkdir -p figures
	cp $< figures/eigenvectors-$(EXPT)$(XTAG)$(STAG).tex

components-tables:
	@mkdir -p figures
	$(MAKE) chinook Q=r VARYING=ac TODO=save-components
	$(MAKE) chinook Q=t VARYING=ac TODO=save-components
	$(MAKE) chinook Q=r VARYING=s TODO=save-components
	$(MAKE) chinook Q=t VARYING=s TODO=save-components
	$(MAKE) coho Q=r VARYING=ac TODO=save-components
	$(MAKE) coho Q=t VARYING=ac TODO=save-components
	$(MAKE) coho Q=r VARYING=s TODO=save-components
	$(MAKE) coho Q=t VARYING=s TODO=save-components

save-components: $(OUTDIR)/variance-components.tex
	@mkdir -p figures
	cp $< figures/variance-components-$(EXPT)-$(QNAME)-$(VARYNAME).tex

########### targets for generating figures for paper 3 ############
## these also put them in figures/

p03-figures:
	$(MAKE) _p03-figures BW=yes

_p03-figures: p03-transfer-figures sweep-a-figures sweep-theta-figures \
 impulse-figures#  sweep-theta-cohort-figures cohort-sweep-sigma-figures \
# extinction-figures cohort-figures

xp03-transfer-figures:
	@mkdir -p figures
	$(MAKE) chinook VARYING=ac Q=t T=1024 TODO=save-transfer
	$(MAKE) chinook VARYING=s  Q=t T=1024 TODO=save-transfer
	$(MAKE) chinook VARYING=s  Q=t EARLY_OCEAN=yes T=1024 TODO=save-transfer
	$(MAKE) chinook VARYING=ac Q=r T=1024 TODO=save-transfer
	$(MAKE) chinook VARYING=s  Q=r T=1024 TODO=save-transfer
	$(MAKE) chinook VARYING=s  Q=r EARLY_OCEAN=yes T=1024 TODO=save-transfer

p03-transfer-figures:
	$(MAKE) grid GRID_TARGET=xsave-transfer VARYING=s Q=t T=1024
	$(MAKE) grid GRID_TARGET=xsave-transfer VARYING=s EARLY_OCEAN=yes Q=t T=1024
	$(MAKE) grid GRID_TARGET=xsave-transfer VARYING=ac Q=t T=1024
	$(MAKE) grid GRID_TARGET=xsave-transfer VARYING=s Q=r T=1024
	$(MAKE) grid GRID_TARGET=xsave-transfer VARYING=s EARLY_OCEAN=yes Q=r T=1024
	$(MAKE) grid GRID_TARGET=xsave-transfer VARYING=ac Q=r T=1024

FIGURES_N_TRIALS=1000

save-sweep-a-figure: sweep-a
	$(MAKE) save-extinction-times SWEEP_VAR=A FIGTAG=-real PLOT_EXT_LOG=yes

sweep-a-figures: sweep-a-s-figure sweep-a-early-s-figure sweep-a-ac-figure
#	$(MAKE) save-sweep-a-figure EXPT=chinook VARYING=s VARIANCE=0.1 EXTINCTION_THRESHOLD=50000 T=100000 TAG=-s-1d EXT_N_TRIALS=FIGURES_N_TRIALS
#	$(MAKE) save-sweep-a-figure EXPT=chinook VARYING=ac VARIANCE=0.5 EXTINCTION_THRESHOLD=500000 T=100000 TAG=-ac-1d EXT_N_TRIALS=FIGURES_N_TRIALS

sweep-a-s-figure:
	$(MAKE) save-sweep-a-figure EXPT=chinook VARYING=s EXTINCTION_THRESHOLD=900000 EXT_N_TRIALS=$(FIGURES_N_TRIALS)

sweep-a-early-s-figure:
	$(MAKE) save-sweep-a-figure EXPT=chinook VARYING=s EARLY_OCEAN=yes EXTINCTION_THRESHOLD=1000000 EXT_N_TRIALS=$(FIGURES_N_TRIALS)

sweep-a-ac-figure:
	$(MAKE) save-sweep-a-figure EXPT=chinook VARYING=ac EXTINCTION_THRESHOLD=1055000 EXT_N_TRIALS=$(FIGURES_N_TRIALS)

sweep-a-weak-figures:
	$(MAKE) save-sweep-a-figure EXPT=chinook VARYING=s VARIANCE=0.0001 EXTINCTION_THRESHOLD=1080000 EXT_N_TRIALS=$(FIGURES_N_TRIALS)
	$(MAKE) save-sweep-a-figure EXPT=chinook VARYING=ac VARIANCE=0.001 EXTINCTION_THRESHOLD=1085000 EXT_N_TRIALS=$(FIGURES_N_TRIALS)

save-sweep-theta-figure: sweep-theta
	$(MAKE) save-extinction-times SWEEP_VAR=THETA

sweep-theta-figures: sweep-theta-s-figures sweep-theta-early-s-figures\
 sweep-theta-ac-figures

sweep-theta-s-figures:
	$(MAKE) save-sweep-theta-figure EXPT=chinook VARYING=s EXTINCTION_THRESHOLD=1000000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex
	$(MAKE) save-sweep-theta-figure EXPT=chinook VARYING=s S_MEAN=0.2828427 RECOMMENDED_S_VARIANCE=0.01 EXTINCTION_THRESHOLD=100000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-small-s
	$(MAKE) save-sweep-theta-figure EXPT=coho VARYING=s EXTINCTION_THRESHOLD=750000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex
	$(MAKE) save-sweep-theta-figure EXPT=coho VARYING=s S_MEAN=0.2828427 RECOMMENDED_S_VARIANCE=0.01 EXTINCTION_THRESHOLD=300000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-small-s
	$(MAKE) save-sweep-theta-figure EXPT=coho VARYING=s S_MEAN=0.2 RECOMMENDED_S_VARIANCE=0.01 EXTINCTION_THRESHOLD=100000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-very-small-s

sweep-theta-early-s-figures:
	$(MAKE) save-sweep-theta-figure EXPT=chinook VARYING=s EARLY_OCEAN=yes EXTINCTION_THRESHOLD=1060000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex
	$(MAKE) save-sweep-theta-figure EXPT=chinook VARYING=s EARLY_OCEAN=yes S_MEAN=0.2828427 RECOMMENDED_S_VARIANCE=0.01 EXTINCTION_THRESHOLD=160000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-small-s
	$(MAKE) save-sweep-theta-figure EXPT=coho VARYING=s EARLY_OCEAN=yes EXTINCTION_THRESHOLD=830000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex
	$(MAKE) save-sweep-theta-figure EXPT=coho VARYING=s EARLY_OCEAN=yes S_MEAN=0.2828427 RECOMMENDED_S_VARIANCE=0.01 EXTINCTION_THRESHOLD=360000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-small-s
	$(MAKE) save-sweep-theta-figure EXPT=coho VARYING=s EARLY_OCEAN=yes S_MEAN=0.2 RECOMMENDED_S_VARIANCE=0.01 EXTINCTION_THRESHOLD=270000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-very-small-s

sweep-theta-ac-figures:
	$(MAKE) save-sweep-theta-figure EXPT=chinook VARYING=ac EXTINCTION_THRESHOLD=1060000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex
	$(MAKE) save-sweep-theta-figure EXPT=chinook VARYING=ac S_MEAN=0.2828427 RECOMMENDED_S_VARIANCE=0.01 EXTINCTION_THRESHOLD=160000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-small-s
	$(MAKE) save-sweep-theta-figure EXPT=coho VARYING=ac EXTINCTION_THRESHOLD=700000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex
	$(MAKE) save-sweep-theta-figure EXPT=coho VARYING=ac S_MEAN=0.2828427 RECOMMENDED_S_VARIANCE=0.01 EXTINCTION_THRESHOLD=380000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-small-s
	$(MAKE) save-sweep-theta-figure EXPT=coho VARYING=ac S_MEAN=0.2 RECOMMENDED_S_VARIANCE=0.01 EXTINCTION_THRESHOLD=270000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-very-small-s

sweep-theta-weak-figures:
	$(MAKE) save-sweep-theta-figure EXPT=chinook VARYING=s VARIANCE=0.00005 EXTINCTION_THRESHOLD=1085000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-weak
	$(MAKE) save-sweep-theta-figure EXPT=chinook VARYING=ac VARIANCE=0.001 EXTINCTION_THRESHOLD=1087000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-weak

sweep-theta-cohort-figures: sweep-theta-s-cohort-figures\
 sweep-theta-early-s-cohort-figures sweep-theta-ac-cohort-figures

sweep-theta-s-cohort-figures:
	$(MAKE) save-sweep-theta-figure EXPT=chinook VARYING=s COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=190000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-cohort
	$(MAKE) save-sweep-theta-figure EXPT=chinook VARYING=s S_MEAN=0.2828427 RECOMMENDED_S_VARIANCE=0.01 COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=1000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-cohort-small-s
	$(MAKE) save-sweep-theta-figure EXPT=coho VARYING=s COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=220000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-cohort
	$(MAKE) save-sweep-theta-figure EXPT=coho VARYING=s S_MEAN=0.2828427 RECOMMENDED_S_VARIANCE=0.01 COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=30000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-cohort-small-s
	$(MAKE) save-sweep-theta-figure EXPT=coho VARYING=s S_MEAN=0.2 RECOMMENDED_S_VARIANCE=0.01 COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=20000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-cohort-very-small-s

sweep-theta-early-s-cohort-figures:
	$(MAKE) save-sweep-theta-figure EXPT=chinook VARYING=s EARLY_OCEAN=yes COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=180000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-cohort
	$(MAKE) save-sweep-theta-figure EXPT=chinook VARYING=s EARLY_OCEAN=yes S_MEAN=0.2828427 RECOMMENDED_S_VARIANCE=0.01 COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=1500 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-cohort-small-s
	$(MAKE) save-sweep-theta-figure EXPT=coho VARYING=s EARLY_OCEAN=yes COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=220000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-cohort
	$(MAKE) save-sweep-theta-figure EXPT=coho VARYING=s EARLY_OCEAN=yes S_MEAN=0.2828427 RECOMMENDED_S_VARIANCE=0.01 COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=30000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-cohort-small-s
	$(MAKE) save-sweep-theta-figure EXPT=coho VARYING=s EARLY_OCEAN=yes S_MEAN=0.2 RECOMMENDED_S_VARIANCE=0.01 COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=8000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-cohort-very-small-s

sweep-theta-ac-cohort-figures:
	$(MAKE) save-sweep-theta-figure EXPT=chinook VARYING=ac COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=175000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-cohort
	$(MAKE) save-sweep-theta-figure EXPT=chinook VARYING=ac S_MEAN=0.2828427 RECOMMENDED_S_VARIANCE=0.01 COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=2500 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-cohort-small-s
	$(MAKE) save-sweep-theta-figure EXPT=coho VARYING=ac COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=220000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-cohort
	$(MAKE) save-sweep-theta-figure EXPT=coho VARYING=ac S_MEAN=0.2828427 RECOMMENDED_S_VARIANCE=0.01 COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=28000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-cohort-small-s
	$(MAKE) save-sweep-theta-figure EXPT=coho VARYING=ac S_MEAN=0.2 RECOMMENDED_S_VARIANCE=0.01 COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=8000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-cohort-very-small-s

xsweep-theta-cohort-no-sigma-figures: sweep-theta-s-cohort-no-sigma-figure\
 sweep-theta-early-s-cohort-no-sigma-figure \
 sweep-theta-ac-cohort-no-sigma-figure

sweep-theta-cohort-no-sigma-figures:
#	$(MAKE) grid GRID_TARGET=save-sweep-theta-figure VARYING=s SIGMA=0 COHORT_EXTINCTION=yes EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes XTAG=-complex-cohort-no-sigma
	$(MAKE) grid GRID_TARGET=save-sweep-theta-figure VARYING=s EARLY_OCEAN=yes SIGMA=0 COHORT_EXTINCTION=yes EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes XTAG=-complex-cohort-no-sigma
#	$(MAKE) grid GRID_TARGET=save-sweep-theta-figure VARYING=ac SIGMA=0 COHORT_EXTINCTION=yes EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes XTAG=-complex-cohort-no-sigma

# cohort extinction thresholds for no-sigma
ifeq ($(SIGMA),0)
ifeq ($(VARYING),s)
ifeq ($(EARLY_OCEAN),yes)
ifeq ($(EXPT),chinook)
ifeq ($(S_MEAN),0.85)
COHORT_EXTINCTION_THRESHOLD=170000
else #S_MEAN
COHORT_EXTINCTION_THRESHOLD=1000
endif #S_MEAN
else  #EXPT=coho
ifeq ($(S_MEAN),0.85)
COHORT_EXTINCTION_THRESHOLD=200000
else #S_MEAN
ifeq ($(S_MEAN),0.2828427)
COHORT_EXTINCTION_THRESHOLD=20000
else  #S_MEAN=0.2
COHORT_EXTINCTION_THRESHOLD=1000
endif #S_MEAN
endif #S_MEAN
endif #EXPT
else  #EARLY_OCEAN != yes
ifeq ($(EXPT),chinook)
ifeq ($(S_MEAN),0.85)
COHORT_EXTINCTION_THRESHOLD=210000
else #S_MEAN
COHORT_EXTINCTION_THRESHOLD=1000
endif #S_MEAN
else  #EXPT=coho
ifeq ($(S_MEAN),0.85)
COHORT_EXTINCTION_THRESHOLD=210000
else #S_MEAN
ifeq ($(S_MEAN),0.2828427)
COHORT_EXTINCTION_THRESHOLD=20000
else  #S_MEAN=0.2
COHORT_EXTINCTION_THRESHOLD=8000
endif #S_MEAN
endif #S_MEAN
endif #EXPT
endif #EARLY_OCEAN
else  #VARYING=ac
ifeq ($(EXPT),chinook)
ifeq ($(S_MEAN),0.85)
COHORT_EXTINCTION_THRESHOLD=188000
else #S_MEAN
COHORT_EXTINCTION_THRESHOLD=1000
endif #S_MEAN
else  #EXPT=coho
ifeq ($(S_MEAN),0.85)
COHORT_EXTINCTION_THRESHOLD=249117
else #S_MEAN
ifeq ($(S_MEAN),0.2828427)
COHORT_EXTINCTION_THRESHOLD=22000
else  #S_MEAN=0.2
COHORT_EXTINCTION_THRESHOLD=3000
endif #S_MEAN
endif #S_MEAN
endif #EXPT
endif #VARYING
endif #SIGMA

sweep-theta-s-cohort-no-sigma-figure:
	$(MAKE) save-sweep-theta-figure EXPT=chinook VARYING=s SIGMA=0 COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=188000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-cohort-no-sigma

sweep-theta-early-s-cohort-no-sigma-figure:
	$(MAKE) save-sweep-theta-figure EXPT=chinook VARYING=s EARLY_OCEAN=yes SIGMA=0 COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=170000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-cohort-no-sigma

sweep-theta-ac-cohort-no-sigma-figure:
	$(MAKE) save-sweep-theta-figure EXPT=chinook VARYING=ac SIGMA=0 COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=210000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-complex-cohort-no-sigma

_sweep-theta-ac-cohort-no-sigma-extra-figures:
	$(MAKE) save-sweep-theta-figure EXPT=chinook VARYING=ac SIGMA=0 COHORT_EXTINCTION=yes EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=yes FIGTAG=-cthresh_$(COHORT_EXTINCTION_THRESHOLD)-complex-cohort-no-sigma

sweep-theta-ac-cohort-no-sigma-extra-figures:
#	17000 is too low
#	$(MAKE) _sweep-theta-ac-cohort-no-sigma-extra-figures COHORT_EXTINCTION_THRESHOLD=180000
##	$(MAKE) _sweep-theta-ac-cohort-no-sigma-extra-figures COHORT_EXTINCTION_THRESHOLD=210500
	$(MAKE) _sweep-theta-ac-cohort-no-sigma-extra-figures COHORT_EXTINCTION_THRESHOLD=210700

sweep-theta-small-s-figures:
	$(MAKE) sweep-theta-figures sweep-theta-cohort-figures RECOMMENDED_S_VARIANCE=0.2828427 XTAG=-small-s
	$(MAKE) sweep-theta-figures sweep-theta-cohort-figures RECOMMENDED_S_VARIANCE=0.2 XTAG=-very-small-s

save-sweep-sigma-figure: extinction-sweep-sigma
	@mkdir -p figures
	$(MAKE) save-extinction-times SWEEP_VAR=SIGMA PLOT_ANALYTIC_EXTINCTION_TIMES=no

extinction-sweep-sigma:
	$(MAKE) ext-sweep SWEEP_VAR=SIGMA SWEEP_FROM=0 SWEEP_TO=$(SIGMA) SWEEP_STEPS=15

cohort-sweep-sigma-figures:
	$(MAKE) save-sweep-sigma-figure EXPT=chinook VARYING=s COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=150000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=no FIGTAG=-cohort-sweep-sigma
	$(MAKE) save-sweep-sigma-figure EXPT=chinook VARYING=s EARLY_OCEAN=yes COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=168000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=no FIGTAG=-cohort-sweep-sigma
	$(MAKE) save-sweep-sigma-figure EXPT=chinook VARYING=ac COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=179000 EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=no FIGTAG=-cohort-sweep-sigma

_extra-ac-sweep-sigma-figures:
	$(MAKE) save-sweep-sigma-figure EXPT=chinook VARYING=ac COHORT_EXTINCTION=yes EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=no FIGTAG=-cohort-sweep-sigma-cthresh_$(COHORT_EXTINCTION_THRESHOLD)
	$(MAKE) save-extinction-path SIGMA=0 EXPT=chinook VARYING=ac COHORT_EXTINCTION=yes EXT_N_TRIALS=$(FIGURES_N_TRIALS) PLOT_EXT_LOG=no PLOT_BY_ANGLE=no XTAG=-cohort-sweep-sigma-cthresh_$(COHORT_EXTINCTION_THRESHOLD)

extra-ac-sweep-sigma-figures:
	$(MAKE) _extra-ac-sweep-sigma-figures COHORT_EXTINCTION_THRESHOLD=210867
	$(MAKE) _extra-ac-sweep-sigma-figures COHORT_EXTINCTION_THRESHOLD=210868

extinction-path-figures: #compare-ac-extinction-path-figures
	$(MAKE) save-extinction-path EXPT=chinook VARYING=s EXTINCTION=yes COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=190000 COLOR=yes COMPLEX_COLOR=yes R=0.95 THETA=1.570795
	$(MAKE) save-extinction-path EXPT=chinook VARYING=s EARLY_OCEAN=yes EXTINCTION=yes COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=180000 COLOR=yes COMPLEX_COLOR=yes R=0.95 THETA=1.570795

compare-ac-extinction-path-figures:
	$(MAKE) save-extinction-path save-u EXPT=chinook VARYING=ac EXTINCTION=yes COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=179000 SIGMA=0 XTAG=-chinook-$(VARYNAME)-sigma_0
	$(MAKE) save-extinction-path save-u EXPT=chinook VARYING=ac EXTINCTION=yes COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=179000 SIGMA=0.133333333333333 XTAG=-chinook-$(VARYNAME)-sigma_0.133333333333333
	$(MAKE) save-extinction-path save-u EXPT=chinook VARYING=ac EXTINCTION=yes COHORT_EXTINCTION=yes COHORT_EXTINCTION_THRESHOLD=179000 SIGMA=0.4 XTAG=-chinook-$(VARYNAME)-sigma_0.4

save-extinction-path: $(OUTDIR)/population.tail.eps
	@mkdir -p figures
	cp $< figures/extinction-path$(XTAG).eps

extinction-figures: chinook-extinction-figures coho-extinction-figures

chinook-extinction-figures:
	@mkdir -p figures
	$(MAKE) chinook VARYING=ac TODO=sweep-variance-small EXT_SWEEP_AFTER=save-extinction-times PLOT_EXT_LOG=yes EXT_N_TRIALS=$(FIGURES_N_TRIALS) EXTINCTION_THRESHOLD=1080000 T=16384 SM_VARIANCE_BEGIN=0.00001 SM_VARIANCE_END=0.00299 SM_VARIANCE_STEPS=20
	$(MAKE) chinook VARYING=s TODO=ext-sweep EXT_SWEEP_AFTER= EXT_N_TRIALS=$(FIGURES_N_TRIALS) EXTINCTION_THRESHOLD=1080000 T=16384 SWEEP_VAR=VARIANCE SWEEP_FROM=0.000002 SWEEP_TO=0.000175 SWEEP_STEPS=6
	$(MAKE) chinook VARYING=s TODO=ext-sweep EXT_SWEEP_AFTER=save-extinction-times PLOT_EXT_LOG=yes EXT_N_TRIALS=$(FIGURES_N_TRIALS) EXTINCTION_THRESHOLD=1080000 T=16384 SWEEP_VAR=VARIANCE SWEEP_FROM=0.0002 SWEEP_TO=0.00299 SWEEP_STEPS=19 SWEEP_TARGETS=sweep-run

coho-extinction-figures:
	@mkdir -p figures
	$(MAKE) coho VARYING=ac TODO=sweep-variance-small EXT_SWEEP_AFTER=save-extinction-times PLOT_EXT_LOG=yes EXT_N_TRIALS=$(FIGURES_N_TRIALS) EXTINCTION_THRESHOLD=880000 T=16384 SM_VARIANCE_END=0.00099 SM_VARIANCE_STEPS=20
	$(MAKE) coho VARYING=s TODO=ext-sweep EXT_SWEEP_AFTER= EXT_N_TRIALS=$(FIGURES_N_TRIALS) EXTINCTION_THRESHOLD=880000 T=16384 SWEEP_VAR=VARIANCE SWEEP_FROM=0.000002 SWEEP_TO=0.00006 SWEEP_STEPS=6
	$(MAKE) coho VARYING=s TODO=ext-sweep EXT_SWEEP_AFTER=save-extinction-times PLOT_EXT_LOG=yes EXT_N_TRIALS=$(FIGURES_N_TRIALS) EXTINCTION_THRESHOLD=880000 T=16384 SWEEP_VAR=VARIANCE SWEEP_FROM=0.000075 SWEEP_TO=0.00099 SWEEP_STEPS=19 SWEEP_TARGETS=sweep-run

cohort-figures: steady-mode-figures extinction-mode-figures impulse-figures

impulse-figures:
	$(MAKE) chinook TODO="save-mode-figures" VARYING=s VARIANCE=0.25 IMPULSE=yes T=32 TAG_MODE=impulse-s-
	$(MAKE) chinook TODO="save-mode-figures" VARYING=s EARLY_OCEAN=yes VARIANCE=0.25 IMPULSE=yes T=32 TAG_MODE=impulse-early-s-
	$(MAKE) chinook TODO="save-mode-figures" VARYING=ac VARIANCE=0.25 IMPULSE=yes T=32 TAG_MODE=impulse-ac-

steady-mode-figures:
#	$(MAKE) chinook TODO="save-mode-figures" VARYING=s VARIANCE=0.1 TAG_MODE=steady-s-
	$(MAKE) chinook TODO="save-mode-figures save-eigenvectors" VARYING=s VARIANCE=0.05 T=1024 TAG_MODE=steady-s-
	$(MAKE) chinook TODO="save-mode-figures" VARYING=ac VARIANCE=0.5 T=1024 TAG_MODE=steady-ac-

extinction-mode-figures:
	$(MAKE) chinook TODO=save-mode-figures EXTINCTION=yes VARYING=s VARIANCE=0.1 EXTINCTION_THRESHOLD=50000 T=100000 TAG_MODE=extinction-s-
	$(MAKE) chinook TODO=save-mode-figures EXTINCTION=yes VARYING=ac VARIANCE=0.5 EXTINCTION_THRESHOLD=500000 T=100000 TAG_MODE=extinction-ac-

save-mode-figures: save-pop save-pop-cohorts $(addprefix save-,$(MODE_NAMES)) \
 $(addprefix save-,$(addsuffix -cohorts,$(MODE_NAMES))) save-u

save-pop : $(OUTDIR)/population.tail.eps
	@mkdir -p figures
	cp $< figures/$(TAG_MODE)$(<F)

plot-cohorts : $(OUTDIR)/population.tail.cum.cohorts.eps

save-pop-cohorts : $(OUTDIR)/population.tail.cum.cohorts.eps
	@mkdir -p figures
	cp $< figures/$(TAG_MODE)$(subst .cum,,$(<F))

save-mode% : $(OUTDIR)/mode%.tail.eps
	@mkdir -p figures
	cp $< figures/$(TAG_MODE)$(<F)

save-mode%-cohorts : $(OUTDIR)/mode%.tail.cum.cohorts.eps
	@mkdir -p figures
	cp $< figures/$(TAG_MODE)$(subst .cum,,$(<F))

save-u : $(OUTDIR)/u.tail.eps
	@mkdir -p figures
	cp $< figures/$(TAG_MODE)u$(XTAG).tail.eps

save-q : $(OUTDIR)/output-quantity.tail.eps
	@mkdir -p figures
	cp $< figures/$(TAG_MODE)$(QNAME).tail.eps

############ section of targets used for 'sweep' ##############

COLLECT_TARGETS = forcing.out transformed-forcing.mag.out\
 transformed-weights.mag.out eigenvalues.out eigenvalues-mag.out deltas.out\
 analytic-variance-gain.out observed-variance-gain.out resonances.mag.out 
# analytic-variance-contributions.out
ifneq ($(EXTINCTION),yes)
COLLECT_TARGETS += transfer-peaks.out
endif
COLLECT_PATHS = $(addprefix $(SWEEP_OUTDIR)/, $(COLLECT_TARGETS))

SWEEP_OUTDIR = out/sweep-$(SWEEP_VAR)$(TAG)

# exceptional collect rules

collect-eigenvalues.out :
	cat $(OUTDIR)/eigenvalues.out >> $(SWEEP_OUTDIR)/eigenvalues.out

# these duplicate commands for print-%-variance-gain
collect-observed-variance-gain.out :
	(echo -n $($(SWEEP_VAR)) ""; perl -e "print `./scripts/variance $(OUTDIR)/output-quantity.out` / $(VARIANCE), \"\\n\"") >> $(SWEEP_OUTDIR)/observed-variance-gain.out

collect-analytic-variance-gain.out :
	(echo -n $($(SWEEP_VAR)) ""; perl -e "print '`cat $(OUTDIR)/analytic-variance.out`' / $(VARIANCE), \"\\n\"") >> $(SWEEP_OUTDIR)/analytic-variance-gain.out

# how to collect in general
collect-%.out :
	(echo -n $($(SWEEP_VAR)) ""; cat $(OUTDIR)/$*.out) >> $(SWEEP_OUTDIR)/$*.out

# new version, uses implicit rules for each collect target
_collect: $(addprefix collect-, $(COLLECT_TARGETS)) collect-eigenvectors

COLLECT_EIGENVECTORS_PREFIX = collect-
COLLECT_EIGENVECTORS_SUFFIX = out
collect-eigenvectors:
	for((i=0;i<$(N);++i)); do \
	 $(MAKE) \
	   $(COLLECT_EIGENVECTORS_PREFIX)u$$i.mag.$(COLLECT_EIGENVECTORS_SUFFIX)  \
	   $(COLLECT_EIGENVECTORS_PREFIX)v$$i.mag.$(COLLECT_EIGENVECTORS_SUFFIX); \
	done;

# new version
plot-collect: $(COLLECT_PATHS:.out=.eps) 

SHOW_CIRCLE=yes
ifeq ($(SHOW_CIRCLE),yes)
CIRCLE_ARG=-c
endif
$(SWEEP_OUTDIR)/eigenvalues.gp: $(SWEEP_OUTDIR)/eigenvalues.out ./scripts/plot-eigenvalues
	./scripts/plot-eigenvalues $< $(CIRCLE_ARG) $(BW_ARG) $(OC_ARG) $(CLOSED_ARG)

# new version
view-collect: $(COLLECT_PATHS:.out=.eps.gv)

# this should only be called recursively
# requires SWEEP_VAR, SWEEP_FROM, SWEEP_TO, SWEEP_STEPS to be set
SWEEP_TARGETS = sweep-prep sweep-run #sweep-collect
sweep: $(SWEEP_TARGETS)

sweep-prep:
	mkdir -p $(SWEEP_OUTDIR)
	$(RM) $(SWEEP_OUTDIR)/*
#	$(RM) $(COMBINED_FORCING) $(COMBINED_TRANSFORMED_FORCING) $(COMBINED_EIGENVALUES) $(COMBINED_EIGENVALUES_MAG) $(COMBINED_DELTAS) $(COMBINED_TRANSFER_PEAKS) $(SWEEP_OUTDIR)/u*.out $(SWEEP_OUTDIR)/v*.out

_sweep:
	for((i=0;i<=$(SWEEP_STEPS);++i)); \
	do $(MAKE) LD_LIBRARY_PATH=$(LD_LIBRARY_PATH) \
	  $(SWEEP_VAR)=`perl -e "print $(SWEEP_FROM)+($(SWEEP_TO) - $(SWEEP_FROM))*$$i/$(SWEEP_STEPS)"` \
	  $(_SWEEP_TARGET);\
	done;

SWEEP_RUN_TARGET = run
sweep-run:
	$(MAKE) _SWEEP_TARGET="clear-outdir $(SWEEP_RUN_TARGET) _collect" _sweep

#sweep-run: sweep-collect

sweep-collect:
	$(MAKE) _SWEEP_TARGET=_collect _sweep

SIGMA_BEGIN=0
SIGMA_END=$(SIGMA)#0.975
SIGMA_STEPS = 39
#INCREMENT = `expr ( $(SIGMA_END) - $(SIGMA_BEGIN) ) / $(STEPS)`

SWEEP_TODO = sweep plot-collect #view-collect
# these are to be called directly
sweep-sigma:
	$(MAKE) SWEEP_VAR=SIGMA SWEEP_FROM=$(SIGMA_BEGIN) SWEEP_TO=$(SIGMA_END) SWEEP_STEPS=$(SIGMA_STEPS) $(SWEEP_TODO)

AC_BEGIN=3
AC_END=4
AC_STEPS = 40

sweep-ac:
	$(MAKE) SWEEP_VAR=AC_MEAN SWEEP_FROM=$(AC_BEGIN) SWEEP_TO=$(AC_END) SWEEP_STEPS=$(AC_STEPS) $(SWEEP_TODO)

################### extinction time experiments ####################

EXT_N_TRIALS = 1 #200
extinction-trials:
	for((i=0;i<$(EXT_N_TRIALS);++i)); \
	do $(MAKE) --no-print-directory LD_LIBRARY_PATH=$(LD_LIBRARY_PATH) \
	  EXTINCTION_THRESHOLD=$(EXTINCTION_THRESHOLD) \
	  run collect-observed-extinction-time; \
	done;

# this happens after each trial
EXT_COLLECT_FILE = $(SWEEP_OUTDIR)/observed-extinction-times.out
collect-observed-extinction-time:
	(echo -n $($(SWEEP_VAR)) ""; cat $(OUTDIR)/extinction-time.out) >> $(EXT_COLLECT_FILE)

# this happens after each set of trials
# note analytic-extinction-time might not be there, not an error
A_EXT_COLLECT_FILE = $(SWEEP_OUTDIR)/analytic-extinction-times.out
collect-analytic-extinction-time:
	if [ -e $(OUTDIR)/analytic-extinction-time.out ]; then \
	  (echo -n $($(SWEEP_VAR)) ""; cat $(OUTDIR)/analytic-extinction-time.out) >> $(A_EXT_COLLECT_FILE); else \
	  touch $(A_EXT_COLLECT_FILE); fi
	cat $(OUTDIR)/extinction-objs.txt >> $(SWEEP_OUTDIR)/extinction-objs.out

EXT_SWEEP_AFTER = plot-extinction-times
ext-sweep:
	$(MAKE) EXTINCTION=yes SWEEP_RUN_TARGET="extinction-trials collect-analytic-extinction-time" sweep $(EXT_SWEEP_AFTER)

#EXT_THRESH_BEGIN=100000
#EXT_THRESH_END=1000000 #500000
EXT_THRESH_BEGIN=800000
EXT_THRESH_END=1089000 #500000
EXT_THRESH_STEPS=20
sweep-extinction-threshold:
	$(MAKE) SWEEP_VAR=EXTINCTION_THRESHOLD SWEEP_FROM=$(EXT_THRESH_BEGIN) SWEEP_TO=$(EXT_THRESH_END) SWEEP_STEPS=$(EXT_THRESH_STEPS) ext-sweep

VARIANCE_EXP_BEGIN=-10
VARIANCE_EXP_END=-1
VARIANCE_EXP_STEPS=9
sweep-variance-exp:
	$(MAKE) SWEEP_VAR=VARIANCE_EXP SWEEP_FROM=$(VARIANCE_EXP_BEGIN) SWEEP_TO=$(VARIANCE_EXP_END) SWEEP_STEPS=$(VARIANCE_EXP_STEPS) ext-sweep

VARIANCE_BEGIN=0
VARIANCE_END=0.1
VARIANCE_STEPS=20
sweep-variance:
	$(MAKE) SWEEP_VAR=VARIANCE SWEEP_FROM=$(VARIANCE_BEGIN) SWEEP_TO=$(VARIANCE_END) SWEEP_STEPS=$(VARIANCE_STEPS) ext-sweep EXTINCTION_THRESHOLD=1000000

SM_VARIANCE_BEGIN=0
SM_VARIANCE_END=0.003
SM_VARIANCE_STEPS=60
sweep-variance-small:
	$(MAKE) SWEEP_VAR=VARIANCE SWEEP_FROM=$(SM_VARIANCE_BEGIN) SWEEP_TO=$(SM_VARIANCE_END) SWEEP_STEPS=$(SM_VARIANCE_STEPS) ext-sweep

A_BEGIN=-0.95
A_END=0.95
A_STEPS=19#38
sweep-a:
	$(MAKE) SWEEP_VAR=A SWEEP_FROM=$(A_BEGIN) SWEEP_TO=$(A_END) SWEEP_STEPS=$(A_STEPS) ext-sweep COLOR=yes

TH_BEGIN=0
#TH_END=6.28318 # 2pi
TH_END=3.14159 # pi
TH_STEPS=16
TH_FIXED_R=0.95
sweep-theta:
	$(MAKE) SWEEP_VAR=THETA SWEEP_FROM=$(TH_BEGIN) SWEEP_TO=$(TH_END) SWEEP_STEPS=$(TH_STEPS) ext-sweep COLOR=yes COMPLEX_COLOR=yes R=$(TH_FIXED_R)

ext_collect:
	$(MAKE) _SWEEP_TARGET="_collect" _sweep

plot-extinction-times: $(SWEEP_OUTDIR)/extinction-times.eps

save-extinction-times: $(SWEEP_OUTDIR)/extinction-times.eps
	@mkdir -p figures
	cp $< figures/extinction-times-$(EXPT)-$(VARYNAME)$(XTAG)$(STAG).eps
	$(MAKE) figures/extinction-times-$(EXPT)-$(VARYNAME)$(XTAG)$(STAG).eps.gv

#save-extinction-times: $(SWEEP_OUTDIR)/extinction-times.eps
#	@mkdir -p figures
#	cp $< figures/extinction-times-$(EXPT)$(FIGTAG)-$(VARYNAME).eps
#	$(MAKE) figures/extinction-times-$(EXPT)$(FIGTAG)-$(VARYNAME)$(XTAG).eps.gv

ifeq ($(PLOT_EXT_LOG),yes)
PE_ARG=-log
endif
ifeq ($(PLOT_BY_ANGLE),yes)
PE_ARG += -angle
endif
NTRIALS_IN_TITLE=yes
ifeq ($(NTRIALS_IN_TITLE),yes)
PE_ARG += -title "$(EXT_N_TRIALS) trials"
endif
ifeq ($(PLOT_ANALYTIC_EXTINCTION_TIME),yes)
%/extinction-times.gp : %/observed-extinction-times.avgstd.out %/analytic-extinction-times.out
	./scripts/plot-extinction-times $^ -g $@ $(PE_ARG)
else
%/extinction-times.gp : %/observed-extinction-times.avgstd.out
	./scripts/plot-extinction-times $^ -g $@ $(PE_ARG)
endif

%/combined-extinction-times.out : %/observed-extinction-times.out %/analytic-extinction-times.out
	paste $^ | sed 's/	/ /g' | cut -d ' ' -f 1,2,4,6 > $@

%/extinction-times-transformed.gp : %/combined-extinction-times.out
	./scripts/plot-extinction-times -t $< $@

################### transform files into other forms of data #############

%.touch :
	touch $*

%.remake :
	$(RM) $*
	$(MAKE) $*

%.more :
	@more $*

more-deltas: $(OUTDIR)/deltas.out $(OUTDIR)/deltas.out.more

%.rm : 
	$(RM) $*

%/population.cum.out : %/population.out
	$(RM) $@
	./scripts/cumulative $< 2 - >$@

%/population.tail.cum.out : %/population.tail.out
	$(RM) $@
	./scripts/cumulative $< 2 - >$@

%.sum.out : %.out
	$(RM) $@
	./scripts/addcols $< 2 - >$@

%.cum.out : %.out
	$(RM) $@
	./scripts/cumulative $< 1 - >$@

%.acv.out : %.out ./scripts/autocov # autocovariance
	./scripts/autocov $<

%.fft.out : %.out ./scripts/fft
	$(RM) $@
	./fft <$< >$@

%.ac.out : %.out ./scripts/ac # remove dc component
	./scripts/ac $<

%.mag.out : %.out ./scripts/mag
	./scripts/mag $<

%.phase.out : %.out ./scripts/phase
	./scripts/phase $<

%.magphase.out : %.out ./scripts/magphase
	./scripts/magphase $<

%.delta-transfer.out : %.x.ac.fft.mag.out %.delta.fft.mag.out # %.delta.ac.fft.mag.out
	./scripts/pwdiv $^ >$@

%.s-transfer.out : %.x.ac.fft.mag.out %.s.fft.mag.out
	./scripts/pwdiv $^ >$@

#%.transfer.out : %.x.ac.fft.mag.out %.eps.fft.mag.out
%.transfer.out : %.x.fft.mag.out %.eps.fft.mag.out
	./scripts/pwdiv $^ >$@

%.etimes.out:# %,thresh=-1.ext.etime.out
	cat $(@:.etimes.out=*.etime.out) >$@

#%.etimes.out: #%.ext.etime.out #fake target but makes it rerun the exe #: %*.etime.out

%.avg.out : %.out
	./scripts/averages $^ >$@

%.logavg.out : %.out
	./scripts/averages -log $^ >$@

%.avgstd.out : %.out
	./scripts/averages -d $^ >$@

BINSIZE=0.05
%.bins.out : %.out
	./scripts/bins $(BINSIZE) <$^ >$@

%.last1000.out : %.out
	tail -1000 $< >$@

%.last100.out : %.out
	tail -100 $< >$@

TAIL_LENGTH = 32
%.tail.out : %.out
	tail -$(TAIL_LENGTH) $< >$@

#################### code for how to plot things ######################

%.eps : %.gp %.out
	gnuplot $<

%.eps : %.m
	octave $<

%.eps : %.gp
	gnuplot $<

.PRECIOUS: %.eps
%.eps.gv : %.eps
	$(GV) $< &

# for convenience
plot-tmf: $(OUTDIR)/transfers.eps.gv $(OUTDIR)/mifi.mag.eps.gv $(OUTDIR)/mfr.mag.eps.gv

plot-transfers: $(OUTDIR)/transfers.eps $(OUTDIR)/transfers.eps.gv

plot-fft: $(OUTDIR)/population-fft.eps.gv

%/population.grey.m : %/population.out
	./scripts/plot-population-octave $< $@

%/population.tail.grey.m : %/population.tail.out
	./scripts/plot-population-octave $< $@

CLIP_FACTOR=1.6
ifeq ($(BW),yes)
TRANSFER_DEST=transfers.unfixed.gp
$(OUTDIR)/transfers.eps : $(OUTDIR)/transfers.unfixed.eps
	scripts/fix-transfer-plots $< >$@
else
TRANSFER_DEST=transfers.gp
endif
$(OUTDIR)/$(TRANSFER_DEST): $(OUTDIR)/observed-transfer.out $(OUTDIR)/analytic-transfer.out $(OUTDIR)/transfer-peaks.out ./scripts/plot-transfers
#	./scripts/plot-transfers $(OUTDIR)/observed-transfer.out $(OUTDIR)/analytic-transfer.out -d -tp $(OUTDIR)/transfer-peaks.out -g $@ -half $(BW_ARG)
	./scripts/plot-transfers $(OUTDIR)/observed-transfer.out $(OUTDIR)/analytic-transfer.out -d -tp $(OUTDIR)/transfer-peaks.out -g $@ -half $(BW_ARG) -clip $(CLIP_FACTOR) $(CLOSED_ARG)
#	./scripts/plot-transfers $(OUTDIR)/observed-transfer.out $(OUTDIR)/analytic-transfer.out -d -g $@ -half $(BW_ARG)
#else
#$(OUTDIR)/transfers.gp: $(OUTDIR)/observed-transfer.out $(OUTDIR)/analytic-transfer.out
#	./scripts/plot-transfers $(OUTDIR)/observed-transfer.out $(OUTDIR)/analytic-transfer.out -d -tp $(OUTDIR)/transfer-peaks.out -g $@ -half $(BW_ARG)
#	./scripts/plot-transfers-with-peaks $(OUTDIR)/observed-transfer.out $(OUTDIR)/analytic-transfer.out $(OUTDIR)/transfer-peaks.out $(OUTDIR)/transfers.gp N=$(N),S_MEAN=$(S_MEAN),S_VAR=$(S_VARIANCE),AC_MEAN=$(AC_MEAN),AC_VAR=$(AC_VARIANCE),SIGMA=$(SIGMA),Q=$(Q),$(VARYING)
#endif

$(OUTDIR)/population-fft.gp : $(OUTDIR)/population.fft.mag.out $(OUTDIR)/eigenvalues.out ./scripts/plot-transfers
	scripts/plot-transfers $(OUTDIR)/population.fft.mag.out -d -x $(OUTDIR)/eigenvalues.out -g $@ $(BW_ARG) -half $(CLOSED_ARG)
#	scripts/plot-fft $< -g $@

#xplot-transfer-components: $(OUTDIR)/analytic-transfer-components.eps
#	$(GV) $(OUTDIR)/transfer-components.eps &

#%/transfer-components.gp: %/analytic-transfer-components.out
#	./scripts/plot-transfer-components $^ $(<:transfer-components=analytic-transfer) $< $(<:transfer-components=transfer-peaks. N=$(N),S_MEAN=$(S_MEAN),S_VAR=$(S_VARIANCE),AC_MEAN=$(AC_MEAN),AC_VAR=$(AC_VARIANCE),SIGMA=$(SIGMA),Q=$(Q),$(VARYING)

$(OUTDIR)/transfer-components.gp: $(OUTDIR)/analytic-transfer.out $(OUTDIR)/analytic-transfer-components.out $(OUTDIR)/transfer-peaks.out
	./scripts/plot-transfer-components $^ $@ N=$(N),S_MEAN=$(S_MEAN),S_VAR=$(S_VARIANCE),AC_MEAN=$(AC_MEAN),AC_VAR=$(AC_VARIANCE),SIGMA=$(SIGMA),Q=$(Q),$(VARYING)

plot-transfer-components: $(OUTDIR)/transfer-components.eps $(OUTDIR)/transfer-components.eps.gv

PLOT_CLOSED=yes
ifeq ($(PLOT_CLOSED),yes)
CLOSED_ARG = -closed
endif
$(OUTDIR)/eigenvalues.gp: $(OUTDIR)/eigenvalues.out ./scripts/plot-eigenvalues
	./scripts/plot-eigenvalues $< -c $(BW_ARG) $(OC_ARG) $(CLOSED_ARG)

plot-eigenvalues: $(OUTDIR)/eigenvalues.eps $(OUTDIR)/eigenvalues.eps.gv

COMPONENT_PLOTS = $(OUTDIR)/transformed-forcing.eps $(OUTDIR)/transformed-weights.eps $(OUTDIR)/mifi.eps
plot-components: $(COMPONENT_PLOTS) $(COMPONENT_PLOTS:eps=eps.gv)

%/transformed-forcing.gp: %/transformed-forcing.out ./scripts/plot-eigenvalues
	./scripts/plot-eigenvalues $< $(BW_ARG)
#	./scripts/plot-eigenvalues $< "VH[i]" $(BW_ARG)
#	gnuplot $(OUTDIR)/transformed-forcing.gp
#	$(GV) $(OUTDIR)/transformed-forcing.eps &

%/transformed-weights.gp: %/transformed-weights.out ./scripts/plot-eigenvalues
	./scripts/plot-eigenvalues $< $(BW_ARG)
#	./scripts/plot-eigenvalues $< "qU[i]" $(BW_ARG)
#	gnuplot $(OUTDIR)/transformed-weights.gp
#	$(GV) $(OUTDIR)/transformed-weights.eps &
#	./scripts/plot-eigenvalues $(OUTDIR)/mifi.out "qU[i] VH[i]"
#	gnuplot $(OUTDIR)/mifi.gp
#	$(GV) $(OUTDIR)/mifi.eps &

%/mifi.gp: %/mifi.out ./scripts/plot-eigenvalues
	./scripts/plot-eigenvalues $< $(BW_ARG)
#	./scripts/plot-eigenvalues $< "qU[i] VH[i]" $(BW_ARG)

# pattern rule appears to require a command even if null
view-%: $(OUTDIR)/%.eps $(OUTDIR)/%.eps.gv
	@#

## some obsolete rules, some not

%/combined-analytic-variance-gain.gp : %/combined-analytic-variance-gain.out
	./scripts/plot-cols -x 1 2 -w l $< $(BW_ARG)

%/combined-observed-variance-gain.gp : %/combined-observed-variance-gain.out
	./scripts/plot-cols -x 1 2 -w l $< $(BW_ARG)

#%/population.gp : %/population.out
#	./scripts/plot-cols -x 1 all -w l $< $(BW_ARG) -nokey

%/population.cum.gp : %/population.cum.out
	./scripts/plot-cols -x 1 all -w l $< $(BW_ARG)

#%/population.tail.gp : %/population.tail.out
#	./scripts/plot-cols -x 1 all -w l $< $(BW_ARG) -nokey

%/population.tail.cum.gp : %/population.tail.cum.out
	./scripts/plot-cols -x 1 all -w l $< $(BW_ARG)

mode%.tail.gp : mode%.tail.out
	./scripts/plot-cols -x 1 all -w l $< $(BW_ARG)

%.cum.gp : %.cum.out
	./scripts/plot-cols -x 1 all -w l $< $(BW_ARG)

%/u.gp : %/u.out
	./scripts/plot-cols -x 1 all -w l $< $(BW_ARG) -nokey

%/u.tail.gp : %/u.tail.out
	./scripts/plot-cols -x 1 all -w l $< $(BW_ARG) -nokey

%/eigenvalues.gp : %/eigenvalues.out ./scripts/plot-eigenvalues
	./scripts/plot-eigenvalues $< $(BW_ARG) $(OC_ARG) $(CLOSED_ARG)

%/deltas.gp : %/deltas.out
	./scripts/plot-cols -x 1 2 3 -w l $< $(BW_ARG)

%/extinction-times.gp : %/extinction-times.out
	./scripts/plot-cols -d -x 1 all -w p $< $(BW_ARG)

%/analytic-variance-components.gp : %/analytic-variance-components.out
	./scripts/plot-cols -x 1 all -w l $< $(BW_ARG)

%.cohorts.gp : %.out
	./scripts/plot-cohorts -c $(AC_MEAN) -x 1 all $< -g $@

%.avg.gp : %.avg.out
	./scripts/plot-2cols 1 2 $<

%.avgstd.gp : %.avgstd.out
	./scripts/plot-avgstd $<

#%.x0.gp : %.x.out plot-x
%.x.gp : %.x.out plot-1stcol-1st200
	./scripts/plot-1stcol-1st200 $<

%.x1x2.gp : %.x.out plot-cov
	./scripts/plot-cov 1 2 $<

%.fft.gp : %.fft.out plot-1stcol
	./scripts/plot-1stcol $< $(BW_ARG) -half

#%.acv.gp : %.acv.out ./plot-acv
#	./plot-acv $<
%.acv.gp : %.acv.out ./plot-1stcol
	./scripts/plot-1stcol $< $(BW_ARG)

%.fft.mag.gp : %.fft.mag.out ./scripts/plot-transfers
	./scripts/plot-transfers $< -g $@ $(BW_ARG) -d

%.mag.gp : %.mag.out
	./scripts/plot-mags $< -bw

#%.transfers.gp : %.s-transfer.out %.an-transfer.mag.out
#%.transfers.gp : %.delta-transfer.out %.an-transfer.mag.out
#%.transfers.gp : %.transfer.out %.an-transfer.mag.out
#	./scripts/plot-transfers $^ -g $@

%.transfers.gp : ./plot-transfers

%.transfer.gp : %.transfer.out ./plot-transfer
	./scripts/plot-transfer $<

%.ratio.eps : %.ratio.out
	./scripts/plot-ratio $<

%.y.gp : %.y.out plot-y
	./scripts/plot-y $<

%.etimes.gp : %.etimes.out %.etimes.avg.out
	./scripts/plot-phase 1 2 $< -lines 1 2 $(<:out=avg.out) \
		$@ $(@:gp=eps)
#	./plot-2cols 1 2 $<

%.etimes.log.gp : %.etimes.out %.etimes.logavg.out
	./scripts/plot-phase -loglin 1 2 $< -lines 1 2 $(<:out=logavg.out) \
		$@ $(@:gp=eps)

## not sure what should be the default plot for just any old .out file
%.gp : %.out
	./scripts/plot-cols all -w l $< $(BW_ARG) -nokey
#	./scripts/plot-1stcol $< $(BW_ARG)
#	./scripts/plot-2cols 1 2 $<
