# Copyright (c) 2020-2023 by Albert Bartók and Miguel Caro

SHELL = /bin/sh

# Include user-modifiable variables from a customizable file.
# Check the makefiles/ directory for a list of tested architectures

include makefiles/Makefile.Ubuntu_gfortran_mpi

# Default locations for various files
BUILD_DIR=build
BIN_DIR=bin
INC_DIR=include
LIB_DIR=lib


# Do not change anything below this line
##########################################################

# Precompiler options can be found here
# -D_CHECK_DEALLOCATE
# -D_DEBUG

F90_OPTS += $(F90_MOD_DIR_OPT) $(INC_DIR) #-D_CHECK_ALLOCATE

PROGRAMS := turbogap


SRC := kinds.f90 \
       constants.f90 \
       error.f90 \
       printing.f90 \
       tg_memory.f90 \
       types.f90 \
       timer.f90 \
       control.f90 \
       timing.f90 \
       functions.f90 \
       splines.f90 \
       calculation.f90 \
       \
       \
       md_types.f90 \
       mc_types.f90 \
       vdw_types.f90 \
       exp_types.f90 \
       \
       soap_turbo_functions.f90 \
       soap_turbo_compress.f90 \
       soap_turbo_radial.f90 \
       soap_turbo_angular.f90 \
       soap_turbo.f90 \
       \
       neighbors_interface.f90 \
       \
       md_utils.f90 \
       \
       write_xyz.f90 \
       \
       \
       state_interface.f90 \
       resamplekin.f90 \
       md_interface.f90 \
       mc_utils.f90 \
       mc_interface.f90 \
       \
       read_utils.f90 \
       read_control.f90 \
       read_mc.f90 \
       read_md.f90 \
       read_xyz.f90 \
       read_gap.f90 \
       read_vdw.f90 \
       read_exp.f90 \
       read_files.f90 \
       \
       local_properties.f90 \
       gap.f90 \
       gap_interface.f90 \
       calculate_gap_soap.f90 \
       \
       vdw_interface.f90 \
       \
       \
       mpi_utils.f90 \
       \
       control_interface.f90 \
       \
       misc.f90 \
       \
       turbogap_main.f90


       # \
       # historylist.f90 \
       # sqnm.f90 \
       # periodic_optimizer.f90 \

       # read.f90 \
       # read_mc.f90 \
       # read_md.f90 \
       # read_control.f90 \
       # read_vdw.f90 \
       # read_exp.f90 \
       # read_stopping.f90 \
       # read_gap.f90 \
       # read_xyz.f90 \
       # read_files.f90 \

OBJ := $(addprefix $(BUILD_DIR)/,$(patsubst %.f90,%.o,$(SRC)))

PROG := $(addprefix $(BIN_DIR)/,$(PROGRAMS))

.SUFFIXES:
.SUFFIXES: .f90 .o
.PHONY: default all programs clean deepclean libturbogap test format format-check

default: libturbogap programs

all: default

clean:
	rm -rf $(OBJ) $(INC_DIR)/*.mod $(PROG)

deepclean:
	rm -rf $(BUILD_DIR) $(BIN_DIR) ${INC_DIR} ${LIB_DIR}

# Smoke-tests every scenario under tests_manual/ (serial + 2-rank MPI) against
# the just-built binary. See tests_manual/run_manual_tests.sh for what "pass"
# means here (no golden output to diff against - this is a crash/regression
# smoke test, not a numerical correctness check). Override the per-run
# timeout with `make test TEST_TIMEOUT=300`.
TEST_TIMEOUT ?= 60
test: programs
	@bash tests_manual/run_manual_tests.sh $(TEST_TIMEOUT)

# Files formatted by fprettify (see .fprettify.rc): explicitly the set of
# hand-written files actually compiled into turbogap/libturbogap (i.e. the
# real source paths behind $(SRC), plus the turbogap.f90 entry point).
# Deliberately NOT a blanket `find src -name '*.f90'`: that also picks up
# the soap_turbo submodule, vendored third_party/ code (keep matching
# upstream), src/allocation/tg_memory.f90 (generated from
# tg_memory_dims.fpp - regenerate it, don't hand-format it), and various
# unused/experimental prototype files scattered under src/ (src/gpu,
# src/allocation/test/, src/allocation/allocation.f90, src/allocation/tg_alloc.f90,
# ...) that were never wired into $(SRC) and aren't real, tested source.
# The list itself lives in .fprettify-files.txt (one path per line) so the
# pre-commit fprettify hook (scripts/run_fprettify.sh) can share the exact
# same scope instead of drifting out of sync with a second copy.
FPRETTIFY_FILES := $(shell cat .fprettify-files.txt)

format:
	fprettify -c .fprettify.rc $(FPRETTIFY_FILES)

# Note: fprettify -d always exits 0, and a fatal parse error on one file
# prints its traceback to stderr while producing no diff at all - so this
# checks fprettify's own exit code AND stderr content (parse errors),
# in addition to whether a diff was produced, rather than trusting the diff
# alone to reflect the real outcome.
format-check:
	@diff_output="$$(fprettify -c .fprettify.rc -d $(FPRETTIFY_FILES) 2>/tmp/fprettify-check-err.$$$$)"; \
	status=$$?; \
	err_output="$$(cat /tmp/fprettify-check-err.$$$$)"; \
	rm -f /tmp/fprettify-check-err.$$$$; \
	ok=1; \
	if [ $$status -ne 0 ]; then \
		echo "fprettify exited with status $$status"; \
		ok=0; \
	fi; \
	if echo "$$err_output" | grep -qi "error\|traceback\|exception"; then \
		echo "$$err_output"; \
		ok=0; \
	fi; \
	if [ -n "$$diff_output" ]; then \
		echo "$$diff_output"; \
		ok=0; \
	fi; \
	if [ "$$ok" -eq 1 ]; then \
		echo "All files formatted correctly."; \
	else \
		echo; \
		echo "Formatting issues found above - run 'make format' to fix (or investigate the fprettify error)."; \
		exit 1; \
	fi

.SECONDEXPANSION:
.SECONDARY: $(OBJS)

programs: $(PROG)

libturbogap: $(OBJ) ${LIB_DIR}
	ar scr $(LIB_DIR)/libturbogap.a $(OBJ)

$(BIN_DIR)/%: src/%.f90 $(OBJ) | $$(@D)
	$(F90) $(PP) $(F90_OPTS) $< -o $@ $(OBJ) $(LIBS)

$(BUILD_DIR)/%.o: src/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@

$(BUILD_DIR)/%.o: src/read/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@

$(BUILD_DIR)/%.o: src/vdw/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@

$(BUILD_DIR)/%.o: src/gap/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@

$(BUILD_DIR)/%.o: src/xyz/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@

$(BUILD_DIR)/%.o: src/neighbors/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@

$(BUILD_DIR)/%.o: src/types/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@

$(BUILD_DIR)/%.o: src/allocation/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@

$(BUILD_DIR)/%.o: src/utils/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@

$(BUILD_DIR)/%.o: src/control/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@

$(BUILD_DIR)/%.o: src/third_party/bussi_thermostat/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@

$(BUILD_DIR)/%.o: src/md/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@

$(BUILD_DIR)/%.o: src/mc/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@

$(BUILD_DIR)/%.o: src/mpi/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@


$(BUILD_DIR)/%.o: src/madlib/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@

# stable quasi newton optimizer
# $(BUILD_DIR)/%.o: src/third_party/vc-sqnm/src/fortran/%.f90 | $$(@D)
# 	$(F90) $(PP) $(F90_OPTS) -c $< -o $@


$(BUILD_DIR)/%.o: src/soap_turbo/src/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@

$(BUILD_DIR): ${INC_DIR}
	mkdir -p $@

$(BIN_DIR):
	mkdir -p $@

$(INC_DIR):
	mkdir -p $@

$(LIB_DIR):
	mkdir -p $@
