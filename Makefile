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

F90_OPTS += $(F90_MOD_DIR_OPT) $(INC_DIR)

PROGRAMS := turbogap


SRC := kinds.f90 \
       error.f90 \
       printing.f90 \
       types.f90 \
       timer.f90 \
       control.f90 \
       timing.f90 \
       functions.f90 \
       splines.f90 \
       calculation.f90 \
       misc.f90 \
       \
       md_types.f90 \
       mc_types.f90 \
       vdw_types.f90 \
       \
       soap_turbo_functions.f90 \
       soap_turbo_compress.f90 \
       soap_turbo_radial.f90 \
       soap_turbo_angular.f90 \
       soap_turbo.f90 \
       \
       \
       neighbors_interface.f90 \
       md_utils.f90 \
       md_interface.f90 \
       mc_interface.f90 \
       state_interface.f90 \
       \
       read_utils.f90 \
       read_control.f90 \
       read_mc.f90 \
       read_md.f90 \
       read_xyz.f90 \
       read_gap.f90 \
       read_vdw.f90 \
       read_files.f90 \
       \
       mpi_utils.f90 \
       \
       control_interface.f90 \
       \
       turbogap_main.f90



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
.PHONY: default all programs clean deepclean libturbogap

default: libturbogap programs

all: default

clean:
	rm -rf $(OBJ) $(INC_DIR)/*.mod $(PROG)

deepclean:
	rm -rf $(BUILD_DIR) $(BIN_DIR) ${INC_DIR} ${LIB_DIR}

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

$(BUILD_DIR)/%.o: src/utils/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@

$(BUILD_DIR)/%.o: src/control/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@

$(BUILD_DIR)/%.o: src/md/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@

$(BUILD_DIR)/%.o: src/mc/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@

$(BUILD_DIR)/%.o: src/mpi/%.f90 | $$(@D)
	$(F90) $(PP) $(F90_OPTS) -c $< -o $@


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
