# List of all available tests
TESTS ?= \
	Flowthrough_trans_periodic \
	Flowthrough_x_inflow_y_outflow \
	acctest_1_maxw_500k_30kms_1deg \
	acctest_2_maxw_500k_100k_20kms_10deg \
	acctest_3_substeps \
	acctest_4_helium \
	acctest_5_proton_antiproton

MAKEFILE_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
include $(MAKEFILE_DIR)/Makefile.test

