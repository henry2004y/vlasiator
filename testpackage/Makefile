# List of all available tests
TESTS ?= \
	advection1d \
	advection2d \
	acctest_1_maxw_500k_30kms_1deg

MAKEFILE_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
include $(MAKEFILE_DIR)/Makefile.test

