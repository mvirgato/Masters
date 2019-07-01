#!/bin/bash

source /cvmfs/belle.cern.ch/sl6/tools/b2setup RELEASE

cd /data/pchan1/B2Kspi0/1-SignalSkim
basf2 B_2Kspi0_skim_udst.py INPUT_GLOB OUTPUT_UDST
