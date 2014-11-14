#!/bin/bash

rm -f log_test_run.txt
./generate_tests.sh
echo "Running test suite for debug flags .." | tee --append log_test_run.txt
echo " (compiling)" | tee --append log_test_run.txt
echo "" | tee --append log_test_run.txt
make >> log_test_run.txt
echo "" | tee --append log_test_run.txt
echo " (running)" | tee --append log_test_run.txt
echo "" | tee --append log_test_run.txt
./bin.debug/tests.x >> log_test_run.txt
echo "Running test suite for release flags .." | tee --append log_test_run.txt
echo " (compiling)" | tee --append log_test_run.txt
echo "" | tee --append log_test_run.txt
make CFG=release >> log_test_run.txt
echo "" | tee --append log_test_run.txt
echo " (running)" | tee --append log_test_run.txt
echo "" | tee --append log_test_run.txt
./bin.release/tests.x | tee --append log_test_run.txt
