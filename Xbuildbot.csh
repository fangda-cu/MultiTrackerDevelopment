#!/bin/tcsh

source build.env

make fresh
make all

wtTestSuite --verbose --project-name=Figaro --email-report=figaro-dev@wetafx.co.nz --email-on=always ./testing/testMaya/testcase_all.py

