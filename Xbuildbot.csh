#!/bin/tcsh

source build.env

make fresh || exit 1
make all -j 8 || exit 1

wtTestSuite --verbose --project-name=Figaro --email-report=figaro-dev@wetafx.co.nz --email-on=always ./testing/testPlayblasts/direct_playblast_optversion_buildbot.py || exit 1



