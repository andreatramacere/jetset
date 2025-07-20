#!/bin/bash
pip install --no-build-isolation --no-deps . 
python -c "import jetset; from jetset.test_data_helper import  test_SEDs; from jetset.ebl_data import *; from jetset.Spectral_Templates_Repo import *; from jetset.jetkernel import mathkernel"
code=$?
if [ $code -ne 0 ]
then
    printf "\n\33[31mError while installing jetset.\33[0m\n"
    exit -63
fi