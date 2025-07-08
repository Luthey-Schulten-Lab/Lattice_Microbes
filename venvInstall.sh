#!/bin/bash

if [ x$VIRTUAL_ENV == x ]; then \
   echo $0: Not in a venv
   exit 1
fi

if [ x$1 == x ]; then \
   echo usage: $0 BUILD-DIR
   exit -1
fi

buildDir=$1

cd $buildDir

cp -r  bin/{lm,lm_python,lm_setdm,lm_setp,lm_setrm} $VIRTUAL_ENV/bin
cp -r python/{pyLM,pySTDLM} lib/{lm.py,_lm.so} $(python -c 'import site; print(site.getsitepackages()[0])')

