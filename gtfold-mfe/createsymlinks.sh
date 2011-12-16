#!/bin/sh

BINDIR=/home/users/pgaurav/git-hub/gtfold_final/bin
echo 'creating symlinks under' $BINDIR
ln -sf gtfold $BINDIR/gtmfe
ln -sf gtfold $BINDIR/gtsubopt
ln -sf gtfold $BINDIR/gtboltzmann

