#!/bin/bash
../../../compile.sh thinFilm.c -DBDF2 -DPRINT 
#For MPI
#-D_MPI=1
if [ $? -eq 0 ]; then
	echo "Compilation Success.";
else
	echo "Compilation error.";
	exit 1;
fi

