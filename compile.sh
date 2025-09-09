# Coloring codes for visuals --- optional 
Color_Off='\033[0m'       # Text Reset
Red='\033[0;31m'          # Red
BGreen='\033[1;32m'       # Green
BRed='\033[1;31m'         # Red

if [[ -z "$1" ]]; then
	echo -e "${Red}Please provide at least one file to compile${Color_Off}";
	echo -e "${BGreen}Usage:${Color_Off} compile.sh + (.c file to compile) + -DFLAGS";
	exit 2
elif [[ "$1" != *"."* ]]; then
	echo -e "${Red}Missing compilation file${Color_Off}";
	echo -e "${BGreen}Usage:${Color_Off} compile.sh + (.c file to compile) + -DFLAGS";
	exit 3
fi

echo -e "${BGreen}Compiling $1 ${Color_Off}";

s=${@:2}
echo -e "${BRed}Additional Flags: $s ${Color_Off}";

if [[ "$*" == *"-D_MPI=1"* ]]
then
	PREFIX=" CC99='mpicc -std=c99' "
	PARFLAGS=" "
else
	PARFLAGS=" -fopenmp "
	PREFIX=" "
fi

GLLIBS=" -L${BASILISK}/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11 "

com=" ${PREFIX} qcc -O2 ${s} -I${PWD}/utils -Wall $PARFLAGS $1 -o run.exe $GLLIBS -lm"

#output compilation command
echo $com
#run the compilation command
eval $com


if [ $? -eq 0 ]; then
	echo -e "${BGreen}Compilation Success ${Color_Off}";
else
	echo -e "${BRed}Compilation Failed ${Color_Off}";
fi
