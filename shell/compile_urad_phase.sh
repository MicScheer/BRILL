
OWD=`pwd`

if test -z $BRILL; then
   echo Shell variable BRILL not defined, trying `pwd`
   BRILL=`pwd`
fi

cd $BRILL/for

rm -f $BRILL/urad_phase.exe

gfortran -c -O2 -cpp \
-ffpe-summary=invalid,zero,overflow \
-fopenmp \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none \
-finit-local-zero -funroll-loops \
urad_modules.f

gfortran -c -O2 -cpp \
-ffpe-summary=invalid,zero,overflow \
-fopenmp \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none \
-finit-local-zero -funroll-loops \
urad_util.f

gfortran -O2 -cpp \
-ffpe-summary=invalid,zero,overflow \
-fopenmp \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none \
-finit-local-zero -funroll-loops \
urad_modules.o urad_util.o \
-o $BRILL/bin/urad_phase.exe \
urad_phase_main.f

cd $OWD
