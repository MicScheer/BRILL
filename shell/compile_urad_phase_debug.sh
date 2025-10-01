
OWD=`pwd`

if test -z $BRILL; then
   echo Shell variable BRILL not defined, trying `pwd`
   BRILL=`pwd`
fi

cd $BRILL/for

rm -f $BRILL/urad_phase.exe

# mshcern.f is created by cat $WAVE_INCL/mshcern/*.f > $BRILL/for/mshcern.f

gfortran -c -g -cpp -w \
-ffpe-summary=invalid,zero,overflow \
-fopenmp \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none \
-finit-local-zero -funroll-loops \
mshcern.f

gfortran -c -g -cpp \
-ffpe-summary=invalid,zero,overflow \
-fopenmp \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none \
-finit-local-zero -funroll-loops \
urad_modules.f

gfortran -c -g -cpp \
-ffpe-summary=invalid,zero,overflow \
-fopenmp \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none \
-finit-local-zero -funroll-loops \
urad_util.f

gfortran -g -cpp \
-ffpe-summary=invalid,zero,overflow \
-fopenmp \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none \
-finit-local-zero -funroll-loops \
urad_modules.o urad_util.o \
mshcern.o \
-o $BRILL/bin/urad_phase_debug.exe \
urad_phase_main_debug.f

cd $OWD
