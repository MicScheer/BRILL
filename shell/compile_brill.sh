
OWD=`pwd`

if test -z $BRILL; then
   echo Shell variable BRILL not defined, trying ~/brill
   BRILL=~/brill
fi

cd $BRILL/for

# +PATCH,//COM/SHELL
# +DECK,gfortran,T=SHELL.
echo
echo
echo
echo
echo '=============================='
echo
echo
echo
echo

rm -f $BRILL/for/mshcern.o 2>/dev/null

gfortran -c -cpp -O2  \
-fd-lines-as-comments \
-fno-automatic \
-ffixed-line-length-none \
-fno-automatic \
 -funroll-loops \
mshcern.f 2>/dev/null \
&&  echo '--- mshcern.f compiled ---' \
|| echo '*** Compilation of mshcern.f failed ***'

rm -f $BRILL/bin/brill.exe 2>/dev/null

gfortran -O2 \
-ffixed-line-length-none \
-fno-automatic \
 -funroll-loops \
mshcern.o \
-o $BRILL/bin/brill.exe $BRILL/for/brill.f 2>/dev/null \
&&  echo '--- brill.f compiled and linked ---' \
|| echo '*** Compilation or linking of brill.f failed ***'

cd $OWD
