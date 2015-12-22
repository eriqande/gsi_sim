bin=gsi_sim-$(uname)

echo "Compiling up executable $bin"
echo


CFLAGS="-I./eca-shared/ecalibs \
	 -I./eca-shared/ranlib/src"


gsisim_SOURCES="./eca-shared/ranlib/linpack/linpack.c \
    ./eca-shared/ranlib/src/com.c \
    ./eca-shared/ranlib/src/ranlib.c \
    ./eca-shared/ecalibs/ECA_MemAlloc.c \
    ./eca-shared/ecalibs/ECA_Opt3.c \
    ./eca-shared/ecalibs/ECA_ReadData.c \
    ./eca-shared/ecalibs/ECA_utilities.c \
    ./eca-shared/ecalibs/MathStatRand.c \
    ./eca-shared/ecalibs/MCTypesEtc.c \
    ./src/gsi_sim.c"


(gcc -O3 $CFLAGS  $gsisim_SOURCES -o $bin -lm) && (echo; echo; echo "Successfully compiled the executable $bin"; echo)

