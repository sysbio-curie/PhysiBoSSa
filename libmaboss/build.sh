cd addons/PhysiBoSS/MaBoSS-env-2.0/engine/src

make install_lib
make MAXNODES=128 install_lib
make MAXNODES=256 install_lib
make MAXNODES=512 install_lib
make MAXNODES=1024 install_lib

make install_alib
make MAXNODES=128 install_alib
make MAXNODES=256 install_alib
make MAXNODES=512 install_alib
make MAXNODES=1024 install_alib

mkdir -p ${PREFIX}/lib
mv ../lib/libMaBoSS* ${PREFIX}/lib

mkdir -p ${PREFIX}/include
mv ../include/*.h ${PREFIX}/include

