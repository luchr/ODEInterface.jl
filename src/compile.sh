#!/bin/bash

FLAGS_I64=(  -fdefault-integer-8  -fdefault-real-8 -fdefault-double-8 )
FLAGS_I32=(  -fdefault-real-8 -fdefault-double-8 )

echo "dopri5"
gfortran -c -fPIC "${FLAGS_I64[@]}"  -o dopri5.o   dopri5.f
gcc  -shared -fPIC -Wl,-soname,libdorpi5.so -lgfortran -o dopri5.so  dopri5.o
rm ./dopri5.o

echo "dopri5_i32"
gfortran -c -fPIC "${FLAGS_I32[@]}"  -o dopri5_i32.o   dopri5.f
gcc  -shared -fPIC -Wl,-soname,libdorpi5_i32.so -lgfortran -o dopri5_i32.so  dopri5_i32.o
rm ./dopri5_i32.o

echo "dopri853"
gfortran -c -fPIC  "${FLAGS_I64[@]}" -o dop853.o   dop853.f
gcc  -shared -fPIC -Wl,-soname,libdop853.so -lgfortran -o dop853.so  dop853.o
rm ./dop853.o

echo "dopri853_i32"
gfortran -c -fPIC  "${FLAGS_I32[@]}" -o dop853_i32.o   dop853.f
gcc  -shared -fPIC -Wl,-soname,libdop853_i32.so -lgfortran -o dop853_i32.so  dop853_i32.o
rm ./dop853_i32.o

echo "odex"
gfortran -c -fPIC  "${FLAGS_I64[@]}" -o odex.o   odex.f
gcc  -shared -fPIC -Wl,-soname,libodex.so -lgfortran -o odex.so  odex.o
rm ./odex.o

echo "odex_i32"
gfortran -c -fPIC  "${FLAGS_I32[@]}" -o odex_i32.o   odex.f
gcc  -shared -fPIC -Wl,-soname,libodex_i32.so -lgfortran -o odex_i32.so  odex_i32.o
rm ./odex_i32.o

echo "dc_lapack, lapack, lapackc"
gfortran -c -fPIC "${FLAGS_I64[@]}"  -o dc_lapack.o   dc_lapack.f
gfortran -c -fPIC "${FLAGS_I64[@]}"  -o lapack.o      lapack.f
gfortran -c -fPIC "${FLAGS_I64[@]}"  -o lapackc.o      lapackc.f

echo "dc_lapack_i32, lapack_i32, lapackc_i32"
gfortran -c -fPIC  "${FLAGS_I32[@]}" -o dc_lapack_i32.o   dc_lapack.f
gfortran -c -fPIC  "${FLAGS_I32[@]}" -o lapack_i32.o      lapack.f
gfortran -c -fPIC  "${FLAGS_I32[@]}" -o lapackc_i32.o      lapackc.f

echo "radau5"
gfortran -c -fPIC "${FLAGS_I64[@]}"  -o radau5.o   radau5.f
gcc -shared -fPIC -Wl,-soname,libradau5.so  -lgfortran  -o radau5.so \
              radau5.o dc_lapack.o lapack.o lapackc.o
rm ./radau5.o

echo "radau5_i32"
gfortran -c -fPIC "${FLAGS_I32[@]}"  -o radau5_i32.o   radau5.f
gcc -shared -fPIC -Wl,-soname,libradau5_i32.so  -lgfortran  -o radau5_i32.so \
              radau5_i32.o dc_lapack_i32.o lapack_i32.o lapackc_i32.o
rm ./radau5_i32.o

echo "radau"
gfortran -c -fPIC "${FLAGS_I64[@]}"  -o radau.o   radau.f
gcc -shared -fPIC -Wl,-soname,libradau5.so  -lgfortran  -o radau.so \
              radau.o dc_lapack.o lapack.o lapackc.o
rm ./radau.o

echo "radau_i32"
gfortran -c -fPIC "${FLAGS_I32[@]}"  -o radau_i32.o   radau.f
gcc -shared -fPIC -Wl,-soname,libradau5_i32.so  -lgfortran  -o radau_i32.so \
              radau_i32.o dc_lapack_i32.o lapack_i32.o lapackc_i32.o
rm ./radau_i32.o

rm ./dc_lapack.o ./lapack.o ./lapackc.o
rm ./dc_lapack_i32.o ./lapack_i32.o ./lapackc_i32.o
