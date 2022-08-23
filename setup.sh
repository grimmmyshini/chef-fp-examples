#!/bin/bash

set -a

ws_dir=$(dirname "$(readlink -f "$0")")
clad_bm_dir=$ws_dir/clad-fp-error-est-examples
pr_mod_dir=$clad_bm_dir/PrintModel

cplusincludepath_add() {
    if [ -d "$1" ] && [[ ":$CPLUS_INCLUDE_PATH:" != *":$1:"* ]]; then
        CPLUS_INCLUDE_PATH="${CPLUS_INCLUDE_PATH:+"$CPLUS_INCLUDE_PATH:"}$1"
    fi
}

ldlibrarypath_add() {
    if [ -d "$1" ] && [[ ":$LD_LIBRARY_PATH:" != *":$1:"* ]]; then
        LD_LIBRARY_PATH="${LD_LIBRARY_PATH:+"$LD_LIBRARY_PATH:"}$1"
    fi
}

librarypath_add() {
    if [ -d "$1" ] && [[ ":$LIBRARY_PATH:" != *":$1:"* ]]; then
        LIBRARY_PATH="${LIBRARY_PATH:+"$LIBRARY_PATH:"}$1"
    fi
}

path_add() {
    if [ -d "$1" ] && [[ ":$PATH:" != *":$1:"* ]]; then
        PATH="${PATH:+"$PATH:"}$1"
    fi
}

cd $ws_dir

cplusincludepath_add $ws_dir/CoDiPack-1.9/include
cplusincludepath_add $ws_dir/benchmark/include
cplusincludepath_add $ws_dir/adapt-fp
cplusincludepath_add $ws_dir/clad/include

ldlibrarypath_add $ws_dir/clang/lib
ldlibrarypath_add $ws_dir/clad/build/lib

librarypath_add $ws_dir/benchmark/build/src

path_add $ws_dir/clang/bin

CODIPACK_HOME=$ws_dir/CoDiPack-1.9

cladbench() {
    $ws_dir/clang/bin/clang++ \
        -Xclang -add-plugin -Xclang clad -Xclang -load -Xclang clad.so \
        -lstdc++ -lm \
        -O2 \
        "$@" \
        -L$ws_dir/benchmark/build/src -lbenchmark -lpthread
}

adaptbench() {
    $ws_dir/clang/bin/clang++ \
        -std=c++11 -O2 \
        "$@" \
        -L$ws_dir/benchmark/build/src -lbenchmark -lpthread
}

benchall() {
    $ws_dir/clang/bin/clang++ \
        -Xclang -add-plugin -Xclang clad -Xclang -load -Xclang clad.so \
        -Xclang -plugin-arg-clad -Xclang -fcustom-estimation-model \
        -Xclang -plugin-arg-clad -Xclang $pr_mod_dir/libPrintModel.so \
        -lstdc++ -lm -std=c++11 \
        "$@" \
        -L$ws_dir/benchmark/build/src \
        -lpthread -lbenchmark -DCODI_ZeroAdjointReverse=0
}

cladgen() {
    $ws_dir/clang/bin/clang++ \
        -Xclang -add-plugin -Xclang clad -Xclang -load -Xclang clad.so \
        -Xclang -plugin-arg-clad -Xclang -fcustom-estimation-model \
        -Xclang -plugin-arg-clad -Xclang $pr_mod_dir/libPrintModel.so \
        -lstdc++ -lm \
        "$@" \
        -Xclang -plugin-arg-clad -Xclang -fgenerate-source-file
}

claddast() {
    $ws_dir/clang/bin/clang++ \
        -Xclang -add-plugin -Xclang clad -Xclang -load -Xclang clad.so \
        -Xclang -plugin-arg-clad -Xclang -fcustom-estimation-model \
        -Xclang -plugin-arg-clad -Xclang $pr_mod_dir/libPrintModel.so \
        -lstdc++ -lm \
        "$@" \
        -Xclang -plugin-arg-clad -Xclang -fdump-derived-fn-ast
}

claddfn() {
    $ws_dir/clang/bin/clang++ \
        -Xclang -add-plugin -Xclang clad -Xclang -load -Xclang clad.so \
        -Xclang -plugin-arg-clad -Xclang -fcustom-estimation-model \
        -Xclang -plugin-arg-clad -Xclang $pr_mod_dir/libPrintModel.so \
        -lstdc++ -lm -v \
        "$@" \
        -Xclang -plugin-arg-clad -Xclang -fdump-derived-fn
}

prmodb() {
    WD=$(pwd)
    cd $pr_mod_dir
    # clang -I$ws_dir/clang/include                                          \
    #       -I$ws_dir/llvm-project/llvm/include                              \
    #       -I$ws_dir/llvm-project/clang/include                             \
    #       -I$ws_dir/clang/tools/clang/include                              \
    #       -fPIC -shared -fno-rtti -Wl,-undefined -Wl,suppress              \
    #       -lm -lstdc++                                                     \
    #       -v ./PrintModel.cpp -o $pr_mod_dir/libPrintModel.so
    $ws_dir/clang/bin/clang++ -v \
        -D_GNU_SOURCE -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS \
        -D__STDC_LIMIT_MACROS -DcladCustomModelPlugin_EXPORTS \
        -I$ws_dir/clad/build/demos/ErrorEstimation/CustomModel \
        -I$ws_dir/clad/demos/ErrorEstimation/CustomModel \
        -I$ws_dir/clad/include \
        -I$ws_dir/llvm-project/clang/include \
        -I$ws_dir/clang/tools/clang/include \
        -I$ws_dir/llvm-project/llvm/include \
        -I$ws_dir/clang/include \
        -I$ws_dir/llvm-project/llvm/tools/clang/include \
        -fPIC -Werror=date-time -std=c++11 \
        -w -fno-common -Woverloaded-virtual -Wcast-qual -fno-strict-aliasing \
        -pedantic -Wno-long-long -Wall -W -Wno-unused-parameter -Wwrite-strings \
        -g -fPIC -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS \
        -D__STDC_LIMIT_MACROS -fno-rtti -MD -MT PrintModel.cpp.o \
        -MF PrintModel.cpp.o.d -o PrintModel.cpp.o -c PrintModel.cpp
    $ws_dir/clang/bin/clang++ -fPIC -fPIC \
        -Werror=date-time -std=c++11 -w -fno-common -Woverloaded-virtual \
        -Wcast-qual -fno-strict-aliasing -pedantic -Wno-long-long -Wall -W \
        -Wno-unused-parameter -Wwrite-strings -g -Wl,-z,defs -Wl,-z,nodelete \
        -shared -Wl,-soname,libPrintModel.so -o libPrintModel.so PrintModel.cpp.o \
        -Wl,--unresolved-symbols=ignore-in-object-files
    cd $WD
}

cladb() {
    WD=$(pwd)
    cd $ws_dir/clad/build
    make -j6
    cd $WD
}

edconf() {
    nvim $ws_dir/conf/setup.sh
}

srcconf() {
    WD=$(pwd)
    srcfpbench
    cd $WD
}

set +a

echo "Loaded configuration!"

setup() {
    # INSTALL LLVM and Clang
    git clone https://github.com/llvm/llvm-project.git --depth=1 --single-branch --branch=llvmorg-9.0.1
    mkdir clang && cd clang
    cmake -DLLVM_ENABLE_PROJECTS=clang -DCMAKE_BUILD_TYPE=Release -G "Unix Makefiles" ../llvm-project/llvm/
    make -j6
    cd ..

    # INSTALL Google Benchmark
    git clone https://github.com/google/benchmark
    cd benchmark/
    cmake -E make_directory "build"
    cmake -E chdir "build" cmake -DBENCHMARK_DOWNLOAD_DEPENDENCIES=on -DCMAKE_BUILD_TYPE=Release ../
    cmake --build "build" --config Release
    cd ..

    # INSTALL ADAPT-FP and CoDiPack
    git clone https://github.com/LLNL/adapt-fp
    cd adapt-fp
    git checkout be4cda8
    cd ..
    wget https://github.com/SciCompKL/CoDiPack/archive/refs/tags/v1.9.1.tar.gz
    tar -xzvf v1.9.1.tar.gz
    rm v1.9.1.tar.gz
    mv CoDiPack-1.9.1/ CoDiPack-1.9/

    # INSTALL Clad
    pip install lit
    git clone https://github.com/sudo-panda/clad --depth=1 --single-branch --branch=ee-bench
    cd clad
    mkdir build && cd build
    cmake .. -DCMAKE_BUILD_TYPE=Debug -DLLVM_DIR=$ws_dir/clang/lib/cmake/llvm/ -DClang_DIR=$ws_dir/clang/lib/cmake/clang -DCMAKE_INSTALL_PREFIX=../inst -DLLVM_EXTERNAL_LIT="$(which lit)"
    make -j6 && make install -j6
    cd ../..

    # INSTALL benchmarks
    git clone https://github.com/grimmmyshini/clad-fp-error-est-examples
    cd clad-fp-error-est-examples/PrintModel
    clang -I$ws_dir/clad/inst/include -I$ws_dir/clang/tools/clang/include -I$ws_dir/llvm-project/clang/include -I$ws_dir/clang/include -I$ws_dir/llvm-project/llvm/include -fPIC -shared -fno-rtti -Wl,-undefined -Wl,suppress PrintModel.cpp -o libPrintModel.so
    cd ../blackscholes/data
    $ws_dir/clang/bin/clang++ inputGen.c -o inputGen.out
    inputGen.out 100 100.txt
    inputGen.out 1000 1000.txt
    inputGen.out 10000 10000.txt
    inputGen.out 100000 100000.txt
    inputGen.out 1000000 1000000.txt
    cd ../../kmeans/data/inpuGen/
    $ws_dir/clang/bin/clang++ datagen.cpp -o datagen.out
    datagen.out 100 ../100.txt -f
    datagen.out 1000 ../1000.txt -f
    datagen.out 10000 ../10000.txt -f
    datagen.out 100000 ../100000.txt -f
    datagen.out 1000000 ../1000000.txt -f
    cd $ws_dir
}

setup
cd clad-fp-error-est-examples
./runner.sh